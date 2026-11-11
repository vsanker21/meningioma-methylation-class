"""Separate a missed clonal arm, a shifted baseline, and an RNA fit residual.

The contrasts are fixed in docs/protocol_locks.md. Stored purity and stored RNA
fractions are not replaced.
"""

import csv
import importlib.util
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.optimize import nnls

ROOT = Path(r"F:\Meningioma Genetics")
WORK = ROOT / "data" / "work"
N_BOOT = 2000
SEED = 1
ARM_CUT = 0.08


def read_csv(path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def as_float(value):
    if value is None or value in ("", "NA", "NaN"):
        return np.nan
    return float(value)


def bootstrap_median(values):
    values = np.asarray(list(values), dtype=float)
    values = values[np.isfinite(values)]
    n = int(values.size)
    if n == 0:
        return n, np.nan, np.nan, np.nan
    median = float(np.median(values))
    if n < 10:
        return n, median, np.nan, np.nan
    rng = np.random.default_rng(SEED)
    draws = values[rng.integers(0, n, size=(N_BOOT, n))]
    lo, hi = np.quantile(np.median(draws, axis=1), [0.025, 0.975])
    return n, median, float(lo), float(hi)


def parse_loss_arms(text):
    if text is None or text in ("", "NA"):
        return []
    levels = []
    for piece in text.split(";"):
        arm, _, level = piece.rpartition(":")
        if arm:
            levels.append(float(level))
    return levels


def add(summary, contrast, cohort, values, extra=None):
    n, median, lo, hi = bootstrap_median(values)
    row = {
        "contrast": contrast,
        "cohort": cohort,
        "n": n,
        "median": median,
        "ci_low": lo,
        "ci_high": hi,
        "fraction_positive": float(np.mean(np.asarray(values)[np.isfinite(values)] > 0)) if n else np.nan,
    }
    if extra:
        row.update(extra)
    summary.append(row)


def load_rna_module():
    path = ROOT / "scripts" / "rna_deconvolution.py"
    spec = importlib.util.spec_from_file_location("rna_deconvolution", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def rna_diagnostics(stored):
    rna = load_rna_module()
    genes, types, signature = rna.load_signature()
    neo = types.index("neoplastic")
    lab_to_sample = rna.bayley_lab_to_sample()
    symbol_of = rna.ensembl_to_symbol()
    jobs = (
        ("GSE183647", rna.TPM, None, False, None),
        ("GSE189521", rna.COUNTS, symbol_of, True, lab_to_sample),
    )
    rows = []
    worst = 0.0
    for cohort, path, symbols, as_cpm, rename in jobs:
        samples, present, values, library = rna.stream_markers(path, genes, symbols)
        index = {gene: i for i, gene in enumerate(present)}
        keep = [gene for gene in genes if gene in index]
        observed = values[[index[gene] for gene in keep], :]
        if as_cpm:
            observed = observed / library * 1e6
        design = signature[[genes.index(gene) for gene in keep], :]
        for i, sample in enumerate(samples):
            sample_id = rename[sample] if rename else sample
            if rename and sample not in rename:
                continue
            coef, rnorm = nnls(design, observed[:, i])
            total = float(coef.sum())
            fraction = float(coef[neo] / total) if total > 0 else np.nan
            previous = stored.get((cohort, sample_id))
            if previous is None or not np.isfinite(previous):
                continue
            worst = max(worst, abs(fraction - previous))
            norm = float(np.linalg.norm(observed[:, i]))
            rows.append({
                "sample_id": sample_id,
                "cohort": cohort,
                "coef_sum": total,
                "relative_residual": rnorm / norm if norm > 0 else np.nan,
                "neoplastic_recomputed": fraction,
            })
    if worst > 1e-4:
        raise SystemExit("recomputed RNA fraction disagrees with the stored fraction by %s" % worst)
    if len(rows) != len(stored):
        raise SystemExit("RNA diagnostic matched %d tumours; the stored file has %d" % (len(rows), len(stored)))
    print("RNA fraction match, max abs difference", worst, "n", len(rows))
    return rows


def main():
    paired = read_csv(WORK / "purity_scale_paired.csv")
    purity = {row["sample_id"]: row for row in read_csv(WORK / "purity_estimates.csv")}
    summary = []

    for cohort in ("GSE183647", "GSE189521"):
        rows = [row for row in paired if row["cohort"] == cohort]
        deepest_gap = [as_float(row["deepest_L"]) - as_float(row["L_rna"]) for row in rows]
        add(summary, "deepest_qualifying_L_minus_L_rna", cohort, deepest_gap)
        spreads = []
        multi = []
        one_arm = 0
        for row in rows:
            levels = parse_loss_arms(purity[row["sample_id"]]["loss_arms"])
            if not levels:
                continue
            spread = max(levels) - min(levels)
            spreads.append(spread)
            if len(levels) == 1:
                one_arm += 1
            else:
                multi.append(spread)
            parsed_min = min(levels)
            if abs(parsed_min - as_float(row["deepest_L"])) > 0.0015:
                raise SystemExit("%s loss_arms minimum disagrees with deepest_L" % row["sample_id"])
        add(summary, "qualifying_arm_spread", cohort, spreads, extra={"n_single_arm": one_arm})
        add(summary, "qualifying_arm_spread_two_or_more", cohort, multi)

    arm_values = defaultdict(list)
    with (WORK / "purity_estimates.arms.csv").open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            value = as_float(row["median_ratio"])
            if np.isfinite(value):
                arm_values[row["sample_id"]].append(value)

    for cohort in ("GSE183647", "GSE189521"):
        gaps = []
        for row in paired:
            if row["cohort"] != cohort:
                continue
            deepest_median = min(arm_values[row["sample_id"]])
            gaps.append(deepest_median - as_float(row["L_rna"]))
        add(summary, "most_negative_arm_median_minus_L_rna", cohort, gaps)

    for cohort in ("GSE183647", "GSE189521", "GSE180061"):
        baseline = []
        loss_side = []
        gain_side = []
        for sample_id, row in purity.items():
            if row["cohort"] != cohort:
                continue
            simple = as_float(row["simple_l_22q"])
            segment = as_float(row["l_22q"])
            if np.isfinite(simple) and np.isfinite(segment):
                baseline.append(simple - segment)
            values = np.asarray(arm_values[sample_id], dtype=float)
            if values.size:
                loss_side.append(float(np.mean(values < -ARM_CUT)))
                gain_side.append(float(np.mean(values > ARM_CUT)))
        add(summary, "simple_l_22q_minus_segment_l_22q", cohort, baseline)
        add(summary, "fraction_arms_below_minus_0.08", cohort, loss_side)
        summary[-1]["fraction_positive"] = np.nan
        add(summary, "fraction_arms_above_plus_0.08", cohort, gain_side)
        summary[-1]["fraction_positive"] = np.nan

    for cohort in ("GSE183647", "GSE189521"):
        gaps = []
        for row in paired:
            if row["cohort"] != cohort:
                continue
            simple = as_float(purity[row["sample_id"]]["simple_l_22q"])
            segment = as_float(purity[row["sample_id"]]["l_22q"])
            if np.isfinite(simple) and np.isfinite(segment):
                gaps.append(simple - as_float(row["L_rna"]))
        add(summary, "simple_l_22q_minus_L_rna", cohort, gaps)

    stored_fraction = {
        (row["cohort"], row["sample_id"]): as_float(row["neoplastic"])
        for row in read_csv(WORK / "rna_neoplastic_fraction.csv")
    }
    diagnostics = rna_diagnostics(stored_fraction)
    diag_path = WORK / "rna_fit_diagnostics.csv"
    with diag_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(diagnostics[0]))
        writer.writeheader()
        writer.writerows(diagnostics)
    for cohort in ("GSE183647", "GSE189521"):
        part = [row for row in diagnostics if row["cohort"] == cohort]
        add(summary, "rna_coefficient_sum", cohort, [row["coef_sum"] for row in part])
        summary[-1]["fraction_positive"] = np.nan
        add(summary, "rna_relative_residual", cohort, [row["relative_residual"] for row in part])
        summary[-1]["fraction_positive"] = np.nan

    out = WORK / "purity_scale_mechanism.csv"
    fieldnames = ["contrast", "cohort", "n", "median", "ci_low", "ci_high", "fraction_positive", "n_single_arm"]
    with out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(summary)
    print("wrote", out)
    for row in summary:
        median = row["median"]
        lo = row["ci_low"]
        hi = row["ci_high"]
        print(
            row["contrast"], row["cohort"], "n", row["n"],
            "median", None if not np.isfinite(median) else round(median, 3),
            "ci", None if not np.isfinite(lo) else (round(lo, 3), round(hi, 3)),
            "frac>0", None if not np.isfinite(row["fraction_positive"]) else round(row["fraction_positive"], 3),
            "single", row.get("n_single_arm", ""),
        )


if __name__ == "__main__":
    main()
