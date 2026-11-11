"""Compare stored arm levels with the clonal level implied by the RNA fraction.

Does not refit intensity purity or the RNA deconvolution. The paired difference
and the bootstrap are the ones fixed in docs/protocol_locks.md.
"""

import csv
from pathlib import Path

import numpy as np

ROOT = Path(r"F:\Meningioma Genetics")
WORK = ROOT / "data" / "work"
CLINICAL = ROOT / "data" / "raw" / "clinical"
N_BOOT = 2000
SEED = 1
TRUNCATION = {"Frameshift", "Nonsense", "Splice"}


def read_csv(path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def as_float(value):
    if value is None or value in ("", "NA", "NaN"):
        return np.nan
    return float(value)


def bootstrap_median(values):
    values = np.asarray(values, dtype=float)
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


def nf2_stratum(cohort, choudhury, bayley_mutation):
    if cohort == "GSE183647":
        if choudhury == "1":
            return "nf2_loss"
        if choudhury == "0":
            return "nf2_not_loss"
        return "nf2_missing"
    if cohort == "GSE189521":
        if bayley_mutation in TRUNCATION:
            return "nf2_loss"
        if bayley_mutation == "In frame deletion":
            return "nf2_inframe_excluded"
        if bayley_mutation == "":
            return "nf2_not_loss"
        return "nf2_missing"
    return "nf2_not_applicable"


def main():
    purity = read_csv(WORK / "purity_estimates.csv")
    rna = {row["sample_id"]: row for row in read_csv(WORK / "rna_neoplastic_fraction.csv")}
    choudhury = {
        row["sample_id"]: row["nf2_loss"]
        for row in read_csv(CLINICAL / "choudhury_2022_natgenet" / "sample_clinical_key.csv")
    }
    bayley = {
        row["Sample"]: row["NF2"]
        for row in read_csv(CLINICAL / "bayley_2022_sciadv" / "sample_clinical_key.csv")
    }

    formula_error = []
    paired = []
    unestimable = []
    ruler_gaps = []
    for row in purity:
        status = row["purity_status"]
        event_l = as_float(row["event_L"])
        stored = as_float(row["purity"])
        if status == "estimated":
            recomputed = 2 * (1 - 2 ** event_l)
            formula_error.append(abs(recomputed - stored))
        rna_row = rna.get(row["sample_id"])
        neoplastic = as_float(rna_row["neoplastic"]) if rna_row else np.nan
        stratum = nf2_stratum(row["cohort"], choudhury.get(row["sample_id"], ""), bayley.get(row["sample_id"], ""))
        if status == "estimated" and np.isfinite(event_l) and np.isfinite(neoplastic):
            if not 0 <= neoplastic <= 1:
                raise SystemExit("%s RNA neoplastic fraction %s is outside [0, 1]" % (row["sample_id"], neoplastic))
            implied = float(np.log2(1 - neoplastic / 2))
            paired.append({
                "sample_id": row["sample_id"],
                "cohort": row["cohort"],
                "group": row["group"],
                "event_arm": row["event_arm"],
                "event_L": event_l,
                "event_coverage": as_float(row["event_coverage"]),
                "deepest_arm": row["deepest_arm"],
                "deepest_L": as_float(row["deepest_L"]),
                "intensity_purity": stored,
                "rna_neoplastic": neoplastic,
                "L_rna": implied,
                "delta_L": event_l - implied,
                "purity_gap": neoplastic - stored,
                "ruler_is_deepest": "yes" if row["event_arm"] == row["deepest_arm"] else "no",
                "nf2_stratum": stratum,
            })
        elif rna_row and status != "estimated":
            unestimable.append({
                "cohort": row["cohort"],
                "group": row["group"],
                "rna_neoplastic": neoplastic,
                "purity_status": status,
            })
        if status == "estimated" and row["event_arm"] == "22q" and row["deepest_arm"] not in ("", "22q"):
            deepest = as_float(row["deepest_L"])
            if np.isfinite(deepest):
                ruler_gaps.append({
                    "cohort": row["cohort"],
                    "event_L": event_l,
                    "deepest_L": deepest,
                    "difference": event_l - deepest,
                })

    errors = np.asarray(formula_error, dtype=float)
    if errors.size and (not np.all(np.isfinite(errors)) or np.max(errors) > 1e-8):
        raise SystemExit("stored purity does not match 2*(1-2^L); max abs error %s" % np.max(errors))
    print("formula check max abs error", 0 if errors.size == 0 else float(np.max(errors)))
    print("paired", len(paired), "unestimable with RNA", len(unestimable), "22q shallower-arm contrasts", len(ruler_gaps))

    paired_path = WORK / "purity_scale_paired.csv"
    with paired_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(paired[0]))
        writer.writeheader()
        writer.writerows(paired)

    summary = []

    def add(contrast, cohort, stratum_type, stratum, observed, reference, difference, extra=None):
        n, median, lo, hi = bootstrap_median(difference)
        _, med_obs, _, _ = bootstrap_median(observed)
        _, med_ref, _, _ = bootstrap_median(reference)
        difference = np.asarray(difference, dtype=float)
        row = {
            "contrast": contrast,
            "cohort": cohort,
            "stratum_type": stratum_type,
            "stratum": stratum,
            "n": n,
            "median_observed": med_obs,
            "median_reference": med_ref,
            "median_difference": median,
            "ci_low": lo,
            "ci_high": hi,
            "fraction_positive": float(np.mean(difference > 0)) if n else np.nan,
        }
        if extra:
            row.update(extra)
        summary.append(row)

    def subset_rows(rows, cohort=None, **equals):
        chosen = []
        for row in rows:
            if cohort is not None and row["cohort"] != cohort:
                continue
            if all(row[key] == value for key, value in equals.items()):
                chosen.append(row)
        return chosen

    cohorts = ["GSE183647", "GSE189521"]
    for cohort in cohorts:
        rows = subset_rows(paired, cohort)
        add(
            "event_L_minus_rna_implied_L", cohort, "all", "all",
            [row["event_L"] for row in rows],
            [row["L_rna"] for row in rows],
            [row["delta_L"] for row in rows],
            extra_gap(rows),
        )
        groups = sorted({row["group"] for row in rows})
        for group in groups:
            part = subset_rows(rows, group=group)
            add(
                "event_L_minus_rna_implied_L", cohort, "group", group,
                [row["event_L"] for row in part],
                [row["L_rna"] for row in part],
                [row["delta_L"] for row in part],
                extra_gap(part),
            )
        for ruler, label in (("22q", "22q"),):
            part = subset_rows(rows, event_arm=ruler)
            add(
                "event_L_minus_rna_implied_L", cohort, "ruler", label,
                [row["event_L"] for row in part],
                [row["L_rna"] for row in part],
                [row["delta_L"] for row in part],
                extra_gap(part),
            )
        other = [row for row in rows if row["event_arm"] != "22q"]
        add(
            "event_L_minus_rna_implied_L", cohort, "ruler", "other",
            [row["event_L"] for row in other],
            [row["L_rna"] for row in other],
            [row["delta_L"] for row in other],
            extra_gap(other),
        )
        for flag in ("yes", "no"):
            part = subset_rows(rows, ruler_is_deepest=flag)
            add(
                "event_L_minus_rna_implied_L", cohort, "ruler_is_deepest", flag,
                [row["event_L"] for row in part],
                [row["L_rna"] for row in part],
                [row["delta_L"] for row in part],
                extra_gap(part),
            )
        for stratum in ("nf2_loss", "nf2_not_loss"):
            part = subset_rows(rows, nf2_stratum=stratum)
            add(
                "event_L_minus_rna_implied_L", cohort, "nf2", stratum,
                [row["event_L"] for row in part],
                [row["L_rna"] for row in part],
                [row["delta_L"] for row in part],
                extra_gap(part),
            )
        for stratum, ruler_name, nf2_name in (
            ("22q_and_nf2_loss", "22q", "nf2_loss"),
            ("22q_and_nf2_not_loss", "22q", "nf2_not_loss"),
        ):
            part = [row for row in rows if row["event_arm"] == "22q" and row["nf2_stratum"] == nf2_name]
            add(
                "event_L_minus_rna_implied_L", cohort, "nf2_and_ruler", stratum,
                [row["event_L"] for row in part],
                [row["L_rna"] for row in part],
                [row["delta_L"] for row in part],
                extra_gap(part),
            )
        held = subset_rows(rows, nf2_stratum="nf2_inframe_excluded")
        if held:
            add(
                "event_L_minus_rna_implied_L", cohort, "nf2", "nf2_inframe_excluded",
                [row["event_L"] for row in held],
                [row["L_rna"] for row in held],
                [row["delta_L"] for row in held],
                extra_gap(held),
            )

    for cohort in cohorts:
        rows = [row for row in unestimable if row["cohort"] == cohort and np.isfinite(row["rna_neoplastic"])]
        add(
            "unestimable_rna_neoplastic", cohort, "all", "all",
            [row["rna_neoplastic"] for row in rows],
            [row["rna_neoplastic"] for row in rows],
            [row["rna_neoplastic"] for row in rows],
        )
        summary[-1]["fraction_positive"] = np.nan
        for group in sorted({row["group"] for row in rows}):
            part = [row for row in rows if row["group"] == group]
            add(
                "unestimable_rna_neoplastic", cohort, "group", group,
                [row["rna_neoplastic"] for row in part],
                [row["rna_neoplastic"] for row in part],
                [row["rna_neoplastic"] for row in part],
            )
            summary[-1]["fraction_positive"] = np.nan

    for cohort in ("GSE183647", "GSE189521", "GSE180061"):
        rows = [row for row in ruler_gaps if row["cohort"] == cohort]
        add(
            "event_L_22q_minus_deepest_other", cohort, "all", "22q_ruler_with_deeper_arm",
            [row["event_L"] for row in rows],
            [row["deepest_L"] for row in rows],
            [row["difference"] for row in rows],
        )

    summary_path = WORK / "purity_scale_summary.csv"
    fieldnames = [
        "contrast", "cohort", "stratum_type", "stratum", "n",
        "median_observed", "median_reference", "median_difference",
        "ci_low", "ci_high", "fraction_positive",
        "median_intensity_purity", "median_rna_neoplastic", "median_purity_gap",
        "purity_gap_ci_low", "purity_gap_ci_high",
    ]
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in summary:
            writer.writerow(row)
    print("wrote", paired_path)
    print("wrote", summary_path)
    for row in summary:
        if row["stratum_type"] in ("all", "ruler", "ruler_is_deepest", "nf2", "nf2_and_ruler") or row["contrast"].startswith("event_L_22q"):
            print(
                row["contrast"], row["cohort"], row["stratum_type"], row["stratum"],
                "n", row["n"],
                "diff", round(row["median_difference"], 3) if np.isfinite(row["median_difference"]) else None,
                "ci", None if not np.isfinite(row["ci_low"]) else (round(row["ci_low"], 3), round(row["ci_high"], 3)),
                "frac", None if not np.isfinite(row["fraction_positive"]) else round(row["fraction_positive"], 3),
            )


def extra_gap(rows):
    n, median, lo, hi = bootstrap_median([row["purity_gap"] for row in rows])
    _, med_intensity, _, _ = bootstrap_median([row["intensity_purity"] for row in rows])
    _, med_rna, _, _ = bootstrap_median([row["rna_neoplastic"] for row in rows])
    return {
        "median_intensity_purity": med_intensity,
        "median_rna_neoplastic": med_rna,
        "median_purity_gap": median,
        "purity_gap_ci_low": lo,
        "purity_gap_ci_high": hi,
    }


if __name__ == "__main__":
    main()
