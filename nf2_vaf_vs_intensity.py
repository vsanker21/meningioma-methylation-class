"""NF2 allele frequency against intensity purity in Nassiri.

The filter and the two-hit formula are fixed in docs/protocol_locks.md.
The intensity estimator is not revised.
"""

import csv
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

ROOT = Path(r"F:\Meningioma Genetics")
MUT = ROOT / "data" / "raw" / "clinical" / "nassiri_2021_cbioportal" / "data_mutations.txt"
PURITY = ROOT / "data" / "work" / "purity_estimates.csv"
OUT_TUMOUR = ROOT / "data" / "work" / "nf2_vaf_tumours.csv"
OUT_SUMMARY = ROOT / "data" / "work" / "nf2_vaf_summary.csv"
N_BOOT = 2000
SEED = 1


def fnum(value):
    if value is None or value in ("", "NA", "NaN", "."):
        return np.nan
    return float(value)


def bootstrap_spearman(x, y):
    ok = np.isfinite(x) & np.isfinite(y)
    x = x[ok]
    y = y[ok]
    n = int(x.size)
    if n < 10:
        return n, np.nan, np.nan, np.nan
    rho = float(spearmanr(x, y).statistic)
    rng = np.random.default_rng(SEED)
    draws = np.empty(N_BOOT)
    for i in range(N_BOOT):
        take = rng.integers(0, n, n)
        draws[i] = spearmanr(x[take], y[take]).statistic
    lo, hi = np.quantile(draws, [0.025, 0.975])
    return n, rho, float(lo), float(hi)


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
    draws = np.median(values[rng.integers(0, n, size=(N_BOOT, n))], axis=1)
    lo, hi = np.quantile(draws, [0.025, 0.975])
    return n, median, float(lo), float(hi)


def main():
    with MUT.open(encoding="utf-8", newline="") as handle:
        calls = [row for row in csv.DictReader(handle, delimiter="\t") if row["Hugo_Symbol"] == "NF2"]
    kept = []
    for row in calls:
        depth = fnum(row["t_depth"])
        alt = fnum(row["t_alt_count"])
        ref = fnum(row["t_ref_count"])
        normal = fnum(row["normal_vaf"])
        if not (depth >= 30 and normal <= 0.05 and alt + ref > 0):
            continue
        vaf = alt / (alt + ref)
        stated = fnum(row["tumor_vaf"])
        if np.isfinite(stated) and abs(stated - vaf) > 0.02:
            raise SystemExit("%s stated VAF %s disagrees with counts %s" % (row["Tumor_Sample_Barcode"], stated, vaf))
        kept.append({
            "sample_id": row["Tumor_Sample_Barcode"],
            "depth": depth,
            "vaf": vaf,
            "rho_twohit": 2 * vaf / (1 + vaf),
            "rho_diploid": 2 * vaf if vaf <= 0.5 else np.nan,
        })
    best = {}
    for row in kept:
        previous = best.get(row["sample_id"])
        if previous is None or row["depth"] > previous["depth"]:
            best[row["sample_id"]] = row
    print("somatic NF2 tumours", len(best), "from calls", len(calls))

    with PURITY.open(encoding="utf-8", newline="") as handle:
        purity = {row["sample_id"]: row for row in csv.DictReader(handle) if row["cohort"] == "GSE180061"}
    joined = []
    missing = []
    for sample_id, row in best.items():
        hit = purity.get(sample_id)
        if hit is None:
            missing.append(sample_id)
            continue
        row = dict(row)
        row["purity_status"] = hit["purity_status"]
        row["intensity_purity"] = fnum(hit["purity"])
        row["l_22q"] = fnum(hit["l_22q"])
        row["event_arm"] = hit["event_arm"]
        row["group"] = hit["group"]
        row["has_22q_loss"] = "yes" if np.isfinite(row["l_22q"]) else "no"
        joined.append(row)
    if missing:
        raise SystemExit("NF2 samples absent from the purity table: %s" % ",".join(missing))

    with OUT_TUMOUR.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(joined[0]))
        writer.writeheader()
        writer.writerows(joined)

    primary = [row for row in joined if row["has_22q_loss"] == "yes" and row["purity_status"] == "estimated"]
    quiet = [row for row in joined if row["purity_status"] != "estimated"]
    n, rho, lo, hi = bootstrap_spearman(
        np.array([row["rho_twohit"] for row in primary]),
        np.array([row["intensity_purity"] for row in primary]),
    )
    _, med_vaf, vaf_lo, vaf_hi = bootstrap_median([row["vaf"] for row in primary])
    _, med_two, two_lo, two_hi = bootstrap_median([row["rho_twohit"] for row in primary])
    _, med_int, int_lo, int_hi = bootstrap_median([row["intensity_purity"] for row in primary])
    _, med_gap, gap_lo, gap_hi = bootstrap_median(
        [row["rho_twohit"] - row["intensity_purity"] for row in primary]
    )
    _, med_quiet, quiet_lo, quiet_hi = bootstrap_median([row["rho_twohit"] for row in quiet])
    n_diploid = sum(np.isfinite(row["rho_diploid"]) for row in primary)
    summary = [
        {"contrast": "spearman_twohit_vs_intensity", "n": n, "estimate": rho, "ci_low": lo, "ci_high": hi},
        {"contrast": "median_vaf_with_22q_loss", "n": len(primary), "estimate": med_vaf, "ci_low": vaf_lo, "ci_high": vaf_hi},
        {"contrast": "median_twohit_fraction_with_22q_loss", "n": len(primary), "estimate": med_two, "ci_low": two_lo, "ci_high": two_hi},
        {"contrast": "median_intensity_with_22q_loss_and_nf2", "n": len(primary), "estimate": med_int, "ci_low": int_lo, "ci_high": int_hi},
        {"contrast": "median_twohit_minus_intensity", "n": len(primary), "estimate": med_gap, "ci_low": gap_lo, "ci_high": gap_hi},
        {"contrast": "median_twohit_fraction_no_intensity_arm", "n": len(quiet), "estimate": med_quiet, "ci_low": quiet_lo, "ci_high": quiet_hi},
        {"contrast": "primary_tumours_with_defined_diploid_reading", "n": n_diploid, "estimate": np.nan, "ci_low": np.nan, "ci_high": np.nan},
    ]
    with OUT_SUMMARY.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary[0]))
        writer.writeheader()
        writer.writerows(summary)
    for row in summary:
        print(row["contrast"], "n", row["n"], "est", row["estimate"], "ci", row["ci_low"], row["ci_high"])
    print("wrote", OUT_SUMMARY)


if __name__ == "__main__":
    main()
