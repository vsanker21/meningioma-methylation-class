"""RNA neoplastic fraction versus intensity purity.

The marker list and the pseudobulk signature are built from the single-cell
object only. This script does not change that list.
"""

import csv
import gzip
from pathlib import Path

import numpy as np
from scipy.optimize import nnls
from scipy.stats import spearmanr

ROOT = Path(r"F:\Meningioma Genetics")
WORK = ROOT / "data" / "work"
MARKERS = WORK / "scrna_markers.csv"
SIGNATURE = WORK / "scrna_pseudobulk_cpm.csv"
PURITY = WORK / "purity_estimates.csv"
TPM = ROOT / "data" / "raw" / "rna" / "GSE183653" / "GSE183653_meningioma_tpm.csv.gz"
COUNTS = ROOT / "data" / "raw" / "rna" / "GSE189672" / "GSE189672_RNAseq_raw_counts.csv.gz"
SHEET = ROOT / "data" / "raw" / "methylation" / "GSE189521" / "GSE189521_Methylation_sample_sheet.csv.gz"
BAYLEY_KEY = ROOT / "data" / "raw" / "clinical" / "bayley_2022_sciadv" / "sample_clinical_key.csv"
ANNOT = ROOT / "data" / "raw" / "GSE" / "GSE252291" / "Human.GRCh38.p13.annot (3).tsv.gz"
OUT = WORK / "rna_neoplastic_fraction.csv"
CONCORD = WORK / "purity_concordance_rna.csv"
RNG = np.random.default_rng(1)


def read_csv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def load_signature():
    markers = read_csv(MARKERS)
    by_type = {}
    for row in markers:
        by_type.setdefault(row["cell_type"], []).append(row["gene"])
    with SIGNATURE.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        types = header[1:]
        cpm = {}
        for row in reader:
            cpm[row[0]] = {types[i]: float(row[i + 1]) for i in range(len(types))}
    genes = []
    seen = set()
    for type_name in types:
        for gene in by_type.get(type_name, []):
            if gene not in seen and gene in cpm:
                genes.append(gene)
                seen.add(gene)
    matrix = np.column_stack([[cpm[gene][type_name] for gene in genes] for type_name in types])
    return genes, types, matrix


def bayley_lab_to_sample():
    with gzip.open(SHEET, "rt", encoding="utf-8", errors="replace") as handle:
        lines = handle.read().splitlines()
    start = next(i for i, line in enumerate(lines) if line.startswith("Sample_Name"))
    rows = list(csv.DictReader(lines[start:]))
    lab_to_sentrix = {}
    for row in rows:
        row = {key.strip(): value.strip() for key, value in row.items() if key}
        lab = row["Sample_Name"].split("#")[-1]
        sentrix = row["Sentrix_ID"] + "_" + row["Sentrix_Position"]
        lab_to_sentrix[lab] = sentrix
    sentrix_to_sample = {row["sentrix"]: row["Sample"] for row in read_csv(BAYLEY_KEY)}
    joined = {}
    missing = []
    for lab, sentrix in lab_to_sentrix.items():
        sample = sentrix_to_sample.get(sentrix)
        if sample is None:
            missing.append(lab)
        else:
            joined[lab] = sample
    print("bayley sheet", len(lab_to_sentrix), "joined", len(joined), "unjoined", len(missing))
    return joined


def ensembl_to_symbol():
    mapping = {}
    with gzip.open(ANNOT, "rt", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            symbol = row["Symbol"].strip()
            ensembl = row["EnsemblGeneID"].strip()
            if symbol and ensembl and ensembl not in mapping:
                mapping[ensembl] = symbol
    return mapping


def stream_markers(path, genes, id_to_symbol=None):
    wanted = set(genes)
    found = {gene: {} for gene in genes}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        header = handle.readline().rstrip("\n").split(",")
        samples = header[1:]
        library = np.zeros(len(samples))
        for line in handle:
            parts = line.rstrip("\n").split(",")
            gene_id = parts[0].strip('"')
            values = np.fromiter((float(x) if x else 0.0 for x in parts[1:]), dtype=float, count=len(samples))
            library += values
            symbol = gene_id if id_to_symbol is None else id_to_symbol.get(gene_id.split(".")[0], "")
            if symbol in wanted:
                previous = found[symbol]
                if not previous:
                    found[symbol] = values
                else:
                    found[symbol] = previous + values
    present = [gene for gene in genes if isinstance(found[gene], np.ndarray)]
    matrix = np.vstack([found[gene] for gene in present])
    return samples, present, matrix, library


def fractions_for(samples, genes, values, library, signature_genes, types, signature, as_cpm):
    index = {gene: i for i, gene in enumerate(genes)}
    keep = [gene for gene in signature_genes if gene in index]
    if len(keep) < 20:
        raise RuntimeError("only %d marker genes in this matrix" % len(keep))
    observed = values[[index[gene] for gene in keep], :]
    if as_cpm:
        observed = observed / library * 1e6
    design = signature[[signature_genes.index(gene) for gene in keep], :]
    out = []
    for i, sample in enumerate(samples):
        coef, _ = nnls(design, observed[:, i])
        total = coef.sum()
        frac = coef / total if total > 0 else np.full(len(types), np.nan)
        row = dict(zip(types, frac))
        row["sample_id"] = sample
        row["n_marker_genes"] = len(keep)
        out.append(row)
    return out


def bootstrap_spearman(x, y):
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]
    n = len(x)
    if n < 5:
        return n, np.nan, np.nan, np.nan
    rho = spearmanr(x, y).statistic
    draws = np.empty(2000)
    for i in range(2000):
        take = RNG.integers(0, n, n)
        draws[i] = spearmanr(x[take], y[take]).statistic
    low, high = np.quantile(draws, [0.025, 0.975])
    return n, rho, low, high


def main():
    genes, types, signature = load_signature()
    if "neoplastic" not in types:
        raise RuntimeError("signature has no neoplastic column")
    print("markers", len(genes), "types", ",".join(types))
    lab_to_sample = bayley_lab_to_sample()
    symbol_of = ensembl_to_symbol()

    tpm_samples, tpm_genes, tpm_values, _ = stream_markers(TPM, genes)
    count_samples, count_genes, count_values, count_library = stream_markers(COUNTS, genes, symbol_of)
    print("choudhury marker genes", len(tpm_genes), "bayley marker genes", len(count_genes))

    choudhury = fractions_for(tpm_samples, tpm_genes, tpm_values, None, genes, types, signature, as_cpm=False)
    for row in choudhury:
        row["cohort"] = "GSE183647"
    bayley_rows = fractions_for(count_samples, count_genes, count_values, count_library, genes, types, signature, as_cpm=True)
    bayley = []
    unmatched = []
    for row in bayley_rows:
        sample = lab_to_sample.get(row["sample_id"])
        if sample is None:
            unmatched.append(row["sample_id"])
            continue
        row["lab_id"] = row["sample_id"]
        row["sample_id"] = sample
        row["cohort"] = "GSE189521"
        bayley.append(row)
    print("bayley counts", len(bayley_rows), "joined", len(bayley), "unmatched", len(unmatched))
    if unmatched:
        print(" unmatched", ",".join(unmatched[:15]))

    purity = {row["sample_id"]: row for row in read_csv(PURITY)}
    combined = choudhury + bayley
    fieldnames = ["sample_id", "cohort", "lab_id", "n_marker_genes", "intensity_purity", "intensity_status"] + types
    with OUT.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in combined:
            hit = purity.get(row["sample_id"])
            row["intensity_purity"] = hit["purity"] if hit else ""
            row["intensity_status"] = hit["purity_status"] if hit else ""
            writer.writerow(row)

    summary = []
    for cohort, rows in (("GSE183647", choudhury), ("GSE189521", bayley)):
        xs, ys, unest = [], [], []
        for row in rows:
            hit = purity.get(row["sample_id"])
            if hit is None or not np.isfinite(row["neoplastic"]):
                continue
            if hit["purity_status"] == "estimated" and hit["purity"] != "":
                xs.append(float(hit["purity"]))
                ys.append(float(row["neoplastic"]))
            elif hit["purity_status"] != "estimated":
                unest.append(float(row["neoplastic"]))
        n, rho, low, high = bootstrap_spearman(np.array(xs), np.array(ys))
        summary.append({
            "cohort": cohort,
            "n_paired": n,
            "spearman": "" if not np.isfinite(rho) else f"{rho:.4f}",
            "ci_low": "" if not np.isfinite(low) else f"{low:.4f}",
            "ci_high": "" if not np.isfinite(high) else f"{high:.4f}",
            "n_unestimable_with_rna": len(unest),
            "median_rna_neoplastic_unestimable": "" if not unest else f"{np.median(unest):.4f}",
            "median_rna_neoplastic_estimable": "" if not ys else f"{np.median(ys):.4f}",
        })
        print(cohort, summary[-1])
    with CONCORD.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary[0]))
        writer.writeheader()
        writer.writerows(summary)
    print("wrote", OUT)
    print("wrote", CONCORD)


if __name__ == "__main__":
    main()
