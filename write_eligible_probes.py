"""Autosomal cg probes that are not on the published SNP or cross-reactive lists."""

import csv
import gzip
from pathlib import Path

root = Path(r"F:\Meningioma Genetics")
manifest = root / "data" / "raw" / "annotation" / "GPL21145_MethylationEPIC_15073387_v-1-0.csv.gz"
filter_dir = root / "data" / "raw" / "github" / "mnp_training" / "filter"
out = root / "data" / "work" / "eligible_autosomal_cg.txt"

blocked = set()
for name in (
    "snp_7998probes.vh20151030.txt",
    "amb_3965probes.vh20151030.txt",
    "xy_11551probes.vh20151030.txt",
):
    for line in (filter_dir / name).read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line:
            blocked.add(line)

kept = 0
with gzip.open(manifest, "rt", encoding="utf-8", errors="replace", newline="") as handle, out.open(
    "w", encoding="utf-8"
) as dest:
    started = False
    for row in csv.reader(handle):
        if not started:
            if row and row[0] == "IlmnID":
                name_i = row.index("Name")
                chr_i = row.index("CHR")
                started = True
            continue
        if not row or row[0].startswith("["):
            break
        if len(row) <= chr_i:
            continue
        probe = row[name_i]
        chrom = row[chr_i]
        if not probe.startswith("cg"):
            continue
        if chrom in {"X", "Y", "x", "y", "chrX", "chrY"}:
            continue
        if probe in blocked:
            continue
        dest.write(probe + "\n")
        kept += 1
print("eligible", kept)
