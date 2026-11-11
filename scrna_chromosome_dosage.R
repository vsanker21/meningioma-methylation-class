# Chromosome dosage in the Maas single-cell object.
#
# The contrast is fixed in docs/protocol_locks.md. Bulk purity is not revised.

suppressPackageStartupMessages(library(Seurat))

root <- "F:/Meningioma Genetics"
annot_path <- file.path(root, "data/raw/GSE/GSE252291/Human.GRCh38.p13.annot (3).tsv.gz")
work <- file.path(root, "data/work")

chrom_of <- c(
  "NC_000001.11" = "1", "NC_000002.12" = "2", "NC_000003.12" = "3",
  "NC_000004.12" = "4", "NC_000005.10" = "5", "NC_000006.12" = "6",
  "NC_000007.14" = "7", "NC_000008.11" = "8", "NC_000009.12" = "9",
  "NC_000010.11" = "10", "NC_000011.10" = "11", "NC_000012.12" = "12",
  "NC_000013.11" = "13", "NC_000014.9" = "14", "NC_000015.10" = "15",
  "NC_000016.10" = "16", "NC_000017.11" = "17", "NC_000018.10" = "18",
  "NC_000019.10" = "19", "NC_000020.11" = "20", "NC_000021.9" = "21",
  "NC_000022.11" = "22"
)

message("reading gene map")
annot <- utils::read.delim(gzfile(annot_path), stringsAsFactors = FALSE, quote = "")
annot <- annot[annot$GeneType == "protein-coding" & annot$ChrAcc %in% names(chrom_of), ]
annot <- annot[!duplicated(annot$Symbol), ]
annot$chrom <- chrom_of[annot$ChrAcc]

message("loading single-cell counts")
obj <- readRDS(file.path(root, "data/raw/scrna/zenodo_18863971/seurat_single_cell.rds"))
meta <- obj@meta.data
labels <- as.character(meta$merged_annot)
keep <- labels %in% c("neoplastic", "macrophage")
assay <- obj[["RNA"]]
counts <- if ("counts" %in% Layers(assay)) LayerData(assay, layer = "counts") else GetAssayData(obj, assay = "RNA", slot = "counts")
counts <- counts[, keep, drop = FALSE]
meta <- meta[keep, , drop = FALSE]
labels <- labels[keep]
rm(obj, assay)
gc()

shared <- intersect(rownames(counts), annot$Symbol)
counts <- counts[shared, , drop = FALSE]
annot <- annot[match(shared, annot$Symbol), ]
keep_gene <- Matrix::rowSums(counts) >= 100
counts <- counts[keep_gene, , drop = FALSE]
annot <- annot[keep_gene, , drop = FALSE]
message("genes ", nrow(counts), " cells ", ncol(counts), " chr22 ", sum(annot$chrom == "22"), " chr2 ", sum(annot$chrom == "2"))

lib <- Matrix::colSums(counts)
counts <- Matrix::t(Matrix::t(counts) / lib) * 10000
counts@x <- log1p(counts@x)

score_of <- function(chrom, drop_gene = NULL) {
  on <- annot$chrom == chrom
  if (!is.null(drop_gene)) on <- on & annot$Symbol != drop_gene
  bg <- annot$chrom != chrom
  Matrix::colMeans(counts[on, , drop = FALSE]) - Matrix::colMeans(counts[bg, , drop = FALSE])
}

scores <- list(
  chr22 = score_of("22"),
  chr22_no_nf2 = score_of("22", "NF2"),
  chr2 = score_of("2")
)

tumours <- split(seq_along(labels), meta$orig.ident)
rows <- list()
for (tumour in names(tumours)) {
  idx <- tumours[[tumour]]
  neo <- idx[labels[idx] == "neoplastic"]
  mac <- idx[labels[idx] == "macrophage"]
  if (length(neo) < 20 || length(mac) < 20) next
  one <- data.frame(tumour = tumour, n_neoplastic = length(neo), n_macrophage = length(mac), stringsAsFactors = FALSE)
  for (name in names(scores)) {
    nscore <- scores[[name]][neo]
    mscore <- scores[[name]][mac]
    mm <- median(mscore)
    md <- mad(mscore)
    one[[paste0(name, "_neo_median")]] <- median(nscore)
    one[[paste0(name, "_mac_median")]] <- mm
    one[[paste0(name, "_difference")]] <- median(nscore) - mm
    one[[paste0(name, "_neo_mad")]] <- mad(nscore)
    one[[paste0(name, "_mac_mad")]] <- md
    one[[paste0(name, "_fraction_below")]] <- if (is.finite(md) && md > 0) mean(nscore < mm - 2 * md) else NA_real_
  }
  rows[[length(rows) + 1]] <- one
}
per_tumour <- do.call(rbind, rows)
utils::write.csv(per_tumour, file.path(work, "scrna_dosage_by_tumour.csv"), row.names = FALSE)
message("tumours ", nrow(per_tumour))

bootstrap_median <- function(values) {
  values <- values[is.finite(values)]
  n <- length(values)
  if (n < 10) return(c(n = n, median = median(values), lo = NA, hi = NA))
  set.seed(1)
  draws <- replicate(2000, median(values[sample.int(n, replace = TRUE)]))
  ci <- stats::quantile(draws, c(0.025, 0.975))
  c(n = n, median = median(values), lo = unname(ci[1]), hi = unname(ci[2]))
}

summary_rows <- list()
for (name in names(scores)) {
  for (measure in c("difference", "fraction_below")) {
    stats <- bootstrap_median(per_tumour[[paste0(name, "_", measure)]])
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      score = name, measure = measure,
      n = unname(stats["n"]), median = unname(stats["median"]),
      ci_low = unname(stats["lo"]), ci_high = unname(stats["hi"]),
      stringsAsFactors = FALSE
    )
  }
}
result <- do.call(rbind, summary_rows)
utils::write.csv(result, file.path(work, "scrna_dosage_summary.csv"), row.names = FALSE)
print(result)
message("wrote dosage summaries")
