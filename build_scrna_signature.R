# Meningioma single-cell reference for the RNA purity check.
#
# Cell types are the merged_annot labels in the Maas Zenodo Seurat object.
# That object has no fibroblast cluster, so none is added from a generic panel.
# The signature is the pseudobulk count per million of each label. Markers are
# chosen here, from this object only: a gene must be at least 4-fold above the
# next-highest label and at least 20 counts per million, and at most 50 genes
# are kept per label, ranked by that fold change. If a label yields fewer than
# 10 genes, the floor is lowered to 5 counts per million and the fold rule is
# kept. The bulk matrices are not read here.

suppressPackageStartupMessages(library(Seurat))

obj <- readRDS("F:/Meningioma Genetics/data/raw/scrna/zenodo_18863971/seurat_single_cell.rds")
assay <- obj[["RNA"]]
counts <- if ("counts" %in% Layers(assay)) LayerData(assay, layer = "counts") else GetAssayData(obj, assay = "RNA", slot = "counts")
labels <- as.character(obj$merged_annot)
message("counts ", nrow(counts), " x ", ncol(counts))

types <- names(sort(table(labels), decreasing = TRUE))
pseudo <- sapply(types, function(type) Matrix::rowSums(counts[, labels == type, drop = FALSE]))
lib <- sapply(types, function(type) sum(obj$nCount_RNA[labels == type]))
cpm <- sweep(pseudo, 2, lib, "/") * 1e6
mode(cpm) <- "numeric"

pick <- function(floor_cpm) {
  chosen <- list()
  for (type in types) {
    own <- cpm[, type]
    other <- do.call(pmax, as.data.frame(cpm[, setdiff(types, type), drop = FALSE]))
    fold <- own / pmax(other, 1e-8)
    ok <- which(own >= floor_cpm & fold >= 4)
    ok <- ok[order(fold[ok], decreasing = TRUE)]
    ok <- head(ok, 50)
    if (!length(ok)) next
    chosen[[type]] <- data.frame(
      gene = rownames(cpm)[ok],
      cell_type = type,
      cpm = own[ok],
      fold = fold[ok],
      floor_cpm = floor_cpm,
      stringsAsFactors = FALSE
    )
  }
  chosen
}

chosen <- pick(20)
n_each <- sapply(types, function(type) if (type %in% names(chosen)) nrow(chosen[[type]]) else 0)
message("markers at floor 20: ", paste(types, n_each, sep = "=", collapse = ", "))
thin <- types[n_each < 10]
if (length(thin)) {
  message("lowering the floor to 5 for: ", paste(thin, collapse = ", "))
  relaxed <- pick(5)
  for (type in thin) chosen[[type]] <- relaxed[[type]]
}
markers <- do.call(rbind, chosen)
out_dir <- "F:/Meningioma Genetics/data/work"
dir.create(out_dir, showWarnings = FALSE)
utils::write.csv(cpm, file.path(out_dir, "scrna_pseudobulk_cpm.csv"), quote = TRUE)
utils::write.csv(markers, file.path(out_dir, "scrna_markers.csv"), row.names = FALSE)
counts_out <- data.frame(cell_type = types, n_cells = as.integer(table(labels)[types]), library_counts = lib)
utils::write.csv(counts_out, file.path(out_dir, "scrna_celltype_counts.csv"), row.names = FALSE)
message("markers ", nrow(markers), " types ", paste(unique(markers$cell_type), collapse = ","))
message("wrote signature")
