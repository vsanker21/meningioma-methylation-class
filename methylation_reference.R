# Reference profiles for the secondary methylation estimator.
#
# Monocyte and lymphocyte profiles are medians from FlowSorted.Blood.EPIC
# purified cells, processed with preprocessNoob, the same normalisation used
# for the tumours. The lymphocyte profile is the mean of the B-cell, CD4 T,
# CD8 T and NK medians, so one subtype is not weighted by how many samples it
# has. The endothelial profile is the median of the two Moss vascular
# endothelial replicates. Neutrophils, artificial mixtures and fibroblasts
# are not components.
#
# A probe qualifies for a lineage only when that lineage's beta differs from
# both of the other two by at least 0.5 in the same direction. Endothelial
# probes also have to put both replicates on the same side of both immune
# medians. At most 100 probes are kept per lineage, ranked by the smaller of
# the two differences. If a lineage has fewer than 20 probes, it is reselected
# at 0.4 with the same-direction rule kept. These probes are later excluded
# from clustering.

suppressPackageStartupMessages(library(minfi))

root <- "F:/Meningioma Genetics"
work <- file.path(root, "data/work")
ref_dir <- file.path(root, "data/raw/methylation/reference")

one_beta <- function(basename) {
  message(basename)
  rg <- read.metharray(basename, verbose = FALSE, force = TRUE)
  rg@annotation <- c(array = "IlluminaHumanMethylationEPIC", annotation = "ilm10b2.hg19")
  getBeta(preprocessNoob(rg))[, 1]
}

median_profile <- function(basenames) {
  columns <- lapply(basenames, one_beta)
  probes <- Reduce(intersect, lapply(columns, names))
  mat <- sapply(columns, function(values) values[probes])
  setNames(apply(mat, 1, median, na.rm = TRUE), probes)
}

message("monocyte and lymphocyte medians")
samples <- utils::read.csv(file.path(work, "methylation_reference_samples.csv"), stringsAsFactors = FALSE)
blood <- function(cell) file.path(ref_dir, "blood_idat", samples$sentrix[samples$cell_type == cell])
mono <- median_profile(blood("Mono"))
bcell <- median_profile(blood("Bcell"))
cd4 <- median_profile(blood("CD4T"))
cd8 <- median_profile(blood("CD8T"))
nk <- median_profile(blood("NK"))
shared <- Reduce(intersect, list(names(bcell), names(cd4), names(cd8), names(nk)))
lymph <- rowMeans(cbind(bcell[shared], cd4[shared], cd8[shared], nk[shared]))

message("endothelial replicates")
endo_files <- file.path(ref_dir, c("GSM3455806", "GSM3455822"))
endo_reps <- lapply(endo_files, one_beta)
endo_probes <- Reduce(intersect, lapply(endo_reps, names))
endo <- setNames(apply(sapply(endo_reps, function(values) values[endo_probes]), 1, median, na.rm = TRUE), endo_probes)

probes <- Reduce(intersect, list(names(mono), names(lymph), names(endo), scan(
  file.path(work, "eligible_autosomal_cg.txt"), what = "", quiet = TRUE
)))
message("shared eligible probes ", length(probes))
profiles <- cbind(monocyte = mono[probes], lymphocyte = lymph[probes], endothelial = endo[probes])
rep_mat <- cbind(endo_reps[[1]][probes], endo_reps[[2]][probes])

same_side <- function(target, other_a, other_b) {
  d1 <- target - other_a
  d2 <- target - other_b
  same <- sign(d1) == sign(d2) & sign(d1) != 0
  score <- pmin(abs(d1), abs(d2))
  score[!same | is.na(same)] <- NA_real_
  score
}

select_lineage <- function(lineage, threshold) {
  others <- setdiff(colnames(profiles), lineage)
  score <- same_side(profiles[, lineage], profiles[, others[1]], profiles[, others[2]])
  if (lineage == "endothelial") {
    immune <- profiles[, c("monocyte", "lymphocyte")]
    for (rep_i in 1:2) {
      rep_score <- same_side(rep_mat[, rep_i], immune[, 1], immune[, 2])
      score[is.na(rep_score)] <- NA_real_
    }
  }
  keep <- which(!is.na(score) & score >= threshold)
  keep <- keep[order(score[keep], decreasing = TRUE)]
  head(keep, 100)
}

chosen <- list()
for (lineage in colnames(profiles)) {
  idx <- select_lineage(lineage, 0.5)
  threshold_used <- 0.5
  if (length(idx) < 20) {
    idx <- select_lineage(lineage, 0.4)
    threshold_used <- 0.4
    message(lineage, " fell back to 0.4, kept ", length(idx))
  }
  if (length(idx) < 20) {
    stop(lineage, " still has only ", length(idx), " probes")
  }
  chosen[[lineage]] <- data.frame(
    probe = probes[idx],
    lineage = lineage,
    threshold = threshold_used,
    stringsAsFactors = FALSE
  )
  message(lineage, " probes ", length(idx), " at ", threshold_used)
}

probe_table <- do.call(rbind, chosen)
utils::write.csv(probe_table, file.path(work, "methylation_reference_probes.csv"), row.names = FALSE)
utils::write.csv(profiles[unique(probe_table$probe), ], file.path(work, "methylation_reference_betas.csv"), row.names = TRUE)
message("wrote reference probes and betas")
