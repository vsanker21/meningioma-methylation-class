# Tumour betas at the reference-selected probes, then the constrained projection.
#
# Resumable: a sample already present in the beta file is skipped. The
# projection and the concordance are written only once every cohort sample
# has a row, so a partial file cannot be mistaken for the result.
#
# Coefficients are non-negative and sum to at most 1. The remainder is
# unassigned. It is not called tumour.

suppressPackageStartupMessages({
  library(minfi)
  library(quadprog)
})

root <- "F:/Meningioma Genetics"
work <- file.path(root, "data/work")
beta_out <- file.path(work, "methylation_tumour_betas.csv")
fraction_out <- file.path(work, "methylation_fractions.csv")
concordance_out <- file.path(work, "purity_concordance_methylation.csv")

probe_table <- utils::read.csv(file.path(work, "methylation_reference_probes.csv"), stringsAsFactors = FALSE)
reference <- utils::read.csv(file.path(work, "methylation_reference_betas.csv"), row.names = 1, stringsAsFactors = FALSE)
probes <- unique(probe_table$probe)
reference <- reference[probes, c("monocyte", "lymphocyte", "endothelial")]

manifest <- utils::read.csv(file.path(work, "purity_run_manifest.csv"), stringsAsFactors = FALSE)
done <- character(0)
if (file.exists(beta_out)) {
  done <- utils::read.csv(beta_out, stringsAsFactors = FALSE)$sample_id
}

one_beta <- function(basename) {
  rg <- read.metharray(basename, verbose = FALSE, force = TRUE)
  rg@annotation <- c(array = "IlluminaHumanMethylationEPIC", annotation = "ilm10b2.hg19")
  getBeta(preprocessNoob(rg))[, 1]
}

for (i in seq_len(nrow(manifest))) {
  sample_id <- manifest$sample_id[i]
  if (sample_id %in% done) next
  message(i, "/", nrow(manifest), " ", sample_id)
  values <- tryCatch(one_beta(manifest$basename[i])[probes], error = function(err) {
    message("failed ", sample_id, ": ", conditionMessage(err))
    setNames(rep(NA_real_, length(probes)), probes)
  })
  row <- data.frame(sample_id = sample_id, t(as.numeric(values[probes])), check.names = FALSE)
  colnames(row) <- c("sample_id", probes)
  utils::write.table(
    row, beta_out, sep = ",", row.names = FALSE, quote = FALSE,
    col.names = !file.exists(beta_out), append = file.exists(beta_out)
  )
  if (i %% 25 == 0) gc(verbose = FALSE)
}

betas <- utils::read.csv(beta_out, stringsAsFactors = FALSE, check.names = FALSE)
if (nrow(betas) < nrow(manifest)) {
  stop("beta file has ", nrow(betas), " rows; concordance is not written until all ", nrow(manifest), " exist")
}

project <- function(y) {
  ok <- is.finite(y) & apply(reference, 1, function(row) all(is.finite(row)))
  if (sum(ok) < 50) return(c(monocyte = NA, lymphocyte = NA, endothelial = NA))
  x <- as.matrix(reference[ok, ])
  y <- y[ok]
  dmat <- crossprod(x)
  diag(dmat) <- diag(dmat) + 1e-8
  amat <- cbind(diag(ncol(x)), rep(-1, ncol(x)))
  solve.QP(dmat, crossprod(x, y), amat, c(0, 0, 0, -1))$solution
}

coeff <- t(apply(betas[, probes, drop = FALSE], 1, project))
colnames(coeff) <- c("monocyte", "lymphocyte", "endothelial")
fractions <- data.frame(
  sample_id = betas$sample_id,
  coeff,
  stromal = rowSums(coeff),
  unassigned = 1 - rowSums(coeff),
  stringsAsFactors = FALSE
)
fractions <- merge(fractions, manifest[, c("sample_id", "cohort", "group")], by = "sample_id", all.x = TRUE)
utils::write.csv(fractions, fraction_out, row.names = FALSE)

intensity <- utils::read.csv(file.path(work, "purity_estimates.csv"), stringsAsFactors = FALSE)
rna <- utils::read.csv(file.path(work, "rna_neoplastic_fraction.csv"), stringsAsFactors = FALSE)
joined <- merge(fractions, intensity[, c("sample_id", "purity", "purity_status")], by = "sample_id")
joined <- merge(joined, rna[, c("sample_id", "neoplastic")], by = "sample_id", all.x = TRUE)

spearman_ci <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  n <- length(x)
  if (n < 10) return(c(n = n, rho = NA, lo = NA, hi = NA))
  rho <- cor(x, y, method = "spearman")
  set.seed(1)
  draws <- replicate(2000, {
    idx <- sample.int(n, replace = TRUE)
    cor(x[idx], y[idx], method = "spearman")
  })
  ci <- stats::quantile(draws, c(0.025, 0.975), na.rm = TRUE)
  c(n = n, rho = rho, lo = unname(ci[1]), hi = unname(ci[2]))
}

rows <- list()
add <- function(cohort, contrast, stats) {
  rows[[length(rows) + 1]] <<- data.frame(
    cohort = cohort, contrast = contrast,
    n = unname(stats["n"]), rho = unname(stats["rho"]),
    ci_low = unname(stats["lo"]), ci_high = unname(stats["hi"]),
    stringsAsFactors = FALSE
  )
}
for (cohort in unique(joined$cohort)) {
  part <- joined[joined$cohort == cohort, ]
  estimated <- part$purity_status == "estimated"
  add(cohort, "stromal_vs_one_minus_intensity", spearman_ci(part$stromal[estimated], 1 - part$purity[estimated]))
  add(cohort, "unassigned_vs_intensity", spearman_ci(part$unassigned[estimated], part$purity[estimated]))
  add(cohort, "stromal_vs_one_minus_rna_neoplastic", spearman_ci(part$stromal, 1 - part$neoplastic))
  message(cohort,
    " stromal median estimable ", round(median(part$stromal[estimated], na.rm = TRUE), 3),
    " unestimable ", round(median(part$stromal[!estimated], na.rm = TRUE), 3))
}
result <- do.call(rbind, rows)
utils::write.csv(result, concordance_out, row.names = FALSE)
print(result)
message("wrote ", fraction_out)
