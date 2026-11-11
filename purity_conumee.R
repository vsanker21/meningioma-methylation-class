# Purity from total array intensity.
#
# conumee fits each tumour to the sex-matched DKFZ normal panel and segments
# the autosomal log2 intensity ratio. Those normals are a reference panel, not
# the methylation classifier. Sex chromosomes are left out of the fit so X/Y
# intensity does not rescale the autosomes.
#
# A single-copy loss in a tumour of purity rho has log2 ratio
# log2(1 - rho/2) once the diploid peak is at zero. The diploid peak is the
# mode of the bin ratios, not the median: a median is pulled down when many
# arms are lost. rho = 2 * (1 - 2^L).
#
# The event is 22q when that arm carries a loss covering at least one third of
# its length. 22q is the founder NF2 event, so it is the single-copy ruler when
# it qualifies. Otherwise the ruler is the deepest remaining autosomal loss
# with the same coverage rule. A subclonal loss cannot be deeper than the
# clonal one, and the shallowest arm is where a small reference shelf appears
# (chromosome 6 in a copy-quiet tumour). A shift deeper than a pure
# single-copy loss is not used. No such arm means purity is unestimable.

suppressPackageStartupMessages({
  library(minfi)
  library(conumee)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: purity_conumee.R manifest.csv output.csv")
manifest_path <- args[[1]]
out_path <- args[[2]]

ref_dir <- "F:/Meningioma Genetics/data/raw/github/mnp_training/CNV_data"
message("loading normal panel")
load(file.path(ref_dir, "CNanalysis4_conumee_ANNO.vh20150715.RData"))
load(file.path(ref_dir, "CNanalysis4_conumee_REF-M.vh20150715.RData"))
load(file.path(ref_dir, "CNanalysis4_conumee_REF-F.vh20150715.RData"))

anno_auto <- anno
skip_parm <- c("13", "14", "15", "21", "22")

chr_key <- function(x) {
  x <- as.character(x)
  sub("^chr", "", x, ignore.case = TRUE)
}

arm_table <- function(anno_obj) {
  # anno@genome columns are chr, size, pq. pq is the centromere coordinate.
  g <- anno_obj@genome
  data.frame(
    chr = chr_key(g$chr),
    centro_start = g$pq,
    centro_end = g$pq,
    chr_start = 1,
    chr_end = g$size,
    stringsAsFactors = FALSE
  )
}

arms <- arm_table(anno_auto)
message("genome columns: ", paste(colnames(anno_auto@genome), collapse = ","))
message("arms: ", nrow(arms))

prepare_anno <- function(anno_obj, probe_ids) {
  keep <- names(anno_obj@probes) %in% probe_ids
  dropped <- sum(!keep)
  if (dropped) {
    anno_obj@probes <- anno_obj@probes[keep]
    message("dropped ", dropped, " annotation probes absent from this array")
  }
  anno_obj
}

purity_from_bins <- function(object) {
  bins <- object@anno@bins
  ratio <- object@bin$ratio[names(bins)]
  chr <- chr_key(seqnames(bins))
  mid <- bins$midpoint
  ok <- is.finite(ratio)
  mode <- density(ratio[ok], n = 1024, bw = "SJ")
  diploid <- mode$x[which.max(mode$y)]
  resid <- ratio - diploid
  noise <- mad(resid[ok & abs(resid) < 0.05], na.rm = TRUE)
  if (!is.finite(noise) || noise <= 0) noise <- mad(resid[ok], na.rm = TRUE)
  threshold <- max(0.08, 3 * noise)
  single_copy_floor <- log2(0.5) - 0.1

  seg <- object@seg$summary
  seg$chr <- chr_key(seg$chrom)
  calls <- list()
  for (i in seq_len(nrow(arms))) {
    info <- arms[i, ]
    for (arm in c("p", "q")) {
      if (arm == "p" && info$chr %in% skip_parm) next
      if (arm == "p") {
        a0 <- info$chr_start
        a1 <- info$centro_start
        in_arm <- chr == info$chr & mid < info$centro_start
      } else {
        a0 <- info$centro_end
        a1 <- info$chr_end
        in_arm <- chr == info$chr & mid > info$centro_end
      }
      arm_bp <- a1 - a0
      if (!is.finite(arm_bp) || arm_bp <= 0 || sum(in_arm, na.rm = TRUE) < 20) next
      on_chr <- seg$chr == info$chr
      covered <- 0
      weighted <- 0
      for (j in which(on_chr)) {
        overlap <- min(seg$loc.end[j], a1) - max(seg$loc.start[j], a0)
        if (!is.finite(overlap) || overlap <= 0) next
        level <- seg$seg.mean[j] - diploid
        if (level < -threshold) {
          covered <- covered + overlap
          weighted <- weighted + level * overlap
        }
      }
      frac <- covered / arm_bp
      if (frac < 1/3) next
      L <- weighted / covered
      bin_resid <- resid[in_arm & is.finite(resid)]
      calls[[length(calls) + 1]] <- data.frame(
        arm = paste0(info$chr, arm),
        L = L,
        coverage = frac,
        mad = mad(bin_resid),
        deeper = L < single_copy_floor,
        stringsAsFactors = FALSE
      )
    }
  }
  calls <- if (length(calls)) do.call(rbind, calls) else data.frame()
  usable <- if (nrow(calls)) calls[!calls$deeper, ] else calls
  status <- "unestimable_no_arm_loss"
  event_arm <- NA_character_
  event_L <- NA_real_
  event_cov <- NA_real_
  purity <- NA_real_
  deepest_arm <- NA_character_
  deepest_L <- NA_real_
  if (nrow(usable)) {
    deepest <- usable[which.min(usable$L), ]
    deepest_arm <- deepest$arm[1]
    deepest_L <- deepest$L[1]
    pick <- if ("22q" %in% usable$arm) usable[usable$arm == "22q", ] else deepest
    event_arm <- pick$arm[1]
    event_L <- pick$L[1]
    event_cov <- pick$coverage[1]
    purity <- 2 * (1 - 2^event_L)
    status <- if (is.finite(purity) && purity > 0 && purity <= 1.05) "estimated" else "unestimable_out_of_range"
  } else if (nrow(calls) && all(calls$deeper)) {
    status <- "unestimable_deeper_than_single_copy"
  }
  q22 <- if (nrow(calls)) calls[calls$arm == "22q", ] else calls
  data.frame(
    diploid_mode = diploid,
    bin_shift = unname(object@bin$shift),
    fit_noise = unname(object@fit$noise),
    noise_mad = noise,
    loss_threshold = threshold,
    event_arm = event_arm,
    event_L = event_L,
    event_coverage = event_cov,
    purity = purity,
    purity_status = status,
    l_22q = if (nrow(q22)) q22$L[1] else NA_real_,
    cov_22q = if (nrow(q22)) q22$coverage[1] else NA_real_,
    deepest_arm = deepest_arm,
    deepest_L = deepest_L,
    n_loss_arms = if (nrow(calls)) nrow(calls) else 0,
    loss_arms = if (nrow(calls)) paste0(calls$arm, ":", sprintf("%.3f", calls$L), collapse = ";") else "",
    stringsAsFactors = FALSE
  )
}

# minfi::getSex writes predictedSex from this cutoff. Its extra kmeans check
# needs at least two samples, so a per-tumour call uses the cutoff alone.
predict_sex <- function(mset) {
  mset@annotation <- c(array = "IlluminaHumanMethylationEPIC", annotation = "ilm10b2.hg19")
  gset <- mapToGenome(mset)
  cn <- getCN(gset)
  chr <- as.character(seqnames(gset))
  x_med <- median(cn[chr == "chrX", 1], na.rm = TRUE)
  y_med <- median(cn[chr == "chrY", 1], na.rm = TRUE)
  sex <- if (is.finite(x_med) && is.finite(y_med) && (y_med - x_med) < -2) "F" else "M"
  list(sex = sex, x_med = x_med, y_med = y_med)
}

arm_median <- function(ratio, bins) {
  chr <- chr_key(seqnames(bins))
  mid <- bins$midpoint
  rows <- list()
  for (i in seq_len(nrow(arms))) {
    info <- arms[i, ]
    for (arm in c("p", "q")) {
      if (arm == "p" && info$chr %in% skip_parm) next
      in_arm <- if (arm == "p") chr == info$chr & mid < info$centro_start else chr == info$chr & mid > info$centro_end
      rows[[length(rows) + 1]] <- data.frame(
        arm = paste0(info$chr, arm),
        median_ratio = median(ratio[in_arm], na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

bin_median_ratio <- function(anno_obj, probe_ratio) {
  hits <- as.matrix(findOverlaps(query = anno_obj@bins, subject = anno_obj@probes))
  grouped <- split(probe_ratio[names(anno_obj@probes)[hits[, "subjectHits"]]],
                   names(anno_obj@bins)[hits[, "queryHits"]])
  sapply(grouped, median, na.rm = TRUE)[names(anno_obj@bins)]
}

one_sample <- function(row) {
  rg <- read.metharray(row$basename, verbose = FALSE, force = TRUE)
  rg@annotation <- c(array = "IlluminaHumanMethylationEPIC", annotation = "ilm10b2.hg19")
  mset <- preprocessRaw(rg)
  sex_info <- predict_sex(mset)
  sex <- sex_info$sex
  ref <- if (identical(sex, "F")) refF.data else refM.data
  ref_name <- if (identical(sex, "F")) "REF-F" else "REF-M"
  loaded <- CNV.load(mset)
  use_anno <- prepare_anno(anno_auto, intersect(rownames(loaded@intensity), rownames(ref@intensity)))
  fit <- CNV.fit(loaded, ref, use_anno)
  fit <- CNV.bin(fit)
  fit <- CNV.segment(fit)
  out <- purity_from_bins(fit)
  probes <- names(use_anno@probes)
  ref_med <- matrixStats::rowMedians(as.matrix(ref@intensity[probes, , drop = FALSE]))
  query_int <- loaded@intensity[probes, 1]
  simple <- log2(pmax(query_int, 1) / pmax(ref_med, 1))
  names(simple) <- probes
  simple_bins <- bin_median_ratio(use_anno, simple)
  simple_ok <- is.finite(simple_bins)
  simple_mode <- density(simple_bins[simple_ok], n = 1024, bw = "SJ")
  simple_diploid <- simple_mode$x[which.max(simple_mode$y)]
  info22 <- arms[arms$chr == "22", ]
  in22 <- chr_key(seqnames(use_anno@bins)) == "22" & use_anno@bins$midpoint > info22$centro_end
  simple_l_22q <- median(simple_bins[in22], na.rm = TRUE) - simple_diploid
  out$simple_diploid_mode <- simple_diploid
  out$simple_l_22q <- simple_l_22q
  out$simple_purity_22q <- 2 * (1 - 2^simple_l_22q)
  out$sample_id <- row$sample_id
  out$cohort <- row$cohort
  out$sex_predicted <- sex
  out$x_med <- sex_info$x_med
  out$y_med <- sex_info$y_med
  out$sex_clinical <- row$sex_clinical
  out$ref_used <- ref_name
  out$group <- row$group
  arms_out <- arm_median(fit@bin$ratio - out$diploid_mode, fit@anno@bins)
  arms_out$sample_id <- row$sample_id
  attr(out, "arms") <- arms_out
  out
}

manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
result_cols <- c(
  "sample_id", "cohort", "group", "sex_clinical", "sex_predicted",
  "x_med", "y_med", "ref_used", "purity", "purity_status",
  "event_arm", "event_L", "event_coverage", "l_22q", "cov_22q",
  "deepest_arm", "deepest_L", "diploid_mode", "bin_shift", "fit_noise",
  "noise_mad", "loss_threshold", "n_loss_arms", "loss_arms",
  "simple_diploid_mode", "simple_l_22q", "simple_purity_22q"
)
arms_path <- sub("\\.csv$", ".arms.csv", out_path)
done <- if (file.exists(out_path)) utils::read.csv(out_path, stringsAsFactors = FALSE)$sample_id else character()
header <- !file.exists(out_path)
arms_header <- !file.exists(arms_path)

blank_result <- function(sample_id, cohort, status) {
  row <- as.list(rep(NA, length(result_cols)))
  names(row) <- result_cols
  row$sample_id <- sample_id
  row$cohort <- cohort
  row$purity_status <- status
  as.data.frame(row, stringsAsFactors = FALSE)
}

append_csv <- function(path, row, write_header) {
  utils::write.table(row, path, sep = ",", row.names = FALSE, col.names = write_header,
                     append = !write_header, quote = TRUE)
}

for (i in seq_len(nrow(manifest))) {
  id <- manifest$sample_id[i]
  if (id %in% done) next
  message("sample ", i, "/", nrow(manifest), " ", id)
  result <- tryCatch(one_sample(manifest[i, ]), error = function(e) {
    message("FAILED ", id, ": ", conditionMessage(e))
    blank_result(id, manifest$cohort[i], paste0("failed: ", conditionMessage(e)))
  })
  arms_out <- attr(result, "arms")
  result <- result[intersect(result_cols, names(result))]
  missing <- setdiff(result_cols, names(result))
  for (nm in missing) result[[nm]] <- NA
  append_csv(out_path, result[result_cols], header)
  header <- FALSE
  if (!is.null(arms_out)) {
    append_csv(arms_path, arms_out, arms_header)
    arms_header <- FALSE
  }
}

out <- utils::read.csv(out_path, stringsAsFactors = FALSE)
message("status counts")
print(table(out$cohort, out$purity_status, useNA = "ifany"))
estimated <- out$purity_status == "estimated"
message("estimated ", sum(estimated, na.rm = TRUE), " of ", nrow(out))
if (any(estimated, na.rm = TRUE)) print(summary(out$purity[estimated]))
message("wrote ", out_path)
