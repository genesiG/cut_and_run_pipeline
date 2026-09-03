#!/usr/bin/env Rscript

# =============================================================================
# compute_lfc_bins.R
#
# Genome-wide spike-in-normalized differential binding (1uM vs DMSO) over
# fixed-width tiled bins.  Produces a shrunken log2FC table per target for
# use in the delta-vs-delta (step_7j) correlation analysis.
#
# Normalization convention
# ------------------------
# spike-in SF (column "SF" in *_spikein_SF.txt):
#   SF = min(raw_SF) / raw_SF    (min-normalised; ranges from 0 < SF <= 1)
#   This is a DIVISIVE factor: samples with more spike-in get a smaller SF,
#   which shrinks their effective library size and therefore inflates their
#   CPM — correctly compensating for the lower input material.
#
# csaw/edgeR normFactors convention:
#   effective_lib_size = lib.size * norm.factors
#   edgeR *multiplies* lib.size by norm.factors when computing CPM.
#   Therefore to apply a spike-in SF directly we want:
#     normFactors(counts) <- 1 / sf_matched$SF
#   which is equivalent to bamCov_SF (already stored in the SF file).
#
# Dispersion handling for single-replicate targets (EZH2, GSTCBX7)
# ------------------------------------------------------------------
# We cannot estimate a per-target dispersion with n=1 per condition.
# Instead we "borrow" the BCV estimate from a replicated target (K27me2).
# If a pre-computed K27me2 dispersion RDS exists we load it; otherwise we
# set a sensible conservative common dispersion (bcv = 0.4 → disp = 0.16).
# =============================================================================

suppressPackageStartupMessages({
  library(reticulate)
  library(GenomicRanges)
  library(csaw)
  library(edgeR)
  library(dplyr)
})

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
args_cli <- commandArgs(trailingOnly = TRUE)

get_flag <- function(flag, args = args_cli, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0) return(default)
  if (idx[1] + 1 > length(args)) return(default)
  args[idx[1] + 1]
}

target       <- get_flag("--target",     default = "K27me2")
bin_size     <- as.integer(get_flag("--bin_size",   default = "10000"))
filter_min   <- as.numeric(get_flag("--filter_min", default = "1"))   # aveLogCPM threshold
regions_file <- get_flag("--regions",    default = "Analysis_Data/peaks/csaw/H3K27me2_DMSO.w150.d50.filt2.w2000.d500.filt1.lfc1.merge100.tmm.bed")
out_dir      <- get_flag("--out_dir",    default = "Analysis_Data/delta_vs_delta")

cat(sprintf("\n=== compute_lfc_bins.R ===\n"))
cat(sprintf("  Target       : %s\n", target))
cat(sprintf("  Bin size     : %d bp\n", bin_size))
cat(sprintf("  Filter min   : aveLogCPM > %.2f\n", filter_min))
cat(sprintf("  Regions BED  : %s\n", ifelse(is.null(regions_file), "None (genome-wide bins)", regions_file)))
cat(sprintf("  Output dir   : %s\n\n", out_dir))

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Load config via reticulate
# ---------------------------------------------------------------------------
use_python(Sys.which("python"), required = TRUE)
py_run_file("Scripts/config.py")

bam_dir       <- py$PROCESSEDBAMDIR
is_paired_end <- py$IS_PAIRED_END
mapq          <- as.integer(py$MAPQ)
spikein_dir   <- file.path(py$SCALINGDIR, "spikein")

# ---------------------------------------------------------------------------
# Targets with only 1 replicate per condition — must deduplicate
# ---------------------------------------------------------------------------
SINGLE_REP_TARGETS <- c("EZH2", "GSTCBX7")
has_replicates <- !(target %in% SINGLE_REP_TARGETS)

# ---------------------------------------------------------------------------
# Load metadata and select BAM files for this target + both conditions
# ---------------------------------------------------------------------------
meta_file <- file.path(py$METADATA,
                       paste0("sample_metadata_", target, "_processed.txt"))
if (!file.exists(meta_file)) {
  stop("Metadata file not found: ", meta_file)
}
meta <- read.table(meta_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat("Metadata loaded:\n"); print(meta); cat("\n")

# Filter to the exact target antibody and the two treatment groups
# This is critical for combined metadata files (e.g. CBX2_processed.txt
# contains both ANTIBODY="CBX2" and ANTIBODY="GSTCBX2" rows).
meta <- meta[meta$ANTIBODY == target & meta$GROUP %in% c("DMSO", "1uM"), ]
meta$GROUP <- factor(meta$GROUP, levels = c("DMSO", "1uM"))
meta <- meta[order(meta$GROUP), ]   # DMSO first

if (nrow(meta) == 0) {
  stop("No rows in metadata for ANTIBODY == '", target, "' and GROUP in {DMSO, 1uM}")
}

bam_suffix <- ".qc.sort.markdup.mapq30.final.bam"
bam_paths  <- file.path(bam_dir, meta$ID)
missing    <- bam_paths[!file.exists(bam_paths)]
if (length(missing) > 0) stop("Missing BAM files:\n", paste(missing, collapse = "\n"))

cat(sprintf("Using %d BAM file(s) for target %s:\n", length(bam_paths), target))
cat(paste0("  ", basename(bam_paths), "\n"), sep = "")

# ---------------------------------------------------------------------------
# csaw readParam
# Targets without replicates: set dedup = TRUE to reduce noise
# ---------------------------------------------------------------------------
pe_str <- if (is_paired_end) "both" else "none"

if (!has_replicates) {
  cat("\nNOTE: Single-replicate target detected — setting dedup = TRUE.\n")
  param <- readParam(minq = mapq, pe = pe_str, dedup = TRUE)
} else {
  param <- readParam(minq = mapq, pe = pe_str, dedup = FALSE)
}

# ---------------------------------------------------------------------------
# 1. Genome-wide bin counting or region counting via csaw
# ---------------------------------------------------------------------------
if (!is.null(regions_file) && nzchar(regions_file) && regions_file != "None" && regions_file != "NULL") {
  if (!file.exists(regions_file)) stop("Specified regions BED file not found: ", regions_file)
  cat(sprintf("\nCounting reads over fixed peakset regions from %s using csaw::regionCounts...\n", basename(regions_file)))
  bed_df <- read.table(regions_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  gr_regions <- GenomicRanges::GRanges(
    seqnames = bed_df$V1,
    ranges   = IRanges::IRanges(start = bed_df$V2 + 1, end = bed_df$V3)
  )
  counts <- csaw::regionCounts(bam_paths, regions = gr_regions, param = param)
  cat(sprintf("  Total regions counted: %d\n", nrow(counts)))
  
  # ---------------------------------------------------------------------------
  # 2. Filtering (Light coverage floor across fixed peakset regions)
  # ---------------------------------------------------------------------------
  floor_val <- if (filter_min == 1.0) 0.0 else filter_min
  cat(sprintf("\nApplying light coverage floor (aveLogCPM > %.2f) to fixed peakset regions...\n", floor_val))
  alc       <- aveLogCPM(asDGEList(counts))
  keep      <- alc > floor_val
  counts_f  <- counts[keep, ]
  cat(sprintf("  Regions before floor: %d\n  Regions after floor : %d (%.1f%% retained)\n",
              nrow(counts), nrow(counts_f), 100 * nrow(counts_f) / nrow(counts)))
  alc       <- aveLogCPM(asDGEList(counts_f))
} else {
  cat(sprintf("\nCounting reads into %d-bp bins genome-wide via csaw::windowCounts...\n", bin_size))
  counts <- windowCounts(bam_paths, bin = TRUE, width = bin_size, param = param)
  cat(sprintf("  Total bins: %d\n", nrow(counts)))
  
  # ---------------------------------------------------------------------------
  # 2. Filter low-abundance bins by aveLogCPM
  # ---------------------------------------------------------------------------
  cat(sprintf("\nFiltering bins with aveLogCPM > %.2f ...\n", filter_min))
  alc       <- aveLogCPM(asDGEList(counts))
  keep      <- alc > filter_min
  counts_f  <- counts[keep, ]
  cat(sprintf("  Bins before filter: %d\n  Bins after  filter: %d\n",
              nrow(counts), nrow(counts_f)))
}

# ---------------------------------------------------------------------------
# 3. Load spike-in scaling factors
# ---------------------------------------------------------------------------
sf_file <- file.path(spikein_dir, paste0(target, "_spikein_SF.txt"))
if (!file.exists(sf_file)) {
  stop("Spike-in SF file not found: ", sf_file)
}
cat("\nLoading spike-in SFs from:", sf_file, "\n")
sf_df <- read.table(sf_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
print(sf_df)

# csaw stores BAM paths in counts$bam.files (not in colnames, which may be NULL).
# Match SFs to the exact BAM file order that windowCounts used.
bam_files_used <- basename(counts_f$bam.files)   # basenames in csaw column order
cat("\ncsaw BAM order (bam.files):\n")
cat(paste0("  ", bam_files_used, "\n"), sep = "")

sf_matched <- sf_df[match(bam_files_used, sf_df$ID), ]

if (any(is.na(sf_matched$SF))) {
  stop("SF matching failed. BAMs not found in SF file:\n",
       paste(bam_files_used[is.na(sf_matched$SF)], collapse = "\n"))
}

# Verify exact 1-to-1 alignment
stopifnot(identical(bam_files_used, sf_matched$ID))

# csaw normFactors convention: normFactors are DIVISIVE.
# SF (min-normalised) is already the right quantity for division.
# Equivalently: normFactors(counts_f) <- 1 / sf_matched$bamCov_SF
# which equals sf_matched$SF by construction (bamCov_SF = 1/SF).
cat("\nAssigning spike-in normFactors (= SF, divisive)...\n")
cat("  SF values: ", paste(round(sf_matched$SF, 4), collapse = ", "), "\n")
counts_f$norm.factors <- sf_matched$SF

# ---------------------------------------------------------------------------
# 4. Convert to DGEList and fit GLM
# ---------------------------------------------------------------------------
y <- asDGEList(counts_f)
y$samples$norm.factors <- sf_matched$SF

n_dmso <- sum(meta$GROUP == "DMSO")
n_1um  <- sum(meta$GROUP == "1uM")
group  <- factor(c(rep("DMSO", n_dmso), rep("1uM", n_1um)), levels = c("DMSO", "1uM"))
y$samples$group <- group
design <- model.matrix(~ group)

cat(sprintf("\nDesign matrix: %d DMSO + %d 1uM replicates\n", n_dmso, n_1um))

if (has_replicates) {
  # Standard path: estimate dispersions from data and use quasi-likelihood F-test
  cat("Estimating dispersions from data (replicated target)...\n")
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  res <- glmQLFTest(fit, coef = 2)
} else {
  # Single-replicate path: no residual df available for quasi-likelihood dispersion estimation.
  # Use exact negative binomial GLM fit with fixed 0.05 dispersion as requested to remove ordering dependencies across jobs.
  disp_val <- 0.05
  cat(sprintf("Single-replicate target detected: setting common.dispersion = %.4f (no across-job borrowing)\n", disp_val))
  y$common.dispersion <- disp_val
  
  fit <- glmFit(y, design, dispersion = y$common.dispersion)
  res <- glmLRT(fit, coef = 2)
}

# ---------------------------------------------------------------------------
# 5. Shrinkage via edgeR::predFC (prior.count = 3)
# ---------------------------------------------------------------------------
cat("\nShrinking log2FC with predFC (prior.count = 3)...\n")
if (has_replicates) {
  lfc_shrunk_mat <- predFC(y, design, prior.count = 3)
} else {
  lfc_shrunk_mat <- predFC(y, design, prior.count = 3, dispersion = y$common.dispersion)
}
# predFC returns a matrix of shrunken coefficients; column 2 = group1uM contrast
lfc_shrunk <- lfc_shrunk_mat[, 2]

# ---------------------------------------------------------------------------
# 6. Assemble output table
# ---------------------------------------------------------------------------
tt       <- res$table
regions  <- rowRanges(counts_f)

out_df <- data.frame(
  seqnames   = as.character(seqnames(regions)),
  start      = start(regions),
  end        = end(regions),
  bin_id     = paste0(as.character(seqnames(regions)), ":",
                      start(regions), "-", end(regions)),
  logFC_raw  = tt$logFC,
  logFC_shrunk = lfc_shrunk,
  aveLogCPM  = aveLogCPM(y),
  PValue     = tt$PValue,
  FDR        = p.adjust(tt$PValue, method = "BH"),
  stringsAsFactors = FALSE
)

# Estimate logFC SE from GLM fit (approximate: SE ≈ logFC / t-stat)
stat_val  <- if (!is.null(res$table$F)) sqrt(res$table$F) else sqrt(res$table$LR)
out_df$lfcSE <- abs(out_df$logFC_raw) / pmax(stat_val, 0.001)

# Re-order columns
out_df <- out_df[, c("seqnames", "start", "end", "bin_id",
                      "logFC_raw", "logFC_shrunk", "lfcSE",
                      "aveLogCPM", "PValue", "FDR")]

# ---------------------------------------------------------------------------
# 7. Save outputs
# ---------------------------------------------------------------------------
target_safe  <- tolower(gsub("-", "_", target))
out_prefix   <- file.path(out_dir, sprintf("lfc_bins_%s_%dbp", target_safe, bin_size))

# RDS (fast binary, used by step_7j)
rds_path <- paste0(out_prefix, ".rds")
saveRDS(out_df, rds_path)
cat("Saved RDS:", rds_path, "\n")

# CSV (human-readable)
csv_path <- paste0(out_prefix, ".csv")
write.csv(out_df, csv_path, row.names = FALSE)
cat("Saved CSV:", csv_path, "\n")

# Quick summary
cat(sprintf("\nSummary for %s (%d-bp bins, after filter):\n", target, bin_size))
cat(sprintf("  N bins       : %d\n", nrow(out_df)))
cat(sprintf("  logFC_shrunk : min=%.3f  median=%.3f  max=%.3f\n",
            min(out_df$logFC_shrunk), median(out_df$logFC_shrunk), max(out_df$logFC_shrunk)))
cat(sprintf("  Bins FDR<0.1 : %d\n", sum(out_df$FDR < 0.1, na.rm = TRUE)))

# ---------------------------------------------------------------------------
# 8. Filter_min Sensitivity Evaluation & Diagnostic Visuals (only when tiling bins)
# ---------------------------------------------------------------------------
if (is.null(regions_file) || !nzchar(regions_file) || regions_file == "None" || regions_file == "NULL") {
  cat("\nComputing filter_min sensitivity metrics across candidate thresholds...\n")

  # Evaluate grid of candidate thresholds
  tau_grid <- sort(unique(c(seq(-1.0, 3.0, by = 0.25), filter_min)))
  total_bins <- length(alc)
  raw_counts_mat <- assay(counts)

  sens_list <- lapply(tau_grid, function(tau) {
    k_tau <- alc > tau
    n_ret <- sum(k_tau)
    n_dis <- total_bins - n_ret
    
    med_alc <- if (n_ret > 0) median(alc[k_tau]) else NA_real_
    mean_cnt <- if (n_ret > 0) mean(raw_counts_mat[k_tau, ]) else NA_real_
    
    data.frame(
      filter_min_threshold = tau,
      bins_retained        = n_ret,
      pct_retained         = round(100 * n_ret / total_bins, 2),
      bins_discarded       = n_dis,
      pct_discarded        = round(100 * n_dis / total_bins, 2),
      median_aveLogCPM     = round(med_alc, 3),
      mean_raw_count       = round(mean_cnt, 2),
      is_current_choice    = (abs(tau - filter_min) < 1e-5),
      stringsAsFactors     = FALSE
    )
  })
  sens_df <- do.call(rbind, sens_list)

  # Save sensitivity metrics table
  sens_tsv <- file.path(out_dir, sprintf("lfc_bins_filter_metrics_%s_%dbp.tsv", target_safe, bin_size))
  write.table(sens_df, sens_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
  cat("Saved filter_min sensitivity metrics TSV:", sens_tsv, "\n")

  # Print clean summary table to log for common cutoffs
  cat("\n--- Filter Threshold Sensitivity Summary ---\n")
  common_taus <- c(-0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5)
  print(sens_df[sens_df$filter_min_threshold %in% common_taus | sens_df$is_current_choice, ], row.names = FALSE)
  cat("--------------------------------------------\n")

  # Helper to generate 4-panel diagnostic plot (used for both PDF and PNG)
  generate_filter_plots <- function() {
    old_par <- par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 1), oma = c(0, 0, 2.5, 0), cex.main = 1.1)
    on.exit(par(old_par))
    
    # Panel 1: aveLogCPM Density across ALL bins genome-wide
    d_alc <- density(alc, na.rm = TRUE)
    plot(d_alc, main = "Genome-Wide aveLogCPM Density", xlab = "aveLogCPM (all bins)",
         ylab = "Density", col = "#1f77b4", lwd = 2.5, type = "l")
    x_shade <- c(d_alc$x[d_alc$x <= filter_min], filter_min, min(d_alc$x))
    y_shade <- c(d_alc$y[d_alc$x <= filter_min], 0, 0)
    polygon(x_shade, y_shade, col = rgb(0.8, 0.2, 0.2, 0.3), border = NA)
    abline(v = filter_min, col = "#d62728", lwd = 2, lty = 2)
    abline(v = c(0.0, 0.5, 2.0), col = c("#7f7f7f", "#2ca02c", "#9467bd"), lty = 3, lwd = 1.5)
    legend("topright", legend = c(sprintf("Current: %.2f", filter_min), "0.0", "0.5", "2.0"),
           col = c("#d62728", "#7f7f7f", "#2ca02c", "#9467bd"), lty = c(2,3,3,3), lwd = 2, cex = 0.8, bty = "n")
    
    # Panel 2: Tradeoff Curve (% Discarded vs % Retained)
    plot(sens_df$filter_min_threshold, sens_df$pct_discarded, type = "o", pch = 16, col = "#d62728",
         lwd = 2, ylim = c(0, 100), xlab = "filter_min Threshold (aveLogCPM)", ylab = "% of Genome-Wide Bins",
         main = "Bins Discarded vs Retained")
    lines(sens_df$filter_min_threshold, sens_df$pct_retained, type = "o", pch = 16, col = "#1f77b4", lwd = 2)
    abline(v = filter_min, col = "#333333", lty = 2, lwd = 1.8)
    curr_row <- sens_df[sens_df$is_current_choice, ][1, ]
    points(filter_min, curr_row$pct_discarded, pch = 19, col = "#d62728", cex = 1.5)
    points(filter_min, curr_row$pct_retained, pch = 19, col = "#1f77b4", cex = 1.5)
    legend("center", legend = c("% Discarded", "% Retained", sprintf("Current (%.1f%% ret)", curr_row$pct_retained)),
           col = c("#d62728", "#1f77b4", "#333333"), lty = c(1, 1, 2), pch = c(16, 16, NA), lwd = 2, cex = 0.8, bty = "n")
    
    # Panel 3: Mean Raw Read Count vs filter_min Threshold
    plot(sens_df$filter_min_threshold, sens_df$mean_raw_count, type = "o", pch = 16, col = "#2ca02c",
         lwd = 2, xlab = "filter_min Threshold (aveLogCPM)", ylab = "Mean Raw Read Count / Sample",
         main = "Signal Quality vs Filtering Threshold")
    abline(v = filter_min, col = "#333333", lty = 2, lwd = 1.8)
    points(filter_min, curr_row$mean_raw_count, pch = 19, col = "#2ca02c", cex = 1.5)
    text(filter_min, curr_row$mean_raw_count, pos = 4, cex = 0.85, font = 2,
         labels = sprintf("%.1f reads", curr_row$mean_raw_count))
    
    # Panel 4: MA-Like Plot of Retained Bins Highlighted by Significance
    n_plot <- min(nrow(out_df), 15000)
    idx_p <- if (nrow(out_df) > n_plot) sample(nrow(out_df), n_plot) else seq_len(nrow(out_df))
    sub_df <- out_df[idx_p, ]
    sig_idx <- which(!is.na(sub_df$FDR) & sub_df$FDR < 0.1)
    non_idx <- setdiff(seq_len(nrow(sub_df)), sig_idx)
    plot(sub_df$aveLogCPM[non_idx], sub_df$logFC_shrunk[non_idx], pch = 16, col = rgb(0.5, 0.5, 0.5, 0.25),
         cex = 0.6, xlab = "aveLogCPM (Retained Bins)", ylab = "Shrunken log2FC",
         main = "MA Plot of Retained Bins (FDR < 0.1)")
    if (length(sig_idx) > 0) {
      points(sub_df$aveLogCPM[sig_idx], sub_df$logFC_shrunk[sig_idx], pch = 16, col = rgb(0.9, 0.2, 0.1, 0.7), cex = 0.8)
    }
    abline(h = 0, col = "#333333", lty = 3)
    abline(v = filter_min, col = "#d62728", lty = 2, lwd = 2)
    legend("topleft", legend = c(sprintf("Significant (n=%d)", sum(!is.na(out_df$FDR) & out_df$FDR < 0.1)),
                                 sprintf("Non-sig (n=%d)", sum(!is.na(out_df$FDR) & out_df$FDR >= 0.1))),
           col = c("#e63946", "#808080"), pch = 16, cex = 0.8, bty = "n")
    mtext(sprintf("LFC Bin Filtering & Threshold Diagnostics: %s (%d bp)", target, bin_size),
          outer = TRUE, cex = 1.3, font = 2)
  }

  # Save PDF plot
  pdf_path <- file.path(out_dir, sprintf("lfc_bins_filter_diagnostics_%s_%dbp.pdf", target_safe, bin_size))
  pdf(pdf_path, width = 10, height = 8)
  generate_filter_plots()
  dev.off()
  cat("Saved filter diagnostics PDF:", pdf_path, "\n")

  # Save PNG plot
  png_path <- file.path(out_dir, sprintf("lfc_bins_filter_diagnostics_%s_%dbp.png", target_safe, bin_size))
  png(png_path, width = 1200, height = 960, res = 130)
  generate_filter_plots()
  dev.off()
  cat("Saved filter diagnostics PNG:", png_path, "\n")
} else {
  cat("\nCounting over fixed pre-specified peak regions: skipping filter sensitivity grid plots.\n")
}

cat("\n=== compute_lfc_bins.R complete ===\n")
