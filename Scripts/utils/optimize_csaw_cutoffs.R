#!/usr/bin/env Rscript

# =============================================================================
# optimize_csaw_cutoffs.R  –  Empirical cut-off optimizer for csaw filtering
#
# Performs a parallelised grid search over (c_small, c_large) enrichment
# thresholds and identifies the combination that maximises the number of
# significant differentially-bound regions (DBRs) at a given FDR target.
#
# CLI arguments (same interface as callpeaks_csaw.R):
#   --chip_bams    Space-separated paths to ChIP/CUT&TAG BAM files
#   --ctrl_bams    Space-separated paths to control (IgG) BAM files
#   --sample_id    Label used for naming output files
#   --workdir      Project root (WORKDIR from config.py)
#   --codedir      Scripts directory (CODEDIR from config.py)
#   --norm_method  CHIPSEQSPIKEINFREE | SPIKEIN | TMM  (default: TMM)
#
# Outputs (all under Analysis_Data/csaw/cutoff_analysis/{sample_clean}/):
#   {sample_clean}_cutoff_results.csv         — full grid results table
#   {sample_clean}_pareto_heatmap.pdf         — 2-D tile Pareto plot
#   {sample_clean}_opt_filter_histograms.pdf  — filter histogram @ optimal cutoff
#   {sample_clean}_opt_bcv.pdf                — BCV plot @ optimal cutoff
#   {sample_clean}_opt.bed                    — peaks @ optimal cutoff only
# =============================================================================

suppressPackageStartupMessages({
  options(repos = c(CRAN = "https://cran.rstudio.com"))
  if (!require("reticulate",   quietly = TRUE)) install.packages("reticulate")
  if (!require("BiocParallel", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install("BiocParallel")
  }
  library(reticulate)
  library(BiocParallel)
})

# ---------------------------------------------------------------------------
# 1. Parse CLI arguments
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

get_flag_values <- function(flag, args, required = TRUE) {
  idx <- which(args == flag)
  if (length(idx) == 0) {
    if (required) stop(paste("Missing required argument:", flag))
    return(character(0))
  }
  start <- idx[1] + 1
  if (start > length(args)) stop(paste("No value supplied for:", flag))
  end <- start
  while (end <= length(args) && !startsWith(args[end], "--")) end <- end + 1
  args[start:(end - 1)]
}
get_flag_value <- function(flag, args, required = TRUE) {
  vals <- get_flag_values(flag, args, required)
  if (length(vals) == 0) return(NULL)
  vals[1]
}

chip_bams   <- get_flag_values("--chip_bams",  args)
ctrl_bams   <- get_flag_values("--ctrl_bams",  args)
sample_id   <- get_flag_value( "--sample_id",  args)
workdir     <- get_flag_value( "--workdir",    args)
codedir     <- get_flag_value( "--codedir",    args)
norm_method <- get_flag_value( "--norm_method", args, required = FALSE)
if (is.null(norm_method)) norm_method <- "TMM"

# "Cleaned" sample name for output files (strip trailing underscores/spaces)
sample_clean <- gsub("[_\\s]+$", "", sample_id)

cat("=== optimize_csaw_cutoffs.R ===\n")
cat("  sample_id    :", sample_id, "\n")
cat("  sample_clean :", sample_clean, "\n")
cat("  norm_method  :", norm_method, "\n")
cat("  ChIP BAMs    :", length(chip_bams), "\n")
cat("  Control BAMs :", length(ctrl_bams), "\n\n")

# ---------------------------------------------------------------------------
# 2. Bootstrap environment
# ---------------------------------------------------------------------------
setwd(workdir)
use_python(Sys.which("python"), required = TRUE)
py_run_file(file.path(codedir, "config.py"))

importabledir <- file.path(workdir, "Importable_Data")
analysisdir   <- file.path(workdir, "Analysis_Data")
peakdir       <- if (!is.null(py$PEAKDIR)) py$PEAKDIR else file.path(analysisdir, "peaks")
csawdir       <- if (!is.null(py$CSAW_OUTDIR)) py$CSAW_OUTDIR else file.path(peakdir, "csaw")
rdsdir        <- file.path(importabledir, "rds")
csawrdsdir    <- file.path(rdsdir, "csaw")
normdir       <- file.path(analysisdir, "normalization")
outdir        <- file.path(csawdir, "cutoff_analysis", sample_clean)
beddir        <- csawdir
for (d in c(csawrdsdir, outdir, beddir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# 3. Load required packages (sourced from csaw_parameters.R)
# ---------------------------------------------------------------------------
source(file.path(codedir, "utils", "csaw_parameters.R"))

load_pkg <- function(pkg, bioc = FALSE) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    if (bioc) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE, quietly = TRUE)
  }
}
for (pkg in packages)    load_pkg(pkg, bioc = FALSE)
for (pkg in biopackages) load_pkg(pkg, bioc = TRUE)
load_pkg("ggplot2",      bioc = FALSE)
load_pkg("ggrepel",      bioc = FALSE)

# ---------------------------------------------------------------------------
# 4. Validate BAMs
# ---------------------------------------------------------------------------
bam.files <- c(chip_bams, ctrl_bams)
for (f in bam.files) if (!file.exists(f)) stop("BAM file not found: ", f)
n_chip <- length(chip_bams)
n_ctrl <- length(ctrl_bams)

# ---------------------------------------------------------------------------
# 5. Load or generate RDS objects
#    Naming convention (mirrors callpeaks_csaw.R):
#      {sample_clean}.w{width}.d{spacing}.rds
#      {sample_clean}.bins{bin_size}.rds
# ---------------------------------------------------------------------------
rds_small <- file.path(csawrdsdir, paste0(sample_clean, ".w", width.small, ".d", spacing.small, ".rds"))
rds_large <- file.path(csawrdsdir, paste0(sample_clean, ".w", width.large, ".d", spacing.large, ".rds"))
rds_bins  <- file.path(csawrdsdir, paste0(sample_clean, ".bins", bin.size, ".rds"))

ext.val <- NA  # paired-end: NA

if (all(file.exists(rds_small, rds_large, rds_bins))) {
  cat("Loading pre-computed RDS objects...\n")
  data.small <- readRDS(rds_small)
  data.large <- readRDS(rds_large)
  bins       <- readRDS(rds_bins)
  cat("  Library sizes:", paste(bins$totals, collapse = ", "), "\n")
} else {
  cat("RDS not found — generating window counts from BAMs...\n")

  if (param$pe == "none") {
    max.delay <- 500
    dedup.on  <- readParam(minq = param$minq, pe = "none", dedup = TRUE)
    x        <- correlateReads(c(chip_bams[1], ctrl_bams[1]), max.delay, param = dedup.on)
    ext.val  <- maximizeCcf(x)
    cat("  Estimated fragment length:", ext.val, "bp\n")
  }

  cat("  Counting small windows (", width.small, "bp)...\n")
  data.small <- windowCounts(bam.files, ext = ext.val, width = width.small,
                             spacing = spacing.small, param = param)
  cat("  Counting large windows (", width.large, "bp)...\n")
  data.large <- windowCounts(bam.files, ext = ext.val, width = width.large,
                             spacing = spacing.large, param = param)
  cat("  Counting background bins (", bin.size, "bp)...\n")
  bins       <- windowCounts(bam.files, bin = TRUE, width = bin.size, param = param)
  cat("  Library sizes:", paste(bins$totals, collapse = ", "), "\n")

  saveRDS(data.small, rds_small)
  saveRDS(data.large, rds_large)
  saveRDS(bins,       rds_bins)
  cat("  RDS saved.\n")
}

# ---------------------------------------------------------------------------
# 6. Compute global background enrichment scores (once — independent of cutoff)
# ---------------------------------------------------------------------------
cat("\nComputing global background filter scores...\n")
filter.Global.small <- filterWindowsGlobal(data.small, bins)
filter.Global.large <- filterWindowsGlobal(data.large, bins)

logFC_small <- filter.Global.small$filter   # log2 enrichment over background
logFC_large <- filter.Global.large$filter

cat("  Small window scores: median =", round(median(logFC_small), 3),
    " range [", round(min(logFC_small), 2), ",", round(max(logFC_small), 2), "]\n")
cat("  Large window scores: median =", round(median(logFC_large), 3),
    " range [", round(min(logFC_large), 2), ",", round(max(logFC_large), 2), "]\n\n")

# ---------------------------------------------------------------------------
# 7. Normalisation  (computed once; injected into each iteration)
# ---------------------------------------------------------------------------
cat("Computing normalisation factors (method:", norm_method, ")...\n")
USE_PRECOMPUTED_SF <- toupper(norm_method) %in% c("CHIPSEQSPIKEINFREE", "SPIKEIN")

norm_factors_vec <- NULL   # will be set below

if (USE_PRECOMPUTED_SF) {
  antibody <- sub(".*_", "", sample_id)
  sf_subdir <- switch(toupper(norm_method),
    "CHIPSEQSPIKEINFREE" = "chipseqspikeinfree",
    "SPIKEIN"           = "spikein",
    "chipseqspikeinfree"
  )
  sf_file <- file.path(normdir, sf_subdir, antibody, paste0(antibody, "_SF.txt"))

  if (!file.exists(sf_file)) {
    cat("  WARNING: SF file not found:", sf_file, "— falling back to TMM.\n")
    USE_PRECOMPUTED_SF <- FALSE
  } else {
    sf_tbl <- read.table(sf_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    igG_sf_file <- file.path(normdir, sf_subdir, "IgG", "IgG_SF.txt")
    if (file.exists(igG_sf_file)) {
      igG_sf_tbl <- read.table(igG_sf_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
      sf_tbl <- rbind(sf_tbl, igG_sf_tbl[!igG_sf_tbl$ID %in% sf_tbl$ID, ])
    }
    sf_lookup      <- setNames(sf_tbl$SF, basename(sf_tbl$ID))
    all_basenames  <- basename(bam.files)
    precomputed_sf <- sapply(all_basenames, function(bn) {
      sf <- sf_lookup[bn]
      if (is.na(sf)) { cat("  NOTE: no SF for", bn, "— defaulting to 1\n"); return(1) }
      sf
    })
    inv_sf           <- 1 / precomputed_sf
    norm_factors_vec <- inv_sf / mean(inv_sf)
    cat("  Precomputed norm factors:", paste(round(norm_factors_vec, 4), collapse = ", "), "\n\n")
  }
}

# Compute TMM-based norm factors on the full bins object (used as fallback or default)
bins_norm <- if (!USE_PRECOMPUTED_SF) normFactors(bins) else bins
if (!USE_PRECOMPUTED_SF) {
  norm_factors_vec <- bins_norm$norm.factors
  cat("  TMM norm factors:", paste(round(norm_factors_vec, 4), collapse = ", "), "\n\n")
}

# ---------------------------------------------------------------------------
# 8. Design matrix + contrast  (fixed across all iterations)
# ---------------------------------------------------------------------------
ctrl_label <- if (py$EXPERIMENT %in% c("cutandrun", "cutandtag")) "IgG" else "Input"
chip_label <- "ChIP"
grouping   <- factor(c(rep(chip_label, n_chip), rep(ctrl_label, n_ctrl)),
                     levels = c(chip_label, ctrl_label))
design     <- model.matrix(~0 + grouping)
colnames(design) <- levels(grouping)
contrast   <- makeContrasts(contrasts = paste0(chip_label, " - ", ctrl_label),
                            levels    = design)
has_reps   <- (nrow(design) - ncol(design)) > 0
cat("  Grouping:", paste(as.character(grouping), collapse = ", "), "\n")
cat("  Replicates available:", has_reps, "\n\n")

# ---------------------------------------------------------------------------
# 9. Parameter grid
# ---------------------------------------------------------------------------
grid_config <- list(
  c_small = seq(2.0, 6.0, by = 0.2),
  c_large = seq(1.0, 5.0, by = 0.2)
)
grid_df <- expand.grid(c_small = grid_config$c_small,
                       c_large = grid_config$c_large)
cat("Grid search: ", nrow(grid_df), "combinations (",
    length(grid_config$c_small), "x", length(grid_config$c_large), ")\n\n")

# ---------------------------------------------------------------------------
# 10. Pre-extract plain serializable objects from SE (avoids bplapply crash)
#     BiocParallel's MulticoreParam serializes all arguments sent to workers.
#     RangedSummarizedExperiment objects contain environments that cannot be
#     serialized — extract the plain matrix and GRanges beforehand.
# ---------------------------------------------------------------------------
counts_small  <- assay(data.small)        # plain integer matrix
counts_large  <- assay(data.large)
ranges_small  <- rowRanges(data.small)    # GRanges — serializable
ranges_large  <- rowRanges(data.large)

# ---------------------------------------------------------------------------
# 11. Worker function — one iteration
#     Receives only plain scalars, matrices, and GRanges (all serializable).
# ---------------------------------------------------------------------------
run_one_iter <- function(i, grid_df,
                         counts_small, counts_large,
                         ranges_small, ranges_large,
                         logFC_small, logFC_large,
                         norm_factors_vec, design, contrast,
                         has_reps, merge.bp, qvalue, lfc) {

  suppressPackageStartupMessages({
    library(GenomicRanges)
    library(edgeR)
    library(csaw)
  })

  cs <- grid_df$c_small[i]
  cl <- grid_df$c_large[i]

  # -- Filter --
  keep_s <- logFC_small >= cs
  keep_l <- logFC_large >= cl

  n_small_kept <- sum(keep_s)
  n_large_kept <- sum(keep_l)

  # Bail out if nothing retained
  if (n_small_kept == 0 && n_large_kept == 0) {
    return(data.frame(c_small = cs, c_large = cl,
                      n_small_kept = 0L, n_large_kept = 0L,
                      n_clusters = 0L, n_sig = 0L,
                      width_min = NA, width_q25 = NA, width_median = NA,
                      width_q75 = NA, width_max = NA))
  }

  # -- Fit edgeR separately on small and large windows --
  # Small
  counts_s <- counts_small[keep_s, , drop = FALSE]
  dge_s    <- DGEList(counts = counts_s)
  dge_s$samples$norm.factors <- norm_factors_vec
  if (has_reps) {
    dge_s <- tryCatch(estimateDisp(dge_s, design, robust = FALSE),
                      error = function(e) { dge_s$common.dispersion <- 0.05; dge_s })
    fit_s <- glmQLFit(dge_s, design, robust = FALSE)
    res_s <- glmQLFTest(fit_s, contrast = contrast)
  } else {
    dge_s$common.dispersion <- 0.05
    fit_s <- glmFit(dge_s, design)
    res_s <- glmLRT(fit_s, contrast = contrast)
  }

  # Large
  counts_l <- counts_large[keep_l, , drop = FALSE]
  dge_l    <- DGEList(counts = counts_l)
  dge_l$samples$norm.factors <- norm_factors_vec
  if (has_reps) {
    dge_l <- tryCatch(estimateDisp(dge_l, design, robust = FALSE),
                      error = function(e) { dge_l$common.dispersion <- 0.05; dge_l })
    fit_l <- glmQLFit(dge_l, design, robust = FALSE)
    res_l <- glmQLFTest(fit_l, contrast = contrast)
  } else {
    dge_l$common.dispersion <- 0.05
    fit_l <- glmFit(dge_l, design)
    res_l <- glmLRT(fit_l, contrast = contrast)
  }

  # -- Multi-scale merge via mergeResultsList (canonical csaw approach) --
  # Wrap ranges back into minimal SummarizedExperiment-like list entries
  # that mergeResultsList expects: just needs rowRanges and the result table.
  # We use a named list approach compatible with csaw >= 1.28.
  gr_s <- ranges_small[keep_s]
  gr_l <- ranges_large[keep_l]

  merged <- mergeResultsList(
    list(SummarizedExperiment::SummarizedExperiment(
           assays = list(counts = counts_s), rowRanges = gr_s),
         SummarizedExperiment::SummarizedExperiment(
           assays = list(counts = counts_l), rowRanges = gr_l)),
    tab.list   = list(res_s$table, res_l$table),
    equiweight = TRUE,
    tol        = merge.bp
  )

  tabcom  <- merged$combined
  tabbest <- merged$best

  n_clusters <- nrow(tabcom)
  sig_mask   <- tabcom$FDR < qvalue & tabbest$rep.logFC > lfc

  sig_regions <- merged$regions[which(sig_mask)]
  n_sig       <- sum(sig_mask)

  wq <- if (n_sig > 0) {
    w <- width(sig_regions)
    c(min = min(w), q25 = unname(quantile(w, 0.25)), median = median(w),
      q75 = unname(quantile(w, 0.75)), max = max(w))
  } else {
    c(min = NA, q25 = NA, median = NA, q75 = NA, max = NA)
  }

  data.frame(c_small      = cs,
             c_large      = cl,
             n_small_kept = n_small_kept,
             n_large_kept = n_large_kept,
             n_clusters   = n_clusters,
             n_sig        = n_sig,
             width_min    = wq["min"],
             width_q25    = wq["q25"],
             width_median = wq["median"],
             width_q75    = wq["q75"],
             width_max    = wq["max"])
}

# ---------------------------------------------------------------------------
# 12. Parallel grid search  (SnowParam — socket-based, avoids edgeR env crash)
# ---------------------------------------------------------------------------
n_cores <- min(8L, parallel::detectCores(logical = FALSE))
cat("Running grid search on", n_cores, "cores (SnowParam)...\n")
bp <- SnowParam(workers = n_cores, type = "SOCK", progressbar = TRUE)

results_list <- bplapply(
  seq_len(nrow(grid_df)),
  run_one_iter,
  grid_df          = grid_df,
  counts_small     = counts_small,
  counts_large     = counts_large,
  ranges_small     = ranges_small,
  ranges_large     = ranges_large,
  logFC_small      = logFC_small,
  logFC_large      = logFC_large,
  norm_factors_vec = norm_factors_vec,
  design           = design,
  contrast         = contrast,
  has_reps         = has_reps,
  merge.bp         = merge.bp,
  qvalue           = qvalue,
  lfc              = lfc,
  BPPARAM          = bp
)

results_df <- do.call(rbind, results_list)
results_df <- results_df[order(-results_df$n_sig), ]
rownames(results_df) <- NULL

cat("Grid search complete. Top 5 results:\n")
print(head(results_df, 5))

# ---------------------------------------------------------------------------
# 13. Save results CSV
# ---------------------------------------------------------------------------
csv_path <- file.path(outdir, paste0(sample_clean, "_cutoff_results.csv"))
write.csv(results_df, csv_path, row.names = FALSE)
cat("\nResults saved:", csv_path, "\n")

# ---------------------------------------------------------------------------
# 14. Pareto heatmap
# ---------------------------------------------------------------------------
heatmap_path <- file.path(outdir, paste0(sample_clean, "_pareto_heatmap.pdf"))
pdf(heatmap_path, width = 10, height = 7)

# Identify the optimal cell
opt_idx   <- which.max(results_df$n_sig)
opt_cs    <- results_df$c_small[opt_idx]
opt_cl    <- results_df$c_large[opt_idx]
opt_n_sig <- results_df$n_sig[opt_idx]

p_heat <- ggplot(results_df, aes(x = factor(round(c_small, 2)),
                                 y = factor(round(c_large, 2)),
                                 fill = n_sig)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = n_sig), size = 2.8, color = "white", fontface = "bold") +
  # Highlight optimal cell with a red border
  geom_tile(data = results_df[results_df$c_small == opt_cs & results_df$c_large == opt_cl, ],
            aes(x = factor(round(c_small, 2)), y = factor(round(c_large, 2))),
            fill = NA, color = "firebrick", linewidth = 1.2) +
  scale_fill_gradientn(
    colours  = c("#1a1a2e", "#16213e", "#0f3460", "#533483", "#e94560", "#f5a623"),
    name     = "Significant\nDBRs"
  ) +
  labs(
    title    = paste0(sample_clean, " - csaw filtering cut-off optimisation"),
    subtitle = sprintf("Optimal: c_small = %.1f, c_large = %.1f  |  max DBRs = %d  (FDR < %.2f, logFC > %.1f)",
                       opt_cs, opt_cl, opt_n_sig, qvalue, lfc),
    x        = "Small-window cut-off (c_small, log2 FC over background)",
    y        = "Large-window cut-off (c_large, log2 FC over background)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold"),
    axis.text.x   = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

print(p_heat)
dev.off()
cat("Pareto heatmap saved:", heatmap_path, "\n")

# ---------------------------------------------------------------------------
# 15. Re-run at optimal cutoff — produce BCV plot + filter histograms + BED
# ---------------------------------------------------------------------------
cat("\n=== Generating outputs for optimal cut-off (c_small =", opt_cs, ", c_large =", opt_cl, ") ===\n")

keep_s_opt <- logFC_small >= opt_cs
keep_l_opt <- logFC_large >= opt_cl

cat("  Small windows retained:", sum(keep_s_opt), "/", length(keep_s_opt), "\n")
cat("  Large windows retained:", sum(keep_l_opt), "/", length(keep_l_opt), "\n")

# --- Filter histograms ---
hist_path <- file.path(outdir, paste0(sample_clean, "_opt_filter_histograms.pdf"))
pdf(hist_path, width = 10, height = 5)
df_hist <- rbind(
  data.frame(score  = logFC_small,
             window = paste0("Small (", width.small, " bp)"),
             cutoff = opt_cs),
  data.frame(score  = logFC_large,
             window = paste0("Large (", width.large, " bp)"),
             cutoff = opt_cl)
)
cutoffs_df <- data.frame(
  window    = c(paste0("Small (", width.small, " bp)"), paste0("Large (", width.large, " bp)")),
  cutoff    = c(opt_cs, opt_cl),
  kept_pct  = c(round(100 * mean(keep_s_opt), 1), round(100 * mean(keep_l_opt), 1))
)
p_hist <- ggplot(df_hist, aes(x = score)) +
  geom_histogram(bins = 100, fill = "#4e79a7", color = "white", linewidth = 0.1) +
  geom_vline(data = cutoffs_df, aes(xintercept = cutoff),
             color = "firebrick", linetype = "dashed", linewidth = 0.8) +
  geom_text(data = cutoffs_df,
            aes(x = cutoff, y = Inf,
                label = sprintf("cutoff = %.2f\n(%s%% kept)", cutoff, kept_pct)),
            hjust = -0.05, vjust = 1.5, color = "firebrick", size = 3) +
  facet_wrap(~window, scales = "free") +
  labs(title    = paste0(sample_clean, " - filter histograms at optimal cut-off"),
       subtitle = sprintf("c_small = %.1f, c_large = %.1f", opt_cs, opt_cl),
       x = "log2 fold-enrichment over background", y = "Number of windows") +
    theme_bw(base_size = 11) +
  theme(strip.background = element_rect(fill = "grey92"),
        plot.title       = element_text(face = "bold"))
print(p_hist)
dev.off()
cat("  Filter histograms saved:", hist_path, "\n")

# --- Full edgeR fit at optimal cut-off for BCV plot + BED ---
# Fit separately on small and large (mirrors grid search worker logic)
counts_s_opt <- counts_small[keep_s_opt, , drop = FALSE]
dge_s_opt    <- DGEList(counts = counts_s_opt)
dge_s_opt$samples$norm.factors <- norm_factors_vec

counts_l_opt <- counts_large[keep_l_opt, , drop = FALSE]
dge_l_opt    <- DGEList(counts = counts_l_opt)
dge_l_opt$samples$norm.factors <- norm_factors_vec

if (has_reps) {
  dge_s_opt <- estimateDisp(dge_s_opt, design)
  dge_l_opt <- estimateDisp(dge_l_opt, design)
  fit_s_opt <- glmQLFit(dge_s_opt, design)
  fit_l_opt <- glmQLFit(dge_l_opt, design)
  res_s_opt <- glmQLFTest(fit_s_opt, contrast = contrast)
  res_l_opt <- glmQLFTest(fit_l_opt, contrast = contrast)
} else {
  dge_s_opt$common.dispersion <- 0.05
  dge_l_opt$common.dispersion <- 0.05
  fit_s_opt <- glmFit(dge_s_opt, design)
  fit_l_opt <- glmFit(dge_l_opt, design)
  res_s_opt <- glmLRT(fit_s_opt, contrast = contrast)
  res_l_opt <- glmLRT(fit_l_opt, contrast = contrast)
}

# BCV plot (use the small-window model which has more windows / better BCV estimate)
bcv_path <- file.path(outdir, paste0(sample_clean, "_opt_bcv.pdf"))
pdf(bcv_path, width = 7, height = 5)
plotBCV(dge_s_opt,
        main = sprintf("%s - BCV at optimal cut-off (c_s=%.1f, c_l=%.1f)",
                       sample_clean, opt_cs, opt_cl))
dev.off()
cat("  BCV plot saved:", bcv_path, "\n")

# --- Cluster-level results via mergeResultsList ---
gr_s_opt <- ranges_small[keep_s_opt]
gr_l_opt <- ranges_large[keep_l_opt]

merged_opt <- mergeResultsList(
  list(SummarizedExperiment::SummarizedExperiment(
         assays = list(counts = counts_s_opt), rowRanges = gr_s_opt),
       SummarizedExperiment::SummarizedExperiment(
         assays = list(counts = counts_l_opt), rowRanges = gr_l_opt)),
  tab.list   = list(res_s_opt$table, res_l_opt$table),
  equiweight = TRUE,
  tol        = merge.bp
)

tabcom_opt  <- merged_opt$combined
tabbest_opt <- merged_opt$best

is_sig_opt <- tabcom_opt$FDR <= qvalue
is_enr_opt <- tabbest_opt$rep.logFC > lfc
keep_peaks <- is_sig_opt & is_enr_opt

ranges_out <- merged_opt$regions[keep_peaks]
fdr_out    <- tabcom_opt$FDR[keep_peaks]
lfc_out    <- tabbest_opt$rep.logFC[keep_peaks]

cat("\n  Significant DBRs (FDR <=", qvalue, ", logFC >", lfc, "):", sum(keep_peaks), "\n")

if (sum(keep_peaks) > 0) {
  # Build BED-compatible GRanges
  mcols(ranges_out)$name  <- paste0("peak_", seq_len(sum(keep_peaks)))
  mcols(ranges_out)$score <- as.integer(round(-10 * log10(pmax(fdr_out, 1e-300))))

  # Harmonise chromosome naming
  new_lvls <- mapSeqlevels(seqlevels(ranges_out), style = "UCSC")
  keep_chr <- !is.na(new_lvls)
  ranges_out <- keepSeqlevels(ranges_out, seqlevels(ranges_out)[keep_chr], pruning.mode = "coarse")
  seqlevels(ranges_out) <- new_lvls[keep_chr]

  bed_name <- paste0(
    sample_clean,
    ".w",   width.small, ".d", spacing.small,
    ".cs",  opt_cs,
    ".w",   width.large, ".d", spacing.large,
    ".cl",  opt_cl,
    ".lfc", lfc,
    ".merge", merge.bp,
    ".", tolower(norm_method),
    ".opt.bed"
  )
  bed_path <- file.path(beddir, bed_name)
  rtracklayer::export.bed(ranges_out, bed_path)
  cat("  Optimal peak BED saved:", bed_path, "\n")
} else {
  cat("  No significant DBRs found at optimal cutoff. Skipping BED export.\n")
}

cat("\n=== optimize_csaw_cutoffs.R complete ===\n")
cat("  Optimal c_small:", opt_cs, "| c_large:", opt_cl,
    "| Significant DBRs:", opt_n_sig, "\n")
