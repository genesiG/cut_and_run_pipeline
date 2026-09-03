#!/usr/bin/env Rscript

# =============================================================================
# callpeaks_csaw.R - Absolute peak calling via csaw differential binding
#
# Strategy: compares one or more ChIP/CUT&RUN replicates against their matched
#           controls (Input for ChIP-seq, IgG for CUT&RUN) to identify absolute
#           enrichment regions. All replicates for both groups are passed in and
#           modelled jointly in the edgeR GLM.
#
# This script is NOT meant to be run interactively — it is driven by the
# Python wrapper step_5c_callpeaks_csaw.py, which builds the bsub batch script.
#
# CLI arguments (all required):
#   --chip_bams   One or more paths to ChIP/CUT&TAG BAM files (space-separated)
#   --ctrl_bams   One or more paths to control (IgG) BAM files
#   --sample_id   Label for output files (e.g. "K27M_K27me3")
#   --workdir     Project root directory (WORKDIR from config)
#   --codedir     Scripts directory (CODEDIR from config)
#   Uses TMM normalization on large background bins.
# =============================================================================

suppressPackageStartupMessages({
  options(repos = c(CRAN = "https://cran.rstudio.com"))
  if (!require("reticulate", quietly = TRUE)) install.packages("reticulate")
  library(reticulate)
})

# ---------------------------------------------------------------------------
# 1. Parse command-line arguments
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
cat("=== callpeaks_csaw.R ===\n")
cat("  sample_id        :", sample_id, "\n")
cat("  ChIP BAMs        :", length(chip_bams), "\n")
for (b in chip_bams) cat("    ", b, "\n")
cat("  Control BAMs     :", length(ctrl_bams), "\n")
for (b in ctrl_bams) cat("    ", b, "\n")
cat("  workdir          :", workdir, "\n")
cat("  codedir          :", codedir, "\n\n")

# ---------------------------------------------------------------------------
# 2. Load config.py to derive workspace paths
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
beddir        <- csawdir
qcdir         <- file.path(csawdir, "qc")
normdir       <- file.path(analysisdir, "normalization")

for (d in c(csawdir, rdsdir, beddir, qcdir, csawrdsdir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# 3. Source csaw_parameters.R
# ---------------------------------------------------------------------------
source(file.path(codedir, "utils", "csaw_parameters.R"))

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE, quietly = TRUE)
}
for (pkg in biopackages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install(pkg)
  }
  library(pkg, character.only = TRUE, quietly = TRUE)
}
suppressPackageStartupMessages(library(ggplot2))

# ---------------------------------------------------------------------------
# 4. Validate BAM files
# ---------------------------------------------------------------------------
bam.files <- c(chip_bams, ctrl_bams)
for (f in bam.files) if (!file.exists(f)) stop(paste("BAM file not found:", f))

n_chip <- length(chip_bams)
n_ctrl <- length(ctrl_bams)
cat("Loaded", n_chip, "ChIP BAM(s) and", n_ctrl, "control BAM(s)\n\n")

# ---------------------------------------------------------------------------
# 5. Estimate fragment length via cross-correlation (PE: NA, SE: auto)
# ---------------------------------------------------------------------------
ext.val <- NA  # paired-end: csaw handles fragment lengths internally

if (param$pe == "none") {
  max.delay <- 500
  dedup.on  <- readParam(minq = param$minq, pe = "none", dedup = TRUE)
  cat("Computing cross-correlation (first ChIP + first control)...\n")
  x        <- correlateReads(c(chip_bams[1], ctrl_bams[1]), max.delay, param = dedup.on)
  frag.len <- maximizeCcf(x)
  cat("  Estimated fragment length:", frag.len, "bp\n\n")
  ext.val <- frag.len
}

# ---------------------------------------------------------------------------
# 6. Count reads into windows across ALL BAMs
# ---------------------------------------------------------------------------
cat("Counting reads into small windows (", width.small, "bp)...\n")
data.small <- windowCounts(bam.files,
                           ext     = ext.val,
                           width   = width.small,
                           spacing = spacing.small,
                           param   = param)

cat("Counting reads into large windows (", width.large, "bp)...\n")
data.large <- windowCounts(bam.files,
                           ext     = ext.val,
                           width   = width.large,
                           spacing = spacing.large,
                           param   = param)

# ---------------------------------------------------------------------------
# 7. Count background bins (all BAMs)
# ---------------------------------------------------------------------------
cat("Counting background bins (", bin.size, "bp)...\n")
bins <- windowCounts(bam.files, bin = TRUE, width = bin.size, param = param)
cat("  Library sizes:", paste(bins$totals, collapse = ", "), "\n\n")

# ---------------------------------------------------------------------------
# 8. Save intermediate RDS objects
# ---------------------------------------------------------------------------
saveRDS(data.small, file.path(csawrdsdir, paste0(sample_id, ".w", width.small,  ".d", spacing.small, ".rds")))
saveRDS(data.large, file.path(csawrdsdir, paste0(sample_id, ".w", width.large,  ".d", spacing.large, ".rds")))
saveRDS(bins,       file.path(csawrdsdir, paste0(sample_id, ".bins", bin.size,  ".rds")))

# ---------------------------------------------------------------------------
# 9. Filter low-abundance windows  +  QC histograms
# ---------------------------------------------------------------------------
cat("Filtering low-abundance windows...\n")

filter.Global.small <- filterWindowsGlobal(data.small, bins)
filter.Global.large <- filterWindowsGlobal(data.large, bins)

keep.small <- filter.Global.small$filter > (small.filt)
keep.large <- filter.Global.large$filter > (large.filt)

cat("  Small windows retained:", sum(keep.small), "/", length(keep.small), "\n")
cat("  Large windows retained:", sum(keep.large), "/", length(keep.large), "\n\n")

# --- QC histograms ---
hist_path <- file.path(qcdir, paste0(sample_id,
                                     ".small.filt", small.filt,
                                     ".large.filt", large.filt,
                                     "_filter_histograms.svg"))
svg(hist_path, width = 10, height = 5)

score_small  <- filter.Global.small$filter
score_large  <- filter.Global.large$filter
cutoff_small <- (small.filt)
cutoff_large <- (large.filt)

df_hist <- rbind(
  data.frame(score = score_small, window = paste0("Small (", width.small, " bp)")),
  data.frame(score = score_large, window = paste0("Large (", width.large, " bp)"))
)
cutoffs <- data.frame(
  window   = c(paste0("Small (", width.small, " bp)"), paste0("Large (", width.large, " bp)")),
  cutoff   = c(cutoff_small, cutoff_large),
  kept_pct = c(round(100 * mean(keep.small), 1), round(100 * mean(keep.large), 1))
)

p <- ggplot(df_hist, aes(x = score)) +
  geom_histogram(bins = 100, fill = "#4e79a7", color = "white", linewidth = 0.1) +
  geom_vline(data = cutoffs,
             aes(xintercept = cutoff),
             color = "firebrick", linetype = "dashed", linewidth = 0.8) +
  geom_text(data = cutoffs,
            aes(x = cutoff, y = Inf,
                label = sprintf("cutoff = %.2f\n(%s%% retained)", cutoff, kept_pct)),
            hjust = -0.05, vjust = 1.5, color = "firebrick", size = 3) +
  facet_wrap(~window, scales = "free") +
  labs(title   = paste0(sample_id, " - global background filter scores"),
       x = "log2 fold-enrichment over background",
       y = "Number of windows") +
  theme_classic(base_size = 14) +
  theme(strip.background = element_rect(fill = "grey92"),
        plot.title = element_text(face = "bold", hjust = 0.5))

print(p)
dev.off()
cat("  Filter histogram saved:", hist_path, "\n\n")

data.small.filt <- data.small[keep.small, ]
data.large.filt <- data.large[keep.large, ]

# ---------------------------------------------------------------------------
# 10. Normalisation
#     use csaw's TMM-on-bins approach for absolute peak calling.
# ---------------------------------------------------------------------------
cat("Normalising (TMM on background bins)...\n")

data.small.filt <- normFactors(bins, se.out = data.small.filt)
data.large.filt <- normFactors(bins, se.out = data.large.filt)
cat("  Small-window TMM factors:", paste(round(data.small.filt$norm.factors, 4), collapse = ", "), "\n")
cat("  Large-window TMM factors:", paste(round(data.large.filt$norm.factors, 4), collapse = ", "), "\n")
cat("\n")

# ---------------------------------------------------------------------------
# 11. Build design matrix
# ---------------------------------------------------------------------------
cat("Fitting edgeR model...\n")

ctrl_label <- if (py$EXPERIMENT %in% c("cutandrun", "cutandtag")) "IgG" else "Input"
chip_label <- "Target"

grouping <- factor(
  c(rep(chip_label, n_chip), rep(ctrl_label, n_ctrl)),
  levels = c(chip_label, ctrl_label)
)

cat("  Design grouping:", paste(as.character(grouping), collapse = ", "), "\n")
design <- model.matrix(~0 + grouping)
colnames(design) <- levels(grouping)

# ---------------------------------------------------------------------------
# 12. Estimate dispersions and fit GLM
# ---------------------------------------------------------------------------
y.small <- asDGEList(data.small.filt)
y.large <- asDGEList(data.large.filt)

has_reps <- (nrow(design) - ncol(design)) > 0

if (has_reps) {
  y.small <- estimateDisp(y.small, design)
  y.large <- estimateDisp(y.large, design)
} else {
  cat("  No replicates detected — using default dispersion 0.05\n")
  y.small$common.dispersion <- 0.05
  y.large$common.dispersion <- 0.05
}

if (has_reps) {
  fit.small <- glmQLFit(y.small, design)
  fit.large <- glmQLFit(y.large, design)
} else {
  fit.small <- glmFit(y.small, design)
  fit.large <- glmFit(y.large, design)
}

# ---------------------------------------------------------------------------
# 13. Examine replicate similarity
# ---------------------------------------------------------------------------
cat("Examining replicate similarity...\n")
mds_path <- file.path(qcdir, paste0(sample_id,
                                    ".small.filt", small.filt,
                                    ".large.filt", large.filt,
                                    "_mds.pdf"))
pdf(mds_path, width = 5, height = 5)
par(mfrow=c(2,2), mar=c(5,4,2,2))
adj.counts <- cpm(y.small, log=TRUE)
for (top in c(100, 500, 1000, 5000)) {
    plotMDS(adj.counts, 
            main=paste("Top ", top, " small windows"), 
            col=c("blue", "blue", "red", "red"),
            labels=c(rep(chip_label, n_chip), rep(ctrl_label, n_ctrl)),
            top=top)
}
adj.counts <- cpm(y.large, log=TRUE)
for (top in c(100, 500, 1000, 5000)) {
    plotMDS(adj.counts, 
            main=paste("Top ", top, " large windows"), 
            col=c("blue", "blue", "red", "red"),
            labels=c(rep(chip_label, n_chip), rep(ctrl_label, n_ctrl)),
            top=top)
}
dev.off()
cat("  MDS plots saved:", mds_path, "\n\n")

# ---------------------------------------------------------------------------
# 14. Test for differential binding
# ---------------------------------------------------------------------------
cat("Testing for differential binding...\n")

contrast <- makeContrasts(
  contrasts = paste0(chip_label, " - ", ctrl_label),
  levels    = design
)

if (has_reps) {
  res.small <- glmQLFTest(fit.small, contrast = contrast)
  res.large <- glmQLFTest(fit.large, contrast = contrast)
} else {
  res.small <- glmLRT(fit.small, contrast = contrast)
  res.large <- glmLRT(fit.large, contrast = contrast)
}


# ---------------------------------------------------------------------------
# 15. Merge windows across scales
# ---------------------------------------------------------------------------
cat("Merging results across window sizes...\n")

merged <- mergeResultsList(
  list(data.small.filt, data.large.filt),
  tab.list   = list(res.small$table, res.large$table),
  equiweight = TRUE,
  tol        = merge.bp
)

tabcom  <- merged$combined
tabbest <- merged$best

is.sig     <- tabcom$FDR <= qvalue
is.sig.pos <- (tabbest$rep.logFC > lfc)[is.sig]

cat("  Total merged regions :", nrow(tabcom), "\n")
cat("  Significant (FDR <=", qvalue, "):", sum(is.sig), "\n")
cat("  Enriched (FDR <=", qvalue, " and logFC >", lfc, "):", sum(is.sig.pos), "\n")
cat("  Direction table:\n")
print(table(tabcom$direction[is.sig]))

# ---------------------------------------------------------------------------
# 16. Build and export BED file
# ---------------------------------------------------------------------------
ranges        <- merged$regions
mcols(ranges) <- DataFrame(tabcom, best.logFC = tabbest$rep.logFC)
genome(ranges) <- build

bed <- ranges %>%
  as.data.frame() %>%
  dplyr::filter(rep.logFC > lfc & FDR <= qvalue)

bed$seqnames <- mapSeqlevels(as.character(bed$seqnames), style = style) %>% as.factor()
bed <- dplyr::filter(bed, !is.na(bed$seqnames))

bed$name  <- paste0("peak_", seq_len(nrow(bed)))
bed$score <- as.integer(round(-10 * log10(pmax(bed$FDR, 1e-300))))
bed$signal <- as.integer(bed$rep.logFC)

bed <- makeGRangesFromDataFrame(bed, keep.extra.columns = TRUE)

bed.name <- paste0(
  sample_id,
  ".w",   width.small, ".d", spacing.small, ".filt", small.filt,
  ".w",   width.large, ".d", spacing.large, ".filt", large.filt,
  ".lfc", lfc, ".merge", merge.bp,
  ".tmm.bed"
)

bed.path <- file.path(beddir, bed.name)
rtracklayer::export.bed(bed, bed.path)

cat("\nPeak BED written to:", bed.path, "\n")
cat("  Peaks exported:   ", length(bed), "\n")
cat("\n=== Done ===\n")
