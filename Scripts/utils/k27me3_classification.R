#!/usr/bin/env Rscript

# =============================================================================
# step_7a_k27me3_classification.R
#
# Defines two mutually exclusive H3K27me3 peak sets from inhibitor treatment
# data (analogous to step_7a_define_nucleation_spreading.R):
#
#   Retained H3K27me3 — peaks from H3K27me3_1uMcpd cells (H3K27me3 persists
#                        despite EZH2 inhibition), after reduce().
#
#   Lost H3K27me3     — peaks from H3K27me3_DMSO cells with any overlap to
#                        retained (1uM) sites subtracted via GRanges::setdiff().
#                        Regions narrower than 200 bp after subtraction are
#                        discarded.
#
# Supports the same peak callers as step_7a (--peak-caller flag):
#   csaw   -> config.PEAKDIR/csaw/   (*.tmm.bed)
#   macs   -> config.PEAKDIR/macs/   (*_peaks.broadPeak or narrowPeak)
#   seacr  -> config.PEAKDIR/seacr/  (*.relaxed.bed)
#
# Usage (run from project root, NOT via bsub):
#   Rscript Scripts/step_7a_k27me3_classification.R --peak-caller csaw
#
# Alternatively, supply explicit BED paths to bypass the peak-dir search:
#   Rscript Scripts/step_7a_k27me3_classification.R \
#     --bed-dmso  Analysis_Data/peaks/csaw/H3K27me3_DMSO.w150.d50.filt3.w2000.d500.filt1.5.lfc1.merge100.tmm.bed \
#     --bed-1um   Analysis_Data/peaks/csaw/H3K27me3_1uMcpd.w150.d50.filt3.w2000.d500.filt1.5.lfc1.merge100.tmm.bed
#
# Outputs (config.ANALYSIS_DATA/k27me3_classification/):
#   retained_H3K27me3.bed  — 0-based BED, 1uMcpd peaks after reduce()
#   lost_H3K27me3.bed      — 0-based BED, DMSO setdiff(retained) >= 200 bp
#   k27me2_in_lost_k27me3.bed — 0-based BED, subset of K27me2 peaks inside lost H3K27me3
# =============================================================================

suppressPackageStartupMessages({
  library(reticulate)
  library(GenomicRanges)
  library(rtracklayer)
})

# =============================================================================
# 1. CLI argument parsing
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)

get_flag <- function(flag, required = FALSE, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0L) {
    if (required) stop(sprintf("Missing required argument: %s", flag))
    return(default)
  }
  val_idx <- idx[1L] + 1L
  if (val_idx > length(args) || startsWith(args[val_idx], "--")) {
    stop(sprintf("No value supplied for %s", flag))
  }
  args[val_idx]
}

CALLER   <- get_flag("--peak-caller", required = FALSE, default = "csaw")
CALLER   <- tolower(trimws(CALLER))
BED_DMSO <- get_flag("--bed-dmso",   required = FALSE, default = NULL)
BED_1UM  <- get_flag("--bed-1um",    required = FALSE, default = NULL)
BED_K27ME2 <- get_flag("--bed-k27me2", required = FALSE, default = NULL)

valid_callers <- c("macs", "csaw", "seacr")
if (!CALLER %in% valid_callers) {
  stop(sprintf(
    "Invalid --peak-caller '%s'. Must be one of: %s",
    CALLER, paste(valid_callers, collapse = ", ")
  ))
}

cat(sprintf("\n========================================\n"))
cat(sprintf("  step_7a_k27me3  |  peak caller: %s\n", toupper(CALLER)))
cat(sprintf("========================================\n\n"))

# =============================================================================
# 2. Load config.py via reticulate
# =============================================================================
use_python(Sys.which("python"), required = TRUE)
py_run_file("Scripts/config.py")

PEAKDIR     <- py$PEAKDIR
ANALYSIS_DATA <- py$ANALYSIS_DATA
MACS_BROAD  <- isTRUE(py$MACS_CALL_BROAD_PEAKS)

# Caller-specific peak directory
caller_peak_dir <- file.path(PEAKDIR, CALLER)

# Output directory
out_dir <- file.path(ANALYSIS_DATA, "k27me3_classification")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("Peak directory  : %s\n", caller_peak_dir))
cat(sprintf("Output directory: %s\n\n", out_dir))

# =============================================================================
# 3. Load peaks
# =============================================================================
cat("=== Loading peaks ===\n")

if (!is.null(BED_DMSO) && !is.null(BED_1UM)) {
  # --- Explicit BED paths supplied -------------------------------------------
  cat(sprintf("Using explicit BED files:\n  DMSO : %s\n  1uM  : %s\n\n",
              BED_DMSO, BED_1UM))

  read_bed6 <- function(path) {
    if (!file.exists(path)) stop("BED file not found: ", path)
    df <- read.table(path, header = FALSE, sep = "\t",
                     col.names = c("chr", "start", "end", "name", "score", "strand"),
                     stringsAsFactors = FALSE)
    makeGRangesFromDataFrame(df,
      seqnames.field          = "chr",
      start.field             = "start",
      end.field               = "end",
      strand.field            = "strand",
      starts.in.df.are.0based = TRUE,
      keep.extra.columns      = FALSE)
  }

  dmso_raw <- read_bed6(BED_DMSO)
  um1_raw  <- read_bed6(BED_1UM)

} else {
  # --- Auto-discover via read_peaks.R ----------------------------------------
  source("Scripts/utils/read_peaks.R")

  dmso_raw <- load_caller_peaks(caller_peak_dir, CALLER, "H3K27me3_DMSO",   macs_broad = MACS_BROAD)
  um1_raw  <- load_caller_peaks(caller_peak_dir, CALLER, "H3K27me3_1uMcpd", macs_broad = MACS_BROAD)
}

cat(sprintf("H3K27me3_DMSO   raw peaks : %d\n", length(dmso_raw)))
cat(sprintf("H3K27me3_1uMcpd raw peaks : %d\n", length(um1_raw)))

# Read K27me2 peaks if supplied or available
if (is.null(BED_K27ME2)) {
  default_k27me2_path <- file.path(PEAKDIR, "csaw", "H3K27me2_1uMcpd.w150.d50.filt2.w2000.d500.filt1.lfc1.merge100.tmm.bed")
  if (file.exists(default_k27me2_path)) {
    BED_K27ME2 <- default_k27me2_path
  }
}

if (!is.null(BED_K27ME2) && file.exists(BED_K27ME2)) {
  cat(sprintf("\nLoading K27me2 peaks from: %s\n", BED_K27ME2))
  read_bed6 <- function(path) {
    df <- read.table(path, header = FALSE, sep = "\t",
                     col.names = c("chr", "start", "end", "name", "score", "strand"),
                     stringsAsFactors = FALSE)
    makeGRangesFromDataFrame(df,
      seqnames.field          = "chr",
      start.field             = "start",
      end.field               = "end",
      strand.field            = "strand",
      starts.in.df.are.0based = TRUE,
      keep.extra.columns      = FALSE)
  }
  k27me2_raw <- read_bed6(BED_K27ME2)
  cat(sprintf("H3K27me2 raw peaks : %d\n", length(k27me2_raw)))
} else {
  k27me2_raw <- NULL
}

# =============================================================================
# 4. Define Retained and Lost sites
# =============================================================================
cat("\n=== Defining retained and lost H3K27me3 sites ===\n")

# reduce() merges overlapping/adjacent peaks (equivalent to bedtools merge)
retained_gr <- reduce(um1_raw,  ignore.strand = TRUE)
dmso_gr     <- reduce(dmso_raw, ignore.strand = TRUE)

cat(sprintf("H3K27me3_1uMcpd after reduce : %d regions\n", length(retained_gr)))
cat(sprintf("H3K27me3_DMSO   after reduce : %d regions\n", length(dmso_gr)))

# setdiff() is the GRanges equivalent of bedtools subtract:
# returns portions of DMSO that do NOT overlap any 1uM peak.
lost_raw <- setdiff(dmso_gr, retained_gr, ignore.strand = TRUE)
cat(sprintf("Lost (DMSO setdiff 1uM, pre-filter): %d regions\n", length(lost_raw)))

# Discard sub-200 bp fragments produced by edge-trimming
MIN_WIDTH <- 200L
lost_gr   <- lost_raw[width(lost_raw) >= MIN_WIDTH]
cat(sprintf("Lost (>= %d bp) : %d regions\n", MIN_WIDTH, length(lost_gr)))

# =============================================================================
# 5. Sanity check: zero overlap expected between the two final sets
# =============================================================================
cat("\n=== Sanity check ===\n")
n_overlap <- length(findOverlaps(retained_gr, lost_gr))
if (n_overlap > 0L) {
  warning(sprintf(
    "UNEXPECTED: %d overlaps found between retained and lost H3K27me3 sites!",
    n_overlap
  ))
} else {
  cat("OK — 0 overlaps between retained and lost H3K27me3 sites.\n")
}

# =============================================================================
# 6. Intersect K27me2 peaks with Lost H3K27me3 regions
# =============================================================================
if (!is.null(k27me2_raw)) {
  cat("\n=== Intersecting K27me2 peaks with Lost H3K27me3 regions ===\n")
  k27me2_in_lost <- subsetByOverlaps(k27me2_raw, lost_gr, ignore.strand = TRUE)
  cat(sprintf("K27me2 peaks occurring in Lost H3K27me3 regions: %d (of %d total K27me2 peaks)\n",
              length(k27me2_in_lost), length(k27me2_raw)))
} else {
  k27me2_in_lost <- NULL
}

# =============================================================================
# 7. Summary statistics
# =============================================================================
summarise_gr <- function(gr, label) {
  w <- width(gr)
  cat(sprintf(
    "\n%s\n  Regions      : %d\n  Total (Mb)   : %.2f\n  Median width : %d bp\n  Width range  : %d - %d bp\n",
    label, length(gr), sum(w) / 1e6, as.integer(median(w)), min(w), max(w)
  ))
}

cat("\n=== Summary ===")
summarise_gr(retained_gr, "Retained H3K27me3 (1uMcpd peaks, reduced)")
summarise_gr(lost_gr,     "Lost H3K27me3 (DMSO setdiff 1uM, >= 200 bp)")
if (!is.null(k27me2_in_lost)) {
  summarise_gr(k27me2_in_lost, "K27me2 peaks in Lost H3K27me3 regions")
}

# =============================================================================
# 8. Export BED files (0-based)
# =============================================================================
cat("\n=== Exporting BED files ===\n")

export_bed_auto <- function(gr, path) {
  df <- as.data.frame(gr)
  if (all(c("name", "score", "strand") %in% colnames(df))) {
    out <- df[, c("seqnames", "start", "end", "name", "score", "strand")]
    colnames(out)[1:3] <- c("chr", "start", "end")
    out$start <- out$start - 1L
  } else {
    out <- df[, c("seqnames", "start", "end")]
    colnames(out) <- c("chr", "start", "end")
    out$start <- out$start - 1L
  }
  write.table(out, file = path, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
  cat(sprintf("  Written: %s  (%d regions, %.2f Mb)\n",
              basename(path), nrow(out), sum(out$end - out$start) / 1e6))
}

retained_bed <- file.path(out_dir, "retained_H3K27me3.bed")
lost_bed     <- file.path(out_dir, "lost_H3K27me3.bed")

export_bed_auto(retained_gr, retained_bed)
export_bed_auto(lost_gr,     lost_bed)

if (!is.null(k27me2_in_lost)) {
  k27me2_intersect_bed <- file.path(out_dir, "k27me2_in_lost_k27me3.bed")
  export_bed_auto(k27me2_in_lost, k27me2_intersect_bed)
}

# =============================================================================
# 9. Done
# =============================================================================
cat(sprintf("\n========================================\n"))
cat(sprintf("  step_7a_k27me3 complete  |  caller: %s\n", toupper(CALLER)))
cat(sprintf("========================================\n"))
cat("\nOutputs:\n")
cat(sprintf("  %s\n", retained_bed))
cat(sprintf("  %s\n", lost_bed))
if (!is.null(k27me2_in_lost)) {
  cat(sprintf("  %s\n", k27me2_intersect_bed))
}
cat(sprintf("\nDone.\n"))
