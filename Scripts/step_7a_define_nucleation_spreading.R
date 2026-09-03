#!/usr/bin/env Rscript

# =============================================================================
# step_7a_define_nucleation_spreading.R
#
# Defines two mutually exclusive H3K27me3 peak sets from CUT&TAG data:
#
#   Nucleation sites — reduced peaks from K27M_DMSO cells (H3K27me3 is retained
#                      despite the dominant-negative H3K27M histone mutation).
#
#   Spreading sites  — reduced peaks from K27MKO_DMSO cells with any overlap to
#                      nucleation sites subtracted out via GRanges::setdiff().
#                      Regions narrower than 200 bp after subtraction are
#                      discarded.
#
# Supports three peak callers via --peak-caller CLI flag:
#   macs   -> config.PEAKDIR/macs/   (*_peaks.broadPeak or narrowPeak)
#   csaw   -> config.PEAKDIR/csaw/   (*.tmm.bed)
#   seacr  -> config.PEAKDIR/seacr/  (*.relaxed.bed)
#
# Full annotation, metrics, and visualizations are handled by the companion
# batch job (step_7b_peak_annotation.py -> utils/peak_annotation.R) which
# runs inside bsub with sufficient memory for Bioconductor annotation packages.
#
# Usage (run from project root, NOT via bsub):
#   Rscript Scripts/step_7a_define_nucleation_spreading.R --peak-caller csaw
#   Rscript Scripts/step_7a_define_nucleation_spreading.R --peak-caller seacr
#   Rscript Scripts/step_7a_define_nucleation_spreading.R --peak-caller macs
#
# Outputs (config.PEAKANNODIR/<caller>/):
#   nucleation_sites.bed   — 0-based BED, K27M_DMSO peaks after reduce()
#   spreading_sites.bed    — 0-based BED, K27MKO_DMSO setdiff nucleation >= 200 bp
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

CALLER <- get_flag("--peak-caller", required = TRUE)
CALLER <- tolower(trimws(CALLER))

valid_callers <- c("macs", "csaw", "seacr")
if (!CALLER %in% valid_callers) {
  stop(sprintf(
    "Invalid --peak-caller '%s'. Must be one of: %s",
    CALLER, paste(valid_callers, collapse = ", ")
  ))
}

cat(sprintf("\n========================================\n"))
cat(sprintf("  step_7a  |  peak caller: %s\n", toupper(CALLER)))
cat(sprintf("========================================\n\n"))

# =============================================================================
# 2. Load config.py via reticulate
# =============================================================================
use_python(Sys.which("python"), required = TRUE)
py_run_file("Scripts/config.py")

PEAKDIR     <- py$PEAKDIR
PEAKANNODIR <- py$PEAKANNODIR
MACS_BROAD  <- isTRUE(py$MACS_CALL_BROAD_PEAKS)

# Caller-specific peak directory
caller_peak_dir <- file.path(PEAKDIR, CALLER)

# Output directory: PEAKANNODIR/<caller>/
out_dir <- file.path(PEAKANNODIR, CALLER)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("Peak directory  : %s\n", caller_peak_dir))
cat(sprintf("Output directory: %s\n\n", out_dir))

# =============================================================================
# 3. Source per-caller peak reader helper
# =============================================================================
source("Scripts/utils/read_peaks.R")

# =============================================================================
# 4. Load peaks — K27M_DMSO (nucleation source) and K27MKO_DMSO (spreading source)
# =============================================================================
cat("=== Loading peaks ===\n")

k27m_raw   <- load_caller_peaks(caller_peak_dir, CALLER, "K27M_DMSO",   macs_broad = MACS_BROAD)
k27mko_raw <- load_caller_peaks(caller_peak_dir, CALLER, "K27MKO_DMSO", macs_broad = MACS_BROAD)

cat(sprintf("K27M_DMSO   raw peaks : %d\n", length(k27m_raw)))
cat(sprintf("K27MKO_DMSO raw peaks : %d\n", length(k27mko_raw)))

# =============================================================================
# 5. Define nucleation and spreading sites
# =============================================================================
cat("\n=== Defining nucleation and spreading sites ===\n")

# reduce() merges overlapping/adjacent peaks (equivalent to bedtools merge)
nucleation_gr <- reduce(k27m_raw,   ignore.strand = TRUE)
k27mko_gr     <- reduce(k27mko_raw, ignore.strand = TRUE)

cat(sprintf("K27M_DMSO   after reduce : %d regions\n", length(nucleation_gr)))
cat(sprintf("K27MKO_DMSO after reduce : %d regions\n", length(k27mko_gr)))

# setdiff() is the GRanges equivalent of bedtools subtract:
# returns portions of K27MKO that do NOT overlap nucleation sites.
spreading_raw <- setdiff(k27mko_gr, nucleation_gr, ignore.strand = TRUE)
cat(sprintf("Spreading (pre-filter): %d regions\n", length(spreading_raw)))

# Discard sub-200 bp fragments produced by edge-trimming
MIN_WIDTH    <- 200L
spreading_gr <- spreading_raw[width(spreading_raw) >= MIN_WIDTH]
cat(sprintf("Spreading (>= %d bp) : %d regions\n", MIN_WIDTH, length(spreading_gr)))

# =============================================================================
# 6. Sanity check: zero overlap expected between the two final sets
# =============================================================================
cat("\n=== Sanity check ===\n")
n_overlap <- length(findOverlaps(nucleation_gr, spreading_gr))
if (n_overlap > 0L) {
  warning(sprintf(
    "UNEXPECTED: %d overlaps found between nucleation and spreading sites!",
    n_overlap
  ))
} else {
  cat("OK — 0 overlaps between nucleation and spreading sites.\n")
}

# =============================================================================
# 7. Summary statistics (GRanges only, no annotation packages needed)
# =============================================================================
summarise_gr <- function(gr, label) {
  w <- width(gr)
  cat(sprintf(
    "\n%s\n  Regions      : %d\n  Total (Mb)   : %.2f\n  Median width : %d bp\n  Width range  : %d – %d bp\n",
    label, length(gr), sum(w) / 1e6, as.integer(median(w)), min(w), max(w)
  ))
}

cat("\n=== Summary ===")
summarise_gr(nucleation_gr, "Nucleation sites (K27M_DMSO H3K27me3)")
summarise_gr(spreading_gr,  "Spreading sites  (K27MKO_DMSO minus nucleation)")

# =============================================================================
# 8. Export BED files (0-based)
# =============================================================================
cat("\n=== Exporting BED files ===\n")

export_bed_0based <- function(gr, path) {
  df <- as.data.frame(gr)[, c("seqnames", "start", "end")]
  colnames(df) <- c("chr", "start", "end")
  df$start <- df$start - 1L          # GRanges is 1-based; BED is 0-based
  write.table(df, file = path, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
  cat(sprintf("  Written: %s  (%d regions, %.2f Mb)\n",
              basename(path), nrow(df), sum(df$end - df$start) / 1e6))
}

nucleation_bed <- file.path(out_dir, "nucleation_sites.bed")
spreading_bed  <- file.path(out_dir, "spreading_sites.bed")

export_bed_0based(nucleation_gr, nucleation_bed)
export_bed_0based(spreading_gr,  spreading_bed)

# =============================================================================
# 9. Done
# =============================================================================
cat(sprintf("\n========================================\n"))
cat(sprintf("  step_7a complete  |  caller: %s\n", toupper(CALLER)))
cat(sprintf("========================================\n"))
cat("\nOutputs:\n")
cat(sprintf("  %s\n", nucleation_bed))
cat(sprintf("  %s\n", spreading_bed))
cat(sprintf("\nDone.\n"))
