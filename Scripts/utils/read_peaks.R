#!/usr/bin/env Rscript

# =============================================================================
# utils/read_peaks.R
#
# Helper functions to import peak files from different callers into a
# normalized GRanges object with a common metadata schema, regardless of the
# format differences between SEACR, MACS3, and csaw.
#
# Public API:
#   load_caller_peaks(peak_dir, caller, group, macs_broad = TRUE)
#     -> GRanges with mcols: $score, $peak_name, $caller
#
# Supported callers: "seacr" | "macs" | "csaw"
#
# Column conventions:
#   SEACR  (0-based 6-col): chr, start, end, total_signal, max_signal, max_signal_region
#   MACS3  broadPeak (0-based 9-col): chr, start, end, name, score, strand,
#                                     signalValue, pValue, qValue
#   MACS3  narrowPeak (0-based 10-col): same + peak (summit offset)
#   csaw   (0-based 6-col BED): chr, start, end, name, score, strand
# =============================================================================

suppressPackageStartupMessages({
  library(GenomicRanges)
})

# ---------------------------------------------------------------------------
# Internal: SEACR BED reader
#   score = total_signal (total AUC of the fragment pile-up)
# ---------------------------------------------------------------------------
.read_seacr_bed <- function(path) {
  if (!file.exists(path)) stop("SEACR file not found: ", path)
  df <- read.table(path, header = FALSE, sep = "\t",
                   col.names = c("chr", "start", "end",
                                 "total_signal", "max_signal", "max_signal_region"),
                   stringsAsFactors = FALSE)
  gr <- makeGRangesFromDataFrame(df,
    seqnames.field          = "chr",
    start.field             = "start",
    end.field               = "end",
    starts.in.df.are.0based = TRUE,
    keep.extra.columns      = FALSE)
  gr$score     <- df$total_signal
  gr$peak_name <- paste0("seacr_", seq_along(gr))
  gr
}

# ---------------------------------------------------------------------------
# Internal: MACS3 broadPeak / narrowPeak reader
#   score = signalValue (continuous fold-enrichment; more informative
#           than the integer 'score' column which is -10*log10(qvalue))
# ---------------------------------------------------------------------------
.read_macs_bed <- function(path, broad = TRUE) {
  if (!file.exists(path)) stop("MACS3 peak file not found: ", path)
  col_names <- c("chr", "start", "end", "name", "score_int", "strand",
                 "signalValue", "pValue", "qValue")
  if (!broad) col_names <- c(col_names, "peak_offset")
  df <- read.table(path, header = FALSE, sep = "\t",
                   col.names   = col_names,
                   fill        = TRUE,   # broadPeak may omit trailing columns
                   stringsAsFactors = FALSE)
  gr <- makeGRangesFromDataFrame(df,
    seqnames.field          = "chr",
    start.field             = "start",
    end.field               = "end",
    strand.field            = "strand",
    starts.in.df.are.0based = TRUE,
    keep.extra.columns      = FALSE)
  gr$score     <- df$signalValue
  gr$peak_name <- df$name
  gr
}

# ---------------------------------------------------------------------------
# Internal: csaw BED reader
#   score = -10*log10(FDR) stored as integer in column 5
#   (exported by callpeaks_csaw.R via rtracklayer::export.bed)
# ---------------------------------------------------------------------------
.read_csaw_bed <- function(path) {
  if (!file.exists(path)) stop("csaw BED file not found: ", path)
  df <- read.table(path, header = FALSE, sep = "\t",
                   col.names = c("chr", "start", "end", "name", "score", "strand"),
                   stringsAsFactors = FALSE)
  gr <- makeGRangesFromDataFrame(df,
    seqnames.field          = "chr",
    start.field             = "start",
    end.field               = "end",
    strand.field            = "strand",
    starts.in.df.are.0based = TRUE,
    keep.extra.columns      = FALSE)
  gr$score     <- df$score
  gr$peak_name <- df$name
  gr
}

# ---------------------------------------------------------------------------
# Internal: locate peak files for a given (caller, group) combination.
#
# Matching strategy:
#   Files must begin with <group> followed immediately by '_' or '.' to
#   avoid matching "K27M" when searching for "K27MKO" and vice-versa.
#   All patterns anchor on the antibody label "K27me3" for specificity.
#
# Returns: character vector of matched file paths (>= 1 expected)
# ---------------------------------------------------------------------------
.find_peak_files <- function(peak_dir, caller, group, macs_broad = TRUE) {
  if (!dir.exists(peak_dir)) {
    stop(sprintf("Peak directory does not exist: %s", peak_dir))
  }

  # Escape any regex metacharacters in the group name
  safe_group <- gsub("([.|()\\^{}+$*?])", "\\\\\\1", group)

  pattern <- switch(caller,
    seacr = paste0("^", safe_group, "[_.].*K27me3.*\\.relaxed\\.bed$"),
    macs  = {
      ext <- if (macs_broad) "broadPeak" else "narrowPeak"
      paste0("^", safe_group, "[_.].*K27me3.*_peaks\\.", ext, "$")
    },
    csaw  = paste0("^", safe_group, "[_.].*K27me3.*\\.tmm\\.bed$"),
    stop("Unknown caller '", caller, "'. Must be one of: macs, csaw, seacr")
  )

  files <- list.files(peak_dir, pattern = pattern, full.names = TRUE)

  if (length(files) == 0L) {
    stop(sprintf(
      paste0(
        "No %s peak files found for group '%s'\n",
        "  Directory : %s\n",
        "  Pattern   : %s\n",
        "  Files present:\n    %s"
      ),
      toupper(caller), group,
      peak_dir, pattern,
      paste(list.files(peak_dir), collapse = "\n    ")
    ))
  }

  message(sprintf("  Found %d file(s) for group '%s' [%s]:",
                  length(files), group, toupper(caller)))
  for (f in files) message("    ", basename(f))
  files
}

# ---------------------------------------------------------------------------
# Public: load_caller_peaks()
#
# Finds all peak files for <group> in <peak_dir> using <caller>-specific
# file patterns, reads each one with the appropriate parser, concatenates
# if multiple replicates are present, and returns a single GRanges.
#
# Args:
#   peak_dir   - path to caller-specific peak directory
#                (e.g. config.PEAKDIR + "/csaw/")
#   caller     - one of "seacr", "macs", "csaw"
#   group      - sample group label (e.g. "K27M" or "K27MKO")
#   macs_broad - if TRUE (default), look for broadPeak; else narrowPeak.
#                Ignored for seacr and csaw.
#
# Returns GRanges with mcols:
#   $score      - numeric peak score (caller-specific units)
#   $peak_name  - character peak identifier
#   $caller     - character caller label
# ---------------------------------------------------------------------------
load_caller_peaks <- function(peak_dir, caller, group, macs_broad = TRUE) {
  caller <- tolower(trimws(caller))
  files  <- .find_peak_files(peak_dir, caller, group, macs_broad)

  gr_list <- lapply(files, function(f) {
    message("  Reading: ", basename(f))
    gr <- switch(caller,
      seacr = .read_seacr_bed(f),
      macs  = .read_macs_bed(f, broad = macs_broad),
      csaw  = .read_csaw_bed(f)
    )
    gr$caller <- caller
    gr
  })

  if (length(gr_list) == 1L) {
    gr_list[[1L]]
  } else {
    message(sprintf("  Concatenating %d replicate GRanges objects ...", length(gr_list)))
    do.call(c, gr_list)
  }
}
