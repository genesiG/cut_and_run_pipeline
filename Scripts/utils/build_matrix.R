#!/usr/bin/env Rscript

# =============================================================================
# utils/build_matrix.R
#
# Generalised helper: build EnrichedHeatmap normalizedMatrix objects from
# bigWig signal over a set of BED regions.
#
# Usable in two ways:
#   1. Standalone CLI:
#        Rscript utils/build_matrix.R \
#          --regions_bed  tss_windows.bed \
#          --bigwig_files a.bw b.bw \
#          --bigwig_labels DMSO UNC1999 \
#          --window_bp    3000 \
#          --out_dir      /path/to/out \
#          --out_prefix   tss
#
#   2. Sourced from another script:
#        source("utils/build_matrix.R")
#        build_matrix_main(regions_bed = ..., ...)
#
# CLI arguments:
#   --regions_bed     0-based 3-col BED defining input regions
#   --bigwig_files    1+ bigWig paths (space-separated)
#   --bigwig_labels   matching labels (same count as --bigwig_files)
#   --window_bp       half-window in bp (e.g. 3000 for ±3 kb)
#   --n_bins          bins per half-window (default: 100)
#   --out_dir         output directory
#   --out_prefix      filename prefix (e.g. "tss" or "genebody")
#   --body_scale      flag (no value) — if present, keep full region widths
#                     for gene-body scaling (include_target = TRUE)
#
# Outputs (in out_dir/):
#   {prefix}_norm_mats.rds        named list of normalizedMatrix objects
#   {prefix}_targets.rds          original GRanges
#   {prefix}_targets_center.rds   1-bp centered GRanges (or original if body_scale)
# =============================================================================

suppressPackageStartupMessages({
  library(rtracklayer)
  library(EnrichedHeatmap)
  library(GenomicRanges)
  library(parallel)     # mclapply — fork-based true parallelism on Linux
})

# =============================================================================
# Core function (called both from CLI and when sourced)
# =============================================================================
build_matrix_main <- function(regions_bed,
                               bigwig_files,
                               bigwig_labels,
                               window_bp,
                               n_bins             = 100L,
                               out_dir,
                               out_prefix,
                               body_scale         = FALSE,
                               clustering_samples = character(0)) {
  # How many parallel worker processes to use.
  # Set MC_CORES in the shell environment (e.g. in the #BSUB batch script)
  # before calling Rscript.  Defaults to 1 (serial) if not set.
  # On Linux, mclapply() uses fork() — true OS-level parallelism.
  # On Windows, mclapply() silently falls back to lapply().
  mc_cores <- as.integer(Sys.getenv("MC_CORES", unset = "1"))

  window_bp <- as.integer(window_bp)
  n_bins    <- as.integer(n_bins)

  if (!file.exists(regions_bed))
    stop("Regions BED not found: ", regions_bed)
  if (length(bigwig_files) != length(bigwig_labels))
    stop("--bigwig_files and --bigwig_labels must have equal length")
  for (bw in bigwig_files)
    if (!file.exists(bw)) stop("bigWig not found: ", bw)

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  cat("=== build_matrix.R ===\n")
  cat(sprintf("  Prefix     : %s\n", out_prefix))
  cat(sprintf("  Regions    : %s\n", regions_bed))
  cat(sprintf("  Window     : +/- %d bp   Bins per side: %d\n", window_bp, n_bins))
  cat(sprintf("  Body scale : %s\n", body_scale))
  cat(sprintf("  Samples    : %s\n", paste(bigwig_labels, collapse = ", ")))

  # ---------------------------------------------------------------------------
  # Load regions
  # ---------------------------------------------------------------------------
  bed_raw <- read.table(regions_bed, header = FALSE, sep = "\t",
                        stringsAsFactors = FALSE, comment.char = "#")
  bed <- bed_raw[, 1:min(4L, ncol(bed_raw)), drop = FALSE]
  if (ncol(bed) >= 4L) {
    colnames(bed)[1:4] <- c("chr", "start", "end", "name")
  } else {
    colnames(bed)[1:3] <- c("chr", "start", "end")
    bed$name <- paste0("region_", seq_len(nrow(bed)))
  }

  targets <- makeGRangesFromDataFrame(bed,
                seqnames.field          = "chr",
                start.field             = "start",
                end.field               = "end",
                keep.extra.columns      = TRUE,
                starts.in.df.are.0based = TRUE)
  names(targets) <- targets$name
  cat(sprintf("  Regions loaded: %d\n", length(targets)))

  # ---------------------------------------------------------------------------
  # Reference points for normalizeToMatrix
  # ---------------------------------------------------------------------------
  if (!body_scale) {
    targets_ref <- resize(targets, width = 1L, fix = "center")
    names(targets_ref) <- names(targets)
    import_windows <- resize(targets_ref, width = 2L * window_bp, fix = "center")
  } else {
    targets_ref    <- targets
    import_windows <- resize(targets, width = width(targets) + 2L * window_bp,
                             fix = "center")
  }

  # ---------------------------------------------------------------------------
  # Build normalizedMatrix per sample
  # ---------------------------------------------------------------------------
  build_nm <- function(label, bw_path) {
    cat(sprintf("  [%s] Importing bigWig ...\n", label))
    bw_signal <- import(bw_path,
                        format    = "BigWig",
                        selection = BigWigSelection(import_windows))
    cat(sprintf("  [%s] Building normalizedMatrix ...\n", label))
    nm <- normalizeToMatrix(
      signal         = bw_signal,
      target         = targets_ref,
      extend         = window_bp,
      w              = window_bp / n_bins,
      value_column   = "score",
      mean_mode      = "w0",
      include_target = body_scale,
      smooth         = FALSE
    )
    cat(sprintf("  [%s] Dimensions: %d x %d\n", label, nrow(nm), ncol(nm)))
    nm
  }

  cat(sprintf("  Parallelism : MC_CORES=%d (mclapply fork)\n", mc_cores))

  # ---------------------------------------------------------------------------
  # Build normalizedMatrix per sample — true parallel on Linux via fork().
  # Each child process independently imports one bigWig and runs
  # normalizeToMatrix(). Shared read-only objects (import_windows, targets_ref,
  # etc.) are inherited via copy-on-write, not duplicated in full.
  # ---------------------------------------------------------------------------
  raw_results <- mclapply(
    seq_along(bigwig_labels),
    function(i) build_nm(bigwig_labels[i], bigwig_files[i]),
    mc.cores  = mc_cores,
    mc.silent = FALSE
  )

  # Detect any child-process failures
  failed <- vapply(raw_results, inherits, logical(1), what = "try-error")
  if (any(failed)) {
    stop("mclapply: the following samples failed:\n",
         paste(bigwig_labels[failed], collapse = "\n"))
  }

  norm_mats <- setNames(raw_results, bigwig_labels)

  # ---------------------------------------------------------------------------
  # Save normalizedMatrix RDS outputs
  # ---------------------------------------------------------------------------
  saveRDS(norm_mats,   file.path(out_dir, paste0(out_prefix, "_norm_mats.rds")))
  saveRDS(targets,     file.path(out_dir, paste0(out_prefix, "_targets.rds")))
  saveRDS(targets_ref, file.path(out_dir, paste0(out_prefix, "_targets_center.rds")))

  cat(sprintf("  Saved: %s_norm_mats.rds / _targets.rds / _targets_center.rds\n",
              out_prefix))

  # ---------------------------------------------------------------------------
  # Optional: export clustering matrix (subset of samples) as gzip TSV.
  # The Python step_8e script reads this with np.loadtxt() for k-means.
  # ---------------------------------------------------------------------------
  if (length(clustering_samples) > 0L) {
    missing_cs <- setdiff(clustering_samples, bigwig_labels)
    if (length(missing_cs) > 0L)
      stop("--clustering_samples not found in bigwig_labels: ",
           paste(missing_cs, collapse = ", "))

    cat(sprintf("  Building cluster matrix from: %s\n",
                paste(clustering_samples, collapse = " + ")))
    cluster_mat <- do.call(
      cbind,
      lapply(clustering_samples, function(s) as.matrix(norm_mats[[s]]))
    )
    cluster_mat[is.na(cluster_mat)] <- 0L
    cat(sprintf("  Cluster matrix dimensions: %d x %d\n",
                nrow(cluster_mat), ncol(cluster_mat)))

    gz_path <- file.path(out_dir, paste0(out_prefix, "_cluster_matrix.tsv.gz"))
    gz_con  <- gzcon(file(gz_path, "wb"))
    write.table(cluster_mat, gz_con,
                sep = "\t", row.names = FALSE, col.names = FALSE)
    close(gz_con)
    cat(sprintf("  Cluster matrix saved: %s\n", gz_path))
  }

  cat(sprintf("=== build_matrix.R complete (%s) ===\n\n", out_prefix))
  invisible(norm_mats)
}

# =============================================================================
# CLI entry-point (only active when run as a standalone script)
# =============================================================================
if (sys.nframe() == 0L && !interactive()) {

  .get_flag_values <- function(flag, args, required = TRUE) {
    idx <- which(args == flag)
    if (length(idx) == 0L) {
      if (required) stop("Missing required argument: ", flag)
      return(character(0))
    }
    start <- idx[1L] + 1L
    if (start > length(args)) stop("No value supplied for: ", flag)
    end <- start
    while (end <= length(args) && !startsWith(args[end], "--")) end <- end + 1L
    args[start:(end - 1L)]
  }

  .get_flag_value <- function(flag, args, required = TRUE, default = NULL) {
    vals <- .get_flag_values(flag, args, required = FALSE)
    if (length(vals) == 0L) {
      if (required) stop("Missing required argument: ", flag)
      return(default)
    }
    vals[1L]
  }

  .args <- commandArgs(trailingOnly = TRUE)

  build_matrix_main(
    regions_bed        = .get_flag_value("--regions_bed",   .args),
    bigwig_files       = .get_flag_values("--bigwig_files",  .args),
    bigwig_labels      = .get_flag_values("--bigwig_labels", .args),
    window_bp          = as.integer(.get_flag_value("--window_bp", .args)),
    n_bins             = as.integer(.get_flag_value("--n_bins", .args, required = FALSE, default = "100")),
    out_dir            = .get_flag_value("--out_dir",   .args),
    out_prefix         = .get_flag_value("--out_prefix", .args),
    body_scale         = "--body_scale" %in% .args,
    clustering_samples = .get_flag_values("--clustering_samples", .args, required = FALSE)
  )
}
