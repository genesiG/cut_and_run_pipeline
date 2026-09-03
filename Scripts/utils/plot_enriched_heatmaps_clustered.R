#!/usr/bin/env Rscript

# =============================================================================
# utils/plot_enriched_heatmaps_clustered.R
#
# Generalised helper: render EnrichedHeatmaps with k-means cluster row-split
# from pre-built normalizedMatrix RDS files (produced by utils/build_matrix.R)
# and cluster assignment CSV (produced by step_8e_cbx2peaks_clustering.py).
#
# All configuration is read from a single JSON file passed on the CLI.
# This avoids shell-quoting issues with multi-line strings and makes the
# parameters fully inspectable and reproducible.
#
# Usable in two ways:
#   1. Standalone CLI:
#        Rscript utils/plot_enriched_heatmaps_clustered.R /path/to/config.json
#
#   2. Sourced from another script:
#        source("utils/plot_enriched_heatmaps_clustered.R")
#        plot_enriched_heatmaps_clustered_main(cfg)
#
# Required JSON fields:
#   norm_mats_rds          path to {prefix}_norm_mats.rds
#   targets_rds            path to {prefix}_targets.rds
#   cluster_assignments_csv path to cluster_assignments_k{K}.csv
#   window_bp              half-window in bp used during build_matrix
#   out_dir                output directory
#   out_prefix             prefix for output files
#   sample_labels          ordered list: internal names (must match norm_mats keys)
#   sample_display_labels  display labels for plot column titles
#   sample_hues            hex colour per sample (same count as sample_labels)
#   cluster_labels         human-readable label per cluster (k elements)
#   n_clusters             integer k
#
# Optional JSON fields (defaults shown):
#   sort_by_sample         sample to sort rows by within each cluster
#                          (default: first element of sample_labels)
#   z_min                  array of per-sample min values, or null (data-driven)
#   z_max                  array of per-sample max values, or null (data-driven)
#   out_format             "svg" or "pdf" (default: "svg")
#   heatmap_width_cm       width per column in cm (default: 3.5)
#   column_gap_cm          gap between columns in cm (default: 1.0)
#   anno_height_cm         profile annotation height in cm (default: 2.5)
#   split_with_gaps        TRUE/FALSE whether to draw gaps between row clusters
#                          (default: TRUE)
#
# Outputs (all written to out_dir/):
#   {out_prefix}_kmeans_k{K}_heatmap.{format}   main heatmap
#   regions/{out_prefix}_cluster{c}.bed         per-cluster BED files
# =============================================================================

suppressPackageStartupMessages({
  library(jsonlite)
  library(EnrichedHeatmap)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(GenomicRanges)
  library(grid)
})

# NULL coalescing operator — defined here so it is available both in the
# main function body and when the script is sourced by another script.
`%||%` <- function(a, b) if (!is.null(a)) a else b

# =============================================================================
# Core function (called both from CLI and when sourced)
# =============================================================================
plot_enriched_heatmaps_clustered_main <- function(cfg) {

  # ---------------------------------------------------------------------------
  # Unpack and validate config
  # ---------------------------------------------------------------------------
  required_fields <- c("norm_mats_rds", "targets_rds", "cluster_assignments_csv",
                       "window_bp", "out_dir", "out_prefix",
                       "sample_labels", "sample_display_labels", "sample_hues",
                       "cluster_labels", "n_clusters")
  missing_fields <- setdiff(required_fields, names(cfg))
  if (length(missing_fields) > 0L)
    stop("Missing required JSON fields: ", paste(missing_fields, collapse = ", "))

  norm_mats_rds          <- cfg$norm_mats_rds
  targets_rds            <- cfg$targets_rds
  cluster_assignments_csv <- cfg$cluster_assignments_csv
  window_bp              <- as.integer(cfg$window_bp)
  out_dir                <- cfg$out_dir
  out_prefix             <- cfg$out_prefix
  sample_labels          <- cfg$sample_labels
  sample_display_labels  <- cfg$sample_display_labels
  sample_hues            <- cfg$sample_hues
  cluster_labels         <- cfg$cluster_labels
  n_clusters             <- as.integer(cfg$n_clusters)

  # Optional fields with defaults
  sort_by_sample    <- cfg$sort_by_sample  %||% sample_labels[[1L]]
  z_min             <- cfg$z_min           # NULL = data-driven
  z_max             <- cfg$z_max           # NULL = data-driven
  out_format        <- tolower(cfg$out_format      %||% "svg")
  heatmap_width     <- as.numeric(cfg$heatmap_width_cm  %||% 3.5)
  column_gap        <- as.numeric(cfg$column_gap_cm     %||% 1.0)
  anno_height       <- as.numeric(cfg$anno_height_cm    %||% 2.5)
  split_with_gaps   <- isTRUE(cfg$split_with_gaps %||% TRUE)
  profile_linewidth <- as.numeric(cfg$profile_linewidth %||% 1.5)

  # Validate files
  for (f in c(norm_mats_rds, targets_rds, cluster_assignments_csv))
    if (!file.exists(f)) stop("File not found: ", f)
  if (length(sample_labels) != length(sample_hues))
    stop("sample_labels and sample_hues must have equal length")
  if (length(sample_labels) != length(sample_display_labels))
    stop("sample_labels and sample_display_labels must have equal length")
  if (length(cluster_labels) != n_clusters)
    stop("cluster_labels must have exactly n_clusters = ", n_clusters, " elements")

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  regions_dir <- file.path(out_dir, "regions")
  if (!dir.exists(regions_dir)) dir.create(regions_dir, recursive = TRUE)

  cat("=== plot_enriched_heatmaps_clustered.R ===\n")
  cat(sprintf("  Prefix     : %s\n", out_prefix))
  cat(sprintf("  Window     : +/- %d bp\n", window_bp))
  cat(sprintf("  Samples    : %s\n", paste(sample_labels, collapse = ", ")))
  cat(sprintf("  k clusters : %d  (%s)\n", n_clusters,
              paste(cluster_labels, collapse = " | ")))
  cat(sprintf("  Sort by    : %s\n", sort_by_sample))
  cat(sprintf("  Format     : %s\n", out_format))

  # ---------------------------------------------------------------------------
  # Load data
  # ---------------------------------------------------------------------------
  cat("Loading data...\n")
  norm_mats_raw   <- readRDS(norm_mats_rds)
  pp_targets      <- readRDS(targets_rds)
  cluster_ids_raw <- as.integer(read.csv(cluster_assignments_csv, header = FALSE)[[1L]])

  # Subset and order norm_mats to the requested sample order
  missing_s <- setdiff(sample_labels, names(norm_mats_raw))
  if (length(missing_s) > 0L)
    stop("sample_labels not found in norm_mats: ", paste(missing_s, collapse = ", "))
  norm_mats <- norm_mats_raw[sample_labels]

  # Validate cluster assignment length matches regions
  n_regions <- nrow(norm_mats[[1L]])
  if (length(cluster_ids_raw) != n_regions)
    stop(sprintf("cluster_assignments has %d rows but norm_mats has %d rows",
                 length(cluster_ids_raw), n_regions))

  # Map cluster integer IDs to human-readable labels
  cluster_ids <- cluster_ids_raw  # 1-indexed integers

  # ---------------------------------------------------------------------------
  # Color scales: fixed (from JSON) or data-driven (max column mean)
  # ---------------------------------------------------------------------------
  window_kb   <- window_bp / 1000
  axis_labels <- c(paste0("-", window_kb, "kb"), "center", paste0("+", window_kb, "kb"))

  get_valid_val <- function(vec, idx) {
    if (is.null(vec) || length(vec) < idx) return(NULL)
    v <- vec[[idx]]
    if (is.null(v) || is.na(v) || is.nan(v)) return(NULL)
    as.numeric(v)
  }

  make_color_fun <- function(i) {
    hue  <- sample_hues[[i]]
    cmin <- get_valid_val(z_min, i) %||% 0
    cmax <- get_valid_val(z_max, i) %||% max(colMeans(as.matrix(norm_mats[[i]]), na.rm = TRUE), 1e-9)
    colorRamp2(c(cmin, cmax), c("white", hue))
  }

  sample_colors <- setNames(
    lapply(seq_along(sample_labels), make_color_fun),
    sample_labels
  )

  # Report actual color ranges
  for (i in seq_along(sample_labels)) {
    rng <- attr(sample_colors[[sample_labels[i]]], "breaks")
    cat(sprintf("  Color scale [%s]: %.3f - %.3f  hue=%s\n",
                sample_labels[[i]], rng[1L], rng[2L], sample_hues[[i]]))
  }

  # ---------------------------------------------------------------------------
  # Profile annotation ylim: per-cluster worst-case max column mean
  # (prevents profile lines from overflowing the annotation box)
  # ---------------------------------------------------------------------------
  cluster_profile_ymax <- function(nm, ids, headroom = 1.08) {
    mat <- as.matrix(nm)
    max(sapply(sort(unique(ids)), function(cl) {
      max(colMeans(mat[ids == cl, , drop = FALSE], na.rm = TRUE))
    })) * headroom
  }

  profile_ylims <- setNames(
    lapply(seq_along(sample_labels), function(i) {
      s <- sample_labels[[i]]
      ymax <- get_valid_val(z_max, i) %||% cluster_profile_ymax(norm_mats[[s]], cluster_ids)
      c(0, max(ymax, 1e-9))
    }),
    sample_labels
  )

  # ---------------------------------------------------------------------------
  # Row split and order
  # Row split: named factor so cluster labels appear on the heatmap
  # Row order within each cluster: decreasing mean of sort_by_sample
  # ---------------------------------------------------------------------------
  group_levels <- paste0("Cluster ", seq_len(n_clusters))
  named_labels <- cluster_labels  # user-provided, length == n_clusters
  level_to_label <- setNames(named_labels, group_levels)

  row_split <- factor(
    level_to_label[paste0("Cluster ", cluster_ids)],
    levels = named_labels
  )

  # Sort by decreasing mean of sort_by_sample within each cluster
  sort_mat <- as.matrix(norm_mats[[sort_by_sample]])
  row_means <- rowMeans(sort_mat, na.rm = TRUE)
  row_order_vec <- order(cluster_ids, -row_means)

  # ---------------------------------------------------------------------------
  # Cluster colours (for the profile lines legend)
  # ---------------------------------------------------------------------------
  if (!is.null(cfg$cluster_colors) && length(cfg$cluster_colors) >= n_clusters) {
    raw_cols <- as.character(cfg$cluster_colors)[seq_len(n_clusters)]
  } else {
    raw_cols <- brewer.pal(max(3L, n_clusters), "Set1")[seq_len(n_clusters)]
  }
  cluster_colors <- setNames(raw_cols, named_labels)

  # ---------------------------------------------------------------------------
  # Build legends (all packed into one horizontal row at the bottom)
  # ---------------------------------------------------------------------------
  font_legend_title <- as.numeric(cfg$font_legend_title %||% 9)
  font_legend_text  <- as.numeric(cfg$font_legend_text  %||% 8)
  font_col_title    <- as.numeric(cfg$font_col_title    %||% 11)
  font_axis_name    <- as.numeric(cfg$font_axis_name    %||% 8)
  font_anno_axis    <- as.numeric(cfg$font_anno_axis    %||% 7)
  font_row_title    <- as.numeric(cfg$font_row_title    %||% 9)

  lgd_samples <- lapply(sample_labels, function(s) {
    Legend(
      col_fun    = sample_colors[[s]],
      title      = sample_display_labels[[which(sample_labels == s)]],
      direction  = "horizontal",
      title_gp   = gpar(fontsize = font_legend_title, fontface = "bold"),
      labels_gp  = gpar(fontsize = font_legend_text)
    )
  })
  lgd_cluster <- Legend(
    title     = "Cluster profile",
    labels    = named_labels,
    legend_gp = gpar(col = cluster_colors, lwd = profile_linewidth),
    type      = "lines",
    direction = "horizontal",
    title_gp  = gpar(fontsize = font_legend_title, fontface = "bold"),
    labels_gp = gpar(fontsize = font_legend_text)
  )
  combined_legend <- do.call(
    packLegend,
    c(lgd_samples, list(lgd_cluster),
      list(direction = "horizontal", gap = unit(4, "mm")))
  )

  # ---------------------------------------------------------------------------
  # Build heatmap list
  # ---------------------------------------------------------------------------
  row_gap      <- if (split_with_gaps) unit(2, "mm") else unit(0, "mm")
  border_style <- split_with_gaps

  ht_list <- NULL
  for (i in seq_along(sample_labels)) {
    s     <- sample_labels[[i]]
    disp  <- sample_display_labels[[i]]

    top_anno <- HeatmapAnnotation(
      enriched = anno_enriched(
        gp         = gpar(col = cluster_colors, lwd = profile_linewidth),
        ylim       = profile_ylims[[s]],
        axis_param = list(gp = gpar(fontsize = font_anno_axis))
      ),
      height = unit(anno_height, "cm")
    )

    ht <- EnrichedHeatmap(
      norm_mats[[s]],
      col                 = sample_colors[[s]],
      name                = s,
      width               = unit(heatmap_width, "cm"),
      column_title        = disp,
      column_title_gp     = gpar(fontsize = font_col_title, fontface = "bold"),
      axis_name           = axis_labels,
      axis_name_gp        = gpar(fontsize = font_axis_name),
      row_order           = row_order_vec,
      row_split           = row_split,
      row_title_side      = "left",
      row_title_rot       = 90,
      row_title_gp        = gpar(fontsize = font_row_title, fontface = "bold"),
      row_gap             = row_gap,
      border              = border_style,
      top_annotation      = top_anno,
      pos_line            = FALSE,
      use_raster          = TRUE,
      show_heatmap_legend = FALSE   # built manually above
    )
    ht_list <- if (is.null(ht_list)) ht else ht_list + ht
  }

  # ---------------------------------------------------------------------------
  # Save heatmap
  # ---------------------------------------------------------------------------
  out_file <- file.path(out_dir,
                        sprintf("%s_kmeans_k%d_heatmap.%s",
                                out_prefix, n_clusters, out_format))
  total_width_cm <- (length(sample_labels) * heatmap_width) +
                    ((max(1L, length(sample_labels)) - 1L) * column_gap) + 5.0
  svg_w <- total_width_cm / 2.54
  svg_h <- 10

  cat(sprintf("  Saving heatmap: %s (width=%.2f in, height=%.2f in)\n",
              basename(out_file), svg_w, svg_h))

  if (out_format == "svg") {
    svg(out_file, width = svg_w, height = svg_h)
  } else {
    pdf(out_file, width = svg_w, height = svg_h)
  }

  tryCatch(
    draw(ht_list,
         gap                    = unit(column_gap, "cm"),
         annotation_legend_list = list(combined_legend),
         annotation_legend_side = "bottom",
         merge_legend           = FALSE,
         padding                = unit(c(5, 5, 5, 5), "mm")),
    error = function(e) message("  WARNING: draw() failed: ", conditionMessage(e))
  )
  dev.off()

  # ---------------------------------------------------------------------------
  # Export per-cluster BED files
  # ---------------------------------------------------------------------------
  bed_meta           <- as.data.frame(pp_targets)[, c("seqnames", "start", "end")]
  colnames(bed_meta) <- c("chr", "start", "end")
  bed_meta$name      <- names(pp_targets)

  if (!is.null(cfg$regions_bed) && nzchar(cfg$regions_bed) && file.exists(cfg$regions_bed)) {
    orig_bed <- read.table(cfg$regions_bed, header = FALSE, sep = "\t",
                           stringsAsFactors = FALSE, comment.char = "#")
    key_orig <- paste(orig_bed[, 1L], orig_bed[, 2L], orig_bed[, 3L], sep = ":")
    key_targ <- paste(bed_meta$chr, bed_meta$start - 1L, bed_meta$end, sep = ":")
    m_idx    <- match(key_targ, key_orig)
    if (all(!is.na(m_idx))) {
      bed_meta <- orig_bed[m_idx, , drop = FALSE]
      cat(sprintf("  Matched %d regions with original peak file: %s\n",
                  nrow(bed_meta), basename(cfg$regions_bed)))
    }
  }

  for (cl in seq_len(n_clusters)) {
    cl_label <- named_labels[[cl]]
    cl_bed   <- bed_meta[cluster_ids == cl, , drop = FALSE]
    out_bed  <- file.path(regions_dir,
                          sprintf("%s_k%d_%s.bed",
                                  out_prefix, n_clusters,
                                  gsub("\\s+", "_", cl_label)))
    write.table(cl_bed, file = out_bed, sep = "\t", quote = FALSE,
                row.names = FALSE, col.names = FALSE)
    cat(sprintf("  Cluster '%s': %d regions -> %s\n",
                cl_label, nrow(cl_bed), basename(out_bed)))
  }

  cat(sprintf("=== plot_enriched_heatmaps_clustered.R complete (%s) ===\n\n",
              out_prefix))
  invisible(out_file)
}

# =============================================================================
# CLI entry-point (only active when run as a standalone script)
# =============================================================================
if (sys.nframe() == 0L && !interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 1L)
    stop("Usage: Rscript utils/plot_enriched_heatmaps_clustered.R <config.json>")

  cfg_file <- args[1L]
  if (!file.exists(cfg_file))
    stop("Config JSON not found: ", cfg_file)

  cfg <- jsonlite::fromJSON(cfg_file, simplifyVector = TRUE)
  plot_enriched_heatmaps_clustered_main(cfg)
}
