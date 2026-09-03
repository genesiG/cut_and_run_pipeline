#!/usr/bin/env Rscript

# =============================================================================
# utils/plot_enriched_heatmaps.R
#
# Generalised helper: render EnrichedHeatmaps from pre-built normalizedMatrix
# RDS files (produced by utils/build_matrix.R).
#
# No clustering row-split — rows are sorted by decreasing mean signal.
# Generalised to accept any number of samples with arbitrary labels/colours.
#
# Usable in two ways:
#   1. Standalone CLI:
#        Rscript utils/plot_enriched_heatmaps.R \
#          --norm_mats_rds tss_norm_mats.rds \
#          --targets_rds   tss_targets.rds \
#          --window_bp     3000 \
#          --out_dir       /path/to/out \
#          --out_prefix    01a_tss_enriched \
#          --sample_labels DMSO UNC1999 K27M \
#          --sample_hues   "#88419d" "#58135e" "#39107b"
#
#   2. Sourced from another script:
#        source("utils/plot_enriched_heatmaps.R")
#        plot_enriched_heatmaps_main(norm_mats_rds = ..., ...)
#
# CLI arguments:
#   --norm_mats_rds   path to {prefix}_norm_mats.rds
#   --targets_rds     path to {prefix}_targets.rds
#   --window_bp       half-window used during build_matrix
#   --out_dir         output directory
#   --out_prefix      prefix for output file
#   --sample_labels   ordered display labels (must match names in norm_mats)
#   --sample_hues     hex color per sample (same count as --sample_labels)
#   --out_format      "pdf" or "svg" (default: pdf)
#   --heatmap_width   cm per heatmap column (default: 3.5)
#   --column_gap      cm gap between columns (default: 1.0)
#   --anno_height     cm profile annotation height (default: 2.5)
#
# Output: {out_dir}/{out_prefix}_heatmap.{format}
# =============================================================================

suppressPackageStartupMessages({
  library(EnrichedHeatmap)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

# =============================================================================
# Core function
# =============================================================================
plot_enriched_heatmaps_main <- function(norm_mats_rds,
                                         targets_rds,
                                         window_bp,
                                         out_dir,
                                         out_prefix,
                                         sample_labels,
                                         sample_hues,
                                         out_format     = "pdf",
                                         heatmap_width  = 3.5,
                                         column_gap     = 1.0,
                                         anno_height    = 2.5) {

  window_bp     <- as.integer(window_bp)
  heatmap_width <- as.numeric(heatmap_width)
  column_gap    <- as.numeric(column_gap)
  anno_height   <- as.numeric(anno_height)
  out_format    <- tolower(out_format)

  if (!file.exists(norm_mats_rds)) stop("norm_mats_rds not found: ", norm_mats_rds)
  if (!file.exists(targets_rds))   stop("targets_rds not found: ",   targets_rds)
  if (length(sample_labels) != length(sample_hues))
    stop("--sample_labels and --sample_hues must have equal length")

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  cat("=== plot_enriched_heatmaps.R ===\n")
  cat(sprintf("  Prefix  : %s\n", out_prefix))
  cat(sprintf("  Window  : +/- %d bp\n", window_bp))
  cat(sprintf("  Samples : %s\n", paste(sample_labels, collapse = ", ")))
  cat(sprintf("  Format  : %s\n", out_format))

  # ---------------------------------------------------------------------------
  # Load data
  # ---------------------------------------------------------------------------
  norm_mats_raw <- readRDS(norm_mats_rds)
  pp_targets    <- readRDS(targets_rds)

  # Subset and order to requested labels
  missing <- setdiff(sample_labels, names(norm_mats_raw))
  if (length(missing) > 0)
    stop("sample_labels not found in norm_mats: ", paste(missing, collapse = ", "))
  norm_mats <- norm_mats_raw[sample_labels]

  # ---------------------------------------------------------------------------
  # Axis labels
  # ---------------------------------------------------------------------------
  window_kb   <- window_bp / 1000
  axis_labels <- c(paste0("-", window_kb, "kb"), "center",
                   paste0("+", window_kb, "kb"))

  # ---------------------------------------------------------------------------
  # Data-driven color scale: cap = max column mean per sample
  # ---------------------------------------------------------------------------
  sample_colors <- setNames(
    lapply(seq_along(sample_labels), function(i) {
      nm  <- norm_mats[[sample_labels[i]]]
      cap <- max(colMeans(as.matrix(nm), na.rm = TRUE))
      cap <- max(cap, 1e-9)
      colorRamp2(c(0, cap), c("white", sample_hues[i]))
    }),
    sample_labels
  )

  # ---------------------------------------------------------------------------
  # Row order: decreasing mean signal across all samples
  # ---------------------------------------------------------------------------
  combined_mat <- do.call(cbind, lapply(norm_mats, as.matrix))
  row_order_vec <- order(-rowMeans(combined_mat, na.rm = TRUE))

  # ---------------------------------------------------------------------------
  # Legends
  # ---------------------------------------------------------------------------
  font_legend_title <- 9
  font_legend_text  <- 8
  font_col_title    <- 11
  font_axis_name    <- 8
  font_anno_axis    <- 7

  lgd_list <- lapply(sample_labels, function(s) {
    Legend(
      col_fun    = sample_colors[[s]],
      title      = sprintf("H3K27me3\n(%s)", s),
      direction  = "horizontal",
      title_gp   = gpar(fontsize = font_legend_title, fontface = "bold"),
      labels_gp  = gpar(fontsize = font_legend_text)
    )
  })
  combined_legend <- do.call(packLegend,
                             c(lgd_list, list(direction = "horizontal",
                                              gap = unit(4, "mm"))))

  # Profile ylim: max column mean with 8% headroom
  profile_ylims <- lapply(sample_labels, function(s) {
    ymax <- max(colMeans(as.matrix(norm_mats[[s]]), na.rm = TRUE)) * 1.08
    c(0, max(ymax, 1e-9))
  })
  names(profile_ylims) <- sample_labels

  # ---------------------------------------------------------------------------
  # Build heatmap list
  # ---------------------------------------------------------------------------
  ht_list <- NULL
  for (s in sample_labels) {
    top_anno <- HeatmapAnnotation(
      enriched = anno_enriched(
        gp         = gpar(col = sample_hues[sample_labels == s], lwd = 1.5),
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
      column_title        = s,
      column_title_gp     = gpar(fontsize = font_col_title, fontface = "bold"),
      axis_name           = axis_labels,
      axis_name_gp        = gpar(fontsize = font_axis_name),
      row_order           = row_order_vec,
      top_annotation      = top_anno,
      pos_line            = FALSE,
      use_raster          = TRUE,
      show_heatmap_legend = FALSE
    )
    ht_list <- if (is.null(ht_list)) ht else ht_list + ht
  }

  # ---------------------------------------------------------------------------
  # Save output
  # ---------------------------------------------------------------------------
  out_file <- file.path(out_dir, paste0(out_prefix, "_heatmap.", out_format))
  svg_w    <- heatmap_width / 2.54 * length(sample_labels) + 2
  svg_h    <- 10

  if (out_format == "svg") {
    svg(out_file, width = svg_w, height = svg_h)
  } else {
    pdf(out_file, width = svg_w, height = svg_h)
  }

  tryCatch({
    draw(ht_list,
         gap                    = unit(column_gap, "cm"),
         annotation_legend_list = list(combined_legend),
         annotation_legend_side = "bottom",
         merge_legend           = FALSE)
  }, error = function(e) {
    message("  WARNING: heatmap draw failed: ", conditionMessage(e))
  })
  dev.off()

  cat(sprintf("  Saved: %s\n", basename(out_file)))
  cat(sprintf("=== plot_enriched_heatmaps.R complete (%s) ===\n\n", out_prefix))

  invisible(out_file)
}

# =============================================================================
# CLI entry-point
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

  plot_enriched_heatmaps_main(
    norm_mats_rds = .get_flag_value("--norm_mats_rds", .args),
    targets_rds   = .get_flag_value("--targets_rds",   .args),
    window_bp     = as.integer(.get_flag_value("--window_bp", .args)),
    out_dir       = .get_flag_value("--out_dir",    .args),
    out_prefix    = .get_flag_value("--out_prefix", .args),
    sample_labels = .get_flag_values("--sample_labels", .args),
    sample_hues   = .get_flag_values("--sample_hues",   .args),
    out_format    = .get_flag_value("--out_format",   .args, required = FALSE, default = "pdf"),
    heatmap_width = as.numeric(.get_flag_value("--heatmap_width", .args, required = FALSE, default = "3.5")),
    column_gap    = as.numeric(.get_flag_value("--column_gap",    .args, required = FALSE, default = "1.0")),
    anno_height   = as.numeric(.get_flag_value("--anno_height",   .args, required = FALSE, default = "2.5"))
  )
}
