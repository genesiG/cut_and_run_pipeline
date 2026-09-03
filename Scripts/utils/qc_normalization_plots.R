#!/usr/bin/env Rscript

# =============================================================================
# qc_normalization_plots.R  –  MA plots, PCA, and Spearman correlation
#                               heatmaps comparing normalization scale factors.
#
# Inputs (all resolved from config.py):
#   - Processed BAMs in config.PROCESSEDBAMDIR
#   - Scale factor TSV files in config.SCALINGDIR/{spikein,chipseqspikeinfree}
#   - config.NORMALIZATION_RUNS (defines targets and parameters)
#
# Outputs (config.SCALINGDIR/qc/):
#   {target}_MA.pdf          — standard log2 CPM MA plots (all methods)
#   {target}_MA_scaled.pdf   — scaled MA plots (spikein + csif only)
#   {target}_PCA.svg         — PCA per normalization method
#   {target}_correlation.svg — Spearman correlation heatmaps
# =============================================================================

#----------------------#
# IMPORT CONFIG MODULE #
#----------------------#
library(reticulate, quietly = TRUE, verbose = FALSE)

use_python(Sys.which("python"), required = TRUE)
py_run_file("Scripts/config.py")

#------------------------#
# LOAD REQUIRED PACKAGES #
#------------------------#
suppressMessages({
  library(dplyr, quietly = TRUE)
  library(readr, quietly = TRUE)
  library(tidyr, quietly = TRUE)
  library(KernSmooth, quietly = TRUE)
  library(edgeR, quietly = TRUE)
  library(csaw, quietly = TRUE)
  library(stringr, quietly = TRUE)
})

# Import parameters from configuration file
SCALING       <- py$SCALINGDIR
WORKDIR       <- file.path(SCALING, "qc")
BAMDIR        <- py$PROCESSEDBAMDIR
METADATA      <- py$METADATA
IS_PAIRED_END <- py$IS_PAIRED_END
MAPQ          <- py$MAPQ
NORM_RUNS     <- py$NORMALIZATION_RUNS  # named list from config.NORMALIZATION_RUNS

if (py$REMOVE_DUPLICATES == TRUE) {
  SUFFIX <- paste0(".qc.sort.rmdup.mapq", MAPQ, ".final.bam")
  DEDUP  <- FALSE
} else {
  SUFFIX <- paste0(".qc.sort.markdup.mapq", MAPQ, ".final.bam")
  DEDUP  <- TRUE
}

#----------

## PREPARE DATA

# Set up working directory
if (!dir.exists(WORKDIR)) {
  dir.create(WORKDIR, recursive = TRUE)
}
setwd(WORKDIR)

if (IS_PAIRED_END) {
  pe <- "both"
} else {
  pe <- "none"
}

# Build a named character vector: run_name -> target regex pattern
TARGETS <- vapply(
  names(NORM_RUNS),
  function(run_name) as.character(NORM_RUNS[[run_name]][["targets"]]),
  character(1)
)

cat(paste0("Generating MA plots for the following runs: ",
           paste(names(TARGETS), collapse = ", "), "\n"))


#-------------------------#
# EXTRACT READS INTO BINS #
#-------------------------#
for (target in TARGETS) {
  cat(paste0("\nCalculating scale factors for target: ", target, "\n"))
  matching_files <- list.files(path = BAMDIR,
                               pattern = target,
                               full.names = TRUE,
                               ignore.case = TRUE)

  matching_files <- matching_files[!grepl("\\.bai$", matching_files)]
  # Exclude spike-in BAMs — E. coli / mm39 reads must not be used for
  # background bin counting or normalization comparisons.
  matching_files <- matching_files[!grepl("spike_in\\.bam$", matching_files)]


  # --- SAFETY CHECKS ---
  if (length(matching_files) == 0) {
    message("ERROR: No BAM files found for target: ", target)
    next
  }

  if (anyNA(matching_files)) {
    message("ERROR:  NA entries detected in matching_files. Check filename filtering.")
    next
  }

  missing <- matching_files[!file.exists(matching_files)]
  if (length(missing) > 0) {
    message("Warning:  These BAM files do not exist:\n",
            paste(missing, collapse="\n"))
    message("Continuing without these files")
  }

  bam.files    <- matching_files
  sample.names <- basename(bam.files)

  cat("\nUsing the following BAM files to compare normalization factors:\n")
  print(sample.names)

  ## Define readParam object with parameters for extracting reads from BAM files
  param <- readParam(minq = MAPQ,
                     pe = "both",
                     max.frag = 500,
                     dedup = DEDUP)

  # Count reads into large background bins
  counts <- windowCounts(bam.files,
                         bin = TRUE,
                         width = 20000,
                         param = param)

  # 1. Dynamically find Spike-in and TMM scale factor files
  spikein_sf <- list.files(path = file.path(SCALING, "spikein"),
                           pattern = paste0(target, ".*_spikein_SF\\.txt$"),
                           full.names = TRUE,
                           ignore.case = TRUE)

  tmm_sf <- list.files(path = SCALING,
                       pattern = paste0(target, ".*_tmm_SF\\.txt$"),
                       full.names = TRUE,
                       ignore.case = TRUE,
                       recursive = TRUE)

  # 2. Determine EXACT ChIPseqSpikeInFree SF file path based on config params
  run_name <- names(TARGETS)[TARGETS == target][1]
  params   <- NORM_RUNS[[run_name]]

  bin_size  <- if (!is.null(params$bin_size))  params$bin_size  else 10000
  cutoff    <- if (!is.null(params$cutoff))    params$cutoff    else 1.2
  max_turns <- if (!is.null(params$max_turns)) params$max_turns else 0.99

  filename_prefix <- run_name
  if (cutoff != 1.2) {
    filename_prefix <- paste0(filename_prefix, "_cutoff_", cutoff)
  }
  if (max_turns != 0.99) {
    filename_prefix <- paste0(filename_prefix, "_max_turns_", max_turns)
  }

  run_dir <- if (bin_size != 10000) paste0(run_name, "_", bin_size, "bp_bins") else run_name
  csif_sf <- file.path(SCALING, "chipseqspikeinfree", run_dir,
                       paste0(filename_prefix, "_SF.txt"))

  if (!file.exists(csif_sf)) {
    csif_sf <- character(0)
  }

  sf_files <- c(spikein_sf, tmm_sf, csif_sf)

  if (length(sf_files) == 0) {
    stop(paste0("No scale factor files found for target: ", target))
  }

  # Read scale factor files and annotate with normalization method
  sf_data <- lapply(sf_files, function(file) {
    df <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

    if (grepl("_spikein_SF", file)) {
      df$norm_method <- "spikein"
      df$ref_color   <- "black"
    } else if (grepl("_tmm_SF", file)) {
      df$norm_method <- "tmm"
      df$ref_color   <- "purple"
    } else {
      df$norm_method <- "chipseqspikeinfree"
      df$ref_color   <- "red"
    }

    df
  })

  # Merge all SF data into one data frame
  df_merged <- bind_rows(sf_data)

  #-------------------#
  # GENERATE MA PLOTS #
  #-------------------#

  ##  Calculate raw count matrix (bins x samples) for PCA / heatmap
  raw_counts_mat <- assay(counts, "counts")

  ##  Calculate read abundances (log2 CPM) for the standard MA plot
  adj.counts <- edgeR::cpm(counts, log = TRUE)

  ## Helper: strip BAM suffix to get a clean sample label
  clean_label <- function(x) str_remove(x, fixed(SUFFIX))

  # Get unique normalization methods present in the SF files
  method_colors <- df_merged %>%
    select(norm_method, ref_color) %>%
    distinct() %>%
    arrange(norm_method)

  # Number of pairwise comparisons (all vs. sample 1)
  n_comparisons <- length(sample.names) - 1

  # ------------------------------------------------------------------
  # Helper: extract aligned SF vector for a given method
  # ------------------------------------------------------------------
  clean_id <- function(x) sub("\\.qc\\.sort\\..*", "", basename(x))
  clean_sample_names <- clean_id(sample.names)

  get_sf <- function(method) {
    df.sf <- df_merged %>%
      filter(norm_method == method) %>%
      select(ID, SF, ref_color) %>%
      mutate(clean_ID = clean_id(as.character(ID)))

    list(
      sf       = df.sf$SF[match(clean_sample_names, df.sf$clean_ID)],
      line_col = df.sf$ref_color[match(clean_sample_names[1], df.sf$clean_ID)]
    )
  }

  # ==================================================================
  # 1. STANDARD MA PLOT  (log2 CPM, one plot per method per comparison)
  # ==================================================================
  cat("\nExporting standard MA plots...\n")

  # Replace the literal ".*" with an underscore
  target_clean <- gsub("\\.\\*", "_", target)
  # Collapse multiple underscores into a single one (e.g. "K27M__" -> "K27M_")
  target_clean <- gsub("_+", "_", target_clean)
  # Remove trailing or leading underscores (e.g. "K27M_" -> "K27M")
  target_clean <- gsub("^_|_$", "", target_clean)
  pdf(file.path(WORKDIR, paste0(target_clean, "_MA.pdf")),
      width = 12.5, height = 4)

  for (method in method_colors$norm_method) {
    sf_info   <- get_sf(method)
    method_sf <- sf_info$sf
    line_col  <- sf_info$line_col
    if (is.na(line_col) || length(line_col) != 1) line_col <- "black"

    if (length(method_sf) == 0 || all(is.na(method_sf))) {
      cat(paste0("  No scale factors for method: ", method, "\n"))
      next
    }

    par(mfrow = c(1, 3), mar = c(7, 9, 5, 5), mgp = c(3.5, 0.65, 0))

    for (i in seq_len(n_comparisons)) {
      cur.x <- adj.counts[, 1]
      cur.y <- adj.counts[, 1 + i]

      sample1 <- clean_label(sample.names[1])
      sample2 <- clean_label(sample.names[i + 1])

      smoothScatter(
        x = (cur.x + cur.y) / 2,
        y = cur.x - cur.y,
        cex.axis = 1.25, cex.lab = 1.5,
        xlab = "Average signal (logCPM)\nacross samples",
        ylab = "Difference in signal\nbetween samples (log ratio)",
        main = paste0(sample1, "\nvs\n", sample2, "\n(", method, ")")
      )

      all.dist <- diff(log2(method_sf[c(i + 1, 1)]))
      if (!is.na(all.dist)) {
        abline(h = all.dist, col = line_col, lwd = 2)
        legend("topright",
               legend = paste("SF difference:", round(all.dist, 3)),
               col = line_col, lty = 1, lwd = 2, cex = 0.8)
      } else {
        legend("topright", legend = "SF missing", text.col = "red", bty = "n", cex = 0.8)
      }
    }
  }

  dev.off()
  cat("  -> Saved:", paste0(target_clean, "_MA.pdf\n"))


  # ==================================================================
  # 3. PCA  (one page per method)
  # ==================================================================
  cat("\nExporting PCA plots...\n")

  svg(file.path(WORKDIR, paste0(target_clean, "_PCA.svg")),
      width = 7, height = 6)

  for (method in method_colors$norm_method) {
    sf_info   <- get_sf(method)
    method_sf <- sf_info$sf

    if (length(method_sf) == 0 || all(is.na(method_sf))) next

    cat(paste0("  PCA for method: ", method, "\n"))

    valid_idx <- which(!is.na(method_sf))
    if (length(valid_idx) < 3) {
      cat("    Not enough valid samples for PCA. Skipping.\n")
      plot(1, type="n", axes=FALSE, xlab="", ylab="")
      text(1, 1, "Not enough valid samples\nfor PCA", col="red", cex=1.5)
      next
    }

    scaled_mat_pca <- sweep(raw_counts_mat[, valid_idx, drop=FALSE],
                            2, method_sf[valid_idx], `*`)
    scaled_log_pca <- log2(scaled_mat_pca + 1)

    pca_res <- prcomp(t(scaled_log_pca), center = TRUE, scale. = FALSE)
    pca_df  <- as.data.frame(pca_res$x[, 1:min(2, ncol(pca_res$x))])
    pca_df$sample <- clean_label(sample.names[valid_idx])
    var_exp <- round(summary(pca_res)$importance[2, 1:2] * 100, 1)

    par(mfrow = c(1, 1), mar = c(6, 6, 5, 5), mgp = c(3.5, 0.65, 0))
    plot(
      pca_df$PC1, pca_df$PC2,
      pch = 19, cex = 1.4,
      col = valid_idx,
      xlab = paste0("PC1 (", var_exp[1], "% variance)"),
      ylab = paste0("PC2 (", var_exp[2], "% variance)"),
      main = paste0(target, " \u2014 PCA (", method, " scaled)")
    )
    text(pca_df$PC1, pca_df$PC2,
         labels = pca_df$sample, pos = 3, cex = 0.75)
  }

  dev.off()
  cat("  -> Saved:", paste0(target_clean, "_PCA.svg\n"))

  # ==================================================================
  # 4. CORRELATION HEATMAP  (one page per method)
  # ==================================================================
  cat("\nExporting correlation heatmaps...\n")
  svg(file.path(WORKDIR, paste0(target_clean, "_correlation.svg")),
      width = 8, height = 7)

  for (method in method_colors$norm_method) {
    sf_info   <- get_sf(method)
    method_sf <- sf_info$sf

    if (length(method_sf) == 0 || all(is.na(method_sf))) next

    cat(paste0("  Correlation heatmap for method: ", method, "\n"))

    valid_idx <- which(!is.na(method_sf))
    if (length(valid_idx) < 2) {
      cat("    Not enough valid samples for correlation. Skipping.\n")
      plot(1, type="n", axes=FALSE, xlab="", ylab="")
      text(1, 1, "Not enough valid samples\nfor correlation", col="red", cex=1.5)
      next
    }

    scaled_mat_cor <- sweep(raw_counts_mat[, valid_idx, drop=FALSE],
                            2, method_sf[valid_idx], `*`)
    scaled_log_cor <- log2(scaled_mat_cor + 1)

    cor_mat <- cor(scaled_log_cor, method = "spearman")
    colnames(cor_mat) <- rownames(cor_mat) <- clean_label(sample.names[valid_idx])

    n_samp  <- ncol(cor_mat)
    lbl_len <- max(nchar(colnames(cor_mat)))
    par(mfrow = c(1, 1),
        mar   = c(max(8, lbl_len * 0.6), max(8, lbl_len * 0.6), 5, 3),
        mgp   = c(3.5, 0.65, 0),
        xpd   = TRUE) # allows legend outside of plotting area

    col_ramp <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
    image(
      1:n_samp, 1:n_samp,
      t(cor_mat[n_samp:1, , drop=FALSE]),
      col  = col_ramp, zlim = c(-1, 1),
      xaxt = "n", yaxt = "n", xlab = "", ylab = "",
      main = paste0(target, " \u2014 Spearman correlation (", method, " scaled)")
    )
    axis(1, at = 1:n_samp, labels = colnames(cor_mat),
         las = 2, cex.axis = 0.75)
    axis(2, at = 1:n_samp, labels = rev(rownames(cor_mat)),
         las = 2, cex.axis = 0.75)

    for (ri in 1:n_samp) {
      for (ci in 1:n_samp) {
        text(ci, n_samp + 1 - ri,
             labels = sprintf("%.2f", cor_mat[ri, ci]),
             cex = max(0.5, 0.9 - n_samp * 0.04))
      }
    }
    legend("bottomright", inset = -0.10,
           legend = c("-1.0", "0.0", "+1.0"),
           fill   = col_ramp[c(1, 50, 100)],
           title  = "r", cex = 0.7, bty = "n")
  }

  dev.off()
  cat("  -> Saved:", paste0(target_clean, "_correlation.svg\n"))

  cat(paste0("\nCompleted all plots for target: ", target, "\n"))
  cat("-----------------------------------\n")
}

cat("\nDone.\n\n")
