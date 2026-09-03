#!/usr/bin/env Rscript

# =============================================================================
# step_7_inhibitor_resistance.R
#
# Quantifies how well H3K27me3 signal is maintained at nucleation sites vs.
# spreading sites upon 500 nM EZH2 inhibitor treatment.
#
# Design
# ------
# BigWig files (chipseqspikeinfree.rpkm.bw) and BAM files are available.
# We compute average signal via csaw+edgeR for the K27MKO cell line
# to output DA metrics (logFC and aveLogCPM), and build heatmaps using
# deepTools (computeMatrix and plotHeatmap).
#
# Outputs (Analysis_Data/inhibitor_resistance/)
# -----------------------------------------------
#   per_region_summary.csv              — full data table
#   01_violin_log2fc_nucleation_vs_spreading.pdf
#   02_dotplot_enrichment_logCPM.pdf
#   03_maplot_log2fc_vs_logCPM.pdf
#   heatmap_4conditions_K27M.pdf        — K27M nucleation + spreading, 4 cols
#   heatmap_4conditions_K27MKO.pdf      — K27MKO nucleation + spreading, 4 cols
# =============================================================================

suppressPackageStartupMessages({
  library(reticulate)
  library(GenomicRanges)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(csaw)
  library(edgeR)
})

# Parse CLI arguments using base R to avoid missing 'argparse' package
args_cli <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  idx <- which(args_cli == flag)
  if (length(idx) == 0) return(default)
  if (idx >= length(args_cli)) stop(paste(flag, "requires a value"))
  args_cli[idx + 1]
}

caller_val <- get_arg("--caller")
if (is.null(caller_val)) {
  stop("Usage: Rscript Scripts/utils/inhibitor_resistance_analysis.R --caller <name> [--top-pct <num>]")
}
args <- list(caller = caller_val)

top_pct_str <- get_arg("--top-pct")
top_pct <- if (!is.null(top_pct_str)) as.numeric(top_pct_str) else NULL
if (!is.null(top_pct) && (top_pct <= 0 || top_pct >= 100)) {
  stop("--top-pct must be a number between 0 and 100 (exclusive)")
}
if (!is.null(top_pct)) {
  cat(sprintf("Top-%%: %.1f%% of regions by DMSO logCPM will be used for subset plots and BEDs.\n", top_pct))
}

get_arg_2vals <- function(flag) {
  idx <- which(args_cli == flag)
  if (length(idx) == 0) return(NULL)
  if (idx + 1 >= length(args_cli)) stop(paste(flag, "requires two numeric values: <min> <max>"))
  val1 <- as.numeric(args_cli[idx + 1])
  val2 <- as.numeric(args_cli[idx + 2])
  if (is.na(val1) || is.na(val2)) stop(paste(flag, "requires two numeric values: <min> <max>"))
  c(val1, val2)
}
violin_y_limits <- get_arg_2vals("--violin_y_limits")
if (!is.null(violin_y_limits)) {
  cat(sprintf("Violin plot Y-axis limits overridden by CLI: [%.2f, %.2f]\n", violin_y_limits[1], violin_y_limits[2]))
}
legend_stat_raw <- get_arg("--legend_stat", default = get_arg("--legend-stat", default = "cohens_d"))
legend_stat <- tolower(legend_stat_raw)
if (!legend_stat %in% c("cohens_d", "summary", "stats", "median_mean_n", "mean_median_n", "median", "median_only", "med", "prop_depleted", "percentage", "prop", "pct", "odds_ratio", "or")) {
  stop("--legend_stat must be 'cohens_d', 'summary', 'median', 'prop_depleted', or 'odds_ratio'")
}
stat_display_name <- if (legend_stat %in% c("cohens_d", "d")) "Cohen's D" else if (legend_stat %in% c("median", "median_only", "med")) "Median" else if (legend_stat %in% c("prop_depleted", "percentage", "prop", "pct")) "% Depleted (LFC<0)" else if (legend_stat %in% c("odds_ratio", "or")) "Odds Ratio (LFC<0)" else "Median, Mean, and N"
cat(sprintf("Violin plot legend will display: %s\n", stat_display_name))

# --- Load config via reticulate -----------------------------------------------
use_python(Sys.which("python"), required = TRUE)
py_run_file("Scripts/config.py")

caller_dir <- file.path(py$PEAKANNODIR, args$caller)
bw_dir     <- file.path(py$IMPORTABLE_DATA, "bigwig")
bam_dir    <- file.path(py$IMPORTABLE_DATA, "bam")
out_dir    <- file.path(py$ANALYSIS_DATA, "inhibitor_resistance", args$caller)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

PSEUDOCOUNT <- 0.01   # CPM-scale pseudocount for log2FC

# --- aesthetics ---
my_theme <- theme_bw(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
    strip.background = element_blank(),
    strip.text = element_text(size = 17),
    plot.title = element_text(face = "bold"),
    legend.position = "none",
    axis.text = element_text(color = "black", size = 14),
    axis.title = element_text(color = "black", size = 15),
    axis.line = element_blank()
  )

# =============================================================================
# 1.  File manifest (BAMs for K27MKO DA analysis)
# =============================================================================
bam_files_k27mko <- list(
  K27MKO_DMSO_A  = file.path(bam_dir, "K27MKO_DMSO_K27me3_A_R1.qc.sort.markdup.mapq30.final.bam"),
  K27MKO_DMSO_B  = file.path(bam_dir, "K27MKO_DMSO_K27me3_B_R1.qc.sort.markdup.mapq30.final.bam"),
  K27MKO_500nM_A = file.path(bam_dir, "K27MKO_500nM_K27me3_A_R1.qc.sort.markdup.mapq30.final.bam"),
  K27MKO_500nM_B = file.path(bam_dir, "K27MKO_500nM_K27me3_B_R1.qc.sort.markdup.mapq30.final.bam")
)

missing_bam <- names(bam_files_k27mko)[!sapply(bam_files_k27mko, file.exists)]
if (length(missing_bam) > 0) stop("Missing BAM files: ", paste(missing_bam, collapse = ", "))

cat("All BAM files found.\n")

# =============================================================================
# 2.  Load nucleation and spreading site BEDs as GRanges
# =============================================================================
read_bed <- function(path) {
  cat("Reading:", path, "\n")
  df <- read.table(path, header = FALSE, sep = "\t",
                   col.names = c("chr", "start", "end"))
  makeGRangesFromDataFrame(df,
                           seqnames.field          = "chr",
                           start.field             = "start",
                           end.field               = "end",
                           starts.in.df.are.0based = TRUE)
}

nucleation_gr <- read_bed(file.path(caller_dir, "nucleation_sites.bed"))
spreading_gr  <- read_bed(file.path(caller_dir, "spreading_sites.bed"))

nucleation_gr$site_type <- "Nucleation"
spreading_gr$site_type  <- "Spreading"

cat(sprintf("Nucleation sites: %d\nSpreading sites:  %d\n",
            length(nucleation_gr), length(spreading_gr)))

all_sites <- c(nucleation_gr, spreading_gr)

# =============================================================================
# 3. csaw regionCounts and edgeR for K27MKO
# =============================================================================
cat("\n=== Counting reads in regions using csaw (K27MKO) ===\n")
pe_param <- readParam(pe = if (py$IS_PAIRED_END) "both" else "none")

# Load spikein-free normalization factors for K27MKO
sf_dir <- file.path(py$IMPORTABLE_DATA, "normalization", "chipseqspikeinfree", "K27MKO_K27me3")
sf_files <- list.files(sf_dir, pattern = "_SF\\.txt$", full.names = TRUE)
if (length(sf_files) == 0) stop("No *_SF.txt files in: ", sf_dir)
primary_sf <- sf_files[!grepl("original_SF\\.txt$", sf_files)]
if (length(primary_sf) == 0) primary_sf <- sf_files
sf_table <- read.table(primary_sf[1], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat(sprintf("Loaded SF table: %s (%d rows)\n", basename(primary_sf[1]), nrow(sf_table)))

get_sf <- function(bam_basename) {
  row <- sf_table[sf_table$ID == bam_basename, ]
  if (nrow(row) == 0) stop("Sample not found in SF table: ", bam_basename)
  sf_val <- suppressWarnings(as.numeric(row$SF[1]))
  if (is.na(sf_val)) {
    warning(sprintf("SF for %s is NA — falling back to 1", bam_basename))
    sf_val <- 1
  }
  sf_val
}

bams_char <- unname(unlist(bam_files_k27mko))
sfs <- sapply(basename(bams_char), get_sf)
cat(sprintf("  SFs: %s\n", paste(round(sfs, 3), collapse = ", ")))

counts <- regionCounts(bams_char, all_sites, param = pe_param)
counts$norm.factors <- sfs

cat("Running edgeR differential analysis...\n")
y <- asDGEList(counts)
group <- factor(c("DMSO", "DMSO", "500nM", "500nM"), levels = c("DMSO", "500nM"))
y$samples$group <- group

design <- model.matrix(~ group)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)
res <- glmQLFTest(fit, coef = 2) # tests 500nM vs DMSO

tt <- res$table
all_sites$logFC <- tt$logFC
all_sites$aveLogCPM <- aveLogCPM(y)
all_sites$logCPM_DMSO <- rowMeans(cpm(y[, 1:2], log = TRUE))
all_sites$logCPM_500nM <- rowMeans(cpm(y[, 3:4], log = TRUE))
all_sites$PValue <- tt$PValue
all_sites$FDR <- p.adjust(tt$PValue, method = "BH")

df_all <- as.data.frame(all_sites)
df_all$site_type <- factor(df_all$site_type, levels = c("Nucleation", "Spreading"))

csv_path <- file.path(out_dir, "per_region_summary.csv")
write.csv(df_all, csv_path, row.names = FALSE)
cat("Saved:", csv_path, "\n")

# =============================================================================
# 4. Plots
# =============================================================================
cat("\n=== Generating Plots ===\n")
nuc_col <- "#FF0066"
spr_col <- "#8c6bb1"

# 4a. Violin plot with effect size
calc_effect_size <- function(data) {
  data %>%
    group_by(site_type) %>%
    summarise(
      median_logFC = median(logFC, na.rm = TRUE),
      mean_logFC   = mean(logFC, na.rm = TRUE),
      median_logCPM = median(aveLogCPM, na.rm = TRUE),
      mean_logCPM   = mean(aveLogCPM, na.rm = TRUE),
      cohens_d = (mean(logCPM_500nM, na.rm = TRUE) - mean(logCPM_DMSO, na.rm = TRUE)) /
                 sqrt((var(logCPM_500nM, na.rm = TRUE) + var(logCPM_DMSO, na.rm = TRUE)) / 2),
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(
      label = sprintf("d = %.2f\nn = %d", cohens_d, n)
    )
}

eff_sz <- calc_effect_size(df_all)
ylims <- quantile(df_all$logFC, probs = c(0.001, 0.999), na.rm = TRUE)

make_plots <- function(df, suffix = "", title_suffix = "") {
  cat(sprintf("\n=== Generating Plots (suffix: '%s') ===\n", suffix))
  
  eff_sz_sub <- calc_effect_size(df)
  if (legend_stat %in% c("cohens_d", "d")) {
    legend_title  <- "Cohen's D"
    ez_map_fc     <- setNames(sprintf("%s = %.2f", eff_sz_sub$site_type, eff_sz_sub$cohens_d), eff_sz_sub$site_type)
    ez_map_enrich <- ez_map_fc
  } else if (legend_stat %in% c("median", "median_only", "med")) {
    legend_title  <- "Median"
    ez_map_fc     <- setNames(sprintf("%s = %.2f", eff_sz_sub$site_type, eff_sz_sub$median_logFC), eff_sz_sub$site_type)
    ez_map_enrich <- setNames(sprintf("%s = %.2f", eff_sz_sub$site_type, eff_sz_sub$median_logCPM), eff_sz_sub$site_type)
  } else if (legend_stat %in% c("prop_depleted", "percentage", "prop", "pct")) {
    legend_title  <- "% Depleted (LFC < 0)"
    pct_map <- df %>% group_by(site_type) %>% summarise(pct = mean(logFC < 0, na.rm = TRUE) * 100, .groups = "drop")
    ez_map_fc <- setNames(sprintf("%s = %.1f%%", pct_map$site_type, pct_map$pct), pct_map$site_type)
    ez_map_enrich <- ez_map_fc
  } else if (legend_stat %in% c("odds_ratio", "or")) {
    legend_title  <- "Odds Ratio (LFC < 0)"
    nuc_dep <- sum(df$site_type == "Nucleation" & df$logFC < 0, na.rm = TRUE)
    nuc_not <- sum(df$site_type == "Nucleation" & df$logFC >= 0, na.rm = TRUE)
    spr_dep <- sum(df$site_type == "Spreading" & df$logFC < 0, na.rm = TRUE)
    spr_not <- sum(df$site_type == "Spreading" & df$logFC >= 0, na.rm = TRUE)
    ft <- fisher.test(matrix(c(spr_dep, spr_not, nuc_dep, nuc_not), nrow = 2, byrow = TRUE))
    or_val <- ft$estimate
    p_val  <- ft$p.value
    p_str  <- ifelse(p_val < 1e-15, "< 1e-15", sprintf("= %.1e", p_val))
    ez_map_fc <- c("Nucleation" = "Nucleation (ref, OR=1.0)",
                   "Spreading"  = sprintf("Spreading (OR=%.2f, p%s)", or_val, p_str))
    ez_map_enrich <- ez_map_fc
  } else {
    legend_title  <- "Median, Mean, N"
    ez_map_fc     <- setNames(sprintf("%s: Med=%.2f, Mean=%.2f, N=%d", eff_sz_sub$site_type, eff_sz_sub$median_logFC, eff_sz_sub$mean_logFC, eff_sz_sub$n), eff_sz_sub$site_type)
    ez_map_enrich <- setNames(sprintf("%s: Med=%.2f, Mean=%.2f, N=%d", eff_sz_sub$site_type, eff_sz_sub$median_logCPM, eff_sz_sub$mean_logCPM, eff_sz_sub$n), eff_sz_sub$site_type)
  }
  y_lim_use <- if (!is.null(violin_y_limits)) violin_y_limits else ylims
  
  # 4a. Violin
  p_violin <- ggplot(df, aes(x = site_type, y = logFC, fill = site_type)) +
    geom_violin(trim = TRUE, alpha = 0.7, color = NA) +
    geom_boxplot(width = 0.12, outlier.shape = NA, color = "grey20", fill = "white") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
    scale_fill_manual(name = legend_title,
                      values = c("Nucleation" = nuc_col, "Spreading" = spr_col),
                      labels = ez_map_fc,
                      guide = guide_legend(override.aes = list(fill = NA, color = NA),
                                           keywidth = unit(0, "pt"), keyheight = unit(0, "pt"))) +
    scale_y_continuous(limits = y_lim_use) +
    my_theme +
    theme(legend.position = "right",
          legend.key      = element_blank(),
          plot.title      = element_text(face = "bold")) +
    labs(title = sprintf("K27MKO: H3K27me3 resistance to EZH2 inhibition%s", title_suffix),
         subtitle = "Effect size directly from csaw/edgeR regionCounts",
         y = expression(atop(textstyle("Change in H3K27me3"), textstyle("("*log[2]*" Fold-Change)"))),
         x = "")

  svg(file.path(out_dir, sprintf("01_violin_log2fc_nucleation_vs_spreading%s.svg", suffix)),
      width = 6.5, height = 5)
  print(p_violin)
  dev.off()
  cat(sprintf("Saved: 01_violin_log2fc_nucleation_vs_spreading%s.svg\n", suffix))

  # 4b. Dotplot: Average H3K27me3 enrichment (aveLogCPM) for Nucleation vs Spreading
  p_dot <- ggplot(df, aes(x = site_type, y = aveLogCPM, color = site_type)) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
    scale_color_manual(values = c("Nucleation" = nuc_col, "Spreading" = spr_col)) +
    my_theme +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold")) +
    labs(title = sprintf("Average H3K27me3 abundance%s", title_suffix),
         y = "Average H3K27me3 abundance (Log2 CPM)",
         x = "")

  svg(file.path(out_dir, sprintf("02_dotplot_enrichment_logCPM%s.svg", suffix)), width = 5, height = 5)
  print(p_dot)
  dev.off()
  cat(sprintf("Saved: 02_dotplot_enrichment_logCPM%s.svg\n", suffix))

  # 4b-violin. Violin: Average H3K27me3 enrichment (aveLogCPM) for Nucleation vs Spreading
  p_dot_violin <- ggplot(df, aes(x = site_type, y = aveLogCPM, fill = site_type)) +
    geom_violin(trim = TRUE, alpha = 0.7, color = NA) +
    geom_boxplot(width = 0.12, outlier.shape = NA, color = "grey20", fill = "white") +
    scale_fill_manual(name = legend_title,
                      values = c("Nucleation" = nuc_col, "Spreading" = spr_col),
                      labels = ez_map_enrich,
                      guide = guide_legend(override.aes = list(fill = NA, color = NA),
                                           keywidth = unit(0, "pt"), keyheight = unit(0, "pt"))) +
    my_theme +
    theme(legend.position = "right",
          legend.key      = element_blank(),
          plot.title      = element_text(face = "bold")) +
    labs(title = sprintf("Average H3K27me3 abundance%s", title_suffix),
         y = "Average H3K27me3 abundance (Log2 CPM)",
         x = "")

  svg(file.path(out_dir, sprintf("02b_violin_enrichment_logCPM%s.svg", suffix)), width = 6.5, height = 5)
  print(p_dot_violin)
  dev.off()
  cat(sprintf("Saved: 02b_violin_enrichment_logCPM%s.svg\n", suffix))

  # 4c. MA Plot
  lfc_thresh <- 0
  df$ma_col <- ifelse(df$logFC >= lfc_thresh, "black",
                      ifelse(df$site_type == "Nucleation", nuc_col, spr_col))

  p_ma <- ggplot(df, aes(x = aveLogCPM, y = logFC, color = ma_col)) +
    geom_point(alpha = 0.4, size = 1) +
    geom_hline(yintercept = lfc_thresh, linetype = "dashed", color = "black") +
    geom_smooth(method = "gam", color = "black", se = TRUE, linetype="dotted") +
    facet_wrap(~ site_type) +
    scale_color_identity() +
    my_theme +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold")) +
    labs(title = sprintf("H3K27me3 (K27MKO)%s", title_suffix),
         y = expression(atop(textstyle("Change in H3K27me3"), textstyle("("*log[2]*" Fold-Change)"))),
         x = "Average H3K27me3 abundance (Log2 CPM)")

  svg(file.path(out_dir, sprintf("03_maplot_log2fc_vs_logCPM%s.svg", suffix)), width = 7, height = 4)
  print(p_ma)
  dev.off()
  cat(sprintf("Saved: 03_maplot_log2fc_vs_logCPM%s.svg\n", suffix))
}

# Generate plots for all regions
make_plots(df_all, suffix = "", title_suffix = "")

# Generate plots and BED files for top-pct subset (if requested)
if (!is.null(top_pct)) {
  cat(sprintf("\n=== Filtering for top %.1f%% of regions by DMSO logCPM ===\n", top_pct))
  
  # Group by site_type so we take top X% of Nucleation and top X% of Spreading separately
  df_top <- do.call(rbind, lapply(split(df_all, df_all$site_type), function(sub_df) {
    pct_thresh <- quantile(sub_df$logCPM_DMSO, probs = 1 - top_pct / 100, na.rm = TRUE)
    sub_df[sub_df$logCPM_DMSO >= pct_thresh, ]
  }))
  rownames(df_top) <- NULL
  
  cat(sprintf("Retained %d / %d Nucleation sites and %d / %d Spreading sites.\n",
              sum(df_top$site_type == "Nucleation"), sum(df_all$site_type == "Nucleation"),
              sum(df_top$site_type == "Spreading"), sum(df_all$site_type == "Spreading")))
  
  make_plots(df_top, suffix = sprintf("_top%.0fpct", top_pct),
             title_suffix = sprintf(" (Top %.0f%%)", top_pct))
  
  # Export BED files for deepTools top-pct heatmaps
  nuc_top <- df_top[df_top$site_type == "Nucleation", ]
  spr_top <- df_top[df_top$site_type == "Spreading", ]
  
  nuc_bed_out <- file.path(out_dir, sprintf("nucleation_sites_top%.0fpct.bed", top_pct))
  spr_bed_out <- file.path(out_dir, sprintf("spreading_sites_top%.0fpct.bed", top_pct))
  
  write.table(data.frame(chr = nuc_top$seqnames, start = nuc_top$start, end = nuc_top$end),
              nuc_bed_out, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(data.frame(chr = spr_top$seqnames, start = spr_top$start, end = spr_top$end),
              spr_bed_out, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  cat(sprintf("Saved top %.0f%% BED files for deepTools:\n  %s\n  %s\n",
              top_pct, nuc_bed_out, spr_bed_out))
}

cat(sprintf("\n=== step_7c differential analysis complete. All outputs in: %s ===\n", out_dir))
