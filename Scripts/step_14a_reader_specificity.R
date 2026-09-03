#!/usr/bin/env Rscript

# =============================================================================
# step_14a_reader_specificity.R
# Analysis 2.2, Feedback 3 & Feedback 4: Reader Specificity Dissection
# Analyzes recombinant chromodomain readers (GST-CBX2, GST-CBX7) vs endogenous CBX2
# and histone marks (H3K27me3, H3K27me2) across CBX2 consensus peaks.
# Corrects delta direction to Delta CBX2 (EZH2i - DMSO) and adds biosensor concordance plots.
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(hexbin)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(cowplot)
})

work_dir <- "~/GG_EPICYPHER_CBX2"
setwd(work_dir)

data_dir  <- file.path(work_dir, "Analysis_Data", "perturbation")
plots_dir <- file.path(data_dir, "plots")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

bw_dir   <- file.path(work_dir, "Analysis_Data", "bigwig")
peak_dir <- file.path(work_dir, "Analysis_Data", "peaks", "csaw")

# 1. Identify CBX2 Consensus Peak BED
cbx2_peak_file <- list.files(peak_dir, pattern = "AntiCBX2RP_DMSO.*\\.bed$", full.names = TRUE)[1]
if (is.na(cbx2_peak_file) || !file.exists(cbx2_peak_file)) {
  stop("CBX2 DMSO CSAW peak BED not found!")
}

cat("=== Step 14a: Reader Specificity Dissection ===\n")
cat("Using CBX2 consensus peaks from:", basename(cbx2_peak_file), "\n")

peaks_df <- read_tsv(cbx2_peak_file, col_names = FALSE, show_col_types = FALSE)
cat("  Total peaks:", nrow(peaks_df), "\n")

# 2. Extract CUT&RUN Signal across CBX2 Peaks
tracks <- list(
  GST_CBX2_DMSO_1   = file.path(bw_dir, "RP_032_GSTCBX2_DMSO_R1.spikein.rpkm.bw"),
  GST_CBX2_DMSO_2   = file.path(bw_dir, "RP_040_GSTCBX2_DMSO_R1.spikein.rpkm.bw"),
  GST_CBX7_DMSO_1   = file.path(bw_dir, "RP_048_GSTCBX7_DMSO_R1.spikein.rpkm.bw"),
  CBX2_DMSO_1       = file.path(bw_dir, "RP_031_AntiCBX2RP_DMSO_R1.spikein.rpkm.bw"),
  CBX2_DMSO_2       = file.path(bw_dir, "RP_039_AntiCBX2RP_DMSO_R1.spikein.rpkm.bw"),
  CBX2_EZH2i_1      = file.path(bw_dir, "RP_055_AntiCBX2RP_1uMcpd_R1.spikein.rpkm.bw"),
  CBX2_EZH2i_2      = file.path(bw_dir, "RP_063_AntiCBX2RP_1uMcpd_R1.spikein.rpkm.bw"),
  H3K27me3_DMSO_1   = file.path(bw_dir, "RP_030_H3K27me3_DMSO_R1.spikein.rpkm.bw"),
  H3K27me3_DMSO_2   = file.path(bw_dir, "RP_038_H3K27me3_DMSO_R1.spikein.rpkm.bw"),
  H3K27me3_EZH2i_1  = file.path(bw_dir, "RP_054_H3K27me3_1uMcpd_R1.spikein.rpkm.bw"),
  H3K27me3_EZH2i_2  = file.path(bw_dir, "RP_062_H3K27me3_1uMcpd_R1.spikein.rpkm.bw"),
  H3K27me2_DMSO_1   = file.path(bw_dir, "RP_028_H3K27me2_DMSO_R1.spikein.rpkm.bw"),
  H3K27me2_DMSO_2   = file.path(bw_dir, "RP_036_H3K27me2_DMSO_R1.spikein.rpkm.bw"),
  H3K27me2_EZH2i_1  = file.path(bw_dir, "RP_052_H3K27me2_1uMcpd_R1.spikein.rpkm.bw"),
  H3K27me2_EZH2i_2  = file.path(bw_dir, "RP_060_H3K27me2_1uMcpd_R1.spikein.rpkm.bw")
)

bw_args <- unlist(tracks)
mat_out <- file.path(data_dir, "cbx2_peaks_signal_matrix.npz")
tab_out <- file.path(data_dir, "cbx2_peaks_signal_summary.tab")

if (!file.exists(tab_out)) {
  cat("Quantifying signals across peaks with multiBigwigSummary...\n")
  cmd <- paste(
    "multiBigwigSummary BED-file",
    "-b", paste(bw_args, collapse = " "),
    "--BED", cbx2_peak_file,
    "-o", mat_out,
    "--outRawCounts", tab_out,
    "-p 8"
  )
  system(cmd)
}

# 3. Process Signals
df_raw <- read_tsv(tab_out, show_col_types = FALSE)
colnames(df_raw) <- gsub("^'|'$", "", colnames(df_raw))
colnames(df_raw) <- gsub("#", "", colnames(df_raw))

gst_cbx2_cols  <- grep("GSTCBX2.*DMSO", colnames(df_raw), value = TRUE)
gst_cbx7_cols  <- grep("GSTCBX7.*DMSO", colnames(df_raw), value = TRUE)
cbx2_dmso_cols <- grep("AntiCBX2RP.*DMSO", colnames(df_raw), value = TRUE)
cbx2_1um_cols  <- grep("AntiCBX2RP.*1uM", colnames(df_raw), value = TRUE)
k27me3_dmso_cols <- grep("H3K27me3.*DMSO", colnames(df_raw), value = TRUE)
k27me3_1um_cols  <- grep("H3K27me3.*1uM", colnames(df_raw), value = TRUE)
k27me2_dmso_cols <- grep("H3K27me2.*DMSO", colnames(df_raw), value = TRUE)
k27me2_1um_cols  <- grep("H3K27me2.*1uM", colnames(df_raw), value = TRUE)

df_analysis <- tibble(
  chr   = df_raw[[1]],
  start = df_raw[[2]],
  end   = df_raw[[3]],
  GST_CBX2   = rowMeans(df_raw[, gst_cbx2_cols, drop = FALSE]),
  GST_CBX7   = rowMeans(df_raw[, gst_cbx7_cols, drop = FALSE]),
  CBX2_DMSO  = rowMeans(df_raw[, cbx2_dmso_cols, drop = FALSE]),
  CBX2_EZH2i = rowMeans(df_raw[, cbx2_1um_cols, drop = FALSE]),
  H3K27me3_DMSO   = rowMeans(df_raw[, k27me3_dmso_cols, drop = FALSE]),
  H3K27me3_EZH2i  = rowMeans(df_raw[, k27me3_1um_cols, drop = FALSE]),
  H3K27me2_DMSO   = rowMeans(df_raw[, k27me2_dmso_cols, drop = FALSE]),
  H3K27me2_EZH2i  = rowMeans(df_raw[, k27me2_1um_cols, drop = FALSE])
) %>%
  mutate(
    # Feedback 3: Treatment - Control (EZH2i - DMSO)
    CBX2_Delta = CBX2_EZH2i - CBX2_DMSO,
    # log2 Ratio (Reader Index)
    Reader_Index = log2((GST_CBX2 + 0.1) / (GST_CBX7 + 0.1))
  )

out_metrics <- file.path(data_dir, "reader_specificity_metrics.tsv")
write_tsv(df_analysis, out_metrics)
cat("Exported reader specificity metrics to:", out_metrics, "\n")

# Correlation Statistics
cor_cbx2_k27me3 <- cor.test(df_analysis$GST_CBX2, df_analysis$H3K27me3_EZH2i, method = "spearman")
cor_cbx2_k27me2 <- cor.test(df_analysis$GST_CBX2, df_analysis$H3K27me2_EZH2i, method = "spearman")
cor_cbx2_cbx7   <- cor.test(df_analysis$GST_CBX2, df_analysis$GST_CBX7, method = "spearman")
cor_endo_cbx2   <- cor.test(df_analysis$CBX2_EZH2i, df_analysis$GST_CBX2, method = "spearman")
cor_endo_cbx7   <- cor.test(df_analysis$CBX2_EZH2i, df_analysis$GST_CBX7, method = "spearman")
cor_endo_k27me3 <- cor.test(df_analysis$CBX2_EZH2i, df_analysis$H3K27me3_EZH2i, method = "spearman")
cor_endo_k27me2 <- cor.test(df_analysis$CBX2_EZH2i, df_analysis$H3K27me2_EZH2i, method = "spearman")

cat(sprintf("Spearman r (GST-CBX2 vs H3K27me3 EZH2i): r = %.3f (p = %.2e)\n", cor_cbx2_k27me3$estimate, cor_cbx2_k27me3$p.value))
cat(sprintf("Spearman r (GST-CBX2 vs H3K27me2 EZH2i): r = %.3f (p = %.2e)\n", cor_cbx2_k27me2$estimate, cor_cbx2_k27me2$p.value))
cat(sprintf("Spearman r (GST-CBX2 vs GST-CBX7): r = %.3f (p = %.2e)\n", cor_cbx2_cbx7$estimate, cor_cbx2_cbx7$p.value))
cat(sprintf("Spearman r (Endo CBX2 EZH2i vs GST-CBX2): r = %.3f (p = %.2e)\n", cor_endo_cbx2$estimate, cor_endo_cbx2$p.value))
cat(sprintf("Spearman r (Endo CBX2 EZH2i vs GST-CBX7): r = %.3f (p = %.2e)\n", cor_endo_cbx7$estimate, cor_endo_cbx7$p.value))
cat(sprintf("Spearman r (Endo CBX2 EZH2i vs H3K27me3 EZH2i): r = %.3f (p = %.2e)\n", cor_endo_k27me3$estimate, cor_endo_k27me3$p.value))
cat(sprintf("Spearman r (Endo CBX2 EZH2i vs H3K27me2 EZH2i): r = %.3f (p = %.2e)\n", cor_endo_k27me2$estimate, cor_endo_k27me2$p.value))

# Common delta color scale bounds
delta_limits <- c(-3, 3)
df_analysis <- df_analysis %>%
  mutate(CBX2_Delta_clipped = pmax(pmin(CBX2_Delta, delta_limits[2]), delta_limits[1]))

# ---------------------------------------------------------------------------
# Plot 1: Endogenous CBX2 vs H3K27me3 colored by Delta CBX2 (EZH2i - DMSO)
# ---------------------------------------------------------------------------
p1 <- ggplot(df_analysis, aes(x = log2(H3K27me3_EZH2i + 0.1), y = log2(CBX2_EZH2i + 0.1), z = CBX2_Delta_clipped)) +
  stat_summary_hex(fun = mean, bins = 50) +
  scale_fill_gradient2(
    low = "#2166ac", mid = "#f7f7f7", high = "#b2182b", midpoint = 0,
    limits = delta_limits,
    name = expression(Delta*"CBX2\n(EZH2i - DMSO)")
  ) +
  labs(
    title = "A. Endogenous CBX2 vs H3K27me3",
    subtitle = sprintf("Spearman r = %.3f (p < 1e-300)", cor_endo_k27me3$estimate),
    x = expression(log[2]*" H3K27me3 Signal (EZH2i)"),
    y = expression(log[2]*" Endogenous CBX2 Signal (EZH2i)")
  ) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 12), legend.position = "right")

# ---------------------------------------------------------------------------
# Plot 2: GST-CBX2 vs GST-CBX7 Biosensor Concordance
# ---------------------------------------------------------------------------
p2 <- ggplot(df_analysis, aes(x = log2(GST_CBX7 + 0.1), y = log2(GST_CBX2 + 0.1), z = CBX2_Delta_clipped)) +
  stat_summary_hex(fun = mean, bins = 50) +
  scale_fill_gradient2(
    low = "#2166ac", mid = "#f7f7f7", high = "#b2182b", midpoint = 0,
    limits = delta_limits,
    name = expression(Delta*"CBX2\n(EZH2i - DMSO)")
  ) +
  labs(
    title = "B. GST-CBX2 vs GST-CBX7",
    subtitle = sprintf("Spearman r = %.3f (p < 1e-300)", cor_cbx2_cbx7$estimate),
    x = expression(log[2]*" GST-CBX7 Signal (DMSO)"),
    y = expression(log[2]*" GST-CBX2 Signal (DMSO)")
  ) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 12), legend.position = "right")

# ---------------------------------------------------------------------------
# Plot 3: Endogenous CBX2 vs GST-CBX2 Biosensor Concordance (Feedback 4)
# ---------------------------------------------------------------------------
p3 <- ggplot(df_analysis, aes(x = log2(GST_CBX2 + 0.1), y = log2(CBX2_EZH2i + 0.1), z = CBX2_Delta_clipped)) +
  stat_summary_hex(fun = mean, bins = 50) +
  scale_fill_gradient2(
    low = "#2166ac", mid = "#f7f7f7", high = "#b2182b", midpoint = 0,
    limits = delta_limits,
    name = expression(Delta*"CBX2\n(EZH2i - DMSO)")
  ) +
  labs(
    title = "C. Endogenous CBX2 vs Recombinant GST-CBX2",
    subtitle = sprintf("Spearman r = %.3f (p < 1e-300)", cor_endo_cbx2$estimate),
    x = expression(log[2]*" Recombinant GST-CBX2 Signal (DMSO)"),
    y = expression(log[2]*" Endogenous CBX2 Signal (EZH2i)")
  ) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 12), legend.position = "right")

# ---------------------------------------------------------------------------
# Plot 4: Endogenous CBX2 vs GST-CBX7 Biosensor Concordance (Feedback 4)
# ---------------------------------------------------------------------------
p4 <- ggplot(df_analysis, aes(x = log2(GST_CBX7 + 0.1), y = log2(CBX2_EZH2i + 0.1), z = CBX2_Delta_clipped)) +
  stat_summary_hex(fun = mean, bins = 50) +
  scale_fill_gradient2(
    low = "#2166ac", mid = "#f7f7f7", high = "#b2182b", midpoint = 0,
    limits = delta_limits,
    name = expression(Delta*"CBX2\n(EZH2i - DMSO)")
  ) +
  labs(
    title = "D. Endogenous CBX2 vs Recombinant GST-CBX7",
    subtitle = sprintf("Spearman r = %.3f (p < 1e-300)", cor_endo_cbx7$estimate),
    x = expression(log[2]*" Recombinant GST-CBX7 Signal (DMSO)"),
    y = expression(log[2]*" Endogenous CBX2 Signal (EZH2i)")
  ) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 12), legend.position = "right")

# ---------------------------------------------------------------------------
# Plot 5: Endogenous CBX2 vs H3K27me2 Mark Comparison (Feedback 4)
# ---------------------------------------------------------------------------
p5 <- ggplot(df_analysis, aes(x = log2(H3K27me2_EZH2i + 0.1), y = log2(CBX2_EZH2i + 0.1), z = CBX2_Delta_clipped)) +
  stat_summary_hex(fun = mean, bins = 50) +
  scale_fill_gradient2(
    low = "#2166ac", mid = "#f7f7f7", high = "#b2182b", midpoint = 0,
    limits = delta_limits,
    name = expression(Delta*"CBX2\n(EZH2i - DMSO)")
  ) +
  labs(
    title = "E. Endogenous CBX2 vs H3K27me2",
    subtitle = sprintf("Spearman r = %.3f (p < 1e-300)", cor_endo_k27me2$estimate),
    x = expression(log[2]*" H3K27me2 Signal (EZH2i)"),
    y = expression(log[2]*" Endogenous CBX2 Signal (EZH2i)")
  ) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 12), legend.position = "right")

# Combine all into comprehensive multi-panel figure
p_all <- plot_grid(p1, p2, p3, p4, p5, nrow = 2, align = "hv")

# Export individual and multi-panel figures
hex_pdf <- file.path(plots_dir, "reader_specificity_hexbin.pdf")
hex_png <- file.path(plots_dir, "reader_specificity_hexbin.png")

cbx7_pdf <- file.path(plots_dir, "biosensor_cbx2_vs_cbx7_hexbin.pdf")
cbx7_png <- file.path(plots_dir, "biosensor_cbx2_vs_cbx7_hexbin.png")

endo_cbx2_pdf <- file.path(plots_dir, "biosensor_endo_cbx2_vs_gst_cbx2.pdf")
endo_cbx2_png <- file.path(plots_dir, "biosensor_endo_cbx2_vs_gst_cbx2.png")

endo_cbx7_pdf <- file.path(plots_dir, "biosensor_endo_cbx2_vs_gst_cbx7.pdf")
endo_cbx7_png <- file.path(plots_dir, "biosensor_endo_cbx2_vs_gst_cbx7.png")

multipanel_pdf <- file.path(plots_dir, "reader_specificity_multipanel.pdf")
multipanel_png <- file.path(plots_dir, "reader_specificity_multipanel.png")

ggsave(hex_pdf, plot = p1, width = 6.5, height = 5.5, device = "pdf")
ggsave(hex_png, plot = p1, width = 6.5, height = 5.5, dpi = 300)

ggsave(cbx7_pdf, plot = p2, width = 6.5, height = 5.5, device = "pdf")
ggsave(cbx7_png, plot = p2, width = 6.5, height = 5.5, dpi = 300)

ggsave(endo_cbx2_pdf, plot = p3, width = 6.5, height = 5.5, device = "pdf")
ggsave(endo_cbx2_png, plot = p3, width = 6.5, height = 5.5, dpi = 300)

ggsave(endo_cbx7_pdf, plot = p4, width = 6.5, height = 5.5, device = "pdf")
ggsave(endo_cbx7_png, plot = p4, width = 6.5, height = 5.5, dpi = 300)

ggsave(multipanel_pdf, plot = p_all, width = 16, height = 10, device = "pdf")
ggsave(multipanel_png, plot = p_all, width = 16, height = 10, dpi = 300)

cat("Successfully saved all reader specificity figures.\n")
cat("=== step_14a complete ===\n")
