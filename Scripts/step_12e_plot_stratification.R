#!/usr/bin/env Rscript

# =============================================================================
# step_12e_plot_stratification.R
# Analysis 1.2, Feedback 2 & Feedback 4: Direct vs. Indirect Target Stratification
# Multi-panel Box/Violin (Linear & Log2FC) & ECDF plots across RNA-seq response tiers
# for CBX2, H3K27me3, H3K27me2, H2AK119ub, and IgG
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(cowplot)
})

work_dir <- "~/GG_EPICYPHER_CBX2"
setwd(work_dir)

data_dir  <- file.path(work_dir, "Analysis_Data", "perturbation")
sig_dir   <- file.path(data_dir, "signal")
plots_dir <- file.path(data_dir, "plots")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

# 1. Load Gene Universe
univ_file <- file.path(data_dir, "gene_universe.tsv")
gene_df <- read_tsv(univ_file, show_col_types = FALSE)

# 2. Load Promoter Signal Summary Table from deepTools
sig_file <- file.path(sig_dir, "promoter_signal_summary.tab")
sig_raw <- read_tsv(sig_file, show_col_types = FALSE)

# Clean column names (strip quotes and paths)
colnames(sig_raw) <- gsub("^'|'$", "", colnames(sig_raw))
colnames(sig_raw) <- gsub("#", "", colnames(sig_raw))

# Columns in raw table: 'chr', 'start', 'end', then BigWig sample paths
# Map back to BED regions: promoters_2kb_all.bed has coordinates: chr, start, end
prom_bed <- file.path(data_dir, "bed", "promoters_2kb_all.bed")
bed_df <- read_tsv(prom_bed, col_names = c("chr", "start", "end", "gene", "score", "strand"), show_col_types = FALSE)

# Add gene names
sig_raw$gene <- bed_df$gene

# Identify columns for each target
cbx2_dmso_cols   <- grep("AntiCBX2RP.*DMSO", colnames(sig_raw), value = TRUE)
cbx2_1um_cols    <- grep("AntiCBX2RP.*1uM", colnames(sig_raw), value = TRUE)
k27me3_dmso_cols <- grep("H3K27me3.*DMSO", colnames(sig_raw), value = TRUE)
k27me3_1um_cols  <- grep("H3K27me3.*1uM", colnames(sig_raw), value = TRUE)
k27me2_dmso_cols <- grep("H3K27me2.*DMSO", colnames(sig_raw), value = TRUE)
k27me2_1um_cols  <- grep("H3K27me2.*1uM", colnames(sig_raw), value = TRUE)
k119ub_dmso_cols <- grep("H2AK119ub.*DMSO", colnames(sig_raw), value = TRUE)
k119ub_1um_cols  <- grep("H2AK119ub.*1uM", colnames(sig_raw), value = TRUE)
igg_dmso_cols    <- grep("RbIgG.*DMSO", colnames(sig_raw), value = TRUE)
igg_1um_cols     <- grep("RbIgG.*1uM", colnames(sig_raw), value = TRUE)

cat("Detected signal columns:\n")
cat("  CBX2 DMSO    :", paste(basename(cbx2_dmso_cols), collapse = ", "), "\n")
cat("  CBX2 EZH2i   :", paste(basename(cbx2_1um_cols), collapse = ", "), "\n")
cat("  H3K27me3 DMSO:", paste(basename(k27me3_dmso_cols), collapse = ", "), "\n")
cat("  H3K27me3 EZH2i:", paste(basename(k27me3_1um_cols), collapse = ", "), "\n")
cat("  H2AK119ub DMSO:", paste(basename(k119ub_dmso_cols), collapse = ", "), "\n")
cat("  H2AK119ub EZH2i:", paste(basename(k119ub_1um_cols), collapse = ", "), "\n")
cat("  IgG DMSO     :", paste(basename(igg_dmso_cols), collapse = ", "), "\n")
cat("  IgG EZH2i    :", paste(basename(igg_1um_cols), collapse = ", "), "\n")

# Compute mean signal per condition and delta (EZH2i - DMSO)
pseudo <- 1
sig_df <- tibble(
  gene = sig_raw$gene,
  
  CBX2_DMSO   = rowMeans(sig_raw[, cbx2_dmso_cols, drop = FALSE]),
  CBX2_EZH2i  = rowMeans(sig_raw[, cbx2_1um_cols, drop = FALSE]),
  Delta_CBX2  = CBX2_EZH2i - CBX2_DMSO,
  Log2FC_CBX2 = log2((CBX2_EZH2i + pseudo) / (CBX2_DMSO + pseudo)),
  
  H3K27me3_DMSO   = rowMeans(sig_raw[, k27me3_dmso_cols, drop = FALSE]),
  H3K27me3_EZH2i  = rowMeans(sig_raw[, k27me3_1um_cols, drop = FALSE]),
  Delta_H3K27me3  = H3K27me3_EZH2i - H3K27me3_DMSO,
  Log2FC_H3K27me3 = log2((H3K27me3_EZH2i + pseudo) / (H3K27me3_DMSO + pseudo)),
  
  H3K27me2_DMSO   = rowMeans(sig_raw[, k27me2_dmso_cols, drop = FALSE]),
  H3K27me2_EZH2i  = rowMeans(sig_raw[, k27me2_1um_cols, drop = FALSE]),
  Delta_H3K27me2  = H3K27me2_EZH2i - H3K27me2_DMSO,
  Log2FC_H3K27me2 = log2((H3K27me2_EZH2i + pseudo) / (H3K27me2_DMSO + pseudo)),
  
  H2AK119ub_DMSO   = rowMeans(sig_raw[, k119ub_dmso_cols, drop = FALSE]),
  H2AK119ub_EZH2i  = rowMeans(sig_raw[, k119ub_1um_cols, drop = FALSE]),
  Delta_H2AK119ub  = H2AK119ub_EZH2i - H2AK119ub_DMSO,
  Log2FC_H2AK119ub = log2((H2AK119ub_EZH2i + pseudo) / (H2AK119ub_DMSO + pseudo)),
  
  IgG_DMSO   = rowMeans(sig_raw[, igg_dmso_cols, drop = FALSE]),
  IgG_EZH2i  = rowMeans(sig_raw[, igg_1um_cols, drop = FALSE]),
  Delta_IgG  = IgG_EZH2i - IgG_DMSO,
  Log2FC_IgG = log2((IgG_EZH2i + pseudo) / (IgG_DMSO + pseudo))
)

# Load csaw CBX2 LFC
csaw_file <- file.path(data_dir, "csaw_cbx2_promoter_lfc.tsv")
csaw_df <- read_tsv(csaw_file, show_col_types = FALSE)

# 3. Merge with Expression Tiers
merged_df <- gene_df %>%
  dplyr::select(gene, tier_CBX2KO, log2FC_CBX2KO, padj_CBX2KO, tier_EZH2i, log2FC_EZH2i, padj_EZH2i) %>%
  inner_join(sig_df, by = "gene") %>%
  left_join(csaw_df %>% dplyr::select(gene, log2FC_CBX2_binding, logCPM_CBX2_EZH2i), by = "gene")

# Ensure tier order: Downregulated, Unchanged, Upregulated
tier_levels <- c("Downregulated", "Unchanged", "Upregulated")
merged_df$tier_EZH2i <- factor(merged_df$tier_EZH2i, levels = tier_levels)

# Use csaw log2FC for CBX2 instead of deepTools bigwig Log2FC
merged_df$Log2FC_CBX2 <- merged_df$log2FC_CBX2_binding

# Pivot long for Linear Delta
long_delta <- merged_df %>%
  pivot_longer(
    cols = starts_with("Delta_"),
    names_to = "Target",
    names_prefix = "Delta_",
    values_to = "Delta_Signal"
  )

# Pivot long for Log2FC
long_lfc <- merged_df %>%
  pivot_longer(
    cols = starts_with("Log2FC_"),
    names_to = "Target",
    names_prefix = "Log2FC_",
    values_to = "Log2FC_Signal"
  )

target_order <- c("CBX2", "H3K27me3", "H3K27me2", "H2AK119ub", "IgG")
long_delta$Target <- factor(long_delta$Target, levels = target_order)
long_lfc$Target   <- factor(long_lfc$Target, levels = target_order)

# 4. Statistical Testing (Wilcoxon rank-sum test: Upregulated vs Unchanged)
stat_delta <- long_delta %>%
  group_by(Target) %>%
  summarise(
    wilcox_p_up = tryCatch(wilcox.test(Delta_Signal[tier_EZH2i == "Upregulated"],
                                       Delta_Signal[tier_EZH2i == "Unchanged"])$p.value, error = function(e) NA),
    median_up = median(Delta_Signal[tier_EZH2i == "Upregulated"], na.rm = TRUE),
    median_unchanged = median(Delta_Signal[tier_EZH2i == "Unchanged"], na.rm = TRUE),
    .groups = "drop"
  )

cat("Statistical comparison (Linear Delta):\n")
print(stat_delta)

# 5. Plotting Functions and Multi-panel Generation
tier_colors <- c("Downregulated" = "#2166ac", "Unchanged" = "#999999", "Upregulated" = "#b2182b")

# Helper function to generate individual plots with DEG-derived zooming
create_strat_plot <- function(df, target, y_col, y_label, title_label) {
  df_sub <- df %>% filter(Target == target)
  
  # Calculate limits based ONLY on DEGs
  deg_df <- df_sub %>% filter(tier_EZH2i %in% c("Upregulated", "Downregulated"))
  ymin <- min(deg_df[[y_col]], na.rm = TRUE)
  ymax <- max(deg_df[[y_col]], na.rm = TRUE)
  margin <- (ymax - ymin) * 0.2
  plot_ymin <- ymin - margin
  plot_ymax <- ymax + margin
  
  # To avoid the huge Unchanged density squishing the violin, we filter Unchanged outliers from the plot dataframe
  # Note: they are still in the data, but removed from visual calculation so it doesn't squish. 
  # Actually coord_cartesian(ylim) zooms in without dropping them from stat calculations, which is exactly what we want!
  
  p <- ggplot(df_sub, aes(x = tier_EZH2i, y = .data[[y_col]], fill = tier_EZH2i)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.5) +
    geom_violin(alpha = 0.5, trim = TRUE, scale = "width", draw_quantiles = c(0.5)) +
    geom_boxplot(width = 0.25, outlier.shape = NA, alpha = 0.75, color = "black", linewidth = 0.4) +
    geom_jitter(data = deg_df,
                aes(color = tier_EZH2i), width = 0.15, size = 1.8, alpha = 0.9) +
    scale_fill_manual(values = tier_colors, guide = "none") +
    scale_color_manual(values = tier_colors, guide = "none") +
    coord_cartesian(ylim = c(plot_ymin, plot_ymax)) +
    labs(title = title_label, x = NULL, y = y_label) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
      axis.text.x = element_text(angle = 30, hjust = 1, face = "bold"),
      panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3)
    )
  return(p)
}

# 5.1 Plot 1: Linear Delta Stratification
plots_delta <- lapply(target_order, function(tgt) {
  create_strat_plot(long_delta, tgt, "Delta_Signal", expression(Delta*" Signal (EZH2i - DMSO)"), tgt)
})
p_violin <- plot_grid(plotlist = plots_delta, nrow = 1, align = "h")
title_violin <- ggdraw() + draw_label("Target Stratification: Promoter CUT&RUN Signal Delta", fontface = 'bold', size = 16)
sub_violin <- ggdraw() + draw_label("EZH2i - DMSO across EZH2i expression response tiers", size = 12, color = "grey30")
p_violin_full <- plot_grid(title_violin, sub_violin, p_violin, ncol = 1, rel_heights = c(0.08, 0.05, 1))

# 5.2 Plot 2: Log2 Fold Change Stratification
plots_lfc <- lapply(target_order, function(tgt) {
  create_strat_plot(long_lfc, tgt, "Log2FC_Signal", expression(log[2]*"FC (EZH2i / DMSO)"), tgt)
})
p_lfc <- plot_grid(plotlist = plots_lfc, nrow = 1, align = "h")
title_lfc <- ggdraw() + draw_label("Target Stratification: Promoter CUT&RUN Log2 Fold Change", fontface = 'bold', size = 16)
sub_lfc <- ggdraw() + draw_label("log2(EZH2i / DMSO) across EZH2i expression response tiers", size = 12, color = "grey30")
p_lfc_full <- plot_grid(title_lfc, sub_lfc, p_lfc, ncol = 1, rel_heights = c(0.08, 0.05, 1))

# 5.3 Plot 3: LogCPM for EZH2i (CBX2 only as requested)
df_cpm <- merged_df %>% filter(!is.na(logCPM_CBX2_EZH2i))
p_cpm <- ggplot(df_cpm, aes(x = tier_EZH2i, y = logCPM_CBX2_EZH2i, fill = tier_EZH2i)) +
  geom_violin(alpha = 0.5, trim = TRUE, scale = "width", draw_quantiles = c(0.5)) +
  geom_boxplot(width = 0.25, outlier.shape = NA, alpha = 0.75, color = "black", linewidth = 0.4) +
  geom_jitter(data = filter(df_cpm, tier_EZH2i %in% c("Upregulated", "Downregulated")),
              aes(color = tier_EZH2i), width = 0.15, size = 1.8, alpha = 0.9) +
  scale_fill_manual(values = tier_colors, guide = "none") +
  scale_color_manual(values = tier_colors, guide = "none") +
  labs(title = "CBX2 Promoter Binding (EZH2i)", x = "EZH2i Transcriptional Response Tier", y = "logCPM (EZH2i samples)") +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    axis.text.x = element_text(angle = 30, hjust = 1, face = "bold"),
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3)
  )

# Export figures
violin_pdf <- file.path(plots_dir, "target_stratification_violin.pdf")
violin_png <- file.path(plots_dir, "target_stratification_violin.png")

lfc_pdf <- file.path(plots_dir, "target_stratification_log2fc.pdf")
lfc_png <- file.path(plots_dir, "target_stratification_log2fc.png")

cpm_pdf <- file.path(plots_dir, "target_stratification_cbx2_logcpm.pdf")
cpm_png <- file.path(plots_dir, "target_stratification_cbx2_logcpm.png")

ggsave(violin_pdf, plot = p_violin_full, width = 16, height = 6.5, device = "pdf", bg = "white")
ggsave(violin_png, plot = p_violin_full, width = 16, height = 6.5, dpi = 300, bg = "white")

ggsave(lfc_pdf, plot = p_lfc_full, width = 16, height = 6.5, device = "pdf", bg = "white")
ggsave(lfc_png, plot = p_lfc_full, width = 16, height = 6.5, dpi = 300, bg = "white")

ggsave(cpm_pdf, plot = p_cpm, width = 6, height = 5.5, device = "pdf")
ggsave(cpm_png, plot = p_cpm, width = 6, height = 5.5, dpi = 300)

cat("Successfully saved:\n")
cat("  Linear Violin/Box :", violin_pdf, "\n")
cat("  Log2FC Violin/Box :", lfc_pdf, "\n")
cat("  CBX2 LogCPM Violin:", cpm_pdf, "\n")
cat("=== step_12e complete ===\n")

