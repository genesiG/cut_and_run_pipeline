#!/usr/bin/env Rscript

# =============================================================================
# step_12c_plot_dual_perturbation.R
# Analysis 1.1 & Feedback 5: Dual-Perturbation Functional Coupling & Combined Perturbation
# - Panel A: CBX2 KO vs EZH2i (colored by Delta CBX2 promoter binding)
# - Panel B: Combined CBX2 KO + EZH2i vs CBX2 KO single perturbation
# - Panel C: Direct correlation between Delta CBX2 promoter binding and Combined expression change
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(readr)
  library(cowplot)
})

work_dir <- "~/GG_EPICYPHER_CBX2"
setwd(work_dir)

data_dir  <- file.path(work_dir, "Analysis_Data", "perturbation")
plots_dir <- file.path(data_dir, "plots")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

# 1. Load Gene Universe and Differential Binding
univ_file <- file.path(data_dir, "gene_universe.tsv")
bind_file <- file.path(data_dir, "csaw_cbx2_promoter_lfc.tsv")

if (!file.exists(univ_file) || !file.exists(bind_file)) {
  stop("Input files not found! Ensure step_12a and step_12b have completed.")
}

gene_df <- read_tsv(univ_file, show_col_types = FALSE)
bind_df <- read_tsv(bind_file, show_col_types = FALSE)

# Merge datasets
df <- gene_df %>%
  inner_join(bind_df %>% dplyr::select(gene, log2FC_CBX2_binding, FDR_CBX2_binding, logCPM_CBX2), by = "gene")

cat("Loaded dataset with", nrow(df), "genes.\n")

# Color scale bounds for visual clarity
lfc_bind_limits <- c(-3, 3)
df_plot <- df %>%
  mutate(
    color_val = pmax(pmin(log2FC_CBX2_binding, lfc_bind_limits[2]), lfc_bind_limits[1])
  )

# Candidate genes to label
top_targets <- df_plot %>%
  filter(
    (log2FC_CBX2KO > 0.8 & log2FC_EZH2i > 0.8) |
    (tier_CBX2KO == "Upregulated" & tier_EZH2i == "Upregulated") |
    (log2FC_CBX2_binding < -1.2 & (log2FC_CBX2KO > 0.5 | log2FC_Dual > 1.0))
  ) %>%
  arrange(log2FC_CBX2_binding) %>%
  head(15)

# ---------------------------------------------------------------------------
# Panel A: CBX2 KO vs EZH2i
# ---------------------------------------------------------------------------
cor_a <- cor.test(df_plot$log2FC_CBX2KO, df_plot$log2FC_EZH2i, method = "spearman")

p_a <- ggplot(df_plot, aes(x = log2FC_CBX2KO, y = log2FC_EZH2i)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_point(aes(color = color_val), alpha = 0.7, size = 1.4) +
  scale_color_gradient2(
    low = "#2166ac", mid = "#f7f7f7", high = "#b2182b", midpoint = 0,
    limits = lfc_bind_limits,
    name = expression(Delta*" CBX2 Binding\n("*log[2]*"FC under EZH2i)")
  ) +
  geom_text_repel(
    data = top_targets, aes(label = gene),
    size = 3.2, fontface = "bold", box.padding = 0.3, max.overlaps = 20, segment.color = "grey40"
  ) +
  labs(
    title = "A. Dual-Perturbation Functional Coupling",
    subtitle = sprintf("CBX2 KO vs. EZH2i (Spearman r = %.3f, p = %.2e)", cor_a$estimate, cor_a$p.value),
    x = expression(log[2]*"FC(Expression) [ CBX2 KO vs sgLUC ]"),
    y = expression(log[2]*"FC(Expression) [ EZH2i vs DMSO ]")
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10, color = "grey30"),
    legend.position = "right",
    panel.grid.major = element_line(color = "grey94", linewidth = 0.3)
  )

# ---------------------------------------------------------------------------
# Panel B: Combined Perturbation (CBX2 KO + EZH2i) vs CBX2 KO Single
# ---------------------------------------------------------------------------
cor_b <- cor.test(df_plot$log2FC_CBX2KO, df_plot$log2FC_Dual, method = "spearman")

p_b <- ggplot(df_plot, aes(x = log2FC_CBX2KO, y = log2FC_Dual)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey30", linewidth = 0.6) +
  geom_point(aes(color = color_val), alpha = 0.7, size = 1.4) +
  scale_color_gradient2(
    low = "#2166ac", mid = "#f7f7f7", high = "#b2182b", midpoint = 0,
    limits = lfc_bind_limits,
    name = expression(Delta*" CBX2 Binding\n("*log[2]*"FC under EZH2i)")
  ) +
  geom_text_repel(
    data = top_targets, aes(label = gene),
    size = 3.2, fontface = "bold", box.padding = 0.3, max.overlaps = 20, segment.color = "grey40"
  ) +
  labs(
    title = "B. Combined (CBX2 KO + EZH2i) vs. Single CBX2 KO",
    subtitle = sprintf("Expression response (Spearman r = %.3f, p = %.2e)", cor_b$estimate, cor_b$p.value),
    x = expression(log[2]*"FC(Expression) [ CBX2 KO vs sgLUC ]"),
    y = expression(log[2]*"FC(Expression) [ CBX2 KO + EZH2i vs sgLUC ]")
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10, color = "grey30"),
    legend.position = "right",
    panel.grid.major = element_line(color = "grey94", linewidth = 0.3)
  )

# ---------------------------------------------------------------------------
# Panel C: Delta CBX2 Binding vs Combined Expression Change (DEGs only)
# ---------------------------------------------------------------------------
df_c <- df_plot %>% filter(tier_EZH2i %in% c("Upregulated", "Downregulated"))

cor_c_up <- cor.test(df_c$log2FC_CBX2_binding[df_c$tier_EZH2i == "Upregulated"], df_c$log2FC_Dual[df_c$tier_EZH2i == "Upregulated"], method = "spearman")
cor_c_dn <- cor.test(df_c$log2FC_CBX2_binding[df_c$tier_EZH2i == "Downregulated"], df_c$log2FC_Dual[df_c$tier_EZH2i == "Downregulated"], method = "spearman")

p_c <- ggplot(df_c, aes(x = log2FC_CBX2_binding, y = log2FC_Dual)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_point(aes(color = color_val), alpha = 0.7, size = 1.4) +
  geom_smooth(method = "lm", color = "black", linewidth = 0.8, linetype = "solid", se = TRUE) +
  scale_color_gradient2(
    low = "#2166ac", mid = "#f7f7f7", high = "#b2182b", midpoint = 0,
    limits = lfc_bind_limits,
    name = expression(Delta*" CBX2 Binding\n("*log[2]*"FC under EZH2i)")
  ) +
  facet_wrap(~tier_EZH2i, scales = "free") +
  labs(
    title = "C. Delta CBX2 Binding vs. Combined Expression Change",
    subtitle = sprintf("DEGs only. Spearman r (Up) = %.2f, r (Dn) = %.2f", cor_c_up$estimate, cor_c_dn$estimate),
    x = expression(Delta*" CBX2 Promoter Binding [ "*log[2]*"FC (EZH2i vs DMSO) ]"),
    y = expression(log[2]*"FC(Expression) [ CBX2 KO + EZH2i vs sgLUC ]")
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10, color = "grey30"),
    legend.position = "right",
    panel.grid.major = element_line(color = "grey94", linewidth = 0.3)
  )

# ---------------------------------------------------------------------------
# Panel D: Delta CBX2 Binding vs EZH2i Expression Change (DEGs only)
# ---------------------------------------------------------------------------
df_d <- df_plot %>% filter(tier_EZH2i %in% c("Upregulated", "Downregulated"))

cor_d_up <- cor.test(df_d$log2FC_CBX2_binding[df_d$tier_EZH2i == "Upregulated"], df_d$log2FC_EZH2i[df_d$tier_EZH2i == "Upregulated"], method = "spearman")
cor_d_dn <- cor.test(df_d$log2FC_CBX2_binding[df_d$tier_EZH2i == "Downregulated"], df_d$log2FC_EZH2i[df_d$tier_EZH2i == "Downregulated"], method = "spearman")

p_d <- ggplot(df_d, aes(x = log2FC_CBX2_binding, y = log2FC_EZH2i)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_point(aes(color = color_val), alpha = 0.7, size = 1.4) +
  geom_smooth(method = "lm", color = "black", linewidth = 0.8, linetype = "solid", se = TRUE) +
  scale_color_gradient2(
    low = "#2166ac", mid = "#f7f7f7", high = "#b2182b", midpoint = 0,
    limits = lfc_bind_limits,
    name = expression(Delta*" CBX2 Binding\n("*log[2]*"FC under EZH2i)")
  ) +
  facet_wrap(~tier_EZH2i, scales = "free") +
  labs(
    title = "D. Delta CBX2 Binding vs. EZH2i Expression Change",
    subtitle = sprintf("DEGs only. Spearman r (Up) = %.2f, r (Dn) = %.2f", cor_d_up$estimate, cor_d_dn$estimate),
    x = expression(Delta*" CBX2 Promoter Binding [ "*log[2]*"FC (EZH2i vs DMSO) ]"),
    y = expression(log[2]*"FC(Expression) [ EZH2i vs DMSO ]")
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10, color = "grey30"),
    legend.position = "right",
    panel.grid.major = element_line(color = "grey94", linewidth = 0.3)
  )

# Combine into multi-panel figure
p_combined <- plot_grid(p_a, p_b, p_c, p_d, nrow = 2, align = "hv", axis = "tb")

# Export figures
scatter_pdf <- file.path(plots_dir, "dual_perturbation_scatter.pdf")
scatter_png <- file.path(plots_dir, "dual_perturbation_scatter.png")

combined_pdf <- file.path(plots_dir, "dual_perturbation_multipanel.pdf")
combined_png <- file.path(plots_dir, "dual_perturbation_multipanel.png")

ggsave(scatter_pdf, plot = p_a, width = 7.5, height = 6.5, device = "pdf")
ggsave(scatter_png, plot = p_a, width = 7.5, height = 6.5, dpi = 300)

ggsave(combined_pdf, plot = p_combined, width = 18, height = 11, device = "pdf")
ggsave(combined_png, plot = p_combined, width = 18, height = 11, dpi = 300)

cat("Successfully saved:\n")
cat("  Single 4-quadrant plot :", scatter_pdf, "\n")
cat("  Multi-panel comparison :", combined_pdf, "\n")
cat("=== step_12c complete ===\n")
