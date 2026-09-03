#!/usr/bin/env Rscript

# cbx2_deg_plots.R
# Approach #2: ΔExpression-centric CBX2 repression analysis — R plots

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(showtext)
})

args <- commandArgs(trailingOnly = TRUE)
get_flag <- function(flag, args, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0L) return(default)
  args[idx[1L] + 1L]
}

out_dir  <- get_flag("--out_dir",  args, "Analysis_Data/cbx2_deg")
plots_dir <- file.path(out_dir, "plots")
rds_dir   <- file.path(out_dir, "rds")
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== cbx2_deg_plots.R ===\n")
cat("Plots dir:", plots_dir, "\n")

# --- Fonts & Theme ---
font_add("Helvetica", regular = "/usr/share/fonts/dejavu-sans-fonts/DejaVuSans.ttf")
showtext_auto()
PLOT_FONT <- "Helvetica"

my_theme <- theme_bw(base_size = 14, base_family = PLOT_FONT) +
  theme(
    panel.grid      = element_blank(),
    panel.border    = element_rect(color = "black", fill = NA, linewidth = 1.2),
    plot.title      = element_text(face = "bold", hjust = 0.5, size = 13),
    axis.text       = element_text(color = "black"),
    legend.position = "right"
  )

# --- Colors ---
# Group colors
col_grp <- c(
  "A" = "#DAA520",   # goldenrod — core CBX2 targets
  "B" = "#800080",   # purple    — CBX2-only
  "C" = "#4169E1",   # royalblue — synergistic
  "D" = "#A0A0A0"    # gray      — NS background
)
grp_labels <- c(
  "A" = "A: CBX2-sensitive at DMSO",
  "B" = "B: EZH2i-sensitive",
  "C" = "C: CBX2-sensitive at EZH2i",
  "D" = "D: NS background"
)

# --- Load master ---
master_rds <- file.path(rds_dir, "cbx2_deg_master.rds")
if (!file.exists(master_rds)) stop("cbx2_deg_master.rds not found! Run cbx2_deg_prep.R first.")
df <- readRDS(master_rds)
df$deg_group <- factor(df$deg_group, levels = c("A", "B", "C", "D"))

cat("Loaded cbx2_deg_master.rds. Priority group sizes:\n")
print(table(df$deg_group))
cat("\nNon-exclusive group memberships:\n")
if ("deg_groups" %in% colnames(df)) print(table(df$deg_groups))

# ============================================================
# Plot 1: Boxplot — CBX2 promoter signal (EZH2i) by group
# This answers: do genes that become derepressed upon CBX2 KO+EZH2i
# carry more CBX2 occupancy at their promoters under EZH2i treatment?
# ============================================================
p1 <- ggplot(df, aes(x = deg_group, y = cbx2_prom_1um, fill = deg_group)) +
  geom_violin(alpha = 0.6, color = NA) +
  geom_boxplot(width = 0.25, color = "black", outlier.shape = NA, fill = "white") +
  scale_fill_manual(values = col_grp, labels = grp_labels, name = "Group") +
  scale_x_discrete(labels = grp_labels) +
  labs(
    title = "CBX2 Promoter Signal (EZH2i) by DEG Group",
    x     = "Gene Group",
    y     = "CBX2 logCPM at Promoter (EZH2i / 1uM)"
  ) +
  my_theme +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "none")

svg(file.path(plots_dir, "boxplot_cbx2_prom_by_group.svg"), width = 7, height = 5)
print(p1)
invisible(dev.off())
cat("Saved boxplot_cbx2_prom_by_group.svg\n")

# ============================================================
# Plot 2: Boxplot — CBX2 gene body signal (EZH2i) by group
# This answers: do genes that become derepressed under EZH2i treatment
# have more CBX2 occupancy across their gene bodies?
# ============================================================
p2 <- ggplot(df, aes(x = deg_group, y = cbx2_gene_1um, fill = deg_group)) +
  geom_violin(alpha = 0.6, color = NA) +
  geom_boxplot(width = 0.25, color = "black", outlier.shape = NA, fill = "white") +
  scale_fill_manual(values = col_grp, labels = grp_labels, name = "Group") +
  scale_x_discrete(labels = grp_labels) +
  labs(
    title = "CBX2 Gene Body Signal (EZH2i) by DEG Group",
    x     = "Gene Group",
    y     = "CBX2 logCPM at Gene Body (EZH2i / 1uM)"
  ) +
  my_theme +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "none")

svg(file.path(plots_dir, "boxplot_cbx2_gene_by_group.svg"), width = 7, height = 5)
print(p2)
invisible(dev.off())
cat("Saved boxplot_cbx2_gene_by_group.svg\n")

# ============================================================
# Plot 3: Side-by-side expression boxplot (WT vs KO+DMSO vs KO+EZH2i) by group
# Visualizes derepression magnitude across all three conditions
# ============================================================
if (!"exp_cbx2ko_1um" %in% colnames(df)) {
  warning("exp_cbx2ko_1um column missing — showing only two conditions")
  df$exp_cbx2ko_1um <- NA_real_
}

expr_long <- df %>%
  select(gene_name, deg_group, exp_control, exp_cbx2ko, exp_cbx2ko_1um) %>%
  pivot_longer(cols = c(exp_control, exp_cbx2ko, exp_cbx2ko_1um),
               names_to  = "condition",
               values_to = "expression") %>%
  mutate(
    condition = factor(condition,
                       levels = c("exp_control", "exp_cbx2ko", "exp_cbx2ko_1um"),
                       labels = c("LUC+DMSO\n(WT)", "CBX2 KO\n+DMSO", "CBX2 KO\n+EZH2i"))
  ) %>%
  filter(!is.na(expression))

p3 <- ggplot(expr_long, aes(x = deg_group, y = expression, fill = condition)) +
  geom_boxplot(width = 0.6, color = "black", outlier.shape = NA,
               position = position_dodge(width = 0.75)) +
  scale_fill_manual(values = c("LUC+DMSO\n(WT)"    = "#2196F3",
                                "CBX2 KO\n+DMSO"    = "#FF7043",
                                "CBX2 KO\n+EZH2i"   = "#8B0000"),
                    name = "Condition") +
  scale_x_discrete(labels = grp_labels) +
  labs(
    title = "Gene Expression by Group: WT vs CBX2 KO (+/- EZH2i)",
    x     = "Gene Group",
    y     = "rlogCPM (RNA-seq)"
  ) +
  my_theme +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

svg(file.path(plots_dir, "boxplot_expression_wt_vs_ko.svg"), width = 9, height = 5)
print(p3)
invisible(dev.off())
cat("Saved boxplot_expression_wt_vs_ko.svg\n")

# ============================================================
# Plot 4a: Volcano — CBX2 KO + DMSO vs LUC + DMSO
# ============================================================
plot_volcano <- function(deg_df, title, out_file, grp_A, grp_B, grp_C) {
  deg_df$group <- "D"
  deg_df$group[deg_df$Gene %in% grp_C] <- "C"
  deg_df$group[deg_df$Gene %in% grp_B] <- "B"
  deg_df$group[deg_df$Gene %in% grp_A] <- "A"
  deg_df$group <- factor(deg_df$group, levels = c("A", "B", "C", "D"))
  deg_df$neg_log10_fdr <- -log10(pmax(deg_df$padj, 1e-100))

  # Label only A and B genes
  label_df <- deg_df[deg_df$group %in% c("A", "B"), ]

  p <- ggplot(deg_df, aes(x = log2FoldChange, y = neg_log10_fdr,
                           color = group, alpha = group, size = group)) +
    geom_point() +
    scale_color_manual(values = col_grp, labels = grp_labels, name = "Group") +
    scale_alpha_manual(values = c("A" = 1, "B" = 1, "C" = 0.7, "D" = 0.3), guide = "none") +
    scale_size_manual(values = c("A" = 2.5, "B" = 2.5, "C" = 1.5, "D" = 0.8), guide = "none") +
    geom_vline(xintercept = c(-0.585, 0.585), linetype = "dashed", color = "gray40") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40") +
    labs(
      title = title,
      x     = "log2 Fold Change (CBX2 KO / LUC)",
      y     = "-log10(FDR)"
    ) +
    my_theme

  if (nrow(label_df) > 0) {
    p <- p + ggrepel::geom_text_repel(
      data = label_df,
      aes(label = Gene),
      size = 3, color = "black", family = PLOT_FONT,
      max.overlaps = 30, box.padding = 0.3
    )
  }
  svg(out_file, width = 8, height = 6)
  print(p)
  invisible(dev.off())
  cat("Saved", out_file, "\n")
}

# Extract exact groups from the master dataframe we built in prep
# Non-exclusive: use per-gene flags if available, otherwise fall back to priority group
if ("in_grp_A" %in% colnames(df)) {
  grp_A_genes <- df$gene_name[df$in_grp_A]
  grp_B_genes <- df$gene_name[df$in_grp_B]
  grp_C_genes <- df$gene_name[df$in_grp_C]
} else {
  grp_A_genes <- df$gene_name[df$deg_group == "A"]
  grp_B_genes <- df$gene_name[df$deg_group == "B"]
  grp_C_genes <- df$gene_name[df$deg_group == "C"]
}

# Try ggrepel, fall back to plain labels if missing
if (!requireNamespace("ggrepel", quietly = TRUE)) {
  ggrepel <- list(geom_text_repel = function(...) ggplot2::geom_text(...))
}

setwd(get_flag("--workdir", args, "~/GG_EPICYPHER_CBX2"))

# Also we need to make sure deg_dmso and deg_1um are loaded for volcano values
deg_dmso <- read.csv("Importable_Data/rnaseq/results_CBX2_DMSO_vs_LUC_DMSO.csv", stringsAsFactors = FALSE)
deg_1um_ko <- read.csv("Importable_Data/rnaseq/results_CBX2_1uM_vs_LUC_1uM.csv", stringsAsFactors = FALSE)

plot_volcano(
  deg_dmso,
  "CBX2 KO + DMSO vs LUC + DMSO",
  file.path(plots_dir, "volcano_CBX2ko_DMSO.svg"),
  grp_A_genes, grp_B_genes, grp_C_genes
)
plot_volcano(
  deg_1um_ko,
  "CBX2 KO + EZH2i vs LUC + EZH2i",
  file.path(plots_dir, "volcano_CBX2ko_1uM.svg"),
  grp_A_genes, grp_B_genes, grp_C_genes
)

cat("=== Plots completed successfully ===\n")
