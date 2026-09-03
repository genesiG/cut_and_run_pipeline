#!/usr/bin/env Rscript

# cbx2_expression_plots.R
# Phase 2 of CBX2 Expression Analysis

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(showtext)
})

args <- commandArgs(trailingOnly = TRUE)
get_flag_value <- function(flag, args, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0L) return(default)
  args[idx[1L] + 1L]
}

out_dir <- get_flag_value("--out_dir", args, "Analysis_Data/cbx2_expression")
csaw_dir <- file.path(out_dir, "csaw")
plots_dir <- file.path(out_dir, "plots")

cat("=== cbx2_expression_plots.R ===\n")
cat("Plots dir :", plots_dir, "\n")

# Fonts
font_add("Helvetica", regular = "/usr/share/fonts/dejavu-sans-fonts/DejaVuSans.ttf")
showtext_auto()
PLOT_FONT <- "Helvetica"

# Colors
col_cbx2 <- c("High" = "#58135e", "Medium" = "#810f7c", "Low" = "#88419d")
col_mark <- c("High" = "#1f78b4", "Medium" = "#33a02c", "Low" = "#b2df8a")
col_expr <- c("High" = "goldenrod", "Medium" = "gray", "Low" = "purple")

my_theme <- theme_bw(base_size = 14, base_family = PLOT_FONT) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
    plot.title = element_text(face = "bold", hjust = 0.5, size=12),
    axis.text = element_text(color = "black")
  )

master_rds <- file.path(csaw_dir, "cbx2_expression_master.rds")
if (!file.exists(master_rds)) stop("master_rds not found!")
df <- readRDS(master_rds)

conditions <- list(
  control = "terc_control",
  ezh2i   = "terc_ezh2i",
  cbx2ko  = "terc_cbx2ko"
)

treatments <- list(
  dmso = c(cbx2_prom="cbx2_prom_dmso", cbx2_gene="cbx2_gene_dmso", 
           me2_prom="me2_prom_dmso", me2_gene="me2_gene_dmso",
           me3_prom="me3_prom_dmso", me3_gene="me3_gene_dmso"),
  `1uM` = c(cbx2_prom="cbx2_prom_1um", cbx2_gene="cbx2_gene_1um",
            me2_prom="me2_prom_1um", me2_gene="me2_gene_1um",
            me3_prom="me3_prom_1um", me3_gene="me3_gene_1um")
)

plot_stacked_bar <- function(data, x_col, fill_col, fill_colors, title, xlab, ylab, legend_title, out_file) {
  bar_df <- data %>%
    group_by(!!sym(x_col), !!sym(fill_col)) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(!!sym(x_col)) %>%
    mutate(pct = count / sum(count) * 100)
    
  p_bar <- ggplot(bar_df, aes(x = !!sym(x_col), y = pct, fill = !!sym(fill_col))) +
    geom_bar(stat = "identity", color = "black", width = 0.7) +
    scale_fill_manual(values = fill_colors) +
    labs(title = title, x = xlab, y = ylab, fill = legend_title) +
    my_theme
    
  svg(out_file, width = 6, height = 5)
  print(p_bar)
  invisible(dev.off())
}

plot_boxplot <- function(data, x_col, y_col, fill_col, fill_colors, title, xlab, ylab, out_file) {
  p_box <- ggplot(data, aes(x = !!sym(x_col), y = !!sym(y_col), fill = !!sym(fill_col))) +
    geom_violin(alpha=0.6, color=NA) +
    geom_boxplot(width=0.2, color="black", outlier.shape=NA, fill="white") +
    scale_fill_manual(values = fill_colors) +
    labs(title = title, x = xlab, y = ylab) +
    my_theme + theme(legend.position = "none")
    
  svg(out_file, width = 5, height = 5)
  print(p_box)
  invisible(dev.off())
}

for (cond_name in names(conditions)) {
  c_dir <- file.path(plots_dir, cond_name)
  dir.create(c_dir, showWarnings = FALSE, recursive = TRUE)
  
  expr_col <- conditions[[cond_name]]
  sub_df <- df[!is.na(df[[expr_col]]), ]
  sub_df[[expr_col]] <- factor(sub_df[[expr_col]], levels = c("Low", "Medium", "High"))
  
  for (trt_name in names(treatments)) {
    for (reg_type in c("prom", "gene")) {
      cbx2_val_col <- treatments[[trt_name]][paste0("cbx2_", reg_type)]
      cbx2_cat_col <- paste0("terc_", cbx2_val_col)
      
      me2_val_col <- treatments[[trt_name]][paste0("me2_", reg_type)]
      me2_cat_col <- paste0("terc_", me2_val_col)
      
      me3_val_col <- treatments[[trt_name]][paste0("me3_", reg_type)]
      me3_cat_col <- paste0("terc_", me3_val_col)
      
      # 1. Original Plot: CBX2 vs Gene Expression
      plot_stacked_bar(sub_df, expr_col, cbx2_cat_col, col_cbx2,
                       sprintf("CBX2 %s (%s) Abundance by %s Expression", toupper(reg_type), trt_name, toupper(cond_name)),
                       "Gene Expression Category", "Percentage of Genes (%)", "CBX2 Abundance",
                       file.path(c_dir, sprintf("stackedbar_cbx2_vs_expr_%s_%s.svg", reg_type, trt_name)))
                       
      plot_boxplot(sub_df, expr_col, cbx2_val_col, expr_col, col_expr,
                   sprintf("CBX2 %s (%s) logFC over %s Expression", toupper(reg_type), trt_name, toupper(cond_name)),
                   "Gene Expression Category", "CBX2 Abundance (logFC over IgG)",
                   file.path(c_dir, sprintf("boxplot_cbx2_vs_expr_%s_%s.svg", reg_type, trt_name)))
                   
      # 2. CBX2 vs H3K27me2
      plot_stacked_bar(sub_df, me2_cat_col, cbx2_cat_col, col_cbx2,
                       sprintf("CBX2 %s (%s) Abundance by H3K27me2", toupper(reg_type), trt_name),
                       "H3K27me2 Abundance Category", "Percentage of Genes (%)", "CBX2 Abundance",
                       file.path(c_dir, sprintf("stackedbar_cbx2_vs_me2_%s_%s.svg", reg_type, trt_name)))
                       
      plot_boxplot(sub_df, me2_cat_col, cbx2_val_col, me2_cat_col, col_mark,
                   sprintf("CBX2 %s (%s) logFC over H3K27me2", toupper(reg_type), trt_name),
                   "H3K27me2 Abundance Category", "CBX2 Abundance (logFC over IgG)",
                   file.path(c_dir, sprintf("boxplot_cbx2_vs_me2_%s_%s.svg", reg_type, trt_name)))
                   
      # 3. CBX2 vs H3K27me3
      plot_stacked_bar(sub_df, me3_cat_col, cbx2_cat_col, col_cbx2,
                       sprintf("CBX2 %s (%s) Abundance by H3K27me3", toupper(reg_type), trt_name),
                       "H3K27me3 Abundance Category", "Percentage of Genes (%)", "CBX2 Abundance",
                       file.path(c_dir, sprintf("stackedbar_cbx2_vs_me3_%s_%s.svg", reg_type, trt_name)))
                       
      plot_boxplot(sub_df, me3_cat_col, cbx2_val_col, me3_cat_col, col_mark,
                   sprintf("CBX2 %s (%s) logFC over H3K27me3", toupper(reg_type), trt_name),
                   "H3K27me3 Abundance Category", "CBX2 Abundance (logFC over IgG)",
                   file.path(c_dir, sprintf("boxplot_cbx2_vs_me3_%s_%s.svg", reg_type, trt_name)))
                   
      # 4. H3K27me2/3 vs Gene Expression
      plot_stacked_bar(sub_df, expr_col, me2_cat_col, col_mark,
                       sprintf("H3K27me2 %s (%s) by %s Expression", toupper(reg_type), trt_name, toupper(cond_name)),
                       "Gene Expression Category", "Percentage of Genes (%)", "H3K27me2 Abundance",
                       file.path(c_dir, sprintf("stackedbar_me2_vs_expr_%s_%s.svg", reg_type, trt_name)))
                       
      plot_stacked_bar(sub_df, expr_col, me3_cat_col, col_mark,
                       sprintf("H3K27me3 %s (%s) by %s Expression", toupper(reg_type), trt_name, toupper(cond_name)),
                       "Gene Expression Category", "Percentage of Genes (%)", "H3K27me3 Abundance",
                       file.path(c_dir, sprintf("stackedbar_me3_vs_expr_%s_%s.svg", reg_type, trt_name)))
    }
  }
}
cat("=== Plots completed successfully ===\n")
