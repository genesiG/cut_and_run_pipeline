#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dplyr))

setwd("~/GG_EPICYPHER_CBX2")

# Load DEG tables
a  <- read.csv("Importable_Data/rnaseq/results_CBX2_DMSO_vs_LUC_DMSO.csv")
b  <- read.csv("Importable_Data/rnaseq/results_LUC_1uM_vs_LUC_DMSO.csv")
c2 <- read.csv("Importable_Data/rnaseq/results_CBX2_1uM_vs_LUC_1uM.csv")

LFC <- 0.585; FDR <- 0.05
is_up <- function(df) !is.na(df$padj) & df$log2FoldChange >= LFC & df$padj <= FDR

uA <- a$Gene[is_up(a)]
uB <- b$Gene[is_up(b)]
uC <- c2$Gene[is_up(c2)]

cat("=== New non-exclusive group sizes ===\n")
cat("A (CBX2 KO+DMSO vs LUC+DMSO):", length(uA), "\n")
cat("B (LUC+1uM vs LUC+DMSO):", length(uB), "\n")
cat("C (CBX2 KO+1uM vs LUC+1uM):", length(uC), "\n")
cat("\n=== Overlaps ===\n")
cat("A intersect C:", length(intersect(uA, uC)), "\n")
cat("B intersect C:", length(intersect(uB, uC)), "\n")
cat("A intersect B:", length(intersect(uA, uB)), "\n")
cat("A intersect B intersect C:", length(Reduce(intersect, list(uA, uB, uC))), "\n")

cat("\nCDKN2A in A:", "CDKN2A" %in% uA, "\n")
cat("CDKN2A in B:", "CDKN2A" %in% uB, "\n")
cat("CDKN2A in C:", "CDKN2A" %in% uC, "\n")

cat("\n=== cbx2_expression_master.rds ===\n")
m <- readRDS("Analysis_Data/cbx2_expression/csaw/cbx2_expression_master.rds")
cat("columns:", paste(colnames(m), collapse=", "), "\n")
cat("nrow:", nrow(m), "\n")

cat("\nCDKN2A row:\n")
print(m[m$gene_name == "CDKN2A", ])

cat("\n=== CBX2 prom signal summary by group (current master) ===\n")
master <- readRDS("Analysis_Data/cbx2_deg/rds/cbx2_deg_master.rds")
cat("Master columns:", paste(colnames(master), collapse=", "), "\n")
cat("\ncbx2_prom_1um by group (median):\n")
print(tapply(master$cbx2_prom_1um, master$deg_group, median, na.rm=TRUE))

cat("\n=== Correlation: CBX2 prom signal vs gene expression (WT) ===\n")
cor_val <- cor(master$cbx2_prom_1um, master$exp_control, use="complete.obs", method="spearman")
cat("Spearman cor(cbx2_prom_1um, exp_control):", round(cor_val, 4), "\n")
cor_val2 <- cor(master$cbx2_prom_1um, master$exp_cbx2ko, use="complete.obs", method="spearman")
cat("Spearman cor(cbx2_prom_1um, exp_cbx2ko):", round(cor_val2, 4), "\n")
