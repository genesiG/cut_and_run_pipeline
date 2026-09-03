#!/usr/bin/env Rscript

# cbx2_deg_prep.R
# Approach #2: ΔExpression-centric CBX2 repression analysis
# Classify genes by derepression upon CBX2 KO (+/- EZH2i) and write BED files

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
get_flag <- function(flag, args, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0L) return(default)
  args[idx[1L] + 1L]
}

workdir  <- get_flag("--workdir",  args, "~/GG_EPICYPHER_CBX2")
out_dir  <- get_flag("--out_dir",  args, "Analysis_Data/cbx2_deg")

cat("=== cbx2_deg_prep.R ===\n")
cat("Workdir :", workdir, "\n")
cat("Out dir :", out_dir, "\n")

setwd(workdir)

# Output directories
bed_dir <- file.path(out_dir, "bed_files")
rds_dir <- file.path(out_dir, "rds")
dir.create(bed_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(rds_dir, showWarnings = FALSE, recursive = TRUE)

# --- Thresholds ---
LFC_THRESH <- 0.585   # log2(1.5): at least 50% expression increase
FDR_THRESH <- 0.05

# --- Load DEG tables ---
deg_ko_dmso    <- read.csv("Importable_Data/rnaseq/results_CBX2_DMSO_vs_LUC_DMSO.csv", stringsAsFactors = FALSE)
deg_ezh2i      <- read.csv("Importable_Data/rnaseq/results_LUC_1uM_vs_LUC_DMSO.csv", stringsAsFactors = FALSE)
deg_ko_1um_vs_luc_1um <- read.csv("Importable_Data/rnaseq/results_CBX2_1uM_vs_LUC_1uM.csv", stringsAsFactors = FALSE)

cat(sprintf("\nDEG tables loaded.\n  KO+DMSO table: %d genes\n  WT+EZH2i table:  %d genes\n  KO+EZH2i (vs LUC+EZH2i) table: %d genes\n",
            nrow(deg_ko_dmso), nrow(deg_ezh2i), nrow(deg_ko_1um_vs_luc_1um)))

# --- Helper functions for classification ---
is_up <- function(df, lfc, fdr) {
  !is.na(df$padj) & df$log2FoldChange >= lfc & df$padj <= fdr
}

up_ko_dmso  <- deg_ko_dmso$Gene[is_up(deg_ko_dmso, LFC_THRESH, FDR_THRESH)]
up_ezh2i    <- deg_ezh2i$Gene[is_up(deg_ezh2i, LFC_THRESH, FDR_THRESH)]
# Group C: CBX2 KO provides significant extra upregulation ON TOP OF EZH2i alone.
# Comparison: CBX2 KO + EZH2i vs LUC + EZH2i (isolates the CBX2 effect under EZH2i).
up_ko_ezh2i <- deg_ko_1um_vs_luc_1um$Gene[is_up(deg_ko_1um_vs_luc_1um, LFC_THRESH, FDR_THRESH)]

# NOTE: Groups are NON-EXCLUSIVE. A gene can belong to A AND C, B AND C, etc.
# Each group is defined independently:
# Group A (CBX2-sensitive at DMSO): Upregulated when CBX2 is lost in standard conditions
grp_A <- up_ko_dmso

# Group B (EZH2i-sensitive): Upregulated by EZH2i alone
grp_B <- up_ezh2i

# Group C (CBX2-sensitive at EZH2i): CBX2 KO provides a significant boost OVER EZH2i alone
grp_C <- up_ko_ezh2i

# Group D (NS Background): NOT in any of A, B, or C
all_genes <- unique(c(deg_ko_dmso$Gene, deg_ezh2i$Gene, deg_ko_1um_vs_luc_1um$Gene))
grp_D <- setdiff(all_genes, c(grp_A, grp_B, grp_C))

cat(sprintf("\nGene classification (LFC >= %.3f, FDR <= %.2f):\n", LFC_THRESH, FDR_THRESH))
cat(sprintf("  Group A (CBX2-sensitive at DMSO):        %d genes\n", length(grp_A)))
cat(sprintf("  Group B (EZH2i-sensitive):               %d genes\n", length(grp_B)))
cat(sprintf("  Group C (CBX2-sensitive at EZH2i):       %d genes\n", length(grp_C)))
cat(sprintf("  Group D (NS background):                 %d genes\n", length(grp_D)))
cat(sprintf("  Overlap A & C:                           %d genes\n", length(intersect(grp_A, grp_C))))
cat(sprintf("  Overlap B & C:                           %d genes\n", length(intersect(grp_B, grp_C))))

# --- Load coordinates from cbx2_expression_master.rds ---
master_path <- "Analysis_Data/cbx2_expression/csaw/cbx2_expression_master.rds"
if (!file.exists(master_path)) stop("cbx2_expression_master.rds not found! Run cbx2_expression_prep.R first.")
master_df <- readRDS(master_path)
cat(sprintf("\nLoaded cbx2_expression_master.rds (%d genes with coordinates)\n", nrow(master_df)))

# --- Augment master with CBX2 KO + EZH2i expression (not stored in original master) ---
rlog_path <- "Importable_Data/rds/rnaseq/rlogData.rds"
rlog_data  <- readRDS(rlog_path)
if (inherits(rlog_data, "SummarizedExperiment")) rlog_data <- SummarizedExperiment::assay(rlog_data)
ko_1um_cols <- c("CBX2KO_Clone1_RepA_1uM", "CBX2KO_Clone1_RepB_1uM",
                  "CBX2KO_Clone3_RepA_1uM", "CBX2KO_Clone3_RepB_1uM")  # replace with your sample IDs
if (all(ko_1um_cols %in% colnames(rlog_data))) {
  exp_cbx2ko_1um_vec <- rowMeans(rlog_data[, ko_1um_cols], na.rm = TRUE)
  exp_cbx2ko_1um_df  <- data.frame(
    gene_name      = rownames(rlog_data),
    exp_cbx2ko_1um = exp_cbx2ko_1um_vec,
    stringsAsFactors = FALSE
  )
  master_df <- left_join(master_df, exp_cbx2ko_1um_df, by = "gene_name")
  cat(sprintf("Added exp_cbx2ko_1um column (CBX2 KO + EZH2i rlogCPM)\n"))
} else {
  warning("CBX2 KO+1uM columns not found in rlogData — exp_cbx2ko_1um will be NA")
  master_df$exp_cbx2ko_1um <- NA_real_
}

# --- Load DEG logFC columns for each gene for the master table ---
deg_info_dmso <- deg_ko_dmso[, c("Gene", "log2FoldChange", "padj", "baseMean", "effect")]
colnames(deg_info_dmso) <- c("gene_name", "lfc_cbx2ko_dmso", "fdr_cbx2ko_dmso",
                              "baseMean_dmso", "effect_dmso")
deg_info_1um <- deg_ko_1um_vs_luc_1um[, c("Gene", "log2FoldChange", "padj", "effect")]
colnames(deg_info_1um) <- c("gene_name", "lfc_cbx2ko_1um", "fdr_cbx2ko_1um", "effect_1um")

# Merge with master
df <- master_df %>%
  left_join(deg_info_dmso, by = "gene_name") %>%
  left_join(deg_info_1um,  by = "gene_name")

# Assign NON-EXCLUSIVE group membership as a comma-separated list
df$in_grp_A <- df$gene_name %in% grp_A
df$in_grp_B <- df$gene_name %in% grp_B
df$in_grp_C <- df$gene_name %in% grp_C
df$in_grp_D <- df$gene_name %in% grp_D

df$deg_groups <- mapply(function(a, b, c2, d) {
  grps <- c(if (a) "A", if (b) "B", if (c2) "C", if (d) "D")
  if (length(grps) == 0) "D" else paste(grps, collapse=",")
}, df$in_grp_A, df$in_grp_B, df$in_grp_C, df$in_grp_D)

# For single-group plotting convenience, assign priority: C > A > B > D
df$deg_group <- "D"
df$deg_group[df$in_grp_B] <- "B"
df$deg_group[df$in_grp_A] <- "A"
df$deg_group[df$in_grp_C] <- "C"  # C takes highest priority

cat("\nFinal group sizes in master dataframe (priority assignment):\n")
print(table(df$deg_group))
cat("\nNon-exclusive group membership breakdown:\n")
print(table(df$deg_groups))

# --- Helper: write BED ---
write_bed <- function(sub_df, out_path, use_promoters = TRUE) {
  if (use_promoters) {
    # Promoter: we stored promoter coords in cbx2_expression_prep.R via promoters()
    # Reconstruct from gene coords: upstream 3000, downstream 3000 of TSS
    tss_pos <- ifelse(sub_df$strand == "+", sub_df$start, sub_df$end)
    prom_start <- ifelse(sub_df$strand == "+", tss_pos - 3000, tss_pos - 3000)
    prom_end   <- ifelse(sub_df$strand == "+", tss_pos + 3000, tss_pos + 3000)
    prom_start <- pmax(prom_start - 1, 0)  # BED is 0-based
    bed <- data.frame(
      chr    = sub_df$chr,
      start  = prom_start,
      end    = prom_end,
      name   = sub_df$gene_name,
      score  = ".",
      strand = sub_df$strand
    )
  } else {
    bed <- data.frame(
      chr    = sub_df$chr,
      start  = sub_df$start - 1,
      end    = sub_df$end,
      name   = sub_df$gene_name,
      score  = ".",
      strand = sub_df$strand
    )
  }
  write.table(bed, out_path, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
  cat(sprintf("  Wrote %d regions -> %s\n", nrow(bed), out_path))
}

# --- Write BED files for each group ---
# Groups are NON-EXCLUSIVE: a gene appears in BED for every group it qualifies for.
cat("\nWriting promoter BED files:\n")
for (grp in c("A", "B", "C", "D")) {
  flag_col <- paste0("in_grp_", grp)
  sub <- df[df[[flag_col]] == TRUE, ]
  write_bed(sub, file.path(bed_dir, sprintf("promoters_group%s.bed", grp)), use_promoters = TRUE)
}

cat("\nWriting gene body BED files:\n")
for (grp in c("A", "B", "C", "D")) {
  flag_col <- paste0("in_grp_", grp)
  sub <- df[df[[flag_col]] == TRUE, ]
  write_bed(sub, file.path(bed_dir, sprintf("genebodies_group%s.bed", grp)), use_promoters = FALSE)
}

# --- Save master RDS ---
out_rds <- file.path(rds_dir, "cbx2_deg_master.rds")
saveRDS(df, out_rds)
cat(sprintf("\nSaved cbx2_deg_master.rds to %s\n", out_rds))
cat("=== Prep completed successfully ===\n")
