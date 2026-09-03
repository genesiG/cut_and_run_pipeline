#!/usr/bin/env Rscript

# =============================================================================
# step_12a_prepare_gene_universe.R
# Builds unified gene reference table & promoter BEDs for CBX2 analysis
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(rtracklayer)
  library(GenomicRanges)
})

work_dir <- "~/GG_EPICYPHER_CBX2"
setwd(work_dir)

out_dir <- file.path(work_dir, "Analysis_Data", "perturbation")
bed_dir <- file.path(out_dir, "bed")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(bed_dir, recursive = TRUE, showWarnings = FALSE)

# 1. Load RNA-seq DESeq2 results
rnaseq_dir <- file.path(work_dir, "Importable_Data", "rnaseq")

cbx2_ko_file <- file.path(rnaseq_dir, "results_CBX2_DMSO_vs_LUC_DMSO.csv")
ezh2i_file   <- file.path(rnaseq_dir, "results_LUC_1uM_vs_LUC_DMSO.csv")
cbx2_ezh2i_file <- file.path(rnaseq_dir, "results_CBX2_1uM_vs_LUC_DMSO.csv")
cbx2_1uM_vs_ezh2i_file <- file.path(rnaseq_dir, "results_CBX2_1uM_vs_LUC_1uM.csv")

cat("Loading RNA-seq DEG results...\n")
res_cbx2_ko <- read_csv(cbx2_ko_file, show_col_types = FALSE) %>%
  dplyr::rename(
    gene = Gene,
    baseMean_cbx2 = baseMean,
    log2FC_CBX2KO = log2FoldChange,
    lfcSE_CBX2KO = lfcSE,
    stat_CBX2KO = stat,
    pvalue_CBX2KO = pvalue,
    padj_CBX2KO = padj,
    shrunkLFC_CBX2KO = shrunkLFC,
    tier_CBX2KO = effect
  )

res_ezh2i <- read_csv(ezh2i_file, show_col_types = FALSE) %>%
  dplyr::rename(
    gene = Gene,
    baseMean_ezh2i = baseMean,
    log2FC_EZH2i = log2FoldChange,
    lfcSE_EZH2i = lfcSE,
    stat_EZH2i = stat,
    pvalue_EZH2i = pvalue,
    padj_EZH2i = padj,
    shrunkLFC_EZH2i = shrunkLFC,
    tier_EZH2i = effect
  )

res_dual <- read_csv(cbx2_ezh2i_file, show_col_types = FALSE) %>%
  dplyr::rename(
    gene = Gene,
    log2FC_Dual = log2FoldChange,
    padj_Dual = padj,
    tier_Dual = effect
  ) %>%
  dplyr::select(gene, log2FC_Dual, padj_Dual, tier_Dual)

# Merge RNA-seq results
rna_df <- res_cbx2_ko %>%
  inner_join(res_ezh2i %>% dplyr::select(gene, baseMean_ezh2i, log2FC_EZH2i, lfcSE_EZH2i, stat_EZH2i, pvalue_EZH2i, padj_EZH2i, shrunkLFC_EZH2i, tier_EZH2i), by = "gene") %>%
  left_join(res_dual, by = "gene")

# Ensure tier classifications are clean and explicit: Upregulated, Downregulated, Unchanged
rna_df <- rna_df %>%
  mutate(
    tier_CBX2KO = case_when(
      !is.na(padj_CBX2KO) & padj_CBX2KO < 0.05 & log2FC_CBX2KO > 1 ~ "Upregulated",
      !is.na(padj_CBX2KO) & padj_CBX2KO < 0.05 & log2FC_CBX2KO < -1 ~ "Downregulated",
      TRUE ~ "Unchanged"
    ),
    tier_EZH2i = case_when(
      !is.na(padj_EZH2i) & padj_EZH2i < 0.05 & log2FC_EZH2i > 1 ~ "Upregulated",
      !is.na(padj_EZH2i) & padj_EZH2i < 0.05 & log2FC_EZH2i < -1 ~ "Downregulated",
      TRUE ~ "Unchanged"
    )
  )

cat("  Loaded", nrow(rna_df), "genes with RNA-seq data.\n")
cat("  CBX2 KO tiers:\n")
print(table(rna_df$tier_CBX2KO))
cat("  EZH2i tiers:\n")
print(table(rna_df$tier_EZH2i))

# 2. Load TSS BED
tss_bed_file <- "/path/to/reference/t2t.curated.norandom.tss.bed"
cat("Loading TSS annotations from:", tss_bed_file, "\n")
tss_df <- read_tsv(tss_bed_file, col_names = c("chr", "start", "end", "strand", "transcript_id", "gene_name", "gene_symbol"), show_col_types = FALSE)

# Format chromosome names to match BAM / BigWig files (chr prefix)
tss_df <- tss_df %>%
  mutate(
    chr = ifelse(grepl("^chr", chr), chr, paste0("chr", chr)),
    # Calculate exact single-base TSS coordinate
    tss = ifelse(strand == "+", start, end)
  )

# Collapse multiple transcripts per gene symbol to canonical / first entry
tss_unique <- tss_df %>%
  distinct(gene_symbol, .keep_all = TRUE) %>%
  dplyr::select(gene_symbol, chr, strand, tss, transcript_id)

cat("  Loaded", nrow(tss_df), "TSS entries (", nrow(tss_unique), "unique gene symbols).\n")

# 3. Merge RNA-seq and TSS data
gene_universe <- rna_df %>%
  inner_join(tss_unique, by = c("gene" = "gene_symbol")) %>%
  mutate(
    # Promoter windows: 2kb and 2.5kb around TSS
    promoter_2kb_start = pmax(0, tss - 2000),
    promoter_2kb_end   = tss + 2000,
    promoter_2.5kb_start = pmax(0, tss - 2500),
    promoter_2.5kb_end   = tss + 2500,
    # Baseline expression: average baseMean from DMSO controls
    baseline_expr = baseMean_cbx2
  ) %>%
  arrange(desc(baseline_expr))

cat("  Final matched gene universe count:", nrow(gene_universe), "genes.\n")

# 4. Export TSV
univ_tsv_path <- file.path(out_dir, "gene_universe.tsv")
write_tsv(gene_universe, univ_tsv_path)
cat("Exported gene universe to:", univ_tsv_path, "\n")

# 5. Export Promoter BED files
export_promoter_bed <- function(df, start_col, end_col, filename, score_col = "baseline_expr") {
  bed_df <- df %>%
    transmute(
      chrom = chr,
      chromStart = as.integer(df[[start_col]]),
      chromEnd = as.integer(df[[end_col]]),
      name = gene,
      score = as.integer(pmin(1000, round(df[[score_col]]))),
      strand = strand
    )
  write.table(bed_df, file.path(bed_dir, filename), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

# 2kb promoter BEDs (for CSAW / signal quantification)
export_promoter_bed(gene_universe, "promoter_2kb_start", "promoter_2kb_end", "promoters_2kb_all.bed")
export_promoter_bed(gene_universe %>% filter(tier_CBX2KO == "Upregulated"), "promoter_2kb_start", "promoter_2kb_end", "promoters_2kb_CBX2KO_upregulated.bed")
export_promoter_bed(gene_universe %>% filter(tier_CBX2KO == "Downregulated"), "promoter_2kb_start", "promoter_2kb_end", "promoters_2kb_CBX2KO_downregulated.bed")
export_promoter_bed(gene_universe %>% filter(tier_CBX2KO == "Unchanged"), "promoter_2kb_start", "promoter_2kb_end", "promoters_2kb_CBX2KO_unchanged.bed")

# 2.5kb promoter BEDs (for deepTools heatmaps)
export_promoter_bed(gene_universe, "promoter_2.5kb_start", "promoter_2.5kb_end", "promoters_2.5kb_all.bed")

cat("Exported promoter BED files to:", bed_dir, "\n")
cat("=== step_12a complete ===\n")
