#!/usr/bin/env Rscript

# cbx2_expression_prep.R
# Phase 1 of CBX2 Expression Analysis

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(GenomicFeatures)
  library(rtracklayer)
  library(csaw)
  library(edgeR)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)

get_flag_value <- function(flag, args, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0L) return(default)
  args[idx[1L] + 1L]
}

out_dir       <- get_flag_value("--out_dir", args, "Analysis_Data/cbx2_expression")
workdir       <- get_flag_value("--workdir", args, "~/GG_EPICYPHER_CBX2")
codedir       <- get_flag_value("--codedir", args, "~/GG_EPICYPHER_CBX2/Scripts")

cat("=== cbx2_expression_prep.R ===\n")
cat("Output dir :", out_dir, "\n")

bed_dir <- file.path(out_dir, "bed_files")
dir.create(bed_dir, showWarnings = FALSE, recursive = TRUE)
csaw_dir <- file.path(out_dir, "csaw")
dir.create(csaw_dir, showWarnings = FALSE, recursive = TRUE)

setwd(workdir)

# Hardcode config values to avoid reticulate bugs during DESeq2 loading
GNM_GTF <- "/path/to/reference/t2t.ncbiRefSeq.curated.norandom.gtf"
PROCESSEDBAMDIR <- "~/GG_EPICYPHER_CBX2/Analysis_Data/bam/processed"
IS_PAIRED_END <- TRUE

cat("GTF File :", GNM_GTF, "\n")

txdb <- suppressWarnings(makeTxDbFromGFF(GNM_GTF, format = "gtf", taxonomyId = 9606L))

# 2. Define Promoters (+/- 3kb) and Gene Bodies (TSS to TES)
genes_gr <- suppressWarnings(genes(txdb))
# filter for standard chromosomes
genes_gr <- keepStandardChromosomes(genes_gr, pruning.mode = "coarse")
cat("Total genes in GTF:", length(genes_gr), "\n")

promoters_gr <- suppressWarnings(promoters(genes_gr, upstream=3000, downstream=3000))
cat("Total promoters:", length(promoters_gr), "\n")

# 3. Load rlogData.rds and Map to GTF
rlog_path <- file.path(workdir, "Importable_Data", "rds", "rnaseq", "rlogData.rds")
rlog_data <- readRDS(rlog_path)
if (inherits(rlog_data, "SummarizedExperiment")) {
  rlog_data <- SummarizedExperiment::assay(rlog_data)
}
cat("Loaded rlogData.rds. Rows:", nrow(rlog_data), "Cols:", ncol(rlog_data), "\n")

gene_symbols <- rownames(rlog_data)

mapped_genes <- gene_symbols[gene_symbols %in% names(genes_gr)]
cat("Mapped", length(mapped_genes), "out of", length(gene_symbols), "genes (", round(100*length(mapped_genes)/length(gene_symbols), 1), "%)\n")

rlog_data <- rlog_data[mapped_genes, ]
mapped_gene_ids <- mapped_genes

valid_idx <- which(names(genes_gr) %in% mapped_gene_ids)
genes_gr <- genes_gr[valid_idx]
promoters_gr <- promoters_gr[valid_idx]
mapped_gene_ids <- names(genes_gr)
mapped_gene_names <- names(genes_gr)

rlog_data <- rlog_data[mapped_gene_names, ]

# 4. Compute Average Expression per Condition
exp_control <- rowMeans(rlog_data[, c("Ctrl_Rep1_DMSO", "Ctrl_Rep2_DMSO")], na.rm=TRUE)  # replace with your sample IDs
exp_ezh2i   <- rowMeans(rlog_data[, c("Ctrl_Rep1_1uM", "Ctrl_Rep2_1uM")], na.rm=TRUE)
exp_cbx2ko  <- rowMeans(rlog_data[, c("CBX2KO_Clone1_RepA_DMSO", "CBX2KO_Clone1_RepB_DMSO", "CBX2KO_Clone3_RepA_DMSO", "CBX2KO_Clone3_RepB_DMSO")], na.rm=TRUE)

# 5. Define Terciles and Write BEDs
process_condition <- function(exp_vec, cond_name) {
  q1 <- quantile(exp_vec, 1/3, na.rm=TRUE)
  q3 <- quantile(exp_vec, 2/3, na.rm=TRUE)
  
  terciles <- rep("Medium", length(exp_vec))
  terciles[exp_vec <= q1] <- "Low"
  terciles[exp_vec >= q3] <- "High"
  
  cat(sprintf("\nCondition %s cutoffs: Q1=%.2f, Q3=%.2f\n", cond_name, q1, q3))
  cat(sprintf("  Low: %d, Medium: %d, High: %d\n", sum(terciles=="Low"), sum(terciles=="Medium"), sum(terciles=="High")))
  
  for(terc in c("Low", "Medium", "High")) {
    idx <- which(terciles == terc)
    
    p_df <- data.frame(
      chr = seqnames(promoters_gr[idx]),
      start = start(promoters_gr[idx]) - 1,
      end = end(promoters_gr[idx]),
      name = mapped_gene_names[idx],
      score = ".",
      strand = strand(promoters_gr[idx])
    )
    write.table(p_df, file.path(bed_dir, sprintf("promoters_%s_%s.bed", cond_name, terc)), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
    
    g_df <- data.frame(
      chr = seqnames(genes_gr[idx]),
      start = start(genes_gr[idx]) - 1,
      end = end(genes_gr[idx]),
      name = mapped_gene_names[idx],
      score = ".",
      strand = strand(genes_gr[idx])
    )
    write.table(g_df, file.path(bed_dir, sprintf("genebodies_%s_%s.bed", cond_name, terc)), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  }
  return(terciles)
}

terc_control <- process_condition(exp_control, "control")
terc_ezh2i   <- process_condition(exp_ezh2i, "ezh2i")
terc_cbx2ko  <- process_condition(exp_cbx2ko, "cbx2ko")

bam_dir <- PROCESSEDBAMDIR
# CBX2
bam_cbx2_dmso <- file.path(bam_dir, c("CBX2_DMSO_Rep1.qc.sort.markdup.mapq30.final.bam", "CBX2_DMSO_Rep2.qc.sort.markdup.mapq30.final.bam"))  # replace with your sample IDs
bam_cbx2_1um  <- file.path(bam_dir, c("CBX2_1uM_Rep1.qc.sort.markdup.mapq30.final.bam", "CBX2_1uM_Rep2.qc.sort.markdup.mapq30.final.bam"))

# IgG (Input)
bam_igg_dmso <- file.path(bam_dir, c("IgG_DMSO_Rep1.qc.sort.markdup.mapq30.final.bam", "IgG_DMSO_Rep2.qc.sort.markdup.mapq30.final.bam"))
bam_igg_1um  <- file.path(bam_dir, c("IgG_1uM_Rep1.qc.sort.markdup.mapq30.final.bam", "IgG_1uM_Rep2.qc.sort.markdup.mapq30.final.bam"))

# H3K27me2
bam_me2_dmso <- file.path(bam_dir, c("H3K27me2_DMSO_Rep1.qc.sort.markdup.mapq30.final.bam", "H3K27me2_DMSO_Rep2.qc.sort.markdup.mapq30.final.bam"))
bam_me2_1um  <- file.path(bam_dir, c("H3K27me2_1uM_Rep1.qc.sort.markdup.mapq30.final.bam", "H3K27me2_1uM_Rep2.qc.sort.markdup.mapq30.final.bam"))

# H3K27me3
bam_me3_dmso <- file.path(bam_dir, c("H3K27me3_DMSO_Rep1.qc.sort.markdup.mapq30.final.bam", "H3K27me3_DMSO_Rep2.qc.sort.markdup.mapq30.final.bam"))
bam_me3_1um  <- file.path(bam_dir, c("H3K27me3_1uM_Rep1.qc.sort.markdup.mapq30.final.bam", "H3K27me3_1uM_Rep2.qc.sort.markdup.mapq30.final.bam"))

pe_param <- readParam(pe = if(isTRUE(IS_PAIRED_END)) "both" else "none")

get_logcpm <- function(bams, regions) {
  counts <- regionCounts(bams, regions, param=pe_param)
  y <- asDGEList(counts)
  # Basic logCPM (no spike-in normalization since we just want relative abundance)
  # log=TRUE gives log2 counts per million
  cpm_mat <- cpm(y, log=TRUE)
  res <- rowMeans(cpm_mat)
  rm(counts, y, cpm_mat)
  gc(verbose=FALSE)
  return(res)
}

cat("\nQuantifying over Promoters...\n")
cpm_prom_cbx2_dmso <- get_logcpm(bam_cbx2_dmso, promoters_gr)
cpm_prom_cbx2_1um  <- get_logcpm(bam_cbx2_1um, promoters_gr)
cpm_prom_igg_dmso  <- get_logcpm(bam_igg_dmso, promoters_gr)
cpm_prom_igg_1um   <- get_logcpm(bam_igg_1um, promoters_gr)
cpm_prom_me2_dmso  <- get_logcpm(bam_me2_dmso, promoters_gr)
cpm_prom_me2_1um   <- get_logcpm(bam_me2_1um, promoters_gr)
cpm_prom_me3_dmso  <- get_logcpm(bam_me3_dmso, promoters_gr)
cpm_prom_me3_1um   <- get_logcpm(bam_me3_1um, promoters_gr)

cat("Quantifying over Gene Bodies...\n")
cpm_gene_cbx2_dmso <- get_logcpm(bam_cbx2_dmso, genes_gr)
cpm_gene_cbx2_1um  <- get_logcpm(bam_cbx2_1um, genes_gr)
cpm_gene_igg_dmso  <- get_logcpm(bam_igg_dmso, genes_gr)
cpm_gene_igg_1um   <- get_logcpm(bam_igg_1um, genes_gr)
cpm_gene_me2_dmso  <- get_logcpm(bam_me2_dmso, genes_gr)
cpm_gene_me2_1um   <- get_logcpm(bam_me2_1um, genes_gr)
cpm_gene_me3_dmso  <- get_logcpm(bam_me3_dmso, genes_gr)
cpm_gene_me3_1um   <- get_logcpm(bam_me3_1um, genes_gr)

# Calculate CBX2 vs IgG logFC
lfc_prom_cbx2_dmso <- cpm_prom_cbx2_dmso - cpm_prom_igg_dmso
lfc_prom_cbx2_1um  <- cpm_prom_cbx2_1um - cpm_prom_igg_1um
lfc_gene_cbx2_dmso <- cpm_gene_cbx2_dmso - cpm_gene_igg_dmso
lfc_gene_cbx2_1um  <- cpm_gene_cbx2_1um - cpm_gene_igg_1um

# Helper function to categorize into true terciles
categorize_terciles <- function(vec) {
  q1 <- quantile(vec, 1/3, na.rm=TRUE)
  q3 <- quantile(vec, 2/3, na.rm=TRUE)
  terc <- rep("Medium", length(vec))
  terc[vec <= q1] <- "Low"
  terc[vec >= q3] <- "High"
  factor(terc, levels = c("High", "Medium", "Low"))
}

# Combine into master dataframe
master_df <- data.frame(
  gene_id = mapped_gene_ids,
  gene_name = mapped_gene_names,
  chr = as.character(seqnames(genes_gr)),
  start = start(genes_gr),
  end = end(genes_gr),
  strand = as.character(strand(genes_gr)),
  
  exp_control = exp_control,
  exp_ezh2i = exp_ezh2i,
  exp_cbx2ko = exp_cbx2ko,
  terc_control = terc_control,
  terc_ezh2i = terc_ezh2i,
  terc_cbx2ko = terc_cbx2ko,
  
  cbx2_prom_dmso = lfc_prom_cbx2_dmso,
  cbx2_prom_1um = lfc_prom_cbx2_1um,
  cbx2_gene_dmso = lfc_gene_cbx2_dmso,
  cbx2_gene_1um = lfc_gene_cbx2_1um,
  
  me2_prom_dmso = cpm_prom_me2_dmso,
  me2_prom_1um = cpm_prom_me2_1um,
  me2_gene_dmso = cpm_gene_me2_dmso,
  me2_gene_1um = cpm_gene_me2_1um,
  
  me3_prom_dmso = cpm_prom_me3_dmso,
  me3_prom_1um = cpm_prom_me3_1um,
  me3_gene_dmso = cpm_gene_me3_dmso,
  me3_gene_1um = cpm_gene_me3_1um,
  
  stringsAsFactors = FALSE
)

# Pre-calculate terciles for CBX2, me2, me3 for faster plotting downstream
master_df$terc_cbx2_prom_dmso <- categorize_terciles(master_df$cbx2_prom_dmso)
master_df$terc_cbx2_prom_1um  <- categorize_terciles(master_df$cbx2_prom_1um)
master_df$terc_cbx2_gene_dmso <- categorize_terciles(master_df$cbx2_gene_dmso)
master_df$terc_cbx2_gene_1um  <- categorize_terciles(master_df$cbx2_gene_1um)

master_df$terc_me2_prom_dmso <- categorize_terciles(master_df$me2_prom_dmso)
master_df$terc_me2_prom_1um  <- categorize_terciles(master_df$me2_prom_1um)
master_df$terc_me2_gene_dmso <- categorize_terciles(master_df$me2_gene_dmso)
master_df$terc_me2_gene_1um  <- categorize_terciles(master_df$me2_gene_1um)

master_df$terc_me3_prom_dmso <- categorize_terciles(master_df$me3_prom_dmso)
master_df$terc_me3_prom_1um  <- categorize_terciles(master_df$me3_prom_1um)
master_df$terc_me3_gene_dmso <- categorize_terciles(master_df$me3_gene_dmso)
master_df$terc_me3_gene_1um  <- categorize_terciles(master_df$me3_gene_1um)

out_rds <- file.path(csaw_dir, "cbx2_expression_master.rds")
saveRDS(master_df, out_rds)
cat("\nSaved master dataframe to", out_rds, "\n")
cat("=== Prep completed successfully ===\n")
