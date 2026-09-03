#!/usr/bin/env Rscript

# =============================================================================
# step_12b_csaw_cbx2_differential.R
# Sliding-window differential CBX2 occupancy across promoters (TSS ± 2kb)
# using csaw + spike-in normalization (EZH2i vs DMSO)
# =============================================================================

suppressPackageStartupMessages({
  library(csaw)
  library(edgeR)
  library(GenomicRanges)
  library(rtracklayer)
  library(dplyr)
  library(readr)
})

work_dir <- "~/GG_EPICYPHER_CBX2"
setwd(work_dir)

out_dir <- file.path(work_dir, "Analysis_Data", "perturbation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 1. Load Gene Universe / Promoters
univ_path <- file.path(out_dir, "gene_universe.tsv")
if (!file.exists(univ_path)) {
  stop("gene_universe.tsv not found! Run step_12a first.")
}
cat("Loading gene universe from:", univ_path, "\n")
gene_df <- read_tsv(univ_path, show_col_types = FALSE)

# Build GRanges for promoters (TSS ± 2kb)
promoter_gr <- GRanges(
  seqnames = gene_df$chr,
  ranges   = IRanges(start = gene_df$promoter_2kb_start, end = gene_df$promoter_2kb_end),
  strand   = gene_df$strand,
  gene     = gene_df$gene
)
names(promoter_gr) <- gene_df$gene

# 2. Define CBX2 BAM files and condition
bam_dir <- file.path(work_dir, "Analysis_Data", "bam", "processed")

bams <- c(
  file.path(bam_dir, "RP_031_AntiCBX2RP_DMSO_R1.qc.sort.markdup.mapq30.final.bam"),
  file.path(bam_dir, "RP_039_AntiCBX2RP_DMSO_R1.qc.sort.markdup.mapq30.final.bam"),
  file.path(bam_dir, "RP_055_AntiCBX2RP_1uMcpd_R1.qc.sort.markdup.mapq30.final.bam"),
  file.path(bam_dir, "RP_063_AntiCBX2RP_1uMcpd_R1.qc.sort.markdup.mapq30.final.bam")
)

sample_ids <- c("RP_031_DMSO", "RP_039_DMSO", "RP_055_EZH2i", "RP_063_EZH2i")
conditions <- factor(c("DMSO", "DMSO", "EZH2i", "EZH2i"), levels = c("DMSO", "EZH2i"))

for (b in bams) {
  if (!file.exists(b)) stop(paste("BAM not found:", b))
}

# 3. Load Spike-in Scale Factors
sf_file <- file.path(work_dir, "Analysis_Data", "normalization", "spikein", "CBX2_spikein_SF.txt")
cat("Loading Spike-in factors from:", sf_file, "\n")
sf_df <- read_tsv(sf_file, show_col_types = FALSE)

# Match scale factors in order of bams
bam_basenames <- basename(bams)
sf_matched <- sf_df %>%
  slice(match(bam_basenames, ID))

cat("Matched spike-in scale factors:\n")
print(sf_matched)

# 4. Count reads in sliding windows within promoter regions
cat("\nCounting reads in sliding windows (150 bp window, 50 bp spacing)...\n")
param <- readParam(
  minq = 30,
  pe = "both",
  dedup = FALSE
)

# Extract sliding windows
win_counts <- windowCounts(
  bams,
  width = 150,
  spacing = 50,
  param = param
)

# Subset windows overlapping promoter regions
olaps <- findOverlaps(win_counts, promoter_gr)
win_promoters <- win_counts[queryHits(olaps), ]
mcols(win_promoters)$gene <- promoter_gr$gene[subjectHits(olaps)]

cat("Total sliding windows overlapping promoters:", length(win_promoters), "\n")

# Also count full promoter region counts
cat("Counting reads across full 4kb promoter windows...\n")
region_counts <- regionCounts(
  bams,
  regions = promoter_gr,
  param = param
)

# 5. Apply Spike-in Normalization
# In edgeR: effective_lib_size = lib.size * norm.factors
# Since SF = min(raw_SF)/raw_SF is divisive, norm.factors = 1 / SF = bamCov_SF
cat("Applying spike-in normalization factors...\n")
spike_norm_factors <- sf_matched$bamCov_SF

y_win <- asDGEList(win_promoters)
y_win$samples$norm.factors <- spike_norm_factors

y_reg <- asDGEList(region_counts)
y_reg$samples$norm.factors <- spike_norm_factors

# 6. Fit edgeR model and test contrast (EZH2i vs DMSO)
cat("Fitting edgeR model and testing differential binding...\n")
design <- model.matrix(~conditions)
colnames(design) <- c("Intercept", "EZH2i_vs_DMSO")

# Estimate dispersion
y_win <- estimateDisp(y_win, design)
fit_win <- glmQLFit(y_win, design, robust = TRUE)
res_win <- glmQLFTest(fit_win, coef = "EZH2i_vs_DMSO")

y_reg <- estimateDisp(y_reg, design)
fit_reg <- glmQLFit(y_reg, design, robust = TRUE)
res_reg <- glmQLFTest(fit_reg, coef = "EZH2i_vs_DMSO")

# 7. Merge sliding-window results per gene promoter
cat("Merging sliding windows per gene promoter...\n")
win_to_promoter <- findOverlaps(rowRanges(win_promoters), promoter_gr)

olap_df <- data.frame(
  prom_id = subjectHits(win_to_promoter),
  win_lfc = res_win$table$logFC[queryHits(win_to_promoter)],
  win_p   = res_win$table$PValue[queryHits(win_to_promoter)]
)

win_summary <- olap_df %>%
  group_by(prom_id) %>%
  summarise(
    win_best_log2FC = win_lfc[which.min(win_p)],
    win_min_pval = min(win_p),
    win_num = n(),
    .groups = "drop"
  )

# Map window summaries to full promoter universe
win_best_lfc <- rep(NA_real_, length(promoter_gr))
win_min_p    <- rep(NA_real_, length(promoter_gr))
win_n        <- rep(0L, length(promoter_gr))

win_best_lfc[win_summary$prom_id] <- win_summary$win_best_log2FC
win_min_p[win_summary$prom_id]    <- win_summary$win_min_pval
win_n[win_summary$prom_id]        <- win_summary$win_num

# Combine window-level summaries with region-level counts
res_table <- data.frame(
  gene = promoter_gr$gene,
  chr = as.character(seqnames(promoter_gr)),
  promoter_start = start(promoter_gr),
  promoter_end = end(promoter_gr),
  log2FC_CBX2_binding = res_reg$table$logFC,
  logCPM_CBX2 = res_reg$table$logCPM,
  logCPM_CBX2_EZH2i = rowMeans(cpm(y_reg, log = TRUE)[, conditions == "EZH2i", drop = FALSE]),
  pvalue_CBX2_binding = res_reg$table$PValue,
  FDR_CBX2_binding = p.adjust(res_reg$table$PValue, method = "BH"),
  win_best_log2FC = win_best_lfc,
  win_min_pval = win_min_p,
  win_num = win_n
)

# 8. Export differential binding table
out_tsv <- file.path(out_dir, "csaw_cbx2_promoter_lfc.tsv")
write_tsv(res_table, out_tsv)
cat("Saved differential CBX2 binding to:", out_tsv, "\n")

out_rds <- file.path(out_dir, "csaw_cbx2_promoter_lfc.rds")
saveRDS(list(res_table = res_table, win_counts = win_promoters, region_counts = region_counts, res_win = res_win, res_reg = res_reg), out_rds)

cat("=== step_12b complete ===\n")
