#!/usr/bin/env Rscript

# =============================================================================
# step_13b_chromatin_states.R
# Analysis 2.1: Partition promoters (TSS ± 2.5kb) into 4 chromatin states:
#   1. Canonical Polycomb: H3K27me3(+) / H2AK119ub(+)
#   2. Non-canonical PRC1: H3K27me3(-) / H2AK119ub(+)
#   3. PRC2-only / Primed: H3K27me3(+) / H2AK119ub(-)
#   4. Unmodified / Active: H3K27me3(-) / H2AK119ub(-)
# =============================================================================

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(dplyr)
  library(readr)
})

work_dir <- "~/GG_EPICYPHER_CBX2"
setwd(work_dir)

data_dir  <- file.path(work_dir, "Analysis_Data", "perturbation")
state_dir <- file.path(data_dir, "chromatin_states")
peak_dir  <- file.path(work_dir, "Analysis_Data", "peaks", "csaw")
dir.create(state_dir, recursive = TRUE, showWarnings = FALSE)

# 1. Load Gene Universe
univ_file <- file.path(data_dir, "gene_universe.tsv")
if (!file.exists(univ_file)) {
  stop("gene_universe.tsv not found! Run step_12a first.")
}
gene_df <- read_tsv(univ_file, show_col_types = FALSE)

# Build promoter GRanges (TSS ± 2.5kb)
promoter_gr <- GRanges(
  seqnames = gene_df$chr,
  ranges   = IRanges(start = gene_df$promoter_2.5kb_start, end = gene_df$promoter_2.5kb_end),
  strand   = gene_df$strand,
  gene     = gene_df$gene,
  baseline_expr = gene_df$baseline_expr
)
names(promoter_gr) <- gene_df$gene

# 2. Find Peak BED files
k27me3_beds <- list.files(peak_dir, pattern = "^H3K27me3_DMSO.*\\.bed$", full.names = TRUE)
k119ub_beds <- list.files(peak_dir, pattern = "^H2AK119ub_DMSO.*\\.bed$", full.names = TRUE)
cbx2_beds   <- list.files(peak_dir, pattern = "^AntiCBX2RP_DMSO.*\\.bed$", full.names = TRUE)

cat("Found peak files in csaw directory:\n")
cat("  H3K27me3 peaks:", ifelse(length(k27me3_beds) > 0, basename(k27me3_beds[1]), "NONE"), "\n")
cat("  H2AK119ub peaks:", ifelse(length(k119ub_beds) > 0, basename(k119ub_beds[1]), "NONE"), "\n")
cat("  CBX2 peaks    :", ifelse(length(cbx2_beds) > 0, basename(cbx2_beds[1]), "NONE"), "\n")

if (length(k27me3_beds) == 0 || length(k119ub_beds) == 0) {
  cat("\nWarning: Waiting for csaw peak calling to finish if still running.\n")
}

# Load peaks
k27me3_gr <- if (length(k27me3_beds) > 0) import(k27me3_beds[1]) else GRanges()
k119ub_gr <- if (length(k119ub_beds) > 0) import(k119ub_beds[1]) else GRanges()
cbx2_gr   <- if (length(cbx2_beds) > 0)   import(cbx2_beds[1])   else GRanges()

# Ensure chromosome styles match
seqlevelsStyle(k27me3_gr) <- "UCSC"
seqlevelsStyle(k119ub_gr) <- "UCSC"
seqlevelsStyle(cbx2_gr)   <- "UCSC"

# 3. Intersect promoters with peaks
has_k27me3 <- if (length(k27me3_gr) > 0) overlapsAny(promoter_gr, k27me3_gr) else rep(FALSE, length(promoter_gr))
has_k119ub <- if (length(k119ub_gr) > 0) overlapsAny(promoter_gr, k119ub_gr) else rep(FALSE, length(promoter_gr))
has_cbx2   <- if (length(cbx2_gr) > 0)   overlapsAny(promoter_gr, cbx2_gr)   else rep(FALSE, length(promoter_gr))

# 4. Assign 4 Chromatin States
state_label <- case_when(
  has_k27me3 & has_k119ub   ~ "Canonical_Polycomb",
  !has_k27me3 & has_k119ub  ~ "Noncanonical_PRC1",
  has_k27me3 & !has_k119ub  ~ "PRC2_only",
  TRUE                      ~ "Unmodified_Active"
)

promoter_df <- gene_df %>%
  mutate(
    has_H3K27me3_peak = has_k27me3,
    has_H2AK119ub_peak = has_k119ub,
    has_CBX2_peak = has_cbx2,
    chromatin_state = state_label
  )

cat("\nChromatin State Summary:\n")
print(table(promoter_df$chromatin_state))

# 5. Export state TSV table
state_tsv <- file.path(data_dir, "chromatin_states.tsv")
write_tsv(promoter_df, state_tsv)
cat("\nSaved chromatin state assignments to:", state_tsv, "\n")

# 6. Export sorted BED file for each state (ordered by baseline expression descending)
states <- c("Canonical_Polycomb", "Noncanonical_PRC1", "PRC2_only", "Unmodified_Active")

for (st in states) {
  st_df <- promoter_df %>%
    filter(chromatin_state == st) %>%
    arrange(desc(baseline_expr)) %>%
    transmute(
      chrom = chr,
      chromStart = as.integer(promoter_2.5kb_start),
      chromEnd = as.integer(promoter_2.5kb_end),
      name = gene,
      score = as.integer(pmin(1000, round(baseline_expr))),
      strand = strand
    )
  
  st_bed <- file.path(state_dir, paste0("state_", st, ".bed"))
  write.table(st_df, st_bed, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  cat(sprintf("  Exported %s: %d promoters -> %s\n", st, nrow(st_df), st_bed))
}

cat("=== step_13b complete ===\n")
