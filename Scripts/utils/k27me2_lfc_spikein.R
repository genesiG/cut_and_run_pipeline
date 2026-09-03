#!/usr/bin/env Rscript

# =============================================================================
# k27me2_lfc_spikein.R
#
# Classifies pre-specified genomic regions into Retained and Lost H3K27me2
# based on log2 fold-change between DMSO and 1uMcpd inhibitor treatment.
#
# Uses csaw::regionCounts over pre-specified peaks and applies pre-computed
# spike-in normalization factors from K27me2_spikein_SF.txt directly in edgeR.
#
# Thresholds (both criteria must be met):
#   Retained H3K27me2: log2FC (1uM vs DMSO) > -1.5  AND  FDR < 0.1
#   Lost H3K27me2:     log2FC (1uM vs DMSO) < -1.5  AND  FDR < 0.1
#   Unclassified:      FDR >= 0.1 (excluded from both BED outputs)
# =============================================================================

suppressPackageStartupMessages({
  library(reticulate)
  library(GenomicRanges)
  library(rtracklayer)
  library(csaw)
  library(edgeR)
  library(dplyr)
})

args_cli <- commandArgs(trailingOnly = TRUE)

get_flag_value <- function(flag, args, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0) {
    return(default)
  }
  if (idx[1] + 1 > length(args)) {
    return(default)
  }
  args[idx[1] + 1]
}

regions_path <- get_flag_value("--regions", args_cli, "Analysis_Data/peaks/csaw/H3K27me2_DMSO.w150.d50.filt2.w2000.d500.filt1.lfc1.merge100.tmm.bed")
sf_path <- get_flag_value("--sf_file", args_cli, "Analysis_Data/normalization/spikein/K27me2_spikein_SF.txt")
out_dir <- get_flag_value("--out_dir", args_cli, "Analysis_Data/k27me2_classification")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 1. Load config via reticulate
use_python(Sys.which("python"), required = TRUE)
py_run_file("Scripts/config.py")

bam_dir <- py$PROCESSEDBAMDIR

# 2. Define all 4 K27me2 BAM files
bam_filenames <- c(
  "RP_028_H3K27me2_DMSO_R1.qc.sort.markdup.mapq30.final.bam",
  "RP_036_H3K27me2_DMSO_R1.qc.sort.markdup.mapq30.final.bam",
  "RP_052_H3K27me2_1uMcpd_R1.qc.sort.markdup.mapq30.final.bam",
  "RP_060_H3K27me2_1uMcpd_R1.qc.sort.markdup.mapq30.final.bam"
)

bam_paths <- file.path(bam_dir, bam_filenames)
missing_bams <- bam_paths[!file.exists(bam_paths)]
if (length(missing_bams) > 0) {
  stop("Missing BAM files:\n", paste(missing_bams, collapse = "\n"))
}

# 3. Load spike-in SFs
cat("Loading spike-in normalization factors from:", sf_path, "\n")
sf_df <- read.table(sf_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
sf_matched <- sf_df[match(bam_filenames, sf_df$ID), ]

if (any(is.na(sf_matched$SF))) {
  stop("Failed to match all BAM filenames in SF file:", sf_path)
}

cat("Matched spike-in SFs:\n")
print(data.frame(ID = bam_filenames, SF = sf_matched$SF, raw_SF = sf_matched$raw_SF))

# 4. Load pre-specified regions
cat("\nLoading regions from:", regions_path, "\n")
if (!file.exists(regions_path)) {
  stop("Regions file not found:", regions_path)
}
regions_gr <- rtracklayer::import(regions_path, format = "BED")
cat("Loaded", length(regions_gr), "regions.\n")

# 5. Count reads in regions using csaw
cat("\nCounting reads in regions across all 4 BAMs...\n")
pe_param <- readParam(pe = if (py$IS_PAIRED_END) "both" else "none")
counts <- regionCounts(bam_paths, regions_gr, param = pe_param)

# 6. Fit edgeR GLM with spike-in normalization
cat("\nFitting edgeR GLM with spike-in normalization...\n")
y <- asDGEList(counts)
group <- factor(c("DMSO", "DMSO", "1uM", "1uM"), levels = c("DMSO", "1uM"))
y$samples$group <- group

# Set lib.size = 1 and norm.factors = SF so that effective library size is exactly SF
y$samples$lib.size <- 1
y$samples$norm.factors <- sf_matched$SF

design <- model.matrix(~group)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)

# Contrast 1uM vs DMSO (coef = 2)
res <- glmQLFTest(fit, coef = 2)
tt <- res$table

regions_gr$logFC <- tt$logFC
regions_gr$aveLogCPM <- aveLogCPM(y)
regions_gr$PValue <- tt$PValue
regions_gr$FDR <- p.adjust(tt$PValue, method = "BH")

# Save full summary table
df_all <- as.data.frame(regions_gr)
csv_out <- file.path(out_dir, "per_region_lfc.csv")
write.csv(df_all, csv_out, row.names = FALSE)
cat("Saved full results table to:", csv_out, "\n")

# 7. Classify into Retained and Lost
# Lost:     logFC (1uM vs DMSO) < -1.5  AND  FDR < 0.1
# Retained: logFC (1uM vs DMSO) > -1.5  (complement of Lost; no FDR filter)
FDR_THRESHOLD <- 0.1

lost_gr <- regions_gr[regions_gr$logFC < -1.5 & regions_gr$FDR < FDR_THRESHOLD]
retained_gr <- regions_gr[regions_gr$logFC > -1.5]

cat(sprintf("\nClassification Summary:\n"))
cat(sprintf("  Total input regions:                            %d\n", length(regions_gr)))
cat(sprintf("  Lost     (logFC < -1.5 & FDR < %.2f):         %d\n", FDR_THRESHOLD, length(lost_gr)))
cat(sprintf("  Retained (logFC > -1.5, no FDR filter):        %d\n", length(retained_gr)))

# Format and export BED files
export_clean_bed <- function(gr, prefix, filepath) {
  if (length(gr) == 0) {
    warning("No regions to export for:", filepath)
    return()
  }
  bed_df <- as.data.frame(gr)[, c("seqnames", "start", "end", "logFC", "FDR")]
  bed_df$name <- paste0(prefix, "_", seq_len(nrow(bed_df)))
  bed_df$score <- as.integer(round(-10 * log10(pmax(bed_df$FDR, 1e-300))))
  bed_df$score <- pmin(bed_df$score, 1000)

  # Ensure standard 6-column BED or 4-column BED format for deepTools
  # deepTools computeMatrix accepts: chr, start, end, name, score, strand
  bed_export <- data.frame(
    chr = bed_df$seqnames,
    start = bed_df$start,
    end = bed_df$end,
    name = bed_df$name,
    score = bed_df$score,
    strand = "."
  )
  write.table(bed_export, filepath, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  cat("Exported BED:", filepath, sprintf("(%d regions)\n", nrow(bed_export)))
}

retained_bed <- file.path(out_dir, "retained_H3K27me2.bed")
lost_bed <- file.path(out_dir, "lost_H3K27me2.bed")

export_clean_bed(retained_gr, "retained", retained_bed)
export_clean_bed(lost_gr, "lost", lost_bed)

cat("\n=== k27me2_lfc_spikein.R complete ===\n")
