#!/usr/bin/env Rscript

# =============================================================================
# classify_regions_lfc.R
#
# Generic classifier: divides pre-specified genomic regions into Retained
# and Lost based on spike-in normalised log2 fold-change (1uMcpd vs DMSO).
#
# Uses csaw::regionCounts over pre-specified peaks and applies pre-computed
# spike-in normalization factors from a {TARGET}_spikein_SF.txt file via edgeR.
#
# Classification criteria
# -----------------------
#   Lost:     logFC (1uM vs DMSO) < lfc_lost   AND  FDR < fdr_threshold
#   Retained: logFC (1uM vs DMSO) > lfc_retained  (complement; no FDR filter)
#   Unclassified: everything else (excluded from BED outputs)
#
# CLI flags (all have defaults; see below)
# ----------------------------------------
#   --target        STRING   Antibody target label (e.g. K27me2, K27me3).
#                            Used to auto-discover metadata and SF files.
#   --metadata      PATH     Tab-delimited file with columns: ID, ANTIBODY, GROUP
#                            (GROUP must contain "DMSO" and "1uM" levels).
#   --regions       PATH     Pre-specified BED regions to classify.
#   --sf_file       PATH     Spike-in SF file (ID, raw_SF, SF, bamCov_SF columns).
#   --out_dir       PATH     Directory for output files.
#   --lfc_retained  NUMBER   logFC threshold above which a region is Retained
#                            (1uM vs DMSO; default -1.0).
#   --lfc_lost      NUMBER   logFC threshold below which a region is Lost
#                            (1uM vs DMSO; default -1.0).
#   --fdr           NUMBER   FDR threshold for Lost classification (default 0.1).
# =============================================================================

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(csaw)
  library(edgeR)
  library(dplyr)
})

# ---------------------------------------------------------------------------
# CLI parsing
# ---------------------------------------------------------------------------
args_cli <- commandArgs(trailingOnly = TRUE)

get_flag <- function(flag, args, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0 || idx[1] + 1 > length(args)) return(default)
  args[idx[1] + 1]
}

target       <- get_flag("--target",       args_cli, "K27me2")
metadata_path<- get_flag("--metadata",     args_cli,
                          file.path("Metadata",
                                    sprintf("sample_metadata_%s_processed.txt", target)))
regions_path <- get_flag("--regions",      args_cli,
                          file.path("Analysis_Data", "peaks", "csaw",
                                    sprintf("H3K27me2_DMSO.w150.d50.filt2.w2000.d500.filt1.lfc1.merge100.tmm.bed")))
sf_path      <- get_flag("--sf_file",      args_cli,
                          file.path("Analysis_Data", "normalization", "spikein",
                                    sprintf("%s_spikein_SF.txt", target)))
out_dir      <- get_flag("--out_dir",      args_cli,
                          file.path("Analysis_Data",
                                    sprintf("%s_classification", tolower(target))))
lfc_retained <- as.numeric(get_flag("--lfc_retained", args_cli, "-1.0"))
lfc_lost     <- as.numeric(get_flag("--lfc_lost",     args_cli, "-1.0"))
fdr_thresh   <- as.numeric(get_flag("--fdr",          args_cli, "0.1"))

cat(sprintf("=== classify_regions_lfc.R ===\n"))
cat(sprintf("  Target        : %s\n", target))
cat(sprintf("  Metadata      : %s\n", metadata_path))
cat(sprintf("  Regions       : %s\n", regions_path))
cat(sprintf("  SF file       : %s\n", sf_path))
cat(sprintf("  Output dir    : %s\n", out_dir))
cat(sprintf("  LFC retained  : > %.2f (no FDR filter)\n", lfc_retained))
cat(sprintf("  LFC lost      : < %.2f AND FDR < %.2f\n", lfc_lost, fdr_thresh))

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# 1. Resolve BAM directory via reticulate → config.py
# ---------------------------------------------------------------------------
# Use the RETICULATE_PYTHON env var if already set (avoids the warning from
# calling use_python() when RETICULATE_PYTHON is already exported by the job).
py_bin <- Sys.getenv("RETICULATE_PYTHON", unset = NA)
if (is.na(py_bin) || !nzchar(py_bin)) {
  py_bin <- Sys.which("python3")
}

suppressPackageStartupMessages(library(reticulate))
# Only call use_python when RETICULATE_PYTHON is NOT set, to avoid the warning.
if (is.na(Sys.getenv("RETICULATE_PYTHON", unset = NA)) ||
    !nzchar(Sys.getenv("RETICULATE_PYTHON", unset = ""))) {
  use_python(py_bin, required = TRUE)
}
py_run_file("Scripts/config.py")

bam_dir      <- py$PROCESSEDBAMDIR
is_pe        <- isTRUE(py$IS_PAIRED_END)

# ---------------------------------------------------------------------------
# 2. Load metadata and discover BAM files
# ---------------------------------------------------------------------------
cat("\nLoading metadata from:", metadata_path, "\n")
if (!file.exists(metadata_path)) stop("Metadata file not found: ", metadata_path)

meta <- read.table(metadata_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# Expect columns: ID, ANTIBODY, GROUP  (GROUP: DMSO | 1uM)

# Reorder so DMSO samples come first, then 1uM
meta <- meta[order(meta$GROUP != "DMSO"), ]

bam_filenames <- meta$ID
bam_paths     <- file.path(bam_dir, bam_filenames)
groups        <- meta$GROUP

missing_bams <- bam_paths[!file.exists(bam_paths)]
if (length(missing_bams) > 0) {
  stop("Missing BAM files:\n", paste(missing_bams, collapse = "\n"))
}
cat(sprintf("Found %d BAM files (%d DMSO, %d 1uM).\n",
            length(bam_filenames),
            sum(groups == "DMSO"),
            sum(groups == "1uM")))

# ---------------------------------------------------------------------------
# 3. Load spike-in scaling factors
# ---------------------------------------------------------------------------
cat("Loading spike-in normalization factors from:", sf_path, "\n")
if (!file.exists(sf_path)) stop("SF file not found: ", sf_path)
sf_df      <- read.table(sf_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
sf_matched <- sf_df[match(bam_filenames, sf_df$ID), ]

if (any(is.na(sf_matched$SF))) {
  stop("Failed to match all BAM filenames to the SF file.\nUnmatched:\n",
       paste(bam_filenames[is.na(sf_matched$SF)], collapse = "\n"))
}
cat("Matched spike-in SFs:\n")
print(data.frame(ID = bam_filenames, GROUP = groups, SF = sf_matched$SF))

# ---------------------------------------------------------------------------
# 4. Load pre-specified regions
# ---------------------------------------------------------------------------
cat("\nLoading regions from:", regions_path, "\n")
if (!file.exists(regions_path)) stop("Regions file not found: ", regions_path)
regions_gr <- rtracklayer::import(regions_path, format = "BED")
cat("Loaded", length(regions_gr), "regions.\n")

# ---------------------------------------------------------------------------
# 5. Count reads in regions (csaw)
# ---------------------------------------------------------------------------
cat("\nCounting reads in regions across all BAMs...\n")
pe_param <- readParam(pe = if (is_pe) "both" else "none")
counts   <- regionCounts(bam_paths, regions_gr, param = pe_param)

# ---------------------------------------------------------------------------
# 6. edgeR GLM with spike-in normalization
# ---------------------------------------------------------------------------
cat("\nFitting edgeR GLM with spike-in normalization...\n")
y         <- asDGEList(counts)
group_fac <- factor(groups, levels = c("DMSO", "1uM"))
y$samples$group        <- group_fac
y$samples$norm.factors <- sf_matched$SF

design <- model.matrix(~ group_fac)
y      <- estimateDisp(y, design)
fit    <- glmQLFit(y, design)

# Contrast: 1uM vs DMSO (coef = 2)
res <- glmQLFTest(fit, coef = 2)
tt  <- res$table

regions_gr$logFC     <- tt$logFC
regions_gr$aveLogCPM <- aveLogCPM(y)
regions_gr$PValue    <- tt$PValue
regions_gr$FDR       <- p.adjust(tt$PValue, method = "BH")

# Save full results CSV
df_all  <- as.data.frame(regions_gr)
csv_out <- file.path(out_dir, "per_region_lfc.csv")
write.csv(df_all, csv_out, row.names = FALSE)
cat("Saved full results to:", csv_out, "\n")

# ---------------------------------------------------------------------------
# 7. Classify into Retained / Lost / Unclassified
# ---------------------------------------------------------------------------
is_lost     <- (regions_gr$logFC <= lfc_lost) & (regions_gr$FDR < fdr_thresh)
is_lost[is.na(is_lost)] <- FALSE

is_retained <- (round(regions_gr$logFC, 1) >= lfc_retained) | (regions_gr$FDR >= fdr_thresh)
is_retained[is.na(is_retained)] <- FALSE

lost_gr         <- regions_gr[is_lost]
retained_gr     <- regions_gr[is_retained]
unclassified_gr <- regions_gr[!is_lost & !is_retained]

cat(sprintf("\nClassification Summary:\n"))
cat(sprintf("  Total input regions                            : %d\n", length(regions_gr)))
cat(sprintf("  Lost     (logFC < %.2f AND FDR < %.2f)        : %d\n",
            lfc_lost, fdr_thresh, length(lost_gr)))
cat(sprintf("  Retained (logFC > %.2f, no FDR filter)        : %d\n",
            lfc_retained, length(retained_gr)))
cat(sprintf("  Unclassified                                   : %d\n",
            length(unclassified_gr)))

# ---------------------------------------------------------------------------
# 8. Export BED files
# ---------------------------------------------------------------------------
export_clean_bed <- function(gr, prefix, filepath) {
  if (length(gr) == 0) {
    warning("No regions to export for: ", filepath)
    return(invisible(NULL))
  }
  bed_df        <- as.data.frame(gr)
  if ("name" %in% colnames(bed_df) && !all(is.na(bed_df$name))) {
    final_names <- bed_df$name
  } else {
    final_names <- paste0(prefix, "_", seq_len(nrow(bed_df)))
  }
  bed_export    <- data.frame(
    chr    = bed_df$seqnames,
    start  = bed_df$start,
    end    = bed_df$end,
    name   = final_names,
    score  = bed_df$logFC,
    strand = "."
  )
  write.table(bed_export, filepath, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
  cat("Exported BED:", filepath, sprintf("(%d regions)\n", nrow(bed_export)))
}

target_lc        <- tolower(target)
retained_bed     <- file.path(out_dir, sprintf("retained_%s.bed",  target_lc))
lost_bed         <- file.path(out_dir, sprintf("lost_%s.bed",      target_lc))
unclassified_bed <- file.path(out_dir, sprintf("unclassified_%s.bed", target_lc))

export_clean_bed(retained_gr,     "retained",     retained_bed)
export_clean_bed(lost_gr,         "lost",         lost_bed)
export_clean_bed(unclassified_gr, "unclassified", unclassified_bed)

cat(sprintf("\n=== classify_regions_lfc.R complete [target: %s] ===\n", target))
