#!/usr/bin/env Rscript

# =============================================================================
# utils/peak_annotation.R
#
# Comprehensive peak annotation, metrics computation, and visualization engine
# for CUT&TAG nucleation/spreading sites.
#
# Designed to run as a batch job via bsub (step_7b_peak_annotation.py),
# which allocates sufficient memory for loading Bioconductor annotation
# packages (ChIPseeker, GenomicFeatures, BSgenome).
#
# CLI arguments (required):
#   --bed_files   Space-separated absolute paths to BED files to annotate.
#                 Each BED must be a 0-based 3-column file (chr start end).
#   --labels      Space-separated labels, one per BED file (e.g. Nucleation Spreading).
#                 Order must match --bed_files.
#   --out_dir     Output directory for tables and plots.
#   --workdir     Project root (for sourcing config.py).
#   --codedir     Scripts directory (for sourcing config.py).
#
# CLI arguments (optional — enable EnrichedHeatmap-based TSS/gene-body plots):
#   --bigwig_files    Space-separated paths to bigWig files.
#   --bigwig_labels   Labels for each bigWig (defaults to basename w/o extension).
#   --window_bp       Half-window for heatmaps in bp (default: 3000).
#   --n_bins          Bins per half-window (default: 100).
#
# Outputs (in --out_dir):
#   region_metrics.tsv             — per-region annotation table
#   summary_stats.tsv              — per-label aggregate statistics
#   overlap_matrix.tsv             — pairwise overlap/Jaccard QC
#   plots/
#     tss_windows.bed                (TSS ±3 kb windows, 0-based BED)
#     gene_body.bed                  (gene body coordinates, 0-based BED)
#     01a_tss_enriched_heatmap.pdf   (EnrichedHeatmap, requires --bigwig_files)
#     01b_genebody_enriched_heatmap.pdf (EnrichedHeatmap, requires --bigwig_files)
#     02_avg_profile_comparison.pdf
#     03_annotation_pie_<label>.pdf  (one per label)
#     04_annotation_bar_comparison.pdf
#     05_dist_to_tss.pdf
#     06_width_distribution.pdf
#     07_polycomb_domain_spectrum.pdf
#     08_gc_cpgoe_distribution.pdf   (if BSgenome available)
#     09_chromosome_distribution.pdf
#     10_tss_dist_ecdf.pdf
#     11_chromosome_coverage.pdf
#     12_chromosome_coverage_<label>.pdf (one per label)
#     13_kegg_compareCluster.pdf     (always written; placeholder if no enrichment)
# =============================================================================

# =============================================================================
# 0. Libraries
# =============================================================================
suppressPackageStartupMessages({
  library(reticulate)
  library(GenomicRanges)
  library(GenomicFeatures)
  library(rtracklayer)
  library(Biostrings)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(RColorBrewer)
  library(paletteer)
  library(gridExtra)
})

if (!requireNamespace("ChIPseeker", quietly = TRUE)) {
  stop(paste0(
    "ChIPseeker is required. Update your conda environment:\n",
    "  conda env update -f environment.yml\n",
    "then re-submit this job."
  ))
}
library(ChIPseeker)
options(mc.cores = 4)

CLUSTERP_AVAILABLE <- requireNamespace("clusterProfiler", quietly = TRUE) &&
                      requireNamespace("org.Hs.eg.db",    quietly = TRUE)
if (CLUSTERP_AVAILABLE) {
  suppressPackageStartupMessages({
    library(clusterProfiler)
    library(org.Hs.eg.db)
  })
}

# =============================================================================
# 1. CLI argument parsing
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)

get_flag_values <- function(flag, args, required = TRUE) {
  idx <- which(args == flag)
  if (length(idx) == 0L) {
    if (required) stop(sprintf("Missing required argument: %s", flag))
    return(character(0))
  }
  start <- idx[1L] + 1L
  if (start > length(args)) stop(sprintf("No value supplied for %s", flag))
  end <- start
  while (end <= length(args) && !startsWith(args[end], "--")) end <- end + 1L
  args[start:(end - 1L)]
}

get_flag_value <- function(flag, args, required = TRUE) {
  vals <- get_flag_values(flag, args, required)
  if (length(vals) == 0L) return(NULL)
  vals[1L]
}

bed_files     <- get_flag_values("--bed_files", args)
labels        <- get_flag_values("--labels",    args)
out_dir       <- get_flag_value( "--out_dir",   args)
workdir       <- get_flag_value( "--workdir",   args)
codedir       <- get_flag_value( "--codedir",   args)
bigwig_files  <- get_flag_values("--bigwig_files",  args, required = FALSE)
bigwig_labels <- get_flag_values("--bigwig_labels", args, required = FALSE)
window_bp_arg <- get_flag_value( "--window_bp", args, required = FALSE)
n_bins_arg    <- get_flag_value( "--n_bins",    args, required = FALSE)

HEATMAP_WINDOW_BP <- if (!is.null(window_bp_arg)) as.integer(window_bp_arg) else 3000L
HEATMAP_N_BINS    <- if (!is.null(n_bins_arg))    as.integer(n_bins_arg)    else 100L

if (length(bigwig_files) > 0 && length(bigwig_labels) == 0) {
  bigwig_labels <- tools::file_path_sans_ext(basename(bigwig_files))
}
USE_BIGWIGS <- length(bigwig_files) > 0

if (length(bed_files) != length(labels)) {
  stop(sprintf(
    "--bed_files (%d) and --labels (%d) must have the same length.",
    length(bed_files), length(labels)
  ))
}
if (USE_BIGWIGS && length(bigwig_files) != length(bigwig_labels)) {
  stop(sprintf(
    "--bigwig_files (%d) and --bigwig_labels (%d) must have the same length.",
    length(bigwig_files), length(bigwig_labels)
  ))
}

cat("=== peak_annotation.R ===\n")
cat(sprintf("  Labels    : %s\n", paste(labels, collapse = ", ")))
cat(sprintf("  BED files : %d\n", length(bed_files)))
for (i in seq_along(bed_files)) {
  cat(sprintf("    [%s] %s\n", labels[i], basename(bed_files[i])))
}
cat(sprintf("  Output    : %s\n\n", out_dir))

plots_dir <- file.path(out_dir, "plots")
dir.create(out_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 2. Load config.py via reticulate
# =============================================================================
setwd(workdir)
use_python(Sys.which("python"), required = TRUE)
py_run_file(file.path(codedir, "config.py"))

GNM_GTF <- trimws(py$GNM_GTF)

# =============================================================================
# 3. Read BED files as GRanges
# =============================================================================
cat("=== Reading BED files ===\n")

read_3col_bed <- function(path, label) {
  if (!file.exists(path)) stop(sprintf("BED file not found: %s", path))
  df <- read.table(path, header = FALSE, sep = "\t",
                   col.names = c("chr", "start", "end"),
                   stringsAsFactors = FALSE)
  gr <- makeGRangesFromDataFrame(df,
    seqnames.field          = "chr",
    start.field             = "start",
    end.field               = "end",
    starts.in.df.are.0based = TRUE)
  gr$label <- label
  cat(sprintf("  [%s] %d regions (%.2f Mb)\n",
              label, length(gr), sum(width(gr)) / 1e6))
  gr
}

gr_list <- mapply(read_3col_bed, bed_files, labels,
                  SIMPLIFY = FALSE, USE.NAMES = TRUE)
names(gr_list) <- labels

# =============================================================================
# 4. Build TxDb from project GTF
# =============================================================================
cat("\n=== Building TxDb from GTF ===\n")
cat(sprintf("  GTF: %s\n", GNM_GTF))

if (!file.exists(GNM_GTF)) {
  stop("GTF file not found: ", GNM_GTF, "\n  Check GNM_GTF in Scripts/config.py")
}

txdb <- suppressWarnings(
  makeTxDbFromGFF(GNM_GTF, format = "gtf", taxonomyId = 9606L)
)
cat("  TxDb built successfully.\n")

TSS_UPSTREAM   <- 3000L
TSS_DOWNSTREAM <- 3000L

# =============================================================================
# 5. BSgenome sequence extraction (GC content, CpG O/E)
# =============================================================================
cat("\n=== BSgenome: GC content and CpG O/E ===\n")

BSgenome_AVAILABLE <- FALSE
bsg <- NULL

for (pkg in c("BSgenome.Hsapiens.UCSC.hs1",
               "BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0")) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    library(pkg, character.only = TRUE)
    bsg <- get(pkg)
    BSgenome_AVAILABLE <- TRUE
    cat(sprintf("  Using BSgenome package: %s\n", pkg))
    break
  }
}

if (!BSgenome_AVAILABLE) {
  cat("  BSgenome not available — GC/CpG metrics will be NA.\n")
}

std_chrs <- paste0("chr", c(1:22, "X", "Y", "M"))

compute_seq_metrics <- function(gr, bsg) {
  if (!BSgenome_AVAILABLE) {
    return(data.frame(gc = rep(NA_real_, length(gr)), cpgoe = rep(NA_real_, length(gr)),
                      row.names = seq_along(gr)))
  }
  gr_std <- keepSeqlevels(gr, intersect(seqlevels(gr), std_chrs),
                          pruning.mode = "coarse")
  result <- data.frame(gc = rep(NA_real_, length(gr)), cpgoe = rep(NA_real_, length(gr)),
                       row.names = seq_along(gr))
  if (length(gr_std) == 0L) return(result)
  seqs <- tryCatch({
    gr_seq <- gr_std
    seqlevelsStyle(gr_seq) <- seqlevelsStyle(bsg)[1]
    getSeq(bsg, gr_seq)
  }, error = function(e) {
    warning("getSeq failed: ", conditionMessage(e))
    NULL
  })
  if (is.null(seqs)) return(result)
  bases  <- letterFrequency(seqs, letters = c("C","G"), as.prob = FALSE)
  cpg_f  <- as.numeric(dinucleotideFrequency(seqs, as.prob = FALSE)[, "CG"])
  len    <- as.numeric(width(gr_std))
  gc     <- as.numeric(rowSums(bases)) / len
  nC     <- as.numeric(bases[, "C"]); nG <- as.numeric(bases[, "G"])
  cpgoe  <- ifelse(nC > 0 & nG > 0, (cpg_f * len) / (nC * nG), NA_real_)

  std_idx <- which(as.character(seqnames(gr)) %in% std_chrs)
  result$gc[std_idx]    <- gc
  result$cpgoe[std_idx] <- cpgoe
  result
}

seq_metrics_list <- lapply(gr_list, compute_seq_metrics, bsg = bsg)

# =============================================================================
# 6. ChIPseeker annotation
# =============================================================================
cat("\n=== ChIPseeker annotation ===\n")

anno_list <- lapply(labels, function(lbl) {
  cat(sprintf("  Annotating [%s] ...\n", lbl))
  annotatePeak(
    gr_list[[lbl]],
    tssRegion = c(-TSS_UPSTREAM, TSS_DOWNSTREAM),
    TxDb      = txdb,
    verbose   = FALSE
  )
})
names(anno_list) <- labels
anno_df_list <- lapply(anno_list, as.data.frame)

# =============================================================================
# 7. Assemble region_metrics.tsv
# =============================================================================
cat("\n=== Assembling region_metrics.tsv ===\n")

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b

metrics_df_list <- mapply(function(lbl, gr, anno_df, seq_df) {
  data.frame(
    label           = lbl,
    seqnames        = as.character(seqnames(gr)),
    start           = start(gr),
    end             = end(gr),
    width_bp        = width(gr),
    gc_content      = seq_df$gc,
    cpg_oe          = seq_df$cpgoe,
    annotation      = anno_df$annotation,
    dist_to_tss_bp  = anno_df$distanceToTSS,
    nearest_gene_id = anno_df$geneId,
    gene_biotype    = anno_df$geneBiotype %||% NA_character_,
    stringsAsFactors = FALSE
  )
}, labels, gr_list, anno_df_list, seq_metrics_list, SIMPLIFY = FALSE)

metrics_all <- bind_rows(metrics_df_list)
metrics_all$label <- factor(metrics_all$label, levels = labels)

metrics_path <- file.path(out_dir, "region_metrics.tsv")
write.table(metrics_all, metrics_path, sep = "\t", quote = FALSE,
            row.names = FALSE)
cat(sprintf("  Written: %s  (%d rows)\n", basename(metrics_path), nrow(metrics_all)))

# =============================================================================
# 8. Assemble summary_stats.tsv
# =============================================================================
cat("\n=== Assembling summary_stats.tsv ===\n")

anno_fractions <- function(ann) {
  list(
    pct_promoter   = mean(grepl("Promoter",         ann)) * 100,
    pct_exon       = mean(grepl("Exon",             ann)) * 100,
    pct_intron     = mean(grepl("Intron",           ann)) * 100,
    pct_downstream = mean(grepl("Downstream",       ann)) * 100,
    pct_intergenic = mean(grepl("Intergenic|Distal",ann)) * 100
  )
}

summary_list <- mapply(function(lbl, gr, metrics_df, anno_df, seq_df) {
  w  <- width(gr)
  af <- anno_fractions(anno_df$annotation)
  data.frame(
    label             = lbl,
    n_regions         = length(gr),
    total_bp          = sum(w),
    total_Mb          = round(sum(w) / 1e6, 3),
    median_width_bp   = as.integer(median(w)),
    mean_width_bp     = round(mean(w), 1),
    sd_width_bp       = round(sd(w), 1),
    width_q25_bp      = as.integer(quantile(w, 0.25)),
    width_q75_bp      = as.integer(quantile(w, 0.75)),
    width_max_bp      = max(w),
    mean_gc           = round(mean(seq_df$gc,    na.rm = TRUE), 4),
    mean_cpg_oe       = round(mean(seq_df$cpgoe, na.rm = TRUE), 4),
    pct_promoter      = round(af$pct_promoter,   2),
    pct_exon          = round(af$pct_exon,        2),
    pct_intron        = round(af$pct_intron,      2),
    pct_downstream    = round(af$pct_downstream,  2),
    pct_intergenic    = round(af$pct_intergenic,  2),
    n_unique_chroms   = length(unique(as.character(seqnames(gr)))),
    n_unique_genes    = length(unique(na.omit(anno_df$geneId))),
    median_tss_dist   = as.integer(median(abs(anno_df$distanceToTSS), na.rm = TRUE)),
    mean_tss_dist     = round(mean(abs(anno_df$distanceToTSS), na.rm = TRUE), 1),
    stringsAsFactors  = FALSE
  )
}, labels, gr_list, metrics_df_list, anno_df_list, seq_metrics_list,
SIMPLIFY = FALSE)

summary_all <- bind_rows(summary_list)
summary_path <- file.path(out_dir, "summary_stats.tsv")
write.table(summary_all, summary_path, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Written: %s\n", basename(summary_path)))
print(summary_all[, c("label","n_regions","total_Mb","median_width_bp",
                      "pct_promoter","pct_intron","pct_intergenic",
                      "median_tss_dist")])

# =============================================================================
# 9. Overlap matrix (pairwise Jaccard + overlap count for all label pairs)
# =============================================================================
cat("\n=== Assembling overlap_matrix.tsv ===\n")

jaccard_bp <- function(A, B) {
  inter_bp <- sum(width(GenomicRanges::intersect(A, B, ignore.strand = TRUE)))
  union_bp  <- sum(width(GenomicRanges::union(A, B, ignore.strand = TRUE)))
  if (union_bp == 0L) return(0)
  inter_bp / union_bp
}

pairs   <- expand.grid(set_A = labels, set_B = labels, stringsAsFactors = FALSE)
overlap_df <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(i) {
  A <- gr_list[[pairs$set_A[i]]]
  B <- gr_list[[pairs$set_B[i]]]
  data.frame(
    set_A               = pairs$set_A[i],
    set_B               = pairs$set_B[i],
    n_overlapping_pairs = length(findOverlaps(A, B, ignore.strand = TRUE)),
    jaccard_bp          = round(jaccard_bp(A, B), 5),
    stringsAsFactors    = FALSE
  )
}))

overlap_path <- file.path(out_dir, "overlap_matrix.tsv")
write.table(overlap_df, overlap_path, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Written: %s\n", basename(overlap_path)))
print(overlap_df)

# =============================================================================
# 10. Visualizations
# =============================================================================
cat("\n=== Generating plots ===\n")

# ---------------------------------------------------------------------------
# Colour palettes
# ---------------------------------------------------------------------------
# Two-group palette: fixed colours for nucleation / spreading sites
# (applied to plots that use the per-label COLS map)
nuclei_col    <- "#FF0066"   # nucleation sites
spread_col    <- "#8c6bb1"   # spreading sites
base_palette  <- c(nuclei_col, spread_col, "#3AA57A", "#A53A6E", "#7A3AA5")
COLS          <- setNames(base_palette[seq_along(labels)], labels)

# Helper: sample N colours from paletteer_c("ggthemes::Purple") at draw time.
# Returns a discrete fill or colour scale; n is resolved automatically from
# the number of unique fill/colour values in the plot.
purple_fill_scale <- function(...) {
  discrete_scale(
    aesthetics = "fill",
    palette    = function(n) as.character(paletteer::paletteer_c("ggthemes::Purple", n = n)),
    ...
  )
}
purple_colour_scale <- function(...) {
  discrete_scale(
    aesthetics = "colour",
    palette    = function(n) as.character(paletteer::paletteer_c("ggthemes::Purple", n = n)),
    ...
  )
}

# ---------------------------------------------------------------------------
# save_pdf: two-variant helper
#   save_pdf()         — expr must return a ggplot (or grid) object
#   save_pdf_base()    — expr uses base-R graphics drawn directly to device
# ---------------------------------------------------------------------------
save_pdf <- function(filename, expr, width = 7, height = 5) {
  path <- file.path(plots_dir, filename)
  pdf(path, width = width, height = height)
  tryCatch({
    p <- force(expr)
    if (inherits(p, c("ggplot", "gg", "HeatmapList", "Heatmap", "trellis"))) {
      print(p)
    } else if (inherits(p, "grob") || inherits(p, "gtable")) {
      grid::grid.draw(p)
    }
  }, error = function(e) {
    message("  WARNING: plot failed [", filename, "]: ", conditionMessage(e))
  })
  dev.off()
  cat(sprintf("  Saved: %s\n", filename))
}

# For base-R graphics functions that draw directly to the device
save_pdf_base <- function(filename, expr, width = 7, height = 5) {
  path <- file.path(plots_dir, filename)
  pdf(path, width = width, height = height)
  tryCatch(
    force(expr),
    error = function(e)
      message("  WARNING: plot failed [", filename, "]: ", conditionMessage(e))
  )
  dev.off()
  cat(sprintf("  Saved: %s\n", filename))
}

# --------------------------------------------------------------------------
# Plots 01: TSS and gene-body heatmaps
# --------------------------------------------------------------------------
promoter_windows <- getPromoters(TxDb = txdb,
                                 upstream   = TSS_UPSTREAM,
                                 downstream = TSS_DOWNSTREAM)

# Export tss_windows.bed (0-based BED)
tss_windows_bed_path <- file.path(plots_dir, "tss_windows.bed")
tss_df <- as.data.frame(promoter_windows)[, c("seqnames", "start", "end")]
tss_df$start <- tss_df$start - 1L   # GRanges 1-based → BED 0-based
write.table(tss_df, tss_windows_bed_path,
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
cat(sprintf("  Saved: tss_windows.bed (%d regions)\n", nrow(tss_df)))

# Export gene_body.bed (0-based BED)
cat("  Generating gene_body.bed ...\n")
gene_body_gr       <- suppressWarnings(genes(txdb))
gene_body_gr       <- keepStandardChromosomes(gene_body_gr, pruning.mode = "coarse")
gene_body_df       <- as.data.frame(gene_body_gr)[, c("seqnames", "start", "end")]
gene_body_df$start <- gene_body_df$start - 1L
gene_body_bed_path <- file.path(plots_dir, "gene_body.bed")
write.table(gene_body_df, gene_body_bed_path,
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
cat(sprintf("  Saved: gene_body.bed (%d genes)\n", nrow(gene_body_df)))

# Always compute tag matrix for plot 02 (average profile comparison)
tagMatrix_list <- lapply(labels, function(lbl) {
  cat(sprintf("  Computing tag matrix [%s] ...\n", lbl))
  getTagMatrix(gr_list[[lbl]], windows = promoter_windows)
})
names(tagMatrix_list) <- labels

# Default hue palette for bigWig samples
DEFAULT_BW_HUES <- c("#88419d", "#58135e", "#39107b", "#6a51a3",
                     "#3f007d", "#54278f", "#7a0177", "#ae017e")

if (USE_BIGWIGS) {
  utils_dir <- file.path(codedir, "utils")
  source(file.path(utils_dir, "build_matrix.R"))
  source(file.path(utils_dir, "plot_enriched_heatmaps.R"))
  bw_hues <- DEFAULT_BW_HUES[seq_along(bigwig_labels)]

  # Plot 01a: TSS EnrichedHeatmap
  cat("\n  Building TSS normalizedMatrix ...\n")
  build_matrix_main(
    regions_bed   = tss_windows_bed_path,
    bigwig_files  = bigwig_files,
    bigwig_labels = bigwig_labels,
    window_bp     = HEATMAP_WINDOW_BP,
    n_bins        = HEATMAP_N_BINS,
    out_dir       = plots_dir,
    out_prefix    = "tss",
    body_scale    = FALSE
  )
  cat("  Plotting TSS EnrichedHeatmap ...\n")
  plot_enriched_heatmaps_main(
    norm_mats_rds = file.path(plots_dir, "tss_norm_mats.rds"),
    targets_rds   = file.path(plots_dir, "tss_targets.rds"),
    window_bp     = HEATMAP_WINDOW_BP,
    out_dir       = plots_dir,
    out_prefix    = "01a_tss_enriched",
    sample_labels = bigwig_labels,
    sample_hues   = bw_hues,
    out_format    = "pdf"
  )

  # Plot 01b: gene-body EnrichedHeatmap (body-scaled, ±2 kb flanks)
  cat("\n  Building gene-body normalizedMatrix ...\n")
  build_matrix_main(
    regions_bed   = gene_body_bed_path,
    bigwig_files  = bigwig_files,
    bigwig_labels = bigwig_labels,
    window_bp     = 2000L,
    n_bins        = 50L,
    out_dir       = plots_dir,
    out_prefix    = "genebody",
    body_scale    = TRUE
  )
  cat("  Plotting gene-body EnrichedHeatmap ...\n")
  plot_enriched_heatmaps_main(
    norm_mats_rds = file.path(plots_dir, "genebody_norm_mats.rds"),
    targets_rds   = file.path(plots_dir, "genebody_targets.rds"),
    window_bp     = 2000L,
    out_dir       = plots_dir,
    out_prefix    = "01b_genebody_enriched",
    sample_labels = bigwig_labels,
    sample_hues   = bw_hues,
    out_format    = "pdf"
  )

} else {
  # Fallback: geom_raster heatmaps from ChIPseeker tag matrices (no bigWigs provided)
  hm_cols <- colorRampPalette(
    as.character(paletteer::paletteer_c("ggthemes::Gold-Purple Diverging", n = 100))
  )(256)

  subsample_tm <- function(tm, max_peaks = 3000L) {
    tm <- tm[rowSums(tm) > 0L, , drop = FALSE]
    if (nrow(tm) <= max_peaks) return(tm)
    tm[sample.int(nrow(tm), max_peaks, replace = FALSE), , drop = FALSE]
  }
  bin_tm_cols <- function(tm, max_bins = 300L) {
    if (ncol(tm) <= max_bins) return(tm)
    tm[, round(seq(1, ncol(tm), length.out = max_bins)), drop = FALSE]
  }
  tm_to_df <- function(tm_list, x_seq_fn) {
    do.call(rbind, lapply(names(tm_list), function(nm) {
      tm   <- bin_tm_cols(subsample_tm(tm_list[[nm]]))
      tm_n <- tm / (rowSums(tm) + 1e-9)
      xs   <- x_seq_fn(ncol(tm_n))
      data.frame(
        sample   = nm,
        peak     = rep(seq_len(nrow(tm_n)), ncol(tm_n)),
        position = rep(xs, each = nrow(tm_n)),
        value    = as.vector(tm_n),
        stringsAsFactors = FALSE
      )
    }))
  }
  hm_theme <- theme_classic(base_size = 13) +
    theme(axis.text.y  = element_blank(),
          axis.ticks.y = element_blank(),
          strip.text   = element_text(face = "bold"))

  # Plot 01a: TSS heatmap (geom_raster fallback)
  save_pdf("01a_tss_profile_heatmap.pdf",
           width = 5 * max(1L, length(labels)), height = 7, {
    df        <- tm_to_df(tagMatrix_list,
                          function(n) seq(-TSS_UPSTREAM, TSS_DOWNSTREAM, length.out = n))
    df$sample <- factor(df$sample, levels = labels)
    q99       <- quantile(df$value, 0.99, na.rm = TRUE)
    ggplot(df, aes(x = position, y = peak, fill = value)) +
      geom_raster(interpolate = FALSE) +
      scale_fill_gradientn(colours = hm_cols,
                           limits  = c(0, max(q99, 1e-9)),
                           oob     = scales::squish,
                           name    = "Freq") +
      scale_x_continuous(name = "Genomic region (bp from TSS)") +
      scale_y_continuous(name = "Peaks", expand = c(0, 0)) +
      facet_wrap(~sample, nrow = 1) +
      hm_theme
  })

  # Plot 01b: gene-body heatmap (geom_raster fallback)
  cat("  Computing gene-body windows ...\n")
  gene_body_windows <- getBioRegion(TxDb       = txdb,
                                    by         = "gene",
                                    type       = "body",
                                    upstream   = 2000L,
                                    downstream = 1000L)
  tagMatrix_body_list <- lapply(labels, function(lbl) {
    cat(sprintf("  Computing gene-body tag matrix [%s] ...\n", lbl))
    getTagMatrix(gr_list[[lbl]], windows = gene_body_windows, nbin = 800L)
  })
  names(tagMatrix_body_list) <- labels

  save_pdf("01b_genebody_profile_heatmap.pdf",
           width = 5 * max(1L, length(labels)), height = 7, {
    df        <- tm_to_df(tagMatrix_body_list, seq_len)
    df$sample <- factor(df$sample, levels = labels)
    n_cols    <- max(df$position)
    q99       <- quantile(df$value, 0.99, na.rm = TRUE)
    ggplot(df, aes(x = position, y = peak, fill = value)) +
      geom_raster(interpolate = FALSE) +
      scale_fill_gradientn(colours = hm_cols,
                           limits  = c(0, max(q99, 1e-9)),
                           oob     = scales::squish,
                           name    = "Freq") +
      scale_x_continuous(name   = "Gene body (5'→3') + flanks",
                         breaks = c(1L, round(n_cols / 2L), n_cols),
                         labels = c("-2 kb", "Body (scaled)", "+1 kb")) +
      scale_y_continuous(name = "Peaks", expand = c(0, 0)) +
      facet_wrap(~sample, nrow = 1) +
      hm_theme
  })
}

# --------------------------------------------------------------------------
# Plot 02: Average profile comparison with 95% CI
# Built manually instead of plotAvgProf(conf=…) because ChIPseeker v1.38.0
# does not forward the conf argument when tagMatrix is a list.
# --------------------------------------------------------------------------
save_pdf("02_avg_profile_comparison.pdf", width = 5, height = 5, {
  df_prof <- do.call(rbind, lapply(labels, function(lbl) {
    tm  <- tagMatrix_list[[lbl]]
    xs  <- seq(-TSS_UPSTREAM, TSS_DOWNSTREAM, length.out = ncol(tm))
    n   <- nrow(tm)
    mu  <- colMeans(tm)
    se  <- apply(tm, 2L, sd) / sqrt(n)
    z   <- qnorm(0.975)
    data.frame(sample = lbl, position = xs,
               mean = mu, lower = mu - z * se, upper = mu + z * se,
               stringsAsFactors = FALSE)
  }))
  df_prof$sample <- factor(df_prof$sample, levels = labels)

  ggplot(df_prof, aes(x = position, colour = sample, fill = sample)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.20, colour = NA) +
    geom_line(aes(y = mean), linewidth = 1) +
    scale_colour_manual(values = COLS) +
    scale_fill_manual(values   = COLS) +
    scale_x_continuous(name   = "Genomic region (5'→3')",
                       breaks = c(-TSS_UPSTREAM, 0L, TSS_DOWNSTREAM),
                       labels = c(sprintf("-%d kb", TSS_UPSTREAM / 1000L),
                                  "TSS",
                                  sprintf("+%d kb", TSS_DOWNSTREAM / 1000L))) +
    labs(title  = "Average peak profile around TSS",
         y      = "Peak count frequency",
         colour = NULL, fill = NULL) +
    theme_classic(base_size = 14) +
    theme(legend.position = "top",
          plot.title      = element_text(face = "bold", hjust = 0.5),
          panel.border    = element_rect(colour = "black", fill = NA, linewidth = 1))
})

# --------------------------------------------------------------------------
# Plots 03: Annotation pie (one per label)
# --------------------------------------------------------------------------
for (lbl in labels) {
  fname <- sprintf("03_annotation_pie_%s.pdf", gsub("[^[:alnum:]_]", "_", lbl))
  # Use ChIPseeker's internal getAnnoStat to find how many slices the pie
  # will have (it merges annotation strings into broad categories internally),
  # so we sample exactly that many colours from the Purple palette.
  n_slices <- tryCatch(
    nrow(ChIPseeker:::getAnnoStat(anno_list[[lbl]])),
    error = function(e) 8L
  )
  n_slices <- max(n_slices, 2L)
  pie_cols  <- as.character(
    paletteer::paletteer_c("ggthemes::Purple", n = n_slices)
  ) %>% rev()
  save_pdf_base(fname, width = 6, height = 6, {
    plotAnnoPie(anno_list[[lbl]], col = pie_cols)
    title(main = sprintf("%s — genomic annotation", lbl), cex.main = 1.1)
  })
}

# --------------------------------------------------------------------------
# Plot 04: Annotation bar comparison
# --------------------------------------------------------------------------
save_pdf("04_annotation_bar_comparison.pdf", width = 9, height = 5, {
  plotAnnoBar(anno_list) +
    purple_fill_scale() +
    labs(title = "Genomic feature distribution by site type",
         x     = "Percentage (%)",
         fill  = NULL) +
    theme_classic(base_size = 14) +
    theme(legend.position = "right",
          plot.title      = element_text(face = "bold"),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))
})

# --------------------------------------------------------------------------
# Plot 05: Distance to TSS
# --------------------------------------------------------------------------
save_pdf("05_dist_to_tss.pdf", width = 9, height = 5, {
  plotDistToTSS(anno_list,
                title = "Distance to nearest TSS") +
    purple_fill_scale() +
    labs(fill = NULL) +
    theme_classic(base_size = 14) +
    theme(legend.position = "right",
          plot.title      = element_text(face = "bold", hjust = 0.5),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))
})

# --------------------------------------------------------------------------
# Plot 06: Width distribution (violin + boxplot, log scale)
# --------------------------------------------------------------------------
save_pdf("06_width_distribution.pdf", width = max(5, 2 * length(labels) + 2), height = 6, {
  df_w <- bind_rows(lapply(labels, function(lbl) {
    data.frame(label = lbl, width_bp = width(gr_list[[lbl]]))
  }))
  df_w$label <- factor(df_w$label, levels = labels)

  ggplot(df_w, aes(x = label, y = width_bp, fill = label)) +
    geom_violin(trim = TRUE, alpha = 0.75, color = NA) +
    geom_boxplot(width = 0.10, outlier.shape = NA,
                 fill = "white", color = "grey30", linewidth = 0.5) +
    scale_y_log10(labels = label_comma(), name = "Peak width (bp, log₁₀ scale)") +
    scale_fill_manual(values = COLS) +
    labs(title = "Peak width distribution",
         x     = "Site type",
         fill  = NULL) +
    theme_bw(base_size = 13) +
    theme(legend.position = "none",
          plot.title      = element_text(face = "bold"))
})

# --------------------------------------------------------------------------
# Plot 07: Polycomb domain size spectrum (ECDF + density)
# --------------------------------------------------------------------------
save_pdf("07_polycomb_domain_spectrum.pdf", width = 10, height = 6, {
  df_w <- bind_rows(lapply(labels, function(lbl) {
    data.frame(label = lbl, width_kb = width(gr_list[[lbl]]) / 1e3)
  }))
  df_w$label <- factor(df_w$label, levels = labels)

  p_ecdf <- ggplot(df_w, aes(x = width_kb, colour = label)) +
    stat_ecdf(geom = "step", linewidth = 1.1, pad = FALSE) +
    scale_x_log10(name = "Peak width (kb, log₁₀)", labels = label_comma(suffix=" kb")) +
    scale_colour_manual(values = COLS) +
    scale_y_continuous(name = "Cumulative fraction", labels = label_percent()) +
    labs(title    = "Polycomb domain size spectrum",
         subtitle = "Cumulative distribution of H3K27me3 peak widths",
         colour   = NULL) +
    theme_bw(base_size = 13) +
    theme(legend.position = "top",
          plot.title      = element_text(face = "bold"))

  p_dens <- ggplot(df_w, aes(x = width_kb, fill = label)) +
    geom_density(alpha = 0.55, linewidth = 0.4) +
    scale_x_log10(name = "Peak width (kb)", labels = label_comma(suffix=" kb")) +
    scale_fill_manual(values = COLS) +
    labs(y = "Density", fill = NULL) +
    theme_bw(base_size = 11) +
    theme(legend.position = "none")

  grid.arrange(p_ecdf, p_dens, ncol = 2, widths = c(1.6, 1))
})

# --------------------------------------------------------------------------
# Plot 08: GC content and CpG O/E
# --------------------------------------------------------------------------
save_pdf("08_gc_cpgoe_distribution.pdf", width = 11, height = 5, {
  df_gc <- metrics_all %>% filter(!is.na(gc_content))
  if (nrow(df_gc) == 0) {
    plot.new()
    text(0.5, 0.5, "BSgenome not available.\nInstall bioconductor-bsgenome.hsapiens.ucsc.hs1",
         cex = 1.2, col = "grey40")
  } else {
    p_gc <- ggplot(df_gc, aes(x = gc_content, fill = label)) +
      geom_density(alpha = 0.5) +
      scale_fill_manual(values = COLS) +
      labs(title = "GC content per peak", x = "GC fraction", y = "Density",
           fill = NULL) +
      theme_classic(base_size = 14) +
      theme(legend.position = "top",
            plot.title = element_text(face = "bold", size = 13),
            panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))

    p_cpg <- ggplot(df_gc %>% filter(!is.na(cpg_oe)),
                    aes(x = cpg_oe, fill = label)) +
      geom_density(alpha = 0.5) +
      coord_cartesian(xlim = c(0, 1.5)) +
      scale_fill_manual(values = COLS) +
      labs(title    = "CpG O/E ratio per peak",
           subtitle = "(nCpG × length) / (nC × nG)",
           x = "CpG O/E", y = "Density",
           fill = NULL) +
      theme_classic(base_size = 14) +
      theme(legend.position = "top",
            plot.title = element_text(face = "bold", size = 13),
            panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))

    grid.arrange(p_gc, p_cpg, ncol = 2)
  }
})

# --------------------------------------------------------------------------
# Plot 09: Chromosome distribution
# --------------------------------------------------------------------------
save_pdf("09_chromosome_distribution.pdf", width = 12, height = 6, {
  df_chr <- bind_rows(lapply(labels, function(lbl) {
    data.frame(label = lbl, chr = as.character(seqnames(gr_list[[lbl]])))
  })) %>%
    filter(grepl("^chr([0-9]+|X|Y)$", chr)) %>%
    mutate(
      chr_num = sub("chr", "", chr),
      chr = factor(chr, levels = paste0("chr", c(as.character(1:22), "X", "Y"))),
      label = factor(label, levels = labels)
    ) %>%
    group_by(label, chr) %>%
    summarise(n_peaks = n(), .groups = "drop")

  ggplot(df_chr, aes(x = chr, y = n_peaks, fill = label)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.85) +
    scale_fill_manual(values = COLS) +
    scale_y_continuous(labels = label_comma()) +
    labs(title = "Peak count by chromosome",
         x = "Chromosome", y = "Number of peaks", fill = NULL) +
    theme_bw(base_size = 12) +
    theme(axis.text.x  = element_text(angle = 60, hjust = 1, size = 9),
          legend.position = "top",
          plot.title    = element_text(face = "bold"))
})

# --------------------------------------------------------------------------
# Plot 10: Distance-to-TSS cumulative distribution (ggplot version,
#          complementing ChIPseeker's plotDistToTSS bar chart)
# --------------------------------------------------------------------------
save_pdf("10_tss_dist_ecdf.pdf", width = 5.25, height = 5, {
  df_tss <- metrics_all %>%
    mutate(abs_dist_kb = abs(dist_to_tss_bp) / 1e3) %>%
    filter(!is.na(abs_dist_kb))
  df_tss$label <- factor(df_tss$label, levels = labels)

  ggplot(df_tss, aes(x = abs_dist_kb, colour = label)) +
    stat_ecdf(geom = "step", linewidth = 1.1, pad = FALSE) +
    scale_x_log10(name = "Distance to nearest TSS (kb, log₁₀)",
                  labels = label_comma(suffix = " kb"),
                  limits = c(0.01, NA)) +
    scale_y_continuous(name = "Cumulative fraction", labels = label_percent()) +
    scale_colour_manual(values = COLS) +
    labs(title    = "Distance to nearest TSS — cumulative distribution",
         subtitle = "Nucleation sites expected closer to TSS (promoter polycomb)",
         colour   = NULL) +
    theme_classic(base_size = 14) +
    theme(legend.position = "top",
          plot.title      = element_text(face = "bold", hjust = 0.5),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))
})

# --------------------------------------------------------------------------
# Plot 11: Genomic coverage per chromosome (total bp per chr, stacked)
# --------------------------------------------------------------------------
save_pdf("11_chromosome_coverage.pdf", width = 12, height = 6, {
  df_cov <- bind_rows(lapply(labels, function(lbl) {
    gr <- gr_list[[lbl]]
    data.frame(
      label    = lbl,
      chr      = as.character(seqnames(gr)),
      width_bp = width(gr)
    )
  })) %>%
    filter(grepl("^chr([0-9]+|X|Y)$", chr)) %>%
    mutate(
      chr   = factor(chr, levels = paste0("chr", c(as.character(1:22), "X", "Y"))),
      label = factor(label, levels = labels)
    ) %>%
    group_by(label, chr) %>%
    summarise(total_Mb = sum(width_bp) / 1e6, .groups = "drop")

  ggplot(df_cov, aes(x = chr, y = total_Mb, fill = label)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.85) +
    scale_fill_manual(values = COLS) +
    labs(title = "Genomic coverage by chromosome",
         x = "Chromosome", y = "Total coverage (Mb)", fill = NULL) +
    theme_bw(base_size = 12) +
    theme(axis.text.x  = element_text(angle = 60, hjust = 1, size = 9),
          legend.position = "top",
          plot.title    = element_text(face = "bold"))
})

# --------------------------------------------------------------------------
# Plots 12 & 13: Per-label chromosome coverage (highest to lowest)
# --------------------------------------------------------------------------
for (lbl in labels) {
  slug  <- gsub("[^[:alnum:]_]", "_", lbl)
  fname <- sprintf("12_chromosome_coverage_%s.pdf", slug)

  save_pdf(fname, width = 10, height = 5, {
    # Build coverage table for this single label
    gr_lbl <- gr_list[[lbl]]
    df_lbl <- data.frame(
      chr      = as.character(seqnames(gr_lbl)),
      width_bp = width(gr_lbl)
    ) %>%
      filter(grepl("^chr([0-9]+|X|Y)$", chr)) %>%
      group_by(chr) %>%
      summarise(total_Mb = sum(width_bp) / 1e6, .groups = "drop") %>%
      arrange(desc(total_Mb)) %>%
      mutate(chr = factor(chr, levels = chr))   # order by coverage

    fill_col <- COLS[[lbl]]

    ggplot(df_lbl, aes(x = chr, y = total_Mb)) +
      geom_bar(stat = "identity", fill = fill_col, alpha = 0.87, width = 0.7) +
      scale_y_continuous(labels = label_comma(suffix = " Mb"),
                         name   = "Total coverage (Mb)") +
      labs(title    = sprintf("%s \u2014 coverage by chromosome", lbl),
           subtitle = "Chromosomes ordered by total coverage (high to low)",
           x        = NULL) +
      theme_bw(base_size = 12) +
      theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 9),
            plot.title  = element_text(face = "bold"),
            panel.grid.major.x = element_blank())
  })
}

# --------------------------------------------------------------------------
# Plot 13: Comparative KEGG pathway enrichment (compareCluster + dotplot)
# Always writes a PDF — placeholder page if enrichment fails or packages absent.
# --------------------------------------------------------------------------
save_pdf("13_kegg_compareCluster.pdf", width = 10, height = 8, {
  if (!CLUSTERP_AVAILABLE) {
    plot.new()
    text(0.5, 0.5,
         paste0("Plot 13: Comparative KEGG Pathway Enrichment\n\n",
                "clusterProfiler and/or org.Hs.eg.db not installed.\n",
                "Install with:\n",
                "  conda install -c bioconda bioconductor-clusterprofiler\n",
                "  conda install -c bioconda bioconductor-org.hs.eg.db"),
         cex = 0.9, col = "grey40", adj = c(0.5, 0.5))
  } else {
    genes_by_label <- lapply(anno_list,
                             function(x) unique(na.omit(as.data.frame(x)$geneId)))
    genes_entrez <- lapply(genes_by_label, function(ids) {
      tryCatch(
        bitr(ids, fromType = "ENSEMBL", toType = "ENTREZID",
             OrgDb = org.Hs.eg.db)$ENTREZID,
        error = function(e) {
          warning("Gene ID conversion failed: ", conditionMessage(e))
          character(0)
        }
      )
    })

    ck_kegg <- tryCatch(
      compareCluster(geneCluster   = genes_entrez,
                     fun           = "enrichKEGG",
                     organism      = "hsa",
                     pvalueCutoff  = 0.05,
                     pAdjustMethod = "BH"),
      error = function(e) {
        message("  WARNING: compareCluster KEGG failed: ", conditionMessage(e))
        NULL
      }
    )

    if (!is.null(ck_kegg) && nrow(as.data.frame(ck_kegg)) > 0L) {
      dotplot(ck_kegg, showCategory = 15,
              title = "Comparative KEGG Pathway Enrichment") +
        theme_classic(base_size = 13) +
        theme(plot.title  = element_text(face = "bold", hjust = 0.5),
              axis.text.x = element_text(angle = 30, hjust = 1))
    } else {
      plot.new()
      text(0.5, 0.5,
           "Plot 13: Comparative KEGG Pathway Enrichment\n\nNo KEGG terms enriched at p < 0.05.",
           cex = 1.1, col = "grey40", adj = c(0.5, 0.5))
    }
  }
})

# =============================================================================
# 11. Done
# =============================================================================
cat(sprintf("\n=== peak_annotation.R complete ===\n"))
cat(sprintf("Output directory: %s\n", out_dir))
cat(sprintf("  Tables  : region_metrics.tsv  summary_stats.tsv  overlap_matrix.tsv\n"))
cat(sprintf("  Plots   : %s/\n\n", plots_dir))
