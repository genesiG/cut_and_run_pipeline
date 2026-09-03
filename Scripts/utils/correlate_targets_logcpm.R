#!/usr/bin/env Rscript

# =============================================================================
# correlate_targets_logcpm.R
#
# Quantifies the correlation between two antibody targets (e.g. H3K27me2 and
# CBX2) within a single condition (e.g. 1 uM EZH2i-treated) by:
#   1. Reading pre-specified regions (e.g. lost_k27me3.bed) as GRanges
#   2. Counting reads in those regions using csaw::regionCounts over the
#      selected group's BAM files for each target
#   3. Computing per-sample logCPM via edgeR::cpm(y, log = TRUE) with
#      standard library-size normalization ONLY (no spike-in SFs).
#
#      Rationale: spike-in normalization factors (SF) are designed to correct
#      for global abundance changes between conditions (DMSO vs 1uM) within
#      the same antibody target. Here we are comparing two *different* antibody
#      targets within a *single* condition, so no cross-condition correction is
#      needed or appropriate. Standard CPM (correcting only for sequencing depth
#      between replicates) is the correct normalization.
#
#   4. Averaging replicates and producing a ggplot2 hexbin scatter plot
#      with Spearman + Pearson correlation annotations, using my_theme style
#      (Helvetica / Arial via showtext)
#
# CLI arguments
# -------------
#   --regions          PATH   BED file defining the regions to count over
#   --target_a         STR    First antibody target label  (e.g. K27me2)
#   --target_b         STR    Second antibody target label (e.g. CBX2)
#   --group            STR    Group to select from metadata (default: 1uM)
#   --metadata_a       PATH   Metadata for target A (auto-discovered)
#   --metadata_b       PATH   Metadata for target B (auto-discovered)
#   --antibody_a       STR    Antibody value in metadata ANTIBODY column for A
#                             (default: same as target_a)
#   --antibody_b       STR    Antibody value in metadata ANTIBODY column for B
#                             (default: same as target_b)
#   --label_a          STR    Display label for axis / title (default: target_a)
#   --label_b          STR    Display label for axis / title (default: target_b)
#   --out_dir          PATH   Output directory (default: Analysis_Data/correlation)
#   --out_prefix       STR    Output file prefix (auto-derived if omitted)
# =============================================================================

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(csaw)
  library(edgeR)
  library(ggplot2)
  library(showtext)
  library(reticulate)
})

# ---------------------------------------------------------------------------
# Font setup (Helvetica → Arial → sans-serif fallback via showtext)
# ---------------------------------------------------------------------------
find_font_path <- function(pattern) {
  res <- tryCatch(
    system(sprintf("fc-list : file family | grep -i '%s' | grep '\\.ttf' | head -1 | awk -F: '{print $1}' | xargs", pattern), intern = TRUE),
    error = function(e) ""
  )
  if (length(res) > 0 && nzchar(res[1]) && file.exists(res[1])) {
    return(res[1])
  }
  return("")
}

helv_path  <- find_font_path("Helvetica")
arial_path <- find_font_path("Arial")
fallback_path <- find_font_path("DejaVuSans.ttf")
if (!nzchar(fallback_path) && file.exists("$HOME/miniconda3/envs/chipseq/fonts/DejaVuSans.ttf")) {
  fallback_path <- "$HOME/miniconda3/envs/chipseq/fonts/DejaVuSans.ttf"
}
if (!nzchar(fallback_path) && file.exists("/usr/share/fonts/dejavu-sans-fonts/DejaVuSans.ttf")) {
  fallback_path <- "/usr/share/fonts/dejavu-sans-fonts/DejaVuSans.ttf"
}

if (nzchar(helv_path)) {
  font_add("Helvetica", regular = helv_path)
  PLOT_FONT <- "Helvetica"
} else if (nzchar(arial_path)) {
  font_add("Arial", regular = arial_path)
  PLOT_FONT <- "Arial"
} else if (nzchar(fallback_path)) {
  font_add("Helvetica", regular = fallback_path)
  font_add("Arial", regular = fallback_path)
  PLOT_FONT <- "Helvetica"
} else {
  PLOT_FONT <- "sans"
}
showtext_auto()
cat(sprintf("Font selected: %s (path: %s)\n", PLOT_FONT, if (nzchar(helv_path)) helv_path else if (nzchar(arial_path)) arial_path else fallback_path))

# ---------------------------------------------------------------------------
# CLI parsing
# ---------------------------------------------------------------------------
args_cli <- commandArgs(trailingOnly = TRUE)

get_flag <- function(flag, args, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0 || idx[1] + 1 > length(args)) return(default)
  args[idx[1] + 1]
}

regions_path <- get_flag("--regions",    args_cli,
  file.path("Analysis_Data", "k27me3_classification", "lost_k27me3.bed"))
target_a     <- get_flag("--target_a",   args_cli, "K27me2")
target_b     <- get_flag("--target_b",   args_cli, "CBX2")
group_use    <- get_flag("--group",      args_cli, "1uM")
antibody_a   <- get_flag("--antibody_a", args_cli, target_a)
antibody_b   <- get_flag("--antibody_b", args_cli, target_b)
label_a      <- get_flag("--label_a",    args_cli, target_a)
label_b      <- get_flag("--label_b",    args_cli, target_b)
out_dir      <- get_flag("--out_dir",    args_cli,
  file.path("Analysis_Data", "correlation"))

out_prefix_default <- paste0(
  tolower(gsub("[^A-Za-z0-9]", "_", label_a)), "_vs_",
  tolower(gsub("[^A-Za-z0-9]", "_", label_b)), "_",
  tolower(gsub("[^A-Za-z0-9]", "_", group_use)),
  "_over_",
  tolower(gsub("[^A-Za-z0-9]", "_",
    tools::file_path_sans_ext(basename(regions_path))))
)
out_prefix <- get_flag("--out_prefix", args_cli, out_prefix_default)

metadata_a <- get_flag("--metadata_a", args_cli,
  file.path("Metadata",
            paste0("sample_metadata_", target_a, "_processed.txt")))
metadata_b <- get_flag("--metadata_b", args_cli,
  file.path("Metadata",
            paste0("sample_metadata_", target_b, "_processed.txt")))

cat("=== correlate_targets_logcpm.R ===\n")
cat(sprintf("  Target A      : %s (%s)\n", target_a, label_a))
cat(sprintf("  Target B      : %s (%s)\n", target_b, label_b))
cat(sprintf("  Group         : %s\n",   group_use))
cat(sprintf("  Normalization : library-size (CPM) only — no spike-in SFs\n"))
cat(sprintf("  Regions       : %s\n",   regions_path))
cat(sprintf("  Metadata A    : %s\n",   metadata_a))
cat(sprintf("  Metadata B    : %s\n",   metadata_b))
cat(sprintf("  Output dir    : %s\n",   out_dir))
cat(sprintf("  Output prefix : %s\n",   out_prefix))

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Resolve BAM directory via reticulate → config.py
# ---------------------------------------------------------------------------
py_bin <- Sys.getenv("RETICULATE_PYTHON", unset = NA)
if (is.na(py_bin) || !nzchar(py_bin)) {
  py_bin <- Sys.which("python3")
}
if (is.na(Sys.getenv("RETICULATE_PYTHON", unset = NA)) ||
    !nzchar(Sys.getenv("RETICULATE_PYTHON", unset = ""))) {
  use_python(py_bin, required = TRUE)
}
py_run_file("Scripts/config.py")

bam_dir <- py$PROCESSEDBAMDIR
is_pe   <- isTRUE(py$IS_PAIRED_END)
cat(sprintf("\nBAM directory : %s\nPaired-end    : %s\n", bam_dir, is_pe))

# ---------------------------------------------------------------------------
# Helper: load metadata for one target, filter to antibody + group
# ---------------------------------------------------------------------------
load_bams_for_target <- function(metadata_path, antibody_col,
                                  group_sel, bam_dir) {
  if (!file.exists(metadata_path))
    stop("Metadata not found: ", metadata_path)

  meta <- read.table(metadata_path, header = TRUE, sep = "\t",
                     stringsAsFactors = FALSE)
  # Keep rows matching antibody and group
  idx  <- meta$ANTIBODY == antibody_col & meta$GROUP == group_sel
  if (sum(idx) == 0)
    stop(sprintf("No samples in metadata matching ANTIBODY='%s' GROUP='%s'",
                 antibody_col, group_sel))
  meta <- meta[idx, ]

  bam_paths <- file.path(bam_dir, meta$ID)
  missing   <- bam_paths[!file.exists(bam_paths)]
  if (length(missing) > 0)
    stop("Missing BAM files:\n", paste(missing, collapse = "\n"))

  cat(sprintf("  %d BAM(s) selected (ANTIBODY=%s, GROUP=%s):\n",
              nrow(meta), antibody_col, group_sel))
  for (i in seq_len(nrow(meta)))
    cat(sprintf("    %s\n", meta$ID[i]))

  list(paths = bam_paths, ids = meta$ID)
}

# ---------------------------------------------------------------------------
# Load regions
# ---------------------------------------------------------------------------
cat("\nLoading regions from:", regions_path, "\n")
if (!file.exists(regions_path))
  stop("Regions file not found: ", regions_path)
regions_gr <- rtracklayer::import(regions_path, format = "BED")
cat(sprintf("Loaded %d regions.\n", length(regions_gr)))

# ---------------------------------------------------------------------------
# Resolve BAMs and SFs for each target
# ---------------------------------------------------------------------------
cat(sprintf("\n--- Target A: %s ---\n", target_a))
info_a <- load_bams_for_target(metadata_a, antibody_a, group_use, bam_dir)
cat(sprintf("\n--- Target B: %s ---\n", target_b))
info_b <- load_bams_for_target(metadata_b, antibody_b, group_use, bam_dir)

# ---------------------------------------------------------------------------
# Count reads with csaw::regionCounts
# ---------------------------------------------------------------------------
pe_param <- readParam(pe = if (is_pe) "both" else "none")

cat(sprintf("\nCounting reads for %s ...\n", target_a))
counts_a <- regionCounts(info_a$paths, regions_gr, param = pe_param)

cat(sprintf("Counting reads for %s ...\n", target_b))
counts_b <- regionCounts(info_b$paths, regions_gr, param = pe_param)

# ---------------------------------------------------------------------------
# logCPM with standard library-size normalization (no spike-in SFs)
#
# Within a single condition we are comparing two *different* antibody targets,
# not the same target across conditions. Spike-in SFs are condition-contrast
# correction factors (DMSO vs 1uM) and would distort the within-condition
# comparison by applying unequal scaling to the two targets. Standard CPM
# (library-size normalization only) is the correct approach here.
# ---------------------------------------------------------------------------
compute_logcpm <- function(counts_se, sample_ids) {
  y            <- asDGEList(counts_se)
  lcpm_mat     <- cpm(y, log = TRUE)   # uses default lib.size from counts
  colnames(lcpm_mat) <- sample_ids
  lcpm_mat
}

cat(sprintf("\nComputing logCPM for %s (library-size norm) ...\n", target_a))
lcpm_a <- compute_logcpm(counts_a, info_a$ids)

cat(sprintf("Computing logCPM for %s (library-size norm) ...\n", target_b))
lcpm_b <- compute_logcpm(counts_b, info_b$ids)

# Average replicates
mean_lcpm_a <- rowMeans(lcpm_a)
mean_lcpm_b <- rowMeans(lcpm_b)

cat(sprintf("\nMean logCPM %s: min=%.2f, max=%.2f, median=%.2f\n",
            target_a, min(mean_lcpm_a), max(mean_lcpm_a),
            median(mean_lcpm_a)))
cat(sprintf("Mean logCPM %s: min=%.2f, max=%.2f, median=%.2f\n",
            target_b, min(mean_lcpm_b), max(mean_lcpm_b),
            median(mean_lcpm_b)))

# ---------------------------------------------------------------------------
# Save per-region data table
# ---------------------------------------------------------------------------
df_out <- data.frame(
  chr    = as.character(seqnames(regions_gr)),
  start  = start(regions_gr),
  end    = end(regions_gr),
  logCPM_A = mean_lcpm_a,
  logCPM_B = mean_lcpm_b
)
colnames(df_out)[4:5] <- c(paste0("logCPM_", target_a),
                            paste0("logCPM_", target_b))
csv_path <- file.path(out_dir, paste0(out_prefix, "_logcpm.csv"))
write.csv(df_out, csv_path, row.names = FALSE)
cat("Saved logCPM table to:", csv_path, "\n")

# ---------------------------------------------------------------------------
# Correlation statistics
# ---------------------------------------------------------------------------
rho   <- cor(mean_lcpm_a, mean_lcpm_b, method = "spearman")
r_val <- cor(mean_lcpm_a, mean_lcpm_b, method = "pearson")
n_reg <- length(mean_lcpm_a)

cat(sprintf("\nSpearman rho = %.3f\n", rho))
cat(sprintf("Pearson  r   = %.3f\n", r_val))
cat(sprintf("N regions    = %d\n",   n_reg))

# ---------------------------------------------------------------------------
# my_theme (mirrors inhibitor_resistance_analysis.R)
# ---------------------------------------------------------------------------
my_theme <- theme_bw(base_size = 13, base_family = PLOT_FONT) +
  theme(
    panel.grid   = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
    strip.background = element_blank(),
    strip.text   = element_text(size = 17),
    plot.title   = element_text(face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 11),
    axis.text    = element_text(color = "black", size = 14),
    axis.title   = element_text(color = "black", size = 15),
    axis.line    = element_blank()
  )

# ---------------------------------------------------------------------------
# Scatterplot: hexbin density + Spearman / Pearson annotations
# ---------------------------------------------------------------------------
df_plot <- data.frame(x = mean_lcpm_a, y = mean_lcpm_b)

corr_label <- sprintf(
  "Spearman \u03c1 = %.3f\nPearson r = %.3f\nN = %s",
  rho, r_val, format(n_reg, big.mark = ",")
)

x_rng <- range(df_plot$x, finite = TRUE)
y_rng <- range(df_plot$y, finite = TRUE)
ann_x <- x_rng[1] + 0.05 * diff(x_rng)
ann_y <- y_rng[2] - 0.05 * diff(y_rng)

p <- ggplot(df_plot, aes(x = x, y = y)) +
  geom_hex(bins = 60) +
  scale_fill_viridis_c(option = "plasma", name = "# regions") +
  geom_smooth(method = "lm", se = FALSE, color = "#39107b",
              linewidth = 0.8, linetype = "dashed") +
  annotate("text",
    x    = ann_x,
    y    = ann_y,
    label = corr_label,
    hjust = 0,
    vjust = 1,
    size  = 4,
    family = PLOT_FONT,
    color  = "black"
  ) +
  labs(
    x     = paste0(label_a, " logCPM (", group_use, ")"),
    y     = paste0(label_b, " logCPM (", group_use, ")"),
    title = paste0(label_a, " vs ", label_b,
                   " abundance\nover ", basename(regions_path))
  ) +
  my_theme

svg_path <- file.path(out_dir, paste0(out_prefix, "_scatter.svg"))
svg(svg_path, width = 5.5, height = 4.5)
print(p)
invisible(dev.off())
cat("Saved scatter plot to:", svg_path, "\n")

cat(sprintf("\n=== correlate_targets_logcpm.R complete ===\n"))
