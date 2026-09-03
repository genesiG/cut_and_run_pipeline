#!/usr/bin/env Rscript

# =============================================================================
# qc_3_alignment_stats.R  –  Parse alignment QC files and produce summary
#                             tables and a comprehensive set of barplots.
#
# Reads from Analysis_Data/qc/alignment/ (populated by qc_3a in step_3_qc_bam.py
# from the PRE-FILTER BAM generated in step_3_process_bam.py — before chrM and
# blacklist removal — so stats correctly include chrM and spike-in reads).
#
# Outputs (in Analysis_Data/qc/alignment/):
#   alignment_stats_summary.tsv       : per-sample summary table
#   alignment_stats_barplot.html      : stacked bar — raw counts (millions)
#   alignment_stats_barplot_pct.html  : stacked bar — proportions (horizontal, %)
#   alignment_stats_pct_spikein_*.html: % spike-in per sample (split by ANTIBODY_PATTERNS)
#   alignment_stats_pct_chrM.html     : % chrM per sample (horizontal)
#   alignment_stats_pct_unmapped.html : % unmapped per sample (horizontal)
#   alignment_stats_n_spikein.html    : spike-in read count (horizontal, 1000 ref line)
#   alignment_stats_bowtie2.html      : bowtie2 alignment score breakdown
#   alignment_stats_duplicates.html   : Picard optical/PCR duplicates (if metrics exist)
# =============================================================================

suppressPackageStartupMessages({
  library(reticulate)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(plotly)
  library(htmlwidgets)
})

# ---------------------------------------------------------------------------
# 1. Load config.py
# ---------------------------------------------------------------------------
use_python(Sys.which("python"), required = TRUE)
py_run_file("Scripts/config.py")

ALIGNDIR      <- file.path(py$QCDIR1, "alignment")   # canonical QC stats dir (pre-filter)
BT2_LOG_DIR   <- file.path(py$ANALYSIS_DATA, "alignment", "log")  # bowtie2 stderr logs
METADATA      <- py$METADATA
ANALYSIS_DATA <- py$ANALYSIS_DATA
PROJECT       <- py$PROJECT_NAME
SPIKE_PREFIX  <- if (!is.null(py$SPIKE_CHR_PREFIX) && nchar(py$SPIKE_CHR_PREFIX) > 0)
                   py$SPIKE_CHR_PREFIX else "ecoli_"

# ---------------------------------------------------------------------------
# 2. Output directory  →  Analysis_Data/qc/alignment/
# ---------------------------------------------------------------------------
out_dir <- ALIGNDIR
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ---------------------------------------------------------------------------
# 3. Read sample list
# ---------------------------------------------------------------------------
samples_file <- if (py$IS_PAIRED_END) "paired_samples.txt" else "samples.txt"
samples_path <- file.path(METADATA, samples_file)
if (!file.exists(samples_path)) stop(paste("Sample file not found:", samples_path))

meta         <- read.delim(samples_path, header = TRUE, sep = "\t",
                           stringsAsFactors = FALSE)
sample_names <- meta$sample.name

# ---------------------------------------------------------------------------
# 4. Helper functions
# ---------------------------------------------------------------------------

parse_flagstat <- function(path) {
  if (!file.exists(path)) {
    warning(paste("flagstat not found:", path))
    return(list(total = NA_real_, mapped = NA_real_))
  }
  lines       <- readLines(path)
  total_line  <- lines[grepl("in total", lines)][1]
  mapped_line <- lines[grepl("^[0-9]+ \\+ [0-9]+ mapped", lines)][1]
  exn <- function(l) as.numeric(strsplit(trimws(l), " ")[[1]][1])
  list(total = exn(total_line), mapped = exn(mapped_line))
}

parse_idxstats <- function(path, spike_prefix) {
  if (!file.exists(path)) {
    warning(paste("idxstats not found:", path))
    return(list(chrM = NA_real_, spike_in = NA_real_, mapped = NA_real_))
  }
  df <- read.delim(path, header = FALSE, sep = "\t",
                   col.names = c("chr", "length", "mapped", "unmapped"),
                   stringsAsFactors = FALSE)
  df <- df[df$chr != "*", ]
  list(
    chrM     = sum(df$mapped[df$chr == "chrM"],                   na.rm = TRUE),
    spike_in = sum(df$mapped[startsWith(df$chr, spike_prefix)],   na.rm = TRUE),
    mapped   = sum(df$mapped, na.rm = TRUE)
  )
}

parse_bowtie2 <- function(path) {
  ## Returns counts in READ PAIRS (bowtie2 PE mode: 1 unit = 1 pair)
  if (!file.exists(path)) return(NULL)
  lines <- readLines(path)

  exn <- function(pattern) {
    l <- lines[grepl(pattern, lines, fixed = FALSE)]
    if (length(l) == 0) return(NA_real_)
    as.numeric(trimws(strsplit(trimws(l[1]), " ")[[1]][1]))
  }

  total     <- exn("reads; of these:")
  unaligned <- exn("aligned concordantly 0 times")
  unique_n  <- exn("aligned concordantly exactly 1 time")
  multi     <- exn("aligned concordantly >1 times")

  if (is.na(total)) return(NULL)
  list(total = total, unique = unique_n, multi = multi, unaligned = unaligned)
}

parse_picard <- function(path) {
  ## Parse Picard MarkDuplicates metrics file.
  ## Returns list(read_pairs, unique, pcr_dup, optical_dup) in PAIRS, or NULL.
  if (!file.exists(path)) return(NULL)
  lines <- readLines(path)

  # Find the header line containing "READ_PAIRS_EXAMINED"
  hdr_idx <- which(grepl("READ_PAIRS_EXAMINED", lines))
  if (length(hdr_idx) == 0) return(NULL)

  header <- strsplit(lines[hdr_idx[1]], "\t")[[1]]
  data_l  <- lines[hdr_idx[1] + 1]
  if (is.na(data_l) || !nzchar(trimws(data_l))) return(NULL)

  vals <- strsplit(data_l, "\t")[[1]]
  if (length(vals) < length(header)) return(NULL)
  names(vals) <- header

  pairs     <- suppressWarnings(as.numeric(vals["READ_PAIRS_EXAMINED"]))
  dup_pairs <- suppressWarnings(as.numeric(vals["READ_PAIR_DUPLICATES"]))
  opt_pairs <- suppressWarnings(as.numeric(vals["READ_PAIR_OPTICAL_DUPLICATES"]))

  if (any(is.na(c(pairs, dup_pairs, opt_pairs)))) return(NULL)
  list(
    read_pairs  = pairs,
    unique      = pairs - dup_pairs,
    pcr_dup     = dup_pairs - opt_pairs,
    optical_dup = opt_pairs
  )
}

# ---------------------------------------------------------------------------
# 5. Build per-sample alignment stats (from flagstat + idxstat)
# ---------------------------------------------------------------------------
cat("Parsing alignment statistics for", length(sample_names), "samples...\n")

stats_list <- lapply(sample_names, function(sn) {
  flag <- parse_flagstat(file.path(ALIGNDIR, paste0(sn, ".flagstat.txt")))
  idx  <- parse_idxstats(file.path(ALIGNDIR, paste0(sn, ".idxstat.txt")), SPIKE_PREFIX)

  # Usable reads = mapped - chrM - spike-in
  mapped <- if (!is.na(flag$mapped)) flag$mapped else idx$mapped
  usable <- pmax(mapped - idx$chrM - idx$spike_in, 0)

  data.frame(
    sample       = sn,
    total_reads  = flag$total,
    mapped_reads = mapped,
    pct_mapped   = round(100 * mapped / flag$total, 2),
    pct_chrM     = round(100 * idx$chrM     / mapped, 2),
    pct_spike_in = round(100 * idx$spike_in / mapped, 2),
    usable_reads = usable,
    pct_unmapped = round(100 * (flag$total - mapped) / flag$total, 2),
    n_chrM       = idx$chrM,
    n_spike_in   = idx$spike_in,
    stringsAsFactors = FALSE
  )
})
stats_df <- bind_rows(stats_list)

# ---------------------------------------------------------------------------
# 6. Export TSV
# ---------------------------------------------------------------------------
tsv_path <- file.path(out_dir, "alignment_stats_summary.tsv")
write.table(stats_df, tsv_path, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)
cat("Summary TSV written to:", tsv_path, "\n")

# ---------------------------------------------------------------------------
# 7. Shared colour palette
# ---------------------------------------------------------------------------
cat_colors <- c(
  "Usable"   = "#88419d",
  "ChrM"     = "gray60",
  "Spike-in" = "#E88A29",
  "Unmapped" = "black"
)

# Common theme additions
theme_qc <- theme_classic(base_size = 11) +
  theme(legend.position = "right")

# Width helper
bar_width <- max(10, 0.25 * nrow(stats_df) + 4)

# ---------------------------------------------------------------------------
# 8a. Plot A — Raw counts stacked barplot (VERTICAL)
# Stacking order from bottom: Usable → ChrM → Spike-in → Unmapped
# ---------------------------------------------------------------------------
plot_df <- stats_df %>%
  mutate(
    n_chrM_p   = round(mapped_reads * pct_chrM    / 100),
    n_spike_p  = round(mapped_reads * pct_spike_in / 100),
    n_unmapped = total_reads - mapped_reads,
    n_usable   = usable_reads
  ) %>%
  select(sample, n_usable, n_chrM_p, n_spike_p, n_unmapped) %>%
  pivot_longer(cols = -sample, names_to = "category", values_to = "reads") %>%
  mutate(
    category = factor(category,
                      levels = c("n_usable", "n_chrM_p", "n_spike_p", "n_unmapped"),
                      labels = c("Usable", "ChrM", "Spike-in", "Unmapped"))
  )

p_counts <- ggplot(plot_df, aes(x = sample, y = reads / 1e6, fill = category)) +
  geom_col(width = 0.75, colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = cat_colors, name = "Read category") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    title   = paste0(PROJECT, " - Alignment Statistics"),
    x       = NULL, y = "Reads (millions)",
    caption = paste0("Spike-in prefix: '", SPIKE_PREFIX, "'")
  ) +
  theme_qc +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

html_out <- file.path(out_dir, "alignment_stats_barplot.html")
saveWidget(ggplotly(p_counts, width = bar_width * 100, height = 500), html_out)
cat("Raw-count barplot written to HTML.\n")

# ---------------------------------------------------------------------------
# 8b. Plot B — Proportions stacked barplot (HORIZONTAL, 10% breaks)
# ---------------------------------------------------------------------------
plot_df_pct <- plot_df %>%
  group_by(sample) %>%
  mutate(pct = 100 * reads / sum(reads, na.rm = TRUE)) %>%
  ungroup()

p_pct <- ggplot(plot_df_pct, aes(y = sample, x = pct, fill = category)) +
  geom_col(width = 0.75, colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = cat_colors, name = "Read category") +
  scale_x_continuous(
    breaks = seq(0, 100, 10),
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.01)),
    limits = c(0, 101)
  ) +
  labs(
    title = paste0(PROJECT, " - Alignment Statistics (Proportions)"),
    y = NULL, x = "Proportion of reads (%)",
    caption = paste0("Spike-in prefix: '", SPIKE_PREFIX, "'")
  ) +
  theme_qc +
  theme(axis.text.y = element_text(size = 7))

html_out <- file.path(out_dir, "alignment_stats_barplot_pct.html")
saveWidget(ggplotly(p_pct, width = 900, height = max(500, (0.25 * nrow(stats_df) + 2) * 100)), html_out)
cat("Proportion barplot written to HTML.\n")

# ---------------------------------------------------------------------------
# Helper: single-metric horizontal bar (reused for plots C, D, E, F)
# ---------------------------------------------------------------------------
make_hbar <- function(df_in, x_col, title_txt, x_label,
                      fill_col = "#88419d", vline = NULL, vline_label = NULL, x_limits = NULL) {
  df_in$y_ord <- reorder(df_in$sample, df_in[[x_col]])
  p <- ggplot(df_in, aes(y = y_ord, x = .data[[x_col]])) +
    geom_col(fill = fill_col, width = 0.7) +
    labs(title = paste0(PROJECT, " - ", title_txt), y = NULL, x = x_label) +
    theme_classic(base_size = 11) +
    theme(axis.text.y = element_text(size = 7))
  if (!is.null(x_limits)) {
    p <- p + scale_x_continuous(limits = x_limits)
  }
  if (!is.null(vline)) {
    p <- p + geom_vline(xintercept = vline, linetype = "dashed",
                        colour = "red", linewidth = 0.7)
    if (!is.null(vline_label))
      p <- p + annotate("text", x = vline, y = 0.5, label = vline_label,
                        angle = 90, vjust = -0.4, hjust = 0, size = 3, colour = "red")
  }
  p
}

plot_h <- max(5, 0.28 * nrow(stats_df) + 2)

# ---------------------------------------------------------------------------
# Plot C — % Spike-in per sample (split by ANTIBODY_PATTERNS)
# ---------------------------------------------------------------------------
# Remove old/stale spike-in HTML files
old_spi_files <- list.files(out_dir, pattern = "^alignment_stats_pct_spikein.*\\.html$", full.names = TRUE)
if (length(old_spi_files) > 0) file.remove(old_spi_files)

ab_patterns <- py$ANTIBODY_PATTERNS
for (pattern in names(ab_patterns)) {
  label <- ab_patterns[[pattern]]
  sub_df <- stats_df[grepl(pattern, stats_df$sample, ignore.case = TRUE), ]
  
  if (nrow(sub_df) > 0) {
    sub_h <- max(5, 0.28 * nrow(sub_df) + 2)
    safe_key <- gsub("[^A-Za-z0-9_-]", "_", pattern)
    title_txt <- if (label == pattern) {
      paste0("% Spike-in Reads (", pattern, ")")
    } else {
      paste0("% Spike-in Reads (", label, " / ", pattern, ")")
    }
    
    p_spi_pct <- make_hbar(sub_df, "pct_spike_in",
                           title_txt,
                           "Spike-in reads (% of mapped)", fill_col = "#E88A29")
    html_out <- file.path(out_dir, paste0("alignment_stats_pct_spikein_", safe_key, ".html"))
    saveWidget(ggplotly(p_spi_pct, width = 700, height = sub_h * 100), html_out)
    cat(paste0("% Spike-in barplot for '", pattern, "' written to HTML: ", html_out, "\n"))
  } else {
    cat(paste0("NOTE: No samples matched pattern '", pattern, "' — skipping spike-in plot.\n"))
  }
}

# ---------------------------------------------------------------------------
# Plot D — % ChrM per sample
# ---------------------------------------------------------------------------
p_chrM_pct <- make_hbar(stats_df, "pct_chrM",
                         "% Mitochondrial Reads per Sample",
                         "ChrM reads (% of mapped)", fill_col = "gray60")
html_out <- file.path(out_dir, "alignment_stats_pct_chrM.html")
saveWidget(ggplotly(p_chrM_pct, width = 700, height = plot_h * 100), html_out)
cat("% ChrM barplot written to HTML.\n")

# ---------------------------------------------------------------------------
# Plot E — % Unmapped per sample
# ---------------------------------------------------------------------------
p_unmap_pct <- make_hbar(stats_df, "pct_unmapped",
                          "% Unmapped Reads per Sample",
                          "Unmapped reads (% of total)", fill_col = "black", x_limits = c(0, 100))
html_out <- file.path(out_dir, "alignment_stats_pct_unmapped.html")
saveWidget(ggplotly(p_unmap_pct, width = 700, height = plot_h * 100), html_out)
cat("% Unmapped barplot written to HTML.\n")

# ---------------------------------------------------------------------------
# Plot F — N Spike-in reads (count) with 1000-read reference line
# ---------------------------------------------------------------------------
p_spi_n <- make_hbar(stats_df, "n_spike_in",
                     "Spike-in Read Counts per Sample",
                     "Spike-in reads (n)",
                     fill_col = "#E88A29",
                     vline = 1000, vline_label = "1,000 reads")
html_out <- file.path(out_dir, "alignment_stats_n_spikein.html")
saveWidget(ggplotly(p_spi_n, width = 700, height = plot_h * 100), html_out)
cat("N Spike-in barplot written to HTML.\n")

# ---------------------------------------------------------------------------
# 9. Plot G — Bowtie2 alignment score breakdown (HORIZONTAL stacked, in pairs)
# ---------------------------------------------------------------------------
bt2_list <- lapply(sample_names, function(sn) {
  log_path <- file.path(BT2_LOG_DIR, paste0(sn, ".error"))
  r <- parse_bowtie2(log_path)
  if (is.null(r)) {
    warning(paste("Bowtie2 log not found or unreadable:", log_path))
    return(NULL)
  }
  data.frame(
    sample   = sn,
    unique   = r$unique,
    multi    = r$multi,
    unaligned = r$unaligned,
    stringsAsFactors = FALSE
  )
})
bt2_df <- bind_rows(bt2_list)

if (nrow(bt2_df) > 0) {
  # Build percentage version for plotting
  bt2_pct <- bt2_df %>%
    rowwise() %>%
    mutate(
      total_     = unique + multi + unaligned,
      pct_unique = 100 * unique    / total_,
      pct_multi  = 100 * multi     / total_,
      pct_unali  = 100 * unaligned / total_
    ) %>%
    ungroup() %>%
    select(sample, pct_unique, pct_multi, pct_unali) %>%
    pivot_longer(cols = -sample, names_to = "category", values_to = "pct") %>%
    mutate(
      category = factor(category,
                        levels = c("pct_unique", "pct_multi", "pct_unali"),
                        labels = c("Uniquely mapped", "Multimapped", "Not aligned"))
    )

  bt2_colors <- c(
    "Uniquely mapped" = "#88419d",
    "Multimapped"     = "#E88A29",
    "Not aligned"     = "black"
  )

  p_bt2 <- ggplot(bt2_pct,
                  aes(y = reorder(sample, pct), x = pct, fill = category)) +
    geom_col(width = 0.75, colour = "white", linewidth = 0.2) +
    scale_fill_manual(values = bt2_colors, name = "Alignment") +
    scale_x_continuous(
      breaks = seq(0, 100, 10),
      labels = function(x) paste0(x, "%"),
      limits = c(0, 101),
      expand = expansion(mult = c(0, 0.01))
    ) +
    labs(
      title = paste0(PROJECT, " - Bowtie2 Alignment Statistics"),
      y = NULL, x = "Proportion of read pairs (%)"
    ) +
    theme_qc +
    theme(axis.text.y = element_text(size = 7))

  html_out <- file.path(out_dir, "alignment_stats_bowtie2.html")
  saveWidget(ggplotly(p_bt2, width = 900, height = plot_h * 100), html_out)
  cat("Bowtie2 barplot written to HTML.\n")
} else {
  cat("WARNING: No bowtie2 logs parsed — skipping bowtie2 barplot.\n")
}

# ---------------------------------------------------------------------------
# 10. Plot H — Picard duplicate metrics (HORIZONTAL stacked, %)
#     Shows: unique / PCR (non-optical) duplicates / optical duplicates.
#     Skipped gracefully if Picard metrics files are absent (run qc_3a first).
# ---------------------------------------------------------------------------
picard_list <- lapply(sample_names, function(sn) {
  p_path <- file.path(ALIGNDIR, paste0(sn, ".picard_metrics.txt"))
  r <- parse_picard(p_path)
  if (is.null(r)) {
    warning(paste("Picard metrics not found or unreadable:", p_path))
    return(NULL)
  }
  data.frame(
    sample      = sn,
    unique      = r$unique,
    pcr_dup     = r$pcr_dup,
    optical_dup = r$optical_dup,
    stringsAsFactors = FALSE
  )
})
picard_df <- bind_rows(picard_list)

if (nrow(picard_df) > 0) {
  dup_pct <- picard_df %>%
    rowwise() %>%
    mutate(
      total_      = unique + pcr_dup + optical_dup,
      pct_unique  = 100 * unique      / total_,
      pct_pcr     = 100 * pcr_dup     / total_,
      pct_optical = 100 * optical_dup / total_
    ) %>%
    ungroup() %>%
    select(sample, pct_unique, pct_pcr, pct_optical) %>%
    pivot_longer(cols = -sample, names_to = "category", values_to = "pct") %>%
    mutate(
      category = factor(category,
                        levels = c("pct_unique", "pct_pcr", "pct_optical"),
                        labels = c("Unique pairs", "PCR duplicates", "Optical duplicates"))
    )

  dup_colors <- c(
    "Unique pairs"       = "#88419d",
    "PCR duplicates"     = "#E88A29",
    "Optical duplicates" = "#D65C5C"
  )

  p_dup <- ggplot(dup_pct,
                  aes(y = reorder(sample, pct), x = pct, fill = category)) +
    geom_col(width = 0.75, colour = "white", linewidth = 0.2) +
    scale_fill_manual(values = dup_colors, name = "Category") +
    scale_x_continuous(
      breaks = seq(0, 100, 10),
      labels = function(x) paste0(x, "%"),
      limits = c(0, 101),
      expand = expansion(mult = c(0, 0.01))
    ) +
    labs(
      title   = paste0(PROJECT, " - Duplicate Statistics (Picard MarkDuplicates)"),
      y       = NULL, x = "Proportion of examined pairs (%)",
      caption = "Optical duplicate pixel distance: 2500 (patterned flow cell)"
    ) +
    theme_qc +
    theme(axis.text.y = element_text(size = 7))

  html_out <- file.path(out_dir, "alignment_stats_duplicates.html")
  saveWidget(ggplotly(p_dup, width = 900, height = plot_h * 100), html_out)
  cat("Picard duplicate barplot written to HTML.\n")
} else {
  cat("NOTE: No Picard metrics found — skipping duplicate barplot.\n")
  cat("      Run: python3 Scripts/step_3_qc_bam.py --only qc_3a\n")
  cat("      Then wait for jobs to finish (bjobs -u $USER) before re-running qc_3.\n")
}

cat("\n=== Alignment stats summary ===\n")
print(stats_df[, c("sample","total_reads","mapped_reads","pct_mapped",
                    "pct_chrM","pct_spike_in","usable_reads")],
      row.names = FALSE)
cat("\n=== Done ===\n")
