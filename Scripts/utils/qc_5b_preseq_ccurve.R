#!/usr/bin/env Rscript

# =============================================================================
# qc_5b_preseq_ccurve.R  –  Plot preseq c_curve library complexity curves
#
# Imports all *.preseq_ccurve.txt files produced by qc_5b (step_3_qc.py) and
# generates a complexity curve plot modelled after the MultiQC preseq panel:
#
#   - X-axis : total molecules sequenced (including duplicates)
#   - Y-axis : distinct molecules (unique library complexity)
#   - One line per sample, each with a distinct shade from the Purple palette
#   - A dashed diagonal reference line (perfect complexity: distinct = total)
#   - Axis limits and breaks derived from the actual data range (no hardcoded scale)
#
# Output:
#   {config.LIB_COMPLEXITY_DIR}/preseq/preseq_ccurve.svg
#
# Inputs:
#   {config.LIB_COMPLEXITY_DIR}/preseq/*.preseq_ccurve.txt
#   {config.LIB_COMPLEXITY_DIR}/preseq/*.preseq_ccurve.txt  (no metadata required)
# =============================================================================

suppressPackageStartupMessages({
  library(reticulate)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(stringr)
  library(paletteer)
  library(plotly)
  library(htmlwidgets)
})

# ---------------------------------------------------------------------------
# 1. Load config.py via reticulate
# ---------------------------------------------------------------------------
use_python(Sys.which("python"), required = TRUE)
py_run_file("Scripts/config.py")

PRESEQ_DIR   <- file.path(py$LIB_COMPLEXITY_DIR, "preseq")
METADATA_DIR <- py$METADATA
PROJECT      <- py$PROJECT_NAME

# ---------------------------------------------------------------------------
# 2. Discover *.preseq_ccurve.txt files
# ---------------------------------------------------------------------------
ccurve_files <- sort(list.files(
  path       = PRESEQ_DIR,
  pattern    = "\\.preseq_ccurve\\.txt$",
  full.names = TRUE
))

if (length(ccurve_files) == 0) {
  stop(
    "No *.preseq_ccurve.txt files found in: ", PRESEQ_DIR, "\n",
    "Run step_3_qc.py --only qc_5b to generate them."
  )
}
cat("Found", length(ccurve_files), "preseq c_curve file(s).\n")

# ---------------------------------------------------------------------------
# 3. Read and stack all c_curve files
# ---------------------------------------------------------------------------
# preseq c_curve output: tab-separated, two columns (no header by default,
# but some versions emit a header line "total_reads\tdistinct_reads").
# We handle both cases robustly.

read_ccurve <- function(path) {
  raw  <- readLines(path)
  raw  <- raw[nzchar(trimws(raw))]          # drop empty lines
  has_header <- grepl("^total", raw[1], ignore.case = TRUE)
  data_lines <- if (has_header) raw[-1] else raw

  vals <- do.call(rbind, strsplit(data_lines, "\\s+"))
  df   <- data.frame(
    total_reads    = as.numeric(vals[, 1]),
    distinct_reads = as.numeric(vals[, 2]),
    stringsAsFactors = FALSE
  )
  df[!is.na(df$total_reads) & !is.na(df$distinct_reads), ]
}

curve_list <- lapply(ccurve_files, function(f) {
  sn <- sub("\\.preseq_ccurve\\.txt$", "", basename(f))
  tryCatch({
    df        <- read_ccurve(f)
    df$sample <- sn
    df
  }, error = function(e) {
    warning("Could not parse: ", f, " (", e$message, ")")
    NULL
  })
})
curve_list <- Filter(Negate(is.null), curve_list)

if (length(curve_list) == 0) {
  stop("All c_curve files failed to parse. Check file format.")
}

all_curves <- bind_rows(curve_list)
cat("Loaded c_curve data for", length(unique(all_curves$sample)), "sample(s).\n")

# ---------------------------------------------------------------------------
# 4. Build a per-sample colour palette  ("ggthemes::Purple")
#    Samples are ordered alphabetically; the darkest shade is assigned to
#    the first sample, lightest to the last.
# ---------------------------------------------------------------------------
samples       <- sort(unique(all_curves$sample))
is_ctrl       <- grepl("IgG|Input", samples, ignore.case = TRUE)
samples       <- c(samples[!is_ctrl], samples[is_ctrl])

all_curves$sample <- factor(all_curves$sample, levels = samples)

n_samples_pal <- max(length(samples), 2L)
smp_pal       <- as.character(paletteer_c("ggthemes::Purple", n_samples_pal))

# Darkest hue → first sample in levels (IPs), lightest -> last (IgG/Input)
sample_colors <- setNames(rev(smp_pal), samples)

# ---------------------------------------------------------------------------
# 6. Axis-label formatter and data-driven axis limits / breaks
# ---------------------------------------------------------------------------
fmt_millions <- function(x) {
  ifelse(is.na(x), NA_character_,
         paste0(formatC(x / 1e6, format = "f", digits = 2), " M"))
}

# Compute range from the actual data — no hardcoded floor
x_max <- max(all_curves$total_reads,    na.rm = TRUE)
y_max <- max(all_curves$distinct_reads, na.rm = TRUE)

# Choose a break interval that yields roughly 5-8 breaks on each axis.
# We pick the smallest "nice" step (0.5, 1, 2, 5, 10, 25, 50 … ×10^k)
# such that ceiling(axis_max / step) <= 8.
nice_step <- function(data_max, max_breaks = 8L) {
  magnitude <- 10 ^ floor(log10(data_max))      # order-of-magnitude anchor
  candidates <- c(0.5, 1, 2, 5, 10, 25, 50) * magnitude
  steps <- candidates[ceiling(data_max / candidates) <= max_breaks]
  if (length(steps) == 0L) steps <- candidates[length(candidates)]
  min(steps)
}

x_step    <- nice_step(x_max)
y_step    <- nice_step(y_max)

x_lim_max <- ceiling(x_max / x_step) * x_step
y_lim_max <- ceiling(y_max / y_step) * y_step

x_breaks  <- seq(0, x_lim_max, by = x_step)
y_breaks  <- seq(0, y_lim_max, by = y_step)

# ---------------------------------------------------------------------------
# 8. Build the ggplot
# ---------------------------------------------------------------------------

# Reference diagonal: perfect-complexity library (distinct == total)
diag_df <- data.frame(
  x = c(0, min(x_lim_max, y_lim_max)),
  y = c(0, min(x_lim_max, y_lim_max))
)

p <- ggplot(all_curves, aes(x = total_reads, y = distinct_reads,
                             colour = sample)) +

  # — Reference diagonal (perfect library: distinct == total) —
  geom_line(data = diag_df, aes(x = x, y = y),
            inherit.aes = FALSE,
            colour = "black", linetype = "dashed",
            linewidth = 1.1, alpha = 0.5) +

  # — Complexity curves (one line per sample) —
  geom_line(linewidth = 1.1, alpha = 0.85) +

  # — Colour scale: one shade per sample —
  scale_colour_manual(values = sample_colors, name = "Sample") +

  # — Axis scales: millions labels —
  scale_x_continuous(
    breaks = x_breaks,
    labels = fmt_millions,
    limits = c(0, x_lim_max),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_y_continuous(
    breaks = y_breaks,
    labels = fmt_millions,
    limits = c(0, y_lim_max),
    expand = expansion(mult = c(0, 0.02))
  ) +

  # — Labels —
  labs(
    title    = "Preseq: Complexity curve",
    subtitle = paste0(PROJECT, " — one line per sample, coloured by sample"),
    x        = "Total molecules (including duplicates)",
    y        = "Unique molecules",
    caption  = paste0("n = ", length(unique(all_curves$sample)),
                      " samples | preseq c_curve")
  ) +

  # — Theme —
  theme_classic(base_size = 13) +
  theme(
    plot.title       = element_text(face = "bold", size = 15),
    plot.subtitle    = element_text(size = 11, colour = "grey40"),
    plot.caption     = element_text(size = 9, colour = "grey50", hjust = 1),
    axis.text.x      = element_text(angle = 30, 
                                    hjust = 1,
                                    size = 11,
                                    colour = "black"),
    axis.text.y      = element_text(size = 11,
                                    colour = "black"),
    legend.position  = "right",
    legend.title     = element_text(face = "bold", size = 11),
    legend.text      = element_text(size = 11),
    panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.35),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

# ---------------------------------------------------------------------------
# 9. Export as SVG
# ---------------------------------------------------------------------------
n_samples  <- length(unique(all_curves$sample))
html_width <- max(900, 600 + n_samples * 18)

out_html <- file.path(PRESEQ_DIR, "preseq_ccurve.html")

saveWidget(ggplotly(p, width = html_width, height = 600), out_html)

cat("\nPreseq complexity curve written to HTML:", out_html, "\n")
cat(sprintf("  Samples plotted : %d\n", n_samples))
cat("=== Done ===\n")
