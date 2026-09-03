#!/usr/bin/env Rscript

# =============================================================================
# qc_4_fragment_sizes.R  –  Plot fragment size distributions from PE BAM files
#
# Uses csaw::getPESizes() on processed BAM files to extract insert-size
# distributions, then produces three outputs:
#
#   (a) Density line plot overlaying all samples, coloured by antibody target
#       with vertical dashed lines at 120, 200, 400, 600 bp
#   (b) Box plot of median fragment sizes grouped by antibody target
#   (c) TSV table: mean, SD, median fragment size per sample
#
# Color scheme: paletteer_c("grDevices::Purple-Blue") for antibody targets
#
# Inputs:
#   - Processed BAM files in config.PROCESSEDBAMDIR
#   - Sample metadata from config.METADATA (samples.txt or paired_samples.txt)
# =============================================================================

suppressPackageStartupMessages({
  library(reticulate)
  library(csaw)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(paletteer)
  library(stringr)
  library(plotly)
  library(htmlwidgets)
})

# ---------------------------------------------------------------------------
# 1. Load config.py
# ---------------------------------------------------------------------------
use_python(Sys.which("python"), required = TRUE)
py_run_file("Scripts/config.py")

PROCESSEDBAMDIR <- py$PROCESSEDBAMDIR
METADATA        <- py$METADATA
ANALYSIS_DATA   <- py$ANALYSIS_DATA
MAPQ            <- as.integer(py$MAPQ)
IS_PE           <- isTRUE(py$IS_PAIRED_END)
REMOVE_DUPS     <- isTRUE(py$REMOVE_DUPLICATES)

# BAM suffix convention (matches step_3_process_bam.py)
SPECIES <- py$SPECIES
BAM_SUFFIX <- if (REMOVE_DUPS) {
  paste0(".qc.sort.rmdup.mapq", MAPQ, ".final.bam")
} else {
  paste0(".qc.sort.markdup.mapq", MAPQ, ".final.bam")
}

# ---------------------------------------------------------------------------
# 2. Output directory
# ---------------------------------------------------------------------------
out_dir <- file.path(ANALYSIS_DATA, "qc", "fragment_sizes")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ---------------------------------------------------------------------------
# 3. Read sample metadata
# ---------------------------------------------------------------------------
samples_file <- if (IS_PE) "paired_samples.txt" else "samples.txt"
samples_path <- file.path(METADATA, samples_file)

if (!file.exists(samples_path)) stop(paste("Sample file not found:", samples_path))

meta <- read.delim(samples_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
sample_names <- meta$sample.name

# Load detailed metadata (ID, ANTIBODY, GROUP columns)
meta_detail_path <- file.path(METADATA, "sample_metadata.txt")
meta_detail <- if (file.exists(meta_detail_path)) {
  read.delim(meta_detail_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
} else {
  # Fallback: derive antibody from sample name using config ANTIBODY_PATTERNS
  data.frame(
    ID      = sample_names,
    ANTIBODY = NA_character_,
    GROUP   = NA_character_,
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# 4. Collect fragment sizes per sample
# ---------------------------------------------------------------------------
cat("Collecting fragment sizes from", length(sample_names), "samples...\n")

sizes_list  <- list()
summary_rows <- list()

for (sn in sample_names) {
  bam_path <- file.path(PROCESSEDBAMDIR, paste0(sn, BAM_SUFFIX))

  if (!file.exists(bam_path)) {
    warning(paste("BAM not found, skipping:", bam_path))
    next
  }

  cat("  Processing:", sn, "\n")
  tryCatch({
    pe_sizes <- getPESizes(bam_path)
    sz <- pe_sizes$sizes
    sz <- sz[sz > 0 & sz < 2000]  # discard implausible fragments

    sizes_list[[sn]] <- data.frame(
      sample = sn,
      size   = sz,
      stringsAsFactors = FALSE
    )

    summary_rows[[sn]] <- data.frame(
      sample     = sn,
      mean_size  = round(mean(sz, na.rm = TRUE), 2),
      sd_size    = round(sd(sz,   na.rm = TRUE), 2),
      median_size = round(median(sz, na.rm = TRUE), 2),
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    warning(paste("Error processing", sn, ":", e$message))
  })
}

if (length(sizes_list) == 0) {
  stop("No fragment size data could be collected. Check BAM files and suffix.")
}

all_sizes   <- bind_rows(sizes_list)
summary_df  <- bind_rows(summary_rows)

# Annotate with antibody target from detailed metadata
all_sizes <- all_sizes %>%
  left_join(meta_detail %>% select(ID, ANTIBODY, GROUP), by = c("sample" = "ID"))

summary_df <- summary_df %>%
  left_join(meta_detail %>% select(ID, ANTIBODY, GROUP), by = c("sample" = "ID"))

# Reorder ANTIBODY factors so IgG/Input is last
ab_levels <- sort(unique(all_sizes$ANTIBODY))
is_ctrl <- grepl("IgG|Input", ab_levels, ignore.case = TRUE)
ab_levels <- c(ab_levels[!is_ctrl], ab_levels[is_ctrl])

all_sizes <- all_sizes %>%
  mutate(ANTIBODY = factor(ANTIBODY, levels = ab_levels))
summary_df <- summary_df %>%
  mutate(ANTIBODY = factor(ANTIBODY, levels = ab_levels))

# ---------------------------------------------------------------------------
# 5. Export TSV table
# ---------------------------------------------------------------------------
tsv_path <- file.path(out_dir, "fragment_sizes_summary.tsv")
write.table(summary_df, tsv_path,
            sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)
cat("Fragment size summary written to:", tsv_path, "\n")

# ---------------------------------------------------------------------------
# 6. Build colour palette — "ggthemes::Purple" per antibody target
# ---------------------------------------------------------------------------
antibodies    <- ab_levels
n_ab          <- max(length(antibodies), 2)
ab_colors     <- setNames(
  as.character(paletteer_c("ggthemes::Purple", n_ab)),
  antibodies
)

# ---------------------------------------------------------------------------
# 7. (a) Density line plot
# ---------------------------------------------------------------------------
# Aggregate density per sample so the plot is not too heavy
density_df <- all_sizes %>%
  group_by(sample, ANTIBODY) %>%
  summarise(
    density_data = list(density(size, n = 512, from = 0, to = 2000)),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    x = list(density_data$x),
    y = list(density_data$y)
  ) %>%
  select(-density_data) %>%
  unnest(cols = c(x, y))

vlines <- c(120, 200, 400, 600)

p_density <- ggplot(density_df,
                    aes(x = x, y = y, colour = ANTIBODY, group = sample)) +
  geom_line(alpha = 0.65, linewidth = 0.55) +
  geom_vline(xintercept = vlines, linetype = "dashed",
             colour = "grey40", linewidth = 0.4) +
  annotate("text", x = vlines + 15, y = Inf,
           label = paste0(vlines, " bp"),
           angle = 90, vjust = 1.5, hjust = 1.1,
           size = 2.8, colour = "grey30") +
  scale_colour_manual(values = ab_colors, name = "Antibody") +
  scale_x_continuous(limits = c(0, 1000)) +
  labs(
    title    = paste0(py$PROJECT_NAME, " — Fragment Size Distributions"),
    subtitle = "One line per sample, coloured by antibody target",
    x        = "Fragment size (bp)",
    y        = "Density"
  ) +
  theme_classic(base_size = 11) +
  theme(legend.position = "right")

density_html <- file.path(out_dir, "fragment_size_density.html")
saveWidget(ggplotly(p_density, width = 900, height = 500), density_html)
cat("Density plot written to HTML:", density_html, "\n")

# ---------------------------------------------------------------------------
# 8. (b) Box plot — median fragment sizes per sample grouped by antibody
# ---------------------------------------------------------------------------
p_box <- ggplot(summary_df,
                aes(x = ANTIBODY, y = median_size,
                    fill = ANTIBODY, colour = ANTIBODY)) +
  geom_boxplot(alpha = 0.4, outlier.shape = NA, linewidth = 0.5) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  scale_fill_manual(values   = ab_colors, name = "Antibody") +
  scale_colour_manual(values = ab_colors, name = "Antibody") +
  labs(
    title = paste0(py$PROJECT_NAME, " — Median Fragment Size by Antibody Target"),
    x     = "Antibody target",
    y     = "Median fragment size (bp)"
  ) +
  theme_classic(base_size = 11) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 35, hjust = 1))

box_html <- file.path(out_dir, "fragment_size_boxplot.html")
saveWidget(ggplotly(p_box, width = 700, height = 500), box_html)
cat("Box plot written to HTML:", box_html, "\n")

cat("\n=== Fragment size summary ===\n")
print(summary_df, row.names = FALSE)
cat("\n=== Done ===\n")
