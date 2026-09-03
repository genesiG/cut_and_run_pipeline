#!/usr/bin/env Rscript
#
# SpikeInFreeInstall.R
# --------------------
# Ensures the ChIPseqSpikeInFree GitHub package is installed into the active
# conda / R library before any parallel normalization jobs are submitted.
#
# Called ONCE by step_4a_chipseqspikeinfree.py in a blocking bsub job.
# Subsequent ChIPseqSpikeInFree.R jobs will find the package already present
# and skip installation entirely.
#

# ── devtools ────────────────────────────────────────────────────────────────
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools", repos = "https://cran.rstudio.com",
                   quiet = TRUE)
}

# ── ChIPseqSpikeInFree ──────────────────────────────────────────────────────
if (!requireNamespace("ChIPseqSpikeInFree", quietly = TRUE)) {
  message("ChIPseqSpikeInFree not found – installing from GitHub ...")
  pak::pak("stjude/ChIPseqSpikeInFree")
} else {
  message("ChIPseqSpikeInFree already installed.")
}

pkg_ver <- packageVersion("ChIPseqSpikeInFree")
cat(paste0("ChIPseqSpikeInFree version: ", pkg_ver, "\n"))

