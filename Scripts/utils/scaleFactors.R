### Import config module
options(repos = c(CRAN = "https://cran.rstudio.com"))

if (!require("reticulate", quietly = TRUE)) {
  install.packages("reticulate", quietly = TRUE)
}
library(reticulate, quietly = TRUE)

use_python(Sys.which("python"), required = TRUE)

py_run_file("Scripts/config.py")

#----------
# LOAD REQUIRED PACKAGES
#----------

if (!require(csaw, quietly = TRUE)) {
  BiocManager::install("csaw", quietly = TRUE)
  library(csaw, quietly = TRUE)
}
if (!require(edgeR)) {
  BiocManager::install("edgeR")
  library(edgeR, quietly = TRUE)
}

if (!require("BiocParallel", quietly = TRUE)) {
  BiocManager::install("BiocParallel", quietly = TRUE)
}

if (!require("KernSmooth", quietly = TRUE)) {
  install.packages("KernSmooth",
                   quietly = TRUE)
}

library(BiocParallel, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(readr, quietly = TRUE)
library(KernSmooth, quietly = TRUE)

# Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 3) {
  TARGETS <- unlist(strsplit(args[1], ","))
  TARGETS <- trimws(TARGETS) # remove trailing spaces
  EFFICIENCY_BIAS <- as.logical(args[2])
  BIN_WIDTH <- as.integer(args[3])
}


WORKDIR <- file.path(py$SCALINGDIR, "tmm")
BAMDIR <- py$PROCESSEDBAMDIR
METADATA <- py$METADATA
IS_PAIRED_END <- py$IS_PAIRED_END
MAPQ <- py$MAPQ

if (IS_PAIRED_END == TRUE) {
  file_path <- file.path(METADATA, "paired_samples.txt")
} else {
  file_path <- file.path(METADATA, "samples.txt")
}

if (py$REMOVE_DUPLICATES == TRUE) {
  SUFFIX <- paste0(".qc.sort.rmdup.mapq", MAPQ, ".final.bam")
  DEDUP <- FALSE
} else {
  SUFFIX <- paste0(".qc.sort.markdup.mapq", MAPQ, ".final.bam")
  DEDUP <- TRUE
}
#----------

#-----------------------------------------------------------------------------#
# CALCULATE SCALE FACTORS FOR NORMALIZATION OF GENOME COVERAGE DATA (bigWigs) #
#-----------------------------------------------------------------------------#

#----------
# PREPARE DATA
#----------

# Set up working directory
if (!dir.exists(WORKDIR)) {
  dir.create(WORKDIR, recursive = TRUE)
}
setwd(WORKDIR)

## Load Object with the maximum fragment length calculated for your data set (from getPESizes)    
rds_path <- file.path(WORKDIR, "max_frag_size.rds")
maxFrag <- readRDS(rds_path)

## Define readParam object with parameters for extracting reads from BAM files
## For consistency, here the parameters sould be the same if used in differential binding analyses


if (IS_PAIRED_END) {
  pe <- "both"
} else {
  pe <- "none"
}

param <- readParam(minq = MAPQ,
                   pe = pe,
                   max.frag = maxFrag,
                   dedup = DEDUP)

#----------
# RUN PIPELINE
#----------
for (target in TARGETS) {
  print("\n-----------------------------------")
  print(paste0("Calculating scale factors for target: ", target))
  print("-----------------------------------\n")

  # Use list.files to get filenames for each target
  matching_files <- list.files(path = BAMDIR,
                               pattern = target,
                               full.names = TRUE,
                               ignore.case = TRUE)

  matching_files <- matching_files[!grepl(".bai", matching_files)]
  # Filter to only keep files ending with the correct SUFFIX
  matching_files <- matching_files[endsWith(matching_files, SUFFIX)]

  ## Prepare vectors to store data
  normfacs <- vector()
  libSizes <- vector()

  ## Prepare vectors with paths to bam files
  print("Using the following BAM files to compute scale factors: ")
  print(matching_files)
  if (length(matching_files) == 0) {
    cat("Could not find matching files for this target, skipping")
    next
  }
  bam.files <- matching_files

  ## Count reads into windows (normalization for efficiency bias, assumes non-DB majority)
  ## or bins (normalization for composition bias, assumes systematic differences accross samples)
  ## Normalization for efficiency bias using TMM normalization from edgeR used to get
  ## normalization factors between replicates in the same conditions
  ## Normalization for composition bias (counting reads into bins) used to get
  ## normalization factors when normalizing in relation to treated samples (i.e. those samples without replicates)

  if (EFFICIENCY_BIAS) {
    BIN_WIDTH <- 150
  }

  msg <- paste0("\nEfficiency bias set to ", EFFICIENCY_BIAS)
  print(msg)

  counts <- windowCounts(bam.files,
                         bin = !EFFICIENCY_BIAS,
                         width = BIN_WIDTH,
                         param = param)

  if (!EFFICIENCY_BIAS == TRUE) {
    print("Counting reads into bins -- normalizing for composition bias")
    print("To normalize for composition bias, use large bins (>= 10kb)")
    msg <- paste0("Bin width set to: ", BIN_WIDTH, " bp")
    print(msg)
  } else {
    print("Counting reads into windows --  normalizing for efficiency bias")
    msg <- paste0("Window size set to: ", BIN_WIDTH, " bp")
    print(msg)
  }
  msg <- paste0("Using a max fragment size of ", maxFrag, " bp")
  print(msg)
  msg <- paste0("Excluding reads with a quality score lower than ", MAPQ)
  print(msg)
  msg <- paste0("Deduplication is set to ", DEDUP)
  print(msg)
  if (pe == "both") {
    print("Extracting reads in paired-end mode")
  } else {
    print("Extracting reads in single-end mode")
  }

  ## Get normalization factors from TMM normalization
  ## Use weigth = FALSE when normalizing for composition bias
  ## (avoids trimming potentially DB bins)
  normfacs <- append(normfacs,
                     normFactors(counts,
                                 se.out = FALSE,
                                 weighted = EFFICIENCY_BIAS))

  print("Getting normalization factors.")
  print("To normalize for composition bias, weighting must be set to FALSE")
  msg <- paste0("Weighting is currently set to ", EFFICIENCY_BIAS)
  print(msg)


  ## Get library sizes:
  libSizes <- append(libSizes,
                     counts$totals)
  print("Calculating size factors as:")
  print("   Normalization factors * Library Size / 1000000")
  ## Calculate size factors:
  SizeFactors <- normfacs * libSizes / 1000000

  # Extract sample names from filenames
  sample.names <- basename(matching_files)
  sample.names <- sub(SUFFIX, "", sample.names)

  ## Data table with sample names and corresponding normalization factors
  normSamples <- data.frame(ID = sample.names,
                            SF = SizeFactors)

  print("Exporting table with scale factors.")
  print("Use this table in the config file when running bamCoverage!")

  target <- gsub("|", "_", target, fixed = TRUE)
  ## Save samples file to be used in the config file
  scale_factors_txt <- file.path(WORKDIR, paste0("tmm_scaling_factors_",
                                                 target,
                                                 ".txt"))
  write_delim(normSamples,
              file = scale_factors_txt,
              delim = "\t",
              col_names = TRUE)

  print("Exporting table with sample order used in MA plot")
  ## Save data table with the order of samples as used in the MA plot
  sample_order <- file.path(WORKDIR, paste0("sample_order_",
                                            target,
                                            ".txt"))

  write_delim(normSamples,
              file = sample_order,
              delim = "\t",
              col_names = TRUE)

  #----------
  # CHECK NORMALIZATION WITH MA PLOT
  #----------

  # The log-ratio of normalization factors should pass through the center of the cloud in the plot
  # Clouds at low A-values represent background and clouds at high A-values represent bound regions
  # Normalization for composition bias: line should cross the center of the "background" cloud
  # Normalization for efficiency bias: line should cross the center of the "bound regions" cloud


  ##  Calculate abundances (log2 CPM)
  adj.counts <- cpm(asDGEList(counts), log=TRUE)

  print("Exporting MA plot")
  pdf_name <- file.path(WORKDIR, paste0(target, "_MA_plot.pdf"))

  pdf(pdf_name, width = 11, height = 7.5)
  par(mfrow = c(2, 3), mar = c(5, 4, 5, 4))
  for (i in seq_len(length(bam.files) -1)) {
    cur.x <- adj.counts[,1] # signal abundances for sample 1
    cur.y <- adj.counts[,1 + i]
    smoothScatter(x = (cur.x + cur.y) / 2,
                  y = cur.x - cur.y,
                  xlab = "Average signal (logCPM) across samples",
                  ylab = "Difference in signal between samples (log ratio)",
                  main = paste(sample.names[1],
                               "\n", "vs\n",
                               sample.names[i + 1]))

    # Ratio between normalization factors for each comparison, log2-normalized
    all.dist <- diff(log2(normfacs[c(i + 1, 1)]))
    # Add reference line marking the log2-ratio of normalization factors
    abline(h = all.dist, col = "red")
  }
  print("Closing PDF")
  dev.off()
  print("-----------------------------------\n")
}

print("\nDone.")
