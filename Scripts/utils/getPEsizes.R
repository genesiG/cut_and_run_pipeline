#!/usr/bin/env Rscript
# Script to get maximum fragment size from CUT&RUN BAM files

# LOAD REQUIRED PACKAGES
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager",
                   repos = "https://cran.rstudio.com",
                   quietly = TRUE)

if (!require(csaw)) {
  BiocManager::install("csaw", quietly = TRUE, update = FALSE)
  library(csaw, quietly = TRUE)
}
if (!require(edgeR)) {
  BiocManager::install("edgeR", quietly = TRUE, update = FALSE)
  library(edgeR, quietly = TRUE)
}
if (!require(readr)) {
  BiocManager::install("readr", quietly = TRUE, update = FALSE)
  library(readr, quietly = TRUE)
}
if (!require(dplyr)) {
  install.packages("dplyr", quietly = TRUE, update = FALSE)
  library(dplyr, quietly = TRUE)
}

### Import config module
if (!require("reticulate", quietly = TRUE)) {
  install.packages("reticulate",
                   repos = "https://cran.rstudio.com",
                   update = FALSE)
}
library(reticulate, quietly = TRUE)

use_python(Sys.which("python"), required = TRUE)
py_run_file("Scripts/config.py")

# CHANGE FOR EACH PROJECT
WORKDIR <- file.path(py$SCALINGDIR, "tmm")
BAMDIR <- py$PROCESSEDBAMDIR
METADATA <- py$METADATA
IS_PAIRED_END <- py$IS_PAIRED_END
MAPQ <- py$MAPQ
file_path <- ifelse(IS_PAIRED_END,
                    file.path(METADATA, "paired_samples.txt"),
                    file.path(METADATA, "samples.txt"))

DEDUP <- as.logical(py$REMOVE_DUPLICATES)

if (DEDUP) {
  SUFFIX <- paste0(".qc.sort.rmdup.mapq", MAPQ, ".final.bam")
} else {
  SUFFIX <- paste0(".qc.sort.markdup.mapq", MAPQ, ".final.bam")
}

# Set up working directory
if (!dir.exists(WORKDIR)) {
  dir.create(WORKDIR, recursive = TRUE)
}
setwd(WORKDIR)

rds_path <- file.path(WORKDIR, "max_frag_size.rds")
txt_path <- file.path(WORKDIR, "max_frag_size.txt")
frag_sizes_path <- file.path(WORKDIR, "frag_sizes.txt")

if (file.exists(rds_path) || file.exists(txt_path)) {
  cat("\nMax frag size file already found. Using this file for normalization\n")
  cat("\nTo recalculate fragment sizes, delete the current max frag size file\n")
} else {
  # Read samples file
  samples <- read.table(file_path,
                        sep = "\t",
                        header = TRUE)
  
  # Generate data frame to store fragment sizes
  samples.names <- vector()
  max.frag <- vector()
  sizesDF <- data.frame(sample_name = character(),
                        meanFragSizes = numeric(),
                        sdFragSizes = numeric(),
                        medianFragSizes = numeric())
  
  # Iterate through each sample to get fragment sizes
  for (i in seq_along(samples$sample.name))  {
    files.names <- samples$sample.name[i]
    cat(paste0("File name: ", files.names))
    samples.names <- append(samples.names, files.names)
  
    bam.files <- file.path(BAMDIR, paste0(files.names, SUFFIX))
  
    cat(paste0("\nGetting fragment sizes from:\n", bam.files))
    PEsize <- getPESizes(bam.files)
  
    max.frag <- append(max.frag, max(PEsize$sizes))
    fragSize <- PEsize$sizes
  
    # Create a new data frame with the current sample information
    new_row <- data.frame(sample_name = files.names,
                          meanFragSizes = mean(fragSize),
                          sdFragSizes = sd(fragSize),
                          medianFragSizes = median(fragSize))
  
    # Bind the new row to the existing data frame
    sizesDF <- rbind(sizesDF, new_row)
  
    print(summary(fragSize))
  }
  
  cat("\nGetting maximum fragment size\n")
  maxFrag <- max(max.frag)
  
  cat("\nExporting results\n")
  
  # Export results
  saveRDS(maxFrag, rds_path)
  write(maxFrag, txt_path)
  write.table(sizesDF,
              frag_sizes_path,
              sep = "\t",
              quote = FALSE,
              row.names = FALSE,
              col.names = TRUE)

}