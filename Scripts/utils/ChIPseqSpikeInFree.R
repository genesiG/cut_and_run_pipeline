#!/usr/bin/env Rscript

# Read command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Set parameters from arguments
bin_size <- as.integer(args[1])
samples_file <- args[2]
suffix <- args[3]
cutoff <- as.numeric(args[4])
max_turns <- as.numeric(args[5])
targets_arg <- args[6]
output_arg <- args[7]
paired_arg <- args[8]
run_name <- args[9]

# Validate parameters
if (length(args) != 9) {
  stop("Requires 9 parameters: \n  bin_size, samples_file, suffix, cutoff, max_turns, targets, output_dir, paired_end, run_name")
}

# Parse TARGETS (split by comma, then by pipe for groups)
# Trim whitespace from each target
targets_list <- trimws(unlist(strsplit(targets_arg, ",")))

cat("Running with parameters:\n")
cat(paste0(
  "Bin size: ", bin_size, "\n",
  "Samples file: ", samples_file, "\n",
  "QC cutoff: ", cutoff, "\n",
  "Max turns: ", max_turns, "\n",
  "Targets: ", targets_arg, "\n",
  "Output directory: ", output_arg, "\n"
))

### Prepare workspace

# Load packages
suppressMessages(library(dplyr, quietly = TRUE))
suppressMessages(library(edgeR, quietly = TRUE))
suppressMessages(library(csaw, quietly = TRUE))
suppressMessages(library(ChIPseqSpikeInFree, quietly = TRUE))
###

library(reticulate, quietly = TRUE, verbose = FALSE)
use_python(Sys.which("python"), required = TRUE)

py_run_file("Scripts/config.py")

###

### Set up file paths for pipeline parameters
# Specify path to metadata file
samples_path <- file.path(py$METADATA, samples_file)

# Read sample metadata
sample_metadata <- read.delim(samples_path, header = TRUE, sep = "\t")

# Assign path to chromosome sizes
chrom_sizes <- file.path(py$GNM_SIZES)

###

### Process each target separately
for (target in targets_list) {
  cat("\n========================================\n")
  cat(paste0("Processing target: ", target, "\n"))
  cat("========================================\n")

  # Pattern matching: split by pipe to get alternative patterns
  target_patterns <- unlist(strsplit(target, "\\|"))

  # Filter samples based on target patterns
  # Search in ID column which contains the full sample information
  sample_matches <- rep(FALSE, nrow(sample_metadata))
  for (pattern in target_patterns) {
    sample_matches <- sample_matches | grepl(pattern, sample_metadata$ID,
                                             ignore.case = TRUE)
    cat(paste0("\nMatches for pattern '", pattern, "':\n"))
    cat(sample_matches)
  }

  if (sum(sample_matches) == 0) {
    cat(paste0("\nWARNING: No samples found matching target '",
               target, "'. Skipping.\n"))
    next
  }

  # Subset samples for this target
  cat(paste0("\nFull sample metadata:\n"))
  print(sample_metadata)
  target_samples <- sample_metadata[sample_matches, ]
  cat(paste0("\nAfter subsetting for target '", target, "':\n"))
  print(target_samples)
  
  # Extract base sample names from ID column (remove .bam extension if present)
  bam_names <- sub("\\.bam$", "", target_samples$ID)

  cat(paste0("\nFound ", length(bam_names), " samples matching target '",
             target, "':\n"))
  print(bam_names)

  # Reconstruct expected BAM file paths based on the suffix parameter
  expected_bams <- file.path(py$PROCESSEDBAMDIR, paste0(bam_names, suffix))
  
  missing_files <- !file.exists(expected_bams)
  
  if (any(missing_files)) {
    cat("WARNING: The following expected BAM files were not found:\n")
    print(expected_bams[missing_files])
    cat("Skipping these files.\n")

    bams <- expected_bams[!missing_files]
    
    # Update base names from the ID metadata column to exclude missing files
    valid_bam_names <- sapply(bam_names, function(x) any(grepl(x, bams)))
    bam_names <- bam_names[valid_bam_names]
    target_samples <- subset(target_samples, ID %in% bam_names)
    
    if (length(bams) == 0) {
      cat(paste0("ERROR: No valid BAM files left for target '", target, "'. Skipping.\n"))
      next
    }
  } else {
    bams <- expected_bams
    cat("\nFound all expected BAM files:\n")
    print(bams)
  }

  # WRITE UPDATED METADATA FILE FOR THIS TARGET
  target_clean <- run_name

  # Create metadata file for this target
  meta_file <- target_samples
  # Update ID column to point to the deduplicated BAM files
  meta_file$ID <- basename(bams)

  # Check the metadata file used for this target
  cat(paste0("\nTarget-specific sample metadata used for ",
             target_clean,
             " normalization:\n"))
  print(meta_file)

  # Export metadata file for this target
  meta_file_path <- file.path(py$METADATA, paste0("sample_metadata_",
                                                  target_clean,
                                                  "_processed.txt"))
  write.table(meta_file,
              file = meta_file_path,
              sep = "\t",
              row.names = FALSE,
              col.names = TRUE,
              quote = FALSE)

  ### Generate output folder
  filename <- target_clean

  target_dir <- file.path(output_arg, filename)
  if (bin_size != 10000) {
    target_dir <- file.path(output_arg,
                            paste0(filename, "_", bin_size, "bp_bins"))
  }

  cat(paste0("\nCreating output directory: ", target_dir, "\n"))

  if (!dir.exists(target_dir)) {
    dir.create(target_dir, recursive = TRUE)
  }

  ### Generate output prefix
  if (cutoff != 1.2) {
    filename <- paste0(filename, "_cutoff_", cutoff)
  }
  if (max_turns != 0.99) {
    filename <- paste0(filename, "_max_turns_", max_turns)
  }

  prefix <- file.path(target_dir, filename)

  # Check whether rawCounts and parsedMatrix files from a previous run exist
  existing_raw <- list.files(target_dir, pattern = "_rawCounts\\.txt$", full.names = TRUE)
  existing_parsed <- list.files(target_dir, pattern = "_parsedMatrix\\.txt$", full.names = TRUE)

  if (length(existing_raw) > 0 && length(existing_parsed) > 0) {
    cat("\nCounts and parsed matrix exist from a previous run. Copying them to reuse...\n")
    src_raw <- existing_raw[1]
    src_parsed <- existing_parsed[1]
    
    dest_raw <- paste0(prefix, "_rawCounts.txt")
    dest_parsed <- paste0(prefix, "_parsedMatrix.txt")
    
    if (src_raw != dest_raw) {
      cat(paste0("Copying ", basename(src_raw), " to ", basename(dest_raw), "\n"))
      file.copy(src_raw, dest_raw, overwrite = TRUE)
    }
    if (src_parsed != dest_parsed) {
      cat(paste0("Copying ", basename(src_parsed), " to ", basename(dest_parsed), "\n"))
      file.copy(src_parsed, dest_parsed, overwrite = TRUE)
    }
  } else {
    cat("\nCounts or parsed matrix do not exist. ChIPseqSpikeInFree will generate them from BAM files.\n")
  }
  ### Use csaw to count reads in bins of specified size
  binned_filename <- paste0(target_clean, "_10kb_bins.txt")
  binned_counts <- file.path(output_arg, binned_filename)
  
  if (file.exists(binned_counts)) {
    cat(paste0("\nBinned counts file already exists for target '",
               target, "'. Skipping read counting step.\n"))
    counts_df <- read.table(binned_counts,
                            sep = "\t",
                            header = TRUE,
                            row.names = 1)
    raw_counts <- as.matrix(counts_df)
  } else {
    cat("\nCounting reads in bins...\n")
    if (paired_arg == "TRUE") {
      pe <- "both"
    } else {
      pe <- "none"
    }
    param <- readParam(minq = 10,
                      pe = pe,
                      max.frag = 500,
                      dedup = FALSE)

    counts <- windowCounts(bams,
                          bin = TRUE,
                          width = 10000,
                          param = param)
    raw_counts <- assay(counts)
  }

  # Add sample names to counts object
  bamnames <- basename(bams)
  colnames(raw_counts) <- bamnames

  # Export binned count table
  write.table(as.data.frame(raw_counts),
              file = binned_counts,
              sep = "\t",
              row.names = TRUE,
              col.names = TRUE,
              quote = FALSE)

  count_table <- file.path(paste0(prefix, "_rawCounts.txt"))
  # Export binned count table for use in ChIPseqSpikeInFree (optional)
  # write.table(as.data.frame(raw_counts),
  #             file = count_table,
  #             sep = "\t",
  #             row.names = TRUE,
  #             col.names = TRUE,
  #             quote = FALSE)

  ### Run pipeline
  cat("\nRunning ChIPseqSpikeInFree pipeline...\n")
  cat(paste0("Output prefix: ", prefix, "\n\n"))

  tryCatch({
    ChIPseqSpikeInFree(
      bamFiles = bams,
      chromFile = chrom_sizes,
      metaFile = meta_file_path,
      prefix = prefix,
      binSize = bin_size,
      cutoff_QC = cutoff,
      maxLastTurn = max_turns,
      ncores = 8
    )
    cat(paste0("\nSuccessfully completed normalization for target: ",
               target, "\n"))

  }, error = function(e) {
    cat(paste0("\nERROR processing target '", target, "': ", e$message, "\n"))
  })


  # Load output table with normalization factors
  output_table <- file.path(paste0(prefix, "_SF.txt"))
  dat <-  read.table(output_table,
                     sep = "\t",
                     header = TRUE,
                     fill = TRUE,
                     stringsAsFactors = FALSE,
                     quote = "",
                     check.names = FALSE)

  # Extract normalization factors
  normfacs <- dat$SF

  # Treat Inf values identically to NA: both represent cases where the
  # algorithm could not compute a meaningful scaling factor.
  normfacs[is.infinite(normfacs)] <- NA

  if (any(is.na(normfacs))) {
    cat(paste0("\nWARNING: NA or Inf values found in SF for target '",
               target, "'. Rescaling normalization factors for this target.\n"))

    # Rename original SF table to keep a copy of the original values
    write.table(dat,
                file.path(paste0(prefix, "_original_SF.txt")),
                sep = "\t",
                quote = FALSE,
                row.names = FALSE,
                col.names = TRUE)

    #' Rescale normalization factors
    #' Samples with complete enrichment loss (NA / Inf) are imputed to the
    #' maximum finite SF observed across the run.
    #' This way samples with complete enrichment loss receive the highest SF.
    max_normfacs <- max(normfacs, na.rm = TRUE)
    if (is.infinite(max_normfacs) || is.na(max_normfacs)) {
      max_normfacs <- 1
    }
    normfacs[is.na(normfacs)] <- max_normfacs
    dat$original_SF <- dat$SF
    dat$SF <- normfacs

    # Export rescaled SF table
    write.table(dat,
                output_table,
                sep = "\t",
                quote = FALSE,
                row.names = FALSE,
                col.names = TRUE)
  }

  #----------
  # CHECK NORMALIZATION WITH MA PLOT
  #----------

  #' The log-ratio of normalization factors should pass through the
  #' center of the cloud in the plot.
  #' Clouds at low A-values represent background
  #' Clouds at high A-values represent bound regions

  # Import raw counts table
  # counts <- read.table(count_table,
  #                      sep = "\t",
  #                      header = TRUE,
  #                      row.names = 1)

  # Reorder columns to match the order of samples in the SF table
  counts <- raw_counts[, dat$ID]

  # Extract sample names from the count table
  sample.names <- colnames(counts)

  ##  Calculate read abundances (log2 CPM)
  adj.counts <- edgeR::cpm(counts, log = TRUE)

  cat("\nExporting MA plot\n")
  pdf_name <- file.path(paste0(prefix, "_MA_plot.pdf"))

  pdf(pdf_name, width = 11, height = 7.5)
  par(mfrow = c(2, 3), mar = c(5, 4, 5, 4))
  for (i in seq_len(length(sample.names) - 1)) {
    cur.x <- adj.counts[, 1] # signal abundances for sample 1
    cur.y <- adj.counts[, 1 + i]

    sample1 <- gsub(suffix, "", sample.names[1])
    sample2 <- gsub(suffix,"", sample.names[i + 1])
    smoothScatter(x = (cur.x + cur.y) / 2,
                  y = cur.x - cur.y,
                  xlab = "Mean signal (logCPM) across samples (A)",
                  ylab = "Log2 Fold-Change (M)",
                  main = paste(sample1,
                               "\n", "vs\n",
                               sample2))
    # Difference between normalization factors, log2-normalized
    all.dist <- diff(log2(normfacs[c(i + 1, 1)]))
    # Reference line marking the log2-difference between normalization factors
    abline(h = all.dist, col = "red")

  }
  cat("Closing pdf")
  dev.off()

  cat(paste0("\nCompleted processing for target: ", target, "\n"))
}



cat("\n========================================\n")
cat("All targets processed.\n")
cat("========================================\n")
