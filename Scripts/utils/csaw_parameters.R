# Requirements
packages <- c("dplyr", "BiocManager")
biopackages <- c("GenomicRanges", "GenomeInfoDb", "edgeR", "csaw")

# Define readParam object
param <- csaw::readParam(
    minq = 30, # minimum alingment quality score
    pe = ifelse(py$IS_PAIRED_END, "both", "none"), # paired-end data
    dedup = FALSE
) # CUT and RUN data

# Bin size for read counting
bin.size <- 10000

# Additional parameters to count reads into windows
extension <- 200
width.small <- 150
spacing.small <- 50
width.large <- 2000
spacing.large <- 500

# Log2 fold-enrichment threshold over background
# for filtering low-abundance windows
small.filt <- 3.0
large.filt <- 1.0

# Additional parameters for genomic analysis
build <- "hs1"
style <- "UCSC"

# Parameters for peak calling / differential binding
merge.bp <- 500
lfc <- 1
qvalue <- 0.05
