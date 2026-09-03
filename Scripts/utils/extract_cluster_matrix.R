#!/usr/bin/env Rscript
# =============================================================================
# extract_cluster_matrix.R
# Extracts binned profiles for a chosen subset of samples from norm_mats.rds
# and writes a gzip TSV file for k-means clustering in step_8e.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
  stop("Usage: Rscript extract_cluster_matrix.R <norm_mats_rds> <out_tsv_gz> <sample1> [sample2 ...]")
}

rds_path <- args[1L]
out_path <- args[2L]
samples  <- args[3L:length(args)]

if (!file.exists(rds_path)) {
  stop("norm_mats file not found: ", rds_path)
}

norm_mats <- readRDS(rds_path)
missing_s <- setdiff(samples, names(norm_mats))
if (length(missing_s) > 0L) {
  stop("Requested samples not found in norm_mats.rds: ", paste(missing_s, collapse = ", "),
       "\nAvailable samples: ", paste(names(norm_mats), collapse = ", "))
}

cat(sprintf("Extracting %d samples for cluster matrix: %s\n",
            length(samples), paste(samples, collapse = " + ")))

mat_list <- lapply(samples, function(s) as.matrix(norm_mats[[s]]))
cluster_mat <- do.call(cbind, mat_list)
cluster_mat[is.na(cluster_mat)] <- 0

gz_con <- gzcon(file(out_path, "wb"))
write.table(cluster_mat, gz_con, sep = "\t", row.names = FALSE, col.names = FALSE)
close(gz_con)

cat(sprintf("Cluster matrix saved (%d regions x %d bins): %s\n",
            nrow(cluster_mat), ncol(cluster_mat), out_path))
