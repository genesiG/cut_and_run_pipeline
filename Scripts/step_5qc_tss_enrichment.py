#!/usr/bin/env python3
"""
step_5qc_tss_enrichment.py - Calculate TSS enrichment for all samples
"""

import os
import subprocess
import glob
import config


# Directories
tss_dir = os.path.join(config.QCDIR1, "tss_enrichment")
os.makedirs(tss_dir, exist_ok=True)
log_dir = os.path.join(tss_dir, "log")
os.makedirs(log_dir, exist_ok=True)

NUM_CORES = 16
MAX_MEM = 64000
MEM_PER_CORE = MAX_MEM / NUM_CORES

norm_method = getattr(config, 'NORMALIZE_USING', 'None').lower()

def get_bigwig_files_for_group(group_pattern):
    """
    Find all BigWig files in config.BIGWIGDIR that match a group pattern.
    Example:
        group_pattern = "DMSO" will match all files containing "DMSO".
    """
    all_bw = glob.glob(os.path.join(config.BIGWIGDIR,
                                    f"*_R1.{norm_method}.bin_size_{config.BIN_SIZE}.bw"))
    return [bw for bw in all_bw if group_pattern in os.path.basename(bw)]

def get_clean_sample_names(bigwig_files):
    """
    Strip the suffix from BigWig file names to get clean sample names.
    Removes: .{norm_method}.bin_size_{config.BIN_SIZE}.bw
    """
    suffix = f".{norm_method}.bin_size_{config.BIN_SIZE}.bw"
    return [os.path.basename(bw).replace(suffix, '') for bw in bigwig_files]

def run_tss_analysis_group(group_name, bigwig_files):
    """
    Run TSS enrichment analysis for a group of samples.
    Steps:
    1. Compute matrix around TSS using all BigWig files in the group
    2. Plot profile and heatmap for the group
    3. Calculate enrichment score
    """
    # Join BigWig files into space-separated string for computeMatrix
    bigwig_input = " ".join(bigwig_files)
    sample_labels = " ".join(get_clean_sample_names(bigwig_files))
    print(f"Sample labels for {group_name}: {sample_labels}")
    sample_labels = sample_labels.replace("_", "\n")  # Escape underscores for plotting

    # Output files (group-based)
    matrix_file = os.path.join(tss_dir, f"{group_name}_tss_matrix.gz")
    plot_file = os.path.join(tss_dir, f"{group_name}_tss_profile.pdf")
    heatmap_file = os.path.join(tss_dir, f"{group_name}_tss_heatmap.pdf")
    data_file = os.path.join(tss_dir, f"{group_name}_tss_data.txt")

    if os.path.exists(matrix_file):
        print("[INFO] TSS matrix already exists, skipping computeMatrix step.")
        matrix_cmd = ""
    else:
        matrix_cmd = f"""
computeMatrix reference-point \\
  -R {config.TSS_BED} \\
  -S {bigwig_input} \\
  --referencePoint TSS \\
  -a 2000 -b 2000 \\
  --binSize 50 \\
  --skipZeros \\
  --missingDataAsZero \\
  -o {matrix_file} \\
  -p {NUM_CORES}
  """

    # Log & batch script
    log_out = os.path.join(log_dir, f"{group_name}.log")
    log_err = os.path.join(log_dir, f"{group_name}.error")
    batch_file = os.path.join(log_dir, f"{group_name}.batch")

    batch_cmd = f"""#!/bin/bash
#BSUB -P {config.PROJECT_NAME}
#BSUB -J tss_{group_name}
#BSUB -oo {log_out}
#BSUB -eo {log_err}
#BSUB -q {config.BSUB_QUEUE}
#BSUB -n {NUM_CORES}
#BSUB -M {MAX_MEM}
#BSUB -R "rusage[mem={MEM_PER_CORE}] span[hosts=1]"

# Activate environment
source {config.HOMEDIR}/miniconda3/etc/profile.d/conda.sh
conda activate chipseq

### 1. Compute matrix around TSS for group
{matrix_cmd}

### 2. Plot profile
plotProfile \\
  -m {matrix_file} \\
  -out {plot_file} \\
  --samplesLabel '{sample_labels}' \\
  --plotTitle "{group_name} TSS Enrichment" \\
  --plotType "se" \\
  --outFileNameData {data_file}

### 3. Plot heatmap
plotHeatmap \\
  -m {matrix_file} \\
  -out {heatmap_file} \\
  --samplesLabel '{sample_labels}' \\
  --colorList "lightblue,yellow,red" \\
  --plotTitle "{group_name} TSS Enrichment" \\
  --whatToShow "heatmap and colorbar"

### 4. Calculate enrichment score
max_signal=$(awk 'NR>1 {{print $NF}}' {data_file} | sort -nr | head -1)
min_signal=$(awk 'NR>1 {{print $NF}}' {data_file} | sort -n | head -1)
tss_enrichment=$(echo "$max_signal / $min_signal" | bc -l)

echo "TSS enrichment score: $tss_enrichment" > {tss_dir}/{group_name}_tss_score.txt

echo "=== Resource usage summary ==="
bjobs -l $LSB_JOBID | grep -A10 "Resource usage"
"""

    with open(batch_file, "w") as f:
        f.write(batch_cmd)

    print(f"Submitting TSS analysis for group: {group_name}")
    subprocess.run(f"bsub < {batch_file}", shell=True)

def main():
    """
    Main function to run TSS enrichment analysis per group.
    """
    print("Starting group-based TSS enrichment analysis")

    # Ensure TSS BED file exists
    if not os.path.exists(config.TSS_BED):
        raise FileNotFoundError(f"TSS BED file not found at {config.TSS_BED}.")

    # Extract group names from config.TARGETS (comma-separated)
    groups = [g.strip() for g in config.TARGETS.split(",")]

    for group in groups:
        bigwig_files = get_bigwig_files_for_group(group)
        if bigwig_files:
            bigwig_files.sort()  # Ensures alphabetical order
            run_tss_analysis_group(group, bigwig_files)
        else:
            print(f"Warning: No BigWig files found for group {group}")

    print(f"Submitted TSS analysis for {len(groups)} groups")

if __name__ == "__main__":
    main()
