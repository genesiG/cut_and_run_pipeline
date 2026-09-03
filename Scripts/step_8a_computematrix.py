#!/usr/bin/env python3

"""
step_5_compute_matrix.py - Generate matrices with signal over genome regions
Steps:
1. Define signal (bigwig) files and genome regions (peaks/BED files)
2. Build batch script with computeMatrix commands
3. Submit jobs to the cluster using bsub
"""

### Import modules
import os
import time
import pandas as pd
import glob
import argparse

import config
###

### Configuration
scaling_file = "scale_UNC_vs_DMSO_SF.txt"
avg_bin = "mean" # mean, median, min, max, std, sum
bin_size = 100 # (bp) Default: 10
base_prefix = "retained_h3k27me2"

# --- Window configurations: (upstream_bp, downstream_bp) ---
# Each entry produces an independent matrix with a unique output name.
WINDOW_CONFIGS = [
    (10000, 10000),
    #(5000, 5000),
    #(2000, 2000),   # +/- 2 kb
]
# ---

regions_dir = config.PEAKDIR
scale_regions = False
scale_to = 1000 # (bp) only used if scale_regions is True
reference_point = "center" # "center", "TSS", "TES" (only used if scale_regions is False)
missing_data_as_zero = True
skip_zeros = False
scale = 1 # Default: 1
suffix = config.NORMALIZE_USING.lower()

# Input files
bw_files = [
    os.path.join(config.BIGWIGDIR, "RP_036_H3K27me2_DMSO_R1.spikein.rpkm.bw"),
    os.path.join(config.BIGWIGDIR, "RP_060_H3K27me2_1uMcpd_R1.spikein.rpkm.bw"),
    os.path.join(config.BIGWIGDIR, "RP_038_H3K27me3_DMSO_R1.spikein.rpkm.bw"),
    os.path.join(config.BIGWIGDIR, "RP_062_H3K27me3_1uMcpd_R1.spikein.rpkm.bw"),
    os.path.join(config.BIGWIGDIR, "RP_039_AntiCBX2RP_DMSO_R1.spikein.rpkm.bw"),
    os.path.join(config.BIGWIGDIR, "RP_063_AntiCBX2RP_1uMcpd_R1.spikein.rpkm.bw"),
    os.path.join(config.BIGWIGDIR, "RP_042_EZH2_DMSO_R1.spikein.rpkm.bw"),
    os.path.join(config.BIGWIGDIR, "RP_066_EZH2_1uMcpd_R1.spikein.rpkm.bw"),
]
region_files = [
    # os.path.join(config.PEAKDIR,
    # "csaw",
    # "H3K27me2_1uMcpd.w150.d50.filt3.w2000.d500.filt1.lfc1.merge100.tmm.bed"),
    os.path.join(config.ANALYSIS_DATA,
    "k27me2_classification",
    "retained_k27me2.bed"
    #"intersections",
    #"lost_h3k27me2_in_lost_h3k27me3.bed"
    )
]

###

def run_compute_matrix(upstream_region: int, downstream_region: int):
    """
    Creates and submits a batch script for computeMatrix operations
    for a single (upstream, downstream) window around the reference point.
    """
    # Derive a compact suffix from the window size, e.g. 2000/2000 -> "2kb"
    half_kb = upstream_region // 1000
    output_suffix = f"{half_kb}kb"
    output_prefix = f"{base_prefix}_{output_suffix}"

    # Workspace configuration
    job_name = f"{output_prefix}_computeMatrix"
    work_dir = config.DEEPTOOLSDIR
    matrix_dir = os.path.join(work_dir, "matrix")

    # Output files
    matrix_file = os.path.join(matrix_dir, f"{output_prefix}.gz")

    # Command line parameters
    command_param = []
    if scale_regions:
        command_param.append("scale-regions")
        command_param.append(f"--regionBodyLength {scale_to}")
    else:
        command_param.append("reference-point")
        command_param.append(f"--referencePoint {reference_point}")

    optional_param = []
    if missing_data_as_zero:
        optional_param.append("--missingDataAsZero")
    if skip_zeros:
        optional_param.append("--skipZeros")

    # Log files
    log_dir = os.path.join(matrix_dir, "log")
    log_out = os.path.join(log_dir, f"{output_prefix}.log")
    log_err = os.path.join(log_dir, f"{output_prefix}.error")
    batch_file = os.path.join(log_dir, f"{output_prefix}_computeMatrix.batch")

    # Create workspace environment
    os.makedirs(work_dir, exist_ok=True)
    os.makedirs(matrix_dir, exist_ok=True)
    os.makedirs(log_dir, exist_ok=True)

    # Build batch command
    batch_cmd = f"""\
#!/bin/bash
#BSUB -P {config.PROJECT_NAME}
#BSUB -J {job_name}
#BSUB -oo {log_out}
#BSUB -eo {log_err}
#BSUB -q {config.BSUB_QUEUE}
#BSUB -n 6  # number of cores to request
#BSUB -M 49152 # total memory limit (MB)
#BSUB -R "rusage [mem=8192] span[hosts=1]" # Reserve the different CPU cores under the same node

source {config.HOMEDIR}/miniconda3/etc/profile.d/conda.sh
conda activate chipseq

computeMatrix {' '.join(command_param)} \\
  -S {' '.join(bw_files)} \\
  -R {' '.join(region_files)} \\
  -b {upstream_region} -a {downstream_region} \\
  -p 8 \\
  -bs {bin_size} \\
  {' '.join(optional_param)} \\
  --averageTypeBins {avg_bin} \\
  --scale {scale} \\
  -o {matrix_file}

# Print resource usage at the end of the job
echo "=== Resource usage summary ==="
bjobs -l $LSB_JOBID | grep -A10 "Resource usage"
"""

    # Write batch file
    with open(batch_file, 'w', encoding='utf-8') as f:
        f.write(batch_cmd)

    print(f"Submitting {job_name}")
    os.system(f"bsub < {batch_file}")
    time.sleep(1)

def main():
    """Main workflow controller — submits one job per window config"""
    parser = argparse.ArgumentParser(description="Generate matrices with signal over genome regions")
    parser.add_argument("--base_prefix", type=str, help="Overwrite base_prefix")
    parser.add_argument("--avg_bin", type=str, choices=["mean", "median", "min", "max", "std", "sum"], help="Overwrite avg_bin")
    parser.add_argument("--bin_size", type=int, help="Overwrite bin_size")
    parser.add_argument("--normfacs", type=str, help="Overwrite scaling_file")
    parser.add_argument("--regions", type=str, nargs="+", help="Overwrite region_files")
    args = parser.parse_args()

    global base_prefix, avg_bin, bin_size, scaling_file, region_files

    if args.avg_bin:
        avg_bin = args.avg_bin

    if args.base_prefix:
        base_prefix = args.base_prefix
    elif args.avg_bin:
        # Re-evaluate base_prefix if avg_bin was provided but base_prefix wasn't
        base_prefix = "k27me2_1uM_peaks"

    if args.bin_size:
        bin_size = args.bin_size

    if args.normfacs:
        scaling_file = args.normfacs
    
    if args.regions:
        region_files = args.regions


    # Log summary of files being used
    print("Found bigWig files:")
    for f in bw_files:
        print(f"   {f}")

    print("\nComputing signal over regions:")
    for f in region_files:
        print(f"   {f}")

    for upstream, downstream in WINDOW_CONFIGS:
        run_compute_matrix(upstream, downstream)

if __name__ == "__main__":
    main()
