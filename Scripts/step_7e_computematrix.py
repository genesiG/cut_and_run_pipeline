#!/usr/bin/env python3

"""
step_7e_computematrix.py - Generate matrices with signal over H3K27me3 DMSO peaks.

Mirrors step_8a_computematrix.py but targets the H3K27me3 DMSO csaw peaks directly.

Same sample order and bigWig files as step_8a (H3K27me2, H3K27me3, CBX2, EZH2).
"""

### Import modules
import os
import time
import argparse

import config
###

### Configuration
avg_bin   = "mean"   # mean, median, min, max, std, sum
bin_size  = 100      # (bp) Default: 10
base_prefix = "k27me3_dmso_peaks"

# --- Window configurations: (upstream_bp, downstream_bp) ---
WINDOW_CONFIGS = [
    (10000, 10000),
    (5000, 5000),
    (2000, 2000),   # +/- 2 kb
]
# ---

scale_regions    = False
scale_to         = 1000   # (bp) only used if scale_regions is True
reference_point  = "center"   # "center", "TSS", "TES"
missing_data_as_zero = True
skip_zeros       = False
scale            = 1     # Default: 1

# Input bigWig files — same order as step_8a
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

# Region file: H3K27me3 DMSO csaw peaks
region_files = [
    os.path.join(
        config.PEAKDIR, "csaw",
        "H3K27me3_DMSO.w150.d50.filt3.w2000.d500.filt2.lfc1.merge500.tmm.bed"
    )
]
###


def run_compute_matrix(upstream_region: int, downstream_region: int):
    """
    Creates and submits a batch script for computeMatrix operations
    for a single (upstream, downstream) window around the reference point.
    """
    half_kb       = upstream_region // 1000
    output_suffix = f"{half_kb}kb"
    output_prefix = f"{base_prefix}_{output_suffix}"

    job_name   = f"{output_prefix}_computeMatrix"
    work_dir   = config.DEEPTOOLSDIR
    matrix_dir = os.path.join(work_dir, "matrix")
    matrix_file = os.path.join(matrix_dir, f"{output_prefix}.gz")

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

    log_dir    = os.path.join(matrix_dir, "log")
    log_out    = os.path.join(log_dir, f"{output_prefix}.log")
    log_err    = os.path.join(log_dir, f"{output_prefix}.error")
    batch_file = os.path.join(log_dir, f"{output_prefix}_computeMatrix.batch")

    os.makedirs(work_dir,   exist_ok=True)
    os.makedirs(matrix_dir, exist_ok=True)
    os.makedirs(log_dir,    exist_ok=True)

    batch_cmd = f"""\
#!/bin/bash
#BSUB -P {config.PROJECT_NAME}
#BSUB -J {job_name}
#BSUB -oo {log_out}
#BSUB -eo {log_err}
#BSUB -q {config.BSUB_QUEUE}
#BSUB -n 6
#BSUB -M 49152
#BSUB -R "rusage [mem=8192] span[hosts=1]"

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

    with open(batch_file, 'w', encoding='utf-8') as f:
        f.write(batch_cmd)

    print(f"Submitting {job_name}")
    os.system(f"bsub < {batch_file}")
    time.sleep(1)


def main():
    """Main workflow controller — submits one job per window config"""
    parser = argparse.ArgumentParser(
        description="Generate matrices with signal over H3K27me3 DMSO peaks"
    )
    parser.add_argument("--base_prefix", type=str,
                        help="Overwrite base_prefix (default: k27me3_dmso_peaks)")
    parser.add_argument("--avg_bin", type=str,
                        choices=["mean", "median", "min", "max", "std", "sum"],
                        help="Overwrite avg_bin")
    parser.add_argument("--bin_size", type=int, help="Overwrite bin_size")
    parser.add_argument("--regions", type=str, nargs="+", help="Overwrite region_files")
    args = parser.parse_args()

    global base_prefix, avg_bin, bin_size, region_files

    if args.avg_bin:
        avg_bin = args.avg_bin

    if args.base_prefix:
        base_prefix = args.base_prefix
    elif args.avg_bin:
        base_prefix = "k27me3_dmso_peaks"

    if args.bin_size:
        bin_size = args.bin_size

    if args.regions:
        region_files = args.regions

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
