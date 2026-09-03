#!/usr/bin/env python3

"""
step_5_generate_bigwig.py - Script to generate scaled bigWig files using bamCoverage
Steps:
1. Determine normalization method from config
2. For each sample, build batch script with appropriate bamCoverage parameters
3. Submit jobs to cluster using bsub
"""

### Import modules
import os
import time
import glob
import argparse
import pandas as pd
import config
###

# Working directory
work_dir = os.path.join(config.BIGWIGDIR)
log_dir = os.path.join(work_dir, "log")
bin_size = config.BIN_SIZE
suffix = (
    f"qc.sort.rmdup.mapq{config.MAPQ}.final"
    if getattr(config, 'REMOVE_DUPLICATES', False)
    else f"qc.sort.markdup.mapq{config.MAPQ}.final"
)

# Create directories if needed
os.makedirs(work_dir, exist_ok=True)
os.makedirs(log_dir, exist_ok=True)

def run_bigwigcompare(sample_name, norm_method, normfactor=None):
    """
    Creates and submits a batch script to generate bigWig files
    """
    # Remove suffix from filename
    out_name = sample_name.replace(f".{suffix}.bam",'')

    # Job name based on sample name
    job_name = f"{out_name}.{bin_size}.bigwig"

    # Input files
    input_bam = os.path.join(config.PROCESSEDBAMDIR, sample_name)

    # Output files
    norm = norm_method.lower()
    output_bw = os.path.join(work_dir, f"{out_name}.{norm}.bw")
    if norm_method.lower() in ["chipseqspikeinfree", "spikein"]:
        output_bw = output_bw.replace(f".{norm}.bw", f".{norm}.rpkm.bw")

    if os.path.exists(output_bw):
        print(f"Skipping {job_name}, {output_bw} already exists.")
        return

    # Log files
    log_out = os.path.join(log_dir, f"{job_name}.log")
    log_err = os.path.join(log_dir, f"{job_name}.error")
    batch_file = os.path.join(log_dir, f"{job_name}.batch")

    # Build bamCoverage command
    bamcoverage_cmd = f"""
bamCoverage -b {input_bam} \\
    -p 8 \\
    --binSize {bin_size} \\
    --skipNonCoveredRegions \\
    --centerReads \\
    --outFileFormat bigwig \\
    -o {output_bw} \\"""

    # Add built-in normalization method if specified
    if norm_method.lower() in ["cpm", "rpkm", "bpm", "rpgc"]:
        bamcoverage_cmd += f"""
    --normalizeUsing {norm_method.upper()} \\"""
    elif norm_method.lower() == "chipseqspikeinfree":
        bamcoverage_cmd += """
    --normalizeUsing RPKM \\"""
    elif norm_method.lower() == "spikein":
        bamcoverage_cmd += """
    --normalizeUsing RPKM \\"""

    # Handle missing or NA scale factors
    if normfactor is None or str(normfactor).strip().lower() in ["na", "nan"]:
        normfactor = 1

    # Add scale factor
    bamcoverage_cmd += f"""
    --scaleFactor {normfactor} \\
    """

    # Build batch script
    batch_content = f"""#!/bin/bash
#BSUB -P {config.PROJECT_NAME}
#BSUB -J {job_name}
#BSUB -oo {log_out}
#BSUB -eo {log_err}
#BSUB -q {config.BSUB_QUEUE}
#BSUB -n 8
#BSUB -M 16384
#BSUB -R "rusage [mem=2048] span[hosts=1]"

# Load environment
source {config.HOMEDIR}/miniconda3/etc/profile.d/conda.sh
conda activate chipseq

# Run bamCoverage
{bamcoverage_cmd}

# Print resource usage at the end of the job
echo "=== Resource usage summary ==="
bjobs -l $LSB_JOBID | grep -A10 "Resource usage"
"""

    # Write batch file
    with open(batch_file, 'w') as f:
        f.write(batch_content)

    print(f"Submitting job: {job_name}")
    os.system(f"bsub < {batch_file}")
    time.sleep(1)

def main():
    """
    Generate scaled bigWig files from alignment BAM files using DeepTools bamCoverage.
    """
    parser = argparse.ArgumentParser(description="Generate scaled bigWig files using bamCoverage")
    parser.add_argument("--patterns",
                        nargs="*",
                        help="Optional list of patterns to match sample names (e.g. IgG K27ac)")
    args = parser.parse_args()

    # Determine normalization method
    norm_method = getattr(config, 'NORMALIZE_USING', 'None').upper()

    # 1. Gather all samples from metadata
    try:
        sample_file = "samples.txt"
        if hasattr(config, 'IS_PAIRED_END') and config.IS_PAIRED_END:
            sample_file = "paired_samples.txt"

        samples_path = os.path.join(config.METADATA, sample_file)
        samples_df = pd.read_csv(samples_path, sep='\t')
        all_sample_names = samples_df['sample.name'].tolist()
    except Exception as e:
        print(f"ERROR reading sample metadata: {str(e)}")
        return

    # Filter by patterns if requested
    if args.patterns:
        all_sample_names = [s for s in all_sample_names
                            if any(p.lower() in s.lower() for p in args.patterns)]
        print(f"Filtered to {len(all_sample_names)} samples matching patterns: {args.patterns}")
    else:
        print(f"Processing {len(all_sample_names)} samples.")

    # 2. Pre-load scaling factors if applicable
    sf_lookup = {}
    if norm_method in ["CHIPSEQSPIKEINFREE", "TMM", "SPIKEIN"]:
        scaling_files = []
        if norm_method == "CHIPSEQSPIKEINFREE":
            for run_name, params in config.NORMALIZATION_RUNS.items():
                b_size = params["bin_size"]
                cutoff = params.get("cutoff", 1.2)
                max_turns = params.get("max_turns", 0.99)

                run_dir = f"{run_name}_{b_size}bp_bins" if b_size != 10000 else run_name

                r_filename = run_name
                if cutoff != 1.2:
                    r_filename += f"_cutoff_{cutoff}"
                if max_turns != 0.99:
                    r_filename += f"_max_turns_{max_turns}"

                sf_path = os.path.join(config.SCALINGDIR,
                                       "chipseqspikeinfree",
                                       run_dir,
                                       f"{r_filename}_SF.txt")
                if os.path.exists(sf_path):
                    scaling_files.append(sf_path)
        elif norm_method == "TMM":
            scaling_files = glob.glob(os.path.join(config.SCALINGDIR, "tmm_scaling_factors*"))
        elif norm_method == "SPIKEIN":
            scaling_files = glob.glob(os.path.join(config.SCALINGDIR, "spikein", "*spikein_SF.txt"))

        if not scaling_files:
            print(f"WARNING: No scaling factor files found in {config.SCALINGDIR}")

        for sf_path in scaling_files:
            try:
                sf_df = pd.read_csv(sf_path, sep='\t')
                print(f"Loaded scaling factors from: {sf_path}")
                for _, row in sf_df.iterrows():
                    sample_id = row['ID']
                    scale_factor = (row['SF']
                                    if norm_method in ["CHIPSEQSPIKEINFREE", "TMM"]
                                    else row['bamCov_SF'])
                    sf_lookup[sample_id] = scale_factor
            except Exception as e:
                print(f"ERROR reading scaling file {sf_path}: {str(e)}")

    # 3. Process each sample
    for base_sample in all_sample_names:
        bam_name = f"{base_sample}.{suffix}.bam"
        bam_path = os.path.join(config.PROCESSEDBAMDIR, bam_name)

        if not os.path.exists(bam_path):
            print(f"WARNING: BAM file not found: {bam_name}")
            continue

        is_control = ("igg" in base_sample.lower()) or ("input" in base_sample.lower())

        if is_control:
            if norm_method == "TMM":
                # Use TMM if global is TMM
                if bam_name in sf_lookup:
                    scale_factor = sf_lookup[bam_name]
                    normfactor = 1 / scale_factor
                    run_bigwigcompare(bam_name, "TMM", normfactor)
                else:
                    print(f"WARNING: No TMM scale factor found for {bam_name}")
            else:
                # Fallback to CPM or RPKM
                fallback_method = "RPKM" if norm_method == "RPKM" else "CPM"
                run_bigwigcompare(bam_name, fallback_method, 1.0)
        else:
            if norm_method in ["CHIPSEQSPIKEINFREE", "TMM", "SPIKEIN"]:
                if bam_name in sf_lookup:
                    scale_factor = sf_lookup[bam_name]
                    # TMM and SpikeInFree scaling factors are meant for division
                    normfactor = (1 / scale_factor
                                  if norm_method in ["CHIPSEQSPIKEINFREE", "TMM"]
                                  else scale_factor)
                    run_bigwigcompare(bam_name, norm_method, normfactor)
                else:
                    print(f"WARNING: No {norm_method} scale factor found for {bam_name}")
            else:
                # Built-in methods
                nf = 1.0 if norm_method == "NONE" else 1.0
                run_bigwigcompare(bam_name, norm_method, nf)

if __name__ == "__main__":
    main()
