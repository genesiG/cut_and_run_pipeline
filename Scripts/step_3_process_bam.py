#!/usr/bin/env python3

"""
step_3_process_bam.py - BAM processing pipeline for CUT&RUN / CUT&TAG / ChIP-seq.

For each sample, submits a batch job that performs the following steps:

  A. QC filtering
     samtools view applies the MAPQ threshold and bitwise FLAG filters.
     Output still contains chrM and spike-in chromosomes -- needed for
     accurate pre-filter statistics and spike-in splitting (Step C).

  B. Pre-filter alignment statistics
     flagstat / idxstats on the QC-filtered BAM, capturing the full
     alignment landscape before any genome-specific filtering.

  C. Spike-in BAM splitting  [spike-in species only]
     The QC-filtered BAM is split into:
       - spike-in BAM  (E. coli or mouse spike-in chromosomes only)
       - discard       (the human split at this stage is NOT kept; the
                        final.bam produced in Step E is the correct human
                        denominator per the EpiCypher protocol)
     Splitting is performed here, before any human-genome filters, so that
     spike-in read counts are unaffected by blacklist or chrM removal.

  D. BAM-native human-genome filtering
     samtools view subsets to a whitelist of canonical human chromosomes
     (excludes chrM, spike-in contigs, and zero-read unplaced scaffolds).
     bedtools intersect -abam removes ENCODE DAC blacklisted regions while
     preserving all SAM auxiliary tags and FLAG values intact.

  E. Duplicate handling with Picard MarkDuplicates
     Operates on the fully filtered, tag-intact BAM from Step D.
     IgG controls always have duplicates removed regardless of the global
     REMOVE_DUPLICATES setting.

  F. Cleanup of intermediate files.

Output BAM naming convention
-----------------------------
  *.qc.mapq{MAPQ}.spike_in.bam
      E. coli (or mouse) spike-in reads only.
      MAPQ-filtered; NOT deduplicated, NOT blacklist-filtered.
      Used as the spike-in numerator in step_4c.

  *.qc.sort.markdup.mapq{MAPQ}.final.bam   [REMOVE_DUPLICATES = False]
  *.qc.sort.rmdup.mapq{MAPQ}.final.bam     [REMOVE_DUPLICATES = True]
      Human reads only: chrM removed, blacklist removed, coordinate-sorted,
      duplicate-handled by Picard.
      Used as the total-uniquely-aligned-reads denominator in step_4c,
      matching the EpiCypher protocol definition: "filter out multi-mapping
      reads, reads assigned to ENCODE DAC exclusion list regions, and
      duplicate reads (as desired) to determine the total number of
      uniquely aligned reads."
      Also the analysis-ready BAMs for bigWig generation, peak calling,
      and differential binding.
"""

import os
import time

import config
from utils import ldsample

# Working directories
work_dir          = config.BAMDIR
original_bam_dir  = config.ORIGINALBAMDIR
processed_bam_dir = config.PROCESSEDBAMDIR
qc_alignment_dir  = os.path.join(config.QCDIR1, "alignment")
log_dir           = os.path.join(processed_bam_dir, "log")

os.makedirs(work_dir, exist_ok=True)
os.makedirs(log_dir, exist_ok=True)
os.makedirs(processed_bam_dir, exist_ok=True)
os.makedirs(qc_alignment_dir, exist_ok=True)

NUM_CORES    = 4
MAX_MEM      = 32000
MEM_PER_CORE = MAX_MEM / NUM_CORES


def process_bam(sample_name):
    """Build and submit a batch script to process one sample."""

    # ------------------------------------------------------------------
    # File paths
    # ------------------------------------------------------------------
    job_name = f'{sample_name}.processbam'

    original_bam  = os.path.join(original_bam_dir,  f"{sample_name}.bam")
    output_sam    = os.path.join(original_bam_dir,  f"{sample_name}.sam")

    # QC-filtered BAM — contains chrM + spike-in; used only in Steps A–C
    output_bam    = os.path.join(processed_bam_dir, f"{sample_name}.qc.bam")

    # Intermediate BAMs from human-genome filtering (Step D)
    chroms_filtered_bam    = os.path.join(processed_bam_dir, f"{sample_name}.chroms_filtered.bam")
    blacklist_filtered_bam = os.path.join(processed_bam_dir, f"{sample_name}.blacklist_filtered.bam")

    # Coordinate-sorted, fully filtered BAM — input to Picard (Step E)
    bam_presort   = os.path.join(processed_bam_dir, f"{sample_name}.sort.final.bam")

    # QC output paths
    qc_flagstat   = os.path.join(qc_alignment_dir, f"{sample_name}.flagstat.txt")
    qc_idxstat    = os.path.join(qc_alignment_dir, f"{sample_name}.idxstat.txt")
    picard_metrics = os.path.join(qc_alignment_dir, f"{sample_name}.picard_metrics.txt")

    mapq = config.MAPQ
    log_out   = os.path.join(log_dir, f"{sample_name}.log")
    log_err   = os.path.join(log_dir, f"{sample_name}.error")
    batch_file = os.path.join(log_dir, f"{sample_name}.batch")

    # ------------------------------------------------------------------
    # samtools FLAG filters
    # -F 4    exclude unmapped reads
    # -F 8    exclude reads with unmapped mate (PE only)
    # -F 256  exclude secondary alignments
    # -F 512  exclude reads failing vendor QC
    # -F 1024 exclude PCR/optical duplicates (kept here; Picard marks them)
    # -F 2048 exclude supplementary alignments
    # ------------------------------------------------------------------
    if config.IS_PAIRED_END:
        samtools_flags = "-f 1 -f 2 -F 4 -F 8 -F 256 -F 512 -F 1024 -F 2048"
    else:
        samtools_flags = "-F 4 -F 256 -F 512 -F 1024 -F 2048"

    # ------------------------------------------------------------------
    # Step C: Spike-in split block
    # ------------------------------------------------------------------
    if config.SPECIES in ["t2t_ecoli", "t2t_mm39"]:
        spiked_chr_prefix = "ecoli_chr" if config.SPECIES == "t2t_ecoli" else "mm39_chr"

        # Spike-in BAM: MAPQ-filtered E. coli reads only.
        # Named *.qc.mapq{MAPQ}.spike_in.bam — the absence of 'markdup' / 'rmdup'
        # in the name is intentional: Picard has not run at this point.
        spiked_bam = os.path.join(
            processed_bam_dir,
            f"{sample_name}.qc.mapq{mapq}.spike_in.bam"
        )

        # The human split at this pre-filter stage is NOT retained as a named
        # output. The final.bam from Step E is the correct human denominator.
        spike_in_split_block = f"""
# =============================================================================
# STEP C: Spike-in BAM splitting (before human-genome filters)
# Spike-in reads are counted from MAPQ-filtered reads only — no deduplication,
# no blacklist filtering. This matches the EpiCypher protocol Step 2, which
# filters only for unique alignment before quantifying spike-in reads.
# =============================================================================
echo '=== Step C: Spike-in BAM splitting ==='

SPIKE_IN_CHROMS=$(samtools idxstats {output_bam} \\
    | grep '^{spiked_chr_prefix}' \\
    | awk '$3 > 0 {{print $1}}' \\
    | xargs)

if [[ -z "$SPIKE_IN_CHROMS" ]]; then
    echo "WARNING: No spike-in chromosomes found with prefix '{spiked_chr_prefix}'. Check alignment."
else
    echo "Spike-in chromosomes: $SPIKE_IN_CHROMS"
    samtools view -@ {NUM_CORES} -b \\
        {output_bam} $SPIKE_IN_CHROMS \\
        -o {spiked_bam}
    samtools index {spiked_bam}
    echo "Spike-in BAM written: {spiked_bam}"
fi

echo '=== Step C complete ==='
"""
        chrm_grep_exclusion = f"grep -v '^{spiked_chr_prefix}'"
        spike_cleanup = f"rm -f {spiked_bam}"  # NOT removed — needed by step_4c

    else:
        spike_in_split_block = "# Step C: Spike-in splitting skipped (non-spike-in species)"
        chrm_grep_exclusion  = "grep -v '^__placeholder_never_matches__'"
        spiked_bam           = None
        spike_cleanup        = ""

    # ------------------------------------------------------------------
    # Step D: BAM-native human-genome filtering
    # ------------------------------------------------------------------
    if config.REMOVE_BLACKLIST:
        blacklist_step = f"""
# Remove ENCODE DAC blacklisted regions.
# -abam preserves all SAM auxiliary tags and correct FLAG values.
bedtools intersect \\
    -v \\
    -abam {chroms_filtered_bam} \\
    -b {config.BLACKLIST} \\
    > {blacklist_filtered_bam}
rm -f {chroms_filtered_bam}
"""
        input_for_sort          = blacklist_filtered_bam
        rm_filter_intermediates = f"rm -f {blacklist_filtered_bam}"
    else:
        blacklist_step          = "# Blacklist removal skipped (REMOVE_BLACKLIST=False)"
        input_for_sort          = chroms_filtered_bam
        rm_filter_intermediates = f"rm -f {chroms_filtered_bam}"

    blacklist_filter_block = f"""
# =============================================================================
# STEP D: BAM-native human-genome filtering
# Removes chrM, spike-in contigs, unplaced scaffolds, and blacklisted regions
# entirely within BAM space — preserving all SAM tags and FLAG values.
# =============================================================================
echo '=== Step D: BAM-native filtering ==='

KEEP_CHROMS=$(samtools idxstats {output_bam} \\
    | grep -v '^chrM' \\
    | grep -v '^\\*' \\
    | {chrm_grep_exclusion} \\
    | awk '$3 > 0 {{print $1}}' \\
    | xargs)

echo "Chromosomes retained (first 10): $(echo $KEEP_CHROMS | tr ' ' '\\n' | head -10 | xargs)"

samtools view -@ {NUM_CORES} -b -h \\
    {output_bam} $KEEP_CHROMS \\
    -o {chroms_filtered_bam}

{blacklist_step}

samtools sort -m 2G -@ {NUM_CORES} \\
    {input_for_sort} \\
    -o {bam_presort}
samtools index {bam_presort}

echo "Filtered, sorted BAM ready: {bam_presort}"
echo '=== Step D complete ==='
"""

    # ------------------------------------------------------------------
    # Step E: Picard MarkDuplicates
    # ------------------------------------------------------------------
    bam_final_markdup = os.path.join(
        processed_bam_dir,
        f"{sample_name}.qc.sort.markdup.mapq{mapq}.final.bam"
    )

    if config.REMOVE_DUPLICATES:
        bam_final = os.path.join(
            processed_bam_dir,
            f"{sample_name}.qc.sort.rmdup.mapq{mapq}.final.bam"
        )
        picard_command = f"""
picard MarkDuplicates \\
    -INPUT {bam_presort} \\
    -OUTPUT {bam_final_markdup} \\
    -METRICS_FILE {picard_metrics} \\
    -OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \\
    -REMOVE_DUPLICATES false \\
    -TMP_DIR {processed_bam_dir}
samtools index {bam_final_markdup}

picard MarkDuplicates \\
    -INPUT {bam_presort} \\
    -OUTPUT {bam_final} \\
    -METRICS_FILE {picard_metrics.replace('.txt', '.rmdup.txt')} \\
    -OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \\
    -REMOVE_DUPLICATES true \\
    -TMP_DIR {processed_bam_dir}
samtools index {bam_final}
"""
    else:
        # Mark duplicates only; IgG controls are always fully deduplicated.
        bam_final = bam_final_markdup
        picard_command = f"""
if [[ "{sample_name}" == *"IgG"* ]]; then
    echo "IgG sample — removing duplicates regardless of REMOVE_DUPLICATES setting"
    picard MarkDuplicates \\
        INPUT={bam_presort} \\
        OUTPUT={bam_final} \\
        METRICS_FILE={picard_metrics} \\
        OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \\
        REMOVE_DUPLICATES=true \\
        TMP_DIR={processed_bam_dir}
else
    picard MarkDuplicates \\
        INPUT={bam_presort} \\
        OUTPUT={bam_final} \\
        METRICS_FILE={picard_metrics} \\
        OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \\
        REMOVE_DUPLICATES=false \\
        TMP_DIR={processed_bam_dir}
fi
samtools index {bam_final}
"""

    # ------------------------------------------------------------------
    # Step F: Cleanup
    # ------------------------------------------------------------------
    # Check if outputs already exist to avoid redundant jobs
    outputs_exist = False
    if config.SPECIES in ["t2t_ecoli", "t2t_mm39"]:
        if (os.path.exists(bam_final) and os.path.exists(f"{bam_final}.bai") and 
            os.path.exists(spiked_bam) and os.path.exists(f"{spiked_bam}.bai")):
            outputs_exist = True
    else:
        if os.path.exists(bam_final) and os.path.exists(f"{bam_final}.bai"):
            outputs_exist = True

    if outputs_exist:
        print(f"Skipping {job_name} - output files already exist.")
        return
    cleanup_block = f"""
rm -f {output_bam} {output_bam}.bai
{rm_filter_intermediates}
rm -f {bam_presort} {bam_presort}.bai
"""
    if config.REMOVE_DUPLICATES:
        # markdup BAM served only as a duplicate-metrics intermediate
        cleanup_block += f"rm -f {bam_final_markdup} {bam_final_markdup}.bai\n"

    # ------------------------------------------------------------------
    # Assemble batch script
    # ------------------------------------------------------------------
    batch_cmd = f"""\
#!/bin/bash
#BSUB -P {config.PROJECT_NAME}
#BSUB -J {job_name}
#BSUB -oo {log_out}
#BSUB -eo {log_err}
#BSUB -q {config.BSUB_QUEUE}
#BSUB -n {NUM_CORES}
#BSUB -M {MAX_MEM}
#BSUB -R "rusage[mem={MEM_PER_CORE}] span[hosts=1]"

source {config.HOMEDIR}/miniconda3/etc/profile.d/conda.sh
conda activate chipseq

# Remove SAM output to save storage space
rm -f {output_sam}

# =============================================================================
# STEP A: QC filtering
# =============================================================================
echo '=== Step A: QC filtering ==='
samtools view \\
    -@ {NUM_CORES} -h -b \\
    -q {mapq} \\
    {samtools_flags} \\
    {original_bam} \\
    -o {output_bam}

samtools sort -m 2G -@ {NUM_CORES} {output_bam} -o {output_bam}.tmp
mv {output_bam}.tmp {output_bam}
samtools index {output_bam}

# =============================================================================
# STEP B: Pre-filter alignment statistics (includes chrM and spike-in reads)
# =============================================================================
echo '=== Step B: Pre-filter alignment statistics ==='
mkdir -p {qc_alignment_dir}
samtools flagstat {output_bam} > {qc_flagstat}
samtools idxstats {output_bam} > {qc_idxstat}

{spike_in_split_block}

{blacklist_filter_block}

# =============================================================================
# STEP E: Duplicate handling with Picard MarkDuplicates
# Input: fully filtered, tag-intact BAM from Step D.
# Output: the analysis-ready final.bam, used as the total-uniquely-aligned-reads
# denominator in spike-in normalization (step_4c), bigWig generation, and
# peak calling.
# =============================================================================
echo '=== Step E: Duplicate handling ==='
{picard_command}

# =============================================================================
# STEP F: Cleanup
# =============================================================================
echo '=== Step F: Cleanup ==='
{cleanup_block}

echo '=== Job complete ==='
echo "=== Resource usage summary ==="
bjobs -l $LSB_JOBID | grep -A10 "Resource usage"
"""

    with open(batch_file, 'w', encoding='utf-8') as fh:
        fh.write(batch_cmd)

    print(f"Submitting job: {job_name}")
    os.system(f"bsub < {batch_file}")
    time.sleep(1)


def main():
    """Load sample metadata and submit one job per sample."""
    if config.IS_PAIRED_END:
        ldsample.load_samples("paired_samples.txt")
        for sample_name, paired_sample in ldsample.SAMPLES_CTL.items():
            if not paired_sample:
                print(f"WARNING: No paired sample found for {sample_name}")
            process_bam(sample_name)
    else:
        ldsample.load_samples("samples.txt")
        for sample_name in ldsample.SAMPLES.keys():
            process_bam(sample_name)


if __name__ == "__main__":
    main()
