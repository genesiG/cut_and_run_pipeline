#!/usr/bin/env python3

"""
step_3_qc.py  –  Pre-peak QC: submit all independent bsub jobs
                  (no inter-step dependencies within this script)

Run this FIRST, then wait for all jobs to finish, then run:
    python3 Scripts/step_3_qc_report.py

Steps (all fire-and-forget bsub submissions):

  fastqc → FastQC on every raw/trimmed FASTQ file
  qc_3a  → Alignment stats from ORIGINALBAMDIR + Picard MarkDuplicates metrics
  qc_4   → Fragment size distributions                                      [R]
  qc_5   → Library complexity (Picard NRF metrics)
  qc_5b  → Preseq c_curve complexity curves on ORIGINALBAMDIR BAMs
  qc_6   → Fingerprint plots (deepTools plotFingerprint, grouped)
  qc_7   → Reproducibility (multiBamSummary + correlation + PCA)

None of these steps depend on outputs from other steps in this script.
All can be submitted in parallel.

Usage:
    python3 Scripts/step_3_qc.py              # run all steps
    python3 Scripts/step_3_qc.py --from qc_4  # start from qc_4
    python3 Scripts/step_3_qc.py --only qc_3a # run only one step
"""

import os
import time
import glob
import argparse

import config
from utils import ldsample

# ---------------------------------------------------------------------------
# Step registry
# Each entry: (step_key, script_path_or_None, is_rscript)
# None script_path = step is handled inline by a dedicated run_*() function
# ---------------------------------------------------------------------------
QC_STEPS = [
    ("fastqc", None,                               False),   # FastQC per sample
    ("qc_3a",  None,                               False),   # alignment stats + Picard
    ("qc_4",   "utils/qc_4_fragment_sizes.R",      True),
    ("qc_5",   "utils/qc_5_library_complexity.py", False),
    ("qc_5b",  None,                               False),   # preseq c_curve
    ("qc_6",   "utils/qc_6_fingerprint.py",        False),
    ("qc_7",   "utils/qc_7_reproducibility.py",    False),
]

SCRIPT_DIR    = config.CODEDIR
ALIGN_QC_DIR  = os.path.join(config.QCDIR1, "alignment")   # canonical alignment stats dir
PRESEQ_DIR    = os.path.join(config.LIB_COMPLEXITY_DIR, "preseq")

# FastQC output base: trimmed or untrimmed depending on config
if config.USE_TRIMMOMATIC:
    FASTQC_DIR = os.path.join(config.QCDIR1, "trimmed")
else:
    FASTQC_DIR = os.path.join(config.QCDIR1, "untrimmed")


# ===========================================================================
# fastqc — FastQC on raw/trimmed FASTQ files  (from step_0e_fastqc.py)
# ===========================================================================

def _submit_fastqc_job(sample_name: str) -> None:
    """Submit a bsub job that runs FastQC on a single sample's FASTQ file."""
    job_name = f"{sample_name}.fastqc"
    log_dir  = os.path.join(FASTQC_DIR, "log")
    os.makedirs(FASTQC_DIR, exist_ok=True)
    os.makedirs(log_dir,    exist_ok=True)

    if config.USE_TRIMMOMATIC:
        in_file = os.path.join(config.TRIMDIR, f"{sample_name}.fastq.gz")
    else:
        in_file = os.path.join(config.DATADIR, f"{sample_name}.fastq.gz")

    log_out = os.path.join(log_dir, f"{sample_name}.fastqc.log")
    log_err = os.path.join(log_dir, f"{sample_name}.fastqc.error")
    batch_f = os.path.join(log_dir, f"{sample_name}.fastqc.batch")

    batch_cmd = f"""\
#!/bin/bash
#BSUB -P {config.PROJECT_NAME}
#BSUB -J {job_name}
#BSUB -oo {log_out}
#BSUB -eo {log_err}
#BSUB -q {config.BSUB_QUEUE}
#BSUB -n 4
#BSUB -M 40960
#BSUB -R "rusage[mem=10240] span[hosts=1]"

source {config.HOMEDIR}/miniconda3/etc/profile.d/conda.sh
conda activate chipseq

fastqc -o {FASTQC_DIR} --extract {in_file}
rm -f {FASTQC_DIR}/*fastqc.zip
echo "FastQC done for: {sample_name}"
"""
    with open(batch_f, 'w', encoding='utf-8') as fh:
        fh.write(batch_cmd)
    print(f"  bsub < {batch_f}")
    os.system(f"bsub < {batch_f}")
    time.sleep(1)


def run_fastqc() -> None:
    """
    Submit one FastQC bsub job per sample.
    Reads FASTQ files from config.DATADIR (or TRIMDIR if USE_TRIMMOMATIC).
    Output: config.QCDIR1/untrimmed/ or config.QCDIR1/trimmed/
    These outputs are consumed by MultiQC in step_3_qc_report.py.
    """
    if config.IS_PAIRED_END:
        ldsample.load_samples("paired_samples.txt")
        all_samples = list(ldsample.SAMPLES_CTL.keys())
    else:
        ldsample.load_samples("samples.txt")
        all_samples = list(ldsample.SAMPLES.keys())

    if not all_samples:
        print("  WARNING: No samples found. Check samples file.")
        return

    print(f"  Found {len(all_samples)} sample(s). Submitting FastQC jobs...")
    for sn in all_samples:
        _submit_fastqc_job(sn)
    print(f"  [fastqc] {len(all_samples)} jobs submitted.")


# ===========================================================================
# qc_3a — Alignment stats (flagstat/idxstats from ORIGINALBAMDIR) + Picard
# ===========================================================================

def _submit_stats_job(sample_name: str, orig_bam: str,
                      flagstat_out: str, idxstat_out: str) -> None:
    """
    Submit a bsub job that coordinate-sorts the original (unfiltered) BAM,
    runs samtools flagstat + idxstats, then deletes the temporary sorted BAM.

    Using the ORIGINALBAMDIR BAM captures unmapped, chrM, and spike-in reads
    that are removed by quality-filtering in step_3_process_bam.py.
    """
    log_out    = os.path.join(ALIGN_QC_DIR, f"{sample_name}.stats.log")
    log_err    = os.path.join(ALIGN_QC_DIR, f"{sample_name}.stats.error")
    batch_f    = os.path.join(ALIGN_QC_DIR, f"{sample_name}.stats.batch")

    batch_cmd = f"""\
#!/bin/bash
#BSUB -P {config.PROJECT_NAME}
#BSUB -J {sample_name}.bamstats
#BSUB -oo {log_out}
#BSUB -eo {log_err}
#BSUB -q {config.BSUB_QUEUE}
#BSUB -n 4
#BSUB -M 16000
#BSUB -R "rusage[mem=4000] span[hosts=1]"

source {config.HOMEDIR}/miniconda3/etc/profile.d/conda.sh
conda activate chipseq

echo "Generating alignment stats for: {sample_name}"
echo "Input (original BAM): {orig_bam}"

# flagstat — total, mapped, unmapped counts
samtools flagstat {orig_bam} > {flagstat_out}
echo "flagstat written to: {flagstat_out}"

# idxstats — per-chromosome read counts (chrM, ecoli_chr, ...)
samtools idxstats {orig_bam} > {idxstat_out}
echo "idxstats written to: {idxstat_out}"

echo "Done: {sample_name}"
"""
    with open(batch_f, 'w', encoding='utf-8') as fh:
        fh.write(batch_cmd)
    os.system(f"bsub < {batch_f}")
    time.sleep(0.5)


def _submit_picard_job(sample_name: str, markdup_bam: str,
                       picard_out: str) -> None:
    """
    Submit a bsub job that:
    1. Adds @RG tags via picard AddOrReplaceReadGroups
       (required by picard MarkDuplicates; sambamba omits RG headers).
    2. Runs picard MarkDuplicates → metrics distinguishing optical from PCR dups.
    3. Deletes all temporary BAMs.
    """
    bam_rg     = os.path.join(config.PROCESSEDBAMDIR, f"{sample_name}.picard_rg.bam")
    picard_tmp = os.path.join(config.PROCESSEDBAMDIR, f"{sample_name}.picard_tmp.bam")
    log_out    = os.path.join(ALIGN_QC_DIR, f"{sample_name}.picard.log")
    log_err    = os.path.join(ALIGN_QC_DIR, f"{sample_name}.picard.error")
    batch_f    = os.path.join(ALIGN_QC_DIR, f"{sample_name}.picard.batch")

    batch_cmd = f"""\
#!/bin/bash
#BSUB -P {config.PROJECT_NAME}
#BSUB -J {sample_name}.picard
#BSUB -oo {log_out}
#BSUB -eo {log_err}
#BSUB -q {config.BSUB_QUEUE}
#BSUB -n 4
#BSUB -M 32000
#BSUB -R "rusage[mem=8000] span[hosts=1]"

source {config.HOMEDIR}/miniconda3/etc/profile.d/conda.sh
conda activate chipseq

echo "Running Picard MarkDuplicates for: {sample_name}"
echo "Input (markdup BAM): {markdup_bam}"

# Step 1: Add @RG read-group tags (required by Picard MarkDuplicates;
#         sambamba does not emit @RG headers).
picard AddOrReplaceReadGroups \\
    INPUT={markdup_bam} \\
    OUTPUT={bam_rg} \\
    RGID={sample_name} \\
    RGLB=lib1 \\
    RGPL=ILLUMINA \\
    RGPU=unit1 \\
    RGSM={sample_name} \\
    TMP_DIR={config.PROCESSEDBAMDIR}

# Step 2: Picard MarkDuplicates — optical vs PCR duplicate metrics
picard MarkDuplicates \\
    INPUT={bam_rg} \\
    OUTPUT={picard_tmp} \\
    METRICS_FILE={picard_out} \\
    OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \\
    REMOVE_DUPLICATES=false \\
    TMP_DIR={config.PROCESSEDBAMDIR}

rm -f {bam_rg} {picard_tmp}
echo "Picard metrics written to: {picard_out}"
"""
    with open(batch_f, 'w', encoding='utf-8') as fh:
        fh.write(batch_cmd)
    os.system(f"bsub < {batch_f}")
    time.sleep(0.5)


def run_qc_3a() -> None:
    """
    Generate pre-peak alignment QC statistics from the RAW bowtie2 output BAMs
    (config.ORIGINALBAMDIR/{sample}.bam, before any quality filtering).

    Using the original BAMs captures the full alignment picture:
      - Unmapped reads (removed by -F 4 in step_3_process_bam.py)
      - ChrM reads    (removed by grep -v chrM)
      - Spike-in reads (removed by BAM splitting)
      - Multi-mapped reads (removed by MAPQ filter)

    Phase 1 — Sort + flagstat + idxstats per sample (bsub; force-regenerate).
    Phase 2 — Picard MarkDuplicates per sample (bsub; skip if metrics exist).

    Output consumed by: qc_3_alignment_stats.R in step_3_qc_report.py.
    """
    os.makedirs(ALIGN_QC_DIR, exist_ok=True)

    markdup_suffix = f"qc.sort.markdup.mapq{config.MAPQ}.final.bam"
    bam_glob       = os.path.join(config.PROCESSEDBAMDIR, f"*.{markdup_suffix}")
    bam_files      = sorted(glob.glob(bam_glob))

    if not bam_files:
        print(f"  WARNING: No markdup BAMs found matching: {bam_glob}")
        return

    print(f"  Found {len(bam_files)} markdup BAM(s).")
    stats_submitted  = []
    picard_submitted = []

    for markdup_bam in bam_files:
        base        = os.path.basename(markdup_bam)
        sample_name = base[: -len(markdup_suffix) - 1]

        # Phase 1: Alignment stats from ORIGINALBAMDIR (always regenerate)
        orig_bam     = os.path.join(config.ORIGINALBAMDIR, f"{sample_name}.bam")
        flagstat_out = os.path.join(ALIGN_QC_DIR, f"{sample_name}.flagstat.txt")
        idxstat_out  = os.path.join(ALIGN_QC_DIR, f"{sample_name}.idxstat.txt")

        if not os.path.exists(orig_bam):
            print(f"  WARNING: Original BAM not found: {orig_bam}")
        else:
            print(f"  [qc_3a] Submitting stats job: {sample_name}")
            _submit_stats_job(sample_name, orig_bam, flagstat_out, idxstat_out)
            stats_submitted.append(sample_name)

        # Phase 2: Picard MarkDuplicates metrics (skip if already present)
        picard_out = os.path.join(ALIGN_QC_DIR, f"{sample_name}.picard_metrics.txt")
        if os.path.exists(picard_out):
            print(f"  [qc_3a] SKIP Picard (exists): {sample_name}")
        else:
            print(f"  [qc_3a] Submitting Picard job: {sample_name}")
            _submit_picard_job(sample_name, markdup_bam, picard_out)
            picard_submitted.append(sample_name)

    total = len(stats_submitted) + len(picard_submitted)
    print(f"\n  [qc_3a] {len(stats_submitted)} stats + "
          f"{len(picard_submitted)} Picard jobs submitted ({total} total).")
    print("  NOTE: Wait for ALL jobs to finish before running step_3_qc_report.py.")
    print("        Monitor with: bjobs -u $USER")


# ===========================================================================
# qc_5b — preseq c_curve on ORIGINALBAMDIR BAMs
# ===========================================================================

def _submit_preseq_job(sample_name: str, orig_bam: str,
                       out_file: str) -> None:
    """Submit a bsub job that runs `preseq c_curve` on a raw original BAM."""
    log_out = os.path.join(PRESEQ_DIR, f"{sample_name}.preseq.log")
    log_err = os.path.join(PRESEQ_DIR, f"{sample_name}.preseq.error")
    batch_f = os.path.join(PRESEQ_DIR, f"{sample_name}.preseq.batch")

    pe_flag = "-P" if config.IS_PAIRED_END else ""

    batch_cmd = f"""\
#!/bin/bash
#BSUB -P {config.PROJECT_NAME}
#BSUB -J {sample_name}.preseq
#BSUB -oo {log_out}
#BSUB -eo {log_err}
#BSUB -q {config.BSUB_QUEUE}
#BSUB -n 1
#BSUB -M 16000
#BSUB -R "rusage[mem=16000] span[hosts=1]"

source {config.HOMEDIR}/miniconda3/etc/profile.d/conda.sh
conda activate chipseq

echo "Running preseq c_curve for: {sample_name}"
echo "Input BAM: {orig_bam}"

# Compute a sensible step size: total_reads / 20, clamped to [10000, 1000000].
# This ensures we always get ~20 data points on the curve even for small IgG
# libraries where the default 1e6 step would exceed the library size entirely.
TOTAL_READS=$(samtools flagstat {orig_bam} | awk 'NR==1{{print $1}}')
STEP_SIZE=$(python3 -c "
total = int('$TOTAL_READS') if '$TOTAL_READS'.isdigit() else 1000000
step  = max(10000, min(1000000, total // 20))
print(step)
")
echo "Library size: $TOTAL_READS reads  →  step size: $STEP_SIZE"

preseq c_curve \\
    -B \\
    {pe_flag} \\
    -v \\
    -s $STEP_SIZE \\
    -o {out_file} \\
    {orig_bam}

echo "preseq c_curve written to: {out_file}"
"""
    with open(batch_f, 'w', encoding='utf-8') as fh:
        fh.write(batch_cmd)
    os.system(f"bsub < {batch_f}")
    time.sleep(0.5)


def run_qc_5b() -> None:
    """
    Submit one `preseq c_curve` bsub job per sample using the raw bowtie2
    output BAMs in config.ORIGINALBAMDIR (before quality filtering or
    deduplication).  Captures the true empirical duplicate distribution.

    Output: config.LIB_COMPLEXITY_DIR/preseq/{sample}.preseq_ccurve.txt
    Existing files are skipped (idempotent).
    """
    os.makedirs(PRESEQ_DIR, exist_ok=True)

    orig_bams = sorted(glob.glob(os.path.join(config.ORIGINALBAMDIR, "*.bam")))
    if not orig_bams:
        print(f"  WARNING: No BAMs found in ORIGINALBAMDIR: {config.ORIGINALBAMDIR}")
        return

    print(f"  Found {len(orig_bams)} original BAM(s).")
    submitted = []

    for orig_bam in orig_bams:
        sample_name = os.path.basename(orig_bam).replace(".bam", "")
        out_file    = os.path.join(PRESEQ_DIR,
                                   f"{sample_name}.preseq_ccurve.txt")
        if os.path.exists(out_file):
            print(f"  [qc_5b] SKIP (exists): {sample_name}")
            continue
        print(f"  [qc_5b] Submitting preseq job: {sample_name}")
        _submit_preseq_job(sample_name, orig_bam, out_file)
        submitted.append(sample_name)

    print(f"\n  [qc_5b] {len(submitted)} preseq c_curve jobs submitted.")
    if submitted:
        print("  Output directory:", PRESEQ_DIR)
        print("  Monitor with: bjobs -u $USER")


# ===========================================================================
# Generic step runners
# ===========================================================================

def run_python_step(script_path: str) -> None:
    """
    Import and run the main() function of a Python QC utility script.
    Executes in the same interpreter; catches SystemExit gracefully.
    """
    import importlib.util
    spec   = importlib.util.spec_from_file_location("_qc_module", script_path)
    module = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(module)
    except SystemExit:
        pass   # some scripts call sys.exit() — ignore

    if hasattr(module, "main"):
        module.main()
    else:
        print(f"  WARNING: {os.path.basename(script_path)} has no main().")


def run_rscript_step(script_path: str, step_key: str) -> None:
    """Submit an R script as a bsub batch job."""
    script_name = os.path.basename(script_path)

    resource_map = {
        "qc_4": (4, 64000),
    }
    num_cores, max_mem = resource_map.get(step_key, (2, 16000))
    mem_per_core = int(max_mem / num_cores)

    log_dir      = os.path.join(config.QCDIR1, "log")
    os.makedirs(log_dir, exist_ok=True)

    job_name     = f"{config.PROJECT_NAME}.{step_key}"
    log_out_file = os.path.join(log_dir, f"{step_key}.log")
    log_err_file = os.path.join(log_dir, f"{step_key}.error")
    batch_file   = os.path.join(log_dir, f"{step_key}.batch")

    batch_cmd = f"""\
#!/bin/bash
#BSUB -P {config.PROJECT_NAME}
#BSUB -J {job_name}
#BSUB -oo {log_out_file}
#BSUB -eo {log_err_file}
#BSUB -q {config.BSUB_QUEUE}
#BSUB -n {num_cores}
#BSUB -M {max_mem}
#BSUB -R "rusage[mem={mem_per_core}] span[hosts=1]"

source {config.HOMEDIR}/miniconda3/etc/profile.d/conda.sh
conda activate chipseq

echo "Running {script_name}"
Rscript {script_path}

echo "=== Resource usage summary ==="
bjobs -l $LSB_JOBID | grep -A10 "Resource usage"
"""
    with open(batch_file, 'w', encoding='utf-8') as fh:
        fh.write(batch_cmd)

    print(f"  Submitting job: {job_name}")
    os.system(f"bsub < {batch_file}")
    time.sleep(1)


# ===========================================================================
# CLI
# ===========================================================================

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "step_3_qc — Submit all independent pre-peak QC bsub jobs.\n"
            "Run step_3_qc_report.py after all jobs finish."
        )
    )
    grp = parser.add_mutually_exclusive_group()
    grp.add_argument(
        "--from", dest="from_step", default=None,
        choices=[s[0] for s in QC_STEPS],
        help="Start from this step (inclusive).",
    )
    grp.add_argument(
        "--only", dest="only_step", default=None,
        choices=[s[0] for s in QC_STEPS],
        help="Run only this step.",
    )
    return parser.parse_args()


def main() -> None:
    """Submit all independent QC bsub jobs in the configured order."""
    args  = parse_args()
    steps = QC_STEPS[:]

    if args.only_step:
        steps = [s for s in steps if s[0] == args.only_step]
    elif args.from_step:
        keys  = [s[0] for s in steps]
        start = keys.index(args.from_step)
        steps = steps[start:]

    print("=" * 60)
    print(f"  {config.PROJECT_NAME} — step_3_qc.py")
    print(f"  Steps: {', '.join(s[0] for s in steps)}")
    print("=" * 60)

    # Inline step dispatch table
    inline_runners = {
        "fastqc": run_fastqc,
        "qc_3a":  run_qc_3a,
        "qc_5b":  run_qc_5b,
    }

    for step_key, script_name, is_rscript in steps:
        print(f"\n{'='*60}")
        label = script_name or f"({step_key})"
        print(f"[{step_key}] Running: {label}")
        print("=" * 60)

        t0 = time.time()
        try:
            if step_key in inline_runners:
                inline_runners[step_key]()
            else:
                script_path = os.path.join(SCRIPT_DIR, script_name)
                if not os.path.exists(script_path):
                    print(f"  [{step_key}] SKIPPED — script not found: {script_path}")
                    continue
                if is_rscript:
                    run_rscript_step(script_path, step_key)
                else:
                    run_python_step(script_path)
        except Exception as exc:   # noqa: BLE001
            print(f"  [{step_key}] ERROR: {exc}")

        elapsed = time.time() - t0
        print(f"[{step_key}] Done in {elapsed:.1f}s")

    print("\n" + "=" * 60)
    print("  All step_3_qc jobs submitted.")
    print("  Wait for all bsub jobs to finish, then run:")
    print("    python3 Scripts/step_3_qc_report.py")
    print("  Monitor: bjobs -u $USER")
    print("=" * 60)


if __name__ == "__main__":
    main()
