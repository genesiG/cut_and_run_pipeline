# CUT&RUN / ChIP-seq Multi-Omic Analysis Pipeline

A config-driven, HPC-batch-scheduled pipeline for CUT&RUN, CUT&Tag, ChIP-seq, ATAC-seq, and
RNA-seq analysis, built around chromatin profiling of histone post-translational modifications
(H3K27me2/me3, H2AK119ub) and a Polycomb reader protein. It takes raw paired-end FASTQ through
alignment, spike-in/TMM normalization, peak calling, differential binding, and multi-omic
integration with RNA-seq, and it includes an agentic orchestration layer that runs the pipeline
end-to-end, submits and polls HPC batch jobs, triages failures, and writes methods sections from
the artifacts of a completed run.

This repository is a sanitized snapshot of a working research pipeline. Absolute paths and the
originating HPC username have been replaced with placeholders, and directories holding raw
sequencing data (FASTQ/BAM/BigWig) and unpublished project-specific metadata have been excluded.

## Architecture

```
Scripts/
├── config.py                  # Single source of truth: paths, genome build, normalization,
│                               # antibody/group naming, MACS/csaw parameters — every step
│                               # reads from here instead of hardcoding values.
├── step_0*.py                 # Setup: rename raw FASTQ, build sample sheets, FastQC
├── step_1_trimmomatic.py      # Adapter/quality trimming
├── step_2_align_bowtie2.py    # Alignment (spike-in-aware genome index)
├── step_3_process_bam.py      # QC filter → spike-in split → blacklist filter → dedup
├── step_4*.py                 # Normalization: ChIPseqSpikeInFree / TMM (csaw) / spike-in
├── step_5*.py                 # BigWig coverage track generation
├── step_6*.py                 # Peak calling: MACS3, SEACR, csaw sliding-window
├── step_7*–9*                 # Peak classification, annotation, differential binding,
│                               # deepTools heatmaps/profiles, genome track rendering
├── step_10–14                 # Multi-omic integration: chromatin state × RNA-seq
│                               # expression stratification, reader-specificity analysis
├── run_R_steps.sh             # Orchestrates the R-based steps in one conda environment
└── utils/                     # ~50 shared helpers (R + Python): csaw parameterization,
                                # QC (FRIP, fingerprint, library complexity, reproducibility),
                                # normalization scaling factors, plotting, peak I/O

Skills/, Specialists/, AGENTS.md   # Agentic orchestration layer (see below)
```

**Design principles**
- Paths are derived from `config.py`.
- Every step is per-sample and submits its own HPC batch job (`bsub`), writing its own
  `.batch`/`.log`/`.error` files so a failed sample can be re-run in isolation.
- `utils/ldsample.py` centralizes sample handling: it reads `Metadata/samples.txt`, resolves each
  sample name to its FASTQ/BAM source file, and tracks sample/input-control pairings.

## Agentic orchestration layer

`AGENTS.md`, `Specialists/`, and `Skills/` define a multi-agent system (built for Claude Code /
Google Antigravity) that runs this pipeline autonomously on an LSF HPC cluster:

- **Orchestrator** — discovers pending steps, submits batch jobs, polls with exponential backoff,
  and hands off to QC after every heavy step.
- **QC** — inspects logs and outputs after each step, issues a PASS/WARN/FAIL verdict against
  defined thresholds (alignment rate, empty outputs, etc.), and routes failures to Engineer.
- **Engineer** — diagnoses failures from `.error` logs, patches the offending script following the
  codebase's script-authoring conventions, and hands the fix back for resubmission.
- **Analyst** — generates and interprets QC figures, flagging anomalous samples by name.

A dedicated skill (`skill_write_methods.md`) reconstructs a publication-ready Materials & Methods
section — with correct tool versions, parameters, and citations — directly from the batch scripts,
logs, and conda environment of a completed run, with an explicit zero-hallucination mandate: every
parameter must be verified against what was actually executed, not assumed from a script default.

## Usage

1. Create the conda environment: `mamba env create -f environment.yml`
2. Edit `Scripts/config.py`: set `USERNAME`, `PROJECT_NAME`, `SPECIES`/reference genome paths,
   and experiment type.
3. Populate `Metadata/samples.txt` (sample name ↔ input-control pairing) and place raw FASTQ under
   `Original_Data/fastq_files/`.
4. Run steps in order (`Scripts/step_0a...` → `step_14...`), or let the orchestrator agent drive
   the run — see `Skills/skill_orchestrator.md`.

## License

MIT — see [LICENSE](LICENSE).
