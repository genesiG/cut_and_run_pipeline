---
name: skill_write_methods
description: >
  Guides an AI agent step-by-step through writing a publication-quality
  Materials and Methods (M&M) section for any genomics pipeline
  (ATAC-seq, ChIP-seq, CUT&RUN/CUT&TAG, RNA-seq, CRISPR screen, etc.).
  The agent reads batch job scripts and the project config.py to reconstruct
  exact parameters, resolves helper/utility scripts for analytical details,
  inspects the active conda environment(s) for precise tool versioning, and
  produces a Microsoft Word-compatible document (.docx preferred) or Markdown (.md fallback) with a properly formatted Reference section.
  Trigger this skill whenever the user asks to "write methods", "write M&M",
  "describe the pipeline", or "write a methods section" for a project.
---

# SKILL: skill_write_methods.md

## Purpose

Reconstruct a complete, publication-ready **Materials and Methods** section
from the artefacts of a completed genomics pipeline — specifically the LSF
batch scripts (`.batch`), pipeline configuration (`config.py`), helper/utility
R or Python scripts called within batch jobs, and conda environment exports.

The output **must** be formatted in a **Microsoft Word `.docx` or other MS Word-compatible format (preferential)**. If `.docx` or Word-compatible generation tools (`python-docx`, `pandoc`) are unavailable across all active environments, fall back cleanly to a structured **Markdown (`.md`)** file. Never produce plain text.

---

## 0 — Pre-flight Checks

Before reading any files, run the following **once** in a persistent terminal:

```bash
hostname      # must be a compute node (nodeXXX); if login node, run: bsub -Is bash
whoami        # must be <your_username>
```

Never read files or run conda commands from a login node. The compute node
constraint from `AGENTS.md` applies here too.

---

## 1 — Discovery Phase (Always perform in this order)

### 1.1  Read `config.py`

Read the project's `Scripts/config.py` (or the equivalent `config.py` at the
project root). Extract and record **every** variable you will need to
parameterize the M&M text:

| Variable to extract | Used in M&M subsection |
|---|---|
| `PROJECT_NAME` | Header / project context |
| `EXPERIMENT` | Determines which subsections to write |
| `SPECIES` / `GNM_IDX` | Reference genome subsection |
| `ALIGNER` | Alignment subsection |
| `USE_TRIMMOMATIC` | Trimming subsection (skip if False) |
| `IS_PAIRED_END` | Read strategy sentence |
| `FORWARD_PATTERN` / `REVERSE_PATTERN` | FASTQ naming convention |
| `REMOVE_DUPLICATES` | BAM processing subsection |
| `REMOVE_BLACKLIST` | BAM processing subsection |
| `MAPQ` | Filtering subsection |
| `IS_RAW_DATA_FASTQ` | Whether to include a trimming/alignment section |
| `BIN_SIZE` | BigWig/coverage track subsection |
| `NORMALIZE_USING` | Normalization subsection |
| `TARGETS` | Normalization grouping logic |
| `MACS_Q_VALUE`, `MACS_FE_CUTOFF`, `MACS_MAX_GAP`, `MACS_MIN_LENGTH`, `MACS_CALL_SUMMITS`, `MACS_CALL_BROAD_PEAKS`, `GENOME_SIZE_FOR_MACS` | Peak calling subsection |
| `REPLICATES` | Peak calling groupings |
| Any spike-in or reference-genome normalization flags | Normalization subsection |

> **TIP:** If `config.py` uses conditionals (e.g. `if SPECIES == "t2t":`),
> resolve the actual paths and values that were active during the run — i.e.
> the branch that matches the recorded `SPECIES` value.

---

### 1.2  Identify which output directories exist

List the standard output directories. Only write M&M subsections for steps
whose output directories exist and are non-empty:

```
Importable_Data/trimmed/                   → Trimming step
Analysis_Data/alignment/                   → Alignment step
Importable_Data/bam/processed/             → BAM processing step
Importable_Data/bam/original/              → Initial SAM→BAM conversion
Analysis_Data/normalization/               → Normalization step
Analysis_Data/bigwig/                      → Coverage track generation
Analysis_Data/peaks/                       → Peak calling
Analysis_Data/peaks/annotation/            → Peak annotation
Analysis_Data/deeptools/                   → deepTools heatmaps/profiles
Analysis_Data/qc/                          → QC (FastQC / MultiQC)
```

For RNA-seq pipelines, also check:
```
Analysis_Data/counts/                      → Feature counting (featureCounts)
Analysis_Data/dge/                         → Differential gene expression
```

For CRISPR screens, also check:
```
Analysis_Data/counts/                      → MAGeCK count
Analysis_Data/mageck/                      → MAGeCK test / MLE
```

---

### 1.3  Read ONE batch file per step

For **each** non-empty output directory identified above, locate the
corresponding **`log/`** subdirectory and read **exactly one** `.batch` file
(pick the first alphabetically, or any representative sample).

**Standard batch file locations and naming patterns:**

| Step | Log directory | Batch file pattern |
|---|---|---|
| Trimming | `Importable_Data/trimmed/log/` | `<sample>.batch` |
| Alignment | `Analysis_Data/alignment/log/` | `<sample>.batch` |
| BAM processing | `Importable_Data/bam/log/` | `<sample>.batch` |
| Normalization | `Analysis_Data/normalization/log/` | `tmm_*.batch` or `normalization.batch` |
| BigWig | `Analysis_Data/bigwig/log/` | `<sample>.<binsize>.bigwig.batch` |
| Peak calling | `Analysis_Data/peaks/log/` | `<sample_or_group>.batch` |
| Peak annotation | `Analysis_Data/peaks/annotation/log/` | `*.batch` |
| DGE | `Analysis_Data/dge/log/` | `*.batch` |
| MAGeCK count | `Analysis_Data/counts/log/` | `*.batch` |
| MAGeCK test | `Analysis_Data/mageck/log/` | `*.batch` |
| deepTools | `Analysis_Data/deeptools/log/` | `*.batch` |

From each batch file, extract:
- The **tool name and command** (e.g., `trimmomatic`, `bowtie2`, `samtools`, `bamCoverage`)
- **All flags and parameters** passed to the tool
- The **input / output file paths** (to understand the workflow topology)
- Whether the job calls any **helper / utility scripts** (see §1.4)

---

### 1.4  Resolve helper and utility scripts

A batch file is "shallow" if it delegates work to a secondary R or Python
script (not named `step_N_*.py`). Common patterns:

```bash
# Examples of helper script calls inside batch files:
Rscript Scripts/scaleFactors.R ...
Rscript Scripts/getPEsizes.R
Rscript Scripts/utils/csaw_parameters.R ...
Rscript Scripts/utils/callpeaks_csaw.R ...
python3 Scripts/utils/count_matrix.py ...
```

**Rule:** Whenever a batch file calls a helper/utility script, you **MUST**
read the full content of that helper script before writing the M&M text for
that step. The helper script contains the true analytical parameters (window
sizes, statistical thresholds, normalization formulas, etc.) that must appear
in the methods.

Steps:
1. Identify the script path from the batch file command line.
2. Resolve relative paths against the project `WORKDIR`.
3. Read the script using `view_file`.
4. **CRITICAL MANDATE — Verify Actual Executed Parameters:** You MUST cross-reference the helper script against the exact command-line arguments passed in the `.batch` file AND check the corresponding `.error`/`.log`/`.out` files in the log directory. Helper scripts (such as `scaleFactors.R` or `csaw_parameters.R`) often have conditional branches (e.g. `if (EFFICIENCY_BIAS == TRUE) width=150` vs `if (EFFICIENCY_BIAS == FALSE) width=10000`) or dynamically set flags (e.g. setting `dedup=FALSE` when reading `.rmdup.` BAMs). Never report default variables from an R/Python script without verifying what was actually printed and executed in the `.batch` and `.error`/`.log` files!
5. Extract all verified values relevant to the method description (e.g., bin width, efficiency bias setting, deduplication flag, maximum fragment size from logs, statistical thresholds).

---

### 1.5  Inspect conda environments for tool versions

The batch files will declare which conda environment to activate
(e.g., `conda activate chipseq`). For **each unique environment name**
encountered across all batch files, export its package list:

```bash
conda run -n <ENV_NAME> conda list --export 2>/dev/null | head -200
# OR
conda env export -n <ENV_NAME> 2>/dev/null | head -200
```

From the output, extract the version string for **every tool** mentioned in
the batch files and helper scripts. Record them in a local table:

| Tool | Package name in conda | Version found |
|---|---|---|
| Trimmomatic | `trimmomatic` | e.g., `0.39` |
| Bowtie2 | `bowtie2` | e.g., `2.5.3` |
| SAMtools | `samtools` | e.g., `1.20` |
| Sambamba | `sambamba` | e.g., `1.0.1` |
| BEDTools | `bedtools` | e.g., `2.31.1` |
| deepTools | `deeptools` | e.g., `3.5.5` |
| MACS3 / MACS2 | `macs3` or `macs2` | e.g., `3.0.1` |
| HISAT2 | `hisat2` | e.g., `2.2.1` |
| featureCounts | `subread` | e.g., `2.0.6` |
| MAGeCK | `mageck` | e.g., `0.5.9.5` |
| R | `r-base` | e.g., `4.4.1` |
| DESeq2 | `bioconductor-deseq2` | e.g., `1.44.0` |
| edgeR | `bioconductor-edger` | e.g., `4.2.1` |
| csaw | `bioconductor-csaw` | e.g., `1.38.0` |
| limma | `bioconductor-limma` | e.g., `3.60.4` |

> **IMPORTANT:** If a version cannot be found in the conda export, note it as
> "version not pinned" and use the version found in the batch file's shebang
> or path (e.g., `/opt/software/Trimmomatic/0.32/`).

---

## 2 — Writing the M&M Text

Write one subsection per pipeline step. Follow the experiment-type templates
below. Use only the parameters you **actually observed** in the batch files and
helper scripts — do not invent or infer values not present in the code.

### 2.1 — General Writing Rules

0. DO NOT mention that a step was performed "using a custom script"; instead, describe what was done in the step succintly and without mentioning the scripts used.
1. **Always include tool version** in parentheses on first mention:
   e.g., "Trimmomatic (v0.39)".
2. **Always include key parameters** as they appear in the command line or
   helper script, written in a readable prose style.
3. **Reference genome/annotation:** Spell out the full genome name and
   assembly version (e.g., "T2T-CHM13v2.0 human reference genome" or
   "GRCm39 mouse reference genome"), not just the short code.
4. **Paired-end vs. single-end:** State this explicitly in the alignment
   subsection. Derive it from `config.IS_PAIRED_END`.
5. **Blacklist:** If `config.REMOVE_BLACKLIST == True`, state which blacklist
   was used (derive from `config.BLACKLIST` path).
6. **Duplicates:** State whether duplicates were "marked and removed" or
   "marked only" based on `config.REMOVE_DUPLICATES`.
7. **Statistical thresholds and cutoffs** from helper scripts must be stated
   numerically (FDR < 0.05, minimum fold-change, etc.).
8. Write in **past tense**, third person (or first person plural — match the
   user's stated preference if given).
9. Do NOT describe steps that were skipped (e.g., if `USE_TRIMMOMATIC = False`,
   omit the trimming subsection entirely).
10. **ABSOLUTE FIDELITY & ZERO-HALLUCINATION MANDATE:** You must copy parameter values (e.g. `dedup=FALSE`, `EFFICIENCY_BIAS=FALSE`, `max.frag=1999`, `bin width=10000`, `NexteraPE-PE.fa`, `--local --very-sensitive-local -N 1 --dovetail`) EXACTLY as executed and logged in the `.batch` and `.error`/`.log` files. NEVER copy placeholder values from this document's templates (e.g. `TruSeq3-PE.fa`, `--end-to-end`, `SLIDINGWINDOW:4:15`, `MINLEN:36`) or default variables from R/Python scripts if the `.batch`/`.error` logs show different executed values! Precision and truthfulness are non-negotiable.

---

### 2.2 — Subsection Templates by Experiment Type

Use the relevant block(s) for the experiment. Sections are modular — compose
only what applies.

---

#### BLOCK A — Read Alignment and Processing (all experiment types except pre-aligned BAM input)

> **Trigger:** `config.IS_RAW_DATA_FASTQ == True` AND `config.USE_TRIMMOMATIC == True`
> **Trigger:** `Importable_Data/bam/log/` exists and contains batch files.

**Template:**

> ## Read Alignment and Processing
> Raw paired-end sequencing reads were processed to remove sequencing adapters and low-quality bases with Trimmomatic (v`<VERSION>`; [citation]) using the ILLUMINACLIP option with the `<ADAPTER_FILE>` adapter file with the following parameters: `<ADAPTER_FILE>:<PARAMETERS>:<TRUE> LEADING:<LEADING> TRAILING:<TRAILING> SLIDINGWINDOW:<WINDOW_SIZE>:<REQUIRED_QUALITY> MINLEN:<MINLEN>`.
> Trimmed reads were aligned to the <human/mouse/species> reference genome (index: `<GNM_IDX basename>`) using Bowtie2 (v`<VERSION>`; [citation]) in `<--local / --end-to-end>` mode with the `--very-sensitive[-local]` preset. Discordant and mixed-orientation alignments were excluded (`--no-discordant --no-mixed`), while dovetailing read pairs were permitted (`--dovetail`).
> Alignment output was converted from SAM to BAM format and then filtered using SAMtools (v`<VERSION>`; [citation]) view (<PARAMETERS> -q <mapq> -f <FLAG_PARAMS> -F <FLAG_PARAMS>) to retain only properly paired reads with a minimum mapping quality score of <mapq>. Filtered BAM files were sorted and indexed using `samtools sort` and `samtools index`. [If REMOVE_DUPLICATES == True:] PCR and optical duplicate reads were identified and removed using Sambamba (<VERSION>; [citation]) `markdup -r` [or Picard (v`<VERSION>`; [citation]) `MarkDuplicates REMOVE_DUPLICATES=true`]. [If REMOVE_BLACKLIST == True — ATAC-seq / ChIP-seq / CUT&RUN path through BED:] To remove reads overlapping known problematic genomic regions, duplicate-removed BAM files were query-name sorted using `samtools sort <SORT_PARAMETERS>` and converted to paired-end BED (BEDPE) format using BEDTools (<VERSION>; [citation]) `bamtobed -bedpe`. Reads overlapping blacklisted regions were excluded using `bedtools intersect -v`, converted back to coordinate-sorted BAM format using `bedtools bedpetobam`, and indexed with `samtools index`. [If blacklist was removed while staying BAM-native:] Blacklisted regions (`<BLACKLIST basename>`) were excluded using `bedtools intersect -v -abam`.

>[NOTE]:Parse the batch file carefully; the exact combination of SAMtools flags
determines this paragraph:

**Standard flags table for SAMtools view:**

| SAMtools FLAG | Meaning |
|---|---|
| `-f 1` | Read is paired |
| `-f 2` | Read mapped in proper pair |
| `-F 4` | Exclude unmapped reads |
| `-F 8` | Exclude reads with unmapped mate |
| `-F 256` | Exclude non-primary alignments |
| `-F 512` | Exclude QC-failed reads |
| `-F 1024` | Exclude PCR/optical duplicates |
| `-F 2048` | Exclude supplementary alignments |
| `-q <N>` | Minimum mapping quality score |

---

**A2 — HISAT2 (RNA-seq)**

> **Trigger:** `config.ALIGNER == "hisat2"` or `EXPERIMENT == "rnaseq"`

**Template:**

> Trimmed reads were aligned to the <human/mouse/species> reference genome (annotation: `<GTF basename>`) using HISAT2 (v`<VERSION>`; [citation]) with parameters `<PARAMETERS>`. Alignment output was converted from SAM to BAM format and then filtered using SAMtools (v`<VERSION>`; [citation]) view (`-h -b -S/<PARAMETERS>` -q <mapq> -f <FLAG_PARAMS> -F <FLAG_PARAMS>) to retain only properly paired reads with a minimum mapping quality score of <mapq>. Filtered BAM files were sorted and indexed using `samtools sort` and `samtools index`.
> [Include any other non-default flags observed in the batch file.]

---

#### BLOCK B — Normalization

**B1 — TMM normalization (ATAC-seq, ChIP-seq without spike-in)**

> **Trigger:** `config.NORMALIZE_USING == "tmm"` and `Analysis_Data/normalization/` exists.

Read the helper scripts called (`getPESizes.R`, `scaleFactors.R`), AND check the `.batch` command-line arguments and `.error`/`.log` files to extract the exact numerical values (`maxFrag` printed in error log, e.g. `1999`), `BIN_WIDTH` (`10000` bp / 10 kb), `EFFICIENCY_BIAS` status (`FALSE`), `dedup` setting passed to `readParam()` (`FALSE` when already `rmdup`), and `weighted` TMM parameter (`FALSE`).

**Template:**

> To normalize chromatin accessibility signals across samples, Trimmed Mean of M-values (TMM) scaling factors were calculated using  edgeR (<VERSION>; [citation]) and csaw (<VERSION>; [citation]). Briefly, fragment size distributions were computed via the `getPESizes()` function and reads were counted over large bins (10 kb) to estimate background noise across samples (`windowCounts` with `param=readParam(minq=<NUMERICAL mapq>, pe="both"[or "none" if single-end], max.frag=<NUMERICAL maxFrag>, dedup=<EXACT DEDUP FLAG FROM LOG, e.g. FALSE>)`). [If EFFICIENCY_BIAS=FALSE:] To normalize for composition bias, unweighted TMM normalization (`normFactors(weighted=FALSE)`) was applied across the background bins and scaling factors were computed as `SizeFactor = normfacs * libSizes / 1e6`.

> Genome coverage tracks were generated in BigWig format using `bamCoverage` from the deepTools suite (<VERSION>; [citation]). Briefly, paired-end reads were centered on the fragment midpoint (`--centerReads`) and coverage was computed at <NUMERICAL_BIN_SIZE_FROM_BATCH_SCRIPT>-bp resolution (`--binSize <NUMERICAL_BIN_SIZE_FROM_BATCH_SCRIPT`) excluding bins with zero read coverage (`--skipNonCoveredRegions`). Each sample was then scaled by the inverse of its corresponding TMM size factor (`--scaleFactor 1/<SizeFactor>`).

**B2 — Spike-in normalization (CUT&RUN, ChIP-seq with exogenous spike)**

> **Trigger:** spike-in step batch file exists.

Read the spike-in normalization script (`step_4c_spikein_normalization.py` or
equivalent) for the exact formula and describe it here. 
>[NOTE:] Write the formula using LaTeX or similar algebraic notation.

**Template:**

> Spike-in normalization was performed by mapping reads to the `<SPIKE_IN_GENOME>` reference genome. A scaling factor for each sample was computed as:

> `<FORMULA extracted from script>`

> Genome coverage tracks were generated in BigWig format using `bamCoverage` from the deepTools suite (<VERSION>; [citation]). Briefly, paired-end reads were centered on the fragment midpoint (`--centerReads`) and coverage was computed at <NUMERICAL_BIN_SIZE_FROM_BATCH_SCRIPT>-bp resolution (`--binSize <NUMERICAL_BIN_SIZE_FROM_BATCH_SCRIPT`) excluding bins with zero read coverage (`--skipNonCoveredRegions`). Each sample was then scaled by its corresponding scaling factor (`--scaleFactor <SF>`).

**B3 — ChIPseqSpikeInFree**

> **Trigger:** `ChIPseqSpikeInFree` appears in a batch file or log or setting as config.py.

**Template:**

> Normalization was performed using ChIPseqSpikeInFree (v`<VERSION>`; [citation]). Coverage was computed across non-overlapping `<BIN_WIDTH from batch file>` bp bins genome-wide [if cutoff == 1.2 and max_turns == 0.99] with default parameters [or with parameters `cutoff = <cutoff>` and `max_turns = <max_turns>`].

> Genome coverage tracks were generated in BigWig format using `bamCoverage` from the deepTools suite (<VERSION>; [citation]). Briefly, paired-end reads were centered on the fragment midpoint (`--centerReads`) and coverage was computed at <NUMERICAL_BIN_SIZE_FROM_BATCH_SCRIPT>-bp resolution (`--binSize <NUMERICAL_BIN_SIZE_FROM_BATCH_SCRIPT`) excluding bins with zero read coverage (`--skipNonCoveredRegions`). Each sample was then scaled by the inverse of its corresponding scaling factor (`--scaleFactor 1/<SF>`).
---

#### BLOCK C — Peak Calling

**C1 — ATAC-seq**

> Chromatin accessible regions were identified using MACS3 (v`<VERSION>`; [citation]) in `--nomodel` mode with shift and extension parameters `--shift -100 --extsize 200`, a q-value cutoff of `<MACS_Q_VALUE>` (`-q <MACS_Q_VALUE>`), and a minimum peak length of `<MACS_MIN_LENGTH>` bp (`--min-length <MACS_MIN_LENGTH>`) with a maximum allowed gap of `<MACS_MAX_GAP>` bp (`--max-gap <MACS_MAX_GAP>`). The effective genome size (`-g`) was set to `<GENOME_SIZE_FOR_MACS>`.

**C2 — ChIP-seq (narrow marks)**

> Peaks were called using MACS3 (v`<VERSION>`; [citation]) with a q-value cutoff of `<MACS_Q_VALUE>` (`-q <MACS_Q_VALUE>`) and a fold-enrichment cutoff of `<MACS_FE_CUTOFF>`. [If paired-end:] The `--format BAMPE` flag was used to leverage paired-end fragment size information alongside effective genome size (`-g <GENOME_SIZE_FOR_MACS>`).

**C3 — ChIP-seq / CUT&RUN (broad marks)**

> Broad chromatin domains were identified using MACS3 (v`<VERSION>`; [citation]) in broad peak mode (`--broad --broad-cutoff <MACS_Q_VALUE>`), with a maximum merge gap of `<MACS_MAX_GAP>` bp (`--max-gap <MACS_MAX_GAP> -g <GENOME_SIZE_FOR_MACS>`).

**C4 — csaw (differential window-based)**

> **CRITICAL:** Read the helper scripts called (e.g., `csaw_parameters.R`, `callpeaks_csaw.R`) before writing this paragraph.

**Template:**

> Differential chromatin accessibility/binding across conditions was assessed using a sliding-window approach implemented in the csaw Bioconductor package (v`<VERSION>`; [citation]; R v`<VERSION>`). Read counts were tallied across `<WINDOW_SIZE>` bp windows sliding across the genome with a step size of `<STEP_SIZE>` bp (`readParam(minq=<minq>, pe="<pe>", max.frag=<max.frag>, dedup=<TRUE/FALSE>)`). Low-abundance windows with average log-CPM below `<FILTER_THRESHOLD>` were filtered prior to testing. TMM [or spike-in normalization] offsets were incorporated (`normOffsets()`), and differential windows between conditions were identified using edgeR (v`<VERSION>`; [citation]) quasi-likelihood F-tests (`glmQLFTest`). P-values were adjusted using the Benjamini-Hochberg procedure, and significant windows were merged (`mergeResults()`, `combineTests()`) at an FDR threshold of `<FDR>`.

---

#### BLOCK D — Peak and Gene Annotation

> **Trigger:** `Analysis_Data/peaks/annotation/` exists, or an annotation batch file is present.

**Template:**

> Genomic feature annotation of peak regions (promoters, exons, introns, and intergenic regions) was performed using ChIPseeker (v`<VERSION>`; [citation]; R v`<VERSION>`) against the `<TxDb package>` transcriptome database. Promoter regions were defined as ± `<PROMOTER_WINDOW>` bp around the transcription start site (TSS).

---

#### BLOCK E — Differential Gene Expression (RNA-seq)

> **Trigger:** `EXPERIMENT == "rnaseq"` and `Analysis_Data/dge/` exists.

Read the DGE helper scripts (e.g., `run_deseq2.R`, `run_edger.R`) along with the batch logs to verify exact formulas and cutoffs before writing this section.

**E1 — DESeq2**

**Template:**

> Gene-level read counts were quantified from aligned BAM files using featureCounts from the Subread package (v`<VERSION>`; [citation]) against the `<GTF basename>` gene annotation (`-s <strandness> -p -Q <minq>`). Differential gene expression analysis across conditions was performed using DESeq2 (v`<VERSION>`; [citation]; R v`<VERSION>`) modeling the design formula `<DESIGN_FORMULA>`. Library size normalization was performed using the median-of-ratios method, and pairwise differential expression between `<GROUP_A>` and `<GROUP_B>` was tested using two-tailed Wald tests. P-values were corrected using the Benjamini-Hochberg procedure (`padj`), with significantly differentially expressed genes defined by an adjusted p-value < `<FDR>` and absolute log₂ fold-change > `<LFC>`.

**E2 — edgeR**

**Template:**

> Gene-level read counts were quantified from aligned BAM files using featureCounts from the Subread package (v`<VERSION>`; [citation]) against the `<GTF basename>` gene annotation (`-s <strandness> -p -Q <minq>`). Differential gene expression analysis across conditions was performed using edgeR (v`<VERSION>`; [citation]; R v`<VERSION>`). Library sizes were normalized using Trimmed Mean of M-values (TMM) scaling (`calcNormFactors`), and empirical Bayes dispersion parameters were estimated (`estimateDisp`) under the design matrix `<DESIGN_MATRIX>`. Differential expression between `<GROUP_A>` and `<GROUP_B>` was evaluated using quasi-likelihood F-tests (`glmQLFTest`). Multiple testing correction was applied via the Benjamini-Hochberg method (`topTags`), defining significant genes by FDR < `<FDR>` and absolute log₂ fold-change > `<LFC>`.

---

#### BLOCK F — CRISPR Screen Analysis

> **Trigger:** `EXPERIMENT == "crisprscreen"` or MAGeCK batch files exist in `Analysis_Data/counts/` or `Analysis_Data/mageck/`.

Read the MAGeCK batch files and helper scripts to verify library definitions, trimming flags, and statistical parameters.

**F1 — MAGeCK Robust Rank Aggregation (RRA)**

**Template:**

> sgRNA read counts were quantified from raw sequencing reads using MAGeCK (v`<VERSION>`; [citation]) `count` against the `<LIBRARY_FILE>` sgRNA library (`--trim-5 <TRIM_5> --mismatch-zero <TRUE/FALSE> --mismatch <MISMATCHES> --norm-method <NORM_METHOD>`). Gene-level essentiality and enrichment scores between conditions (`-c <CONTROL_SAMPLES> -t <TREATMENT_SAMPLES>`) were computed using MAGeCK `test`, which applies modified robust rank aggregation (RRA) (`--gene-lfc-method <LFC_METHOD>`) to aggregate sgRNA-level p-values by gene. Significantly enriched or depleted genes were identified at a false discovery rate (FDR) < `<FDR>`.

**F2 — MAGeCK Maximum Likelihood Estimation (MLE)**

**Template:**

> sgRNA read counts were quantified from raw sequencing reads using MAGeCK (v`<VERSION>`; [citation]) `count` against the `<LIBRARY_FILE>` sgRNA library (`--trim-5 <TRIM_5> --mismatch <MISMATCHES> --norm-method <NORM_METHOD>`). For multi-condition experimental designs, gene-level beta scores and statistical significance were estimated using the MAGeCK `mle` module across the design matrix `<MATRIX_FILE>` (`--design-matrix <MATRIX_FILE> --norm-method <NORM_METHOD>`). Genes showing significant fitness effects across conditions were identified at FDR < `<FDR>`.

---

#### BLOCK G — Visualization and Signal Profiles

> **Trigger:** `Analysis_Data/deeptools/` exists.

Read the deepTools batch files (`computeMatrix`, `plotHeatmap`, `plotProfile`) to extract exact genomic windows and clustering parameters.

**Template:**

> Genomic signal matrices across target features (`<peak summits / TSS / gene bodies>`) were computed from normalized BigWig coverage tracks using `computeMatrix` from the deepTools suite (v`<VERSION>`; [citation]) in `<reference-point / scale-regions>` mode (`--referencePoint <TSS/center> --upstream <UPSTREAM> --downstream <DOWNSTREAM> --binSize <BIN_SIZE> --missingDataAsZero`). [If plotHeatmap:] Signal enrichment heatmaps across features were generated using `plotHeatmap` (`--colorMap <COLORMAP> --kmeans <N>` [or `--hclust <N>`]). [If plotProfile:] Metagene average signal profiles across conditions were generated using `plotProfile` (`--perGroup --plotType <lines/se>`).

---

## 3 — Reference Genome and Annotation Details

Always include a paragraph describing the reference genome and annotation source used:

**T2T (CHM13v2.0):**
> All analyses were performed using the Telomere-to-Telomere (T2T) human reference genome assembly CHM13v2.0.

**hg19 / hg38 / GRCh38:**
> All analyses were performed using the GRCh38 (hg38) human reference genome.

**mm39 / GRCm39:**
> All analyses were performed using the GRCm39 mouse reference genome.

Supplement with the annotation database source if applicable (e.g., NCBI RefSeq, GENCODE, Ensembl) as observed in the `GNM_GTF` config variable.

---

## 4 — Output Format: Microsoft Word (.docx) or MS Word-Compatible Format

The generated Materials and Methods document **MUST** be output in a **Microsoft Word `.docx` or other MS Word-compatible format (preferential)**. If and only if Word-compatible generation tools (`python-docx`, `pandoc`, etc.) cannot be executed across any active environment, output a cleanly structured **Markdown (`.md`) file (fallback)**.

### 4.1 — Preference Order and Tool Availability Check

Before writing the output document, verify which tools are available across your active conda environments or system Python:

1. **Python `python-docx` (`.docx` — Preferred):** Check if `python-docx` is available or can be installed:
   ```bash
   python3 -c "import docx; print('available')" 2>/dev/null \
     || conda run -n chipseq python -c "import docx; print('available')" 2>/dev/null \
     || conda run -n rnaseq python -c "import docx; print('available')" 2>/dev/null \
     || echo "not_available"
   ```
   *If `python-docx` is not installed, attempt to install it via `pip install --user python-docx` or `conda install -y python-docx` in a writable environment.*

2. **Pandoc (`.docx` / `.rtf` — Preferential Alternative):** If `python-docx` cannot be used, check if `pandoc` is available:
   ```bash
   which pandoc 2>/dev/null && echo "pandoc_available" || echo "not_available"
   ```
   *If available, first write the complete M&M text to `Materials_and_Methods.md`, then convert it directly to `.docx`: `pandoc -s Materials_and_Methods.md -o Materials_and_Methods.docx`.*

3. **Markdown (`.md` — Fallback):** If neither `python-docx` nor `pandoc` can generate a Word-compatible file, output a well-structured `Materials_and_Methods.md` file that the user can open and import directly into Microsoft Word.

**Primary Save Location:** `<PROJECT_ROOT>/Materials_and_Methods.docx` (or `.rtf` / `.md` fallback)

---

### 4.2 — Python Script Template for `.docx` Generation (`python-docx`)

Write a Python script to `Scripts/utils/generate_methods_docx.py` and execute it:

```python
import datetime
from docx import Document
from docx.shared import Pt, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH

doc = Document()

# --- Document Title ---
title = doc.add_heading("Materials and Methods", level=1)
title.alignment = WD_ALIGN_PARAGRAPH.CENTER

# --- Subsections ---
# For each analytical step, add a level-2 heading and paragraph body:
# doc.add_heading("<SUBSECTION TITLE>", level=2)
# p = doc.add_paragraph("<SUBSECTION TEXT>")
# p.paragraph_format.first_line_indent = Pt(18)

doc.add_heading("Read Alignment and Processing", level=2)
p1 = doc.add_paragraph(
    "Raw paired-end sequencing reads were processed to remove sequencing adapters..."
)
p1.paragraph_format.first_line_indent = Pt(18)

# --- References Section ---
doc.add_heading("References", level=1)
ref_list = [
    "1. Bolger AM, Lohse M, Usadel B. Trimmomatic: a flexible trimmer for "
    "Illumina sequence data. Bioinformatics. 2014;30(15):2114-2120. "
    "doi:10.1093/bioinformatics/btu170",
    # ... add all references mentioned in the text in numerical order
]
for ref in ref_list:
    doc.add_paragraph(ref, style="List Number")

# --- Footer ---
section = doc.sections[0]
footer = section.footer
footer.paragraphs[0].text = (
    f"Generated by AI agent on {datetime.date.today().isoformat()} | "
    "Project: <PROJECT_NAME>"
)

doc.save("Materials_and_Methods.docx")
print("Successfully saved Materials_and_Methods.docx")
```

Run the generator:
```bash
python3 Scripts/utils/generate_methods_docx.py
```

---

### 4.3 — Markdown Fallback Template (`.md`)

If MS Word-compatible generation is unavailable (`not_available`), save the output to `<PROJECT_ROOT>/Materials_and_Methods.md`:

```markdown
# Materials and Methods

## Read Alignment and Processing

<Text>

## Normalization

<Text>

## Peak Calling

<Text>

...

## References

1. Bolger AM, Lohse M, Usadel B. Trimmomatic: a flexible trimmer for Illumina sequence data. Bioinformatics. 2014;30(15):2114-2120. doi:10.1093/bioinformatics/btu170
2. ...
```

---

## 5 — References Section

The References section **must** appear at the end of every M&M document. Include a numbered citation for **every tool** mentioned in the text, ordered chronologically by first appearance. Use the following curated citations (add more as needed):

```
[Trimmomatic]
Bolger AM, Lohse M, Usadel B. Trimmomatic: a flexible trimmer for Illumina
sequence data. Bioinformatics. 2014;30(15):2114-2120.
doi:10.1093/bioinformatics/btu170

[Bowtie2]
Langmead B, Salzberg SL. Fast gapped-read alignment with Bowtie 2.
Nat Methods. 2012;9(4):357-359. doi:10.1038/nmeth.1923

[HISAT2]
Kim D, Paggi JM, Park C, Bennett C, Salzberg SL. Graph-based genome alignment
and genotyping with HISAT2 and HISAT-genotype. Nat Biotechnol. 2019;37(8):907-915.
doi:10.1038/s41587-019-0201-4

[SAMtools]
Danecek P, Bonfield JK, Liddle J, et al. Twelve years of SAMtools and BCFtools.
GigaScience. 2021;10(2):giab008. doi:10.1093/gigascience/giab008

[Sambamba]
Tarasov A, Vilella AJ, Cuppen E, Nijman IJ, Prins P. Sambamba: fast processing
of NGS alignment formats. Bioinformatics. 2015;31(12):2032-2034.
doi:10.1093/bioinformatics/btv098

[BEDTools]
Quinlan AR, Hall IM. BEDTools: a flexible suite of utilities for comparing
genomic features. Bioinformatics. 2010;26(6):841-842.
doi:10.1093/bioinformatics/btq033

[deepTools / bamCoverage]
Ramírez F, Ryan DP, Grüning B, et al. deepTools2: a next generation web server
for deep-sequencing data analysis. Nucleic Acids Res. 2016;44(W1):W160-165.
doi:10.1093/nar/gkw257

[MACS2 / MACS3]
Zhang Y, Liu T, Meyer CA, et al. Model-based analysis of ChIP-Seq (MACS).
Genome Biol. 2008;9(9):R137. doi:10.1186/gb-2008-9-9-r137

[SEACR]
Meers MP, Tenenbaum D, Bhanu Bhanu K. Peak calling by Sparse Enrichment
Analysis for CUT&RUN chromatin profiling. Epigenetics Chromatin.
2019;12(1):42. doi:10.1186/s13072-019-0287-4

[featureCounts / Subread]
Liao Y, Smyth GK, Shi W. featureCounts: an efficient general purpose program
for assigning sequence reads to genomic features. Bioinformatics.
2014;30(7):923-930. doi:10.1093/bioinformatics/btt656

[DESeq2]
Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion
for RNA-seq data with DESeq2. Genome Biol. 2014;15(12):550.
doi:10.1186/s13059-014-0550-8

[edgeR]
Robinson MD, McCarthy DJ, Smyth GK. edgeR: a Bioconductor package for
differential expression analysis of digital gene expression data.
Bioinformatics. 2010;26(1):139-140. doi:10.1093/bioinformatics/btp616

[csaw]
Lun AT, Smyth GK. csaw: a Bioconductor package for differential binding
analysis of ChIP-seq data using sliding windows. Nucleic Acids Res.
2016;44(5):e45. doi:10.1093/nar/gkv1191

[limma / voom]
Ritchie ME, Phipson B, Wu D, et al. limma powers differential expression
analyses for RNA-sequencing and microarray studies. Nucleic Acids Res.
2015;43(7):e47. doi:10.1093/nar/gkv007

[ChIPseeker]
Yu G, Wang LG, He QY. ChIPseeker: an R/Bioconductor package for ChIP peak
annotation, comparison and visualization. Bioinformatics.
2015;31(14):2382-2383. doi:10.1093/bioinformatics/btv145

[MAGeCK]
Li W, Xu H, Xiao T, et al. MAGeCK enables robust identification of essential
genes from genome-scale CRISPR/Cas9 knockout screens. Genome Biol.
2014;15(12):554. doi:10.1186/s13059-014-0554-4

[ChIPseqSpikeInFree]
Jin H, Kasper LH, Czerwinska JM, et al. ChIPseqSpikeInFree: A ChIP-seq
normalization approach to reveal global changes in histone modifications
without spike-in. Bioinformatics. 2020;36(4):1270-1272.
doi:10.1093/bioinformatics/btz720

[T2T-CHM13v2.0]
Nurk S, Koren S, Rhie A, et al. The complete sequence of a human genome.
Science. 2022;376(6588):44-53. doi:10.1126/science.abj6987
```

---

## 6 — Quality Control Checks Before Finalizing

Before saving and reporting the output document, verify:

- [ ] Every tool mentioned in the M&M text has a version number on first mention (`vX.Y.Z`).
- [ ] Every tool mentioned in the text has a numbered `[citation]` placeholder corresponding to the References section.
- [ ] No steps or tools are described that were skipped (cross-checked against `config.py` and log folders).
- [ ] Key parameters (cutoffs, window sizes, flags) match verbatim what was observed in the `.batch` scripts and helper scripts (`.R`/`.py`).
- [ ] Reference genome and annotation are spelled out fully (e.g., T2T-CHM13v2.0, GRCh38, GRCm39).
- [ ] Output is saved in **MS Word-compatible format (`.docx` preferred)** or structured `.md` (fallback) at `<PROJECT_ROOT>/Materials_and_Methods.<docx|md>`.

---

## 7 — Experiment-Type Checklist

Use this checklist to confirm which M&M BLOCKs apply for each genomics pipeline type:

| Block | ATAC-seq | ChIP-seq | CUT&RUN / CUT&TAG | RNA-seq | CRISPR Screen |
|---|---|---|---|---|---|
| **A — Read Alignment & Processing** | ✓ (Bowtie2) | ✓ (Bowtie2) | ✓ (Bowtie2) | ✓ (HISAT2) | — |
| **B — Normalization & Coverage** | B1 (TMM) | B1 / B2 / B3 | B2 (Spike-in) | — | — |
| **C — Peak Calling** | C1 (ATAC) | C2 / C3 | C3 / C4 | — | — |
| **D — Peak / Gene Annotation** | optional | optional | optional | — | — |
| **E — Differential Gene Expression** | — | — | — | ✓ (DESeq2 / edgeR) | — |
| **F — CRISPR Screen Analysis** | — | — | — | — | ✓ (MAGeCK) |
| **G — Visualization & Profiles** | optional | optional | optional | optional | — |

---

## 8 — Example Invocation Workflow (Summary)

```
1.  Read config.py                         → record all parameters
2.  List non-empty output directories      → determine which blocks apply
3.  For each block:
    a. Read ONE batch file from log/       → extract tool + flags
    b. If helper scripts are called →
       read each helper script             → extract analytical parameters
4.  For each conda env in batch files →
    run: conda env export -n <ENV>         → extract tool versions
5.  Write M&M text block by block          → follow templates in §2
6.  Compile References section             → cite every tool used
7.  Check Word-compatible tools (§4.1)     → python-docx or pandoc
8.  Generate Materials_and_Methods.docx
    (or .md if Word-compatible generation unavailable)
9.  Run quality control checklist (§6)
10. Report output file path to user
```
