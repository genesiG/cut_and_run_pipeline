#!/usr/bin/env python3

"""
config.py - Configuration module for the genomic analysis pipeline
          - Change file paths for your specific project
"""

import os

### SETTINGS (Change as needed)
# Your HPC username
USERNAME = "<your_username>"

# Project name (no spaces allowed)
PROJECT_NAME = "GG_EPICYPHER_CBX2"

BSUB_QUEUE = "rhel9"

# What type of experiment is being analyzed?
# 'rnaseq' | 'chipseq' | 'atacseq' | 'cutandrun' | 'cutandtag'
EXPERIMENT = "cutandrun"

# Genome build (e.g. "hg19" or "t2t" for human genome build)
SPECIES = "t2t_ecoli"

# Alignment tool used in the pipeline
ALIGNER = "bowtie2"

# Whether the data should be processed as PE or SE
IS_PAIRED_END = True
FORWARD_PATTERN = ('_R1',)
REVERSE_PATTERN = ('_R2',)

# BIGWIG
BIN_SIZE = 100
# 'cpm' | 'rpkm' | 'spikein' | 'chipseqspikeinfree' | 'tmm'
NORMALIZE_USING = "SPIKEIN"

# Are you starting the analysis from raw sequencing files (True), or
# from BAM files downloaded from the literature (False)?
IS_RAW_DATA_FASTQ = True

# Minimum read mapping quality (bowtie2)
MAPQ = 30

### METADATA PATTERNS
# Used by step_0c to build samples.txt and sample_metadata.txt
# Maps pattern found in filename → antibody label
# Order matters: more specific patterns must come first
ANTIBODY_PATTERNS = {
    'K27me3': 'K27me3',
    'K27me2': 'K27me2',
    'RbIgG':    'IgG',
    'AntiCBX2RP': 'CBX2',
    'GSTCBX2': 'GSTCBX2',
    'GSTControl': 'GSTControl',
    'EZH2': 'EZH2',
    'GSTCBX7': 'GSTCBX7',
    'K4me3': 'K4me3',
    'H2AK119ub': 'H2AK119ub',
}

# Maps pattern found in filename → group label
GROUP_PATTERNS = {
    'DMSO'    : 'DMSO',
    '1uMcpd'    : '1uM',
    '5uMcpd'  : '5uM'
}

# NORMALIZATION PARAMETERS
# Dictionary where keys are run names (one per antibody target), and values
# are dicts with parameters.
#
# Key design:
#   - 'targets'   : Regex passed to select BAM files for the run.
#                   Pattern matches all samples (e.g. DMSO, EEDi, EPZ)
#                   for a given antibody so that relative changes are captured.
#   - 'bin_size'  : (step_4a only) genomic bin size in bp for ChIPseqSpikeInFree
#   - 'cutoff'    : (step_4a only) signal cutoff for ChIPseqSpikeInFree
#   - 'max_turns' : (step_4a only) fraction of bins used for fitting
NORMALIZATION_RUNS = {
    "IgG": {
        "targets": "IgG",
        "bin_size": 10000,
        "cutoff": 1.2,
        "max_turns": 0.999
    },
    "K27me2": {
        "targets": "K27me2",
        "bin_size": 10000,
        "cutoff": 1.2,
        "max_turns": 0.999
    },
    "K27me3": {
        "targets": "K27me3",
        "bin_size": 10000,
        "cutoff": 1.2,
        "max_turns": 0.99
    },
    "K4me3": {
        "targets": "K4me3",
        "bin_size": 10000,
        "cutoff": 1.2,
        "max_turns": 0.99
    },
    "H2AK119ub": {
        "targets": "H2AK119ub",
        "bin_size": 10000,
        "cutoff": 1.2,
        "max_turns": 0.99
    },  
    "CBX2": {
        "targets": "AntiCBX2",
        "bin_size": 10000,
        "cutoff": 1.2,
        "max_turns": 0.99
    },
    "EZH2": {
        "targets": "EZH2",
        "bin_size": 10000,
        "cutoff": 1.2,
        "max_turns": 0.99
    },
    "GSTCBX2": {
        "targets": "GSTCBX2",
        "bin_size": 10000,
        "cutoff": 1.2,
        "max_turns": 0.99
    },
    "GSTControl": {
        "targets": "GSTControl",
        "bin_size": 10000,
        "cutoff": 1.2,
        "max_turns": 0.99
    },
    "GSTCBX7": {
        "targets": "GSTCBX7",
        "bin_size": 10000,
        "cutoff": 1.2,
        "max_turns": 0.99
    }
}
###

# ---------------------------------------------------------------------------
# Auto-configure run parameters from experiment type.
# These values are derived programmatically and should not be overridden
# manually below — change EXPERIMENT instead.
# ---------------------------------------------------------------------------
if EXPERIMENT in ("cutandrun", "cutandtag"):
    USE_TRIMMOMATIC   = True
    REMOVE_DUPLICATES = False   # Mark duplicates only (IgG samples are always deduped)
    REMOVE_BLACKLIST  = True
elif EXPERIMENT in ("chipseq", "atacseq"):
    USE_TRIMMOMATIC   = True
    REMOVE_DUPLICATES = True
    REMOVE_BLACKLIST  = True
elif EXPERIMENT == "rnaseq":
    USE_TRIMMOMATIC   = True
    REMOVE_DUPLICATES = False
    REMOVE_BLACKLIST  = False
else:
    # Fallback defaults
    USE_TRIMMOMATIC   = False
    REMOVE_DUPLICATES = True
    REMOVE_BLACKLIST  = True


### Set up work directory
HOMEDIR         = os.path.join("/home", USERNAME)
WORKDIR         = os.path.join(HOMEDIR, PROJECT_NAME)
CODEDIR         = os.path.join(WORKDIR, "Scripts")
ORIGINAL_DATA   = os.path.join(WORKDIR, "Original_Data")
METADATA        = os.path.join(WORKDIR, "Metadata")
ANALYSIS_DATA   = os.path.join(WORKDIR, "Analysis_Data")
IMPORTABLE_DATA = os.path.join(WORKDIR, "Importable_Data")
# Conda environment binaries path
CONDA_ENV_BIN = os.path.join(HOMEDIR, "miniconda3", "envs", "chipseq", "bin")
###

### Set up reference files
# Reference genome
if SPECIES == "t2t":
    # bowtie2 indexes (e.g. CHM13v2.0)
    GNM_IDX = f"{HOMEDIR}/mylibrary/{ALIGNER}/t2t/chm13v2.0"
    # Annotation
    GNM_GTF = f"{HOMEDIR}/mylibrary/genomes/t2t/t2t.ncbiRefSeq.curated.norandom.gtf"
    # TSS
    TSS_BED = f"{HOMEDIR}/mylibrary/genomes/t2t/t2t.curated.norandom.tss.bed"
elif SPECIES in ("t2t_ecoli", "t2t_mm39"):
    # bowtie2 index
    GNM_IDX = f"{HOMEDIR}/mylibrary/{ALIGNER}/{SPECIES}/{SPECIES}"
    # Annotation
    GNM_GTF = f"{HOMEDIR}/mylibrary/genomes/t2t/t2t.ncbiRefSeq.curated.norandom.gtf"
    # TSS
    TSS_BED = f"{HOMEDIR}/mylibrary/genomes/t2t/t2t.curated.norandom.tss.bed"
else:
    # bowtie2 index
    GNM_IDX = f"{HOMEDIR}/mylibrary/{ALIGNER}/{SPECIES}/{SPECIES}"
    # Annotation
    GNM_GTF = f"{HOMEDIR}/mylibrary/genomes/{SPECIES}/{SPECIES}.ncbiRefSeq.curated.norandom.gtf "
    # TSS
    TSS_BED = f"{HOMEDIR}/mylibrary/genomes/{SPECIES}/{SPECIES}.curated.norandom.tss.bed"

# Chromosome sizes
GNM_SIZES = f"{HOMEDIR}/mylibrary/chr_size/{SPECIES}.chrom.sizes"

# Blacklisted genome regions
if SPECIES == "mm39":
    # BLACKLIST = f"{HOMEDIR}/mylibrary/blacklist/mm10-blacklist.v2.Liftover.mm39.bed"
    BLACKLIST = (f"{HOMEDIR}/mylibrary/blacklist/{SPECIES}.blacklist.cut.and.run.bed"
                 if EXPERIMENT == "cutandrun"
                 else f"{HOMEDIR}/mylibrary/blacklist/{SPECIES}.blacklist.bed")
elif SPECIES == "t2t_ecoli" or SPECIES == "t2t_mm39":
    BLACKLIST = (f"{HOMEDIR}/mylibrary/blacklist/t2t.blacklist.cut.and.run.bed"
                 if EXPERIMENT == "cutandrun"
                 else f"{HOMEDIR}/mylibrary/blacklist/t2t.blacklist.bed")
else:
    BLACKLIST = (f"{HOMEDIR}/mylibrary/blacklist/{SPECIES}.blacklist.cut.and.run.bed"
                 if EXPERIMENT == "cutandrun"
                 else f"{HOMEDIR}/mylibrary/blacklist/{SPECIES}.blacklist.bed")

###

### Set up data folders
# Output FASTQC reports
QCDIR1 = os.path.join(ANALYSIS_DATA, "qc")

# Output trimmed sequencing data
TRIMDIR = os.path.join(IMPORTABLE_DATA, "trimmed")

# Output alignment data
ALIGNDIR        = os.path.join(ANALYSIS_DATA, "alignment")
BAMDIR          = os.path.join(ANALYSIS_DATA, "bam")
ORIGINALBAMDIR = os.path.join(BAMDIR, "original")

# Output for processed aligned data
PROCESSEDBAMDIR = os.path.join(BAMDIR, "processed")
BEDDIR          = os.path.join(IMPORTABLE_DATA, "bed")
SCALINGDIR      = os.path.join(ANALYSIS_DATA, "normalization")
BEDGRAPHDIR     = os.path.join(ANALYSIS_DATA, "bedgraph")
BIGWIGDIR       = os.path.join(ANALYSIS_DATA, "bigwig")

# Spike-in consistency guard
# If spike-in normalisation is requested, SPECIES must encode the spike-in
# organism so that (a) Bowtie2 aligns to the combined genome in step 2, and
# (b) step 3 knows which contigs to split out.
if NORMALIZE_USING == "SPIKEIN" and not any(
    tag in SPECIES for tag in ("_ecoli", "_mm39")
):
    raise ValueError(
        f"\n\n[config.py] Inconsistent settings detected:\n"
        f"  NORMALIZE_USING = '{NORMALIZE_USING}' requires a spike-in genome,\n"
        f"  but SPECIES = '{SPECIES}' does not include a spike-in suffix.\n\n"
        f"  Fix: set SPECIES to the combined build before running any pipeline step.\n"
        f"  Examples:\n"
        f"    SPECIES = 't2t_ecoli'   # human T2T + E. coli spike-in\n"
        f"    SPECIES = 't2t_mm39'    # human T2T + mouse mm39 spike-in\n"
    )

# PEAK CALLING
# Output peak calling files
PEAKDIR = os.path.join(ANALYSIS_DATA, "peaks")

# Specify replicates for peak calling
# Samples will be searched with format: <ANTIBODY_TARGET>_<SAMPLE>
ANTIBODY_TARGETS = ["DMSO", "1uMcpd"]
SAMPLES          = ["H3K27me3",
                    "H3K27me2",
                    "AntiCBX2RP",
                    "H2AK119ub"
                    ]

# Output peak annotation files
PEAKANNODIR = os.path.join(PEAKDIR, "annotation")

# Parameters for MACS
MACS_Q_VALUE = 0.01
MACS_FE_CUTOFF = 3.0
MACS_MAX_GAP = 500
MACS_CALL_SUMMITS = False
MACS_CALL_BROAD_PEAKS = True
if SPECIES in ["t2t", "t2t_ecoli", "t2t_mm39"]:
    GENOME_SIZE_FOR_MACS = 3.1e9
elif SPECIES == "mm39":
    GENOME_SIZE_FOR_MACS = 2.7e9
else:
    GENOME_SIZE_FOR_MACS = "hs"

# Step 5d — SEACR peak calling output
SEACR_OUTDIR      = os.path.join(PEAKDIR, "seacr")
CSAW_OUTDIR       = os.path.join(PEAKDIR, "csaw")
BEDGRAPHDIR = os.path.join(IMPORTABLE_DATA, "bedgraph")

# Output motif discovery files
MOTIFDIR = os.path.join(PEAKDIR, "motif_analysis")

# Parameters for heatmap generation
BUILD_MATRIX_USING = "deepTools" # options: "EnrichedHeatmap", "deepTools"
DEEPTOOLSDIR       = os.path.join(ANALYSIS_DATA, "deeptools")
HEATMAPDIR         = os.path.join(ANALYSIS_DATA, "heatmaps")
###

### QC sub-directories
# qc_2: MultiQC report
MULTIQCDIR          = os.path.join(QCDIR1, "multiqc")
# qc_3: Alignment statistics
ALIGN_STATS_DIR     = os.path.join(QCDIR1, "alignment_stats")
# qc_4: Fragment size distributions
FRAG_SIZE_DIR       = os.path.join(QCDIR1, "fragment_sizes")
# qc_5: Library complexity
LIB_COMPLEXITY_DIR  = os.path.join(QCDIR1, "library_complexity")
# qc_6: deepTools fingerprint
FINGERPRINT_DIR     = os.path.join(QCDIR1, "fingerprint")
# qc_7: Reproducibility (correlation / PCA)
REPRODUCIBILITY_DIR = os.path.join(QCDIR1, "reproducibility")
# qc_8: FRiP
FRIPDIR             = os.path.join(QCDIR1, "frip")
# qc_9: Peak statistics
PEAK_STATS_DIR      = os.path.join(QCDIR1, "peak_stats")
# Spike-in chromosome prefix used in qc_3 to separate spike-in reads
SPIKE_CHR_PREFIX    = "ecoli_chr" if SPECIES in ["t2t_ecoli"] else (
                       "mm39_chr"  if SPECIES in ["t2t_mm39"]  else "")
###

### BAM FILES AS RAW DATA
# Path to folder with raw data to start analysis
if IS_RAW_DATA_FASTQ:
    DATADIR = os.path.join(ORIGINAL_DATA, "fastq_files")
else:
    DATADIR = os.path.join(BAMDIR, "original")

# Define file extensions from starting data in the analysis
if IS_RAW_DATA_FASTQ:
    EXTENSIONS = ['.fastq.gz', '.fastq']
else:
    EXTENSIONS = ['.bam']
