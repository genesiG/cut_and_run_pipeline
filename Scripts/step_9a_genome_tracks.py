#!/usr/bin/env python3

"""
step_9a_genome_tracks.py  —  Submit bsub jobs that render genome track views
                              using pyGenomeTracks and the .ini config in Scripts/utils/.

Usage (single gene, original behaviour):
    python3 Scripts/step_9a_genome_tracks.py --gene CBX2 --flank 125000

Usage (multiple genes with per-gene flanks):
    python3 Scripts/step_9a_genome_tracks.py \\
        --genes CBX2 CBX4 PSMD8 \\
        --flanks 125000 150000 200000

One bsub job is submitted for each gene : flank pair.

Options:
    --gene    Single gene locus (mutually exclusive with --genes)  [default: CDKN2A]
    --flank   Flanking bp for --gene                               [default: 125000]
    --genes   Space-separated list of gene loci
    --flanks  Space-separated list of flanking bp (must match --genes length)
    --region  Explicit region in UCSC format (overrides gene/flank) [default: None]
    --ini     Path to the .ini track configuration   [default: Scripts/utils/cdkn2a.ini]
    --out     Output SVG path (only valid for a single gene/region) [default: auto]

Outputs: Analysis_Data/genome_tracks/
Logs:    Analysis_Data/genome_tracks/log/
"""

### Import modules
import os
import sys
import time
import argparse
import config
###

### Constants
NUM_CORES    = 4
MAX_MEM      = 16000
MEM_PER_CORE = MAX_MEM // NUM_CORES

# Default gene : flank pairs (used when neither --gene/--genes nor --flanks are supplied)
DEFAULT_GENES  = ["CDKN2A"]
DEFAULT_FLANKS = [125000]
###

### Paths
TRACKS_DIR  = os.path.join(config.ANALYSIS_DATA, "genome_tracks")
LOG_DIR     = os.path.join(TRACKS_DIR, "log")
INI_DEFAULT = os.path.join(config.CODEDIR, "utils", "cdkn2a.ini")

os.makedirs(TRACKS_DIR, exist_ok=True)
os.makedirs(LOG_DIR,    exist_ok=True)
###


def sanitize_region(region: str) -> str:
    """Convert a region string like 'chr9:21,927,246-22,127,450' to a filename-safe slug."""
    return region.replace(":", "_").replace(",", "").replace("-", "_")


def get_gene_coordinates(gtf_path: str, gene: str) -> tuple:
    """Reads the GTF and returns (chrom, start, end) for the given gene."""
    min_pos = float('inf')
    max_pos = float('-inf')
    chrom = None

    if not os.path.exists(gtf_path):
        print(f"ERROR: GTF not found at {gtf_path}", file=sys.stderr)
        sys.exit(1)

    print(f"Scanning {gtf_path} for '{gene}'...")
    with open(gtf_path, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith("#"):
                continue
            if f'"{gene}"' in line:
                parts = line.split('\t')
                min_pos = min(min_pos, int(parts[3]))
                max_pos = max(max_pos, int(parts[4]))
                chrom = parts[0]

    if chrom is None:
        print(f"ERROR: Gene '{gene}' not found in {gtf_path}", file=sys.stderr)
        sys.exit(1)

    return chrom, min_pos, max_pos


def submit_tracks_job(region: str, slug: str, ini_path: str, out_file: str) -> None:
    """Build and submit a single bsub job that runs pyGenomeTracks."""
    job_name    = f"{config.PROJECT_NAME}.genome_tracks.{slug}"
    log_out     = os.path.join(LOG_DIR, f"{slug}.log")
    log_err     = os.path.join(LOG_DIR, f"{slug}.error")
    batch_file  = os.path.join(LOG_DIR, f"{slug}.batch")

    if not os.path.exists(ini_path):
        print(f"ERROR: .ini file not found: {ini_path}", file=sys.stderr)
        sys.exit(1)

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

# Run from project root so config.py resolves all relative paths correctly.
cd {config.WORKDIR}

echo "=== step_9a: genome track view ==="
echo "Project : {config.PROJECT_NAME}"
echo "Region  : {region}"
echo "INI     : {ini_path}"
echo "Output  : {out_file}"
echo "Started : $(date)"
echo ""

pyGenomeTracks \\
    --tracks {ini_path} \\
    --region {region} \\
    --outFileName {out_file} \\
    --width 15 \\
    --dpi 600

echo ""
echo "Finished: $(date)"
echo "=== Resource usage summary ==="
bjobs -l $LSB_JOBID | grep -A10 "Resource usage"
"""

    with open(batch_file, "w", encoding="utf-8") as fh:
        fh.write(batch_cmd)

    print(f"Submitting job : {job_name}")
    print(f"  Region  : {region}")
    print(f"  INI     : {ini_path}")
    print(f"  Output  : {out_file}")
    print(f"  Log     : {log_out}")
    print(f"  Error   : {log_err}")
    os.system(f"bsub < {batch_file}")
    time.sleep(1)


def resolve_region(gene: str, flank: int) -> tuple[str, str]:
    """Return (region_string, slug) for a given gene + flank."""
    chrom, start, end = get_gene_coordinates(config.GNM_GTF, gene)
    pad_start = max(0, start - flank)
    pad_end   = end + flank
    region    = f"{chrom}:{pad_start}-{pad_end}"
    flank_kb  = flank // 1000
    slug      = f"{gene}_{flank_kb}kb"
    return region, slug


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate genome track view plots using pyGenomeTracks."
    )

    # ── Single-gene arguments (original interface) ─────────────────────────────
    parser.add_argument(
        "--gene",
        default=None,
        help="Single gene locus of interest. Mutually exclusive with --genes."
    )
    parser.add_argument(
        "--flank",
        type=int,
        default=None,
        help="Base pairs to pad around --gene. Ignored when --genes is used."
    )
    parser.add_argument(
        "--region",
        default=None,
        help="Explicit genomic region in UCSC format. Overrides gene/flank if provided."
    )

    # ── Multi-gene arguments ────────────────────────────────────────────────────
    parser.add_argument(
        "--genes",
        nargs="+",
        default=None,
        metavar="GENE",
        help=(
            "Space-separated list of gene loci. "
            "Must be paired 1-to-1 with --flanks."
        )
    )
    parser.add_argument(
        "--flanks",
        nargs="+",
        type=int,
        default=None,
        metavar="BP",
        help=(
            "Space-separated list of flanking base pairs, one per gene in --genes."
        )
    )

    # ── Shared arguments ────────────────────────────────────────────────────────
    parser.add_argument(
        "--ini",
        default=INI_DEFAULT,
        help=f"Path to the pyGenomeTracks .ini config file. [default: {INI_DEFAULT}]"
    )
    parser.add_argument(
        "--out",
        default=None,
        help="Output SVG path. Only valid for a single gene/region run. [default: auto]"
    )
    args = parser.parse_args()

    # ── Validate and build the list of (region, slug, out_file) jobs ───────────

    # 1. Explicit region overrides everything
    if args.region:
        if args.genes or args.gene:
            print("WARNING: --region overrides --gene/--genes; gene arguments are ignored.")
        region   = args.region.replace(",", "")
        slug     = sanitize_region(region)
        out_file = args.out or os.path.join(TRACKS_DIR, f"{slug}.svg")
        jobs = [(region, slug, out_file)]

    # 2. Multi-gene mode  (--genes + --flanks)
    elif args.genes:
        if args.flanks is None:
            parser.error("--genes requires --flanks (one value, or one per gene).")
        # Broadcast a single flank value across all genes
        if len(args.flanks) == 1:
            args.flanks = args.flanks * len(args.genes)
        if len(args.genes) != len(args.flanks):
            parser.error(
                f"--genes has {len(args.genes)} entries but --flanks has "
                f"{len(args.flanks)}. Provide either one value (applied to all) "
                f"or one value per gene."
            )
        if args.out:
            parser.error("--out cannot be used with --genes; output paths are auto-generated.")

        jobs = []
        for gene, flank in zip(args.genes, args.flanks):
            region, slug = resolve_region(gene, flank)
            out_file     = os.path.join(TRACKS_DIR, f"{slug}.svg")
            jobs.append((region, slug, out_file))

    # 3. Single-gene mode  (--gene / default)
    else:
        gene  = args.gene  or DEFAULT_GENES[0]
        flank = args.flank or DEFAULT_FLANKS[0]
        region, slug = resolve_region(gene, flank)
        out_file = args.out or os.path.join(TRACKS_DIR, f"{slug}.svg")
        jobs = [(region, slug, out_file)]

    # ── Submit one job per entry ────────────────────────────────────────────────
    print(f"\nSubmitting {len(jobs)} genome-track job(s)...\n")
    for region, slug, out_file in jobs:
        submit_tracks_job(
            region=region,
            slug=slug,
            ini_path=args.ini,
            out_file=out_file,
        )

    print(f"\nAll {len(jobs)} job(s) submitted.")


if __name__ == "__main__":
    main()
