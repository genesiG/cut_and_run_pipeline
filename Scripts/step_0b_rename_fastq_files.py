#!/usr/bin/env python3

"""
step_0b_rename_fastq_files.py - Strip sequencing-core suffixes from fastq filenames.

Naming convention handled (standard Illumina sequencing-core naming, e.g. bcl2fastq/Illumina DRAGEN):
    <SampleName>_S<1-2 digits>_L<3 digits>_R<1|2>_001.fastq.gz
                 ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                 stripped portions

Result:
    <SampleName>_R<1|2>.fastq.gz

Examples
--------
  BT54_K27me3_B_S35_L006_R2_001.fastq.gz  →  BT54_K27me3_B_R2.fastq.gz
  BT88_K27me3_A_S26_L006_R1_001.fastq.gz  →  BT88_K27me3_A_R1.fastq.gz

Usage
-----
  python3 Scripts/step_0b_rename_fastq_files.py [--dry-run]

  --dry-run   Print what would be renamed without touching the files.
"""

import argparse
import os
import re
import sys

# ---------------------------------------------------------------------------
# Import project config
# ---------------------------------------------------------------------------
sys.path.insert(0, os.path.dirname(__file__))
import config

# ---------------------------------------------------------------------------
# Regex that matches the core-specific suffix to be stripped.
#
#   _S<1-2 digits>   – sample index assigned by the sequencer
#   _L<3 digits>     – lane number
#   _001             – read-chunk suffix (always 001 for a single merged file)
#
# The _R<1|2> part is intentionally NOT matched so it is preserved.
# ---------------------------------------------------------------------------
def strip_core_suffixes(filename: str) -> str:
    """Return filename with sequencing-core suffixes removed, or original if no match."""
    new_name = filename
    
    # Remove _500K
    new_name = re.sub(r"_500K", "", new_name)
    
    # Remove _Rpcells_Sxx
    new_name = re.sub(r"_Rpcells_S\d{1,2}", "", new_name)
    
    # Change .R1 / .R2 to _R1 / _R2 for config.py pattern matching
    new_name = new_name.replace(".R1.fastq", "_R1.fastq")
    new_name = new_name.replace(".R2.fastq", "_R2.fastq")
    
    # Standard cleanup just in case
    new_name = re.sub(r"_001(?=\.fastq\.gz$|\.fastq$)", "", new_name)
    new_name = re.sub(r"_S\d{1,2}_L\d{3}(?=_R[12])", "", new_name)

    return new_name


def rename_file(src_path: str, dst_path: str, dry_run: bool) -> None:
    if dry_run:
        print(f"  [DRY-RUN]  {os.path.basename(src_path)}  →  {os.path.basename(dst_path)}")
    else:
        os.rename(src_path, dst_path)
        print(f"  Renamed:   {os.path.basename(src_path)}  →  {os.path.basename(dst_path)}")



def main() -> None:
    parser = argparse.ArgumentParser(
        description="Strip sequencing-core suffixes (_S##_L###, _001) from fastq filenames."
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be renamed without making any changes.",
    )
    args = parser.parse_args()

    data_dir = config.DATADIR          # e.g. Original_Data/fastq_files
    dry_run  = args.dry_run

    if not os.path.isdir(data_dir):
        print(f"ERROR: Data directory not found: {data_dir}")
        sys.exit(1)

    fastq_files = [
        f for f in os.listdir(data_dir)
        if f.endswith(".fastq.gz") or f.endswith(".fastq")
    ]

    if not fastq_files:
        print(f"No fastq files found in {data_dir}")
        sys.exit(0)

    rename_count  = 0
    skip_count    = 0
    conflict_count = 0

    print(f"\n{'[DRY-RUN] ' if dry_run else ''}Scanning {len(fastq_files)} file(s) in {data_dir}\n")

    for fname in sorted(fastq_files):
        new_fname = strip_core_suffixes(fname)

        if new_fname == fname:
            # Nothing to strip – already clean or unexpected format
            print(f"  Skipped (no match):  {fname}")
            skip_count += 1
            continue

        src = os.path.join(data_dir, fname)
        dst = os.path.join(data_dir, new_fname)

        if os.path.exists(dst) and not dry_run:
            print(f"  CONFLICT (target exists, skipping):  {fname}  →  {new_fname}")
            conflict_count += 1
            continue

        rename_file(src, dst, dry_run)
        rename_count += 1

    print(
        f"\n{'[DRY-RUN] ' if dry_run else ''}Done. "
        f"{rename_count} renamed | {skip_count} skipped | {conflict_count} conflicts."
    )


if __name__ == "__main__":
    main()
