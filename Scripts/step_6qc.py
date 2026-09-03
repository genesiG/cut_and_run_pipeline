#!/usr/bin/env python3

"""
step_5_qc_peaks.py - QC steps that require called peaks

Calls each post-peak QC step in order:

  qc_8  → FRiP (Fraction of Reads in Peaks)
  qc_9  → Peak statistics and replicate overlaps

Both steps consume peak files / .bed files as input and must therefore
be run AFTER peak calling (e.g. after step_5 / MACS2).

Usage:
    python3 step_5_qc_peaks.py              # run all steps
    python3 step_5_qc_peaks.py --from qc_9  # start from qc_9
    python3 step_5_qc_peaks.py --only qc_8  # run only one step

Each Python QC step is imported and its main() is called directly.
"""

### Import modules
import os
import sys
import time
import subprocess
import argparse

import config
###

# Ordered list of post-peak QC steps with metadata
# Each entry: (step_key, script_path, is_rscript)
QC_STEPS = [
    ("qc_8", "utils/qc_8_frip.py",        False),
    ("qc_9", "utils/qc_9_peak_stats.py",  False),
]

SCRIPT_DIR = config.CODEDIR


def run_python_step(script_path: str):
    """
    Import and run the main() function of a Python QC script.
    Executes in the same interpreter but catches exceptions gracefully.
    """
    import importlib.util
    spec   = importlib.util.spec_from_file_location("_qc_module", script_path)
    module = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(module)
    except SystemExit:
        pass  # Some scripts may call sys.exit() — ignore

    if hasattr(module, "main"):
        module.main()
    else:
        print(f"  WARNING: {os.path.basename(script_path)} has no main() function.")


def run_rscript_step(script_path: str):
    """
    Call an R script via subprocess, running from the project WORKDIR
    so that py_run_file('Scripts/config.py') resolves correctly.
    """
    cmd = ["Rscript", script_path]
    print(f"  Running: Rscript {os.path.basename(script_path)}")
    result = subprocess.run(cmd, cwd=config.WORKDIR)
    if result.returncode != 0:
        print(f"  WARNING: {os.path.basename(script_path)} exited with "
              f"return code {result.returncode}")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Post-peak QC pipeline — runs qc_8 and qc_9 (peaks required)."
    )
    grp = parser.add_mutually_exclusive_group()
    grp.add_argument(
        "--from", dest="from_step", default=None,
        choices=[s[0] for s in QC_STEPS],
        help="Start from this step (inclusive)."
    )
    grp.add_argument(
        "--only", dest="only_step", default=None,
        choices=[s[0] for s in QC_STEPS],
        help="Run only this step."
    )
    return parser.parse_args()


def main():
    """
    Parses CLI arguments and runs the requested post-peak QC steps in order.
    """
    args  = parse_args()
    steps = QC_STEPS[:]

    # Filter steps based on CLI flags
    if args.only_step:
        steps = [s for s in steps if s[0] == args.only_step]
    elif args.from_step:
        keys  = [s[0] for s in steps]
        start = keys.index(args.from_step)
        steps = steps[start:]

    print("=" * 60)
    print(f"  {config.PROJECT_NAME} — Peak QC (step_5_qc_peaks.py)")
    print(f"  Steps to run: {', '.join(s[0] for s in steps)}")
    print("=" * 60)

    for step_key, script_name, is_rscript in steps:
        script_path = os.path.join(SCRIPT_DIR, script_name)

        if not os.path.exists(script_path):
            print(f"\n[{step_key}] SKIPPED — script not found: {script_path}")
            continue

        print(f"\n{'='*60}")
        print(f"[{step_key}] Running: {script_name}")
        print(f"{'='*60}")

        t0 = time.time()
        try:
            if is_rscript:
                run_rscript_step(script_path)
            else:
                run_python_step(script_path)
        except Exception as exc:
            print(f"[{step_key}] ERROR: {exc}")
        elapsed = time.time() - t0
        print(f"[{step_key}] Completed in {elapsed:.1f}s")

    print("\n" + "=" * 60)
    print("  All peak QC steps submitted / completed.")
    print("  NOTE: Steps that submit bsub jobs may still be running.")
    print("        Check job status with: bjobs -u $USER")
    print("=" * 60)


if __name__ == "__main__":
    main()
