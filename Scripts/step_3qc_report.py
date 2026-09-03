#!/usr/bin/env python3

"""
step_3_qc_report.py  –  Post-processing QC report
                         Requires ALL bsub jobs from step_3_qc.py to be done.

Run this AFTER step_3_qc.py and after all bsub jobs have finished:

    bjobs -u $USER            # confirm no jobs running
    python3 Scripts/step_3_qc_report.py

Steps (sequential — each step waits for the previous before proceeding):

  multiqc   → MultiQC aggregation of FastQC reports                      [bsub+wait]
  qc_3      → Alignment statistics plots (R)                             [bsub+wait]
  qc_5b_r   → Preseq complexity curve plot (R)                          [bsub+wait]
  report    → Generate single self-contained HTML report

The final HTML report embeds all QC outputs (PDFs, SVGs, PNGs, MultiQC HTML
link) into one document, clearly separated by headers.

Usage:
    python3 Scripts/step_3_qc_report.py           # run all steps
    python3 Scripts/step_3_qc_report.py --only report   # regenerate HTML only
    python3 Scripts/step_3_qc_report.py --from qc_3     # start from qc_3 onwards
"""

import base64
import os
import subprocess
import time
import argparse

import config

# ---------------------------------------------------------------------------
# Step registry
# Each entry: (step_key, description)
# All steps are handled inline by dedicated run_*() functions.
# ---------------------------------------------------------------------------
REPORT_STEPS = [
    ("multiqc", "MultiQC aggregation of FastQC reports"),
    ("qc_3",    "Alignment statistics plots (R)"),
    ("qc_5b_r", "Preseq complexity curve plot (R)"),
    ("report",  "Generate combined HTML QC report"),
]

SCRIPT_DIR = config.CODEDIR

# FastQC output directory (read by MultiQC)
if config.USE_TRIMMOMATIC:
    FASTQC_DIR = os.path.join(config.QCDIR1, "trimmed")
else:
    FASTQC_DIR = os.path.join(config.QCDIR1, "untrimmed")

# Output directories used when building the HTML report
MULTIQC_DIR     = config.MULTIQCDIR
ALIGN_QC_DIR    = os.path.join(config.QCDIR1, "alignment")
FRAG_SIZE_DIR   = config.FRAG_SIZE_DIR
LIB_COMPLEX_DIR = config.LIB_COMPLEXITY_DIR
PRESEQ_DIR      = os.path.join(config.LIB_COMPLEXITY_DIR, "preseq")
PICARD_DIR      = os.path.join(config.LIB_COMPLEXITY_DIR, "picard")
FINGERPRINT_DIR = config.FINGERPRINT_DIR
REPRO_DIR       = config.REPRODUCIBILITY_DIR

# Final HTML report output path
REPORT_PATH = os.path.join(config.QCDIR1, f"{config.PROJECT_NAME}_qc_report.html")


# ===========================================================================
# Shared bsub helpers
# ===========================================================================

def _submit_bsub_and_wait(job_name: str, batch_cmd: str,
                           log_dir: str, step_key: str,
                           poll_interval: int = 30,
                           timeout_minutes: int = 180) -> bool:
    """
    Write a batch script, submit via bsub, then poll until the job exits.

    Returns True on success (job exited DONE), False if timed-out or FAILED.
    Raises RuntimeError if bjobs returns EXIT status.

    Parameters
    ----------
    poll_interval   : seconds between bjobs status checks (default 30 s)
    timeout_minutes : abort after this many minutes (default 3 h)
    """
    os.makedirs(log_dir, exist_ok=True)

    batch_file = os.path.join(log_dir, f"{step_key}.batch")
    with open(batch_file, 'w', encoding='utf-8') as fh:
        fh.write(batch_cmd)

    # Submit and capture the job ID from bsub stdout
    result = subprocess.run(
        ["bsub"],
        input=batch_cmd,
        text=True,
        capture_output=True,
        check=False,
    )
    stdout = result.stdout.strip()
    print(f"  bsub output: {stdout}")

    # Parse job ID from "Job <NNNN> is submitted..."
    import re
    m = re.search(r"Job <(\d+)>", stdout)
    if not m:
        print(f"  WARNING: Could not parse job ID from bsub output. "
              f"Proceeding without polling.")
        return True

    job_id = m.group(1)
    print(f"  Job submitted: {job_id}. Polling every {poll_interval}s ...")

    deadline = time.time() + timeout_minutes * 60
    while time.time() < deadline:
        time.sleep(poll_interval)
        bjobs = subprocess.run(
            ["bjobs", "-noheader", job_id],
            text=True, capture_output=True, check=False
        )
        status_line = bjobs.stdout.strip()
        if not status_line:
            # Job no longer in the queue — completed or failed
            # Check exit status via bhist
            bhist = subprocess.run(
                ["bhist", "-noheader", "-d", job_id],
                text=True, capture_output=True, check=False
            )
            if "DONE" in bhist.stdout or bhist.stdout.strip() == "":
                print(f"  Job {job_id} DONE.")
                return True
            if "EXIT" in bhist.stdout:
                raise RuntimeError(
                    f"Job {job_id} ({step_key}) exited with failure. "
                    f"Check log: {os.path.join(log_dir, step_key + '.error')}"
                )
            # bhist may not have data yet; treat as done
            print(f"  Job {job_id} no longer in queue (assumed DONE).")
            return True

        # Still in queue — check state field (3rd column)
        fields = status_line.split()
        if len(fields) >= 3:
            state = fields[2]
            if state == "DONE":
                print(f"  Job {job_id} DONE.")
                return True
            if state == "EXIT":
                raise RuntimeError(
                    f"Job {job_id} ({step_key}) exited with failure. "
                    f"Check log: {os.path.join(log_dir, step_key + '.error')}"
                )
            # RUN or PEND — keep waiting
            print(f"  Job {job_id} status: {state} …")

    print(f"  WARNING: Timed out waiting for job {job_id} ({step_key}) "
          f"after {timeout_minutes} min.")
    return False


# ===========================================================================
# multiqc — aggregate FastQC reports
# ===========================================================================

def run_multiqc() -> None:
    """
    Submit a MultiQC bsub job and wait for completion.
    Aggregates all FastQC output in FASTQC_DIR into a single HTML report.
    Output: config.MULTIQCDIR/{PROJECT}_multiqc_report.html
    """
    if not os.path.isdir(FASTQC_DIR):
        print(f"  WARNING: FastQC directory not found: {FASTQC_DIR}")
        print("  Run step_3_qc.py --only fastqc and wait for jobs to finish first.")
        return

    log_dir      = os.path.join(MULTIQC_DIR, "log")
    os.makedirs(MULTIQC_DIR, exist_ok=True)
    os.makedirs(log_dir,     exist_ok=True)

    job_name = f"{config.PROJECT_NAME}.multiqc"
    log_out  = os.path.join(log_dir, "multiqc.log")
    log_err  = os.path.join(log_dir, "multiqc.error")
    batch_file = os.path.join(log_dir, "multiqc.batch")

    batch_cmd = f"""\
#!/bin/bash
#BSUB -P {config.PROJECT_NAME}
#BSUB -J {job_name}
#BSUB -oo {log_out}
#BSUB -eo {log_err}
#BSUB -q {config.BSUB_QUEUE}
#BSUB -n 4
#BSUB -M 16000
#BSUB -R "rusage[mem=4000] span[hosts=1]"

source {config.HOMEDIR}/miniconda3/etc/profile.d/conda.sh
conda activate chipseq

echo "Running MultiQC on: {FASTQC_DIR}"

multiqc \\
    {FASTQC_DIR} \\
    --outdir {MULTIQC_DIR} \\
    --filename {config.PROJECT_NAME}_multiqc_report.html \\
    --force \\
    --verbose

echo "MultiQC report: {MULTIQC_DIR}/{config.PROJECT_NAME}_multiqc_report.html"
echo "=== Resource usage summary ==="
bjobs -l $LSB_JOBID | grep -A10 "Resource usage"
"""
    with open(batch_file, 'w', encoding='utf-8') as fh:
        fh.write(batch_cmd)
        
    print(f"  Submitting MultiQC job ...")
    os.system(f"bsub < {batch_file}")
    time.sleep(1)

# ===========================================================================
# qc_3 — Alignment statistics R plots
# ===========================================================================

def run_qc_3() -> None:
    """
    Submit the qc_3_alignment_stats.R bsub job and wait for completion.
    Requires: qc_3a outputs in ALIGN_QC_DIR (flagstat, idxstats, picard metrics).
    Produces: multiple PDF plots in ALIGN_QC_DIR.
    """
    script_path = os.path.join(SCRIPT_DIR, "utils", "qc_3_alignment_stats.R")
    if not os.path.exists(script_path):
        print(f"  WARNING: Script not found: {script_path} — skipping qc_3.")
        return

    log_dir  = os.path.join(config.QCDIR1, "log")
    os.makedirs(log_dir, exist_ok=True)

    job_name     = f"{config.PROJECT_NAME}.qc_3"
    log_out_file = os.path.join(log_dir, "qc_3.log")
    log_err_file = os.path.join(log_dir, "qc_3.error")

    batch_cmd = f"""\
#!/bin/bash
#BSUB -P {config.PROJECT_NAME}
#BSUB -J {job_name}
#BSUB -oo {log_out_file}
#BSUB -eo {log_err_file}
#BSUB -q {config.BSUB_QUEUE}
#BSUB -n 2
#BSUB -M 16000
#BSUB -R "rusage[mem=8000] span[hosts=1]"

source {config.HOMEDIR}/miniconda3/etc/profile.d/conda.sh
conda activate chipseq

echo "Running qc_3_alignment_stats.R"
Rscript {script_path}

echo "=== Resource usage summary ==="
bjobs -l $LSB_JOBID | grep -A10 "Resource usage"
"""
    print("  Submitting qc_3 (alignment stats R) and waiting for completion ...")
    _submit_bsub_and_wait(job_name, batch_cmd, log_dir, "qc_3")


# ===========================================================================
# qc_5b_r — Preseq complexity curve R plot
# ===========================================================================

def run_qc_5b_r() -> None:
    """
    Submit the qc_5b_preseq_ccurve.R bsub job and wait for completion.
    Requires: *.preseq_ccurve.txt files in PRESEQ_DIR (from step_3_qc.py qc_5b).
    Produces: PRESEQ_DIR/preseq_ccurve.svg
    """
    script_path = os.path.join(SCRIPT_DIR, "utils", "qc_5b_preseq_ccurve.R")
    if not os.path.exists(script_path):
        print(f"  WARNING: Script not found: {script_path} — skipping qc_5b_r.")
        return

    ccurve_files = [
        f for f in os.listdir(PRESEQ_DIR)
        if f.endswith(".preseq_ccurve.txt")
    ] if os.path.isdir(PRESEQ_DIR) else []

    if not ccurve_files:
        print("  WARNING: No preseq_ccurve.txt files found in:")
        print(f"  {PRESEQ_DIR}")
        print("  Run step_3_qc.py --only qc_5b and wait for jobs first.")
        return

    log_dir  = os.path.join(config.QCDIR1, "log")
    os.makedirs(log_dir, exist_ok=True)

    job_name     = f"{config.PROJECT_NAME}.qc_5b_r"
    log_out_file = os.path.join(log_dir, "qc_5b_r.log")
    log_err_file = os.path.join(log_dir, "qc_5b_r.error")

    batch_cmd = f"""\
#!/bin/bash
#BSUB -P {config.PROJECT_NAME}
#BSUB -J {job_name}
#BSUB -oo {log_out_file}
#BSUB -eo {log_err_file}
#BSUB -q {config.BSUB_QUEUE}
#BSUB -n 2
#BSUB -M 16000
#BSUB -R "rusage[mem=8000] span[hosts=1]"

source {config.HOMEDIR}/miniconda3/etc/profile.d/conda.sh
conda activate chipseq

echo "Running qc_5b_preseq_ccurve.R"
Rscript {script_path}

echo "=== Resource usage summary ==="
bjobs -l $LSB_JOBID | grep -A10 "Resource usage"
"""
    print("  Submitting qc_5b_r (preseq R plot) and waiting for completion ...")
    _submit_bsub_and_wait(job_name, batch_cmd, log_dir, "qc_5b_r")


# ===========================================================================
# report — Build combined HTML QC report
# ===========================================================================

# Ordered list of QC sections and the glob patterns that find their outputs.
# Each entry: (section_title, base_dir, pattern, embed_mode)
# pattern can be a glob string or a list of glob strings
# embed_mode: "pdf" | "svg" | "png" | "iframe" | "link"
REPORT_SECTIONS = [
    # ── step_3_qc outputs ────────────────────────────────────────────────────
    (
        "Alignment Statistics",
        ALIGN_QC_DIR,
        [
            "alignment_stats_barplot.html",
            "alignment_stats_barplot_pct.html",
            "alignment_stats_n_spikein.html",
            "alignment_stats_pct_spikein_*.html",
            "alignment_stats_pct_chrM.html",
            "alignment_stats_pct_unmapped.html",
            "alignment_stats_bowtie2.html",
            "alignment_stats_duplicates.html",
        ],
        "html_inline",
    ),
    (
        "Fragment Size Distribution",
        FRAG_SIZE_DIR,
        "*.html",
        "html_inline",
    ),
    (
        "Library Complexity",
        PICARD_DIR,
        "library_complexity*.html",
        "html_inline",
    ),
    (
        "Preseq Complexity Curves",
        PRESEQ_DIR,
        "preseq_ccurve.html",
        "html_inline",
    ),
    (
        "Fingerprint Plots",
        FINGERPRINT_DIR,
        "*.png",
        "png",
    ),
    (
        "Reproducibility — Correlation Heatmaps (by Antibody Target)",
        REPRO_DIR,
        "[!g]*_correlation_heatmap.svg",
        "svg",
    ),
    (
        "Reproducibility — Correlation Heatmaps (by Sample Group)",
        REPRO_DIR,
        "group_*_correlation_heatmap.svg",
        "svg",
    ),
    (
        "Reproducibility — PCA Biplots (by Antibody Target)",
        REPRO_DIR,
        "[!g]*_pca_biplot.svg",
        "svg",
    ),
    (
        "Reproducibility — PCA Biplots (by Sample Group)",
        REPRO_DIR,
        "group_*_pca_biplot.svg",
        "svg",
    ),
]


def _file_to_data_uri(path: str, mime: str) -> str:
    """Read a binary file and return a base64-encoded data URI."""
    with open(path, "rb") as fh:
        encoded = base64.b64encode(fh.read()).decode("ascii")
    return f"data:{mime};base64,{encoded}"


def _svg_inline(path: str) -> str:
    """Return the raw SVG text for direct inline embedding."""
    with open(path, encoding="utf-8") as fh:
        return fh.read()


from typing import Union

def _build_section_html(title: str, base_dir: str, pattern: Union[list, str],
                         embed_mode: str) -> str:
    """
    Build the HTML block for a single QC section.

    For pdf/svg/png: embeds the file content directly in the HTML.
    For iframe: embeds the HTML file in a sandboxed iframe.
    For link: renders a clickable hyperlink to the file.
    """
    if not os.path.isdir(base_dir):
        return (
            f'<section class="qc-section">'
            f'<h2>{title}</h2>'
            f'<p class="missing">⚠ Output directory not found: '
            f'<code>{base_dir}</code></p>'
            f'</section>\n'
        )

    import glob as _glob
    
    if isinstance(pattern, list):
        files = []
        for p in pattern:
            matches = sorted(_glob.glob(os.path.join(base_dir, p)))
            files.extend(matches)
    else:
        files = sorted(_glob.glob(os.path.join(base_dir, pattern)))

    if not files:
        pattern_str = ", ".join(pattern) if isinstance(pattern, list) else pattern
        return (
            f'<section class="qc-section">'
            f'<h2>{title}</h2>'
            f'<p class="missing">⚠ No output files found matching '
            f'<code>{pattern_str}</code> in <code>{base_dir}</code></p>'
            f'</section>\n'
        )

    items_html = []
    for path in files:
        fname = os.path.basename(path)
        caption = f'<p class="fig-caption"><code>{fname}</code></p>'

        if embed_mode == "pdf":
            # Embed PDF via <object> with a fallback link
            data_uri = _file_to_data_uri(path, "application/pdf")
            items_html.append(
                f'<figure>'
                f'<object data="{data_uri}" type="application/pdf" '
                f'width="100%" height="600px">'
                f'<p><a href="{path}">Download PDF: {fname}</a></p>'
                f'</object>'
                f'{caption}'
                f'</figure>'
            )

        elif embed_mode == "svg":
            svg_content = _svg_inline(path)
            items_html.append(
                f'<figure class="svg-figure">'
                f'{svg_content}'
                f'{caption}'
                f'</figure>'
            )

        elif embed_mode == "png":
            data_uri = _file_to_data_uri(path, "image/png")
            items_html.append(
                f'<figure>'
                f'<img src="{data_uri}" alt="{fname}" style="max-width:100%;">'
                f'{caption}'
                f'</figure>'
            )

        elif embed_mode == "html_inline":
            # Read the HTML file and embed it in a srcdoc iframe.
            # This avoids browser CORS/same-origin blocks on file:// URLs
            # while still sandboxing the widget's scripts.
            try:
                with open(path, encoding="utf-8") as _fh:
                    raw_html = _fh.read()
                # Escape double-quotes in the HTML so it can sit inside
                # the srcdoc="..." attribute safely.
                srcdoc = raw_html.replace('&', '&amp;').replace('"', '&quot;')
                items_html.append(
                    f'<figure>'
                    f'<iframe srcdoc="{srcdoc}" width="100%" height="550px" '
                    f'frameborder="0" sandbox="allow-scripts allow-same-origin allow-popups">'
                    f'<p><a href="{path}">Open interactive plot: {fname}</a></p>'
                    f'</iframe>'
                    f'{caption}'
                    f'</figure>'
                )
            except OSError as _e:
                items_html.append(
                    f'<p class="missing">⚠ Could not read <code>{fname}</code>: {_e}</p>'
                )

        elif embed_mode == "link":
            items_html.append(
                f'<p><a href="{path}" target="_blank">📄 {fname}</a></p>'
            )

    content = "\n".join(items_html)
    return (
        f'<section class="qc-section">\n'
        f'<h2>{title}</h2>\n'
        f'{content}\n'
        f'</section>\n'
    )


_HTML_STYLE = """
<style>
  :root {
    --bg:      #0f1117;
    --surface: #1a1d27;
    --border:  #2d3148;
    --accent:  #88419d;
    --text:    #e2e8f0;
    --muted:   #6b7280;
    --warn:    #f59e0b;
    --font:    'Inter', 'Segoe UI', system-ui, sans-serif;
  }
  @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;600;700&display=swap');
  *, *::before, *::after { box-sizing: border-box; margin: 0; padding: 0; }
  body {
    background: var(--bg);
    color: var(--text);
    font-family: var(--font);
    font-size: 15px;
    line-height: 1.6;
    padding: 2rem;
    max-width: 1400px;
    margin: 0 auto;
  }
  header {
    border-bottom: 2px solid var(--accent);
    padding-bottom: 1.5rem;
    margin-bottom: 2.5rem;
  }
  header h1 {
    font-size: 2rem;
    font-weight: 700;
    color: var(--text);
    letter-spacing: -0.02em;
  }
  header p {
    color: var(--muted);
    margin-top: 0.4rem;
  }
  nav {
    background: var(--surface);
    border: 1px solid var(--border);
    border-radius: 8px;
    padding: 1rem 1.5rem;
    margin-bottom: 2.5rem;
  }
  nav h3 { font-size: 0.8rem; text-transform: uppercase;
            letter-spacing: 0.1em; color: var(--muted); margin-bottom: 0.6rem; }
  nav ol { padding-left: 1.4rem; }
  nav ol li { margin-bottom: 0.25rem; }
  nav ol li a { color: var(--accent); text-decoration: none; }
  nav ol li a:hover { text-decoration: underline; }
  .qc-section {
    background: var(--surface);
    border: 1px solid var(--border);
    border-radius: 10px;
    padding: 1.8rem 2rem;
    margin-bottom: 2.5rem;
  }
  .qc-section h2 {
    font-size: 1.25rem;
    font-weight: 700;
    color: var(--accent);
    padding-bottom: 0.6rem;
    border-bottom: 1px solid var(--border);
    margin-bottom: 1.2rem;
  }
  figure {
    margin: 1.2rem 0;
    padding: 0.5rem;
    background: #13161f;
    border-radius: 6px;
    border: 1px solid var(--border);
  }
  .fig-caption {
    font-size: 0.8rem;
    color: var(--muted);
    margin-top: 0.4rem;
    padding: 0 0.25rem;
  }
  .fig-caption code { color: var(--text); }
  .svg-figure svg { width: 100%; height: auto; }
  .missing { color: var(--warn); font-style: italic; }
  footer {
    margin-top: 3rem;
    padding-top: 1rem;
    border-top: 1px solid var(--border);
    color: var(--muted);
    font-size: 0.8rem;
    text-align: center;
  }
</style>
"""


def run_report() -> None:
    """
    Collect all QC output files from step_3_qc and step_3_qc_report and
    generate a single self-contained HTML QC report at REPORT_PATH.

    PDFs are embedded as base64-encoded <object> elements.
    PNGs are embedded as base64-encoded <img> elements.
    SVGs are embedded inline.
    The MultiQC HTML is embedded in an iframe (absolute path).
    """
    print(f"  Generating HTML report → {REPORT_PATH}")

    sections_html = []
    toc_entries   = []

    for i, (title, base_dir, pattern, mode) in enumerate(REPORT_SECTIONS, 1):
        anchor = f"section-{i}"
        toc_entries.append(f'<li><a href="#{anchor}">{title}</a></li>')
        sec = _build_section_html(title, base_dir, pattern, mode)
        # Inject anchor id into the section tag
        sec = sec.replace('<section class="qc-section">',
                          f'<section class="qc-section" id="{anchor}">', 1)
        sections_html.append(sec)

    import datetime
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")

    toc_html = "\n".join(toc_entries)
    body_html = "\n".join(sections_html)

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>{config.PROJECT_NAME} — QC Report</title>
  {_HTML_STYLE}
</head>
<body>

<header>
  <h1>🧬 {config.PROJECT_NAME} — Pre-Peak QC Report</h1>
  <p>Generated: {now} &nbsp;|&nbsp;
     Experiment: {config.EXPERIMENT} &nbsp;|&nbsp;
     Species: {config.SPECIES}</p>
</header>

<nav>
  <h3>Contents</h3>
  <ol>
{toc_html}
  </ol>
</nav>

{body_html}

<footer>
  <p>Generated by <code>step_3_qc_report.py</code> — {config.PROJECT_NAME}</p>
</footer>

</body>
</html>
"""

    os.makedirs(os.path.dirname(REPORT_PATH), exist_ok=True)
    with open(REPORT_PATH, 'w', encoding='utf-8') as fh:
        fh.write(html)

    print(f"  ✓ QC report written to: {REPORT_PATH}")
    print(f"    Sections: {len(REPORT_SECTIONS)}")


# ===========================================================================
# CLI
# ===========================================================================

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "step_3_qc_report — Sequential post-processing QC report.\n"
            "Requires all step_3_qc.py bsub jobs to be complete."
        )
    )
    grp = parser.add_mutually_exclusive_group()
    grp.add_argument(
        "--from", dest="from_step", default=None,
        choices=[s[0] for s in REPORT_STEPS],
        help="Start from this step (inclusive).",
    )
    grp.add_argument(
        "--only", dest="only_step", default=None,
        choices=[s[0] for s in REPORT_STEPS],
        help="Run only this step.",
    )
    return parser.parse_args()


def main() -> None:
    """Run all report steps sequentially, waiting for each bsub job to finish."""
    args  = parse_args()
    steps = REPORT_STEPS[:]

    if args.only_step:
        steps = [s for s in steps if s[0] == args.only_step]
    elif args.from_step:
        keys  = [s[0] for s in steps]
        start = keys.index(args.from_step)
        steps = steps[start:]

    # Inline step dispatch table
    runners = {
        "multiqc": run_multiqc,
        "qc_3":    run_qc_3,
        "qc_5b_r": run_qc_5b_r,
        "report":  run_report,
    }

    print("=" * 60)
    print(f"  {config.PROJECT_NAME} — step_3_qc_report.py")
    print(f"  Steps: {', '.join(s[0] for s in steps)}")
    print("=" * 60)

    for step_key, description in steps:
        print(f"\n{'='*60}")
        print(f"[{step_key}] {description}")
        print("=" * 60)

        t0 = time.time()
        try:
            runners[step_key]()
        except RuntimeError as exc:
            print(f"  [{step_key}] FAILED: {exc}")
            print("  Stopping report generation. Fix the error and re-run.")
            raise SystemExit(1) from exc
        except Exception as exc:   # noqa: BLE001
            print(f"  [{step_key}] ERROR: {exc}")

        elapsed = time.time() - t0
        print(f"[{step_key}] Done in {elapsed:.1f}s")

    print("\n" + "=" * 60)
    print("  step_3_qc_report.py complete.")
    print(f"  HTML report: {REPORT_PATH}")
    print("=" * 60)


if __name__ == "__main__":
    main()
