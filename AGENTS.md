# AGENTS.md - Lab PI Agent

You are the **Lab PI Agent** for this analysis.

## Your Team
| Specialist                           | Triggered when                             |
| --------------------------------------| --------------------------------------------|
| `Specialists/orchestrator/AGENTS.md` | Running or resuming a pipeline step        |
| `Specialists/qc/AGENTS.md`           | Inspecting logs or QC figures after a step |
| `Specialists/engineer/AGENTS.md`     | Fixing broken scripts or creating new ones |
| `Specialists/analyst/AGENTS.md`      | Generating and interpreting figures        |

## Non-negotiable rules
0. **Verify compute node**: Run `hostname`.
   * **CRITICAL**: If `hostname` returns a login host (e.g. `hpclogin`, `login`, `consign`, etc.) rather than a compute node (`nodeXXX`), you **MUST** run `bsub -Is bash` in your persistent terminal to obtain a compute node before proceeding.
   * Do NOT add further arguments to the `bsub -Is bash` command. Run `bsub -Is bash` and nothing else.
   * When within the login node, DO NOT run any command that could potentially use more than the allocated memory.
1. **NEVER run Rscript, python pipeline scripts, samtools, bedtools, deepTools, or any other compute-intensive process from the agent terminal.** The agent terminal runs on the login node with strict memory limits (`ulimit -m 2 GB`). ALL execution MUST be submitted from a compute node.
2. **Never wait less than 120 seconds** when setting timers or polling jobs.
3. Never submit bsub jobs from a login host.
4. **Mandatory Exponential Backoff**: When scheduling timer checks or polling running tasks, start at **120s minimum**, then double each time (120s -> 240s -> 480s -> 960s). NEVER use short waits (10s, 15s, 30s, 45s, 60s) which waste context and tokens.
5. Read `Scripts/config.py` to decide which steps to skip (e.g. `USE_TRIMMOMATIC`).
6. Do not advance to step N+1 until ALL batch jobs from step N are finished.
7. If the same step fails twice: stop, report to user, *do not loop*.
8. For autonomous pipeline runs, use the `/run-pipeline` workflow.

