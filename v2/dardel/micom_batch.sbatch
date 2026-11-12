#!/bin/bash -l
# The 277-sample MICOM batch (v2/micom/micom_batch.py), one array task per sample.
#
# Submit from v2/dardel/ AFTER the build job has written $MICOM_SCRATCH/manifest.csv and
# AFTER v2/micom/out/alpha_sweep_choice.json exists (the alpha is read from it):
#
#   python v2/micom/micom_batch.py arrays    # prints the two tiered sbatch lines to run
#
# Sizing (measured 2026-11-12, v2/micom/memory_probe.py, one fresh solve at alpha 0.5):
# ONE process, ONE core; peak 0.84 GB at 20 taxa, 2.22 GB at 61 taxa (0.034 GB/taxon;
# 2.76 GB predicted at the cohort max of 77). The alpha sweep's 6-13 GB workers were
# accumulation over ten solves in one process, which a batch task never does.
# Two memory tiers: `python v2/micom/micom_batch.py arrays` prints both sbatch lines
# (<= 61 taxa at --mem=2650M, billed 3 cores; the 13 larger at 3500M, billed 4).
# Time: the largest test community took 41 min fresh; highs_qp now caps every LP at
# 900 s and the polish at 2400 s, so a stalled sample ends inside 90 min with its
# certified gap recorded. A task that still times out leaves no status.json and is
# re-submitted once by index, then reported as failed by `collect`.
#SBATCH -J ki-micom
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=2650M     # overridden per tier on the sbatch line
#SBATCH -t 01:30:00
#SBATCH -o slurm-ki-micom-%A_%a.out
set -euo pipefail
: "${PDC_STORAGE:=/cfs/klemming/projects/supr/naiss2024-6-169/andnair}"
export MICOM_SCRATCH="${MICOM_SCRATCH:-/cfs/klemming/scratch/a/andnair/micom}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export TMPDIR="${MICOM_SCRATCH}/tmp"; mkdir -p "$TMPDIR"   # not the memory-backed /tmp
source "${PDC_STORAGE}/venv-micom/bin/activate"
cd "${PDC_STORAGE}/ki"
python v2/micom/micom_batch.py solve --index "${SLURM_ARRAY_TASK_ID}"
