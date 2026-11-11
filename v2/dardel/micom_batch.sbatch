#!/bin/bash -l
# The 277-sample MICOM batch (v2/micom/micom_batch.py), one array task per sample.
#
# Submit from v2/dardel/ AFTER the build job has written $MICOM_SCRATCH/manifest.csv and
# AFTER v2/micom/out/alpha_sweep_choice.json exists (the alpha is read from it):
#
#   sbatch -A naiss2026-3-773 --array=0-276%32 micom_batch.sbatch
#
# Sizing (measured on the laptop, see v2/micom/README.md finding 7 and the FVA benchmark):
# one sample is ONE process using ONE core, peak ~2.4 GB. The shared partition bills
# max(cores, mem / 888 MB), so --mem=2650M is billed as 3 cores; 3 GB would bill 4.
# 90 min covers the slowest sample measured (57 min). A task that times out leaves no
# status.json and is re-submitted once by index, then reported as failed by `collect`.
#SBATCH -J ki-micom
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=2650M
#SBATCH -t 01:30:00
#SBATCH -o slurm-ki-micom-%A_%a.out
set -euo pipefail
: "${PDC_STORAGE:=/cfs/klemming/projects/supr/naiss2024-6-169/andnair}"
export MICOM_SCRATCH="${MICOM_SCRATCH:-/cfs/klemming/scratch/a/andnair/micom}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
source "${PDC_STORAGE}/venv-micom/bin/activate"
cd "${PDC_STORAGE}/ki"
python v2/micom/micom_batch.py solve --index "${SLURM_ARRAY_TASK_ID}"
