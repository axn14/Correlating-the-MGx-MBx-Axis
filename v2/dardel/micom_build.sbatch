#!/bin/bash -l
# Build the 277 community models for the MICOM batch, onto scratch (not project storage).
#
#   sbatch -A naiss2026-3-773 micom_build.sbatch
#
# Needs Data/micom/agora2_gtdb_species.qza (0.23 GB) staged next to the code, and
# v2/micom/out/batch/taxonomy.tsv (from `micom_batch.py prepare` on the laptop).
#SBATCH -J ki-micom-build
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=7000M
#SBATCH -t 06:00:00
#SBATCH -o slurm-ki-micom-build-%j.out
set -euo pipefail
: "${PDC_STORAGE:=/cfs/klemming/projects/supr/naiss2024-6-169/andnair}"
export MICOM_SCRATCH="${MICOM_SCRATCH:-/cfs/klemming/scratch/a/andnair/micom}"
source "${PDC_STORAGE}/venv-micom/bin/activate"
cd "${PDC_STORAGE}/ki"
python v2/micom/micom_batch.py build --model-db Data/micom/agora2_gtdb_species.qza --threads 8
