#!/bin/bash -l
# Build the 277 community models for the MICOM batch, onto scratch (not project storage).
#
#   sbatch -A naiss2026-3-773 micom_build.sbatch
#
# Needs Data/micom/agora2_gtdb_species.qza (0.23 GB) staged next to the code, and
# v2/micom/out/batch/taxonomy.tsv (from `micom_batch.py prepare` on the laptop).
#
# 2026-11-12, after job 27465147 died OUT_OF_MEMORY in 4 min at 22 GB (11 oom kills):
# given a .qza, micom's build_and_save unzips the WHOLE 1.7 GB database into a temp dir
# in EVERY worker, and on Dardel's compute nodes the default /tmp is memory-backed, so
# 8 workers x 1.7 GB of unzipped JSON landed in the job's RAM on top of each build's own
# ~2.4 GB (measured locally, 61 taxa). Now: the database is extracted ONCE to scratch
# and the build points at the extracted folder (what the successful local builds did),
# and TMPDIR is on scratch as in run_pipeline.sbatch. Workers = cores requested.
# Memory: 6 workers x ~2.4 GB (3 GB at the 77-taxon maximum) + the parent -> 20 GB.
#SBATCH -J ki-micom-build
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 6
#SBATCH --mem=20000M
#SBATCH -t 06:00:00
#SBATCH -o slurm-ki-micom-build-%j.out
set -euo pipefail
: "${PDC_STORAGE:=/cfs/klemming/projects/supr/naiss2024-6-169/andnair}"
export MICOM_SCRATCH="${MICOM_SCRATCH:-/cfs/klemming/scratch/a/andnair/micom}"
export TMPDIR="${MICOM_SCRATCH}/tmp"
mkdir -p "$TMPDIR"
source "${PDC_STORAGE}/venv-micom/bin/activate"
cd "${PDC_STORAGE}/ki"

DB_DIR="${MICOM_SCRATCH}/agora2_extracted"
if ! find "$DB_DIR" -name manifest.csv -print -quit 2>/dev/null | grep -q .; then
    echo "extracting the AGORA2 database once to ${DB_DIR}"
    python - "$DB_DIR" <<'EOF'
import sys, zipfile
zipfile.ZipFile("Data/micom/agora2_gtdb_species.qza").extractall(sys.argv[1])
EOF
fi
MODEL_DB="$(dirname "$(find "$DB_DIR" -name manifest.csv | head -1)")"
echo "model db: ${MODEL_DB} ($(ls "$MODEL_DB" | wc -l) files)"
python v2/micom/micom_batch.py build --model-db "$MODEL_DB" --threads "${SLURM_CPUS_PER_TASK:-6}"
