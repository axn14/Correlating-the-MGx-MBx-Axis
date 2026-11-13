#!/bin/bash -l
# Class G for the stage-associated targets (CONTRACTS 5d-v). Submit from anywhere on Dardel:
#
#   sbatch -A naiss2026-3-773 extract_stage_clusters.sbatch
#
# Sized from extract_gene_families.sbatch: the same 3.26 GB load measured MaxRSS 5.5 GiB and
# 5 min, so 12G and one hour. One core: a load, a match and a row count. Nothing threads.
#SBATCH -J ki-genes-stage
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=12G
#SBATCH -t 01:00:00
#SBATCH -o slurm-ki-genes-stage-%j.out
set -euo pipefail

ROOT="${KI_GENES_ROOT:-/cfs/klemming/projects/supr/naiss2024-6-169/andnair}"
R_MODULE="${R_MODULE:-cray-R/4.4.0}"
echo "== job ${SLURM_JOB_ID:-interactive} on $(hostname) at $(date -u '+%F %T')Z"

for f in profiles/gene_families.rda mappings/map_level4ec_uniref90.txt.gz \
         genes/stage_enzymes.tsv genes/stage_accessions.txt genes/extract_stage_clusters.R; do
    [ -f "${ROOT}/${f}" ] || { echo "missing ${ROOT}/${f}" >&2; exit 1; }
done

# Never pipe `ml`: it is a shell function and a pipe runs it in a subshell.
ml "$R_MODULE"
export TMPDIR="${PDC_TMP:-/cfs/klemming/scratch/${USER:0:1}/${USER}}/tmp"
mkdir -p "$TMPDIR"

start=$(date +%s)
set +e
Rscript --vanilla "${ROOT}/genes/extract_stage_clusters.R" "$ROOT"
rc=$?
set -e
echo "== Rscript exited ${rc} after $(( ($(date +%s) - start) / 60 )) min"
exit $rc
