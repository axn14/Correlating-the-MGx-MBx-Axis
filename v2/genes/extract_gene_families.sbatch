#!/bin/bash -l
# Subset the stratified gene families to the target enzyme clusters. Submit from
# anywhere on Dardel:
#
#   sbatch -A <project> extract_gene_families.sbatch
#
#SBATCH -J ki-genes
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=12G
#SBATCH -t 01:00:00
#SBATCH -o /cfs/klemming/home/a/andnair/.ki-scratch/genes-%j.out
#
# Why 12G. 8,930,907 x 616 dense doubles would be 41.0 GiB, which is what the row and
# column counts imply and what a first guess would have asked for. Streaming the
# serialization measured 4.58 GiB, 0.112 of that, so the first run asked for 32G -- and
# `sacct` then reported MaxRSS 5,777,524K, 5.5 GiB, on a dgCMatrix holding 379,029,747
# nonzeros at 6.89% density. 12G is roughly twice the measured peak, which leaves room
# for the dense copy of the subset and for R's GC having more slack than it did under a
# 32G request. The measurement cost 0.2 core-hours and took this from roughly 90 cores'
# worth of memory in `shared` -- where memory is allocated at about 0.93 GB per core --
# down to 13. If this is ever OOM-killed, read the reported dim and class from the log
# before raising the number.
#
# Wall time: the first run did everything in 5m26s, of which 288 s was load() alone
# (bzip2 decompression plus deserialization of 3.26 GB). One hour is ample.
#
# One core, because this is a load, a vector match and a write. Nothing here threads.
set -euo pipefail

ROOT="${KI_GENES_ROOT:-/cfs/klemming/projects/supr/naiss2024-6-169/andnair}"
R_MODULE="${R_MODULE:-cray-R/4.4.0}"

echo "== job ${SLURM_JOB_ID:-interactive} on $(hostname) at $(date -u '+%F %T')Z"
echo "   root   : ${ROOT}"
echo "   cores  : ${SLURM_CPUS_PER_TASK:-unset}   mem: ${SLURM_MEM_PER_NODE:-unset} MB"

for f in profiles/gene_families.rda genes/target_clusters.tsv genes/curated_clusters.tsv \
         genes/extract_gene_families.R; do
    [ -f "${ROOT}/${f}" ] || { echo "missing ${ROOT}/${f}" >&2; exit 1; }
done

# Never pipe `ml`: it is a shell function and a pipe runs it in a subshell, which
# silently discards the environment changes.
ml "$R_MODULE"
echo "   R      : $(command -v Rscript)"

export TMPDIR="${PDC_TMP:-/cfs/klemming/scratch/${USER:0:1}/${USER}}/tmp"
mkdir -p "$TMPDIR"

start=$(date +%s)
set +e
Rscript --vanilla "${ROOT}/genes/extract_gene_families.R" "$ROOT"
rc=$?
set -e
echo "== Rscript exited ${rc} after $(( ($(date +%s) - start) / 60 )) min"

if [ "$rc" -eq 0 ]; then
    echo "   outputs in ${ROOT}/genes/out:"
    ls -la "${ROOT}/genes/out"
fi
exit $rc
