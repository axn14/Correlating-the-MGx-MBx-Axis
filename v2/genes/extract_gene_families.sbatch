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
#SBATCH --mem=32G
#SBATCH -t 02:00:00
#SBATCH -o /cfs/klemming/home/a/andnair/.ki-scratch/genes-%j.out
#
# Why 32G and not 80G. 8,930,907 x 616 dense doubles would be 41.0 GiB, which is what
# the row and column counts imply and what a first guess would have asked for. The
# measured uncompressed serialization is 4.58 GiB -- 0.112 of that -- so the object is
# not a dense double matrix of that shape, and 32G leaves room for the load plus a
# dense copy of the subset plus R's own overhead. That measurement cost 0.2 core-hours
# on one core and saved asking for roughly 90 cores' worth of memory in `shared`,
# where memory is allocated at about 0.93 GB per core. If this job is OOM-killed, read
# the reported dim and class from the log before raising the number.
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
