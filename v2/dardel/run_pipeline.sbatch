#!/bin/bash -l
# The v2 pipeline as a Dardel batch job. Submit from v2/dardel/ on Dardel:
#
#   sbatch -A <project> run_pipeline.sbatch --dry-run   # print the reuse plan, seconds
#   sbatch -A <project> run_pipeline.sbatch --smoke      # shake-out, NOT publication output
#   sbatch -A <project> run_pipeline.sbatch              # full run
#   sbatch -A <project> run_pipeline.sbatch --from 06    # this notebook onward
#
# Arguments after the script name are passed straight to run_pipeline.py.
#
# -A is deliberately not set here: the allocation is yours, not the repository's, and
# a wrong one baked into a committed file bills the wrong project.
#
#SBATCH -J ki-v2
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=8G
#SBATCH -t 10:00:00
#SBATCH -o slurm-ki-v2-%j.out
#
# Why -p shared and not -p main: main allocates whole nodes and bills 128 cores
# whatever you use, with a 24-hour ceiling. This job is one sequential chain of
# notebooks, so on main a five-hour run would cost roughly 640 core-hours against
# about 40 here. shared allocates by core and allows 7 days.
#
# Why -c 8 and not more: the bill is cores x wall-clock, and this pipeline does not
# scale linearly, so MORE cores cost MORE core-hours for less than proportional time
# saved. 8 cores for about 5 hours is roughly 40 core-hours; 32 cores would finish
# nearer 90 minutes but bill about 48. Measured: the full run took 4h37m at -c 8.
# Raise it with `sbatch -c 16` when wall time matters more than allocation, and keep
# an eye on `projinfo`. Matching the baseline's 8 does NOT make a Dardel run
# numerically comparable to it -- the CPU decides that, not the core count. See
# CONTRACTS.md section 12.
#
# Memory in shared is proportional to the core count (roughly 0.85 GB per core). The
# baseline ran inside a 1500 MB cap, so --mem=8G is generous. Check the real peak
# afterwards with `seff <jobid>`; the notebooks' own peaks are recorded per notebook
# in the run summary.
#
# Expect `squeue` to show 10 CPUs for this job rather than 8. That is not an error and
# does NOT cost more: 8G needs 10 cores' worth of memory, so SLURM allocates 10 while
# `scontrol show job` reports billing=8. SLURM_CPUS_PER_TASK is then 10, so the worker
# budget picks up the two spare cores for free. Do not "fix" it by cutting --mem.
set -euo pipefail

# Confirm with `ml --terse avail cray-python`. PDC/23.12 is retired (the old trees are
# PDCOLD/*); cray-python lives in the Cray PE core modules and needs no PDC module.
PY_MODULE="${PY_MODULE:-cray-python/3.12.12}"

# Finding the project root. SLURM COPIES the batch script into its spool directory and
# runs it from there, so inside a job BASH_SOURCE[0] is /var/spool/slurm/..., not the
# script in the repository -- which on the first attempt silently resolved the root to
# /var/spool and failed looking for a venv there. SLURM_SUBMIT_DIR is the directory
# the job was submitted FROM, which is the one we want. Order: an explicit
# KI_PROJECT_ROOT wins, then SLURM_SUBMIT_DIR, then BASH_SOURCE for interactive use.
if [ -n "${KI_PROJECT_ROOT:-}" ]; then
    PROJECT_ROOT="$KI_PROJECT_ROOT"
elif [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
    PROJECT_ROOT="$(cd "${SLURM_SUBMIT_DIR}/../.." && pwd)"
else
    PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
fi
VENV_PREFIX="${VENV_PREFIX:-}"

echo "== job ${SLURM_JOB_ID:-interactive} on $(hostname) at $(date -u '+%F %T')Z"
echo "   project root : ${PROJECT_ROOT}"
echo "   cores        : ${SLURM_CPUS_PER_TASK:-unset}"

# Assert it rather than trusting it: a wrong root otherwise surfaces as a confusing
# error about a missing venv, three steps later.
if [ ! -f "${PROJECT_ROOT}/v2/run_pipeline.py" ]; then
    echo "project root looks wrong: no v2/run_pipeline.py under ${PROJECT_ROOT}" >&2
    echo "  submit from <root>/v2/dardel, or pass the root explicitly:" >&2
    echo "  sbatch -A <project> --export=ALL,KI_PROJECT_ROOT=<root> run_pipeline.sbatch" >&2
    exit 1
fi

# Never pipe `ml` into another command: it is a shell function, and a pipe runs it in
# a subshell, so the environment changes are silently discarded.
ml "$PY_MODULE"

if [ -z "$VENV_PREFIX" ]; then
    # Default to where 01_setup_env.sh put it. PDC_STORAGE is not set inside a job,
    # so derive it from the project root: the venv is its sibling.
    VENV_PREFIX="$(dirname "$PROJECT_ROOT")/venv-ki-v2"
fi
[ -f "${VENV_PREFIX}/bin/activate" ] \
    || { echo "no venv at ${VENV_PREFIX}; run 01_setup_env.sh, or set VENV_PREFIX" >&2; exit 1; }
# shellcheck disable=SC1091
source "${VENV_PREFIX}/bin/activate"
echo "   python       : $(command -v python) ($(python -V 2>&1))"

# Where the tree is. utils.py resolves DATA_DIR and RESULTS_DIR from this.
export KI_PROJECT_ROOT="$PROJECT_ROOT"

# Worker budget follows the allocation. Every estimator pins n_jobs=1 and random_state,
# so a single fit does not depend on this -- but NB03's permutation loop stops early on
# the exceedance count, checked at joblib batch boundaries, and those boundaries move
# with the worker count. Measured against the baseline, N_perm for one target differed
# by 557 draws. Note this is 10 here, not 8: SLURM allocates 10 cores to carry 8G.
# Set it to 8 explicitly if you want the permutation budget to track the baseline's.
export KI_N_CORES="${SLURM_CPUS_PER_TASK:-8}"

# KI_BLAS_THREADS is deliberately NOT set, so it stays at the baseline 8 that the
# frozen run used. Raising it reorders floating-point reductions inside BLAS and can
# move the last digits of every number downstream of NB01, which would put the
# manuscript's 271 assertions out of step with the run that produced them.

# Headless, and keep scratch rather than the 25 GB / 100 K-file home for temporaries.
export MPLBACKEND=Agg
export TMPDIR="${PDC_TMP:-/cfs/klemming/scratch/${USER:0:1}/${USER}}/tmp"
export NUMBA_CACHE_DIR="${TMPDIR}/numba"
mkdir -p "$TMPDIR" "$NUMBA_CACHE_DIR"

cd "${PROJECT_ROOT}/v2"

# Run python directly rather than through srun: this is a single-process job with no
# MPI, and SLURM 22.05 onward stopped propagating --cpus-per-task to srun steps, so
# an srun wrapper here is a way to silently get one core instead of 32.
start=$(date +%s)
set +e
python run_pipeline.py "$@"
rc=$?
set -e
echo "== run_pipeline.py exited ${rc} after $(( ($(date +%s) - start) / 60 )) min"

if [ -d "${PROJECT_ROOT}/Results/runs" ]; then
    latest="$(ls -1 "${PROJECT_ROOT}/Results/runs" | tail -1)"
    echo "   run directory : Results/runs/${latest}"
    [ -f "${PROJECT_ROOT}/Results/SMOKE_RUN.txt" ] \
        && echo "   !! SMOKE_RUN.txt present: Results/ is NOT publication output"
    echo "   pull it to the laptop with: bash pull_results.sh ${latest}"
fi
exit $rc
