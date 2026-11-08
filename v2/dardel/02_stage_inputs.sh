#!/usr/bin/env bash
# Push code and inputs from the laptop to Dardel. Runs on the LAPTOP, under WSL
# (rsync, ssh and git live there; Windows has none of them on PATH for this tree).
#
#   bash 02_stage_inputs.sh                 # code at HEAD + the cohort tables
#   bash 02_stage_inputs.sh --with-micom    # also the 53.5 GB of prebuilt MICOM models
#   bash 02_stage_inputs.sh --rsync-code    # overlay the working tree (see below)
#
# The pipeline needs far less input than the tree suggests. run_pipeline.py hashes
# Data/*.tsv (flat) plus Data/reference/agora2_species_ec_index.json as its shared
# inputs: 164 MB in total. Data/micom is 53.5 GB across 2,756 files and is read only
# by NB10, which run_pipeline.py does not run, so --with-micom is for the MICOM work
# later and not for a first port.
#
# Code travels as a git bundle rather than an rsync of v2/, because run_pipeline.py
# derives the run id from .git: a real repository on the far side gives a clean
# commit and no "-dirty" tag, where rsyncing a subset of the tree over staged git
# metadata would mark every run dirty and lose the provenance. --rsync-code overlays
# the working tree afterwards for iterating on uncommitted edits; runs made that way
# are correctly tagged "-dirty" and must not be used for publication output.
#
# Results/ is NEVER pushed. run_pipeline.py:283 moves every output directory into
# Results/_archive/<run_id> at the start of a run and prunes to the newest three by
# rmtree, so uploading the frozen baseline would relocate it and then delete it.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_FILE="${ENV_FILE:-${HERE}/dardel.env}"
[ -f "$ENV_FILE" ] || { echo "missing ${ENV_FILE} -- copy dardel.env.example and fill it in" >&2; exit 1; }
# shellcheck disable=SC1090
source "$ENV_FILE"
# shellcheck source=_paths.sh
source "${HERE}/_paths.sh"

: "${PDC_HOST:?}" "${PDC_USER:?}" "${REMOTE_ROOT:?}" "${LOCAL_ROOT:?}"
normalise_local_paths
configure_ssh

WITH_MICOM=0
RSYNC_CODE=0
for arg in "$@"; do
    case "$arg" in
        --with-micom) WITH_MICOM=1 ;;
        --rsync-code) RSYNC_CODE=1 ;;
        *) echo "unknown argument: $arg" >&2; exit 2 ;;
    esac
done

REMOTE="${PDC_USER}@${PDC_HOST}"
BRANCH="${BRANCH:-main}"
echo "== local root: ${LOCAL_ROOT}"

# No --delete anywhere in this script, deliberately: a mistyped path with --delete
# against a project directory that is not backed up has no undo.
RSYNC_OPTS=(-rltvz --human-readable --partial --info=stats1)

echo "== local state"
git -C "$LOCAL_ROOT" rev-parse --short HEAD
if [ -n "$(git -C "$LOCAL_ROOT" status --porcelain)" ]; then
    echo "   WARNING: working tree is dirty. The bundle carries HEAD only." >&2
    if [ "$RSYNC_CODE" -eq 0 ]; then
        echo "   Uncommitted edits will NOT reach Dardel. Commit them, or pass --rsync-code." >&2
    fi
else
    echo "   clean"
fi

echo "== creating the remote tree"
ssh "$REMOTE" "mkdir -p '${REMOTE_ROOT}/Data/reference' '${REMOTE_ROOT}/.stage'"

echo "== bundling the repository"
BUNDLE="$(mktemp -t ki-bundle-XXXXXX.bundle)"
trap 'rm -f "$BUNDLE"' EXIT
git -C "$LOCAL_ROOT" bundle create "$BUNDLE" --branches --tags
rsync "${RSYNC_OPTS[@]}" "$BUNDLE" "${REMOTE}:${REMOTE_ROOT}/.stage/ki.bundle"

echo "== unpacking it into a working repository on Dardel"
# git init is safe in a directory that already holds Data/. The checkout refuses
# rather than clobbering if tracked files were modified on the far side, which is
# the behaviour we want for a machine-managed tree.
ssh "$REMOTE" "set -e
    cd '${REMOTE_ROOT}'
    [ -d .git ] || git init -q
    git fetch -q '${REMOTE_ROOT}/.stage/ki.bundle' '+refs/heads/*:refs/remotes/bundle/*'
    git checkout -q -B '${BRANCH}' 'refs/remotes/bundle/${BRANCH}'
    git rev-parse --short HEAD"

if [ "$RSYNC_CODE" -eq 1 ]; then
    echo "== overlaying the working tree (runs will be tagged -dirty)"
    rsync "${RSYNC_OPTS[@]}" \
        --exclude '__pycache__/' \
        --exclude '.ipynb_checkpoints/' \
        --exclude '*.pyc' \
        --exclude 'dardel/dardel.env' \
        "${LOCAL_ROOT}/v2/" "${REMOTE}:${REMOTE_ROOT}/v2/"
fi

echo "== cohort tables (Data/*.tsv)"
# Globbed, not -r: a recursive sync of Data/ would pull in the 53.5 GB of micom models.
rsync -ltvz --human-readable --partial --info=stats1 \
    "${LOCAL_ROOT}"/Data/*.tsv "${REMOTE}:${REMOTE_ROOT}/Data/"

echo "== AGORA2 EC index (Data/reference/)"
rsync "${RSYNC_OPTS[@]}" \
    "${LOCAL_ROOT}/Data/reference/" "${REMOTE}:${REMOTE_ROOT}/Data/reference/"

if [ "$WITH_MICOM" -eq 1 ]; then
    echo "== MICOM models (Data/micom/) -- 53.5 GB, 2,756 files"
    echo "   resumable: re-run this script if the transfer is interrupted."
    rsync "${RSYNC_OPTS[@]}" --append-verify \
        "${LOCAL_ROOT}/Data/micom/" "${REMOTE}:${REMOTE_ROOT}/Data/micom/"
fi

echo "== verifying the inputs run_pipeline.py will hash"
ssh "$REMOTE" "cd '${REMOTE_ROOT}'
    echo \"  commit          : \$(git rev-parse --short HEAD) dirty=\$([ -n \"\$(git status --porcelain)\" ] && echo yes || echo no)\"
    echo \"  Data/*.tsv      : \$(find Data -maxdepth 1 -name '*.tsv' | wc -l) files\"
    echo \"  Data/reference  : \$(ls -1 Data/reference 2>/dev/null | wc -l) files\"
    echo \"  v2/*.ipynb      : \$(ls -1 v2/*.ipynb 2>/dev/null | wc -l) notebooks\"
    echo \"  Results/        : \$([ -d Results ] && echo present || echo 'absent (correct before the first run)')\""

cat <<EOF

== staged
Next, on Dardel:
  PDC_STORAGE=<your project dir> bash ${REMOTE_ROOT}/v2/dardel/01_setup_env.sh
then
  cd ${REMOTE_ROOT}/v2/dardel && sbatch -A <your-project> run_pipeline.sbatch --dry-run
EOF
