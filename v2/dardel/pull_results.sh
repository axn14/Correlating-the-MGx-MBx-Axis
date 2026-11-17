#!/usr/bin/env bash
# Pull one Dardel run down to the laptop for evaluation. Runs on the LAPTOP, under WSL.
#
#   bash pull_results.sh                       # the newest run on Dardel
#   bash pull_results.sh 20261108T091200Z_abc1234
#   bash pull_results.sh --list                # what runs exist on Dardel
#   bash pull_results.sh --with-intermediate   # also Results/intermediate (the pickle)
#
# One-way, additive, and into a per-run directory of its own. Three reasons it is
# built this way rather than as a sync of Results/:
#
#   1. run_pipeline.py:283 moves every output directory into Results/_archive/<run_id>
#      at the start of a run and prunes to the newest three by rmtree. A bidirectional
#      sync would propagate those moves and deletions to the laptop.
#   2. The frozen baseline at f14c7f9 is what the manuscript's 271 assertions are
#      pinned to. A Dardel run is a NEW run, not a replacement for it.
#   3. rsync --delete against a tree that is not backed up has no undo.
#
# So: no --delete, and the destination is refused if it resolves inside the project's
# own Results/ tree.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_FILE="${ENV_FILE:-${HERE}/dardel.env}"
[ -f "$ENV_FILE" ] || { echo "missing ${ENV_FILE} -- copy dardel.env.example and fill it in" >&2; exit 1; }
# shellcheck disable=SC1090
source "$ENV_FILE"
# shellcheck source=_paths.sh
source "${HERE}/_paths.sh"

: "${PDC_HOST:?}" "${PDC_USER:?}" "${REMOTE_ROOT:?}" "${LOCAL_ROOT:?}" "${LOCAL_PULL:?}"
normalise_local_paths
configure_ssh
REMOTE="${PDC_USER}@${PDC_HOST}"

RUN_ID=""
WITH_INTERMEDIATE=0
for arg in "$@"; do
    case "$arg" in
        --list)
            echo "runs on ${REMOTE}:${REMOTE_ROOT}/Results/runs"
            ssh "$REMOTE" "ls -1 '${REMOTE_ROOT}/Results/runs' 2>/dev/null \
                | while read -r r; do \
                    printf '  %-34s %s\n' \"\$r\" \"\$(du -sh '${REMOTE_ROOT}/Results/runs/'\$r 2>/dev/null | cut -f1)\"; \
                  done" || echo "  (none)"
            exit 0 ;;
        --with-intermediate) WITH_INTERMEDIATE=1 ;;
        -*) echo "unknown option: $arg" >&2; exit 2 ;;
        *)  [ -z "$RUN_ID" ] || { echo "give at most one run id" >&2; exit 2; }
            RUN_ID="$arg" ;;
    esac
done

if [ -z "$RUN_ID" ]; then
    RUN_ID="$(ssh "$REMOTE" "ls -1 '${REMOTE_ROOT}/Results/runs' 2>/dev/null | tail -1")"
    [ -n "$RUN_ID" ] || { echo "no runs found under ${REMOTE_ROOT}/Results/runs" >&2; exit 1; }
    echo "== newest run on Dardel: ${RUN_ID}"
fi

DEST="${LOCAL_PULL%/}/${RUN_ID}"

# The guard. readlink -m resolves without requiring the path to exist yet.
PROTECTED="$(readlink -m "${LOCAL_ROOT}/Results")"
DEST_ABS="$(readlink -m "$DEST")"
case "$DEST_ABS" in
    "$PROTECTED"|"$PROTECTED"/*)
        echo "REFUSED: ${DEST_ABS} is inside ${PROTECTED}." >&2
        echo "That tree is managed by run_pipeline.py, which moves and prunes it." >&2
        echo "Point LOCAL_PULL somewhere outside the project." >&2
        exit 1 ;;
esac
case "$DEST_ABS" in
    "$(readlink -m "$LOCAL_ROOT")"/*)
        echo "WARNING: pulling inside the project tree (${DEST_ABS})." >&2
        echo "It is outside Results/, so nothing will clobber it, but a stray copy" >&2
        echo "there can be mistaken for pipeline output." >&2 ;;
esac

echo "== pulling into ${DEST}"
# Laid out as a project root, not a loose folder: utils.py resolves DATA_DIR and
# RESULTS_DIR from PROJECT_ROOT, so <dest>/Results plus a Data symlink makes the
# pulled run directly readable by utils.py and by verify_manuscript_numbers.py
# under KI_PROJECT_ROOT, with no copying into the real Results/ tree.
mkdir -p "${DEST}/Results"

# No --delete, no --inplace. -z because these are CSVs and logs; --partial so an
# interrupted pull resumes rather than restarting.
RSYNC_OPTS=(-rltvz --human-readable --partial --info=stats1)

SUBDIRS=(tables figures extensions logs)
[ "$WITH_INTERMEDIATE" -eq 1 ] && SUBDIRS+=(intermediate)

for sub in "${SUBDIRS[@]}"; do
    if ssh "$REMOTE" "[ -d '${REMOTE_ROOT}/Results/${sub}' ]"; then
        echo "-- Results/${sub}"
        mkdir -p "${DEST}/Results/${sub}"
        rsync "${RSYNC_OPTS[@]}" \
            "${REMOTE}:${REMOTE_ROOT}/Results/${sub}/" "${DEST}/Results/${sub}/"
    else
        echo "-- Results/${sub} absent on Dardel, skipped"
    fi
done

echo "-- Results/runs/${RUN_ID} (executed notebooks, HTML, logs, provenance)"
mkdir -p "${DEST}/Results/runs/${RUN_ID}"
rsync "${RSYNC_OPTS[@]}" \
    "${REMOTE}:${REMOTE_ROOT}/Results/runs/${RUN_ID}/" \
    "${DEST}/Results/runs/${RUN_ID}/"

# The smoke marker is the one file that decides whether any of this is citable.
SMOKE=0
if ssh "$REMOTE" "[ -f '${REMOTE_ROOT}/Results/SMOKE_RUN.txt' ]"; then
    rsync "${RSYNC_OPTS[@]}" \
        "${REMOTE}:${REMOTE_ROOT}/Results/SMOKE_RUN.txt" "${DEST}/Results/"
    SMOKE=1
fi

# Data is read-only here and 53.7 GB, so it is linked, never copied. utils.py only
# asserts that DATA_DIR is a directory, which a symlink to the real one satisfies.
if [ ! -e "${DEST}/Data" ]; then
    ln -s "$(readlink -m "${LOCAL_ROOT}/Data")" "${DEST}/Data"
    echo "-- linked Data -> ${LOCAL_ROOT}/Data"
fi

echo
echo "== pulled $(du -sh --exclude=Data "$DEST" | cut -f1) to ${DEST}"
find "${DEST}/Results" -maxdepth 1 -mindepth 1 -type d -printf '   Results/%f\n' | sort
if [ "$SMOKE" -eq 1 ]; then
    echo
    echo "!! SMOKE_RUN.txt was present: this run is a shake-out, NOT publication output."
fi

cat <<EOF

Rendered notebooks: ${DEST}/Results/runs/${RUN_ID}/
Open the .html files in a browser -- no Python environment needed to read them.
Tables are CSV under ${DEST}/Results/tables/.

To see how this run compares with the frozen baseline, target by target:

  python "${LOCAL_ROOT}/v2/dardel/compare_to_baseline.py" "${DEST}"

Expect differences. A run on other hardware is a new run: OpenBLAS picks different
kernels per microarchitecture, and nested CV turns last-digit differences into
percent-level ones through hyperparameter argmax. CONTRACTS.md section 12 has the
measured size of it. That is information about robustness, not a problem to fix by
copying files over Results/.

Do NOT reach for writeup/verify_manuscript_numbers.py here. It reads
D:\\KI Research\\Results\\tables by a hard-coded path, so KI_PROJECT_ROOT does not
redirect it: pointed at a pulled tree it checks the laptop's own tables, silently
skips whatever it resolves through utils, and still prints "no mismatches".
EOF
