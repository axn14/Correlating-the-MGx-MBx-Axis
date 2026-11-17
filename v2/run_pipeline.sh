#!/usr/bin/env bash
# Execute the v2 pipeline headless on the server, in dependency order (see CONTRACTS.md).
#   tmux new -s ki 'bash ~/ki/v2/run_pipeline.sh; exec bash'
# 07_mediation is markdown-only and is skipped. Logs + summary go to <project root>/Results/logs/.
# NB03_N_PERM (default 100) sets NB03's permutation count.
set -uo pipefail

cd "$(dirname "$0")"
VENV="${VENV:-$HOME/ki/.venv}"
PY="$VENV/bin/python"
LOG_DIR="../Results/logs"   # v2/ sits directly under the project root
SUMMARY="$LOG_DIR/summary_v2.txt"
mkdir -p "$LOG_DIR"

# user ki is capped at 2 cores / 1500 MB: never let BLAS/OpenMP/joblib spawn more workers.
export MPLBACKEND=Agg
export OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 NUMEXPR_MAX_THREADS=2 LOKY_MAX_CPU_COUNT=2

echo "v2.1 pipeline run started $(date '+%F %T')" > "$SUMMARY"

# ---------------------------------------------------------------------------
# Pre-flight checks. These gate the whole run: a broken group definition, a
# leaked covariate or an inverted effect sign must not be able to produce
# results that look finished. checks.py is fast (synthetic data + metadata).
# ---------------------------------------------------------------------------
echo "[$(date '+%F %T')] pre-flight checks"
if ! "$PY" checks.py > "$LOG_DIR/checks.log" 2>&1; then
    {
        echo "PRE-FLIGHT CHECKS FAILED -- pipeline not started."
        echo "See $LOG_DIR/checks.log. Failing checks:"
        grep -E '^(FAIL|  - )' "$LOG_DIR/checks.log" || true
    } | tee -a "$SUMMARY"
    exit 1
fi
printf '%-42s %-4s  %s\n' "checks.py" "OK" \
    "$(grep -c '^ok' "$LOG_DIR/checks.log") checks passed" | tee -a "$SUMMARY"

# run <notebook> <fatal|nonfatal>
run() {
    local nb="$1" mode="$2" base="${1%.ipynb}" start rc secs peak
    echo "[$(date '+%F %T')] start $nb"
    start=$(date +%s)
    "$PY" pipeline_tools.py memrun "$LOG_DIR/$base.peak" \
        "$VENV/bin/jupyter" nbconvert --to notebook --execute --inplace \
        --ExecutePreprocessor.timeout=-1 --ExecutePreprocessor.kernel_name=python3 \
        "$nb" > "$LOG_DIR/$base.log" 2>&1
    rc=$?
    secs=$(( $(date +%s) - start ))
    peak=$(cat "$LOG_DIR/$base.peak" 2>/dev/null || echo "?")
    if [ $rc -eq 0 ]; then status=OK; else status=FAIL; fi
    printf '%-42s %-4s %6ss  peak %5s MB  (%s)\n' "$nb" "$status" "$secs" "$peak" "$mode" | tee -a "$SUMMARY"
    "$PY" pipeline_tools.py scan "$nb" | sed 's/^/    /' >> "$SUMMARY"
    echo "[$(date '+%F %T')] done  $nb ($status)"
    if [ $rc -ne 0 ] && [ "$mode" = fatal ]; then
        echo "FATAL: $nb failed -- see $LOG_DIR/$base.log; stopping." | tee -a "$SUMMARY"
        cat "$SUMMARY"
        exit 1
    fi
}

run 01_preprocessing.ipynb                  fatal
run 01b_yachida_significant_features.ipynb  nonfatal   # nothing downstream reads its output
run 02_analysis.ipynb                       fatal
run 03_ml_benchmarking.ipynb                fatal

# Sequential on purpose (2-core / 1500 MB cap).
run 04_crc_stage_tumour_site.ipynb          nonfatal   # NB08 reads its stage-trend table,
                                                       # which NB08 treats as optional
run 05_two_stage_model.ipynb                nonfatal   # leaf analysis
run 06_multiomics_factors.ipynb             fatal      # NB08/NB09 need its tables

run 08_pathway_analysis.ipynb               fatal      # writes pair_evidence_nb08_v2.csv
run 09_pathway_network.ipynb                fatal      # reads it

echo "pipeline finished $(date '+%F %T')" | tee -a "$SUMMARY"
{
    echo
    echo "NB10 (MICOM) is NOT part of this run: bash run_micom.sh in its own tmux session."
    echo "After NB10 finishes, re-run 08 and 09 to fill the mechanistic evidence column:"
    echo "  for nb in 08_pathway_analysis 09_pathway_network; do"
    echo "    \$VENV/bin/jupyter nbconvert --to notebook --execute --inplace \$nb.ipynb; done"
} | tee -a "$SUMMARY"
cat "$SUMMARY"
