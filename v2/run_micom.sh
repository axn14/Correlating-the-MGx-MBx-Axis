#!/usr/bin/env bash
# NB10 — MICOM community flux. DECOUPLED from run_pipeline.sh (hours of compute).
# Run AFTER the main pipeline, in its own tmux session:
#   tmux new -s micom 'bash ~/ki/v2/run_micom.sh batch; exec bash'
#
# Modes:
#   bash run_micom.sh smoke   # build models + tradeoff sweep + ONE-sample smoke gate, then stop
#   bash run_micom.sh batch   # the full per-sample batch (hours) + the three output tables
#   bash run_micom.sh         # defaults to smoke (safe)
#
# After a successful batch, re-run NB08 then NB09 so the mechanistic column fills.
set -uo pipefail
cd "$(dirname "$0")"

MODE="${1:-smoke}"
VENV="${VENV:-$HOME/ki/.venv}"
PY="$VENV/bin/python"
LOG_DIR="../Results/logs"   # v2/ sits directly under the project root
SUMMARY="$LOG_DIR/summary_micom_v2.txt"
mkdir -p "$LOG_DIR"

# user ki is capped at 2 cores / 1500 MB. MICOM's own threads come from MICOM_THREADS.
export MPLBACKEND=Agg
export OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 NUMEXPR_MAX_THREADS=2 LOKY_MAX_CPU_COUNT=2
export MICOM_THREADS="${MICOM_THREADS:-2}"
export MICOM_CHUNK="${MICOM_CHUNK:-25}"

if [ "$MODE" = batch ]; then
    export MICOM_RUN_BATCH=1
    echo "MODE=batch : full per-sample solve (this takes hours)."
else
    unset MICOM_RUN_BATCH 2>/dev/null || true
    echo "MODE=smoke : stops after the one-sample smoke gate. Use 'batch' for the full run."
fi

echo "micom run ($MODE) started $(date '+%F %T')" > "$SUMMARY"

# 1. (Re)generate the notebook from its edit script — single source of truth.
echo "[$(date '+%F %T')] building 10_micom_flux.ipynb from edits_nb10.py"
if ! "$PY" edits_nb10.py > "$LOG_DIR/edits_nb10.log" 2>&1; then
    echo "FAILED to build notebook — see $LOG_DIR/edits_nb10.log" | tee -a "$SUMMARY"; exit 1
fi

# 2. Execute it (memory-tracked, no cell timeout). The notebook itself stops after
#    the smoke gate unless MICOM_RUN_BATCH=1, so 'smoke' mode is safe to --execute.
echo "[$(date '+%F %T')] executing NB10 ($MODE)"
start=$(date +%s)
"$PY" pipeline_tools.py memrun "$LOG_DIR/10_micom_flux.peak" \
    "$VENV/bin/jupyter" nbconvert --to notebook --execute --inplace \
    --ExecutePreprocessor.timeout=-1 --ExecutePreprocessor.kernel_name=python3 \
    10_micom_flux.ipynb > "$LOG_DIR/10_micom_flux.log" 2>&1
rc=$?
secs=$(( $(date +%s) - start ))
peak=$(cat "$LOG_DIR/10_micom_flux.peak" 2>/dev/null || echo "?")
if [ $rc -eq 0 ]; then status=OK; else status=FAIL; fi
printf '%-28s %-4s %7ss  peak %5s MB  (%s)\n' "10_micom_flux.ipynb" "$status" "$secs" "$peak" "$MODE" | tee -a "$SUMMARY"
if [ $rc -ne 0 ]; then
    echo "NB10 failed — see $LOG_DIR/10_micom_flux.log" | tee -a "$SUMMARY"; cat "$SUMMARY"; exit 1
fi

echo "micom run ($MODE) finished $(date '+%F %T')" | tee -a "$SUMMARY"
if [ "$MODE" = batch ]; then
    {
        echo
        echo "MICOM batch done. Now re-run NB08 + NB09 to fill the mechanistic column:"
        echo "  for nb in 08_pathway_analysis 09_pathway_network; do"
        echo "    \$VENV/bin/jupyter nbconvert --to notebook --execute --inplace \$nb.ipynb; done"
    } | tee -a "$SUMMARY"
fi
cat "$SUMMARY"
