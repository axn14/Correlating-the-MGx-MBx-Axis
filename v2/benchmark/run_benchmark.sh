#!/usr/bin/env bash
# Benchmark the v2 framework against MIMOSA2 + MelonnPan on the Yachida cohort.
# Decoupled from the main pipeline. Run AFTER NB01-09 have produced
# final_evidence_synthesis_v2.csv.  (NB10/MICOM is not required for the benchmark.)
#
#   bash run_benchmark.sh
#
# Prereqs (one-time):
#   - R with: data.table, devtools, melonnpan (biobakery/melonnpan), and MIMOSA2 (../../MIMOSA2)
#   - MIMOSA2 reference data (AGORA/KEGG). Set MIMOSA2_KEGG_PREFIX / MIMOSA2_DATA_PREFIX,
#     or let MIMOSA2::download_reference_data() fetch them.
set -uo pipefail
cd "$(dirname "$0")"

VENV="${VENV:-$HOME/ki/.venv}"
PY="$VENV/bin/python"
RSCRIPT="${RSCRIPT:-Rscript}"
export OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MPLBACKEND=Agg

echo "[1/3] exporting matrices (python)"
"$PY" export_for_benchmark.py || { echo "export failed"; exit 1; }

echo "[2/3] running MelonnPan + MIMOSA2 (R)"
if ! "$RSCRIPT" run_mimosa2_melonnpan.R; then
    echo "R step reported problems — inspect benchmark/outputs/. Continuing to compare what exists."
fi

echo "[3/3] comparing against the framework (python)"
"$PY" compare_frameworks.py || { echo "compare failed"; exit 1; }

echo "Done. See benchmark/outputs/: benchmark_comparison_v2.csv, benchmark_overlap_v2.csv, benchmark_concordance_v2.png"
