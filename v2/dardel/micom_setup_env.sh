#!/bin/bash -l
# Build the MICOM environment on Dardel. Run ONCE, on a login node.
#
#   PDC_STORAGE=/cfs/klemming/projects/supr/<proj>/<user> bash micom_setup_env.sh
#
# A SEPARATE venv from 01_setup_env.sh's venv-ki-v2: micom and its solvers are pinned
# outside the pipeline lock (CONTRACTS section 8), and installing them into the pipeline
# venv would change the environment every certified pipeline run is pinned to. This one
# takes the pipeline lock first (micom_batch.py imports utils, which needs it) and then
# requirements-micom.txt on top -- the versions the laptop solver check ran with.
set -euo pipefail
: "${PDC_STORAGE:?set PDC_STORAGE to your Klemming project directory}"
VENV_PREFIX="${VENV_PREFIX:-${PDC_STORAGE}/venv-micom}"
REPO="${REPO:-${PDC_STORAGE}/ki}"

ml cray-python/3.12.12
[ -d "$VENV_PREFIX" ] && { echo "exists: $VENV_PREFIX (delete it to rebuild)"; exit 1; }
python3 -m venv "$VENV_PREFIX"
# shellcheck disable=SC1091
source "${VENV_PREFIX}/bin/activate"
pip install --upgrade pip==25.3
pip install --no-cache-dir -r "${REPO}/v2/requirements.lock"
pip install --no-cache-dir -r "${REPO}/v2/requirements-micom.txt"
python - <<'EOF'
import micom, optlang, highspy, osqp, cobra
print("micom", micom.__version__, "optlang", optlang.__version__,
      "highspy", highspy.Highs().version(), "osqp", osqp.__version__, "cobra", cobra.__version__)
EOF
du -sh "$VENV_PREFIX"
