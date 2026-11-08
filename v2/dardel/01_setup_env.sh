#!/bin/bash -l
# Build the pinned pipeline environment on Dardel. Run ONCE, on a login node.
#
#   PDC_STORAGE=/cfs/klemming/projects/supr/<proj>/<user> bash 01_setup_env.sh
#
# Why cray-python + venv rather than miniconda3, having checked both on the machine
# (2026-11-07): `ml cray-python/3.12.12` gives Python 3.12.12 with venv and pip 25.3,
# and miniconda3/26.3.2-2-cpeGNU-26.03 exists under PDC/26.03. venv wins here on
# three counts.
#
#   1. requirements.lock is a `pip freeze` taken on Python 3.12.10, so every wheel
#      carries a cp312 tag and installs unchanged on 3.12.12. conda would supply only
#      the interpreter and pip would do the same work anyway.
#   2. BLAS comparability. The frozen baseline used pip wheels, which bundle
#      OpenBLAS; pip wheels on linux-x86_64 bundle OpenBLAS too, so numpy and scipy
#      link the same library family. conda-forge would substitute its own build, and
#      a different BLAS is exactly the thing that moves the last digits.
#   3. Fewer moving parts: no .condarc relocation, no channel resolution.
#
# The environment does NOT go in $HOME: home is 25 GB and capped at 100 K files, and
# this venv is tens of thousands on its own. It goes in project storage, which has
# the room (checked: ~340 K files of quota free).
set -euo pipefail

: "${PDC_STORAGE:?set PDC_STORAGE to your Klemming project directory (check with projinfo)}"
VENV_PREFIX="${VENV_PREFIX:-${PDC_STORAGE}/venv-ki-v2}"

# Confirm with `ml --terse avail cray-python`. Note PDC/23.12 is retired: the live
# stack is PDC/26.03 and the old trees are renamed PDCOLD/*. cray-python sits in the
# Cray PE core modules and needs no PDC module loaded.
PY_MODULE="${PY_MODULE:-cray-python/3.12.12}"

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOCK="${HERE}/../requirements.lock"
[ -f "$LOCK" ] || { echo "missing $LOCK -- stage the code first (02_stage_inputs.sh)" >&2; exit 1; }

echo "== target"
echo "   storage : ${PDC_STORAGE}"
echo "   venv    : ${VENV_PREFIX}"
echo "   lock    : ${LOCK}"
[ -d "$PDC_STORAGE" ] || { echo "PDC_STORAGE does not exist or is not readable: ${PDC_STORAGE}" >&2; exit 1; }

echo "== loading ${PY_MODULE}"
# Never pipe `ml` into head/grep: it is a shell function, and a pipe runs it in a
# subshell so the environment changes are silently discarded.
ml "$PY_MODULE"
echo "   python  : $(command -v python3) ($(python3 -V 2>&1))"

case "$(python3 -V 2>&1)" in
    *" 3.12."*) ;;
    *) echo "WARNING: requirements.lock was frozen on 3.12; this is $(python3 -V 2>&1)." >&2
       echo "         Wheel tags will not match and pip may build from source." >&2 ;;
esac

if [ -d "$VENV_PREFIX" ]; then
    echo "== venv already exists; installing into it"
else
    echo "== creating the venv"
    python3 -m venv "$VENV_PREFIX"
fi
# shellcheck disable=SC1091
source "${VENV_PREFIX}/bin/activate"

echo "== python in use"
python -V
python -c "import sys; print(sys.executable)"

echo "== installing the pinned environment"
# No --system-site-packages on the venv, deliberately: cray-python ships numpy 2.3.5
# and we need the pinned 2.3.1, so the venv must not see the module's packages.
python -m pip install --no-input --upgrade pip
python -m pip install --no-input -r "$LOCK"

echo "== verifying the imports the pipeline needs"
python - <<'PY'
import importlib
mods = ["numpy", "pandas", "scipy", "sklearn", "xgboost", "lightgbm", "shap",
        "optuna", "mofapy2", "pingouin", "statsmodels", "networkx", "community",
        "matplotlib", "seaborn", "psutil", "nbformat", "nbconvert", "threadpoolctl"]
bad = []
for m in mods:
    try:
        importlib.import_module(m)
    except Exception as exc:
        bad.append(f"{m}: {exc!r}")
print(f"{len(mods) - len(bad)}/{len(mods)} imports ok")
if bad:
    raise SystemExit("FAILED imports:\n  " + "\n  ".join(bad))
PY

echo "== pinned versions actually installed (against the baseline)"
python - <<'PY'
import importlib.metadata as md
want = {"numpy": "2.3.1", "scipy": "1.16.0", "pandas": "2.3.0",
        "scikit-learn": "1.7.0", "xgboost": "3.0.2", "shap": "0.50.0"}
for name, exp in want.items():
    got = md.version(name)
    print(f"  {name:<14} {got:<10} {'ok' if got == exp else 'EXPECTED ' + exp}")
PY

echo "== BLAS the installed numpy actually links"
python - <<'PY'
from threadpoolctl import threadpool_info
for d in threadpool_info():
    print(f"  {d.get('internal_api')} {d.get('version')} {d.get('filepath')}")
PY

echo "== file count in the venv (project quota is 1,000,000 shared)"
find "$VIRTUAL_ENV" -type f 2>/dev/null | wc -l

cat <<EOF

== done
Activate it in a job with:
  ml ${PY_MODULE}
  source ${VENV_PREFIX}/bin/activate

Next:
  cd $(cd "${HERE}" && pwd) && sbatch -A <your-project> run_pipeline.sbatch --dry-run
EOF
