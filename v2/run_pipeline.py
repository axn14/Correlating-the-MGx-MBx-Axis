#!/usr/bin/env python3
"""Execute the v2 pipeline headless, in dependency order (see CONTRACTS.md).

    python run_pipeline.py                  # full run, reusing still-valid results
    python run_pipeline.py --fresh          # recompute the selected notebooks
    python run_pipeline.py --only 02 03     # just these, by stem prefix
    python run_pipeline.py --from 06        # this notebook onward
    python run_pipeline.py --dry-run        # print the plan, with each reuse decision
    python run_pipeline.py --smoke          # fast shake-out pass, NOT publication output

Execution is OUT OF PLACE. Each notebook runs with v2/ as its working directory, but
the executed copy (with outputs) goes to Results/runs/<run_id>/, alongside an HTML
export, the per-notebook log, the run summary and provenance.json. Tracked notebooks
under v2/ are never modified, so git only ever sees source.

Reuse rule. A previous result is reused only if ALL hold:
  - it succeeded;
  - the notebook's CODE cells are byte-identical (outputs and markdown are not part
    of the key);
  - the shared inputs are identical: utils.py, viz.py, requirements.txt, every
    Data/*.tsv and the prebuilt Data/reference index, by sha256;
  - the result-affecting profile (NB03_* knobs) is identical;
  - it finished after every upstream notebook's latest run.
The first version of this runner keyed reuse on the notebook file's size+mtime only,
so editing utils.py -- where the masks, CLR and DA live -- left every notebook
"reusable" and a resume would have served results computed by the old code.

Freshness. When NB01 is about to run, Results/{tables,intermediate,figures} are moved
to Results/_archive/<run_id>/ first, so a file whose writer was removed can never
satisfy a downstream reader. Any other notebook that runs has the outputs its contract
declares moved there first, for the same reason. The newest three archives are kept.

Output contract (output_contract.py). After a notebook exits 0, each output it
declares must exist, be written by this execution, and (CSV) carry a header, the
required columns and enough rows, or (pickle) the required keys; numbers CONTRACTS.md
quotes must match what it produced. A violation turns OK into FAIL. The same check
runs before a result is reused, so a deleted or hand-edited output forces a re-run.
Files a notebook writes that its contract does not declare are listed, not failed.

Thread policy. Oversubscription is the real risk on 8 physical cores: BLAS threads
multiply with joblib workers, so 8 x 8 = 64 threads thrash. NB03 parallelises at the
joblib level and gets BLAS=1 with the cores spent on workers via KI_N_JOBS; every
other notebook gets the cores through BLAS. PYTHONHASHSEED and SOURCE_DATE_EPOCH are
pinned so set iteration order and PDF timestamps cannot vary between runs.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import output_contract as oc

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RESULTS = ROOT / "Results"
LOG_DIR = RESULTS / "logs"
RUNS_DIR = RESULTS / "runs"
ARCHIVE_DIR = RESULTS / "_archive"
LATEST_SUMMARY = LOG_DIR / "summary_v2.txt"
MANIFEST = LOG_DIR / "pipeline_manifest.json"
ARCHIVES_KEPT = 3

# Inputs every notebook depends on. A change to any of them invalidates every result.
SHARED_CODE = ("utils.py", "viz.py", "requirements.txt", "requirements.lock")
DATA_GLOB = "*.tsv"                      # flat cohort tables in Data/ (not Data/micom)
# Prebuilt reference inputs (tools/build_agora2_index.py writes the AGORA2 EC index).
DATA_REFERENCE = ("reference/agora2_species_ec_index.json",)

# 8 physical cores / 16 logical on the Ryzen 7 6800H. Compute-bound tree building
# does not benefit much from SMT, so the budget is physical cores.
N_CORES = 8

_DETERMINISM = {
    "PYTHONHASHSEED": "0",
    "SOURCE_DATE_EPOCH": "1767225600",   # 2026-01-01T00:00:00Z
}
# Linear-algebra notebooks: give the cores to BLAS.
ENV_BLAS = {
    **_DETERMINISM,
    "OMP_NUM_THREADS": str(N_CORES),
    "MKL_NUM_THREADS": str(N_CORES),
    "OPENBLAS_NUM_THREADS": str(N_CORES),
    "NUMEXPR_MAX_THREADS": str(N_CORES),
    "LOKY_MAX_CPU_COUNT": str(N_CORES),
}
# joblib-parallel notebooks: one thread per fit, cores spent on workers.
ENV_JOBLIB = {
    **_DETERMINISM,
    "OMP_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "NUMEXPR_MAX_THREADS": "1",
    "LOKY_MAX_CPU_COUNT": str(N_CORES),
    "KI_N_JOBS": str(N_CORES),
}

# (notebook, fatal?, thread env, note)
# NB04 is absent by design (discarded). NB07 is markdown-only. NB10 (MICOM) is not
# part of this sequence: it needs micom/cobra, which are not in the pinned
# environment (CONTRACTS.md section 8).
NOTEBOOKS: list[tuple[str, bool, dict[str, str], str]] = [
    ("01_preprocessing.ipynb", True, ENV_BLAS, "builds the reduced pickle"),
    ("01b_yachida_significant_features.ipynb", False, ENV_BLAS, "nothing downstream reads it"),
    ("02_analysis.ipynb", True, ENV_BLAS, "differential abundance, correlations"),
    ("03_ml_benchmarking.ipynb", True, ENV_JOBLIB, "nested CV + permutation"),
    # NB05 does not read KI_N_JOBS (its LightGBM fits pin n_jobs themselves), so the
    # joblib profile would only starve its sklearn/BLAS work of threads.
    ("05_two_stage_model.ipynb", False, ENV_BLAS, "two-stage model"),
    ("06_multiomics_factors.ipynb", True, ENV_BLAS, "NB08/NB09 need its tables"),
    ("08_pathway_analysis.ipynb", True, ENV_BLAS, "evidence assembly"),
    ("09_pathway_network.ipynb", True, ENV_BLAS, "network + evidence synthesis"),
]
ORDER = [n for n, *_ in NOTEBOOKS]

# Environment knobs that change RESULTS, not merely speed. KI_N_JOBS is excluded:
# with fixed seeds and single-threaded fits inside the parallel loop it changes wall
# time, not output (verified byte-identical for XGBoost under 1 vs 8 workers).
PROFILE_KEYS = ("NB03_N_PERM", "NB03_N_ITER", "NB03_MAX_TARGETS")
# First notebook whose outputs depend on the profile; everything upstream is immune.
PROFILE_FIRST = "03_ml_benchmarking.ipynb"

# Shake-out profile: exercises every code path end to end in minutes. NOT
# publication output -- N_PERM=19 floors the permutation p at 0.05, and only three
# targets are benchmarked.
SMOKE_PROFILE = {"NB03_N_PERM": "19", "NB03_N_ITER": "2", "NB03_MAX_TARGETS": "3"}
SMOKE_MARKER = RESULTS / "SMOKE_RUN.txt"


# ---------------------------------------------------------------------------
# Fingerprints
# ---------------------------------------------------------------------------

def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def code_digest(nb_path: Path) -> str:
    """sha256 of the notebook's code-cell sources only. Outputs, execution counts and
    markdown are excluded, so executing a notebook never changes its own key."""
    nb = json.loads(nb_path.read_text(encoding="utf-8"))
    sources = ["".join(c.get("source", "")) for c in nb.get("cells", [])
               if c.get("cell_type") == "code"]
    return hashlib.sha256("\x1e".join(sources).encode("utf-8")).hexdigest()


def shared_inputs() -> dict[str, str]:
    """sha256 of every input shared by all notebooks."""
    out = {f"v2/{name}": _sha256_file(HERE / name) for name in SHARED_CODE}
    for p in sorted((ROOT / "Data").glob(DATA_GLOB)):
        out[f"Data/{p.name}"] = _sha256_file(p)
    for rel in DATA_REFERENCE:
        p = ROOT / "Data" / rel
        out[f"Data/{rel}"] = _sha256_file(p) if p.exists() else "absent"
    return out


def inputs_digest(inputs: dict[str, str]) -> str:
    blob = "\n".join(f"{k}={v}" for k, v in sorted(inputs.items()))
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()


def current_profile() -> dict:
    return {k: os.environ.get(k, "") for k in PROFILE_KEYS}


# ---------------------------------------------------------------------------
# Provenance
# ---------------------------------------------------------------------------

def git_info() -> dict:
    """HEAD commit read from .git directly (git is not on the Windows PATH here), and
    the dirty state via `git` if available, else `wsl git`, else 'unknown'."""
    gd = ROOT / ".git"
    info: dict = {"commit": None, "branch": None, "dirty": "unknown"}
    try:
        head = (gd / "HEAD").read_text(encoding="utf-8").strip()
        if head.startswith("ref: "):
            ref = head[5:]
            info["branch"] = ref.rsplit("/", 1)[-1]
            loose = gd / ref
            if loose.exists():
                info["commit"] = loose.read_text(encoding="utf-8").strip()
            elif (gd / "packed-refs").exists():
                for ln in (gd / "packed-refs").read_text(encoding="utf-8").splitlines():
                    if ln.endswith(" " + ref):
                        info["commit"] = ln.split()[0]
        else:
            info["commit"] = head
    except OSError:
        pass

    wsl_root = f"/mnt/{str(ROOT)[0].lower()}{str(ROOT)[2:]}".replace("\\", "/")
    for cmd in (["git", "-C", str(ROOT), "status", "--porcelain"],
                ["wsl", "-e", "git", "-C", wsl_root, "status", "--porcelain"]):
        try:
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        except (OSError, subprocess.SubprocessError):
            continue
        if r.returncode == 0:
            changed = [ln for ln in r.stdout.splitlines() if ln.strip()]
            info["dirty"] = bool(changed)
            info["dirty_files"] = changed[:50]
            break
    return info


def write_provenance(run_dir: Path, run_id: str, mode: str, profile: dict,
                     inputs: dict, git: dict) -> None:
    from importlib import metadata
    pkgs = {}
    for dist in metadata.distributions():
        name = dist.metadata.get("Name")
        if name:
            pkgs[name.lower()] = dist.version
    blas = None
    try:
        from threadpoolctl import threadpool_info
        blas = [{k: d.get(k) for k in ("internal_api", "version", "filepath")}
                for d in threadpool_info()]
    except Exception as exc:          # provenance must never stop a run
        blas = f"unavailable: {exc!r}"
    prov = {
        "run_id": run_id,
        "mode": mode,
        "profile": profile,
        "started_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "git": git,
        "python": sys.version,
        "platform": platform.platform(),
        "machine": platform.machine(),
        "processor": platform.processor(),
        "cpu_count_logical": os.cpu_count(),
        "n_cores_budget": N_CORES,
        "env_blas": ENV_BLAS,
        "env_joblib": ENV_JOBLIB,
        "inputs_sha256": inputs,
        "inputs_digest": inputs_digest(inputs),
        "packages": dict(sorted(pkgs.items())),
        "blas_libraries": blas,
    }
    (run_dir / "provenance.json").write_text(json.dumps(prov, indent=2), encoding="utf-8")


# ---------------------------------------------------------------------------
# Manifest and freshness
# ---------------------------------------------------------------------------

def load_manifest() -> dict:
    if MANIFEST.exists():
        try:
            return json.loads(MANIFEST.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, OSError):
            print("manifest unreadable; treating as empty", file=sys.stderr)
    return {}


def save_manifest(man: dict) -> None:
    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    tmp = MANIFEST.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(man, indent=2), encoding="utf-8")
    tmp.replace(MANIFEST)  # atomic, so a kill mid-write cannot corrupt it


def prune_archives() -> None:
    if ARCHIVE_DIR.exists():
        for old in sorted(p for p in ARCHIVE_DIR.iterdir() if p.is_dir())[:-ARCHIVES_KEPT]:
            shutil.rmtree(old, ignore_errors=True)


def archive_results(run_id: str) -> Path | None:
    """Move generated outputs aside so stale files cannot satisfy readers."""
    dest = ARCHIVE_DIR / run_id
    moved = False
    for sub in oc.OUTPUT_DIRS:
        src = RESULTS / sub
        if src.exists() and any(src.iterdir()):
            dest.mkdir(parents=True, exist_ok=True)
            shutil.move(str(src), str(dest / sub))
            moved = True
    prune_archives()
    return dest if moved else None


def started_epoch(rec: dict) -> float:
    """Start of the execution a manifest record describes (records written before
    started_epoch was kept carry only the finish time and the duration)."""
    if "started_epoch" in rec:
        return rec["started_epoch"]
    # Those records stamped the finish after the warning scan and HTML export, so
    # finish - duration lands a few seconds after the real start; allow a minute.
    return rec.get("finished_epoch", 0.0) - rec.get("seconds", 0.0) - 60.0


def contract_problems(nb: str, since: float, until: float | None) -> list[str]:
    try:
        return oc.verify(nb, since=since, until=until) + oc.verify_doc_facts(nb)
    except Exception as exc:          # a contract that cannot be checked is not satisfied
        return [f"contract check raised {exc!r}"]


def reuse_decision(nb: str, man: dict, code_fp: str, inp_fp: str, profile: dict,
                   upstream_will_run: bool) -> tuple[bool, str]:
    prev = man.get(nb)
    if prev is None:
        return False, "no previous result"
    if prev.get("status") != "OK":
        return False, f"previous status {prev.get('status')}"
    if "code_digest" not in prev:
        return False, "recorded before content hashing (size+mtime key)"
    if prev.get("code_digest") != code_fp:
        return False, "code cells changed"
    if prev.get("inputs_digest") != inp_fp:
        return False, "shared inputs changed (utils/viz/requirements/Data)"
    # The profile only steers NB03, so it can only change NB03 and what consumes it.
    if ORDER.index(nb) >= ORDER.index(PROFILE_FIRST) and prev.get("profile") != profile:
        return False, "profile changed"
    if upstream_will_run:
        return False, "an upstream notebook re-runs"
    i = ORDER.index(nb)
    up = [man[u].get("finished_epoch", 0.0) for u in ORDER[:i] if u in man]
    if up and prev.get("finished_epoch", 0.0) < max(up):
        return False, "an upstream notebook finished more recently"
    # Its outputs must still be the ones it wrote, and still satisfy the contract.
    problems = contract_problems(nb, started_epoch(prev), prev.get("finished_epoch"))
    if problems:
        more = f" (+{len(problems) - 1} more)" if len(problems) > 1 else ""
        return False, f"outputs fail the contract: {problems[0]}{more}"
    return True, f"already OK at {prev.get('finished', '?')}"


# ---------------------------------------------------------------------------
# Execution
# ---------------------------------------------------------------------------

def run_one(nb: str, env_extra: dict[str, str], run_dir: Path) -> tuple[int, float, Path]:
    """Execute one notebook out of place; return (exit code, peak RSS MB, executed copy)."""
    import psutil

    env = os.environ.copy()
    # No MPLBACKEND: the kernel uses the inline backend, so executed notebooks (and their
    # HTML exports) carry their figures. Forcing Agg -- inherited from the bash runner --
    # left them figure-less and emitted ~40 "FigureCanvasAgg is non-interactive" warnings
    # per run. Verified safe first: none of the 19 pyplot-level plt.savefig() calls comes
    # after a plt.show() or in a different cell from its figure's creation.
    env.pop("MPLBACKEND", None)
    env.update(env_extra)
    stem = Path(nb).stem
    executed = run_dir / f"{stem}.ipynb"
    cmd = [
        sys.executable, "-m", "nbconvert",
        "--to", "notebook", "--execute",
        "--output-dir", str(run_dir), "--output", stem,
        "--ExecutePreprocessor.timeout=-1",
        "--ExecutePreprocessor.kernel_name=python3",
        nb,
    ]
    # nbconvert executes with the INPUT notebook's directory as the kernel's cwd, so
    # the notebooks' relative imports and paths behave exactly as before.
    with open(run_dir / f"{stem}.log", "w", encoding="utf-8", errors="replace") as fh:
        proc = subprocess.Popen(cmd, cwd=HERE, env=env, stdout=fh, stderr=subprocess.STDOUT)
        ps = psutil.Process(proc.pid)
        peak = 0.0
        while proc.poll() is None:
            try:
                procs = [ps] + ps.children(recursive=True)
                peak = max(peak, sum(p.memory_info().rss for p in procs if p.is_running()) / 1e6)
            except psutil.Error:
                pass
            time.sleep(0.5)
    return proc.returncode, peak, executed


def export_html(executed: Path) -> None:
    """Best-effort HTML copy for deposit; never fails the run."""
    if not executed.exists():
        return
    subprocess.run([sys.executable, "-m", "nbconvert", "--to", "html",
                    "--output-dir", str(executed.parent), str(executed)],
                   cwd=HERE, capture_output=True, text=True)


def scan_notebook(nb: Path) -> str:
    """Surface warnings/tracebacks from cell output; nbconvert's log omits them."""
    if not nb.exists():
        return "(no executed copy -- nbconvert failed before writing; see the .log)"
    # UTF-8 end to end: decoding the child's output with the Windows locale codepage
    # turned every em dash in a notebook message into mojibake.
    out = subprocess.run([sys.executable, str(HERE / "pipeline_tools.py"), "scan", str(nb)],
                         capture_output=True, text=True, encoding="utf-8", errors="replace",
                         cwd=HERE, env={**os.environ, "PYTHONIOENCODING": "utf-8"})
    return out.stdout.strip() or "(no warnings)"


# Scan content that marks a result as numerically suspect. The notebook still counts as
# having run (its result is reusable), but the summary says WARN, never a clean OK:
# NB09 printed "OK" over 37 non-converged GLASSO permutation fits.
WARN_MARKERS = ("ConvergenceWarning", "RuntimeWarning", "!!!")


def warn_flag(scan: str) -> bool:
    return any(m in scan for m in WARN_MARKERS)


# ---------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fresh", action="store_true",
                    help="recompute the selected notebooks (other records are kept)")
    ap.add_argument("--only", nargs="+", metavar="STEM", help="run only these, by stem prefix")
    ap.add_argument("--from", dest="from_", metavar="STEM", help="start at this notebook")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--smoke", action="store_true",
                    help="shake-out run: NB03 capped to 3 targets, N_ITER=2, N_PERM=19. "
                         "Exercises every notebook fast; NOT publication output.")
    ap.add_argument("--skip-checks", action="store_true",
                    help="skip pre-flight checks.py (not advised)")
    args = ap.parse_args()

    selected = list(NOTEBOOKS)
    if args.only:
        selected = [n for n in selected if any(n[0].startswith(s) for s in args.only)]
    if args.from_:
        idx = next((i for i, n in enumerate(selected) if n[0].startswith(args.from_)), None)
        if idx is None:
            return print(f"--from {args.from_}: no such notebook") or 2
        selected = selected[idx:]
    if not selected:
        return print("nothing selected") or 2

    missing = [n for n, *_ in selected if not (HERE / n).exists()]
    if missing:
        print("MISSING notebooks, not started:")
        for m in missing:
            print(f"  - {m}")
        return 1

    if args.smoke:
        os.environ.update(SMOKE_PROFILE)         # run_one copies os.environ per notebook
    profile = current_profile()
    mode = "SMOKE" if args.smoke else ("CUSTOM" if any(profile.values()) else "FULL")

    man = load_manifest()
    if args.fresh:
        # Clear ONLY the selected notebooks. The first version replaced the whole
        # manifest, so `--only 08 09 --fresh` erased the records for 01-06.
        for nb, *_ in selected:
            man.pop(nb, None)

    inputs = shared_inputs()
    inp_fp = inputs_digest(inputs)

    # Plan every reuse decision up front, so --dry-run shows exactly what will run.
    plan, upstream_will_run = [], False
    for nb, fatal, env_extra, note in selected:
        ok, why = reuse_decision(nb, man, code_digest(HERE / nb), inp_fp, profile,
                                 upstream_will_run)
        plan.append((nb, fatal, env_extra, note, ok, why))
        upstream_will_run = upstream_will_run or not ok

    if args.dry_run:
        print(f"mode {mode}; inputs digest {inp_fp[:12]}")
        print(f"{'notebook':<44} {'action':<7} {'OMP':<4} reason")
        for nb, fatal, env_extra, note, ok, why in plan:
            print(f"{nb:<44} {'reuse' if ok else 'RUN':<7} "
                  f"{env_extra['OMP_NUM_THREADS']:<4} {why}")
        return 0

    git = git_info()
    sha = (git.get("commit") or "nocommit")[:7]
    dirty_tag = {True: "-dirty", False: "", "unknown": "-unknown"}[git["dirty"]]
    run_id = f"{datetime.now(timezone.utc):%Y%m%dT%H%M%SZ}_{sha}{dirty_tag}"
    run_dir = RUNS_DIR / run_id
    run_dir.mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    summary = run_dir / "summary.txt"
    write_provenance(run_dir, run_id, mode, profile, inputs, git)

    def emit(line: str) -> None:
        print(line)
        with open(summary, "a", encoding="utf-8") as fh:
            fh.write(line + "\n")

    def publish_summary() -> None:
        shutil.copyfile(summary, LATEST_SUMMARY)

    emit(f"v2 pipeline run {run_id}  started {datetime.now():%F %T}")
    emit(f"host: {N_CORES} physical cores, python {sys.version.split()[0]}; "
         f"commit {git.get('commit') or '?'}  dirty={git['dirty']}")
    emit(f"mode: {mode}   profile: "
         + (", ".join(f"{k}={v}" for k, v in profile.items() if v) or "defaults")
         + f"   inputs {inp_fp[:12]}")
    if mode == "SMOKE":
        emit("!! SMOKE RUN -- shake-out only. Outputs in Results/ are NOT publication "
             "results (NB03 capped to 3 targets, N_PERM=19). !!")
        SMOKE_MARKER.parent.mkdir(parents=True, exist_ok=True)
        SMOKE_MARKER.write_text(
            f"Results/ was last written by a SMOKE run ({run_id}).\n"
            "NB03 benchmarked only 3 targets with N_ITER=2 and N_PERM=19, and every\n"
            "downstream notebook consumed those outputs. Nothing here is publication\n"
            "output. A full run on a clean tree recomputes everything and deletes this\n"
            "file on success.\n", encoding="utf-8")

    # Pre-flight. A broken group definition, a leaked covariate or an inverted effect
    # sign must not be able to produce results that look finished.
    if not args.skip_checks:
        print(f"[{datetime.now():%F %T}] pre-flight checks")
        chk = run_dir / "checks.log"
        with open(chk, "w", encoding="utf-8", errors="replace") as fh:
            rc = subprocess.run([sys.executable, "checks.py"], cwd=HERE,
                                stdout=fh, stderr=subprocess.STDOUT).returncode
        text = chk.read_text(encoding="utf-8", errors="replace")
        if rc != 0:
            emit("PRE-FLIGHT CHECKS FAILED -- pipeline not started.")
            emit(f"See {chk}. Failing checks:")
            for ln in text.splitlines():
                if ln.startswith("FAIL") or ln.startswith("  - "):
                    emit("  " + ln.strip())
            publish_summary()
            return 1
        emit(f"{'checks.py':<44} {'OK':<4} "
             f"{sum(1 for ln in text.splitlines() if ln.startswith('ok'))} checks passed")

    # Freshness: a run that recomputes NB01 starts from an empty Results/.
    if any(nb == ORDER[0] and not ok for nb, *_r, ok, _w in plan):
        dest = archive_results(run_id)
        if dest:
            emit(f"previous outputs archived to {dest.relative_to(ROOT)}")

    for nb, fatal, env_extra, _note, ok, why in plan:
        path = HERE / nb
        if ok:
            prev = man[nb]
            tag = "; result carries WARN" if prev.get("warn") else ""
            emit(f"{nb:<44} SKIP  ({why}; run {prev.get('run_id', '?')}{tag})")
            for ln in (prev.get("scan") or "(no warnings)").splitlines():
                emit("    " + ln)
            continue

        print(f"[{datetime.now():%F %T}] start {nb}  ({why})")
        # Whatever its contract declares must come from THIS execution: move the
        # previous copies aside (a no-op after NB01 archived the whole tree).
        try:
            n_prev = oc.clear_outputs(nb, ARCHIVE_DIR / run_id)
        except oc.ContractError as exc:
            emit(f"{nb:<44} FAIL  contract is ambiguous: {exc}")
            publish_summary()
            return 1
        if n_prev:
            prune_archives()
        start = time.time()
        rc, peak, executed = run_one(nb, env_extra, run_dir)
        end = time.time()
        secs = end - start
        contract = contract_problems(nb, start, end) if rc == 0 else []
        extra = oc.undeclared(nb, start, end) if rc == 0 else []
        status = "OK" if rc == 0 and not contract else "FAIL"
        scan = scan_notebook(executed)
        export_html(executed)

        warn = status == "OK" and warn_flag(scan)
        man[nb] = {"status": status, "warn": warn, "seconds": round(secs, 1),
                   "peak_mb": round(peak), "finished": f"{datetime.now():%F %T}",
                   "started_epoch": start, "finished_epoch": end,
                   "code_digest": code_digest(path),
                   "inputs_digest": inp_fp, "profile": profile, "mode": mode,
                   "run_id": run_id, "scan": scan, "thread_env": env_extra,
                   "contract": contract, "undeclared_outputs": extra}
        save_manifest(man)

        shown = "WARN" if warn else status
        emit(f"{nb:<44} {shown:<4} {secs:7.0f}s  peak {peak:6.0f} MB"
             + (f"  ({n_prev} previous output(s) archived)" if n_prev else ""))
        for ln in scan.splitlines():
            emit("    " + ln)
        for p in contract:
            emit("    CONTRACT: " + p)
        if extra:
            emit(f"    contract: {len(extra)} output(s) not declared in output_contract.py: "
                 + ", ".join(extra))

        if status != "OK" and fatal:
            cause = ("its outputs break the contract (CONTRACT lines above)" if rc == 0
                     else f"see {run_dir / (path.stem + '.log')}")
            emit(f"FATAL: {nb} failed -- {cause}; stopping.")
            emit("Fix the cause, then re-run; valid results will be reused.")
            publish_summary()
            return 1

    # Anything downstream of what just ran, but not selected, is now stale.
    last = max(ORDER.index(nb) for nb, *_ in selected)
    stale = [nb for nb in ORDER[last + 1:] if nb in man
             and man[nb].get("finished_epoch", 0.0) < max(man[u].get("finished_epoch", 0.0)
                                                         for u in ORDER[:ORDER.index(nb)]
                                                         if u in man)]
    if stale:
        emit("STALE (an upstream notebook re-ran after them): " + ", ".join(stale))
        emit(f"  re-run with: python run_pipeline.py --from {Path(stale[0]).stem[:2]}")

    emit(f"pipeline finished {datetime.now():%F %T}  (mode: {mode})")
    first = ORDER.index(PROFILE_FIRST)
    all_ok = all(man.get(n, {}).get("status") == "OK"
                 and (i < first or man.get(n, {}).get("profile") == profile)
                 for i, n in enumerate(ORDER))
    if mode == "FULL" and all_ok and not stale and SMOKE_MARKER.exists():
        if git["dirty"] is False:
            SMOKE_MARKER.unlink()
            emit("Full run complete on a clean tree: removed Results/SMOKE_RUN.txt.")
        else:
            emit(f"Full run complete, but the tree is dirty={git['dirty']}: "
                 "Results/SMOKE_RUN.txt kept. Commit, then re-run for a baseline.")
    elif mode == "SMOKE":
        emit("SMOKE RUN: Results/ holds shake-out outputs only -- see Results/SMOKE_RUN.txt.")
    emit(f"executed notebooks, HTML, logs and provenance: {run_dir.relative_to(ROOT)}")
    publish_summary()
    return 0


if __name__ == "__main__":
    sys.exit(main())
