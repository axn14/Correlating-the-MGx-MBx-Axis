#!/usr/bin/env python3
"""Execute the v2 pipeline headless, in dependency order (see CONTRACTS.md).

Replaces run_pipeline.sh. The analysis server this project was built for is no
longer available, so the pipeline now runs on the local Windows workstation and
the orchestrator has to be cross-platform. Only stdlib + psutil are used, which
is the same footprint as pipeline_tools.py.

    python run_pipeline.py                  # full run, resuming where it left off
    python run_pipeline.py --fresh          # ignore the manifest, re-run everything
    python run_pipeline.py --only 02 03     # just these, by stem prefix
    python run_pipeline.py --from 06        # this notebook onward
    python run_pipeline.py --dry-run        # print the plan and exit
    python run_pipeline.py --smoke          # fast shake-out pass, NOT publication output

Why a manifest: a notebook that has already succeeded is skipped on the next
invocation. The previous bash runner had no such record, so any failure meant
re-running every completed notebook from scratch -- the stale/looped failure mode
NB03 made expensive. Delete Results/logs/pipeline_manifest.json, or pass --fresh,
to force a clean pass.

Thread policy. Oversubscription is the real risk on 8 physical cores: BLAS
threads multiply with joblib workers, so 8 x 8 = 64 threads thrash. Each notebook
therefore gets its own thread environment below. Notebooks whose cost is linear
algebra (preprocessing, MOFA) get the cores through BLAS; notebooks that
parallelise at the joblib level (NB03's nested CV) get BLAS=1 and spend the cores
on workers instead, via KI_N_JOBS.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
LOG_DIR = HERE.parent / "Results" / "logs"
SUMMARY = LOG_DIR / "summary_v2.txt"
MANIFEST = LOG_DIR / "pipeline_manifest.json"

# 8 physical cores / 16 logical on the Ryzen 7 6800H. Compute-bound tree building
# does not benefit much from SMT, so the budget is physical cores.
N_CORES = 8

# Linear-algebra notebooks: give the cores to BLAS.
ENV_BLAS = {
    "OMP_NUM_THREADS": str(N_CORES),
    "MKL_NUM_THREADS": str(N_CORES),
    "OPENBLAS_NUM_THREADS": str(N_CORES),
    "NUMEXPR_MAX_THREADS": str(N_CORES),
    "LOKY_MAX_CPU_COUNT": str(N_CORES),
}
# joblib-parallel notebooks: one thread per fit, cores spent on workers.
ENV_JOBLIB = {
    "OMP_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "NUMEXPR_MAX_THREADS": "1",
    "LOKY_MAX_CPU_COUNT": str(N_CORES),
    "KI_N_JOBS": str(N_CORES),
}

# (notebook, fatal?, thread env, note)
# NB04 is absent by design: discarded, its tumour-site exercise moved to NB05
# stage 2 and its stage trends to NB02. NB07 is markdown-only. NB10 (MICOM) runs
# separately via run_micom.py -- it is a long solve, not part of this sequence.
NOTEBOOKS: list[tuple[str, bool, dict[str, str], str]] = [
    ("01_preprocessing.ipynb", True, ENV_BLAS, "builds the reduced pickle"),
    ("01b_yachida_significant_features.ipynb", False, ENV_BLAS, "nothing downstream reads it"),
    # As on disk NB02 writes per-comparison DA CSVs, which is what NB03 reads. The tidy
    # differential_abundance_v2.csv exists only in edits_nb02b.py, not yet applied.
    ("02_analysis.ipynb", True, ENV_BLAS, "writes per-comparison DA tables"),
    ("03_ml_benchmarking.ipynb", True, ENV_JOBLIB, "nested CV + permutation"),
    # NB05 does not read KI_N_JOBS (its LightGBM fits pin n_jobs themselves), so the
    # joblib profile would only starve its sklearn/BLAS work of threads.
    ("05_two_stage_model.ipynb", False, ENV_BLAS, "leaf analysis + transfer"),
    ("06_multiomics_factors.ipynb", True, ENV_BLAS, "NB08/NB09 need its tables"),
    ("08_pathway_analysis.ipynb", True, ENV_BLAS, "writes the evidence table"),
    ("09_pathway_network.ipynb", True, ENV_BLAS, "reads it"),
]


# Environment knobs that change RESULTS, not merely speed. A manifest entry is
# reused only if it was produced under the same profile, so a smoke run's capped
# NB03 can never satisfy a full run. (KI_N_JOBS is excluded: with fixed seeds and
# single-threaded fits it changes wall time, not output.)
PROFILE_KEYS = ("NB03_N_PERM", "NB03_N_ITER", "NB03_MAX_TARGETS")

# Shake-out profile: exercises every code path end to end in minutes. NOT
# publication output -- N_PERM=19 floors the permutation p at 0.05, and only three
# targets are benchmarked.
SMOKE_PROFILE = {"NB03_N_PERM": "19", "NB03_N_ITER": "2", "NB03_MAX_TARGETS": "3"}
SMOKE_MARKER = HERE.parent / "Results" / "SMOKE_RUN.txt"


def current_profile() -> dict:
    return {k: os.environ.get(k, "") for k in PROFILE_KEYS}


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


def notebook_fingerprint(path: Path) -> str:
    """Source size + mtime. If a notebook is edited its manifest entry is stale,
    so a resume re-runs it rather than trusting an output from older code."""
    st = path.stat()
    return f"{st.st_size}:{int(st.st_mtime)}"


def run_one(nb: str, env_extra: dict[str, str], peak_file: Path, log_file: Path) -> tuple[int, float]:
    """Execute one notebook in place, returning (exit code, peak RSS in MB)."""
    import psutil

    env = os.environ.copy()
    env["MPLBACKEND"] = "Agg"
    env.update(env_extra)

    cmd = [
        sys.executable, "-m", "nbconvert",
        "--to", "notebook", "--execute", "--inplace",
        "--ExecutePreprocessor.timeout=-1",
        "--ExecutePreprocessor.kernel_name=python3",
        nb,
    ]
    with open(log_file, "w", encoding="utf-8", errors="replace") as fh:
        proc = subprocess.Popen(cmd, cwd=HERE, env=env, stdout=fh,
                                stderr=subprocess.STDOUT)
        ps = psutil.Process(proc.pid)
        peak = 0.0
        while proc.poll() is None:
            try:
                procs = [ps] + ps.children(recursive=True)
                peak = max(peak, sum(p.memory_info().rss
                                     for p in procs if p.is_running()) / 1e6)
            except psutil.Error:
                pass
            time.sleep(0.5)
    peak_file.write_text(f"{peak:.0f}\n", encoding="utf-8")
    return proc.returncode, peak


def scan_notebook(nb: Path) -> str:
    """Surface warnings/tracebacks from cell output; nbconvert's log omits them."""
    out = subprocess.run([sys.executable, str(HERE / "pipeline_tools.py"), "scan", str(nb)],
                         capture_output=True, text=True, cwd=HERE)
    return out.stdout.strip() or "(no warnings)"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fresh", action="store_true", help="ignore the manifest")
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

    if args.dry_run:
        print(f"{'notebook':<44} {'fatal':<6} {'OMP':<5} note")
        for nb, fatal, env, note in selected:
            print(f"{nb:<44} {str(fatal):<6} {env['OMP_NUM_THREADS']:<5} {note}")
        return 0

    if args.smoke:
        os.environ.update(SMOKE_PROFILE)       # run_one copies os.environ per notebook
    profile = current_profile()
    mode = "SMOKE" if args.smoke else ("CUSTOM" if any(profile.values()) else "FULL")

    LOG_DIR.mkdir(parents=True, exist_ok=True)
    man = {} if args.fresh else load_manifest()

    lines = [f"v2 pipeline run started {datetime.now():%F %T}",
             f"host: {N_CORES} physical cores, python {sys.version.split()[0]}",
             f"mode: {mode}   profile: "
             + (", ".join(f"{k}={v}" for k, v in profile.items() if v) or "defaults")]
    if mode == "SMOKE":
        lines.append("!! SMOKE RUN -- shake-out only. Outputs in Results/ are NOT publication "
                     "results (NB03 capped to 3 targets, N_PERM=19). !!")
    SUMMARY.write_text("\n".join(lines) + "\n", encoding="utf-8")
    if mode == "SMOKE":
        SMOKE_MARKER.parent.mkdir(parents=True, exist_ok=True)
        SMOKE_MARKER.write_text(
            f"Results/ was last written by a SMOKE run ({datetime.now():%F %T}).\n"
            "NB03 benchmarked only 3 targets with N_ITER=2 and N_PERM=19, and every\n"
            "downstream notebook consumed those outputs. Nothing here is publication\n"
            "output. A full run (python run_pipeline.py) recomputes everything and\n"
            "deletes this file on success.\n", encoding="utf-8")

    order = [n for n, *_ in NOTEBOOKS]

    def upstream_latest(nb: str) -> float:
        """Most recent finish time of any notebook upstream of `nb` in the manifest."""
        i = order.index(nb)
        t = [man[u].get("finished_epoch", 0.0) for u in order[:i] if u in man]
        return max(t) if t else 0.0

    def emit(line: str) -> None:
        print(line)
        with open(SUMMARY, "a", encoding="utf-8") as fh:
            fh.write(line + "\n")

    # Pre-flight. A broken group definition, a leaked covariate or an inverted
    # effect sign must not be able to produce results that look finished.
    if not args.skip_checks:
        print(f"[{datetime.now():%F %T}] pre-flight checks")
        chk = LOG_DIR / "checks.log"
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
            return 1
        emit(f"{'checks.py':<44} {'OK':<4} "
             f"{sum(1 for ln in text.splitlines() if ln.startswith('ok'))} checks passed")

    for nb, fatal, env_extra, _note in selected:
        path = HERE / nb
        stem = path.stem
        fp = notebook_fingerprint(path)
        prev = man.get(nb)
        # Reuse a previous result only if ALL of these hold:
        #   - it succeeded, from this exact notebook source (fingerprint);
        #   - under the same result-affecting profile (a smoke run never satisfies
        #     a full one, or vice versa);
        #   - it finished AFTER every upstream notebook's latest run. Each notebook
        #     consumes upstream outputs, so if anything upstream was re-run since,
        #     this result was computed from inputs that no longer exist. Checking
        #     only the notebook's own fingerprint -- the first version of this
        #     runner -- would have reused NB03 on top of a freshly re-run NB02.
        reusable = (prev is not None and prev.get("status") == "OK"
                    and prev.get("fingerprint") == fp
                    and prev.get("profile") == profile
                    and prev.get("finished_epoch", 0.0) >= upstream_latest(nb))
        if reusable:
            emit(f"{nb:<44} SKIP  (already OK at {prev.get('finished', '?')})")
            continue

        print(f"[{datetime.now():%F %T}] start {nb}")
        start = time.time()
        rc, peak = run_one(nb, env_extra, LOG_DIR / f"{stem}.peak", LOG_DIR / f"{stem}.log")
        secs = time.time() - start
        status = "OK" if rc == 0 else "FAIL"

        # Fingerprint again: --inplace rewrites the file, so the pre-run value is
        # stale and would make the next resume re-run a notebook that succeeded.
        man[nb] = {"status": status, "seconds": round(secs, 1), "peak_mb": round(peak),
                   "finished": f"{datetime.now():%F %T}", "finished_epoch": time.time(),
                   "fingerprint": notebook_fingerprint(path), "profile": profile,
                   "mode": mode}
        save_manifest(man)

        emit(f"{nb:<44} {status:<4} {secs:7.0f}s  peak {peak:6.0f} MB")
        for ln in scan_notebook(path).splitlines():
            emit("    " + ln)

        if rc != 0 and fatal:
            emit(f"FATAL: {nb} failed -- see {LOG_DIR / (stem + '.log')}; stopping.")
            emit("Fix the cause, then re-run; completed notebooks will be skipped.")
            return 1

    emit(f"pipeline finished {datetime.now():%F %T}  (mode: {mode})")
    all_ok = all(man.get(n, {}).get("status") == "OK"
                 and man.get(n, {}).get("profile") == profile for n, *_ in selected)
    if mode == "FULL" and all_ok and len(selected) == len(NOTEBOOKS) and SMOKE_MARKER.exists():
        SMOKE_MARKER.unlink()
        emit("Full run complete: removed Results/SMOKE_RUN.txt.")
    elif mode == "SMOKE":
        emit("SMOKE RUN: Results/ holds shake-out outputs only -- see Results/SMOKE_RUN.txt.")
    emit("")
    emit("NB10 (MICOM) is not part of this run: it is a long solve against the")
    emit("pre-built AGORA2 community models. Run it separately, then re-run 08 and 09")
    emit("to fill the mechanistic evidence column:")
    emit("  python run_pipeline.py --only 08 09 --fresh")
    return 0


if __name__ == "__main__":
    sys.exit(main())
