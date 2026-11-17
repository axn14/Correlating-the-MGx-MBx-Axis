"""Peak memory of ONE fresh batch solve, to size Dardel's --mem (micom_batch.sbatch).

The alpha sweep's workers reached 6-9 GB, but each ran ten solves in one process; a batch
task runs exactly one. This runs the batch's own step -- load the community model, set
the medium, `forced_secretion.forced_exports` at the chosen alpha -- in a fresh process and
records the peak working set (Windows `peak_wset`: everything the process touched,
including HiGHS's native allocations, which tracemalloc cannot see).

Run:  python v2/micom/memory_probe.py <sample_id> [alpha]
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import pandas as pd
import psutil

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))

import highs_qp  # noqa: E402

highs_qp.install("repair")


def main() -> None:
    from micom import load_pickle
    import forced_secretion as fs
    sid = sys.argv[1]
    alpha = float(sys.argv[2]) if len(sys.argv) > 2 else \
        json.loads((HERE / "out" / "alpha_sweep_choice.json").read_text())["chosen_alpha"]
    proc = psutil.Process()
    t0 = time.time()
    reg = pd.read_csv(HERE.parent / "targets" / "out" / "target_registry.tsv", sep="\t")
    targets = set(reg.loc[reg["export_gate"] == True, "kegg"].dropna())  # noqa: E712
    com = load_pickle(str(HERE / "out" / "community" / "models" / f"{sid}.pickle"))
    after_load = proc.memory_info().peak_wset
    com.solver.configuration.presolve = True
    d = pd.read_csv(HERE / "gut_diet_v2.csv", comment="#")
    med = d.set_index(d["reaction"].str.replace(r"_e$", "_m", regex=True))["flux"]
    com.medium = med[med.index.isin([r.id for r in com.exchanges])]
    df = fs.forced_exports(com, alpha, targets)
    res = {"sample_id": sid, "alpha": alpha, "n_taxa": len(com.taxa),
           "peak_gb_after_load": round(after_load / 1e9, 3),
           "peak_gb_total": round(proc.memory_info().peak_wset / 1e9, 3),
           "seconds": round(time.time() - t0, 1),
           "gap_rel": highs_qp.LAST_RUN.get("gap_rel"),
           "n_candidates": df.attrs.get("n_candidates"), "n_lp": df.attrs.get("n_lp"),
           "n_forced": int(df["forced"].sum()) if "forced" in df else None}
    out = HERE / "out" / "memory_probe.jsonl"
    with open(out, "a") as f:
        f.write(json.dumps(res) + "\n")
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
