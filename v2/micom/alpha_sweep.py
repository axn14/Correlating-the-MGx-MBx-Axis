"""Phase E2: the cooperative-tradeoff sweep, with the certified solver.

The 0.5 tradeoff in NB10 was chosen on 10 pilot samples solved with OSQP, which the
solver check (README finding 7) showed was not converged even where "valid". This re-asks
the question with the certified solver (`highs_qp.install("repair")`, gap <= 1%) and adds
the review's criterion: the per-taxon export calls for the targets must be STABLE in alpha,
or the attribution depends on a free parameter.

Design, fixed before the first run (2026-11-11):
- Samples: the 8 community models `diet_screen.py --community 8` built (out/community/).
- Grid: alpha = 0.1, 0.2, ..., 1.0. Medium: gut_diet_v2 (target_medium.py withdrew nothing).
- Per (sample, alpha): validity (`micom_tools.solution_quality`), fraction of taxa growing,
  certified gap, and the set of taxa with net export > tolerance for each of the 26
  export-gated targets (`v2/targets/out/target_registry.tsv`).
- Stability of alpha: mean over (sample, target) of the Jaccard index between the exporter
  set at alpha and at each neighbour (alpha +/- 0.1), over pairs where either set is
  non-empty.
- RULE: the chosen alpha is the LARGEST alpha with (i) median fraction of taxa growing
  >= 0.95 and (ii) stability >= 0.80. If none qualifies, 0.5 is kept and the table says so.

Run:  python v2/micom/alpha_sweep.py --run       (~10 h wall on the laptop, 4 workers)
      python v2/micom/alpha_sweep.py --analyse   (seconds; writes the choice)
Outputs: out/alpha_sweep/<sample>_<alpha>.pkl, out/alpha_sweep_summary.tsv,
         out/alpha_sweep_stability.tsv, out/alpha_sweep_choice.json
"""
from __future__ import annotations

import argparse
import json
import pickle
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))

import highs_qp  # noqa: E402

highs_qp.install("repair")          # module level: applies in spawned workers too

import micom_tools as mt            # noqa: E402

WORK = HERE / "out" / "community"
SWEEP = HERE / "out" / "alpha_sweep"
REGISTRY = HERE.parent / "targets" / "out" / "target_registry.tsv"
GRID = [round(a / 10, 1) for a in range(1, 11)]
EXPORT_TOL = 1e-6
MIN_GROWING, MIN_STABILITY = 0.95, 0.80


def medium() -> pd.DataFrame:
    d = pd.read_csv(HERE / "gut_diet_v2.csv", comment="#")
    m = d[["reaction", "flux"]].copy()
    m["reaction"] = m["reaction"].str.replace(r"_e$", "_m", regex=True)
    return m


def solve_sample(sid: str) -> str:
    from micom.workflows import grow
    man = pd.read_csv(WORK / "manifest.csv", dtype={"sample_id": str})
    row = man[man["sample_id"] == sid]
    med = medium()
    SWEEP.mkdir(parents=True, exist_ok=True)
    for a in GRID:
        out = SWEEP / f"{sid}_{a:.1f}.pkl"
        if out.exists():
            continue
        t0 = time.time()
        try:
            res = grow(row, model_folder=str(WORK / "models"), medium=med, tradeoff=a,
                       threads=1, strategy="pFBA", presolve=True)
            rec = {"exchanges": res.exchanges, "growth_rates": res.growth_rates,
                   "annotations": res.annotations, "ok": True}
        except Exception as err:
            rec = {"ok": False, "error": str(err)}
        rec.update({"sample_id": sid, "alpha": a, "seconds": time.time() - t0,
                    "solver": dict(highs_qp.LAST_RUN)})
        with open(out, "wb") as f:
            pickle.dump(rec, f)
        print(f"  {sid} alpha={a:.1f}: {rec['seconds']:.0f}s ok={rec['ok']} "
              f"gap={rec['solver'].get('gap_rel', float('nan')):.3g}", flush=True)
    return sid


def run() -> None:
    man = pd.read_csv(WORK / "manifest.csv", dtype={"sample_id": str})
    sids = man["sample_id"].tolist()
    print(f"sweeping {len(sids)} samples x {len(GRID)} alphas with "
          f"{highs_qp.describe()['qp_solver']}", flush=True)
    with ProcessPoolExecutor(max_workers=4) as ex:
        for sid in ex.map(solve_sample, sids):
            print(f"done {sid}", flush=True)


def exporters(rec: dict, kegg_targets: set[str]) -> dict[str, set[str]]:
    """{target KEGG -> taxa with net export > tol} from one solve."""
    ex = rec["exchanges"].copy()
    ann = rec["annotations"]
    kcol = next((c for c in ann.columns if c.lower() in ("kegg", "kegg.compound")), None)
    met_kegg = {}
    if kcol:
        for met, k in zip(ann["metabolite"], ann[kcol].astype(str)):
            for x in k.split(";"):
                if x.strip() in kegg_targets:
                    met_kegg[met] = x.strip()
    ex = ex[ex["taxon"].astype(str) != "medium"]
    ex = ex[ex["flux"] > EXPORT_TOL]
    out: dict[str, set[str]] = {k: set() for k in kegg_targets}
    for met, taxon in zip(ex["metabolite"], ex["taxon"].astype(str)):
        k = met_kegg.get(met)
        if k:
            out[k].add(taxon)
    return out


def analyse() -> None:
    reg = pd.read_csv(REGISTRY, sep="\t")
    targets = set(reg.loc[reg["export_gate"] == True, "kegg"].dropna())  # noqa: E712
    recs = []
    for p in sorted(SWEEP.glob("*.pkl")):
        with open(p, "rb") as f:
            recs.append(pickle.load(f))
    rows, sets = [], {}
    for r in recs:
        base = {"sample_id": r["sample_id"], "alpha": r["alpha"], "seconds": r["seconds"],
                "ok": r["ok"], "gap_rel": r["solver"].get("gap_rel", np.nan)}
        if not r["ok"]:
            rows.append(base); continue
        q = mt.solution_quality(r["exchanges"], r["growth_rates"]).iloc[0]
        g = r["growth_rates"]
        g = g[g["taxon"].astype(str) != "medium"]
        # Clarification, 2026-11-11, written after the first two solves (alpha 0.1:
        # 10091 gap 4.4%, 10335 gap 16.9%): a solve counts as valid only if its certified
        # gap is within the 5% bound solver_check.py fixed as P3. An uncertified optimum
        # cannot qualify an alpha. Reported, not hidden: `certified` is a column.
        certified = bool(base["gap_rel"] <= 0.05) if pd.notna(base["gap_rel"]) else False
        base.update({"certified": certified,
                     "valid": bool(q["valid"]) and certified,
                     "community_growth": float(q["community_growth"]),
                     "frac_taxa_growing": float((g["growth_rate"] > mt.GROWTH_FLOOR).mean())})
        rows.append(base)
        sets[(r["sample_id"], r["alpha"])] = exporters(r, targets)
    summ = pd.DataFrame(rows).sort_values(["sample_id", "alpha"])
    summ.to_csv(HERE / "out" / "alpha_sweep_summary.tsv", sep="\t", index=False)

    stab = []
    for a in GRID:
        js = []
        for (sid, aa), s in sets.items():
            if aa != a:
                continue
            for nb in (round(a - 0.1, 1), round(a + 0.1, 1)):
                t = sets.get((sid, nb))
                if t is None:
                    continue
                for k in targets:
                    u = s[k] | t[k]
                    if u:
                        js.append(len(s[k] & t[k]) / len(u))
        sub = summ[(summ["alpha"] == a) & summ.get("valid", pd.Series(dtype=bool)).fillna(False)]
        stab.append({"alpha": a, "stability_jaccard": float(np.mean(js)) if js else np.nan,
                     "n_comparisons": len(js),
                     "median_frac_growing": float(sub["frac_taxa_growing"].median()) if len(sub) else np.nan,
                     "n_valid": int(len(sub))})
    st = pd.DataFrame(stab)
    st.to_csv(HERE / "out" / "alpha_sweep_stability.tsv", sep="\t", index=False)
    ok = st[(st["median_frac_growing"] >= MIN_GROWING) & (st["stability_jaccard"] >= MIN_STABILITY)]
    chosen = float(ok["alpha"].max()) if len(ok) else 0.5
    choice = {"run_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
              "rule": f"largest alpha with median frac growing >= {MIN_GROWING} and "
                      f"stability >= {MIN_STABILITY}; else 0.5",
              "qualifying": ok["alpha"].tolist(), "chosen_alpha": chosen,
              "fallback_used": len(ok) == 0, "solver": highs_qp.describe()}
    (HERE / "out" / "alpha_sweep_choice.json").write_text(json.dumps(choice, indent=1))
    print(st.round(3).to_string(index=False))
    print(json.dumps(choice, indent=1))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--analyse", action="store_true")
    a = ap.parse_args()
    if a.run:
        run()
    if a.analyse:
        analyse()
