"""The Gurobi-free solver, checked against OSQP on the same community models.

Context: README finding 4 -- under optlang's OSQP/HiGHS hybrid, micom's tradeoff QP ends
in a point the pFBA step cannot use, and micom returns it without raising. A Gurobi
licence was not obtainable (2026-11-11). `highs_qp.py` keeps OSQP for the QP and then
REPAIRS its point with one HiGHS LP -- the exactly feasible flux state whose taxon growth
rates are nearest (L1) to OSQP's -- then polishes it with Wolfe's minimum-norm-point
method (HiGHS LPs only) until the CERTIFIED relative optimality gap is <= 1%. The
certificate is the Frank-Wolfe gap: an upper bound on suboptimality, not an estimate.

Caveat for reading P2: OSQP stopped at its iteration limit on every sample probed so far
(10091, 10723), including ones it "solved" validly, so a difference from OSQP there is
not necessarily an error of the new solver. P2 is kept exactly as fixed below anyway.

This script re-solves the 8 community models that `diet_screen.py --community 8` built
(out/community/, same manifest, same gut_diet_v2 medium, tradeoff 0.5, pFBA, presolve),
one sample at a time so each sample's certificate is recorded, and compares with the
stored OSQP solve (out/community/grow_gut_diet_v2.pkl).

PASS CRITERIA -- fixed 2026-11-11, before this script was first run to completion:
  P1  every sample returns a solution that passes micom_tools.solution_quality.
  P2  on the samples OSQP already solved validly, the solver change moves nothing that
      matters: median per-taxon relative growth difference <= 1%, community growth ratio
      within 1e-3 of 1, and for every probe compound the set of secreting taxa is the
      same (Jaccard = 1) in at least 90% of (sample, compound) pairs with any secretor.
  P3  every sample carries a certificate, with relative optimality gap <= 5%.
  P4  median wall time <= 5 min per sample on the laptop (8 logical cores, 1 used).
A criterion that fails is reported as failed; none is relaxed after the fact.

Run:  python v2/micom/solver_check.py            (~1 h on the laptop, no HPC)
"""
from __future__ import annotations

import json
import pickle
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))

import highs_qp  # noqa: E402

highs_qp.install("repair")          # module level: also applies in spawned workers

import micom_tools as mt            # noqa: E402
from diet_screen import PROBE, _secretors   # noqa: E402

OUT = HERE / "out"
WORK = OUT / "community"
OSQP_PKL = WORK / "grow_gut_diet_v2.pkl"
NEW_PKL = WORK / "grow_gut_diet_v2_highs_repair.pkl"
DIET = HERE / "gut_diet_v2.csv"
TRADEOFF = 0.5


def medium() -> pd.DataFrame:
    d = pd.read_csv(DIET, comment="#")
    med = d[["reaction", "flux"]].copy()
    med["reaction"] = med["reaction"].str.replace(r"_e$", "_m", regex=True)
    return med


def solve_all(manifest: pd.DataFrame) -> tuple[dict, pd.DataFrame]:
    """Solve each sample separately; returns (combined frames, per-sample diagnostics)."""
    from micom.workflows import grow
    med = medium()
    ex, gr, diag = [], [], []
    for sid in manifest["sample_id"]:
        t0 = time.time()
        try:
            res = grow(manifest[manifest["sample_id"] == sid], model_folder=str(WORK / "models"),
                       medium=med, tradeoff=TRADEOFF, threads=1, strategy="pFBA",
                       presolve=True)
            ex.append(res.exchanges)
            gr.append(res.growth_rates)
            ok = True
        except Exception as err:        # micom raises when every sample failed
            ok = False
            print(f"  {sid}: grow failed: {err}", flush=True)
        d = {"sample_id": str(sid), "returned": ok, "seconds": time.time() - t0,
             **{k: v for k, v in highs_qp.LAST_RUN.items() if k != "mode"}}
        diag.append(d)
        print(f"  {sid}: {d['seconds']:.0f}s  osqp={d.get('osqp_status')}  "
              f"gap_rel={d.get('gap_rel', float('nan')):.3g}", flush=True)
    frames = {"exchanges": pd.concat(ex, ignore_index=True),
              "growth_rates": pd.concat(gr, ignore_index=True)}
    return frames, pd.DataFrame(diag)


def per_taxon(gr: pd.DataFrame) -> pd.Series:
    g = gr[gr["taxon"].astype(str) != "medium"].copy()
    g["sample_id"] = g["sample_id"].astype(str)
    return g.set_index(["sample_id", "taxon"])["growth_rate"]


def main() -> None:
    manifest = pd.read_csv(WORK / "manifest.csv", dtype={"sample_id": str})
    with open(OSQP_PKL, "rb") as f:
        old = pickle.load(f)
    if NEW_PKL.exists():
        with open(NEW_PKL, "rb") as f:
            new, diag = pickle.load(f)
        print(f"reusing {NEW_PKL.name}")
    else:
        print(f"solving {len(manifest)} samples with {highs_qp.describe()['qp_solver']}")
        new, diag = solve_all(manifest)
        with open(NEW_PKL, "wb") as f:
            pickle.dump((new, diag), f)

    q_old = mt.solution_quality(old.exchanges, old.growth_rates).set_index("sample_id")
    q_new = mt.solution_quality(new["exchanges"], new["growth_rates"]).set_index("sample_id")
    per = q_old.add_suffix("_osqp").join(q_new.add_suffix("_new"), how="outer")
    per = per.join(diag.set_index("sample_id"), how="left")
    per["growth_ratio"] = per["community_growth_new"] / per["community_growth_osqp"]

    g_old, g_new = per_taxon(old.growth_rates), per_taxon(new["growth_rates"])
    both_g = pd.concat([g_old.rename("osqp"), g_new.rename("new")], axis=1)
    both_g["abs_diff"] = (both_g["new"] - both_g["osqp"]).abs()
    both_g["rel_diff"] = both_g["abs_diff"] / both_g["osqp"].where(both_g["osqp"] > mt.GROWTH_FLOOR)
    agg = both_g.groupby(level="sample_id").agg(
        n_taxa=("new", "size"),
        n_growing_osqp=("osqp", lambda s: int((s > mt.GROWTH_FLOOR).sum())),
        n_growing_new=("new", lambda s: int((s > mt.GROWTH_FLOOR).sum())),
        max_abs_growth_diff=("abs_diff", "max"),
        median_rel_growth_diff=("rel_diff", "median"))
    per = per.join(agg)
    per.index.name = "sample_id"
    per.to_csv(OUT / "solver_check_per_sample.tsv", sep="\t")
    both_g.to_csv(OUT / "solver_check_growth.tsv", sep="\t")

    valid_osqp = sorted(per.index[per["valid_osqp"].fillna(False).astype(bool)])
    ex_old = old.exchanges.copy()
    ex_old["sample_id"] = ex_old["sample_id"].astype(str)
    ex_new = new["exchanges"].copy()
    ex_new["sample_id"] = ex_new["sample_id"].astype(str)
    srows = []
    for sid in per.index:
        for b in PROBE:
            so = set(_secretors(ex_old[ex_old["sample_id"] == sid], b)["taxon"])
            sn = set(_secretors(ex_new[ex_new["sample_id"] == sid], b)["taxon"])
            union = so | sn
            srows.append({"sample_id": sid, "compound": b, "osqp_valid": sid in valid_osqp,
                          "n_osqp": len(so), "n_new": len(sn),
                          "jaccard": (len(so & sn) / len(union)) if union else np.nan,
                          "only_osqp": ";".join(sorted(so - sn)),
                          "only_new": ";".join(sorted(sn - so))})
    sec = pd.DataFrame(srows)
    sec.to_csv(OUT / "solver_check_secretion.tsv", sep="\t", index=False)

    # ---- the pre-registered criteria --------------------------------------------------
    v = per.loc[valid_osqp]
    sv = sec[sec["osqp_valid"] & sec["jaccard"].notna()]
    crit = {
        "P1_all_valid": bool(per["valid_new"].fillna(False).all()),
        "P2_median_rel_growth_diff": float(v["median_rel_growth_diff"].median()),
        "P2_growth_ratio_max_dev": float((v["growth_ratio"] - 1).abs().max()),
        "P2_secretion_identical_share": float((sv["jaccard"] == 1).mean()) if len(sv) else np.nan,
        "P3_max_gap_rel": float(per["gap_rel"].max()),
        "P4_median_seconds": float(per["seconds"].median()),
    }
    passed = {
        "P1": crit["P1_all_valid"],
        "P2": (crit["P2_median_rel_growth_diff"] <= 0.01
               and crit["P2_growth_ratio_max_dev"] <= 1e-3
               and crit["P2_secretion_identical_share"] >= 0.90),
        "P3": crit["P3_max_gap_rel"] <= 0.05 and per["gap_rel"].notna().all(),
        "P4": crit["P4_median_seconds"] <= 300,
    }
    print("\nper sample:")
    cols = ["valid_osqp", "valid_new", "community_growth_osqp", "community_growth_new",
            "max_balance_residual_new", "osqp_status", "max_shift", "gap_rel",
            "n_growing_osqp", "n_growing_new", "median_rel_growth_diff", "seconds"]
    print(per[[c for c in cols if c in per]].round(4).to_string())
    print("\ncriteria:", json.dumps(crit, indent=1))
    print("passed:", passed)

    prov = {"script": "v2/micom/solver_check.py",
            "run_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "solver": highs_qp.describe(), "tradeoff": TRADEOFF, "strategy": "pFBA",
            "samples": list(per.index), "osqp_reference": OSQP_PKL.name,
            "criteria": crit, "passed": passed, "python": platform.python_version(),
            "platform": platform.platform()}
    with open(OUT / "solver_check_provenance.json", "w") as f:
        json.dump(prov, f, indent=1, default=str)


if __name__ == "__main__":
    main()
