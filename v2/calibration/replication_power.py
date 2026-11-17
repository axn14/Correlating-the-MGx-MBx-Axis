"""CONTRACTS 5d-vi(2): could Erawijantari (n = 42) have replicated the stage pairs?

Truth is taken from Yachida: rho_true = m * Rho_partial for every X-available pair, with
m in {1.0, 0.75, 0.5} (no, moderate and strong winner's-curse shrinkage; all reported).
Each replicate draws every pair's Fisher z from N(atanh rho_true, 1.06 / (n - k - 3))
(Bonett-Wright variance for Spearman, partial on k covariates), takes p from t with
n - 2 - k df as utils.partial_spearman_grid does, applies BH across the 3,535-pair family,
and calls support on the same sign and q < 0.05. Pairs are drawn independently: the
expected count is unbiased, its spread understated.

Verdict (fixed in 5d-vi before running): the observed stage-associated count among
O-supported X-available pairs at or above the 2.5th percentile of its expectation at
m = 0.5 => the shortfall is power; below it => the stage pairs replicate less than their
own effect sizes predict.

Run:  python v2/calibration/replication_power.py --results "D:/KI Research-dardel/<run>"
Out:  v2/calibration/out/replication_power_{summary,per_metabolite}.tsv, _provenance.json
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

HERE = Path(__file__).resolve().parent
V2 = HERE.parent
sys.path.insert(0, str(V2))
import utils  # noqa: E402

OUT = HERE / "out"
XDIR = V2 / "replication" / "out"
B = 1000
SEED = 20261113
SHRINK = (1.0, 0.75, 0.5)
LABELS = ("a_priori_class", "stage_associated")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", type=Path, required=True)
    ap.add_argument("--b", type=int, default=B)
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    prov_x = json.loads((XDIR / "erawijantari_provenance.json").read_text())
    n, k = prov_x["n_samples"], len(prov_x["covariates"])
    dof = n - 2 - k
    se = np.sqrt(1.06 / (n - k - 3))

    x = pd.read_csv(XDIR / "erawijantari_species_replication.tsv", sep="\t")
    x = x[x["Source_replication_species_available"]].reset_index(drop=True)
    assert len(x) == prov_x["measurable_pairs"]
    cls = pd.read_csv(a.results / "Results" / "tables" / "final_evidence_classes_v2.csv",
                      usecols=["Species", "Metabolite", "Target_label", "Classes_supporting"],
                      low_memory=False)
    x = x.merge(cls, on=["Species", "Metabolite"], how="left", validate="one_to_one")
    assert x["Target_label"].isin(LABELS).all()
    x["O_supported"] = x["Classes_supporting"].fillna("").str.contains("O")
    focus = x["O_supported"].to_numpy()
    lab = x["Target_label"].to_numpy()
    rho_y = x["Rho_partial"].to_numpy()
    observed = x["Source_replication_species_support"].to_numpy(bool)

    rng = np.random.default_rng(SEED)
    rows, per_met = [], []
    for m in SHRINK:
        mu = np.arctanh(np.clip(m * rho_y, -0.999, 0.999))
        hits = np.zeros((a.b, len(x)), bool)
        for b in range(a.b):
            r = np.tanh(rng.normal(mu, se))
            t = r * np.sqrt(dof / np.clip(1 - r ** 2, 1e-12, None))
            p = 2 * stats.t.sf(np.abs(t), dof)
            q = multipletests(p, method="fdr_bh")[1]
            hits[b] = (q < utils.FDR_THRESHOLD) & (np.sign(r) == np.sign(rho_y))
        for name, sel in [("all_O_supported", focus)] + [(l, focus & (lab == l)) for l in LABELS]:
            cnt = hits[:, sel].sum(axis=1)
            o = int(observed[sel].sum())
            rows.append({"shrinkage_m": m, "set": name, "pairs": int(sel.sum()),
                         "observed_replicated": o,
                         "expected_mean": round(float(cnt.mean()), 1),
                         "expected_lo95": float(np.percentile(cnt, 2.5)),
                         "expected_hi95": float(np.percentile(cnt, 97.5)),
                         "mean_power_per_pair": round(float(hits[:, sel].mean()), 3),
                         "p_obs_or_fewer": round((1 + int((cnt <= o).sum())) / (1 + a.b), 4)})
        pw = hits.mean(axis=0)
        g = pd.DataFrame({"Metabolite": x["Metabolite"], "Target_label": lab, "focus": focus,
                          "power": pw, "observed": observed})
        g = g[g["focus"]].groupby(["Metabolite", "Target_label"]).agg(
            pairs=("power", "size"), expected=("power", "sum"),
            observed=("observed", "sum")).reset_index()
        g["shrinkage_m"] = m
        per_met.append(g)

    summ = pd.DataFrame(rows)
    summ.to_csv(OUT / "replication_power_summary.tsv", sep="\t", index=False)
    pm = pd.concat(per_met)
    pm["expected"] = pm["expected"].round(2)
    pm.to_csv(OUT / "replication_power_per_metabolite.tsv", sep="\t", index=False)

    st = summ[(summ["shrinkage_m"] == 0.5) & (summ["set"] == "stage_associated")].iloc[0]
    verdict = ("power" if st["observed_replicated"] >= st["expected_lo95"]
               else "replicates less than its own effect sizes predict")
    prov = {"script": "v2/calibration/replication_power.py", "contract": "CONTRACTS 5d-vi(2)",
            "run_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "results_root": str(a.results), "n": n, "k_covariates": k, "dof": dof,
            "fisher_se": round(float(se), 4), "family": int(len(x)), "B": a.b, "seed": SEED,
            "shrinkage": list(SHRINK), "stage_verdict_at_m0.5": verdict}
    (OUT / "replication_power_provenance.json").write_text(json.dumps(prov, indent=1))
    print(summ.to_string(index=False))
    print("\nstage verdict at m = 0.5:", verdict)


if __name__ == "__main__":
    main()
