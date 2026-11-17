"""CONTRACTS 5d-vi(1): how many Supported links would the tier rule produce by chance?

Within each target metabolite, every non-O class is permuted independently across the pairs
where it is available. Class O, every gate and every availability flag stay where they
are, so the null keeps each class's coverage and support rate and breaks only the
AGREEMENT between classes, which is the thing a tier claims.

  R  its member columns (genomic_producer, genomic_capability), moved as one block
  G  support
  M  support
  X  the Erawijantari (rho, q), with support and contradiction recomputed against the
     row's own Yachida sign (the rule of replication/erawijantari_species.py)

Tiers come from evidence_tools.class_tier itself; an identity permutation must reproduce
every observed Tier_class before any null is drawn.

Run:  python v2/calibration/class_null.py --results "D:/KI Research-dardel/<run>"
Out:  v2/calibration/out/class_null_{summary.tsv,draws.tsv.gz,provenance.json}
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
V2 = HERE.parent
sys.path.insert(0, str(V2))
import evidence_tools as et  # noqa: E402
import utils  # noqa: E402

OUT = HERE / "out"
XTAB = V2 / "replication" / "out" / "erawijantari_species_replication.tsv"
B = 1000
SEED = 20261113
R_COLS = ["Source_genomic_producer_available", "Source_genomic_producer_support",
          "Source_genomic_capability_available", "Source_genomic_capability_support"]
LABELS = ("a_priori_class", "stage_associated")


def tier_word(t: str) -> str:
    return t.split(" (")[0]


def tiers(rows: list[dict]) -> np.ndarray:
    return np.array([tier_word(et.class_tier(r)) for r in rows])


def x_calls(rho_e, q_e, rho_y):
    sig = q_e < utils.FDR_THRESHOLD
    same = np.sign(rho_e) == np.sign(rho_y)
    return sig & same, sig & ~same


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", type=Path, required=True)
    ap.add_argument("--b", type=int, default=B)
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    d = pd.read_csv(a.results / "Results" / "tables" / "final_evidence_classes_v2.csv",
                    low_memory=False)
    d = d[d["Target_label"] != "not_a_target"].reset_index(drop=True)
    x = pd.read_csv(XTAB, sep="\t")[["Species", "Metabolite", "Rho_erawijantari", "Q_erawijantari"]]
    d = d.merge(x, on=["Species", "Metabolite"], how="left", validate="one_to_one")
    xa = d["Source_replication_species_available"].to_numpy(bool)
    assert d.loc[xa, ["Rho_erawijantari", "Q_erawijantari"]].notna().all().all()

    # The X rule, recomputed from (rho, q), must agree with the stored calls.
    sup, con = x_calls(d["Rho_erawijantari"].to_numpy(), d["Q_erawijantari"].to_numpy(),
                       d["Rho_partial"].to_numpy())
    assert (sup[xa] == d.loc[xa, "Source_replication_species_support"].to_numpy(bool)).all()
    assert (con[xa] == d.loc[xa, "Replication_species_contradicted"].to_numpy(bool)).all()

    rows = d.to_dict("records")
    observed = tiers(rows)
    assert (observed == d["Tier_class"].map(tier_word).to_numpy()).all(), "identity check failed"
    labels = d["Target_label"].to_numpy()

    # Permutation groups: within metabolite, among the rows where the class is available.
    r_av = d[["Source_genomic_producer_available",
              "Source_genomic_capability_available"]].any(axis=1).to_numpy()
    avail = {"R": r_av,
             "G": d["Source_gene_level_available"].to_numpy(bool),
             "M": d["Source_mechanistic_available"].to_numpy(bool),
             "X": xa}
    groups = {c: [np.flatnonzero((d["Metabolite"] == m).to_numpy() & av)
                  for m in d["Metabolite"].unique()]
              for c, av in avail.items()}
    groups = {c: [g for g in gs if len(g) > 1] for c, gs in groups.items()}

    r_block = d[R_COLS].to_numpy()
    g_sup = d["Source_gene_level_support"].to_numpy(bool)
    m_sup = d["Source_mechanistic_support"].to_numpy(bool)
    rho_e = d["Rho_erawijantari"].to_numpy()
    q_e = d["Q_erawijantari"].to_numpy()
    rho_y = d["Rho_partial"].to_numpy()

    def perm_index(cls, rng):
        idx = np.arange(len(d))
        for g in groups[cls]:
            idx[g] = rng.permutation(g)
        return idx

    def stats(t):
        good = np.isin(t, ["Supported", "Conclusive"])
        s = {"supported_plus_all": int(good.sum()),
             "conclusive_all": int((t == "Conclusive").sum())}
        for lab in LABELS:
            s[f"supported_plus_{lab}"] = int((good & (labels == lab)).sum())
            s[f"conclusive_{lab}"] = int(((t == "Conclusive") & (labels == lab)).sum())
        return s

    obs = stats(observed)
    rng = np.random.default_rng(SEED)
    draws = []
    for b in range(a.b):
        ir, ig, im, ix = (perm_index(c, rng) for c in ("R", "G", "M", "X"))
        rb, gs, ms = r_block[ir], g_sup[ig], m_sup[im]
        xs, xc = x_calls(rho_e[ix], q_e[ix], rho_y)
        for i, r in enumerate(rows):
            for j, c in enumerate(R_COLS):
                r[c] = rb[i, j]
            r["Source_gene_level_support"] = gs[i]
            r["Source_mechanistic_support"] = ms[i]
            if xa[i]:
                r["Source_replication_species_support"] = bool(xs[i])
                r["Replication_species_contradicted"] = bool(xc[i])
        draws.append(stats(tiers(rows)))
        if (b + 1) % 100 == 0:
            print(f"  {b + 1}/{a.b}  [{time.time() - t0:.0f} s]", flush=True)

    dr = pd.DataFrame(draws)
    dr.to_csv(OUT / "class_null_draws.tsv.gz", sep="\t", index_label="draw")
    summ = []
    for k, o in obs.items():
        v = dr[k].to_numpy()
        summ.append({"statistic": k, "observed": o, "null_mean": round(float(v.mean()), 2),
                     "null_lo95": float(np.percentile(v, 2.5)),
                     "null_hi95": float(np.percentile(v, 97.5)),
                     "p_empirical": round((1 + int((v >= o).sum())) / (1 + a.b), 4),
                     "expected_false_share": round(float(v.mean()) / o, 3) if o else None})
    summ = pd.DataFrame(summ)
    summ.to_csv(OUT / "class_null_summary.tsv", sep="\t", index=False)
    prov = {"script": "v2/calibration/class_null.py", "contract": "CONTRACTS 5d-vi(1)",
            "run_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "results_root": str(a.results), "target_pairs": len(d), "B": a.b, "seed": SEED,
            "permuted_groups": {c: len(g) for c, g in groups.items()},
            "identity_check": "passed", "seconds": round(time.time() - t0)}
    (OUT / "class_null_provenance.json").write_text(json.dumps(prov, indent=1))
    print(summ.to_string(index=False))


if __name__ == "__main__":
    main()
