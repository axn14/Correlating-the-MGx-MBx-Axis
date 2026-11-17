"""Would losing ~40% of samples to solver failure bias the progression analysis?

WHY THIS SCRIPT EXISTS. CONTRACTS section 8 records that 5 of 12 samples end infeasible
at the pFBA step under HiGHS + OSQP, and `diet_screen.py --community 8` established that
the failure is numerical rather than dietary or biological (the same sample fails under
both media, with every taxon at zero growth and a lumen imbalance seven orders of
magnitude over tolerance). That leaves a choice: buy a QP that solves them, or run the
batch and drop the failures.

Dropping them is only safe if the dropout is INDEPENDENT of stage. If it is not -- if
whatever makes a community model hard to solve also tracks disease stage -- then the
surviving 60% is a stage-biased subsample and every progression comparison computed on it
is confounded by solver behaviour. That is not a small methodological caveat; it is the
kind of error that invalidates the result while leaving every number looking reasonable.

WHAT THIS CAN AND CANNOT ESTABLISH. It cannot test the dropout directly: the 2026-10-23
pilot recorded only a count, not which 5 of the 12 failed, and the 8-sample run has a
single failure, which supports no correlation at all. So this script does the two things
that ARE possible without a working QP:

  1. It measures whether the plausible DRIVERS of solver difficulty -- how many taxa a
     sample's community model contains, how large that model is, how even the abundances
     are -- track the stage axis across all 277 primary samples. If they do not track
     stage, then a dropout linked to them cannot be stage-biased, and dropping failures
     is defensible. If they do, the risk is real and conditional on the link.
  2. It quantifies the consequence under an explicit worst case: if the 42% that fail
     were exactly the hardest 42% by each candidate driver, what would the retained
     cohort's stage composition be?

Point 2 is a conditional, and is reported as one. The premise -- that failure is driven
by these properties -- is UNTESTED here. What makes the conditional worth computing is
that its antecedent cannot be ruled out, and the cost of being wrong is asymmetric: a
free academic licence against a confounded result.

Reads the same inputs as diet_screen.py and writes out/dropout_risk.tsv plus a provenance
JSON. No HPC, no solver, about a minute.
"""
from __future__ import annotations

import json
import pickle
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

HERE = Path(__file__).resolve().parent
V2 = HERE.parent
# Both explicitly: `python -I` implies -P, which drops the script's own directory from
# sys.path, so `import diet_screen` fails under the isolated interpreter without this.
sys.path.insert(0, str(V2))
sys.path.insert(0, str(HERE))

import utils  # noqa: E402
from diet_screen import (AGORA_EXTRACTED, OUT, PREPROCESSED, YACHIDA,  # noqa: E402
                         agora_index, cohort_species)

# ONE axis, not two. `cohort_species()` applies `utils.primary_mask`, which from v2.2
# excludes PostSurgery AND Adenoma, so these 277 samples are exactly Healthy + Early_CRC
# + Advanced_CRC -- and they are exactly the set NB10's batch solves. Running the
# four-level early axis here would re-test the same 277 rows with different numeric
# spacing and no Adenoma rows at all, which is not a second family, just a second
# p-value for the same comparison. The adenoma arm is outside the batch, so the dropout
# question does not reach it.
PRIMARY_ORDER = utils.STAGE_GROUP_ORDER


def shannon(row: np.ndarray) -> float:
    p = row[row > 0]
    if p.size == 0:
        return 0.0
    p = p / p.sum()
    return float(-(p * np.log(p)).sum())


def model_sizes() -> pd.DataFrame | None:
    """Per-taxon reaction and metabolite counts from the 8 solved community models.

    growth_rates carries `reactions` and `metabolites` per taxon, so summing them gives
    the size of the LP/QP that actually went to the solver -- the quantity that stresses
    it, which taxon count only proxies. Returns None when --community has not been run.
    """
    p = OUT / "community" / "grow_gut_diet_v2.pkl"
    if not p.exists():
        return None
    with open(p, "rb") as f:
        res = pickle.load(f)
    gr = res.growth_rates.copy()
    gr["sample_id"] = gr["sample_id"].astype(str)
    out = gr.groupby("sample_id").agg(
        solved_n_taxa=("taxon", "nunique"),
        solved_reactions=("reactions", "sum"),
        solved_metabolites=("metabolites", "sum"),
        solved_n_growing=("growth_rate", lambda s: int((s > 1e-6).sum())))

    # Community growth and target-compound production, so the question "does what MICOM
    # REPORTS scale with taxon count?" can be asked. If it does, a richness-linked
    # dropout biases the reported fluxes themselves, not merely the group sizes.
    ex = res.exchanges.copy()
    ex["sample_id"] = ex["sample_id"].astype(str)
    import micom_tools as _mt
    q = _mt.solution_quality(ex, res.growth_rates)
    q["sample_id"] = q["sample_id"].astype(str)
    out = out.join(q.set_index("sample_id")[["community_growth", "valid"]])
    for b in ("ptrc", "spmd"):
        m = (ex["reaction"].str.contains(f"{b}_m|{b}\\(e\\)|{b}\\[e\\]", regex=True)
             & (ex["taxon"] != "medium") & (ex["flux"] > 0))
        out[f"solved_{b}_production"] = (ex.loc[m].groupby("sample_id")["flux"].sum()
                                         .reindex(out.index).fillna(0.0))
    return out


def per_sample_table() -> pd.DataFrame:
    """One row per primary sample: what its community model would contain, plus stage."""
    _, abund = cohort_species()
    idx = agora_index()
    modelled = set(idx["species_plain"])

    cols_model = [c for c in abund.columns if c in modelled]
    a_all = abund.to_numpy()
    a_mod = abund[cols_model].to_numpy()

    t = pd.DataFrame(index=abund.index)
    # n_taxa is the count that enters `micom.workflows.build`: a species above the 0.1%
    # cutoff AND present in AGORA2. Taxa failing either are not in the model at all, so
    # counting them would describe a community nobody solves.
    t["n_taxa_modelled"] = (a_mod > 0).sum(axis=1)
    t["n_taxa_detected"] = (a_all > 0).sum(axis=1)
    t["agora2_coverage"] = a_mod.sum(axis=1) / np.where(a_all.sum(axis=1) > 0,
                                                        a_all.sum(axis=1), np.nan)
    t["shannon_modelled"] = [shannon(r) for r in a_mod]
    # Evenness is the conditioning proxy: a community of 60 taxa at comparable abundance
    # couples 60 growth variables through the lumen at comparable scale, where one
    # dominant taxon plus a tail does not.
    t["evenness_modelled"] = t["shannon_modelled"] / np.log(
        t["n_taxa_modelled"].where(t["n_taxa_modelled"] > 1))
    t["max_abundance_modelled"] = a_mod.max(axis=1)

    with open(PREPROCESSED, "rb") as f:
        pdata = pickle.load(f)
    meta = pdata["harmonized_meta"][YACHIDA].copy()
    if "Sample" not in meta.columns:
        meta.reset_index(inplace=True)
        meta.rename(columns={meta.columns[0]: "Sample"}, inplace=True)
    meta["Sample"] = meta["Sample"].astype(str)
    meta = meta.set_index("Sample")

    t = t.join(meta[["Stage.Group"]], how="left")
    t["stage_numeric"] = t["Stage.Group"].map(utils.STAGE_GROUP_NUMERIC)

    sizes = model_sizes()
    if sizes is not None:
        t = t.join(sizes, how="left")
    return t


DRIVERS = ["n_taxa_modelled", "shannon_modelled", "evenness_modelled",
           "agora2_coverage", "max_abundance_modelled"]


def stage_association(t: pd.DataFrame) -> pd.DataFrame:
    """Does each candidate driver track stage? Unadjusted and confounder-adjusted.

    One BH family over the five drivers on the three-level primary axis. A driver that
    does not track stage cannot carry a stage bias into the dropout, so these five
    p-values are the whole test and correcting them together is the whole correction.
    """
    with open(PREPROCESSED, "rb") as f:
        pdata = pickle.load(f)
    meta = pdata["harmonized_meta"][YACHIDA].copy()
    if "Sample" not in meta.columns:
        meta.reset_index(inplace=True)
        meta.rename(columns={meta.columns[0]: "Sample"}, inplace=True)
    meta["Sample"] = meta["Sample"].astype(str)
    # build_confounder_matrix wants "Sample" as a COLUMN and returns it as one; the
    # partial-correlation call needs it as the index to align with t.
    conf = (utils.build_confounder_matrix(meta)
            .set_index("Sample")
            .apply(pd.to_numeric, errors="coerce"))

    rows = []
    col = "stage_numeric"
    sub = t[t[col].notna()]
    for d in DRIVERS:
        s = sub[[d, col]].dropna()
        if len(s) < 10:
            continue
        rho, p = stats.spearmanr(s[d], s[col])
        # Kruskal across the named groups, which does not assume the axis is linear --
        # a U-shaped driver would track stage while giving rho near zero.
        groups = [g[d].dropna().to_numpy()
                  for _, g in sub.groupby("Stage.Group") if len(g) >= 3]
        kw_p = stats.kruskal(*groups).pvalue if len(groups) >= 2 else np.nan

        # Confounder-adjusted, on complete cases only (partial_spearman_grid requires
        # them). Age/BMI/Brinkman/alcohol/sex is the matrix NB02 uses and is 100%
        # complete for Yachida, so this should not lose samples.
        common = s.index.intersection(conf.dropna().index)
        if len(common) >= 20:
            # Returns (r, p) as (p, q) arrays, so a 1-vs-1 call is the [0, 0] cell.
            a_r, a_pv = utils.partial_spearman_grid(
                t.loc[common, [d]], t.loc[common, [col]], conf.loc[common])
            a_rho, a_p = float(a_r[0, 0]), float(a_pv[0, 0])
        else:
            a_rho = a_p = np.nan
        rows.append({"driver": d, "n": len(s),
                     "spearman_rho": rho, "p": p,
                     "kruskal_p": kw_p,
                     "adj_rho": a_rho, "adj_p": a_p, "n_adjusted": len(common),
                     **{f"median_{g}": float(sub.loc[sub['Stage.Group'] == g, d]
                                             .median())
                        for g in PRIMARY_ORDER if (sub["Stage.Group"] == g).any()}})
    out = pd.DataFrame(rows)
    out["q"] = multipletests(out["p"], method="fdr_bh")[1]
    out["adj_q"] = multipletests(out["adj_p"], method="fdr_bh")[1]
    return out


def worst_case_dropout(t: pd.DataFrame, fail_rate: float = 5 / 12) -> pd.DataFrame:
    """If the samples that fail were exactly the hardest by each driver, who is left?

    THIS IS A CONDITIONAL, not a measurement. Its premise -- that solver failure is
    monotone in the driver -- is untested, and deliberately extreme: it is the worst the
    bias can be if the link exists at all, which is the number a methods decision needs.
    The null row is the same dropout taken at random, which is what "safe to drop" means.
    """
    n_drop = int(round(fail_rate * len(t)))
    base = (t["Stage.Group"].value_counts(normalize=True) * 100)
    groups = list(base.index)
    rows = [{"scenario": "all 277 samples", "driver": "-", "n_retained": len(t),
             **{f"pct_{k}": round(base[k], 1) for k in groups},
             "healthy_vs_adv_taxa_gap": round(
                 t.loc[t["Stage.Group"] == "Advanced_CRC", "n_taxa_modelled"].median()
                 - t.loc[t["Stage.Group"] == "Healthy", "n_taxa_modelled"].median(), 2)}]

    # THE NULL IS A DISTRIBUTION, NOT ONE DRAW. A single seeded permutation of this
    # cohort moves pct_Healthy by about 5 points on its own, which is the same order as
    # the driver-led scenarios below -- so comparing them to one random draw would
    # invent a difference or hide one. 2,000 draws give the band that a driver-led shift
    # has to clear before it means anything.
    rng = np.random.default_rng(20261108)
    sim = np.empty((2000, len(groups)))
    gap = np.empty(2000)
    codes = t["Stage.Group"].to_numpy()
    taxa = t["n_taxa_modelled"].to_numpy()
    for i in range(2000):
        keep = rng.permutation(len(t))[n_drop:]
        g = codes[keep]
        sim[i] = [(g == k).mean() * 100 for k in groups]
        gap[i] = (np.median(taxa[keep][g == "Advanced_CRC"])
                  - np.median(taxa[keep][g == "Healthy"]))
    rows.append({"scenario": f"random {n_drop} dropped: null mean", "driver": "-",
                 "n_retained": len(t) - n_drop,
                 **{f"pct_{k}": round(sim[:, j].mean(), 1)
                    for j, k in enumerate(groups)},
                 "healthy_vs_adv_taxa_gap": round(float(gap.mean()), 2)})
    rows.append({"scenario": f"random {n_drop} dropped: null 95% band", "driver": "-",
                 "n_retained": len(t) - n_drop,
                 **{f"pct_{k}": (f"{np.percentile(sim[:, j], 2.5):.1f} to "
                                 f"{np.percentile(sim[:, j], 97.5):.1f}")
                    for j, k in enumerate(groups)},
                 "healthy_vs_adv_taxa_gap": (f"{np.percentile(gap, 2.5):.2f} to "
                                             f"{np.percentile(gap, 97.5):.2f}")})

    for d in DRIVERS:
        for direction in ("highest", "lowest"):
            s = t[d].dropna().sort_values(ascending=(direction == "highest"))
            keep = s.index[n_drop:]
            k_g = t.loc[keep, "Stage.Group"]
            v = (k_g.value_counts(normalize=True) * 100)
            kt = t.loc[keep, "n_taxa_modelled"]
            rows.append({"scenario": f"{direction} {d} dropped", "driver": d,
                         "n_retained": len(keep),
                         **{f"pct_{k}": round(v.get(k, 0.0), 1) for k in groups},
                         **{f"n_{k}": int((k_g == k).sum()) for k in groups},
                         "healthy_vs_adv_taxa_gap": round(
                             kt[k_g == "Advanced_CRC"].median()
                             - kt[k_g == "Healthy"].median(), 2)})
    out = pd.DataFrame(rows)
    # The retained per-group n is the argument that survives: the composition barely
    # moves, but 115 samples is 115 samples.
    for k in groups:
        out.loc[out["scenario"] == "all 277 samples", f"n_{k}"] = int(
            (t["Stage.Group"] == k).sum())
    return out


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    t = per_sample_table()
    print(f"per-sample table: {len(t)} primary samples, "
          f"{t['Stage.Group'].isna().sum()} without a Stage.Group")
    print(t.groupby("Stage.Group")[["n_taxa_modelled", "shannon_modelled",
                                    "agora2_coverage"]].median().to_string())

    assoc = stage_association(t)
    print("\n=== does any candidate driver of solver difficulty track stage?")
    show = ["driver", "n", "spearman_rho", "p", "q", "kruskal_p",
            "adj_rho", "adj_q"]
    print(assoc[show].to_string(index=False, float_format=lambda x: f"{x:.4g}"))

    drop = worst_case_dropout(t)
    print("\n=== stage composition if the hardest 42% were the ones that fail "
          "(CONDITIONAL -- premise untested)")
    print(drop.to_string(index=False))

    sizes = [c for c in ("solved_n_taxa", "solved_reactions", "solved_metabolites")
             if c in t.columns]
    if sizes:
        s = t[t["solved_n_taxa"].notna()]
        print(f"\n=== the 8 solved samples: model size, and where the FAILING one sits")
        print(s[["n_taxa_modelled", "solved_n_taxa", "solved_reactions",
                 "solved_metabolites", "shannon_modelled", "Stage.Group"]]
              .to_string())
        if "10091" in s.index:
            for c in ("solved_reactions", "solved_n_taxa", "shannon_modelled",
                      "evenness_modelled"):
                pct = float((s[c] < s.loc["10091", c]).mean() * 100)
                print(f"         10091 {c}: {s.loc['10091', c]:.4g}, "
                      f"larger than {pct:.0f}% of the 8 solved")

        # Does what MICOM REPORTS scale with taxon count? If it does, a richness-linked
        # dropout biases the fluxes themselves and not just the group sizes. n=7, so
        # this is an indication and is reported with its n, never as an estimate.
        ok = s[s["valid"].eq(True)]      # .eq avoids fillna on an object column
        print(f"\n=== does the REPORTED quantity scale with taxon count? (n={len(ok)} "
              f"valid samples -- indicative only)")
        for c in ("community_growth", "solved_ptrc_production",
                  "solved_spmd_production"):
            if c in ok.columns and ok[c].notna().sum() >= 4:
                rho, p = stats.spearmanr(ok["n_taxa_modelled"], ok[c])
                print(f"         n_taxa vs {c:24s} rho {rho:+.3f}  p {p:.3f}")

    t.to_csv(OUT / "dropout_risk.tsv", sep="\t")
    assoc.to_csv(OUT / "dropout_risk_association.tsv", sep="\t", index=False)
    drop.to_csv(OUT / "dropout_risk_scenarios.tsv", sep="\t", index=False)
    prov = {
        "built_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "n_samples": int(len(t)),
        "fail_rate_assumed": 5 / 12,
        "note": ("The dropout itself is NOT measured here: the 2026-10-23 pilot recorded "
                 "a count and not which samples failed, and the 8-sample community run "
                 "has one failure. This measures whether the plausible drivers of solver "
                 "difficulty track stage, and what the bias would be if they drove it."),
        "association": json.loads(assoc.to_json(orient="records")),
        "scenarios": json.loads(drop.to_json(orient="records")),
        "inputs": {"preprocessed": str(PREPROCESSED),
                   "agora2": str(AGORA_EXTRACTED)},
        "software": {"python": sys.version.split()[0], "pandas": pd.__version__},
    }
    (OUT / "dropout_risk_provenance.json").write_text(
        json.dumps(prov, indent=2), encoding="utf-8")
    print(f"\nwrote {OUT / 'dropout_risk.tsv'} and two summary tables")


if __name__ == "__main__":
    main()
