#!/usr/bin/env python3
"""Why the MICOM batch is on hold, answered at the level where it can be answered.

THE TWO BLOCKERS (CONTRACTS.md section 8, measured 2026-10-23 on 12 samples)
---------------------------------------------------------------------------
1. 5 of 12 samples end INFEASIBLE at the pFBA step under the open-source HiGHS+OSQP
   hybrid. Looser tolerances did not help; the same samples solve with
   strategy="none". Read at the time as a solver problem.
2. `western_diet_gut.csv` supplies putrescine and spermidine at 1.0 mmol/gDW/h, so no
   taxon secretes either and the batch would produce no species-level evidence for
   them. Removing both from the diet made 1 of 4 test samples show secretion of both
   -- and 3 of 4 become infeasible.

Those two observations are usually read as separate, and the second is the clue to the
first. A community model is infeasible when its constraints cannot all be met at once,
and at tradeoff 0.5 micom REQUIRES every taxon to reach half its own maximum growth
rate. One taxon that cannot grow at all on the medium therefore makes the whole sample
infeasible, with no error and nothing naming the taxon. That is a medium question, and
it can be settled on single-species models in minutes instead of on community models in
hours.

WHAT THIS SCRIPT DOES
---------------------
    stage 1   the species the cohort actually contains, matched to AGORA2 by the same
              GTDB-native rule NB10 uses (no suffix stripping, no strain collapsing)
    stage 2   every matched model, grown on the diet as written and on the diet with
              each target compound withdrawn. A model that grows on the full diet and
              not without a compound is AUXOTROPHIC for it, and is the mechanism by
              which withdrawing that compound makes a community infeasible
    stage 3   for each auxotroph, the minimum uptake flux that restores growth,
              bisected. This is the number the reformulated diet needs
    stage 4   the reformulated medium, written out with the reasoning in its header

WHAT IT CANNOT SETTLE
---------------------
Whether the REMAINING infeasibility is the solver. Single-model LPs are not the
community QP, and a sample can be infeasible for reasons no single model shows (the
tradeoff constraint couples taxa through a shared lumen). This script bounds the
medium's contribution; whatever is left after a feasible medium is the solver's, and
`--community` measures that directly on real samples.

Run:
    python v2/micom/diet_screen.py                 # stages 1-4, ~15 min
    python v2/micom/diet_screen.py --limit 50      # quick pass over 50 models
    python v2/micom/diet_screen.py --community 6   # also re-solve 6 samples, old vs new
"""

from __future__ import annotations

import argparse
import json
import pickle
import platform
import re
import sys
import time
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent          # v2/micom
V2 = HERE.parent                                 # v2
ROOT = V2.parent                                 # repo root
sys.path.insert(0, str(V2))

import utils        # noqa: E402
import micom_tools as mt  # noqa: E402

OUT = HERE / "out"
DATA = ROOT / "Data"
MICOM_DIR = DATA / "micom"
AGORA_DB = MICOM_DIR / "agora2_gtdb_species.qza"
AGORA_EXTRACTED = MICOM_DIR / "agora2_extracted"
DIET_IN = MICOM_DIR / "western_diet_gut.csv"
DIET_OUT = HERE / "gut_diet_v2.csv"
PREPROCESSED = ROOT / "Results" / "intermediate" / "preprocessed_data_v2.pkl"

YACHIDA = utils.DATASET_DISCOVERY
ABUND_CUTOFF = 1e-3          # NB10's cutoff: taxa below 0.1% never enter a community

# What the chapter wants to OBSERVE being secreted. A compound the diet supplies cannot
# be observed as a product, because every community takes it up instead. These are the
# base VMH ids; the diet file writes them as EX_<id>_e and single models as EX_<id>(e).
#
# Only the first two are in the diet at all. Everything below the rule is already
# absent, which is the part of the 2026-10-23 note that turned out to be wrong: the
# medium does NOT need rebuilding around tryptamine, histamine or N-acetylglutamate,
# because it never supplied them. The decision is confined to putrescine and spermidine.
TARGETS_IN_DIET = ["ptrc", "spmd"]
TARGETS_NOT_IN_DIET = ["agm", "15dap", "sprm", "4abut", "trypta", "hista", "tym",
                       "acglu"]
# Substrates, which must STAY. Withdrawing the substrate of the reaction you want to
# observe does not make the product easier to see, it makes it impossible.
SUBSTRATES = ["orn", "arg_L", "arg__L", "lys_L", "lys__L", "his_L", "his__L",
              "trp_L", "trp__L", "tyr_L", "tyr__L", "glu_L", "glu__L"]

GROWTH_FLOOR = mt.GROWTH_FLOOR   # 1e-6 /h, micom's own "this taxon grows"
SAFETY = 2.0                     # the reformulated supply is this times the measured need


# ===========================================================================
# stage 1 -- the species the cohort contains
# ===========================================================================

def gtdb_binomial(col: str):
    """Species binomial from a column name, keeping GTDB suffixes.

    Copied in behaviour from NB10 cell 3 on purpose. The screen has to match on the
    same names the build matches on, or it would characterise models the communities
    never contain.
    """
    s = str(col)
    if ";s__" in s or s.startswith("s__"):
        m = re.search(r"s__([^;]+)", s)
        name = m.group(1).strip() if m else ""
    else:
        name = s.strip()
    toks = name.split()
    return name if len(toks) >= 2 else None


def cohort_species():
    """Species above NB10's abundance cutoff in the primary samples.

    Returns (per-species summary, per-sample abundance matrix). The matrix is what the
    community-demand calculation needs: a taxon's draw on the lumen is its abundance
    times its per-gram flux, so a per-species average cannot stand in for it.
    """
    with open(PREPROCESSED, "rb") as f:
        pdata = pickle.load(f)
    filt = pdata["filtered_species"][YACHIDA].copy()
    meta = pdata["harmonized_meta"][YACHIDA].copy()
    for df in (filt, meta):
        if "Sample" not in df.columns:
            df.reset_index(inplace=True)
            df.rename(columns={df.columns[0]: "Sample"}, inplace=True)
        df["Sample"] = df["Sample"].astype(str)

    prim = utils.primary_mask(meta.set_index("Sample"))
    keep = meta.set_index("Sample").index[prim.values]
    filt = filt[filt["Sample"].isin(keep)].set_index("Sample")
    cols = list(filt.columns)
    filt = filt[cols].div(filt[cols].sum(axis=1), axis=0).fillna(0.0)

    long = filt.reset_index().melt(id_vars="Sample", var_name="_raw",
                                   value_name="abundance")
    long = long[long["abundance"] > ABUND_CUTOFF].copy()
    long["species"] = long["_raw"].map(gtdb_binomial)
    long = long.dropna(subset=["species"])

    mat = (long.pivot_table(index="Sample", columns="species", values="abundance",
                            aggfunc="sum", fill_value=0.0)
           .reindex(filt.index, fill_value=0.0))
    agg = (long.groupby("species")
           .agg(samples=("Sample", "nunique"),
                mean_abundance=("abundance", "mean"),
                max_abundance=("abundance", "max"))
           .reset_index())
    agg["prevalence"] = agg["samples"] / filt.shape[0]
    agg.attrs["n_samples"] = int(filt.shape[0])
    return agg.sort_values("mean_abundance", ascending=False), mat


def agora_index() -> pd.DataFrame:
    """species -> model json path, from the pre-extracted AGORA2 store."""
    man = next(AGORA_EXTRACTED.rglob("manifest.csv"), None)
    if man is None:
        with ZipFile(str(AGORA_DB)) as z:
            raise FileNotFoundError(
                f"{AGORA_EXTRACTED} holds no manifest.csv; NB10 cell 5 extracts the "
                f".qza once ({len(z.namelist())} entries in the archive).")
    df = pd.read_csv(man)
    df["species_plain"] = df["species"].astype(str).str.replace(r"^s__", "",
                                                                regex=True)
    df["path"] = df["file"].map(lambda f: str(Path(man).parent / f))
    return df


# ===========================================================================
# stage 2 -- the auxotrophy screen
# ===========================================================================

def read_diet() -> pd.DataFrame:
    """The diet as written, plus the base id and the single-model reaction name.

    The file uses `EX_x_e`; NB10 rewrites that to `EX_x_m` for the community lumen, and
    single AGORA2 models spell the same reaction `EX_x(e)`. All three are the same
    exchange and the screen has to speak the single-model dialect.
    """
    d = pd.read_csv(DIET_IN)
    if not {"reaction", "flux"} <= set(d.columns):
        raise ValueError(f"{DIET_IN.name} needs reaction,flux columns")
    dup = d[d["reaction"].duplicated(keep=False)]
    if dup.groupby("reaction")["flux"].nunique().gt(1).any():
        raise ValueError(f"{DIET_IN.name} lists a reaction twice with different "
                         f"fluxes:\n{dup}")
    d = d.drop_duplicates("reaction").copy()
    d["base"] = d["reaction"].str.replace(r"^EX_", "", regex=True) \
                             .str.replace(r"_e$", "", regex=True)
    d["model_rxn"] = "EX_" + d["base"] + "(e)"
    return d


def exchange_census(paths) -> pd.Series:
    """How many of the cohort's models carry each exchange reaction.

    One pass of plain JSON parsing over the model files, no cobra, because all this
    needs is the reaction ids. Used to check that the diet's reaction names are names
    the models actually use.
    """
    import collections
    c = collections.Counter()
    for p in paths:
        with open(p, encoding="utf-8") as f:
            m = json.load(f)
        for r in m["reactions"]:
            rid = r["id"]
            if rid.startswith("EX_"):
                c[rid] += 1
    return pd.Series(c, dtype=int).sort_values(ascending=False)


def audit_diet(diet: pd.DataFrame, census: pd.Series, n_models: int):
    """Which diet lines reach the models, and which are names the models never use.

    A medium is a dictionary lookup. A reaction id the models do not have is not an
    error, a warning or a zero-flux entry -- it is simply absent from every model's
    medium, and the nutrient is silently not supplied. `EX_glc_e` is exactly that case
    here: the AGORA2 models all spell glucose `EX_glc_D(e)`, so the single most
    important carbon source in the file has been reaching none of them.

    An unmatched id is renamed ONLY when exactly one stereoisomer-suffixed form exists
    in the models, which is unambiguous. Anything else is reported unmatched rather
    than guessed at: the seven remaining lines are long-chain fatty acids, AMP,
    adenosylcobalamin and Fe3-dicitrate, which no cohort model exchanges at all.
    """
    d = diet.copy()
    d["n_models"] = d["model_rxn"].map(lambda r: int(census.get(r, 0)))
    rename = {}
    for r in d.loc[d["n_models"] == 0, "model_rxn"]:
        base = r[3:-3]
        cands = [c for c in (f"EX_{base}_D(e)", f"EX_{base}_L(e)")
                 if int(census.get(c, 0)) > 0]
        if len(cands) == 1:
            rename[r] = cands[0]
    d["renamed_to"] = d["model_rxn"].map(rename).fillna("")
    d.loc[d["renamed_to"] != "", "model_rxn"] = d.loc[d["renamed_to"] != "",
                                                      "renamed_to"]
    d["n_models"] = d["model_rxn"].map(lambda r: int(census.get(r, 0)))
    unmatched = d.loc[d["n_models"] == 0, "reaction"].tolist()
    return d, rename, unmatched


def screen_model(path: str, diet: dict[str, float], withdraw: dict[str, list[str]]):
    """Growth on the full diet and on each withdrawal of it, for one model.

    `withdraw` maps a label to the reaction ids to remove. Returns a dict of growth
    rates; np.nan where the model has no exchange for a withdrawn compound (nothing was
    actually withdrawn, so the number would be the full-diet one and would read as
    "not auxotrophic" for a model that simply cannot transport the compound).
    """
    from cobra.io import load_json_model

    model = load_json_model(path)
    have = {r.id for r in model.exchanges}
    medium = {r: f for r, f in diet.items() if r in have}
    if not medium:
        return {"error": "no diet exchange is in this model"}

    out = {"n_exchanges": len(have), "n_medium": len(medium)}
    model.medium = medium
    out["growth_full"] = float(model.slim_optimize() or 0.0)

    for label, rxns in withdraw.items():
        present = [r for r in rxns if r in medium]
        out[f"has_{label}"] = bool([r for r in rxns if r in have])
        if not present:
            out[f"growth_without_{label}"] = np.nan
            continue
        sub = {k: v for k, v in medium.items() if k not in present}
        model.medium = sub
        out[f"growth_without_{label}"] = float(model.slim_optimize() or 0.0)
        model.medium = medium
    return out


SATURATION = 0.99   # "fully supplied" = this fraction of the growth the full diet gives


def saturating_supply(path: str, diet: dict[str, float], rxn: str,
                      hi: float = 1.0, tol: float = 1e-5):
    """Smallest uptake flux of `rxn` at which this model reaches SATURATION of its
    growth on the full diet -- and the dose-response that justifies the number.

    WHY NOT "THE SMALLEST FLUX AT WHICH IT GROWS AT ALL".
    That was the first version of this function and it measured the wrong quantity.
    Spermidine is a biomass precursor in these models, so the uptake a taxon needs is
    PROPORTIONAL to the growth rate it achieves: Bacteroides uniformis grows at 0.0014,
    0.0142 and 0.1424 /h on 1e-5, 1e-4 and 1e-3 mmol/gDW/h, a straight line, and only
    levels off at 0.859 /h somewhere below 0.01. A bisection for "grows at all" returns
    6e-05 and a diet built on it would hold every auxotroph at about 1% of its growth
    rate -- feasible, and meaningless. The quantity that matters is where the curve
    flattens.

    Bisection between 0 and the diet's own value on the predicate "reaches SATURATION",
    which is monotone in the supply, plus a coarse dose-response grid kept for the
    record so the saturation claim can be checked rather than believed.
    """
    from cobra.io import load_json_model

    model = load_json_model(path)
    have = {r.id for r in model.exchanges}
    base = {r: f for r, f in diet.items() if r in have}
    if rxn not in base:
        return np.nan, {}

    def growth(level: float) -> float:
        m = dict(base)
        m[rxn] = level
        model.medium = m
        return float(model.slim_optimize() or 0.0)

    gmax = growth(base[rxn])
    curve = {f"{L:g}": round(growth(L), 6)
             for L in (0.0, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, base[rxn])}
    if gmax <= GROWTH_FLOOR:
        return np.nan, curve
    target = SATURATION * gmax
    if growth(0.0) >= target:
        return 0.0, curve

    lo, hi = 0.0, float(base[rxn])
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        if growth(mid) >= target:
            hi = mid
        else:
            lo = mid
    return hi, curve


def complete_for(path: str, diet: dict[str, float], growth: float = 0.01,
                 max_added: float = 1.0):
    """The smallest set of extra exchanges that lets a starved model grow.

    `micom.media.complete_medium` solves for the minimum-weight addition to a medium
    that supports a target growth rate. A species that cannot grow on the western diet
    at all makes every community containing it infeasible at any nonzero tradeoff, with
    no error and nothing naming the species, so what it is missing is worth knowing by
    name rather than by elimination.

    Returns (additions as {reaction: flux}, error string or "").
    """
    from cobra.io import load_json_model
    from micom.media import complete_medium

    try:
        model = load_json_model(path)
        have = {r.id for r in model.exchanges}
        base = {r: f for r, f in diet.items() if r in have}
        full = complete_medium(model, pd.Series(base), growth=growth,
                               max_import=max_added, minimize_components=True)
        extra = {k: float(v) for k, v in full.items()
                 if k not in base and float(v) > 1e-9}
        return extra, ""
    except Exception as e:
        return {}, f"{type(e).__name__}: {e}"


def community_demand(sat: pd.Series, abundance: pd.DataFrame) -> pd.Series:
    """Per-sample community-level uptake needed, as sum(abundance x per-gram need).

    micom's medium bound is on the LUMEN exchange, so it caps what the whole community
    may import, while `sat` is what one gram of one taxon needs. A taxon's draw on the
    lumen is its abundance times its own flux (the same unit rule `micom_tools` encodes
    and NB10 asserts as the lumen balance), so the community needs the abundance-weighted
    sum. Setting the diet below that starves somebody and brings the infeasibility back
    by a different route; this is the number that must not be undershot.

    Taxa with no measured requirement contribute nothing, which makes this a lower bound
    on demand in the strict sense -- but they are the taxa that grow without the compound
    at all, so their demand is zero by construction.
    """
    cols = [c for c in abundance.columns if c in sat.index]
    w = sat.reindex(cols).fillna(0.0).to_numpy()
    return pd.Series(abundance[cols].to_numpy() @ w, index=abundance.index)


# ===========================================================================
# stage 4 -- the reformulated medium
# ===========================================================================

def write_diet(diet: pd.DataFrame, supply: dict[str, float], notes: dict,
               renamed: dict, unmatched: list, census: pd.Series,
               n_models: int) -> None:
    """The diet as written, with the target compounds reduced to the measured need.

    NOT removed. Removal is what made 3 of 4 samples infeasible on 2026-10-23, and this
    script says why: some taxa cannot grow without them at all. Setting the supply to
    the largest minimum any cohort auxotroph needs, times a safety factor, keeps every
    taxon feasible while leaving the lumen concentration low enough that secretion is
    visible against it. Both numbers are in the header so nobody has to trust the choice
    without seeing it.
    """
    d = diet.copy()
    for base, level in supply.items():
        m = d["base"] == base
        if m.any():
            d.loc[m, "flux"] = level
    # Write the id the models actually use, back in the file's own `EX_x_e` dialect.
    d["reaction"] = d["model_rxn"].str.replace(r"\(e\)$", "_e", regex=True)

    header = [
        "# gut_diet_v2.csv -- western_diet_gut.csv with one reaction id corrected and",
        "# the two observable products reduced, not removed.",
        "# Generated by v2/micom/diet_screen.py; do not hand-edit.",
        "#",
        f"# source       {DIET_IN.name}, {len(diet)} reactions",
        f"# generated    {pd.Timestamp.now(tz='UTC').strftime('%Y-%m-%dT%H:%M:%SZ')}",
        "#",
    ]
    for old, new in renamed.items():
        header.append(
            f"# RENAMED  {old} -> {new}. The old id is in no AGORA2 model; the new one "
            f"is in {int(census.get(new, 0))} of {n_models} cohort models. A medium is "
            f"a dictionary lookup, so the old line supplied this nutrient to nobody and "
            f"said nothing about it.")
    if unmatched:
        header.append(
            f"# NO MODEL exchanges these {len(unmatched)}, so they supply nothing: "
            f"{', '.join(unmatched)}. Left in place; this is AGORA2's coverage, not an "
            f"error in the file.")
    header.append("#")
    for base, level in supply.items():
        was = float(diet.loc[diet["base"] == base, "flux"].iloc[0])
        n = notes.get(base, {})
        header.append(
            f"# {base:6s} {was:g} -> {level:g} mmol/gDW/h, a factor of "
            f"{was / level:.0f} less. {n.get('n_auxotroph', 0)} of "
            f"{n.get('n_tested', 0)} cohort species in AGORA2 cannot grow without it "
            f"at all, which is why it is reduced and not removed; they saturate at up "
            f"to {n.get('max_saturating_supply', float('nan')):.4g} mmol/gDW/h per gram, "
            f"the worst sample's abundance-weighted community demand is "
            f"{n.get('demand_max', float('nan')):.4g}, and this is {SAFETY:g}x that.")
    header.append("#")
    DIET_OUT.write_text("\n".join(header) + "\n" +
                        d[["reaction", "flux", "dilution"]].to_csv(index=False),
                        encoding="utf-8")


# ===========================================================================
# main
# ===========================================================================

def community_check(n_samples: int, diet_old: pd.DataFrame, diet_new: pd.DataFrame,
                    abund: pd.DataFrame, matched: pd.DataFrame, seed: int = 42):
    """Solve the same samples under both media and compare what actually happens.

    This is the only part of the script that can separate the medium from the solver.
    The screen above is single-model LPs; a sample fails or succeeds as a community QP,
    and nothing about one implies the other. So: build a few community models once, run
    micom's own `grow` on each twice, and count.

    Reported per medium: how many samples returned a solution at all, how many returned
    a solution that passes `micom_tools.solution_quality` (community growth > 0 and a
    balanced lumen -- micom returns infeasible points without raising), the community
    growth rate, and how many taxa secrete each target compound. That last column is the
    point of the diet change and the one number the 2026-10-23 note could not produce.
    """
    from micom.workflows import build, grow

    work = OUT / "community"
    work.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    samples = sorted(abund.index)
    pick = sorted(rng.choice(samples, size=min(n_samples, len(samples)),
                             replace=False).tolist())

    keep = set(matched["species"])
    long = (abund.loc[pick].stack().rename("abundance").reset_index())
    long.columns = ["sample_id", "species", "abundance"]
    long = long[(long["abundance"] > ABUND_CUTOFF) & long["species"].isin(keep)].copy()
    long["sample_id"] = long["sample_id"].astype(str)
    long["genus"] = long["species"].str.split().str[0]
    long["id"] = long["species"]
    tax = long[["sample_id", "id", "genus", "species", "abundance"]]

    model_db = str(next(AGORA_EXTRACTED.rglob("manifest.csv")).parent)
    models = work / "models"
    models.mkdir(parents=True, exist_ok=True)
    man_path = work / "manifest.csv"
    if man_path.exists():
        manifest = pd.read_csv(man_path, dtype={"sample_id": str})
        print(f"         reusing {len(manifest)} community models")
    else:
        t0 = time.time()
        print(f"         building {len(pick)} community models "
              f"({tax['species'].nunique()} distinct taxa) ...")
        manifest = build(tax, out_folder=str(models), model_db=model_db,
                         cutoff=ABUND_CUTOFF, threads=2)
        manifest.to_csv(man_path, index=False)
        print(f"         built in {time.time() - t0:.0f}s")
    manifest["sample_id"] = manifest["sample_id"].astype(str)

    rows = []
    for label, d in (("western_diet_gut", diet_old), ("gut_diet_v2", diet_new)):
        med = d[["reaction", "flux"]].copy()
        # NB10's convention: the community lumen exchange is `_m`, not `_e`.
        med["reaction"] = med["reaction"].str.replace(r"_e$", "_m", regex=True)
        res_path = work / f"grow_{label}.pkl"
        if res_path.exists():
            with open(res_path, "rb") as f:
                res = pickle.load(f)
            print(f"         reusing {label} solve")
        else:
            t0 = time.time()
            print(f"         solving {len(manifest)} samples on {label} ...")
            res = grow(manifest, model_folder=str(models), medium=med,
                       tradeoff=0.5, threads=2, strategy="pFBA", presolve=True)
            with open(res_path, "wb") as f:
                pickle.dump(res, f)
            print(f"         {label}: {time.time() - t0:.0f}s")

        ex = res.exchanges.copy()
        ex["sample_id"] = ex["sample_id"].astype(str)
        returned = set(ex["sample_id"])
        q = mt.solution_quality(ex, res.growth_rates)
        valid = set(q.loc[q["valid"], "sample_id"].astype(str))
        secreting = {}
        for b in TARGETS_IN_DIET + ["agm", "15dap", "orn"]:
            m = (ex["reaction"].str.contains(f"{b}_m|{b}\\(e\\)|{b}\\[e\\]", regex=True)
                 & (ex["taxon"] != "medium") & (ex["flux"] > 0))
            secreting[b] = int(ex.loc[m, "taxon"].nunique())
        rows.append({
            "medium": label,
            "samples": len(manifest),
            "returned_a_solution": len(returned),
            "passed_quality_check": len(valid),
            "median_community_growth": float(q["community_growth"].median()),
            **{f"taxa_secreting_{k}": v for k, v in secreting.items()},
        })
        print(f"         {label}: {len(returned)}/{len(manifest)} returned, "
              f"{len(valid)} valid, median community growth "
              f"{q['community_growth'].median():.4f}/h, secreting taxa {secreting}")
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "community_check.tsv", sep="\t", index=False)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0,
                    help="screen only the N most abundant species (a quick pass)")
    ap.add_argument("--community", type=int, default=0,
                    help="also re-solve N samples under both media (slow)")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    # ---- stage 1 -------------------------------------------------------
    sp, abund = cohort_species()
    n_samples = sp.attrs["n_samples"]
    idx = agora_index()
    have = set(idx["species_plain"])
    sp["in_agora2"] = sp["species"].isin(have)
    sp = sp.merge(idx[["species_plain", "path"]], left_on="species",
                  right_on="species_plain", how="left").drop(columns="species_plain")
    matched = sp[sp["in_agora2"]].copy()
    # Coverage PER SAMPLE, as NB10 reports it: what fraction of each sample's community
    # AGORA2 can model. A pooled "share of mean abundance" is a different and lower
    # number, because it averages over species rather than over samples.
    mcols = [c for c in abund.columns if c in set(matched["species"])]
    per_sample_cov = abund[mcols].sum(axis=1) / abund.sum(axis=1)
    cov = float(per_sample_cov.median())
    print(f"stage 1  cohort species above {ABUND_CUTOFF:g} in {n_samples} primary "
          f"samples: {len(sp)}; in AGORA2: {len(matched)}; median per-sample "
          f"abundance coverage {cov:.1%}")

    # ---- stage 1b: does the diet reach the models at all? ----------------------
    paths_all = matched.set_index("species")["path"]
    diet_raw = read_diet()
    census = exchange_census(paths_all.tolist())
    diet_df, renamed, unmatched = audit_diet(diet_raw, census, len(matched))
    print(f"stage 1b diet names vs model names: {len(diet_df)} lines, "
          f"{len(renamed)} renamed, {len(unmatched)} reaching no model")
    for old, new in renamed.items():
        line = diet_df.loc[diet_df["renamed_to"] == new]
        print(f"         RENAMED {old} -> {new} "
              f"({int(census.get(new, 0))} of {len(matched)} models carry it; "
              f"flux {float(line['flux'].iloc[0]):g}). The old id reached none, so this "
              f"nutrient was not being supplied.")
    if unmatched:
        print(f"         no cohort model exchanges: {unmatched} -- left in the file, "
              f"supplying nothing, which is a property of AGORA2 rather than an error")

    # ---- stage 2 -------------------------------------------------------
    diet = dict(zip(diet_df["model_rxn"], diet_df["flux"]))
    in_diet = [b for b in TARGETS_IN_DIET if b in set(diet_df["base"])]
    absent = [b for b in TARGETS_NOT_IN_DIET if b not in set(diet_df["base"])]
    print(f"stage 2  diet: {len(diet_df)} exchanges. Target compounds it supplies: "
          f"{in_diet}. Already absent, so observable as products without any change: "
          f"{absent}")

    withdraw = {b: [f"EX_{b}(e)"] for b in in_diet}
    withdraw["both"] = [f"EX_{b}(e)" for b in in_diet]

    todo = matched if not args.limit else matched.head(args.limit)
    rows, t0 = [], time.time()
    for i, r in enumerate(todo.itertuples(), 1):
        rec = {"species": r.species, "samples": r.samples,
               "mean_abundance": r.mean_abundance}
        try:
            rec.update(screen_model(r.path, diet, withdraw))
        except Exception as e:                      # a model that will not load is data
            rec["error"] = f"{type(e).__name__}: {e}"
        rows.append(rec)
        if i % 50 == 0 or i == len(todo):
            print(f"         {i}/{len(todo)} models, {time.time() - t0:.0f}s")
    scr = pd.DataFrame(rows)

    bad = scr[scr.get("error").notna()] if "error" in scr else scr.iloc[:0]
    ok = scr[scr["growth_full"] > GROWTH_FLOOR] if "growth_full" in scr else scr.iloc[:0]
    dead = (scr[(scr.get("growth_full", pd.Series(dtype=float)) <= GROWTH_FLOOR)]
            if "growth_full" in scr else scr.iloc[:0])
    print(f"stage 2  {len(scr)} screened: {len(ok)} grow on the diet as written, "
          f"{len(dead)} do not grow on it at all, {len(bad)} failed to load")

    # ---- stage 2b: taxa the medium cannot feed ----------------------------------
    #
    # WHAT THIS IS AND IS NOT. An earlier version of this block asserted that a taxon
    # which cannot grow makes its sample infeasible at any tradeoff. That is wrong, and
    # micom's source says so: `Community.cooperative_tradeoff` takes `min_growth=0.0`
    # and constrains `fraction` of the maximal COMMUNITY growth rate, not each taxon's
    # own maximum. A taxon that cannot grow is therefore feasible; it simply grows at
    # zero.
    #
    # It is still worth naming, for two reasons. First, it is a validity problem on its
    # own terms: a community model in which a member contributes no growth and no
    # metabolism is not modelling that member, and the sample's fluxes are those of a
    # smaller community than the one sampled. Second, it lowers the community growth
    # maximum that the tradeoff constraint is a fraction of, which is the regime where
    # the QP gets hard. Whether it is the cause of the observed failures is measured by
    # `--community`, not argued here.
    dead_sp = sorted(dead["species"].tolist())
    present = [c for c in abund.columns if c in set(dead_sp)]
    at_risk = (abund[present] > 0).any(axis=1) if present else pd.Series(
        False, index=abund.index)
    per_sample_dead = (abund[present] > 0).sum(axis=1) if present else 0
    lost = abund[present].sum(axis=1) if present else pd.Series(0.0,
                                                               index=abund.index)
    print(f"stage 2b {len(dead_sp)} species cannot grow on this diet at all. "
          f"{int(at_risk.sum())} of {len(abund)} samples ({at_risk.mean():.1%}) "
          f"contain at least one; they carry a median {lost.median():.2%} and up to "
          f"{lost.max():.2%} of a sample's abundance, and contribute no growth and no "
          f"metabolism to it")
    if len(dead_sp):
        print(f"         they are: {dead_sp}")

    # What each of them is missing, by name.
    gaps = []
    for s in dead_sp:
        extra, err = complete_for(paths_all[s], diet)
        gaps.append({"species": s, "n_added": len(extra), "error": err,
                     "added": ";".join(sorted(extra)) if extra else ""})
        print(f"         {s}: " + (err if err else
                                   f"needs {len(extra)} further exchange(s): "
                                   f"{sorted(extra)[:8]}"))
    gaps_df = pd.DataFrame(gaps)
    gaps_df.to_csv(OUT / "diet_gaps_nongrowing.tsv", sep="\t", index=False)

    notes = {}
    for b in in_diet:
        col = f"growth_without_{b}"
        aux = ok[(ok[col].notna()) & (ok[col] <= GROWTH_FLOOR)]
        notes[b] = {"n_tested": int(len(ok)), "n_auxotroph": int(len(aux)),
                    "auxotrophs": sorted(aux["species"].tolist())}
        print(f"         withdraw {b}: {len(aux)} of {len(ok)} cohort species stop "
              f"growing" + (f" -- {sorted(aux['species'])[:4]}" if len(aux) else ""))
    if "growth_without_both" in ok:
        both = ok[(ok["growth_without_both"].notna()) &
                  (ok["growth_without_both"] <= GROWTH_FLOOR)]
        print(f"         withdraw both: {len(both)} of {len(ok)} stop growing")

    # ---- stage 3 -------------------------------------------------------
    supply, curves = {}, {}
    paths = paths_all
    for b in in_diet:
        aux = notes[b]["auxotrophs"]
        was = float(diet_df.loc[diet_df["base"] == b, "flux"].iloc[0])
        if not aux:
            supply[b] = 0.0
            notes[b].update({"max_saturating_supply": 0.0, "demand_max": 0.0})
            print(f"stage 3  {b}: no auxotroph, so it can be withdrawn outright")
            continue
        sat, cur = {}, {}
        for s in aux:
            sat[s], cur[s] = saturating_supply(paths[s], diet, f"EX_{b}(e)")
        sat = pd.Series(sat).dropna()
        curves[b] = cur
        # The per-gram requirement is not the diet bound. micom caps what the WHOLE
        # community may import from the lumen, so the bound has to cover the
        # abundance-weighted sum over every auxotroph in a sample, in the worst sample.
        dem = community_demand(sat, abund)
        level = float(dem.max()) * SAFETY
        supply[b] = level
        notes[b].update({
            "max_saturating_supply": float(sat.max()),
            "median_saturating_supply": float(sat.median()),
            "saturating_supply_per_species": {k: float(v) for k, v in sat.items()},
            "demand_median": float(dem.median()),
            "demand_max": float(dem.max()),
            "samples_demanding": int((dem > 0).sum()),
        })
        print(f"stage 3  {b}: {len(aux)} auxotrophs saturate between {sat.min():.4g} "
              f"and {sat.max():.4g} mmol/gDW/h per gram. Abundance-weighted community "
              f"demand: median {dem.median():.4g}, max {dem.max():.4g} over "
              f"{len(dem)} samples. Diet {was:g} -> {level:.4g} "
              f"({SAFETY:g}x the worst sample, {was / level:.0f}x less than now)")

    # ---- stage 4 -------------------------------------------------------
    write_diet(diet_df, supply, notes, renamed, unmatched, census, len(matched))
    print(f"stage 4  wrote {DIET_OUT}")

    comm = None
    if args.community:
        print(f"stage 5  community check on {args.community} samples, both media")
        new_df = pd.read_csv(DIET_OUT, comment="#")
        comm = community_check(args.community, diet_raw, new_df, abund, matched)

    scr.to_csv(OUT / "diet_auxotrophy_screen.tsv", sep="\t", index=False)
    sp.to_csv(OUT / "cohort_species_agora2.tsv", sep="\t", index=False)
    prov = {
        "built_utc": pd.Timestamp.now(tz="UTC").strftime("%Y-%m-%dT%H:%M:%SZ"),
        "diet_in": {"file": DIET_IN.name, "reactions": int(len(diet_df)),
                    "renamed": renamed,
                    "reaching_no_model": unmatched},
        "cohort": {"primary_samples": n_samples, "species": int(len(sp)),
                   "in_agora2": int(len(matched)),
                   "abundance_coverage": round(cov, 4),
                   "screened": int(len(scr))},
        "growth_floor": GROWTH_FLOOR,
        "safety_factor": SAFETY,
        "non_growing": {
            "species": dead_sp,
            "samples_containing_one": int(at_risk.sum()),
            "samples_total": int(len(abund)),
            "fraction_of_samples": round(float(at_risk.mean()), 4),
            "median_per_sample": (float(np.median(per_sample_dead))
                                  if present else 0.0),
            "gaps": gaps,
        },
        "targets_supplied_by_diet": in_diet,
        "targets_already_absent": absent,
        "withdrawal": {b: {k: v for k, v in n.items()
                           if k != "saturating_supply_per_species"}
                       for b, n in notes.items()},
        "dose_response": curves,
        "new_supply": supply,
        "community_check": (comm.to_dict("records") if comm is not None
                            else "not run (--community N)"),
        "software": {"python": platform.python_version(),
                     "pandas": pd.__version__},
    }
    (OUT / "diet_screen_provenance.json").write_text(json.dumps(prov, indent=2))
    print(json.dumps({k: prov[k] for k in ("cohort", "withdrawal", "new_supply")},
                     indent=2, default=str))


if __name__ == "__main__":
    main()
