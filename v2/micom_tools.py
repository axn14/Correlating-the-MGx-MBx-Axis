"""Pure helpers for NB10 (MICOM community flux).

Nothing here imports micom, so `checks.py` (D8) tests every rule on synthetic tables in the
pinned environment. Each helper encodes a fact about micom 0.39.1 output that NB10 once got
wrong:

- micom rewrites taxon ids: `Community` replaces every run of characters outside
  [A-Za-z0-9_] with "_" (micom/community.py, `taxonomy.id.str.replace`), so
  "Dialister hominis" comes back from `grow()` as "Dialister_hominis".
  `taxon_name_map` maps the ids back to the names every other notebook uses.
- A taxon's exchange flux is per gram of THAT taxon (mmol/gDW_taxon/h). In the community
  model each taxon exchange feeds the shared lumen with coefficient = its abundance, so the
  lumen balance is  medium flux = sum_t abundance_t * flux_t. A species' contribution is
  abundance * flux; adding raw taxon fluxes mixes units.
- `grow()` drops a sample whose solve failed and raises only if every sample fails.
  `missing_samples` makes the loss visible; downstream a missing sample is "not solved"
  (NaN), never "produced nothing" (0).
- Taxon metabolites are written "ptrc[e]", community-level ones "ptrc_m", exchange
  reactions "EX_ptrc(e)" / "EX_ptrc_m". `classify_metabolite` accepts all of them.
"""
import re

import numpy as np
import pandas as pd

# Polyamine KEGG id -> (name, VMH base ids). Matching is on the bare base id.
POLYAMINE_VMH = {
    "C00134": ("Putrescine", ["ptrc"]),
    "C00315": ("Spermidine", ["spmd"]),
    "C00750": ("Spermine", ["sprm"]),
    "C00179": ("Agmatine", ["agm"]),
    "C01672": ("Cadaverine", ["15dap", "cadav"]),
    "C00077": ("Ornithine", ["orn", "orndp"]),
    "C00062": ("Arginine", ["arg__L", "arg_L"]),
}
NAME_TO_KEGG = {name: k for k, (name, _ids) in POLYAMINE_VMH.items()}
_VMH_BASE = {b.lower(): k for k, (_n, ids) in POLYAMINE_VMH.items() for b in ids}
# compartment tag at the end of a metabolite / exchange id: [e], (e), _e, _m, _lu, _u, _c
_COMPARTMENT = re.compile(r"(\[(e|m|lu|u|c)\]|\((e|m|lu|u|c)\)|_(e|m|lu|u|c))$", re.IGNORECASE)

GROWTH_FLOOR = 1e-6  # micom's growth threshold for "this taxon grows"


def metabolite_base_id(met_id) -> str:
    """Bare, lower-case VMH id: "EX_ptrc(e)", "ptrc[e]", "ptrc_m" -> "ptrc"."""
    s = re.sub(r"^EX_", "", str(met_id))
    return _COMPARTMENT.sub("", s).strip().lower()


def classify_metabolite(met_id):
    """(KEGG, name) when the id is a tracked polyamine, else (None, None)."""
    k = _VMH_BASE.get(metabolite_base_id(met_id))
    return (k, POLYAMINE_VMH[k][0]) if k else (None, None)


def micom_taxon_id(name) -> str:
    """The id micom gives a taxon named `name` (micom 0.39.1 community.py)."""
    return re.sub(r"[^A-Za-z0-9_]+", "_", str(name))


def taxon_name_map(species) -> dict:
    """micom taxon id -> original species name. Raises if two names collapse to one id,
    because their fluxes could then not be told apart."""
    out = {}
    for s in pd.unique(pd.Series(list(species), dtype=object)):
        tid = micom_taxon_id(s)
        if tid in out and out[tid] != s:
            raise ValueError(f"species {out[tid]!r} and {s!r} both become micom id {tid!r}")
        out[tid] = s
    return out


def choose_tradeoff(sweep: pd.DataFrame, floor: float = GROWTH_FLOOR):
    """micom's rule: the LARGEST tradeoff that still lets the largest share of taxa grow.

    Low tradeoffs tie at "every taxon grows" because they cap community growth at that
    fraction of its maximum; taking the first tie (idxmax) picks the most restrictive
    value. The sweep's tradeoff-NaN rows (plain optimisation) are not candidates.
    Returns (tradeoff, per-tradeoff stats).
    """
    d = sweep.dropna(subset=["tradeoff"]).assign(alive=lambda x: x["growth_rate"] > floor)
    stats = d.groupby("tradeoff").agg(frac_alive=("alive", "mean"),
                                      median_growth=("growth_rate", "median"))
    best = stats["frac_alive"].max()
    return float(stats.index[np.isclose(stats["frac_alive"], best)].max()), stats


def missing_samples(requested, exchanges: pd.DataFrame) -> list:
    """Requested sample ids that `grow()` returned no solution for."""
    got = set(exchanges["sample_id"].astype(str))
    return sorted(set(map(str, requested)) - got)


BALANCE_TOL_FACTOR = 10   # residual allowed, in units of the solve tolerance


def solution_quality(exchanges: pd.DataFrame, growth: pd.DataFrame,
                     floor: float = GROWTH_FLOOR) -> pd.DataFrame:
    """Per sample: is the solution micom returned a real one?

    micom 0.39.1 `grow()` keeps whatever the last solve left behind: with the hybrid
    OSQP/HiGHS solver the cooperative-tradeoff + pFBA step ends INFEASIBLE for some
    communities, and grow() returns that point without an error -- community growth ~0,
    hundreds of spurious exchanges, and a lumen that does not balance. Valid means
    community growth > floor AND, over every metabolite, |medium flux - sum_t a_t v_t|
    <= BALANCE_TOL_FACTOR * the solve tolerance (micom drops rows below the tolerance,
    so a feasible solution misses the balance by at most about twice it).
    Columns: sample_id, community_growth, max_balance_residual, tolerance, valid.
    """
    ex = exchanges.copy()
    ex["sample_id"] = ex["sample_id"].astype(str)
    ex["base"] = ex["metabolite"].map(metabolite_base_id)
    tx = ex[ex["taxon"].astype(str) != "medium"]
    md = ex[ex["taxon"].astype(str) == "medium"]
    av = (tx["flux"] * tx["abundance"]).groupby([tx["sample_id"], tx["base"]]).sum()
    med = md.groupby(["sample_id", "base"])["flux"].sum()
    resid = av.sub(med, fill_value=0.0).abs().groupby(level="sample_id").max()
    g = growth.copy()
    g["sample_id"] = g["sample_id"].astype(str)
    g = g[g["taxon"].astype(str) != "medium"]
    cg = (g["growth_rate"] * g["abundance"]).groupby(g["sample_id"]).sum()
    tol = ex.groupby("sample_id")["tolerance"].max()
    out = pd.DataFrame({"community_growth": cg, "max_balance_residual": resid, "tolerance": tol})
    out.index.name = "sample_id"
    out = out.reset_index()
    out["valid"] = ((out["community_growth"] > floor)
                    & (out["max_balance_residual"] <= BALANCE_TOL_FACTOR * out["tolerance"]))
    return out


def _taxon_rows(exchanges: pd.DataFrame, name_map: dict) -> pd.DataFrame:
    tx = exchanges[exchanges["taxon"].astype(str) != "medium"].copy()
    unknown = sorted(set(tx["taxon"].astype(str)) - set(name_map))
    if unknown:
        raise KeyError(f"{len(unknown)} micom taxon id(s) have no species name, "
                       f"e.g. {unknown[:3]}; was the taxonomy table changed after the build?")
    tx["Species"] = tx["taxon"].astype(str).map(name_map)
    return tx


def polyamine_exchanges(exchanges: pd.DataFrame, name_map: dict) -> pd.DataFrame:
    """Taxon-level polyamine exchange rows with species names and community-scale flux.

    Columns: sample_id, Species, Polyamine, KEGG, direction, flux_specific
    (mmol/gDW_taxon/h, positive = secreted into the lumen), abundance, flux_community
    (= abundance * flux_specific, mmol/gDW_community/h).
    """
    tx = _taxon_rows(exchanges, name_map)
    kn = tx["metabolite"].map(classify_metabolite)
    tx["KEGG"] = kn.map(lambda t: t[0])
    tx["Polyamine"] = kn.map(lambda t: t[1])
    tx = tx.dropna(subset=["KEGG"])
    if tx["abundance"].isna().any():
        raise ValueError("taxon exchange rows without abundance; cannot scale to the community")
    tx = tx.rename(columns={"flux": "flux_specific"})
    tx["flux_community"] = tx["flux_specific"] * tx["abundance"]
    tx["sample_id"] = tx["sample_id"].astype(str)
    return tx[["sample_id", "Species", "Polyamine", "KEGG", "direction",
               "flux_specific", "abundance", "flux_community"]].reset_index(drop=True)


def species_export(pa_ex: pd.DataFrame) -> pd.DataFrame:
    """Per (sample, species, polyamine): secretion into the lumen, community scale.

    Rows exist only where the species secretes; a species that is absent, takes the
    compound up, or does not exchange it contributes no row (0 export downstream).
    """
    ex = pa_ex[pa_ex["flux_specific"] > 0]
    return (ex.groupby(["sample_id", "Species", "Polyamine"], as_index=False)
              .agg(export_community=("flux_community", "sum"),
                   export_specific=("flux_specific", "sum")))


def community_balance(exchanges: pd.DataFrame, pa_ex: pd.DataFrame) -> pd.DataFrame:
    """Per (sample, polyamine): the lumen balance.

    Net_medium_flux is the community's net exchange with the environment (positive = net
    production leaving the lumen, negative = net uptake from the diet). Gross_export /
    Gross_uptake are the abundance-weighted taxon sides; Net_medium_flux equals
    Gross_export - Gross_uptake up to solver tolerance (checked in D8 and in NB10).
    """
    md = exchanges[exchanges["taxon"].astype(str) == "medium"].copy()
    kn = md["metabolite"].map(classify_metabolite)
    md["Polyamine"] = kn.map(lambda t: t[1])
    md = md.dropna(subset=["Polyamine"])
    md["sample_id"] = md["sample_id"].astype(str)
    net = md.groupby(["sample_id", "Polyamine"])["flux"].sum().rename("Net_medium_flux")
    v = pa_ex["flux_community"]
    side = (pa_ex.assign(Gross_export=v.clip(lower=0), Gross_uptake=(-v).clip(lower=0),
                         N_exporting_taxa=(v > 0).astype(int),
                         N_uptaking_taxa=(v < 0).astype(int))
                 .groupby(["sample_id", "Polyamine"])
                 [["Gross_export", "Gross_uptake", "N_exporting_taxa", "N_uptaking_taxa"]]
                 .sum())
    out = side.join(net, how="outer").reset_index()
    for c in ("Gross_export", "Gross_uptake", "Net_medium_flux"):
        out[c] = out[c].fillna(0.0)
    for c in ("N_exporting_taxa", "N_uptaking_taxa"):
        out[c] = out[c].fillna(0).astype(int)
    return out
