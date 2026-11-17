"""Phase E1: can each diet-supplied target be withdrawn from the medium?

The review's point: if the medium supplies a target, a positive community-level exchange
can be the supply not taken up rather than production. A TAXON's net export cannot be
pass-through (its own network must make the compound), so class M already counts only
per-taxon net export with FVA minimum > 0 (CONTRACTS 5d). This screen decides the medium
itself, target by target, for the 14 exportable targets the diet supplies (the 12 others
are not in the diet at all, so their export is production by construction).

Rule, fixed before the first run (2026-11-11): withdraw a target's supply (set it to 0) iff
NO cohort species with an AGORA2 model that grows on the full diet stops growing without
it (single-species LP, growth < 1e-6). Otherwise keep the diet's supply -- dietary amino
acids do reach the colon -- and record which species depend on it; class M for that
target then rests on per-taxon net export only, as it does anyway.

Output: out/target_medium.tsv (per target) and out/gut_diet_v2_targets.csv (the medium the
batch uses: gut_diet_v2 minus the withdrawable targets).
Run:    python v2/micom/target_medium.py      (~5 min, no HPC)
"""
from __future__ import annotations

import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))

import diet_screen as ds  # noqa: E402

REGISTRY = HERE.parent / "targets" / "out" / "target_registry.tsv"
OUT = HERE / "out"
GROWTH_FLOOR = 1e-6


def kegg_by_base(paths: list[str]) -> dict[str, set[str]]:
    base2k: dict[str, set[str]] = {}
    for p in paths:
        with open(p, encoding="utf-8") as f:
            for met in json.load(f)["metabolites"]:
                k = str(met.get("annotation", {}).get("kegg.compound", "NaN"))
                if k != "NaN":
                    base2k.setdefault(met["id"].split("[")[0], set()).update(
                        x.strip() for x in k.split(";"))
    return base2k


def _screen(args):
    path, diet, withdraw = args
    try:
        return path, ds.screen_model(path, diet, withdraw)
    except Exception as err:            # one unreadable model is reported, not fatal
        return path, {"error": str(err)}


def main() -> None:
    reg = pd.read_csv(REGISTRY, sep="\t")
    cs = pd.read_csv(HERE / "out" / "cohort_species_agora2.tsv", sep="\t")
    cs = cs[cs["in_agora2"] == True].dropna(subset=["path"])  # noqa: E712
    paths = cs["path"].tolist()
    # The batch's medium is gut_diet_v2, not the western diet ds.read_diet() reads.
    d = pd.read_csv(ds.DIET_OUT, comment="#")
    d["base"] = d["reaction"].str.replace(r"^EX_", "", regex=True).str.replace(r"_e$", "", regex=True)
    d["model_rxn"] = "EX_" + d["base"] + "(e)"
    b2k = kegg_by_base(paths[:300])
    d["kegg"] = d["base"].map(lambda b: b2k.get(b, set()))

    supplied = {}
    for _, t in reg[reg["export_gate"] == True].iterrows():  # noqa: E712
        rows = d[d["kegg"].map(lambda ks: t["kegg"] in ks)]
        if len(rows):
            supplied[t["target"]] = rows["model_rxn"].tolist()
    diet = dict(zip(d["model_rxn"], d["flux"]))

    with ProcessPoolExecutor(max_workers=6) as ex:
        res = dict(ex.map(_screen, [(p, diet, supplied) for p in paths]))

    abund = cs.set_index("path")["mean_abundance"]
    rows = []
    for tgt, rxns in supplied.items():
        lost, can = [], 0
        for p, r in res.items():
            if "error" in r or r.get("growth_full", 0) <= GROWTH_FLOOR:
                continue
            g = r.get(f"growth_without_{tgt}", np.nan)
            if np.isnan(g):
                continue
            can += 1
            if g <= GROWTH_FLOOR:
                lost.append(p)
        names = cs.set_index("path").loc[lost, "species"].tolist() if lost else []
        rows.append({"target": tgt, "diet_reactions": ";".join(rxns),
                     "diet_flux": float(d.loc[d["model_rxn"].isin(rxns), "flux"].sum()),
                     "species_growing_full_diet": can, "species_lost_without": len(lost),
                     "mean_abundance_lost": float(abund.reindex(lost).sum()),
                     "withdraw": len(lost) == 0,
                     "species_lost": ";".join(sorted(names)[:30])})
    tab = pd.DataFrame(rows)
    OUT.mkdir(exist_ok=True)
    tab.to_csv(OUT / "target_medium.tsv", sep="\t", index=False)

    drop = {r for t, rs in supplied.items() if tab.set_index("target").at[t, "withdraw"] for r in rs}
    med = pd.read_csv(ds.DIET_OUT, comment="#")
    med["model_rxn"] = "EX_" + med["reaction"].str.replace(r"^EX_", "", regex=True) \
                                              .str.replace(r"_e$", "", regex=True) + "(e)"
    out = med[~med["model_rxn"].isin(drop)][["reaction", "flux"]]
    with open(OUT / "gut_diet_v2_targets.csv", "w", newline="") as f:
        f.write("# gut_diet_v2 minus the diet-supplied targets that no cohort species needs\n"
                "# (v2/micom/target_medium.py, rule in its docstring). Do not hand-edit.\n")
        out.to_csv(f, index=False, lineterminator="\n")
    print(tab.drop(columns="species_lost").to_string(index=False))
    print(f"withdrawn: {sorted(tab.loc[tab.withdraw, 'target'])}")


if __name__ == "__main__":
    main()
