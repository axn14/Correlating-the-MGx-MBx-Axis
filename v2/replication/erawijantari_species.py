"""Class X (CONTRACTS 5d): species-level replication of Yachida's species-metabolite pairs.

Erawijantari 2020 is the only other cohort on disk with shotgun species AND a CE-TOFMS
metabolome from the same platform and group as Yachida, so it is the only place a
SPECIES-level association can be replicated at all (Sinha and Kim are 16S, genus level,
on other metabolomics platforms).

Independence. 53 of Erawijantari's 54 healthy subjects are Yachida subjects, and their
data are the SAME samples: all 53 metabolome profiles are identical and the species
profiles agree at median r = 1.00 (checked 2026-11-11, `Shared.w.YACHIDA_2019`, matched on
`Subject.Num`). They are excluded. What remains is 42 post-gastrectomy patients and one
healthy subject; the healthy one is dropped too, so the replication sample is one
homogeneous group. **n = 42, all post-gastrectomy** -- an independent cohort, but a small
one with a surgically altered gut. A failure to replicate here is weak evidence against a
pair; a success is real.

Fixed before the first run (2026-11-11):
- Pairs: every Yachida pair in NB08's `pair_evidence_nb08_v2.csv` whose metabolite is one of
  the registry's targets (`v2/targets/out/target_registry.tsv`).
- Measurable: the species passes `PREVALENCE_SPECIES` in the 42 samples and the metabolite
  passes `PREVALENCE_MTB`, both matched by name (species by its `s__` label, metabolite by
  KEGG id). Unmeasurable pairs are `available = False`, never failures.
- Transforms: NB01's own -- `filter_by_prevalence` + `clr_transform` for species,
  `transform_metabolites` (log10 + PQN) for the metabolome.
- Statistic: `partial_spearman_grid`, adjusting for Age, BMI, sex, ever-smoked and alcohol
  (yes/no). Yachida used Brinkman index and continuous alcohol; Erawijantari records
  status only, so these are the closest available.
- Family: all measurable target pairs; BH across them.
- Support: same sign as Yachida's `Rho_partial` and q < 0.05. Contradiction: opposite sign
  and q < 0.05.

Output: out/erawijantari_species_replication.tsv, out/erawijantari_provenance.json.
Run:    python v2/replication/erawijantari_species.py --results "D:/KI Research-dardel/<run>"
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
COHORT = "ERAWIJANTARI-GASTRIC-CANCER-2020"
REGISTRY = V2 / "targets" / "out" / "target_registry.tsv"


def short_species(col: str) -> str:
    return col.split("s__")[-1] if "s__" in col else col


def load_erawijantari():
    d = utils.DATA_DIR
    meta = pd.read_csv(d / f"{COHORT} metadata.tsv", sep="\t", dtype={"Sample": str})
    keep = meta[(meta["Shared.w.YACHIDA_2019"].astype(str).str.upper() != "TRUE")
                & (meta["Study.Group"] == "Gastrectomy")].copy()
    assert not keep["Shared.w.YACHIDA_2019"].astype(str).str.upper().eq("TRUE").any()
    sp = pd.read_csv(d / f"{COHORT} species.tsv", sep="\t", dtype={"Sample": str})
    mt = pd.read_csv(d / f"{COHORT} mtb.tsv", sep="\t", dtype={"Sample": str})
    sp = sp[sp["Sample"].isin(keep["Sample"])].reset_index(drop=True)
    mt = mt[mt["Sample"].isin(keep["Sample"])].reset_index(drop=True)
    return keep, sp, mt


def covariates(meta: pd.DataFrame) -> pd.DataFrame:
    z = pd.DataFrame({"Sample": meta["Sample"].values})
    z["Age"] = pd.to_numeric(meta["Age"], errors="coerce").values
    z["BMI"] = pd.to_numeric(meta["BMI"], errors="coerce").values
    z["Gender_num"] = (meta["Gender"].astype(str) == "Male").astype(float).values
    smk = meta["SmokingStatus"].astype(str)
    z["Ever_smoked"] = np.where(smk.isin(["Smoking", "StoppedSmoking"]), 1.0,
                                np.where(smk == "NotSmoking", 0.0, np.nan))
    alc = meta["AlcoholConsumption"].astype(str)
    z["Alcohol_yes"] = np.where(alc == "Yes", 1.0, np.where(alc == "No", 0.0, np.nan))
    return z


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", type=Path, required=True)
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    reg = pd.read_csv(REGISTRY, sep="\t")
    targets = set(reg["target"])
    pairs = pd.read_csv(a.results / "Results" / "tables" / "pair_evidence_nb08_v2.csv",
                        usecols=["Species", "Metabolite", "Rho_partial", "QValue_partial"])
    pairs = pairs[pairs["Metabolite"].isin(targets)].copy()

    meta, sp, mt = load_erawijantari()
    z = covariates(meta).set_index("Sample")
    complete = z.dropna().index
    sp = sp[sp["Sample"].isin(complete)].reset_index(drop=True)
    mt = mt[mt["Sample"].isin(complete)].reset_index(drop=True)
    assert list(sp["Sample"]) == list(mt["Sample"]) or set(sp["Sample"]) == set(mt["Sample"])
    mt = mt.set_index("Sample").loc[sp["Sample"]].reset_index()
    n = len(sp)

    sp_f = utils.filter_by_prevalence(sp, utils.PREVALENCE_SPECIES)
    sp_c = utils.clr_transform(sp_f)
    mt_f = utils.filter_by_prevalence(mt, utils.PREVALENCE_MTB)
    mt_t, _ = utils.transform_metabolites(mt_f)
    sp_c = sp_c.set_index("Sample")
    mt_t = mt_t.set_index("Sample").loc[sp_c.index]
    zz = z.loc[sp_c.index]

    sp_cols = {short_species(c): c for c in sp_c.columns}
    kegg_e = {c.split("_", 1)[0]: c for c in mt_t.columns}
    pairs["kegg"] = pairs["Metabolite"].str.split("_").str[0]
    pairs["e_species"] = pairs["Species"].map(sp_cols)
    pairs["e_metabolite"] = pairs["kegg"].map(kegg_e)
    meas = pairs.dropna(subset=["e_species", "e_metabolite"]).copy()

    xs = sorted(meas["e_species"].unique())
    ys = sorted(meas["e_metabolite"].unique())
    r, p = utils.partial_spearman_grid(sp_c[xs].values, mt_t[ys].values, zz.values)
    xi = {c: i for i, c in enumerate(xs)}
    yi = {c: j for j, c in enumerate(ys)}
    meas["Rho_erawijantari"] = [r[xi[s], yi[m]] for s, m in zip(meas["e_species"], meas["e_metabolite"])]
    meas["P_erawijantari"] = [p[xi[s], yi[m]] for s, m in zip(meas["e_species"], meas["e_metabolite"])]
    ok = meas["P_erawijantari"].notna()
    meas["Q_erawijantari"] = np.nan
    meas.loc[ok, "Q_erawijantari"] = multipletests(meas.loc[ok, "P_erawijantari"], method="fdr_bh")[1]
    same = np.sign(meas["Rho_erawijantari"]) == np.sign(meas["Rho_partial"])
    sig = meas["Q_erawijantari"] < utils.FDR_THRESHOLD
    meas["Source_replication_species_available"] = ok
    meas["Source_replication_species_support"] = ok & same & sig
    meas["Replication_species_contradicted"] = ok & (~same) & sig

    out = pairs.merge(meas[["Species", "Metabolite", "Rho_erawijantari", "P_erawijantari",
                            "Q_erawijantari", "Source_replication_species_available",
                            "Source_replication_species_support",
                            "Replication_species_contradicted"]],
                      on=["Species", "Metabolite"], how="left")
    for c in ("Source_replication_species_available", "Source_replication_species_support",
              "Replication_species_contradicted"):
        out[c] = out[c].fillna(False).astype(bool)
    out.to_csv(OUT / "erawijantari_species_replication.tsv", sep="\t", index=False)

    # The smallest |rho| that could reach p < 0.05 two-sided here, before any FDR.
    df_t = n - 2 - zz.shape[1]
    t_crit = stats.t.ppf(0.975, df_t)
    rho_min = float(t_crit / np.sqrt(df_t + t_crit ** 2))
    prov = {"script": "v2/replication/erawijantari_species.py",
            "run_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "results_root": str(a.results), "n_samples": int(n),
            "excluded": "53 Yachida-shared subjects (identical samples) + 1 healthy",
            "covariates": list(zz.columns),
            "target_pairs_in_yachida": int(len(pairs)),
            "measurable_pairs": int(ok.sum()),
            "species_measurable": len(xs), "metabolites_measurable": len(ys),
            "min_detectable_abs_rho_uncorrected": round(rho_min, 3),
            "supported": int(out["Source_replication_species_support"].sum()),
            "contradicted": int(out["Replication_species_contradicted"].sum()),
            "nominal_p05_same_sign": int((ok & same & (meas["P_erawijantari"] < 0.05)).sum()),
            "sign_agreement_rate": float(same[ok].mean()) if ok.any() else None}
    (OUT / "erawijantari_provenance.json").write_text(json.dumps(prov, indent=1))
    print(json.dumps(prov, indent=1))


if __name__ == "__main__":
    main()
