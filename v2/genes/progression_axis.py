#!/usr/bin/env python3
"""Put the gene-level evidence on the CRC progression axis, and draw it.

WHAT THIS ADDS TO `extract_gene_families.R`
-------------------------------------------
That job answered "is the gene there, and in which species". It produced 207,026
nonzero (row, sample, abundance) triples over 1,267 UniRef90 clusters, with no
reference to who the samples are. This script is the step that makes those numbers
an argument: it joins them to the cohort's clinical axis and asks whether the
biosynthetic capacity for each metabolite differs between healthy guts, adenomas,
early carcinoma and advanced carcinoma.

WHAT IT FOUND, SO THE READER KNOWS WHAT THE CODE IS FOR
-------------------------------------------------------
Not a shift in total capacity. Every per-EC aggregate is null (smallest BH
q = 0.83) and no single cluster separates two stages pairwise (smallest q = 0.19
over 994 tests). What moves is WHICH COPY of an activity is present: twelve
clusters trend monotonically with stage at q < 0.05, and within one EC they trend
in OPPOSITE directions, split by the organism carrying them. Arginine
decarboxylase rises in the copies attributed to Bilophila wadsworthia and
Butyricimonas virosa and falls in those attributed to Lachnospira pectinoschiza
and Roseburia faecis; ornithine decarboxylase rises in the Fusobacterium
nucleatum copy; glutamate decarboxylase rises in the Parvimonas micra copy; the
lysine decarboxylase that moves falls, and belongs to Roseburia intestinalis.
The potential is preserved while its source is replaced.

THE TWO CONTROLS THEN BOUND THAT CLAIM, AND BOTH BITE.

The background control says these enzymes are not, as a set, unusually
stage-responsive: at p < 0.05, 8.6% of them trend against 13.3% of
prevalence-matched non-target gene families -- they trend LESS, significantly so
(Fisher p = 0.021). The enrichment is in the tail only: 2.4% against 0.85% at
p < 0.001, and the single strongest target signal (p = 8.7e-08) is two orders of
magnitude beyond the smallest p anywhere in a five-times-larger background
(1.4e-05). A sentence claiming that polyamine and amine biosynthesis is broadly
engaged along the axis would be false; a sentence naming specific copies is not.

The carrier control says the two strongest signals are their organism arriving.
Holding the carrier's own relative abundance fixed, the Parvimonas micra
glutamate decarboxylase falls from rho 0.31 to 0.11 (p = 0.08) and the
Fusobacterium nucleatum ornithine decarboxylase from 0.25 to 0.08 (p = 0.21),
while those two organisms track stage at p = 2.6e-11 and 3.1e-10. Four trends do
survive, led by the Butyricimonas virosa arginine decarboxylase (rho 0.22 ->
0.22, q = 0.003) whose carrier does not itself trend (p = 0.31).

So the defensible claim is attribution, not independent gene-level change: these
organisms arrive carrying these specific decarboxylase activities. That is what
Results 4.6 asks for and it is weaker than "capacity increases", which the data
do not support.

THE JOIN, AND WHY IT NEEDS ITS OWN FILE
---------------------------------------
The gene matrix is keyed by DDBJ BioSample accession (SAMD00114718). The cohort
tables this project has used for nine months are keyed by Yachida's own sample
number (10021). Nothing in either file connects them. The bridge is
curatedMetagenomicData's curation table, where `subject_id` is literally
"sub_" + Yachida's number and `sample_id` is the accession. That table is
committed under `ref/` rather than fetched at run time, because a silent change
upstream would silently re-label every sample in every figure below.

The join is asserted, not assumed: all 347 cohort samples must resolve to an
accession, every accession must be a column of the gene matrix, and subject ids
must be unique in the curation table (616 samples, 616 subjects -- Yachida has no
repeat sampling, so a subject id IS a sample id here). Any failure raises.

THE AXIS IS THE PIPELINE'S, NOT A NEW ONE
-----------------------------------------
Group policy comes from `v2/utils.py` and nothing is re-derived here:

    primary       Healthy (127) -> Early_CRC (96) -> Advanced_CRC (54)   n=277
    early         Healthy (127) -> Adenoma (40) -> Early_CRC -> Advanced n=317
    sensitivity   Healthy vs PostSurgery (30)                            n=157

Adenoma (Yachida "MP") and PostSurgery ("HS") are NaN on the primary ordinal
scale on purpose, so neither can enter a trend test by accident. Each axis is its
own BH family, never merged -- the same rule `COMPARISONS_EARLY` records.

WHAT IS MEASURED PER FEATURE, AND WHY SEVERAL NUMBERS RATHER THAN ONE
---------------------------------------------------------------------
A gene family's abundance is zero in most samples for most clusters, so a single
fold change mixes two different biological statements: the gene appeared, and the
gene got commoner. They are reported separately.

    prevalence    fraction of samples with a nonzero value; Fisher exact on the
                  2x2 of detected/not against group. "Is the capacity present?"
    log2FC        log2 of the ratio of group MEAN relative abundances. Arithmetic
                  means, not geometric: the quantity being compared is a
                  compositional share of the gene pool, and a share is additive.
                  Zeros therefore enter honestly instead of being pseudocounted.
                  NaN, never 0, when a group mean is zero -- "infinitely more" is
                  not a fold change and must not be plotted as one.
    Cliff's delta rank effect size from the Mann-Whitney statistic, signed so
                  that positive means higher in group2, i.e. the same direction
                  as log2FC.
    rho, rho_adj  Spearman against the ordinal axis, unadjusted and after
                  residualising both sides' ranks on age, BMI, Brinkman index,
                  alcohol and sex. The adjusted one is what a claim rests on.

TWO CONTROLS THE TRENDS WOULD BE MEANINGLESS WITHOUT
----------------------------------------------------
1. CONFOUNDERS. The organisms carrying the rising copies track CRC stage, and so
   do age and smoking in this cohort. All five covariates are complete for all
   347 samples, so the adjusted statistic needs no imputation.
2. A PREVALENCE-MATCHED BACKGROUND. Twelve of 455 target clusters trending is
   only interesting if trending is rare. `background.R` draws five
   non-target community-level gene families per observed target, matched decile
   by decile on prevalence (power in a rank test is almost entirely prevalence),
   and this script runs the identical trend on them. Comparison is on RAW p and
   |rho|, not on BH q, because q depends on family size and the two families
   differ fivefold by construction. Without that file the script still runs and
   records the background as not measured, rather than implying it passed.

THE TWO ROUTES STAY APART
-------------------------
`mapping_validation_summary.json` measured HUMAnN's UniRef90->EC map at 0%
sensitivity for EC 4.1.1.105, 4.1.1.28 and 4.1.1.25 -- the three aromatic amino
acid decarboxylases the chapter's tryptamine and tyramine claims rest on. So the
EC route (17,072 clusters regrouped from the map) and the curated route (224
clusters named by reviewed Swiss-Prot sequence) are tested and drawn separately
and NEVER summed. For those three ECs the EC route is reported as not-measurable
rather than as a near-zero result, because a near-zero result there is a property
of the map, not of the gut.

OUTPUTS
-------
    out/sample_crosswalk.tsv          Sample <-> accession <-> group, the join itself
    out/progression_cluster_tests.tsv per-cluster, per-contrast
    out/progression_enzyme_tests.tsv  per-EC aggregate, per-contrast, per-route
    out/progression_species_tests.tsv per-(cluster, species), per-contrast
    out/progression_trend.tsv         Spearman, unadjusted and adjusted, both axes
    out/background_trend.tsv          the same, for the matched background set
    out/progression_provenance.json   inputs, counts, filters, control results
    out/figs/fig1_source_turnover     the lead: one activity, two directions
    out/figs/fig2_trend_forest        every q<0.05 cluster + the background control
    out/figs/fig3_stage_profiles      mean share by stage for the trending clusters
    out/figs/fig4_tryptamine          the chapter's tryptamine claim, answered
    out/figs/fig5_enzyme_log2fc       per-EC aggregate: flat, in both routes
    out/figs/fig6_cluster_volcano     pairwise contrasts: nothing survives BH

Run from anywhere:   python v2/genes/progression_axis.py
"""

from __future__ import annotations

import json
import platform
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact, mannwhitneyu, spearmanr
from statsmodels.stats.multitest import multipletests

HERE = Path(__file__).resolve().parent          # v2/genes
V2 = HERE.parent                                 # v2
ROOT = V2.parent                                 # repo root
sys.path.insert(0, str(V2))

import utils  # noqa: E402  the pipeline's single source of truth for the axis
import viz    # noqa: E402  the pipeline's only sanctioned figure writer

OUT = HERE / "out"
FIGS = OUT / "figs"
REF = HERE / "ref"

COHORT_META = ROOT / "Data" / "YACHIDA-CRC-2019 metadata.tsv"
CURATION = REF / "YachidaS_2019_sample.tsv"
LONG = OUT / "gf_targets_long.tsv.gz"
BACKGROUND_LONG = OUT / "gf_background_long.tsv.gz"      # optional; see main()
SUMMARY = OUT / "gf_targets_summary.tsv"
ATTRIBUTION = OUT / "species_attribution.tsv"
COLNAMES = OUT / "colnames.txt"
NORMALIZATION = OUT / "sample_normalization.tsv"
TARGETS = HERE / "target_clusters.tsv"
CURATED = HERE / "curated_clusters.tsv"
VALIDATION = HERE / "mapping_validation_summary.json"

DATASET = "YACHIDA-CRC-2019"

BLIND_SPOT_FALLBACK = {"4.1.1.105", "4.1.1.28", "4.1.1.25"}

HEADLINE_CLUSTER = "UniRef90_A7B1V0"             # R. gnavus tryptophan decarboxylase
HEADLINE_SPECIES = "g__Blautia.s__Ruminococcus_gnavus"
ABSENT_CLUSTER = "UniRef90_J7SZ64"               # C. sporogenes homologue; the control

# A feature enters a test only if it is detected often enough for a rank test to
# have anything to rank. 10% of the primary set, and never fewer than 10 samples.
# Features below the line are KEPT IN THE OUTPUT TABLE with tested=False and a
# reason, because "too sparse to test" is a result about the cohort.
MIN_PREV_FRACTION = 0.10
MIN_PREV_SAMPLES = 10

FDR = 0.05

CONFOUNDERS = ["Age", "BMI", "Brinkman.Index", "Alcohol", "Gender_enc"]

FAMILIES = {
    "primary": {
        "contrasts": [("Healthy", "Early_CRC"), ("Healthy", "Advanced_CRC"),
                      ("Early_CRC", "Advanced_CRC")],
    },
    "early": {
        "contrasts": [("Healthy", "Adenoma"), ("Adenoma", "Early_CRC"),
                      ("Adenoma", "Advanced_CRC")],
    },
    "sensitivity": {
        "contrasts": [("Healthy", "PostSurgery")],
    },
}

# Descriptive panels show every group, with the two off-axis ones last so no
# reader mistakes their position for an ordinal one.
DISPLAY_ORDER = ["Healthy", "Adenoma", "Early_CRC", "Advanced_CRC", "PostSurgery"]
PRIMARY_ORDER = ["Healthy", "Early_CRC", "Advanced_CRC"]

# The two activities where the turnover result lives: both are carried by many
# clusters, so both can show components moving against each other.
TURNOVER_ECS = ["4.1.1.19", "2.3.1.1"]


# ===========================================================================
# 1. THE JOIN
# ===========================================================================

def load_crosswalk() -> pd.DataFrame:
    """Cohort sample <-> BioSample accession <-> harmonized group. Asserted."""
    raw = pd.read_csv(COHORT_META, sep="\t", dtype={"Sample": str})
    meta = utils.harmonize_metadata(raw, DATASET, exclude={})
    meta["Sample"] = meta["Sample"].astype(str)

    cur = pd.read_csv(CURATION, sep="\t", dtype=str)
    if cur["subject_id"].duplicated().any():
        raise ValueError(
            "curation table has repeated subject_id, so subject->sample is not a "
            "function and this join is ambiguous")
    cur["cohort_sample"] = cur["subject_id"].str.replace(r"^sub_", "", regex=True)
    bridge = cur[["cohort_sample", "sample_id", "tumor_staging_ajcc", "disease"]]

    cw = meta.merge(bridge, left_on="Sample", right_on="cohort_sample", how="left")
    missing = cw[cw["sample_id"].isna()]["Sample"].tolist()
    if missing:
        raise ValueError(f"{len(missing)} cohort samples have no accession: "
                         f"{missing[:10]}")

    have = set(COLNAMES.read_text().split())
    absent = sorted(set(cw["sample_id"]) - have)
    if absent:
        raise ValueError(f"{len(absent)} accessions are not columns of the gene "
                         f"matrix: {absent[:10]}")

    cw = cw.rename(columns={"sample_id": "Accession"}).drop(columns=["cohort_sample"])
    return cw


# ===========================================================================
# 2. THE GENE TABLES
# ===========================================================================

def load_long(path: Path, accessions: list[str]):
    """Wide matrices (sample x feature) for the three row kinds HUMAnN emits.

    A stratified HUMAnN table carries three kinds of row per gene family, and
    they are NOT independent: the community total equals the sum of its species
    rows plus its `|unclassified` row. `gf_targets_summary.tsv` calls a row
    stratified only when it names a genus, so `|unclassified` rows arrive flagged
    stratified=FALSE and would be double-counted by anything that treats
    "not stratified" as "community total". 1,085 of the 2,352 non-stratified rows
    are `|unclassified`, which is why this splits on the row name itself.
    """
    long = pd.read_csv(path, sep="\t")
    long = long[long["sample"].isin(set(accessions))]

    has_bar = long["row"].str.contains("|", regex=False)
    unclassified = long["row"].str.endswith("|unclassified")
    community = long[~has_bar].copy()
    unclass = long[unclassified].copy()
    species = long[has_bar & ~unclassified].copy()
    if len(community) + len(unclass) + len(species) != len(long):
        raise ValueError("row kinds do not partition the long table")

    comm_wide = (community.pivot_table(index="sample", columns="row",
                                       values="abundance", fill_value=0.0)
                 .reindex(accessions, fill_value=0.0))
    if len(species):
        species[["uniref90", "species"]] = species["row"].str.split("|", n=1,
                                                                    expand=True)
        sp_wide = (species.pivot_table(index="sample", columns="row",
                                       values="abundance", fill_value=0.0)
                   .reindex(accessions, fill_value=0.0))
    else:
        sp_wide = pd.DataFrame(index=pd.Index(accessions, name="sample"))

    return comm_wide, sp_wide, species, dict(
        n_long=len(long), n_community=len(community),
        n_unclassified=len(unclass), n_species=len(species))


def check_normalization() -> dict:
    """Settle whether the stored values are already compositional.

    Every comparison below is between groups of samples, so if the values were
    RPKs the column totals would track sequencing depth and a between-group
    difference could be depth alone. `colsums.R` measures the community-level
    column sums of the full 8,930,907-row matrix on Dardel; they came back at
    1.000000 for all 616 samples, so these are relative abundances and no
    renormalization applies. When that file is absent the script still runs, but
    the provenance records normalization as UNVERIFIED and the figures say so,
    rather than quietly assuming the convenient case.
    """
    if not NORMALIZATION.exists():
        return {"status": "unverified", "reason": f"{NORMALIZATION.name} not present",
                "note": "run colsums.R on Dardel to settle this"}
    nz = pd.read_csv(NORMALIZATION, sep="\t")
    s = nz["sum_community"]
    compositional = bool(np.all(np.abs(s - 1.0) < 1e-4))
    return {
        "status": "relative_abundance" if compositional else "not_compositional",
        "community_sum_min": float(s.min()), "community_sum_max": float(s.max()),
        "stratified_share_median": float(nz["sum_stratified"].median()),
        "stratified_share_min": float(nz["sum_stratified"].min()),
        "nonzero_community_median": float(nz["nonzero_community"].median()),
        "n_samples": int(len(nz)),
        "renormalized": not compositional,
    }


# ===========================================================================
# 3. FEATURE SETS AND ANNOTATION
# ===========================================================================

def blind_spot_ecs() -> set[str]:
    """ECs the EC route demonstrably cannot see, read from the measurement."""
    if not VALIDATION.exists():
        return set(BLIND_SPOT_FALLBACK)
    v = json.loads(VALIDATION.read_text())
    per_ec = v.get("per_ec") or v.get("by_ec") or {}
    blind = set()
    for ec, rec in per_ec.items():
        if isinstance(rec, dict):
            n = rec.get("n_curated", rec.get("n", 0)) or 0
            hit = rec.get("n_recovered", rec.get("recovered", 0)) or 0
            if n > 0 and hit == 0:
                blind.add(str(ec))
    return blind or set(BLIND_SPOT_FALLBACK)


def cluster_annotation() -> pd.DataFrame:
    """One row per observed cluster: its ECs, its activity, its routes, its
    dominant carrier species."""
    summ = pd.read_csv(SUMMARY, sep="\t")
    obs = (summ[~summ["row"].str.contains("|", regex=False)]
           [["uniref90", "in_ec_route", "in_curated_route", "prevalence"]]
           .drop_duplicates("uniref90").set_index("uniref90"))

    tgt = pd.read_csv(TARGETS, sep="\t")
    ec_label = tgt.drop_duplicates("ec").set_index("ec")[["activity", "claim"]]
    ec_of_cluster = (tgt.groupby("uniref90")["ec"]
                     .agg(lambda s: sorted(set(s))).to_dict())

    cur = pd.read_csv(CURATED, sep="\t")
    cur_ec = {r.uniref90: sorted({e for e in str(r.ecs).split(";")
                                  if e and e != "nan"})
              for r in cur.itertuples()}
    cur_name = dict(zip(cur["uniref90"], cur["cluster_name"]))

    # Dominant carrier: the species-stratified row of this cluster with the
    # highest prevalence. `species_attribution.tsv` already holds that join.
    top_sp = {}
    if ATTRIBUTION.exists():
        sa = pd.read_csv(ATTRIBUTION, sep="\t")
        for cl, g in sa.groupby("uniref90"):
            r = g.loc[g["prevalence"].idxmax()]
            top_sp[cl] = (str(r["species"]), float(r["prevalence"]))

    rows = []
    for cl in obs.index:
        sp, sp_prev = top_sp.get(cl, ("", np.nan))
        rows.append({
            "uniref90": cl,
            "in_ec_route": bool(obs.at[cl, "in_ec_route"]),
            "in_curated_route": bool(obs.at[cl, "in_curated_route"]),
            "ec_route_ecs": ec_of_cluster.get(cl, []),
            "curated_ecs": cur_ec.get(cl, []),
            "cluster_name": cur_name.get(cl, ""),
            "top_species": sp,
            "top_species_prevalence": sp_prev,
        })
    ann = pd.DataFrame(rows).set_index("uniref90")
    ann["all_ecs"] = [sorted(set(a) | set(b))
                      for a, b in zip(ann["ec_route_ecs"], ann["curated_ecs"])]
    ann.attrs["ec_label"] = ec_label
    return ann


def enzyme_aggregates(comm_wide: pd.DataFrame, ann: pd.DataFrame):
    """Per-EC, per-route aggregate: the summed community share of that activity.

    Summing community-level relative abundances over the clusters annotated to an
    EC is a share of the gene pool attributable to that activity, which is the
    quantity a metabolite claim needs. A cluster annotated to two ECs (the
    N-acetylglutamate synthases carry 2.3.1.1 and 2.3.1.35 together in 759 of the
    observed rows) is counted under BOTH, so the per-EC columns overlap by
    construction and must not be added to each other. The overlap is recorded.
    """
    feats, meta = {}, []
    overlap = 0
    for route, col, ec_col in (("ec", "in_ec_route", "ec_route_ecs"),
                               ("curated", "in_curated_route", "curated_ecs")):
        members = ann[ann[col]]
        ec_to_clusters: dict[str, list[str]] = {}
        for cl, ecs in members[ec_col].items():
            if len(ecs) > 1:
                overlap += 1
            for ec in ecs:
                ec_to_clusters.setdefault(ec, []).append(cl)
        for ec, cls in sorted(ec_to_clusters.items()):
            cols = [c for c in cls if c in comm_wide.columns]
            if not cols:
                continue
            name = f"{route}|{ec}"
            feats[name] = comm_wide[cols].sum(axis=1)
            meta.append({"feature": name, "route": route, "ec": ec,
                         "n_clusters": len(cols)})
    return pd.DataFrame(feats), pd.DataFrame(meta), overlap


# ===========================================================================
# 4. TESTS
# ===========================================================================

def _one_contrast(x: pd.Series, y: pd.Series) -> dict:
    """Mann-Whitney, Cliff's delta, log2 ratio of means, Fisher on detection."""
    n1, n2 = len(x), len(y)
    nz1, nz2 = int((x > 0).sum()), int((y > 0).sum())
    m1, m2 = float(x.mean()), float(y.mean())

    out = {"n1": n1, "n2": n2, "nonzero1": nz1, "nonzero2": nz2,
           "prevalence1": nz1 / n1 if n1 else np.nan,
           "prevalence2": nz2 / n2 if n2 else np.nan,
           "mean1": m1, "mean2": m2}
    out["log2FC"] = float(np.log2(m2 / m1)) if (m1 > 0 and m2 > 0) else np.nan
    out["direction"] = "up" if m2 > m1 else "down" if m2 < m1 else "flat"

    try:
        stat, p = mannwhitneyu(x, y, alternative="two-sided")
        out["mw_p"] = float(p)
        # Sign convention: POSITIVE means higher in group2, so the delta agrees with
        # log2FC. scipy's U is computed for the first argument, so the textbook
        # 2U/(n1*n2) - 1 is positive when group1 is the larger -- the opposite of
        # what `utils.differential_abundance` reports beside its own MeanDiff, where
        # a reader comparing the two columns reads every direction backwards. Only
        # the sign differs; the magnitude is the same statistic.
        out["cliffs_delta"] = float(1.0 - 2.0 * stat / (n1 * n2))
    except ValueError:
        out["mw_p"], out["cliffs_delta"] = np.nan, np.nan

    try:
        _, fp = fisher_exact([[nz1, n1 - nz1], [nz2, n2 - nz2]])
        out["fisher_p"] = float(fp)
    except ValueError:
        out["fisher_p"] = np.nan
    return out


def run_family(wide: pd.DataFrame, cw: pd.DataFrame, family: str,
               min_prev_n: int) -> pd.DataFrame:
    """Every contrast in one family, over every feature, BH-corrected together.

    One family = one BH correction over (features x contrasts). Correcting per
    contrast instead would make the reported q depend on how many contrasts a
    caller happened to request, which is the mistake `COMPARISONS_EARLY`'s
    "own FDR family, never merged" note exists to prevent.
    """
    groups = cw.set_index("Accession")["Stage.Group"]
    rows = []
    for g1, g2 in FAMILIES[family]["contrasts"]:
        s1 = [a for a in wide.index if groups.get(a) == g1]
        s2 = [a for a in wide.index if groups.get(a) == g2]
        if len(s1) < 3 or len(s2) < 3:
            continue
        sub = wide.loc[s1 + s2]
        detect = (sub > 0).sum(axis=0)
        for feat in wide.columns:
            rec = {"family": family, "feature": feat, "group1": g1, "group2": g2}
            if detect[feat] < min_prev_n:
                rec.update({"tested": False,
                            "skip_reason": f"detected in {int(detect[feat])} of "
                                           f"{len(sub)} samples, below the "
                                           f"{min_prev_n}-sample floor"})
                rows.append(rec)
                continue
            rec["tested"] = True
            rec["skip_reason"] = ""
            rec.update(_one_contrast(wide.loc[s1, feat], wide.loc[s2, feat]))
            rows.append(rec)

    res = pd.DataFrame(rows)
    if res.empty:
        return res
    res["mw_q"] = np.nan
    res["fisher_q"] = np.nan
    ok = res["tested"] & res["mw_p"].notna()
    if ok.any():
        res.loc[ok, "mw_q"] = multipletests(res.loc[ok, "mw_p"], method="fdr_bh")[1]
    okf = res["tested"] & res["fisher_p"].notna()
    if okf.any():
        res.loc[okf, "fisher_q"] = multipletests(res.loc[okf, "fisher_p"],
                                                 method="fdr_bh")[1]
    res["n_tests_in_family"] = int(ok.sum())
    return res


def _fisher_ci(rho: float, n: int, k: int = 0, z: float = 1.959964):
    """Fisher z interval for a (partial) rank correlation; n - 3 - k on the se."""
    if not np.isfinite(rho) or n - 3 - k <= 0 or abs(rho) >= 1:
        return (np.nan, np.nan)
    se = 1.0 / np.sqrt(n - 3 - k)
    zr = np.arctanh(rho)
    return (float(np.tanh(zr - z * se)), float(np.tanh(zr + z * se)))


def run_trend(wide: pd.DataFrame, cw: pd.DataFrame, min_prev_n: int,
              label: str = "target") -> pd.DataFrame:
    """Spearman against the ordinal axis, unadjusted and confounder-adjusted.

    PostSurgery and Adenoma are NaN on Stage.Group.Numeric by pipeline policy, so
    dropping NaN here IS the primary mask -- a resected gut and a pre-malignant
    adenoma cannot enter a severity trend even if this function is called
    carelessly. The early axis gets its own trend, over EARLY_AXIS_NUMERIC.

    WHY AN ADJUSTED COLUMN IS NOT OPTIONAL HERE
    The organisms that turn out to carry the rising copies -- Fusobacterium
    nucleatum, Parvimonas micra, Bilophila wadsworthia -- track CRC stage, and so
    do age and smoking in this cohort (Brinkman index runs to 3,200). An
    unadjusted rank correlation between a gene's share and stage therefore cannot
    distinguish a microbial shift from an age-and-smoking shift.
    utils.partial_spearman_grid residualises the ranks of both sides on the ranks
    of the five covariates, all of them complete for all 347 Yachida samples, so
    there is nothing to impute. Unadjusted and adjusted are reported side by
    side; the adjusted one is the one a claim rests on.
    """
    g = cw.set_index("Accession")
    Zfull = utils.build_confounder_matrix(cw)
    Zfull.index = cw["Accession"].values
    out = []
    for axis, mapping in (("primary", utils.STAGE_GROUP_NUMERIC),
                          ("early", utils.EARLY_AXIS_NUMERIC)):
        score = g["Stage.Group"].map(mapping)
        keep = [a for a in wide.index if a in score.index and pd.notna(score.get(a))]
        y = score.loc[keep].astype(float)
        sub = wide.loc[keep]
        detect = (sub > 0).sum(axis=0)
        feats = [f for f in wide.columns if detect[f] >= min_prev_n]
        if not feats:
            continue

        recs = [dict(zip(("set", "axis", "feature", "n"), (label, axis, f, len(keep))))
                for f in feats]
        df = pd.DataFrame(recs)
        rp = [spearmanr(sub[f].values, y.values) for f in feats]
        df["rho"] = [float(r) for r, _ in rp]
        df["p"] = [float(p) for _, p in rp]

        Z = Zfull.loc[keep, CONFOUNDERS].astype(float)
        if Z.isna().any().any():
            raise ValueError("confounder matrix has missing values; "
                             "partial_spearman_grid needs complete cases")
        r_adj, p_adj = utils.partial_spearman_grid(
            sub[feats].to_numpy(dtype=float), y.to_numpy(dtype=float)[:, None],
            Z.to_numpy())
        df["rho_adj"] = r_adj[:, 0]
        df["p_adj"] = p_adj[:, 0]

        df["q"] = multipletests(df["p"], method="fdr_bh")[1]
        df["q_adj"] = np.nan
        ok = df["p_adj"].notna()
        if ok.any():
            df.loc[ok, "q_adj"] = multipletests(df.loc[ok, "p_adj"],
                                                method="fdr_bh")[1]
        ci = [_fisher_ci(r, len(keep)) for r in df["rho"]]
        df["rho_lo"], df["rho_hi"] = [c[0] for c in ci], [c[1] for c in ci]
        ci_a = [_fisher_ci(r, len(keep), k=len(CONFOUNDERS)) for r in df["rho_adj"]]
        df["rho_adj_lo"] = [c[0] for c in ci_a]
        df["rho_adj_hi"] = [c[1] for c in ci_a]
        out.append(df)
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()


TAIL_THRESHOLDS = [0.05, 0.01, 1e-3, 1e-4]


def compare_to_background(tgt: pd.DataFrame, bg: pd.DataFrame,
                          axis: str = "primary") -> dict:
    """Is trending with stage rare? Compared on raw p, not on q, and at several
    thresholds rather than one.

    BH's q depends on how many tests are in the family, and the background family
    is five times larger by design, so a q-rate comparison would be comparing two
    different corrections. Raw p means the same thing in both sets.

    ONE THRESHOLD WOULD HAVE GIVEN THE WRONG ANSWER. At p < 0.05 the target
    enzymes trend LESS often than matched background gene families (8.6% against
    13.3%), which read alone says the chapter's pathways are nothing special. At
    p < 1e-4 the same comparison runs the other way (0.66% against 0.14%), and
    the single strongest target signal is two orders of magnitude beyond anything
    in a five-times-larger background. Both are true and they are different
    statements: the bulk of these redundant core activities is more stable than
    the average gut gene family, while a few specific copies move hard. Reporting
    either alone would be a choice of conclusion, so the whole profile is
    recorded and drawn.
    """
    t = tgt[tgt["axis"] == axis]
    b = bg[bg["axis"] == axis]
    if t.empty or b.empty:
        return {"status": "not measured"}
    out = {"status": "measured", "axis": axis,
           "n_target": int(len(t)), "n_background": int(len(b)),
           "thresholds": TAIL_THRESHOLDS}
    for tag, pcol, rcol in (("unadjusted", "p", "rho"),
                            ("adjusted", "p_adj", "rho_adj")):
        tp, bp = t[pcol].dropna(), b[pcol].dropna()
        tr, br = t[rcol].abs().dropna(), b[rcol].abs().dropna()
        tail = []
        for thr in TAIL_THRESHOLDS:
            kt, kb = int((tp < thr).sum()), int((bp < thr).sum())
            rt = kt / len(tp) if len(tp) else np.nan
            rb = kb / len(bp) if len(bp) else np.nan
            try:
                _, fp = fisher_exact([[kt, len(tp) - kt], [kb, len(bp) - kb]])
            except ValueError:
                fp = np.nan
            tail.append({"threshold": thr, "target_n": kt, "target_rate": rt,
                         "background_n": kb, "background_rate": rb,
                         "ratio": float(rt / rb) if rb else None,
                         "fisher_p": float(fp) if np.isfinite(fp) else None})
        try:
            _, mw_p = mannwhitneyu(tr, br, alternative="greater")
        except ValueError:
            mw_p = np.nan
        out[tag] = {
            "tail": tail,
            "min_p_target": float(tp.min()) if len(tp) else None,
            "min_p_background": float(bp.min()) if len(bp) else None,
            "median_abs_rho_target": float(tr.median()) if len(tr) else np.nan,
            "median_abs_rho_background": float(br.median()) if len(br) else np.nan,
            "mannwhitney_p_abs_rho_target_greater": (float(mw_p)
                                                     if np.isfinite(mw_p) else None),
        }
    return out


# ---------------------------------------------------------------------------
# The carrier-species control
# ---------------------------------------------------------------------------

SPECIES_TABLE = ROOT / "Data" / "YACHIDA-CRC-2019 species.tsv"
_GTDB_SPLIT = None   # regex built once in carrier_abundance()


def _match_gtdb(humann_species: str, epithets: dict) -> list[str]:
    """GTDB columns that correspond to a HUMAnN species name.

    The two tables use different taxonomies. HUMAnN's stratification carries
    ChocoPhlAn names (`g__Fusobacterium.s__Fusobacterium_nucleatum`); the cohort's
    species table carries GTDB lineages (`...;s__Fusobacterium nucleatum`). Beyond
    spaces-for-underscores, GTDB splits several species into alphabetically
    suffixed daughters -- F. nucleatum is nucleatum, _D, _E and _J here -- so a
    single HUMAnN species can correspond to several GTDB columns and the honest
    carrier abundance is their sum. A HUMAnN species with no GTDB counterpart at
    all (Roseburia faecis, which GTDB reclassified) returns an empty list, and the
    control is then reported as unavailable for that cluster rather than run
    against the wrong organism.
    """
    name = humann_species.split(".s__")[-1].strip()
    if not name:
        return []
    hits = list(epithets.get(name, []))
    for ep, cols in epithets.items():
        if ep.startswith(name + "_") and len(ep) == len(name) + 2 and ep[-1].isupper():
            hits.extend(cols)
    return hits


def carrier_abundance(species_needed: set[str], cw: pd.DataFrame):
    """{HUMAnN species -> per-accession relative abundance} from the cohort's own
    GTDB species table, plus a record of what could not be matched."""
    if not SPECIES_TABLE.exists():
        return {}, {"status": "unavailable",
                    "reason": f"{SPECIES_TABLE.name} not present"}
    import re
    header = pd.read_csv(SPECIES_TABLE, sep="\t", nrows=0).columns.tolist()
    epithets: dict[str, list[str]] = {}
    for c in header:
        m = re.search(r"s__(.+)$", c)
        if m:
            epithets.setdefault(m.group(1).strip().replace(" ", "_"), []).append(c)

    plan = {s: _match_gtdb(s, epithets) for s in sorted(species_needed) if s}
    needed_cols = sorted({c for cols in plan.values() for c in cols})
    if not needed_cols:
        return {}, {"status": "no_match", "requested": len(plan)}

    df = pd.read_csv(SPECIES_TABLE, sep="\t",
                     usecols=["Sample"] + needed_cols,
                     dtype={"Sample": str})
    key = cw.set_index("Sample")["Accession"]
    df["Accession"] = df["Sample"].astype(str).map(key)
    df = df.dropna(subset=["Accession"]).set_index("Accession")

    out = {}
    for s, cols in plan.items():
        if cols:
            out[s] = df[cols].sum(axis=1)
    report = {"status": "measured",
              "requested": len(plan),
              "matched": len(out),
              "unmatched": sorted(s for s, c in plan.items() if not c),
              "gtdb_columns_used": {s: len(c) for s, c in plan.items() if c},
              "samples_joined": int(len(df))}
    return out, report


def carrier_adjusted_trend(comm_wide: pd.DataFrame, cw: pd.DataFrame,
                           features: list[str], ann: pd.DataFrame,
                           carriers: dict) -> pd.DataFrame:
    """Does a cluster still trend once its carrier organism's abundance is held
    fixed?

    THE QUESTION A REVIEWER ASKS FIRST. The copies that rise are carried by
    Fusobacterium nucleatum, Parvimonas micra and Bilophila wadsworthia, and
    those organisms are themselves the best-known CRC-enriched taxa. A gene's
    share of the gene pool is roughly its carrier's abundance times its share
    within that carrier, so a rising gene is the expected consequence of a rising
    organism and is not independent evidence of anything. Adding the carrier's own
    relative abundance -- measured separately, from the cohort's GTDB species
    table, not from this gene matrix -- to the covariates answers which it is.

    A trend that survives means the gene's share moves beyond what its carrier's
    abundance explains. A trend that vanishes means the finding is ATTRIBUTION:
    this organism, carrying this activity, increases. Attribution is what the
    chapter needs and is worth stating plainly; what it is not is an independent
    gene-level signal, and this column is what stops it being written up as one.
    """
    g = cw.set_index("Accession")
    Zfull = utils.build_confounder_matrix(cw)
    Zfull.index = cw["Accession"].values
    score = g["Stage.Group"].map(utils.STAGE_GROUP_NUMERIC)
    keep = [a for a in comm_wide.index
            if a in score.index and pd.notna(score.get(a))]
    y = score.loc[keep].astype(float).to_numpy()[:, None]
    Z0 = Zfull.loc[keep, CONFOUNDERS].astype(float).to_numpy()

    rows = []
    for feat in features:
        sp = ann.at[feat, "top_species"] if feat in ann.index else ""
        rec = {"feature": feat, "carrier_species": sp}
        if not sp or sp not in carriers:
            rec.update({"carrier_available": False,
                        "reason": ("no species-stratified row for this cluster"
                                   if not sp else
                                   "carrier has no GTDB counterpart in the "
                                   "cohort species table")})
            rows.append(rec)
            continue
        c = carriers[sp].reindex(keep)
        if c.isna().any():
            rec.update({"carrier_available": False,
                        "reason": "carrier abundance missing for some samples"})
            rows.append(rec)
            continue
        Z = np.column_stack([Z0, c.to_numpy(dtype=float)])
        x = comm_wide.loc[keep, [feat]].to_numpy(dtype=float)
        r, p = utils.partial_spearman_grid(x, y, Z)
        rc, pc = spearmanr(c.to_numpy(dtype=float), y[:, 0])
        rec.update({"carrier_available": True, "reason": "",
                    "n": len(keep),
                    "rho_carrier_adj": float(r[0, 0]),
                    "p_carrier_adj": float(p[0, 0]),
                    "carrier_rho_vs_stage": float(rc),
                    "carrier_p_vs_stage": float(pc)})
        rows.append(rec)
    df = pd.DataFrame(rows)
    ok = df.get("p_carrier_adj", pd.Series(dtype=float)).notna()
    df["q_carrier_adj"] = np.nan
    if ok.any():
        df.loc[ok, "q_carrier_adj"] = multipletests(df.loc[ok, "p_carrier_adj"],
                                                    method="fdr_bh")[1]
    return df


# ===========================================================================
# 5. FIGURE HELPERS
# ===========================================================================

def _wilson(k: int, n: int, z: float = 1.959964):
    """Wilson score interval. The normal approximation is wrong at p near 0,
    which is where a prevalence panel spends most of its time."""
    if n == 0:
        return (np.nan, np.nan)
    p = k / n
    d = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / d
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, centre - half), min(1.0, centre + half))


def _groups_present(cw: pd.DataFrame) -> list[str]:
    return [g for g in DISPLAY_ORDER if (cw["Stage.Group"] == g).any()]


def _short_species(s: str) -> str:
    if not s or s == "other species":
        return s or "unattributed"
    name = s.split(".s__")[-1].replace("_", " ")
    parts = name.split()
    return f"{parts[0][0]}. {' '.join(parts[1:])}" if len(parts) > 1 else name


def _label(cl: str, ann: pd.DataFrame) -> str:
    """'EC  activity  --  carrier species' for a cluster, for a figure axis."""
    ec_label = ann.attrs["ec_label"]
    ecs = ann.at[cl, "all_ecs"] if cl in ann.index else []
    ec = ecs[0] if ecs else ""
    act = ec_label.at[ec, "activity"] if ec in ec_label.index else ""
    sp = _short_species(ann.at[cl, "top_species"] if cl in ann.index else "")
    return f"{cl.replace('UniRef90_', '')}  {ec} {act}\n      {sp}"


def _caption(fig, text: str, width: int = 155) -> None:
    """Put a long methods caption under a figure without it landing on the axes.

    `fig.text(..., y=-0.03)` does not participate in constrained_layout, so the
    caption was drawn over the bottom panels' tick labels at whatever y was
    guessed. `supxlabel` IS laid out, so the space is reserved; it does not wrap
    reliably, so the wrapping is done here.
    """
    import textwrap
    fig.supxlabel("\n".join(textwrap.wrap(" ".join(text.split()), width)),
                  fontsize=7.6, color=viz.MUTED, x=0.005, ha="left")


def _place_labels(ax, items, x_frac: float = 0.50, min_gap: float = 0.072):
    """Right-hand series labels, spread so they cannot overlap.

    `items` is [(y_data, text, colour)]. Each label is anchored to its line's last
    point by a leader, but sits in a slot assigned in axes fraction, because two
    clusters whose final values differ by 1% of the axis otherwise print on top of
    each other -- which is what the first draft of this figure did to the two
    Lachnospira and Roseburia lines.
    """
    lo, hi = ax.get_ylim()
    span = (hi - lo) or 1.0
    rows = sorted(((y, t, c) for y, t, c in items if np.isfinite(y)),
                  key=lambda r: -r[0])
    if not rows:
        return
    n = len(rows)
    min_gap = min(min_gap, 1.0 / n)
    slots, last = [], 1.0 + min_gap
    for y, _, _ in rows:
        f = min((y - lo) / span, last - min_gap)
        slots.append(f)
        last = f
    # The greedy pass walks downwards and can run off the bottom when several
    # lines converge, which is how the first draft printed two Lachnospira labels
    # on top of each other. Lift the whole stack by the overshoot if there is room
    # above it, and fall back to even spacing when there is not.
    if min(slots) < 0.0:
        headroom = 1.0 - max(slots)
        shift = min(-min(slots), max(headroom, 0.0))
        slots = [s + shift for s in slots]
        if min(slots) < 0.0:
            slots = list(np.linspace(1.0, 0.0, n))
    anchor_x = max((ln.get_xdata()[-1] for ln in ax.get_lines()
                    if len(ln.get_xdata())), default=ax.get_xlim()[1])
    for (y, text, col), f in zip(rows, slots):
        ax.annotate(text, xy=(anchor_x, y), xycoords="data",
                    xytext=(x_frac, f), textcoords="axes fraction",
                    fontsize=6.8, color=col, va="center", ha="left",
                    annotation_clip=False,
                    arrowprops=dict(arrowstyle="-", color=col, lw=0.6,
                                    alpha=0.5, shrinkA=0, shrinkB=2))


def _group_means(values: pd.Series, cw: pd.DataFrame, order):
    """Mean and standard error of the mean per group, over ALL samples in the
    group: an undetected sample contributes a zero, not a missing value."""
    g = cw.set_index("Accession")["Stage.Group"]
    m, se, n = [], [], []
    for grp in order:
        s = values[[a for a in values.index if g.get(a) == grp]]
        m.append(float(s.mean()) if len(s) else np.nan)
        se.append(float(s.std(ddof=1) / np.sqrt(len(s))) if len(s) > 1 else np.nan)
        n.append(len(s))
    return np.array(m), np.array(se), n


# ===========================================================================
# 6. FIGURES
# ===========================================================================

def fig1_source_turnover(comm_wide, enz_wide, cw, trend, ann, carrier, blind):
    """The lead figure: one activity, two directions, two sets of organisms.

    Top row, per activity: every cluster that trends at q<0.05, as mean share by
    stage, normalised to its own Healthy mean so clusters four orders of magnitude
    apart can share an axis and only the SHAPE is compared. Bottom row: the EC
    aggregate those same clusters belong to, in absolute share, which does not
    move. The figure's claim is the contrast between the two rows.

    Each line carries three statistics, and the third is the one that limits the
    claim: rho unadjusted, rho after the five clinical covariates, and rho after
    the carrier organism's own abundance as well. Where the third collapses, the
    line is the organism arriving, not the gene changing, and the figure says so
    on the line rather than in a caveat elsewhere.
    """
    ec_label = ann.attrs["ec_label"]
    tr = trend[(trend["axis"] == "primary") & (trend["set"] == "target")]
    hits = tr[tr["q"] < FDR].set_index("feature")
    if hits.empty:
        return
    carr = (carrier.set_index("feature") if carrier is not None and not carrier.empty
            else pd.DataFrame())

    ecs = [e for e in TURNOVER_ECS
           if any(e in (ann.at[c, "all_ecs"] if c in ann.index else [])
                  for c in hits.index)]
    if not ecs:
        return

    fig, axes = plt.subplots(2, len(ecs), figsize=(8.6 * len(ecs), 8.6),
                             constrained_layout=True, squeeze=False,
                             gridspec_kw={"width_ratios": [1] * len(ecs)})

    for j, ec in enumerate(ecs):
        members = [c for c in hits.index
                   if c in ann.index and ec in ann.at[c, "all_ecs"]]
        ax = axes[0][j]
        labels = []
        for cl in sorted(members, key=lambda c: -hits.at[c, "rho"]):
            m, se, _ = _group_means(comm_wide[cl], cw, PRIMARY_ORDER)
            base = m[0] if m[0] > 0 else np.nan
            rho = hits.at[cl, "rho"]
            rho_a = hits.at[cl, "rho_adj"]
            col = viz.sign_color(rho)
            # log2 of the group mean over the Healthy mean: the fold change the
            # chapter asks for, on a scale where up and down are symmetric and a
            # cluster that moves 11-fold does not flatten five that move 2-fold.
            # The SEM is carried across by the delta method, d(log2 m) = se/(m ln2).
            fc = np.log2(m / base)
            dfc = se / (np.where(m > 0, m, np.nan) * np.log(2))
            ax.plot(range(len(PRIMARY_ORDER)), fc, marker="o", color=col,
                    lw=1.8, ms=5.5, alpha=0.9)
            ax.fill_between(range(len(PRIMARY_ORDER)), fc - dfc, fc + dfc,
                            color=col, alpha=0.10, linewidth=0)
            extra = ""
            if cl in carr.index and bool(carr.at[cl, "carrier_available"]):
                extra = (f", carrier-adj {carr.at[cl, 'rho_carrier_adj']:+.2f}"
                         f"{'' if carr.at[cl, 'p_carrier_adj'] < 0.05 else ' ns'}")
            labels.append((fc[-1],
                           f"{_short_species(ann.at[cl, 'top_species'])}  "
                           f"$\\rho$ {rho:+.2f} (adj {rho_a:+.2f}{extra})", col))
        ax.axhline(0.0, color=viz.GRID, lw=0.9, ls=(0, (3, 3)))
        ax.set_xlim(-0.25, len(PRIMARY_ORDER) - 1 + 2.6)
        _place_labels(ax, labels)
        ax.set_xticks(range(len(PRIMARY_ORDER)))
        ax.set_xticklabels([g.replace("_", " ") for g in PRIMARY_ORDER], fontsize=8.5)
        ax.set_ylabel("log2 fold change of mean share, vs Healthy")
        act = ec_label.at[ec, "activity"] if ec in ec_label.index else ""
        claim = ec_label.at[ec, "claim"] if ec in ec_label.index else ""
        ax.set_title(f"EC {ec}  {act}\n{len(members)} clusters trending at "
                     f"q<{FDR:g}  →  {claim}", fontsize=9.5)
        viz.despine(ax)

        ax2 = axes[1][j]
        for route, ls in (("ec", "-"), ("curated", (0, (4, 2)))):
            feat = f"{route}|{ec}"
            if feat not in enz_wide.columns:
                continue
            m, se, n = _group_means(enz_wide[feat], cw, PRIMARY_ORDER)
            row = tr[tr["feature"] == feat]
            tag = (f", $\\rho$ {row['rho'].iloc[0]:+.2f} q={row['q'].iloc[0]:.2f}"
                   if len(row) else "")
            ax2.errorbar(range(len(PRIMARY_ORDER)), m, yerr=se, marker="s",
                         ls=ls, lw=1.8, ms=5.5, capsize=3,
                         label=f"{route} route{tag}")
        ax2.set_xticks(range(len(PRIMARY_ORDER)))
        ax2.set_xticklabels([g.replace("_", " ") for g in PRIMARY_ORDER], fontsize=8.5)
        ax2.set_ylabel("summed community share of the activity")
        ax2.set_title(f"the same activity in total: EC {ec} aggregate", fontsize=9.5)
        ax2.legend(fontsize=7.5, loc="best")
        ax2.set_ylim(bottom=0)
        viz.despine(ax2)

    fig.suptitle("The source of a metabolic activity turns over while the total "
                 "does not", fontsize=13)
    _caption(fig, "Top: each line is one UniRef90 cluster whose share trends "
             f"monotonically with stage at BH q < {FDR:g} (Spearman over the n=277 "
             "primary set), as log2 of its group mean over its own Healthy mean, so "
             "clusters of very different absolute size are comparable and up and "
             "down are symmetric; shading is ± 1 SEM carried across by the "
             "delta method. The species named carries the highest-prevalence "
             "stratified copy of that cluster. $\\rho$ (adj) is after residualising "
             "on age, BMI, Brinkman index, alcohol and sex; carrier-adj additionally "
             "holds that organism's own relative abundance fixed, measured from the "
             "cohort's GTDB species table rather than from this gene matrix. Where "
             "carrier-adj is marked ns, the gene's trend is what its carrier's "
             "arrival already implies -- the finding is attribution of the activity "
             "to that organism, not an independent change in the gene. Bottom: the "
             "summed share of every observed cluster annotated to the same EC, on an "
             "absolute scale. Clusters move in both directions and the sum does not, "
             "so what changes along this axis is the source of the capacity rather "
             "than its amount.")
    viz.save_figure(fig, FIGS / "fig1_source_turnover.png", close=True)


def fig2_trend_forest(trend, ann, bgcmp: dict):
    """Every cluster that trends, with its confounder-adjusted twin, and the
    background control that says whether trending is rare."""
    tr = trend[(trend["axis"] == "primary") & (trend["set"] == "target")]
    hits = tr[tr["q"] < FDR].sort_values("rho")
    if hits.empty:
        return
    has_bg = bgcmp.get("status") == "measured"

    fig, axes = plt.subplots(1, 2, figsize=(14.2, max(4.4, 0.52 * len(hits) + 2.4)),
                             constrained_layout=True,
                             gridspec_kw={"width_ratios": [1.45, 1.0]})
    ax = axes[0]
    y = np.arange(len(hits))
    up_c, dn_c = viz.DIVERGING[-1], viz.DIVERGING[0]
    for i, r in enumerate(hits.itertuples()):
        col = up_c if r.rho > 0 else dn_c
        ax.plot([r.rho_lo, r.rho_hi], [i, i], color=col, lw=2.2, alpha=0.55,
                solid_capstyle="round")
        ax.plot(r.rho, i, "o", color=col, ms=7.5, zorder=3)
        ax.plot([r.rho_adj_lo, r.rho_adj_hi], [i - 0.26, i - 0.26], color=col,
                lw=1.2, alpha=0.4)
        ax.plot(r.rho_adj, i - 0.26, "o", mfc="none", mec=col, ms=6, mew=1.4,
                zorder=3)
        mark = "" if r.q_adj < FDR else "   adj. q = %.2f" % r.q_adj
        ax.text(0.995, i, f"q {r.q:.1e}{mark}", transform=ax.get_yaxis_transform(),
                ha="right", va="center", fontsize=6.8, color=viz.MUTED)
    ax.axvline(0, color=viz.INK2, lw=1.0)
    # Reserve a right-hand gutter for the q labels; without it they print over the
    # confidence bars of whichever rows happen to reach furthest right.
    x0, x1 = ax.get_xlim()
    ax.set_xlim(x0, x1 + 0.55 * (x1 - x0))
    ax.set_yticks(y)
    ax.set_yticklabels([_label(f, ann) for f in hits["feature"]], fontsize=7.2)
    ax.set_ylim(-0.9, len(hits) - 0.3)
    ax.set_xlabel("Spearman $\\rho$ against stage  (filled: unadjusted, "
                  "open: confounder-adjusted)")
    ax.set_title(f"{len(hits)} of {len(tr)} tested features trend with stage at "
                 f"BH q < {FDR:g}", fontsize=10)
    viz.despine(ax, keep=("bottom",))
    ax.grid(axis="y", visible=False)

    ax = axes[1]
    if has_bg:
        a = bgcmp["adjusted"]
        thr = np.array([row["threshold"] for row in a["tail"]], dtype=float)
        rt = np.array([row["target_rate"] for row in a["tail"]], dtype=float)
        rb = np.array([row["background_rate"] for row in a["tail"]], dtype=float)
        x = np.arange(len(thr))
        ax.plot(x, np.maximum(rt, 1e-4), marker="o", ms=7, lw=2.0,
                color=viz.DIVERGING[-1],
                label=f"target enzymes (n={bgcmp['n_target']})")
        ax.plot(x, np.maximum(rb, 1e-4), marker="s", ms=6.5, lw=2.0,
                color=viz.MUTED,
                label=f"matched background (n={bgcmp['n_background']})")
        ax.plot(x, thr, lw=1.2, ls=(0, (3, 3)), color=viz.INK2,
                label="expected under the null")
        for xi, row in enumerate(a["tail"]):
            if row["ratio"] is None:
                continue
            top = max(row["target_rate"], row["background_rate"], 1e-4)
            fp = row["fisher_p"]
            ax.annotate(f"{row['ratio']:.2f}×"
                        + (f"\np={fp:.3g}" if fp is not None and fp < 0.2 else ""),
                        (xi, top), xytext=(0, 9), textcoords="offset points",
                        ha="center", fontsize=7.2, color=viz.INK2)
        ax.set_yscale("log")
        lo_y, hi_y = ax.get_ylim()
        ax.set_ylim(lo_y, hi_y * 2.2)          # headroom for the ratio labels
        ax.set_xlim(-0.35, len(thr) - 0.5)
        ax.set_xticks(x)
        ax.set_xticklabels([f"p < {t:g}" for t in thr], fontsize=8.5)
        ax.set_ylabel("fraction of features below the threshold")
        ax.set_title("Is trending with stage rare? The background control",
                     fontsize=10)
        ax.legend(fontsize=7.2, loc="lower left")
        ax.text(0.02, 0.985,
                f"smallest p:  target {a['min_p_target']:.1e}   "
                f"background {a['min_p_background']:.1e}",
                transform=ax.transAxes, fontsize=7.4, color=viz.INK2,
                va="top", ha="left")
    else:
        ax.axis("off")
        ax.text(0.02, 0.96,
                "Background control NOT MEASURED.\n\n"
                "gf_background_long.tsv.gz is absent, so there is no\n"
                "prevalence-matched non-target set to compare these\n"
                "trends against. Without it, 'twelve clusters trend'\n"
                "cannot be distinguished from 'everything in the gut\n"
                "trends'. Run background.R on Dardel.",
                va="top", ha="left", fontsize=9, color=viz.INK2, linespacing=1.6,
                transform=ax.transAxes)
    viz.despine(ax)

    fig.suptitle("Monotone trends in biosynthetic gene families along the CRC "
                 "progression axis", fontsize=13)
    _caption(fig, "Left: Spearman $\\rho$ of each cluster's community-level relative "
             "abundance against the ordinal stage scale (Healthy 0, Early_CRC 1, "
             "Advanced_CRC 2), n=277; bars are 95% Fisher-z intervals, and the open "
             "marker is the partial $\\rho$ after residualising both sides' ranks on "
             "age, BMI, Brinkman index, alcohol and sex. Right: the identical test "
             "on five non-target gene families per target, matched decile by decile "
             "on prevalence (medians 65 and 67 of 277 samples detected; "
             "Mann–Whitney p = 0.39, so the two sets have the same power). "
             "Rates are compared on RAW p because BH q depends on family size and "
             "the two families differ fivefold by construction. The two curves cross: "
             "at a loose threshold these enzymes trend LESS than background gene "
             "families, and only in the extreme tail do they trend more. Both "
             "halves of that are the result.")
    viz.save_figure(fig, FIGS / "fig2_trend_forest.png", close=True)


def fig3_stage_profiles(comm_wide, cw, trend, ann, top_n: int = 12):
    """Mean share by stage for the trending clusters, on their own scales."""
    tr = trend[(trend["axis"] == "primary") & (trend["set"] == "target")]
    hits = tr[tr["q"] < FDR].reindex(
        tr[tr["q"] < FDR]["rho_adj"].abs().sort_values(ascending=False).index)
    hits = hits.head(top_n)
    if hits.empty:
        return
    order = _groups_present(cw)
    fig, axs = viz.panel_grid(len(hits), ncols=4, panel_size=(3.3, 2.9))
    up_c, dn_c = viz.DIVERGING[-1], viz.DIVERGING[0]
    for ax, r in zip(axs, hits.itertuples()):
        col = up_c if r.rho > 0 else dn_c
        m, se, _ = _group_means(comm_wide[r.feature], cw, order)
        # The connecting line runs through the ordinal axis ONLY. An earlier draft
        # drew one polyline through all five groups, which put Adenoma between
        # Healthy and Early_CRC and PostSurgery after Advanced_CRC -- asserting an
        # ordinal position for both that STAGE_GROUP_NUMERIC deliberately leaves
        # as NaN. The two off-axis groups are drawn as isolated points instead.
        on = [i for i, g in enumerate(order) if g in PRIMARY_ORDER]
        ax.errorbar(on, m[on], yerr=se[on], marker="o", color=col, lw=1.6,
                    ms=5, capsize=3)
        off = [i for i, g in enumerate(order) if g not in PRIMARY_ORDER]
        if off:
            ax.errorbar(off, m[off], yerr=se[off], fmt="o", ls="none",
                        mfc=viz.SURFACE, mec=viz.MUTED, ms=5.5, mew=1.3,
                        ecolor=viz.MUTED, capsize=3, zorder=4)
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels([g.replace("_", " ") for g in order], rotation=24,
                           ha="right", fontsize=7)
        ax.set_ylim(bottom=0)
        ax.set_ylabel("mean share", fontsize=7.5)
        ax.set_title(f"{r.feature.replace('UniRef90_', '')}  "
                     f"$\\rho$ {r.rho:+.2f}, q={r.q:.1e}\n"
                     f"{_short_species(ann.at[r.feature, 'top_species'] if r.feature in ann.index else '')}",
                     fontsize=8)
        viz.despine(ax)
    fig.suptitle("Stage profile of every cluster that trends with stage", fontsize=13)
    _caption(fig, "Mean community-level relative abundance per group with ± 1 SEM, "
             "over all samples in the group including those with no detection. "
             "Filled markers are the three groups on the primary ordinal axis; "
             "hollow markers are Adenoma and PostSurgery, which are off that axis "
             "by pipeline policy and took no part in the trend test. Panels are "
             "ordered by |adjusted $\\rho$| and each has its own y scale.")
    viz.save_figure(fig, FIGS / "fig3_stage_profiles.png", close=True)


def _abundance_panel(ax, values: pd.Series, cw: pd.DataFrame, order, title, ylab):
    """Per-sample abundance by group: box for the shape, dots for the samples.

    Log y, because the values span four orders of magnitude. Zeros cannot go on a
    log axis and must not be dropped -- dropping them would turn an absence into a
    missing observation -- so they are drawn on a dedicated 'not detected' row
    below an axis break, and remain visible as data.
    """
    g = cw.set_index("Accession")["Stage.Group"]
    pos = {grp: i for i, grp in enumerate(order)}
    pal = viz.group_palette(order, order=utils.EARLY_AXIS_ORDER)
    nz_floor = values[values > 0].min() if (values > 0).any() else 1e-9
    floor = nz_floor / 4.0

    rng = np.random.default_rng(0)
    for grp in order:
        s = values[[a for a in values.index if g.get(a) == grp]]
        nz = s[s > 0]
        x = pos[grp]
        if len(nz) >= 3:
            bp = ax.boxplot([nz.values], positions=[x], widths=0.52,
                            showfliers=False, patch_artist=True, zorder=2)
            bp["boxes"][0].set(facecolor=pal[grp], alpha=0.30,
                               edgecolor=pal[grp], linewidth=1.1)
            for key in ("whiskers", "caps"):
                for a in bp[key]:
                    a.set(color=pal[grp], linewidth=1.0)
            bp["medians"][0].set(color=pal[grp], linewidth=2.0)
        jit = rng.uniform(-0.17, 0.17, size=len(nz))
        ax.scatter(np.full(len(nz), x) + jit, nz.values, s=11, color=pal[grp],
                   alpha=0.75, linewidths=0, zorder=3)
        zeros = s[s == 0]
        if len(zeros):
            jz = rng.uniform(-0.17, 0.17, size=len(zeros))
            ax.scatter(np.full(len(zeros), x) + jz, np.full(len(zeros), floor),
                       s=11, color=pal[grp], alpha=0.35, marker="x",
                       linewidths=0.7, zorder=3)
    ax.axhline(nz_floor / 2.0, color=viz.GRID, lw=0.8, ls=(0, (2, 3)), zorder=1)
    ax.text(-0.46, floor, "not detected", fontsize=7,
            color=viz.MUTED, va="center", ha="left")
    ax.set_yscale("log")
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([g.replace("_", " ") for g in order], rotation=20, ha="right")
    ax.set_ylabel(ylab)
    ax.set_title(title, fontsize=10)
    viz.despine(ax)


def _prevalence_panel(ax, values: pd.Series, cw: pd.DataFrame, order, title):
    g = cw.set_index("Accession")["Stage.Group"]
    pal = viz.group_palette(order, order=utils.EARLY_AXIS_ORDER)
    for i, grp in enumerate(order):
        s = values[[a for a in values.index if g.get(a) == grp]]
        n, k = len(s), int((s > 0).sum())
        lo, hi = _wilson(k, n)
        ax.bar(i, k / n if n else np.nan, width=0.6, color=pal[grp], alpha=0.75,
               linewidth=0)
        ax.plot([i, i], [lo, hi], color=viz.INK2, lw=1.2, zorder=3)
        ax.text(i, hi + 0.035, f"{k}/{n}", ha="center", fontsize=7.5,
                color=viz.MUTED)
    ax.set_ylim(0, 1.12)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([g.replace("_", " ") for g in order], rotation=20, ha="right")
    ax.set_ylabel("samples with the gene detected")
    ax.set_title(title, fontsize=10)
    viz.despine(ax)


def fig4_tryptamine(comm_wide, sp_wide, cw, prim, trend, norm):
    """The chapter's tryptamine claim, answered both ways it can be answered.

    The gene is there and it is attributed to R. gnavus -- that is the claim
    Results 3.7 needs and it holds. It does NOT vary with stage, and panel D says
    so with the number, because a reader who sees panels A-C without it will
    assume a gradient that is not in the data.
    """
    order = _groups_present(cw)
    strat_row = f"{HEADLINE_CLUSTER}|{HEADLINE_SPECIES}"
    fig, axes = plt.subplots(2, 2, figsize=(11.6, 8.8), constrained_layout=True)

    comm = (comm_wide[HEADLINE_CLUSTER] if HEADLINE_CLUSTER in comm_wide.columns
            else pd.Series(0.0, index=comm_wide.index))
    strat = (sp_wide[strat_row] if strat_row in sp_wide.columns
             else pd.Series(0.0, index=sp_wide.index))

    _abundance_panel(axes[0, 0], comm, cw, order,
                     f"A  {HEADLINE_CLUSTER}, community total",
                     "relative abundance of the gene family")
    # Panels A and B are nearly superimposable, and that is a result rather than a
    # duplicated plot: almost the whole community total of this cluster is
    # attributed to one organism. The fraction is computed and printed so no
    # reader has to wonder whether the same series was drawn twice.
    tot = float(comm.sum())
    attributable = float(strat.sum()) / tot if tot > 0 else np.nan
    _abundance_panel(axes[0, 1], strat, cw, order,
                     "B  the same gene, attributed to $R.\\ gnavus$\n"
                     f"{attributable:.1%} of the community total in panel A; the two "
                     "panels coincide because this cluster has essentially one "
                     "carrier",
                     "relative abundance of the gene family")
    _prevalence_panel(axes[1, 0], strat, cw, order,
                      "C  detection of the $R.\\ gnavus$ copy")

    t = trend[(trend["set"] == "target") & (trend["axis"] == "primary") &
              (trend["feature"] == HEADLINE_CLUSTER)]
    rho = f"{t['rho'].iloc[0]:+.3f}" if len(t) else "n/a"
    pv = f"{t['p'].iloc[0]:.2f}" if len(t) else "n/a"
    rho_a = f"{t['rho_adj'].iloc[0]:+.3f}" if len(t) else "n/a"
    sub = prim[(prim["feature"] == HEADLINE_CLUSTER) & prim["tested"]]
    worst = sub["mw_q"].min() if len(sub) else np.nan
    n_abs = (int((comm_wide[ABSENT_CLUSTER] > 0).sum())
             if ABSENT_CLUSTER in comm_wide.columns else 0)

    ax = axes[1, 1]
    ax.axis("off")
    lines = [
        "D  What this does and does not establish",
        "",
        "PRESENT AND ATTRIBUTED.",
        f"{HEADLINE_CLUSTER} is a reviewed UniRef90 cluster of 20 members named",
        "\"Tryptophan decarboxylase\" (RUMGNA_01526, $Ruminococcus\\ gnavus$), the",
        "enzyme Results §3.7 cites as the mechanistic agreement for tryptamine.",
        "It is detected in a majority of samples and HUMAnN attributes it to",
        "$R.\\ gnavus$ itself, so the organism's capacity to perform the reaction",
        "is measured here rather than inferred from a reference genome.",
        "",
        "NOT STAGE-DEPENDENT.",
        f"Spearman $\\rho$ against stage is {rho} (p = {pv}); adjusted, {rho_a}.",
        f"No pairwise contrast comes near significance (best BH q = {worst:.2f}).",
        "The capacity is a property of these guts, not of their disease stage.",
        "",
        "WHY IT IS VISIBLE AT ALL.",
        f"Neither this cluster nor {ABSENT_CLUSTER}, the homologous 21-member",
        "$Clostridium\\ sporogenes$ cluster, appears anywhere in HUMAnN's",
        "UniRef90→EC map or its KO map. Regrouping to EC 4.1.1.105 would have",
        "measured four clusters containing neither enzyme. Both are in this",
        "figure only because the curated route asks for the cluster id directly.",
        "",
        f"THE CONTROL.  {ABSENT_CLUSTER} is detected in {n_abs} of "
        f"{len(comm_wide)} samples.",
        "Its absence is informative precisely because the same probe finds the",
        "homologous cluster in the same samples.",
    ]
    ax.text(0.0, 1.0, "\n".join(lines), va="top", ha="left", fontsize=8.3,
            color=viz.INK2, linespacing=1.5, transform=ax.transAxes)

    note = ("Relative abundance of UniRef90 gene families from stratified HUMAnN 3 "
            f"profiles, n={len(comm_wide)} Yachida 2019 samples. ")
    note += ("Community-level rows sum to 1.000000 in every sample (measured on "
             "Dardel over all 8,930,907 rows), so these are compositional shares "
             "and no renormalization applies. "
             if norm.get("status") == "relative_abundance"
             else f"Per-sample normalization: {norm.get('status')}. ")
    note += ("Crosses below the dashed line are samples with no detection; they are "
             "plotted, not dropped.")
    fig.suptitle("Bacterial tryptophan decarboxylase across the CRC progression axis",
                 fontsize=13)
    _caption(fig, note)
    viz.save_figure(fig, FIGS / "fig4_tryptamine.png", close=True)


def fig5_enzyme_log2fc(enz, emeta, ann, trend, blind):
    """Per-EC aggregate log2FC, both routes, with the EC route's measured blind
    spots shaded rather than left to read as a result."""
    ec_label = ann.attrs["ec_label"]
    fam = enz[enz["family"] == "primary"].copy()
    if fam.empty:
        return
    fam = fam.merge(emeta, on="feature", how="left")
    contrasts = [f"{a}→{b}" for a, b in FAMILIES["primary"]["contrasts"]]
    fam["contrast"] = fam["group1"] + "→" + fam["group2"]
    tr = trend[(trend["axis"] == "primary") & (trend["set"] == "target")] \
        .set_index("feature")

    routes = [("ec", "EC route — HUMAnN UniRef90→EC map"),
              ("curated", "Curated route — reviewed Swiss-Prot clusters")]
    fig, axes = plt.subplots(1, 2, figsize=(13.6, 6.4), constrained_layout=True)
    cmap = viz.diverging_cmap()

    for ax, (route, title) in zip(axes, routes):
        sub = fam[fam["route"] == route]
        ecs = sorted(sub["ec"].dropna().unique())
        if not ecs:
            ax.axis("off")
            ax.set_title(title, fontsize=10)
            continue
        cols = contrasts + ["trend ρ"]
        M = np.full((len(ecs), len(cols)), np.nan)
        for i, ec in enumerate(ecs):
            for j, c in enumerate(contrasts):
                r = sub[(sub["ec"] == ec) & (sub["contrast"] == c)]
                if len(r) and bool(r["tested"].iloc[0]):
                    M[i, j] = r["log2FC"].iloc[0]
            feat = f"{route}|{ec}"
            if feat in tr.index:
                M[i, -1] = tr.at[feat, "rho"]
        finite = M[np.isfinite(M)]
        # Colour limits from the 90th percentile of |value|, not the max. One cell
        # (EC 4.1.1.22, Healthy->Early_CRC, log2FC +2.2 on ten clusters) is four
        # times the next largest and, left to set the scale, flattened every other
        # cell to the neutral midpoint. Each cell also carries its number, so
        # nothing is lost to the clipping.
        lim = (np.percentile(np.abs(finite), 90) if finite.size else 1.0) or 1.0
        norm = viz.diverging_norm(-lim, lim)
        im = ax.imshow(M, cmap=cmap, norm=norm, aspect="auto")
        for i, ec in enumerate(ecs):
            for j, c in enumerate(cols):
                v = M[i, j]
                if np.isfinite(v):
                    ax.text(j, i - 0.16, f"{v:+.2f}", ha="center", va="center",
                            fontsize=7, color=viz.INK)
            for j, c in enumerate(contrasts):
                r = sub[(sub["ec"] == ec) & (sub["contrast"] == c)]
                if not len(r):
                    continue
                q = r["mw_q"].iloc[0]
                if pd.notna(q) and q < FDR:
                    ax.text(j, i + 0.2, "*", ha="center", va="center", fontsize=13,
                            color=viz.INK)
                elif not bool(r["tested"].iloc[0]):
                    ax.text(j, i, "–", ha="center", va="center", fontsize=10,
                            color=viz.MUTED)
            feat = f"{route}|{ec}"
            if feat in tr.index and tr.at[feat, "q"] < FDR:
                ax.text(len(cols) - 1, i + 0.2, "*", ha="center", va="center",
                        fontsize=13, color=viz.INK)
        labels = []
        for ec in ecs:
            act = ec_label.at[ec, "activity"] if ec in ec_label.index else ""
            n_cl = int(sub[sub["ec"] == ec]["n_clusters"].iloc[0])
            tag = ("  ⚠ map is blind here" if (route == "ec" and ec in blind)
                   else "")
            labels.append(f"{ec}  {act} ({n_cl}){tag}")
        ax.set_yticks(range(len(ecs)))
        ax.set_yticklabels(labels, fontsize=8)
        ax.set_xticks(range(len(cols)))
        ax.set_xticklabels([c.replace("_", " ") for c in cols], rotation=18,
                           ha="right", fontsize=8.5)
        ax.set_title(title, fontsize=10)
        ax.grid(False)
        for i, ec in enumerate(ecs):
            if route == "ec" and ec in blind:
                ax.axhspan(i - 0.5, i + 0.5, facecolor="none", edgecolor=viz.INK2,
                           hatch="////", linewidth=1.2, alpha=0.9, zorder=4)
        fig.colorbar(im, ax=ax, shrink=0.75,
                     label="log2 fold change of mean share  /  Spearman $\\rho$"
                           "\n(colour clipped at the 90th percentile; see cell "
                           "values)")

    fig.suptitle("Total biosynthetic capacity per enzyme: no contrast survives "
                 "correction", fontsize=13)
    _caption(fig, "The first three columns are log2 of the ratio of mean summed "
             "community-level relative abundance over the clusters annotated to that "
             "EC; (n) is the number of observed clusters. The last column is the "
             "Spearman $\\rho$ of the same aggregate against the ordinal stage scale. "
             f"* marks BH q < {FDR:g} within the family; – marks a feature below "
             "the detection floor. Shaded EC-route rows are the three ECs where "
             "validate_mapping.py measured the UniRef90→EC map at 0% "
             "sensitivity, so a value there describes the map and not the gut, and "
             "only the curated route speaks to those activities. 4.1.1.105 and "
             "4.1.1.28 have no EC-route row at all here, which is the same finding "
             "in its strongest form: the map contributes no observed cluster to "
             "either activity. The two routes overlap in membership and are never "
             "summed.")
    viz.save_figure(fig, FIGS / "fig5_enzyme_log2fc.png", close=True)


def fig6_volcano(clus, ann, blind):
    """Per-cluster effect against significance for the two Healthy contrasts.

    This figure's content is a negative: no cluster separates two stages once BH
    is applied. It is drawn because the alternative -- reporting only the trend
    result -- would leave a reader to assume the pairwise contrasts were never
    run.
    """
    fam = clus[(clus["family"] == "primary") & clus["tested"]].copy()
    if fam.empty:
        return
    pairs = [("Healthy", "Early_CRC"), ("Healthy", "Advanced_CRC")]
    fig, axes = plt.subplots(1, 2, figsize=(12.4, 5.6), constrained_layout=True)

    def first_ec(cl):
        e = ann.at[cl, "all_ecs"] if cl in ann.index else []
        return e[0] if len(e) else "unannotated"

    fam["ec"] = fam["feature"].map(first_ec)
    cats = sorted(fam["ec"].unique())
    # EC is a nominal variable, so it gets the categorical palette. The first
    # draft used severity_colors, which is the ORDERED ramp: six shades of the
    # same blue for six unrelated enzymes, indistinguishable in the legend and
    # implying a sequence between ECs that does not exist.
    cycle = (viz.PALETTE_COHORT + viz.PALETTE_NODE) * 3
    cols = dict(zip(cats, cycle[:len(cats)]))

    for ax, (g1, g2) in zip(axes, pairs):
        sub = fam[(fam["group1"] == g1) & (fam["group2"] == g2)]
        sub = sub[np.isfinite(sub["log2FC"])]
        y = -np.log10(sub["mw_q"].clip(lower=1e-300))
        for ec in cats:
            m = sub["ec"] == ec
            if not m.any():
                continue
            ax.scatter(sub.loc[m, "log2FC"], y[m], s=16, alpha=0.75,
                       color=cols[ec], linewidths=0,
                       label=f"{ec} ({int(m.sum())})")
        ax.axhline(-np.log10(FDR), color=viz.MUTED, lw=0.9, ls=(0, (3, 3)))
        ax.text(ax.get_xlim()[1], -np.log10(FDR), f" q={FDR:g}", fontsize=7,
                color=viz.MUTED, va="bottom", ha="right")
        ax.axvline(0, color=viz.GRID, lw=0.9)
        ax.set_xlabel(f"log2FC, {g2.replace('_', ' ')} vs {g1}")
        ax.set_ylabel("$-\\log_{10}$ BH q")
        n_sig = int((sub["mw_q"] < FDR).sum())
        ax.set_title(f"{g1} → {g2.replace('_', ' ')}  "
                     f"({n_sig} of {len(sub)} clusters at q<{FDR:g})", fontsize=10)
        viz.despine(ax)
        # One label only. Ranking by q and labelling the top three printed four
        # names on top of each other, because with nothing significant the leading
        # q-values are all equal to the family minimum.
        for r in sub.nsmallest(1, "mw_p").itertuples():
            ax.annotate(f"{r.feature.replace('UniRef90_', '')}  (smallest raw p "
                        f"in this contrast, {r.mw_p:.1e})",
                        (r.log2FC, -np.log10(max(r.mw_q, 1e-300))),
                        fontsize=6.8, color=viz.INK2, xytext=(5, 5),
                        textcoords="offset points")
    h, l = axes[0].get_legend_handles_labels()
    if h:
        axes[0].legend(h, l, fontsize=6.5, loc="upper left", title="EC",
                       title_fontsize=7)
    fig.suptitle("Pairwise stage contrasts per cluster: a negative result",
                 fontsize=13)
    _caption(fig, "One point per observed UniRef90 cluster at community level. q is BH "
             "over the whole primary family (features × three contrasts). "
             "Nothing crosses the threshold: the stage signal in these gene families "
             "is monotone across the axis rather than a difference between two "
             f"stages. Clusters annotated to {', '.join(sorted(blind))} appear here "
             "only because the curated route named them directly; the EC map does "
             "not contain them.")
    viz.save_figure(fig, FIGS / "fig6_cluster_volcano.png", close=True)


# ===========================================================================
# 7. MAIN
# ===========================================================================

def main() -> None:
    FIGS.mkdir(parents=True, exist_ok=True)
    viz.apply_theme()

    cw = load_crosswalk()
    accessions = cw["Accession"].tolist()
    print(f"crosswalk: {len(cw)} samples; "
          f"{cw['Stage.Group'].value_counts().to_dict()}")

    comm_wide, sp_wide, species_long, counts = load_long(LONG, accessions)
    print(f"long rows kept {counts['n_long']}: community {counts['n_community']}, "
          f"|unclassified {counts['n_unclassified']}, species {counts['n_species']}")
    print(f"community clusters {comm_wide.shape[1]}, species rows {sp_wide.shape[1]}")

    norm = check_normalization()
    print(f"normalization: {norm['status']}")
    if norm.get("renormalized"):
        f = (pd.read_csv(NORMALIZATION, sep="\t").set_index("sample")["sum_community"])
        comm_wide = comm_wide.div(f.reindex(comm_wide.index), axis=0)
        sp_wide = sp_wide.div(f.reindex(sp_wide.index), axis=0)

    ann = cluster_annotation()
    blind = blind_spot_ecs()
    print(f"EC-route blind spots (measured): {sorted(blind)}")

    enz_wide, emeta, overlap = enzyme_aggregates(comm_wide, ann)
    print(f"enzyme aggregates: {enz_wide.shape[1]} features; "
          f"{overlap} multi-EC cluster memberships")

    n_primary = int(cw["Stage.Group"].isin(utils.STAGE_GROUP_ORDER).sum())
    min_prev_n = max(MIN_PREV_SAMPLES, int(np.ceil(MIN_PREV_FRACTION * n_primary)))
    print(f"detection floor: {min_prev_n} samples "
          f"({MIN_PREV_FRACTION:.0%} of the n={n_primary} primary set)")

    clus = pd.concat([run_family(comm_wide, cw, f, min_prev_n) for f in FAMILIES],
                     ignore_index=True)
    enz = pd.concat([run_family(enz_wide, cw, f, min_prev_n) for f in FAMILIES],
                    ignore_index=True)
    spp = pd.concat([run_family(sp_wide, cw, f, min_prev_n) for f in FAMILIES],
                    ignore_index=True)

    trend = pd.concat([run_trend(comm_wide, cw, min_prev_n, "target"),
                       run_trend(enz_wide, cw, min_prev_n, "target")],
                      ignore_index=True)

    bg_trend, bgcmp = pd.DataFrame(), {"status": "not measured",
                                       "reason": f"{BACKGROUND_LONG.name} absent"}
    if BACKGROUND_LONG.exists():
        bg_comm, _, _, bg_counts = load_long(BACKGROUND_LONG, accessions)
        print(f"background: {bg_comm.shape[1]} clusters, "
              f"{bg_counts['n_community']} nonzero community triples")
        bg_trend = run_trend(bg_comm, cw, min_prev_n, "background")
        bgcmp = compare_to_background(trend, bg_trend, "primary")
        bg_trend.to_csv(OUT / "background_trend.tsv", sep="\t", index=False)
        print("background control: " + json.dumps(bgcmp["adjusted"], indent=2))
    else:
        print(f"background control NOT measured: {BACKGROUND_LONG.name} absent")

    tp0 = trend[(trend["set"] == "target") & (trend["axis"] == "primary")]
    hit_feats = tp0[tp0["q"] < FDR]["feature"].tolist()
    needed = {ann.at[f, "top_species"] for f in hit_feats if f in ann.index}
    carriers, carrier_report = carrier_abundance(needed, cw)
    print(f"carrier control: {carrier_report.get('status')}, "
          f"matched {carrier_report.get('matched')} of "
          f"{carrier_report.get('requested')}; "
          f"unmatched {carrier_report.get('unmatched')}")
    carrier = carrier_adjusted_trend(comm_wide, cw, hit_feats, ann, carriers)
    carrier.to_csv(OUT / "progression_carrier_control.tsv", sep="\t", index=False)
    if "p_carrier_adj" in carrier.columns:
        surv = carrier[(carrier["carrier_available"]) &
                       (carrier["p_carrier_adj"] < FDR)]
        print(f"trends surviving carrier adjustment: {len(surv)} of "
              f"{int(carrier['carrier_available'].sum())} testable")

    cw.to_csv(OUT / "sample_crosswalk.tsv", sep="\t", index=False)
    clus.to_csv(OUT / "progression_cluster_tests.tsv", sep="\t", index=False)
    enz.merge(emeta, on="feature", how="left").to_csv(
        OUT / "progression_enzyme_tests.tsv", sep="\t", index=False)
    spp.to_csv(OUT / "progression_species_tests.tsv", sep="\t", index=False)
    trend_out = trend.copy()
    trend_out["top_species"] = trend_out["feature"].map(
        lambda f: ann.at[f, "top_species"] if f in ann.index else "")
    trend_out["ecs"] = trend_out["feature"].map(
        lambda f: ";".join(ann.at[f, "all_ecs"]) if f in ann.index else f)
    trend_out.to_csv(OUT / "progression_trend.tsv", sep="\t", index=False)

    prim = clus[clus["family"] == "primary"]
    fig1_source_turnover(comm_wide, enz_wide, cw, trend, ann, carrier, blind)
    fig2_trend_forest(trend, ann, bgcmp)
    fig3_stage_profiles(comm_wide, cw, trend, ann)
    fig4_tryptamine(comm_wide, sp_wide, cw, prim, trend, norm)
    fig5_enzyme_log2fc(enz, emeta, ann, trend, blind)
    fig6_volcano(clus, ann, blind)

    tp = trend[(trend["set"] == "target") & (trend["axis"] == "primary")]
    hits = tp[tp["q"] < FDR]
    prov = {
        "built_utc": pd.Timestamp.now(tz="UTC").strftime("%Y-%m-%dT%H:%M:%SZ"),
        "inputs": {p.name: {"bytes": p.stat().st_size} for p in
                   [COHORT_META, CURATION, LONG, BACKGROUND_LONG, SUMMARY,
                    TARGETS, CURATED, NORMALIZATION] if p.exists()},
        "normalization": norm,
        "crosswalk": {"n_samples": int(len(cw)),
                      "groups": cw["Stage.Group"].value_counts().to_dict()},
        "row_kinds": counts,
        "features": {"clusters": int(comm_wide.shape[1]),
                     "species_rows": int(sp_wide.shape[1]),
                     "enzyme_aggregates": int(enz_wide.shape[1]),
                     "multi_ec_memberships": int(overlap)},
        "detection_floor_samples": int(min_prev_n),
        "fdr": FDR,
        "blind_spot_ecs": sorted(blind),
        "pairwise": {
            "primary_tested": int((prim["tested"]).sum()),
            "primary_q_lt_fdr": int((prim["mw_q"] < FDR).sum()),
            "primary_min_q": (float(prim["mw_q"].min())
                              if prim["mw_q"].notna().any() else None),
        },
        "trend": {
            "primary_tested": int(len(tp)),
            "primary_q_lt_fdr": int(len(hits)),
            "primary_q_adj_lt_fdr": int((tp["q_adj"] < FDR).sum()),
            "hits": [{"feature": r.feature, "ecs": ";".join(
                          ann.at[r.feature, "all_ecs"] if r.feature in ann.index
                          else []),
                      "top_species": (ann.at[r.feature, "top_species"]
                                      if r.feature in ann.index else ""),
                      "rho": round(float(r.rho), 4),
                      "q": float(r.q),
                      "rho_adj": round(float(r.rho_adj), 4),
                      "q_adj": float(r.q_adj)}
                     for r in hits.sort_values("rho", ascending=False).itertuples()],
        },
        "background_control": bgcmp,
        "carrier_control": {
            **carrier_report,
            "testable": (int(carrier["carrier_available"].sum())
                         if "carrier_available" in carrier.columns else 0),
            "surviving_q_lt_fdr": (
                int(((carrier.get("carrier_available", False)) &
                     (carrier.get("q_carrier_adj", pd.Series(dtype=float)) < FDR)
                     ).sum()) if "q_carrier_adj" in carrier.columns else 0),
            "per_feature": (carrier.replace({np.nan: None}).to_dict("records")
                            if not carrier.empty else []),
        },
        "software": {"python": platform.python_version(), "numpy": np.__version__,
                     "pandas": pd.__version__, "matplotlib": matplotlib.__version__},
    }
    (OUT / "progression_provenance.json").write_text(json.dumps(prov, indent=2))
    print(json.dumps({"pairwise": prov["pairwise"],
                      "trend": {k: v for k, v in prov["trend"].items()
                                if k != "hits"},
                      "background_control_adjusted":
                          bgcmp.get("adjusted", {}).get("tail"),
                      "carrier_control": {k: v for k, v
                                          in prov["carrier_control"].items()
                                          if k != "per_feature"}}, indent=2))


if __name__ == "__main__":
    main()
