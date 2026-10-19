"""NB02 edits, part A: intro, loading/primary masking, feature reduction, dual DA."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from nbedit import Notebook

nb = Notebook(Path(__file__).resolve().parent / "02_analysis.ipynb")

# ---------------------------------------------------------------- intro
nb.set("cell-001", r'''# Notebook 02: Statistical Analysis, Correlation, Network & Progression (v2.1)

**Cohort scope:** `YACHIDA-CRC-2019` (discovery) plus `SINHA_CRC_2016` and
`Kim_adenomas_2020` (external validation) — `utils.DATASETS_ANALYSIS`. The IBD/IBS cohorts
appear only in NB01's exploratory sections. `ERAWIJANTARI-GASTRIC-CANCER-2020` is never used:
53 of its samples are the same Healthy subjects as Yachida's, so it is not independent.

## Sample scope: the primary analysis set

Everything in this notebook runs on the **primary** set, `utils.primary_mask()` —
Healthy + Early_CRC + Advanced_CRC, **n=277** for Yachida. The two groups that are not on
that axis are handled once, here, and never leak downstream:

- **MP** (multiple polypoid adenomas) was dropped at load time in NB01 and does not exist
  in any table this notebook reads.
- **PostSurgery** (HS — cancer resected, currently cancer-free) is split out of the
  matrices in the first section below. The reduced pickle that NB03+ consume therefore
  contains the **primary samples only** under its usual keys, with the PostSurgery rows
  stored separately as `sensitivity_species` / `sensitivity_mtb`. A downstream notebook that
  simply loads `transformed_species` gets the correct sample set by default rather than by
  remembering to filter — the failure mode is absence, not contamination.
  NB01b already quantified how different that group is (its own FDR family).

## Differential abundance is run two ways, and only agreement counts

A rank test and a covariate-adjusted compositional model answer different questions, and the
earlier pipeline reported only the first:

| | Test | Adjusts for confounders | Compositional |
|---|---|---|---|
| `utils.differential_abundance` | Mann-Whitney U + Cliff's δ | no | no |
| `utils.clr_lm_da` | linear model on CLR values | **yes** | **yes** (LinDA bias correction) |

A feature is a **robust DA hit** only when *both* reach q<0.05 **and agree in direction**.
This matters in this cohort specifically: smoking is strongly imbalanced across groups
(median Brinkman index 100 in Healthy vs 196–390 in the CRC groups), so an unadjusted test
cannot distinguish a cancer effect from a smoking effect. Features that pass the rank test
but fail the adjusted model are reported, not hidden — they are candidates explained by
confounding.

**Metabolites** are tested with the same linear model but with `bias_correct=False` and
`already_clr=True`: the metabolite matrix is log10-transformed, not compositional, so the
LinDA mode-subtraction does not apply to it. Only the species model is compositional.

## Steps
- **Primary split + global feature reduction** — top-variance species/metabolites retained
  for this notebook and all of NB03+, saved to `preprocessed_data_reduced_v2.pkl`.
- **Step 6A** — confound / residual diagnostics (the old version's OLS fits all failed
  silently and always reported R²=0.000; fixed via `utils.build_formula`).
- **Step 6** — DA on the fine-grained `Study.Group` scheme (`MeanDiff`, not the old,
  frequently-NaN `Log2FC`).
- **Step 6B** — DA on the grouped `Stage.Group` scheme (Healthy / Early_CRC / Advanced_CRC),
  **both tests**, producing the robust-hit tables NB03 selects its targets from.
- **Step 7** — correlations: species–polyamine, full species–metabolite, partial
  correlations adjusted for the full confounder set (saved for NB09), and **three-cohort**
  replication (Yachida vs Sinha vs Kim) at **genus** level with metabolites matched by
  KEGG ID.
- **Step 8** — networks: species–metabolite co-occurrence, polyamine-centric subnetworks,
  and stage-stratified hub networks.
- **Step 10** — progression: ordinal Spearman against CRC stage (10A) and binary
  early-vs-late DA (10B).

There is **no ML step in this notebook.** The old version's intro advertised a "Step 9: ML
Pipeline" for which no code existed anywhere in it. All ML lives in `03_ml_benchmarking.ipynb`,
which reads this notebook's `preprocessed_data_reduced_v2.pkl`.''')

# ---------------------------------------------------------------- imports
nb.set("cell-002", r'''import sys
sys.path.insert(0, '.')

import os
import pickle
import warnings
warnings.filterwarnings('ignore', category=FutureWarning)

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
from statsmodels.formula.api import ols as sm_ols
from scipy.stats import probplot

import utils
import viz
from utils import (
    RESULTS_DIR,
    DATASETS_ANALYSIS, DATASET_DISCOVERY, DATASETS_VALIDATION,
    POLYAMINE_DATASETS_ANALYSIS, POLYAMINE_KEGG,
    COMPARISONS, COMPARISONS_GROUP,
    CRC_STAGE_ORDER, CRC_STAGE_NUMERIC, CRC_PROGRESSION_COMPARISONS,
    STAGE_GROUP_MAP, STAGE_GROUP_ORDER,
    YACHIDA_COMPLETE_CONFOUNDERS,
    MAX_SPECIES_REDUCED, MAX_MTB_REDUCED, FDR_THRESHOLD,
    differential_abundance, clr_lm_da,
    compute_correlations, find_consistent_correlations,
    volcano_plot, plot_correlation_heatmap,
    progression_correlation_analysis, plot_stage_trend,
    build_formula, primary_mask, sensitivity_mask, apply_mask,
)

for sub in ['figures/qc', 'figures/differential', 'figures/correlation',
            'figures/network', 'figures/progression', 'tables', 'intermediate']:
    (RESULTS_DIR / sub).mkdir(parents=True, exist_ok=True)

viz.apply_theme()   # one shared theme: see viz.py (dpi used to be 100 here, 110 in NB03/05)

print('Libraries loaded.')
print(f'Discovery : {DATASET_DISCOVERY}')
print(f'Validation: {DATASETS_VALIDATION}')''')

# ---------------------------------------------------------------- load + primary split
nb.set("cell-003", r'''# Load preprocessed data from Notebook 01
with open(RESULTS_DIR / 'intermediate' / 'preprocessed_data_v2.pkl', 'rb') as f:
    preprocessed = pickle.load(f)

transformed_species = {ds: preprocessed['transformed_species'][ds] for ds in DATASETS_ANALYSIS}
transformed_mtb     = {ds: preprocessed['transformed_mtb'][ds]     for ds in DATASETS_ANALYSIS}
harmonized_meta     = {ds: preprocessed['harmonized_meta'][ds]     for ds in DATASETS_ANALYSIS}
polyamine_columns   = {ds: preprocessed['polyamine_columns'].get(ds, {})
                       for ds in DATASETS_ANALYSIS}
polyamine_proxies   = {ds: preprocessed.get('polyamine_proxies', {}).get(ds, {})
                       for ds in DATASETS_ANALYSIS}
# mtb_kegg_map: {dataset: {mtb_column_name: kegg_id}} -- built in NB01 from each
# dataset's own mtb.map file. Needed below for KEGG-ID-based cross-cohort matching
# (Step 7.4), since the three cohorts use different metabolomics naming conventions.
mtb_kegg_map  = {ds: preprocessed['mtb_kegg_map'].get(ds, {}) for ds in DATASETS_ANALYSIS}
MTB_PRESCALED = preprocessed.get('mtb_prescaled', {})

# MP must be gone already (dropped at load in NB01).
assert not harmonized_meta[DATASET_DISCOVERY]['Study.Group'].eq('MP').any(), \
    'MP must have been excluded in NB01'

# ---------------------------------------------------------------------------
# PRIMARY SPLIT -- done once, here, for every cohort.
# The sensitivity group (PostSurgery / HS) is removed from the working matrices
# and kept aside. Everything after this cell operates on the primary set only.
# ---------------------------------------------------------------------------
sensitivity_species, sensitivity_mtb, sensitivity_meta = {}, {}, {}

for ds in DATASETS_ANALYSIS:
    meta = harmonized_meta[ds]
    pm, sm = primary_mask(meta), sensitivity_mask(meta)

    if sm.any():
        s_spc, s_mtb = apply_mask(sm, transformed_species[ds], transformed_mtb[ds], meta_df=meta)
        sensitivity_species[ds] = s_spc
        sensitivity_mtb[ds]     = s_mtb
        sensitivity_meta[ds]    = meta[sm].reset_index(drop=True)

    p_spc, p_mtb = apply_mask(pm, transformed_species[ds], transformed_mtb[ds], meta_df=meta)
    transformed_species[ds] = p_spc
    transformed_mtb[ds]     = p_mtb
    harmonized_meta[ds]     = meta[pm].reset_index(drop=True)

    print(f'{ds}: primary n={int(pm.sum())}'
          + (f', held aside as sensitivity n={int(sm.sum())}' if sm.any() else ''))

# Hard guard: no non-primary group may survive in the working tables.
for ds in DATASETS_ANALYSIS:
    bad = set(harmonized_meta[ds]['Stage.Group'].dropna()) - set(STAGE_GROUP_ORDER)
    assert not bad, f'{ds}: non-primary group(s) still present: {bad}'
    for name, tbl in (('species', transformed_species[ds]), ('mtb', transformed_mtb[ds])):
        assert set(tbl['Sample']) == set(harmonized_meta[ds]['Sample']), \
            f'{ds}: {name} samples do not match the primary metadata'

print('\nPreprocessed data loaded; primary analysis set in force.')
for ds in DATASETS_ANALYSIS:
    n_spc  = len(transformed_species[ds].columns) - 1
    n_mtb  = len(transformed_mtb[ds].columns) - 1
    n_pa   = len(polyamine_columns[ds])
    n_px   = len(polyamine_proxies[ds])
    n_kegg = len(mtb_kegg_map[ds])
    print(f'  {ds}: {len(transformed_species[ds])} samples, {n_spc} species, '
          f'{n_mtb} metabolites, {n_pa} polyamines (+{n_px} proxy), {n_kegg} KEGG mappings')''')

# ---------------------------------------------------------------- feature reduction
nb.set("cell-004", r'''---
## Primary split + Global Feature Reduction

Applied immediately after loading. All subsequent cells here (DA, correlations, hub
networks, progression) and NB03+ (ML) use this reduced feature set. Sizes come from the
shared constants `utils.MAX_SPECIES_REDUCED` (500) and `utils.MAX_MTB_REDUCED` (150), not
local magic numbers.

**Caveat (named, not resolved):** this is a pooled-variance pre-filter, computed across all
primary samples before any group split. A feature whose variance is driven by a real
between-group difference belongs in the top-variance set for a good reason, but a feature
with high *within-Healthy* variance for technical reasons can also make the cut ahead of a
smaller, perfectly real group difference. That is a mild selection bias on every downstream
comparison in this notebook. Polyamines and their proxies are exempt — they are always
retained, so the study's pre-specified targets never depend on this filter.

The reduction is computed on the **primary** samples only, so the sensitivity group cannot
influence which features the whole pipeline keeps.''')

nb.set("cell-005", r'''for ds in DATASETS_ANALYSIS:
    spc_cols = [c for c in transformed_species[ds].columns if c != "Sample"]
    mtb_cols = [c for c in transformed_mtb[ds].columns     if c != "Sample"]

    # Polyamines AND their pre-specified proxies are protected from both filters.
    pa_cols = (list(polyamine_columns.get(ds, {}).values())
               + [col for col, _via in polyamine_proxies.get(ds, {}).values()])

    # Species: top-variance selection, polyamine-pathway-agnostic (no target leakage)
    top_spc = (transformed_species[ds][spc_cols]
               .var().nlargest(MAX_SPECIES_REDUCED).index.tolist())
    transformed_species[ds] = transformed_species[ds][["Sample"] + top_spc]

    # Metabolites: top-variance selection + polyamine/proxy protection
    top_mtb = (transformed_mtb[ds][mtb_cols]
               .var().nlargest(MAX_MTB_REDUCED).index.tolist())
    top_mtb = list(dict.fromkeys(top_mtb + [c for c in pa_cols if c in mtb_cols]))
    transformed_mtb[ds] = transformed_mtb[ds][["Sample"] + top_mtb]

    # Keep the sensitivity tables on the SAME feature set, so NB06 can project
    # those samples onto factors learned from the primary set.
    if ds in sensitivity_species:
        sensitivity_species[ds] = sensitivity_species[ds][["Sample"] + top_spc]
        sensitivity_mtb[ds]     = sensitivity_mtb[ds][["Sample"] + top_mtb]

    n_kept_pa = sum(1 for c in pa_cols if c in top_mtb)
    print(f"{ds}: species {len(spc_cols)} -> {len(top_spc)} | "
          f"metabolites {len(mtb_cols)} -> {len(top_mtb)} "
          f"({n_kept_pa}/{len(pa_cols)} polyamine+proxy columns protected)")

# Undetected-metabolite mask (pre-imputation value NaN or <= 0) for the retained
# metabolites -- NB06 uses it to give MOFA+ genuine missing values instead of the
# half-min/log/centred floor. Not defined for pre-scaled tables, where "<= 0" has
# no detection meaning (a z-score below zero is simply below average).
undetected_mtb = {}
for ds in DATASETS_ANALYSIS:
    if MTB_PRESCALED.get(ds):
        print(f"{ds}: pre-scaled metabolites -- no undetected mask (<=0 is not 'undetected')")
        continue
    fm = preprocessed['filtered_mtb'][ds].set_index('Sample')
    keep_cols = [c for c in transformed_mtb[ds].columns if c != 'Sample' and c in fm.columns]
    vals = fm.loc[transformed_mtb[ds]['Sample'], keep_cols]
    und = (vals.isna() | (vals <= 0)).reset_index()
    undetected_mtb[ds] = und
    print(f"{ds}: undetected metabolite cells = {und[keep_cols].values.mean():.1%}")

# Save reduced data so NB03+ load the same feature set AND the same sample set.
reduced_pkl_path = RESULTS_DIR / "intermediate" / "preprocessed_data_reduced_v2.pkl"
with open(reduced_pkl_path, "wb") as f:
    pickle.dump({
        # Primary analysis set -- the default keys every downstream notebook reads.
        "transformed_species": transformed_species,
        "transformed_mtb":     transformed_mtb,
        "harmonized_meta":     harmonized_meta,
        # Sensitivity group, same features, separate keys: opt-in only (NB06 projects it).
        "sensitivity_species": sensitivity_species,
        "sensitivity_mtb":     sensitivity_mtb,
        "sensitivity_meta":    sensitivity_meta,
        "polyamine_columns":   polyamine_columns,
        "polyamine_proxies":   polyamine_proxies,
        "mtb_kegg_map":        mtb_kegg_map,
        "datasets":            DATASETS_ANALYSIS,
        "dataset_discovery":   DATASET_DISCOVERY,
        "datasets_validation": DATASETS_VALIDATION,
        "undetected_mtb":      undetected_mtb,
        "mtb_prescaled":       MTB_PRESCALED,
        "stage_group_order":   STAGE_GROUP_ORDER,
        "sample_scope":        "primary (utils.primary_mask); sensitivity group under sensitivity_* keys",
        "feature_reduction":   {"method": "variance_topN",
                                 "MAX_SPECIES_REDUCED": MAX_SPECIES_REDUCED,
                                 "MAX_MTB_REDUCED": MAX_MTB_REDUCED,
                                 "protected": "polyamines + pre-specified proxies"},
    }, f)
print(f"\nReduced pkl saved -> {reduced_pkl_path} "
      f"({reduced_pkl_path.stat().st_size / 1e6:.1f} MB)")
print(f"  primary samples per cohort     : "
      f"{ {ds: len(transformed_species[ds]) for ds in DATASETS_ANALYSIS} }")
print(f"  sensitivity samples held aside : "
      f"{ {ds: len(v) for ds, v in sensitivity_species.items()} }")''')

# ---------------------------------------------------------------- Step 6A confounders
nb.set("cell-006", r'''---
## Step 6A: Residual / Confound Diagnostics

How much variance in the transformed species data is explained by sample-level covariates?
A high per-feature R² is what motivates adjusting for them, which Step 6B now does
(`utils.clr_lm_da`). Purely diagnostic — it does not modify `transformed_species` or the
pickle.

Covariates: the full Yachida set (`utils.YACHIDA_COMPLETE_CONFOUNDERS` — Age, BMI, Gender,
Brinkman.Index, Alcohol) plus `Stage.Group`, so the diagnostic reports on the same variables
the adjusted DA model uses. Tumour location is *not* among them: it is "-" for every
non-cancer sample, so on this mixed set it would simply re-encode disease status.

**Fix applied:** the old notebook quoted a formula term in `Q("...")` only if the column name
contained a *space*. `Study.Group` contains a `.`, not a space, so patsy tried to parse it as
`Study` attribute-accessed `.Group` — a parse error on *every* feature. A bare `except: pass`
swallowed them all, so the notebook reported a fabricated R²=0.000 for every dataset without
ever fitting a regression. This version builds every formula with `utils.build_formula()`
(which quotes any name containing anything other than letters/digits/underscore) and counts
and prints fit failures instead of swallowing them.''')

nb.set("cell-007", r'''CONFOUNDERS    = YACHIDA_COMPLETE_CONFOUNDERS + ['Stage.Group']
MAX_FEAT_RESID = 100   # test top-variance species (speed)

for ds in DATASETS_ANALYSIS:
    meta     = harmonized_meta[ds].copy()
    spc      = transformed_species[ds]
    spc_cols = [c for c in spc.columns if c != 'Sample']
    top_spc  = spc[spc_cols].var().nlargest(MAX_FEAT_RESID).index.tolist()
    merged   = meta.merge(spc[['Sample'] + top_spc], on='Sample')

    avail = [c for c in CONFOUNDERS if c in merged.columns
             and merged[c].notna().mean() > 0.3
             and merged[c].nunique(dropna=True) > 1]
    if len(avail) < 1:
        print(f'{ds}: insufficient confounder data -- skip'); continue

    # NOTE: categoricals (Gender, Stage.Group) are left as strings on purpose.
    # patsy auto-detects object-dtype columns and dummy-encodes them inside the
    # formula; the old code pd.factorize()'d them into arbitrary integer codes and
    # fit them as continuous, which is wrong for a nominal variable.

    r2_vals, n_failed, fail_msgs = [], 0, []
    for feat in top_spc:
        sub = merged[avail + [feat]].dropna()
        if len(sub) < 10:
            continue
        formula = build_formula(feat, avail)
        try:
            r2_vals.append(sm_ols(formula, data=sub).fit().rsquared)
        except Exception as e:
            n_failed += 1
            if len(fail_msgs) < 3:
                fail_msgs.append(f'{feat[:40]}: {e}')

    status = 'none -- all fits succeeded' if n_failed == 0 else 'see examples below'
    print(f'{ds}: {n_failed} / {len(top_spc)} OLS fits failed ({status})')
    for msg in fail_msgs:
        print(f'    FAILED: {msg[:120]}')

    r2_arr = np.array(r2_vals) if r2_vals else np.array([0.0])
    print(f'{ds} -- Confounder R2 (top-{len(top_spc)} species, adjusting for {avail}):')
    print(f'  mean={r2_arr.mean():.3f}  median={np.median(r2_arr):.3f}  '
          f'max={r2_arr.max():.3f}  >0.10: {(r2_arr > 0.1).mean():.1%}')

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    axes[0].hist(r2_arr, bins=20, color='#3498db', edgecolor='white', alpha=0.85)
    axes[0].axvline(0.1, color='red', lw=1.5, ls='--', label='R2=0.10')
    axes[0].set_xlabel('Confounder R2 per species', fontsize=11)
    axes[0].set_ylabel('Count', fontsize=11)
    axes[0].set_title(f'{ds} Variance explained by covariates', fontsize=10, fontweight='bold')
    axes[0].legend()

    probplot(spc[top_spc[0]].dropna(), plot=axes[1])
    axes[1].set_title(f'Q-Q: {top_spc[0][:40]} (highest-variance species)', fontsize=10, fontweight='bold')

    plt.suptitle(f'{ds} -- Residual / Confound Diagnostics', fontsize=12, fontweight='bold')
    plt.tight_layout()
    # v2.1: was RESULTS_DIR/'figures'/... i.e. the figures ROOT, the only figure in the
    # pipeline that bypassed the figures/<topic>/ convention NB01 sets up.
    out_r = RESULTS_DIR / 'figures' / 'qc' / f'residual_diagnostics_{ds}_v2.png'
    viz.save_figure(plt.gcf(), out_r)
    print(f'  Saved: {out_r}')''')

nb.save()
