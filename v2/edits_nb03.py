"""NB03 edits: robust-hit target selection, grouped stage scheme, leak guard."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from nbedit import Notebook

nb = Notebook(Path(__file__).resolve().parent / "03_ml_benchmarking.ipynb")

nb.set("c00-title", r'''# Notebook 03: ML Benchmarking — CRC Metabolome & Microbial Producer Analysis (v2.1)

**Scope: `YACHIDA-CRC-2019`, primary analysis set only** (n=277: Healthy + Early_CRC +
Advanced_CRC). Sinha and Kim are genus-level and are used for external validation in NB05,
not here. MP is excluded from the study; PostSurgery (HS) was split out in NB02, so the
pickle this notebook reads contains primary samples only.

## Pipeline

1. **Targets** = the union of
   - **robust DA metabolites** — FDR-significant in *both* the unadjusted rank test and the
     confounder-adjusted model, with the same direction (NB02's `Robust` column), from
     `Healthy vs Early_CRC` and `Healthy vs Advanced_CRC`; up to 8 elevated + 7 depleted per
     comparison. There is no nominal-p fallback: a comparison with few robust hits
     contributes few targets, and that is reported rather than padded.
   - **the 7 pre-specified polyamines** — the study's a-priori hypothesis targets, analysed
     whether or not they are differentially abundant, labelled `prespecified_polyamine` and
     never described as "dysregulated".
2. Visualise abundance differences (boxplots + cross-stage comparison).
3. Benchmark one XGBoost classifier per target (species CLR + clinical confounders →
   above/below-median metabolite level):
   - **Performance** = `AUC_nested`: hyperparameters tuned by `RandomizedSearchCV` *inside*
     each outer fold (5 outer folds, 10 candidates, 3 inner folds), scored on the untouched
     outer test fold.
   - **Significance** = permutation test with **pre-specified, data-independent**
     hyperparameters (`FIXED_PARAMS`) for both observed and label-shuffled scores on the
     same folds, BH-corrected across targets (`Perm_qvalue`).
4. Producer analysis for the pre-specified polyamines: out-of-fold SHAP from the same
   fixed-parameter model, with a **direction** per species (Spearman between a species'
   value and its SHAP value) and bootstrap rank CIs.

**Prerequisite:** run `01_preprocessing.ipynb` and `02_analysis.ipynb` first.

---

## Design notes (what was wrong before, and what this does instead)

1. **No disease-status proxy among the features.** Tumour location is "-" for every
   non-cancer sample, so on a Healthy+CRC set a tumour-site dummy is an exact disease
   indicator. `utils.build_confounder_matrix()` returns Age, BMI, Brinkman.Index, Alcohol and
   Gender only, and *raises* if tumour location is requested for a set containing non-cancer
   samples. Measured cost of the old leak: NB05's stage-1 accuracy fell 0.63 → 0.39 once it
   was removed.
2. **Confounders are imputed INSIDE each CV fold** (`IterativeImputer` inside a
   `ColumnTransformer` pipeline step), never on the whole cohort.
3. **No oversampling.** A median split gives balanced classes by construction, so
   synthetic-minority oversampling only added noise. Targets whose minority class has
   < 20 samples (ties at the median) are skipped and logged.
4. **Tuning is nested.** The previous version tuned on all samples and then cross-validated
   on those same samples, which inflates AUC. The permutation test uses fixed hyperparameters
   so observed and null scores come from an identical procedure (the test is exact); with
   `N_PERM = 100` the smallest attainable p is 1/101 ≈ 0.0099.
5. **OOF SHAP is the producer ranking, and it is signed.** `Direction` = sign of the Spearman
   ρ between a species' abundance and its SHAP contribution to "high metabolite". Only
   positively associated species are candidate producers downstream (NB08/NB09) — several
   top-|SHAP| species are *negatively* associated with their polyamine, and the old pipeline
   would have called those producers.
6. **No cross-notebook evidence triangulation here.** That lives at the end of
   `09_pathway_network.ipynb`, and it starts from NB02's partial correlations rather than
   from this notebook's SHAP, so SHAP is not both the candidate generator and its own
   evidence.
7. **Sign conventions** (NB02 `_v2.csv` tables): `MeanDiff = mean(group2) − mean(group1)` and
   `Coef` (adjusted model) share the same sign — **> 0 means higher in the case group** —
   while `CliffsD` carries the **opposite** sign by construction.''')

nb.set("c01-imports", r'''import sys, pickle, warnings
sys.path.insert(0, '.')
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from collections import Counter

from sklearn.experimental import enable_iterative_imputer  # noqa: F401 (must precede IterativeImputer import)
from sklearn.impute import IterativeImputer
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from joblib import Parallel, delayed
from sklearn.model_selection import (
    StratifiedKFold, RepeatedStratifiedKFold, RandomizedSearchCV,
)
from sklearn.metrics import roc_curve, auc
from scipy.stats import spearmanr, rankdata
from statsmodels.stats.multitest import multipletests
from joblib import parallel_backend
import os, time

try:
    from xgboost import XGBClassifier
except ImportError as e:
    raise ImportError(
        "xgboost is required: the v2 pipeline collapses the old 5-classifier "
        "zoo to a single tuned XGBoost model (see markdown above). "
        "Install with `pip install xgboost`."
    ) from e

try:
    import shap
except ImportError as e:
    raise ImportError(
        "shap is required for the producer (Section 3) analysis. "
        "Install with `pip install shap`."
    ) from e

import utils
import viz
from utils import RESULTS_DIR, DATA_DIR, POLYAMINE_KEGG, FDR_THRESHOLD, STAGE_GROUP_ORDER

for sub in ['figures/ml', 'tables', 'intermediate', 'logs']:
    (RESULTS_DIR / sub).mkdir(parents=True, exist_ok=True)

# One shared theme (see viz.py). This notebook used dpi=110 and
# set_context('paper', font_scale=1.2) while NB01/02/06 used dpi=100 and no context, so
# identical fontsize= arguments rendered at different physical sizes between notebooks.
viz.apply_theme()

CRC_DS         = 'YACHIDA-CRC-2019'
# Grouped stage scheme (utils.STAGE_GROUP_MAP): the same groups NB02's robust-hit
# tables are built on, so target selection and modelling use one definition.
GROUP_EARLY    = 'Early_CRC'
GROUP_ADVANCED = 'Advanced_CRC'

print('xgboost and shap available.')
print(f'Stage groups in use: {STAGE_GROUP_ORDER}')
print('Setup complete.')''')

nb.set("c03-load-pickle", r'''# -- Load the feature-reduced pickle (YACHIDA, primary analysis set) -----------
with open(RESULTS_DIR / 'intermediate' / 'preprocessed_data_reduced_v2.pkl', 'rb') as f:
    _prep = pickle.load(f)

spc = _prep['transformed_species'][CRC_DS].copy()
mtb = _prep['transformed_mtb'][CRC_DS].copy()
meta_h = _prep['harmonized_meta'][CRC_DS].copy()
polyamine_columns_ds = _prep['polyamine_columns'][CRC_DS]

print(f"Pickle sample scope: {_prep.get('sample_scope', 'unspecified')}")
print(f'{CRC_DS}: {len(spc)} samples | {spc.shape[1]-1} species (reduced) | '
      f'{mtb.shape[1]-1} metabolites (reduced)')

# The pickle must already be primary-only; assert rather than re-filter, so a
# scope change upstream surfaces here instead of being silently absorbed.
_groups = set(meta_h['Stage.Group'].dropna())
assert _groups <= set(STAGE_GROUP_ORDER), \
    f'non-primary group(s) in the ML input: {_groups - set(STAGE_GROUP_ORDER)}'
assert not meta_h['Study.Group'].isin(['MP', 'HS']).any(), \
    'MP/HS must not reach the ML notebook'

print(f'\nStage.Group distribution:\n{meta_h["Stage.Group"].value_counts().reindex(STAGE_GROUP_ORDER).to_string()}')
print(f'\nFine-grained Study.Group:\n{meta_h["Study.Group"].value_counts().to_string()}')''')

nb.set("c04-confounders", r'''# -- Confounder matrix: Age, BMI, Brinkman.Index, Alcohol, Gender --------------
# From the reduced pickle's harmonized metadata (same samples as spc/mtb).
# Tumour location is NOT requested: this is a Healthy+CRC set, and
# utils.build_confounder_matrix() refuses include_tumor_location=True on any set
# containing a non-cancer sample precisely so this cannot be done by accident.
conf_df = utils.build_confounder_matrix(meta_h)

IMPUTE_COLS = ['Age', 'BMI', 'Brinkman.Index', 'Alcohol']  # genuinely-missing-capable only
CONF_COLS   = [c for c in conf_df.columns if c != 'Sample']

# Guard against any disease-status proxy reaching the feature matrix.
_banned = [c for c in CONF_COLS
           if c.startswith(('TumorSite_', 'TumorLoc_'))
           or 'Stage' in c or c.startswith('Tumor')]
assert not _banned, f'disease-status proxy in features: {_banned}'

# And prove the guard is live rather than trusting it.
try:
    utils.build_confounder_matrix(meta_h, include_tumor_location=True)
    raise AssertionError('tumour-location leak guard is NOT active')
except ValueError:
    print('Leak guard verified: tumour location is refused on this Healthy+CRC set.')

print('Confounder matrix columns:', CONF_COLS)
print('\nMissingness (imputed per-fold below, not here):')
for col in IMPUTE_COLS:
    n_na = conf_df[col].isna().sum()
    print(f'  {col}: {n_na} / {len(conf_df)} missing')
print(conf_df[CONF_COLS].describe().round(2))''')

nb.set("c05-section1-hdr", r'''---
## Section 1: Target metabolites

### Sign conventions (NB02 tables)

| Column | Sign | Meaning |
|---|---|---|
| `MeanDiff` = mean(case) − mean(Healthy) | **> 0** | metabolite **higher in CRC** |
| `Coef` (confounder-adjusted model) | **> 0** | metabolite **higher in CRC** |
| `CliffsD` | **< 0** | metabolite **higher in CRC** (elevated) |
| `CliffsD` | **> 0** | metabolite **lower in CRC** (depleted) |

`CliffsD` therefore has the **opposite** sign to `MeanDiff`/`Coef`. That is a property of
`mannwhitneyu(Healthy, CRC)`: scipy's U counts (Healthy > CRC) pairs, so
`CliffsD = 2U/(n₁n₂) − 1` is positive when Healthy values rank higher.

### Selection

- **Robust DA targets.** Only metabolites with `Robust == True` in NB02's grouped tables —
  q < 0.05 in **both** the Mann-Whitney test and the confounder-adjusted linear model, with
  **matching direction**. Up to 8 most elevated (most negative `CliffsD`) + 7 most depleted
  per comparison, from `Healthy vs Early_CRC` and `Healthy vs Advanced_CRC`.

  Requiring both tests matters here specifically: smoking is strongly imbalanced across
  these groups (median Brinkman index 100 in Healthy vs 196–390 in CRC), so a metabolite
  that passes only the unadjusted test may be tracking smoking rather than cancer. Those
  `mw_only` features are counted and reported below, not silently dropped — and not
  modelled.

- **Pre-specified polyamines.** Every polyamine column present, analysed as an a-priori
  hypothesis target regardless of DA significance (`Target_type = prespecified_polyamine`),
  with its DA q-values reported alongside. This is what keeps the study's central question
  answerable even if no polyamine survives FDR.''')

nb.set("c06-top15-select", r'''# -- Load NB02's grouped, robust-annotated DA table ---------------------------
# v2.1: one tidy table (Dataset | Comparison | Scheme | Assay | Feature | ...) replaces
# the six per-comparison CSVs. Sliced by Comparison + Assay rather than by filename.
_DA_TIDY = RESULTS_DIR / 'tables' / 'differential_abundance_v2.csv'


def _load_da(group):
    if not _DA_TIDY.exists():
        raise FileNotFoundError(
            f'{_DA_TIDY.name} not found. NB02 Step 6B writes it; run 02_analysis.ipynb '
            'first. (Before v2.1 this read per-comparison '
            f'{CRC_DS}_Healthy_vs_<group>_metabolites_DA_v2.csv files.)')
    da = pd.read_csv(_DA_TIDY)
    for col in ('Robust', 'Agreement', 'QValue_adj', 'Coef_adj'):
        if col not in da.columns:
            raise KeyError(
                f'{_DA_TIDY.name} has no {col!r} column -- it was written by an older NB02 '
                'that did not run the confounder-adjusted model. Re-run 02_analysis.ipynb.')
    df = da[(da['Dataset'] == CRC_DS)
            & (da['Comparison'] == f'Healthy_vs_{group}')
            & (da['Assay'] == 'metabolites')].reset_index(drop=True)
    if df.empty:
        raise ValueError(
            f'{_DA_TIDY.name} has no rows for Healthy_vs_{group} metabolites in {CRC_DS}. '
            f'Comparisons present: {sorted(da["Comparison"].unique())}')
    return df

early_da    = _load_da(GROUP_EARLY)
advanced_da = _load_da(GROUP_ADVANCED)

POLYAMINE_COL_NAMES = set(polyamine_columns_ds.values())
KEGG_TO_NAME = {v: POLYAMINE_KEGG.get(k, k) for k, v in polyamine_columns_ds.items()}


def short_name(feat):
    """Strip the KEGG prefix for display."""
    parts = feat.split('_', 1)
    return parts[1] if len(parts) > 1 else feat


def select_robust_balanced(da_df, n_up=8, n_down=7, label=''):
    """Robust-only selection: among features significant in BOTH tests with the same
    direction, take the n_up most elevated and n_down most depleted (by CliffsD).

    No nominal-p and no MW-only fallback: a comparison with few robust hits
    contributes few targets.
    """
    rob = da_df[da_df['Robust'] == True].copy()        # noqa: E712 (CSV may read as object)
    up   = rob[rob['CliffsD'] < 0].nsmallest(n_up, 'CliffsD')
    down = rob[rob['CliffsD'] > 0].nlargest(n_down, 'CliffsD')
    sel = pd.concat([up, down]).reset_index(drop=True)
    sel['is_polyamine'] = sel['Feature'].isin(POLYAMINE_COL_NAMES)
    sel['direction'] = sel['CliffsD'].apply(
        lambda d: 'Elevated in CRC' if d < 0 else 'Depleted in CRC')

    if label:
        agree = da_df['Agreement'].value_counts().to_dict()
        n_mw  = int((da_df['QValue'] < FDR_THRESHOLD).sum())
        n_rob = int((da_df['Robust'] == True).sum())     # noqa: E712
        print(f'{label}:')
        print(f'    unadjusted q<0.05 : {n_mw}')
        print(f'    robust (both+dir) : {n_rob}')
        print(f'    dropped by adjustment (mw_only): {agree.get("mw_only", 0)}'
              f'  | revealed by adjustment (lm_only): {agree.get("lm_only", 0)}'
              f'  | direction conflict: {agree.get("both_q_opposite_dir", 0)}')
        print(f'    selected as targets: {len(sel)} '
              f'({int((sel["CliffsD"] < 0).sum())} elevated + '
              f'{int((sel["CliffsD"] > 0).sum())} depleted), '
              f'{int(sel["is_polyamine"].sum())} polyamine(s)')
    return sel


print('=== Robust DA target selection (grouped stage scheme) ===')
top15_early    = select_robust_balanced(early_da,    label=f'Healthy vs {GROUP_EARLY}')
top15_advanced = select_robust_balanced(advanced_da, label=f'Healthy vs {GROUP_ADVANCED}')

top15_early_set    = set(top15_early['Feature'])
top15_advanced_set = set(top15_advanced['Feature'])

# -- Pre-specified polyamine targets (a-priori hypothesis; not DA-selected) ----
def _da_row(da_df, feat):
    r = da_df[da_df['Feature'] == feat]
    if r.empty:
        return (np.nan, np.nan, np.nan, 'not_tested')
    return (float(r['QValue'].iloc[0]), float(r['CliffsD'].iloc[0]),
            float(r['QValue_adj'].iloc[0]) if pd.notna(r['QValue_adj'].iloc[0]) else np.nan,
            str(r['Agreement'].iloc[0]))

pa_target_rows = []
for kegg_id, col in polyamine_columns_ds.items():
    if col not in mtb.columns:
        print(f'  NOTE: polyamine {POLYAMINE_KEGG.get(kegg_id)} ({col}) absent from the '
              f'reduced matrix -- cannot be modelled')
        continue
    q_e, d_e, qadj_e, ag_e = _da_row(early_da, col)
    q_a, d_a, qadj_a, ag_a = _da_row(advanced_da, col)
    pa_target_rows.append({
        'Feature': col, 'KEGG': kegg_id, 'PA_name': POLYAMINE_KEGG.get(kegg_id, kegg_id),
        'QValue_early': q_e, 'CliffsD_early': d_e, 'QValue_adj_early': qadj_e,
        'Agreement_early': ag_e,
        'QValue_advanced': q_a, 'CliffsD_advanced': d_a, 'QValue_adj_advanced': qadj_a,
        'Agreement_advanced': ag_a,
        'Target_type': 'prespecified_polyamine'})
pa_targets = pd.DataFrame(pa_target_rows)

print(f'\nPre-specified polyamine targets: {len(pa_targets)}')
if not pa_targets.empty:
    print(pa_targets[['PA_name', 'QValue_early', 'QValue_adj_early', 'Agreement_early',
                      'QValue_advanced', 'QValue_adj_advanced', 'Agreement_advanced']]
          .round(4).to_string(index=False))

print(f'\n=== Robust targets: Healthy vs {GROUP_EARLY} ===')
print(top15_early[['Feature', 'CliffsD', 'QValue', 'QValue_adj', 'direction',
                   'is_polyamine']].to_string(index=False) if not top15_early.empty
      else '  (none -- this comparison contributes no robust DA targets)')
print(f'\n=== Robust targets: Healthy vs {GROUP_ADVANCED} ===')
print(top15_advanced[['Feature', 'CliffsD', 'QValue', 'QValue_adj', 'direction',
                      'is_polyamine']].to_string(index=False) if not top15_advanced.empty
      else '  (none -- this comparison contributes no robust DA targets)')''')

# ---- boxplots: switch to the grouped scheme ----
nb.set("c07-boxplot-early", r'''# -- Boxplots: robust targets, Healthy vs Early CRC ---------------------------
stage_mask_early = meta_h['Stage.Group'].isin(['Healthy', GROUP_EARLY])
meta_e   = meta_h[stage_mask_early].copy()
mtb_e    = mtb[mtb['Sample'].isin(meta_e['Sample'])]
merged_e = meta_e.merge(mtb_e, on='Sample')

if top15_early.empty:
    print(f'No robust DA targets for Healthy vs {GROUP_EARLY} -- no boxplot panel.')
else:
    n = len(top15_early)
    ncols = min(5, n)
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 4, nrows * 3.5), squeeze=False)
    axes_flat = axes.flatten()

    for i, (_, row) in enumerate(top15_early.iterrows()):
        ax = axes_flat[i]
        feat = row['Feature']
        vals_h = merged_e[merged_e['Stage.Group'] == 'Healthy'][feat].dropna()
        vals_c = merged_e[merged_e['Stage.Group'] == GROUP_EARLY][feat].dropna()

        crc_color = '#d63031' if row['CliffsD'] < 0 else '#0984e3'
        bplot = ax.boxplot([vals_h, vals_c], patch_artist=True, widths=0.5,
                            showfliers=False, medianprops=dict(color='black', linewidth=1.5))
        bplot['boxes'][0].set_facecolor('#dfe6e9'); bplot['boxes'][0].set_alpha(0.85)
        bplot['boxes'][1].set_facecolor(crc_color); bplot['boxes'][1].set_alpha(0.85)

        ax.set_xticklabels(['Healthy', 'Early\nCRC'], fontsize=9)
        ax.set_title(short_name(feat), fontsize=8.5, fontweight='bold',
                     color='#c0392b' if row['is_polyamine'] else 'black')
        ax.set_ylabel('log₁₀(abund.)', fontsize=7.5)
        ax.text(0.97, 0.97, f"q={row['QValue']:.3f}\nq_adj={row['QValue_adj']:.3f}",
                transform=ax.transAxes, ha='right', va='top', fontsize=6.5, style='italic')
        dir_label = 'Elevated' if row['CliffsD'] < 0 else 'Depleted'
        ax.text(0.03, 0.97, dir_label, transform=ax.transAxes,
                ha='left', va='top', fontsize=7, fontweight='bold', color=crc_color)
        if row['is_polyamine']:
            for spine in ax.spines.values():
                spine.set_edgecolor('#c0392b'); spine.set_linewidth(2)

    for j in range(n, len(axes_flat)):
        axes_flat[j].set_visible(False)

    fig.suptitle('Robust DA metabolites — Early CRC vs Healthy\n'
                 '(significant in BOTH the rank test and the confounder-adjusted model)\n'
                 'Red box = elevated in CRC | Blue box = depleted | Red border = polyamine',
                 fontsize=11, fontweight='bold')
    plt.tight_layout()
    out = RESULTS_DIR / 'figures' / 'ml' / 'top15_metabolites_early_CRC_v2.png'
    plt.savefig(out, dpi=150, bbox_inches='tight')
    plt.show(); plt.close()
    print(f'Saved: {out}')''')

nb.set("c08-interp-early", r'''### Interpretation: early-stage metabolite changes

Only metabolites robust to confounder adjustment are used as DA targets. If few or none
survive for `Healthy vs Early_CRC`, this comparison contributes few or no targets — reported
as is, never padded with nominally significant metabolites. The counts printed above show
how many unadjusted hits were lost to adjustment (`mw_only`), which is the quantity of
interest: in this cohort smoking differs 2–4× across groups, so unadjusted early-stage
"signal" is the most likely place for confounding to masquerade as biology.

Polyamines are analysed regardless, as pre-specified targets.''')

nb.set("c09-boxplot-advanced", r'''# -- Boxplots: robust targets, Healthy vs Advanced CRC ------------------------
stage_mask_adv = meta_h['Stage.Group'].isin(['Healthy', GROUP_ADVANCED])
meta_a   = meta_h[stage_mask_adv].copy()
mtb_a    = mtb[mtb['Sample'].isin(meta_a['Sample'])]
merged_a = meta_a.merge(mtb_a, on='Sample')

if top15_advanced.empty:
    print(f'No robust DA targets for Healthy vs {GROUP_ADVANCED} -- no boxplot panel.')
else:
    n = len(top15_advanced)
    ncols = min(5, n)
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 4, nrows * 3.5), squeeze=False)
    axes_flat = axes.flatten()

    for i, (_, row) in enumerate(top15_advanced.iterrows()):
        ax = axes_flat[i]
        feat = row['Feature']
        vals_h = merged_a[merged_a['Stage.Group'] == 'Healthy'][feat].dropna()
        vals_c = merged_a[merged_a['Stage.Group'] == GROUP_ADVANCED][feat].dropna()

        crc_color = '#d63031' if row['CliffsD'] < 0 else '#0984e3'
        bplot = ax.boxplot([vals_h, vals_c], patch_artist=True, widths=0.5,
                            showfliers=False, medianprops=dict(color='black', linewidth=1.5))
        bplot['boxes'][0].set_facecolor('#dfe6e9'); bplot['boxes'][0].set_alpha(0.85)
        bplot['boxes'][1].set_facecolor(crc_color); bplot['boxes'][1].set_alpha(0.85)

        ax.set_xticklabels(['Healthy', 'Adv.\nCRC'], fontsize=9)
        ax.set_title(short_name(feat), fontsize=8.5, fontweight='bold',
                     color='#c0392b' if row['is_polyamine'] else 'black')
        ax.set_ylabel('log₁₀(abund.)', fontsize=7.5)
        ax.text(0.97, 0.97, f"q={row['QValue']:.3f}\nq_adj={row['QValue_adj']:.3f}",
                transform=ax.transAxes, ha='right', va='top', fontsize=6.5, style='italic')
        dir_label = 'Elevated' if row['CliffsD'] < 0 else 'Depleted'
        ax.text(0.03, 0.97, dir_label, transform=ax.transAxes,
                ha='left', va='top', fontsize=7, fontweight='bold', color=crc_color)
        if row['is_polyamine']:
            for spine in ax.spines.values():
                spine.set_edgecolor('#c0392b'); spine.set_linewidth(2)

    for j in range(n, len(axes_flat)):
        axes_flat[j].set_visible(False)

    fig.suptitle('Robust DA metabolites — Advanced CRC vs Healthy\n'
                 '(significant in BOTH the rank test and the confounder-adjusted model)\n'
                 'Red box = elevated in CRC | Blue box = depleted | Red border = polyamine',
                 fontsize=11, fontweight='bold')
    plt.tight_layout()
    out = RESULTS_DIR / 'figures' / 'ml' / 'top15_metabolites_advanced_CRC_v2.png'
    plt.savefig(out, dpi=150, bbox_inches='tight')
    plt.show(); plt.close()
    print(f'Saved: {out}')''')

nb.set("c10-crossstage-scatter", r'''# -- Cross-stage comparison scatter -------------------------------------------
compare = (early_da[['Feature', 'CliffsD', 'Robust']]
           .rename(columns={'CliffsD': 'CD_early', 'Robust': 'Rob_early'})
           .merge(advanced_da[['Feature', 'CliffsD', 'Robust']]
                  .rename(columns={'CliffsD': 'CD_advanced', 'Robust': 'Rob_advanced'}),
                  on='Feature'))
compare['is_polyamine'] = compare['Feature'].isin(POLYAMINE_COL_NAMES)
compare['in_top15_e']   = compare['Feature'].isin(top15_early_set)
compare['in_top15_a']   = compare['Feature'].isin(top15_advanced_set)

fig, ax = plt.subplots(figsize=(8, 7))
ax.scatter(compare['CD_early'], compare['CD_advanced'],
           c='#bdc3c7', s=25, alpha=0.6, zorder=1, label='Other metabolites')
mask_e = compare['in_top15_e'] & ~compare['is_polyamine']
ax.scatter(compare.loc[mask_e, 'CD_early'], compare.loc[mask_e, 'CD_advanced'],
           c='#3498db', s=60, alpha=0.9, zorder=2, label='Robust target (Early)')
mask_a = compare['in_top15_a'] & ~compare['is_polyamine']
ax.scatter(compare.loc[mask_a, 'CD_early'], compare.loc[mask_a, 'CD_advanced'],
           c='#e67e22', s=60, alpha=0.9, zorder=2, label='Robust target (Advanced)')
mask_pa = compare['is_polyamine']
ax.scatter(compare.loc[mask_pa, 'CD_early'], compare.loc[mask_pa, 'CD_advanced'],
           c='#c0392b', marker='*', s=200, zorder=3, label='Polyamines (pre-specified)')
for _, r in compare[mask_pa].iterrows():
    name = KEGG_TO_NAME.get(r['Feature'], short_name(r['Feature']))
    ax.annotate(name, (r['CD_early'], r['CD_advanced']), fontsize=8, color='#c0392b',
                fontweight='bold', xytext=(6, 4), textcoords='offset points')

ax.axhline(0, color='gray', lw=0.8, ls='--')
ax.axvline(0, color='gray', lw=0.8, ls='--')
ax.set_xlabel("Cliff's δ  (Early CRC vs Healthy)\n"
              "negative = elevated in CRC  |  positive = depleted in CRC", fontsize=11)
ax.set_ylabel("Cliff's δ  (Advanced CRC vs Healthy)\n"
              "negative = elevated in CRC  |  positive = depleted in CRC", fontsize=11)
ax.set_title("Metabolite dysregulation: Early vs Advanced CRC\n"
             "(coloured points are robust in that comparison)", fontsize=11, fontweight='bold')
ax.legend(loc='lower right', fontsize=9, framealpha=0.9)
plt.tight_layout()
out = RESULTS_DIR / 'figures' / 'ml' / 'metabolite_dysregulation_comparison_v2.png'
plt.savefig(out, dpi=150, bbox_inches='tight')
plt.show(); plt.close()
print(f'Saved: {out}')''')

nb.set("c11-section2-hdr", r'''---
## Section 2: ML Benchmark — Microbial Source Prediction

For each target metabolite (robust DA targets + pre-specified polyamines):
- **Features**: species CLR abundances + Age, BMI, Gender, Brinkman Index, Alcohol.
  No tumour location, no stage indicator — see the leak guard in Section 0.
- **Target**: binary high/low metabolite level, split at the cohort median.
- **Model**: XGBoost; confounder imputation (`IterativeImputer`) and `StandardScaler` are
  pipeline steps refit on each training fold.
- **`AUC_nested`** (performance): 5-fold outer CV; inside each training fold a
  `RandomizedSearchCV` (10 candidates, 3 inner folds, ROC-AUC) selects hyperparameters and
  the refit model is scored on the untouched outer test fold.
- **`Perm_pvalue` / `Perm_qvalue`** (significance): observed AUC of the fixed-parameter
  model (`AUC_fixed`) against `N_PERM` label permutations on the same 5 folds, then BH
  across targets. Fixed parameters are what make the comparison exact — the null and the
  observed score come from an identical procedure.
- A timing probe prints the projected runtime before the loop (`NB03_N_PERM` overrides
  `N_PERM`).

**What an AUC here does and does not mean.** The target is a median split of a *measured*
metabolite, so a high AUC says the species composition carries information about that
metabolite's level in this cohort. It is not evidence that those species produce it: that
requires direction (Section 3), genomic capability (NB08), and ideally flux (NB10), which is
why the final synthesis in NB09 counts those as separate evidence sources.''')

nb.set("c13-benchmark-loop", r'''# -- Benchmark: nested-CV AUC + fixed-parameter permutation test -------------
da_features = list(dict.fromkeys(
    top15_early['Feature'].tolist() + top15_advanced['Feature'].tolist()))
pa_features = pa_targets['Feature'].tolist() if not pa_targets.empty else []
target_metabolites = list(dict.fromkeys(da_features + pa_features))
print(f'Target metabolites: {len(target_metabolites)} '
      f'({len(da_features)} robust DA, {len(pa_features)} pre-specified polyamines, '
      f'{len(set(da_features) & set(pa_features))} both)')
if not da_features:
    print('NOTE: no robust DA targets survived confounder adjustment in either comparison.\n'
          '      The benchmark therefore runs on the pre-specified polyamines alone, which\n'
          '      is a valid a-priori analysis -- not a fallback that was tuned to find hits.')

direction_map = {}
_sel_all = [d for d in (top15_early, top15_advanced) if not d.empty]
if _sel_all:
    for _, r in pd.concat(_sel_all).drop_duplicates('Feature').iterrows():
        direction_map[r['Feature']] = 'Elevated' if r['CliffsD'] < 0 else 'Depleted'

N_OUTER     = 5
N_ITER      = 10
N_INNER     = 3
N_PERM      = int(os.environ.get('NB03_N_PERM', 100))
MIN_CLASS_N = 20

# Live progress (notebook output is only written when the notebook finishes).
PROGRESS = RESULTS_DIR / 'logs' / 'nb03_progress_v2.txt'
PROGRESS.parent.mkdir(parents=True, exist_ok=True)
def progress(msg):
    with open(PROGRESS, 'a') as fh:
        fh.write(f"{time.strftime('%H:%M:%S')}  {msg}\n")
PROGRESS.write_text('')

# -- timing probe: project the total runtime before committing to it --
if target_metabolites:
    _X, _y, _fc, _ = build_feature_matrix(spc, mtb, target_metabolites[0], conf_df)
    _yb = binarize_target(_y)
    _cache = build_fold_cache(_X, _yb, n_splits=N_OUTER)
    _t0 = time.time(); run_cv_auc_cached(_cache, _yb); _t_fit = (time.time() - _t0) / N_OUTER
    _fits_per_target = N_OUTER * (N_ITER * N_INNER + 1) + N_OUTER * (N_PERM + 1)
    _proj_h = len(target_metabolites) * _fits_per_target * _t_fit / 2 / 3600   # 2 workers
    print(f'Timing probe: {_t_fit:.2f} s/fit -> projected {_proj_h:.2f} h for '
          f'{len(target_metabolites)} targets ({_fits_per_target} fits each, 2 workers)')
    progress(f'probe {_t_fit:.2f} s/fit (fixed params), projected {_proj_h:.2f} h, '
             f'{len(target_metabolites)} targets, N_PERM={N_PERM}')
    if _proj_h > 2.75:
        print('WARNING: projected NB03 runtime exceeds 2.75 h')

bench_rows = []
roc_store  = {}
tuned_params_by_target = {}
skipped = []

for col_name in target_metabolites:
    if col_name not in mtb.columns:
        print(f'  SKIP (not in MTB): {col_name}')
        skipped.append((col_name, 'absent from reduced matrix')); continue

    X, y, feat_cols, _ = build_feature_matrix(spc, mtb, col_name, conf_df)
    y_bin = binarize_target(y)
    pos_n, neg_n = int(y_bin.sum()), int(len(y_bin) - y_bin.sum())
    if min(pos_n, neg_n) < MIN_CLASS_N:
        print(f'  SKIP (minority class < {MIN_CLASS_N}, ties at median): {col_name}  '
              f'pos={pos_n}, neg={neg_n}')
        skipped.append((col_name, f'minority class {min(pos_n, neg_n)} < {MIN_CLASS_N}'))
        continue

    # Performance: nested CV
    nested_aucs, chosen = nested_cv_auc(X, y_bin, n_outer=N_OUTER, n_iter=N_ITER,
                                        inner=N_INNER, random_state=42, n_jobs=2)
    tuned_params_by_target[col_name] = chosen

    # Significance: fixed-parameter observed vs permuted, identical folds
    fold_cache = build_fold_cache(X, y_bin, n_splits=N_OUTER, random_state=42)
    mean_fpr, mean_tpr, fixed_aucs = run_cv_auc_cached(fold_cache, y_bin, None, random_state=42)
    roc_store[col_name] = (mean_fpr, mean_tpr, fixed_aucs)
    fixed_auc = float(fixed_aucs.mean())
    perm_aucs = run_permutation_null_cached(fold_cache, y_bin, None, n_perm=N_PERM,
                                            random_state=42, base_seed=1000, n_jobs=2)
    perm_p = (np.sum(perm_aucs >= fixed_auc) + 1) / (N_PERM + 1)

    in_e, in_a = col_name in top15_early_set, col_name in top15_advanced_set
    stage_tag = 'Both' if in_e and in_a else 'Early' if in_e else 'Advanced' if in_a else 'none'
    is_pa_target = col_name in pa_features
    target_type = ('robust_DA+prespecified' if (is_pa_target and (in_e or in_a)) else
                   'prespecified_polyamine' if is_pa_target else 'robust_DA')

    bench_rows.append({
        'Metabolite':          col_name,
        'Short_name':          short_name(col_name),
        'Classifier':          'XGBoost',
        'Target_type':         target_type,
        'AUC_nested':          float(nested_aucs.mean()),
        'AUC_nested_std':      float(nested_aucs.std()),
        'AUC_fixed':           fixed_auc,
        'N':                   len(y_bin),
        'N_high':              pos_n,
        'N_low':               neg_n,
        'is_polyamine':        col_name in POLYAMINE_COL_NAMES,
        'CRC_stage':           stage_tag,
        'Direction':           direction_map.get(col_name, 'not DA-selected'),
        'Perm_pvalue':         float(perm_p),
        'Perm_null_AUC_mean':  float(perm_aucs.mean()),
        'N_perm':              N_PERM,
        'Nested_params':       chosen,
    })
    progress(f'{len(bench_rows)}/{len(target_metabolites)} {short_name(col_name)[:40]} '
             f'AUC_nested={nested_aucs.mean():.3f} AUC_fixed={fixed_auc:.3f} perm_p={perm_p:.4f}')
    print(f'  {short_name(col_name):40s} n={len(y_bin):3d}  AUC_nested={nested_aucs.mean():.3f}'
          f'±{nested_aucs.std():.3f}  AUC_fixed={fixed_auc:.3f}  perm_p={perm_p:.4f}  '
          f'(null AUC={perm_aucs.mean():.3f})')

bench_df = pd.DataFrame(bench_rows)
if not bench_df.empty:
    bench_df['Perm_qvalue'] = multipletests(bench_df['Perm_pvalue'], method='fdr_bh')[1]
    # backwards-compatible aliases used by the plots/tables below
    bench_df['AUC_mean'] = bench_df['AUC_nested']
    bench_df['AUC_std']  = bench_df['AUC_nested_std']
    n_q05 = int((bench_df['Perm_qvalue'] < 0.05).sum())
    print(f'\nBenchmark done: {len(bench_df)} metabolites modelled, {n_q05} with Perm_qvalue < 0.05')
else:
    print('\nBenchmark produced no rows -- every target was skipped.')
if skipped:
    print(f'Skipped {len(skipped)} target(s):')
    for name, why in skipped:
        print(f'  {short_name(name)[:45]:45s} {why}')''')

nb.save()
