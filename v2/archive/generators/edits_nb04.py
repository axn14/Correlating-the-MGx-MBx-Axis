"""Rebuild NB04 as the within-cancer stage / tumour-site analysis.

The old notebook's entire subject was the MP (multiple polypoid adenomas) group,
which v2.1 excludes from the study, so there is nothing left to edit -- the
notebook is restated. Its new job is the one analysis in the whole pipeline where
tumour location is a legitimate covariate: every sample is a carcinoma.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from nbedit import Notebook

HERE = Path(__file__).resolve().parent
OLD = HERE / "04_mp_trajectory.ipynb"
NEW = HERE / "04_crc_stage_tumour_site.ipynb"

CELLS = [
("markdown", "nb04-title", r'''# 04 — Within-Cancer Stage Progression & Tumour Site

**Scope: the 150 Yachida samples with a diagnosed carcinoma** (`Stage_0` n=27,
`Stage_I_II` n=69, `Stage_III_IV` n=54), selected by `utils.crc_only_mask()`.
Healthy and PostSurgery samples are deliberately absent.

## Why this notebook exists, and why it is the only place tumour site appears

`Tumor location` is recorded as `"-"` for **every** non-cancer sample in Yachida (all 127
Healthy and all 30 HS). So in any comparison that includes non-cancer samples, a tumour-site
covariate is an exact disease-status indicator: knowing "site = NoTumor" *is* knowing the
patient has no cancer. Using it there inflates every model that sees it — measured directly
in this pipeline, NB05's stage-1 accuracy fell from 0.63 to 0.39 once it was removed.

Within this notebook's sample set that objection disappears, because **every** sample has a
located tumour. Site is therefore a genuine clinical covariate here: rectal, left-colon and
right-colon cancers differ in embryonic origin, blood supply, microbial environment and
clinical behaviour, so adjusting for site is the *correct* thing to do when asking how
polyamines change with stage.

`utils.build_confounder_matrix(meta, include_tumor_location=True)` enforces this rule rather
than documenting it: it raises `ValueError` if any row of the metadata it is given is not a
carcinoma. The first code cell below deliberately triggers that error on the full cohort to
demonstrate the guard is live, then uses the CRC-only subset.

## Analyses

1. **Stage trend per polyamine**, across the ordinal cancer stages (Stage_0 → I/II → III/IV):
   - unadjusted Spearman ρ;
   - **partial** Spearman adjusting for tumour site plus Age, BMI, Gender, Brinkman index
     and Alcohol — i.e. does the stage association survive adjustment for where the tumour
     is and who the patient is?
   - BH-corrected across the 7 polyamines.
2. **Tumour-site differences** in each polyamine (Kruskal–Wallis across Rectum / Left colon /
   Right colon), which is also the check on whether site adjustment was necessary.
3. **Per-site stage trends** (exploratory): the stage trend repeated within each site.
   Sites are small (Rectum 69, Right 44, Left 32, MultiSite 5), so these are reported as
   exploratory and MultiSite is excluded from per-site tests.
4. **Species stage trend**, same adjusted design, as the species-side counterpart.

**Sample-size honesty.** n=150 split three ways, and the site subgroups are smaller still.
Nothing here is powered for a small effect; the point is that the adjusted estimate is the
defensible one, not that this notebook settles the question.

**Downstream:** NB08 reads `nb04_polyamine_stage_trend_v2.csv` for its `Stage_trend_*`
evidence columns. Nothing else consumes these outputs.'''),

("code", "nb04-imports", r'''import sys, pickle, warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.stats import spearmanr, kruskal
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings('ignore')

sys.path.insert(0, '.')
import utils
import viz
from utils import (
    RESULTS_DIR, POLYAMINE_KEGG, CRC_STAGE_NUMERIC,
    CRC_ONLY_STUDY_GROUPS, TUMOR_SITE_ORDER, YACHIDA_COMPLETE_CONFOUNDERS,
    FDR_THRESHOLD, differential_abundance, volcano_plot,
    crc_only_mask, apply_mask,
)

try:
    import pingouin as pg
    HAS_PINGOUIN = True
except ImportError:
    HAS_PINGOUIN = False
    print('pingouin missing -- partial correlations will be skipped (pip install pingouin)')

PREP_PKL = RESULTS_DIR / 'intermediate' / 'preprocessed_data_reduced_v2.pkl'
CRC_DS   = 'YACHIDA-CRC-2019'

TABLES_DIR = RESULTS_DIR / 'tables'
FIG_DIR    = RESULTS_DIR / 'figures' / 'progression'
TABLES_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

viz.apply_theme()   # one shared theme: see viz.py (this notebook set no dpi at all)

print('Imports OK')
print(f'Cancer stages analysed: {CRC_ONLY_STUDY_GROUPS}')
print(f'Tumour sites: {TUMOR_SITE_ORDER}')'''),

("code", "nb04-load", r'''if not PREP_PKL.exists():
    raise FileNotFoundError(
        f'{PREP_PKL} not found -- run 02_analysis.ipynb first. '
        'This notebook depends on its reduced-feature output pickle.')

with open(PREP_PKL, 'rb') as f:
    _prep = pickle.load(f)

spc_all  = _prep['transformed_species'][CRC_DS].copy()
mtb_all  = _prep['transformed_mtb'][CRC_DS].copy()
meta_all = _prep['harmonized_meta'][CRC_DS].copy()
polyamine_cols = _prep['polyamine_columns'][CRC_DS]          # {KEGG_ID: col_name}

print(f"Pickle sample scope: {_prep.get('sample_scope', 'unspecified')}")
print(f'Loaded {len(meta_all)} primary samples.')

# --- Demonstrate that the tumour-location leak guard is live -----------------
# Asking for tumour location on the full primary set (which contains Healthy
# samples) must fail. This is not decoration: it is the check that the rule is
# enforced by code rather than by a comment somebody can ignore.
try:
    utils.build_confounder_matrix(meta_all, include_tumor_location=True)
    raise AssertionError('LEAK GUARD INACTIVE: tumour location was allowed on a '
                         'sample set containing non-cancer samples.')
except ValueError as e:
    print(f'\nLeak guard verified on the full cohort:\n  ValueError: {str(e)[:150]}...')

# --- Restrict to carcinomas, where site IS a legitimate covariate -----------
crc = crc_only_mask(meta_all)
spc, mtb = apply_mask(crc, spc_all, mtb_all, meta_df=meta_all)
meta = meta_all[crc].reset_index(drop=True)

assert len(meta) == int(crc.sum())
assert set(meta['Study.Group']) <= set(CRC_ONLY_STUDY_GROUPS), \
    f'non-carcinoma group in the CRC-only subset: {set(meta["Study.Group"])}'
assert meta['Tumor.Location.Site'].notna().all(), \
    'a carcinoma sample has no tumour site recorded'

# Now the SAME call succeeds, and carries the site dummies.
conf_crc = utils.build_confounder_matrix(meta, include_tumor_location=True)
site_cols = [c for c in conf_crc.columns if c.startswith('TumorSite_')]

pa_cols = [c for c in polyamine_cols.values() if c in mtb.columns]
KEGG_TO_NAME = {v: POLYAMINE_KEGG[k] for k, v in polyamine_cols.items()}

print(f'\nCRC-only analysis set: n={len(meta)}')
print(meta['Study.Group'].value_counts().reindex(CRC_ONLY_STUDY_GROUPS).to_string())
print(f'\nTumour site distribution:')
print(meta['Tumor.Location.Site'].value_counts().reindex(TUMOR_SITE_ORDER).to_string())
print(f'\nConfounder columns now available: {[c for c in conf_crc.columns if c != "Sample"]}')
print(f'Polyamines available: {len(pa_cols)}')'''),

("markdown", "nb04-trend-hdr", r'''## 1. Polyamine trend across cancer stage, adjusted for tumour site

For each polyamine, `CRC.Stage.Numeric` (Stage_0=1, Stage_I_II=2, Stage_III_IV=3 within this
subset) is correlated with abundance:

- `Rho_unadjusted` — plain Spearman, i.e. what the earlier pipeline would have reported.
- `Rho_adjusted` — **partial** Spearman controlling for tumour-site dummies plus Age, BMI,
  Gender, Brinkman index and Alcohol.

Both are BH-corrected across the polyamines tested. The comparison between the two columns
is the result worth reading: a polyamine whose unadjusted trend disappears after adjustment
was tracking site or patient characteristics, not stage.

Note the ordinal coding treats the three stages as equally spaced, which is a modelling
assumption, not a fact about tumour biology. A rank-based correlation only needs the
*ordering* to be right, which it is.'''),

("code", "nb04-trend", r'''# Stage trend per polyamine: unadjusted vs site+confounder-adjusted.
merged = (mtb.merge(meta[['Sample', 'Study.Group', 'CRC.Stage.Numeric',
                          'Tumor.Location.Site']], on='Sample')
             .merge(conf_crc, on='Sample'))

ADJ_COVARS = ['Age', 'BMI', 'Brinkman.Index', 'Alcohol', 'Gender_enc'] + site_cols
usable_covars = [c for c in ADJ_COVARS
                 if c in merged.columns and merged[c].nunique(dropna=True) > 1]
dropped_covars = [c for c in ADJ_COVARS if c not in usable_covars]
if dropped_covars:
    print(f'Constant/absent covariates dropped: {dropped_covars}')
print(f'Adjusting for: {usable_covars}')

rows = []
for col in pa_cols:
    sub = merged[[col, 'CRC.Stage.Numeric'] + usable_covars].dropna()
    if len(sub) < 30 or sub[col].nunique() < 5:
        print(f'  SKIP {KEGG_TO_NAME.get(col, col)}: n={len(sub)} after dropping missing covariates')
        continue

    rho_u, p_u = spearmanr(sub[col], sub['CRC.Stage.Numeric'])

    rho_a, p_a, n_a = np.nan, np.nan, len(sub)
    if HAS_PINGOUIN:
        try:
            r = pg.partial_corr(data=sub, x=col, y='CRC.Stage.Numeric',
                                covar=usable_covars, method='spearman')
            rho_a, p_a = float(r['r'].iloc[0]), float(r['p-val'].iloc[0])
            n_a = int(r['n'].iloc[0])
        except Exception as e:
            print(f'  partial corr failed for {col}: {e}')

    rows.append({'Polyamine': KEGG_TO_NAME.get(col, col), 'Feature': col, 'N': n_a,
                 'Rho_unadjusted': rho_u, 'PValue_unadjusted': p_u,
                 'Rho_adjusted': rho_a, 'PValue_adjusted': p_a})

trend = pd.DataFrame(rows)
if not trend.empty:
    trend['QValue_unadjusted'] = multipletests(trend['PValue_unadjusted'], method='fdr_bh')[1]
    ok = trend['PValue_adjusted'].notna()
    trend['QValue_adjusted'] = np.nan
    if ok.any():
        trend.loc[ok, 'QValue_adjusted'] = multipletests(
            trend.loc[ok, 'PValue_adjusted'], method='fdr_bh')[1]
    trend['Survives_adjustment'] = (trend['QValue_adjusted'] < FDR_THRESHOLD)
    trend['Direction'] = np.where(trend['Rho_adjusted'] > 0,
                                  'increases with stage', 'decreases with stage')
    trend = trend.sort_values('QValue_adjusted', na_position='last')

    print('\n=== Polyamine trend across cancer stage (n=%d carcinomas) ===' % len(meta))
    print(trend[['Polyamine', 'N', 'Rho_unadjusted', 'QValue_unadjusted',
                 'Rho_adjusted', 'QValue_adjusted', 'Survives_adjustment']]
          .round(4).to_string(index=False))

    n_u = int((trend['QValue_unadjusted'] < FDR_THRESHOLD).sum())
    n_a = int(trend['Survives_adjustment'].sum())
    print(f'\nSignificant unadjusted: {n_u} | after site+confounder adjustment: {n_a}')
    lost = trend[(trend['QValue_unadjusted'] < FDR_THRESHOLD) & ~trend['Survives_adjustment']]
    if not lost.empty:
        print('Lost to adjustment (was tracking site or patient characteristics, not stage):')
        for _, r in lost.iterrows():
            print(f"  {r['Polyamine']}: rho {r['Rho_unadjusted']:+.3f} (q={r['QValue_unadjusted']:.3g}) "
                  f"-> {r['Rho_adjusted']:+.3f} (q={r['QValue_adjusted']:.3g})")
    if n_u == 0 and n_a == 0:
        print('No polyamine tracks stage within diagnosed cancer at FDR 0.05. With n=150 '
              'split across three stages this is an underpowered test, so this is an '
              'absence of evidence, not evidence of no stage effect.')

    trend.to_csv(TABLES_DIR / 'nb04_polyamine_stage_trend_v2.csv', index=False)
    print(f"\nSaved: nb04_polyamine_stage_trend_v2.csv")
else:
    trend = pd.DataFrame()
    print('No polyamine could be tested.')'''),

("code", "nb04-trend-fig", r'''# Figure: polyamine distribution by cancer stage, with the adjusted trend annotated.
stage_order  = CRC_ONLY_STUDY_GROUPS
stage_labels = {'Stage_0': 'Stage 0\n(in situ)', 'Stage_I_II': 'Stage I/II',
                'Stage_III_IV': 'Stage III/IV'}

if pa_cols:
    n_pa  = len(pa_cols)
    ncols = min(4, n_pa)
    nrows = (n_pa + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.6 * ncols, 4.0 * nrows), squeeze=False)
    axes_flat = axes.flatten()
    ramp = ['#9ec5f4', '#2a78d6', '#104281']       # sequential: ordered severity

    trend_idx = trend.set_index('Feature') if not trend.empty else pd.DataFrame()

    for i, col in enumerate(pa_cols):
        ax = axes_flat[i]
        plot_df = merged[['Study.Group', col]].dropna()
        plot_df = plot_df[plot_df['Study.Group'].isin(stage_order)]

        sns.boxplot(data=plot_df, x='Study.Group', y=col, order=stage_order,
                    palette=ramp, ax=ax, linewidth=0.8, fliersize=2)
        sns.stripplot(data=plot_df, x='Study.Group', y=col, order=stage_order,
                      color='black', alpha=0.25, size=3, ax=ax)

        meds = plot_df.groupby('Study.Group')[col].median().reindex(stage_order)
        xs = [i2 for i2, s in enumerate(stage_order) if pd.notna(meds.get(s))]
        ys = [meds[s] for s in stage_order if pd.notna(meds.get(s))]
        if len(xs) >= 2:
            ax.plot(xs, ys, '-', color='darkred', lw=1.8, marker='D', ms=5, zorder=4)

        ttl = KEGG_TO_NAME.get(col, col)
        if col in getattr(trend_idx, 'index', []):
            r = trend_idx.loc[col]
            ttl += (f"\nadj. rho={r['Rho_adjusted']:+.2f}, q={r['QValue_adjusted']:.3f}"
                    if pd.notna(r['Rho_adjusted']) else '')
        ax.set_title(ttl, fontsize=9.5, fontweight='bold')
        ax.set_xticks(range(len(stage_order)))
        ax.set_xticklabels([stage_labels[s] for s in stage_order], fontsize=8)
        ax.set_xlabel(''); ax.set_ylabel('log10 abundance', fontsize=8.5)
        for xi, s in enumerate(stage_order):
            ax.text(xi, ax.get_ylim()[0], f'n={int((plot_df["Study.Group"] == s).sum())}',
                    ha='center', va='bottom', fontsize=6, color='grey')

    for j in range(n_pa, len(axes_flat)):
        axes_flat[j].set_visible(False)

    plt.suptitle('Polyamine abundance across diagnosed cancer stage (n=%d carcinomas)\n'
                 'Red line = median trend | annotated rho/q are adjusted for tumour site, '
                 'age, BMI, sex, smoking and alcohol' % len(meta),
                 fontsize=11.5, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    out = FIG_DIR / 'nb04_polyamine_stage_trend_v2.png'
    fig.savefig(out, dpi=150, bbox_inches='tight')
    plt.show(); plt.close()
    print(f'Saved: {out}')'''),

("markdown", "nb04-site-hdr", r'''## 2. Tumour site differences

Do the polyamines differ by where the tumour is? This is both a result in its own right and
the justification for adjusting for site in Section 1: if site matters for a compound, then
an unadjusted stage trend for that compound is confounded by the site mix of each stage
group.

Kruskal–Wallis across Rectum / Left colon / Right colon, BH-corrected across polyamines.
`MultiSite` (n=5) is excluded — it is not a location so much as a statement that there was
more than one.'''),

("code", "nb04-site", r'''SITES_TESTED = ['Rectum', 'Left_colon', 'Right_colon']   # MultiSite (n=5) excluded

site_rows = []
site_df = merged[merged['Tumor.Location.Site'].isin(SITES_TESTED)]
print('Samples per site tested:')
print(site_df['Tumor.Location.Site'].value_counts().reindex(SITES_TESTED).to_string())
print(f"(excluded: {int((merged['Tumor.Location.Site'] == 'MultiSite').sum())} MultiSite)")

for col in pa_cols:
    groups = [g[col].dropna().values for _, g in site_df.groupby('Tumor.Location.Site')
              if len(g[col].dropna()) >= 5]
    if len(groups) < 2:
        continue
    try:
        H, p = kruskal(*groups)
    except ValueError:
        continue
    meds = site_df.groupby('Tumor.Location.Site')[col].median().reindex(SITES_TESTED)
    site_rows.append({'Polyamine': KEGG_TO_NAME.get(col, col), 'Feature': col,
                      'H': H, 'PValue': p,
                      **{f'Median_{s}': meds.get(s, np.nan) for s in SITES_TESTED}})

site_res = pd.DataFrame(site_rows)
if not site_res.empty:
    site_res['QValue'] = multipletests(site_res['PValue'], method='fdr_bh')[1]
    site_res['Significant'] = site_res['QValue'] < FDR_THRESHOLD
    site_res = site_res.sort_values('QValue')
    print('\n=== Polyamine differences by tumour site (Kruskal-Wallis) ===')
    print(site_res[['Polyamine', 'H', 'PValue', 'QValue', 'Significant']
                   + [f'Median_{s}' for s in SITES_TESTED]].round(4).to_string(index=False))
    n_sig = int(site_res['Significant'].sum())
    print(f'\n{n_sig} of {len(site_res)} polyamines differ by tumour site at FDR 0.05.')
    if n_sig:
        print('=> Site is a real covariate for these compounds, so the site-adjusted stage '
              'estimates in Section 1 are the ones to report.')
    else:
        print('=> No polyamine differs detectably by site in this sample. Adjusting for site '
              'was still the correct design choice (it cannot bias the estimate), but it '
              'does not change the conclusions here.')
    site_res.to_csv(TABLES_DIR / 'nb04_polyamine_tumour_site_v2.csv', index=False)
    print('Saved: nb04_polyamine_tumour_site_v2.csv')
else:
    site_res = pd.DataFrame()
    print('No polyamine could be tested across sites.')'''),

("markdown", "nb04-persite-hdr", r'''## 3. Per-site stage trends (exploratory)

The Section 1 stage trend, repeated separately within each tumour site. This asks whether a
stage association is consistent across sites or driven by one of them.

**Labelled exploratory and not corrected across sites**, because the subgroups are small
(Rectum 69, Right colon 44, Left colon 32 — split further across three stages) and the test
is run once per site per polyamine. Nothing downstream consumes these; they exist to show
whether a headline trend is site-specific, which is a caveat on the main result rather than
a result itself.'''),

("code", "nb04-persite", r'''persite_rows = []
for site in SITES_TESTED:
    sub_site = merged[merged['Tumor.Location.Site'] == site]
    n_stages = sub_site['Study.Group'].nunique()
    if len(sub_site) < 20 or n_stages < 2:
        print(f'{site}: n={len(sub_site)}, {n_stages} stage(s) -- too small, skipped')
        continue
    for col in pa_cols:
        s = sub_site[[col, 'CRC.Stage.Numeric']].dropna()
        if len(s) < 15:
            continue
        rho, p = spearmanr(s[col], s['CRC.Stage.Numeric'])
        persite_rows.append({'Site': site, 'N': len(s),
                             'Polyamine': KEGG_TO_NAME.get(col, col), 'Feature': col,
                             'Rho': rho, 'PValue': p})

persite = pd.DataFrame(persite_rows)
if not persite.empty:
    # Nominal p only -- see the header: these are exploratory subgroup tests.
    persite['Nominal_p_lt_05'] = persite['PValue'] < 0.05
    persite = persite.sort_values(['Polyamine', 'Site'])
    print('=== Per-site stage trends (exploratory, nominal p, NOT FDR-corrected) ===')
    print(persite[['Polyamine', 'Site', 'N', 'Rho', 'PValue', 'Nominal_p_lt_05']]
          .round(4).to_string(index=False))

    # Consistency: does a polyamine point the same way in every site?
    cons = (persite.groupby('Polyamine')['Rho']
            .agg(n_sites='size', n_positive=lambda s: int((s > 0).sum()))
            .reset_index())
    cons['Consistent_direction'] = cons['n_positive'].isin([0]) | (cons['n_positive'] == cons['n_sites'])
    print('\nDirection consistency across sites:')
    print(cons.to_string(index=False))
    print('\nA polyamine with Consistent_direction=False points different ways in different '
          'sites, so any pooled stage estimate for it is an average over heterogeneous '
          'subgroups -- read it with that in mind.')

    persite.to_csv(TABLES_DIR / 'nb04_polyamine_stage_trend_per_site_v2.csv', index=False)
    print('\nSaved: nb04_polyamine_stage_trend_per_site_v2.csv')
else:
    persite = pd.DataFrame()
    print('No per-site trend could be computed.')'''),

("markdown", "nb04-species-hdr", r'''## 4. Species stage trend within cancer (site-adjusted)

The species-side counterpart to Section 1: which species track cancer stage once tumour site
and patient covariates are accounted for? BH-corrected across the reduced species set.

This is the within-cancer analogue of NB02's Healthy-vs-CRC contrasts, and it is a harder
test: it asks for a gradient *among* cancers rather than a difference from healthy controls,
with a third of the sample size.'''),

("code", "nb04-species", r'''spc_merged = (spc.merge(meta[['Sample', 'CRC.Stage.Numeric']], on='Sample')
                 .merge(conf_crc, on='Sample'))
spc_feats = [c for c in spc.columns if c != 'Sample']

srows = []
for col in spc_feats:
    s = spc_merged[[col, 'CRC.Stage.Numeric'] + usable_covars].dropna()
    if len(s) < 30 or s[col].nunique() < 5:
        continue
    rho_u, p_u = spearmanr(s[col], s['CRC.Stage.Numeric'])
    rho_a, p_a = np.nan, np.nan
    if HAS_PINGOUIN:
        try:
            r = pg.partial_corr(data=s, x=col, y='CRC.Stage.Numeric',
                                covar=usable_covars, method='spearman')
            rho_a, p_a = float(r['r'].iloc[0]), float(r['p-val'].iloc[0])
        except Exception:
            pass
    srows.append({'Species': col, 'N': len(s),
                  'Rho_unadjusted': rho_u, 'PValue_unadjusted': p_u,
                  'Rho_adjusted': rho_a, 'PValue_adjusted': p_a})

spc_trend = pd.DataFrame(srows)
if not spc_trend.empty:
    spc_trend['QValue_unadjusted'] = multipletests(spc_trend['PValue_unadjusted'],
                                                   method='fdr_bh')[1]
    ok = spc_trend['PValue_adjusted'].notna()
    spc_trend['QValue_adjusted'] = np.nan
    if ok.any():
        spc_trend.loc[ok, 'QValue_adjusted'] = multipletests(
            spc_trend.loc[ok, 'PValue_adjusted'], method='fdr_bh')[1]
    spc_trend['Survives_adjustment'] = spc_trend['QValue_adjusted'] < FDR_THRESHOLD
    spc_trend = spc_trend.sort_values('QValue_adjusted', na_position='last')

    n_u = int((spc_trend['QValue_unadjusted'] < FDR_THRESHOLD).sum())
    n_a = int(spc_trend['Survives_adjustment'].sum())
    print(f'Species tested: {len(spc_trend)}')
    print(f'  track stage unadjusted (q<0.05)          : {n_u}')
    print(f'  survive site+confounder adjustment       : {n_a}')
    if n_a:
        print('\nTop site-adjusted species stage trends:')
        print(spc_trend.head(15)[['Species', 'Rho_adjusted', 'QValue_adjusted',
                                   'Rho_unadjusted', 'QValue_unadjusted']]
              .round(4).to_string(index=False))
    else:
        print('\nNo species tracks stage within diagnosed cancer after adjustment. Given '
              'n=150 across three stages and BH across the full species set, this test has '
              'little power -- reported as a negative result, not as evidence of absence.')

    spc_trend.to_csv(TABLES_DIR / 'nb04_species_stage_trend_v2.csv', index=False)
    print('\nSaved: nb04_species_stage_trend_v2.csv')
else:
    spc_trend = pd.DataFrame()
    print('No species could be tested.')'''),

("markdown", "nb04-summary-hdr", r'''## Summary'''),

("code", "nb04-summary", r'''print('=' * 78)
print('NB04 -- WITHIN-CANCER STAGE & TUMOUR SITE: SUMMARY')
print('=' * 78)
print(f'\nSample set: {len(meta)} carcinomas '
      f'({dict(meta["Study.Group"].value_counts().reindex(CRC_ONLY_STUDY_GROUPS))})')
print(f'Tumour sites: {dict(meta["Tumor.Location.Site"].value_counts())}')
print(f'Adjustment set: {usable_covars}')

print('\n1. Polyamine stage trend (site + confounder adjusted):')
if not trend.empty:
    n_u = int((trend['QValue_unadjusted'] < FDR_THRESHOLD).sum())
    n_a = int(trend['Survives_adjustment'].sum())
    print(f'   {n_u} significant unadjusted -> {n_a} survive adjustment (of {len(trend)} tested)')
else:
    print('   not computed')

print('\n2. Tumour-site differences:')
if not site_res.empty:
    print(f'   {int(site_res["Significant"].sum())} of {len(site_res)} polyamines differ by site (FDR 0.05)')
else:
    print('   not computed')

print('\n3. Per-site stage trends (exploratory):')
if not persite.empty:
    print(f'   {len(persite)} site x polyamine tests across {persite["Site"].nunique()} sites')
else:
    print('   not computed')

print('\n4. Species stage trend (site + confounder adjusted):')
if not spc_trend.empty:
    print(f'   {int(spc_trend["Survives_adjustment"].sum())} of {len(spc_trend)} species '
          f'survive adjustment')
else:
    print('   not computed')

print('\nOutputs:')
for f in ['nb04_polyamine_stage_trend_v2.csv',
          'nb04_polyamine_tumour_site_v2.csv',
          'nb04_polyamine_stage_trend_per_site_v2.csv',
          'nb04_species_stage_trend_v2.csv']:
    print(f'  tables/{f}' + ('  [read by NB08]' if 'stage_trend_v2' in f and 'polyamine' in f else ''))
print('  figures/progression/nb04_polyamine_stage_trend_v2.png')

print('\nCaveats carried forward:')
print('  - n=150 across three stages: underpowered for small effects.')
print('  - Site subgroups are small (69/44/32); per-site results are exploratory.')
print('  - Ordinal stage coding assumes equal spacing of stage ranks.')
print('  - Tumour site is used ONLY here and in NB05 stage 2 / the workbench CRC scope, '
      'where every sample is a carcinoma.')'''),
]

nb = Notebook.create(NEW, CELLS)
nb.save()

if OLD.exists():
    OLD.unlink()
    print(f"removed {OLD.name} (its subject, the MP group, is excluded from v2.1)")
