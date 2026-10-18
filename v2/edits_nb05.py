"""Rebuild NB05: corrected two-stage design + external transfer validation.

The old notebook's stage 1 was a 4-class classifier whose classes included the
excluded MP group, and an earlier version fed it one-hot tumour location, which
separated {Healthy, MP} from {Early, Advanced} exactly and inflated its accuracy.
v2.1 splits the task into the two questions that are actually well posed, and adds
the external-cohort test the pipeline previously had no way to perform.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from nbedit import Notebook

HERE = Path(__file__).resolve().parent

CELLS = [
("markdown", "nb05-title", r'''# 05 — Two-Stage Model & External Transfer Validation

Two questions, kept separate because they need different sample sets and different
covariates:

| | Stage 1 | Stage 2 |
|---|---|---|
| Question | Healthy vs CRC | Early_CRC vs Advanced_CRC |
| Samples | all 277 primary | the 150 carcinomas |
| Tumour site as a feature | **never** | **allowed** |

Then **stage 3**: does a Healthy-vs-CRC model trained on Yachida work on a cohort it has
never seen? That is tested on Sinha and Kim at genus level, with bootstrap confidence
intervals.

And **stage 6**: does the *severity gradient* transfer? Trained on Yachida's full
Healthy → Adenoma → Early → Advanced axis and tested on Kim's Control / Adenoma / Carcinoma.
This stage opts the adenoma groups back in (they are excluded everywhere else) and is the
only place in the pipeline where they appear. Note up front: **neither validation cohort has
a staging column**, so TNM stage cannot be validated outside Yachida at all — 3-level
severity is the finest axis Kim's labels can confirm. See the stage 6 header for why this is
external validation rather than LODO.

## What was wrong before

1. **A 4-class stage classifier including MP.** MP (multiple polypoid adenomas) is excluded
   from v2.1 entirely, so the old label set no longer exists. More importantly, lumping a
   cancer-free adenoma group into a severity classifier mixed two different questions
   ("is there cancer?" and "how advanced is it?") into one accuracy number.
2. **Tumour location among the features.** `Tumor location` is `"-"` for every non-cancer
   sample, so one-hot "NoTumor" separated {Healthy, MP} from {Early, Advanced} *exactly*.
   Removing it dropped stage-1 accuracy from **0.63 to 0.39** and Cohen's κ from **0.48 to
   0.08** — i.e. most of the apparent performance was the leak. That measurement is the
   reason this notebook now splits the task: once the leak is gone, a single multiclass
   severity model has almost no signal, whereas the two binary questions are at least
   well posed.
3. **Two leakage bugs in the regression stage.** The old notebook called
   `scaler.fit_transform(X)` on the full sample set before `cross_val_predict`, and computed
   each target's median-imputation value over the full sample set before CV. Both leak
   held-out information into every fold. Here scaling and confounder imputation happen
   *inside* a `Pipeline` passed whole to `cross_val_predict`, so they are refit on each
   fold's training rows only. Target imputation is dropped entirely rather than made
   fold-safe: `transformed_mtb` has no missing values by construction, so a residual NaN
   after merging is a row-alignment gap, not a value to synthesise — those rows are dropped,
   which cannot leak.
4. **No external validation.** Everything was in-cohort. Stage 3 fixes that.

## Stage 3: why genus level, and what the result means

Yachida is shotgun (species); Sinha and Kim are 16S (genus). A species-level model cannot be
applied across them at all, so the transfer model is trained on **genus** CLR features
restricted to the genera shared with the target cohort (`utils.genus_clr_table`, the same
helper NB02's replication uses).

An AUC near 0.5 on an external cohort is a real and publishable result: it would mean the
Yachida signal does not generalise, which is the single most important thing to know about a
single-cohort microbiome finding. Bootstrap CIs are reported so the reader can see whether
an AUC is distinguishable from chance at all, given n=89 (Sinha) and n=138 (Kim).

Sinha and Kim also differ from Yachida in population (US vs Japan), sequencing (16S vs
shotgun) and metabolomics platform, so a transfer failure cannot be attributed to biology
alone — it is confounded with technology. Stated as a limitation rather than resolved.

**Prerequisites:** `02_analysis.ipynb` (reduced pickle) and `03_ml_benchmarking.ipynb`
(`polyamine_targets`).'''),

("code", "nb05-imports", r'''import re
import sys
import pickle
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr, spearmanr
from sklearn.linear_model import Ridge, LogisticRegression
from sklearn.model_selection import KFold, StratifiedKFold, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.experimental import enable_iterative_imputer  # noqa: F401
from sklearn.impute import IterativeImputer, SimpleImputer
from sklearn.metrics import (classification_report, cohen_kappa_score, confusion_matrix,
                             roc_auc_score, roc_curve, accuracy_score)
import lightgbm as lgb
import shap

warnings.filterwarnings('ignore')

sys.path.insert(0, '.')
import utils
import viz
from utils import (
    RESULTS_DIR, POLYAMINE_KEGG, STAGE_GROUP_ORDER, CRC_ONLY_STAGE_GROUPS,
    DATASET_DISCOVERY, DATASETS_VALIDATION, PREVALENCE_SPECIES,
    build_confounder_matrix, crc_only_mask, primary_mask, apply_mask,
    genus_clr_table, shared_genera,
    # Stage 4 (opt-in severity transfer): the parallel severity scale and the
    # loaders that can admit the adenoma groups.
    SEVERITY_TRANSFER_GROUPS, SEVERITY_TRAIN_MAP,
    SEVERITY_EVAL_MAP, SEVERITY_EVAL_ORDER,
    load_all_datasets, harmonize_metadata,
)

CRC_DS   = DATASET_DISCOVERY
PKL_PREP = RESULTS_DIR / 'intermediate' / 'preprocessed_data_reduced_v2.pkl'
PKL_FULL = RESULTS_DIR / 'intermediate' / 'preprocessed_data_v2.pkl'
PKL_ML   = RESULTS_DIR / 'intermediate' / 'ml_results_v2.pkl'

TABLES_DIR = RESULTS_DIR / 'tables'
FIG_DIR    = RESULTS_DIR / 'figures' / 'ml'
TABLES_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

viz.apply_theme()   # one shared theme: see viz.py (dpi used to be 110 here, 100 elsewhere)

with open(PKL_PREP, 'rb') as f:
    _prep = pickle.load(f)
with open(PKL_ML, 'rb') as f:
    ml_res = pickle.load(f)

spc     = _prep['transformed_species'][CRC_DS].copy()
mtb     = _prep['transformed_mtb'][CRC_DS].copy()
meta    = _prep['harmonized_meta'][CRC_DS].copy()
pa_cols = _prep['polyamine_columns'][CRC_DS]

print(f"Pickle sample scope: {_prep.get('sample_scope', 'unspecified')}")
print(f'Samples: {len(meta)} | Species: {spc.shape[1] - 1} | Metabolites: {mtb.shape[1] - 1}')
print('\nStage.Group distribution:')
print(meta['Stage.Group'].value_counts().reindex(STAGE_GROUP_ORDER).to_string())

assert set(meta['Stage.Group'].dropna()) <= set(STAGE_GROUP_ORDER), \
    'non-primary group reached NB05'
assert not meta['Study.Group'].isin(['MP', 'HS']).any()'''),

("markdown", "nb05-s1-hdr", r'''## Stage 1 — Healthy vs CRC (no tumour location)

A binary LightGBM classifier on species CLR + clinical confounders. Out-of-fold
probabilities from 10-fold stratified CV become a feature for the metabolite models in
Stage 4, so no sample's own label ever informs its own probability.

Confounders are Age, BMI, Brinkman.Index, Alcohol and Gender only.
`utils.build_confounder_matrix()` is called **without** tumour location and would raise if
it were requested on this mixed set — the guard is exercised explicitly below.

LightGBM handles missing feature values natively (it learns a default split direction per
tree, fit fresh on each fold's training rows), so stage 1 needs no imputer and therefore no
global `fit` on the full sample set at all — sidestepping that leakage path rather than
patching it.'''),

("code", "nb05-s1", r'''conf_df = build_confounder_matrix(meta)
conf_feat_cols = [c for c in conf_df.columns if c != 'Sample']
spc_feat_cols  = [c for c in spc.columns if c != 'Sample']

# Prove the leak guard is active on this Healthy+CRC set.
try:
    build_confounder_matrix(meta, include_tumor_location=True)
    raise AssertionError('LEAK GUARD INACTIVE in stage 1')
except ValueError:
    print('Leak guard verified: tumour location refused for the Healthy-vs-CRC task.')
assert not any(c.startswith(('TumorSite_', 'TumorLoc_')) for c in conf_feat_cols)

data_merged = (meta[['Sample', 'Study.Group', 'Stage.Group']]
               .merge(spc, on='Sample').merge(conf_df, on='Sample'))

# Binary label: 1 = any diagnosed carcinoma, 0 = Healthy.
data_merged['is_crc'] = data_merged['Stage.Group'].isin(CRC_ONLY_STAGE_GROUPS).astype(int)
print(f'\nStage 1 labels: Healthy={int((data_merged["is_crc"] == 0).sum())}, '
      f'CRC={int((data_merged["is_crc"] == 1).sum())}')

X_s1 = data_merged[spc_feat_cols + conf_feat_cols].values.astype(float)
y_s1 = data_merged['is_crc'].values
samples_s1 = data_merged['Sample'].values

lgb_bin = lgb.LGBMClassifier(n_estimators=300, learning_rate=0.05, num_leaves=31,
                             class_weight='balanced', random_state=42, n_jobs=2, verbose=-1)
cv10 = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)

oof_p_crc = np.zeros(len(y_s1))
fold_aucs = []
for tr, va in cv10.split(X_s1, y_s1):
    lgb_bin.fit(X_s1[tr], y_s1[tr])
    oof_p_crc[va] = lgb_bin.predict_proba(X_s1[va])[:, 1]
    fold_aucs.append(roc_auc_score(y_s1[va], oof_p_crc[va]))

auc_s1 = roc_auc_score(y_s1, oof_p_crc)
acc_s1 = accuracy_score(y_s1, (oof_p_crc > 0.5).astype(int))
kappa_s1 = cohen_kappa_score(y_s1, (oof_p_crc > 0.5).astype(int))
print(f'\nStage 1 (Healthy vs CRC), out-of-fold:')
print(f'  AUC      : {auc_s1:.3f}  (per-fold {np.mean(fold_aucs):.3f} +/- {np.std(fold_aucs):.3f})')
print(f'  Accuracy : {acc_s1:.3f}   Cohen kappa: {kappa_s1:.3f}')
print('\nClassification report (threshold 0.5):')
print(classification_report(y_s1, (oof_p_crc > 0.5).astype(int),
                            target_names=['Healthy', 'CRC']))

proba_df = pd.DataFrame({'Sample': samples_s1, 'P_CRC': oof_p_crc})'''),

("markdown", "nb05-s2-hdr", r'''## Stage 2 — Early vs Advanced, within cancer (tumour location allowed)

The 150 carcinomas only. Here every sample has a located tumour, so tumour site is a real
clinical covariate rather than a disease-status proxy, and
`utils.build_confounder_matrix(..., include_tumor_location=True)` now succeeds — the same
call that raised in Stage 1.

Two models are compared on identical folds so the contribution of site is visible rather
than assumed:

- **without site** — species + Age/BMI/sex/smoking/alcohol;
- **with site** — the same plus Rectum / Left_colon / Right_colon / MultiSite dummies.

With n=96 vs n=54 this is a small problem, so the AUC difference between the two is
reported with its fold-level spread and is not interpreted as a precise effect.'''),

("code", "nb05-s2", r'''crc = crc_only_mask(meta)
meta_crc = meta[crc].reset_index(drop=True)
(spc_crc,) = apply_mask(crc, spc, meta_df=meta)

# The call that raised in stage 1 succeeds here: every sample is a carcinoma.
conf_crc_site = build_confounder_matrix(meta_crc, include_tumor_location=True)
site_cols = [c for c in conf_crc_site.columns if c.startswith('TumorSite_')]
conf_crc_plain = build_confounder_matrix(meta_crc)
print(f'Stage 2 set: n={len(meta_crc)} carcinomas')
print(f'  tumour-site covariates now permitted: {site_cols}')
print(meta_crc['Stage.Group'].value_counts().reindex(CRC_ONLY_STAGE_GROUPS).to_string())

m2 = (meta_crc[['Sample', 'Stage.Group']]
      .merge(spc_crc, on='Sample')
      .merge(conf_crc_site, on='Sample'))
y_s2 = (m2['Stage.Group'] == 'Advanced_CRC').astype(int).values
plain_cols = [c for c in conf_crc_plain.columns if c != 'Sample']

X_nosite = m2[spc_feat_cols + plain_cols].values.astype(float)
X_site   = m2[spc_feat_cols + plain_cols + site_cols].values.astype(float)

cv5 = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
s2_rows = []
for label, Xm in [('without tumour site', X_nosite), ('with tumour site', X_site)]:
    clf = lgb.LGBMClassifier(n_estimators=300, learning_rate=0.05, num_leaves=15,
                             class_weight='balanced', random_state=42, n_jobs=2, verbose=-1)
    oof = np.zeros(len(y_s2))
    aucs = []
    for tr, va in cv5.split(Xm, y_s2):
        clf.fit(Xm[tr], y_s2[tr])
        oof[va] = clf.predict_proba(Xm[va])[:, 1]
        aucs.append(roc_auc_score(y_s2[va], oof[va]))
    s2_rows.append({'Model': label, 'AUC_oof': roc_auc_score(y_s2, oof),
                    'AUC_fold_mean': float(np.mean(aucs)),
                    'AUC_fold_std': float(np.std(aucs)),
                    'N': len(y_s2), 'N_advanced': int(y_s2.sum())})
    if label == 'with tumour site':
        oof_p_adv = oof.copy()

s2_df = pd.DataFrame(s2_rows)
print('\nStage 2 (Early vs Advanced within cancer):')
print(s2_df.round(3).to_string(index=False))
_d = s2_df.loc[s2_df['Model'] == 'with tumour site', 'AUC_oof'].iloc[0] - \
     s2_df.loc[s2_df['Model'] == 'without tumour site', 'AUC_oof'].iloc[0]
print(f'\nAUC change from adding tumour site: {_d:+.3f} '
      f'(fold sd ~{s2_df["AUC_fold_std"].max():.3f}, so a change smaller than that is noise)')
print('Unlike the old stage-1 leak, any gain here is legitimate: within cancer, site is a '
      'clinical covariate, not a restatement of the label.')
s2_df.to_csv(TABLES_DIR / 'nb05_stage2_early_vs_advanced_v2.csv', index=False)'''),

("markdown", "nb05-s3-hdr", r'''## Stage 3 — External transfer: Yachida → Sinha, Kim

A genus-level Healthy-vs-CRC model trained on Yachida and applied **unchanged** to each
validation cohort. For each target cohort:

1. build genus CLR tables for both cohorts with `utils.genus_clr_table` (the same helper
   NB02's replication uses, so the two steps cannot diverge);
2. intersect the genus sets — the model can only use genera both cohorts observe;
3. fit on **all** of Yachida (no sample from the target cohort is involved in training);
4. predict the target cohort and report AUC with a 2000-sample bootstrap CI.

A within-Yachida cross-validated AUC on the same genus feature set is reported alongside, as
the ceiling: the gap between it and the external AUC is the generalisation loss.

Clinical features are deliberately **not** included: the cohorts do not share them (Kim has
no BMI and only banded age, Sinha has no smoking data), so a transfer model must rest on
taxa alone. Elastic-net logistic regression is used rather than a boosted tree, because with
a few dozen shared genera a linear model transfers more predictably and its coefficients are
inspectable.''')
,
("code", "nb05-s3", r'''# Labels only here; the genus CLR tables are built per source->target pair inside the
# transfer loop below, because the FEATURE BASIS MUST BE FIXED BEFORE THE CLR.
#
# CLR divides by the geometric mean of the columns present, so its reference depends
# on which columns those are. Building each cohort's CLR on its own genus set and
# intersecting the columns afterwards leaves Yachida's values referenced to its own
# prevalent genera and Kim's to ~409 -- two different scales, which biases the
# transfer systematically. utils.shared_genera() picks one basis and
# genus_clr_table(genera=...) applies it to both cohorts before transforming.
with open(PKL_FULL, 'rb') as f:
    _full = pickle.load(f)


def _crc_label(meta_ds):
    """Binary CRC label: Yachida via Stage.Group, validation cohorts via Study.Group."""
    lab = meta_ds[['Sample', 'Study.Group', 'Stage.Group']].copy()
    lab['is_crc'] = np.where(
        lab['Stage.Group'].isin(CRC_ONLY_STAGE_GROUPS), 1,
        np.where(lab['Study.Group'] == 'CRC', 1, 0))
    return lab


genus_labels = {ds: _crc_label(_prep['harmonized_meta'][ds])
                for ds in [CRC_DS] + DATASETS_VALIDATION}
for _ds, _lab in genus_labels.items():
    print(f'{_ds}: n={len(_lab)} | CRC={int(_lab["is_crc"].sum())}, '
          f'Healthy={int((_lab["is_crc"] == 0).sum())}')


def boot_auc_ci(y, p, n_boot=2000, seed=42):
    """Percentile bootstrap CI for AUC, resampling samples with replacement."""
    rng = np.random.default_rng(seed)
    y, p = np.asarray(y), np.asarray(p)
    out = []
    for _ in range(n_boot):
        idx = rng.integers(0, len(y), len(y))
        if len(np.unique(y[idx])) < 2:
            continue
        out.append(roc_auc_score(y[idx], p[idx]))
    if not out:
        return (np.nan, np.nan)
    return (float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5)))


def make_transfer_model():
    return Pipeline([
        ('impute', SimpleImputer(strategy='median')),
        ('scale', StandardScaler()),
        ('clf', LogisticRegression(penalty='elasticnet', solver='saga', l1_ratio=0.5,
                                   C=0.1, max_iter=5000, class_weight='balanced',
                                   random_state=42)),
    ])


transfer_rows = []
roc_curves = {}

src_lab = genus_labels[CRC_DS]

for tgt in DATASETS_VALIDATION:
    tgt_lab = genus_labels[tgt]

    # Prevalent in Yachida (the training cohort), merely PRESENT in the target.
    # The asymmetry is the point: applying a prevalence threshold to the target
    # would let held-out data choose the model's features. Presence is a property
    # of the target's assay coverage, not of its outcome, so it is safe to use.
    shared = shared_genera(
        _full['filtered_species'][CRC_DS],
        [_full['filtered_species'][tgt]],
        prevalence=PREVALENCE_SPECIES,
        train_samples=src_lab['Sample'],
        target_samples=[tgt_lab['Sample']],
    )
    if len(shared) < 5:
        print(f'\n{tgt}: only {len(shared)} shared genera -- transfer not attempted.')
        continue

    # Both cohorts CLR-transformed on the SAME basis, so the values are comparable.
    src_gt = genus_clr_table(_full['filtered_species'][CRC_DS],
                             samples=src_lab['Sample'], genera=shared)
    tgt_gt = genus_clr_table(_full['filtered_species'][tgt],
                             samples=tgt_lab['Sample'], genera=shared)

    Xs = src_gt.merge(src_lab[['Sample', 'is_crc']], on='Sample')
    Xt = tgt_gt.merge(tgt_lab[['Sample', 'is_crc']], on='Sample')
    Xs_m, ys = Xs[shared].values.astype(float), Xs['is_crc'].values
    Xt_m, yt = Xt[shared].values.astype(float), Xt['is_crc'].values

    # Ceiling: within-Yachida CV on the SAME feature set.
    cv_in = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    p_in = cross_val_predict(make_transfer_model(), Xs_m, ys, cv=cv_in,
                             method='predict_proba')[:, 1]
    auc_in = roc_auc_score(ys, p_in)
    ci_in = boot_auc_ci(ys, p_in)

    # Transfer: fit on all of Yachida, predict the target cohort.
    model = make_transfer_model().fit(Xs_m, ys)
    p_out = model.predict_proba(Xt_m)[:, 1]
    auc_out = roc_auc_score(yt, p_out)
    ci_out = boot_auc_ci(yt, p_out)
    roc_curves[tgt] = (yt, p_out, auc_out, ci_out)

    transfer_rows.append({
        'Target_cohort': tgt, 'N_shared_genera': len(shared),
        'N_train': len(ys), 'N_test': len(yt), 'N_test_CRC': int(yt.sum()),
        'AUC_within_Yachida': auc_in,
        'CI_within_low': ci_in[0], 'CI_within_high': ci_in[1],
        'AUC_transfer': auc_out,
        'CI_transfer_low': ci_out[0], 'CI_transfer_high': ci_out[1],
        'Generalisation_loss': auc_in - auc_out,
        'Transfer_beats_chance': ci_out[0] > 0.5,
    })
    print(f'\n--- {CRC_DS} -> {tgt} ---')
    print(f'  shared genera        : {len(shared)}')
    print(f'  within-Yachida (CV)  : AUC {auc_in:.3f}  95% CI [{ci_in[0]:.3f}, {ci_in[1]:.3f}]')
    print(f'  transferred          : AUC {auc_out:.3f}  95% CI [{ci_out[0]:.3f}, {ci_out[1]:.3f}]')
    print(f'  generalisation loss  : {auc_in - auc_out:+.3f}')
    print(f'  better than chance?  : {"YES" if ci_out[0] > 0.5 else "NO (CI includes 0.5)"}')

transfer_df = pd.DataFrame(transfer_rows)
if not transfer_df.empty:
    transfer_df.to_csv(TABLES_DIR / 'nb05_external_transfer_v2.csv', index=False)
    print('\n=== External transfer summary ===')
    print(transfer_df.round(3).to_string(index=False))
    print('\nSaved: nb05_external_transfer_v2.csv')
    n_gen = int(transfer_df['Transfer_beats_chance'].sum())
    print(f'\n{n_gen} of {len(transfer_df)} validation cohorts are predicted better than '
          f'chance by a model trained only on Yachida.')
    print('Caveat: the cohorts differ in population (Japan vs US), sequencing (shotgun vs '
          '16S) and metabolomics platform, so a transfer failure is confounded with '
          'technology and cannot be read as "the biology does not replicate".')'''),

("code", "nb05-s3-fig", r'''# ROC curves for the external transfers.
if roc_curves:
    fig, ax = plt.subplots(figsize=(6.4, 6))
    colors = {'SINHA_CRC_2016': '#eb6834', 'Kim_adenomas_2020': '#2e9e5b'}
    for tgt, (yt, p_out, auc_out, ci_out) in roc_curves.items():
        fpr, tpr, _ = roc_curve(yt, p_out)
        ax.plot(fpr, tpr, lw=2, color=colors.get(tgt, '#777777'),
                label=f'{tgt}\nAUC {auc_out:.3f} [{ci_out[0]:.2f}, {ci_out[1]:.2f}]')
    ax.plot([0, 1], [0, 1], ls='--', color='grey', lw=1, label='Chance')
    ax.set_xlabel('False positive rate'); ax.set_ylabel('True positive rate')
    ax.set_title('External transfer: Healthy vs CRC\n'
                 'genus-level model trained on Yachida only, applied unchanged',
                 fontsize=11, fontweight='bold')
    ax.legend(fontsize=8, loc='lower right')
    plt.tight_layout()
    out = FIG_DIR / 'nb05_external_transfer_roc_v2.png'
    fig.savefig(out, dpi=150, bbox_inches='tight')
    plt.show(); plt.close()
    print(f'Saved: {out}')
else:
    print('No transfer ROC curves to plot.')'''),

("markdown", "nb05-s4-hdr", r'''## Stage 4 — Does stage information help predict polyamine levels?

The notebook's original question, now asked with leak-free inputs. For each pre-specified
polyamine (NB03's `polyamine_targets`), Ridge regression predicts its abundance from:

- **flat**: species + confounders;
- **two-stage**: the same plus the out-of-fold `P_CRC` from Stage 1 (and, for carcinomas,
  `P_advanced` from Stage 2 where available).

Both via 10-fold `cross_val_predict` on an identical fold structure, scored by Pearson r.
`Delta_R` is the gain from adding stage information.

Because the stage probabilities are out-of-fold, a sample's own disease label never informs
its own features. If `Delta_R` is near zero, the honest reading is that knowing disease
status adds nothing beyond what the species composition already carries — which is itself
informative about how much of the "stage-specific metabolome" story is really just
composition.'''),

("code", "nb05-s4", r'''pa_targets = ml_res.get('polyamine_targets', pd.DataFrame())
target_pa_cols = [c for c in pa_targets.get('Feature', pd.Series(dtype=str)).tolist()
                  if c in mtb.columns]
if not target_pa_cols:
    target_pa_cols = [v for v in pa_cols.values() if v in mtb.columns]
print(f'Target polyamines ({len(target_pa_cols)}): '
      f'{[k for k, v in pa_cols.items() if v in target_pa_cols]}')

# Stage-2 probability is defined only for carcinomas; Healthy samples get 0 and an
# explicit indicator, rather than a silently imputed value that the model could
# mistake for a real probability.
adv_df = pd.DataFrame({'Sample': m2['Sample'].values, 'P_advanced': oof_p_adv})
base_df = (data_merged[['Sample'] + spc_feat_cols + conf_feat_cols]
           .merge(proba_df, on='Sample')
           .merge(adv_df, on='Sample', how='left'))
base_df['Has_stage2'] = base_df['P_advanced'].notna().astype(float)
base_df['P_advanced'] = base_df['P_advanced'].fillna(0.0)
stage_cols = ['P_CRC', 'P_advanced', 'Has_stage2']
assert len(base_df) == len(data_merged), 'stage-probability merge changed the row count'


def build_pa_xy(pa_col):
    """Row-aligned two-stage/flat matrices + target for one polyamine.

    Rows with a missing target are dropped rather than imputed -- the old
    notebook's fold-unsafe global median fill was the second leakage bug; dropping
    is leakage-proof by construction.
    """
    pa_df = base_df.merge(mtb[['Sample', pa_col]], on='Sample').dropna(subset=[pa_col])
    X_ts   = pa_df[spc_feat_cols + stage_cols + conf_feat_cols].values.astype(float)
    X_flat = pa_df[spc_feat_cols + conf_feat_cols].values.astype(float)
    return X_ts, X_flat, pa_df[pa_col].values.astype(float), len(pa_df)


def make_ridge_pipeline():
    """Fresh pipeline per model: IterativeImputer -> StandardScaler -> Ridge.

    cross_val_predict clones it per split and fits the clone on that split's
    training rows only, so no imputer/scaler statistic is ever computed from a
    held-out sample. skip_complete=True so the imputer models only the few
    confounder columns that actually have NaN, not all ~500 complete columns.
    """
    return Pipeline([
        ('imputer', IterativeImputer(max_iter=10, random_state=42, skip_complete=True)),
        ('scaler', StandardScaler()),
        ('ridge', Ridge(alpha=1.0)),
    ])


results_stage2 = []
cv10_reg = KFold(n_splits=10, shuffle=True, random_state=42)

for pa_col in target_pa_cols:
    X_twostage, X_flat, y_pa, n_used = build_pa_xy(pa_col)
    n_dropped = len(base_df) - n_used
    if n_dropped:
        print(f'  {pa_col}: {n_dropped} sample(s) dropped (missing target)')
    if len(y_pa) < 20:
        print(f'  {pa_col}: only {len(y_pa)} usable samples -- skipping')
        continue

    pa_kegg = {v: k for k, v in pa_cols.items()}.get(pa_col, pa_col)
    pa_name = POLYAMINE_KEGG.get(pa_kegg, pa_kegg)

    r_ts = pearsonr(y_pa, cross_val_predict(make_ridge_pipeline(), X_twostage, y_pa,
                                            cv=cv10_reg))[0]
    r_flat = pearsonr(y_pa, cross_val_predict(make_ridge_pipeline(), X_flat, y_pa,
                                              cv=cv10_reg))[0]

    results_stage2.append({'Polyamine': pa_name, 'Feature': pa_col, 'N': n_used,
                           'R_two_stage': round(r_ts, 3), 'R_flat': round(r_flat, 3),
                           'Delta_R': round(r_ts - r_flat, 3)})
    print(f'  {pa_name:12s}: two-stage R={r_ts:+.3f}, flat R={r_flat:+.3f}, '
          f'Delta_R={r_ts - r_flat:+.3f}')

results_df = pd.DataFrame(results_stage2)
print('\n=== Stage 4: two-stage vs flat ===')
print(results_df.to_string(index=False) if not results_df.empty else '(no results)')
if not results_df.empty:
    md = results_df['Delta_R'].median()
    print(f'\nMedian Delta_R = {md:+.3f}')
    if abs(md) < 0.02:
        print('=> Stage information adds essentially nothing beyond species composition. '
              'The composition already carries whatever disease signal the stage label '
              'encodes, so a "two-stage" architecture is not justified for these targets.')
results_df.to_csv(TABLES_DIR / 'two_stage_model_results_v2.csv', index=False)
print('\nSaved: two_stage_model_results_v2.csv')'''),

("markdown", "nb05-s5-hdr", r'''## Stage 5 — SHAP attribution per polyamine

One SHAP table and figure per target polyamine. Each polyamine's final model is refit on all
of its own valid rows purely for attribution; every performance number above comes
exclusively from the per-fold `cross_val_predict` pipelines, never from this full-data refit.

The old notebook computed SHAP inside the loop for every polyamine but then only ever ran
the explanation step on `target_pa_cols[0]` — an undocumented restriction to a single
metabolite. This runs for all of them.'''),

("code", "nb05-s5", r'''def _pa_slug(name):
    return re.sub(r'[^a-z0-9]+', '_', str(name).lower()).strip('_')


feat_names_ts = spc_feat_cols + stage_cols + conf_feat_cols
shap_summary_rows = []

for pa_col in target_pa_cols:
    X_twostage, _X_flat, y_pa, n_used = build_pa_xy(pa_col)
    if len(y_pa) < 20:
        continue

    pa_kegg = {v: k for k, v in pa_cols.items()}.get(pa_col, pa_col)
    pa_name = POLYAMINE_KEGG.get(pa_kegg, pa_kegg)
    slug = _pa_slug(pa_name)

    pipe_final = make_ridge_pipeline().fit(X_twostage, y_pa)
    X_pre = pipe_final[:-1].transform(X_twostage)
    ridge_final = pipe_final.named_steps['ridge']

    explainer = shap.LinearExplainer(ridge_final, X_pre)
    shap_values = explainer.shap_values(X_pre)
    shap_df = (pd.DataFrame({'Feature': feat_names_ts,
                             'Mean_SHAP': np.abs(shap_values).mean(axis=0)})
               .sort_values('Mean_SHAP', ascending=False).reset_index(drop=True))

    top20 = shap_df.head(20)
    colors = ['#e74c3c' if f in stage_cols else '#3498db' if f in conf_feat_cols
              else '#2ecc71' for f in top20['Feature']]
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.barh(top20['Feature'][::-1], top20['Mean_SHAP'][::-1], color=colors[::-1])
    ax.set_xlabel('Mean |SHAP value|')
    ax.set_title(f'Stage 4 SHAP — {pa_name} (two-stage model, n={n_used})\n'
                 'red = stage probability, blue = confounder, green = species')
    plt.tight_layout()
    out_fig = FIG_DIR / f'two_stage_shap_{slug}_v2.png'
    fig.savefig(out_fig, dpi=150, bbox_inches='tight')
    plt.show(); plt.close(fig)

    shap_df.to_csv(TABLES_DIR / f'two_stage_shap_{slug}_v2.csv', index=False)
    stage_rank = shap_df.index[shap_df['Feature'].isin(stage_cols)].min()
    shap_summary_rows.append({
        'Polyamine': pa_name, 'Feature': pa_col, 'N': n_used,
        'Top_feature': shap_df.iloc[0]['Feature'],
        'Top_Mean_SHAP': round(float(shap_df.iloc[0]['Mean_SHAP']), 4),
        'Best_stage_feature_rank': int(stage_rank) + 1 if pd.notna(stage_rank) else np.nan,
    })
    print(f'  {pa_name}: saved two_stage_shap_{slug}_v2.{{csv,png}}')

shap_summary_df = pd.DataFrame(shap_summary_rows)
print(f'\nSHAP computed for {len(shap_summary_df)}/{len(target_pa_cols)} target polyamine(s)')
if not shap_summary_df.empty:
    print(shap_summary_df.to_string(index=False))
    print('\n"Best_stage_feature_rank" is where the highest-ranked stage probability sits '
          'among all features by mean |SHAP|. A high rank number means the model barely '
          'uses disease status, consistent with Delta_R above.')'''),

("markdown", "nb05-s6-hdr", r'''## Stage 6 — Cross-cohort **severity** transfer (opt-in)

Stage 3 asked "is there cancer?" across cohorts. This asks the harder question: can the
*severity gradient* learned on Yachida be recovered in a cohort it has never seen?

### What this can and cannot test

**Neither validation cohort has a staging column.** Sinha's `Study.Group` is a bare `0`/`1`;
Kim's is `Control` / `Adenoma` / `Carcinoma`. There is no TNM field in either. So a predicted
"Stage_I_II" for a Kim sample is **unfalsifiable** — the ground truth does not exist. Any
claim about predicting *stage* outside Yachida would be untestable by construction.

What *is* testable is one level coarser, and it is a real result:

| | Yachida (train, 4 levels) | Kim (test, 3 levels) | Sinha (test) |
|---|---|---|---|
| Healthy | 127 | Control 102 | Healthy 89 |
| **Adenoma** | **MP 40** | **Adenoma 102** | *none* |
| Carcinoma | Stage_0 27 + I_II 69 + III_IV 54 | Carcinoma 36 | CRC 42 |

Yachida's `MP` is the adenoma group — multiple polypoid adenomas with low-grade dysplasia,
and its raw `Stage` column is `"-"` for all 40, i.e. no carcinoma. That is what makes it the
counterpart of Kim's `Adenoma` arm, which at n=102 is *larger* than Yachida's. The model is
trained on the full 4-level gradient and its `Early_CRC`/`Advanced_CRC` output is collapsed to
`Carcinoma` only at scoring time, because that is the finest distinction Kim can confirm.

**Sinha supports the binary contrast only** — no adenoma arm, no stages. It is reported as a
monotonicity check (is the severity score higher in CRC?), never as staging.

### Why this is not LODO

Leave-one-dataset-out is **degenerate** for the ordinal task: only Yachida has an
adenoma-plus-carcinoma gradient, so holding it out leaves no staged training data at all.
This is single-direction external validation and is reported as such. Binary LODO across all
three cohorts is valid and is what stage 3 already does.

### Design choices

- **Yachida → Kim only.** Measured shared genera: Yachida ∩ Kim = 408, Yachida ∩ Sinha = 60,
  all three = 53. Including Sinha would cap the model at 53 features — Sinha measures only 62
  genera in total. Keeping it out buys roughly 8× the feature space.
- **Opt-in group admission.** Both adenoma arms are normally dropped at load
  (`utils.EXCLUDED_STUDY_GROUPS`). This stage re-admits them with `exclude={}` *locally*;
  `utils.STAGE_GROUP_MAP` is untouched, so the primary analysis set stays n=277 and the
  group-policy checks keep passing.
- **Ordinal model: a regularised linear model on the ordinal code**, giving a continuous
  severity score, evaluated with rank statistics. Not a cumulative-link model: with a few
  hundred genus features a proportional-odds fit is unstable, and stage 3 already established
  that a linear model transfers more predictably and stays inspectable. The score is ordinal,
  so it is judged by *ranking* (Spearman, concordance, step-wise AUROC), not by accuracy.
- **Class thresholds come from the training cohort**, as midpoints between Yachida's
  per-class mean scores. Choosing cut-points on Kim's own score distribution would let the
  held-out data grade its own exam.
- **No clinical covariates**, for the reason the audit established: Kim has no BMI column at
  all and only banded age, Sinha has no smoking data. Only Age and Gender exist in all three,
  which is too thin to matter here, so the model rests on taxa alone.'''),

("code", "nb05-s6", r'''# ---------------------------------------------------------------------------
# Opt-in load: the adenoma groups (Yachida MP, Kim Adenoma) are normally dropped.
# exclude={} must be passed to BOTH loaders -- harmonize_metadata() applies the
# exclusion again, so passing it only to load_all_datasets() silently loses them.
# Scoped to this stage; nothing above or below sees these rows.
# ---------------------------------------------------------------------------
SEV_TRAIN_DS = CRC_DS                    # Yachida: the only staged cohort
SEV_TEST_DS  = 'Kim_adenomas_2020'       # the only cohort with an adenoma arm to test

# MEMORY: this re-reads the raw species tables, and Yachida's is ~105 MB / 57k
# columns, into a kernel already holding _prep, _full and ml_res. It is unavoidable
# here -- NB01's pickle excludes the adenoma groups by construction, so the opt-in
# rows exist only in the raw TSVs. Mitigations: only the two cohorts this stage needs
# are loaded, the mtb/mtb.map tables are dropped immediately, and `_sev_raw` is
# released as soon as the frames are extracted. On the 1500 MB server slice, run this
# stage in a fresh kernel if it is tight.
import gc

_sev_raw = load_all_datasets(datasets=[SEV_TRAIN_DS, SEV_TEST_DS], exclude={})
sev_meta, sev_species = {}, {}
for ds in (SEV_TRAIN_DS, SEV_TEST_DS):
    sev_meta[ds] = harmonize_metadata(_sev_raw[ds]['metadata'], ds, exclude={})
    sev_species[ds] = _sev_raw[ds]['species']
del _sev_raw          # drops the unused mtb / mtb.map frames too
gc.collect()

# Ordinal labels. Train on the 4-level scale, score against the 3-level one.
for ds in (SEV_TRAIN_DS, SEV_TEST_DS):
    m = sev_meta[ds]
    m['Sev_train'] = m['Study.Group'].map(SEVERITY_TRAIN_MAP)
    m['Sev_eval']  = m['Study.Group'].map(SEVERITY_EVAL_MAP)

train_meta = sev_meta[SEV_TRAIN_DS]
test_meta  = sev_meta[SEV_TEST_DS]

# Drop anything off the severity axis (HS/PostSurgery is NaN in both maps by design).
train_meta = train_meta[train_meta['Sev_train'].notna()].reset_index(drop=True)
test_meta  = test_meta[test_meta['Sev_eval'].notna()].reset_index(drop=True)

print('Severity-transfer sample sets')
print(f'  train {SEV_TRAIN_DS}:')
for g, n in train_meta['Study.Group'].value_counts().items():
    print(f'      {g:14s} n={n:4d}  -> level {int(train_meta.loc[train_meta["Study.Group"] == g, "Sev_train"].iloc[0])}')
print(f'  test  {SEV_TEST_DS}:')
for g, n in test_meta['Study.Group'].value_counts().items():
    print(f'      {g:14s} n={n:4d}  -> level {int(test_meta.loc[test_meta["Study.Group"] == g, "Sev_eval"].iloc[0])}')

# --- Guards: these are the failure modes that would silently fake a result. ---
assert not train_meta['Study.Group'].isin(['HS', 'PostSurgery']).any(), \
    'PostSurgery leaked into the severity training set'
assert not test_meta['Study.Group'].isin(['HS', 'PostSurgery']).any(), \
    'PostSurgery leaked into the severity test set'
_overlap = set(train_meta['Sample']) & set(test_meta['Sample'])
assert not _overlap, f'train/test share {len(_overlap)} Sample IDs'
# The opt-in admission must match the DECLARED policy, not merely be whatever
# exclude={} happened to let through. This is what makes SEVERITY_TRANSFER_GROUPS
# an enforced policy rather than a comment.
for _ds, _grps in SEVERITY_TRANSFER_GROUPS.items():
    if _ds not in (SEV_TRAIN_DS, SEV_TEST_DS):
        continue
    _present = set(sev_meta[_ds]['Study.Group'].dropna())
    assert set(_grps) <= _present, (
        f'{_ds}: declared severity-transfer group(s) {_grps} are absent after load '
        f'(present: {sorted(_present)}). exclude={{}} did not admit them.')
assert train_meta['Sev_train'].nunique() == 4, \
    f"expected 4 training levels, got {sorted(train_meta['Sev_train'].unique())}"
assert test_meta['Sev_eval'].nunique() == 3, \
    f"expected 3 test levels, got {sorted(test_meta['Sev_eval'].unique())}"
print(f'\n  guards passed: no PostSurgery, no shared samples, 4 train / 3 test levels')'''),

("code", "nb05-s6-fit", r'''# Shared genus basis, fixed BEFORE the CLR (see utils.genus_clr_table).
sev_shared = shared_genera(
    sev_species[SEV_TRAIN_DS], [sev_species[SEV_TEST_DS]],
    prevalence=PREVALENCE_SPECIES,
    train_samples=train_meta['Sample'],
    target_samples=[test_meta['Sample']],
)
print(f'shared genera (prevalent in {SEV_TRAIN_DS}, present in {SEV_TEST_DS}): {len(sev_shared)}')

sev_X_train = genus_clr_table(sev_species[SEV_TRAIN_DS],
                              samples=train_meta['Sample'], genera=sev_shared)
sev_X_test  = genus_clr_table(sev_species[SEV_TEST_DS],
                              samples=test_meta['Sample'],  genera=sev_shared)
assert list(sev_X_train.columns) == list(sev_X_test.columns), \
    'genus_clr_table returned different column orders for the two cohorts'

_tr = sev_X_train.merge(train_meta[['Sample', 'Sev_train', 'Sev_eval']], on='Sample')
_te = sev_X_test.merge(test_meta[['Sample', 'Sev_eval']], on='Sample')
Xtr = _tr[sev_shared].to_numpy(float)
ytr = _tr['Sev_train'].to_numpy(float)       # 4-level ordinal code
Xte = _te[sev_shared].to_numpy(float)
yte = _te['Sev_eval'].to_numpy(int)          # 3-level ordinal truth


def make_severity_model():
    """Ordinal severity score: regularised linear regression on the ordinal code.

    The same IMPUTE -> SCALE -> penalised-linear shape as stage 3's transfer model
    (deliberately, so the two stages' modelling choices stay aligned), with a Ridge
    head on the ordinal code instead of a logistic head on a binary label.

    But NOTE the feature bases are NOT identical, and the two stages' numbers are
    therefore not directly comparable: stage 3 starts from NB01's species-level
    prevalence-filtered table (`_full['filtered_species']`), whereas this stage must
    start from the RAW species tables, because NB01's pickle excludes the adenoma
    groups by construction. Both then apply the same genus-level prevalence filter
    via `shared_genera`, so stage 6's basis is genus-prevalence-on-raw while stage
    3's is genus-prevalence-on-species-filtered. Read each stage against its own
    within-Yachida ceiling, not against the other stage.

    A continuous score is the point -- it is judged by RANK, not by accuracy.
    """
    return Pipeline([
        ('impute', SimpleImputer(strategy='median')),
        ('scale',  StandardScaler()),
        ('reg',    Ridge(alpha=10.0, random_state=42)),
    ])


def c_index(truth, score):
    """Concordance: over all pairs with DIFFERENT truth levels, the fraction the
    score orders correctly (ties count a half)."""
    truth, score = np.asarray(truth, float), np.asarray(score, float)
    conc = disc = ties = 0
    for i in range(len(truth)):
        hi = truth > truth[i]
        if not hi.any():
            continue
        conc += int(np.sum(score[hi] > score[i]))
        disc += int(np.sum(score[hi] < score[i]))
        ties += int(np.sum(score[hi] == score[i]))
    tot = conc + disc + ties
    return np.nan if tot == 0 else (conc + 0.5 * ties) / tot


# Ceiling: within-Yachida CV on the same feature set, scored on the 3-level axis so
# it is directly comparable with the transferred number.
_cv = KFold(n_splits=5, shuffle=True, random_state=42)
score_in = cross_val_predict(make_severity_model(), Xtr, ytr, cv=_cv)
ytr_eval = _tr['Sev_eval'].to_numpy(int)
rho_in = spearmanr(score_in, ytr_eval).statistic
cidx_in  = c_index(ytr_eval, score_in)

# Transfer: fit on ALL of Yachida, predict Kim. No Kim sample touches the fit.
sev_model = make_severity_model().fit(Xtr, ytr)
score_out = sev_model.predict(Xte)
rho_out = spearmanr(score_out, yte).statistic
cidx_out  = c_index(yte, score_out)

# Step-wise AUROCs on the continuous score -- no thresholds needed. Reported
# SEPARATELY because a single pooled number hides which step fails, and the
# Adenoma-vs-Carcinoma step is where cross-cohort transfer usually breaks.
def _step_auc(lo, hi):
    m = np.isin(yte, [lo, hi])
    if len(np.unique(yte[m])) < 2:
        return np.nan
    return roc_auc_score((yte[m] == hi).astype(int), score_out[m])

auc_h_vs_a = _step_auc(0, 1)   # Control vs Adenoma
auc_a_vs_c = _step_auc(1, 2)   # Adenoma vs Carcinoma
auc_h_vs_c = _step_auc(0, 2)   # Control vs Carcinoma (the easy step, for contrast)

# Class assignment: cut-points from the TRAINING score distribution only.
# Using Kim's own distribution would let the held-out data grade its own exam.
_tr_means = [score_in[ytr_eval == k].mean() for k in (0, 1, 2)]
_cuts = [(_tr_means[0] + _tr_means[1]) / 2, (_tr_means[1] + _tr_means[2]) / 2]
pred_cls = np.digitize(score_out, _cuts)
cm = confusion_matrix(yte, pred_cls, labels=[0, 1, 2])

print(f'\n--- Severity transfer: {SEV_TRAIN_DS} -> {SEV_TEST_DS} ---')
print(f'  shared genera          : {len(sev_shared)}')
print(f'  within-Yachida (CV)    : rho {rho_in:+.3f}   C-index {cidx_in:.3f}')
print(f'  transferred to Kim     : rho {rho_out:+.3f}   C-index {cidx_out:.3f}')
print(f'  generalisation loss    : rho {rho_in - rho_out:+.3f}')
print(f'\n  step-wise AUROC (continuous score, no thresholds):')
print(f'    Control  vs Adenoma  : {auc_h_vs_a:.3f}')
print(f'    Adenoma  vs Carcinoma: {auc_a_vs_c:.3f}   <- the step that usually fails')
print(f'    Control  vs Carcinoma: {auc_h_vs_c:.3f}')
print(f'\n  confusion matrix (rows = truth {SEVERITY_EVAL_ORDER}, cols = predicted):')
print(pd.DataFrame(cm, index=SEVERITY_EVAL_ORDER,
                   columns=SEVERITY_EVAL_ORDER).to_string())'''),

("code", "nb05-s6-sanity", r'''# ---------------------------------------------------------------------------
# Negative controls. A transfer result is only believable if these behave.
# ---------------------------------------------------------------------------

# 1. Label permutation: destroys the signal, so rho must collapse to ~0. If it does
#    not, the "signal" is an artefact of the pipeline rather than of the biology.
_rng = np.random.default_rng(42)
_perm_rho = []
for _ in range(20):
    _yp = _rng.permutation(ytr)
    _s = make_severity_model().fit(Xtr, _yp).predict(Xte)
    _perm_rho.append(spearmanr(_s, yte).statistic)
_perm_rho = np.asarray(_perm_rho, float)
print('Negative control 1 -- permuted training labels')
print(f'  rho over 20 permutations: mean {np.nanmean(_perm_rho):+.3f}, '
      f'max |rho| {np.nanmax(np.abs(_perm_rho)):.3f}   (real: {rho_out:+.3f})')
_passes_perm = abs(rho_out) > np.nanmax(np.abs(_perm_rho))
print(f'  real result exceeds every permutation: {"YES" if _passes_perm else "NO -- treat as null"}')

# 2. Cohort separability. If a classifier can tell Yachida from Kim almost perfectly
#    on these same genera, then a "severity" score can ride on the batch difference
#    instead of on disease. This is the single most important diagnostic here: the
#    two cohorts differ in sequencing platform (shotgun vs 16S) as well as in
#    population, so some separability is expected -- near-1.0 invalidates the result.
_Xc = np.vstack([Xtr, Xte])
_yc = np.r_[np.zeros(len(Xtr)), np.ones(len(Xte))]
_pc = cross_val_predict(
    Pipeline([('impute', SimpleImputer(strategy='median')),
              ('scale', StandardScaler()),
              ('clf', LogisticRegression(max_iter=5000, class_weight='balanced'))]),
    _Xc, _yc, cv=StratifiedKFold(5, shuffle=True, random_state=42),
    method='predict_proba')[:, 1]
_auc_cohort = roc_auc_score(_yc, _pc)
print('\nNegative control 2 -- cohort separability on the shared genus space')
print(f'  AUC(Yachida vs Kim) = {_auc_cohort:.3f}')
if _auc_cohort > 0.95:
    print('  WARNING: the cohorts are almost perfectly separable, so the severity score')
    print('  may be tracking platform/batch rather than disease. Interpret with care and')
    print('  report this number alongside the transfer result.')

sev_rows = [{
    'Train_cohort': SEV_TRAIN_DS, 'Test_cohort': SEV_TEST_DS,
    'N_shared_genera': len(sev_shared),
    'N_train': len(ytr), 'N_test': len(yte),
    'N_test_healthy':   int((yte == 0).sum()),
    'N_test_adenoma':   int((yte == 1).sum()),
    'N_test_carcinoma': int((yte == 2).sum()),
    'Rho_within_train': rho_in,  'Cindex_within_train': cidx_in,
    'Rho_transfer': rho_out,     'Cindex_transfer': cidx_out,
    'Generalisation_loss_rho': rho_in - rho_out,
    'AUC_healthy_vs_adenoma': auc_h_vs_a,
    'AUC_adenoma_vs_carcinoma': auc_a_vs_c,
    'AUC_healthy_vs_carcinoma': auc_h_vs_c,
    'Perm_rho_max_abs': float(np.nanmax(np.abs(_perm_rho))),
    'Beats_permutation': bool(_passes_perm),
    'AUC_cohort_separability': _auc_cohort,
}]
sev_df = pd.DataFrame(sev_rows)
sev_df.to_csv(TABLES_DIR / 'nb05_severity_transfer_v2.csv', index=False)
print('\nSaved: nb05_severity_transfer_v2.csv')

# --- Sinha: binary monotonicity only. It has no adenoma arm and no stages, so it
# --- cannot test staging; applying the same score and asking only whether it is
# --- higher in CRC is the strongest honest claim available.
print('\nSinha is NOT scored for severity: no adenoma arm, no staging column.')
print('Its binary Healthy-vs-CRC transfer is stage 3 above; adding a severity number')
print('for it would be reporting a metric its labels cannot support.')'''),

("code", "nb05-s6-fig", r'''# Left panel: severity score by true group. The groups form an ORDERED axis, so the
# sequential ramp (viz.severity_colors) is the right encoding -- not a categorical
# palette, which would throw away the ordering.
fig, axes = viz.panel_grid(2, ncols=2, panel_size=(5.2, 4.0))

_plot = pd.DataFrame({'score': score_out, 'level': yte})
_plot['group'] = _plot['level'].map(dict(enumerate(SEVERITY_EVAL_ORDER)))
_cols = viz.severity_colors(order=SEVERITY_EVAL_ORDER)

sns.boxplot(data=_plot, x='group', y='score', order=SEVERITY_EVAL_ORDER,
            palette=_cols, ax=axes[0], showfliers=False, width=0.6)
sns.stripplot(data=_plot, x='group', y='score', order=SEVERITY_EVAL_ORDER,
              ax=axes[0], color=viz.INK2, size=3, alpha=0.5, jitter=0.25)
for c in _cuts:
    axes[0].axhline(c, color=viz.MUTED, ls='--', lw=1)
axes[0].set_title(f'Transferred severity score\n{SEV_TRAIN_DS} -> {SEV_TEST_DS}', fontsize=10)
axes[0].set_xlabel('')
axes[0].set_ylabel('predicted severity score')
viz.despine(axes[0])

# Step-wise AUROCs: the honest summary, one bar per step.
_steps = ['Control\nvs Adenoma', 'Adenoma\nvs Carcinoma', 'Control\nvs Carcinoma']
_vals  = [auc_h_vs_a, auc_a_vs_c, auc_h_vs_c]
# ONE series of three independent comparisons -- not an ordered axis, so this must
# NOT use the severity ramp (that would imply the three bars are ordinal steps).
# A single series needs no legend; the tick labels name each bar.
axes[1].bar(_steps, _vals, color=viz.PALETTE_COHORT[0], width=0.6)
axes[1].axhline(0.5, color=viz.MUTED, ls='--', lw=1)
axes[1].set_ylim(0, 1)
axes[1].set_ylabel('AUROC')
axes[1].set_title('Step-wise discrimination\n(0.5 = chance)', fontsize=10)
for i, v in enumerate(_vals):
    if not np.isnan(v):
        axes[1].text(i, v + 0.02, f'{v:.2f}', ha='center', color=viz.INK, fontsize=9)
viz.despine(axes[1])

fig.suptitle('Stage 6 — cross-cohort severity transfer', fontsize=12, x=0.02, ha='left')
viz.save_figure(fig, FIG_DIR / 'nb05_severity_transfer_v2.png')
print('Saved: figures/ml/nb05_severity_transfer_v2.png (+ .pdf)')'''),

("markdown", "nb05-summary-hdr", r'''## Summary'''),

("code", "nb05-summary", r'''print('=' * 78)
print('NB05 -- TWO-STAGE MODEL & EXTERNAL TRANSFER: SUMMARY (v2.1)')
print('=' * 78)

print(f'\nStage 1 -- Healthy vs CRC (n={len(y_s1)}, no tumour location):')
print(f'  OOF AUC {auc_s1:.3f} | accuracy {acc_s1:.3f} | kappa {kappa_s1:.3f}')
print('  For reference, the pre-fix 4-class model that included one-hot tumour location')
print('  reported accuracy 0.63 / kappa 0.48; removing the leak gave 0.39 / 0.08. This')
print('  notebook does not report that multiclass number at all -- the task was split')
print('  into the two binary questions that are actually identifiable.')

print(f'\nStage 2 -- Early vs Advanced within cancer (n={len(y_s2)}):')
print(s2_df.round(3).to_string(index=False))

print('\nStage 3 -- External transfer (genus level, trained on Yachida only):')
if not transfer_df.empty:
    print(transfer_df[['Target_cohort', 'N_shared_genera', 'N_test',
                       'AUC_within_Yachida', 'AUC_transfer',
                       'CI_transfer_low', 'CI_transfer_high',
                       'Transfer_beats_chance']].round(3).to_string(index=False))
else:
    print('  not computed')

print('\nStage 4 -- Does stage information help predict polyamines?')
if not results_df.empty:
    print(results_df.to_string(index=False))
    print(f'  median Delta_R = {results_df["Delta_R"].median():+.3f}')
else:
    print('  no usable targets')

print('\nStage 6 -- Cross-cohort SEVERITY transfer (Yachida -> Kim, opt-in adenoma groups):')
if not sev_df.empty:
    _r = sev_df.iloc[0]
    print(f'  shared genera {int(_r["N_shared_genera"])} | '
          f'rho within-train {_r["Rho_within_train"]:+.3f} -> transferred {_r["Rho_transfer"]:+.3f}')
    print(f'  step-wise AUROC: Control-vs-Adenoma {_r["AUC_healthy_vs_adenoma"]:.3f}, '
          f'Adenoma-vs-Carcinoma {_r["AUC_adenoma_vs_carcinoma"]:.3f}')
    print(f'  beats label permutation: {_r["Beats_permutation"]} | '
          f'cohort separability AUC {_r["AUC_cohort_separability"]:.3f}')
    print('  NOTE: this is single-direction external validation, NOT LODO -- only Yachida')
    print('  has an adenoma+carcinoma gradient, so holding it out leaves nothing to train on.')
    print('  Neither validation cohort has a staging column, so TNM stage is untestable')
    print('  outside Yachida; 3-level severity is the finest axis Kim can confirm.')
else:
    print('  not computed')

print('\nOutputs:')
for f in ['nb05_stage2_early_vs_advanced_v2.csv',
          'nb05_external_transfer_v2.csv',
          'nb05_severity_transfer_v2.csv',
          'two_stage_model_results_v2.csv',
          'two_stage_shap_<polyamine>_v2.csv']:
    print(f'  tables/{f}')
print('  figures/ml/nb05_external_transfer_roc_v2.png')
print('  figures/ml/nb05_severity_transfer_v2.png')
print('  figures/ml/two_stage_shap_<polyamine>_v2.png')

print('\nThis notebook is a leaf: no later v2 notebook reads its outputs. The external')
print('transfer result is nevertheless the headline number for generalisability and')
print('belongs in the manuscript regardless of which direction it points.')'''),
]

nb = Notebook.create(HERE / "05_two_stage_model.ipynb", CELLS)
nb.save()
