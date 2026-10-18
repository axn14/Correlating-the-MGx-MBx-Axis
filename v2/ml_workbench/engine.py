"""Evaluation engine for the v2 ML workbench.

No Streamlit imports: loky worker processes import this module, so every task
function below must stay importable and picklable. Everything learned from
data (feature selection, scaling, class-balance weights, tuning) happens
inside a training fold only.
"""
from __future__ import annotations

import functools
import hashlib
import json
import pickle
import sys
import time
import warnings
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Memory
from scipy.stats import mannwhitneyu
from sklearn.ensemble import RandomForestClassifier
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import (RepeatedStratifiedKFold, StratifiedKFold,
                                     StratifiedShuffleSplit)
from sklearn.preprocessing import StandardScaler
from sklearn.utils.class_weight import compute_sample_weight

V2_DIR = Path(__file__).resolve().parent.parent
if str(V2_DIR) not in sys.path:
    sys.path.insert(0, str(V2_DIR))
import utils  # noqa: E402
from utils import DATA_DIR, RESULTS_DIR  # noqa: E402

warnings.filterwarnings('ignore', category=ConvergenceWarning)
# LightGBM 4.6 + scikit-learn 1.7 warn on every numpy predict even though fit also used numpy.
warnings.filterwarnings('ignore', message='X does not have valid feature names')

# Bump when pipeline logic changes, so the on-disk fold cache is not reused.
# v4 (2026-10-08): candidate_targets() now screens on the confounder-adjusted Robust
#   flag from the tidy DA table instead of a nominal PValue<0.05 on the unadjusted
#   fine-grained family, so the target set itself changed and every cached fold
#   computed under the old screen must be discarded.
#   NOTE this was 2 while CONTRACTS section 9 documented 3 -- so the cache was not
#   being invalidated as documented and could serve folds from a superseded group
#   policy. Keep this constant and CONTRACTS in step.
# v2: HS excluded, Adenoma group, no tumour-location features.
ENGINE_VERSION = 4

CRC_DS = 'YACHIDA-CRC-2019'
INTERMEDIATE = RESULTS_DIR / 'intermediate'
CACHE_DIR = INTERMEDIATE / 'ml_workbench_cache_v2'
RUN_LOG = INTERMEDIATE / 'ml_workbench_runs_v2.jsonl'
OPTUNA_JOURNAL = INTERMEDIATE / 'ml_workbench_optuna_v2.log'

# No tumour-location group: "NoTumor" exactly identifies the non-cancer groups, so it
# leaks disease status and was removed from the feature set (utils.build_confounder_matrix).
GROUPS = ['species', 'polyamine_species', 'clinical']
GROUP_LABELS = {
    'species': 'Species (other)',
    'polyamine_species': 'Polyamine-pathway species (AGORA2)',
    'clinical': 'Clinical confounders',
}
CLINICAL_COLS = ['Age', 'BMI', 'Brinkman.Index', 'Alcohol', 'Gender_enc']

MODELS = ['xgboost', 'lightgbm', 'elasticnet', 'random_forest']
MODEL_LABELS = {'xgboost': 'XGBoost', 'lightgbm': 'LightGBM',
                'elasticnet': 'Elastic-net logistic', 'random_forest': 'Random forest'}
DEFAULT_PARAMS = {
    # max_bin=64: measured 4.4 vs 12.8 CPU-s per fit at the default 256 on this data.
    'xgboost': dict(n_estimators=300, max_depth=4, learning_rate=0.05, subsample=0.8,
                    colsample_bytree=0.8, min_child_weight=1.0, reg_lambda=1.0, max_bin=64),
    'lightgbm': dict(n_estimators=300, num_leaves=15, learning_rate=0.05, min_child_samples=10,
                     subsample=0.8, subsample_freq=1, colsample_bytree=0.8, reg_lambda=1.0,
                     max_bin=63),
    'elasticnet': dict(C=0.1, l1_ratio=0.5),
    'random_forest': dict(n_estimators=300, max_features=0.2, min_samples_leaf=2, max_depth=12),
}

FS_METHODS = ['none', 'univariate', 'l1', 'shap', 'stability']
FS_LABELS = {'none': 'None (all features)', 'univariate': 'Univariate (Mann-Whitney) top-k',
             'l1': 'L1-logistic', 'shap': 'SHAP-ranked top-k', 'stability': 'Stability selection'}

memory = Memory(CACHE_DIR, verbose=0)


# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------

def polyamine_pathway_species(species_cols) -> set:
    """Species whose AGORA2 reconstruction encodes a biosynthetic enzyme for a
    polyamine proper (utils.POLYAMINE_PATHWAY_KEGG; ornithine/arginine synthesis
    ECs are near-universal and would put almost every species in this group)."""
    idx = utils.load_agora2_ec_index()
    return {s for s in species_cols
            if any(utils.species_encodes_enzyme(s, k, idx) for k in utils.POLYAMINE_PATHWAY_KEGG)}


def candidate_targets(mtb: pd.DataFrame, polyamine_cols: dict) -> pd.DataFrame:
    """Metabolites that are ROBUST DA hits in either Healthy-vs-CRC comparison, plus
    every polyamine column (always kept, as pre-specified hypothesis targets).

    v2.1 -- two fixes, both of which changed which targets this returns:

    1. It read the FINE-GRAINED DA family (Healthy vs Stage_I_II / Stage_III_IV),
       which NB02 wrote with an unadjusted Mann-Whitney only. NB03 reads the GROUPED,
       confounder-adjusted family (Healthy vs Early_CRC / Advanced_CRC) instead, so
       the workbench and NB03 were screening the same cohort with different
       statistics and disagreeing about which was authoritative. Both now read the
       one tidy table NB02 writes.
    2. It screened on a nominal `PValue < 0.05`. CONTRACTS section 9 states candidate
       targets come from the robust-hit tables, not a nominal-p screen -- an
       unadjusted p-value here lets smoking-associated features in as "findings"
       (Brinkman index differs 2-4x across these groups in this cohort). The screen
       is now the `Robust` flag: significant under BOTH tests, in the same direction.
    """
    path = RESULTS_DIR / 'tables' / 'differential_abundance_v2.csv'
    if not path.exists():
        raise FileNotFoundError(
            f'{path.name} not found. NB02 Step 6B writes it; run 02_analysis.ipynb. '
            '(Before v2.1 this read the per-comparison *_DA_v2.csv files, which are '
            'no longer written.)'
        )
    da = pd.read_csv(path)
    missing = {'Dataset', 'Comparison', 'Assay', 'Feature', 'Robust'} - set(da.columns)
    if missing:
        raise KeyError(
            f'{path.name} is missing {sorted(missing)} -- it predates the v2.1 tidy DA '
            'table. Re-run 02_analysis.ipynb.'
        )

    tabs = []
    for label, group in (('early', 'Early_CRC'), ('advanced', 'Advanced_CRC')):
        sub = da[(da['Dataset'] == CRC_DS)
                 & (da['Comparison'] == f'Healthy_vs_{group}')
                 & (da['Assay'] == 'metabolites')]
        if sub.empty:
            raise ValueError(
                f'{path.name} has no Healthy_vs_{group} metabolite rows for {CRC_DS}. '
                f'Comparisons present: {sorted(da["Comparison"].unique())}'
            )
        cols = ['CliffsD', 'PValue', 'QValue', 'QValue_adj', 'Robust']
        tabs.append(sub.set_index('Feature')[cols].add_prefix(f'{label}_'))

    t = tabs[0].join(tabs[1], how='outer')
    pa = set(polyamine_cols.values())
    # Robust in either comparison, or a pre-specified polyamine.
    robust_either = (t['early_Robust'].fillna(False).astype(bool)
                     | t['advanced_Robust'].fillna(False).astype(bool))
    keep = robust_either | t.index.isin(pa)
    t = t[keep & t.index.isin(mtb.columns)].copy()
    t['is_polyamine'] = t.index.isin(pa)
    t['is_robust'] = robust_either.reindex(t.index).fillna(False).astype(bool)
    # Rank on the ADJUSTED q-value, consistent with the screen above.
    t['min_q'] = t[['early_QValue_adj', 'advanced_QValue_adj']].min(axis=1)
    t['name'] = [f.split('_', 1)[1] if '_' in f else f for f in t.index]
    return t.sort_values('min_q').rename_axis('Feature').reset_index()


@functools.lru_cache(maxsize=2)
def load_workbench_data(holdout_frac: float = 0.2, split_seed: int = 42) -> dict:
    """Load once per process. The holdout split is stratified on Stage.Group,
    fixed across targets, and fully determined by (holdout_frac, split_seed)."""
    with open(INTERMEDIATE / 'preprocessed_data_reduced_v2.pkl', 'rb') as f:
        prep = pickle.load(f)
    spc = prep['transformed_species'][CRC_DS].set_index('Sample')
    mtb = prep['transformed_mtb'][CRC_DS].set_index('Sample')
    meta = prep['harmonized_meta'][CRC_DS].set_index('Sample')
    # Confounders from the pickle's harmonized metadata (same samples; HS excluded).
    conf = utils.build_confounder_matrix(prep['harmonized_meta'][CRC_DS]).set_index('Sample')

    samples = pd.Index(sorted(set(spc.index) & set(mtb.index) & set(conf.index) & set(meta.index)))
    X = pd.concat([spc.loc[samples], conf.loc[samples, CLINICAL_COLS]], axis=1).astype(float)

    pa_species = polyamine_pathway_species(list(spc.columns))
    groups = {c: ('polyamine_species' if c in pa_species else 'species') for c in spc.columns}
    groups.update({c: 'clinical' for c in CLINICAL_COLS})

    stage = meta.loc[samples, 'Stage.Group'].astype(str)
    sss = StratifiedShuffleSplit(n_splits=1, test_size=holdout_frac, random_state=split_seed)
    dev_pos, hold_pos = next(sss.split(np.zeros(len(samples)), stage))

    return dict(
        X=X, groups=groups, mtb=mtb.loc[samples], stage=stage.to_numpy(),
        samples=samples, dev_pos=np.sort(dev_pos), hold_pos=np.sort(hold_pos),
        targets=candidate_targets(mtb, prep['polyamine_columns'][CRC_DS]),
        holdout_frac=holdout_frac, split_seed=split_seed,
    )


def target_labels(data: dict, target: str):
    """High/low label split at the DEVELOPMENT-set median, applied unchanged
    to the holdout."""
    y = data['mtb'][target].to_numpy(float)
    thr = float(np.median(y[data['dev_pos']]))
    return (y > thr).astype(int), thr


def active_features(data: dict, group_weights: dict):
    """Drop weight-0 groups; return (feature names, per-feature weights)."""
    feats = [f for f in data['X'].columns if group_weights.get(data['groups'][f], 1.0) > 0]
    fw = np.array([group_weights.get(data['groups'][f], 1.0) for f in feats], dtype=float)
    return feats, fw


# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

def default_cfg(target: str) -> dict:
    return dict(
        target=target, model='xgboost', params={},
        group_weights={'species': 1.0, 'polyamine_species': 1.0, 'clinical': 1.0},
        class_balance=False,
        stage_weights={g: 1.0 for g in utils.STAGE_GROUP_ORDER},
        fs_method='none', fs_k=50,
        fs_params={'l1_C': 0.1, 'stability_B': 30, 'stability_pi': 0.6},
        n_splits=5, n_repeats=1, cv_seed=42,
        holdout_frac=0.2, split_seed=42,
    )


def canon(cfg: dict) -> str:
    return json.dumps({**cfg, 'engine_version': ENGINE_VERSION}, sort_keys=True, default=float)


def cfg_hash(cfg: dict) -> str:
    return hashlib.sha1(canon(cfg).encode()).hexdigest()[:12]


def n_folds(cfg: dict) -> int:
    return cfg['n_splits'] * cfg['n_repeats']


# ---------------------------------------------------------------------------
# Pipeline pieces (always called on a training fold)
# ---------------------------------------------------------------------------

def sample_weights(y_tr, stage_tr, cfg) -> np.ndarray:
    w = np.ones(len(y_tr))
    if cfg['class_balance']:
        w *= compute_sample_weight('balanced', y_tr)
    sw = cfg['stage_weights']
    w *= np.array([sw.get(s, 1.0) for s in stage_tr])
    return w


class _Fitted:
    def __init__(self, name, est, scaler=None, col_scale=None):
        self.name, self.est, self.scaler, self.col_scale = name, est, scaler, col_scale

    def _prep(self, X):
        return self.scaler.transform(X) * self.col_scale if self.scaler is not None else X

    def predict_proba(self, X) -> np.ndarray:
        return self.est.predict_proba(self._prep(X))[:, 1]

    def importance(self, X) -> np.ndarray:
        """Mean |contribution| per feature (SHAP for trees, |coef| for elastic-net)."""
        if self.name == 'xgboost':
            from xgboost import DMatrix
            contrib = self.est.get_booster().predict(DMatrix(X), pred_contribs=True)[:, :-1]
            return np.abs(contrib).mean(axis=0)
        if self.name == 'lightgbm':
            return np.abs(self.est.predict(X, pred_contrib=True)[:, :-1]).mean(axis=0)
        if self.name == 'random_forest':
            import shap
            sv = np.asarray(shap.TreeExplainer(self.est).shap_values(X, check_additivity=False))
            sv = sv[..., 1] if sv.ndim == 3 else sv
            return np.abs(sv).mean(axis=0)
        return np.abs(self.est.coef_.ravel())


def fit_model(name, params, X_tr, y_tr, w_tr, fw) -> _Fitted:
    """Group weights act model-natively where a mechanism exists: XGBoost
    feature_weights (column-sampling probability), LightGBM feature_contri
    (split-gain multiplier), elastic-net column scaling after standardisation
    (per-feature penalty relief). Random forest has none."""
    p = {**DEFAULT_PARAMS[name], **(params or {})}
    uniform = fw is None or np.allclose(fw, fw[0])
    if name == 'xgboost':
        from xgboost import XGBClassifier
        est = XGBClassifier(n_jobs=1, random_state=0, eval_metric='logloss', tree_method='hist',
                            feature_weights=None if uniform else fw, **p)
        return _Fitted(name, est.fit(X_tr, y_tr, sample_weight=w_tr))
    if name == 'lightgbm':
        from lightgbm import LGBMClassifier
        extra = {} if uniform else {'feature_contri': fw.tolist()}
        est = LGBMClassifier(n_jobs=1, random_state=0, verbose=-1, **p, **extra)
        return _Fitted(name, est.fit(X_tr, y_tr, sample_weight=w_tr))
    if name == 'elasticnet':
        scaler = StandardScaler().fit(X_tr)
        col_scale = np.ones(X_tr.shape[1]) if uniform else fw
        est = LogisticRegression(penalty='elasticnet', solver='saga', max_iter=3000, **p)
        est.fit(scaler.transform(X_tr) * col_scale, y_tr, sample_weight=w_tr)
        return _Fitted(name, est, scaler, col_scale)
    if name == 'random_forest':
        est = RandomForestClassifier(n_jobs=1, random_state=0, **p)
        return _Fitted(name, est.fit(X_tr, y_tr, sample_weight=w_tr))
    raise ValueError(f'unknown model {name!r}')


def _l1_coefs(Z, y, w, C):
    lr = LogisticRegression(penalty='l1', solver='liblinear', C=C)
    return np.abs(lr.fit(Z, y, sample_weight=w).coef_.ravel())


def select_features(cfg, params, X_tr, y_tr, w_tr, fw):
    """Return (selected column indices, note). Scores are multiplied by the
    feature-group weights before ranking."""
    method, p = cfg['fs_method'], X_tr.shape[1]
    k = min(int(cfg['fs_k']), p)
    if method == 'none':
        return np.arange(p), ''

    def univariate():
        pv = mannwhitneyu(X_tr[y_tr == 1], X_tr[y_tr == 0], axis=0).pvalue
        score = -np.log10(np.clip(np.nan_to_num(pv, nan=1.0), 1e-300, 1.0)) * fw
        return np.argsort(-score, kind='stable')[:k]

    if method == 'univariate':
        return univariate(), ''
    if method == 'shap':
        imp = fit_model(cfg['model'], params, X_tr, y_tr, w_tr, fw).importance(X_tr) * fw
        return np.argsort(-imp, kind='stable')[:k], ''

    Z = StandardScaler().fit_transform(X_tr) * fw
    C = cfg['fs_params']['l1_C']
    if method == 'l1':
        coef = _l1_coefs(Z, y_tr, w_tr, C)
        nz = np.flatnonzero(coef)
        if len(nz) == 0:
            return univariate(), 'L1 kept no features at this C; fell back to univariate'
        return nz[np.argsort(-coef[nz], kind='stable')][:k], ''
    if method == 'stability':
        B, pi = int(cfg['fs_params']['stability_B']), float(cfg['fs_params']['stability_pi'])
        rng = np.random.default_rng(0)
        n = len(y_tr)
        counts = np.zeros(p)
        for _ in range(B):
            sub = rng.choice(n, n // 2, replace=False)
            if len(np.unique(y_tr[sub])) == 2:
                counts += _l1_coefs(Z[sub], y_tr[sub], w_tr[sub], C) > 0
        freq = counts / B
        sel = np.flatnonzero(freq >= pi)
        if len(sel) == 0:
            return np.argsort(-freq, kind='stable')[:k], f'no feature reached pi={pi}; kept top-{k} by frequency'
        return sel[np.argsort(-freq[sel], kind='stable')][:k], ''
    raise ValueError(f'unknown feature-selection method {method!r}')


def run_pipeline(cfg, params, X_tr, y_tr, stage_tr, X_te, fw):
    """weights -> in-fold selection -> fit -> predict. Returns (proba, selected idx, note)."""
    w_tr = sample_weights(y_tr, stage_tr, cfg)
    sel, note = select_features(cfg, params, X_tr, y_tr, w_tr, fw)
    fm = fit_model(cfg['model'], params, X_tr[:, sel], y_tr, w_tr, fw[sel])
    return fm.predict_proba(X_te[:, sel]), sel, note


def cv_folds(y, n_splits, n_repeats, seed):
    return list(RepeatedStratifiedKFold(n_splits=n_splits, n_repeats=n_repeats,
                                        random_state=seed).split(np.zeros(len(y)), y))


def _dev_arrays(cfg):
    data = load_workbench_data(cfg['holdout_frac'], cfg['split_seed'])
    feats, fw = active_features(data, cfg['group_weights'])
    y_all, thr = target_labels(data, cfg['target'])
    dev = data['dev_pos']
    X = data['X'][feats].to_numpy(float)
    return data, feats, fw, X[dev], y_all[dev], data['stage'][dev], dev


# ---------------------------------------------------------------------------
# Worker tasks
# ---------------------------------------------------------------------------

@memory.cache
def _run_fold_cached(cfg_json: str, fold_i: int) -> dict:
    cfg = json.loads(cfg_json)
    data, feats, fw, X, y, stage, dev = _dev_arrays(cfg)
    tr, te = cv_folds(y, cfg['n_splits'], cfg['n_repeats'], cfg['cv_seed'])[fold_i]
    t0 = time.perf_counter()
    proba, sel, note = run_pipeline(cfg, cfg['params'], X[tr], y[tr], stage[tr], X[te], fw)

    # Baseline on identical folds: clinical confounders only,
    # same model family at defaults, same sample weights, no selection.
    conf_cols = [c for c in data['X'].columns if data['groups'][c] == 'clinical']
    Xc = data['X'][conf_cols].to_numpy(float)[dev]
    base = fit_model(cfg['model'], None, Xc[tr], y[tr], sample_weights(y[tr], stage[tr], cfg), None)
    base_proba = base.predict_proba(Xc[te])

    return dict(fold=fold_i, auc=float(roc_auc_score(y[te], proba)),
                baseline_auc=float(roc_auc_score(y[te], base_proba)),
                test_idx=te.tolist(), proba=proba.tolist(), baseline_proba=base_proba.tolist(),
                selected=[feats[i] for i in sel], note=note,
                seconds=time.perf_counter() - t0)


def run_fold(cfg_json: str, fold_i: int) -> dict:
    """One CV fold of the development set (memoised on disk by config)."""
    return _run_fold_cached(cfg_json, fold_i)


def assemble_run(cfg: dict, fold_results: list) -> dict:
    data, feats, fw, X, y, stage, dev = _dev_arrays(cfg)
    fold_results = sorted(fold_results, key=lambda r: r['fold'])
    oof_sum, base_sum, oof_cnt = np.zeros(len(y)), np.zeros(len(y)), np.zeros(len(y))
    freq = Counter()
    for r in fold_results:
        oof_sum[r['test_idx']] += r['proba']
        base_sum[r['test_idx']] += r['baseline_proba']
        oof_cnt[r['test_idx']] += 1
        freq.update(r['selected'])
    nf = len(fold_results)
    aucs = np.array([r['auc'] for r in fold_results])
    base = np.array([r['baseline_auc'] for r in fold_results])
    return dict(
        run_id=cfg_hash(cfg), kind='run', created=time.time(), cfg=cfg,
        target=cfg['target'], model=cfg['model'], fs_method=cfg['fs_method'],
        auc_mean=float(aucs.mean()), auc_sd=float(aucs.std()), fold_aucs=aucs.tolist(),
        baseline_auc_mean=float(base.mean()), baseline_fold_aucs=base.tolist(),
        gain=float(aucs.mean() - base.mean()),
        n_features_mean=float(np.mean([len(r['selected']) for r in fold_results])),
        n_features_total=len(feats),
        selection_freq={f: c / nf for f, c in freq.items()},
        oof=(oof_sum / np.maximum(oof_cnt, 1)).tolist(),
        baseline_oof=(base_sum / np.maximum(oof_cnt, 1)).tolist(), y_dev=y.tolist(),
        notes=sorted({r['note'] for r in fold_results if r['note']}),
        seconds=float(sum(r['seconds'] for r in fold_results)),
    )


def permutation_chunk(cfg_json: str, seeds: list) -> list:
    """Null CV-AUCs: the full pipeline (weights, in-fold selection, model with
    fixed params) re-run on shuffled labels. Folds are re-stratified on each
    shuffled label vector, as sklearn's permutation_test_score does, so the
    statistic is the same function the observed run used."""
    cfg = json.loads(cfg_json)
    _, _, fw, X, y, stage, _ = _dev_arrays(cfg)
    out = []
    for s in seeds:
        yp = np.random.default_rng(int(s)).permutation(y)
        aucs = [roc_auc_score(yp[te], run_pipeline(cfg, cfg['params'], X[tr], yp[tr], stage[tr], X[te], fw)[0])
                for tr, te in cv_folds(yp, cfg['n_splits'], cfg['n_repeats'], cfg['cv_seed'])]
        out.append(float(np.mean(aucs)))
    return out


def holdout_eval(cfg_json: str, n_boot: int = 2000) -> dict:
    """Fit on the whole development set, score the locked holdout once."""
    cfg = json.loads(cfg_json)
    data = load_workbench_data(cfg['holdout_frac'], cfg['split_seed'])
    feats, fw = active_features(data, cfg['group_weights'])
    y_all, thr = target_labels(data, cfg['target'])
    dev, hold = data['dev_pos'], data['hold_pos']
    X = data['X'][feats].to_numpy(float)
    proba, sel, note = run_pipeline(cfg, cfg['params'], X[dev], y_all[dev], data['stage'][dev], X[hold], fw)
    yh = y_all[hold]

    conf_cols = [c for c in data['X'].columns if data['groups'][c] == 'clinical']
    Xc = data['X'][conf_cols].to_numpy(float)
    base = fit_model(cfg['model'], None, Xc[dev], y_all[dev],
                     sample_weights(y_all[dev], data['stage'][dev], cfg), None)
    base_proba = base.predict_proba(Xc[hold])

    rng = np.random.default_rng(0)
    boots = []
    for _ in range(n_boot):
        b = rng.integers(0, len(hold), len(hold))
        if len(np.unique(yh[b])) == 2:
            boots.append(roc_auc_score(yh[b], proba[b]))
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return dict(auc=float(roc_auc_score(yh, proba)), ci_low=float(lo), ci_high=float(hi),
                baseline_auc=float(roc_auc_score(yh, base_proba)),
                n_holdout=int(len(hold)), n_pos=int(yh.sum()), threshold=thr,
                n_selected=int(len(sel)), note=note, proba=proba.tolist(), y=yh.tolist())


# ---------------------------------------------------------------------------
# Tuning (Optuna)
# ---------------------------------------------------------------------------

def suggest_params(trial, model: str) -> dict:
    """Search spaces. Keys are the estimator's own parameter names, so
    study.best_params can be used directly as `params`."""
    if model == 'xgboost':
        return dict(n_estimators=trial.suggest_int('n_estimators', 100, 500, step=50),
                    max_depth=trial.suggest_int('max_depth', 2, 6),
                    learning_rate=trial.suggest_float('learning_rate', 0.01, 0.2, log=True),
                    subsample=trial.suggest_float('subsample', 0.6, 1.0),
                    colsample_bytree=trial.suggest_float('colsample_bytree', 0.4, 1.0),
                    min_child_weight=trial.suggest_float('min_child_weight', 1.0, 10.0, log=True),
                    reg_lambda=trial.suggest_float('reg_lambda', 1e-3, 10.0, log=True))
    if model == 'lightgbm':
        return dict(n_estimators=trial.suggest_int('n_estimators', 100, 500, step=50),
                    num_leaves=trial.suggest_int('num_leaves', 4, 31),
                    learning_rate=trial.suggest_float('learning_rate', 0.01, 0.2, log=True),
                    min_child_samples=trial.suggest_int('min_child_samples', 5, 40),
                    subsample=trial.suggest_float('subsample', 0.6, 1.0),
                    colsample_bytree=trial.suggest_float('colsample_bytree', 0.4, 1.0),
                    reg_lambda=trial.suggest_float('reg_lambda', 1e-3, 10.0, log=True))
    if model == 'elasticnet':
        return dict(C=trial.suggest_float('C', 1e-3, 10.0, log=True),
                    l1_ratio=trial.suggest_float('l1_ratio', 0.0, 1.0))
    if model == 'random_forest':
        return dict(n_estimators=trial.suggest_int('n_estimators', 100, 500, step=50),
                    max_features=trial.suggest_float('max_features', 0.05, 0.5),
                    min_samples_leaf=trial.suggest_int('min_samples_leaf', 1, 10),
                    max_depth=trial.suggest_int('max_depth', 3, 20))
    raise ValueError(f'unknown model {model!r}')


def optuna_storage():
    from optuna.storages import JournalStorage
    from optuna.storages.journal import JournalFileBackend, JournalFileOpenLock
    OPTUNA_JOURNAL.parent.mkdir(parents=True, exist_ok=True)
    path = str(OPTUNA_JOURNAL)
    # OpenLock: the default symlink lock needs admin rights on Windows.
    return JournalStorage(JournalFileBackend(path, lock_obj=JournalFileOpenLock(path)))


def create_study(cfg: dict) -> str:
    import optuna
    name = f"{cfg['target']}|{cfg['model']}|{cfg['fs_method']}|{int(time.time())}"
    study = optuna.create_study(study_name=name, storage=optuna_storage(), direction='maximize')
    study.set_user_attr('cfg', cfg)
    return name


def _objective_factory(cfg, X, y, stage, fw, folds):
    import optuna

    def objective(trial):
        params = suggest_params(trial, cfg['model'])
        aucs = []
        for i, (tr, te) in enumerate(folds):
            proba = run_pipeline(cfg, params, X[tr], y[tr], stage[tr], X[te], fw)[0]
            aucs.append(roc_auc_score(y[te], proba))
            trial.report(float(np.mean(aucs)), i)
            if trial.should_prune():
                raise optuna.TrialPruned()
        return float(np.mean(aucs))
    return objective


def _pruner():
    import optuna
    return optuna.pruners.MedianPruner(n_startup_trials=5, n_warmup_steps=1)


def optuna_worker(study_name: str, cfg_json: str, n_trials: int, worker_seed: int) -> int:
    """Several of these share one journal-backed study; the UI reads it live.
    Tuning CV runs on the development set, so its best score is optimistic;
    the holdout or nested_audit_fold gives the honest estimate."""
    import optuna
    optuna.logging.set_verbosity(optuna.logging.WARNING)
    cfg = json.loads(cfg_json)
    _, _, fw, X, y, stage, _ = _dev_arrays(cfg)
    folds = cv_folds(y, cfg['n_splits'], 1, cfg['cv_seed'])
    study = optuna.load_study(study_name=study_name, storage=optuna_storage(),
                              sampler=optuna.samplers.TPESampler(seed=worker_seed), pruner=_pruner())
    study.optimize(_objective_factory(cfg, X, y, stage, fw, folds), n_trials=n_trials)
    return n_trials


def nested_audit_fold(cfg_json: str, outer_i: int, n_trials: int) -> dict:
    """Leak-free estimate for a tuned model: tune on the outer-training part
    (3-fold inner CV), score the untouched outer-test fold."""
    import optuna
    optuna.logging.set_verbosity(optuna.logging.WARNING)
    cfg = json.loads(cfg_json)
    _, _, fw, X, y, stage, _ = _dev_arrays(cfg)
    tr, te = cv_folds(y, cfg['n_splits'], 1, cfg['cv_seed'])[outer_i]
    inner = list(StratifiedKFold(3, shuffle=True, random_state=cfg['cv_seed']).split(X[tr], y[tr]))
    study = optuna.create_study(direction='maximize', pruner=_pruner(),
                                sampler=optuna.samplers.TPESampler(seed=outer_i))
    study.optimize(_objective_factory(cfg, X[tr], y[tr], stage[tr], fw, inner), n_trials=n_trials)
    proba = run_pipeline(cfg, study.best_params, X[tr], y[tr], stage[tr], X[te], fw)[0]
    return dict(fold=outer_i, auc=float(roc_auc_score(y[te], proba)),
                inner_best=float(study.best_value), best_params=study.best_params)


# ---------------------------------------------------------------------------
# Scoring and bookkeeping (main process)
# ---------------------------------------------------------------------------

COMPOSITE_PARTS = ['auc', 'stability', 'parsimony', 'gain']


def composite_scores(runs: pd.DataFrame, weights: dict) -> pd.DataFrame:
    """Min-max normalise each part across the given runs, then weight. A part
    with no spread across runs scores 1 for every run (neutral)."""
    parts = pd.DataFrame({
        'auc': runs['auc_mean'],
        'stability': 1 - runs['auc_sd'],
        'parsimony': 1 - runs['n_features_mean'] / runs['n_features_total'],
        'gain': runs['gain'],
    }, index=runs.index)
    spread = (parts.max() - parts.min()).replace(0, np.nan)
    norm = ((parts - parts.min()) / spread).fillna(1.0)
    w = pd.Series({k: float(weights.get(k, 0.0)) for k in COMPOSITE_PARTS})
    score = (norm * w).sum(axis=1) / w.sum() if w.sum() > 0 else pd.Series(0.0, index=runs.index)
    return norm.add_prefix('norm_').assign(composite=score)


def perm_pvalue(observed: float, null: list):
    """(b+1)/(N+1) with a Clopper-Pearson 95% CI on the Monte Carlo estimate."""
    from scipy.stats import binomtest
    null = np.asarray(null)
    b, n = int((null >= observed).sum()), len(null)
    ci = binomtest(b, n).proportion_ci(method='exact') if n else None
    return (b + 1) / (n + 1), (ci.low, ci.high) if ci else (np.nan, np.nan), b, n


def append_log(record: dict) -> None:
    RUN_LOG.parent.mkdir(parents=True, exist_ok=True)
    with open(RUN_LOG, 'a', encoding='utf-8') as f:
        f.write(json.dumps(record, default=float) + '\n')


def read_log() -> list:
    if not RUN_LOG.exists():
        return []
    with open(RUN_LOG, encoding='utf-8') as f:
        return [json.loads(line) for line in f if line.strip()]
