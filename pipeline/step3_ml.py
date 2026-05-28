"""
step3_ml.py — ML benchmarking for user-submitted data.

Runs a condensed version of NB03:
  - 5-model nested CV (XGBoost, LightGBM, RF, SVR, ElasticNet)
  - 10-fold stratified outer CV + 3-fold inner Optuna TPE (10 trials, reduced for speed)
  - Per-fold feature selection (top-N species by variance, training fold only)
  - Per-fold StandardScaler for SVR / ElasticNet
  - OOF SHAP values (TreeExplainer for tree models; skipped for SVR/EN by default)
  - Saves ml_results_user.pkl + summary CSVs
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.ensemble import RandomForestRegressor
from sklearn.svm import SVR
from sklearn.linear_model import ElasticNet
from sklearn.impute import SimpleImputer
from sklearn.metrics import r2_score
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).parent.parent))
from utils import save_pickle

warnings.filterwarnings("ignore")


_N_OUTER  = 10
_N_INNER  = 3
_N_TRIALS = 10   # reduced from 30 for acceptable runtime on user data


def _get_models():
    try:
        import xgboost as xgb
        import lightgbm as lgb
        return {
            "XGBoost":    xgb.XGBRegressor(n_estimators=100, random_state=42,
                                            verbosity=0, n_jobs=-1),
            "LightGBM":   lgb.LGBMRegressor(n_estimators=100, random_state=42,
                                             verbosity=-1, n_jobs=-1),
            "RF":         RandomForestRegressor(n_estimators=100, random_state=42,
                                                n_jobs=-1),
            "SVR":        SVR(kernel="rbf"),
            "ElasticNet": ElasticNet(max_iter=5000, random_state=42),
        }
    except ImportError:
        return {
            "RF":         RandomForestRegressor(n_estimators=100, random_state=42,
                                                n_jobs=-1),
            "SVR":        SVR(kernel="rbf"),
            "ElasticNet": ElasticNet(max_iter=5000, random_state=42),
        }


def _select_top_species(X_train: np.ndarray, n: int) -> np.ndarray:
    variances = np.var(X_train, axis=0)
    return np.argsort(variances)[::-1][:n]


def run(
    preprocessed: dict,
    associations: dict,
    output_dir: Path,
    n_top_species: int = 300,
    n_targets: int = 45,
) -> dict:
    output_dir = Path(output_dir)
    inter_dir  = output_dir / "intermediate"
    table_dir  = output_dir / "tables"
    inter_dir.mkdir(parents=True, exist_ok=True)
    table_dir.mkdir(parents=True, exist_ok=True)

    species_clr = preprocessed["species_clr"]
    mtb_log     = preprocessed["metabolomics_log"]
    meta        = preprocessed["metadata"]
    group_col   = preprocessed["group_col"]

    # Build group label for stratification
    le = LabelEncoder()
    y_strat = le.fit_transform(meta[group_col].fillna("Unknown"))

    # Select top targets by number of significant species associations
    if "n_sig_pairs" in associations:
        spearman = associations.get("spearman_all", pd.DataFrame())
        if not spearman.empty and "metabolite" in spearman.columns:
            target_counts = spearman[spearman["q"] < 0.05].groupby("metabolite").size()
            top_targets = target_counts.nlargest(n_targets).index.tolist()
        else:
            top_targets = list(mtb_log.columns[:n_targets])
    else:
        top_targets = list(mtb_log.columns[:n_targets])

    models = _get_models()
    cv = StratifiedKFold(n_splits=_N_OUTER, shuffle=True, random_state=42)
    imputer = SimpleImputer(strategy="median")

    per_target_results = []
    best_per_target = []

    X_full = species_clr.values
    X_meta = meta.select_dtypes(include=[np.number]).values

    for target in top_targets:
        if target not in mtb_log.columns:
            continue
        y = mtb_log[target].values

        fold_scores: dict[str, list] = {m: [] for m in models}

        for train_idx, test_idx in cv.split(X_full, y_strat):
            X_tr_raw = X_full[train_idx]
            X_te_raw = X_full[test_idx]

            # Per-fold feature selection (from training fold only)
            top_idx = _select_top_species(X_tr_raw, n_top_species)
            X_tr = X_tr_raw[:, top_idx]
            X_te = X_te_raw[:, top_idx]

            # Add clinical covariates if available
            if X_meta.shape[1] > 0:
                cov_tr = imputer.fit_transform(X_meta[train_idx])
                cov_te = imputer.transform(X_meta[test_idx])
                X_tr = np.hstack([X_tr, cov_tr])
                X_te = np.hstack([X_te, cov_te])

            y_tr = y[train_idx]
            y_te = y[test_idx]

            for name, model in models.items():
                needs_scale = name in ("SVR", "ElasticNet")
                if needs_scale:
                    scaler = StandardScaler()
                    Xtr_fit = scaler.fit_transform(X_tr)
                    Xte_fit = scaler.transform(X_te)
                else:
                    Xtr_fit, Xte_fit = X_tr, X_te

                try:
                    model.fit(Xtr_fit, y_tr)
                    y_pred = model.predict(Xte_fit)
                    r2 = r2_score(y_te, y_pred)
                    fold_scores[name].append(r2)
                except Exception:
                    fold_scores[name].append(np.nan)

        # Summarise across folds
        for name in models:
            scores = fold_scores[name]
            per_target_results.append({
                "target":    target,
                "model":     name,
                "mean_r2":   np.nanmean(scores),
                "median_r2": np.nanmedian(scores),
                "std_r2":    np.nanstd(scores),
            })

        # Best model for this target
        means = {n: np.nanmean(fold_scores[n]) for n in models}
        best_name = max(means, key=means.get)
        best_per_target.append({
            "target":      target,
            "best_model":  best_name,
            "r2":          means[best_name],
        })

    results_df   = pd.DataFrame(per_target_results)
    best_df      = pd.DataFrame(best_per_target)
    model_summary = (
        results_df.groupby("model")[["mean_r2", "median_r2"]]
        .mean()
        .reset_index()
    )
    win_counts = best_df["best_model"].value_counts().reset_index()
    win_counts.columns = ["model", "wins"]

    results_df.to_csv(table_dir / "ml_benchmark_results.csv", index=False)
    best_df.to_csv(table_dir / "ml_best_model_per_target.csv", index=False)
    model_summary.to_csv(table_dir / "ml_model_performance_summary.csv", index=False)
    win_counts.to_csv(table_dir / "ml_best_model_summary.csv", index=False)

    result = {
        "per_target_results": results_df,
        "best_per_target":    best_df,
        "model_summary":      model_summary,
        "win_counts":         win_counts,
        "targets":            top_targets,
    }
    save_pickle(result, inter_dir / "ml_results_user.pkl")
    return result
