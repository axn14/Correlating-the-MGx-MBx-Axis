"""Calibration of NB03's permutation test when confounders drive the target.

Question: when a metabolite depends on the confounders only, which test wrongly
declares the species predictive? H0 is "species independent of the target given the
confounders"; every rejection below is a false positive.

Data: the real Yachida species (CLR, the 500 NB02 keeps) and the real confounders
(Age, BMI, Brinkman index, alcohol, sex; imputed from the confounders alone), read
from Results/intermediate/preprocessed_data_reduced_v2.pkl -- run NB01-NB02 first.
Synthetic targets:
  lin   y = 0.8 z(Age) + 0.5 z(Brinkman) + 0.4 z(sex) + N(0, 1)
  zinf  the same latent with the lowest 60% set to a floor: zero-inflated, so the
        split is detection status, and detection depends on the confounders
Procedures (fixed-parameter XGBoost as in NB03, stratified 5-fold CV AUC):
  old  X = [confounders, species], label = median split of y (detection for zinf)
       -- NB03 before 2026-10-20
  a    X = species, label = median split of the OLS residual of y on the confounders
       -- NB03 since (utils.residualise_on_confounders)
  c    X = species residualised on the confounders, label as `old`

Design. For a fixed X and class balance the label-permutation null does not depend on
the target (a permutation of a 0/1 vector with k ones is a uniform random k-subset), so
one pooled null per (X, k) serves every seed exactly and each seed costs only its
observed fit. Folds are stratified on each labelling (observed and permuted alike), a
valid variant of NB03's fixed-partition design. p = (#null >= obs + 1) / (B + 1).

Result (2026-10-20, 100 seeds, B = 200; ~14 min on 8 cores):
  old  rejects 100% (lin and zinf), mean AUC 0.81
  a    rejects 6.0% (lin) and 5.0% (zinf), mean AUC 0.49; p slightly conservative
  c    rejects 99%, mean AUC 0.69: a zero-abundance sample's species residual is a
       linear function of the confounders, and the trees recover it

usage: python studies/permutation_null_calibration.py [n_seeds=100] [n_null=200]
"""
import pickle
import sys
import time
from pathlib import Path

import numpy as np
from joblib import Parallel, delayed, parallel_backend
from scipy.stats import kstest
from sklearn.experimental import enable_iterative_imputer  # noqa: F401
from sklearn.impute import IterativeImputer
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from xgboost import XGBClassifier

V2 = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(V2))
import utils  # noqa: E402

N_SEEDS = int(sys.argv[1]) if len(sys.argv) > 1 else 100
N_NULL = int(sys.argv[2]) if len(sys.argv) > 2 else 200
N_JOBS = 8
FIXED = dict(n_estimators=200, max_depth=3, learning_rate=0.05, subsample=0.8,
             colsample_bytree=0.8, min_child_weight=1)          # NB03's FIXED_PARAMS


def cv_auc(X, yb):
    aucs = []
    for tr, te in StratifiedKFold(5, shuffle=True, random_state=42).split(X, yb):
        m = XGBClassifier(random_state=42, n_jobs=1, eval_metric="logloss",
                          tree_method="hist", max_bin=64, **FIXED)
        m.fit(X[tr], yb[tr])
        aucs.append(roc_auc_score(yb[te], m.predict_proba(X[te])[:, 1]))
    return float(np.mean(aucs))


def median_split(v):
    return (v > np.median(v)).astype(int)


def main():
    with open(V2.parent / "Results" / "intermediate" / "preprocessed_data_reduced_v2.pkl",
              "rb") as fh:
        prep = pickle.load(fh)
    ds = utils.DATASET_DISCOVERY
    spc, meta = prep["transformed_species"][ds], prep["harmonized_meta"][ds]
    conf = utils.build_confounder_matrix(meta)
    d = spc.merge(conf, on="Sample")
    zcols = [c for c in conf.columns if c != "Sample"]
    Z = IterativeImputer(max_iter=10, random_state=42).fit_transform(d[zcols].to_numpy(float))
    S = d[[c for c in spc.columns if c != "Sample"]].to_numpy(float)
    n = len(S)
    Z1 = np.column_stack([np.ones(n), Z])
    S_res = S - Z1 @ np.linalg.lstsq(Z1, S, rcond=None)[0]
    zs = (Z - Z.mean(0)) / Z.std(0)
    ix = {c: i for i, c in enumerate(zcols)}
    sex = ix[next(c for c in zcols if c.startswith("Gender"))]
    X_SETS = {"old": np.column_stack([Z, S]), "a": S, "c": S_res}
    k_zinf = int(round(0.4 * n))
    t0 = time.time()

    nulls = {}
    with parallel_backend("threading", n_jobs=N_JOBS):
        for key in sorted({("old", n // 2), ("old", k_zinf), ("a", n // 2),
                           ("c", n // 2), ("c", k_zinf)}):
            rng = np.random.default_rng(7)
            labs = []
            for _ in range(N_NULL):
                y = np.zeros(n, int)
                y[rng.choice(n, key[1], replace=False)] = 1
                labs.append(y)
            nulls[key] = np.asarray(Parallel()(delayed(cv_auc)(X_SETS[key[0]], y)
                                               for y in labs))
            print(f"null {key}: mean {nulls[key].mean():.3f}, 95th pct "
                  f"{np.quantile(nulls[key], .95):.3f} ({(time.time() - t0) / 60:.1f} min)",
                  flush=True)

    jobs, rows = [], []
    for seed in range(N_SEEDS):
        rng = np.random.default_rng(seed)
        latent = (0.8 * zs[:, ix["Age"]] + 0.5 * zs[:, ix["Brinkman.Index"]]
                  + 0.4 * zs[:, sex] + rng.normal(0, 1, n))
        cut = np.sort(latent)[n - k_zinf - 1]
        y_zinf = np.where(latent > cut, latent, cut - 1.0)
        for tname, y, yb in (("lin", latent, median_split(latent)),
                             ("zinf", y_zinf, (y_zinf > y_zinf.min()).astype(int))):
            r = utils.residualise_on_confounders(y, Z)
            for pname, lab in (("old", yb), ("a", median_split(r)), ("c", yb)):
                jobs.append((X_SETS[pname], lab))
                rows.append((tname, pname, int(lab.sum())))
    with parallel_backend("threading", n_jobs=N_JOBS):
        obs = Parallel()(delayed(cv_auc)(X, lab) for X, lab in jobs)

    print(f"\n{N_SEEDS} seeds, pooled nulls of {N_NULL} ({(time.time() - t0) / 60:.1f} min). "
          "A valid test rejects ~0.05.")
    for tname in ("lin", "zinf"):
        for pname in ("old", "a", "c"):
            sel = [(o, k) for o, (t, p, k) in zip(obs, rows) if t == tname and p == pname]
            ps = np.asarray([(np.sum(nulls[(pname, k)] >= o) + 1) / (N_NULL + 1)
                             for o, k in sel])
            print(f"  {tname:4s} {pname:3s} reject@.05 {np.mean(ps < 0.05):.3f}  mean AUC "
                  f"{np.mean([o for o, _ in sel]):.3f}  KS-vs-uniform p "
                  f"{kstest(ps, 'uniform').pvalue:.3f}  p quartiles "
                  f"{np.quantile(ps, [.25, .5, .75]).round(2)}")


if __name__ == "__main__":
    main()
