"""Does sparsity set the scaling exponent? A two-sided test.

Run from v2/ with the project venv:

    python test_sparsity_scaling.py

Writes Results/extensions/early_axis_sparsity_{sweep,fits}_v2.csv.

## Why this is the remaining candidate

Fitting MDE = k * D^p across designs, four explanations for the species assay's shallow
exponent have been tested and excluded: the balanced-design family, the cohort, the
multiplicity burden and the compositional bias correction. A fifth survives, and it orders
the four analysed tables perfectly.

Sparsity here means the fraction of the ANALYSED table sitting at the transform's floor,
which is recoverable exactly because each transform ties its floor in a known way:
multiplicative replacement gives every zero in a sample the same delta, so after CLR
former zeros are tied within a ROW at that row's minimum, while the metabolome places
non-detects at one floor per feature, tied within a COLUMN. Measured that way:

    Kim species          108 features   floor 0.504   p = 0.562
    Yachida metabolites  277 features   floor 0.382   p = 0.770
    Yachida species      500 features   floor 0.271   p = 0.851
    Kim metabolites      462 features   floor 0.042   p = 0.926

Spearman -1.000. Four points order perfectly by chance 1 in 24 for a directional
hypothesis, so that is a lead, not a result — and across those points sparsity is still
confounded with feature count, assay type and taxonomic level. This separates them.

## The two directions

**REMOVE.** Within Yachida species alone, two fixed-size subsets of 108 features chosen by
per-feature prevalence: the least sparse 108 and the most sparse 108. Same cohort, same
assay, same transform, same taxonomic level, same samples, same feature count — only the
floor fraction moves (0.04 against 0.54 on the full table). If p falls as sparsity rises,
sparsity is doing the work, and this also breaks the genus-versus-species-level confound
that the four-point comparison cannot.

**INDUCE.** Kim metabolites, natively the densest table at 0.042, censored at rising
per-feature quantiles to manufacture floor fractions near 0.25 and 0.50 — matching Yachida
species and Kim species respectively. If p falls from 0.93 toward 0.56, that is the forward
direction, and the two together are much stronger than either alone.

Censoring destroys real signal and is **a mechanistic probe, not an analysis**: none of its
differential abundance results is used anywhere. The same applies to the prevalence-split
subsets, which are not a pre-specified feature selection.

## A limitation this inherits

The earlier exclusion of the multiplicity burden was established at ~4% sparsity, because
that experiment held sparsity constant by construction. Whether multiplicity matters at 27%
or 50% is untested, and a saturation mechanism would operate through that interaction rather
than through count alone. Discussion 4.8 records the exclusion as measured in a low-sparsity
regime.

Nothing here touches the pipeline: it reads NB02's frozen pickle read-only and writes to
Results/extensions/, which the runner does not archive.
"""
import pickle
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

sys.path.insert(0, str(Path(__file__).resolve().parent))
import utils  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")

Y, K = "YACHIDA-CRC-2019", "Kim_adenomas_2020"
GROUP1, GROUP2 = "Healthy", "Adenoma"
DELTAS = [0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00]
DESIGNS_Y = [(15, 15), (25, 25), (40, 40), (100, 30), (80, 40), (127, 40)]
DESIGNS_K = [(15, 15), (25, 25), (40, 40), (102, 30), (80, 40), (60, 60), (102, 102)]
SUBSET_SIZE = 108
CENSOR_Q = [0.0, 0.25, 0.50]
N_REPS = 2
SEED = 20261107
OUT = utils.RESULTS_DIR / "extensions"
OUT.mkdir(parents=True, exist_ok=True)


def design_factor(n1, n2):
    return float(np.sqrt((n1 + n2 + 1.0) / (n1 * n2)))


def mde(effects, powers, level=0.80):
    x, y = np.asarray(effects, float), np.asarray(powers, float)
    o = np.argsort(x)
    x, y = x[o], y[o]
    hit = np.flatnonzero(y >= level)
    if hit.size == 0:
        return np.nan
    i = hit[0]
    if i == 0:
        return float(x[0])
    return float(x[i - 1] + (level - y[i - 1]) * (x[i] - x[i - 1]) / (y[i] - y[i - 1]))


def m1(D, k, p):
    return k * D ** p


def floor_fraction(df, by_row):
    M = df.drop(columns=["Sample"]).to_numpy(float)
    ref = M.min(axis=1, keepdims=True) if by_row else M.min(axis=0, keepdims=True)
    return float(np.mean(M == ref))


def sweep(tab, meta, arms, designs, covars, gcol, bias, rng, label, extra):
    rows, n_proof = [], 0
    for n1, n2 in designs:
        full = (n1 == len(arms[GROUP1]) and n2 == len(arms[GROUP2]))
        reps = 1 if full else N_REPS
        for rep in range(reps):
            t0 = time.perf_counter()
            keep = np.concatenate([rng.choice(arms[GROUP1], n1, replace=False),
                                   rng.choice(arms[GROUP2], n2, replace=False)])
            meta_s = meta[meta["Sample"].isin(keep)]
            tab_s = tab[tab["Sample"].isin(keep)]
            try:
                S = utils.RobustDASensitivity(tab_s, meta_s, GROUP1, GROUP2, gcol,
                                              covars, bias_correct=bias)
            except Exception as exc:                            # noqa: BLE001
                print(f"    ({n1},{n2}) rep {rep}: skipped -- {exc}")
                continue
            n_proof += S.verify_against_pipeline(
                tab_s, meta_s, covars,
                shift_cases=[(str(rng.choice(S.features)), -0.45)])
            powers = []
            for d in DELTAS:
                ev = rob = 0
                for f in S.features:
                    for sgn in (-1.0, 1.0):
                        r = S.evaluate_shift(f, sgn * d)
                        if r is None:
                            continue
                        ev += 1
                        rob += r["Robust"]
                powers.append(rob / ev if ev else np.nan)
            rows.append({**extra, "Arm": label, "N1": n1, "N2": n2, "Rep": rep,
                         "Design_factor": design_factor(n1, n2),
                         "MDE": mde(DELTAS, powers)})
            print(f"    ({n1:3d},{n2:3d}) D={design_factor(n1, n2):.4f} "
                  f"rep {rep + 1}/{reps}  MDE {rows[-1]['MDE']:.3f}  "
                  f"({time.perf_counter() - t0:.0f}s)")
    return rows, n_proof


def main():
    with open(utils.RESULTS_DIR / "intermediate"
              / "preprocessed_data_reduced_v2.pkl", "rb") as fh:
        D = pickle.load(fh)

    rng = np.random.default_rng(SEED)
    all_rows, n_proof = [], 0

    # ---------------- REMOVE: prevalence-stratified subsets of Yachida species --------
    gcol_y = utils.EARLY_AXIS_GROUP_COL[Y]
    meta_y = pd.concat([D["harmonized_meta"][Y], D["adenoma_meta"][Y]],
                       ignore_index=True)
    cov_y = [c for c in utils.COHORT_CONFOUNDERS.get(Y, [])
             if c in meta_y.columns and meta_y[c].notna().mean() > 0.5
             and meta_y[c].nunique(dropna=True) > 1]
    tab_y = pd.concat([D["transformed_species"][Y], D["adenoma_species"][Y]],
                      ignore_index=True)
    arms_y = {g: meta_y.loc[meta_y[gcol_y] == g, "Sample"].to_numpy()
              for g in (GROUP1, GROUP2)}
    feats_y = [c for c in tab_y.columns if c != "Sample"]
    My = tab_y[feats_y].to_numpy(float)
    per_feat = (My == My.min(axis=1, keepdims=True)).mean(axis=0)
    order = np.argsort(per_feat)
    subsets = {"dense": [feats_y[i] for i in order[:SUBSET_SIZE]],
               "sparse": [feats_y[i] for i in order[-SUBSET_SIZE:]]}

    print(f"=== REMOVE: Yachida species, two {SUBSET_SIZE}-feature subsets "
          f"({len(arms_y[GROUP1])}/{len(arms_y[GROUP2])} samples, covars {cov_y}) ===")
    for name, feats in subsets.items():
        sub = tab_y[["Sample"] + feats]
        ff = floor_fraction(sub, by_row=True)
        print(f"\n  -- {name} subset: {len(feats)} features, "
              f"floor fraction {ff:.3f} --")
        r, np_ = sweep(sub, meta_y, arms_y, DESIGNS_Y, cov_y, gcol_y, True, rng,
                       f"remove_{name}",
                       {"Test": "REMOVE", "Cohort": "Yachida", "Assay": "species",
                        "N_features": len(feats), "Floor_fraction": ff})
        all_rows += r
        n_proof += np_

    # ---------------- INDUCE: censor Kim metabolites upward ---------------------------
    gcol_k = utils.EARLY_AXIS_GROUP_COL[K]
    meta_k = pd.concat([D["harmonized_meta"][K], D["adenoma_meta"][K]],
                       ignore_index=True)
    cov_k = [c for c in utils.COHORT_CONFOUNDERS.get(K, [])
             if c in meta_k.columns and meta_k[c].notna().mean() > 0.5
             and meta_k[c].nunique(dropna=True) > 1]
    tab_k = pd.concat([D["transformed_mtb"][K], D["adenoma_mtb"][K]], ignore_index=True)
    arms_k = {g: meta_k.loc[meta_k[gcol_k] == g, "Sample"].to_numpy()
              for g in (GROUP1, GROUP2)}
    feats_k = [c for c in tab_k.columns if c != "Sample"]

    print(f"\n=== INDUCE: Kim metabolites censored upward "
          f"({len(arms_k[GROUP1])}/{len(arms_k[GROUP2])} samples, covars {cov_k}) ===")
    for q in CENSOR_Q:
        cens = tab_k.copy()
        if q > 0:
            Mk = cens[feats_k].to_numpy(float)
            thr = np.quantile(Mk, q, axis=0, keepdims=True)
            fl = np.broadcast_to(Mk.min(axis=0, keepdims=True), Mk.shape)
            cens[feats_k] = np.where(Mk <= thr, fl, Mk)
        ff = floor_fraction(cens, by_row=False)
        print(f"\n  -- censor quantile {q:.2f}: floor fraction {ff:.3f} --")
        r, np_ = sweep(cens, meta_k, arms_k, DESIGNS_K, cov_k, gcol_k, False, rng,
                       f"induce_q{q:.2f}",
                       {"Test": "INDUCE", "Cohort": "Kim", "Assay": "metabolites",
                        "N_features": len(feats_k), "Floor_fraction": ff})
        all_rows += r
        n_proof += np_

    out = pd.DataFrame(all_rows)
    out.to_csv(OUT / "early_axis_sparsity_sweep_v2.csv", index=False)
    print(f"\nSaved: {OUT / 'early_axis_sparsity_sweep_v2.csv'} ({len(out)} rows); "
          f"{n_proof} perturbed cases proven")

    print("\n=== exponent against sparsity, within each test ===")
    fits = []
    for (test, arm), g in out.dropna(subset=["MDE"]).groupby(["Test", "Arm"]):
        per = g.groupby("Design_factor")["MDE"].mean().reset_index()
        popt, _ = curve_fit(m1, per["Design_factor"].to_numpy(float),
                            per["MDE"].to_numpy(float), p0=(2.5, 1.0), maxfev=40000)
        pred = m1(per["Design_factor"].to_numpy(float), *popt)
        yy = per["MDE"].to_numpy(float)
        r2 = 1 - ((yy - pred) ** 2).sum() / ((yy - yy.mean()) ** 2).sum()
        ff = float(g["Floor_fraction"].iloc[0])
        print(f"  {test:7s} {arm:14s} floor {ff:.3f}  p = {popt[1]:.3f}  "
              f"R2 = {r2:.4f}  (equiv. n^-{popt[1] / 2:.3f})")
        fits.append({"Test": test, "Arm": arm, "Floor_fraction": ff,
                     "p_in_D": float(popt[1]), "R2": r2,
                     "N_features": int(g["N_features"].iloc[0])})

    fdf = pd.DataFrame(fits)
    fdf.to_csv(OUT / "early_axis_sparsity_fits_v2.csv", index=False)
    print(f"\nSaved: {OUT / 'early_axis_sparsity_fits_v2.csv'}")

    print("\n=== verdict ===")
    for test in ("REMOVE", "INDUCE"):
        t = fdf[fdf.Test == test].sort_values("Floor_fraction")
        if len(t) < 2:
            continue
        lo, hi = t.iloc[0], t.iloc[-1]
        d = hi.p_in_D - lo.p_in_D
        print(f"  {test}: floor {lo.Floor_fraction:.3f} -> {hi.Floor_fraction:.3f}, "
              f"p {lo.p_in_D:.3f} -> {hi.p_in_D:.3f}  (change {d:+.3f})")
        print(f"    -> {'SPARSITY LOWERS THE EXPONENT' if d < -0.15 else 'no material effect'}")


if __name__ == "__main__":
    main()
