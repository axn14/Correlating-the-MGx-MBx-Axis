"""Does the multiplicity burden set the scaling exponent, holding everything else fixed?

Run from v2/ with the project venv:

    python test_multiplicity_scaling.py

Writes Results/extensions/early_axis_multiplicity_scaling_v2.csv.

## The question

Fitting MDE = k * D^p across designs gives exponents that line up with the number of
features tested rather than with the cohort:

    Kim species          108 features   p = 0.56
    Yachida metabolites  277 features   p = 0.77
    Kim metabolites      462 features   p = 0.93
    Yachida species      500 features   p = 0.85

Kim's two assays sit at opposite ends, so this is not a cohort effect. But across those
four points feature count is confounded with assay type, with platform, and with sparsity
-- Kim species is genus-level 16S, far denser and less zero-inflated than either
metabolome -- so any of those could drive it.

This separates them. One cohort, one assay, one set of samples: Kim metabolites, with the
feature set randomly subsampled to several sizes. Cohort, assay, platform and the sparsity
profile of the features are held fixed by construction, and only the number of tests
entering the Benjamini-Hochberg correction varies. If p falls toward 0.56 as the feature
count falls toward 108, the multiplicity burden is doing the work.

A proposed mechanism, which this design tests rather than assumes: with few features the
BH threshold is loose, many features are already near detectable, the power curve
saturates and extra samples buy less -- flattening the scaling. With many features the
threshold is strict, the design sits far from saturation, and the classical rank-statistic
rate dominates.

## What it would mean for a targeted panel

Methods 2.15 and Discussion 4.8 note that a targeted panel of a few dozen metabolites
would need materially fewer samples than the untargeted assay. That stays true -- a looser
threshold lowers the minimum detectable effect at any given n. What this analysis would ADD
is that such a panel would also scale more shallowly, so additional samples buy less and
its required-n extrapolation is less constrained, not more. Cheaper to reach a given effect,
harder to size reliably.

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

COHORT = "Kim_adenomas_2020"
GROUP1, GROUP2 = "Healthy", "Adenoma"
ASSAY = "metabolites"
DELTAS = [0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00]
DESIGNS = [(15, 15), (25, 25), (40, 40), (102, 30), (80, 40), (60, 60), (102, 102)]
FEATURE_COUNTS = [108, 231, 462]      # 108 matches Kim species; 462 is the full set
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


def main():
    with open(utils.RESULTS_DIR / "intermediate"
              / "preprocessed_data_reduced_v2.pkl", "rb") as fh:
        D = pickle.load(fh)

    gcol = utils.EARLY_AXIS_GROUP_COL[COHORT]
    meta_full = pd.concat([D["harmonized_meta"][COHORT], D["adenoma_meta"][COHORT]],
                          ignore_index=True)
    covars = [c for c in utils.COHORT_CONFOUNDERS.get(COHORT, [])
              if c in meta_full.columns and meta_full[c].notna().mean() > 0.5
              and meta_full[c].nunique(dropna=True) > 1]
    tab_full = pd.concat([D["transformed_mtb"][COHORT], D["adenoma_mtb"][COHORT]],
                         ignore_index=True)
    all_feats = [c for c in tab_full.columns if c != "Sample"]
    arms = {g: meta_full.loc[meta_full[gcol] == g, "Sample"].to_numpy()
            for g in (GROUP1, GROUP2)}
    print(f"{COHORT} {ASSAY}: {len(all_feats)} features, "
          f"{len(arms[GROUP1])}/{len(arms[GROUP2])} samples; covariates {covars}")

    rng = np.random.default_rng(SEED)
    # One fixed feature subset per count, shared across every design, so that the only
    # thing varying within an arm is the design -- not which features were drawn.
    subsets = {}
    for k in FEATURE_COUNTS:
        subsets[k] = (all_feats if k >= len(all_feats)
                      else list(rng.choice(all_feats, k, replace=False)))

    rows, n_proof = [], 0
    for k in FEATURE_COUNTS:
        feats = subsets[k]
        tab_k = tab_full[["Sample"] + list(feats)]
        # Sparsity of the drawn subset, to show it tracks the full set rather than
        # drifting with the draw -- otherwise count and sparsity would still be confounded.
        floor_frac = float(np.mean([
            (tab_k[f] == tab_k[f].min()).mean() for f in feats]))
        print(f"\n--- {len(feats)} features (mean fraction at the detection floor "
              f"{floor_frac:.3f}) ---")
        for n1, n2 in DESIGNS:
            full = (n1 == len(arms[GROUP1]) and n2 == len(arms[GROUP2]))
            reps = 1 if full else N_REPS
            for rep in range(reps):
                t0 = time.perf_counter()
                keep = np.concatenate([rng.choice(arms[GROUP1], n1, replace=False),
                                       rng.choice(arms[GROUP2], n2, replace=False)])
                meta_s = meta_full[meta_full["Sample"].isin(keep)]
                tab_s = tab_k[tab_k["Sample"].isin(keep)]
                try:
                    S = utils.RobustDASensitivity(tab_s, meta_s, GROUP1, GROUP2, gcol,
                                                  covars, bias_correct=False)
                except Exception as exc:                       # noqa: BLE001
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
                rows.append({"N_features": len(feats), "Floor_fraction": floor_frac,
                             "N1": n1, "N2": n2, "Rep": rep,
                             "Design_factor": design_factor(n1, n2),
                             "MDE": mde(DELTAS, powers)})
                print(f"    ({n1:3d},{n2:3d}) D={design_factor(n1, n2):.4f} "
                      f"rep {rep + 1}/{reps}  MDE {rows[-1]['MDE']:.3f}  "
                      f"({time.perf_counter() - t0:.0f}s)")

    out = pd.DataFrame(rows)
    out.to_csv(OUT / "early_axis_multiplicity_scaling_v2.csv", index=False)
    print(f"\nSaved: {OUT / 'early_axis_multiplicity_scaling_v2.csv'} "
          f"({len(out)} rows); {n_proof} perturbed cases proven")

    print("\n=== exponent against feature count, everything else held fixed ===")
    fits = []
    for k, g in out.dropna(subset=["MDE"]).groupby("N_features"):
        per = g.groupby("Design_factor")["MDE"].mean().reset_index()
        popt, _ = curve_fit(m1, per["Design_factor"].to_numpy(float),
                            per["MDE"].to_numpy(float), p0=(2.5, 1.0), maxfev=40000)
        pred = m1(per["Design_factor"].to_numpy(float), *popt)
        yy = per["MDE"].to_numpy(float)
        r2 = 1 - ((yy - pred) ** 2).sum() / ((yy - yy.mean()) ** 2).sum()
        mde_102 = float(g[(g.N1 == 102) & (g.N2 == 102)]["MDE"].mean())
        print(f"  {k:4d} features  p = {popt[1]:.3f}  k = {popt[0]:.3f}  "
              f"R2 = {r2:.4f}  (equiv. n^-{popt[1] / 2:.3f})   "
              f"MDE at 102/arm = {mde_102:.3f}")
        fits.append({"N_features": int(k), "p_in_D": float(popt[1]),
                     "k": float(popt[0]), "R2": r2, "MDE_at_102": mde_102,
                     "Floor_fraction": float(g["Floor_fraction"].iloc[0])})

    fit_df = pd.DataFrame(fits).sort_values("N_features")
    fit_df.to_csv(OUT / "early_axis_multiplicity_fits_v2.csv", index=False)
    print(f"\nSaved: {OUT / 'early_axis_multiplicity_fits_v2.csv'}")

    lo, hi = fit_df.iloc[0], fit_df.iloc[-1]
    print(f"\n  p rises from {lo.p_in_D:.3f} at {int(lo.N_features)} features "
          f"to {hi.p_in_D:.3f} at {int(hi.N_features)}")
    print(f"  MDE at 102/arm falls from {lo.MDE_at_102:.3f} to {hi.MDE_at_102:.3f} "
          "as features are REMOVED"
          if lo.MDE_at_102 > hi.MDE_at_102 else
          f"  MDE at 102/arm rises from {lo.MDE_at_102:.3f} to {hi.MDE_at_102:.3f} "
          "as features are added")
    print(f"  floor fraction across arms: "
          f"{fit_df.Floor_fraction.min():.3f}-{fit_df.Floor_fraction.max():.3f} "
          "(sparsity held approximately constant by construction)")
    if hi.p_in_D > lo.p_in_D:
        print("\n  -> MULTIPLICITY sets the exponent: within one assay, on one set of")
        print("     samples, with sparsity held fixed, removing features flattens the")
        print("     scaling. Feature count, not cohort and not assay type.")
    else:
        print("\n  -> NOT reproduced within-assay: the four-point pattern across cohorts")
        print("     must be driven by assay type or sparsity rather than by count.")


if __name__ == "__main__":
    main()
