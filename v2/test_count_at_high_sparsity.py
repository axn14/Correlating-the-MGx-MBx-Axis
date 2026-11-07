"""Does feature count matter at HIGH sparsity? The last cell of the 2x2.

Run from v2/ with the project venv:

    python test_count_at_high_sparsity.py

Writes Results/extensions/early_axis_count_highsparse_{sweep,fits}_v2.csv.

## The cell

Two things were shown to set, or not set, the scaling exponent p in MDE = k * D^p:

    feature count, at ~4% sparsity   p flat at 0.96 across 462 -> 231 -> 108
    sparsity, at fixed count         p 0.914 -> 0.534 as the floor fraction rises

That leaves one cell: does count matter at HIGH sparsity? It is not idle. A saturation
mechanism -- few tests, loose threshold, many features already near detectable, so extra
samples buy less -- would plausibly operate through the interaction rather than through
count alone, and the low-sparsity row cannot exclude it.

## Why this design rather than feature-band matching

A concurrent analysis filled this cell by drawing features from a narrow floor-fraction
band inside Yachida species, and reported p flat at 0.743-0.755 across a 4-fold count
range. That is the natural approach, but it couples the thing being varied to the thing
being held fixed: sparsity is controlled BY selecting features, so each count arm is a
different selection, and that analysis had to widen its band to fill the largest arm
because only 91 of 500 features fell inside it.

Here sparsity is set by CENSORING and count by RANDOM subsampling, which are independent.
Every arm is the same censored table, so the floor fraction is identical by construction
rather than matched after the fact, and the features in each arm are an unbiased sample of
it. If count mattered at high sparsity this design would show it without a selection
confound to argue about.

Censoring destroys real signal and is **a mechanistic probe, not an analysis**: no
differential abundance result from it is used anywhere.

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
DELTAS = [0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00]
DESIGNS = [(15, 15), (25, 25), (40, 40), (102, 30), (80, 40), (60, 60), (102, 102)]
CENSOR_Q = 0.50                     # targets a floor fraction near 0.50
FEATURE_COUNTS = [108, 231, 462]
N_REPS = 2
SEED = 20261108
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
    meta = pd.concat([D["harmonized_meta"][COHORT], D["adenoma_meta"][COHORT]],
                     ignore_index=True)
    covars = [c for c in utils.COHORT_CONFOUNDERS.get(COHORT, [])
              if c in meta.columns and meta[c].notna().mean() > 0.5
              and meta[c].nunique(dropna=True) > 1]
    tab = pd.concat([D["transformed_mtb"][COHORT], D["adenoma_mtb"][COHORT]],
                    ignore_index=True)
    feats_all = [c for c in tab.columns if c != "Sample"]
    arms = {g: meta.loc[meta[gcol] == g, "Sample"].to_numpy()
            for g in (GROUP1, GROUP2)}

    # Censor ONCE, so every count arm is a subsample of the identical sparse table.
    M = tab[feats_all].to_numpy(float)
    thr = np.quantile(M, CENSOR_Q, axis=0, keepdims=True)
    floor = np.broadcast_to(M.min(axis=0, keepdims=True), M.shape)
    cens = tab.copy()
    cens[feats_all] = np.where(M <= thr, floor, M)
    Mc = cens[feats_all].to_numpy(float)
    ff_full = float(np.mean(Mc == Mc.min(axis=0, keepdims=True)))
    print(f"{COHORT} metabolites censored at q={CENSOR_Q:.2f}: "
          f"floor fraction {ff_full:.3f} over all {len(feats_all)} features; "
          f"covariates {covars}")

    rng = np.random.default_rng(SEED)
    rows, n_proof = [], 0
    for k in FEATURE_COUNTS:
        feats = (feats_all if k >= len(feats_all)
                 else list(rng.choice(feats_all, k, replace=False)))
        sub = cens[["Sample"] + list(feats)]
        Ms = sub[list(feats)].to_numpy(float)
        ff = float(np.mean(Ms == Ms.min(axis=0, keepdims=True)))
        print(f"\n--- {len(feats)} features, floor fraction {ff:.3f} ---")
        for n1, n2 in DESIGNS:
            full = (n1 == len(arms[GROUP1]) and n2 == len(arms[GROUP2]))
            reps = 1 if full else N_REPS
            for rep in range(reps):
                t0 = time.perf_counter()
                keep = np.concatenate([rng.choice(arms[GROUP1], n1, replace=False),
                                       rng.choice(arms[GROUP2], n2, replace=False)])
                meta_s = meta[meta["Sample"].isin(keep)]
                tab_s = sub[sub["Sample"].isin(keep)]
                try:
                    S = utils.RobustDASensitivity(tab_s, meta_s, GROUP1, GROUP2, gcol,
                                                  covars, bias_correct=False)
                except Exception as exc:                        # noqa: BLE001
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
                rows.append({"N_features": len(feats), "Floor_fraction": ff,
                             "N1": n1, "N2": n2, "Rep": rep,
                             "Design_factor": design_factor(n1, n2),
                             "MDE": mde(DELTAS, powers)})
                print(f"    ({n1:3d},{n2:3d}) D={design_factor(n1, n2):.4f} "
                      f"rep {rep + 1}/{reps}  MDE {rows[-1]['MDE']:.3f}  "
                      f"({time.perf_counter() - t0:.0f}s)")

    out = pd.DataFrame(rows)
    out.to_csv(OUT / "early_axis_count_highsparse_sweep_v2.csv", index=False)
    print(f"\nSaved: {OUT / 'early_axis_count_highsparse_sweep_v2.csv'} "
          f"({len(out)} rows); {n_proof} perturbed cases proven")

    print("\n=== exponent against count, at matched high sparsity ===")
    fits = []
    for k, g in out.dropna(subset=["MDE"]).groupby("N_features"):
        per = g.groupby("Design_factor")["MDE"].mean().reset_index()
        popt, _ = curve_fit(m1, per["Design_factor"].to_numpy(float),
                            per["MDE"].to_numpy(float), p0=(2.5, 1.0), maxfev=40000)
        pred = m1(per["Design_factor"].to_numpy(float), *popt)
        yy = per["MDE"].to_numpy(float)
        r2 = 1 - ((yy - pred) ** 2).sum() / ((yy - yy.mean()) ** 2).sum()
        mde102 = float(g[(g.N1 == 102) & (g.N2 == 102)]["MDE"].mean())
        print(f"  {k:4d} features  floor {g['Floor_fraction'].iloc[0]:.3f}  "
              f"p = {popt[1]:.3f}  R2 = {r2:.4f}   MDE at 102/arm = {mde102:.3f}")
        fits.append({"N_features": int(k),
                     "Floor_fraction": float(g["Floor_fraction"].iloc[0]),
                     "p_in_D": float(popt[1]), "R2": r2, "MDE_at_102": mde102})

    fdf = pd.DataFrame(fits).sort_values("N_features")
    fdf.to_csv(OUT / "early_axis_count_highsparse_fits_v2.csv", index=False)
    print(f"\nSaved: {OUT / 'early_axis_count_highsparse_fits_v2.csv'}")

    spread = float(fdf.p_in_D.max() - fdf.p_in_D.min())
    ff_spread = float(fdf.Floor_fraction.max() - fdf.Floor_fraction.min())
    print(f"\n  p spread across a {max(FEATURE_COUNTS) / min(FEATURE_COUNTS):.1f}-fold "
          f"count range: {spread:.3f}")
    print(f"  floor fraction spread across arms: {ff_spread:.4f} "
          "(identical by construction, not matched after the fact)")
    if spread < 0.05:
        print("\n  -> COUNT DOES NOT MATTER AT HIGH SPARSITY EITHER. With the")
        print("     low-sparsity row, the 2x2 is complete: sparsity sets the exponent,")
        print("     count never does, and there is no interaction between them.")
    else:
        print("\n  -> count DOES matter at high sparsity: an interaction exists and the")
        print("     low-sparsity exclusion cannot be generalised.")


if __name__ == "__main__":
    main()
