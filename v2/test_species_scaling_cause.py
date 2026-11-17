"""Why does the species assay scale so much more shallowly than the metabolome?

Run from v2/ with the project venv:

    python test_species_scaling_cause.py

Writes Results/extensions/early_axis_species_cause_v2.csv.

## What is already excluded

Fitting MDE = k * D^p across designs, the Kim species assay gives p = 0.56 while the Kim
metabolome gives p = 0.93 -- equivalent to n^-0.28 against n^-0.46, so species sensitivity
improves with sample size at barely half the classical rate. Three explanations have been
tested and ruled out:

  * a balanced-design artefact -- excluded by the mixed-design sweep, where D does real
    work and the exponent is unchanged (required_n_designs.py);
  * the cohort -- excluded because both assays come from the SAME cohort and samples;
  * the multiplicity burden -- excluded by subsampling the metabolome's features to 108,
    231 and 462, where p stayed flat at 0.96 (test_multiplicity_scaling.py).

So it is a property of the species assay itself. This isolates the most likely component.

## The candidate

The two assays are not tested identically, by design. Species are compositional, so the
adjusted model applies the LinDA bias correction: each feature's group coefficient has the
MODE of all features' coefficients subtracted, because a centred-log-ratio value is a
ratio to the sample's geometric mean and a genuine shift in a subset of taxa imparts a
shared offset to every coefficient (Methods 2.5). The metabolite matrix is not a
composition and gets no such correction.

That correction is estimated from the data, by kernel density over the coefficients. Its
precision therefore improves with sample size too -- but it is a single shared quantity
estimated once per fit, not a per-feature statistic, so there is no reason for it to
improve at the per-feature rate. If the shared term is what limits the species assay, then
removing it should restore something closer to the classical exponent.

**This is a diagnostic, not an alternative analysis.** Running the species assay without
the bias correction is not a valid compositional test and none of its differential
abundance results are used anywhere; the comparison exists only to attribute the scaling.
Both arms are reported so the attribution can be judged.

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
    tab = pd.concat([D["transformed_species"][COHORT], D["adenoma_species"][COHORT]],
                    ignore_index=True)
    arms = {g: meta_full.loc[meta_full[gcol] == g, "Sample"].to_numpy()
            for g in (GROUP1, GROUP2)}
    print(f"{COHORT} species: {tab.shape[1] - 1} features, "
          f"{len(arms[GROUP1])}/{len(arms[GROUP2])}; covariates {covars}")

    rng = np.random.default_rng(SEED)
    rows, n_proof = [], 0
    for bias in (True, False):
        label = ("LinDA bias correction ON (the pipeline's criterion)" if bias
                 else "bias correction OFF (DIAGNOSTIC ONLY -- not a valid "
                      "compositional test)")
        print(f"\n--- {label} ---")
        for n1, n2 in DESIGNS:
            full = (n1 == len(arms[GROUP1]) and n2 == len(arms[GROUP2]))
            reps = 1 if full else N_REPS
            for rep in range(reps):
                t0 = time.perf_counter()
                keep = np.concatenate([rng.choice(arms[GROUP1], n1, replace=False),
                                       rng.choice(arms[GROUP2], n2, replace=False)])
                meta_s = meta_full[meta_full["Sample"].isin(keep)]
                tab_s = tab[tab["Sample"].isin(keep)]
                try:
                    S = utils.RobustDASensitivity(tab_s, meta_s, GROUP1, GROUP2, gcol,
                                                  covars, bias_correct=bias)
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
                rows.append({"Bias_correction": bias, "N1": n1, "N2": n2, "Rep": rep,
                             "Design_factor": design_factor(n1, n2),
                             "MDE": mde(DELTAS, powers)})
                print(f"    ({n1:3d},{n2:3d}) D={design_factor(n1, n2):.4f} "
                      f"rep {rep + 1}/{reps}  MDE {rows[-1]['MDE']:.3f}  "
                      f"({time.perf_counter() - t0:.0f}s)")

    out = pd.DataFrame(rows)
    out.to_csv(OUT / "early_axis_species_cause_v2.csv", index=False)
    print(f"\nSaved: {OUT / 'early_axis_species_cause_v2.csv'} ({len(out)} rows); "
          f"{n_proof} perturbed cases proven")

    print("\n=== exponent with and without the shared bias-correction term ===")
    fits = []
    for bias, g in out.dropna(subset=["MDE"]).groupby("Bias_correction"):
        per = g.groupby("Design_factor")["MDE"].mean().reset_index()
        popt, _ = curve_fit(m1, per["Design_factor"].to_numpy(float),
                            per["MDE"].to_numpy(float), p0=(2.5, 1.0), maxfev=40000)
        pred = m1(per["Design_factor"].to_numpy(float), *popt)
        yy = per["MDE"].to_numpy(float)
        r2 = 1 - ((yy - pred) ** 2).sum() / ((yy - yy.mean()) ** 2).sum()
        mde102 = float(g[(g.N1 == 102) & (g.N2 == 102)]["MDE"].mean())
        print(f"  bias correction {'ON ' if bias else 'OFF'}  p = {popt[1]:.3f}  "
              f"R2 = {r2:.4f}  (equiv. n^-{popt[1] / 2:.3f})   "
              f"MDE at 102/arm = {mde102:.3f}")
        fits.append({"Bias_correction": bool(bias), "p_in_D": float(popt[1]),
                     "k": float(popt[0]), "R2": r2, "MDE_at_102": mde102})

    fdf = pd.DataFrame(fits)
    fdf.to_csv(OUT / "early_axis_species_cause_fits_v2.csv", index=False)
    print(f"\nSaved: {OUT / 'early_axis_species_cause_fits_v2.csv'}")

    on = fdf[fdf.Bias_correction].iloc[0]
    off = fdf[~fdf.Bias_correction].iloc[0]
    print(f"\n  p: {on.p_in_D:.3f} with the correction, {off.p_in_D:.3f} without "
          f"(change {off.p_in_D - on.p_in_D:+.3f})")
    if off.p_in_D > on.p_in_D + 0.15:
        print("  -> the shared bias-correction term accounts for much of the shallow")
        print("     species scaling: removing it moves the exponent toward the")
        print("     classical rate. The correction is still REQUIRED for a")
        print("     compositional assay -- this names a cost of it, not a fix.")
    else:
        print("  -> the bias correction does NOT explain it; the shallow species scaling")
        print("     lies elsewhere in the assay (sparsity, the CLR reference itself, or")
        print("     the genus-level aggregation).")


if __name__ == "__main__":
    main()
