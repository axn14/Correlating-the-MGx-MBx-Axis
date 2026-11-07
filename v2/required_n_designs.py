"""Is the MDE intercept identifiable on Kim once the designs are not all balanced?

Run from v2/ with the project venv:

    python required_n_designs.py

Writes Results/extensions/early_axis_design_sweep_v2.csv.

## Why this exists

required_n_curve.py swept BALANCED Kim designs only. Over a balanced family the
Mann-Whitney design factor D = sqrt((n1+n2+1)/(n1*n2)) is in one-to-one correspondence
with n per arm, so a power law in n and a linear-in-D model with an intercept are both
two-parameter curves through the same points and can trade places on noise alone. A peer
session's sweep mixed balanced and unbalanced Yachida designs and found a positive
intercept that beat the power law out of sample; that is a real criticism of the balanced
design, not of the arithmetic.

It also correctly identified a pathology in the three-parameter form a + c*n^p that
Discussion 4.8 relied on: as p -> 0 the constant and c*n^p become collinear, so the fit
runs away along that ridge. The observed a=-516, c=517, p=-0.000 was that ridge, and
reading it as "the intercept is not identifiable" confused a parametrisation failure with
a property of the data. The two-parameter a + b*D has no such failure mode.

This script therefore repeats the sweep on Kim with designs that are deliberately NOT all
balanced, including two structurally different designs at nearly the same D -- (40,40) at
0.2250 and (60,30) at 0.2249 -- so D has to do real work. The models are then compared out
of sample on the four largest designs, which is the extrapolation direction the
manuscript's recommendation depends on.

Nothing here touches the pipeline: it reads NB02's frozen pickle read-only and writes to
Results/extensions/, which the runner does not archive.
"""
import pickle
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import utils  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")

COHORT = "Kim_adenomas_2020"
GROUP1, GROUP2 = "Healthy", "Adenoma"
DELTAS = [0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00]
#          n1,  n2   -- deliberately mixed; (40,40) and (60,30) sit at the same D
DESIGNS = [(15, 15), (25, 25), (102, 20), (40, 40), (60, 30), (102, 30),
           (80, 40), (60, 60), (102, 45), (102, 102)]
N_REPS = 3
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
    arms = {g: meta_full.loc[meta_full[gcol] == g, "Sample"].to_numpy()
            for g in (GROUP1, GROUP2)}
    print(f"{COHORT}: {len(arms[GROUP1])}/{len(arms[GROUP2])}; covariates {covars}")
    for n1, n2 in DESIGNS:
        if n1 > len(arms[GROUP1]) or n2 > len(arms[GROUP2]):
            raise SystemExit(f"design ({n1},{n2}) exceeds the available arms")

    rng = np.random.default_rng(SEED)
    rows, n_proof = [], 0
    for assay, store, aden, bias in (("species", "transformed_species",
                                      "adenoma_species", True),
                                     ("metabolites", "transformed_mtb",
                                      "adenoma_mtb", False)):
        tab_full = pd.concat([D[store][COHORT], D[aden][COHORT]], ignore_index=True)
        for n1, n2 in DESIGNS:
            full = (n1 == len(arms[GROUP1]) and n2 == len(arms[GROUP2]))
            reps = 1 if full else N_REPS
            for rep in range(reps):
                t0 = time.perf_counter()
                keep = np.concatenate([rng.choice(arms[GROUP1], n1, replace=False),
                                       rng.choice(arms[GROUP2], n2, replace=False)])
                meta_s = meta_full[meta_full["Sample"].isin(keep)]
                tab_s = tab_full[tab_full["Sample"].isin(keep)]
                try:
                    S = utils.RobustDASensitivity(tab_s, meta_s, GROUP1, GROUP2, gcol,
                                                  covars, bias_correct=bias)
                except Exception as exc:                       # noqa: BLE001
                    print(f"  {assay} ({n1},{n2}) rep {rep}: skipped -- {exc}")
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
                rows.append({"Assay": assay, "N1": n1, "N2": n2, "Rep": rep,
                             "Design_factor": design_factor(n1, n2),
                             "Balanced": n1 == n2, "N_features": len(S.features),
                             "MDE": mde(DELTAS, powers)})
                print(f"  {assay:11s} ({n1:3d},{n2:3d}) D={design_factor(n1, n2):.4f} "
                      f"rep {rep + 1}/{reps}  MDE {rows[-1]['MDE']:.3f}  "
                      f"({time.perf_counter() - t0:.0f}s)")

    out = pd.DataFrame(rows)
    out.to_csv(OUT / "early_axis_design_sweep_v2.csv", index=False)
    print(f"\nSaved: {OUT / 'early_axis_design_sweep_v2.csv'} ({len(out)} rows); "
          f"{n_proof} perturbed cases proven")

    print("\n=== matched-D check: structurally different designs at the same D ===")
    for assay, g in out.groupby("Assay"):
        a = g[(g.N1 == 40) & (g.N2 == 40)]["MDE"]
        b = g[(g.N1 == 60) & (g.N2 == 30)]["MDE"]
        if len(a) and len(b):
            print(f"  {assay:11s} (40,40) D=0.2250 {a.mean():.3f} +- {a.std():.3f}  |  "
                  f"(60,30) D=0.2249 {b.mean():.3f} +- {b.std():.3f}  |  "
                  f"difference {100 * abs(a.mean() - b.mean()) / a.mean():.1f}%")


if __name__ == "__main__":
    main()
