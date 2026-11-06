"""How many pre-malignant samples per arm would the robust-DA criterion actually need?

Run from v2/ with the project venv:

    python required_n_curve.py

Writes Results/extensions/early_axis_required_n_v2.csv (the measured curve),
Results/extensions/early_axis_required_n_fit_v2.csv (the fitted scaling and the inverted
sample sizes) and a figure beside them.

## Why this is a standalone script and not a notebook cell

The analysis is frozen at commit f14c7f9 (run 20261031T082418Z_f14c7f9). Adding a cell to
NB02 would change its code and so invalidate NB03's reuse key, costing a five-hour re-run
and a fresh reconciliation of every number in the manuscript. This script imports the same
`utils.RobustDASensitivity` the pipeline used, reads NB02's frozen pickle read-only, and
writes to Results/extensions/ -- which the runner does not archive, unlike
Results/{tables,intermediate,figures}. Nothing here can alter a published result.

## What it measures

CONTRACTS.md section 2c reports a minimum detectable effect of |Cliff's delta| 0.36-0.51
against a largest observed effect of 0.33, and concludes that a prospective cohort needs
"more than 102 adenomas per arm". That is a floor, not a design. This turns it into one.

For a grid of sample sizes, the Kim arms are subsampled to n per group and the sensitivity
analysis of section 2c is re-run in full on the subsample: every feature, both effect
directions, the same two-test agreement criterion, the same multiplicity. The minimum
detectable effect at 80% power is read off each curve by interpolation, a power-law scaling
is fitted across sample sizes, and that fit is inverted to give the n per arm needed to
detect a target effect.

Kim is used because its arms are balanced (102 against 102), so a subsample is a smaller
balanced design and the result reads directly as a prospective recommendation. Yachida
cannot serve: its adenoma arm is 40, which is below most of the grid and unbalanced
against 127 controls.

## What it does not establish

The fit is extrapolated beyond the measured range to reach the smaller target effects, and
an extrapolation is an assumption about the scaling, not a measurement. The exponent is
reported so a reader can judge it against the theoretical -1/2 for a rank statistic.
Multiplicity is held fixed at the cohort's own feature count, which is the right choice
for a study using the same assay and the wrong one for a targeted panel -- a panel of
twenty metabolites would need far fewer samples, and that is a different calculation.
Subsampling also cannot create heterogeneity the cohort does not have, so a prospective
cohort recruited across more sites should expect to need more, not fewer, samples.
"""
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import utils  # noqa: E402
import viz  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
viz.apply_theme()

COHORT = "Kim_adenomas_2020"
GROUP1, GROUP2 = "Healthy", "Adenoma"
DELTAS = [0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.60, 0.70, 0.80, 0.90]
N_GRID = [30, 45, 60, 80, 102]
N_REPS = 3                     # 1 at the full size, where a "subsample" is the whole arm
POWER = 0.80
TARGETS = [0.35, 0.30, 0.25, 0.20]
SEED = 20261106

OUT = utils.RESULTS_DIR / "extensions"
OUT.mkdir(parents=True, exist_ok=True)


def mde(effect_sizes, powers, level=POWER):
    """Smallest effect reaching `level`, linearly interpolated; NaN if never reached."""
    x = np.asarray(effect_sizes, float)
    y = np.asarray(powers, float)
    order = np.argsort(x)
    x, y = x[order], y[order]
    hit = np.flatnonzero(y >= level)
    if hit.size == 0:
        return np.nan
    i = hit[0]
    if i == 0:
        return float(x[0])
    x0, x1, y0, y1 = x[i - 1], x[i], y[i - 1], y[i]
    return float(x0 if y1 == y0 else x0 + (level - y0) * (x1 - x0) / (y1 - y0))


def main():
    with open(utils.RESULTS_DIR / "intermediate"
              / "preprocessed_data_reduced_v2.pkl", "rb") as fh:
        import pickle
        D = pickle.load(fh)

    gcol = utils.EARLY_AXIS_GROUP_COL[COHORT]
    meta_full = pd.concat([D["harmonized_meta"][COHORT], D["adenoma_meta"][COHORT]],
                          ignore_index=True)
    covars = [c for c in utils.COHORT_CONFOUNDERS.get(COHORT, [])
              if c in meta_full.columns and meta_full[c].notna().mean() > 0.5
              and meta_full[c].nunique(dropna=True) > 1]
    arms = {g: meta_full.loc[meta_full[gcol] == g, "Sample"].to_numpy()
            for g in (GROUP1, GROUP2)}
    n_full = min(len(arms[GROUP1]), len(arms[GROUP2]))
    print(f"{COHORT}: {len(arms[GROUP1])} {GROUP1} / {len(arms[GROUP2])} {GROUP2}; "
          f"covariates {covars}")
    if max(N_GRID) > n_full:
        raise SystemExit(f"grid asks for n={max(N_GRID)} but only {n_full} available")

    assays = (("species", "transformed_species", "adenoma_species", True),
              ("metabolites", "transformed_mtb", "adenoma_mtb", False))
    rng = np.random.default_rng(SEED)
    rows, n_proof_total = [], 0

    for assay, store, aden, bias in assays:
        tab_full = pd.concat([D[store][COHORT], D[aden][COHORT]], ignore_index=True)
        for n in N_GRID:
            reps = 1 if n >= n_full else N_REPS
            for rep in range(reps):
                t0 = time.perf_counter()
                keep = np.concatenate([
                    rng.choice(arms[GROUP1], n, replace=False),
                    rng.choice(arms[GROUP2], n, replace=False)])
                meta_s = meta_full[meta_full["Sample"].isin(keep)]
                tab_s = tab_full[tab_full["Sample"].isin(keep)]
                S = utils.RobustDASensitivity(tab_s, meta_s, GROUP1, GROUP2, gcol,
                                              covars, bias_correct=bias)
                # Prove the shortcut on this subsample too: the engine is re-built per
                # draw, so its baseline equivalence must hold per draw.
                probe = [(str(rng.choice(S.features)), -0.45)]
                n_proof_total += S.verify_against_pipeline(tab_s, meta_s, covars,
                                                          shift_cases=probe)
                powers = []
                for d in DELTAS:
                    n_eval = n_rob = 0
                    for f in S.features:
                        for sgn in (-1.0, 1.0):
                            r = S.evaluate_shift(f, sgn * d)
                            if r is None:
                                continue
                            n_eval += 1
                            n_rob += r["Robust"]
                    powers.append(n_rob / n_eval if n_eval else np.nan)
                    rows.append({"Cohort": COHORT, "Assay": assay, "N_per_arm": n,
                                 "Replicate": rep, "Effect_size": d,
                                 "N_features": len(S.features),
                                 "N_evaluations": n_eval,
                                 "Power_robust": powers[-1]})
                print(f"  {assay:11s} n={n:3d} rep {rep + 1}/{reps}  "
                      f"MDE {mde(DELTAS, powers):.3f}  ({time.perf_counter() - t0:.0f}s)")

    curve = pd.DataFrame(rows)
    curve.to_csv(OUT / "early_axis_required_n_v2.csv", index=False)
    print(f"\nSaved: {OUT / 'early_axis_required_n_v2.csv'} ({len(curve)} rows); "
          f"{n_proof_total} perturbed cases proven across subsamples")

    # ---- per-draw MDE, then the scaling ------------------------------------
    per_draw = (curve.groupby(["Assay", "N_per_arm", "Replicate"])
                     .apply(lambda g: mde(g["Effect_size"], g["Power_robust"]),
                            include_groups=False)
                     .rename("MDE").reset_index())
    print("\n=== minimum detectable effect per draw ===")
    print(per_draw.pivot_table(index="N_per_arm", columns="Assay", values="MDE",
                               aggfunc=["mean", "std"]).to_string())

    fits = []
    for assay, g in per_draw.groupby("Assay"):
        g = g.dropna(subset=["MDE"])
        if g["N_per_arm"].nunique() < 3:
            print(f"{assay}: too few sample sizes to fit")
            continue
        # MDE = c * n**p, fitted in logs; p is -1/2 for an ideal rank statistic.
        p, logc = np.polyfit(np.log(g["N_per_arm"]), np.log(g["MDE"]), 1)
        c = float(np.exp(logc))
        pred = c * g["N_per_arm"].to_numpy(float) ** p
        ss_res = float(((g["MDE"] - pred) ** 2).sum())
        ss_tot = float(((g["MDE"] - g["MDE"].mean()) ** 2).sum())
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
        for target in TARGETS:
            need = float((target / c) ** (1.0 / p))
            fits.append({"Cohort": COHORT, "Assay": assay, "Exponent": float(p),
                         "Coefficient": c, "R2": r2,
                         "N_measured_max": int(g["N_per_arm"].max()),
                         "Target_abs_CliffsD": target,
                         "Required_N_per_arm": need,
                         "Extrapolated": bool(need > g["N_per_arm"].max())})
        print(f"\n{assay}: MDE = {c:.3f} * n^({p:.3f})   R2 = {r2:.3f}")
        for f in [x for x in fits if x["Assay"] == assay]:
            flag = "  (EXTRAPOLATED)" if f["Extrapolated"] else ""
            print(f"    |d| = {f['Target_abs_CliffsD']:.2f} needs "
                  f"n = {f['Required_N_per_arm']:.0f} per arm{flag}")

    fit_df = pd.DataFrame(fits)
    fit_df.to_csv(OUT / "early_axis_required_n_fit_v2.csv", index=False)
    print(f"\nSaved: {OUT / 'early_axis_required_n_fit_v2.csv'}")

    # ---- figure ------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(7.5, 5))
    for assay, g in per_draw.groupby("Assay"):
        m = g.groupby("N_per_arm")["MDE"].agg(["mean", "std"]).reset_index()
        ax.errorbar(m["N_per_arm"], m["mean"], yerr=m["std"].fillna(0),
                    marker="o", ms=5, capsize=3, label=f"{assay} (measured)")
        row = fit_df[fit_df.Assay == assay]
        if not row.empty:
            c, p = row["Coefficient"].iloc[0], row["Exponent"].iloc[0]
            xs = np.linspace(min(N_GRID), 600, 200)
            ax.plot(xs, c * xs ** p, ls="--", lw=1,
                    label=f"{assay} fit: n^{p:.2f}")
    ax.axvspan(max(N_GRID), 600, color="grey", alpha=0.12)
    ax.text(max(N_GRID) * 1.05, 0.92, "extrapolated", fontsize=8, color="grey",
            transform=ax.get_xaxis_transform())
    for t in TARGETS:
        ax.axhline(t, ls=":", lw=0.8, color="crimson", alpha=0.6)
        ax.text(560, t, f"|δ|={t:.2f}", fontsize=7, color="crimson",
                va="bottom", ha="right")
    ax.set_xlabel("samples per arm")
    ax.set_ylabel("minimum detectable |Cliff's δ| at 80% power")
    ax.set_title("What a prospective pre-malignant cohort would need\n"
                 f"{COHORT.split('_')[0]}, {GROUP1} vs {GROUP2}, "
                 "robust-DA criterion at the cohort's own multiplicity", fontsize=10)
    ax.legend(fontsize=8)
    fig.tight_layout()
    viz.save_figure(fig, OUT / "early_axis_required_n_v2.png")
    print(f"Saved: {OUT / 'early_axis_required_n_v2.png'}")


if __name__ == "__main__":
    main()
