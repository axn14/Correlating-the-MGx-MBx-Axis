"""Does the MDE/sample-size relation have a positive intercept?

A peer session reports, from Yachida species at five designs, that MDE/D is not constant
(2.39 -> 2.68 as the Mann-Whitney design factor D falls) and that a two-point fit gives
MDE = 2.07*D + 0.12. If a positive intercept is real, a pure power law in n understates
the n required for small target effects, and worse the smaller the target -- which would
make the Discussion 4.8 table optimistic.

Their evidence is a different cohort, assay, covariate set and design family from the
Kim-balanced fit in 4.8, so it is tested here against the Kim curve data directly. Three
models, and an honest extrapolation test: fit on the SMALL sizes only and predict the
measured value at the largest, which is the thing the table actually relies on.

  M1  MDE = c * n^p              (what 4.8 currently uses)
  M2  MDE = a + c * n^p          (free intercept)
  M3  MDE = a + b * D            (their parametrisation; D = sqrt((2n+1)/n^2) when balanced)
"""
import sys

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

sys.path.insert(0, r"D:\KI Research\v2")
import utils  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
pd.set_option("display.width", 200)

EXT = utils.RESULTS_DIR / "extensions"
TARGETS = [0.35, 0.30, 0.25, 0.20]
curve = pd.read_csv(EXT / "early_axis_required_n_v2.csv")


def mde(g, level=0.80):
    g = g.sort_values("Effect_size")
    x, y = g["Effect_size"].to_numpy(float), g["Power_robust"].to_numpy(float)
    hit = np.flatnonzero(y >= level)
    if hit.size == 0:
        return np.nan
    i = hit[0]
    if i == 0:
        return float(x[0])
    return float(x[i - 1] + (level - y[i - 1]) * (x[i] - x[i - 1]) / (y[i] - y[i - 1]))


per_draw = (curve.groupby(["Assay", "N_per_arm", "Replicate"])
                 .apply(mde, include_groups=False).rename("MDE").reset_index())


def design_factor(n):           # balanced arms
    return np.sqrt((2.0 * n + 1.0) / (n * n))


def m1(n, c, p):
    return c * n ** p


def m2(n, a, c, p):
    return a + c * n ** p


def m3(n, a, b):
    return a + b * design_factor(n)


def fit_all(nn, yy):
    out = {}
    p, logc = np.polyfit(np.log(nn), np.log(yy), 1)
    out["M1"] = (m1, (float(np.exp(logc)), float(p)))
    try:
        popt, _ = curve_fit(m2, nn, yy, p0=(0.1, 3.0, -0.5), maxfev=40000)
        out["M2"] = (m2, tuple(float(v) for v in popt))
    except Exception as e:                                    # noqa: BLE001
        print(f"      M2 did not converge: {e}")
    b, a = np.polyfit(design_factor(nn), yy, 1)
    out["M3"] = (m3, (float(a), float(b)))
    return out


def r2(f, params, nn, yy):
    pred = f(nn, *params)
    ss_res = float(((yy - pred) ** 2).sum())
    ss_tot = float(((yy - yy.mean()) ** 2).sum())
    return 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan


def required_n(name, params, target):
    """Invert a model for the n per arm reaching `target`; NaN if unreachable."""
    if name == "M1":
        c, p = params
        return float((target / c) ** (1.0 / p))
    if name == "M2":
        a, c, p = params
        if target <= a:
            return np.nan
        return float(((target - a) / c) ** (1.0 / p))
    a, b = params
    D = (target - a) / b
    if D <= 0:
        return np.nan
    return float((1.0 + np.sqrt(1.0 + D * D)) / (D * D))


print("=== peer's ratio diagnostic, applied to Kim ===")
print("If MDE/D rises as D falls, a positive intercept is indicated.\n")
for assay, g in per_draw.groupby("Assay"):
    m = g.groupby("N_per_arm")["MDE"].mean()
    print(f"  {assay}")
    for n, v in m.items():
        D = design_factor(n)
        print(f"    n={n:3d}/arm  D={D:.4f}  MDE={v:.3f}  MDE/D={v / D:6.2f}")
    ratios = [m[n] / design_factor(n) for n in m.index]
    trend = "RISES as D falls (intercept indicated)" if ratios[-1] > ratios[0] \
        else "falls or flat (no intercept indicated)"
    print(f"    -> MDE/D {trend}: {ratios[0]:.2f} -> {ratios[-1]:.2f}\n")

print("=== model comparison on all sizes ===")
rows = []
for assay, g in per_draw.groupby("Assay"):
    g = g.dropna(subset=["MDE"])
    nn, yy = g["N_per_arm"].to_numpy(float), g["MDE"].to_numpy(float)
    fits = fit_all(nn, yy)
    print(f"\n  {assay}")
    for name, (f, params) in fits.items():
        k = len(params)
        n_obs = len(nn)
        rss = float(((yy - f(nn, *params)) ** 2).sum())
        aic = n_obs * np.log(rss / n_obs) + 2 * k
        if name == "M1":
            desc = f"c={params[0]:.3f}, p={params[1]:+.3f}"
        elif name == "M2":
            desc = f"a={params[0]:+.3f}, c={params[1]:.3f}, p={params[2]:+.3f}"
        else:
            desc = f"a={params[0]:+.3f}, b={params[1]:.3f}"
        print(f"    {name}  {desc:46s} R2={r2(f, params, nn, yy):.4f}  AIC={aic:7.1f}")
        for t in TARGETS:
            rows.append({"Assay": assay, "Model": name, "Target": t,
                         "Required_N": required_n(name, params, t)})

print("\n=== required n per arm, by model ===")
req = pd.DataFrame(rows).pivot_table(index=["Assay", "Target"], columns="Model",
                                     values="Required_N")
print(req.round(0).to_string())

print("\n=== THE TEST THAT MATTERS: fit on n<=80, predict n=102 ===")
print("Extrapolation is what 4.8 relies on, so each model is judged on it.\n")
for assay, g in per_draw.groupby("Assay"):
    g = g.dropna(subset=["MDE"])
    train = g[g["N_per_arm"] <= 80]
    held = g[g["N_per_arm"] == 102]
    if train.empty or held.empty:
        continue
    truth = float(held["MDE"].mean())
    nn, yy = train["N_per_arm"].to_numpy(float), train["MDE"].to_numpy(float)
    print(f"  {assay}: measured MDE at n=102 = {truth:.4f}")
    for name, (f, params) in fit_all(nn, yy).items():
        pred = float(f(np.array([102.0]), *params)[0])
        print(f"    {name} predicts {pred:.4f}   error {pred - truth:+.4f} "
              f"({100 * (pred - truth) / truth:+.1f}%)")
    print()
