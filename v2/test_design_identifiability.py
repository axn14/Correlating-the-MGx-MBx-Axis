"""Is the MDE intercept identified on Kim once designs are not all balanced?

Run from v2/ after required_n_designs.py:

    python test_design_identifiability.py

Writes Results/extensions/early_axis_design_models_v2.csv.

Three two-parameter models, all in the design factor D so that balanced and unbalanced
designs are on one scale, plus the one-parameter through-origin model as the weak
benchmark it is:

  M0  c * D            1 parameter, through origin
  M1  k * D^p          2 parameters, free exponent, through origin
  M2  a + b * D        2 parameters, intercept

The honest comparison is M1 against M2: both have two parameters, and a model with an
intercept beating a one-parameter constrained model establishes very little. Judged out of
sample on the largest designs, which is the extrapolation direction the manuscript's
recommendation depends on.

The three-parameter a + c*n^p form is deliberately absent. As p -> 0 the constant and
c*n^p become collinear and the fit runs away along that ridge, which is what produced the
a=-516, c=517, p=-0.000 result an earlier revision of Discussion 4.8 mistook for evidence
that the intercept was unidentifiable. It was evidence about the parametrisation.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

sys.path.insert(0, str(Path(__file__).resolve().parent))
import utils  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
pd.set_option("display.width", 220)

EXT = utils.RESULTS_DIR / "extensions"
TARGETS = [0.35, 0.30, 0.25, 0.20]
FIT_MIN_D = 0.2080          # fit on the smaller designs, predict the larger ones


def m0(D, c):
    return c * D


def m1(D, k, p):
    return k * D ** p


def m2(D, a, b):
    return a + b * D


MODELS = {"M0 c*D": (m0, (2.5,)), "M1 k*D^p": (m1, (2.5, 1.0)),
          "M2 a+b*D": (m2, (0.1, 2.0))}


def n_from_D(D):
    if not np.isfinite(D) or D <= 0:
        return np.nan
    return float((1.0 + np.sqrt(1.0 + D * D)) / (D * D))


def invert(name, popt, t):
    if name.startswith("M0"):
        return n_from_D(t / popt[0])
    if name.startswith("M1"):
        k, p = popt
        return n_from_D((t / k) ** (1.0 / p))
    a, b = popt
    return n_from_D((t - a) / b) if t > a else np.nan


def r2_of(f, popt, DD, yy):
    pred = f(DD, *popt)
    return 1 - ((yy - pred) ** 2).sum() / ((yy - yy.mean()) ** 2).sum()


def main():
    path = EXT / "early_axis_design_sweep_v2.csv"
    if not path.exists():
        raise SystemExit("run required_n_designs.py first")
    sweep = pd.read_csv(path).dropna(subset=["MDE"])

    rows = []
    for assay, g in sweep.groupby("Assay"):
        per = (g.groupby(["N1", "N2", "Design_factor", "Balanced"])["MDE"]
                .agg(["mean", "std", "count"]).reset_index()
                .sort_values("Design_factor"))
        print(f"\n{'=' * 76}\n{assay}: {len(per)} designs "
              f"({int(per.Balanced.sum())} balanced, {int((~per.Balanced).sum())} not)")
        print(per.round(4).to_string(index=False))

        train = per[per["Design_factor"] >= FIT_MIN_D]
        held = per[per["Design_factor"] < FIT_MIN_D]
        DD = train["Design_factor"].to_numpy(float)
        yy = train["mean"].to_numpy(float)
        print(f"\n  fit on {len(train)} designs (D >= {FIT_MIN_D}), "
              f"predict {len(held)}")

        fits = {}
        for name, (f, p0) in MODELS.items():
            popt, _ = curve_fit(f, DD, yy, p0=p0, maxfev=40000)
            fits[name] = (f, popt)
            print(f"    {name:10s} params {np.round(popt, 4)}  "
                  f"in-sample R2 {r2_of(f, popt, DD, yy):.4f}")

        print()
        errs = {k: [] for k in fits}
        for _, h in held.iterrows():
            line = (f"    ({int(h.N1):3d},{int(h.N2):3d}) D={h.Design_factor:.4f} "
                    f"{'bal' if h.Balanced else 'unb'} measured {h['mean']:.4f} | ")
            for name, (f, popt) in fits.items():
                pred = float(f(np.array([h.Design_factor]), *popt)[0])
                e = 100 * (pred - h["mean"]) / h["mean"]
                errs[name].append(e)
                line += f"{name.split()[0]} {e:+.1f}%  "
            print(line)
        print("\n    mean |held-out error|: " + "  ".join(
            f"{k.split()[0]} {np.mean(np.abs(v)):.2f}%" for k, v in errs.items()))

        # Fit on everything for the recommendation itself.
        DD_all = per["Design_factor"].to_numpy(float)
        yy_all = per["mean"].to_numpy(float)
        print("\n    required n per arm (balanced equivalent), fitted on all designs:")
        for t in TARGETS:
            vals = {}
            for name, (f, p0) in MODELS.items():
                popt, _ = curve_fit(f, DD_all, yy_all, p0=p0, maxfev=40000)
                vals[name] = invert(name, popt, t)
            finite = [v for v in vals.values() if np.isfinite(v)]
            fold = (max(finite) / min(finite)) if len(finite) > 1 else np.nan
            print(f"      |d|={t:.2f}  " + "  ".join(
                f"{k.split()[0]} {('%7.0f' % v) if np.isfinite(v) else '   n/a'}"
                for k, v in vals.items()) + f"   fold {fold:.1f}x")
            rows.append({"Assay": assay, "Target": t,
                         **{k.split()[0]: v for k, v in vals.items()},
                         "Fold_spread": fold})

        for name, (f, p0) in MODELS.items():
            popt, _ = curve_fit(f, DD_all, yy_all, p0=p0, maxfev=40000)
            rows.append({"Assay": assay, "Target": np.nan, "Model": name.split()[0],
                         "Intercept": (popt[0] if name.startswith("M2") else np.nan),
                         "R2_all": r2_of(f, popt, DD_all, yy_all),
                         "Mean_holdout_err_pct": float(np.mean(np.abs(errs[name])))})

    out = pd.DataFrame(rows)
    out.to_csv(EXT / "early_axis_design_models_v2.csv", index=False)
    print(f"\nSaved: {EXT / 'early_axis_design_models_v2.csv'}")


if __name__ == "__main__":
    main()
