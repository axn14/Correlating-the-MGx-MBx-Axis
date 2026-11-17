#!/usr/bin/env python3
"""Compare a pulled Dardel run against the laptop's frozen baseline.

    python compare_to_baseline.py <pulled run dir>            # the summary
    python compare_to_baseline.py <pulled run dir> --full      # every target
    python compare_to_baseline.py <pulled run dir> --baseline <other root>

<pulled run dir> is what pull_results.sh writes: a directory holding Results/tables.

Why this exists rather than writeup/verify_manuscript_numbers.py: that script reads
the baseline's tables by a hard-coded absolute path, so it cannot be pointed at
another tree. Aimed at a pulled run it quietly checks the laptop's own tables and
reports no mismatches, which is worse than reporting nothing.

What to expect: differences. A run on other hardware is a new run. OpenBLAS dispatches
microarchitecture-specific kernels, so the same version reduces in a different order on
Zen 2 and Zen 3+, and nested CV amplifies that into percent-level AUC moves because
choosing hyperparameters is an argmax over near-tied fold scores. CONTRACTS.md section
12 records the measured size. The question this script answers is not "do the numbers
match" -- they will not -- but "does anything the chapter claims depend on which
machine ran it".
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

BASELINE_DEFAULT = Path(r"D:\KI Research")
BENCH = "tables/ml_benchmark_dysregulated_metabolites_v2.csv"
KEY = ["Metabolite", "Classifier"]
FDR = 0.05

# Rounded to three places, these are the AUCs the chapter quotes. A target whose
# rounded AUC moves is one whose sentence in the draft would have to change.
QUOTED = {
    "N-Acetylglutamate", "Glu", "Phe", "Putrescine", "Phenethylamine", "Cadaverine",
    "Tryptamine", "Histamine", "GABA", "Trp", "Tyr", "Lys", "His",
    "Spermine", "N1-Acetylspermine",
}


def load(root: Path, what: str) -> pd.DataFrame:
    p = root / "Results" / BENCH
    if not p.is_file():
        sys.exit(f"no {BENCH} under {what} ({root})")
    return pd.read_csv(p).set_index(KEY).sort_index()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("run", type=Path, help="the pulled run directory")
    ap.add_argument("--baseline", type=Path, default=BASELINE_DEFAULT)
    ap.add_argument("--full", action="store_true", help="print every target, not just movers")
    args = ap.parse_args()

    a = load(args.baseline, "baseline")
    b = load(args.run, "pulled run")
    print(f"baseline  {args.baseline}   {len(a)} targets")
    print(f"this run  {args.run}   {len(b)} targets")

    for lbl, idx in (("only in the baseline", a.index.difference(b.index)),
                     ("only in this run", b.index.difference(a.index))):
        if len(idx):
            print(f"  !! {lbl}: {[i[0] for i in idx]}")
    common = a.index.intersection(b.index)
    a, b = a.loc[common], b.loc[common]

    # 1. The gate. This is the only difference that changes what the chapter can claim.
    pa, pb = a.Perm_qvalue < FDR, b.Perm_qvalue < FDR
    print(f"\n== targets beating the permutation null (q<{FDR})")
    print(f"   baseline {int(pa.sum())}   this run {int(pb.sum())}")
    flip = common[pa.values != pb.values]
    if len(flip):
        print("   !! GATE FLIPPED -- these targets change side, so the figures and the")
        print("      producer loop downstream of NB03 differ from the baseline:")
        for k in flip:
            print(f"        {k[0]:42s} q {a.loc[k, 'Perm_qvalue']:.4f} -> "
                  f"{b.loc[k, 'Perm_qvalue']:.4f}   "
                  f"({'now passes' if pb.loc[k] else 'now fails'})")
    else:
        print("   no target changed side of the gate")

    # 2. The quoted AUCs, at the precision the draft quotes them.
    print("\n== AUCs the chapter quotes, rounded to 3 places as the draft rounds them")
    moved = []
    for k in common:
        nm = str(a.loc[k, "Short_name"])
        if nm not in QUOTED:
            continue
        x, y = round(float(a.loc[k, "AUC_nested"]), 3), round(float(b.loc[k, "AUC_nested"]), 3)
        if x != y:
            moved.append((nm, x, y))
    if moved:
        for nm, x, y in sorted(moved, key=lambda t: -abs(t[2] - t[1])):
            print(f"   {nm:42s} {x:.3f} -> {y:.3f}   ({y - x:+.3f})")
        print(f"   {len(moved)} of {len(QUOTED)} quoted AUCs would have to be requoted"
              f" for this run")
    else:
        print("   none moved at 3 decimal places")

    # 3. The magnitude, so the reader can judge whether it is noise or a real shift.
    print("\n== spread across all targets")
    for col in ("AUC_nested", "AUC_confounders_adjusted", "Perm_qvalue",
                "Perm_null_AUC_mean", "R2_confounders"):
        if col not in a.columns or col not in b.columns:
            continue
        d = (b[col] - a[col]).abs()
        print(f"   {col:26s} max |diff| {d.max():.3g}   median {d.median():.3g}"
              f"   moved {int((d > 0).sum())}/{len(d)}")
    if "N_perm" in a.columns:
        d = (b.N_perm - a.N_perm).abs()
        print(f"   {'N_perm':26s} max |diff| {int(d.max()):d} draws"
              f"   moved {int((d > 0).sum())}/{len(d)}  "
              f"(adaptive stopping checks at joblib batch boundaries,"
              f" which move with the worker count)")

    if args.full:
        print("\n== every target")
        t = pd.DataFrame({
            "auc_base": a.AUC_nested, "auc_run": b.AUC_nested,
            "q_base": a.Perm_qvalue, "q_run": b.Perm_qvalue,
            "nperm_base": a.N_perm, "nperm_run": b.N_perm,
        })
        t["d_auc"] = t.auc_run - t.auc_base
        with pd.option_context("display.width", 200, "display.max_rows", 200):
            print(t.reindex(t.d_auc.abs().sort_values(ascending=False).index).to_string())

    print("\nThe frozen baseline stays the baseline. Differences here are evidence about"
          "\nhow robust each claim is to the machine, not output to substitute for it.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
