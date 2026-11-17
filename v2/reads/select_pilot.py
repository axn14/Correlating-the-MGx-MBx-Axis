#!/usr/bin/env python3
"""Choose the 30 pilot samples, deterministically and outcome-blind.

    python select_pilot.py                 # writes pilot_30.tsv next to this file
    python select_pilot.py --n-per-group 6

The selection rule, declared before it is run:

  Six samples from each of the five groups on the progression axis -- Healthy,
  MP (adenoma), Stage_0, Stage_I_II, Stage_III_IV -- taken at evenly spaced
  quantiles of total FASTQ volume within each group.

Two reasons for that rule rather than a random draw.

First, the pilot's job is partly to measure cost, and HUMAnN's runtime scales with
read depth. Sampling the depth range inside each group returns a cost distribution
that extrapolates to 347; sampling at random returns six draws that may all sit near
the median and understate the tail. Evenly spaced quantiles are also reproducible
without a seed, so there is no hidden choice in which samples these are.

Second, the selection must not look at any metabolite or any species abundance. The
pilot is meant to answer whether a decarboxylase gene is present and how it tracks
the disease axis; choosing samples on cadaverine, say, to make the Allisonella
question easier to answer would guarantee the answer. Volume and group label are the
only inputs, and neither is an outcome.

HS (post-surgery) is excluded: it is not a point on the progression axis, and the
chapter treats it separately throughout.
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
META = ROOT / "Data" / "YACHIDA-CRC-2019 metadata.tsv"
MANIFEST = HERE / "manifest_347.tsv"
OUT = HERE / "pilot_30.tsv"

AXIS = ["Healthy", "MP", "Stage_0", "Stage_I_II", "Stage_III_IV"]


def quantile_picks(n_items: int, k: int) -> list[int]:
    """k indices spread evenly across 0..n_items-1, endpoints included."""
    if k >= n_items:
        return list(range(n_items))
    if k == 1:
        return [n_items // 2]
    return [round(i * (n_items - 1) / (k - 1)) for i in range(k)]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-per-group", type=int, default=6)
    args = ap.parse_args()

    if not META.is_file():
        sys.exit(f"metadata not found: {META}")
    if not MANIFEST.is_file():
        sys.exit(f"manifest not found: {MANIFEST} -- build it on Dardel first")

    with open(META, newline="", encoding="utf-8-sig") as fh:
        meta = {r["Sample"].strip(): r for r in csv.DictReader(fh, delimiter="\t")}

    vol: dict[str, int] = defaultdict(int)
    runs: dict[str, str] = {}
    with open(MANIFEST, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            vol[r["sample"]] += int(r["bytes"])
            runs[r["sample"]] = r["run"]

    missing = sorted(set(vol) - set(meta))
    if missing:
        sys.exit(f"manifest has {len(missing)} samples absent from the metadata: {missing[:5]}")

    by_group: dict[str, list[str]] = defaultdict(list)
    for s in vol:
        by_group[meta[s]["Study.Group"].strip()].append(s)

    print(f"{'group':<14} {'available':>9} {'picked':>6}   depth range of the picks")
    chosen: list[dict] = []
    for g in AXIS:
        pool = sorted(by_group.get(g, []), key=lambda s: vol[s])
        if not pool:
            sys.exit(f"no samples in group {g!r}")
        idx = quantile_picks(len(pool), args.n_per_group)
        picks = [pool[i] for i in idx]
        lo, hi = vol[picks[0]] / 1e9, vol[picks[-1]] / 1e9
        print(f"{g:<14} {len(pool):>9} {len(picks):>6}   {lo:.2f}-{hi:.2f} GB")
        for rank, s in zip(idx, picks):
            m = meta[s]
            chosen.append({
                "sample": s,
                "run": runs[s],
                "group": g,
                "volume_bytes": vol[s],
                "volume_gb": f"{vol[s] / 1e9:.2f}",
                "depth_rank_in_group": f"{rank + 1}/{len(pool)}",
                "age": m.get("Age", ""),
                "gender": m.get("Gender", ""),
                "bmi": m.get("BMI", ""),
                "stage": m.get("Stage", ""),
            })

    total = sum(c["volume_bytes"] for c in chosen)
    print(f"\n{len(chosen)} samples, {total / 1e9:.1f} GB to transfer in total")
    print(f"peak storage under stream-and-discard: one sample at a time, "
          f"{max(c['volume_bytes'] for c in chosen) / 1e9:.1f} GB plus its intermediates")

    # Sex and age balance are not selection criteria, but a wildly skewed draw would
    # be worth knowing about before the figure is built on it.
    for key in ("gender",):
        counts: dict[str, int] = defaultdict(int)
        for c in chosen:
            counts[c[key] or "?"] += 1
        print(f"incidental {key} balance: {dict(counts)}")
    ages = [float(c["age"]) for c in chosen if c["age"]]
    if ages:
        print(f"incidental age range: {min(ages):.0f}-{max(ages):.0f}, "
              f"median {sorted(ages)[len(ages) // 2]:.0f}")

    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, delimiter="\t", lineterminator="\n",
                           fieldnames=list(chosen[0]))
        w.writeheader()
        w.writerows(chosen)
    print(f"\nwrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
