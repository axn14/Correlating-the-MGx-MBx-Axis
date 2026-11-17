"""Does sequencing depth track the stage axis? (CONTRACTS 5d-ii, item 4)

Yachida's metadata carries no read counts. The ENA manifest (`v2/reads/manifest_347.tsv`)
carries the size of each sample's gzipped FASTQ pair, which is proportional to the number
of reads to within compression variation -- a depth proxy. Depth can only confound a stage
result if it differs by stage, so this screens that first. Re-running NB03 with depth as a
covariate (~4.5 h) is warranted only if it does.

Output: out/sample_depth_proxy.tsv (per sample) and out/depth_screen.json (the test).
Run:    python v2/targets/depth_screen.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
V2 = HERE.parent
sys.path.insert(0, str(V2))
import utils  # noqa: E402

OUT = HERE / "out"


def main() -> None:
    man = pd.read_csv(V2 / "reads" / "manifest_347.tsv", sep="\t", dtype={"sample": str})
    depth = man.groupby("sample")["bytes"].sum().rename("fastq_gz_bytes")
    meta = pd.read_csv(utils.DATA_DIR / "YACHIDA-CRC-2019 metadata.tsv", sep="\t",
                       dtype={"Sample": str})
    meta["Sample"] = meta["Sample"].astype(str)
    df = meta.set_index("Sample").join(depth, how="inner")
    df["stage_axis"] = df["Study.Group"].map(utils.CRC_STAGE_NUMERIC) \
        if hasattr(utils, "CRC_STAGE_NUMERIC") else np.nan
    df["log10_bytes"] = np.log10(df["fastq_gz_bytes"])
    OUT.mkdir(parents=True, exist_ok=True)
    df[["Study.Group", "stage_axis", "fastq_gz_bytes", "log10_bytes"]].to_csv(
        OUT / "sample_depth_proxy.tsv", sep="\t")

    on_axis = df.dropna(subset=["stage_axis"])
    rho, p = stats.spearmanr(on_axis["stage_axis"], on_axis["log10_bytes"])
    groups = [g["log10_bytes"].values for _, g in on_axis.groupby("Study.Group")]
    kw = stats.kruskal(*groups)
    by_group = on_axis.groupby("Study.Group")["fastq_gz_bytes"].median() / 1e9
    res = {"n_samples": int(len(on_axis)), "spearman_rho_stage": float(rho),
           "spearman_p": float(p), "kruskal_H": float(kw.statistic),
           "kruskal_p": float(kw.pvalue),
           "median_gb_by_group": {k: round(float(v), 3) for k, v in by_group.items()},
           "verdict": ("depth tracks stage: add it to NB03's confounders"
                       if p < 0.05 else
                       "depth does not track stage at p<0.05: not a stage confounder here")}
    (OUT / "depth_screen.json").write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
