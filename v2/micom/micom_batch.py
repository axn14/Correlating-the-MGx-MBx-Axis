"""The 277-sample MICOM batch for class M (CONTRACTS 5d), built to run on Dardel.

Four steps, each a subcommand:

  prepare            (laptop) the per-sample taxonomy table micom builds from -- the same
                     cohort, cutoff and AGORA2 matching as `diet_screen.py` -- written to
                     out/batch/taxonomy.tsv, small and auditable, so Dardel needs no pickle.
  build              (Dardel) micom.workflows.build for every sample, models on scratch.
  solve --index i    (Dardel, one SLURM array task per sample) the certified cooperative
                     tradeoff at the chosen alpha and the forced-secretion test
                     (`forced_secretion.forced_exports`), one directory per sample.
  collect            (either) concatenate, and call class M per (species, target).

Class M rule (CONTRACTS 5d, fixed 2026-11-11): a (species, target) pair is SUPPORTED when
the taxon's export is forced (FVA minimum > tolerance) in >= 20% of the samples that
contain it, with >= 10 such samples, AND the target's forced-export base rate among taxa
that can export it is <= 0.2 (`utils.EVIDENCE_MAX_BASE_RATE`). AVAILABLE when the taxon is
modelled in >= 10 samples and the target passed the export gate. A sample whose QP did not
certify (gap > `highs_qp.CERT_TOL`) still contributes, flagged; `collect` reports how many.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))

import highs_qp  # noqa: E402

highs_qp.install("repair")

BATCH = HERE / "out" / "batch"
REGISTRY = HERE.parent / "targets" / "out" / "target_registry.tsv"
CHOICE = HERE / "out" / "alpha_sweep_choice.json"
MIN_FRACTION, MIN_SAMPLES = 0.20, 10


def _scratch() -> Path:
    return Path(os.environ.get("MICOM_SCRATCH", BATCH))


def prepare() -> None:
    import diet_screen as ds
    sp, abund = ds.cohort_species()
    idx = ds.agora_index()
    sp = sp.merge(idx[["species_plain", "path"]], left_on="species",
                  right_on="species_plain", how="left")
    matched = set(sp.loc[sp["path"].notna(), "species"])
    long = abund.stack().rename("abundance").reset_index()
    long.columns = ["sample_id", "species", "abundance"]
    long = long[(long["abundance"] > ds.ABUND_CUTOFF) & long["species"].isin(matched)].copy()
    long["sample_id"] = long["sample_id"].astype(str)
    long["genus"] = long["species"].str.split().str[0]
    long["id"] = long["species"]
    BATCH.mkdir(parents=True, exist_ok=True)
    long[["sample_id", "id", "genus", "species", "abundance"]].to_csv(
        BATCH / "taxonomy.tsv", sep="\t", index=False)
    print(f"{long['sample_id'].nunique()} samples, {long['species'].nunique()} taxa, "
          f"{len(long)} rows -> {BATCH / 'taxonomy.tsv'}")


def build(model_db: str, threads: int) -> None:
    from micom.workflows import build as mbuild
    import diet_screen as ds
    tax = pd.read_csv(BATCH / "taxonomy.tsv", sep="\t", dtype={"sample_id": str})
    models = _scratch() / "models"
    models.mkdir(parents=True, exist_ok=True)
    man = mbuild(tax, out_folder=str(models), model_db=model_db,
                 cutoff=ds.ABUND_CUTOFF, threads=threads)
    man.to_csv(_scratch() / "manifest.csv", index=False)
    print(f"built {len(man)} community models in {models}")


def _medium() -> pd.Series:
    d = pd.read_csv(HERE / "gut_diet_v2.csv", comment="#")
    return d.set_index(d["reaction"].str.replace(r"_e$", "_m", regex=True))["flux"]


def solve(index: int, alpha: float | None) -> None:
    from micom import load_pickle
    import forced_secretion as fs
    if alpha is None:
        alpha = json.loads(CHOICE.read_text())["chosen_alpha"]
    man = pd.read_csv(_scratch() / "manifest.csv", dtype={"sample_id": str})
    row = man.iloc[index]
    sid = str(row["sample_id"])
    out = _scratch() / "results" / sid
    if (out / "status.json").exists():
        print(f"{sid}: already solved"); return
    out.mkdir(parents=True, exist_ok=True)
    reg = pd.read_csv(REGISTRY, sep="\t")
    targets = set(reg.loc[reg["export_gate"] == True, "kegg"].dropna())  # noqa: E712
    t0 = time.time()
    status = {"sample_id": sid, "alpha": alpha, "solver": highs_qp.describe()}
    try:
        com = load_pickle(str(_scratch() / "models" / row["file"]))
        com.solver.configuration.presolve = True
        med = _medium()
        com.medium = med[med.index.isin([r.id for r in com.exchanges])]
        df = fs.forced_exports(com, alpha, targets)
        df.insert(0, "sample_id", sid)
        df.to_csv(out / "exports.tsv", sep="\t", index=False)
        g = pd.Series(df.attrs["growth_rates"], name="growth_rate").rename_axis("taxon")
        g.reset_index().assign(sample_id=sid).to_csv(out / "growth.tsv", sep="\t", index=False)
        status.update({k: v for k, v in df.attrs.items() if k != "growth_rates"})
        status.update({"ok": True, "qp": dict(highs_qp.LAST_RUN)})
    except Exception as err:
        status.update({"ok": False, "error": repr(err)})
    status["seconds"] = time.time() - t0
    (out / "status.json").write_text(json.dumps(status, indent=1, default=str))
    print(json.dumps({k: status.get(k) for k in ("sample_id", "ok", "seconds", "n_candidates")}))


def arrays(threshold: int = 61) -> None:
    """SLURM --array index lists for the two memory tiers, from the built manifest.

    Measured (memory_probe.py, one fresh solve at alpha 0.5): 0.84 GB at 20 taxa, 2.22 GB
    at 61; linear, 0.034 GB per taxon, so the cohort maximum of 77 taxa predicts 2.76 GB.
    Communities with <= `threshold` taxa run at --mem=2650M (billed 3 cores), the rest at
    --mem=3500M (billed 4).
    """
    man = pd.read_csv(_scratch() / "manifest.csv", dtype={"sample_id": str})
    n = man["found_taxa"].astype(int)
    small = man.index[n <= threshold].tolist()
    large = man.index[n > threshold].tolist()

    def compact(ix):
        out, i = [], 0
        while i < len(ix):
            j = i
            while j + 1 < len(ix) and ix[j + 1] == ix[j] + 1:
                j += 1
            out.append(f"{ix[i]}-{ix[j]}" if j > i else f"{ix[i]}")
            i = j + 1
        return ",".join(out)
    print(f"# {len(small)} samples <= {threshold} taxa, {len(large)} larger")
    print(f"sbatch -A naiss2026-3-773 --mem=2650M --array={compact(small)}%32 micom_batch.sbatch")
    if large:
        print(f"sbatch -A naiss2026-3-773 --mem=3500M --array={compact(large)} micom_batch.sbatch")


def collect() -> None:
    import utils
    res = _scratch() / "results"
    stat = pd.DataFrame([json.loads(p.read_text()) for p in res.glob("*/status.json")])
    ex = pd.concat([pd.read_csv(p, sep="\t", dtype={"sample_id": str})
                    for p in res.glob("*/exports.tsv")], ignore_index=True)
    gr = pd.concat([pd.read_csv(p, sep="\t", dtype={"sample_id": str})
                    for p in res.glob("*/growth.tsv")], ignore_index=True)
    tax = pd.read_csv(BATCH / "taxonomy.tsv", sep="\t", dtype={"sample_id": str})
    id2sp = dict(zip(tax["id"].str.replace(" ", "_"), tax["species"]))

    present = tax.groupby("species")["sample_id"].nunique()
    ex["species"] = ex["taxon"].map(id2sp)
    per = (ex.groupby(["species", "kegg"])
             .agg(samples_capable=("sample_id", "nunique"),
                  samples_forced=("forced", "sum"),
                  median_export_pfba=("export_pfba", "median"))
             .reset_index())
    per["samples_present"] = per["species"].map(present)
    per["frac_forced"] = per["samples_forced"] / per["samples_present"]
    per["meets_sample_rule"] = (per["frac_forced"] >= MIN_FRACTION) & (per["samples_forced"] >= MIN_SAMPLES)
    elig = per[per["samples_present"] >= MIN_SAMPLES]
    base = elig.groupby("kegg")["meets_sample_rule"].mean().rename("forced_base_rate")
    per = per.merge(base, on="kegg", how="left")
    per["Source_mechanistic_available"] = per["samples_present"] >= MIN_SAMPLES
    per["Source_mechanistic_support"] = (per["Source_mechanistic_available"]
                                         & per["meets_sample_rule"]
                                         & (per["forced_base_rate"] <= utils.EVIDENCE_MAX_BASE_RATE))
    BATCH.mkdir(parents=True, exist_ok=True)
    ex.to_csv(BATCH / "micom_forced_exports_v2.tsv.gz", sep="\t", index=False)
    gr.to_csv(BATCH / "micom_growth_v2.tsv", sep="\t", index=False)
    stat.to_csv(BATCH / "micom_solve_status_v2.tsv", sep="\t", index=False)
    per.to_csv(BATCH / "class_m.tsv", sep="\t", index=False)
    gap = stat.get("qp", pd.Series(dtype=object)).map(
        lambda q: q.get("gap_rel", np.nan) if isinstance(q, dict) else np.nan)
    print(json.dumps({"samples_solved": int(stat["ok"].sum()), "samples_failed": int((~stat["ok"]).sum()),
                      "uncertified_gt_1pct": int((gap > highs_qp.CERT_TOL).sum()),
                      "pairs_available": int(per["Source_mechanistic_available"].sum()),
                      "pairs_supported": int(per["Source_mechanistic_support"].sum())}, indent=1))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("step", choices=["prepare", "build", "solve", "arrays", "collect"])
    ap.add_argument("--index", type=int, help="solve: row of the manifest (SLURM_ARRAY_TASK_ID)")
    ap.add_argument("--alpha", type=float, help="solve: override the sweep's chosen alpha")
    ap.add_argument("--model-db", help="build: AGORA2 .qza or extracted folder")
    ap.add_argument("--threads", type=int, default=8)
    a = ap.parse_args()
    if a.step == "prepare":
        prepare()
    elif a.step == "build":
        build(a.model_db, a.threads)
    elif a.step == "arrays":
        arrays()
    elif a.step == "solve":
        solve(a.index if a.index is not None else int(os.environ["SLURM_ARRAY_TASK_ID"]), a.alpha)
    else:
        collect()
