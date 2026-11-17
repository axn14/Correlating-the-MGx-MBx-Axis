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


def _samples() -> list[str]:
    """The batch's sample order: sorted ids from taxonomy.tsv, so an array index means
    the same sample on every machine and in every resubmission."""
    tax = pd.read_csv(BATCH / "taxonomy.tsv", sep="\t", dtype={"sample_id": str})
    return sorted(tax["sample_id"].unique())


def _build_one(sid: str) -> tuple[Path, float]:
    """Build ONE sample's community model in this process, unless it already exists.

    2026-11-12: the separate 6-worker build job (27475925) sat at 19.3 of 20 GB for 93 min
    without writing a model -- workers thrashing at the memory ceiling. Building inside the
    per-sample task needs one model in memory at a time (2.38 GB peak at 61 taxa, measured
    locally) and isolates a failure to its sample. Each sample gets its own output folder
    because micom.workflows.build also writes a manifest.csv there, and concurrent tasks
    must not share one.
    """
    from micom.workflows import build as mbuild
    import diet_screen as ds
    folder = _scratch() / "models" / sid
    model = folder / f"{sid}.pickle"
    if model.exists():
        return model, 0.0
    db = os.environ.get("MICOM_MODEL_DB")
    if not db:
        raise RuntimeError("set MICOM_MODEL_DB to the extracted AGORA2 folder (manifest.csv)")
    tax = pd.read_csv(BATCH / "taxonomy.tsv", sep="\t", dtype={"sample_id": str})
    folder.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    mbuild(tax[tax["sample_id"] == sid], out_folder=str(folder), model_db=db,
           cutoff=ds.ABUND_CUTOFF, threads=1)
    return model, time.time() - t0


def solve(index: int, alpha: float | None) -> None:
    from micom import load_pickle
    import forced_secretion as fs
    if alpha is None:
        alpha = json.loads(CHOICE.read_text())["chosen_alpha"]
    sid = _samples()[index]
    out = _scratch() / "results" / sid
    if (out / "status.json").exists():
        print(f"{sid}: already solved"); return
    out.mkdir(parents=True, exist_ok=True)
    reg = pd.read_csv(REGISTRY, sep="\t")
    targets = set(reg.loc[reg["export_gate"] == True, "kegg"].dropna())  # noqa: E712
    t0 = time.time()
    status = {"sample_id": sid, "alpha": alpha, "solver": highs_qp.describe()}
    try:
        model, build_s = _build_one(sid)
        status["build_seconds"] = build_s
        com = load_pickle(str(model))
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
    try:
        import psutil
        mi = psutil.Process().memory_info()
        status["peak_gb"] = round(getattr(mi, "peak_wset", 0) / 1e9, 3) or None
    except Exception:
        pass
    if status.get("peak_gb") is None:
        import resource
        status["peak_gb"] = round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6, 3)
    (out / "status.json").write_text(json.dumps(status, indent=1, default=str))
    print(json.dumps({k: status.get(k) for k in ("sample_id", "ok", "seconds", "n_candidates")}))


def arrays(threshold: int = 59) -> None:
    """SLURM --array index lists for the two memory tiers, from taxonomy.tsv's sample order.

    Measured 2026-11-12, build and solve as SEPARATE processes (the task's peak is the
    larger of the two): build 1.02 GB at 20 taxa, 2.38 at 61; solve 0.84 and 2.22. With
    15% headroom, <= 59 taxa fit --mem=2664M (billed 3 cores: 254 samples) and the rest
    --mem=3552M (billed 4: 23 samples, up to 77 taxa). Estimated 307 billed core-h.
    """
    tax = pd.read_csv(BATCH / "taxonomy.tsv", sep="\t", dtype={"sample_id": str})
    counts = tax.groupby("sample_id").size()
    n = pd.Series([int(counts[s]) for s in _samples()])
    small = n.index[n <= threshold].tolist()
    large = n.index[n > threshold].tolist()

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
    print(f"sbatch -A naiss2026-3-773 --mem=2664M --array={compact(small)}%32 micom_batch.sbatch")
    if large:
        print(f"sbatch -A naiss2026-3-773 --mem=3552M --array={compact(large)} micom_batch.sbatch")


def _class_m(ex: pd.DataFrame, tax: pd.DataFrame, samples: set[str], max_base_rate: float) -> pd.DataFrame:
    """Class M per (species, target) over the given samples only (CONTRACTS 5d)."""
    e = ex[ex["sample_id"].isin(samples)]
    t = tax[tax["sample_id"].isin(samples)]
    present = t.groupby("species")["sample_id"].nunique()
    per = (e.groupby(["species", "kegg"])
             .agg(samples_capable=("sample_id", "nunique"),
                  samples_forced=("forced", "sum"),
                  median_export_pfba=("export_pfba", "median"))
             .reset_index())
    per["samples_present"] = per["species"].map(present).fillna(0).astype(int)
    per["frac_forced"] = per["samples_forced"] / per["samples_present"].where(per["samples_present"] > 0)
    per["meets_sample_rule"] = (per["frac_forced"] >= MIN_FRACTION) & (per["samples_forced"] >= MIN_SAMPLES)
    elig = per[per["samples_present"] >= MIN_SAMPLES]
    base = elig.groupby("kegg")["meets_sample_rule"].mean().rename("forced_base_rate")
    per = per.merge(base, on="kegg", how="left")
    per["available"] = per["samples_present"] >= MIN_SAMPLES
    per["support"] = (per["available"] & per["meets_sample_rule"]
                      & (per["forced_base_rate"] <= max_base_rate))
    return per


def collect() -> None:
    """Concatenate the batch and call class M.

    CONTRACTS 5d requires the certified gap <= 1% (highs_qp.CERT_TOL) for a sample to count
    toward class M. On the 277-sample batch (2026-11-13) 137 samples met it and 139 did not
    (79 above 5%), so class M is called on the certified samples only, with every count's
    denominator restricted to them. The same call at the 5% bound (solver_check P3) is
    written alongside as a NON-VOTING sensitivity column, so the cost of the rule is visible.
    """
    import utils
    res = _scratch() / "results"
    stat = pd.DataFrame([json.loads(p.read_text()) for p in res.glob("*/status.json")])
    stat["sample_id"] = stat["sample_id"].astype(str)
    stat["gap_rel"] = stat.get("qp", pd.Series(dtype=object)).map(
        lambda q: q.get("gap_rel", np.nan) if isinstance(q, dict) else np.nan).astype(float)
    stat["certified_1pct"] = stat["ok"].astype(bool) & (stat["gap_rel"] <= highs_qp.CERT_TOL)
    stat["certified_5pct"] = stat["ok"].astype(bool) & (stat["gap_rel"] <= 0.05)
    ex = pd.concat([pd.read_csv(p, sep="\t", dtype={"sample_id": str})
                    for p in res.glob("*/exports.tsv")], ignore_index=True)
    gr = pd.concat([pd.read_csv(p, sep="\t", dtype={"sample_id": str})
                    for p in res.glob("*/growth.tsv")], ignore_index=True)
    tax = pd.read_csv(BATCH / "taxonomy.tsv", sep="\t", dtype={"sample_id": str})
    id2sp = dict(zip(tax["id"].str.replace(" ", "_"), tax["species"]))
    ex["species"] = ex["taxon"].map(id2sp)
    unmapped = int(ex["species"].isna().sum())

    s1 = set(stat.loc[stat["certified_1pct"], "sample_id"])
    s5 = set(stat.loc[stat["certified_5pct"], "sample_id"])
    m1 = _class_m(ex, tax, s1, utils.EVIDENCE_MAX_BASE_RATE)
    m5 = _class_m(ex, tax, s5, utils.EVIDENCE_MAX_BASE_RATE)
    per = m1.rename(columns={"available": "Source_mechanistic_available",
                             "support": "Source_mechanistic_support"})
    per = per.merge(m5[["species", "kegg", "available", "support"]].rename(
        columns={"available": "sens5pct_available", "support": "sens5pct_support"}),
        on=["species", "kegg"], how="outer")
    for c in ("Source_mechanistic_available", "Source_mechanistic_support",
              "sens5pct_available", "sens5pct_support"):
        per[c] = per[c].eq(True)

    BATCH.mkdir(parents=True, exist_ok=True)
    ex.to_csv(BATCH / "micom_forced_exports_v2.tsv.gz", sep="\t", index=False)
    gr.to_csv(BATCH / "micom_growth_v2.tsv", sep="\t", index=False)
    stat.drop(columns=[c for c in ("qp", "solver") if c in stat]).to_csv(
        BATCH / "micom_solve_status_v2.tsv", sep="\t", index=False)
    per.to_csv(BATCH / "class_m.tsv", sep="\t", index=False)
    summary = {"samples": int(len(stat)), "solved": int(stat["ok"].sum()),
               "failed": stat.loc[~stat["ok"].astype(bool), "sample_id"].tolist(),
               "certified_1pct": len(s1), "certified_5pct": len(s5),
               "exports_unmapped_taxa": unmapped,
               "pairs_available": int(per["Source_mechanistic_available"].sum()),
               "pairs_supported": int(per["Source_mechanistic_support"].sum()),
               "sens5pct_pairs_available": int(per["sens5pct_available"].sum()),
               "sens5pct_pairs_supported": int(per["sens5pct_support"].sum())}
    (BATCH / "class_m_summary.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary, indent=1))

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("step", choices=["prepare", "build", "build-one", "solve", "arrays", "collect"])
    ap.add_argument("--index", type=int, help="solve: row of the manifest (SLURM_ARRAY_TASK_ID)")
    ap.add_argument("--alpha", type=float, help="solve: override the sweep's chosen alpha")
    ap.add_argument("--model-db", help="build: AGORA2 .qza or extracted folder")
    ap.add_argument("--threads", type=int, default=8)
    a = ap.parse_args()
    if a.step == "prepare":
        prepare()
    elif a.step == "build":
        build(a.model_db, a.threads)
    elif a.step == "build-one":
        # its own process, so the solve that follows starts clean: measured locally, a
        # build+solve in ONE process peaks at 3.20 GB (61 taxa) against 2.22 for the
        # solve alone -- the build's objects are not all released.
        idx = a.index if a.index is not None else int(os.environ["SLURM_ARRAY_TASK_ID"])
        sid = _samples()[idx]
        model, secs = _build_one(sid)
        try:
            import psutil
            pk = getattr(psutil.Process().memory_info(), "peak_wset", 0) / 1e9
        except Exception:
            pk = 0
        if not pk:
            import resource
            pk = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6
        print(json.dumps({"sample_id": sid, "build_seconds": round(secs, 1),
                          "build_peak_gb": round(pk, 3), "model": str(model)}))
    elif a.step == "arrays":
        arrays()
    elif a.step == "solve":
        solve(a.index if a.index is not None else int(os.environ["SLURM_ARRAY_TASK_ID"]), a.alpha)
    else:
        collect()
