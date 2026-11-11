"""The target registry: one row per trinity target, with the gates CONTRACTS 5d-ii needs.

For each of the targets NB03 ran (the section 5c set, recomputed per run), this records:

  label            `stage_associated` (a Robust DA hit, Yachida, Healthy vs Early or
                   Advanced CRC, NB02) or `a_priori_class` (in the pre-specified
                   decarboxylation / polyamine class only). The 33 vs 49 question, answered
                   per row.
  direction        sign of NB02's MeanDiff (CRC minus Healthy) in each robust contrast
  confounder gate  `evidence_tools.confounder_unattributable` on NB03's own AUCs
  export gate      does any cohort species' AGORA2 model carry an EXCHANGE reaction for the
                   compound (matched by the model's own `kegg.compound` annotation)? A
                   compound no model can export has no modelled microbial export route, so
                   MICOM (class M, and R's MICOM part) is unavailable for it.
  capability       how many cohort models touch the compound at all (any compartment)
  host enzymes     KEGG: EC numbers linked to the compound that have human (hsa) genes.
                   Direction-agnostic, exactly like AGORA2 capability: the enzyme may
                   consume the compound. It answers "can the host act on this", not "does
                   the host make it".
  origin_flag      microbial_only     no human enzyme linked to the compound
                   host_and_microbial human enzyme(s) and >= 1 cohort model touches it
                   host_dominant      human enzyme(s), no cohort model touches it
                   unannotated        neither KEGG-human nor AGORA2 annotates it
                   unassessable       no KEGG id
  note             the explicit biological caveats CONTRACTS 5d-ii names

Run:  python v2/targets/build_registry.py --results "D:/KI Research-dardel/<run_id>"
      (the directory holding Results/tables). KEGG REST answers are cached in
      v2/targets/out/kegg_cache.json, so a re-run is offline and reproducible.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
import urllib.request
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
V2 = HERE.parent
sys.path.insert(0, str(V2))

import utils            # noqa: E402
import evidence_tools as et  # noqa: E402

OUT = HERE / "out"
COHORT_MODELS = V2 / "micom" / "out" / "cohort_species_agora2.tsv"
KEGG_CACHE = OUT / "kegg_cache.json"
YACHIDA = "YACHIDA-CRC-2019"
CONTRASTS = {"Healthy_vs_Early_CRC": "early", "Healthy_vs_Advanced_CRC": "advanced"}

# The caveats CONTRACTS 5d-ii requires to be stated, and why.
NOTES = {
    "C00624": "N-acetylglutamate: in bacteria an intracellular intermediate of arginine "
              "biosynthesis (argA), not a recognised export product; host enterocytes also "
              "make it (NAGS). Faecal levels may reflect lysis or host origin.",
    "C05135": "N-acetylhistamine (Yachida label '4-(beta-acetylaminoethyl)imidazole'): a "
              "histamine catabolite; host mast-cell and epithelial histamine turnover can "
              "supply it.",
}
POLYAMINE_NOTE = ("polyamine: host ODC1/SAT1/SMS/SRM metabolism and tumour tissue are "
                  "major sources; failed or near-failed the confounder screen in NB03.")


def kegg_of(target: str) -> str | None:
    head = str(target).split("_", 1)[0]
    return head if head.startswith("C") and head[1:].isdigit() else None


def load_tables(results: Path):
    t = results / "Results" / "tables"
    bench = pd.read_csv(t / "ml_benchmark_dysregulated_metabolites_v2.csv")
    da = pd.read_csv(t / "differential_abundance_v2.csv")
    return bench, da


def label_targets(bench: pd.DataFrame, da: pd.DataFrame) -> pd.DataFrame:
    targets = sorted(bench["Metabolite"].unique())
    y = da[(da["Dataset"] == YACHIDA) & (da["Robust"] == True)  # noqa: E712
           & (da["Assay"] == "metabolites") & (da["Scheme"] == "grouped")
           & (da["Comparison"].isin(CONTRASTS))]
    rows = []
    class_kegg = set(utils.DECARBOXYLATION_TARGET_KEGG) | set(utils.POLYAMINE_KEGG)
    for tgt in targets:
        k = kegg_of(tgt)
        hits = y[y["Feature"] == tgt]
        robust_in = {}
        for _, h in hits.iterrows():
            robust_in[CONTRASTS[h["Comparison"]]] = "up" if float(h["MeanDiff"]) > 0 else "down"
        b = bench[bench["Metabolite"] == tgt].iloc[0]
        rows.append({
            "target": tgt, "kegg": k, "short_name": b.get("Short_name", tgt),
            "label": "stage_associated" if robust_in else "a_priori_class",
            "in_prespecified_class": k in class_kegg,
            "robust_early": robust_in.get("early", ""),
            "robust_advanced": robust_in.get("advanced", ""),
            "auc_model": float(b["AUC_nested"]),
            "auc_confounders_adjusted": float(b["AUC_confounders_adjusted"]),
            "perm_q": float(b["Perm_qvalue"]),
            "confounder_unattributable": et.confounder_unattributable(
                b["AUC_nested"], b["AUC_confounders_adjusted"]),
        })
    return pd.DataFrame(rows)


def agora2_census(keggs: set[str]) -> tuple[Counter, Counter, int]:
    """(models exporting, models touching, models read) per KEGG id, cohort species only."""
    cs = pd.read_csv(COHORT_MODELS, sep="\t")
    paths = [p for p in cs.loc[cs["in_agora2"] == True, "path"].dropna()]  # noqa: E712
    export, touch = Counter(), Counter()
    for p in paths:
        with open(p, encoding="utf-8") as f:
            m = json.load(f)
        # An annotation can hold several ids ("C00025;C00217" for glutamate), so it is
        # split, and a metabolite maps to every target id it carries.
        # Annotations are not always repeated on every compartment's copy, so they are
        # pooled by base id ("citr_L" for citr_L[c] and citr_L[e]).
        by_base: dict[str, set[str]] = {}
        for met in m["metabolites"]:
            raw = str(met.get("annotation", {}).get("kegg.compound", "NaN"))
            ks = {x.strip() for x in raw.split(";")} & keggs
            if ks:
                by_base.setdefault(met["id"].split("[")[0], set()).update(ks)
        kegg_by_met = {met["id"]: by_base[met["id"].split("[")[0]]
                       for met in m["metabolites"] if met["id"].split("[")[0] in by_base}
        touched = set().union(*kegg_by_met.values()) if kegg_by_met else set()
        exported = set()
        for r in m["reactions"]:
            if r["id"].startswith("EX_"):
                for met_id in r["metabolites"]:
                    exported |= kegg_by_met.get(met_id, set())
        for k in touched:
            touch[k] += 1
        for k in exported:
            export[k] += 1
    return export, touch, len(paths)


def _kegg_get(path: str, cache: dict) -> str:
    if path in cache:
        return cache[path]
    with urllib.request.urlopen(f"https://rest.kegg.jp/{path}", timeout=30) as r:
        txt = r.read().decode()
    cache[path] = txt
    time.sleep(0.4)                      # KEGG asks for <= 3 requests per second
    return txt


def human_enzymes(keggs: list[str]) -> dict[str, list[str]]:
    cache = json.loads(KEGG_CACHE.read_text()) if KEGG_CACHE.exists() else {}
    out = {}
    for k in keggs:
        ecs = sorted({ln.split("\t")[1].replace("ec:", "")
                      for ln in _kegg_get(f"link/enzyme/cpd:{k}", cache).splitlines()
                      if "\t" in ln})
        human = []
        for i in range(0, len(ecs), 10):
            chunk = ecs[i:i + 10]
            txt = _kegg_get("link/hsa/" + "+".join(f"ec:{e}" for e in chunk), cache)
            human += sorted({ln.split("\t")[0].replace("ec:", "")
                             for ln in txt.splitlines() if "\thsa:" in ln})
        out[k] = sorted(set(human))
    cache["_fetched_utc"] = cache.get("_fetched_utc") or \
        datetime.now(timezone.utc).isoformat(timespec="seconds")
    KEGG_CACHE.write_text(json.dumps(cache, indent=0, sort_keys=True))
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", type=Path, required=True,
                    help="run root holding Results/tables (e.g. a pulled Dardel run)")
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    bench, da = load_tables(a.results)
    reg = label_targets(bench, da)
    keggs = sorted({k for k in reg["kegg"] if k})
    export, touch, n_models = agora2_census(set(keggs))
    human = human_enzymes(keggs)

    reg["agora2_models_exporting"] = reg["kegg"].map(lambda k: export.get(k, 0) if k else pd.NA)
    reg["agora2_models_touching"] = reg["kegg"].map(lambda k: touch.get(k, 0) if k else pd.NA)
    reg["export_gate"] = reg["agora2_models_exporting"].map(
        lambda n: (bool(n > 0)) if pd.notna(n) else pd.NA)
    reg["human_ec"] = reg["kegg"].map(lambda k: ";".join(human.get(k, [])) if k else "")

    def origin(r):
        if not r["kegg"]:
            return "unassessable"
        has_h = bool(r["human_ec"])
        has_m = (r["agora2_models_touching"] or 0) > 0
        if not has_h and not has_m:
            return "unannotated"          # neither side annotates it: no call either way
        if not has_h:
            return "microbial_only"
        return "host_and_microbial" if has_m else "host_dominant"
    reg["origin_flag"] = reg.apply(origin, axis=1)
    reg["attribution_requires_lysis"] = reg["export_gate"].map(
        lambda g: (not g) if pd.notna(g) else pd.NA)

    def note(r):
        parts = []
        if r["kegg"] in NOTES:
            parts.append(NOTES[r["kegg"]])
        if r["kegg"] in utils.POLYAMINE_KEGG:
            parts.append(POLYAMINE_NOTE)
        if r["confounder_unattributable"]:
            parts.append("confounders alone match or beat the species model (NB03).")
        return " ".join(parts)
    reg["note"] = reg.apply(note, axis=1)

    reg.to_csv(OUT / "target_registry.tsv", sep="\t", index=False)
    prov = {"script": "v2/targets/build_registry.py",
            "run_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "results_root": str(a.results), "n_targets": int(len(reg)),
            "labels": reg["label"].value_counts().to_dict(),
            "origin_flags": reg["origin_flag"].value_counts().to_dict(),
            "export_gate_failed": reg.loc[reg["export_gate"] == False, "target"].tolist(),  # noqa: E712
            "agora2_cohort_models_read": n_models,
            "kegg_cache": str(KEGG_CACHE.name)}
    (OUT / "target_registry_provenance.json").write_text(json.dumps(prov, indent=1, default=str))
    print(json.dumps(prov, indent=1, default=str))


if __name__ == "__main__":
    main()
