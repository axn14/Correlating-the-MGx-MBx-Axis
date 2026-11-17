"""Producing enzymes for the stage-associated targets, from AGORA2's directional reactions.

Class G (CONTRACTS 5d) asks whether reads of a *producing* enzyme for the target are
attributed to the species in Yachida. Until 2026-11-13 it covered only the a-priori class
targets (curated decarboxylase / arginine-route enzymes, `compositional.EC_PRODUCTS`), so it
was silent on every one of the 30 stage-associated targets. This script supplies their
producing enzymes, under a rule fixed before any read is counted (CONTRACTS 5d-v):

  An EC produces target T if, in at least one of the 1,440 AGORA2 reconstructions, a
  reaction annotated with that complete (four-level) EC can NET-PRODUCE T in an allowed
  direction: T has a positive coefficient and upper_bound > 0, or a negative coefficient
  and lower_bound < 0.

Excluded: exchange / demand / sink reactions (one metabolite, or ids EX_, DM_, sink_),
transport (T on both sides of the reaction, in any compartment), and partial ECs
("2.6.1.-"), which HUMAnN's level-4 EC map cannot carry.

Why AGORA2 and not KEGG: KEGG's compound-enzyme links and reaction equations carry no
direction, so a consuming enzyme would be counted as a producer. AGORA2's bounds encode
curated reversibility. A reversible reaction counts in both directions, which is true of
the enzyme. Ubiquitous enzymes (aminotransferases, core biosynthesis) are expected; class G's
base-rate gate (<= utils.EVIDENCE_MAX_BASE_RATE, 5d-v) is what stops them voting.

Independence: the reference only defines WHICH enzymes count. Whether a species carries
one is read from HUMAnN reads in these patients, as for every other class-G call.

Output: stage_enzymes.tsv (kegg, target, ec, n_models, n_reactions, example_reactions)
        stage_enzymes_provenance.json
Run:    python v2/genes/stage_enzymes.py
"""
from __future__ import annotations

import json
import re
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
REGISTRY = HERE.parent / "targets" / "out" / "target_registry.tsv"
MODEL_ROOT = REPO / "Data" / "micom" / "agora2_extracted"
FULL_EC = re.compile(r"^\d+\.\d+\.\d+\.\d+$")


def model_dir() -> Path:
    hits = sorted(MODEL_ROOT.glob("*/data/manifest.csv"))
    if len(hits) != 1:
        sys.exit(f"expected one AGORA2 manifest under {MODEL_ROOT}, found {len(hits)}")
    return hits[0].parent


def ecs_of(reaction: dict) -> set[str]:
    raw = str(reaction.get("annotation", {}).get("ec-code", ""))
    parts = re.split(r"[;,\s]+|\bor\b|\band\b", raw)
    return {p.strip() for p in parts if FULL_EC.match(p.strip())}


def main() -> None:
    reg = pd.read_csv(REGISTRY, sep="\t", dtype=str)
    stage = reg[(reg["label"] == "stage_associated") & reg["kegg"].notna()]
    targets = dict(zip(stage["kegg"], stage["short_name"]))
    keggs = set(targets)
    print(f"{len(targets)} stage-associated targets with a KEGG id "
          f"({(reg['label'] == 'stage_associated').sum() - len(targets)} without)")

    d = model_dir()
    paths = sorted(p for p in d.glob("*.json"))
    models_by: dict[tuple[str, str], set[str]] = defaultdict(set)
    rxns_by: dict[tuple[str, str], set[str]] = defaultdict(set)
    for i, p in enumerate(paths, 1):
        m = json.loads(p.read_text(encoding="utf-8"))
        # pooled by base id, exactly as build_registry.agora2_census does
        by_base: dict[str, set[str]] = {}
        for met in m["metabolites"]:
            raw = str(met.get("annotation", {}).get("kegg.compound", "NaN"))
            ks = {x.strip() for x in raw.split(";")} & keggs
            if ks:
                by_base.setdefault(met["id"].split("[")[0], set()).update(ks)
        if not by_base:
            continue
        for r in m["reactions"]:
            rid = r["id"]
            mets = r["metabolites"]
            if len(mets) < 2 or rid.startswith(("EX_", "DM_", "sink_")):
                continue
            ecs = ecs_of(r)
            if not ecs:
                continue
            lb, ub = float(r.get("lower_bound", 0)), float(r.get("upper_bound", 0))
            net: dict[str, float] = defaultdict(float)
            sides: dict[str, set[int]] = defaultdict(set)
            for met_id, coef in mets.items():
                for k in by_base.get(met_id.split("[")[0], ()):
                    net[k] += coef
                    sides[k].add(1 if coef > 0 else -1)
            for k, c in net.items():
                if len(sides[k]) > 1:          # on both sides: transport, not synthesis
                    continue
                if (c > 0 and ub > 0) or (c < 0 and lb < 0):
                    for ec in ecs:
                        models_by[(k, ec)].add(p.stem)
                        rxns_by[(k, ec)].add(rid)
        if i % 200 == 0:
            print(f"  {i}/{len(paths)} models")

    rows = [{"kegg": k, "target": targets[k], "ec": ec,
             "n_models": len(models_by[(k, ec)]), "n_reactions": len(rxns_by[(k, ec)]),
             "example_reactions": ";".join(sorted(rxns_by[(k, ec)])[:5])}
            for (k, ec) in models_by]
    out = pd.DataFrame(rows).sort_values(["target", "n_models"], ascending=[True, False])
    out.to_csv(HERE / "stage_enzymes.tsv", sep="\t", index=False)

    per = out.groupby("target")["ec"].nunique()
    none = sorted(set(targets.values()) - set(per.index))
    prov = {"script": "v2/genes/stage_enzymes.py",
            "run_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "models_scanned": len(paths), "model_dir": str(d.relative_to(REPO)),
            "targets": len(targets), "targets_with_producing_ec": int(per.size),
            "targets_without": none, "ec_rows": len(out),
            "distinct_ecs": int(out["ec"].nunique()),
            "per_target": per.to_dict()}
    (HERE / "stage_enzymes_provenance.json").write_text(json.dumps(prov, indent=1))
    print(json.dumps(prov, indent=1))


if __name__ == "__main__":
    main()
