"""Post-NB09 step: tier every pair by evidence CLASS (CONTRACTS 5d).

NB09 still writes its source-count `Tier_final`, kept for comparison. This step joins the
class inputs that are computed from the SAME run's outputs and so cannot live inside NB09:

  v2/targets/out/target_registry.tsv          (build_registry.py, from NB02/NB03 tables)
  v2/genes/out/class_g.tsv                    (genes/compositional.py)
  v2/replication/out/erawijantari_species_replication.tsv  (from NB08's pair table)
  v2/micom/out/batch/class_m.tsv              (micom_batch.py collect; optional)

Order on a fresh run: run_pipeline.py -> targets/build_registry.py --results <run> ->
replication/erawijantari_species.py --results <run> -> [MICOM batch -> collect] -> this.

Writes <run>/Results/tables/final_evidence_classes_v2.csv and prints the tier table and
the old-vs-new cross-tabulation.

Run:  python v2/assemble_classes.py --results "D:/KI Research-dardel/<run_id>"
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import evidence_tools as et  # noqa: E402

INPUTS = {
    "registry": HERE / "targets" / "out" / "target_registry.tsv",
    "class_g": HERE / "genes" / "out" / "class_g.tsv",
    "class_x": HERE / "replication" / "out" / "erawijantari_species_replication.tsv",
    "class_m": HERE / "micom" / "out" / "batch" / "class_m.tsv",
}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", type=Path, required=True)
    a = ap.parse_args()
    tables = a.results / "Results" / "tables"
    pairs = pd.read_csv(tables / "final_evidence_synthesis_metabolome_v2.csv", low_memory=False)
    got = {k: (pd.read_csv(p, sep="\t") if p.exists() else None) for k, p in INPUTS.items()}
    if got["registry"] is None:
        raise SystemExit("no target registry: run v2/targets/build_registry.py first")
    out = et.attach_class_evidence(pairs, got["registry"], got["class_g"],
                                   got["class_x"], got["class_m"])
    out.to_csv(tables / "final_evidence_classes_v2.csv", index=False)

    tgt = out[out["Target_label"] != "not_a_target"]
    tier = tgt["Tier_class"].str.extract(r"^([A-Za-z -]+?) \(")[0]
    print("inputs present:", {k: v is not None for k, v in got.items()})
    print(f"\n{len(out):,} pairs; {len(tgt):,} on the {tgt['Metabolite'].nunique()} targets")
    print("\nclass tier, target pairs:")
    print(tier.value_counts().to_string())
    print("\nold source tier (rows) x class tier (columns), target pairs:")
    print(pd.crosstab(tgt["Tier_final"].str.split(" ").str[0], tier).to_string())
    concl = tgt[tier == "Conclusive"]
    if len(concl):
        print("\nConclusive pairs:")
        print(concl[["Species", "Metabolite", "Rho_partial", "Classes_supporting",
                     "Target_label", "Target_origin_flag"]].to_string(index=False))
    summary = {"pairs": int(len(out)), "target_pairs": int(len(tgt)),
               "class_tiers": tier.value_counts().to_dict(),
               "inputs": {k: str(p) if got[k] is not None else None for k, p in INPUTS.items()}}
    (tables / "final_evidence_classes_summary_v2.json").write_text(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
