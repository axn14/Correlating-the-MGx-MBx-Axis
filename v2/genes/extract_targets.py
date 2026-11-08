#!/usr/bin/env python3
"""Turn the ten target EC numbers into the set of UniRef90 clusters to look for.

    python extract_targets.py                       # writes target_clusters.tsv here
    python extract_targets.py --map /path/to/map_level4ec_uniref90.txt.gz

The gene families in `profiles/gene_families.rda` are named `UniRef90_XXXXXX` for the
community total and `UniRef90_XXXXXX|g__Genus.s__Species` for the species-stratified
rows. This script produces the cluster ids that correspond to our enzyme activities, so
the 44 GB matrix can be subset to a few thousand rows in one pass and never loaded again.

The EC list is IMPORTED from v2/reads/build_enzyme_db.py rather than restated. That is
deliberate: the curated-sequence arm and the mapped-gene-family arm are then provably
asking about the same ten reactions, and the Swiss-Prot members become a usable check on
this mapping instead of a parallel set of claims about a slightly different list.

A note on what the mapping is and is not. HUMAnN's map assigns an EC to a UniRef90
cluster from UniProt annotation of the cluster's members, so it inherits UniProt's
annotation quality -- including the automatic TrEMBL annotations that the curated arm
deliberately holds at arm's length. A cluster mapped to EC 4.1.1.18 is therefore
"annotated as lysine decarboxylase", not "a characterised lysine decarboxylase". That
distinction is why the Swiss-Prot cross-check below exists and why both numbers get
reported, never merged.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import importlib.util
import json
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
READS = HERE.parent / "reads"
DEFAULT_MAP = Path(
    "/cfs/klemming/projects/supr/naiss2024-6-169/andnair/mappings/"
    "map_level4ec_uniref90.txt.gz"
)


def load_targets() -> dict[str, tuple[str, str]]:
    """Import TARGETS from the enzyme-probe script, so there is one declaration."""
    src = READS / "build_enzyme_db.py"
    if not src.is_file():
        sys.exit(f"cannot find {src}: the EC list lives there and is not restated here")
    spec = importlib.util.spec_from_file_location("_build_enzyme_db", src)
    mod = importlib.util.module_from_spec(spec)
    # build_enzyme_db.py does its network work inside main(), so importing is safe.
    spec.loader.exec_module(mod)
    return mod.TARGETS


def read_map(path: Path, wanted: set[str]) -> dict[str, list[str]]:
    """One line per EC: the EC, then tab-separated UniRef90 cluster ids."""
    found: dict[str, list[str]] = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            ec, _, rest = line.rstrip("\n").partition("\t")
            if ec in wanted:
                found[ec] = [c for c in rest.split("\t") if c]
                if len(found) == len(wanted):
                    break
    return found


def sprot_accessions() -> dict[str, set[str]]:
    """EC -> UniProt accessions from the curated arm, for the cross-check."""
    path = READS / "enzymes" / "members_sprot.tsv"
    if not path.is_file():
        return {}
    out: dict[str, set[str]] = defaultdict(set)
    with path.open(newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            acc = (row.get("accession") or "").strip()
            if not acc:
                continue
            # the EC column may hold several, semicolon- or comma-separated
            raw = (row.get("ec") or row.get("ec_number") or "").replace(",", ";")
            for ec in (e.strip() for e in raw.split(";")):
                if ec:
                    out[ec].add(acc)
    return dict(out)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--map", type=Path, default=DEFAULT_MAP)
    ap.add_argument("--out", type=Path, default=HERE / "target_clusters.tsv")
    args = ap.parse_args()

    if not args.map.is_file():
        sys.exit(f"no mapping at {args.map} -- run fetch_mappings.sh first")

    targets = load_targets()
    found = read_map(args.map, set(targets))

    missing = [ec for ec in targets if ec not in found]
    sprot = sprot_accessions()

    rows: list[tuple[str, str, str, str]] = []
    print(f"{'EC':<10} {'clusters':>9} {'sprot':>6}  {'repr?':>5}  activity")
    print(f"{'-' * 10} {'-' * 9:>9} {'-' * 6:>6}  {'-' * 5:>5}  {'-' * 8}")
    for ec, (label, claim) in targets.items():
        clusters = found.get(ec, [])
        accs = sprot.get(ec, set())
        # A UniRef90 cluster is named for its representative member, so a curated
        # accession that IS a representative shows up as UniRef90_<accession>. This
        # counts how many of the curated sequences are directly recoverable that way.
        # It is a lower bound on agreement: a curated sequence inside a cluster named
        # for a different representative is not counted here.
        direct = sum(1 for a in accs if f"UniRef90_{a}" in set(clusters))
        print(f"{ec:<10} {len(clusters):>9} {len(accs):>6}  {direct:>5}  {label}")
        for c in clusters:
            rows.append((ec, c, label, claim))

    if missing:
        print(f"\nNOT IN THE MAP: {', '.join(missing)}")
        print("  An EC absent from the map cannot be measured by this route at all.")
        print("  That is a statement about the mapping, not about the samples.")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["ec", "uniref90", "activity", "claim"])
        w.writerows(rows)

    uniq = {c for _, c, _, _ in rows}
    prov = {
        "built_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "map": str(args.map),
        "ec_source": "v2/reads/build_enzyme_db.py TARGETS (imported, not restated)",
        "targets": len(targets),
        "ecs_found": len(found),
        "ecs_missing": missing,
        "cluster_rows": len(rows),
        "clusters_unique": len(uniq),
        "per_ec": {ec: len(cs) for ec, cs in found.items()},
    }
    (args.out.parent / "target_clusters_provenance.json").write_text(
        json.dumps(prov, indent=2) + "\n", encoding="utf-8"
    )

    print(f"\n{len(rows)} (ec, cluster) pairs over {len(uniq)} distinct clusters")
    print(f"wrote {args.out}")
    # Clusters are shared between ECs where a protein carries more than one activity,
    # which is why the distinct count matters for the subset and the pair count does not.
    if len(uniq) != len(rows):
        print(f"  note: {len(rows) - len(uniq)} pairs are clusters shared across ECs")

    # The two named loci are the chapter's mechanistic argument for tryptamine, so
    # report explicitly whether this route can see them at all.
    print("\nthe named tryptophan decarboxylase loci, under EC 4.1.1.105:")
    for c in found.get("4.1.1.105", []):
        print(f"  {c}")
    for acc, who in (("A7B1V0", "RUMGNA_01526, R. gnavus"),
                     ("J7SZ64", "CLOSPO_02083, C. sporogenes")):
        tag = "present as a cluster representative" \
            if f"UniRef90_{acc}" in uniq else "NOT a cluster representative"
        print(f"  {acc} ({who}): {tag}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
