#!/usr/bin/env python3
"""Measure how much of the curated enzyme set HUMAnN's UniRef90->EC map can actually see.

    python validate_mapping.py --map /path/to/map_level4ec_uniref90.txt.gz

This is the check that decides whether the mapped-gene-family route can carry the
chapter's specific claims, and it has to run before any conclusion is drawn from the
44 GB matrix.

The question it answers. `extract_targets.py` found four UniRef90 clusters mapped to
EC 4.1.1.105, and neither of the two named bacterial tryptophan decarboxylases --
A7B1V0 (RUMGNA_01526, R. gnavus) and J7SZ64 (CLOSPO_02083, C. sporogenes) -- is one of
their representatives. Results 3.7 rests its mechanistic agreement for tryptamine on
exactly those enzymes. So either they are non-representative MEMBERS of one of those
four clusters, in which case the route sees them, or their clusters are not annotated
4.1.1.105 in the map, in which case EC 4.1.1.105 by this route does not measure the
bacterial activity the chapter is about.

A cluster is named for its representative member only, so representativeness is the
wrong test. The right one is membership, and UniProt will state it directly: every
UniProtKB entry carries its UniRef90 cluster as a cross-reference. This script asks for
that, for every accession in the curated arm, and reports per EC:

    recovered  curated sequences whose UniRef90 cluster the map assigns to the SAME EC
    reassigned curated sequences whose cluster is in the map, but under a different EC
    unmapped   curated sequences whose cluster the map does not carry at all

`recovered` is the sensitivity of this route against characterised enzymes. It is a
number worth reporting in its own right: it is the first thing a reviewer should ask of
any EC-level HUMAnN result, and it is rarely stated.

Nothing here is a claim about the samples. It is a claim about the instrument.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import importlib.util
import io
import json
import sys
import time
import urllib.parse
import urllib.request
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
READS = HERE.parent / "reads"
# UniProtKB does NOT expose a UniRef cross-reference as a return field -- every spelling
# of it (uniref_cluster_90, xref_uniref90, uniref90, cc_uniref) is rejected with
# "Invalid fields parameter value". The UniRef database has to be queried instead, where
# `uniprot_id:<acc> AND identity:0.9` returns the UniRef90 cluster containing that
# member. One accession per query, because a batched query returns the clusters but not
# which input produced which, and the `members` field is truncated at ten so it cannot
# be used to re-attribute them. 363 queries is a few minutes and leaves no ambiguity.
API = "https://rest.uniprot.org/uniref/search"

NAMED = {
    "A7B1V0": "RUMGNA_01526, R. gnavus tryptophan decarboxylase",
    "J7SZ64": "CLOSPO_02083, C. sporogenes tryptophan decarboxylase",
}


def load_targets() -> dict[str, tuple[str, str]]:
    src = READS / "build_enzyme_db.py"
    spec = importlib.util.spec_from_file_location("_build_enzyme_db", src)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.TARGETS


def fetch(query: str, fields: str, size: int = 10) -> str:
    params = {"query": query, "format": "tsv", "fields": fields, "size": str(size)}
    url = f"{API}?{urllib.parse.urlencode(params)}"
    req = urllib.request.Request(url, headers={"User-Agent": "ki-v2-genes/1.0"})
    for attempt in (1, 2, 3):
        try:
            with urllib.request.urlopen(req, timeout=180) as r:
                return r.read().decode("utf-8")
        except Exception as exc:
            if attempt == 3:
                raise
            print(f"      retry {attempt} after {exc!r}")
            time.sleep(attempt * 5)
    return ""


def curated() -> dict[str, set[str]]:
    """accession -> set of target ECs UniProt annotates it with."""
    path = READS / "enzymes" / "members_sprot.tsv"
    if not path.is_file():
        sys.exit(f"no curated arm at {path} -- run v2/reads/build_enzyme_db.py first")
    out: dict[str, set[str]] = defaultdict(set)
    with path.open(newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            acc = (row.get("accession") or "").strip()
            if not acc:
                continue
            raw = (row.get("ec") or "").replace(",", ";")
            for ec in (e.strip() for e in raw.split(";")):
                if ec:
                    out[acc].add(ec)
    return dict(out)


def cluster_of(accessions: list[str], cache: Path) -> dict[str, dict[str, str]]:
    """accession -> {cluster, name, size}, from UniRef, one query per accession.

    Cached to disk because this is ~363 network round trips and the rest of the script
    is iterated on far more often than UniProt changes.
    """
    found: dict[str, dict[str, str]] = {}
    if cache.is_file():
        with cache.open(newline="", encoding="utf-8") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                found[row["accession"]] = {"cluster": row["cluster"],
                                           "name": row.get("name", ""),
                                           "size": row.get("size", "")}
        print(f"  {len(found)} from cache {cache.name}")

    todo = [a for a in accessions if a not in found]
    for n, acc in enumerate(todo, 1):
        if n % 25 == 1 or n == len(todo):
            print(f"  UniRef: {n}/{len(todo)}")
        tsv = fetch(f"uniprot_id:{acc} AND identity:0.9", "id,name,count", size=5)
        rows = list(csv.DictReader(io.StringIO(tsv), delimiter="\t"))
        if rows:
            r = rows[0]
            found[acc] = {
                "cluster": (r.get("Cluster ID") or "").strip(),
                "name": (r.get("Cluster Name") or "").strip(),
                "size": (r.get("Size") or "").strip(),
            }
        else:
            found[acc] = {"cluster": "", "name": "", "size": ""}

    if todo:
        with cache.open("w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["accession", "cluster", "name", "size"])
            for acc in sorted(found):
                d = found[acc]
                w.writerow([acc, d["cluster"], d["name"], d["size"]])
        print(f"  cached {len(found)} to {cache.name}")
    return found


def ec_of_cluster(map_path: Path) -> dict[str, set[str]]:
    """cluster -> every EC the map assigns it. Full pass; the file is 26 MB."""
    out: dict[str, set[str]] = defaultdict(set)
    with gzip.open(map_path, "rt") as fh:
        for line in fh:
            ec, _, rest = line.rstrip("\n").partition("\t")
            for c in rest.split("\t"):
                if c:
                    out[c].add(ec)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--map", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=HERE / "mapping_validation.tsv")
    args = ap.parse_args()

    targets = load_targets()
    cur = curated()
    print(f"=== curated arm: {len(cur)} accessions carrying a target EC ===")

    print("=== resolving each accession's UniRef90 cluster at UniProt ===")
    info = cluster_of(sorted(cur), HERE / "uniref90_cache.tsv")
    clusters = {a: d["cluster"] for a, d in info.items() if d["cluster"]}
    print(f"  resolved {len(clusters)} of {len(cur)}")

    print(f"=== reading {args.map.name} ===")
    cluster_ec = ec_of_cluster(args.map)
    print(f"  map carries {len(cluster_ec)} distinct clusters")

    rows: list[dict] = []
    per_ec: dict[str, dict[str, int]] = {ec: {"curated": 0, "recovered": 0,
                                              "reassigned": 0, "unmapped": 0,
                                              "no_cluster": 0}
                                         for ec in targets}
    for acc, ecs in sorted(cur.items()):
        cl = clusters.get(acc, "")
        mapped = cluster_ec.get(cl, set()) if cl else set()
        for ec in ecs:
            if ec not in per_ec:
                continue          # a secondary EC outside our ten
            per_ec[ec]["curated"] += 1
            if not cl:
                verdict = "no_cluster"
            elif ec in mapped:
                verdict = "recovered"
            elif mapped:
                verdict = "reassigned"
            else:
                verdict = "unmapped"
            per_ec[ec][verdict] += 1
            rows.append({"accession": acc, "ec": ec, "uniref90": cl,
                         "map_ecs": ";".join(sorted(mapped)), "verdict": verdict})

    print()
    hdr = f"{'EC':<10} {'curated':>7} {'recov':>6} {'reasg':>6} {'unmap':>6} {'none':>5} {'sens':>7}  activity"
    print(hdr)
    print("-" * len(hdr))
    for ec, (label, _claim) in targets.items():
        s = per_ec[ec]
        sens = (100.0 * s["recovered"] / s["curated"]) if s["curated"] else float("nan")
        print(f"{ec:<10} {s['curated']:>7} {s['recovered']:>6} {s['reassigned']:>6} "
              f"{s['unmapped']:>6} {s['no_cluster']:>5} {sens:>6.0f}%  {label}")

    tot = {k: sum(s[k] for s in per_ec.values())
           for k in ("curated", "recovered", "reassigned", "unmapped", "no_cluster")}
    overall = 100.0 * tot["recovered"] / tot["curated"] if tot["curated"] else 0.0
    print("-" * len(hdr))
    print(f"{'ALL':<10} {tot['curated']:>7} {tot['recovered']:>6} {tot['reassigned']:>6} "
          f"{tot['unmapped']:>6} {tot['no_cluster']:>5} {overall:>6.0f}%")

    # The curated route's target list: the UniRef90 clusters that actually contain
    # characterised enzymes. This is what the matrix gets queried for alongside the EC
    # route, and for the claims the EC map cannot see it is the ONLY route.
    cur_out = HERE / "curated_clusters.tsv"
    with cur_out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["uniref90", "cluster_name", "cluster_size",
                    "ecs", "accessions", "in_ec_map"])
        by_cluster: dict[str, dict] = {}
        for acc, ecs in cur.items():
            cl = clusters.get(acc, "")
            if not cl:
                continue
            d = by_cluster.setdefault(cl, {"ecs": set(), "accs": set(),
                                           "name": info[acc]["name"],
                                           "size": info[acc]["size"]})
            d["ecs"] |= {e for e in ecs if e in targets}
            d["accs"].add(acc)
        for cl in sorted(by_cluster):
            d = by_cluster[cl]
            w.writerow([cl, d["name"], d["size"], ";".join(sorted(d["ecs"])),
                        ";".join(sorted(d["accs"])),
                        "yes" if cl in cluster_ec else "no"])
    n_only = sum(1 for cl in by_cluster if cl not in cluster_ec)
    print(f"\ncurated route: {len(by_cluster)} distinct UniRef90 clusters hold the "
          f"{len(clusters)} characterised sequences")
    print(f"  of those, {n_only} carry NO EC in HUMAnN's map and are therefore "
          f"invisible to EC regrouping")
    print(f"wrote {cur_out}")

    print("\n=== the two named tryptophan decarboxylase loci ===")
    for acc, who in NAMED.items():
        cl = clusters.get(acc, "")
        mapped = sorted(cluster_ec.get(cl, set())) if cl else []
        print(f"  {acc}  {who}")
        print(f"    UniRef90 cluster : {cl or '(UniProt returned none)'}"
              f"  [{info.get(acc, {}).get('name', '')}, "
              f"{info.get(acc, {}).get('size', '?')} members]")
        print(f"    map assigns it   : {', '.join(mapped) if mapped else 'NOTHING'}")
        if "4.1.1.105" in mapped:
            print("    -> the mapped route CAN see this enzyme under EC 4.1.1.105")
        elif mapped:
            print("    -> visible, but under a DIFFERENT EC: measuring 4.1.1.105 "
                  "alone would miss it")
        else:
            print("    -> the mapped route CANNOT see this enzyme at all")

    with args.out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, delimiter="\t", lineterminator="\n",
                           fieldnames=["accession", "ec", "uniref90", "map_ecs", "verdict"])
        w.writeheader()
        w.writerows(rows)

    (HERE / "mapping_validation_summary.json").write_text(json.dumps({
        "built_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "map": str(args.map),
        "curated_source": "v2/reads/enzymes/members_sprot.tsv",
        "definitions": {
            "recovered": "curated sequence's UniRef90 cluster is mapped to the same EC",
            "reassigned": "cluster is in the map but under other ECs only",
            "unmapped": "cluster carries no EC in the map",
            "no_cluster": "UniProt returned no UniRef90 cross-reference",
        },
        "per_ec": per_ec,
        "overall": tot,
        "sensitivity_pct": round(overall, 1),
        "curated_clusters": len(by_cluster),
        "curated_clusters_absent_from_ec_map": n_only,
        "named_loci": {a: {"cluster": clusters.get(a, ""),
                           "cluster_name": info.get(a, {}).get("name", ""),
                           "cluster_size": info.get(a, {}).get("size", ""),
                           "map_ecs": sorted(cluster_ec.get(clusters.get(a, ""), set()))}
                       for a in NAMED},
    }, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
