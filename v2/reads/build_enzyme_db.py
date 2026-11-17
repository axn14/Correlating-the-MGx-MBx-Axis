#!/usr/bin/env python3
"""Build the curated enzyme reference set the targeted read search runs against.

    python build_enzyme_db.py              # writes enzymes/ next to this file
    python build_enzyme_db.py --unreviewed # also take TrEMBL, for coverage

This is the scientific core of the read programme, so membership is decided by
annotation and not by name matching. Every sequence in the output carries one of the
target EC numbers in UniProt's own `ec` field, or is one of a short list of named
accessions admitted deliberately. Nothing enters because its description happened to
contain a string.

That rule exists because the first version of this script did match on strings, and
81% of what it produced was wrong. `cadA` is the gene symbol for inducible lysine
decarboxylase AND for the cadmium-transporting ATPase, so a free-text query returned
eleven cadmium pumps and six CadC regulators; `gadB` dragged in GadX transcriptional
regulators; and a cadaverine/lysine antiporter came along too -- a transporter, not a
decarboxylase. None of that fails loudly. It quietly counts reads against the wrong
protein, which is the one error this probe cannot be allowed to make.

Two N-acetylglutamate routes are kept apart on purpose. ArgA (EC 2.3.1.1) forms
N-acetylglutamate de novo; ArgJ (EC 2.3.1.35) forms it by recycling the acetyl group
from N-acetylornithine. Both bear on the N-acetylglutamate finding in Discussion 4.4,
they are different mechanisms, and summing them into one "synthase" number would
misreport which one the data support. They are reported as separate classes.

What this set is NOT: a claim that these are the only organisms able to perform these
reactions. It is a probe. A read aligning to a characterised lysine decarboxylase says
lysine decarboxylase is in the sample; it does not say whose. Attributing a gene to an
organism needs assembly, which is the follow-on work.
"""
from __future__ import annotations

import argparse
import csv
import io
import json
import sys
import time
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / "enzymes"
API = "https://rest.uniprot.org/uniprotkb/stream"
BACTERIA = 2  # NCBI taxonomy id

# EC -> (class label, the metabolite claim it speaks to). These are the reactions the
# chapter makes claims about: the amino-acid decarboxylase class pre-specified in
# Methods, plus the polyamine arm, which the MICOM diet reformulation still needs even
# though spermine and N1-acetylspermine failed confounder control in NB03.
TARGETS: dict[str, tuple[str, str]] = {
    "4.1.1.105": ("tryptophan/aromatic decarboxylase", "tryptamine, phenethylamine"),
    "4.1.1.28":  ("aromatic-L-amino-acid decarboxylase", "tyramine, phenethylamine"),
    "4.1.1.25":  ("tyrosine decarboxylase", "tyramine, phenethylamine"),
    "4.1.1.18":  ("lysine decarboxylase", "cadaverine -- the Allisonella question"),
    "4.1.1.22":  ("histidine decarboxylase", "histamine"),
    "4.1.1.15":  ("glutamate decarboxylase", "GABA"),
    "2.3.1.1":   ("N-acetylglutamate synthase, de novo (ArgA)", "N-acetylglutamate"),
    "2.3.1.35":  ("N-acetylglutamate synthase, recycling (ArgJ)", "N-acetylglutamate"),
    "4.1.1.17":  ("ornithine decarboxylase", "putrescine"),
    "4.1.1.19":  ("arginine decarboxylase", "agmatine"),
}

# Gene symbols queried against UniProt's `gene` FIELD -- never as free text. A hit here
# still has to carry a target EC to be admitted, so this only widens recall within the
# target activities; it cannot introduce an off-target protein.
GENES = ["ldcC", "cadA", "hdcA", "hdcB", "tyrDC", "gadA", "gadB",
         "argA", "argJ", "speA", "speC", "speF"]

# Admitted by accession regardless of EC annotation, because the chapter's mechanistic
# argument names them specifically and an EC gate would drop one that is annotated only
# by activity description. Each is justified, and each is checked for presence below.
NAMED_LOCI = {
    "RUMGNA_01526": "R. gnavus tryptophan decarboxylase; Results 3.7's mechanistic "
                    "agreement for tryptamine, also active on phenylalanine",
    "CLOSPO_02083": "C. sporogenes tryptophan decarboxylase, the other characterised "
                    "member of the same family",
}

FIELDS = "accession,id,ec,gene_names,protein_name,organism_name,reviewed,lineage"


def fetch(query: str, fmt: str, fields: str | None = None) -> str:
    params = {"query": query, "format": fmt}
    if fields:
        params["fields"] = fields
    url = f"{API}?{urllib.parse.urlencode(params)}"
    req = urllib.request.Request(url, headers={"User-Agent": "ki-v2-reads/1.0"})
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


def read_fasta(text: str) -> dict[str, str]:
    """accession -> full fasta record."""
    out: dict[str, str] = {}
    for block in text.split("\n>"):
        block = block.strip()
        if not block:
            continue
        if not block.startswith(">"):
            block = ">" + block
        head = block.split("\n", 1)[0]
        parts = head[1:].split("|")
        acc = parts[1] if len(parts) >= 2 else head[1:].split()[0]
        out[acc] = block.rstrip() + "\n"
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--unreviewed", action="store_true",
                    help="include TrEMBL as well as Swiss-Prot (far larger, noisier)")
    args = ap.parse_args()
    OUT.mkdir(exist_ok=True)

    rev_clause = "" if args.unreviewed else " AND (reviewed:true)"
    ec_clause = " OR ".join(f"(ec:{ec})" for ec in TARGETS)
    gene_clause = " OR ".join(f"(gene:{g})" for g in GENES)
    query = f"(({ec_clause}) OR ({gene_clause})) AND (taxonomy_id:{BACTERIA}){rev_clause}"

    print("=== querying UniProt once, then gating on the EC annotation ===")
    print(f"  taxonomy : NCBI {BACTERIA} (Bacteria)")
    print(f"  reviewed : {'Swiss-Prot + TrEMBL' if args.unreviewed else 'Swiss-Prot only'}")
    meta_tsv = fetch(query, "tsv", FIELDS)
    rows = list(csv.DictReader(io.StringIO(meta_tsv), delimiter="\t"))
    print(f"  returned : {len(rows)} entries")

    fasta = read_fasta(fetch(query, "fasta"))
    print(f"  sequences: {len(fasta)}")

    # The gate. An entry is admitted only if UniProt annotates it with a target EC, or
    # it is one of the named loci. Everything else is recorded as rejected, with why.
    keep: dict[str, dict] = {}
    rejected: list[dict] = []
    per_ec: dict[str, list[str]] = defaultdict(list)
    named_found: dict[str, str] = {}

    for r in rows:
        acc = (r.get("Entry") or r.get("accession") or "").strip()
        ecs = {e.strip() for e in (r.get("EC number") or "").replace(";", " ").split() if e.strip()}
        hit = sorted(ecs & set(TARGETS))
        genes = (r.get("Gene Names") or "")
        locus = next((l for l in NAMED_LOCI if l.lower() in genes.lower()), None)
        if locus:
            named_found[locus] = acc
        if hit or locus:
            keep[acc] = r
            for ec in (hit or ["(named locus)"]):
                per_ec[ec].append(acc)
        else:
            rejected.append({"accession": acc,
                             "protein": r.get("Protein names", "")[:70],
                             "genes": genes[:40],
                             "ec": r.get("EC number", "")[:40]})

    print(f"\n=== the EC gate ===")
    print(f"  admitted : {len(keep)}")
    print(f"  rejected : {len(rejected)}  (no target EC, not a named locus)")
    if rejected:
        worst = Counter(x["protein"] for x in rejected).most_common(8)
        print("  the rejections, by protein name -- these are what string matching let in:")
        for name, c in worst:
            print(f"     {c:>4}  {name}")

    print(f"\n=== admitted sequences by activity ===")
    for ec in list(TARGETS) + ["(named locus)"]:
        accs = per_ec.get(ec, [])
        if not accs:
            # An empty "(named locus)" bucket is the good case, not a problem: it means
            # every named locus also carried a target EC and was counted under it.
            note = ("-- empty is correct here: the named loci carry target ECs"
                    if ec == "(named locus)" else "-- NONE, investigate")
            print(f"  EC {ec:<10} {0:>5}   {note}")
            continue
        label, claim = TARGETS.get(ec, ("named characterised loci", "see NAMED_LOCI"))
        print(f"  EC {ec:<10} {len(accs):>5}   {label}")
        print(f"  {'':<14}          -> {claim}")

    print(f"\n=== the named enzymes the chapter's claims rest on ===")
    ok = True
    for locus, why in NAMED_LOCI.items():
        acc = named_found.get(locus)
        if acc:
            print(f"  {locus:<14} PRESENT as {acc}")
        else:
            print(f"  {locus:<14} ABSENT   <-- the tryptamine claim has no probe without it")
            ok = False
        print(f"  {'':<14} {why}")

    # Write only the admitted sequences.
    missing_seq = [a for a in keep if a not in fasta]
    if missing_seq:
        print(f"\n  note: {len(missing_seq)} admitted entries had no sequence returned; "
              f"excluded from the FASTA")
    # Every output is suffixed by arm. The two arms are reported separately and never
    # merged -- the specific arm supports presence, the sensitivity arm says whether an
    # absence is real -- so sharing a filename between them silently destroys the
    # record of one. It did exactly that once.
    arm = "trembl" if args.unreviewed else "sprot"
    out_fa = OUT / f"enzymes_{arm}.faa"
    out_fa.write_text("".join(fasta[a] for a in keep if a in fasta), encoding="utf-8")
    n_written = sum(1 for a in keep if a in fasta)
    print(f"\nwrote {out_fa.name}: {n_written} sequences")

    with open(OUT / f"members_{arm}.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["accession", "ec", "gene_names", "protein_name", "organism"])
        for acc, r in sorted(keep.items()):
            w.writerow([acc, r.get("EC number", ""), r.get("Gene Names", ""),
                        r.get("Protein names", ""), r.get("Organism", "")])
    with open(OUT / f"rejected_{arm}.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, delimiter="\t", lineterminator="\n",
                          fieldnames=["accession", "protein", "genes", "ec"])
        w.writeheader()
        w.writerows(rejected)

    (OUT / f"provenance_{arm}.json").write_text(json.dumps({
        "built_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "api": API,
        "query": query,
        "reviewed_only": not args.unreviewed,
        "returned": len(rows),
        "admitted": len(keep),
        "rejected": len(rejected),
        "written": n_written,
        "per_ec": {k: len(v) for k, v in sorted(per_ec.items())},
        "named_loci": named_found,
        "targets": {k: list(v) for k, v in TARGETS.items()},
    }, indent=2), encoding="utf-8")
    print(f"wrote members_{arm}.tsv, rejected_{arm}.tsv, provenance_{arm}.json")

    if not ok:
        print("\nThe named loci are the probe for the chapter's central mechanistic")
        print("claim. Do not run the search until they are in the set.")
        return 1
    if n_written < 50:
        print("\nWARNING: fewer than 50 sequences; re-run with --unreviewed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
