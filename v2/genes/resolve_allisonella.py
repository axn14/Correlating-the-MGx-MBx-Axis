#!/usr/bin/env python3
"""Settle the Allisonella question, with its own positive control.

    python resolve_allisonella.py

Results §3.7 ranks *Allisonella* first for cadaverine while noting that the genus is
characterised by histidine decarboxylation. Whether it carries a lysine decarboxylase at
all is therefore a claim the chapter makes and cannot currently support.

The extraction answered the easy half. *Allisonella histaminiformans* is profiled in this
cohort with **1,487 stratified gene families**, so it is not a coverage gap -- and none of
those 1,487 is in any of the 416 target lysine decarboxylase clusters.

That is not yet a finding, for a reason the data itself supplies. *Allisonella* is also
absent from `species_attribution.tsv` under **EC 4.1.1.22**, histidine decarboxylase --
the activity the organism is literally named for. The only species the EC route finds for
4.1.1.22 are Enterobacteriaceae carrying the P28577/P28578-type enzyme. *Allisonella*'s
histidine decarboxylase is the pyruvoyl-dependent HdcA family, a different protein
entirely, and HUMAnN's EC map does not connect it. So an EC-route absence for this
organism demonstrably means nothing on its own.

Which gives the negative claim a control that is hard to argue with:

    POSITIVE CONTROL  a histidine decarboxylase must be among these 1,487 clusters,
                      because the organism is defined by that activity
    THE CLAIM         is a lysine decarboxylase among them?

If the positive control fails, the probe cannot detect decarboxylases in this organism
and no absence claim is reportable. If it passes and lysine decarboxylase is still
absent, the absence is evidence.

Method, and why it is the reverse of `build_enzyme_db.py`. There, names were forbidden
and EC annotation was the gate, because the task was precision -- keeping cadmium pumps
out of a lysine decarboxylase set. Here the task is **recall**: finding any candidate the
EC map missed. So names are used to nominate candidates and EC annotation then confirms
each one individually. Names never decide the answer; they only decide what gets checked.
"""
from __future__ import annotations

import csv
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
OUT = HERE / "out"
UNIREF = "https://rest.uniprot.org/uniref/search"
UNIPROT = "https://rest.uniprot.org/uniprotkb/search"
BATCH = 40

# Words that nominate a cluster for an EC check. Deliberately wide: a false nomination
# costs one confirmation query, a missed one costs the claim.
NOMINATE = ("decarboxylase", "decarboxylating", "carboxy-lyase")

# The activities we care about, and the EC that confirms each.
LYSINE = {"4.1.1.18"}
HISTIDINE = {"4.1.1.22"}


def get(url_base: str, params: dict) -> str:
    url = f"{url_base}?{urllib.parse.urlencode(params)}"
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


def resolve(ids: list[str], cache: Path) -> dict[str, dict[str, str]]:
    """UniRef90 cluster id -> {protein, ec, organism, reviewed} for its representative.

    Straight to UniProtKB on the representative accession, rather than asking UniRef for
    the cluster name first. The first version did it the other way and 253 of 1,487
    clusters came back nameless: these profiles were built against a 2021 ChocoPhlAn
    release and UniRef has since merged or retired those cluster ids, so `id:UniRef90_X`
    finds nothing while the underlying accession is still there. Going at the accession
    recovers them, and returns the EC in the same pass -- so nomination and confirmation
    come from one query instead of two.

    `UniRef90_UPI...` clusters are representative-of-UniParc and have no UniProtKB entry;
    they are reported as unresolvable rather than silently counted as negative.
    """
    found: dict[str, dict[str, str]] = {}
    if cache.is_file():
        with cache.open(newline="", encoding="utf-8") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                found[row["uniref90"]] = {k: row.get(k, "") for k in
                                          ("representative", "protein", "ec",
                                           "organism", "reviewed", "uniref_name")}
        print(f"  {len(found)} from cache")

    todo = [i for i in ids if i not in found]
    accs = {c: c.split("_", 1)[1] for c in todo if "_" in c}
    kb = [(c, a) for c, a in sorted(accs.items()) if not a.startswith("UPI")]
    for k in range(0, len(kb), BATCH):
        chunk = kb[k:k + BATCH]
        q = " OR ".join(f"accession:{a}" for _, a in chunk)
        print(f"  UniProtKB: {k + 1}-{k + len(chunk)} of {len(kb)}")
        tsv = get(UNIPROT, {"query": f"({q})", "format": "tsv", "size": str(len(chunk)),
                            "fields": "accession,ec,protein_name,organism_name,reviewed"})
        by_acc = {(r.get("Entry") or "").strip(): r
                  for r in csv.DictReader(io.StringIO(tsv), delimiter="\t")}
        for cid, acc in chunk:
            r = by_acc.get(acc) or {}
            found[cid] = {
                "representative": acc,
                "protein": (r.get("Protein names") or "").strip(),
                "ec": (r.get("EC number") or "").strip(),
                "organism": (r.get("Organism") or "").strip(),
                "reviewed": (r.get("Reviewed") or "").strip(),
                "uniref_name": "",
            }
    for c in todo:
        found.setdefault(c, {"representative": c.split("_", 1)[1] if "_" in c else "",
                             "protein": "", "ec": "", "organism": "",
                             "reviewed": "", "uniref_name": ""})

    if todo:
        with cache.open("w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["uniref90", "representative", "protein", "ec",
                        "organism", "reviewed", "uniref_name"])
            for c in sorted(found):
                d = found[c]
                w.writerow([c, d["representative"], d["protein"], d["ec"],
                            d["organism"], d["reviewed"], d["uniref_name"]])
        print(f"  cached {len(found)}")
    return found


def reference_has(activity_ec: str, name: str, organism_substr: str) -> dict[str, int]:
    """Does UniProt hold this activity for this organism AT ALL?

    This is what separates the two ways a control can fail. If the probe misses an
    enzyme UniProt knows about, that is a detection failure and the depth or the
    pangenome is to blame. If UniProt has no such enzyme for the organism, then no
    UniRef90 cluster for it exists, so HUMAnN could never report one and the control
    cannot pass by construction -- which makes the whole question unanswerable by any
    UniRef90-based route, rather than answered in the negative.
    """
    out = {}
    for label, query in (
        ("by_ec", f"(ec:{activity_ec}) AND (organism_name:{organism_substr})"),
        ("by_name", f'(protein_name:"{name}") AND (organism_name:{organism_substr})'),
        ("by_gene", f"(gene:hdcA) AND (organism_name:{organism_substr})"),
        ("any_entry", f"organism_name:{organism_substr}"),
    ):
        tsv = get(UNIPROT, {"query": query, "format": "tsv",
                            "fields": "accession", "size": "500"})
        out[label] = max(0, len(tsv.strip().splitlines()) - 1)
    return out


def main() -> int:
    src = OUT / "allisonella_genes.tsv"
    if not src.is_file():
        sys.exit(f"no {src} -- run extract_gene_families.sbatch on Dardel first")

    rows = list(csv.DictReader(src.open(newline="", encoding="utf-8"), delimiter="\t"))
    ids = sorted({r["uniref90"] for r in rows})
    prev = {r["uniref90"]: (int(r["n_samples_nz"]), float(r["max"])) for r in rows}
    species = sorted({r["species"] for r in rows})
    print(f"=== {len(rows)} stratified rows, {len(ids)} distinct clusters ===")
    for s in species:
        print(f"  {s}")

    print("=== resolving every cluster's representative at UniProtKB ===")
    conf = resolve(ids, HERE / "allisonella_resolved.tsv")
    named = sum(1 for c in ids if conf.get(c, {}).get("protein"))
    unres = [c for c in ids if not conf.get(c, {}).get("protein")]
    print(f"  resolved {named} of {len(ids)}; {len(unres)} unresolvable")
    if unres:
        upi = sum(1 for c in unres if "_UPI" in c)
        print(f"    of those, {upi} are UniParc-representative clusters with no "
              f"UniProtKB entry")
        print(f"    and {len(unres) - upi} have an accession UniProtKB did not return")
        print("    -- these are a genuine blind spot in what follows, not a negative")

    names = {c: conf.get(c, {}).get("protein", "") for c in ids}
    nominated = [c for c in ids
                 if any(w in names.get(c, "").lower() for w in NOMINATE)]
    print(f"\n=== {len(nominated)} clusters nominated by protein name ===")
    for c in nominated:
        n, _mx = prev.get(c, (0, 0.0))
        print(f"  {c:<24} {names[c][:64]:<64} {n:>3}/616")

    # Separately: anything mentioning histidine at all, whether or not the name says
    # decarboxylase. The positive control is too important to rest on one keyword --
    # A. histaminiformans carries the pyruvoyl-dependent HdcA family, which UniProt
    # sometimes names by its proenzyme or by a related activity.
    hist_any = [c for c in ids if "histidine" in names.get(c, "").lower()]
    print(f"\n=== {len(hist_any)} clusters whose protein name mentions histidine ===")
    for c in hist_any:
        n, _mx = prev.get(c, (0, 0.0))
        print(f"  {c:<24} EC {conf[c]['ec'] or '-':<12} {names[c][:52]:<52} {n:>3}/616")

    if not nominated and not hist_any:
        print("  nothing nominated at all.")

    lys, his, other = [], [], []
    for c in set(nominated) | set(hist_any):
        d = conf.get(c, {})
        ecs = {e.strip() for e in d.get("ec", "").replace(",", ";").split(";") if e.strip()}
        nm = d.get("protein", "").lower()
        if ecs & LYSINE:
            lys.append(c)
        elif (ecs & HISTIDINE) or ("histidine decarboxylase" in nm):
            # Name admitted here only for the control, and only for this exact phrase.
            # The control's job is to show the probe CAN see a decarboxylase in this
            # organism; a pyruvoyl-dependent HdcA frequently carries the name without
            # a 4.1.1.22 annotation in TrEMBL, and refusing it on that basis would make
            # the control fail for a reason that has nothing to do with the organism.
            his.append(c)
        else:
            other.append(c)
    lys, his, other = sorted(lys), sorted(his), sorted(other)

    def show(title: str, cs: list[str]) -> None:
        print(f"\n  {title}: {len(cs)}")
        for c in cs:
            d = conf.get(c, {})
            n, mx = prev.get(c, (0, 0.0))
            print(f"    {c:<24} EC {d.get('ec', '') or '-':<18} {n:>3}/616  max {mx:.3g}")
            print(f"      {d.get('protein', '')[:92]}")
            print(f"      {d.get('organism', '')[:70]}  [{d.get('reviewed', '')}]")

    show("histidine decarboxylase, EC 4.1.1.22 -- THE POSITIVE CONTROL", his)
    show("lysine decarboxylase, EC 4.1.1.18 -- THE CLAIM", lys)
    show("other decarboxylases", other)

    ref = {}
    if not his:
        print("\n=== the control failed, so: does the REFERENCE even have the enzyme? ===")
        # MetaPhlAn calls the organism Allisonella histaminiformans; UniProt's entries
        # for it are filed under Dialister histaminiformans, the same organism after the
        # genus was reassigned. Matching on the species epithet covers both.
        ref = reference_has("4.1.1.22", "histidine decarboxylase", "histaminiformans")
        print(f"  UniProt entries for *histaminiformans*        : {ref['any_entry']}")
        print(f"  ... annotated EC 4.1.1.22                     : {ref['by_ec']}")
        print(f"  ... named 'histidine decarboxylase'           : {ref['by_name']}")
        print(f"  ... carrying gene hdcA                        : {ref['by_gene']}")

    print("\n" + "=" * 78)
    if not his and ref and ref["any_entry"] > 0 and ref["by_ec"] == 0 \
            and ref["by_name"] == 0 and ref["by_gene"] == 0:
        print("THE QUESTION IS NOT ANSWERABLE BY THIS ROUTE, and the reason is clean.")
        print()
        print(f"UniProt holds {ref['any_entry']}+ entries for this organism and NOT ONE is")
        print("annotated as a histidine decarboxylase -- not by EC, not by protein name,")
        print("not by gene symbol -- for the organism whose species epithet means")
        print("'histamine-forming'. So no UniRef90 cluster for that enzyme exists, HUMAnN")
        print("could never report one, and the positive control cannot pass by")
        print("construction.")
        print()
        print("That makes the absence of a lysine decarboxylase among all",
              len(ids), "of this")
        print("organism's gene families UNINTERPRETABLE. Results 3.7's caveat about")
        print("Allisonella and cadaverine must stand as a caveat: this evidence neither")
        print("supports nor refutes it, and reporting the absence as a finding would be")
        print("reporting a gap in UniProt as a fact about the gut.")
        print()
        print("What would answer it: assembly and binning of Allisonella from the reads")
        print("to recover its actual gene complement -- the Paper 2 work in")
        print("writeup/dardel_plan.md section 2 -- or a targeted experiment. Not this.")
        verdict = "unanswerable_reference_lacks_control_enzyme"
    elif not his:
        print("POSITIVE CONTROL FAILED. No histidine decarboxylase was confirmed in")
        print("A. histaminiformans, an organism defined by that activity, although")
        print("UniProt does hold such an enzyme for it. That is a detection failure --")
        print("depth, or the ChocoPhlAn pangenome -- so the absence of a lysine")
        print("decarboxylase is uninformative and must NOT be reported as a finding.")
        verdict = "control_failed_detection"
    elif lys:
        print("Allisonella DOES carry a lysine decarboxylase in this cohort. Results")
        print("3.7's cadaverine ranking is consistent with a direct route, and the")
        print("caveat about the genus being a histidine decarboxylase specialist")
        print("should be revised rather than retained.")
        verdict = "lysine_dc_present"
    else:
        print("POSITIVE CONTROL PASSED and NO lysine decarboxylase was found.")
        print("The probe detects decarboxylases in this organism -- it found the")
        print("histidine decarboxylase the genus is named for -- and it finds no")
        print("lysine decarboxylase among all", len(ids), "of its gene families.")
        print("Results 3.7 ranks Allisonella first for cadaverine; this says the")
        print("organism most likely does not make cadaverine directly, so that")
        print("ranking reflects correlation rather than a producer relationship.")
        print("That is a correction to the chapter, supported by a control.")
        verdict = "lysine_dc_absent_control_passed"
    print("=" * 78)

    print("\nThe limits of this, stated rather than left implied:")
    print("  - it measures DNA. A gene is not expression and not flux.")
    print("  - the positive control licenses an absence claim about THIS organism in")
    print("    THIS cohort, not about the genus in general.")
    print("  - the representative's EC stands in for the cluster. A cluster can hold")
    print("    members of mixed annotation; UniRef90 is 90% identity, which is tight")
    print("    but not a guarantee of shared function.")
    print("  - only Yachida has shotgun data, so nothing here replicates.")

    with (OUT / "allisonella_verdict.tsv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["uniref90", "representative", "ec", "protein",
                    "organism", "reviewed", "n_samples_nz", "max", "role"])
        for role, cs in (("positive_control_histidine_dc", his),
                         ("claim_lysine_dc", lys), ("other_decarboxylase", other)):
            for c in cs:
                d = conf.get(c, {})
                n, mx = prev.get(c, (0, 0.0))
                w.writerow([c, d.get("representative", ""), d.get("ec", ""),
                            d.get("protein", ""), d.get("organism", ""),
                            d.get("reviewed", ""), n, mx, role])

    (OUT / "allisonella_verdict.json").write_text(json.dumps({
        "built_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "species": species,
        "stratified_rows": len(rows),
        "distinct_clusters": len(ids),
        "clusters_resolved": named,
        "clusters_unresolvable": len(unres),
        "nominated_by_name": len(nominated),
        "name_mentions_histidine": len(hist_any),
        "confirmed_histidine_dc": len(his),
        "confirmed_lysine_dc": len(lys),
        "other_decarboxylases": len(other),
        "reference_check": ref,
        "verdict": verdict,
    }, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {OUT / 'allisonella_verdict.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
