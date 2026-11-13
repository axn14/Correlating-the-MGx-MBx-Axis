"""Nofinska et al. 2026: which deposited subjects are new to this study? (goals doc, step 1)

Public metadata only; no reads are downloaded. Nofinska's metagenomes are in DDBJ
PRJDB19991 and PRJDB37508. Every BioSample there is titled "... patient NNNNN", the
National Cancer Center numbering that Yachida 2019 (PRJDB4176, "human feces_NNNNN") and
our own Yachida/Erawijantari tables (`Subject`) also use. A shared number is therefore
read as a shared person, and that reading is checked where it can be:

  1. Within Yachida, the numbering must be a stable person ID. Some PRJDB4176 BioSamples
     state sex and age (most leave them blank); for every one of those also in our Yachida
     table, sex must agree and age must agree within 1 year. Asserted at >= 95% of >= 50.
  2. The Nofinska BioSamples carry no age or sex, so their overlap rests on the
     number alone. Chance agreement is reported as the overlap expected if Nofinska's
     numbers were drawn at random from the same numeric range.

Independence follows the Erawijantari precedent (`erawijantari_species.py`): a Nofinska
subject counts as independent only if its number is in NONE of PRJDB4176, our Yachida
table or our Erawijantari table. The data request asks the authors to confirm this.

Run:  python v2/replication/nofinska_overlap.py [--fetch]
      --fetch re-downloads the four ENA reports (about 0.3 MB) into nofinska_raw/;
      otherwise the committed copies are used.
Out:  out/nofinska_runs.tsv           one row per run, with the overlap flags
      out/nofinska_subjects.tsv       one row per patient number
      out/nofinska_overlap_summary.json
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
V2 = HERE.parent
sys.path.insert(0, str(V2))
import utils  # noqa: E402

RAW = HERE / "nofinska_raw"
OUT = HERE / "out"
ENA = "https://www.ebi.ac.uk/ena/portal/api"
RUN_FIELDS = ("run_accession,sample_accession,sample_title,library_strategy,instrument_model,"
              "read_count,base_count,fastq_bytes,first_public")
REPORTS = {
    "PRJDB19991_runs.tsv": f"{ENA}/filereport?accession=PRJDB19991&result=read_run&fields={RUN_FIELDS}&format=tsv&limit=0",
    "PRJDB37508_runs.tsv": f"{ENA}/filereport?accession=PRJDB37508&result=read_run&fields={RUN_FIELDS}&format=tsv&limit=0",
    "PRJDB4176_runs.tsv": f"{ENA}/filereport?accession=PRJDB4176&result=read_run&fields={RUN_FIELDS}&format=tsv&limit=0",
    "PRJDB4176_samples.tsv": f"{ENA}/search?result=sample&query=study_accession=PRJDB4176&fields=sample_accession,sample_title,sex,age&format=tsv&limit=0",
}
NOFINSKA = ("PRJDB19991", "PRJDB37508")
# NCC numbers, with PRJDB4176's repeat samples suffixed ".1"; other sites keep their own
# label (PRJDB37508 has two "KAGAWA-n" samples), which cannot collide with an NCC number.
PATIENT = re.compile(r"(?:patient|feces_)\s*(\d{4,6})(?:\.\d+)?\s*$")
OTHER_SITE = re.compile(r"patient\s+([A-Z]+-\d+)\s*$")


def fetch() -> None:
    RAW.mkdir(parents=True, exist_ok=True)
    for name, url in REPORTS.items():
        req = urllib.request.Request(url, headers={"User-Agent": "generic-client"})
        with urllib.request.urlopen(req, timeout=120) as r:
            (RAW / name).write_bytes(r.read())
        print(f"fetched {name}")


def patient_no(title: str) -> str | None:
    m = PATIENT.search(str(title)) or OTHER_SITE.search(str(title))
    return m.group(1) if m else None


def gb(fastq_bytes: str) -> float:
    return sum(int(x) for x in str(fastq_bytes).split(";") if x.strip().isdigit()) / 1e9


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fetch", action="store_true")
    a = ap.parse_args()
    if a.fetch or not all((RAW / n).exists() for n in REPORTS):
        fetch()
    OUT.mkdir(parents=True, exist_ok=True)

    runs = {p: pd.read_csv(RAW / f"{p}_runs.tsv", sep="\t", dtype=str)
            for p in (*NOFINSKA, "PRJDB4176")}
    for p, df in runs.items():
        df["patient"] = df["sample_title"].map(patient_no)
        df["project"] = p
        assert df["patient"].notna().all(), f"{p}: a sample title has no patient number"
    old = runs["PRJDB4176"]
    old = old[old["library_strategy"] == "WGS"]
    in_4176 = set(old["patient"])

    yach = pd.read_csv(utils.DATA_DIR / "YACHIDA-CRC-2019 metadata.tsv", sep="\t", dtype=str)
    eraw = pd.read_csv(utils.DATA_DIR / "ERAWIJANTARI-GASTRIC-CANCER-2020 metadata.tsv",
                       sep="\t", dtype=str)
    yach["patient"] = yach["Subject"].str.extract(r"^(\d+)")[0]
    eraw["patient"] = eraw["Subject"].str.extract(r"^(\d+)")[0]
    y_group = yach.set_index("patient")["Study.Group"]
    e_group = eraw.groupby("patient")["Study.Group"].agg(lambda s: "/".join(sorted(set(s))))

    # 1. The numbering is a stable person ID within the Yachida deposits.
    smp = pd.read_csv(RAW / "PRJDB4176_samples.tsv", sep="\t", dtype=str)
    smp["patient"] = smp["sample_title"].map(patient_no)
    chk = smp.merge(yach[["patient", "Gender", "Age"]], on="patient", how="inner")
    # Most PRJDB4176 BioSamples leave sex and age blank; compare only those that state both.
    chk = chk[chk["sex"].notna() & chk["age"].notna()]
    sex_ok = chk["sex"].str.lower() == chk["Gender"].str.lower()
    age_ok = (pd.to_numeric(chk["age"], errors="coerce")
              - pd.to_numeric(chk["Age"], errors="coerce")).abs() <= 1
    id_agreement = float((sex_ok & age_ok).mean())
    assert len(chk) >= 50 and id_agreement >= 0.95, (
        f"patient numbers are not a stable ID: {id_agreement:.1%} of {len(chk)} agree")

    # 2. Flag every Nofinska run.
    new = pd.concat([runs[p] for p in NOFINSKA], ignore_index=True)
    new["fastq_gb"] = new["fastq_bytes"].map(gb)
    new["in_PRJDB4176"] = new["patient"].isin(in_4176)
    new["in_yachida_table"] = new["patient"].isin(set(yach["patient"]))
    new["yachida_group"] = new["patient"].map(y_group)
    new["in_erawijantari_table"] = new["patient"].isin(set(eraw["patient"]))
    new["erawijantari_group"] = new["patient"].map(e_group)
    new["independent"] = ~(new["in_PRJDB4176"] | new["in_yachida_table"]
                           | new["in_erawijantari_table"])
    cols = ["project", "run_accession", "sample_accession", "patient", "read_count",
            "fastq_gb", "in_PRJDB4176", "in_yachida_table", "yachida_group",
            "in_erawijantari_table", "erawijantari_group", "independent", "first_public"]
    new[cols].to_csv(OUT / "nofinska_runs.tsv", sep="\t", index=False)

    subj = new.groupby("patient").agg(
        projects=("project", lambda s: "/".join(sorted(set(s)))), runs=("run_accession", "size"),
        fastq_gb=("fastq_gb", "sum"), in_PRJDB4176=("in_PRJDB4176", "first"),
        in_yachida_table=("in_yachida_table", "first"), yachida_group=("yachida_group", "first"),
        in_erawijantari_table=("in_erawijantari_table", "first"),
        independent=("independent", "first")).reset_index()
    subj.to_csv(OUT / "nofinska_subjects.tsv", sep="\t", index=False)

    # Chance overlap: Nofinska numbers drawn at random from the same numeric range.
    nums = pd.to_numeric(subj["patient"], errors="coerce").dropna()
    lo, hi = int(nums.min()), int(nums.max())
    ref = {int(x) for x in in_4176 | set(yach["patient"]) | set(eraw["patient"].dropna())}
    ref_in_range = sum(lo <= x <= hi for x in ref)
    expected_chance = len(subj) * ref_in_range / (hi - lo + 1)

    ind = subj[subj["independent"]]
    summary = {
        "script": "v2/replication/nofinska_overlap.py",
        "run_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "sources": {k: v for k, v in REPORTS.items()},
        "id_check_PRJDB4176_vs_yachida_table": {
            "subjects_compared": int(len(chk)), "sex_and_age_agree": round(id_agreement, 4)},
        "nofinska_runs": {p: int((new["project"] == p).sum()) for p in NOFINSKA},
        "nofinska_subjects": int(len(subj)),
        "fastq_gb_total": round(float(new["fastq_gb"].sum()), 1),
        "subjects_in_PRJDB4176": int(subj["in_PRJDB4176"].sum()),
        "subjects_in_our_yachida_table": int(subj["in_yachida_table"].sum()),
        "our_yachida_overlap_by_group": subj.loc[subj["in_yachida_table"], "yachida_group"]
                                        .value_counts().to_dict(),
        "subjects_in_our_erawijantari_table": int(subj["in_erawijantari_table"].sum()),
        "independent_subjects": int(len(ind)),
        "independent_by_project": ind["projects"].value_counts().to_dict(),
        "independent_fastq_gb": round(float(ind["fastq_gb"].sum()), 1),
        "non_ncc_subjects": sorted(subj.loc[~subj["patient"].str.isdigit(), "patient"]),
        "number_range": [lo, hi],
        "overlap_expected_by_chance": round(expected_chance, 1),
        "overlap_observed_any": int((~subj["independent"]).sum()),
        "caveat": ("Nofinska BioSamples carry no age or sex, so their overlap rests on the "
                   "shared patient number. Independence is provisional until the authors' "
                   "metadata confirms it."),
    }
    (OUT / "nofinska_overlap_summary.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
