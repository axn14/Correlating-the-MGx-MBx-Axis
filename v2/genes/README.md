# `v2/genes` — gene-level evidence from the stratified profiles

Discussion §4.6 says why this directory exists: *"of the six declared sources, the two
that speak most directly to whether an organism **can** perform a reaction were either
uninformative by base rate or unavailable."* That is the reason the eight three-source
pairs are prioritised hypotheses rather than producer relationships.

`curatedMetagenomicData` distributes bioBakery 3 HUMAnN output for `YachidaS_2019`, and
it is **species-stratified**: row names are `UniRef90_XXXXXX` for the community total and
`UniRef90_XXXXXX|g__Genus.s__Species` for the per-species contribution. The attribution
is already in the row names, which is what the read programme in `../reads` would have
needed a separate assembly step to reach. See `../reads/README.md` for why that route is
superseded, and `writeup/read_programme_alternative.md` for the measurements behind it.

## Read this first: the EC map has a blind spot on the chapter's lead claims

The obvious way to turn 8.9 M UniRef90 gene families into enzyme activities is HUMAnN's
own `map_level4ec_uniref90`, which is what `humann_regroup_table` uses. We measured that
map against the 363 curated Swiss-Prot sequences in `../reads/enzymes/` before trusting
it, and it fails exactly where the thesis makes its strongest claims:

| EC | activity | curated seqs | recoverable via the EC map | sensitivity |
|---|---|---|---|---|
| 4.1.1.105 | tryptophan decarboxylase → **tryptamine** | 2 | 0 | **0%** |
| 4.1.1.28 | aromatic-L-amino-acid dc | 1 | 0 | **0%** |
| 4.1.1.25 | tyrosine decarboxylase → **tyramine** | 4 | 0 | **0%** |
| 4.1.1.18 | lysine dc → **cadaverine** | 8 | 3 | 38% |
| 4.1.1.15 | glutamate dc → GABA | 16 | 8 | 50% |
| 4.1.1.22 | histidine dc → histamine | 19 | 14 | 74% |
| 2.3.1.1 | ArgA, N-acetylglutamate de novo | 194 | 145 | 75% |
| 4.1.1.19 | arginine dc → agmatine | 114 | 92 | 81% |
| 2.3.1.35 | ArgJ, N-acetylglutamate recycling | 130 | 112 | 86% |
| 4.1.1.17 | ornithine dc → putrescine | 5 | 5 | 100% |
| | **overall** | **493** | **379** | **77%** |

The two characterised gut bacterial tryptophan decarboxylases are the clearest case:

- **`UniRef90_A7B1V0`** — A7B1V0 = `RUMGNA_01526`, *R. gnavus*; UniRef names the cluster
  "Tryptophan decarboxylase", 20 members, reviewed.
- **`UniRef90_J7SZ64`** — J7SZ64 = `CLOSPO_02083`, *C. sporogenes*; same cluster name,
  21 members, reviewed.

**Neither appears anywhere in `map_level4ec_uniref90`, nor in `map_ko_uniref90`.** The
four clusters the EC map *does* assign to 4.1.1.105 are `UniRef90_A0A286LEZ8`,
`UniRef90_P0DPA6`, `UniRef90_P93082` and `UniRef90_P93083` — none of them these enzymes.

So regrouping gene families to EC 4.1.1.105 and reporting the result would have measured
a set that excludes both enzymes Results §3.7 names, and a near-zero answer would have
read as *no bacterial tryptamine synthesis* when it actually means *the mapping does not
carry these enzymes*. Nothing about that failure is loud. It is the same class of error
as the `cadA` collision documented in `../reads/README.md`: a wrong number, not an error
message.

Of the 224 distinct UniRef90 clusters holding the 363 curated sequences, **28 carry no EC
in HUMAnN's map at all.**

## Therefore: two routes, reported separately, never summed

| route | target list | what a number from it means |
|---|---|---|
| **EC** | `target_clusters.tsv` — 17,072 clusters from HUMAnN's map across the ten ECs | community-level activity, directly comparable with published HUMAnN EC output |
| **curated** | `curated_clusters.tsv` — 224 clusters containing characterised Swiss-Prot sequences | abundance of genes in clusters holding a *characterised* enzyme |

Summing them would rebuild the base-rate flaw §4.6 identifies in the AGORA2 producer
call, and would also mix "annotated as" with "characterised as". A cluster mapped to
EC 4.1.1.18 is annotated as a lysine decarboxylase, which inherits UniProt's automatic
TrEMBL annotation quality; a cluster in the curated list contains a sequence a curator
has read the literature for. Both are useful. They are not the same claim.

This is the same availability-versus-support separation §4.7 claims as a methodological
contribution, applied to gene families.

## What came out

The extraction ran on 2026-11-08 (job 27319804, 5m09s, 14 cores, MaxRSS 5.5 GiB). The
object asserts its own name — `2021-10-14.YachidaS_2019.gene_families` — and is a
`dgCMatrix`, 8,930,907 × 616, 379,029,747 nonzeros at 6.89% density. **2,371,507 rows
(26.6%) are species-stratified**, not the 5,351,296 quoted in
`writeup/read_programme_alternative.md`; that figure should be corrected there.

3,516 rows matched, over **1,267 of the 17,100 target clusters**, 1,164 of them
species-stratified.

### The tryptamine mechanism is now attributed, and only the curated route found it

```
UniRef90_A7B1V0                                      nonzero in 379/616 samples
UniRef90_A7B1V0|g__Blautia.s__Ruminococcus_gnavus    nonzero in 378/616 samples
UniRef90_A7B1V0|unclassified                         nonzero in   8/616 samples
```

The *R. gnavus* tryptophan decarboxylase that Results §3.7 cites as the mechanistic
agreement for tryptamine is present in **378 of 616 samples and attributed to
*R. gnavus* itself**. That is the "can this organism perform the reaction" evidence
§4.6 says is missing, as a per-sample per-species gene abundance rather than an
inference from a reference genome.

The EC route would have returned nothing for it. This row exists in the output only
because `curated_clusters.tsv` asked for the cluster id directly.

`UniRef90_J7SZ64`, the *C. sporogenes* enzyme, is **not present in this cohort** — a
real negative, since the probe demonstrably detects the homologous cluster.

### Species attribution per activity

| EC | stratified rows | species | most prevalent carrier |
|---|---|---|---|
| 2.3.1.1 | 369 | 283 | *F. prausnitzii* (535/616) |
| 2.3.1.35 | 261 | 230 | *F. prausnitzii* (535/616) |
| 4.1.1.19 | 333 | 222 | *B. uniformis* (549/616) |
| 4.1.1.15 | 148 | 103 | *B. uniformis* (560/616) |
| 4.1.1.17 | 156 | 55 | *H. parainfluenzae* (264/616) |
| 4.1.1.18 | 142 | 46 | *P. succinatutens* (237/616) |
| 4.1.1.22 | 12 | 10 | *R. ornithinolytica* (30/616) |
| 4.1.1.25 | 11 | 10 | *M. smithii* (31/616) |

Note how the base rate differs by two orders of magnitude between the N-acetylglutamate
routes (283 species) and the aromatic decarboxylases (10). The §4.6 base-rate problem
applies here too and any producer call must be read against the column of species counts,
not in isolation.

### The *Allisonella* question: not answerable by this route, and we can prove why

Results §3.7 ranks *Allisonella* first for cadaverine while noting the genus is
characterised by histidine decarboxylation. The extraction found
*A. histaminiformans* profiled with **1,487 stratified gene families**, and **none** of
them in any of the 416 target lysine decarboxylase clusters. That looks like an answer.
It is not, and the reason is the control.

`resolve_allisonella.py` resolved all 1,487 cluster representatives at UniProtKB
(they are filed under *Dialister histaminiformans*, the same organism after the genus
was reassigned). It then asked the positive-control question: a **histidine**
decarboxylase must be among them, because the organism is named for that activity.

There is none. Nor any lysine decarboxylase. The 11 decarboxylases present are
ornithine, phosphatidylserine, uroporphyrinogen, oxaloacetate, diaminopimelate,
S-adenosylmethionine and orotidine-5'-phosphate; the six "histidine" names are tRNA
ligases, histidine kinases and a HIT-family protein.

So the probe cannot see the one enzyme this organism is defined by, and the diagnosis
goes one level deeper: **UniProt holds 500+ entries for *histaminiformans* and not one is
annotated as a histidine decarboxylase** — not by EC 4.1.1.22, not by protein name, not
by gene `hdcA`. No UniRef90 cluster for that enzyme exists, so HUMAnN could never report
one, and the control could not have passed by construction.

**The absence of a lysine decarboxylase in *Allisonella* is therefore uninterpretable.**
Reporting it as a finding would be reporting a gap in UniProt as a fact about the gut.
§3.7's caveat stands as a caveat — this evidence neither supports nor refutes it.
Answering it needs assembly and binning of *Allisonella* to recover its actual gene
complement (§2 of `writeup/dardel_plan.md`, Paper 2), or a targeted experiment.

This is the one result here that is a negative about the method rather than about the
data, and it is the reason the control was built in rather than bolted on.

## Files

| | |
|---|---|
| `fetch_mappings.sh` | pulls the three mapping files from Zenodo with md5 verification |
| `extract_targets.py` | ten ECs → 17,072 UniRef90 clusters. Imports the EC list from `../reads/build_enzyme_db.py` rather than restating it |
| `validate_mapping.py` | the table above. Resolves every curated accession's UniRef90 cluster and checks it against the map |
| `extract_gene_families.R` | the one job that loads the 3.26 GB object, with the self-naming assertion |
| `extract_gene_families.sbatch` | 1 core, 32 GB, ~2 h ceiling |
| `target_clusters.tsv`, `curated_clusters.tsv` | the two target lists |
| `uniref90_cache.tsv` | accession → UniRef90 cluster, cached so 363 round trips happen once |
| `mapping_validation.tsv` | per-accession verdicts behind the table |
| `resolve_allisonella.py` | the *Allisonella* question and its positive control |
| `allisonella_resolved.tsv` | all 1,487 cluster representatives, named and EC-annotated |
| `out/species_attribution.tsv` | which species carries which activity, per sample count |
| `out/gf_targets_long.tsv.gz` | 207,026 nonzero (row, sample, abundance) triples — the join key for the metabolite data |
| `out/gf_targets_summary.tsv` | per-row prevalence and maximum, with route and EC labels |
| `out/allisonella_genes.tsv` | the full 1,487-gene complement |
| `out/rownames_index.txt.gz` | all 8.9 M row names. **Gitignored** at 42 MB; regenerate with the batch job. Keeping it means no later question needs a second five-minute load |

## Provenance of the mappings

From **Zenodo record 17100034**, `doi:10.5281/zenodo.17100034` ("ChocoPhlAn mapping
files"), md5-verified on fetch:

| file | bytes | md5 |
|---|---|---|
| `map_level4ec_uniref90.txt.gz` | 27,591,865 | `3ace48af4aed6a2334959ed64d1f76e8` |
| `map_ko_uniref90.txt.gz` | 17,547,123 | `44c7c7b6ce27bd4600f2e99de3c0aa43` |
| `map_ec_name.txt.gz` | 75,265 | `d324852a5092d38bd3a7035f3506e8f1` |

Not from HUMAnN's own downloader, which takes the whole `full_mapping_v201901b.tar.gz`
bundle — **2,734,293,864 B**, measured by HEAD on 2026-11-08 — to deliver the 26 MB file
we need, carries no checksum, and is served over a redirect chain ending at a Globus
endpoint. Zenodo serves the files individually under a citable DOI. 45 MB instead of
2.73 GB.

One path trap, recorded because several tutorials get it wrong:
`…/humann_data/full_mapping/<file>` is a 404, while `…/humann_data/full_mapping_v201901b.tar.gz`
with no subdirectory is a 200.

## Two measurements that sized the work

**The object is not what its shape implies.** 8,930,907 × 616 dense doubles would be
41.0 GiB. The measured uncompressed serialization is **4.58 GiB**, 0.112 of that, so the
`.rda` is not a dense double matrix of that shape. One core and 0.2 core-hours of
streaming `bzcat | wc -c` established it, which is why the job asks for 32 GB rather than
roughly 90 cores' worth of memory in `shared`.

**The `.rda` files are bzip2, not gzip** (`42 5a 68` = `BZh`). `zcat` fails on them with
"not in gzip format"; use `bzcat`, or just `load()`.

## What this can and cannot establish

It measures **DNA**. A stratified gene family says an organism carrying a gene annotated
for the reaction is present, and how abundant that gene is. It is not expression and not
flux, so it is not evidence the reaction ran. That is exactly the "can this organism
perform the reaction" source §4.6 says is missing — no more.

Read-based and profile-based results are also **single-cohort by construction**: only
Yachida has shotgun data. Kim and Sinha are genus-level amplicon cohorts, so nothing here
replicates across cohorts.

The *Allisonella* question has a third possible answer that must not be collapsed into
the other two. Results §3.7 ranks the genus first for cadaverine while noting it is
characterised by histidine decarboxylation. The stratified rows can show a lysine
decarboxylase under `g__Allisonella`, or show none — but "none" is only informative if
*Allisonella* is stratified in this matrix at all, so the extraction reports that count
separately rather than reading an absent genus as an absent gene.
