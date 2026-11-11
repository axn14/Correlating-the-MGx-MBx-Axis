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
`dgCMatrix`, 8,930,907 × 616, 379,029,747 nonzeros at 6.89% density.

The row counts need stating precisely, because two different true numbers get called
"the stratified rows" and an earlier draft of this file called one of them wrong:

| row kind | count | share |
|---|---|---|
| community total, `UniRef90_X` | 3,579,611 | 40.1% |
| species-stratified, `UniRef90_X\|g__G.s__S` | 2,371,507 | 26.6% |
| `UniRef90_X\|unclassified` | 2,979,789 | 33.4% |
| **rows containing a pipe** | **5,351,296** | **59.9%** |

So the 5,351,296 in `writeup/read_programme_alternative.md` is correct for "rows below
the community level", and 2,371,507 is correct for "rows naming a species". They are not
interchangeable: anything that treats *not community-level* as *attributed to a species*
double-counts, because the community total is the sum of the species rows **plus** the
unclassified row. `progression_axis.py` splits on the row name for exactly this reason.
Measured by `colsums.R`, job 27324295.

Those column sums also settle the units: the community-level rows sum to **1.000000 in
every one of the 616 samples**, so these are relative abundances and no renormalization
applies to any comparison below. Between 61% and 93% of each sample's gene pool
(median 83%) is attributable to a named species.

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

## Reformulated, 2026-11-11: `compositional.py` (supersedes the "source replacement" wording)

A methodological review raised two problems with the section below, and both were correct.
First, the gene families are **relative**, so a flat EC total cannot show that capacity
is "conserved" in absolute terms. Second, a gene is collinear with its carrier, so
partialling the carrier out is not a clean test. `compositional.py` re-asks both
questions and adds log10 sequencing depth (`v2/targets/depth_screen.py`) to the confounders.

**1. EC totals as an additive log-ratio to a fixed reference** (5,839 prevalence-matched
background families, `gf_background_long`). This is scale-invariant, unlike a share.
**No EC trends with stage** (9 tests, all q > 0.6). Arginine decarboxylase, EC 4.1.1.19:
ρ −0.015, p 0.80. The flat result holds in compositional form. It is still a statement
about **relative metagenomic space**: Yachida has no load measurement, so absolute
capacity per gram of stool is not estimable.

**2. log(gene / clade): the dose per carrier genome**, on samples where both are detected.
**None of 279 testable gene × species rows changes with stage** (all q > 0.44).
- **The *Butyricimonas virosa* arginine-decarboxylase result does not survive.** Its gene
  share trends (ρ 0.15, p 0.013), but per genome it is flat (ρ −0.07, p 0.50). The earlier
  ρ = 0.22 "carrier-adjusted" figure came from the partial-correlation formulation and
  **must not be cited**.
- *Parvimonas micra* glutamate decarboxylase: the share trends strongly (ρ 0.28), the
  ratio is flat (ρ 0.004).

**What is defensible is taxonomic turnover within a functional group, in relative
metagenomic space.** The organisms carrying these enzymes change with stage, while each
organism's own dose of the gene does not. "Source replacement of a conserved capacity" is
withdrawn as wording.

**3. Class G** (`out/class_g.tsv`, CONTRACTS §5d), for NB09. A (GTDB species, target) pair
is supported when reads of a curated *producing* enzyme for the target are attributed to
the species in at least 10% and at least 10 of the Yachida samples.
- Taxonomy matching: HUMAnN (ChocoPhlAn) names are matched to GTDB after stripping GTDB's
  suffixes from both genus and epithet. `_match_gtdb` missed genus suffixes, so it never
  matched *Ruminococcus_B gnavus*.
- 158 GTDB species match uniquely. 352 are GTDB splits of one HUMAnN species; those are
  **unavailable**, not credited to every daughter species.
- 115 (species, target) pairs are supported.
- **Tryptamine:** *R. gnavus* is supported (UniRef90_A7B1V0 in 214 of the Yachida
  samples). *E. ramosum* is assessable and has **no** tryptophan-decarboxylase reads
  attributed to it.

## On the progression axis: `progression_axis.py`

> **Superseded in part (2026-11-11):** see the section above. The *B. virosa* "survives carrier control" result and the "capacity conserved / source replaced" wording below are withdrawn.


The extraction above says whether a gene is there and whose it is. It says nothing about
disease, because the gene matrix carries no clinical labels. `progression_axis.py` joins
the two and asks whether biosynthetic capacity moves along the CRC axis.

**The join needs its own file.** The gene matrix is keyed by DDBJ BioSample accession
(`SAMD00114718`); every cohort table in this project is keyed by Yachida's own sample
number (`10021`). The bridge is curatedMetagenomicData's curation table, where
`subject_id` is literally `sub_` + the cohort number — committed here as
`ref/YachidaS_2019_sample.tsv` rather than fetched at run time, because a silent change
upstream would silently re-label every sample in every figure. All 347 cohort samples
resolve, all 347 accessions are columns of the matrix, and the 616 subject ids are
unique, so the mapping is a function. Each of those is asserted; any failure raises.

Groups come from `v2/utils.py` unchanged: primary Healthy 127 → Early_CRC 96 →
Advanced_CRC 54 (n=277), early axis adding Adenoma 40 (n=317), PostSurgery 30 as a
sensitivity group that is never pooled. Each axis is its own BH family.

### Nothing separates two stages

| level | features tested | q < 0.05 | smallest q |
|---|---|---|---|
| per cluster, community level | 994 | **0** | 0.185 |
| per species-stratified row | 748 | **0** | 0.204 |
| per-EC aggregate, both routes | 39 | **0** | 0.825 |

The per-EC aggregates are flat on every contrast and on the trend. Whatever happens
along this axis, it is not a change in the total gene pool devoted to these activities.

### What does move is which copy is present

Spearman against the ordinal scale finds **12 of 455 features at BH q < 0.05**, and
**4 survive adjustment for age, BMI, Brinkman index, alcohol and sex** (all five complete
for all 347 samples, so nothing is imputed). Within a single EC they move in **opposite
directions**, split by the organism carrying them:

| cluster | EC | carrier | ρ | q | ρ adj | ρ carrier-adj |
|---|---|---|---|---|---|---|
| `UniRef90_A8SL35` | 4.1.1.15 glutamate dc | *P. micra* | +0.31 | 3.8e-05 | +0.30 | +0.11 (p 0.08) |
| `UniRef90_Q8RG18` | 4.1.1.17 ornithine dc | *F. nucleatum* | +0.25 | 7.9e-03 | +0.26 | +0.08 (p 0.21) |
| `UniRef90_R5HT10` | 2.3.1.1/35 ArgA | unattributed | +0.24 | 9.6e-03 | +0.23 | — |
| `UniRef90_A0A174QH87` | 4.1.1.19 arginine dc | unattributed | +0.23 | 1.1e-02 | +0.21 | — |
| `UniRef90_R5P0Y5` | 4.1.1.19 arginine dc | *B. virosa* | +0.22 | 1.5e-02 | +0.21 | **+0.22 (q 0.003)** |
| `UniRef90_E5Y860` | 4.1.1.19 arginine dc | *B. wadsworthia* | +0.20 | 2.8e-02 | +0.18 | +0.09 (ns) |
| `UniRef90_A0A174M373` | 4.1.1.19 arginine dc | *L. pectinoschiza* | −0.20 | 4.0e-02 | −0.19 | **−0.16 (q 0.039)** |
| `UniRef90_A0A0M6WRA6` | 2.3.1.1/35 ArgA | *R. faecis* | −0.21 | 2.6e-02 | −0.17 | — |
| `UniRef90_A0A173T3V0` | 4.1.1.19 arginine dc | *R. faecis* | −0.21 | 2.6e-02 | −0.17 | — |
| `UniRef90_A0A174GRE2` | 4.1.1.19 arginine dc | *L. pectinoschiza* | −0.21 | 1.8e-02 | −0.19 | **−0.14 (q 0.048)** |
| `UniRef90_R6K290` | 2.3.1.1/35 ArgA | *L. pectinoschiza* | −0.22 | 1.7e-02 | −0.20 | **−0.14 (q 0.048)** |
| `UniRef90_D4L1D2` | 4.1.1.18 lysine dc | *R. intestinalis* | −0.22 | 1.5e-02 | −0.20 | −0.10 (ns) |

Six clusters of arginine decarboxylase trend, three up and three down. The EC 4.1.1.19
aggregate over all 295 observed clusters does not move at all (ρ −0.03, q 0.91). The
same holds for ArgA. **The capacity is conserved while its source is replaced** — commensal
Lachnospiraceae copies fall, *Fusobacterium*, *Parvimonas*, *Bilophila* and
*Butyricimonas* copies rise. `figs/fig1_source_turnover.png` is that result in one
picture.

### Two controls, and both of them bite

**Is trending with stage rare?** `background.R` draws five non-target community-level
gene families per observed target, matched decile by decile on prevalence — the thing a
rank test's power actually depends on. The matching holds: medians 65 and 67 detected of
277, Mann–Whitney p = 0.39. Comparison is on raw p, because BH's q depends on family
size and the two families differ fivefold by construction.

| threshold | target (n=455) | background (n=2,107) | ratio | Fisher p |
|---|---|---|---|---|
| p < 0.05 | 7.3% | 10.9% | **0.66×** | 0.021 |
| p < 0.01 | 3.1% | 2.5% | 1.25× | 0.42 |
| p < 0.001 | 1.3% | 0.38% | **3.47×** | 0.025 |
| p < 1e-4 | 0.44% | 0.09% | 4.63× | 0.15 |
| smallest p | **3.1e-07** | 3.6e-05 | | |

(confounder-adjusted; unadjusted is the same shape, 8.6% against 13.3% at p < 0.05.)

Read at one threshold this says opposite things, so it is reported at four. **These
enzymes are, as a set, LESS stage-responsive than random gut gene families** — which is
what you would expect of redundant core metabolism, and which rules out any sentence
claiming that amine biosynthesis is broadly engaged along the axis. The enrichment is in
the tail alone: a handful of specific copies move harder than anything in a five-times
larger background.

**Is the gene moving, or its organism arriving?** A gene's share of the gene pool is
roughly its carrier's abundance times its share within that carrier, so a rising gene in
a rising organism is not independent evidence. `progression_axis.py` adds the carrier's
own relative abundance — from the cohort's GTDB species table, an independent
measurement, not from this gene matrix — as a sixth covariate.

The two strongest findings collapse. *P. micra*'s glutamate decarboxylase goes from
ρ 0.31 to 0.11 (p 0.08); *F. nucleatum*'s ornithine decarboxylase from 0.25 to 0.08
(p 0.21). Both organisms themselves track stage at p = 2.6e-11 and 3.1e-10. Four trends
survive, led by `UniRef90_R5P0Y5` (*B. virosa* arginine decarboxylase, ρ 0.22 → 0.22,
q 0.003) **whose carrier does not itself trend** (ρ 0.06, p 0.31).

So the defensible claim is **attribution, not independent gene-level change**: these
organisms arrive carrying these specific decarboxylase activities. That is the
"can this organism perform the reaction" evidence §4.6 asks for. It is weaker than
"capacity increases", which these data do not support.

One carrier, *R. faecis*, has no GTDB counterpart in the cohort's species table (GTDB
reclassified it), so two clusters are marked untestable rather than tested against the
wrong organism. *F. nucleatum* matches four GTDB columns, which are summed.

### The chapter's tryptamine claim, answered both ways

`UniRef90_A7B1V0` is **flat**: ρ +0.013 (p 0.83), adjusted +0.008, no pairwise contrast
near significance. Prevalence runs 78/127 Healthy, 65/96 Early, 30/54 Advanced. The
presence-and-attribution result above stands; the stage-dependence does not exist.
`figs/fig4_tryptamine.png` states both on the figure, because a reader shown the
attribution panels alone will assume a gradient that is not in the data.

### Figures

| | |
|---|---|
| `figs/fig1_source_turnover` | one activity, two directions, two sets of organisms; the aggregate underneath, flat |
| `figs/fig2_trend_forest` | every q<0.05 cluster with its adjusted twin, and the background control |
| `figs/fig3_stage_profiles` | mean share by stage for all twelve, off-axis groups drawn unconnected |
| `figs/fig4_tryptamine` | the chapter's lead claim: present, attributed, not stage-dependent |
| `figs/fig5_enzyme_log2fc` | per-EC aggregate, both routes, blind-spot rows hatched |
| `figs/fig6_cluster_volcano` | the pairwise negative, drawn so it is not assumed unrun |

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
| `extract_gene_families.sbatch` | 1 core, 12 GB, 1 h ceiling — sized from the measured 5.5 GiB peak, not from the matrix shape |
| `colsums.R` | per-sample community/stratified column sums of the full matrix. The job that settled the units and the row-kind counts (27324295) |
| `background.R` | the prevalence-matched background set: 6,335 non-target clusters, 5 per observed target, matched decile by decile (27324874) |
| `progression_axis.py` | the join to the clinical axis, every test, both controls, all six figures |
| `ref/YachidaS_2019_sample.tsv` | curatedMetagenomicData's curation table — the only thing connecting `SAMD…` accessions to Yachida's sample numbers. Committed so the join cannot change under us |
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
| `out/sample_normalization.tsv` | per-sample community and stratified column sums; the evidence that these are relative abundances |
| `out/sample_crosswalk.tsv` | the asserted join: cohort sample ↔ accession ↔ group |
| `out/progression_cluster_tests.tsv` | per cluster × contrast × family, with the untested ones kept and the reason recorded |
| `out/progression_enzyme_tests.tsv`, `out/progression_species_tests.tsv` | the same at the aggregate and species-row levels |
| `out/progression_trend.tsv` | Spearman unadjusted and adjusted, both axes, with carrier species and ECs |
| `out/background_trend.tsv` | the identical trend on the matched background set |
| `out/progression_carrier_control.tsv` | per hit: trend after holding the carrier organism's abundance fixed |
| `out/gf_background_long.tsv.gz` | the background set's abundances, 442,566 triples over 6,335 clusters |

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
replicates across cohorts. Every progression-axis number above is therefore a discovery
result with no external validation available, and should be written up as one.

Four further limits on the progression analysis specifically:

- **The pairwise contrasts are null and the trends are not a substitute for them.** A
  monotone rank correlation over 277 samples is a different and more powerful test than
  Healthy-vs-Advanced over 181, but it also assumes the three stages are equally spaced
  on the scale `STAGE_GROUP_NUMERIC` gives them. That coding is a modelling choice, not
  a measurement.
- **The carrier control is only as good as the carrier assignment.** "Carrier" here means
  the species-stratified row of that cluster with the highest prevalence, which for a
  cluster spread across many organisms is a simplification. Two clusters have no
  stratified row at all and two have a carrier GTDB does not recognise, so eight of the
  twelve hits are testable and four are not.
- **`unattributed` is not `absent`.** Two of the twelve trending clusters have no species
  row, meaning HUMAnN could not assign their reads to a pangenome. They are shown and
  labelled, not dropped, but nothing here says which organism they belong to.
- **The background control compares gene families, not pathways.** It establishes that
  these clusters are not unusually stage-responsive as a set. It does not establish that
  no amine-biosynthesis pathway responds — a pathway-level test over a matched pathway
  background is a different experiment and has not been run.

The *Allisonella* question has a third possible answer that must not be collapsed into
the other two. Results §3.7 ranks the genus first for cadaverine while noting it is
characterised by histidine decarboxylation. The stratified rows can show a lysine
decarboxylase under `g__Allisonella`, or show none — but "none" is only informative if
*Allisonella* is stratified in this matrix at all, so the extraction reports that count
separately rather than reading an absent genus as an absent gene.
