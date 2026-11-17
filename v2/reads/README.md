# `v2/reads` — the gene-level evidence the chapter is missing

Discussion §4.6 names three intended evidence sources that were never delivered, and
says plainly why it matters: *"of the six declared sources, the two that speak most
directly to whether an organism **can** perform a reaction were either uninformative by
base rate or unavailable."* That is why the eight three-source pairs are prioritised
hypotheses rather than established producer relationships. Everything here exists to
change that.

Nothing in this directory moves bulk data. The reads are fetched, profiled and deleted
per sample; see **Stream and discard** below.

## Status — read this before running anything here

**The FASTQ route is superseded.** `writeup/read_programme_alternative.md` establishes
that `curatedMetagenomicData` already distributes bioBakery 3 HUMAnN output for
`YachidaS_2019`, and the profiles are in project storage at `../../profiles/`
(3.1 GB). Verified independently on 2026-11-08: both small matrices load naming
themselves `2021-10-14.YachidaS_2019.*`, and **347 of 347** of our samples are columns
in each. `pathway_abundance` is 31,291 × 616 with 96.9% of rows species-stratified;
`gene_families` is 8.9 M UniRef90 rows, of which 5.35 M carry a
`|g__Genus.s__Species` suffix.

That is strictly more than a read programme would have produced — species-resolved
rather than needing a separate attribution step, for all 347 samples rather than a
pilot's 30, at zero HUMAnN core-hours instead of ~2,500 — so **do not download reads to
get gene abundances.**

| | |
|---|---|
| `manifest_347.tsv` | FASTQ URLs, sizes and MD5s. **347 of 347**, one WGS run each, 1.713 TB. Kept for the designed-subsample fallback, not for a full download |
| `pilot_30.tsv` | the 30-sample pilot selection. **Shelved**, superseded as above. Kept because the rule is declared and reproducible if a read pilot is ever wanted |
| `enzymes/` | the two-arm probe, 363 Swiss-Prot / 32,348 with TrEMBL. **Role changed**: now the validation reference for the UniRef90→EC mapping, and the fallback probe if the mapped ECs prove too thin |
| profiling | not run here; the profiles were obtained rather than computed |

The fallback, if UniRef90→EC mapping cannot carry the specific claims: a **designed
subsample** of 40 samples from each tail of the relevant metabolite distribution,
stage-balanced — roughly 400 GB, 4–11 h at the measured 4.2 MB/s. Extreme-phenotype
sampling is the efficient design for a gene–metabolite association anyway. That is what
`manifest_347.tsv` and `enzymes/` are retained for.

## The sample selection rule, declared before use

Six samples from each of the five groups on the progression axis — Healthy, MP
(adenoma), Stage_0, Stage_I_II, Stage_III_IV — taken at **evenly spaced quantiles of
FASTQ volume within each group**. HS (post-surgery) is excluded; it is not a point on
the progression axis and the chapter treats it separately.

Quantiles rather than a random draw, for two reasons. HUMAnN's runtime scales with read
depth, and the pilot's second job is to measure cost: sampling the depth range returns a
distribution that extrapolates to 347, where six random draws may cluster near the
median and understate the tail. The realised spread is 1.25–9.96 GB, an 8× range. And
quantiles need no seed, so there is no hidden choice in which samples these are.

The rule looks at group label and file size only. **No metabolite and no species
abundance enters it.** Choosing samples on cadaverine to make the *Allisonella* question
easier to answer would guarantee the answer.

## Two arms, reported separately, never summed

`build_enzyme_db.py` builds two reference sets for the same ten activities:

- **specific** (`enzymes_sprot.faa`, 363 seqs) — Swiss-Prot only. A curator has read the
  literature for every entry, so a hit means a read aligned to a *characterised*
  enzyme. This arm supports presence.
- **sensitivity** (`enzymes_trembl.faa`, 32,348 seqs) — adds TrEMBL. Annotations are
  automatic and propagated. This arm says whether an absence in the specific arm is a
  real absence or a thin reference set.

Merging them would rebuild the exact flaw §4.6 identifies in the AGORA2 producer call:
nearly every modelled species "encodes the enzyme", the base rate goes too high to
discriminate, and the pre-declared 0.20 gate rejects it for six of seven compounds. The
separation is the same availability-versus-support discipline §4.7 claims as a
methodological contribution, applied to reads.

The *Allisonella* question is what makes two arms necessary rather than tidy. Results
§3.7 wants to say something about whether the genus carries a lysine decarboxylase.
Swiss-Prot has **8** bacterial sequences for EC 4.1.1.18; finding nothing against 8
sequences is evidence about the reference set, not about *Allisonella*. TrEMBL has 1007.

| EC | specific | sensitivity | activity → claim |
|---|---|---|---|
| 4.1.1.105 | 2 | 30 | tryptophan dc → tryptamine, phenethylamine |
| 4.1.1.28 | 1 | 163 | aromatic-L-amino-acid dc |
| 4.1.1.25 | 4 | 414 | tyrosine dc → tyramine |
| 4.1.1.18 | 8 | 1007 | lysine dc → **cadaverine, the Allisonella question** |
| 4.1.1.22 | 19 | 424 | histidine dc → histamine |
| 4.1.1.15 | 16 | 3368 | glutamate dc → GABA |
| 2.3.1.1 | 194 | 18810 | ArgA, N-acetylglutamate de novo |
| 2.3.1.35 | 130 | 12808 | ArgJ, N-acetylglutamate by recycling |
| 4.1.1.17 | 5 | 2576 | ornithine dc → putrescine |
| 4.1.1.19 | 114 | 5678 | arginine dc → agmatine |

Note the asymmetry: 324 of the specific arm's 363 sequences serve N-acetylglutamate,
while the two claims the chapter leads with get the thinnest probes. That is a property
of what has been curated, not a choice, and it belongs in any write-up of these results.

**ArgA and ArgJ are kept apart.** N-acetylglutamate is formed de novo by ArgA
(EC 2.3.1.1) or by recycling the acetyl group from N-acetylornithine by ArgJ
(EC 2.3.1.35). Both bear on §4.4's N-acetylglutamate finding, they are different
mechanisms, and one number covering both would misreport which the data support.

## Membership is decided by annotation, not by name

Every sequence carries a target EC in UniProt's own `ec` field, or is one of two named
loci admitted deliberately. The first version of the script matched on strings and 81%
of its output was wrong, because `cadA` is the gene symbol for inducible lysine
decarboxylase **and** for the cadmium-translocating ATPase. At TrEMBL scale that one
collision contributes **5,191** cadmium pumps — the lysine decarboxylase probe for the
cadaverine question would have been roughly 85% cation transporters. `gadB` similarly
pulls in GadX regulators, and a cadaverine/lysine *antiporter* arrives too: a
transporter, not an enzyme.

None of that fails loudly. It quietly counts reads against the wrong protein. The EC
gate makes it structurally impossible; `rejected_*.tsv` records what was turned away.

Two loci are admitted by name because the chapter's mechanistic argument names them,
and both are verified present on every build:

- `RUMGNA_01526` → **A7B1V0**, the *R. gnavus* tryptophan decarboxylase that Results
  §3.7 cites as the mechanistic agreement for tryptamine, also active on phenylalanine
  — which is why *R. gnavus* ranks second for phenethylamine.
- `CLOSPO_02083` → **J7SZ64**, the *C. sporogenes* enzyme of the same family.

Both turn out to be the only two reviewed bacterial EC 4.1.1.105 entries in existence.

## Stream and discard

Per sample: fetch → `fastp` QC → host removal → MetaPhlAn + HUMAnN + targeted DIAMOND
→ keep the profile tables, which are a few MB → **delete the FASTQ**. Peak footprint is
whatever is in flight, about 10–25 GB, in scratch, which is not migrated and purges
itself.

This is not an optimisation. Klemming is shared, 86% full, and being migrated to
Arrhenius under a quota the project answers for; the FASTQs are public and
re-downloadable from ENA at any time, so storing them would spend part of the project's
single migration request on data that needs no preserving. Downloading 1.7 TB before
processing anything would also have been 630× the project footprint this work actually
occupies (2.7 GB).

Bandwidth does not go away: ENA gives ~4.2 MB/s per stream from the login node, and PDC
names `dardel-ftn01.pdc.kth.se` for bulk transfers, with login nodes for small ones.
Compute nodes are not documented for transfers and measured 366 KB/s, an order of
magnitude worse.

## What this can and cannot establish

A read aligning to a characterised lysine decarboxylase says lysine decarboxylase is in
the sample. **It does not say whose.** Attributing a gene to an organism needs assembly,
binning and a dereplicated catalogue — §2 of `writeup/dardel_plan.md`, 20,000–30,000
core-hours, and the item the plan flags as most likely to overrun. The targeted search
is the confirmatory check on the eight pairs, not the attribution itself.

Read-based results are also **single-cohort by construction**: only Yachida has shotgun
data. Kim and Sinha are genus-level amplicon cohorts, so nothing here replicates across
cohorts, and that limitation belongs wherever these results are reported.
