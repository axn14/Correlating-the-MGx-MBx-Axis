# `v2/replication` — what can and cannot replicate a species claim (CONTRACTS §5d)

## Two kinds of replication, kept apart

| | Cohorts | Resolution | Role in the tiers |
|---|---|---|---|
| **Ecological consistency** | Sinha 2016, Kim 2020 | 16S, genus level; other metabolomics platforms | Reported only; **never votes** for a species claim (`evidence_tools.NON_VOTING_SOURCES`). |
| **Species-level replication (class X)** | Erawijantari 2020 | shotgun species; CE-TOFMS from the same group and platform as Yachida | Votes as class X. |

Metabolites are matched across cohorts by identifier, not by name. NB02 already does this
through `utils.metabolite_identity` (KEGG *or* HMDB).

## Erawijantari: what was usable

- **53 of its 54 healthy subjects are Yachida subjects, and they are the same samples.**
  All 53 metabolome profiles are identical, and the species profiles agree at median
  r = 1.00. These samples were excluded, and so was the one remaining healthy subject.
- **The replication sample is 42 post-gastrectomy patients.** It is independent, but small
  and surgically altered. The smallest detectable |ρ| is 0.33 before FDR correction. A
  failure to replicate is weak evidence; a success is real.

## Result (run `20261108T133844Z_1720a81`, 2026-11-11)

`python v2/replication/erawijantari_species.py --results "<pulled run>"`

| | |
|---|---|
| Yachida target pairs (NB08, the 49 targets) | 4,812 |
| measurable in Erawijantari | 4,037 (457 species × 42 metabolites) |
| **replicated** (same sign, BH q < 0.05) | **256** |
| contradicted (opposite sign, q < 0.05) | 0 |
| sign agreement | **75%** |

**Null check.** The metabolome was shuffled across samples 5 times. That gave 0 replicated
pairs and 46–59% sign agreement. The signal is real, not an artefact of the pipeline.

**By target.** Most replications are for glutamate (96), putrescine (72), N-acetylglutamate
(43) and N-methylglutamate (21). Tryptamine is not measured in Erawijantari, so the
*E. ramosum* and *R. gnavus* pairs are *not available* rather than failed. Two of the five
N-acetylglutamate pairs the old tiers called Strong replicate: CAG-83 sp000435555
(ρ 0.47 → 0.56) and *Ruminiclostridium_E* sp002437415 (0.42 → 0.55).

Columns written for NB09: `Source_replication_species_available`,
`Source_replication_species_support` and `Replication_species_contradicted`.
