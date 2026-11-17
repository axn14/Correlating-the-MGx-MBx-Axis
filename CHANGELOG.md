# Changelog

## v2.0 (November 2026)

v2.0 is a rewrite of the analysis, not an update to it. It changes how samples are grouped,
how targets are chosen and what counts as evidence. Several results reported with
`v1.0-thesis` do not survive these changes. They are listed here so that nobody relies on
them.

### Corrections to v1.0-thesis

**Sample groups.** v1 grouped Yachida's samples as Early CRC = stage 0 + HS (n = 57) and
Advanced CRC = stage I–II + stage III–IV + MP (n = 163). HS are patients sampled after their
cancer was surgically removed, and MP are people with multiple polypoid adenomas. Neither
belongs in a cancer stage. The v1 wiki also described HS as "high-grade adenoma" and MP as
"metachronous polyp", which is incorrect. v2 uses Early CRC = stage 0 + I–II (n = 96) and
Advanced CRC = stage III–IV (n = 54). HS is analysed as a separate sensitivity arm and MP as
an adenoma arm. Much of v1's "early CRC" signal (for example 38 differentially abundant
species) compared healthy guts with resected ones. Under the corrected groups, 0 species are
robustly different between healthy people and early CRC.

**Differential abundance.** v1 accepted features at uncorrected p < 0.05 in some gates. v2
requires q < 0.05 in both a rank test and a covariate-adjusted compositional model, with
sequencing depth among the covariates.

**Target selection.** v1 chose targets partly by how well species predicted them in the same
cohort, so later tests of whether species account for those targets were partly answered in
advance. v2 admits a metabolite only if it is robustly differential or belongs to an enzyme
class fixed from the literature, and a check stops the pipeline if target selection reads any
species–metabolite statistic.

**Evidence counting.** v1's nine-stream score and "Holy Trinity" counted several readings of
the same association (SHAP, MOFA, mediation, within-stage correlation) as separate evidence,
used genus-level community-model results for species-level claims, and allowed tumour location
into models where it encodes the diagnosis. Its mediation result (30 of 30 forward paths
significant) followed from how the pairs were selected. v2 groups evidence into five
independence classes and counts each class once, never lets a genus vote for a species, and
records tumour location only within cancers.

**Calibration.** v1 did not ask how often its tiers would be met by chance. v2 measures it
with a within-metabolite permutation null and a replication power analysis.

### Effect on reported results

| v1 claim | v2 status |
|---|---|
| Three species-level Trinity pairs: *F. animalis* and phenylalanine, *B. wadsworthia* and isoleucine, *S. wadsworthensis* and leucine | *F. animalis* and phenylalanine: association present (ρ 0.34) but no gene evidence and no replication; a hypothesis only. The other two associations are not significant. |
| Ten genus-level Trinity associations | Not supported; a genus cannot stand in for a species. |
| *Escherichia* and lactate as a validated driver (thesis abstract) | Already withdrawn by the May 2026 review (0 Trinity pairs); not supported. |
| Early CRC dysbiosis (38 species) | Artefact of pooling post-surgery patients; 0 robust species. |
| *Veillonella* and putrescine | Holds. *V. atypica* and putrescine is the one link that meets every v2 criterion. |
| *R. gnavus* and tryptamine | Holds, with patient-level gene evidence for the characterised enzyme. |

### New in v2.0

The changes above led to these additions:
- species-level replication in Erawijantari et al. 2020, with the 53 subjects shared with
  Yachida removed;
- gene evidence read from the patients' own metagenomes;
- community models solved to a certified optimality gap without a commercial solver;
- export, origin and confounder gates per target;
- the calibration analyses;
- a dated rule book (`v2/CONTRACTS.md`) in which each amendment was committed before the
  results it governs.

## v1.0-thesis (May 2026)

This is the pipeline behind the MSc thesis, kept unchanged under the tag `v1.0-thesis` for
reference. Do not use its results; see above.
