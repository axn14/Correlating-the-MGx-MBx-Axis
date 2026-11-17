# What changed in v2.0

The thesis pipeline (tag `v1.0-thesis`) reported several producer pairs that v2.0 does not
support. This page explains why, so that the earlier numbers are not reused. The full list is
in `CHANGELOG.md` in the repository.

## The sample groups were wrong

Yachida's cohort contains six groups. v1 placed post-surgery patients (HS) in Early CRC and
people with multiple adenomas (MP) in Advanced CRC, and put stage I–II cancers in Advanced
CRC as well. A bowel after tumour resection is neither a healthy control nor a tumour-bearing
gut, so the v1 comparison "healthy vs early CRC" was mostly healthy vs resected. That is
where most of v1's early-cancer species signal came from.

| Group | n | v1 label | v2 label |
|---|---|---|---|
| Healthy | 127 | Healthy | Healthy |
| Stage 0 | 27 | Early CRC | Early CRC |
| HS (after resection, cancer-free) | 30 | Early CRC | separate sensitivity arm |
| Stage I–II | 69 | Advanced CRC | Early CRC |
| Stage III–IV | 54 | Advanced CRC | Advanced CRC |
| MP (multiple polypoid adenomas) | 40 | Advanced CRC | separate adenoma arm |

## The evidence was counted more than once

v1 scored pairs across nine streams. Several of them (SHAP importance, MOFA co-loading,
mediation, within-stage correlation) are different readings of the same association in the
same samples, so a single strong correlation could collect several points. The "Trinity"
community-model pillar was evaluated at genus level and then applied to species. v2 assigns
each source to one of five independence classes and counts each class once. See
[[Evidence classes]].

## Targets were partly chosen by the result

v1 picked target metabolites partly by how well species predicted them. Testing afterwards
whether species account for those metabolites then answers itself. v2 admits a target only if
it differs between healthy people and cancer, or belongs to an enzyme class fixed from the
literature. A check stops the pipeline if target selection reads any species–metabolite
statistic.

## Nothing measured the chance rate

v2 permutes each class within its metabolite and recomputes the tiers 1,000 times. About 45 of
the 82 Supported links would be expected by chance. v1 had no equivalent, so its tiers could
not be interpreted.

## What still holds

*Veillonella* and putrescine, and *R. gnavus* and tryptamine, hold up under the stricter rules.
*Fusobacterium animalis* and phenylalanine remains a reasonable hypothesis: the association is
there, but there is no gene evidence and no replication.
