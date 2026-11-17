# Evidence classes

Every source of evidence belongs to one class. A class supports a species–metabolite pair
when any of its members does, and it counts once however many members it has. A class that
cannot assess a pair is unavailable, which is not a vote against.

| Class | What it asks | Source |
|---|---|---|
| O, observational | Does the species track the metabolite in Yachida? | Prediction from species beyond the confounders (XGBoost, permutation test); MOFA+ same-sign co-loading |
| R, reference genome | Can the species' reference genome make it? | AGORA2 reconstructions, base-rate gated |
| G, patient gene reads | Do the patients' own reads carry a producing enzyme in this species? | HUMAnN 3 gene families; reads in at least 10% and 10 samples; enzyme carried by at most 20% of species |
| M, community model | Is the species forced to export it in a community model? | MICOM, certified solver, flux-variability minimum above zero |
| X, replication | Does the association replicate in an independent cohort? | Erawijantari 2020, same sign and q < 0.05, shared subjects removed |

## Tiers

| Tier | Rule |
|---|---|
| Conclusive | At least three classes, one of them G or M, no contradiction, target passes the gates |
| Supported | Two classes, at least one not O |
| Single-class | One class other than O |
| Observational only | O alone |
| Unsupported | None |

## Gates per target

- Export: does any of 1,440 AGORA2 models have an exchange reaction for the metabolite?
- Origin: do human enzymes also make it?
- Confounders: do age, BMI, sex, smoking, alcohol and depth alone predict it as well as the
  species do?

Genus-level agreement in the 16S cohorts (Sinha 2016, Kim 2020) is reported as ecological
context and never counted for a species.
