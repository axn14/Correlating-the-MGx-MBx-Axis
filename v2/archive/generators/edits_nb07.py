"""NB07 edits: update the exclusion rationale for the v2.1 group policy."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from nbedit import Notebook

nb = Notebook(Path(__file__).resolve().parent / "07_mediation.ipynb")

nb.set("2dfe4d39", r'''# 07 — Mediation Analysis (excluded from v2.1)

This notebook is intentionally **not** a working analysis. It records what the original
exploratory mediation analysis (`exploratory/07_mediation.ipynb`) attempted and why it was
dropped rather than repaired, mirroring the same decision made earlier for
`CRC_Gut_Analysis_Git/05_mediation_network.ipynb`.

It is kept as a file so the exclusion is a documented decision rather than a gap somebody
later "fixes" by reinstating the original code.

## What it tried

Species → polyamine → CRC-stage mediation via `statsmodels.stats.mediation.Mediation`
(Imai–Keele–Tingley potential-outcomes simulation) for each (species, polyamine) pair drawn
from the top-5 SHAP-ranked "producer" species per polyamine.

## Why it was dropped, not fixed

**1. Structural circularity.** Mediator selection (top-5 SHAP producers per polyamine) and
the mediation test itself (does that species → polyamine path mediate the stage outcome) drew
on the *same model fit to the same cohort*. The species most likely to pass were selected
precisely because they already looked like strong producers in that dataset, so the a-path is
close to guaranteed to look strong regardless of the true causal structure. That is not a bug
to patch; it is the same circular reasoning that removed mediation from the git-tracked
pipeline's deliverables.

**2. A real implementation bug compounded it.** `Mediation.fit()` hardcodes a treatment=0 vs
treatment=1 contrast (confirmed from the installed `statsmodels` source — this is the
Imai–Keele–Tingley design, not the classical Baron & Kenny approach the original notebook's
markdown claimed). The exposure (raw species abundance) was never rescaled into that 0–1
range, unlike the z-scored covariates in the same model, so even setting circularity aside
the reported ACME values were not per-unit effects and were not comparable across species —
yet they were ranked against each other in a forest plot.

**3. Its outcome variable inherited the stage-coding bug.** The outcome was the old,
incorrectly ordered `CRC.Stage.Numeric`, in which MP (multiple polypoid adenomas, a
cancer-free lesion) was coded as *more* severe than Stage III/IV cancer. That was the
highest-stakes place the bug reached, because this was a formal causal-effect estimate rather
than a descriptive correlation.

v2.1 removes that encoding problem at the source — `utils.py` excludes MP from the study
entirely and places HS ("healthy after surgery") in its own off-axis `PostSurgery` group with
a NaN ordinal position, so no severity scale contains a non-cancer group. But fixing the
labels does not touch findings 1–2, which are the reasons mediation stays excluded.

## A further reason, specific to v2.1

Mediation needs a defensible causal ordering, and this design cannot supply one. All three
cohorts are **cross-sectional**, with one sample per subject (verified in NB04's predecessor:
unique `Subject` per row). Species abundance, metabolite level and disease stage are measured
at the same moment, so "species → metabolite → stage" is an assumed direction, not an observed
one. The reverse path — a tumour altering the gut environment, which alters both the community
and its metabolites — is equally consistent with every number in this pipeline.

## What v2.1 does instead

Evidence for species–metabolite relationships is assembled, **without** any causal claim, as
explicitly separate and independently computed sources:

| Source | Notebook | Claim it supports |
|---|---|---|
| confounder-adjusted partial correlation | `02_analysis.ipynb` | the association exists net of age/BMI/sex/smoking/alcohol |
| out-of-fold SHAP + permutation-tested model | `03_ml_benchmarking.ipynb` | the species carries predictive information about the metabolite |
| AGORA2 genome-scale enzyme presence | `08_pathway_analysis.ipynb` | the genome *can* produce the compound |
| MICOM community flux | `10_micom_flux.ipynb` | community modelling predicts production flux |
| external-cohort replication | `02_analysis.ipynb` | the association appears in an independent population |

`09_pathway_network.ipynb` counts how many of those independent sources support each pair and
assigns a pre-declared tier. A `Strong` pair there is a prioritised hypothesis for
experimental testing — which is what this data can support — and is never described as a
demonstrated causal pathway.''')

nb.save()
