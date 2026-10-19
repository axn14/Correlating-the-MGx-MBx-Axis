"""NB08 edits: restructure into (A) a pair-level evidence table that NB09 consumes,
(B) the enrichment tests, (C) clearly-labelled supplementary descriptive indices."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from nbedit import Notebook

nb = Notebook(Path(__file__).resolve().parent / "08_pathway_analysis.ipynb")

# ------------------------------------------------------------------ title
nb.set("nb08-00-title", r'''# 08 — Pathway-Level Evidence Assembly

**This notebook's job changed in v2.1.** It used to be a bag of loosely related
pathway analyses (R1–R8) whose outputs NB09 then partly recomputed for itself. It now has
one primary output:

> **`tables/pair_evidence_nb08_v2.csv`** — one row per (species, polyamine) pair, with every
> independent line of evidence for that pair as its own column.

NB09's final synthesis **reads that file** instead of re-deriving the same quantities, so a
pair's evidence is computed once, in one place, with one definition.

## Structure

| Part | Content | Role |
|---|---|---|
| **A** | Per-species / per-pair facts: AGORA2 enzyme evidence (R1), SCFA guild membership (R3), *Fusobacterium* context (R6) | feed Part B |
| **B** | **`pair_evidence_nb08_v2.csv`** — the assembled pair-level table | **primary output, read by NB09** |
| **C** | Enrichment tests: are SHAP producers enriched for the biosynthetic enzyme (R1), are polyamine producers enriched for SCFA capability (R3) | standalone results |
| **D** | *Supplementary* descriptive indices: oxidative stress (R4), oral invasion (R5), polyamine/dicarboxylic-acid overlap (R8) | reported, not evidence |

## Where the candidate pairs come from, and why it matters

The candidate pool is **NB02's partial correlations** (`{cohort}_polyamine_partial_correlations_v2.csv`)
— species–polyamine associations adjusted for age, BMI, sex, smoking and alcohol — **not**
NB03's SHAP ranking.

This is the fix for the central circularity of the earlier pipeline: SHAP chose which pairs
to look at, and SHAP-derived lists were then counted as evidence *for those same pairs*. Here
the statistical association nominates a pair, and SHAP importance is one independent column
among several. A pair can therefore be nominated and then fail on every other axis, which is
a result the old design could not express.

## Evidence columns, and what each one is worth

| Column group | Source | What it establishes | What it does not |
|---|---|---|---|
| `Rho_partial`, `QValue_partial` | NB02 | the association exists after confounder adjustment | nothing about direction of causation |
| `SHAP_*`, `Model_perm_q` | NB03 | the species carries predictive information about the metabolite | not production — and only if `Model_perm_q < 0.05` |
| `AGORA2_encodes_enzyme` | AGORA2 GENRE | the genome encodes a **producing** enzyme (`utils.POLYAMINE_EC`) | not that the gene is expressed here |
| `Polyamine_stage_q`, `Species_stage_q` | NB04 | the compound/species tracks stage within cancer | not a species–metabolite link |
| `Replicated_*`, `Replication_same_direction` | NB02 genus level | the association appears in an independent cohort | only at genus resolution, sometimes via a proxy compound |
| `MICOM_flux_*` | NB10 | community metabolic modelling predicts production flux | `NaN` / "pending" until NB10 has run |

`Evidence_count` sums the **independent** sources that support a pair. It deliberately counts
*sources*, not columns: the three SHAP columns are one source, because they come from one
model.

**`AGORA2_encodes_enzyme` is tri-state.** `True` / `False` only when AGORA2 has a
reconstruction for that species; `None` means *no evidence either way* and is never counted
as a refutation. The share of species with no model is reported, because it caps how much
genomic evidence this design can ever supply.

## Fixes carried over from the v2 review

1. **Genus-level KEGG organism lookups replaced by species-level AGORA2 evidence.** The old
   approach inferred a study species' enzyme content from *any* organism sharing its genus in
   KEGG — a congener, not the species — over a live, uncached REST dependency. There are now
   no network calls anywhere in this notebook.
2. **R2 dropped.** It tested five pre-fetched KOs against an unrelated module's gene set and
   was non-functional; AGORA2 species-level evidence subsumes it.
3. **R7 dropped.** Its mediation design took its mediator list from the same SHAP ranking it
   then re-tested (see `07_mediation.ipynb` for the full reasoning).
4. **Dual-producer labelling bug fixed.** The old R1 built a combined list and then labelled
   every entry by membership in the first list, so a genuine dual ornithine/putrescine
   producer silently lost one of its two rows. Labels are now attached at construction.
5. **Only positive-direction SHAP species are candidate producers** (NB03's `Direction`):
   a species whose higher abundance pushes the prediction towards *low* metabolite is not a
   producer, and several top-|SHAP| species are exactly that.
6. **The polyamine EC map was corrected** — it previously listed enzymes that *consume*
   ornithine and arginine as their producers, which inverted the producer call for those
   compounds.
7. **The SCFA guild is EC-defined** (`utils.SCFA_EC`), not a curated genus list assembled
   from this same cohort's SHAP output, which had made that enrichment circular.
8. **BH correction within every table** (`Q_*` columns), and `r8_*` is always written because
   NB09 depends on its existence.
9. **MOFA factors are identified dynamically**, never by hardcoded number: the old notebook
   asserted "Factor 5 = oxidative stress" and "Factor 3 = oral invasion", but NB06's factor
   count and ordering vary by run, so each index is computed from abundances and the
   best-correlating fitted factor is then *reported*.''')

# ------------------------------------------------------------------ imports / inputs
nb.set("nb08-01-imports", r'''import sys
sys.path.insert(0, '.')

import pickle
import warnings

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm
from scipy.stats import spearmanr, kruskal, fisher_exact
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.multitest import multipletests

import utils

RESULTS_DIR = utils.RESULTS_DIR
TABLES      = RESULTS_DIR / 'tables'
FIGURES     = RESULTS_DIR / 'figures' / 'pathway'
TABLES.mkdir(parents=True, exist_ok=True)
FIGURES.mkdir(parents=True, exist_ok=True)

YACHIDA = utils.DATASET_DISCOVERY

with open(RESULTS_DIR / 'intermediate' / 'preprocessed_data_reduced_v2.pkl', 'rb') as f:
    pdata = pickle.load(f)
with open(RESULTS_DIR / 'intermediate' / 'ml_results_v2.pkl', 'rb') as f:
    ml_results = pickle.load(f)

fac_scores   = pd.read_csv(TABLES / 'mofa_factor_scores_v2.csv')
fac_loadings = pd.read_csv(TABLES / 'mofa_factor_loadings_v2.csv')

spc_df  = pdata['transformed_species'][YACHIDA].copy()
mtb_df  = pdata['transformed_mtb'][YACHIDA].copy()
meta_df = pdata['harmonized_meta'][YACHIDA].copy()
polyamine_columns = pdata.get('polyamine_columns', {}).get(YACHIDA, {})   # {kegg_id: column}
mtb_kegg_map      = pdata.get('mtb_kegg_map', {}).get(YACHIDA, {})        # {column: kegg_id}

print(f"Pickle sample scope: {pdata.get('sample_scope', 'unspecified')}")
assert set(meta_df['Stage.Group'].dropna()) <= set(utils.STAGE_GROUP_ORDER), \
    'non-primary group reached NB08'

# ---------------------------------------------------------------------------
# Inputs for the pair-evidence table (Part B). Each is REQUIRED except MICOM,
# which is pending until NB10 runs; a missing required input is a hard error so
# the evidence table can never be silently built with a column of all-NaN.
# ---------------------------------------------------------------------------
def _require_csv(path, produced_by):
    if not Path(path).exists():
        raise FileNotFoundError(
            f'{Path(path).name} not found. It is produced by {produced_by}; run that '
            'notebook before NB08. (Set V2_ALLOW_MISSING=1 only for debugging.)')
    return pd.read_csv(path)


import os
ALLOW_MISSING = os.environ.get('V2_ALLOW_MISSING') == '1'


def _optional_csv(path, produced_by):
    p = Path(path)
    if p.exists():
        return pd.read_csv(p)
    print(f'  [pending] {p.name} absent -- produced by {produced_by}; '
          f'its columns will be NaN and marked pending.')
    return None


# NB02: the SHAP-independent candidate pool.
partial_corr = _require_csv(
    TABLES / f'{YACHIDA}_polyamine_partial_correlations_v2.csv', '02_analysis.ipynb')
# NB02: genus-level cross-cohort replication.
replication = _optional_csv(
    TABLES / 'consistent_polyamine_correlations_genus_v2.csv', '02_analysis.ipynb')
# NB04: within-cancer stage trends.
pa_stage_trend  = _optional_csv(TABLES / 'nb04_polyamine_stage_trend_v2.csv',
                                '04_crc_stage_tumour_site.ipynb')
spc_stage_trend = _optional_csv(TABLES / 'nb04_species_stage_trend_v2.csv',
                                '04_crc_stage_tumour_site.ipynb')
# NB10: MICOM community fluxes (runs separately, after the main pipeline).
micom_pairs = _optional_csv(TABLES / 'nb10_micom_species_flux_pairs_v2.csv',
                            '10_micom_flux.ipynb (run_micom.sh)')

pa_producers_all = ml_results.get(
    'pa_producers', pd.DataFrame(columns=['Species', 'Polyamine', 'Rank']))
# Candidate producers = species whose higher abundance pushes the model towards HIGH
# metabolite levels (NB03 'Direction'); negatively associated species are not producers.
if 'Direction' in pa_producers_all.columns:
    pa_producers = pa_producers_all[pa_producers_all['Direction'] == 'positive'].copy()
else:
    warnings.warn('pa_producers has no Direction column (old NB03 output) -- using all rows')
    pa_producers = pa_producers_all.copy()
print(f"\nProducer rows: {len(pa_producers_all)} total, "
      f"{len(pa_producers)} positive-direction")
print(f"Partial-correlation candidate pool: {len(partial_corr)} pairs "
      f"({int((partial_corr['QValue_partial'] < 0.05).sum())} at q<0.05, "
      f"{int(((partial_corr['QValue_partial'] < 0.05) & (partial_corr['Rho_partial'] > 0)).sum())} positive)")

print(f"\nLoaded: spc={spc_df.shape}, mtb={mtb_df.shape}, meta={meta_df.shape}, "
      f"fac_scores={fac_scores.shape}, fac_loadings={fac_loadings.shape}")

# -- Align every table to the same sample set (never assume matching row order) --
common_samples = sorted(set(spc_df['Sample']) & set(mtb_df['Sample'])
                        & set(meta_df['Sample']) & set(fac_scores['Sample']))
print(f"Common samples across spc/mtb/meta/factor-scores: {len(common_samples)}")

spc  = spc_df.set_index('Sample').loc[common_samples]
mtb  = mtb_df.set_index('Sample').loc[common_samples]
meta = meta_df.set_index('Sample').loc[common_samples]
fac  = fac_scores.set_index('Sample').loc[common_samples]

CRC_STAGE   = meta['CRC.Stage.Numeric'].astype(float)   # utils' corrected tiers
STAGE_GROUP = meta['Stage.Group']                        # Healthy/Early_CRC/Advanced_CRC
STUDY_GROUP = meta['Study.Group']''')

# ------------------------------------------------------------------ R1 header
nb.set("nb08-05-r1-md", r'''## Part A · R1 — AGORA2 enzyme evidence for SHAP producer species

For every (species, polyamine) pair in NB03's positive-direction SHAP producer list, does
that species' own AGORA2 reconstruction encode an enzyme that **produces** that polyamine
(`utils.species_encodes_enzyme`, EC map `utils.POLYAMINE_EC`)?

This cell produces the per-pair genomic fact that Part B joins into the evidence table, and
the enrichment test in Part C uses the same evidence.

**Tri-state, not boolean.** `True`/`False` only where AGORA2 has a model for the species;
`None` means no evidence either way and is excluded from the enrichment test rather than
counted as a negative. Treating "no model" as "no enzyme" would manufacture refutations for
poorly-characterised taxa — which are exactly the taxa a discovery analysis surfaces.''')

# ------------------------------------------------------------------ R3 header
nb.set("nb08-08-r3-md", r'''## Part C · R3 — SCFA / polyamine-producer enrichment

Fisher's exact test for enrichment of short-chain-fatty-acid (SCFA) producers among the
positive-direction SHAP polyamine-producer species — overall and per polyamine, BH-corrected
as one family.

**SCFA producer = the AGORA2 model encodes a butyrate (2.8.3.8, 2.7.2.7) or propionate
(2.8.3.1) synthesis enzyme** (`utils.SCFA_EC`, matched anywhere in the model via
`utils.species_has_any_ec`).

The previous version defined the SCFA guild from a curated genus list that had itself been
assembled from this cohort's SHAP top genera — so "SCFA producers are enriched among
polyamine producers" was partly a restatement of how the list was built. An EC-based
definition is independent of this cohort entirely.

Test universe: reduced species with an AGORA2 model. Species without one carry no evidence
and are excluded rather than treated as non-SCFA.''')

# ------------------------------------------------------------------ R6 header
nb.set("nb08-14-r6-md", r'''## Part A · R6 — *Fusobacterium* context

Direct Spearman correlations between each *Fusobacterium* species and ornithine, putrescine
and CRC stage, plus AGORA2 species-level evidence for whether each one encodes the terminal
ornithine/putrescine biosynthesis enzyme.

This replaces the old notebook's hardcoded, unverified `FUSO_CODES = {'fnn', 'fpe', ...}`
KEGG organism-code list with the same species-level AGORA2 lookup R1 uses.

*Fusobacterium* gets its own section because it is the most-reported CRC-associated taxon in
the literature, so whether it appears in this cohort's polyamine associations is worth
stating explicitly either way. Part B carries a `Fusobacterium` flag per pair so the final
synthesis can report it without this section having privileged status.''')

nb.save()
