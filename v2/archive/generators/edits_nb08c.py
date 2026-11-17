"""NB08: remove all AGORA2 genome-scale evidence (R1, R3, R6's enzyme columns)
and rebuild the pair-evidence table without the genomic/mechanistic sources."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from nbedit import Notebook

nb = Notebook(Path(__file__).resolve().parent / "08_pathway_analysis.ipynb")

# ------------------------------------------------------------------ title
nb.set("nb08-00-title", r'''# 08 — Pathway-Level Evidence Assembly

**Primary output:** `tables/pair_evidence_nb08_v2.csv` — one row per (species, polyamine)
pair, with each independent line of evidence as its own column. NB09's final synthesis
**reads that file** rather than re-deriving the same quantities, so each column is computed
once, in one place, with one definition.

## Genome-scale evidence is removed in v2.1 — read this first

The original design had two AGORA2-based sources: **enzyme presence** (can this species make
the compound?) and **MICOM community flux** (does the modelled community export it?). Both
are gone, for an infrastructure reason, not a scientific one:

- the EC index needs 1.7 GB of AGORA2 reconstructions resident on disk to build;
- MICOM needs roughly 100 MB of model **per sample** (307 samples; an earlier attempt
  produced 36 GB and took ~4.6 h on 4 threads).

On the available 2-core / 1500 MiB machine the MICOM step saturated the host until it stopped
responding. Rather than keep a half-working dependency, v2.1 removes both.

**What this costs, stated plainly.** The remaining evidence is entirely statistical and
predictive. Nothing in this pipeline now establishes that a species *can* synthesise a
polyamine — only that its abundance co-varies with the metabolite, that a model finds it
predictive, that the association reproduces in another cohort, or that the two co-load on a
shared latent factor. A top-tier pair here is a **prioritised hypothesis for experimental or
genomic follow-up**, and the manuscript should say so rather than implying mechanism.

The curated EC map (`utils.POLYAMINE_EC`) is kept and still guarded by `checks.py` D1,
because it fixes a map that previously listed *consuming* enzymes as producers. Re-enabling
genomic evidence later therefore starts from a corrected map; `utils.py`'s AGORA2 section
documents the one-off procedure.

| Removed | Was |
|---|---|
| **R1** | AGORA2 enzyme evidence per (species, polyamine) + Fisher enrichment |
| **R3** | SCFA-guild enrichment (guild defined by AGORA2 EC numbers) |
| **R6** enzyme columns | per-*Fusobacterium* enzyme lookups (correlations retained) |
| **NB10** | MICOM community flux |
| NB09 **Module B** | AGORA2 producer enzyme profile |
| NB09 **Module G** | AGORA2-defined functional guild scores × MOFA factors |

## Structure

| Part | Content | Role |
|---|---|---|
| **A** | *Fusobacterium* context (R6, correlations only) | feeds Part B's context column |
| **B** | **`pair_evidence_nb08_v2.csv`** | **primary output, read by NB09** |
| **D** | *Supplementary* descriptive indices: oxidative stress (R4), oral invasion (R5), polyamine/dicarboxylic-acid overlap (R8) | reported, not evidence |

## Where the candidate pairs come from, and why it matters

The candidate pool is **NB02's partial correlations**
(`{cohort}_polyamine_partial_correlations_v2.csv`) — species–polyamine associations adjusted
for age, BMI, sex, smoking and alcohol — **not** NB03's SHAP ranking.

This is the fix for the central circularity of the earlier pipeline: SHAP chose which pairs
to examine, and SHAP-derived lists were then counted as evidence *for those same pairs*. Here
the adjusted association nominates a pair, and SHAP importance is one independent column
among several. A pair can be nominated and then fail every other axis — a result the old
design could not express.

## Evidence columns and what each is worth

| Column group | Source | Establishes | Does not establish |
|---|---|---|---|
| `Rho_partial`, `QValue_partial` | NB02 | the association survives confounder adjustment | direction of causation |
| `SHAP_*`, `Model_perm_q` | NB03 | the species carries predictive information about the metabolite | production |
| `Replicated_*`, `Replication_same_direction` | NB02 genus level | the association appears in an independent cohort | species-level identity (genus resolution; sometimes a proxy compound) |
| `Polyamine_stage_q`, `Species_stage_q` | NB04 | the compound/species tracks stage within cancer | a species–metabolite link |

`Evidence_count` sums **independent sources**, not columns: the SHAP columns are one source
because they come from one model. The statistical association is **not** counted — it is how
a pair was nominated, so counting it would give every row a free point.''')

# ------------------------------------------------- setup cell: drop AGORA2 index
nb.set("nb08-04-agora2", r'''# ---------------------------------------------------------------------------
# Shared helper. (The AGORA2 enzyme index that used to be built here is gone --
# see the title cell: genome-scale evidence is removed in v2.1. utils still
# carries the corrected EC map and the loader, but the loader now REFUSES to
# scan Data/micom from inside a notebook, which is what overloaded the server.)
# ---------------------------------------------------------------------------
assert utils.AGORA2_ENABLED is False, \
    'utils.AGORA2_ENABLED is True but this notebook was written for the disabled path'


def bh_family(df, pcols):
    """BH-correct every p-value column in `pcols` JOINTLY (one family per table)
    and add a matching Q_<col> for each. NaN p-values stay NaN."""
    pcols = [c for c in pcols if c in df.columns]
    if df.empty or not pcols:
        return df
    stacked = df[pcols].stack()                       # drops NaN
    if stacked.empty:
        return df
    q = pd.Series(multipletests(stacked.values, method='fdr_bh')[1], index=stacked.index)
    q = q.unstack()
    for c in pcols:
        df[f'Q_{c}'] = q[c].reindex(df.index) if c in q.columns else np.nan
    return df


print('Genome-scale (AGORA2) evidence: DISABLED in v2.1.')
print('  no enzyme-presence column, no MICOM flux column')
print('  evidence sources available here: statistical (nomination), predictive, replication')''')

# ------------------------------------------------- R1 removed
nb.set("nb08-05-r1-md", r'''## R1 — removed (AGORA2 enzyme evidence)

R1 asked, for every (species, polyamine) pair in NB03's positive-direction SHAP list, whether
that species' own AGORA2 reconstruction encodes an enzyme that **produces** that polyamine,
and then tested with Fisher's exact test whether SHAP producers were enriched for it.

It is removed because the enzyme index cannot be built on the available hardware (see the
title cell). This is the single biggest loss in v2.1: it was the only source that spoke to
metabolic *capability* rather than co-variation.

Two things are preserved so it can be restored cheaply:

- `utils.POLYAMINE_EC` — the **corrected** producer EC map. The pre-v2 version listed enzymes
  that *consume* ornithine and arginine as their producers, which inverted the producer call
  for those compounds; `checks.py` D1 still asserts the required producers are present and
  the known consumers are absent.
- `utils.build_agora2_ec_index()` and the species lookups, which need only the small cached
  JSON once it has been built once.

No R1 code cell follows.''')
nb.delete("nb08-06-r1")

# ------------------------------------------------- R3 removed
nb.set("nb08-08-r3-md", r'''## R3 — removed (SCFA-guild enrichment)

R3 tested whether short-chain-fatty-acid producers were enriched among the polyamine-producer
species, with the SCFA guild defined by AGORA2 EC numbers (`utils.SCFA_EC`: butyrate 2.8.3.8 /
2.7.2.7, propionate 2.8.3.1).

It is removed with the rest of the AGORA2 dependency. Note what would be required to keep it
honestly: the *previous* version of this test defined the SCFA guild from a curated genus list
that had itself been assembled from this cohort's SHAP top genera, which made "SCFA producers
are enriched among polyamine producers" partly a restatement of how the list was built. An
EC-based definition was the fix. Reverting to a literature genus list to keep the test alive
would reintroduce that circularity, so the test is dropped rather than weakened.

`utils.SCFA_EC` is retained for a future re-enable (`checks.py` D6).

No R3 code cell follows.''')
nb.delete("nb08-09-r3")

# ------------------------------------------------- R6: correlations only
nb.set("nb08-14-r6-md", r'''## Part A · R6 — *Fusobacterium* context (correlations only)

Spearman correlations between each *Fusobacterium* species and ornithine, putrescine and CRC
stage, BH-corrected within the table.

The AGORA2 enzyme columns this section used to carry
(`AGORA2_Encodes_Ornithine_Enzyme` / `..._Putrescine_Enzyme`) are gone, so what remains is
association only — it says whether *Fusobacterium* abundance tracks these compounds in this
cohort, not whether it can make them.

*Fusobacterium* gets its own section because it is the most-reported CRC-associated taxon in
the literature, so whether it appears in this cohort's polyamine associations is worth stating
either way. Part B carries an `Is_Fusobacterium` flag per pair, so the synthesis can report it
without this section having privileged status.''')

nb.set("nb08-15-r6", r'''orn_col = pa_col_map.get('C00077')
put_col = pa_col_map.get('C00134')

fuso_rows = []
for fsp in FUSO_COLS:
    x = spc[fsp]
    rho_o, p_o = spearmanr(x, mtb[orn_col]) if orn_col is not None else (np.nan, np.nan)
    rho_p, p_p = spearmanr(x, mtb[put_col]) if put_col is not None else (np.nan, np.nan)
    rho_s, p_s = spearmanr(x, CRC_STAGE)
    fuso_rows.append({
        'Species': fsp,
        'Rho_ornithine': rho_o, 'Pval_ornithine': p_o,
        'Rho_putrescine': rho_p, 'Pval_putrescine': p_p,
        'Rho_stage': rho_s, 'Pval_stage': p_s,
    })

fuso_df = pd.DataFrame(fuso_rows)
fuso_df = bh_family(fuso_df, ['Pval_ornithine', 'Pval_putrescine', 'Pval_stage'])
fuso_df.to_csv(TABLES / 'r6_fusobacterium_polyamine_correlations_v2.csv', index=False)

print(f"R6 -- {len(fuso_df)} Fusobacterium species evaluated (association only; "
      f"no enzyme evidence in v2.1).")
if not fuso_df.empty:
    print(fuso_df.round(4).to_string(index=False))
else:
    warnings.warn('R6: no Fusobacterium species found in the reduced species table.')
print("Saved: r6_fusobacterium_polyamine_correlations_v2.csv")''')

nb.save()
