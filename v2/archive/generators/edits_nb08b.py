"""NB08 edits part 2: insert the pair-evidence builder (Part B) and relabel the
descriptive indices as supplementary (Part D)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from nbedit import Notebook

nb = Notebook(Path(__file__).resolve().parent / "08_pathway_analysis.ipynb")

# ============================================================ Part B (new)
nb.insert_after("nb08-15-r6", "nb08-16-pairev-md", "markdown", r'''---
# Part B · Pair-level evidence table — **the primary output**

One row per (species, polyamine) pair, assembled from every independent source. NB09's final
synthesis reads this file; it does not recompute any of these quantities.

**Candidate pool.** Pairs with `QValue_partial < 0.05` in NB02's confounder-adjusted partial
correlations. Both directions are retained: a negative association is evidence *against*
production and is informative, so it is kept and flagged rather than filtered out at this
stage. `Candidate_producer = (Rho_partial > 0)` marks the subset that could plausibly be
producers.

**`Evidence_count` counts independent sources, not columns:**

1. `statistical` — always true for rows in this table (that is how they were nominated), so
   it is **not** counted; counting it would give every row a free point.
2. `predictive` — positive-direction SHAP **and** the target's model beat chance
   (`Model_perm_q < 0.05`). A SHAP ranking from a model no better than chance is not
   evidence, which is why the model's own q-value gates this column.
3. `genomic` — AGORA2 encodes a producing enzyme (`True`, not `None`).
4. `mechanistic` — MICOM predicts a positive species-level flux contribution (NB10).
5. `replication` — the genus-level association reproduces in ≥1 external cohort with the
   same sign.

**Tiers** are a pre-declared function of that count, fixed before the numbers were seen:
`Strong` ≥ 3 independent sources, `Moderate` = 2, `Weak` = 1, `Statistical only` = 0.

**Proxy and resolution caveats travel with the row.** `Replication_proxy_used` records when
an external cohort only measured an N-acetyl derivative, and replication is genus-level
throughout — so a replicated row is evidence about a genus, not about the species named in
`Species`. Those columns exist so the limitation cannot be lost when the table is read.''')

nb.insert_after("nb08-16-pairev-md", "nb08-16-pairev", "code", r'''# ---------------------------------------------------------------------------
# Part B -- assemble the pair-level evidence table.
# ---------------------------------------------------------------------------
Q_CAND = utils.FDR_THRESHOLD

cand = partial_corr[partial_corr['QValue_partial'] < Q_CAND].copy()
cand = cand[['Species', 'Polyamine', 'Metabolite', 'Rho_partial', 'QValue_partial',
             'PValue_partial', 'KEGG', 'Measured_as', 'Covariates']].copy()
cand['Candidate_producer'] = cand['Rho_partial'] > 0
print(f'Candidate pairs (q<{Q_CAND}): {len(cand)}  '
      f'({int(cand["Candidate_producer"].sum())} positive / candidate producers)')

# ---- 1. predictive evidence (NB03 SHAP, gated on the model beating chance) ----
if not pa_producers.empty:
    shap_cols = ['Species', 'Polyamine', 'Mean_SHAP_OOF', 'Mean_SHAP_signed',
                 'SHAP_dir_rho', 'Direction', 'Rank', 'Rank_stable',
                 'Model_AUC_nested', 'Model_perm_q']
    shap_src = (pa_producers[[c for c in shap_cols if c in pa_producers.columns]]
                .drop_duplicates(subset=['Species', 'Polyamine']))
    ev = cand.merge(shap_src, on=['Species', 'Polyamine'], how='left',
                    suffixes=('', '_shap'))
else:
    ev = cand.copy()
    for c in ['Mean_SHAP_OOF', 'Mean_SHAP_signed', 'SHAP_dir_rho', 'Rank',
              'Model_AUC_nested', 'Model_perm_q']:
        ev[c] = np.nan
    ev['Direction'] = None
    ev['Rank_stable'] = False

ev['SHAP_positive'] = ev.get('Direction').eq('positive') if 'Direction' in ev else False
ev['Model_beats_chance'] = ev['Model_perm_q'] < utils.FDR_THRESHOLD
ev['Evidence_predictive'] = ev['SHAP_positive'] & ev['Model_beats_chance']

# ---- 2. genomic evidence (tri-state; None = no genomic reference available) ----
# genomic_evidence() is defined in NB08's evidence-setup cell. While no genomic
# reference index is built it returns None throughout, so AGORA2_model_available is
# False everywhere and Evidence_genomic is never counted as a negative.
kegg_by_name = {v: k for k, v in utils.POLYAMINE_KEGG.items()}
ev['KEGG_parent'] = ev['Polyamine'].map(kegg_by_name)
ev['AGORA2_encodes_enzyme'] = [
    genomic_evidence(sp, kid) if pd.notna(kid) else None
    for sp, kid in zip(ev['Species'], ev['KEGG_parent'])
]
ev['AGORA2_model_available'] = ev['AGORA2_encodes_enzyme'].notna()
ev['Evidence_genomic'] = ev['AGORA2_encodes_enzyme'].eq(True)

# ---- 3. SCFA guild + Fusobacterium context (per species) ----
ev['SCFA_producer'] = [scfa_ev.get(sp) for sp in ev['Species']]
ev['Is_Fusobacterium'] = ev['Species'].str.contains('Fusobact', case=False, na=False)

# ---- 4. stage trends (NB04): per-polyamine and per-species, joined separately ----
if pa_stage_trend is not None and not pa_stage_trend.empty:
    t = pa_stage_trend[['Polyamine', 'Rho_adjusted', 'QValue_adjusted']].rename(
        columns={'Rho_adjusted': 'Polyamine_stage_rho',
                 'QValue_adjusted': 'Polyamine_stage_q'})
    ev = ev.merge(t, on='Polyamine', how='left')
else:
    ev['Polyamine_stage_rho'] = np.nan
    ev['Polyamine_stage_q'] = np.nan

if spc_stage_trend is not None and not spc_stage_trend.empty:
    t = spc_stage_trend[['Species', 'Rho_adjusted', 'QValue_adjusted']].rename(
        columns={'Rho_adjusted': 'Species_stage_rho',
                 'QValue_adjusted': 'Species_stage_q'})
    ev = ev.merge(t, on='Species', how='left')
else:
    ev['Species_stage_rho'] = np.nan
    ev['Species_stage_q'] = np.nan

# ---- 5. cross-cohort replication (genus level) ----
# Species -> genus, normalised the same way NB02's replication step normalised it,
# so "Blautia_A wexlerae" and genus "Blautia" match.
ev['Genus'] = [utils.normalise_genus_name(str(s).split()[0]) if str(s).strip() else ''
               for s in ev['Species']]

if replication is not None and not replication.empty:
    rep = replication.copy()
    gcol = 'Genus' if 'Genus' in rep.columns else 'Species'
    rep['_genus'] = [utils.normalise_genus_name(g) for g in rep[gcol]]
    rep['_pa'] = rep['Polyamine'] if 'Polyamine' in rep.columns else rep['Metabolite_key']
    keep = ['_genus', '_pa', 'n_datasets', 'mean_rho', 'Same_direction',
            'Cohorts_tested', 'Proxy_used_in']
    rep = rep[[c for c in keep if c in rep.columns]].drop_duplicates(subset=['_genus', '_pa'])
    ev = ev.merge(rep, left_on=['Genus', 'Polyamine'], right_on=['_genus', '_pa'], how='left')
    ev = ev.drop(columns=[c for c in ('_genus', '_pa') if c in ev.columns])
    ev = ev.rename(columns={'n_datasets': 'Replication_n_cohorts',
                            'mean_rho': 'Replication_mean_rho',
                            'Same_direction': 'Replication_same_direction',
                            'Cohorts_tested': 'Replication_cohorts',
                            'Proxy_used_in': 'Replication_proxy_used'})
else:
    for c in ['Replication_n_cohorts', 'Replication_mean_rho',
              'Replication_same_direction', 'Replication_cohorts',
              'Replication_proxy_used']:
        ev[c] = np.nan

ev['Replication_same_direction'] = ev['Replication_same_direction'].fillna(False).astype(bool)
ev['Evidence_replication'] = ev['Replication_same_direction']

# ---- 6. MICOM mechanistic evidence (pending until NB10 has run) ----
if micom_pairs is not None and not micom_pairs.empty:
    m = micom_pairs.rename(columns={'Rho': 'MICOM_flux_rho', 'QValue': 'MICOM_flux_q'})
    mcols = [c for c in ['Species', 'Polyamine', 'MICOM_flux_rho', 'MICOM_flux_q'] if c in m.columns]
    ev = ev.merge(m[mcols].drop_duplicates(subset=['Species', 'Polyamine']),
                  on=['Species', 'Polyamine'], how='left')
    ev['Evidence_mechanistic'] = (ev['MICOM_flux_q'] < utils.FDR_THRESHOLD) & (ev['MICOM_flux_rho'] > 0)
    ev['MICOM_status'] = np.where(ev['MICOM_flux_q'].notna(), 'computed', 'not_modelled')
else:
    ev['MICOM_flux_rho'] = np.nan
    ev['MICOM_flux_q'] = np.nan
    ev['Evidence_mechanistic'] = False
    ev['MICOM_status'] = 'pending'

# ---- Evidence count + pre-declared tiers ----
EVIDENCE_SOURCES = ['Evidence_predictive', 'Evidence_genomic',
                    'Evidence_mechanistic', 'Evidence_replication']
for c in EVIDENCE_SOURCES:
    ev[c] = ev[c].fillna(False).astype(bool)
ev['Evidence_count'] = ev[EVIDENCE_SOURCES].sum(axis=1).astype(int)
ev['Evidence_sources'] = [
    ', '.join(s.replace('Evidence_', '') for s in EVIDENCE_SOURCES if row[s]) or 'none'
    for _, row in ev[EVIDENCE_SOURCES].iterrows()
]
ev['Tier'] = pd.cut(ev['Evidence_count'], bins=[-1, 0, 1, 2, 99],
                    labels=['Statistical only', 'Weak', 'Moderate', 'Strong'])

ev = ev.sort_values(['Evidence_count', 'QValue_partial'],
                    ascending=[False, True]).reset_index(drop=True)

out_pair = TABLES / 'pair_evidence_nb08_v2.csv'
ev.to_csv(out_pair, index=False)

print(f'\n=== pair_evidence_nb08_v2.csv: {len(ev)} pairs ===')
print(f'Saved: {out_pair.name}   (read by 09_pathway_network.ipynb)')
print('\nEvidence source support (independent sources; statistical is how a pair was '
      'nominated and is not counted):')
for c in EVIDENCE_SOURCES:
    print(f'  {c.replace("Evidence_", ""):13s}: {int(ev[c].sum()):4d} pairs')
print(f'  MICOM status : {ev["MICOM_status"].iloc[0] if len(ev) else "n/a"}')

print('\nTier distribution:')
print(ev['Tier'].value_counts().reindex(
    ['Strong', 'Moderate', 'Weak', 'Statistical only']).to_string())

print('\nAGORA2 coverage of candidate species (caps how much genomic evidence is possible):')
n_cov = int(ev['AGORA2_model_available'].sum())
print(f'  {n_cov}/{len(ev)} pairs ({100 * n_cov / max(len(ev), 1):.0f}%) have a '
      f'reconstruction for the species; the rest are "no evidence", not "no enzyme".')

if (ev['Evidence_count'] >= 2).any():
    print('\nPairs with >=2 independent evidence sources:')
    show = ['Species', 'Polyamine', 'Rho_partial', 'QValue_partial',
            'Evidence_sources', 'Evidence_count', 'Tier', 'Replication_proxy_used']
    print(ev[ev['Evidence_count'] >= 2][[c for c in show if c in ev.columns]]
          .head(30).round(4).to_string(index=False))
else:
    print('\nNo pair reaches 2 independent evidence sources yet. With MICOM pending '
          '(NB10) the mechanistic column is unavailable, so the attainable maximum is '
          'currently 3.')

pair_evidence = ev   # name used by the summary cell''')

# ============================================================ Part D relabels
nb.set("nb08-10-r4-md", r'''---
# Part D · Supplementary descriptive indices

**These are descriptive, not evidence.** None of R4, R5 or R8 contributes a column to
`pair_evidence_nb08_v2.csv`, and NB09's synthesis does not read them. They are retained
because they describe the cohort and were part of the original analysis plan, and because
removing them silently would hide what was looked at.

Each index is a composite of a handful of metabolites or species, so it inherits whichever of
those columns survived NB02's top-variance reduction — the counts printed below say which.
An index built from two of four intended markers is reported as such rather than presented as
"the" index.

---

## R4 · Oxidative-stress metabolite index vs CRC stage / MOFA factors

An oxidative-stress index built directly from urate + glucuronate abundance, correlated with
CRC stage, with the best-correlating fitted MOFA factor reported as a data-driven candidate
"oxidative-stress axis".

The factor is **identified, not assumed**: the old notebook hardcoded "Factor 5 = oxidative
stress", but NB06's factor count and ordering vary between runs, so a fixed number is
meaningless. The index never depends on a factor number.''')

nb.set("nb08-12-r5-md", r'''## R5 · Oral-invasion index + dicarboxylic-acid correlations  *(supplementary)*

A composite oral-invasion score from the canonical oral-invader species, correlated with CRC
stage (plus an Age/BMI-adjusted partial correlation via OLS residualisation), the fitted MOFA
factor it aligns with best, and the individual dicarboxylic acids against stage / oral
invasion / that factor.

Interpretation is limited by construction: the index is a mean of up to three species that
may or may not have survived feature reduction, and "oral invasion" is inferred from the
presence of oral-associated taxa in stool, not observed directly.''')

nb.set("nb08-17-r8-md", r'''## R8 · Polyamine / dicarboxylic-acid overlap  *(supplementary)*

Correlations between polyamines and dicarboxylic acids, BH-corrected, written
unconditionally because NB09 checks for this file's existence.

**Known limitation, left as-is.** This section reports every pair passing its threshold, which
on a cohort this size is a large and not very selective list — it describes co-variation in
the metabolome rather than identifying anything specific. It is kept for continuity with the
original analysis plan and is explicitly *not* used as evidence for any species–metabolite
pair. A threshold that would make it selective would be chosen after seeing the data, which
is precisely what the rest of this pipeline avoids.''')

nb.save()
