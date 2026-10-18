"""NB09 edits: the final synthesis now READS NB08's pair-evidence table and adds
only the one source NB08 does not compute (MOFA same-sign co-loading)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from nbedit import Notebook

nb = Notebook(Path(__file__).resolve().parent / "09_pathway_network.ipynb")

nb.set("synthesis-md", r'''## Final Synthesis — Cross-Notebook Evidence Triangulation

**This section no longer recomputes evidence.** It reads
`tables/pair_evidence_nb08_v2.csv`, which NB08 assembles, and adds exactly one source that
NB08 cannot: MOFA **same-sign co-loading**, which needs NB06's loadings that this notebook
already has open.

That split is the point. Previously NB08 and NB09 each derived overlapping versions of the
same per-pair facts, so a pair's "evidence" depended on which notebook you read. Now there is
one table, one definition per column, and one place each column is computed.

### The five independent sources

| Source | Computed in | A pair is supported when |
|---|---|---|
| **statistical** | NB02 | *how the pair was nominated* — adjusted partial correlation q<0.05. Not counted, or every row would get a free point. |
| **predictive** | NB03 | positive-direction OOF SHAP **and** that target's model beat its permutation null (`Model_perm_q < 0.05`) |
| **genomic** | AGORA2 via NB08 | the species' reconstruction encodes an enzyme that **produces** the compound |
| **mechanistic** | NB10 | MICOM predicts a positive species-level flux contribution (`NaN`/pending until NB10 runs) |
| **replication** | NB02 | the genus-level association reproduces in ≥1 external cohort **with the same sign** |
| **co-loading** | *here*, from NB06 | species and metabolite load on a shared MOFA factor with the **same sign** (\|weight\| > 0.10) |

### Pre-declared tiers

Fixed before the numbers were seen, and counted over **independent sources**:
`Strong` ≥ 3, `Moderate` = 2, `Weak` = 1, `Statistical only` = 0.

### What this table is not

- **Co-loading is the weakest of the five.** Two features loading on a shared factor is
  consistent with a direct relationship but equally with both responding to a third thing —
  disease stage being the obvious candidate, since the factors are sign-oriented to stage.
  It is counted as one source, never as confirmation on its own.
- **Replication is genus-level.** A replicated row is evidence about the genus, not about the
  named species, and sometimes rests on an N-acetyl proxy compound. Both caveats travel in
  `Replication_cohorts` / `Replication_proxy_used`.
- **Missing coverage is not negative evidence.** A species with no AGORA2 reconstruction, or a
  polyamine NB03 could not model, yields `None` for that source — never a refutation. The
  availability counts printed below say how much evidence was *obtainable* per pair, which
  bounds what any tier can mean.
- **No causal claim.** Every source here is observational or genome-based; even a `Strong`
  pair is a prioritised hypothesis for experimental testing, not a demonstrated
  producer–product relationship.''')

nb.set("synthesis-build", r'''# ---------------------------------------------------------------------------
# Read NB08's assembled pair-level evidence table. This is the single source for
# the statistical / predictive / genomic / mechanistic / replication columns.
# ---------------------------------------------------------------------------
PAIR_EV = TABLES / 'pair_evidence_nb08_v2.csv'
if not PAIR_EV.exists():
    raise FileNotFoundError(
        f'{PAIR_EV.name} not found. It is produced by 08_pathway_analysis.ipynb (Part B); '
        'run NB08 before NB09. The synthesis deliberately does not fall back to '
        'recomputing it -- that duplication is what v2.1 removed.')

pair_ev = pd.read_csv(PAIR_EV)
print(f'Loaded {PAIR_EV.name}: {len(pair_ev)} pairs, {pair_ev.shape[1]} columns')

NB08_SOURCES = ['Evidence_predictive', 'Evidence_genomic',
                'Evidence_mechanistic', 'Evidence_replication']
missing = [c for c in NB08_SOURCES + ['Species', 'Polyamine', 'Rho_partial',
                                      'QValue_partial'] if c not in pair_ev.columns]
if missing:
    raise KeyError(f'{PAIR_EV.name} is missing expected column(s): {missing}. '
                   'Re-run 08_pathway_analysis.ipynb with the current version.')

for c in NB08_SOURCES:
    pair_ev[c] = pair_ev[c].fillna(False).astype(bool)

micom_status = pair_ev['MICOM_status'].iloc[0] if 'MICOM_status' in pair_ev.columns and len(pair_ev) else 'unknown'
print(f'MICOM (mechanistic) status: {micom_status}')
if micom_status == 'pending':
    print('  -> NB10 has not run yet, so the mechanistic column is unavailable and the')
    print('     attainable maximum is 4 of 5 sources. Re-run NB08 and NB09 after NB10.')

# ---------------------------------------------------------------------------
# The one source NB08 cannot compute: MOFA same-sign co-loading (needs NB06).
# ---------------------------------------------------------------------------
kegg_by_name = {v: k for k, v in utils.POLYAMINE_KEGG.items()}


def resolve_polyamine_column(kegg_id):
    if kegg_id is None or pd.isna(kegg_id):
        return None
    if kegg_id in polyamine_columns and polyamine_columns[kegg_id] in mtb_df.columns:
        return polyamine_columns[kegg_id]
    inv_map = {v: k for k, v in (mtb_kegg_map or {}).items()}
    if kegg_id in inv_map and inv_map[kegg_id] in mtb_df.columns:
        return inv_map[kegg_id]
    for col in mtb_df.columns:
        if col.startswith(str(kegg_id) + '_'):
            return col
    return None


MOFA_W_THRESH = 0.10
species_features_all, mtb_features_all = set(), set()
species_strong_factors, mtb_strong_factors = {}, {}

if HAS_MOFA_LOADINGS and not mofa_loadings.empty:
    cols_lower = {c.lower(): c for c in mofa_loadings.columns}
    view_col = cols_lower.get('view', 'View')
    feat_col = cols_lower.get('feature', 'Feature')
    fac_col  = cols_lower.get('factor', 'Factor')
    wt_col   = cols_lower.get('weight', 'Weight')

    species_features_all = set(mofa_loadings.loc[mofa_loadings[view_col] == 'species', feat_col])
    mtb_features_all     = set(mofa_loadings.loc[mofa_loadings[view_col] == 'metabolites', feat_col])

    strong = mofa_loadings[mofa_loadings[wt_col].abs() > MOFA_W_THRESH]
    # {feature: {factor: sign}} so co-loading can require the SAME sign, i.e. positive
    # co-variation, not merely "both matter on this factor".
    species_strong_factors = {
        f: dict(zip(g[fac_col], np.sign(g[wt_col])))
        for f, g in strong[strong[view_col] == 'species'].groupby(feat_col)
    }
    mtb_strong_factors = {
        f: dict(zip(g[fac_col], np.sign(g[wt_col])))
        for f, g in strong[strong[view_col] == 'metabolites'].groupby(feat_col)
    }
    print(f'MOFA loadings: {len(species_features_all)} species, '
          f'{len(mtb_features_all)} metabolites with |weight|>{MOFA_W_THRESH} available')
else:
    print('MOFA loadings unavailable -- the co-loading source will be marked unavailable '
          'for every pair (not counted as evidence against).')

co_avail, co_support, co_factors = [], [], []
for sp, pa in zip(pair_ev['Species'], pair_ev['Polyamine']):
    kegg_id = kegg_by_name.get(pa)
    mtb_col = resolve_polyamine_column(kegg_id)
    avail = bool(HAS_MOFA_LOADINGS and sp in species_features_all
                 and mtb_col is not None and mtb_col in mtb_features_all)
    if avail:
        sf = species_strong_factors.get(sp, {})
        mf = mtb_strong_factors.get(mtb_col, {})
        shared_same_sign = [f for f in set(sf) & set(mf) if sf[f] == mf[f]]
        co_avail.append(True)
        co_support.append(bool(shared_same_sign))
        co_factors.append(', '.join(sorted(shared_same_sign)) or '')
    else:
        co_avail.append(False)
        co_support.append(False)
        co_factors.append('')

pair_ev['Coloading_available']  = co_avail
pair_ev['Evidence_coloading']   = co_support
pair_ev['Coloading_factors']    = co_factors
print(f'\nCo-loading: available for {int(pd.Series(co_avail).sum())}/{len(pair_ev)} pairs, '
      f'supporting {int(pd.Series(co_support).sum())}')''')

nb.set("synthesis-score", r'''# ---------------------------------------------------------------------------
# Recount over all five independent sources and apply the pre-declared tiers.
# ---------------------------------------------------------------------------
ALL_SOURCES = NB08_SOURCES + ['Evidence_coloading']
SOURCE_LABEL = {
    'Evidence_predictive':  'predictive (NB03 SHAP + model q)',
    'Evidence_genomic':     'genomic (AGORA2 enzyme)',
    'Evidence_mechanistic': 'mechanistic (NB10 MICOM flux)',
    'Evidence_replication': 'replication (external cohort, same sign)',
    'Evidence_coloading':   'co-loading (NB06 MOFA, same sign)',
}

# Availability per source: what evidence was obtainable for this pair at all.
synth_df = pair_ev.copy()
synth_df['Avail_predictive']  = synth_df.get('Model_perm_q', pd.Series(np.nan, index=synth_df.index)).notna()
synth_df['Avail_genomic']     = synth_df.get('AGORA2_model_available', False)
synth_df['Avail_mechanistic'] = synth_df.get('MICOM_flux_q', pd.Series(np.nan, index=synth_df.index)).notna()
synth_df['Avail_replication'] = synth_df.get('Replication_n_cohorts', pd.Series(np.nan, index=synth_df.index)).notna()
synth_df['Avail_coloading']   = synth_df['Coloading_available']
AVAIL_COLS = ['Avail_predictive', 'Avail_genomic', 'Avail_mechanistic',
              'Avail_replication', 'Avail_coloading']
for c in AVAIL_COLS:
    synth_df[c] = synth_df[c].fillna(False).astype(bool)

synth_df['N_Sources_Available']  = synth_df[AVAIL_COLS].sum(axis=1).astype(int)
synth_df['N_Sources_Supporting'] = synth_df[ALL_SOURCES].sum(axis=1).astype(int)
synth_df['Evidence_sources'] = [
    ', '.join(s.replace('Evidence_', '') for s in ALL_SOURCES if row[s]) or 'none'
    for _, row in synth_df[ALL_SOURCES].iterrows()
]


def evidence_level(n):
    if n >= 3:
        return 'Strong (>=3 sources)'
    if n == 2:
        return 'Moderate (2 sources)'
    if n == 1:
        return 'Weak (1 source)'
    return 'Statistical only (0 sources)'


synth_df['Evidence_Level'] = synth_df['N_Sources_Supporting'].map(evidence_level)

LEAD_COLS = ['Species', 'Genus', 'Polyamine', 'KEGG_parent', 'Measured_as',
             'Rho_partial', 'QValue_partial', 'Candidate_producer',
             'Evidence_predictive', 'Evidence_genomic', 'Evidence_mechanistic',
             'Evidence_replication', 'Evidence_coloading',
             'N_Sources_Supporting', 'N_Sources_Available', 'Evidence_sources',
             'Evidence_Level',
             'Mean_SHAP_OOF', 'Direction', 'Model_perm_q',
             'AGORA2_encodes_enzyme', 'AGORA2_model_available',
             'MICOM_flux_rho', 'MICOM_flux_q', 'MICOM_status',
             'Replication_n_cohorts', 'Replication_same_direction',
             'Replication_cohorts', 'Replication_proxy_used',
             'Coloading_factors',
             'Polyamine_stage_rho', 'Polyamine_stage_q',
             'Species_stage_rho', 'Species_stage_q',
             'SCFA_producer', 'Is_Fusobacterium']
ordered = [c for c in LEAD_COLS if c in synth_df.columns]
ordered += [c for c in synth_df.columns if c not in ordered]
synth_df = synth_df[ordered].sort_values(
    ['N_Sources_Supporting', 'QValue_partial'], ascending=[False, True]).reset_index(drop=True)

synth_df.to_csv(TABLES / 'final_evidence_synthesis_v2.csv', index=False)
print(f'Saved: final_evidence_synthesis_v2.csv  ({len(synth_df)} pairs, '
      f'{synth_df.shape[1]} columns)')

if not synth_df.empty:
    print('\n--- Evidence obtainable per pair (bounds what any tier can mean) ---')
    av = synth_df['N_Sources_Available'].value_counts()
    for k in range(5, -1, -1):
        if av.get(k):
            print(f'  {k} of 5 source(s) available: {int(av.get(k))} pairs')

    print('\n--- Per-source availability and support ---')
    for s in ALL_SOURCES:
        a = s.replace('Evidence_', 'Avail_')
        print(f'  {SOURCE_LABEL[s]:42s} available {int(synth_df[a].sum()):4d}/'
              f'{len(synth_df)}, supporting {int(synth_df[s].sum()):4d}')

    print('\n--- Tier distribution (pre-declared thresholds) ---')
    print(synth_df['Evidence_Level'].value_counts().reindex(
        ['Strong (>=3 sources)', 'Moderate (2 sources)', 'Weak (1 source)',
         'Statistical only (0 sources)']).to_string())

    top = synth_df[synth_df['N_Sources_Supporting'] >= 2]
    print(f'\n--- Pairs with >=2 independent sources ({len(top)}) ---')
    if not top.empty:
        show = ['Species', 'Polyamine', 'Rho_partial', 'QValue_partial',
                'Evidence_sources', 'N_Sources_Supporting', 'N_Sources_Available',
                'Replication_proxy_used']
        print(top[[c for c in show if c in top.columns]].head(40).round(4).to_string(index=False))
    else:
        print('  none')
        print('  With MICOM pending this is 4 obtainable sources at most; the honest '
              'reading is that no pair is yet corroborated beyond its statistical\n'
              '  association plus at most one other line of evidence.')

    # Positive-direction candidates specifically: the plausible producers.
    prod = synth_df[synth_df.get('Candidate_producer', False) == True]   # noqa: E712
    print(f'\nOf {len(synth_df)} pairs, {len(prod)} are positive-direction '
          f'(candidate producers); {int((prod["N_Sources_Supporting"] >= 2).sum())} of those '
          f'have >=2 supporting sources.')''')

nb.set("synthesis-fig", r'''if not synth_df.empty and synth_df['N_Sources_Supporting'].max() >= 2:
    high_conf = synth_df[synth_df['N_Sources_Supporting'] >= 2]
    pivot_e = high_conf.pivot_table(index='Species', columns='Polyamine',
                                     values='N_Sources_Supporting', fill_value=0, aggfunc='max')
    fig, ax = plt.subplots(figsize=(max(6, len(pivot_e.columns) * 2.0),
                                     max(5, len(pivot_e) * 0.35 + 2)))
    sns.heatmap(pivot_e, cmap='YlOrRd', vmin=0, vmax=5, annot=True, fmt='d',
                linewidths=0.4, ax=ax,
                cbar_kws={'label': 'Independent sources supporting (of 5)'})
    ax.set_title('Final evidence synthesis — pairs with >=2 independent sources\n'
                 'sources: predictive (SHAP+model q) | genomic (AGORA2) | mechanistic (MICOM) | '
                 'replication (external, same sign) | co-loading (MOFA, same sign)\n'
                 'candidate pool = NB02 confounder-adjusted partial correlations, not SHAP',
                 fontsize=9.5, fontweight='bold')
    plt.tight_layout()
    out_synth = FIGURES / 'final_evidence_synthesis_heatmap_v2.png'
    fig.savefig(out_synth, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'Saved: {out_synth}')

    # Second panel: how many pairs each source supports, so the reader can see
    # which source is doing the work rather than inferring it from the heatmap.
    counts = {SOURCE_LABEL[s]: int(synth_df[s].sum()) for s in ALL_SOURCES}
    avails = {SOURCE_LABEL[s]: int(synth_df[s.replace('Evidence_', 'Avail_')].sum())
              for s in ALL_SOURCES}
    fig, ax = plt.subplots(figsize=(9, 4.2))
    ypos = range(len(counts))
    ax.barh(list(ypos), list(avails.values()), color='#dfe6e9',
            label='available (evidence obtainable)')
    ax.barh(list(ypos), list(counts.values()), color='#2a78d6', label='supporting')
    ax.set_yticks(list(ypos)); ax.set_yticklabels(list(counts), fontsize=8.5)
    ax.invert_yaxis()
    ax.set_xlabel(f'Pairs (of {len(synth_df)})')
    ax.set_title('Which evidence source supports how many pairs', fontsize=11, fontweight='bold')
    ax.legend(fontsize=8, loc='lower right')
    for i, (k, v) in enumerate(counts.items()):
        ax.text(v + len(synth_df) * 0.01, i, str(v), va='center', fontsize=8)
    plt.tight_layout()
    out_src = FIGURES / 'final_evidence_source_support_v2.png'
    fig.savefig(out_src, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'Saved: {out_src}')
else:
    print('Synthesis figure skipped -- no pairs with >=2 supporting sources.')''')

nb.save()
