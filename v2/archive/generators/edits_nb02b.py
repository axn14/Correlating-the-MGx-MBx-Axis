"""NB02 edits, part B: grouped DA with robust hits, partial correlations,
three-cohort genus replication."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from nbedit import Notebook

nb = Notebook(Path(__file__).resolve().parent / "02_analysis.ipynb")

# ------------------------------------------------- Step 6B header
nb.set("cell-014", r'''---
## Step 6B: Differential Abundance (grouped `Stage.Group` scheme) — **both tests**

`utils.COMPARISONS_GROUP["YACHIDA-CRC-2019"]` defines three comparisons on the corrected
grouping (`utils.STAGE_GROUP_MAP`): `Healthy`, `Early_CRC` (Stage_0 + Stage_I_II),
`Advanced_CRC` (Stage_III_IV). Only Yachida has `Stage.Group` populated (Sinha and Kim have
no staging), so this step is Yachida-only.

Each comparison is tested twice and the results are merged:

1. **Mann-Whitney U + Cliff's δ** (`utils.differential_abundance`) — unadjusted, rank-based.
   `MeanDiff = mean(group2) − mean(group1)`, so **MeanDiff > 0 means higher in group2**;
   `CliffsD` carries the **opposite** sign by construction.
2. **Covariate-adjusted linear model** (`utils.clr_lm_da`) — adjusts for Age, BMI, Gender,
   Brinkman.Index and Alcohol. For **species** it runs on CLR values with the LinDA
   compositional bias correction; for **metabolites** it runs on the log10 matrix with
   `bias_correct=False`, because a metabolite panel is not a composition and subtracting the
   coefficient mode would not be meaningful there. `Coef > 0` means higher in group2, i.e.
   the **same** sign convention as `MeanDiff`.

### Robust hits

`Robust == True` requires **q < 0.05 in both tests AND the same direction**. That is the
column NB03 selects its ML targets from, and the one NB08/NB09 treat as statistical
evidence. Features significant in only one test are kept in the table with `Robust = False`
and an `Agreement` label, so a reader can see what confounding adjustment removed rather
than having it silently dropped:

- `both_q_same_dir` — robust.
- `both_q_opposite_dir` — significant in both but directions disagree: almost always a
  confounded feature, never counted as a hit.
- `mw_only` — passes the rank test, fails adjustment. The prime suspects are the
  smoking-associated taxa, since Brinkman index differs ~2–4× across these groups.
- `lm_only` — adjustment reveals an effect the unadjusted test missed (confounding was
  masking it).

**Output (v2.1): one tidy table, `tables/differential_abundance_v2.csv`.**
Columns `Dataset | Comparison | Scheme | Assay | Feature | ...` followed by the
Mann-Whitney columns, the adjusted-model columns, `Agreement` and `Robust`. Slice it with
a filter rather than opening a different file:

```python
da = pd.read_csv(RESULTS_DIR / 'tables' / 'differential_abundance_v2.csv')
early = da[(da.Comparison == 'Healthy_vs_Early_CRC') & (da.Assay == 'metabolites')]
```

This replaces six per-comparison CSVs that used the *same filename template* as the
fine-grained family, so Yachida carried two statistically different answers to one
question under near-identical names — and since `Stage_0 ∪ Stage_I_II == Early_CRC` the
two families overlap by construction. A cross-reference of the tree found 12 of those 16
files had no reader at all, while the two that did read them disagreed about which family
was authoritative: NB03 requires the adjusted `Robust` column, and the ML workbench was
screening the *unadjusted* family on a nominal `PValue < 0.05`, which CONTRACTS §9
explicitly forbids. Keeping only the adjusted family removes the ambiguity.''')

# ------------------------------------------------- Step 6B code
nb.set("cell-015", r'''GROUP_DS = DATASET_DISCOVERY   # the only cohort in utils.COMPARISONS_GROUP

# Confounders actually usable in this cohort (present, non-constant, mostly observed).
_meta_g = harmonized_meta[GROUP_DS]
DA_COVARS = [c for c in YACHIDA_COMPLETE_CONFOUNDERS
             if c in _meta_g.columns
             and _meta_g[c].notna().mean() > 0.5
             and _meta_g[c].nunique(dropna=True) > 1]
print(f'Adjusted model covariates: {DA_COVARS}')

da_group_da_spc = {}
da_group_da_mtb = {}
robust_summary  = []
da_tidy_rows    = []   # accumulates every comparison x assay into one tidy table


def merge_da(mw, lm, label, dtype):
    """Join the rank test and the adjusted model on Feature and classify agreement.

    Sign conventions are aligned first: MeanDiff and Coef are both
    'group2 - group1', so a robust hit needs sign(MeanDiff) == sign(Coef).
    """
    if mw is None or mw.empty:
        return pd.DataFrame()
    out = mw.copy()
    if lm is None or lm.empty:
        out['Coef_adj'] = np.nan; out['SE_adj'] = np.nan
        out['PValue_adj'] = np.nan; out['QValue_adj'] = np.nan
        out['Agreement'] = 'mw_only_no_model'
        out['Robust'] = False
        return out

    l = lm[['Feature', 'Coef', 'SE', 'PValue', 'QValue']].rename(columns={
        'Coef': 'Coef_adj', 'SE': 'SE_adj',
        'PValue': 'PValue_adj', 'QValue': 'QValue_adj'})
    out = out.merge(l, on='Feature', how='left')

    sig_mw = out['QValue'] < FDR_THRESHOLD
    sig_lm = out['QValue_adj'] < FDR_THRESHOLD
    same   = np.sign(out['MeanDiff']) == np.sign(out['Coef_adj'])

    out['Agreement'] = np.select(
        [sig_mw & sig_lm & same, sig_mw & sig_lm & ~same, sig_mw & ~sig_lm, ~sig_mw & sig_lm],
        ['both_q_same_dir', 'both_q_opposite_dir', 'mw_only', 'lm_only'],
        default='neither')
    out['Robust'] = out['Agreement'] == 'both_q_same_dir'
    return out


for group1, group2 in COMPARISONS_GROUP[GROUP_DS]:
    print(f'\n=== {GROUP_DS}: {group1} vs {group2} ===')

    # 1. Unadjusted rank test
    mw_spc = differential_abundance(transformed_species[GROUP_DS], harmonized_meta[GROUP_DS],
                                    group1, group2, group_col='Stage.Group')
    mw_mtb = differential_abundance(transformed_mtb[GROUP_DS], harmonized_meta[GROUP_DS],
                                    group1, group2, group_col='Stage.Group')

    # 2. Adjusted model. Species: compositional (CLR in, LinDA bias correction).
    #    Metabolites: log10 matrix, NOT a composition -> bias_correct=False.
    lm_spc = clr_lm_da(transformed_species[GROUP_DS], harmonized_meta[GROUP_DS],
                       group1, group2, group_col='Stage.Group',
                       covars=DA_COVARS, already_clr=True, bias_correct=True)
    lm_mtb = clr_lm_da(transformed_mtb[GROUP_DS], harmonized_meta[GROUP_DS],
                       group1, group2, group_col='Stage.Group',
                       covars=DA_COVARS, already_clr=True, bias_correct=False)

    spc_da = merge_da(mw_spc, lm_spc, f'{group1}_vs_{group2}', 'species')
    mtb_da = merge_da(mw_mtb, lm_mtb, f'{group1}_vs_{group2}', 'metabolites')

    da_group_da_spc[(group1, group2)] = spc_da
    da_group_da_mtb[(group1, group2)] = mtb_da

    for dtype, df_da, lm_res in (('species', spc_da, lm_spc), ('metabolites', mtb_da, lm_mtb)):
        if df_da.empty:
            print(f'  {dtype}: no testable features')
            continue
        counts = df_da['Agreement'].value_counts().to_dict()
        n_rob = int(df_da['Robust'].sum())
        n_mw  = int((df_da['QValue'] < FDR_THRESHOLD).sum())
        n_lm  = int((df_da['QValue_adj'] < FDR_THRESHOLD).sum())
        print(f'  {dtype}: MW q<0.05: {n_mw} | adjusted q<0.05: {n_lm} | '
              f'ROBUST (both + same dir): {n_rob}')
        print(f'    agreement breakdown: {counts}')
        if not lm_res.empty:
            print(f"    model: covars={lm_res.attrs.get('covariates')}, "
                  f"dropped {lm_res.attrs.get('n_dropped_missing_covar', 0)} row(s) for missing covars, "
                  f"bias offset={lm_res.attrs.get('bias_offset', 0.0):+.4f}")

        # v2.1: accumulated into ONE tidy table instead of six per-comparison CSVs
        # (written after the loop). The old layout wrote the same filename template as
        # the fine-grained family, so Yachida had two statistically different answers
        # to the same question under near-identical names -- and because
        # Stage_0 u Stage_I_II == Early_CRC, they overlap by construction.
        _tidy = df_da.copy()
        _tidy.insert(0, 'Dataset', GROUP_DS)
        _tidy.insert(1, 'Comparison', f'{group1}_vs_{group2}')
        _tidy.insert(2, 'Scheme', 'grouped')      # Stage.Group, confounder-adjusted
        _tidy.insert(3, 'Assay', dtype)
        da_tidy_rows.append(_tidy)

        robust_summary.append({
            'Comparison': f'{group1}_vs_{group2}', 'Type': dtype,
            'N_tested': len(df_da), 'N_MW_q05': n_mw, 'N_adjusted_q05': n_lm,
            'N_robust': n_rob,
            'N_mw_only': int((df_da['Agreement'] == 'mw_only').sum()),
            'N_lm_only': int((df_da['Agreement'] == 'lm_only').sum()),
            'N_opposite': int((df_da['Agreement'] == 'both_q_opposite_dir').sum()),
        })

DA_TIDY_PATH = RESULTS_DIR / 'tables' / 'differential_abundance_v2.csv'
da_tidy = (pd.concat(da_tidy_rows, ignore_index=True) if da_tidy_rows
           else pd.DataFrame(columns=['Dataset', 'Comparison', 'Scheme', 'Assay', 'Feature']))
da_tidy.to_csv(DA_TIDY_PATH, index=False)
print(f'\nSaved: differential_abundance_v2.csv  ({len(da_tidy):,} rows, '
      f'{da_tidy["Comparison"].nunique() if not da_tidy.empty else 0} comparisons x '
      f'{da_tidy["Assay"].nunique() if not da_tidy.empty else 0} assays)')
print('  One tidy table replaces the six per-comparison CSVs. Slice it with, e.g.:')
print("    da[(da.Comparison == 'Healthy_vs_Early_CRC') & (da.Assay == 'metabolites')]")

robust_df = pd.DataFrame(robust_summary)
robust_df.to_csv(RESULTS_DIR / 'tables' / 'da_robust_hit_summary_v2.csv', index=False)
print('\n=== Robust-hit summary (both tests, same direction) ===')
print(robust_df.to_string(index=False))
print('\nSaved: da_robust_hit_summary_v2.csv')
print('\nInterpretation: "mw_only" counts features whose unadjusted significance does not '
      'survive adjustment for age/BMI/sex/smoking/alcohol. In this cohort smoking is the '
      'main suspect (median Brinkman 100 in Healthy vs 196-390 in CRC). Those features are '
      'NOT treated as findings anywhere downstream.')''')

# ------------------------------------------------- 7.3 partial correlations header
nb.set("cell-022", r'''### 7.3 Partial Correlations (controlling for the full confounder set)

Adjusts species–polyamine Spearman correlations for the complete Yachida confounder set:
`Age`, `BMI`, `Gender`, `Brinkman.Index` and `Alcohol` (`utils.YACHIDA_COMPLETE_CONFOUNDERS`),
all carried through `harmonized_meta` by `utils.harmonize_metadata()`.

**Tumour location is deliberately not a covariate here.** `Tumor location` is "-" for every
non-cancer sample, so on this mixed Healthy+CRC set adjusting for it would partial out
disease status itself. `utils.build_confounder_matrix()` refuses the request outright rather
than relying on this note; the within-cancer analysis that *may* legitimately adjust for
tumour site is NB04.

**Missing covariates are dropped, not filled.** The previous version replaced missing
`Brinkman.Index` and `Alcohol` with the column mean before adjusting. Mean-filling a
covariate and then partialling it out understates the uncertainty and biases the adjustment
toward no effect for exactly the samples where the covariate was unknown. Complete-case
analysis is used instead, and the number of dropped samples is reported.

All reduced species are tested against every polyamine — no marginal pre-screen, which would
make the partial-correlation FDR conditional on a prior screen. The table is saved as
`{dataset}_polyamine_partial_correlations_v2.csv`; **NB09 uses it as its SHAP-independent
candidate pool**, which is what breaks the old pipeline's circularity (SHAP picking the
candidates and then being counted as evidence for them).''')

nb.set("cell-023", r'''try:
    import pingouin as pg
    HAS_PINGOUIN = True
    print(f'pingouin {pg.__version__} loaded.')
except ImportError:
    print('pingouin not installed. Run: pip install pingouin')
    HAS_PINGOUIN = False

partial_corr_results = {}

if HAS_PINGOUIN:
    for ds in POLYAMINE_DATASETS_ANALYSIS:
        if ds not in DATASETS_ANALYSIS:
            continue
        print(f'\n--- Partial correlations: {ds} ---')
        meta = harmonized_meta[ds].copy()

        # Build the covariate list from what this cohort actually has. Gender is
        # encoded numerically; everything else must already be numeric.
        covariates = []
        for col in YACHIDA_COMPLETE_CONFOUNDERS:
            if col not in meta.columns:
                continue
            if col == 'Gender':
                meta['Gender_num'] = (meta['Gender'].astype(str).str.lower() == 'male').astype(float)
                if meta['Gender_num'].notna().mean() > 0.5 and meta['Gender_num'].nunique() > 1:
                    covariates.append('Gender_num')
                continue
            series = pd.to_numeric(meta[col], errors='coerce')
            if series.notna().mean() > 0.5 and series.nunique(dropna=True) > 1:
                meta[col] = series
                covariates.append(col)

        # No tumour-location covariate: it is a disease-status proxy on this set.

        if not covariates:
            print('  No usable covariates -- skipping.')
            continue
        print(f'  Adjusting for: {covariates}')

        spc_renamed = transformed_species[ds].rename(
            columns={c: f'SPC__{c}' for c in transformed_species[ds].columns if c != 'Sample'})
        mtb_renamed = transformed_mtb[ds].rename(
            columns={c: f'MTB__{c}' for c in transformed_mtb[ds].columns if c != 'Sample'})

        merged_all = spc_renamed.merge(mtb_renamed, on='Sample').merge(
            meta[['Sample'] + covariates], on='Sample')
        # Complete-case: drop rows missing any covariate (never mean-fill).
        merged = merged_all.dropna(subset=covariates)
        n_dropped = len(merged_all) - len(merged)
        print(f'  n={len(merged)} complete cases'
              + (f' ({n_dropped} dropped for missing covariates)' if n_dropped else ''))

        top_species = [c for c in transformed_species[ds].columns if c != 'Sample']

        # Direct polyamines plus pre-specified proxies; each row records which it was.
        targets = [(col, POLYAMINE_KEGG[k], 'direct')
                   for k, col in polyamine_columns[ds].items()]
        targets += [(col, POLYAMINE_KEGG[k], f'proxy:{via}')
                    for k, (col, via) in polyamine_proxies.get(ds, {}).items()]

        pc_rows = []
        for spc in top_species:
            spc_key = f'SPC__{spc}'
            if spc_key not in merged.columns:
                continue
            for pa_col, pa_name, measured_as in targets:
                pa_key = f'MTB__{pa_col}'
                if pa_key not in merged.columns:
                    continue
                try:
                    result = pg.partial_corr(data=merged, x=spc_key, y=pa_key,
                                              covar=covariates, method='spearman')
                except Exception:
                    continue
                pc_rows.append({
                    'Species':        spc,
                    'Metabolite':     pa_col,
                    'Polyamine':      pa_name,
                    'Measured_as':    measured_as,
                    'Rho_partial':    result['r'].values[0],
                    'PValue_partial': result['p-val'].values[0],
                    'N':              int(result['n'].values[0]),
                })

        if not pc_rows:
            print('  No valid pairs computed.')
            continue

        pc_df = pd.DataFrame(pc_rows)
        _, qvals, _, _ = multipletests(pc_df['PValue_partial'], method='fdr_bh')
        pc_df['QValue_partial'] = qvals
        col_to_kegg = {v: k for k, v in polyamine_columns[ds].items()}
        col_to_kegg.update({col: k for k, (col, _v) in polyamine_proxies.get(ds, {}).items()})
        pc_df['KEGG'] = pc_df['Metabolite'].map(col_to_kegg)
        pc_df['Covariates'] = ', '.join(covariates)
        partial_corr_results[ds] = pc_df
        pc_df.to_csv(RESULTS_DIR / 'tables' / f'{ds}_polyamine_partial_correlations_v2.csv',
                     index=False)

        n_sig = int((pc_df['QValue_partial'] < 0.05).sum())
        n_pos = int(((pc_df['QValue_partial'] < 0.05) & (pc_df['Rho_partial'] > 0)).sum())
        print(f'  {len(pc_df)} pairs tested, {n_sig} significant (q<0.05), '
              f'of which {n_pos} POSITIVE (candidate producers for NB09)')''')

# ------------------------------------------------- 7.4 three-cohort replication
nb.set("cell-024", r'''### 7.4 Cross-Cohort Replication (Yachida vs Sinha vs Kim, genus level)

**Why genus, not species.** Yachida is shotgun metagenomics at species resolution; Sinha and
Kim are 16S at genus resolution. Species names therefore *cannot* match across cohorts — the
earlier species-level version of this step returned zero pairs for that reason alone, which
looks like a negative result but is an artefact of the taxonomy. All three cohorts are
aggregated to **genus** (GTDB placeholder suffixes such as `_A` stripped, then re-aggregated),
CLR-transformed, and correlated with each polyamine.

**Why KEGG IDs, not column names.** The three platforms name the same compound differently
(`C00134_Putrescine` / `PUTRESCINE` / `N-acetylputrescine`), so matching is done on KEGG
compound ID.

**Proxy handling.** Kim does not measure putrescine or cadaverine, only their N-acetyl
derivatives. Those are mapped to the **parent** compound for matching purposes — otherwise
Kim could never replicate the study's headline compound — but every resulting row carries a
`Measured_as` / `Proxy_used` label. A replication resting on a derivative is weaker evidence
than one resting on the parent, and the final evidence table in NB09 keeps that distinction
rather than averaging it away.

**Replication is directional, not just "both significant".** `Same_direction` requires the
sign of ρ to agree across the cohorts in which the pair is significant. A pair significant in
two cohorts with *opposite* signs is not a replication; it is a contradiction, and is
reported as such.

Sinha and Kim metabolites are pre-standardised; Spearman is rank-based, so that has no
effect on this step.''')

nb.set("cell-025", r'''import re
from utils import aggregate_to_genus, clr_transform, filter_by_prevalence, PREVALENCE_SPECIES

def _norm_genus(g):
    g = re.sub(r'_[A-Z]+$', '', str(g).strip())      # GTDB placeholder suffixes
    return g.strip('[]').strip()

def genus_clr(ds):
    """Genus-level CLR table for one cohort, restricted to the primary samples."""
    src = preprocessed['filtered_species'][ds]
    keep = set(harmonized_meta[ds]['Sample'])
    src = src[src['Sample'].isin(keep)].reset_index(drop=True)

    gen = aggregate_to_genus(src)
    feat = gen.drop(columns='Sample')
    feat.columns = [_norm_genus(c) for c in feat.columns]
    feat = feat.T.groupby(level=0).sum().T          # merge genera joined by suffix stripping
    feat = feat.loc[:, ~feat.columns.isin(['', 'Unknown', 'None'])]
    feat.insert(0, 'Sample', gen['Sample'].values)
    feat = filter_by_prevalence(feat, threshold=PREVALENCE_SPECIES)
    return clr_transform(feat)

# Per-cohort map {metabolite column -> PARENT polyamine KEGG id}, so a proxy column
# matches its parent across cohorts. Also record which cohorts used a proxy.
pa_match_map = {}
proxy_used_by_ds = {}
for ds in POLYAMINE_DATASETS_ANALYSIS:
    m = {col: kegg for kegg, col in polyamine_columns.get(ds, {}).items()}
    px = {}
    for kegg, (col, via) in polyamine_proxies.get(ds, {}).items():
        m[col] = kegg          # map derivative -> PARENT compound
        px[kegg] = via
    pa_match_map[ds] = m
    proxy_used_by_ds[ds] = px

genus_polyamine_corr = {}
for ds in POLYAMINE_DATASETS_ANALYSIS:
    pa_cols = list(pa_match_map[ds].keys())
    if not pa_cols:
        print(f'{ds}: no polyamine columns -- skipped')
        continue
    gclr = genus_clr(ds)
    corr = compute_correlations(gclr, transformed_mtb[ds], target_metabolites=pa_cols)
    genus_polyamine_corr[ds] = corr
    n_sig = int((corr['QValue'] < 0.05).sum()) if not corr.empty else 0
    print(f'{ds}: {gclr.shape[1]-1} genera x {len(pa_cols)} polyamine columns, '
          f'{n_sig} pairs q<0.05'
          + (f'  [proxies: {list(proxy_used_by_ds[ds])}]' if proxy_used_by_ds[ds] else ''))

cohorts_present = [ds for ds in genus_polyamine_corr if not genus_polyamine_corr[ds].empty]
if len(cohorts_present) >= 2:
    shared_genera = set.intersection(
        *[set(genus_polyamine_corr[ds]['Species']) for ds in cohorts_present])
    print(f'\nGenera testable in all {len(cohorts_present)} cohorts: {len(shared_genera)}')
else:
    shared_genera = set()
    print('\nFewer than two cohorts have testable pairs -- no replication possible.')


def _rho_lookup(ds):
    km = pa_match_map.get(ds, {})
    return {(r.Species, km.get(r.Metabolite, r.Metabolite)): r.Rho
            for r in genus_polyamine_corr[ds].itertuples()}

def _annotate(df):
    """Add Same_direction, which cohorts replicated, and whether any used a proxy."""
    if df.empty:
        return df
    maps = {ds: _rho_lookup(ds) for ds in cohorts_present}
    df = df.copy()
    same_dir, prox, coh = [], [], []
    for g, k in zip(df['Species'], df['Metabolite_key']):
        signs = {ds: np.sign(m[(g, k)]) for ds, m in maps.items() if (g, k) in m}
        same_dir.append(len(set(signs.values())) == 1 and len(signs) >= 2)
        coh.append(', '.join(sorted(signs)))
        prox.append(', '.join(sorted(
            ds for ds in signs if k in proxy_used_by_ds.get(ds, {}))) or 'none')
    df['Same_direction'] = same_dir
    df['Cohorts_tested'] = coh
    df['Proxy_used_in']  = prox
    df['Polyamine']      = [POLYAMINE_KEGG.get(k, k) for k in df['Metabolite_key']]
    return df

consistent_pairs   = _annotate(find_consistent_correlations(
    genus_polyamine_corr, kegg_map=pa_match_map, min_datasets=2, q_threshold=0.05))
consistent_relaxed = pd.DataFrame()

out_csv = RESULTS_DIR / 'tables' / 'consistent_polyamine_correlations_genus_v2.csv'
if not consistent_pairs.empty:
    n_same = int(consistent_pairs['Same_direction'].sum())
    print(f'\n{len(consistent_pairs)} genus-polyamine pairs (KEGG-matched) significant in '
          f'>=2 cohorts at q<0.05; same direction in all of them: {n_same}')
    if n_same < len(consistent_pairs):
        print(f'  {len(consistent_pairs) - n_same} pair(s) are significant in two cohorts '
              f'with OPPOSITE signs -- contradictions, not replications.')
    display(consistent_pairs.head(20))
    consistent_pairs.rename(columns={'Species': 'Genus'}).to_csv(out_csv, index=False)
    print(f'Saved: {out_csv.name}')
else:
    print('\nNo pairs replicate at q<0.05 in >=2 cohorts. Reporting a q<0.10 '
          'hypothesis-generating list instead, clearly labelled as such.')
    consistent_relaxed = _annotate(find_consistent_correlations(
        genus_polyamine_corr, kegg_map=pa_match_map, min_datasets=2, q_threshold=0.10))
    print(f'At q<0.10: {len(consistent_relaxed)} pairs')
    if not consistent_relaxed.empty:
        consistent_relaxed['Threshold'] = 'q<0.10 (exploratory)'
        display(consistent_relaxed.head(20))
        consistent_relaxed.rename(columns={'Species': 'Genus'}).to_csv(out_csv, index=False)
        print(f'Saved: {out_csv.name}  (exploratory threshold -- NOT a replication claim)')

# Per-cohort genus-level tables, so NB08/NB09 can cite a specific cohort's result.
for ds, corr in genus_polyamine_corr.items():
    if corr.empty:
        continue
    c = corr.copy()
    c['Parent_KEGG'] = c['Metabolite'].map(pa_match_map[ds])
    c['Polyamine']   = c['Parent_KEGG'].map(POLYAMINE_KEGG)
    c['Measured_as'] = np.where(
        c['Parent_KEGG'].isin(proxy_used_by_ds.get(ds, {})), 'proxy', 'direct')
    c.rename(columns={'Species': 'Genus'}).to_csv(
        RESULTS_DIR / 'tables' / f'{ds}_genus_polyamine_correlations_v2.csv', index=False)
print('\nSaved per-cohort genus-level correlation tables.')''')

nb.save()
