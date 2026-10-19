"""NB02 edits, part C: replication plot, stratified networks, dashboard,
progression scope, final save/summary."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from nbedit import Notebook

nb = Notebook(Path(__file__).resolve().parent / "02_analysis.ipynb")

# ------------------------------------------------- replication dot plot
nb.set("cell-026", r'''# Cross-cohort dot plot. consistent_pairs carries a KEGG-matched 'Metabolite_key'
# (a PARENT compound id), so for each cohort we look up which of ITS OWN columns map
# to that parent via pa_match_map -- which is also how a proxy column is found.

def _lookup_corr(ds, species, metabolite_key):
    corr_df = genus_polyamine_corr.get(ds)
    if corr_df is None or corr_df.empty:
        return None
    ds_map = pa_match_map.get(ds, {})
    candidate_cols = [col for col, kid in ds_map.items() if kid == metabolite_key]
    if not candidate_cols:
        return None
    match = corr_df[(corr_df['Species'] == species) &
                    (corr_df['Metabolite'].isin(candidate_cols))]
    return None if match.empty else match.iloc[0]


plot_consistent = consistent_pairs if not consistent_pairs.empty else consistent_relaxed

if not plot_consistent.empty:
    top_pairs = plot_consistent.head(15)
    plot_data = []
    for _, pair in top_pairs.iterrows():
        for ds in POLYAMINE_DATASETS_ANALYSIS:
            if ds not in genus_polyamine_corr:
                continue
            row = _lookup_corr(ds, pair['Species'], pair['Metabolite_key'])
            if row is None:
                continue
            is_proxy = pair['Metabolite_key'] in proxy_used_by_ds.get(ds, {})
            plot_data.append({
                'Pair':        f"{pair['Species'][:24]}\n<-> {pair.get('Polyamine', pair['Metabolite_key'])}",
                'Cohort':      ds.split('-')[0].split('_')[0],
                'Rho':         row['Rho'],
                'Significant': row['QValue'] < 0.05,
                'Proxy':       is_proxy,
            })

    if plot_data:
        plot_df = pd.DataFrame(plot_data)
        cohort_colors = {'YACHIDA': '#2a78d6', 'SINHA': '#eb6834', 'Kim': '#2e9e5b'}
        fig, ax = plt.subplots(figsize=(13, max(6, len(top_pairs) * 0.62)))

        for cohort, sub_c in plot_df.groupby('Cohort'):
            col = cohort_colors.get(cohort, '#777777')
            for sig, sub in sub_c.groupby('Significant'):
                # Proxy measurements get a hollow marker: visibly weaker evidence.
                for proxy, s in sub.groupby('Proxy'):
                    ax.scatter(
                        s['Rho'], s['Pair'],
                        facecolors='none' if proxy else col,
                        edgecolors=col, linewidths=1.6,
                        marker='o' if sig else 'x',
                        s=110 if sig else 55, alpha=0.9,
                        label=f"{cohort}{' (proxy)' if proxy else ''}"
                              f"{'' if sig else ' ns'}")

        ax.axvline(0, color='grey', ls='--', alpha=0.5)
        ax.set_xlabel('Spearman rho (genus vs polyamine, CLR)')
        handles, labels = ax.get_legend_handles_labels()
        seen = dict(zip(labels, handles))
        ax.legend(seen.values(), seen.keys(), fontsize=7, loc='best', framealpha=0.9)
        ax.set_title('Cross-cohort replication: genus-polyamine correlations\n'
                     'genus level, matched by KEGG parent compound | filled=q<0.05, '
                     'hollow=measured via N-acetyl proxy, x=not significant',
                     fontsize=10, fontweight='bold')
        plt.tight_layout()
        plt.savefig(RESULTS_DIR / 'figures' / 'correlation' / 'cross_cohort_consistency_v2.png',
                    dpi=150, bbox_inches='tight')
        plt.show(); plt.close()
    else:
        print('No plottable cross-cohort pairs found.')
else:
    print('No consistent pairs (strict or exploratory) to plot.')''')

# ------------------------------------------------- stratified hub networks
nb.set("cell-035", r'''---
## Stratified Hub Species Network Analysis

Separate polyamine-centric correlation networks per CRC stage group
(`YACHIDA-CRC-2019`, the only cohort with staging), using `utils.STAGE_GROUP_ORDER`:
**Healthy, Early_CRC, Advanced_CRC**. MP is excluded from the study and PostSurgery is not
on the severity axis, so neither appears here.

**Caveat — this comparison is confounded by sample size and must not be read as biology.**
The groups differ substantially in n (127 / 96 / 54) while every network uses the same
q<0.05 and |rho|>=0.3 thresholds. Statistical power scales with n, so the Healthy network
will contain more edges than the Advanced_CRC network *by construction*, whether or not its
true correlation structure is any denser. Edge counts across groups are therefore reported
for completeness but are not interpreted as "the network densifies/sparsifies with stage";
no downstream notebook consumes these counts as evidence.''')

nb.set("cell-036", r'''CRC_DS_HUB = DATASET_DISCOVERY

# Stage.Group comes from utils.harmonize_metadata() (NB01) via utils.STAGE_GROUP_MAP --
# nothing is re-derived locally.
print('Primary stage groups (utils.STAGE_GROUP_ORDER):', STAGE_GROUP_ORDER)
print('utils.STAGE_GROUP_MAP (fine-grained Study.Group -> grouped Stage.Group):')
for sg, grp in STAGE_GROUP_MAP.items():
    note = '   [sensitivity: not analysed here]' if grp not in STAGE_GROUP_ORDER else ''
    print(f'  {sg:14s} -> {grp}{note}')


def filter_by_group(ds, group_name):
    """Species and metabolite frames for samples whose Stage.Group == group_name."""
    meta    = harmonized_meta[ds]
    samples = meta[meta['Stage.Group'] == group_name]['Sample'].values
    spc_sub = transformed_species[ds][
        transformed_species[ds]['Sample'].isin(samples)].reset_index(drop=True)
    mtb_sub = transformed_mtb[ds][
        transformed_mtb[ds]['Sample'].isin(samples)].reset_index(drop=True)
    return spc_sub, mtb_sub


_counts = (harmonized_meta[CRC_DS_HUB].groupby('Stage.Group')['Sample'].count()
           .reindex(STAGE_GROUP_ORDER))
print('\nSamples per Stage.Group (primary set only):')
print(_counts.to_string())
print(f'\nTotal primary: {int(_counts.sum())}')
print('Group sizes differ, so edge counts below are power-limited -- see the caveat above.')''')

nb.set("cell-040", r'''# Non-Healthy multi-panel hub species networks (Early CRC | Advanced CRC)
from matplotlib.lines import Line2D

NON_HEALTHY = [g for g in STAGE_GROUP_ORDER if g != 'Healthy']   # Early_CRC, Advanced_CRC
PANEL_TITLES = {
    'Early_CRC':    'Early CRC (Stage 0 + I/II)',
    'Advanced_CRC': 'Advanced CRC (Stage III/IV)',
}

fig, axes = plt.subplots(1, len(NON_HEALTHY), figsize=(14 * len(NON_HEALTHY), 14),
                          squeeze=False)
axes = axes.flatten()
fig.suptitle(
    'Non-Healthy Hub Species Networks by CRC Stage Group\n'
    '(YACHIDA-CRC-2019, primary set; MP excluded from the study, PostSurgery off-axis)\n'
    'Edge counts are power-limited by group size and are not interpreted as biology',
    fontsize=15, fontweight='bold', y=1.02
)

for ax, group_name in zip(axes, NON_HEALTHY):
    G_g   = stratified_networks.get(group_name, nx.Graph())
    n_grp = int((harmonized_meta[CRC_DS_HUB]['Stage.Group'] == group_name).sum())
    title = f"{PANEL_TITLES[group_name]}  (n={n_grp})"

    if G_g.number_of_nodes() == 0:
        ax.text(0.5, 0.5, 'No significant\ncorrelations found',
                ha='center', va='center', transform=ax.transAxes, fontsize=13)
        ax.set_title(title, fontsize=14, fontweight='bold')
        ax.axis('off')
        continue

    pos = nx.spring_layout(G_g, k=2.0, seed=42)

    node_colors, node_sizes = [], []
    for n in G_g.nodes():
        if n in pa_node_set_hub:
            node_colors.append('#e74c3c'); node_sizes.append(600)
        elif G_g.nodes[n].get('node_type') == 'species':
            node_colors.append('#3498db'); node_sizes.append(220)
        else:
            node_colors.append('#2ecc71'); node_sizes.append(160)

    edge_colors = ['#e74c3c' if G_g[u][v].get('sign') == 'negative' else '#27ae60'
                   for u, v in G_g.edges()]
    nx.draw(G_g, pos, node_color=node_colors, node_size=node_sizes,
            edge_color=edge_colors, alpha=0.7, ax=ax, with_labels=False, width=0.8)

    labels = {}
    for n in G_g.nodes():
        if n in pa_node_set_hub:
            labels[n] = (n.split('_', 1)[-1] if '_' in n else n)[:28]
        elif G_g.nodes[n].get('node_type') == 'species' and G_g.degree(n) >= 3:
            labels[n] = n[:26]
    nx.draw_networkx_labels(G_g, pos, labels, font_size=6, ax=ax)

    n_spc = sum(1 for _, d in G_g.nodes(data=True) if d.get('node_type') == 'species')
    ax.set_title(f'{title}\n{G_g.number_of_nodes()} nodes ({n_spc} spc) | '
                 f'{G_g.number_of_edges()} edges', fontsize=13, fontweight='bold')

legend_elements = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#e74c3c', markersize=14,
           label='Polyamine'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#3498db', markersize=14,
           label='Species'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#2ecc71', markersize=12,
           label='Other metabolite'),
    Line2D([0], [0], color='#e74c3c', lw=2.5, label='Negative correlation'),
    Line2D([0], [0], color='#27ae60', lw=2.5, label='Positive correlation'),
]
fig.legend(handles=legend_elements, loc='lower center', ncol=5, fontsize=12,
           bbox_to_anchor=(0.5, -0.04), framealpha=0.9)

plt.tight_layout()
out_multi = RESULTS_DIR / 'figures' / 'network' / 'hub_network_non_healthy_multifig_v2.png'
plt.savefig(out_multi, dpi=150, bbox_inches='tight')
plt.show(); plt.close()
print(f'Saved: {out_multi}')''')

# ------------------------------------------------- Step 10 header
nb.set("cell-041", r'''---
## Step 10: CRC Progression Analysis

**Scope: Yachida only** (Sinha and Kim are binary Healthy/CRC cohorts with no staging), and
within Yachida the **primary** set only.

Two complementary approaches:
- **10A — ordinal Spearman**: correlate species and metabolites with ordered stage
  (`utils.CRC_STAGE_NUMERIC`: Healthy=0, Stage_0=1, Stage_I_II=2, Stage_III_IV=3).
  `HS`/`PostSurgery` is **NaN** on this scale by construction, so it cannot enter the trend
  test even if a caller forgets to filter — and it has already been removed from these
  matrices anyway. MP does not exist in the data.
- **10B — binary early vs late**: `Stage_I_II` vs `Stage_III_IV`, excluding Healthy and
  Stage_0, giving a direct test of what separates early from advanced disease.''')

nb.set("cell-041b", r'''### 8.5 Stage-group dashboard

Cohort composition and hub-network size across the three primary `Stage.Group` tiers, in
clinical order (`Healthy -> Early_CRC -> Advanced_CRC`). The left panel uses a single-hue
sequential ramp stepped by severity rather than unrelated categorical colours, because these
are an *ordered* progression; the darkening blue reinforces the same ordering
`utils.STAGE_GROUP_ORDER` encodes everywhere else.

The right panel repeats the caveat that matters: network size tracks group size here, so the
two series are shown together with n annotated and are not read as a biological trend.''')

nb.set("cell-041c", r'''# Stage-group dashboard: sample composition (left) + hub-network size (right).
# Sequential blue ramp for the ordered severity tiers (never unrelated categorical
# hues for what is one ordinal axis); categorical blue/orange for the two
# independent network-size series on the right.
STAGE_RAMP = ['#9ec5f4', '#2a78d6', '#104281']          # one step per primary tier
SERIES_BLUE, SERIES_ORANGE = '#2a78d6', '#eb6834'
INK, INK_SECONDARY, MUTED, GRID = '#0b0b0b', '#52514e', '#898781', '#e1e0d9'

_counts_dash = (harmonized_meta[CRC_DS_HUB]
                .groupby('Stage.Group')['Sample'].count()
                .reindex(STAGE_GROUP_ORDER, fill_value=0))
_net_sizes = [
    (g, stratified_networks.get(g, nx.Graph()).number_of_nodes(),
        stratified_networks.get(g, nx.Graph()).number_of_edges())
    for g in STAGE_GROUP_ORDER
]

fig, (ax_n, ax_g) = plt.subplots(1, 2, figsize=(13, 4.6))
fig.patch.set_facecolor('#fcfcfb')

# -- left: sample N per stage group, sequential-by-severity --
y_pos = range(len(STAGE_GROUP_ORDER))
ax_n.barh(y_pos, _counts_dash.values, color=STAGE_RAMP, height=0.62,
          edgecolor='white', linewidth=0.6)
for i, v in enumerate(_counts_dash.values):
    ax_n.text(v + max(_counts_dash.values) * 0.015, i, f'n={v}',
              va='center', fontsize=10, color=INK)
ax_n.set_yticks(list(y_pos))
ax_n.set_yticklabels(STAGE_GROUP_ORDER, fontsize=10.5, color=INK)
ax_n.invert_yaxis()
ax_n.set_xlabel('Samples', fontsize=10.5, color=INK_SECONDARY)
ax_n.set_title(f'Cohort composition (primary set, n={int(_counts_dash.sum())})',
               fontsize=12, fontweight='bold', color=INK, loc='left')
for spine in ['top', 'right']:
    ax_n.spines[spine].set_visible(False)
for spine in ['left', 'bottom']:
    ax_n.spines[spine].set_color(GRID)
ax_n.tick_params(colors=MUTED)
ax_n.set_facecolor('#fcfcfb')

# -- right: hub-network nodes vs edges per stage group --
labels_g = [g for g, _, _ in _net_sizes]
nodes_g  = [n for _, n, _ in _net_sizes]
edges_g  = [e for _, _, e in _net_sizes]
x = range(len(labels_g))
w = 0.36
ax_g.bar([i - w/2 for i in x], nodes_g, width=w, color=SERIES_BLUE,
         label='Nodes', edgecolor='white', linewidth=0.6)
ax_g.bar([i + w/2 for i in x], edges_g, width=w, color=SERIES_ORANGE,
         label='Edges', edgecolor='white', linewidth=0.6)
_ymax = max(edges_g + nodes_g + [1])
for i, (n, e) in enumerate(zip(nodes_g, edges_g)):
    if n: ax_g.text(i - w/2, n + _ymax * 0.02, str(n), ha='center', fontsize=9, color=INK_SECONDARY)
    if e: ax_g.text(i + w/2, e + _ymax * 0.02, str(e), ha='center', fontsize=9, color=INK_SECONDARY)
ax_g.set_xticks(list(x))
ax_g.set_xticklabels([f'{g}\n(n={_counts_dash[g]})' for g in labels_g],
                     fontsize=9, color=INK)
ax_g.set_ylabel('Count', fontsize=10.5, color=INK_SECONDARY)
ax_g.set_title('Hub-network size — tracks group size, not biology',
               fontsize=12, fontweight='bold', color=INK, loc='left')
ax_g.legend(frameon=False, fontsize=9.5, loc='upper right')
for spine in ['top', 'right']:
    ax_g.spines[spine].set_visible(False)
for spine in ['left', 'bottom']:
    ax_g.spines[spine].set_color(GRID)
ax_g.tick_params(colors=MUTED)
ax_g.set_facecolor('#fcfcfb')
ax_g.grid(axis='y', color=GRID, linewidth=0.7, zorder=0)
ax_g.set_axisbelow(True)

plt.tight_layout()
out_dash = RESULTS_DIR / 'figures' / 'network' / 'stage_group_dashboard_v2.png'
plt.savefig(out_dash, dpi=150, bbox_inches='tight', facecolor=fig.get_facecolor())
plt.show(); plt.close()
print(f'Saved: {out_dash}')''')

nb.set("cell-042", r'''CRC_DS = DATASET_DISCOVERY

print('Yachida stage distribution (primary set):')
print(harmonized_meta[CRC_DS]['Study.Group'].value_counts()
      .reindex(CRC_STAGE_ORDER, fill_value=0).to_string())

print('\nOrdinal stage coding actually in use (utils.CRC_STAGE_NUMERIC):')
print(harmonized_meta[CRC_DS][['Study.Group', 'CRC.Stage.Numeric']]
      .drop_duplicates().sort_values('CRC.Stage.Numeric').to_string(index=False))

# The trend tests below depend on this: every primary sample must have a finite
# ordinal position, and the sensitivity group must be absent.
assert harmonized_meta[CRC_DS]['CRC.Stage.Numeric'].notna().all(), \
    'a primary sample has no ordinal stage'
assert not harmonized_meta[CRC_DS]['Stage.Group'].eq('PostSurgery').any()
print(f'\nAll {len(harmonized_meta[CRC_DS])} primary samples carry a finite ordinal stage.')''')

nb.set("cell-048", r'''### 10B: Binary Early vs Late CRC

Comparison: `Stage_I_II` vs `Stage_III_IV`. Healthy and Stage_0 samples are excluded, so
this is a within-cancer contrast: what separates early-invasive from advanced disease.
MP is not in the study and PostSurgery is not in these matrices.

Note this is a *different* contrast from Step 6B's `Early_CRC vs Advanced_CRC`, which pools
Stage_0 into the early arm. Both are reported because they answer slightly different
questions, and the stage-group version is the one NB03 draws targets from.

Tumour location is a legitimate covariate for a contrast like this one (both arms are
carcinomas) — that analysis is NB04's, which adjusts for site and reports
rectum/left/right subgroups. It is not repeated here.''')

# ------------------------------------------------- save + summary
nb.set("cell-054", r'''---
## Results are already on disk as CSVs — no results pickle

`analysis_results_v2.pkl` was **removed in v2.1**. A full read-vs-write cross-reference of
the tree found **no reader anywhere** for it, and every meaningful key was already written
as a CSV by the cells above:

| former pickle key | the CSV that already holds it |
|---|---|
| `da_group_da_spc` / `da_group_da_mtb` | `differential_abundance_v2.csv` (tidy, adjusted) |
| `da_results` | the fine-grained, unadjusted DA family (see note below) |
| `partial_corr` | `{ds}_polyamine_partial_correlations_v2.csv` |
| `genus_polyamine_corr` | `{ds}_genus_polyamine_correlations_v2.csv` |
| `consistent_pairs` + `consistent_threshold` | `consistent_polyamine_correlations_genus_v2.csv` |
| `robust_summary` | `da_robust_hit_summary_v2.csv` |
| `prog_corr_ordinal` | `CRC_progression_ordinal_correlations_v2.csv` |
| `network_metrics` | `polyamine_hub_species_v2.csv`, `stratified_hub_species_v2.csv` |

The remaining keys (`full_corr`, `polyamine_corr`, `pa_match_map`, `proxy_used_by_ds`,
`da_covariates`, `sample_scope`) had no reader *and* no CSV — they were written and never
looked at. `full_corr` in particular is every species × every metabolite; NB01's own notes
record an earlier pipeline where that same object grew the pickle to ~800 MB with nothing
downstream reading it, and it had been reintroduced here under a new filename.

Anything genuinely needed downstream should be added as a named CSV with a stated reader,
not as a key in a catch-all pickle.''')

nb.set("cell-055", r'''# No results pickle: see the markdown above. This cell deliberately writes nothing and
# only reports where the results actually live, so a reader of the executed notebook is
# not left looking for a file that no longer exists.
print('NB02 results are on disk as CSVs under tables/ -- no analysis_results_v2.pkl.')
print('Removed in v2.1: it had zero readers and duplicated the CSVs written above.')''')

nb.set("cell-056", r'''print('\n' + '=' * 80)
print('ANALYSIS COMPLETE -- SUMMARY (v2.1, primary analysis set)')
print('=' * 80)

print(f'\nSample scope: primary set only (utils.primary_mask). '
      f'MP excluded from the study; PostSurgery held aside.')
for ds in DATASETS_ANALYSIS:
    print(f'  {ds}: n={len(harmonized_meta[ds])}')

print('\n1. Differential Abundance (fine-grained Study.Group scheme, unadjusted):')
for ds in DATASETS_ANALYSIS:
    for (g1, g2), results in da_results[ds].items():
        n_spc = results['species']['Significant'].sum()     if not results['species'].empty     else 0
        n_mtb = results['metabolites']['Significant'].sum() if not results['metabolites'].empty else 0
        print(f'  {ds} ({g1} vs {g2}): {n_spc} species, {n_mtb} metabolites sig.')

print('\n1B. Differential Abundance (grouped scheme, BOTH tests -- robust = both + same direction):')
print(robust_df.to_string(index=False))

print('\n2. Polyamine correlations (species level, per cohort):')
for ds in POLYAMINE_DATASETS_ANALYSIS:
    if ds in polyamine_corr and not polyamine_corr[ds].empty:
        n_sig = int((polyamine_corr[ds]['QValue'] < 0.05).sum())
        print(f'  {ds}: {n_sig} significant species-polyamine pairs (q<0.05)')

print('\n2B. Partial correlations (adjusted for the full confounder set; NB09 candidate pool):')
for ds, pc in partial_corr_results.items():
    n_sig = int((pc['QValue_partial'] < 0.05).sum())
    n_pos = int(((pc['QValue_partial'] < 0.05) & (pc['Rho_partial'] > 0)).sum())
    print(f'  {ds}: {len(pc)} pairs | {n_sig} q<0.05 | {n_pos} positive '
          f'(covariates: {pc["Covariates"].iloc[0]})')

# Read the replication frame directly rather than via the removed results pickle.
# `consistent_pairs` is the q<0.05 claim; `consistent_relaxed` is the q<0.10
# exploratory fallback -- the same if/else that decides which one is written to
# consistent_polyamine_correlations_genus_v2.csv (step 7.4).
_final_consistent = consistent_pairs if not consistent_pairs.empty else consistent_relaxed
_consistent_threshold = 'q<0.05' if not consistent_pairs.empty else 'q<0.10 (exploratory)'
print(f'\n3. Cross-cohort replication (genus level, KEGG-parent matched, '
      f'{_consistent_threshold}):')
print(f'  {len(_final_consistent)} genus-polyamine pairs significant in >=2 cohorts')
if not _final_consistent.empty and 'Same_direction' in _final_consistent.columns:
    n_same = int(_final_consistent['Same_direction'].sum())
    print(f'    of which same direction in every cohort tested: {n_same}')
    n_prox = int((_final_consistent['Proxy_used_in'] != 'none').sum())
    print(f'    resting on an N-acetyl proxy in >=1 cohort: {n_prox} '
          f'(weaker evidence; labelled in the CSV)')

print('\n4. Networks (full species-metabolite):')
for ds, G in networks.items():
    print(f'  {ds}: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges')

print('\n5. Stratified hub networks (power-limited by group size -- not a biological trend):')
for group_name in STAGE_GROUP_ORDER:
    G_g = stratified_networks.get(group_name, nx.Graph())
    n_grp = int((harmonized_meta[CRC_DS_HUB]['Stage.Group'] == group_name).sum())
    print(f'  {group_name} (n={n_grp}): {G_g.number_of_nodes()} nodes, '
          f'{G_g.number_of_edges()} edges')

print('\n6. CRC progression 10A (ordinal, FDR fixed inside utils):')
if not prog_corr.empty:
    n_sig_species_u     = prog_corr.loc[prog_corr['QVal_species_stage'] < 0.05, 'Species'].nunique()
    n_sig_metabolites_u = prog_corr.loc[prog_corr['QVal_metabolite_stage'] < 0.05, 'Metabolite'].nunique()
    n_sig_all3 = int(((prog_corr['QVal_species_stage'] < 0.05) &
                      (prog_corr['QVal_metabolite_stage'] < 0.05) &
                      (prog_corr['QVal_spc_mtb'] < 0.05)).sum())
    print(f'  {n_sig_species_u} unique species, {n_sig_metabolites_u} unique metabolites '
          f'track ordinal stage (q<0.05)')
    print(f'  {n_sig_all3} species-metabolite pairs significant in all 3 progression tests')

print('\n7. CRC progression 10B (binary Stage_I_II vs Stage_III_IV):')
if not da_early_late_spc.empty:
    print(f'  {int(da_early_late_spc["Significant"].sum())} species, '
          f'{int(da_early_late_mtb["Significant"].sum())} metabolites significantly different')

print('\n' + '=' * 80)
print(f'Results saved to: {RESULTS_DIR}')
print('  Tables  : Results/tables/*_v2.csv  (differential_abundance_v2.csv is the DA table)')
print('  Figures : Results/figures/**/*_v2.png')
print('  Pickles : Results/intermediate/preprocessed_data_reduced_v2.pkl')
print('            (no analysis_results_v2.pkl -- removed in v2.1, it had no readers)')''')

nb.save()
