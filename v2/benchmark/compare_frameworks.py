"""Benchmark step 3/3 — compare MIMOSA2 and MelonnPan against the v2 framework.

For each polyamine, takes the top candidate producer species from each method and asks
how much they agree, and what the framework surfaces that the single-method tools miss
(and vice versa). Writes:
  outputs/benchmark_comparison_v2.csv      per (polyamine, species) union with per-method flags
  outputs/benchmark_overlap_v2.csv         per-polyamine overlap summary (counts + Jaccard)
  outputs/benchmark_concordance_v2.png     overlap heatmap

The honest question this answers for reviewers: does the multi-evidence framework recover
what established tools recover, and does its extra evidence change the shortlist?

Run after the R step:  python compare_frameworks.py
"""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import utils

HERE = Path(__file__).resolve().parent
OUT  = HERE / 'outputs'; OUT.mkdir(exist_ok=True)
INP  = HERE / 'inputs'
TABLES = utils.RESULTS_DIR / 'tables'
TOP_K = 15   # per-method shortlist depth per polyamine

PA_NAMES = list(utils.POLYAMINE_KEGG.values())          # Putrescine, ... Arginine
PA_BY_KEGG = dict(utils.POLYAMINE_KEGG)                 # C00134 -> Putrescine
kegg_names = (pd.read_csv(INP / 'metabolite_kegg_names.tsv', sep='\t')
              if (INP / 'metabolite_kegg_names.tsv').exists() else pd.DataFrame())

def _norm(s):   # genus+species, GTDB-suffix-insensitive, for cross-tool name matching
    return (utils.normalise_genus_name(str(s).split()[0]) + ' ' + ' '.join(str(s).split()[1:])).strip().lower() \
        if str(s).split() else str(s).lower()

# ---- 1. Framework: supporting species per polyamine (>=1 non-statistical source) ----
fw = pd.read_csv(TABLES / 'final_evidence_synthesis_v2.csv')
sup_col = 'N_Sources_Supporting'
fw_top = {}
for pa in PA_NAMES:
    sub = fw[fw['Polyamine'] == pa].copy()
    if sub.empty:
        fw_top[pa] = set(); continue
    sub = sub.sort_values([sup_col, 'QValue_partial'], ascending=[False, True]).head(TOP_K)
    fw_top[pa] = {_norm(s) for s in sub['Species']}

# ---- 2. MelonnPan: top |weight| taxa per polyamine metabolite ----
mp_top = {pa: set() for pa in PA_NAMES}
mp_file = OUT / 'melonnpan_weights.tsv'
if mp_file.exists():
    w = pd.read_csv(mp_file, sep='\t')
    feat_col = w.columns[0]
    w = w.set_index(feat_col)
    # columns are metabolites (KEGG or names). Map polyamine KEGG -> column.
    for kegg, pa in PA_BY_KEGG.items():
        col = next((c for c in w.columns if c == kegg or str(c).startswith(kegg)), None)
        if col is None and not kegg_names.empty:  # try display name
            nm = kegg_names.loc[kegg_names['KEGG'] == kegg, 'Name']
            col = next((c for c in w.columns if nm.notna().any() and c in set(nm)), None)
        if col is not None:
            s = w[col].astype(float).abs().sort_values(ascending=False)
            mp_top[pa] = {_norm(i) for i in s[s > 0].head(TOP_K).index}
else:
    print('MelonnPan weights not found — its column will be empty.  # >>> run R step')

# ---- 3. MIMOSA2: top |var_share| taxa per polyamine compound ----
m2_top = {pa: set() for pa in PA_NAMES}
m2_file = OUT / 'mimosa2_contributions.tsv'
if m2_file.exists():
    c = pd.read_csv(m2_file, sep='\t')
    # expected columns ~ compound/metabolite, taxon/species, var_share/contribution
    comp_col = next((x for x in c.columns if x.lower() in ('compound', 'metabolite', 'met')), c.columns[0])
    tax_col  = next((x for x in c.columns if x.lower() in ('taxon', 'species', 'otu', 'feature')), c.columns[1])
    val_col  = next((x for x in c.columns if 'share' in x.lower() or 'contrib' in x.lower()), c.columns[-1])
    for kegg, pa in PA_BY_KEGG.items():
        sub = c[c[comp_col].astype(str).str.startswith(kegg)]
        if sub.empty:
            continue
        sub = sub.reindex(sub[val_col].abs().sort_values(ascending=False).index).head(TOP_K)
        m2_top[pa] = {_norm(t) for t in sub[tax_col]}
else:
    print('MIMOSA2 contributions not found — its column will be empty.  # >>> run R step')

# ---- 4. Build the union comparison + overlap summary ----
rows, ov = [], []
for pa in PA_NAMES:
    F, M, P = fw_top[pa], m2_top[pa], mp_top[pa]
    allsp = F | M | P
    for sp in sorted(allsp):
        rows.append({'Polyamine': pa, 'Species': sp,
                     'Framework': sp in F, 'MIMOSA2': sp in M, 'MelonnPan': sp in P,
                     'N_methods': int(sp in F) + int(sp in M) + int(sp in P)})
    def jac(a, b):
        return (len(a & b) / len(a | b)) if (a | b) else np.nan
    ov.append({'Polyamine': pa, 'Framework_n': len(F), 'MIMOSA2_n': len(M), 'MelonnPan_n': len(P),
               'FW_MIMOSA2_shared': len(F & M), 'FW_MelonnPan_shared': len(F & P),
               'MIMOSA2_MelonnPan_shared': len(M & P),
               'Jaccard_FW_MIMOSA2': jac(F, M), 'Jaccard_FW_MelonnPan': jac(F, P),
               'FW_unique': len(F - M - P)})

comp = pd.DataFrame(rows).sort_values(['Polyamine', 'N_methods'], ascending=[True, False])
comp.to_csv(OUT / 'benchmark_comparison_v2.csv', index=False)
ovdf = pd.DataFrame(ov)
ovdf.to_csv(OUT / 'benchmark_overlap_v2.csv', index=False)
print('Wrote benchmark_comparison_v2.csv and benchmark_overlap_v2.csv')
print(ovdf.to_string(index=False))

# ---- 5. Concordance figure: shared-counts heatmap ----
mat = ovdf.set_index('Polyamine')[['FW_MIMOSA2_shared', 'FW_MelonnPan_shared', 'MIMOSA2_MelonnPan_shared']]
if mat.notna().any().any():
    fig, ax = plt.subplots(figsize=(7, 0.5 * len(mat) + 2))
    im = ax.imshow(mat.values, cmap='YlGnBu', aspect='auto')
    ax.set_xticks(range(mat.shape[1]))
    ax.set_xticklabels(['FW∩MIMOSA2', 'FW∩MelonnPan', 'MIMOSA2∩MelonnPan'], rotation=20, ha='right')
    ax.set_yticks(range(len(mat))); ax.set_yticklabels(mat.index)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            ax.text(j, i, int(mat.values[i, j]), ha='center', va='center', fontsize=9)
    ax.set_title(f'Shared top-{TOP_K} producer species across methods\n(per polyamine)', fontsize=11)
    fig.colorbar(im, ax=ax, label='species in common')
    plt.tight_layout()
    fig.savefig(OUT / 'benchmark_concordance_v2.png', dpi=150, bbox_inches='tight')
    print('Wrote benchmark_concordance_v2.png')
else:
    print('No tool outputs to plot yet — run the R step first.')
