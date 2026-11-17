"""NB09 edit: DEMOTE Module D (GLASSO network) to SUPPLEMENTARY.

Why: GLASSO is not one of the five scored evidence sources (predictive, genomic,
mechanistic, replication, co-loading) — it is a descriptive partial-correlation
network. In the Oct-1 run it also did not converge (dual gap negative at 300/400
iterations; see Results/logs/summary_v2.txt). This edit:

  1. relabels Module D as Supplementary in the markdown + figure title,
  2. raises GraphicalLassoCV / GraphicalLasso max_iter and tightens tol so it
     converges before anything is shown,
  3. nudges the alpha-grid floor up slightly (a larger minimum penalty is where
     the CV fit was failing on p>n),
  4. moves its outputs into a `supp_` namespace (figures/pathway_network/supplementary/,
     supp_d_glasso_*.csv/png). Nothing downstream reads these — GLASSO is not a
     scored source — so the rename is safe.

Idempotent and additive: it only `set()`s the four Module-D cells, so re-running
edits_nb09.py (which sets the synthesis-* cells) does not collide with it.
Run it after edits_nb09.py, before executing the notebook.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from nbedit import Notebook

nb = Notebook(Path(__file__).resolve().parent / "09_pathway_network.ipynb")

# ------------------------------------------------------------------ module-d-md
nb.set("module-d-md", r'''## Supplementary — Module D: GLASSO Partial-Correlation Network

**Supplementary, not a scored evidence source.** The five evidence sources counted in
the final synthesis are predictive (NB03), genomic (AGORA2), mechanistic (NB10 MICOM),
replication (NB02, external) and co-loading (NB06 MOFA). This GLASSO network is a
*descriptive* view of direct species-metabolite partial correlations (transitive
correlations through stage confounding removed); it informs none of the tiers and is
reported as a supplementary figure only.

A graphical lasso is fit on the joint species-metabolite matrix; candidate edges are
given an empirical, Westfall-Young max-T FWER p-value by a permutation null. Convergence
tolerance and iteration caps were raised relative to the first run (which did not
converge at 300/400 iterations on this p>n problem); if the CV fit still fails to
converge the module is skipped and marked unavailable rather than shown.''')

# ------------------------------------------------------------------ module-d-fit
nb.set("module-d-fit", r'''# SUPPLEMENTARY (not a scored evidence source). Feature selection for the GLASSO
# network: with all ~650 features (p > n = 317) the CV-selected alpha was so large that
# no edge survived. Use ~110 features instead:
#   top-60 species by variance + species significant in NB02's partial-correlation pool
#   (cap 80 species total) + all polyamines + top-40 other metabolites by variance.
pc_sig_species = (partial_corr.loc[partial_corr['QValue_partial'] < utils.FDR_THRESHOLD]
                  .sort_values('QValue_partial')['Species'].drop_duplicates().tolist())
sel_spc = list(spc_imp.var().nlargest(60).index)
for s in pc_sig_species:
    if len(sel_spc) >= 80:
        break
    if s in spc_imp.columns and s not in sel_spc:
        sel_spc.append(s)
pa_in_mtb = [c for c in polyamine_columns.values() if c in mtb_imp.columns]
sel_mtb = pa_in_mtb + [c for c in mtb_imp.var().sort_values(ascending=False).index
                       if c not in pa_in_mtb][:40]

spc_d = spc_imp[sel_spc]
mtb_d = mtb_imp[sel_mtb]
joint = np.hstack([
    StandardScaler().fit_transform(spc_d.values),
    StandardScaler().fit_transform(mtb_d.values),
])
joint_cols = list(spc_d.columns) + list(mtb_d.columns)
n_spc = spc_d.shape[1]
n_mtb = mtb_d.shape[1]
print(f"[Supplementary Module D] Joint matrix: {joint.shape} "
      f"({n_spc} species + {n_mtb} metabolites, {len(pa_in_mtb)} polyamines)")

# Alpha grid floor raised (-2.3 -> -2.0): the smallest penalties are where the CV fit
# failed to converge on this p>n problem. max_iter 300 -> 1000, tol 1e-3 -> 1e-4.
ALPHAS = np.logspace(-2.0, -0.3, 10)
print(f"Fitting GraphicalLassoCV over alphas {np.round(ALPHAS, 4).tolist()} "
      f"(max_iter=1000, tol=1e-4) ...")
glasso = GraphicalLassoCV(alphas=ALPHAS, cv=3, max_iter=1000, tol=1e-4, n_jobs=1)  # 1: memory cap
HAS_GLASSO = True
try:
    glasso.fit(joint)
except Exception as fit_exc:
    print(f"!!! Supplementary Module D: GraphicalLassoCV failed to converge: {fit_exc} "
          f"— Module D is SKIPPED (it is supplementary and not a scored source).")
    HAS_GLASSO = False

if HAS_GLASSO:
    prec = glasso.precision_
    print(f"Alpha selected: {glasso.alpha_:.4f}"
          f"{'  (at the edge of the grid)' if glasso.alpha_ in (ALPHAS.min(), ALPHAS.max()) else ''}")
    diag = np.sqrt(np.maximum(np.diag(prec), 1e-10))
    pcorr = -prec / np.outer(diag, diag)
    np.fill_diagonal(pcorr, 0)
    pcorr_cross = pcorr[:n_spc, n_spc:]   # rows=species, cols=metabolites
    print(f"Partial correlation matrix: {pcorr.shape}, cross-block: {pcorr_cross.shape}, "
          f"non-zero cross edges: {int((np.abs(pcorr_cross) > 1e-8).sum())}")''')

# ------------------------------------------------------------------ module-d-permutation
nb.set("module-d-permutation", r'''# ── SUPPLEMENTARY Module D: permutation null model for GLASSO edges ──────────
# Shuffles the sample correspondence of the metabolite block only (species block stays
# in its original sample order), breaking species-metabolite association while
# preserving each block's own within-block covariance structure, then refits
# GraphicalLasso at the ALREADY-SELECTED alpha (no re-CV — far too slow per permutation).
N_PERMUTATIONS = 100
MAGNITUDE_PREFILTER = 0.08   # candidate-edge magnitude floor, for reporting only —
                              # final "Significant" status is set by the permutation
                              # test below, not by this threshold alone.

HAS_GLASSO_NULL = False
glasso_edges = pd.DataFrame(columns=['Species', 'Metabolite', 'PartialCorr', 'AbsPartialCorr',
                                      'Empirical_PValue', 'QValue', 'Significant', 'Validated'])
SUPP_EDGES_CSV = 'supp_d_glasso_partial_correlation_edges_v2.csv'

if HAS_GLASSO:
    rng = np.random.default_rng(42)
    fixed_alpha = glasso.alpha_
    spc_scaled = StandardScaler().fit_transform(spc_d.values)
    mtb_scaled = StandardScaler().fit_transform(mtb_d.values)

    null_abs_list = []
    for perm_i in range(N_PERMUTATIONS):
        perm_idx = rng.permutation(len(common))
        joint_perm = np.hstack([spc_scaled, mtb_scaled[perm_idx]])
        gl_perm = GraphicalLasso(alpha=fixed_alpha, max_iter=1000, tol=1e-4)
        try:
            gl_perm.fit(joint_perm)
        except Exception as perm_exc:
            print(f"  Permutation {perm_i}: GraphicalLasso failed to converge — skipped "
                  f"({perm_exc}).")
            continue
        prec_perm = gl_perm.precision_
        diag_perm = np.sqrt(np.maximum(np.diag(prec_perm), 1e-10))
        pcorr_perm = -prec_perm / np.outer(diag_perm, diag_perm)
        null_abs_list.append(np.abs(pcorr_perm[:n_spc, n_spc:]))

    n_valid_perms = len(null_abs_list)
    print(f"Permutations completed: {n_valid_perms}/{N_PERMUTATIONS} converged.")

    observed_abs = np.abs(pcorr_cross)
    rows, cols = np.where(observed_abs >= MAGNITUDE_PREFILTER)
    print(f"Candidate edges (|partial_corr| >= {MAGNITUDE_PREFILTER}): {len(rows)}")

    if n_valid_perms < 5:
        print(f"WARNING: only {n_valid_perms} permutations converged — permutation-based "
              f"significance is unreliable. Saving magnitude-threshold candidate edges "
              f"with Validated=False (UNVALIDATED — no statistical null-model support) "
              f"rather than silently claiming FDR-corrected significance.")
        glasso_edges = pd.DataFrame({
            'Species':          [spc_d.columns[r] for r in rows],
            'Metabolite':       [mtb_d.columns[c] for c in cols],
            'PartialCorr':      [pcorr_cross[r, c] for r, c in zip(rows, cols)],
        })
        glasso_edges['AbsPartialCorr']   = glasso_edges['PartialCorr'].abs()
        glasso_edges['Empirical_PValue'] = np.nan
        glasso_edges['QValue']           = np.nan
        glasso_edges['Significant']      = False
        glasso_edges['Validated']        = False
    else:
        null_stack = np.stack(null_abs_list, axis=0)   # (n_valid_perms, n_spc, n_mtb)
        counts_ge = (null_stack >= observed_abs[None, :, :]).sum(axis=0)
        emp_p = (counts_ge + 1) / (n_valid_perms + 1)

        all_rows, all_cols = np.where(np.ones_like(observed_abs, dtype=bool))
        full_df = pd.DataFrame({
            'Species':          [spc_d.columns[r] for r in all_rows],
            'Metabolite':       [mtb_d.columns[c] for c in all_cols],
            'PartialCorr':      [pcorr_cross[r, c]  for r, c in zip(all_rows, all_cols)],
            'Empirical_PValue': [emp_p[r, c]        for r, c in zip(all_rows, all_cols)],
        })
        full_df['AbsPartialCorr'] = full_df['PartialCorr'].abs()
        _, qvals, _, _ = multipletests(full_df['Empirical_PValue'], method='fdr_bh')
        full_df['QValue'] = qvals
        # Westfall-Young max-statistic p-value (FWER over ALL cross edges): per-edge BH
        # cannot work here -- with ~n_spc*n_mtb edges and N_PERMUTATIONS permutations the
        # smallest attainable BH q is far above 0.05. The max-T null compares each
        # observed |partial corr| with the largest |partial corr| of every permuted fit.
        max_null = null_stack.reshape(n_valid_perms, -1).max(axis=1)
        p_maxT = (np.sum(max_null[:, None] >= observed_abs.ravel()[None, :], axis=0) + 1) / (n_valid_perms + 1)
        p_maxT = p_maxT.reshape(observed_abs.shape)
        full_df['P_maxT_FWER'] = [p_maxT[r, c] for r, c in zip(all_rows, all_cols)]
        full_df['Significant'] = full_df['P_maxT_FWER'] < 0.05
        print(f"max-T null: 95th percentile of max |partial corr| = {np.quantile(max_null, 0.95):.3f}")
        full_df['Validated'] = True
        HAS_GLASSO_NULL = True

        glasso_edges = full_df[
            (full_df['AbsPartialCorr'] >= MAGNITUDE_PREFILTER) | full_df['Significant']
        ].sort_values('AbsPartialCorr', ascending=False).reset_index(drop=True)

    glasso_edges.to_csv(TABLES / SUPP_EDGES_CSV, index=False)
    n_sig = int(glasso_edges['Significant'].sum()) if 'Significant' in glasso_edges else 0
    print(f"[Supplementary] GLASSO edges saved: {len(glasso_edges)} "
          f"(max-T permutation FWER-significant: {n_sig})")
    print(f"Saved: {SUPP_EDGES_CSV}")
else:
    print("Supplementary Module D skipped (GraphicalLassoCV did not converge) — no edges.")
    glasso_edges.to_csv(TABLES / SUPP_EDGES_CSV, index=False)
    print(f"Saved empty {SUPP_EDGES_CSV} (Module D unavailable).")''')

# ------------------------------------------------------------------ module-d-network-fig
nb.set("module-d-network-fig", r'''PA_KEGG_IDS = set(utils.POLYAMINE_KEGG.keys())
PA_COLS = {c for c in mtb_imp.columns
           if any(c.startswith(k + '_') for k in PA_KEGG_IDS) or c in polyamine_columns.values()}

# Supplementary figures live in their own subfolder and carry a supp_ prefix.
SUPP_FIG_DIR = FIGURES / 'supplementary'
SUPP_FIG_DIR.mkdir(parents=True, exist_ok=True)

if HAS_GLASSO and not glasso_edges.empty:
    import networkx as nx

    top_edges = glasso_edges.head(60)
    G = nx.Graph()
    for _, row in top_edges.iterrows():
        G.add_node(row['Species'], node_type='species')
        G.add_node(row['Metabolite'], node_type='metabolite')
        G.add_edge(row['Species'], row['Metabolite'],
                   weight=abs(row['PartialCorr']),
                   sign='pos' if row['PartialCorr'] > 0 else 'neg')

    color_map, size_map = {}, {}
    for n in G.nodes():
        if G.nodes[n].get('node_type') == 'species':
            color_map[n] = '#3498db'; size_map[n] = 160
        elif n in PA_COLS:
            color_map[n] = '#e74c3c'; size_map[n] = 400
        else:
            color_map[n] = '#2ecc71'; size_map[n] = 200
    node_colors = [color_map[n] for n in G.nodes()]
    node_sizes  = [size_map[n] for n in G.nodes()]
    edge_colors = ['#e74c3c' if G[u][v].get('sign') == 'neg' else '#27ae60'
                   for u, v in G.edges()]

    fig, ax = plt.subplots(figsize=(18, 12))
    pos = nx.spring_layout(G, k=1.8, seed=42)
    nx.draw(G, pos, node_color=node_colors, node_size=node_sizes,
            edge_color=edge_colors, alpha=0.75, ax=ax, with_labels=False, width=0.9)
    labels = {n: (n.split('_', 1)[-1] if '_' in n else n)[:28]
              for n in G.nodes()
              if G.nodes[n].get('node_type') == 'metabolite' or G.degree(n) >= 4}
    nx.draw_networkx_labels(G, pos, labels, font_size=7, ax=ax)

    legend = [
        mpatches.Patch(color='#e74c3c', label='Polyamine metabolite'),
        mpatches.Patch(color='#2ecc71', label='Other metabolite'),
        mpatches.Patch(color='#3498db', label='Species'),
    ]
    ax.legend(handles=legend, fontsize=9)
    validated_note = 'permutation-FDR validated' if HAS_GLASSO_NULL else 'UNVALIDATED (magnitude only)'
    ax.set_title(f'Supplementary Figure — GLASSO Partial Correlation Network '
                 f'(top {len(top_edges)} edges, {validated_note})\n'
                 'Descriptive only; not a scored evidence source. '
                 'Red edges = negative partial correlation, green = positive', fontsize=12)
    plt.tight_layout()
    out_d = SUPP_FIG_DIR / 'supp_d_glasso_network_v2.png'
    fig.savefig(out_d, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {out_d}")
else:
    print("Supplementary Module D network figure skipped — no GLASSO edges available.")''')

nb.save()
