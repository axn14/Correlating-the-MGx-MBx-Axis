"""NB06 edits: v2.1 group labels, PostSurgery projection onto primary-set factors."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from nbedit import Notebook

nb = Notebook(Path(__file__).resolve().parent / "06_multiomics_factors.ipynb")

nb.set("cell-001", r'''# 06 — MOFA+ Joint Multi-Omic Factor Analysis

MOFA+ jointly decomposes microbiome and metabolome data into shared latent factors, each an
axis of co-variation across the two views (species abundances, metabolite concentrations).
This notebook fits a **single-group** MOFA+ model on the **primary analysis set**, prunes to
genuinely active factors, canonicalises factor sign, tests each active factor against CRC
stage (Kruskal–Wallis + Spearman, BH-FDR corrected), **projects the held-aside PostSurgery
samples onto the fitted factors**, and saves the full loadings/scores/association tables for
`08_pathway_analysis.ipynb` and `09_pathway_network.ipynb`.

**Reads:** `RESULTS_DIR/intermediate/preprocessed_data_reduced_v2.pkl` (NB02 output), whose
default keys are primary-only; the PostSurgery rows come from its `sensitivity_*` keys.

**Writes** (exact filenames NB08/NB09 depend on):
- `tables/mofa_factor_stage_association_v2.csv` — KW + Spearman vs `Stage.Group` /
  `Stage.Group.Numeric`, BH-FDR corrected across all tested factors.
- `tables/mofa_factor_loadings_v2.csv` — **every** active factor × **every** feature in both
  views (not a top-N subset).
- `tables/mofa_factor_scores_v2.csv` — `Sample` + every `FactorN` + `Study.Group` +
  `CRC.Stage.Numeric` + `Stage.Group` (+ `Stage.Group.Numeric`).
- `tables/mofa_postsurgery_projection_v2.csv` — the sensitivity group's projected scores.

## Sample scope

Factors are learned from the **277 primary samples** (Healthy / Early_CRC / Advanced_CRC)
only. The 30 `PostSurgery` samples are **not** used to fit the model; they are projected
afterwards (Section 14). That ordering matters: if they had been included in the fit, they
would have influenced the axes against which they are then compared, and any apparent
"PostSurgery sits between Healthy and CRC" statement would be partly circular.

## Sign-canonicalisation convention (read this if you are NB08/NB09)

MOFA+'s factor sign and scale are arbitrary — the model is invariant to flipping a factor's
sign and negating its loadings together. Before anything is saved, every factor's sign is
fixed so that, consistently:

> **A higher `FactorN` score means a more advanced CRC disease state**
> (positively correlated with `CRC.Stage.Numeric`, the fine-grained ordinal encoding from
> `utils.CRC_STAGE_NUMERIC`: Healthy=0, Stage_0=1, Stage_I_II=2, Stage_III_IV=3).

`PostSurgery` is deliberately **NaN** on that scale and therefore plays no part in sign
canonicalisation — it is not a point on a severity axis.

So a **positive** `Weight` for a feature on `FactorN` in `mofa_factor_loadings_v2.csv` means
that feature trends with the more-advanced-disease direction of that factor, uniformly for
every factor.

## Why a single MOFA+ group, not one group per stage

MOFA+ centres every feature **within each group**. With groups = stage categories, each
stage's mean is removed before the factors are learned, so between-stage differences —
exactly what the stage-association test looks for — cannot appear in any factor. The earlier
pipeline did precisely that and then tested for stage effects, which was self-defeating.
A single group keeps those differences and makes the stage test meaningful.

## Bugs fixed vs. the old `exploratory/06_multiomics_factors.ipynb`

1. Hardcoded dead path and wrong-drive input pickle — replaced with `sys.path.insert(0, '.')`
   and `RESULTS_DIR` from `utils`.
2. `calculate_variance_explained()` indexing — verified against the installed `mofapy2`
   source (`mofapy2/core/BayesNet.py`): it returns a **list of length `len(GROUPS)`**, each
   element an `(n_views, K)` array, indexed `[group][view]`. The old dict-style
   `.get('r2_per_factor', ...)` raised `AttributeError`, swallowed by a broad `try/except`,
   so the variance-explained figures never rendered on a real fit.
3. Undetected metabolite values become `NaN` (MOFA+'s native missing-value handling) using
   NB02's `undetected_mtb` mask (pre-imputation value missing or ≤ 0). The earlier
   "tied-for-column-minimum" detection changed 0 cells, because per-sample median centring
   breaks the ties.
4. Factor count is not hardcoded: fit with a generous ceiling (15), then prune via
   `removeInactiveFactors()` — looped, because that method removes only ONE factor per call.
5. BH-FDR correction across all tested factors before anything is called significant.
6. **Every** active factor's loadings are saved (both views, no top-N filter) — a hard NB08
   contract requirement, since NB08 addresses specific `FactorN` columns.
7. Factor sign canonicalised against `CRC.Stage.Numeric`.
8. All stage encoding comes from `utils` and is asserted against the pickle, never
   redefined locally.
9. Every output directory is created immediately before the `savefig()`/`to_csv()` using it.''')

# Stage-association header: 3 primary groups, no Adenoma.
nb.set("cell-017", r'''## 8. Factor ↔ CRC stage association (bug fixes #5, #8)

Kruskal–Wallis across `Stage.Group` (Healthy / Early_CRC / Advanced_CRC, per
`utils.STAGE_GROUP_MAP` / `STAGE_GROUP_ORDER`) plus Spearman correlation against the ordinal
`Stage.Group.Numeric` (primary) and the fine-grained `CRC.Stage.Numeric` (secondary), each
BH-FDR corrected across all tested factors.

`PostSurgery` is absent from the fit and from this test; it is projected and described
separately in Section 14.

The two tests answer different questions and can disagree legitimately: Kruskal–Wallis asks
whether *any* group differs, while Spearman asks for a monotonic gradient across the ordered
groups. A factor that rises steadily with stage but with wide within-group spread can reach
Spearman significance without reaching KW significance. Both are reported; neither is
presented as the single arbiter.''')

# Palette: primary groups + PostSurgery for the projection figure.
nb.set("cell-020", r'''# Primary stage groups get the ordered severity ramp; PostSurgery (projected,
# off-axis) is grey so it can never be mistaken for a point on that axis.
stage_palette = {'Healthy': '#2ecc71', 'Early_CRC': '#f39c12',
                 'Advanced_CRC': '#c0392b', 'PostSurgery': '#95a5a6'}
top_factors = kw_df.head(min(3, len(kw_df)))['Factor'].tolist()

if top_factors:
    fig, axes = plt.subplots(1, len(top_factors), figsize=(7 * len(top_factors), 5), squeeze=False)
    axes = axes[0]
    for ax, fc in zip(axes, top_factors):
        row = kw_df[kw_df['Factor'] == fc].iloc[0]
        data_plot = [Z_df.loc[Z_df['Stage.Group'] == sg, fc].dropna().to_numpy()
                     for sg in stage_groups_present]
        bp = ax.boxplot(data_plot, patch_artist=True,
                         medianprops=dict(color='black', linewidth=2),
                         flierprops=dict(marker='o', markersize=3, alpha=0.4))
        for patch, sg in zip(bp['boxes'], stage_groups_present):
            patch.set_facecolor(stage_palette.get(sg, '#cccccc'))
            patch.set_alpha(0.75)
        for j, (sg, vals) in enumerate(zip(stage_groups_present, data_plot)):
            ax.scatter(np.random.normal(j + 1, 0.07, size=len(vals)), vals,
                       color=stage_palette.get(sg, '#cccccc'), s=12, alpha=0.5, zorder=3)
        ax.set_xticks(range(1, len(stage_groups_present) + 1))
        ax.set_xticklabels([f'{sg}\n(n={len(v)})' for sg, v in
                            zip(stage_groups_present, data_plot)], fontsize=9)
        ax.set_ylabel(f'{fc} score (sign-canonicalized)', fontsize=10)
        ax.set_title(f"{fc}   KW q={row['KW_qval']:.4f}\n"
                     f"Spearman rho={row['Spearman_rho']:+.3f} (q={row['Spearman_qval']:.4f})",
                     fontsize=10)
        ax.axhline(0, color='grey', ls='--', lw=0.8, alpha=0.5)
        ax.grid(True, axis='y', alpha=0.3)

    plt.suptitle('MOFA+ factor scores by CRC Stage.Group (primary set, n=%d)\n'
                 f'top {len(top_factors)} of {len(kw_df)} tested factors by KW FDR-q'
                 % len(Z_df), fontsize=13, fontweight='bold')
    plt.tight_layout()

    figs_dir = RESULTS_DIR / 'figures' / 'multiomics'
    figs_dir.mkdir(parents=True, exist_ok=True)
    out_fig = figs_dir / 'mofa_factor_scores_by_stage_v2.png'
    fig.savefig(out_fig, dpi=150, bbox_inches='tight')
    plt.show(); plt.close()
    print(f'Saved: {out_fig}')
else:
    print('No factor had >=2 stage groups with n>3 samples -- boxplot skipped.')''')

# ---- new section 14: project the PostSurgery samples ----
nb.insert_after("cell-029", "cell-029b", "markdown", r'''## 14. Projecting the PostSurgery (HS) samples onto the fitted factors

The 30 `PostSurgery` samples took no part in fitting the model. Here they are projected onto
the already-learned factor axes and compared with the primary groups.

**Method.** For each held-out sample, its factor scores solve a least-squares problem against
the fitted loading matrix `W` (species and metabolite views stacked), after centring each
feature by the **primary-set** mean — the same reference the model was fit against. Features
that are missing for a sample (undetected metabolites) are dropped from that sample's system
rather than imputed, so a projection uses only observed evidence.

This is an approximation to a full variational projection: it ignores the view-specific noise
precisions MOFA+ learned, and so treats both views as equally reliable. It is adequate for
placing samples on an existing axis, which is all it is used for, and it is **not** used for
any significance test.

**What this can and cannot show.** A projected group sitting between Healthy and CRC is
consistent with "a resected gut is intermediate", but it cannot establish that: these
patients are post-treatment and formerly had cancer, so disease history, surgery and
treatment are inseparable here. NB01b already quantified how far this group sits from Healthy
in a direct test. This section is descriptive, and is reported as such.''')

nb.insert_after("cell-029b", "cell-029c", "code", r'''# Project the held-aside PostSurgery samples onto the fitted MOFA+ factors.
sens_spc_all = _prep.get('sensitivity_species', {}).get(CRC_DS)
sens_mtb_all = _prep.get('sensitivity_mtb', {}).get(CRC_DS)
sens_meta    = _prep.get('sensitivity_meta', {}).get(CRC_DS)

proj_df = pd.DataFrame()
if sens_spc_all is None or sens_meta is None or len(sens_meta) == 0:
    print('No sensitivity (PostSurgery) samples in the pickle -- projection skipped.')
else:
    # Same feature order as the fitted views.
    s_spc = sens_spc_all.set_index('Sample').reindex(columns=spc_names)
    s_mtb = sens_mtb_all.set_index('Sample').reindex(columns=mtb_names)
    common_sens = [s for s in sens_meta['Sample'] if s in s_spc.index and s in s_mtb.index]
    s_spc, s_mtb = s_spc.loc[common_sens], s_mtb.loc[common_sens]

    # Centre on the PRIMARY-set means: the reference the model was fitted against.
    mu_spc = np.nanmean(spc_aligned.to_numpy(dtype=float), axis=0)
    mu_mtb = np.nanmean(mtb_aligned.to_numpy(dtype=float), axis=0)

    W_all = np.vstack([W_spc, W_mtb])                 # (n_features, K)
    Xc = np.hstack([s_spc.to_numpy(dtype=float) - mu_spc,
                    s_mtb.to_numpy(dtype=float) - mu_mtb])   # (n_sens, n_features)

    Z_proj = np.full((len(common_sens), n_factors_active), np.nan)
    for i in range(len(common_sens)):
        obs = np.isfinite(Xc[i])
        if obs.sum() < n_factors_active * 3:          # too little evidence to place it
            continue
        sol, *_ = np.linalg.lstsq(W_all[obs], Xc[i][obs], rcond=None)
        Z_proj[i] = sol

    proj_df = pd.DataFrame(Z_proj, columns=factor_cols)
    proj_df.insert(0, 'Sample', common_sens)
    proj_df = proj_df.merge(
        sens_meta[['Sample', 'Study.Group', 'Stage.Group']], on='Sample', how='left')
    proj_df['Projected'] = True

    n_ok = int(np.isfinite(Z_proj).all(axis=1).sum())
    print(f'Projected {n_ok}/{len(common_sens)} PostSurgery samples onto '
          f'{n_factors_active} factors.')

    out_proj = RESULTS_DIR / 'tables' / 'mofa_postsurgery_projection_v2.csv'
    proj_df.to_csv(out_proj, index=False)
    print(f'Saved: {out_proj.name}')

    # Where does PostSurgery sit relative to the primary groups, per factor?
    comp_rows = []
    for fc in factor_cols:
        row = {'Factor': fc, 'PostSurgery_median': float(np.nanmedian(proj_df[fc]))}
        for sg in stage_groups_present:
            row[f'{sg}_median'] = float(Z_df.loc[Z_df['Stage.Group'] == sg, fc].median())
        comp_rows.append(row)
    comp = pd.DataFrame(comp_rows)
    print('\nFactor medians: primary groups (fitted) vs PostSurgery (projected)')
    print(comp.round(3).to_string(index=False))
    print('\nDescriptive only -- no test is run on the projected scores, and the projection\n'
          'ignores MOFA+\'s learned per-view noise precisions.')

    # Figure: primary groups as boxes, projected samples overlaid.
    if top_factors:
        fig, axes = plt.subplots(1, len(top_factors),
                                 figsize=(6.5 * len(top_factors), 5), squeeze=False)
        axes = axes[0]
        for ax, fc in zip(axes, top_factors):
            groups_plot = stage_groups_present + ['PostSurgery']
            data_plot = [Z_df.loc[Z_df['Stage.Group'] == sg, fc].dropna().to_numpy()
                         for sg in stage_groups_present]
            data_plot.append(proj_df[fc].dropna().to_numpy())
            bp = ax.boxplot(data_plot, patch_artist=True,
                            medianprops=dict(color='black', linewidth=2),
                            flierprops=dict(marker='o', markersize=3, alpha=0.4))
            for patch, sg in zip(bp['boxes'], groups_plot):
                patch.set_facecolor(stage_palette.get(sg, '#cccccc'))
                patch.set_alpha(0.75)
                if sg == 'PostSurgery':
                    patch.set_hatch('//')         # hatched = projected, not fitted
            for j, (sg, vals) in enumerate(zip(groups_plot, data_plot)):
                ax.scatter(np.random.normal(j + 1, 0.07, size=len(vals)), vals,
                           color=stage_palette.get(sg, '#cccccc'), s=12, alpha=0.5, zorder=3)
            ax.set_xticks(range(1, len(groups_plot) + 1))
            ax.set_xticklabels([f'{sg}\n(n={len(v)})' for sg, v in zip(groups_plot, data_plot)],
                               fontsize=8.5)
            ax.axhline(0, color='grey', ls='--', lw=0.8, alpha=0.5)
            ax.set_ylabel(f'{fc} score', fontsize=10)
            ax.set_title(fc, fontsize=10, fontweight='bold')
            ax.grid(True, axis='y', alpha=0.3)
        plt.suptitle('MOFA+ factors: primary groups (fitted) vs PostSurgery (projected, hatched)\n'
                     'PostSurgery took no part in fitting these axes',
                     fontsize=12, fontweight='bold')
        plt.tight_layout()
        figs_dir = RESULTS_DIR / 'figures' / 'multiomics'
        figs_dir.mkdir(parents=True, exist_ok=True)
        _pfig = figs_dir / 'mofa_postsurgery_projection_v2.png'
        fig.savefig(_pfig, dpi=150, bbox_inches='tight')
        plt.show(); plt.close()
        print(f'Saved: {_pfig}')''')

# Summary: mention the projection output.
nb.set("cell-030", r'''print('=== NB06 MOFA+ Analysis Summary ===')
print(f'\nModel: {"MOFA+ (mofapy2)" if MOFA_SUCCESS else "PCA fallback"}')
print(f'Active factors: {n_factors_active}')
print(f'Fitted on     : {len(common_samples)} PRIMARY samples '
      f'({dict(Z_df["Stage.Group"].value_counts())})')
if not proj_df.empty:
    print(f'Projected     : {len(proj_df)} PostSurgery samples (not used in fitting)')

print('\nFactors significant vs Stage.Group after BH-FDR (q<0.05):')
_sig_kw = kw_df[kw_df['Significant_KW_FDR05']]
print(_sig_kw[['Factor', 'KW_qval', 'Spearman_rho', 'Spearman_qval']].to_string(index=False)
      if not _sig_kw.empty else '  (none by Kruskal-Wallis)')
_sig_sp = kw_df[kw_df['Significant_Spearman_FDR05']]
print('\nFactors with a significant monotonic stage gradient (Spearman q<0.05):')
print(_sig_sp[['Factor', 'Spearman_rho', 'Spearman_qval', 'KW_qval']].to_string(index=False)
      if not _sig_sp.empty else '  (none by Spearman)')
if _sig_kw.empty and not _sig_sp.empty:
    print('\nNote: a Spearman-only hit means a monotonic trend across ordered groups that')
    print('Kruskal-Wallis (which ignores the ordering) does not detect. That is a real')
    print('result, not a contradiction -- the ordered test is simply more powerful for a')
    print('gradient.')

print('\nOutputs:')
print(f'  - tables/mofa_factor_stage_association_v2.csv  ({len(kw_df)} factors)')
print(f'  - tables/mofa_factor_loadings_v2.csv            ({len(loadings_df)} rows, '
      f'{n_factors_active} factors x {len(spc_names)+len(mtb_names)} features)')
print(f'  - tables/mofa_factor_scores_v2.csv              ({len(Z_df)} primary samples)')
if not proj_df.empty:
    print(f'  - tables/mofa_postsurgery_projection_v2.csv     ({len(proj_df)} samples)')
print('  - figures/multiomics/mofa_variance_explained_v2.png')
print('  - figures/multiomics/mofa_factor_scores_by_stage_v2.png')
print('  - figures/multiomics/mofa_factor_biplot_v2.png')
print('  - figures/multiomics/mofa_loading_heatmap_v2.png')
if not proj_df.empty:
    print('  - figures/multiomics/mofa_postsurgery_projection_v2.png')''')

nb.save()
