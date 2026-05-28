"""
step2_associations.py — Statistical associations for user-submitted data.

Runs the same pipeline as NB02:
  - Differential abundance (Mann-Whitney U + BH FDR) for each group pair
  - Pairwise Spearman correlations (pre-filtered |r| >= min_corr)
  - Partial correlations adjusted for available confounders
  - Species co-abundance matrix
  - Saves analysis_results_user.pkl + key CSVs to output_dir
"""

from __future__ import annotations

import sys
from pathlib import Path
from itertools import combinations

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent))
from utils import (
    differential_abundance,
    spearman_correlation_matrix,
    partial_corr_residuals,
    species_coabundance_matrix,
    save_pickle,
    FDR_THRESHOLD,
    MIN_CORR,
)


_CONFOUNDER_CANDIDATES = ["age", "bmi", "sex", "gender", "alcohol", "smoking"]


def _available_confounders(metadata_df: pd.DataFrame) -> list[str]:
    """Return confounders present as columns in metadata (case-insensitive)."""
    meta_lower = {c.lower(): c for c in metadata_df.columns}
    return [meta_lower[c] for c in _CONFOUNDER_CANDIDATES if c in meta_lower]


def run(preprocessed: dict, output_dir: Path) -> dict:
    """
    Returns associations dict with keys:
      da_results, spearman_sig, partial_corr_sig, coabundance
    """
    output_dir = Path(output_dir)
    inter_dir  = output_dir / "intermediate"
    table_dir  = output_dir / "tables"
    inter_dir.mkdir(parents=True, exist_ok=True)
    table_dir.mkdir(parents=True, exist_ok=True)

    species_clr = preprocessed["species_clr"]
    mtb_log     = preprocessed["metabolomics_log"]
    meta        = preprocessed["metadata"]
    group_col   = preprocessed["group_col"]

    groups = sorted(meta[group_col].dropna().unique())

    # Differential abundance for all group pairs
    da_results: dict[str, pd.DataFrame] = {}
    for g1, g2 in combinations(groups, 2):
        key = f"{g1}_vs_{g2}"
        for label, df, transform in [
            ("species",      species_clr, "clr"),
            ("metabolomics", mtb_log,     "log"),
        ]:
            da = differential_abundance(df, meta, group_col, g1, g2, transform)
            da_results[f"{label}_{key}"] = da
            da.to_csv(table_dir / f"da_{label}_{key}.csv")

    # Spearman species–metabolite correlations
    spearman_all = spearman_correlation_matrix(
        species_clr, mtb_log, fdr=FDR_THRESHOLD, min_r=MIN_CORR
    )
    spearman_all.to_csv(table_dir / "spearman_correlations_all.csv", index=False)

    # Partial correlations (adjust for available confounders)
    confounders = _available_confounders(meta)
    if confounders:
        partial_rows = []
        sig = spearman_all[spearman_all["q"] < FDR_THRESHOLD]
        for _, row in sig.iterrows():
            cov = meta[confounders].apply(pd.to_numeric, errors="coerce")
            rp = partial_corr_residuals(
                species_clr[row["species"]].values,
                mtb_log[row["metabolite"]].values,
                cov.values,
            )
            partial_rows.append({**row.to_dict(), "partial_rho": rp["r"], "partial_p": rp["p-val"]})
        partial_df = pd.DataFrame(partial_rows)
        partial_df.to_csv(table_dir / "spearman_correlations_partial.csv", index=False)
    else:
        partial_df = spearman_all.copy()
        partial_df["partial_rho"] = partial_df.get("rho", 0.0)

    # Species co-abundance
    coabundance = species_coabundance_matrix(species_clr, fdr=FDR_THRESHOLD, min_r=MIN_CORR)
    coabundance.to_csv(table_dir / "species_coabundance.csv", index=False)

    result = {
        "da_results":       da_results,
        "spearman_all":     spearman_all,
        "spearman_partial": partial_df,
        "coabundance":      coabundance,
        "confounders_used": confounders,
        "groups":           groups,
        "n_sig_pairs":      int((spearman_all["q"] < FDR_THRESHOLD).sum()),
    }

    save_pickle(result, inter_dir / "analysis_results_user.pkl")
    return result
