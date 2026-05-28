"""
step1_preprocess.py — Preprocessing wrapper for user-submitted data.

Applies the same pipeline as NB01:
  - Prevalence filtering (species ≥ prev_spe, metabolites ≥ prev_mtb)
  - CLR transform (pseudocount 1e-4) for species
  - log10(x+1) + sample-wise median centering for metabolites
  - Saves preprocessed_user.pkl to output_dir/intermediate/
"""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
import pandas as pd

# Allow importing utils from repo root
sys.path.insert(0, str(Path(__file__).parent.parent))
from utils import (
    clr_transform,
    log10_transform,
    qc_filter_species,
    qc_filter_metabolites,
    compute_sample_qc,
    compute_feature_qc,
    save_pickle,
)


def run(
    species_df: pd.DataFrame,
    metabolomics_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    group_col: str,
    output_dir: Path,
    prev_spe: float = 0.10,
    prev_mtb: float = 0.15,
) -> dict:
    """
    Returns preprocessed data dict with keys:
      species_raw, species_clr, metabolomics_raw, metabolomics_log,
      metadata, group_col, qc_species, qc_metabolomics
    """
    output_dir = Path(output_dir)
    inter_dir  = output_dir / "intermediate"
    inter_dir.mkdir(parents=True, exist_ok=True)

    # Prevalence filter
    species_filt = qc_filter_species(species_df, prevalence_min=prev_spe)
    mtb_filt     = qc_filter_metabolites(metabolomics_df, prevalence_min=prev_mtb)

    # Transforms
    species_clr = clr_transform(species_filt)
    mtb_log     = log10_transform(mtb_filt)

    # QC summaries
    qc_spe = compute_sample_qc(species_filt, "species")
    qc_mtb = compute_sample_qc(mtb_filt, "metabolomics")

    result = {
        "species_raw":        species_filt,
        "species_clr":        species_clr,
        "metabolomics_raw":   mtb_filt,
        "metabolomics_log":   mtb_log,
        "metadata":           metadata_df,
        "group_col":          group_col,
        "qc_species":         qc_spe,
        "qc_metabolomics":    qc_mtb,
        "n_species":          species_clr.shape[1],
        "n_metabolites":      mtb_log.shape[1],
        "n_samples":          len(metadata_df),
        "prev_spe_threshold": prev_spe,
        "prev_mtb_threshold": prev_mtb,
    }

    save_pickle(result, inter_dir / "preprocessed_user.pkl")
    return result
