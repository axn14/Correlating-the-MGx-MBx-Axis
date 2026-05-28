"""
step4_mofa.py — MOFA+ multi-omic integration for user-submitted data.

Wraps mofapy2 to run with n_factors latent factors (default 10, reduced from 15
for smaller user datasets). Falls back gracefully if mofapy2 is not installed.
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent))
from utils import save_pickle

warnings.filterwarnings("ignore")


def run(preprocessed: dict, output_dir: Path, n_factors: int = 10) -> dict:
    output_dir = Path(output_dir)
    inter_dir  = output_dir / "intermediate"
    table_dir  = output_dir / "tables"
    inter_dir.mkdir(parents=True, exist_ok=True)
    table_dir.mkdir(parents=True, exist_ok=True)

    try:
        from mofapy2.run.entry_point import entry_point
    except ImportError:
        msg = "mofapy2 not installed — skipping MOFA+ step"
        warnings.warn(msg)
        return {"skipped": True, "reason": msg}

    species_clr = preprocessed["species_clr"]
    mtb_log     = preprocessed["metabolomics_log"]

    # MOFA+ expects list of views: each view is (samples × features)
    ent = entry_point()
    ent.set_data_options(scale_groups=False, scale_views=False)
    ent.set_data_matrix(
        [[species_clr.values, mtb_log.values]],
        likelihoods=["gaussian", "gaussian"],
        views_names=["species", "metabolomics"],
        groups_names=["all_samples"],
        samples_names=[list(species_clr.index)],
        features_names=[list(species_clr.columns), list(mtb_log.columns)],
    )
    ent.set_model_options(factors=n_factors)
    ent.set_train_options(
        iter=500, convergence_mode="fast",
        startELBO=1, freqELBO=5,
        seed=42, verbose=False,
    )
    ent.build()
    ent.run()

    # Extract results
    Z = pd.DataFrame(
        ent.model.getExpectations()["Z"]["E"][0],
        index=species_clr.index,
        columns=[f"Factor{i+1}" for i in range(n_factors)],
    )

    factor_loadings_spe = pd.DataFrame(
        ent.model.getExpectations()["W"]["E"][0],
        index=species_clr.columns,
        columns=[f"Factor{i+1}" for i in range(n_factors)],
    )
    factor_loadings_mtb = pd.DataFrame(
        ent.model.getExpectations()["W"]["E"][1],
        index=mtb_log.columns,
        columns=[f"Factor{i+1}" for i in range(n_factors)],
    )

    Z.to_csv(table_dir / "mofa_factor_scores.csv")
    factor_loadings_spe.to_csv(table_dir / "mofa_factor_loadings_species.csv")
    factor_loadings_mtb.to_csv(table_dir / "mofa_factor_loadings_metabolites.csv")

    result = {
        "factor_scores":         Z,
        "loadings_species":      factor_loadings_spe,
        "loadings_metabolomics": factor_loadings_mtb,
        "n_factors":             n_factors,
        "skipped":               False,
    }
    save_pickle(result, inter_dir / "mofa_results_user.pkl")
    return result
