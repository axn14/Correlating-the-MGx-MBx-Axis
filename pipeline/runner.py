"""
runner.py — Orchestrates the full 5-step pipeline for user-submitted data.

Usage (from Python or Streamlit):
    from pipeline.runner import run_full_pipeline
    results = run_full_pipeline(
        species_df, metabolomics_df, metadata_df,
        group_col="condition",
        output_dir=Path("results/user_run"),
        progress_callback=lambda step, name, pct: print(f"[{pct}%] {name}"),
    )
"""

from __future__ import annotations

import threading
import traceback
from pathlib import Path
from typing import Callable

import pandas as pd

from pipeline import step1_preprocess, step2_associations, step3_ml, step4_mofa, step5_network


def run_full_pipeline(
    species_df:       pd.DataFrame,
    metabolomics_df:  pd.DataFrame,
    metadata_df:      pd.DataFrame,
    group_col:        str,
    output_dir:       Path,
    prev_spe:         float = 0.10,
    prev_mtb:         float = 0.15,
    n_top_species:    int   = 300,
    n_ml_targets:     int   = 45,
    mofa_factors:     int   = 10,
    progress_callback: Callable | None = None,
) -> dict:
    """
    Run all pipeline steps sequentially and return a results dict.
    progress_callback(step_num, step_name, pct_done) is called after each step.
    Raises on error; caller should catch and report to UI.
    """
    output_dir = Path(output_dir)
    results    = {}

    def _cb(step, name, pct):
        if progress_callback:
            progress_callback(step, name, pct)

    _cb(0, "Starting pipeline…", 0)

    results["preprocessed"] = step1_preprocess.run(
        species_df, metabolomics_df, metadata_df,
        group_col=group_col,
        output_dir=output_dir,
        prev_spe=prev_spe,
        prev_mtb=prev_mtb,
    )
    _cb(1, "Preprocessing complete", 20)

    results["associations"] = step2_associations.run(
        results["preprocessed"], output_dir=output_dir
    )
    _cb(2, "Associations complete", 40)

    results["ml"] = step3_ml.run(
        results["preprocessed"], results["associations"],
        output_dir=output_dir,
        n_top_species=n_top_species,
        n_targets=n_ml_targets,
    )
    _cb(3, "ML benchmarking complete", 65)

    results["mofa"] = step4_mofa.run(
        results["preprocessed"], output_dir=output_dir, n_factors=mofa_factors
    )
    _cb(4, "MOFA+ complete", 80)

    results["network"] = step5_network.run(
        results["preprocessed"], results["associations"], results["ml"],
        output_dir=output_dir,
    )
    _cb(5, "Network & mediation complete", 100)

    return results


def run_in_background(
    species_df:        pd.DataFrame,
    metabolomics_df:   pd.DataFrame,
    metadata_df:       pd.DataFrame,
    group_col:         str,
    output_dir:        Path,
    status_dict:       dict,
    **kwargs,
) -> threading.Thread:
    """
    Run the pipeline in a background thread. Update status_dict in-place.

    status_dict keys written:
      running (bool), step (int), step_name (str), pct (int),
      error (str|None), done (bool)
    """
    status_dict.update({"running": True, "step": 0, "step_name": "Starting…",
                         "pct": 0, "error": None, "done": False})

    def _progress(step, name, pct):
        status_dict.update({"step": step, "step_name": name, "pct": pct})

    def _worker():
        try:
            run_full_pipeline(
                species_df, metabolomics_df, metadata_df,
                group_col=group_col,
                output_dir=output_dir,
                progress_callback=_progress,
                **kwargs,
            )
            status_dict.update({"running": False, "done": True, "pct": 100})
        except Exception as e:
            status_dict.update({
                "running": False, "done": False,
                "error": f"{type(e).__name__}: {e}\n{traceback.format_exc()}",
            })

    t = threading.Thread(target=_worker, daemon=True)
    t.start()
    return t
