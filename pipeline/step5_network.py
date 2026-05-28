"""
step5_network.py — Partial-correlation network and bootstrap mediation.

Builds a networkx graph from significant species–metabolite Spearman pairs,
computes centrality, and runs bootstrap mediation for top SHAP-ranked pairs.
Evidence streams E1–E7 are scored here.
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent))
from utils import (
    bootstrap_mediation,
    annotate_pathway,
    extract_genus,
    save_pickle,
    FDR_THRESHOLD,
)

warnings.filterwarnings("ignore")


def _build_network(spearman_sig: pd.DataFrame):
    try:
        import networkx as nx
    except ImportError:
        return None, pd.DataFrame()

    G = nx.Graph()
    for _, row in spearman_sig.iterrows():
        G.add_edge(row["species"], row["metabolite"], weight=abs(row.get("rho", 0)))

    centrality = nx.betweenness_centrality(G, weight="weight", normalized=True)
    cent_df = pd.DataFrame([
        {"node": n, "betweenness": c, "degree": G.degree(n)}
        for n, c in centrality.items()
    ]).sort_values("betweenness", ascending=False)

    return G, cent_df


def run(
    preprocessed: dict,
    associations: dict,
    ml_results: dict,
    output_dir: Path,
    n_mediation_pairs: int = 20,
    n_boot: int = 500,
) -> dict:
    output_dir = Path(output_dir)
    inter_dir  = output_dir / "intermediate"
    table_dir  = output_dir / "tables"
    inter_dir.mkdir(parents=True, exist_ok=True)
    table_dir.mkdir(parents=True, exist_ok=True)

    meta      = preprocessed["metadata"]
    group_col = preprocessed["group_col"]
    species   = preprocessed["species_clr"]
    mtb       = preprocessed["metabolomics_log"]

    spearman_all = associations.get("spearman_all", pd.DataFrame())
    spearman_sig = (
        spearman_all[spearman_all["q"] < FDR_THRESHOLD]
        if not spearman_all.empty and "q" in spearman_all.columns
        else pd.DataFrame()
    )

    # Network
    G, centrality_df = _build_network(spearman_sig)
    if not centrality_df.empty:
        centrality_df.to_csv(table_dir / "network_centrality.csv", index=False)

    hub_metabolites = pd.DataFrame()
    hub_species     = pd.DataFrame()
    if not centrality_df.empty:
        all_mtb_cols  = set(mtb.columns)
        all_spe_cols  = set(species.columns)
        hub_metabolites = centrality_df[
            centrality_df["node"].isin(all_mtb_cols)
        ].head(30)
        hub_species = centrality_df[
            centrality_df["node"].isin(all_spe_cols)
        ].head(30)
        hub_metabolites.to_csv(table_dir / "hub_metabolites.csv", index=False)
        hub_species.to_csv(table_dir / "hub_species.csv", index=False)

    # Bootstrap mediation — reverse direction (Species → Metabolite → Outcome)
    # Select top pairs from ML results
    best_df = ml_results.get("best_per_target", pd.DataFrame())
    if best_df.empty or spearman_sig.empty:
        mediation_df = pd.DataFrame()
    else:
        top_targets = best_df.nlargest(n_mediation_pairs, "r2")["target"].tolist()
        # Encode group as numeric outcome
        le_map = {g: i for i, g in enumerate(sorted(meta[group_col].dropna().unique()))}
        y_outcome = meta[group_col].map(le_map).fillna(0).values

        confounders = associations.get("confounders_used", [])
        cov = meta[confounders].apply(pd.to_numeric, errors="coerce").values if confounders else None

        med_rows = []
        for target in top_targets:
            if target not in mtb.columns:
                continue
            # Pick best-correlated species for this metabolite
            pairs = spearman_sig[spearman_sig["metabolite"] == target]
            if pairs.empty:
                continue
            best_spe = pairs.iloc[0]["species"]
            if best_spe not in species.columns:
                continue
            x = species[best_spe].values
            m = mtb[target].values
            try:
                med = bootstrap_mediation(
                    x, m, y_outcome,
                    n_boot=n_boot, ci=0.95,
                    random_state=42,
                    covariates=cov,
                )
                med_rows.append({
                    "species":    best_spe,
                    "metabolite": target,
                    "acme":       med.get("acme", np.nan),
                    "acme_lo":    med.get("acme_ci_lo", np.nan),
                    "acme_hi":    med.get("acme_ci_hi", np.nan),
                    "p_acme":     med.get("p_acme", np.nan),
                    "direction":  "reverse",
                })
            except Exception:
                pass

        mediation_df = pd.DataFrame(med_rows)
        if not mediation_df.empty:
            mediation_df.to_csv(table_dir / "bootstrap_mediation_results.csv", index=False)

    result = {
        "network_nodes":    G.number_of_nodes() if G else 0,
        "network_edges":    G.number_of_edges() if G else 0,
        "centrality":       centrality_df,
        "hub_metabolites":  hub_metabolites,
        "hub_species":      hub_species,
        "mediation":        mediation_df,
    }
    save_pickle(result, inter_dir / "network_results_user.pkl")
    return result
