import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import streamlit as st
import pandas as pd
from utils_st import show_figure, load_csv, metric_row

st.set_page_config(page_title="Evidence Integration", layout="wide")
st.title("NB07–NB09 — Evidence Integration & TRINITY Validation")

metric_row([
    ("Evidence streams", "9  (E1–E9)"),
    ("TRINITY-validated genus pairs", "11"),
    ("Species-level trinity pairs", "3"),
])
st.divider()

tab_matrix, tab_trinity, tab_ref, tab_figs = st.tabs([
    "Evidence Matrix", "TRINITY Pairs", "E1–E9 Reference", "Figures"
])

# ── Tab 1: Evidence Matrix ────────────────────────────────────────────────────
with tab_matrix:
    st.subheader("General Metabolite Evidence Matrix (NB08)")
    df = load_csv("general_metabolite_evidence_matrix.csv")
    if df is not None:
        # Highlight trinity-validated rows if column present
        trinity_col = next((c for c in df.columns if "trinity" in c.lower()), None)
        if trinity_col:
            n_trinity = df[trinity_col].sum() if df[trinity_col].dtype == bool else (df[trinity_col] == True).sum()
            st.caption(f"{n_trinity} TRINITY-validated pairs in this table.")
        st.dataframe(df, use_container_width=True, hide_index=True)
    else:
        st.info("general_metabolite_evidence_matrix.csv not found.")

    df2 = load_csv("producer_evidence_matrix_v2.csv")
    if df2 is not None:
        st.subheader("Polyamine Producer Evidence Matrix")
        st.dataframe(df2, use_container_width=True, hide_index=True)

# ── Tab 2: TRINITY Pairs ──────────────────────────────────────────────────────
with tab_trinity:
    st.subheader("Confirmed TRINITY Pairs")

    species_trinity = pd.DataFrame({
        "Level":       ["Species", "Species", "Species", "Genus"],
        "Taxon":       [
            "Bilophila wadsworthia",
            "Sutterella wadsworthensis_A",
            "Fusobacterium animalis",
            "Alistipes faecigallinarum",
        ],
        "Metabolite":  ["Isoleucine", "Leucine", "Phenylalanine", "Succinate"],
        "Enrichment":  ["Advanced CRC", "Advanced CRC", "Advanced CRC", "Early + Advanced CRC"],
        "Adv. p-value": ["0.039", "<0.001", "0.001", "0.001"],
        "Notes": [
            "", "", "",
            "Genus-level only — not a species-level trinity pair",
        ],
    })
    st.dataframe(species_trinity, use_container_width=True, hide_index=True)

    st.markdown(
        """
        > **Note:** *Alistipes shahii* → GABA is **not** in the dataset.
        > The genus-level pair is *Alistipes faecigallinarum* → Succinate.
        """
    )

    df = load_csv("general_metabolite_top_candidates.csv")
    if df is not None:
        st.subheader("Top Novel Candidate Producers (NB08)")
        st.dataframe(df, use_container_width=True, hide_index=True)

# ── Tab 3: E1–E9 Reference ────────────────────────────────────────────────────
with tab_ref:
    st.subheader("Evidence Streams E1–E9")
    st.markdown(
        """
        | Stream | Code | Description | Cluster |
        |--------|------|-------------|---------|
        | SHAP stability | E1 | Mean \|SHAP\| across 10 outer folds; top-20 cutoff | STAT |
        | Spearman significance | E2 | Partial-corr q<0.05 (adj. age/BMI/sex/alcohol/Brinkman) | STAT |
        | KEGG enzyme | E3 | Known KEGG pathway EC entry for species→metabolite | ORTHOGONAL |
        | GLASSO co-abundance | E4 | Edge in GLASSO network (q<0.05) | ORTHOGONAL |
        | MOFA factor loading | E5 | Co-load on same stage-associated MOFA+ factor | ORTHOGONAL |
        | Mediation significance | E6 | Significant ACME, n_boot=1,000 (q<0.05) | ORTHOGONAL |
        | Within-stage Spearman | E7 | Significant in ≥1 CRC stage (q<0.05) | ORTHOGONAL |
        | GutSMASH BGC | E8 | BGC classification for metabolite class | TRINITY |
        | MICOM flux | E9 | Non-zero community FBA flux (tradeoff=0.5) | TRINITY |

        **TRINITY criterion:** STAT (E1 or E2 ≥ 1) **AND** E8 = 1 **AND** E9 = 1
        """
    )

# ── Tab 4: Figures ────────────────────────────────────────────────────────────
with tab_figs:
    show_figure("Trinity Validation", "TRINITY validation overview")
    col1, col2 = st.columns(2)
    with col1:
        show_figure("nb07_glasso_top_edges", "GLASSO top co-abundance edges (E4)")
    with col2:
        show_figure("nb08_shap_stage_bubble_heatmap", "SHAP × stage bubble heatmap")
    show_figure("nb08_directional_stage_bubble_heatmap", "Directional species–metabolite associations")
