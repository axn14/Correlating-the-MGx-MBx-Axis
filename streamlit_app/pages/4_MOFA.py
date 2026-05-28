import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import streamlit as st
from utils_st import show_figure, load_csv, metric_row

st.set_page_config(page_title="MOFA+", layout="wide")
st.title("NB04 — MOFA+ Multi-Omic Factor Analysis")

metric_row([
    ("Latent factors", "15"),
    ("Convergence (iterations)", "86"),
    ("Stage-associated factors", "3"),
    ("Max McFadden R² (stage)", "0.010"),
])
st.divider()

# Figures
st.subheader("Factor Analysis Plots")
col1, col2 = st.columns(2)
with col1:
    show_figure("nb04_mofa_factor_boxplots", "Factor scores by CRC stage")
    show_figure("nb04_mofa_variance_per_modality", "Variance explained per modality")
with col2:
    show_figure("nb04_mofa_combined_loadings", "Top feature loadings (species + metabolites)")
    show_figure("nb04_mofa_ordinal_variance_explained", "Stage variance explained (ordinal R²)")

st.caption(
    "All 3 stage-associated factors (q < 0.05) explain < 1% of stage variance each, "
    "consistent with heterogeneous CRC-stage microbiome signal."
)
st.divider()

# Tables
tab1, tab2, tab3 = st.tabs(["Factor–stage associations", "Top loadings", "Variance explained"])
with tab1:
    df = load_csv("mofa_factor_stage_association.csv")
    if df is None:
        df = load_csv("mofa_ordinal_regression_results.csv")
    if df is not None:
        st.dataframe(df, use_container_width=True, hide_index=True)

with tab2:
    df = load_csv("mofa_top_loadings.csv")
    if df is None:
        df = load_csv("mofa_factor_loadings.csv")
    if df is not None:
        st.dataframe(df.head(100), use_container_width=True, hide_index=True)

with tab3:
    df = load_csv("mofa_factor_variance_explained.csv")
    if df is not None:
        st.dataframe(df, use_container_width=True, hide_index=True)
