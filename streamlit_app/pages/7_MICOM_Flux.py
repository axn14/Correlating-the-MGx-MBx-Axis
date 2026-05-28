import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import streamlit as st
from utils_st import show_figure, load_csv, metric_row

st.set_page_config(page_title="MICOM Flux", layout="wide")
st.title("NB09 — MICOM Community Flux Balance Analysis")

metric_row([
    ("Reference models", "AGORA103"),
    ("Tradeoff parameter", "0.5"),
    ("Initialisation", "Equal-abundance"),
    ("Evidence stream", "E9"),
])
st.divider()

# Figures
st.subheader("Flux Analysis Figures")
col1, col2 = st.columns(2)
with col1:
    show_figure("AGORA Coverage", "AGORA103 model coverage")
    show_figure("Polyamine detection rates", "Polyamine flux detection rates")
with col2:
    show_figure("Distribution", "Flux distribution across species")

st.info(
    "MICOM provides the E9 evidence stream: a non-zero predicted community metabolic flux "
    "for a species–metabolite pair is required for TRINITY status, alongside E8 (GutSMASH BGC) "
    "and the statistical streams (E1/E2)."
)
st.divider()

# Tables
tab1, tab2, tab3 = st.tabs(["Flux summary", "SHAP × MICOM cross-reference", "All-metabolite flux"])
with tab1:
    df = load_csv("micom_flux_summary.csv")
    if df is None:
        df = load_csv("micom_polyamine_flux_summary.csv")
    if df is not None:
        st.dataframe(df, use_container_width=True, hide_index=True)
    else:
        st.info("micom_flux_summary.csv not found.")

with tab2:
    df = load_csv("micom_shap_trinity_crossref.csv")
    if df is not None:
        st.dataframe(df, use_container_width=True, hide_index=True)

with tab3:
    df = load_csv("micom_all_metabolite_flux_summary.csv")
    if df is not None:
        st.dataframe(df, use_container_width=True, hide_index=True)
