import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import streamlit as st
from utils_st import show_figure, load_csv, metric_row

st.set_page_config(page_title="Mediation & Network", layout="wide")
st.title("NB05 — Partial-Correlation Network & Bootstrap Mediation")

metric_row([
    ("Network edges (q<0.05)", "1,992"),
    ("Species nodes", "406"),
    ("Metabolite nodes", "120"),
    ("Putrescine betweenness", "0.213"),
    ("Reverse ACME significant", "0 / 30"),
])
st.divider()

# Figures
col1, col2 = st.columns(2)
with col1:
    show_figure("nb05_centrality_scatter", "Node centrality — species vs metabolites")
with col2:
    show_figure("nb05_mediation_forest", "Bootstrap mediation — dual-direction forest plot")

st.info(
    "**Interpretation:** In the reverse direction (Species → Metabolite → CRC stage), "
    "zero of 30 top SHAP-ranked pairs showed a significant ACME (all 95% CIs crossed zero). "
    "This supports **direct microbial metabolite production** rather than stage-mediated metabolic cascades."
)
st.divider()

# Tables
tab1, tab2, tab3 = st.tabs(["Mediation results", "Hub metabolites", "Hub species"])
with tab1:
    df = load_csv("bootstrap_mediation_results.csv")
    if df is None:
        df = load_csv("mediation_results.csv")
    if df is not None:
        st.dataframe(df, use_container_width=True, hide_index=True)
    else:
        st.info("bootstrap_mediation_results.csv not found.")

with tab2:
    df = load_csv("hub_metabolites.csv")
    if df is not None:
        st.dataframe(df, use_container_width=True, hide_index=True)

with tab3:
    df = load_csv("hub_species.csv")
    if df is not None:
        st.dataframe(df, use_container_width=True, hide_index=True)

# Network centrality
df = load_csv("network_centrality.csv")
if df is not None:
    st.subheader("Full Network Centrality Table")
    st.dataframe(df.head(50), use_container_width=True, hide_index=True)
    st.caption("Showing top 50 rows by betweenness centrality.")
