import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import streamlit as st
from utils_st import show_figure, load_csv, metric_row

st.set_page_config(page_title="Preprocessing", layout="wide")
st.title("NB01 — Data Preprocessing & Quality Control")

metric_row([
    ("Raw species", "57,702"),
    ("Species retained (≥10% prev.)", "4,392"),
    ("Raw metabolites", "450"),
    ("Metabolites retained (≥15% prev.)", "249"),
    ("Final samples", "347"),
])
st.divider()

st.subheader("Sample QC")
col1, col2 = st.columns(2)
with col1:
    show_figure("nb01_sample_qc", "Per-sample QC metrics")
with col2:
    show_figure("nb01_shannon_diversity", "Shannon diversity by CRC stage")

st.subheader("Feature Detection & Distributions")
col1, col2 = st.columns(2)
with col1:
    show_figure("nb01_feature_detection_rates", "Feature detection rates")
with col2:
    show_figure("nb01_feature_distributions_post_transform", "Post-transform distributions")

st.subheader("Ordination")
col1, col2 = st.columns(2)
with col1:
    show_figure("nb01_pca_species", "PCA — CLR-transformed species")
with col2:
    show_figure("nb01_pca_metabolites", "PCA — log₁₀ metabolites")

# Shannon diversity table
df = load_csv("shannon_diversity_test.csv")
if df is not None:
    st.subheader("Shannon Diversity Test Results")
    st.dataframe(df, use_container_width=True, hide_index=True)

# PCA variance
col1, col2 = st.columns(2)
with col1:
    df = load_csv("pca_variance_species.csv")
    if df is not None:
        st.subheader("PCA Variance — Species")
        st.dataframe(df, use_container_width=True, hide_index=True)
with col2:
    df = load_csv("pca_variance_metabolites.csv")
    if df is not None:
        st.subheader("PCA Variance — Metabolites")
        st.dataframe(df, use_container_width=True, hide_index=True)
