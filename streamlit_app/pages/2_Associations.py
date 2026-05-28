import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import streamlit as st
import pandas as pd
from utils_st import show_figure, load_csv, metric_row

st.set_page_config(page_title="Associations", layout="wide")
st.title("NB02 — Statistical Associations")

metric_row([
    ("Sig. associations (q<0.05)", "6,465"),
    ("Surviving partial-corr adj.", "6,454"),
    ("PERMANOVA pseudo-F", "1.61"),
    ("PERMANOVA p-value", "0.011"),
    ("PERMANOVA R²", "0.009"),
])
st.divider()

# Volcano plots
st.subheader("Differential Abundance — YACHIDA (Healthy vs. Advanced CRC)")
col1, col2 = st.columns(2)
with col1:
    show_figure("nb02_volcano_metabolites", "Volcano plot — Metabolites")
with col2:
    show_figure("nb02_volcano_species", "Volcano plot — Species")

st.divider()

# DA tables with comparison selector
st.subheader("Differential Abundance Tables")
comparison = st.selectbox(
    "Select comparison",
    ["Healthy_vs_Advanced_CRC", "Healthy_vs_Early_CRC", "Early_CRC_vs_Advanced_CRC"],
)
modality = st.radio("Modality", ["Metabolites", "Species"], horizontal=True)

prefix = "metabolites" if modality == "Metabolites" else "species"
csv_name = f"YACHIDA-CRC-2019_Healthy_vs_{comparison.split('_vs_')[1]}_3grp_{prefix}_DA.csv"
if comparison == "Early_CRC_vs_Advanced_CRC":
    csv_name = f"YACHIDA-CRC-2019_Early_CRC_vs_Advanced_CRC_3grp_{prefix}_DA.csv"
elif comparison == "Healthy_vs_Advanced_CRC":
    csv_name = f"YACHIDA-CRC-2019_Healthy_vs_Advanced_CRC_3grp_{prefix}_DA.csv"
else:
    csv_name = f"YACHIDA-CRC-2019_Healthy_vs_Early_CRC_3grp_{prefix}_DA.csv"

df = load_csv(csv_name)
if df is None:
    # fallback to generic da_ files
    df = load_csv(f"da_{prefix.lower()}_{comparison.lower()}.csv")

if df is not None:
    sig = df[df.get("q", df.get("padj", pd.Series(dtype=float))) < 0.05] if any(
        c in df.columns for c in ["q", "padj"]
    ) else df
    st.caption(f"{len(sig)} significant rows shown (q < 0.05)")
    st.dataframe(sig, use_container_width=True, hide_index=True)
else:
    st.info("DA table not found for this comparison.")

st.divider()

# PERMANOVA and Mantel
col1, col2 = st.columns(2)
with col1:
    df = load_csv("permanova_results.csv")
    if df is not None:
        st.subheader("PERMANOVA Results")
        st.dataframe(df, use_container_width=True, hide_index=True)
with col2:
    df = load_csv("mantel_test_results.csv")
    if df is not None:
        st.subheader("Mantel Test (microbiome–metabolome co-structure)")
        st.dataframe(df, use_container_width=True, hide_index=True)

# Top correlated pairs
df = load_csv("CRC_progression_top_pairs.csv")
if df is not None:
    st.subheader("Top Species–Metabolite Correlation Pairs")
    st.dataframe(df.head(50), use_container_width=True, hide_index=True)
