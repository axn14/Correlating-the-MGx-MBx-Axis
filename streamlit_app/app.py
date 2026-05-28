import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

import streamlit as st
import pandas as pd
from utils_st import metric_row, load_csv

st.set_page_config(
    page_title="CRC Multi-Omic Pipeline",
    page_icon="🧬",
    layout="wide",
)

st.title("CRC Multi-Omic Pipeline — Results Explorer")
st.markdown(
    """
    This dashboard presents results from an integrated gut metagenomics and serum metabolomics
    study of colorectal cancer (CRC) progression using the **YACHIDA-CRC-2019** cohort (n=347).
    Nine analysis notebooks identified microbial producers of cancer-associated metabolites through
    a nine-stream evidence integration framework.

    Use the **sidebar** to navigate between analysis sections.
    """
)

st.divider()

# ── Key findings ──────────────────────────────────────────────────────────────
st.subheader("Key Findings")
metric_row([
    ("PERMANOVA p-value", "0.011"),
    ("Sig. species–metabolite pairs", "6,465"),
    ("TRINITY-validated genus pairs", "11"),
    ("Species-level trinity pairs", "3"),
])

st.divider()

# ── Cohort demographics ───────────────────────────────────────────────────────
st.subheader("Cohort Demographics — YACHIDA-CRC-2019")

cohort_df = load_csv("cohort_summary_yachida.csv")
if cohort_df is not None:
    st.dataframe(cohort_df, use_container_width=True, hide_index=True)
else:
    cohort = pd.DataFrame({
        "Stage": ["Healthy", "Stage_0", "High-risk adenoma (HS)",
                  "Stage I–II CRC", "Stage III–IV CRC", "Multiple polyp (MP)", "Total"],
        "n": [127, 27, 30, 69, 54, 40, 347],
    })
    st.dataframe(cohort, use_container_width=True, hide_index=True)

st.divider()

# ── Pipeline overview ─────────────────────────────────────────────────────────
st.subheader("Analysis Pipeline")
st.markdown(
    """
    | Notebook | Section | Key output |
    |----------|---------|------------|
    | NB01 | Preprocessing & QC | 4,392 species, 249 metabolites retained |
    | NB02 | Statistical associations | 6,465 sig. pairs; PERMANOVA R²=0.009 |
    | NB03 | ML benchmarking | RF best overall (mean R²=0.052) |
    | NB04 | MOFA+ integration | 3 stage-associated factors (q<0.05) |
    | NB05 | Mediation & network | 0/30 reverse ACME significant |
    | NB06 | GutSMASH BGC mining | E8 stream: biosynthetic gene clusters |
    | NB07 | Advanced evidence integration | E1–E8 evidence matrix (polyamine focus) |
    | NB08 | General source attribution | All-metabolite evidence matrix |
    | NB09 | MICOM flux | E9 stream: community FBA; 3 species trinity pairs |
    """
)

st.info("Navigate using the sidebar to explore each analysis section.")
