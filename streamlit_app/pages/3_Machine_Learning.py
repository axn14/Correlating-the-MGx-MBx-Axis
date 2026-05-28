import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import streamlit as st
import pandas as pd
from utils_st import show_figure, load_csv, metric_row

st.set_page_config(page_title="Machine Learning", layout="wide")
st.title("NB03 — Machine Learning Benchmarking")

metric_row([
    ("Metabolite targets", "45"),
    ("Outer CV folds", "10"),
    ("Inner CV folds (Optuna)", "3"),
    ("Optuna trials / model / fold", "30"),
    ("RF mean R²", "0.052"),
])
st.divider()

# Model comparison (verified from NB03 outputs)
st.subheader("Model Comparison — Win Rate & Performance")

model_df = pd.DataFrame({
    "Model":     ["ElasticNet", "XGBoost", "SVR (RBF)", "Random Forest", "LightGBM"],
    "Wins / 45": [13, 11, 11, 8, 2],
    "Mean R²":   [-0.522, 0.002, 0.015, 0.052, -0.001],
    "Median R²": ["—", -0.005, "—", 0.032, "—"],
})
st.dataframe(model_df, use_container_width=True, hide_index=True)

st.markdown(
    """
    **Best-predicted metabolites:**
    - Sebacate (C08277): Random Forest, R² = 0.392, ρ = 0.675
    - N-Acetylputrescine (C02714): XGBoost, R² = 0.356, ρ = 0.645
    """
)
st.divider()

# Figures
st.subheader("Benchmark Figures")
col1, col2 = st.columns(2)
with col1:
    show_figure("nb03_benchmark_heatmap", "Per-target R² heatmap")
    show_figure("nb03_model_radar", "Radar chart — model profile")
    show_figure("nb03_power_curves", "Power / learning curves")
with col2:
    show_figure("nb03_model_comparison_violin", "R² distribution by model")
    show_figure("nb03_model_win_rate", "Win-rate bar chart")
    show_figure("nb03_r2_distribution", "R² histogram across targets")

st.divider()

# CSV tables
tab1, tab2, tab3 = st.tabs(["Best model per target", "Model performance summary", "SHAP results"])
with tab1:
    df = load_csv("ml_best_model_per_target.csv")
    if df is None:
        df = load_csv("ml_best_model_summary.csv")
    if df is not None:
        st.dataframe(df, use_container_width=True, hide_index=True)
    else:
        st.info("ml_best_model_per_target.csv not found.")

with tab2:
    df = load_csv("ml_model_performance_summary.csv")
    if df is None:
        df = load_csv("ml_benchmark_results.csv")
    if df is not None:
        st.dataframe(df, use_container_width=True, hide_index=True)

with tab3:
    df = load_csv("advanced_crc_ml_shap_results.csv")
    if df is not None:
        st.dataframe(df.head(100), use_container_width=True, hide_index=True)
        st.caption("Showing top 100 rows. Full table has more entries.")
