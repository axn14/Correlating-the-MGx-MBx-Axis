"""
0_Upload_and_Run.py — Upload your own data and run the pipeline.

Streamlit page that: accepts 3 file uploads → validates format →
lets user configure parameters → runs the 5-step pipeline in a
background thread → shows real-time progress.
"""

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
sys.path.insert(0, str(Path(__file__).parent.parent))

import streamlit as st
import pandas as pd

st.set_page_config(page_title="Upload & Run", layout="wide")
st.title("Upload Your Own Data & Run Pipeline")

st.info(
    "Upload a species abundance table, metabolomics table, and sample metadata. "
    "The pipeline will run steps 1–5 (Preprocessing → Associations → ML → MOFA+ → Network). "
    "**E8 (GutSMASH) and E9 (MICOM) are not available for user-submitted data** — "
    "they require assembled MAGs and AGORA103 model coverage respectively."
)

# ── File uploaders ────────────────────────────────────────────────────────────
st.subheader("1 — Upload files")
col1, col2, col3 = st.columns(3)

with col1:
    st.markdown("**Species / Microbiome table**")
    species_file = st.file_uploader(
        "TSV or CSV — rows = taxa, columns = samples (or transposed)",
        type=["tsv", "csv", "txt"],
        key="species_upload",
    )
    if species_file:
        st.caption("Accepted formats: MetaPhlAn3/4 profile, generic OTU/ASV table")

with col2:
    st.markdown("**Metabolomics table**")
    mtb_file = st.file_uploader(
        "TSV or CSV — feature IDs (KEGG/HMDB/names) as rows or columns",
        type=["tsv", "csv", "txt"],
        key="mtb_upload",
    )
    if mtb_file:
        st.caption("Feature IDs: KEGG (C12345), HMDB IDs, or compound names")

with col3:
    st.markdown("**Sample metadata**")
    meta_file = st.file_uploader(
        "TSV or CSV — must include sample_id column + at least one grouping column",
        type=["tsv", "csv", "txt"],
        key="meta_upload",
    )

# ── Format preview & validation ───────────────────────────────────────────────
if species_file and mtb_file and meta_file:
    try:
        from pipeline.data_intake import (
            detect_and_load_species,
            detect_and_load_metabolomics,
            detect_and_load_metadata,
            validate_overlap,
        )

        with st.spinner("Detecting file formats…"):
            species_file.seek(0)
            mtb_file.seek(0)
            meta_file.seek(0)

            species_df, spe_fmt = detect_and_load_species(species_file.read())
            mtb_df, mtb_fmt, kegg_map = detect_and_load_metabolomics(mtb_file.read())
            meta_df = detect_and_load_metadata(meta_file.read())

        st.subheader("2 — Format preview")
        col1, col2, col3 = st.columns(3)
        with col1:
            st.caption(f"Format: {spe_fmt}")
            st.caption(f"{species_df.shape[0]} samples × {species_df.shape[1]} species")
            st.dataframe(species_df.head(3), use_container_width=True)
        with col2:
            st.caption(f"Format: {mtb_fmt}")
            st.caption(f"{mtb_df.shape[0]} samples × {mtb_df.shape[1]} metabolites")
            st.dataframe(mtb_df.head(3), use_container_width=True)
        with col3:
            st.caption(f"{meta_df.shape[0]} samples × {meta_df.shape[1]} columns")
            st.dataframe(meta_df.head(3), use_container_width=True)

        # Overlap validation
        overlap = validate_overlap(species_df, mtb_df, meta_df)
        n_common = overlap["n_common_samples"]
        for w in overlap["warnings"]:
            st.warning(w)
        st.success(f"{n_common} samples shared across all three files.")

        if n_common < 10:
            st.error("Too few shared samples to run the pipeline reliably (< 10).")
            st.stop()

        # ── Configuration ─────────────────────────────────────────────────────
        st.subheader("3 — Configure parameters")
        col1, col2, col3 = st.columns(3)
        with col1:
            group_col = st.selectbox(
                "Grouping / condition column",
                options=[c for c in meta_df.columns],
            )
        with col2:
            prev_spe = st.slider("Species prevalence threshold", 0.05, 0.30, 0.10, 0.01)
            prev_mtb = st.slider("Metabolite prevalence threshold", 0.05, 0.40, 0.15, 0.01)
        with col3:
            n_targets = st.slider("Max ML targets", 10, 50, 45, 5)
            mofa_factors = st.slider("MOFA+ latent factors", 5, 20, 10, 1)

        # Show group distribution
        if group_col:
            group_counts = meta_df[group_col].value_counts()
            st.caption("Group distribution: " + ", ".join(
                f"{g}={n}" for g, n in group_counts.items()
            ))

        # ── Run button ────────────────────────────────────────────────────────
        st.subheader("4 — Run pipeline")
        st.markdown(
            "**Estimated runtime:** 5–30 min (depends on dataset size and hardware). "
            "The ML step (step 3) is the longest."
        )

        if "pipeline_status" not in st.session_state:
            st.session_state.pipeline_status = {}

        status = st.session_state.pipeline_status
        run_dir = Path("results") / "user_run"

        if not status.get("running") and not status.get("done"):
            if st.button("Run pipeline", type="primary", use_container_width=True):
                from pipeline.data_intake import align_to_common_samples
                from pipeline.runner import run_in_background

                common = overlap["common_samples"]
                spe_aligned, mtb_aligned, meta_aligned = align_to_common_samples(
                    species_df, mtb_df, meta_df, common
                )
                run_dir.mkdir(parents=True, exist_ok=True)

                run_in_background(
                    spe_aligned, mtb_aligned, meta_aligned,
                    group_col=group_col,
                    output_dir=run_dir,
                    status_dict=st.session_state.pipeline_status,
                    prev_spe=prev_spe,
                    prev_mtb=prev_mtb,
                    n_ml_targets=n_targets,
                    mofa_factors=mofa_factors,
                )
                st.rerun()

        # ── Progress display ──────────────────────────────────────────────────
        if status.get("running") or status.get("done"):
            pct  = status.get("pct", 0)
            name = status.get("step_name", "")
            err  = status.get("error")
            done = status.get("done", False)

            if err:
                st.error(f"Pipeline failed:\n\n```\n{err}\n```")
            elif done:
                st.success("Pipeline complete! Navigate to pages 1–5 to explore your results.")
                st.session_state.run_output_dir = run_dir
                if st.button("Reset — upload new data"):
                    st.session_state.pipeline_status = {}
                    st.rerun()
            else:
                st.progress(pct / 100, text=f"Step {status.get('step', 0)}/5 — {name}")
                time.sleep(2)
                st.rerun()

    except Exception as e:
        st.error(f"Error reading files: {e}")
        st.stop()

else:
    st.info("Upload all three files to continue.")
