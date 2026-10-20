"""Executable output contract: what each pipeline notebook must leave in Results/.

CONTRACTS.md section 6 describes these outputs in prose; this module is the version a
machine checks. run_pipeline.py uses it three ways:

  * before a notebook runs, every file matching its declared outputs is moved to the
    run's archive, so whatever exists afterwards was written by THIS execution;
  * after the notebook runs, verify() checks each declared output: it exists, is newer
    than the notebook's start, a CSV has a header row, the required columns and enough
    rows, a pickle holds the required keys. Any violation turns OK into FAIL;
  * before a previous result is reused, verify() runs against it again, so a deleted
    or hand-edited output forces a re-run instead of being served.

DOC_FACTS ties numbers quoted in CONTRACTS.md to the tables that produce them, so the
document cannot drift from the results. checks.py E8 keeps the reverse direction
honest: every output CONTRACTS.md names is declared here or listed in NOT_PRODUCED.

Why this exists: the reproducibility review found CONTRACTS.md listing outputs that no
code wrote (pair_evidence_nb08_v2.csv, NB05's transfer tables, NB06's PostSurgery
projection), NB09 "OK" while documented to raise without one of them, and two
2-byte CSVs with no header that pandas cannot read. Every run had reported success.

Paths are relative to Results/. A '*' declares a family (one file per cohort, target or
comparison); min_files sets how many it must match. min_rows=0 and min_files=0 are
allowed only with a `note` saying when the output is legitimately empty or absent.
"""
from __future__ import annotations

import csv
import fnmatch
import pickle
import re
import shutil
from dataclasses import dataclass, field
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent / "Results"
CONTRACTS_MD = HERE / "CONTRACTS.md"
OUTPUT_DIRS = ("tables", "intermediate", "figures")

# Filesystem timestamps vs time.time(): allow this much clock slack.
MTIME_SLACK_S = 2.0

csv.field_size_limit(1 << 30)


@dataclass(frozen=True)
class Out:
    path: str                                   # relative to Results/; '*' = family
    cols: tuple[str, ...] = ()                  # CSV: required header columns
    keys: tuple[str, ...] = ()                  # pickle: required top-level keys
    min_rows: int = 1                           # CSV: data rows required per file
    min_files: int = 1                          # files the path must match
    values: dict = field(default_factory=dict)  # CSV: {column: values that must occur}
    allowed: dict = field(default_factory=dict) # CSV: {column: the only values permitted}
    note: str = ""                              # why min_rows/min_files may be 0

    def __post_init__(self):
        if (self.min_rows == 0 or self.min_files == 0) and not self.note:
            raise ValueError(f"{self.path}: an output allowed to be empty or absent "
                             "must say when and why (note=...)")


VALIDATION = ("SINHA_CRC_2016", "Kim_adenomas_2020")

# Selected-hit tables (significance-filtered) can be empty in a valid run; test tables
# (one row per feature or pair tested) cannot.
_HITS = "a significance-filtered hit table; empty is a valid result"
_NO_GENOMIC = ("needs a genomic reference index, which is unavailable "
               "(utils.AGORA2_ENABLED is False; gutSMASH pending)")
_MOFA_ONLY = {"Method": ("MOFA+",)}

CONTRACT_OUTPUTS: dict[str, tuple[Out, ...]] = {
    "01_preprocessing.ipynb": (
        Out("intermediate/preprocessed_data_v2.pkl", keys=(
            "transformed_species", "transformed_mtb", "filtered_species", "filtered_mtb",
            "harmonized_meta", "polyamine_columns", "polyamine_proxies", "mtb_kegg_map",
            "da_results_spc", "da_results_mtb", "mtb_prescaled", "mtb_normalisation",
            "excluded_study_groups", "sensitivity_groups", "stage_group_order",
            "sensitivity_stage_groups", "datasets", "dataset_discovery",
            "datasets_validation")),
        Out("tables/polyamine_cohort_coverage_v2.csv",
            cols=("KEGG", "Name", "YACHIDA", "SINHA", "Kim")),
        Out("figures/qc/*_pca_v2.png", min_files=9),       # 3 cohorts x 3 views
        Out("figures/qc/mtb_total_signal_by_group_v2.png"),
        Out("figures/qc/species_shannon_by_group_v2.png"),
        Out("figures/qc/species_detection_rates_v2.png"),
        Out("figures/qc/rda_species_metabolites_v2.png"),
        Out("figures/differential/*_da_subplots_v2.png", min_files=3),
    ),
    "01b_yachida_significant_features.ipynb": (
        Out("intermediate/yachida_significant_subset_v2.pkl", keys=(
            "species_filtered", "metabolites_filtered", "samples_meta",
            "sig_species_union", "sig_metabolites_union", "thresholds",
            "comparisons_primary", "sensitivity_da")),
        Out("tables/YACHIDA_v2_significant_species_per_comparison.csv",
            cols=("Feature", "Comparison", "QValue", "CliffsD", "Direction"),
            min_rows=0, note=_HITS),
        Out("tables/YACHIDA_v2_significant_metabolites_per_comparison.csv",
            cols=("Feature", "Comparison", "QValue", "CliffsD", "Direction"),
            min_rows=0, note=_HITS),
        Out("tables/YACHIDA_v2_significant_species_union.csv",
            cols=("Feature", "N_comparisons"), min_rows=0, note=_HITS),
        Out("tables/YACHIDA_v2_significant_metabolites_union.csv",
            cols=("Feature", "N_comparisons"), min_rows=0, note=_HITS),
        Out("tables/YACHIDA_v2_significant_samples.csv", cols=("Sample", "Stage.Group")),
        Out("tables/YACHIDA_v2_sensitivity_healthy_vs_postsurgery_*.csv", min_files=2,
            cols=("Feature", "PValue", "QValue", "CliffsD")),
        Out("figures/differential/YACHIDA_v2_*_annotated_volcano.png", min_files=6),
        Out("figures/differential/YACHIDA_v2_significant_*_heatmap.png", min_files=2),
    ),
    "02_analysis.ipynb": (
        Out("intermediate/preprocessed_data_reduced_v2.pkl", keys=(
            "transformed_species", "transformed_mtb", "harmonized_meta",
            "sensitivity_species", "sensitivity_mtb", "sensitivity_meta",
            "polyamine_columns", "polyamine_proxies", "mtb_kegg_map", "undetected_mtb",
            "mtb_prescaled", "datasets", "dataset_discovery", "datasets_validation",
            "stage_group_order", "sample_scope", "feature_reduction")),
        Out("tables/differential_abundance_v2.csv",
            cols=("Dataset", "Comparison", "Scheme", "Assay", "Feature", "MeanDiff",
                  "PValue", "CliffsD", "QValue", "Coef_adj", "QValue_adj", "Agreement",
                  "Robust"),
            # CONTRACTS.md section 2: the grouped, confounder-adjusted scheme only, which
            # exists for Yachida alone (the validation cohorts have no Stage.Group).
            values={"Dataset": ("YACHIDA-CRC-2019",), "Scheme": ("grouped",),
                    "Assay": ("species", "metabolites"),
                    "Comparison": ("Healthy_vs_Early_CRC", "Healthy_vs_Advanced_CRC",
                                   "Early_CRC_vs_Advanced_CRC")}),
        Out("tables/da_robust_hit_summary_v2.csv",
            cols=("Comparison", "Type", "N_tested", "N_robust")),
        Out("tables/polyamine_DA_summary_v2.csv",
            cols=("Dataset", "Comparison", "Polyamine", "QValue")),
        Out("tables/*_polyamine_partial_correlations_v2.csv", min_files=3,
            cols=("Species", "Metabolite", "Polyamine", "Rho_partial", "PValue_partial",
                  "QValue_partial", "N", "Covariates")),
        Out("tables/*_genus_polyamine_correlations_v2.csv", min_files=3,
            cols=("Genus", "Metabolite", "Rho", "PValue", "QValue", "Parent_KEGG")),
        Out("tables/consistent_polyamine_correlations_genus_v2.csv",
            cols=("Genus", "Metabolite_key", "n_datasets", "Same_direction"),
            min_rows=0, note=_HITS),
        Out("tables/replication_genus_metabolome_v2.csv",
            cols=("Genus", "Compound", "Rho_discovery", "Q_discovery", "Canonical",
                  "Replication_tier")
                 + tuple(f"{c}_{ds}" for ds in VALIDATION
                         for c in ("In", "Rho", "Q", "Replicated", "Contradicted"))),
        Out("tables/polyamine_hub_species_v2.csv",
            cols=("Dataset", "Species", "N_polyamine_connections"), min_rows=0, note=_HITS),
        Out("tables/stratified_hub_species_v2.csv",
            cols=("Group", "Species", "N_polyamine_connections"), min_rows=0, note=_HITS),
        Out("tables/CRC_progression_ordinal_correlations_v2.csv",
            cols=("Species", "Metabolite", "Rho_spc_mtb", "QVal_spc_mtb",
                  "Rho_species_stage", "Rho_metabolite_stage")),
        Out("tables/CRC_progression_top_pairs_v2.csv", cols=("Species", "Metabolite")),
        Out("tables/CRC_StageI_II_vs_III_IV_*_DA_v2.csv", min_files=2,
            cols=("Feature", "PValue", "CliffsD", "QValue")),
        Out("tables/CRC_binary_stage_correlations_v2.csv",
            cols=("Species", "Metabolite", "Rho", "PValue", "QValue")),
        Out("tables/CRC_progression_species_overlap_10A_10B_v2.csv",
            cols=("Species_overlap",), min_rows=0, note=_HITS),
        Out("tables/CRC_progression_metabolite_overlap_10A_10B_v2.csv",
            cols=("Metabolite_overlap",), min_rows=0, note=_HITS),
        Out("figures/qc/residual_diagnostics_*_v2.png", min_files=3),
        Out("figures/qc/residual_diagnostics_*_v2.pdf", min_files=3),
        Out("figures/differential/*_mtb_volcano_v2.png", min_files=3),
        Out("figures/differential/*_spc_volcano_v2.png", min_files=3),
        Out("figures/differential/da_results_dashboard_v2.png"),
        Out("figures/correlation/*_polyamine_heatmap_v2.png", min_files=3),
        Out("figures/correlation/cross_cohort_consistency_v2.png"),
        Out("figures/network/*_polyamine_network_v2.png", min_files=3),
        Out("figures/network/hub_network_*_v2.png", min_files=1),
        Out("figures/network/stage_group_dashboard_v2.png"),
        Out("figures/progression/ordinal_progression_scatter_v2.png"),
        Out("figures/progression/binary_StageI_II_vs_III_IV_*_volcano_v2.png", min_files=2),
        Out("figures/progression/trend_*_v2.png", min_files=0,
            note="one panel per top species-metabolite stage-trend pair; none if no "
                 "pair qualifies"),
    ),
    "03_ml_benchmarking.ipynb": (
        Out("intermediate/ml_results_v2.pkl", keys=(
            "top15_early", "top15_advanced", "benchmark_df", "pa_producers",
            "oof_shap_store", "polyamine_targets", "n_perm")),
        Out("tables/ml_benchmark_dysregulated_metabolites_v2.csv",
            cols=("Metabolite", "Classifier", "Target_type", "Target_definition",
                  "AUC_nested", "AUC_fixed", "N", "Perm_statistic", "Perm_pvalue",
                  "Perm_qvalue", "N_perm", "R2_confounders", "AUC_confounders_raw",
                  "AUC_confounders_adjusted")),
        Out("tables/polyamine_producer_candidates_v2.csv",
            cols=("Polyamine", "Species", "Mean_SHAP_OOF", "Direction", "Rank",
                  "Rank_stable", "Model_AUC_nested", "Model_perm_q"),
            min_rows=0, note="no producers when no polyamine model is benchmarked "
                             "(class imbalance) or none has a stable positive species"),
        Out("figures/ml/top15_metabolites_early_CRC_v2.png"),
        Out("figures/ml/top15_metabolites_advanced_CRC_v2.png"),
        Out("figures/ml/metabolite_dysregulation_comparison_v2.png"),
        Out("figures/ml/benchmark_auc_bars_v2.png"),
        Out("figures/ml/polyamine_producer_heatmap_v2.png", min_files=0,
            note="drawn only when producer candidates exist"),
        Out("figures/ml/polyamine_*_oof_shap_producers_v2.png", min_files=0,
            note="one per benchmarked polyamine with producer candidates"),
    ),
    "05_two_stage_model.ipynb": (
        Out("tables/nb05_stage2_early_vs_advanced_v2.csv",
            cols=("Model", "AUC_oof", "N", "N_advanced")),
        Out("tables/nb05_external_transfer_v2.csv",
            cols=("Target_cohort", "N_shared_genera", "N_test", "AUC_within_Yachida",
                  "AUC_transfer", "CI_transfer_low", "CI_transfer_high",
                  "P_transfer_vs_chance", "AUC_transfer_presence_only",
                  "Train_zero_frac_harmonised", "Target_zero_frac"),
            values={"Target_cohort": VALIDATION}),
        Out("tables/nb05_severity_transfer_v2.csv",
            cols=("Train_cohort", "Test_cohort", "N_shared_genera", "Rho_within_train",
                  "Rho_transfer", "Perm_p", "N_permutations", "AUC_cohort_separability")),
        Out("tables/two_stage_model_results_v2.csv",
            cols=("Polyamine", "R_two_stage", "R_flat", "Delta_R")),
        Out("tables/two_stage_shap_*_v2.csv", cols=("Feature", "Mean_SHAP")),
        Out("figures/ml/nb05_external_transfer_roc_v2.png"),
        Out("figures/ml/nb05_severity_transfer_v2.png"),
        Out("figures/ml/nb05_severity_transfer_v2.pdf"),
        Out("figures/ml/two_stage_shap_*_v2.png"),
    ),
    # Method must be MOFA+ on every row: NB06's PCA fallback (opt-in only) labels its
    # outputs Method=PCA, and PCA axes are never MOFA evidence.
    "06_multiomics_factors.ipynb": (
        Out("tables/mofa_factor_loadings_v2.csv",
            cols=("Factor", "View", "Feature", "Weight", "Method"), allowed=_MOFA_ONLY),
        Out("tables/mofa_factor_scores_v2.csv",
            cols=("Sample", "Factor1", "Stage.Group", "Method"), allowed=_MOFA_ONLY),
        Out("tables/mofa_factor_stage_association_v2.csv",
            cols=("Factor", "KW_pval", "KW_qval", "Spearman_rho", "Spearman_qval", "Method"),
            allowed=_MOFA_ONLY),
        Out("tables/mofa_postsurgery_projection_v2.csv",
            cols=("Sample", "Factor1", "Projected", "Method"), allowed=_MOFA_ONLY),
        Out("figures/multiomics/mofa_variance_explained_v2.png"),
        Out("figures/multiomics/mofa_factor_scores_by_stage_v2.png"),
        Out("figures/multiomics/mofa_factor_biplot_v2.png"),
        Out("figures/multiomics/mofa_loading_heatmap_v2.png"),
        Out("figures/multiomics/mofa_postsurgery_projection_v2.png"),
    ),
    "08_pathway_analysis.ipynb": (
        Out("tables/r1_agora2_producer_enrichment_v2.csv",
            cols=("Species", "Polyamine", "Rank", "KEGG_ID", "Encodes_Enzyme",
                  "AGORA2_Model_Available"),
            min_rows=0, note="one row per NB03 positive-direction producer; none if NB03 "
                             "found none"),
        Out("tables/r1_agora2_enrichment_test_v2.csv",
            cols=("Polyamine", "KEGG_ID", "n_universe", "OddsRatio", "PValue", "QValue"),
            min_rows=0, note="no polyamine is testable: the test " + _NO_GENOMIC),
        Out("tables/r3_scfa_polyamine_enrichment_v2.csv",
            cols=("Test", "OddsRatio", "PValue", "QValue", "Significant"),
            min_rows=0, note="no producers to test"),
        Out("tables/r4_factor_oxidative_index_correlations_v2.csv",
            cols=("Factor", "Rho_oxidative_index", "PValue", "Q_PValue")),
        Out("tables/r4_blautia_oxidative_correlations_v2.csv",
            cols=("Species", "Rho_stage", "Pval_stage")),
        Out("tables/r5_factor_oral_invasion_correlations_v2.csv",
            cols=("Factor", "Rho_oral_invasion", "PValue", "Q_PValue")),
        Out("tables/r5_dicarboxylic_acid_correlations_v2.csv",
            cols=("Metabolite", "KEGG_ID", "Rho_stage", "P_stage")),
        Out("tables/r6_fusobacterium_polyamine_correlations_v2.csv",
            cols=("Species", "Rho_ornithine", "Rho_putrescine", "Rho_stage")),
        Out("tables/r8_polyamine_dicarb_overlap_v2.csv",
            cols=("Species", "Polyamine", "KEGG_ID", "Rho_polyamine", "QValue_polyamine",
                  "Dicarb_Metabolite", "Rho_dicarb", "QValue_dicarb"),
            min_rows=0, note=_HITS + " (always written: NB09 reads it)"),
        Out("figures/pathway/r4_oxidative_stress_index_by_stage_v2.png"),
        Out("figures/pathway/r5_oral_invasion_dicarb_v2.png"),
        Out("figures/pathway/r8_dicarboxylic_acid_index_v2.png"),
    ),
    "09_pathway_network.ipynb": (
        Out("tables/d_glasso_partial_correlation_edges_v2.csv",
            cols=("Species", "Metabolite", "PartialCorr", "P_maxT_FWER", "Significant",
                  "Alpha", "Alpha_at_grid_edge", "N_permutations", "N_perm_failed"),
            min_rows=0, note="GLASSO may select no species-metabolite edge"),
        Out("tables/final_evidence_synthesis_v2.csv",
            cols=("Species", "Polyamine", "KEGG_ID", "Rho_partial", "QValue_partial",
                  "N_Sources_Available", "N_Sources_Supporting", "Evidence_Level")),
        Out("tables/c_pathway_mapping_v2.csv", cols=("Feature", "KEGG_ID", "Reactome_ID")),
        Out("tables/c_species_pathway_correlations_v2.csv",
            cols=("Pathway", "Species", "Rho", "PValue", "QValue")),
        Out("tables/g_guild_members_v2.csv", cols=("Guild", "Species"),
            min_rows=0, note="guild membership " + _NO_GENOMIC),
        Out("tables/g_pathway_enrichment_scores_v2.csv", cols=("Sample", "Stage.Group")),
        Out("tables/g_pathway_factor_correlations_v2.csv",
            cols=("Pathway", "Factor", "Rho", "PValue", "QValue", "Significant"),
            min_rows=0, note="no guild scores to correlate: guild membership " + _NO_GENOMIC),
        Out("tables/b_agora2_producer_enzyme_profile_v2.csv", min_files=0,
            note="Module B " + _NO_GENOMIC),
        Out("figures/pathway_network/d_glasso_network_v2.png"),
        Out("figures/pathway_network/c_species_pathway_heatmap_v2.png"),
        Out("figures/pathway_network/final_evidence_synthesis_heatmap_v2.png", min_files=0,
            note="drawn only when some pair has >= 2 supporting sources"),
        Out("figures/pathway_network/g_pathway_mofa_factor_heatmap_v2.png", min_files=0,
            note="drawn only when guild scores exist"),
        Out("figures/pathway_network/g_pathway_stage_boxplots_v2.png", min_files=0,
            note="drawn only when guild scores exist"),
    ),
}

# Outputs CONTRACTS.md names that no executed notebook produces. checks.py E8 fails on
# any other name the document mentions without a declaration above.
NOT_PRODUCED = {
    "pair_evidence_nb08_v2.csv": "design target of the NB08 rebuild; the executed NB08 "
                                 "does not write it (CONTRACTS.md NB08 status banner)",
    "analysis_results_v2.pkl": "retired: zero readers, duplicated the CSVs",
    "*_vs_*_DA_v2.csv": "retired per-comparison DA family, replaced by "
                        "differential_abundance_v2.csv (the Stage I-II vs III-IV "
                        "contrast is declared under NB02)",
    "nb10_*_v2.csv": "NB10 (MICOM) is not runnable in the pinned environment",
    "nb04_*": "NB04 is discarded",
}


@dataclass(frozen=True)
class DocFact:
    label: str
    pattern: str            # regex on CONTRACTS.md; group 1 = the quoted number
    table: str              # CSV relative to Results/
    column: str
    where: tuple = ()       # (column, value) selecting exactly one row


# Numbers CONTRACTS.md quotes that a run must reproduce. A quoted value matches when the
# produced value, rounded to the quoted precision, equals it.
DOC_FACTS: dict[str, tuple[DocFact, ...]] = {
    "05_two_stage_model.ipynb": (
        DocFact("stage 3 shared genera, Yachida-Sinha",
                r"stage 3 Yachida ∩ Sinha = (\d+)",
                "tables/nb05_external_transfer_v2.csv", "N_shared_genera",
                ("Target_cohort", "SINHA_CRC_2016")),
        DocFact("stage 3 shared genera, Yachida-Kim",
                r"stage 3 Yachida ∩ Sinha = \d+,\s+Yachida ∩ Kim = (\d+)",
                "tables/nb05_external_transfer_v2.csv", "N_shared_genera",
                ("Target_cohort", "Kim_adenomas_2020")),
        DocFact("stage 6 shared genera, Yachida-Kim",
                r"stage 6 \(adenomas admitted\) Yachida ∩ Kim = (\d+)",
                "tables/nb05_severity_transfer_v2.csv", "N_shared_genera"),
        DocFact("stage 3 transfer AUC, Sinha", r"Yachida → Sinha AUC\s+([\d.]+)",
                "tables/nb05_external_transfer_v2.csv", "AUC_transfer",
                ("Target_cohort", "SINHA_CRC_2016")),
        DocFact("stage 3 transfer p, Sinha", r"Yachida → Sinha AUC\s+[\d.]+\s+\(p = ([\d.]+)",
                "tables/nb05_external_transfer_v2.csv", "P_transfer_vs_chance",
                ("Target_cohort", "SINHA_CRC_2016")),
        DocFact("stage 3 presence-only AUC, Sinha", r"presence-only ([\d.]+)\)",
                "tables/nb05_external_transfer_v2.csv", "AUC_transfer_presence_only",
                ("Target_cohort", "SINHA_CRC_2016")),
        DocFact("stage 3 transfer AUC, Kim", r"Yachida → Kim AUC ([\d.]+)",
                "tables/nb05_external_transfer_v2.csv", "AUC_transfer",
                ("Target_cohort", "Kim_adenomas_2020")),
        DocFact("stage 3 transfer p, Kim", r"Yachida → Kim AUC [\d.]+ \(p = ([\d.]+)",
                "tables/nb05_external_transfer_v2.csv", "P_transfer_vs_chance",
                ("Target_cohort", "Kim_adenomas_2020")),
        DocFact("stage 6 within-Yachida rho", r"within-Yachida ρ ([+−-]?[\d.]+)",
                "tables/nb05_severity_transfer_v2.csv", "Rho_within_train"),
        DocFact("stage 6 transferred rho", r"transferred to Kim ρ ([+−-]?[\d.]+)",
                "tables/nb05_severity_transfer_v2.csv", "Rho_transfer"),
        DocFact("stage 6 permutation p", r"transferred to Kim ρ [+−-]?[\d.]+, permutation "
                                         r"p = ([\d.]+)",
                "tables/nb05_severity_transfer_v2.csv", "Perm_p"),
        DocFact("stage 6 cohort separability AUC", r"separability AUC ([\d.]+)",
                "tables/nb05_severity_transfer_v2.csv", "AUC_cohort_separability"),
    ),
}


class ContractError(RuntimeError):
    pass


# ---------------------------------------------------------------------------

def _matches(rel: str, pattern: str) -> bool:
    # fnmatch's '*' crosses '/', so anchor the directory part explicitly.
    return (rel.count("/") == pattern.count("/")) and fnmatch.fnmatchcase(rel, pattern)


def owners(rel: str) -> list[str]:
    """Notebooks whose contract declares this Results-relative path."""
    return [nb for nb, outs in CONTRACT_OUTPUTS.items()
            if any(_matches(rel, o.path) for o in outs)]


def _files(o: Out, results: Path) -> list[Path]:
    if "*" not in o.path:
        p = results / o.path
        return [p] if p.is_file() else []
    return sorted(p for p in results.glob(o.path) if p.is_file()
                  and _matches(p.relative_to(results).as_posix(), o.path))


def _csv_shape(path: Path, want: dict) -> tuple[list[str], int, dict]:
    """(header, data rows, {column: set of values seen}) streamed with the csv module."""
    seen = {c: set() for c in want}
    with open(path, newline="", encoding="utf-8") as fh:
        rd = csv.reader(fh)
        header = next(rd, None) or []
        idx = {c: header.index(c) for c in want if c in header}
        n = 0
        for row in rd:
            n += 1
            for c, i in idx.items():
                if i < len(row):
                    seen[c].add(row[i])
    return header, n, seen


def _check_file(path: Path, o: Out, rel: str, since: float | None,
                until: float | None) -> list[str]:
    st = path.stat()
    out = []
    if since is not None and st.st_mtime < since - MTIME_SLACK_S:
        out.append(f"{rel}: stale -- last written before this notebook started")
    if until is not None and st.st_mtime > until + MTIME_SLACK_S:
        out.append(f"{rel}: modified after its notebook finished (edited by hand, or "
                   "written by another notebook)")
    if st.st_size == 0:
        return out + [f"{rel}: empty file"]
    suffix = path.suffix.lower()
    if suffix == ".csv":
        header, n, seen = _csv_shape(path, {**o.values, **o.allowed})
        if not header or all(not h.strip() for h in header):
            return out + [f"{rel}: no header row (a column-less DataFrame was written; "
                          "build empty tables with explicit columns)"]
        missing = [c for c in o.cols if c not in header]
        if missing:
            out.append(f"{rel}: missing column(s) {missing}")
        if n < o.min_rows:
            out.append(f"{rel}: {n} data row(s), contract requires >= {o.min_rows}")
        for col, required in o.values.items():
            absent = [v for v in required if v not in seen.get(col, set())]
            if absent:
                out.append(f"{rel}: column {col!r} lacks {absent}")
        for col, permitted in o.allowed.items():
            extra = sorted(seen.get(col, set()) - set(permitted))
            if extra:
                out.append(f"{rel}: column {col!r} holds {extra}; only {list(permitted)} "
                           "is permitted")
    elif suffix == ".pkl":
        if o.keys:
            with open(path, "rb") as fh:
                obj = pickle.load(fh)          # the pipeline's own output, not untrusted
            if not isinstance(obj, dict):
                out.append(f"{rel}: expected a dict pickle, got {type(obj).__name__}")
            else:
                missing = [k for k in o.keys if k not in obj]
                if missing:
                    out.append(f"{rel}: missing key(s) {missing}")
    elif suffix == ".png":
        with open(path, "rb") as fh:
            if fh.read(8) != b"\x89PNG\r\n\x1a\n":
                out.append(f"{rel}: not a PNG file")
    elif suffix == ".pdf":
        with open(path, "rb") as fh:
            if fh.read(5) != b"%PDF-":
                out.append(f"{rel}: not a PDF file")
    return out


def verify(nb: str, since: float | None = None, until: float | None = None,
           results: Path = RESULTS) -> list[str]:
    """Contract violations for one notebook's outputs ([] = satisfied). `since`/`until`
    bound the execution that is supposed to have written them."""
    problems = []
    for o in CONTRACT_OUTPUTS.get(nb, ()):
        files = _files(o, results)
        if len(files) < o.min_files:
            problems.append(f"{o.path}: missing" if "*" not in o.path else
                            f"{o.path}: {len(files)} file(s), contract requires "
                            f">= {o.min_files}")
            continue
        for f in files:
            problems += _check_file(f, o, f.relative_to(results).as_posix(), since, until)
    return problems


def _quoted_precision(text: str) -> int:
    return len(text.split(".", 1)[1]) if "." in text else 0


def verify_doc_facts(nb: str, results: Path = RESULTS,
                     contracts_md: Path = CONTRACTS_MD) -> list[str]:
    """Numbers CONTRACTS.md quotes for this notebook that its outputs do not reproduce."""
    facts = DOC_FACTS.get(nb, ())
    if not facts:
        return []
    import pandas as pd

    text = contracts_md.read_text(encoding="utf-8")
    problems, tables = [], {}
    for f in facts:
        m = re.search(f.pattern, text)
        if not m:
            problems.append(f"CONTRACTS.md no longer quotes {f.label} "
                            f"(pattern {f.pattern!r}); update DOC_FACTS with the text")
            continue
        quoted = m.group(1).replace("−", "-").lstrip("+")
        if f.table not in tables:
            tables[f.table] = pd.read_csv(results / f.table)
        df = tables[f.table]
        if f.where:
            df = df[df[f.where[0]].astype(str) == str(f.where[1])]
        if len(df) != 1:
            problems.append(f"{f.label}: {f.table} has {len(df)} matching rows, expected 1")
            continue
        value = float(df[f.column].iloc[0])
        nd = _quoted_precision(quoted)
        if round(value, nd) != round(float(quoted), nd):
            problems.append(f"{f.label}: CONTRACTS.md quotes {m.group(1)}, this run produced "
                            f"{value:.{max(nd, 3)}f} ({f.table}:{f.column}) -- update the "
                            "document")
    return problems


def clear_outputs(nb: str, archive_dir: Path, results: Path = RESULTS) -> int:
    """Move every existing file this notebook's contract declares into archive_dir, so
    what exists after the notebook runs was written by that run. Refuses if a file is
    also declared by another notebook (an ambiguous contract)."""
    moved = 0
    for o in CONTRACT_OUTPUTS.get(nb, ()):
        for f in _files(o, results):
            rel = f.relative_to(results).as_posix()
            who = owners(rel)
            if who != [nb]:
                raise ContractError(f"{rel} is declared by {who}; contract patterns must "
                                    "not overlap")
            dest = archive_dir / rel
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(f), str(dest))
            moved += 1
    return moved


def undeclared(nb: str, since: float, until: float, results: Path = RESULTS) -> list[str]:
    """Files written during this notebook's window that its contract does not declare."""
    out = []
    for sub in OUTPUT_DIRS:
        base = results / sub
        if not base.exists():
            continue
        for p in base.rglob("*"):
            if not p.is_file():
                continue
            mt = p.stat().st_mtime
            if since - MTIME_SLACK_S <= mt <= until + MTIME_SLACK_S:
                rel = p.relative_to(results).as_posix()
                if nb not in owners(rel):
                    out.append(rel)
    return sorted(out)
