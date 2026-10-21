"""
Shared utility functions for the v2 microbiome-metabolome analysis pipeline.

v2 vs the `exploratory/` + `Python Patches/utils.py` pipeline it replaces:
see v2/CHANGES.md for the full list of fixes. The short version:
  - RESULTS_DIR now actually points at the shared results tree (the old
    module pointed at a near-empty stub directory).
  - The CRC stage ordering is fixed at the source (see "STAGE ENCODING"
    below) instead of being wrong in utils.py and then hand-re-broken
    locally in half the notebooks that used it.
  - differential_abundance() no longer takes log2() of an already-
    log/CLR-transformed ratio (which was frequently NaN); it reports a
    linear-space MeanDiff instead.
  - progression_correlation_analysis() no longer FDR-corrects a table
    that silently duplicates each unique test dozens to hundreds of times.
  - Confounder handling treats "no tumor" as its own category instead of
    something MICE imputes as if it were missing data.
  - Adds the AGORA2 species-encodes-enzyme lookup (ported from
    CRC_Gut_Analysis_Git/utils.py) to replace genus-level KEGG lookups
    used as production/enrichment evidence.
"""

import os
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import entropy, mannwhitneyu, spearmanr, rankdata
from statsmodels.stats.multitest import multipletests
from sklearn.decomposition import PCA
from sklearn.exceptions import ConvergenceWarning

# ============================================================================
# NUMERICAL HEALTH
# ============================================================================

# What fit_strict raises when a fit is not numerically sound. Catch this tuple where a
# failure has a defined treatment (e.g. counted as a permutation exceedance).
NUMERICAL_FAILURES = (ConvergenceWarning, RuntimeWarning, FloatingPointError)


def fit_strict(estimator, *args, **kwargs):
    """Fit `estimator` with convergence and numerical warnings promoted to errors.

    A warning is not an exception. NB09's GLASSO permutation loop caught exceptions
    only, so 37 of 100 refits that raised ConvergenceWarning were counted as
    converged and the cell printed "100/100 converged". Use this for every fit whose
    convergence a reported result depends on, and catch NUMERICAL_FAILURES wherever
    a failure has a defined, conservative treatment. Returns the fitted estimator.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("error", ConvergenceWarning)
        warnings.simplefilter("error", RuntimeWarning)
        return estimator.fit(*args, **kwargs)


# ============================================================================
# PATHS
# ============================================================================

# Resolved from this file's location (<project root>/v2/utils.py), so the same
# code works wherever the tree is checked out. Set KI_PROJECT_ROOT to override.
#
# This was parents[2] while the pipeline lived at <root>/Notebooks/v2/. Moving it
# to <root>/v2/ made parents[2] the DRIVE root, so DATA_DIR/RESULTS_DIR silently
# pointed at nonexistent directories and every loader returned empty tables --
# hence the explicit check below instead of 32 per-file "not found" warnings.
_ENV_ROOT = os.environ.get("KI_PROJECT_ROOT")
PROJECT_ROOT = (Path(_ENV_ROOT).expanduser().resolve() if _ENV_ROOT
                else Path(__file__).resolve().parents[1])
DATA_DIR    = PROJECT_ROOT / "Data"
RESULTS_DIR = PROJECT_ROOT / "Results"

if not DATA_DIR.is_dir():
    warnings.warn(
        f"DATA_DIR does not exist: {DATA_DIR} (PROJECT_ROOT={PROJECT_ROOT}). "
        "load_all_datasets() will return nothing but None tables. If the project "
        "tree moved, set the KI_PROJECT_ROOT environment variable or correct "
        "PROJECT_ROOT in utils.py.",
        stacklevel=2,
    )

# ============================================================================
# DATASETS
# ============================================================================

# The IBD/IBS/ESRD cohorts (FRANZOSA-IBD-2019, iHMP-IBDMDB-2019, MARS-IBS-2020,
# WANG_ESRD_2020) were REMOVED from the project on 2026-10-07: nothing downstream
# of NB01 used them, and loading all 8 raw tables exceeded the server's 1.5 GB
# cap. Their files are no longer on disk, so they are gone from DATASETS_ALL,
# COMPARISONS, POLYAMINE_DATASETS and STUDY_GROUP_MAP rather than left as
# entries that resolve to missing files.
DATASETS_ALL = [
    "ERAWIJANTARI-GASTRIC-CANCER-2020",
    "YACHIDA-CRC-2019",
    "SINHA_CRC_2016",
    "Kim_adenomas_2020",
]
DATASETS = DATASETS_ALL

# Yachida (CRC) + Sinha (CRC) are the two cohorts used for primary/
# confirmatory analysis in NB02 onward. (The old exploratory pipeline's
# markdown wrongly called the second cohort "Erawijantari" in several
# places — Erawijantari is gastric cancer, not CRC, and is not in this set.)
DATASETS_ANALYSIS = [
    "YACHIDA-CRC-2019",
    "SINHA_CRC_2016",
    "Kim_adenomas_2020",
]

# Discovery cohort vs external validation cohorts. Yachida is the only cohort with
# CRC staging, metabolite-level shotgun species data and full confounders; Sinha and
# Kim are genus-level 16S with pre-scaled (Metabolon-style) metabolites, used only to
# test whether a Yachida finding reproduces in an independent population.
DATASET_DISCOVERY  = "YACHIDA-CRC-2019"
DATASETS_VALIDATION = ["SINHA_CRC_2016", "Kim_adenomas_2020"]

# Kept as an empty list, not deleted: NB01 still imports it and tags figures with
# `if ds in DATASETS_IBD`. Those cohorts are gone, so every such test is now False
# and the exploratory tagging is simply inert.
DATASETS_IBD: list[str] = []

FILE_TYPES = ["metadata", "mtb", "mtb.map", "species"]

POLYAMINE_DATASETS = [
    "ERAWIJANTARI-GASTRIC-CANCER-2020",
    "YACHIDA-CRC-2019",
    "SINHA_CRC_2016",
    "Kim_adenomas_2020",
]
POLYAMINE_DATASETS_ANALYSIS = ["YACHIDA-CRC-2019", "SINHA_CRC_2016", "Kim_adenomas_2020"]

# Which polyamines each cohort actually measures (checked against each mtb.map):
#   Yachida : all 7 (putrescine, spermidine, spermine, agmatine, cadaverine,
#             ornithine, arginine) + N-acetyl derivatives
#   Sinha   : putrescine, cadaverine, ornithine, arginine (no spermidine/spermine/
#             agmatine); N-acetylputrescine present
#   Kim     : spermidine, ornithine, arginine directly; putrescine and cadaverine
#             ONLY as the N-acetyl derivatives (C02714 N-acetylputrescine,
#             N-acetyl-cadaverine) -- see POLYAMINE_PROXY_KEGG.
#
# Acetylated derivatives are NOT the parent compound: N-acetylation is a host/
# bacterial inactivation step, so the derivative tracks the parent only loosely.
# They are therefore declared as explicit, pre-specified PROXIES, used only when a
# cohort lacks the parent, and every result built on one is labelled as such.
POLYAMINE_PROXY_KEGG = {
    "C00134": ["C02714"],   # Putrescine  <- N-acetylputrescine
    "C00315": ["C00612"],   # Spermidine  <- N1-acetylspermidine
    "C00750": ["C02567"],   # Spermine    <- N1-acetylspermine
    "C00077": ["C00437"],   # Ornithine   <- N-acetylornithine
}

# Some platforms report an acetyl derivative that has no KEGG ID at all (Kim's
# "N-acetyl-cadaverine" is annotated with neither KEGG nor HMDB). These names are
# matched case-insensitively, after punctuation is stripped, and ONLY for a parent
# with no direct measurement and no KEGG proxy. Listed explicitly here so the
# proxy set stays pre-specified rather than discovered per cohort.
POLYAMINE_PROXY_NAMES = {
    "C01672": ["n acetyl cadaverine", "n acetylcadaverine"],   # Cadaverine
    "C00179": ["n acetylagmatine", "n acetyl agmatine"],        # Agmatine
}

POLYAMINE_KEGG = {
    "C00134": "Putrescine",
    "C00315": "Spermidine",
    "C00750": "Spermine",
    "C00179": "Agmatine",
    "C01672": "Cadaverine",
    "C00077": "Ornithine",
    "C00062": "Arginine",
}

# Curated EC numbers for reactions that PRODUCE each polyamine KEGG compound
# — used as an override by species_encodes_enzyme() below, since it's more
# specific than "any EC-annotated reaction that makes it".
#
# v2 fix: the previous map listed enzymes that CONSUME the compound for
# ornithine (ODC 4.1.1.17 decarboxylates ornithine -> putrescine) and
# arginine (ADC 4.1.1.19 / agmatinase 3.5.3.11 act on arginine/agmatine),
# and gave putrescine LDC/ADC (whose products are cadaverine/agmatine).
# Every EC below has the compound as its product.
POLYAMINE_EC: dict[str, list[str]] = {
    "C00134": ["4.1.1.17", "3.5.3.11", "3.5.1.53"],  # Putrescine: ODC, agmatinase, N-carbamoylputrescine amidohydrolase
    "C00315": ["2.5.1.16", "4.1.1.96"],              # Spermidine: spermidine synthase, carboxynorspermidine decarboxylase
    "C00750": ["2.5.1.22"],                          # Spermine: spermine synthase
    "C00179": ["4.1.1.19"],                          # Agmatine: arginine decarboxylase (ADC)
    "C01672": ["4.1.1.18"],                          # Cadaverine: lysine decarboxylase (LDC)
    "C00077": ["3.5.3.1", "3.5.1.16", "2.3.1.35"],   # Ornithine: arginase, acetylornithine deacetylase, ornithine acetyltransferase
    "C00062": ["4.3.2.1"],                           # Arginine: argininosuccinate lyase
}

# Polyamines proper (excluding the near-universal amino-acid precursors
# ornithine/arginine, whose synthesis ECs are present in ~90% of AGORA2
# species and therefore carry no pathway-membership information).
POLYAMINE_PATHWAY_KEGG = ["C00134", "C00315", "C00750", "C00179", "C01672"]

# EC-defined functional guilds (used instead of curated genus lists, which
# in the old pipeline were copied from the same cohort's SHAP top genera).
# Note 2.8.3.8 (butyryl-CoA:acetate CoA-transferase) also acts on
# acetoacetate, so it is a strong but not exclusive butyrate marker.
SCFA_EC: dict[str, list[str]] = {
    "butyrate":   ["2.8.3.8", "2.7.2.7"],   # butyryl-CoA:acetate CoA-transferase, butyrate kinase
    "propionate": ["2.8.3.1"],              # propionate CoA-transferase
}
PURINE_CATAB_EC = ["1.17.1.4", "3.5.2.5"]   # xanthine dehydrogenase, allantoinase

COMPARISONS = {
    "ERAWIJANTARI-GASTRIC-CANCER-2020": [("Healthy", "Gastrectomy")],
    "SINHA_CRC_2016":                  [("Healthy", "CRC")],
    "Kim_adenomas_2020":               [("Healthy", "CRC")],
    "YACHIDA-CRC-2019": [
        ("Healthy", "Stage_I_II"),
        ("Healthy", "Stage_III_IV"),
        ("Healthy", "Stage_0"),
    ],
}

# ============================================================================
# STAGE ENCODING (single source of truth — see CHANGES.md finding #4)
# ----------------------------------------------------------------------
# Yachida et al. 2019 (DOI 10.1038/s41591-019-0458-7) study groups, as
# confirmed by the user who owns this analysis:
#   Healthy      n=127  colonoscopy-confirmed controls
#   HS           n=30   "healthy after surgery" — patients who had colorectal
#                       cancer resected and are currently cancer-free
#   MP           n=40   multiple polypoid adenomas with low-grade dysplasia
#   Stage_0      n=27   intramucosal carcinoma / carcinoma in situ
#   Stage_I_II   n=69   invasive CRC, TNM I–II
#   Stage_III_IV n=54   invasive CRC, TNM III–IV
# Crosstabulating Study.Group against the raw `Stage` column confirms
# Healthy, HS and MP carry Stage == "-" (no current carcinoma).
#
# v2.1 group policy (user decision, 2026-10-01):
#   * MP is EXCLUDED from the study entirely (EXCLUDED_STUDY_GROUPS).
#   * HS is kept as its own group, "PostSurgery", and is NEVER pooled with
#     Healthy or with any CRC group. A resected gut is neither a clean control
#     (the anatomy and transit are altered) nor an active tumour, so pooling it
#     either way biases the comparison. It is used only where a step asks for it
#     explicitly — the Healthy-vs-PostSurgery sensitivity analysis.
#   * PostSurgery is deliberately OUTSIDE the ordinal severity scale: its
#     STAGE_GROUP_NUMERIC / CRC_STAGE_NUMERIC entries are NaN, so it cannot
#     silently enter a Spearman/trend test even if a caller forgets to mask it.
#
# For comparison, the earlier pipelines were wrong in different ways: the
# publication pipeline put HS in Early_CRC and MP in Advanced_CRC; the first v2
# run dropped HS and called MP "Adenoma".
# ----------------------------------------------------------------------

EXCLUDED_STUDY_GROUPS = {
    "YACHIDA-CRC-2019": ["MP"],
    # Kim's adenoma arm is dropped for the same reason MP is: adenomas are out of
    # scope for v2.1, so the validation cohort is Control vs Carcinoma only.
    "Kim_adenomas_2020": ["Adenoma"],
}

# Groups that are LOADED and kept in the tables, but excluded from every primary
# analysis; a step must opt in explicitly (see primary_mask / sensitivity_mask).
SENSITIVITY_GROUPS = {"YACHIDA-CRC-2019": ["HS"]}

CRC_STAGE_ORDER = ["Healthy", "Stage_0", "Stage_I_II", "Stage_III_IV"]
CRC_STAGE_NUMERIC = {
    "Healthy":      0,
    "Stage_0":      1,   # carcinoma in situ (earliest actual cancer diagnosis)
    "Stage_I_II":   2,
    "Stage_III_IV": 3,
    # HS has no position on a severity axis — see the policy note above.
    "HS":           np.nan,
}

# Fine-grained order INCLUDING the sensitivity group, for descriptive plots and
# tables that legitimately want to show it alongside the primary groups.
CRC_STAGE_ORDER_WITH_SENSITIVITY = ["Healthy", "PostSurgery", "Stage_0", "Stage_I_II", "Stage_III_IV"]

# Progression comparisons for ordinal/trend analyses (diagnosed cancer only)
CRC_PROGRESSION_COMPARISONS = [
    ("Stage_0",    "Stage_I_II"),
    ("Stage_0",    "Stage_III_IV"),
    ("Stage_I_II", "Stage_III_IV"),
]

# -----------------------------------------------------------------------
# GROUPED STAGING SCHEME (confirmatory analysis — Yachida-only)
# -----------------------------------------------------------------------
#   Healthy      : n=127  colonoscopy-confirmed controls
#   Early_CRC    : n=96   Stage_0 (n=27) + Stage_I_II (n=69) — diagnosed,
#                          localized/early-invasive cancer
#   Advanced_CRC : n=54   Stage_III_IV only — regionally/distantly advanced
#   --> PRIMARY analysis set: n=277
#
#   PostSurgery  : n=30   HS — sensitivity analyses only, never pooled
#   (MP, n=40: excluded from the study entirely)
#   --> loaded total: n=307
# -----------------------------------------------------------------------

STAGE_GROUP_MAP = {
    "Healthy":      "Healthy",
    "Stage_0":      "Early_CRC",
    "Stage_I_II":   "Early_CRC",
    "Stage_III_IV": "Advanced_CRC",
    "HS":           "PostSurgery",
}
# Primary, ordered severity axis. PostSurgery is intentionally absent.
STAGE_GROUP_ORDER = ["Healthy", "Early_CRC", "Advanced_CRC"]
STAGE_GROUP_NUMERIC = {"Healthy": 0, "Early_CRC": 1, "Advanced_CRC": 2, "PostSurgery": np.nan}

# Groups that exist in the data but are not on the primary axis.
SENSITIVITY_STAGE_GROUPS = ["PostSurgery"]
STAGE_GROUP_ORDER_WITH_SENSITIVITY = STAGE_GROUP_ORDER + SENSITIVITY_STAGE_GROUPS

# Convenience: fine-grained Study.Group -> grouped numeric directly
CRC_STAGE_GROUP_NUMERIC = {sg: STAGE_GROUP_NUMERIC[grp] for sg, grp in STAGE_GROUP_MAP.items()}

# Comparisons under the grouped scheme (pass group_col="Stage.Group" to
# differential_abundance()).
COMPARISONS_GROUP = {
    "YACHIDA-CRC-2019": [
        ("Healthy",   "Early_CRC"),      # early dysregulation vs baseline
        ("Healthy",   "Advanced_CRC"),   # advanced dysregulation vs baseline
        ("Early_CRC", "Advanced_CRC"),   # progression signal (diagnosed cancer only)
    ]
}

# Sensitivity-only comparisons: run and reported separately, never merged into the
# primary FDR family (that is the point of keeping them apart).
COMPARISONS_SENSITIVITY = {
    "YACHIDA-CRC-2019": [
        ("Healthy", "PostSurgery"),      # does a resected gut look like a control?
    ]
}

# CRC-only groups: the one context where tumour location is a legitimate covariate
# (every sample in these groups has, or had, a located tumour).
CRC_ONLY_STAGE_GROUPS = ["Early_CRC", "Advanced_CRC"]
CRC_ONLY_STUDY_GROUPS = ["Stage_0", "Stage_I_II", "Stage_III_IV"]

# -----------------------------------------------------------------------
# OPT-IN SEVERITY-TRANSFER SCALE (NB05 stage 6 only)
# -----------------------------------------------------------------------
# A SECOND, PARALLEL severity axis, used only by the cross-cohort severity-transfer
# analysis. It is deliberately separate from STAGE_GROUP_MAP / STAGE_GROUP_ORDER:
# adding MP to those would break the group policy (checks.py A1 asserts MP has no
# Stage.Group mapping) and, now that primary_mask() excludes rather than enumerates,
# would let adenomas walk into every primary analysis.
#
# WHAT THIS CAN AND CANNOT TEST (audited 2026-10-07, see CONTRACTS):
#   Neither validation cohort has ANY staging column. Sinha's Study.Group is a bare
#   0/1; Kim's is Control/Adenoma/Carcinoma. So TNM stage is UNFALSIFIABLE outside
#   Yachida -- a model that emits "Stage_I_II" for a Kim sample cannot be checked.
#   The coarsest axis Kim can actually verify is 3-level ordinal severity, and that
#   is what SEVERITY_EVAL_* encodes. Sinha supports the binary contrast only.
#
# Yachida MP (n=40) is the adenoma group: multiple polypoid adenomas with low-grade
# dysplasia, and its raw `Stage` column is "-" for all 40 (no carcinoma), which is
# what makes it the counterpart of Kim's Adenoma arm (n=102).
SEVERITY_TRANSFER_GROUPS = {"YACHIDA-CRC-2019": ["MP"], "Kim_adenomas_2020": ["Adenoma"]}

# 4-level TRAIN scale -- Yachida only, the full gradient. HS/PostSurgery is absent on
# purpose: a resected gut has no position on a severity axis.
SEVERITY_TRAIN_ORDER = ["Healthy", "Adenoma", "Early_CRC", "Advanced_CRC"]
SEVERITY_TRAIN_MAP = {
    "Healthy":      0,
    "MP":           1,   # adenoma
    "Stage_0":      2,   # carcinoma in situ
    "Stage_I_II":   2,
    "Stage_III_IV": 3,
}

# 3-level EVAL scale -- the cross-cohort axis. The model's Early_CRC/Advanced_CRC
# output collapses to "Carcinoma" at scoring time, because that is the finest
# distinction Kim's labels can confirm.
SEVERITY_EVAL_ORDER = ["Healthy", "Adenoma", "Carcinoma"]
SEVERITY_EVAL_MAP = {
    # Yachida (fine-grained Study.Group)
    "Healthy":      0,
    "MP":           1,
    "Stage_0":      2,
    "Stage_I_II":   2,
    "Stage_III_IV": 2,
    # Kim / Sinha, post-harmonisation
    "Adenoma":      1,
    "CRC":          2,
    # Kim raw labels, in case a caller reads its metadata without harmonising
    "Control":      0,
    "Carcinoma":    2,
    # Never on a severity axis.
    "HS":           np.nan,
    "PostSurgery":  np.nan,
}

# Collapsing the 4-level train scale onto the 3-level eval scale.
SEVERITY_TRAIN_TO_EVAL = {0: 0, 1: 1, 2: 2, 3: 2}

# -----------------------------------------------------------------------
# ANALYSIS-SET MASKS
# -----------------------------------------------------------------------
# Every primary analysis selects its samples through primary_mask(). This is the
# single mechanism keeping the sensitivity group (PostSurgery/HS) out of the main
# results: a notebook cannot include it by forgetting to filter, only by calling
# sensitivity_mask() or with_sensitivity() on purpose.


def primary_mask(meta_df, group_col="Stage.Group"):
    """Boolean mask for the PRIMARY analysis set: everything except the
    sensitivity group(s) (PostSurgery / HS).

    For Yachida, whose Stage.Group is fully mapped, this is exactly
    Healthy + Early_CRC + Advanced_CRC (n=277 of 307).

    For cohorts with no staging (Sinha, Kim) every sample is primary.

    v2.1 FIX (2026-10-07): this was an INCLUDE-list,
    `meta_df[group_col].isin(STAGE_GROUP_ORDER)`, guarded by an all-NaN
    short-circuit. Sinha and Kim harmonize to Healthy/"CRC" and "CRC" is not a
    STAGE_GROUP_MAP key, so their Stage.Group is PARTIALLY NaN: controls mapped to
    "Healthy", every case to NaN. The short-circuit therefore never fired and the
    include-list silently dropped every case from both validation cohorts
    (Sinha 131 -> 89, Kim 138 -> 102). Downstream that was invisible:
    differential_abundance() returns an empty frame when a group has <3 samples,
    and NB02 reports "0 significant" for it. Excluding the sensitivity groups
    instead of enumerating the primary ones cannot fail that way.

    Because unmapped rows now count as primary, a staged cohort must have a fully
    mapped Stage.Group; checks.py B1 asserts that for Yachida.
    """
    if group_col not in meta_df.columns:
        return pd.Series(True, index=meta_df.index)
    return ~meta_df[group_col].isin(SENSITIVITY_STAGE_GROUPS)


def sensitivity_mask(meta_df, group_col="Stage.Group"):
    """Boolean mask selecting ONLY the sensitivity group(s) (PostSurgery / HS)."""
    if group_col not in meta_df.columns:
        return pd.Series(False, index=meta_df.index)
    return meta_df[group_col].isin(SENSITIVITY_STAGE_GROUPS)


def crc_only_mask(meta_df, group_col="Stage.Group"):
    """Boolean mask for samples with a diagnosed carcinoma (Early_CRC +
    Advanced_CRC). This is the ONLY set in which tumour location is a legitimate
    covariate — see build_confounder_matrix()."""
    if group_col not in meta_df.columns:
        return pd.Series(False, index=meta_df.index)
    return meta_df[group_col].isin(CRC_ONLY_STAGE_GROUPS)


def apply_mask(mask, *frames, sample_col="Sample", meta_df=None):
    """Filter a metadata mask across aligned feature tables by Sample ID.

    `mask` is positional/index-aligned to `meta_df` (default: the first frame).
    Returns a list of row-filtered copies, in the order given, each reindexed.
    Using Sample IDs rather than positions means a feature table whose row order
    differs from the metadata's cannot be silently mis-sliced.
    """
    base = meta_df if meta_df is not None else frames[0]
    keep = set(base.loc[mask, sample_col])
    return [f[f[sample_col].isin(keep)].reset_index(drop=True) for f in frames]


# -----------------------------------------------------------------------
# YACHIDA CLINICAL CONFOUNDERS
# -----------------------------------------------------------------------
# v2.1: the full set is adjusted for, not just Age/BMI/Gender. Smoking differs
# strongly by group in this cohort (median Brinkman index 100 in Healthy vs
# 196–540 in the disease groups), so leaving it out would let a smoking effect
# masquerade as a cancer effect.
YACHIDA_COMPLETE_CONFOUNDERS = ["Age", "BMI", "Gender", "Brinkman.Index", "Alcohol"]

# Candidate confounders per cohort -- what each cohort actually records, measured on the
# harmonised metadata (2026-10-20). Yachida records smoking as the Brinkman index; Kim
# only as a binary Smoker flag (99% populated) and age in three bands; Sinha records
# neither smoking nor alcohol. Applying Yachida's list everywhere left Kim adjusted for
# age and sex only. Callers still drop any candidate that is <50% populated or constant.
COHORT_CONFOUNDERS = {
    "YACHIDA-CRC-2019":  ["Age", "BMI", "Gender", "Brinkman.Index", "Alcohol"],
    "SINHA_CRC_2016":    ["Age", "BMI", "Gender"],
    "Kim_adenomas_2020": ["Age", "Gender", "Smoker"],
}
YACHIDA_ALL_CONFOUNDERS = YACHIDA_COMPLETE_CONFOUNDERS
# Tumour location is a disease-status proxy WHENEVER non-cancer samples are in
# scope: it is "-" for every Healthy/HS/MP sample, so "NoTumor" alone identifies
# the non-cancer groups. It is admissible ONLY within crc_only_mask() samples,
# which build_confounder_matrix() enforces rather than documents.
YACHIDA_DISEASE_PROXY_COLS = ["Tumor.Location.Category"]

# Global Study.Group -> numeric ordering (cross-dataset; kept consistent
# with the Yachida-specific CRC_STAGE_NUMERIC tiers for the groups they share).
#
# Only labels the remaining four cohorts actually use. The IBD/IBS/ESRD removal
# took "nonIBD" (iHMP), "H"/"D"/"C" (MARS), "UC"/"CD" (Franzosa, iHMP) and "ESRD"
# (Wang) with it. "Control" and "Carcinoma" stay because they are Kim's RAW
# labels: harmonize_metadata() remaps them to Healthy/CRC before this map is
# applied, so they are a safety net for anything reading Kim's metadata directly.
STUDY_GROUP_MAP = {
    "Healthy": 0, "HP": 0, "Control": 0,
    # Pre-malignant, between healthy and the earliest carcinoma. 0.5 rather than a
    # renumbering so the existing integer tiers keep their meaning.
    #
    # These exist so that Study.Group.Numeric is not silently NaN for an admitted
    # adenoma row. Nothing currently READS Study.Group.Numeric (NB05 stage 6 labels
    # via SEVERITY_TRAIN_MAP / SEVERITY_EVAL_MAP, not this map), so today they are
    # defensive only -- but an unmapped label here is exactly the silent-NaN trap
    # that voided Kim's adenomas once already.
    #
    # NOTE this is STUDY_GROUP_MAP, not STAGE_GROUP_MAP -- the latter deliberately
    # has no MP entry (checks.py A1), which is what keeps adenomas out of the
    # primary analysis set now that primary_mask() excludes rather than enumerates.
    "MP": 0.5, "Adenoma": 0.5,
    "Stage_0": 1,
    "Stage_I_II": 2,
    "Stage_III_IV": 3, "Carcinoma": 3, "CRC": 3,
    "Gastrectomy": 7,
    # Not on the severity axis (see STAGE ENCODING): kept NaN on purpose.
    "HS": np.nan, "PostSurgery": np.nan,
}

# ============================================================================
# QUALITY-CONTROL THRESHOLDS (single source — see CHANGES.md finding #3)
# ============================================================================
PREVALENCE_SPECIES = 0.20   # species: >=20% of samples non-zero
PREVALENCE_MTB      = 0.15   # metabolites: >=15% (sparser; stricter reduces FDR inflation)
FDR_THRESHOLD        = 0.05
MAX_SPECIES_REDUCED  = 500   # top-variance species retained for NB02+ correlation/ML analyses
MAX_MTB_REDUCED       = 150   # top-variance metabolites retained for NB02+ correlation/ML analyses

# ============================================================================
# DATA LOADING
# ============================================================================

def _read_wide_numeric_tsv(fpath):
    """Memory-lean reader for very wide numeric TSVs ("Sample" + numeric
    columns, e.g. Yachida species: 57,703 columns). pandas.read_csv peaks at
    ~1.1-1.3 GB on that file because of parser buffers and allocator
    fragmentation; parsing each line straight into one float64 matrix needs
    ~2x the data size (~300 MB). Values are parsed with correct rounding, so
    they can differ from read_csv's fast float parser in the last bit
    (relative difference <= 1e-12). Falls back to pandas for anything that
    does not look like a plain Sample + numeric table."""
    with open(fpath, "r", encoding="utf-8") as f:
        header = f.readline().rstrip("\r\n").split("\t")
        if header[0] != "Sample":
            return pd.read_csv(fpath, sep="\t")
        n_feat = len(header) - 1
        samples, rows = [], []
        for line in f:
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) != n_feat + 1:
                return pd.read_csv(fpath, sep="\t")
            samples.append(parts[0])
            try:
                rows.append(np.array(parts[1:], dtype=np.float64))
            except ValueError:  # "NA", "" etc.
                rows.append(pd.to_numeric(pd.Series(parts[1:]), errors="coerce").to_numpy(np.float64))
    X = np.vstack(rows)
    del rows
    out = pd.DataFrame(X, columns=header[1:], copy=False)
    # Infer the Sample dtype the way read_csv would (Yachida IDs are integers;
    # keeping them as strings would silently break every Sample-based join
    # with the metadata table, which read_csv parses as int64).
    sample_col = pd.Series(samples)
    num = pd.to_numeric(sample_col, errors="coerce")
    if num.notna().all():
        sample_col = num.astype("int64") if (num == num.round()).all() else num
    out.insert(0, "Sample", sample_col.values)
    return out


def load_mtb_map(dataset, data_dir=DATA_DIR):
    """Read one cohort's metabolite map without loading its (large) abundance tables.
    Kim's ships as "mtb_map.tsv" rather than "mtb.map.tsv"; both names are tried."""
    for p in (data_dir / f"{dataset} mtb.map.tsv", data_dir / f"{dataset} mtb_map.tsv"):
        if p.exists():
            return pd.read_csv(p, sep="\t")
    raise FileNotFoundError(f"no mtb.map file for {dataset} in {data_dir}")


def load_all_datasets(data_dir=DATA_DIR, datasets=DATASETS, exclude=None):
    """Load all datasets into a nested dict: data[dataset_name][file_type] -> DataFrame.

    `exclude` ({dataset: [Study.Group, ...]}, default EXCLUDED_STUDY_GROUPS)
    drops those groups' samples from metadata, mtb and species right after
    loading -- i.e. before any QC/prevalence step -- so no downstream
    notebook ever sees them. Pass exclude={} to load everything.
    """
    if exclude is None:
        exclude = EXCLUDED_STUDY_GROUPS
    data = {}
    for ds in datasets:
        data[ds] = {}
        for ft in FILE_TYPES:
            # Kim's metabolite map ships as "mtb_map.tsv", not "mtb.map.tsv";
            # try the documented name first, then the known variant, so a missing
            # map is a real absence rather than a filename mismatch.
            candidates = [data_dir / f"{ds} {ft}.tsv"]
            if ft == "mtb.map":
                candidates.append(data_dir / f"{ds} mtb_map.tsv")
            fpath = next((p for p in candidates if p.exists()), None)
            if fpath is not None:
                if ft == "species" and fpath.stat().st_size > 20e6:
                    data[ds][ft] = _read_wide_numeric_tsv(fpath)
                else:
                    data[ds][ft] = pd.read_csv(fpath, sep="\t")
            else:
                warnings.warn(f"File not found: {candidates[0]}")
                data[ds][ft] = None

        drop_groups = exclude.get(ds, [])
        meta = data[ds].get("metadata")
        if drop_groups and meta is not None and "Study.Group" in meta.columns:
            drop_samples = set(meta.loc[meta["Study.Group"].isin(drop_groups), "Sample"])
            dropped = {}
            for ft in ("metadata", "mtb", "species"):
                df = data[ds].get(ft)
                if df is not None and "Sample" in df.columns:
                    keep = ~df["Sample"].isin(drop_samples)
                    dropped[ft] = int((~keep).sum())
                    data[ds][ft] = df[keep].reset_index(drop=True)
            print(f"[load_all_datasets] {ds}: excluded Study.Group {drop_groups} "
                  f"({len(drop_samples)} samples); rows dropped per table: {dropped}")
            if len(set(dropped.values())) > 1:
                raise ValueError(f"{ds}: exclusion dropped different row counts per table "
                                 f"{dropped} -- Sample IDs do not match across tables")
    return data


# ============================================================================
# METADATA HARMONIZATION
# ============================================================================

def encode_tumor_location_category(loc_series):
    """Categorical tumor-location encoding that treats "no tumor" as a real
    category rather than a value to later impute.

    Categories: "NoTumor" (raw "-", i.e. the non-cancer groups — Healthy, HS and
    the excluded MP; every Stage_0 patient has a recorded tumor location),
    "Distal" (Rectum/Left colon), "Proximal" (Right colon), "MultiSite" (more
    than one site recorded), "Unknown" (anything else / NaN in the source).

    This coarse Distal/Proximal split is kept for backwards compatibility and for
    descriptive tables. Analyses use encode_tumor_location(), which keeps rectum
    separate from left colon (they differ in blood supply, microbiota and
    clinical behaviour) and returns NaN — not a category — for non-cancer
    samples, so a tumour-location covariate cannot encode "has no tumour".
    """
    def _encode(v):
        if pd.isna(v):
            return "Unknown"
        s = str(v).strip()
        if s in ("-", ""):
            return "NoTumor"
        v_lower = s.lower()
        n_sites = v_lower.count("colon") + v_lower.count("rectum")
        if n_sites > 1:
            return "MultiSite"
        if "rectum" in v_lower or "left" in v_lower:
            return "Distal"
        if "right" in v_lower:
            return "Proximal"
        return "Unknown"

    return loc_series.apply(_encode).astype("category")


# Analysis-grade tumour sites. Yachida counts among the 150 carcinoma samples:
# Rectum 69, Right_colon 44, Left_colon 32, MultiSite 5.
TUMOR_SITE_ORDER = ["Rectum", "Left_colon", "Right_colon", "MultiSite"]


def encode_tumor_location(loc_series):
    """Tumour site for WITHIN-CANCER analyses: Rectum / Left_colon / Right_colon /
    MultiSite, and NaN for samples with no tumour recorded.

    Returning NaN (rather than a "NoTumor" level) for non-cancer samples is
    deliberate: it means a model fed this column on a mixed Healthy/CRC set gets
    missing values instead of a perfect disease-status indicator, and
    build_confounder_matrix() refuses the request outright.
    """
    def _encode(v):
        if pd.isna(v):
            return np.nan
        s = str(v).strip()
        if s in ("-", ""):
            return np.nan
        v_lower = s.lower()
        n_sites = v_lower.count("colon") + v_lower.count("rectum")
        if n_sites > 1:
            return "MultiSite"
        if "rectum" in v_lower:
            return "Rectum"
        if "left" in v_lower:
            return "Left_colon"
        if "right" in v_lower:
            return "Right_colon"
        return np.nan

    return loc_series.apply(_encode)


def harmonize_metadata(meta_df, dataset_name, exclude=None):
    """Map dataset-specific column names to a common schema.

    Returned columns
    ----------------
    Core (all datasets):
        Dataset, Sample, Subject, Study.Group, Study.Group.Numeric,
        Age, Gender, BMI

    Yachida staging (YACHIDA-CRC-2019 only; NaN elsewhere):
        CRC.Stage.Numeric   - fine-grained ordinal (Healthy 0, Stage_0 1,
                              Stage_I_II 2, Stage_III_IV 3; HS = NaN — see
                              STAGE ENCODING above)
        Stage.Group          - "Healthy" / "Early_CRC" / "Advanced_CRC" /
                               "PostSurgery" (HS; sensitivity only)
        Stage.Group.Numeric   - grouped ordinal (0-2; PostSurgery = NaN)

    Yachida confounders (100% complete; NaN for other datasets):
        Brinkman.Index, Alcohol

    Smoking (Kim: binary from "Chem ID / Smoking history"; NaN elsewhere):
        Smoker

    Yachida tumour site:
        Tumor.Location.Raw       - original string
        Tumor.Location.Category  - NoTumor / Distal / Proximal / MultiSite / Unknown
                                   (descriptive only)
        Tumor.Location.Site       - Rectum / Left_colon / Right_colon / MultiSite,
                                   NaN for non-cancer samples. Usable as a
                                   covariate ONLY within crc_only_mask() samples;
                                   build_confounder_matrix() enforces this.

    Rows in EXCLUDED_STUDY_GROUPS (Yachida MP, Kim Adenoma) are dropped
    defensively. Rows in SENSITIVITY_GROUPS (Yachida HS) are KEPT and labelled
    PostSurgery; primary analyses exclude them via primary_mask().

    `exclude` ({dataset: [Study.Group, ...]}, default EXCLUDED_STUDY_GROUPS)
    overrides which groups are dropped. Pass exclude={} to keep everything --
    needed by the opt-in severity-transfer analysis, which trains on the adenoma
    groups (Yachida MP, Kim Adenoma). Without this parameter the module constant
    was hardcoded here, so load_all_datasets(exclude={}) was silently undone one
    step later and the adenoma rows vanished anyway.
    """
    if exclude is None:
        exclude = EXCLUDED_STUDY_GROUPS
    drop_groups = exclude.get(dataset_name, [])
    if drop_groups and "Study.Group" in meta_df.columns:
        n_drop = int(meta_df["Study.Group"].isin(drop_groups).sum())
        if n_drop:
            warnings.warn(f"harmonize_metadata: dropping {n_drop} {drop_groups} rows "
                          f"from {dataset_name} (should already be excluded at load)")
            meta_df = meta_df[~meta_df["Study.Group"].isin(drop_groups)]
    meta_df = meta_df.reset_index(drop=True)

    # Build from Sample first: assigning a scalar to an empty DataFrame
    # creates a 0-row column that later array assignments reindex to NaN.
    hm = pd.DataFrame({"Sample": meta_df["Sample"].values})
    hm.insert(0, "Dataset", dataset_name)
    hm["Subject"] = meta_df["Subject"].values

    sg = meta_df["Study.Group"].copy()
    if dataset_name == "SINHA_CRC_2016":
        sg = sg.map({0: "Healthy", 1: "CRC"})
    elif dataset_name == "Kim_adenomas_2020":
        # Control -> Healthy, Carcinoma -> CRC (matching Sinha's labels so the two
        # validation cohorts share one comparison definition).
        #
        # "Adenoma" MUST be carried through rather than left out of this map. It is
        # normally dropped at load (EXCLUDED_STUDY_GROUPS), but the severity-transfer
        # analysis admits it via exclude={} -- and an unmapped label becomes NaN, so
        # omitting it here would silently void all 102 adenoma labels in exactly the
        # way the primary_mask include-list voided every validation case.
        sg = sg.map({"Control": "Healthy", "Carcinoma": "CRC", "Adenoma": "Adenoma"})

    hm["Study.Group"]         = sg.values
    hm["Study.Group.Numeric"] = sg.map(STUDY_GROUP_MAP).values

    # .values everywhere: index-aligned assignment silently yields NaN if
    # meta_df's index is not a clean RangeIndex.
    hm["CRC.Stage.Numeric"]   = sg.map(CRC_STAGE_NUMERIC).values
    hm["Stage.Group"]         = sg.map(STAGE_GROUP_MAP).values
    hm["Stage.Group.Numeric"] = hm["Stage.Group"].map(STAGE_GROUP_NUMERIC).values

    _AGE_MIDPOINTS = {
        'Under 20': 17.5, '20-29': 24.5, '30-39': 34.5,
        '40-49': 44.5, '50-59': 54.5, '60-69': 64.5, '70+': 75.0,
    }
    if "consent_age" in meta_df.columns:
        hm["Age"] = pd.to_numeric(meta_df["consent_age"], errors="coerce").values
    elif "Age" in meta_df.columns:
        if dataset_name == "Kim_adenomas_2020":
            hm["Age"] = meta_df["Age"].map(_AGE_MIDPOINTS).values
        else:
            hm["Age"] = pd.to_numeric(meta_df["Age"], errors="coerce").values
    else:
        hm["Age"] = np.nan

    hm["Gender"] = meta_df["Gender"].values if "Gender" in meta_df.columns else np.nan

    if "BMI_at_baseline" in meta_df.columns:
        hm["BMI"] = pd.to_numeric(meta_df["BMI_at_baseline"], errors="coerce").values
    elif "BMI" in meta_df.columns:
        hm["BMI"] = pd.to_numeric(meta_df["BMI"], errors="coerce").values
    else:
        hm["BMI"] = np.nan

    if dataset_name == "YACHIDA-CRC-2019":
        hm["Brinkman.Index"] = pd.to_numeric(
            meta_df.get("Brinkman Index", pd.Series(np.nan, index=meta_df.index)),
            errors="coerce"
        ).values
        hm["Alcohol"] = pd.to_numeric(
            meta_df.get("Alcohol", pd.Series(np.nan, index=meta_df.index)),
            errors="coerce"
        ).values
        raw_loc = meta_df.get(
            "Tumor location", pd.Series("-", index=meta_df.index)
        ).fillna("-").astype(str)
        hm["Tumor.Location.Raw"]      = raw_loc.replace("-", np.nan).values
        hm["Tumor.Location.Category"] = encode_tumor_location_category(raw_loc).values
        hm["Tumor.Location.Site"]     = encode_tumor_location(raw_loc).values
    else:
        hm["Brinkman.Index"]          = np.nan
        hm["Alcohol"]                 = np.nan
        hm["Tumor.Location.Raw"]      = np.nan
        hm["Tumor.Location.Category"] = np.nan
        hm["Tumor.Location.Site"]     = np.nan

    # Smoking: Kim records it categorically; Yachida's Brinkman.Index already
    # carries smoking dose, so a redundant binary is not derived from it.
    if dataset_name == "Kim_adenomas_2020" and "Chem ID / Smoking history" in meta_df.columns:
        hm["Smoker"] = (meta_df["Chem ID / Smoking history"].astype(str).str.lower()
                        .map({"smoke": 1.0, "non-smoke": 0.0}).values)
    else:
        hm["Smoker"] = np.nan

    return hm


def build_confounder_matrix(meta_df, include_tumor_location=False):
    """Build a numeric-ready confounder matrix for ML/regression use:
    Age, BMI, Brinkman.Index, Alcohol, Gender_enc.

    Tumour location (`include_tumor_location=True`) is admissible ONLY when every
    row of `meta_df` is a diagnosed carcinoma — Early_CRC or Advanced_CRC. The raw
    column is "-" for all non-cancer samples, so on a mixed set a tumour-site
    dummy is an exact disease-status indicator and inflates every model that sees
    it (measured: NB05 stage-1 accuracy 0.63 -> 0.39 once it was removed).

    This function therefore does not merely document that rule, it enforces it:
    asking for tumour location on a set containing any non-cancer sample raises
    ValueError. Callers that want it must first subset with
    utils.crc_only_mask(). Columns are TumorSite_<Rectum|Left_colon|Right_colon|
    MultiSite>; a site with no samples in the subset is omitted, and
    drop_first is NOT applied so each site stays individually interpretable.

    Only Age/BMI/Brinkman.Index/Alcohol may genuinely be missing; callers
    should impute those INSIDE each CV fold, not on this matrix as a whole.
    """
    out = meta_df[["Sample", "Age", "BMI", "Brinkman.Index", "Alcohol"]].copy()
    out["Gender_enc"] = (meta_df["Gender"].astype(str).str.lower() == "male").astype(float)

    if not include_tumor_location:
        return out

    if "Stage.Group" in meta_df.columns:
        non_cancer = meta_df.loc[~crc_only_mask(meta_df), "Stage.Group"]
        if len(non_cancer):
            offenders = sorted(set(non_cancer.dropna().astype(str))) or ["<missing Stage.Group>"]
            raise ValueError(
                "include_tumor_location=True requires a CRC-only sample set "
                f"({CRC_ONLY_STAGE_GROUPS}), but {len(non_cancer)} of {len(meta_df)} rows are "
                f"{offenders}. Tumour location is '-' for every non-cancer sample, so as a "
                "covariate it would encode disease status itself. Subset with "
                "utils.crc_only_mask(meta) first."
            )
    else:
        raise ValueError(
            "include_tumor_location=True needs a 'Stage.Group' column to verify that every "
            "sample is a carcinoma; this metadata has none."
        )

    if "Tumor.Location.Raw" not in meta_df.columns:
        raise ValueError("include_tumor_location=True but 'Tumor.Location.Raw' is absent.")

    site = encode_tumor_location(meta_df["Tumor.Location.Raw"])
    n_missing = int(site.isna().sum())
    if n_missing:
        raise ValueError(
            f"{n_missing} CRC sample(s) have no usable tumour location; cannot build a "
            "tumour-site covariate without imputing disease-adjacent information."
        )
    dummies = pd.get_dummies(pd.Categorical(site, categories=TUMOR_SITE_ORDER),
                             prefix="TumorSite", dtype=float)
    dummies = dummies.loc[:, dummies.sum(axis=0) > 0]          # drop empty sites
    dummies.index = out.index
    return pd.concat([out, dummies], axis=1)


# ============================================================================
# METABOLITE IDENTITY & CROSS-COHORT MATCHING  (v2.2)
# ============================================================================
# The v2.2 pivot is metabolome-wide, so metabolites need stable cross-cohort
# identifiers rather than the 7 hand-curated KEGG ids the polyamine era used.
#
# Measured on the three analysis cohorts' mtb.map files:
#   Yachida   450 compounds (KEGG 432, HMDB 421); 445 carry at least one id
#   Sinha     530 compounds (HMDB 354, KEGG 189)   <- KEGG alone loses most of Sinha
#   Kim       462 compounds (HMDB 366, KEGG 265); only cohort with Superpathway/Subpathway
#
# Matching on KEGG *or* HMDB rather than KEGG alone is what recovers the larger
# shared space -- HMDB alone is no better than KEGG, but the union is:
#                       KEGG only      KEGG or HMDB
#   all three cohorts         86              90
#   Yachida n Kim            147             152
#   Yachida n Sinha          107             113

_ID_SEP = re.compile(r"[;,|]")
_NULLISH = {"", "NA", "NAN", "NONE", "-", "N/A"}


def _clean_ids(value):
    """Every usable identifier in one mtb.map cell, as a list ([] if none)."""
    if value is None:
        return []
    if isinstance(value, float) and np.isnan(value):
        return []
    out = []
    for part in _ID_SEP.split(str(value)):
        p = part.strip()
        if p and p.upper() not in _NULLISH:
            out.append(p)
    return out


def metabolite_identity(mtb_map_df, dataset_name=None):
    """Resolve each metabolite to stable identifiers for cross-cohort matching.

    Returns one row per compound in `mtb_map_df`:
        Compound       the column name as it appears in the cohort's mtb matrix
        KEGG, HMDB     first id of each type, or NaN
        Canonical      KEGG if present, else HMDB, else "name:<compound>"
        Ids            frozenset of ALL ids (KEGG and HMDB) for this compound
        HighConfidence the cohort's own annotation-quality flag, where present
        Superpathway,  Metabolon pathway classes -- Kim only, NaN elsewhere
        Subpathway

    `Ids` is the matching key, not `Canonical`: two cohorts may annotate the same
    compound with different id types, so membership is tested by set intersection.
    A compound with no id at all is cohort-local and gets a "name:" canonical, so it
    can still be analysed within its own cohort but can never spuriously "match".
    """
    if mtb_map_df is None or len(mtb_map_df) == 0:
        return pd.DataFrame(columns=["Compound", "KEGG", "HMDB", "Canonical", "Ids",
                                     "HighConfidence", "Superpathway", "Subpathway",
                                     "Dataset"])

    cols = set(mtb_map_df.columns)
    if "Compound" not in cols:
        raise ValueError(
            f"mtb.map for {dataset_name or '<unknown>'} has no 'Compound' column; "
            f"got {sorted(cols)}"
        )

    rows = []
    for _, r in mtb_map_df.iterrows():
        kegg = _clean_ids(r.get("KEGG")) if "KEGG" in cols else []
        hmdb = _clean_ids(r.get("HMDB")) if "HMDB" in cols else []
        ids = frozenset(kegg) | frozenset(hmdb)
        compound = str(r["Compound"])
        canonical = (kegg[0] if kegg else hmdb[0] if hmdb else f"name:{compound}")
        rows.append({
            "Compound":       compound,
            "KEGG":           kegg[0] if kegg else np.nan,
            "HMDB":           hmdb[0] if hmdb else np.nan,
            "Canonical":      canonical,
            "Ids":            ids,
            "HighConfidence": r.get("High.Confidence.Annotation", np.nan),
            "Superpathway":   r.get("Superpathway", np.nan),
            "Subpathway":     r.get("Subpathway", np.nan),
            "Dataset":        dataset_name,
        })
    return pd.DataFrame(rows)


# Replication eligibility, most to least informative. "discovery_only" means the
# compound is not MEASURED elsewhere -- which is categorically different from having
# been tested and not reproduced.
REPLICATION_TIERS = ["all_three", "yachida_kim", "yachida_sinha", "discovery_only"]

# Measured partition of Yachida's 445 id-carrying compounds (see the table above).
# all_three 90 + kim-only 62 + sinha-only 23 + neither 270 = 445.
EXPECTED_MTB_OVERLAP = {"kim": 152, "sinha": 113, "both": 90, "annotated": 445}


def replication_tier(discovery_identity, target_identities):
    """Tag each discovery compound with which cohorts could ever replicate it.

    `discovery_identity` is a metabolite_identity() frame (Yachida).
    `target_identities`  is {cohort_name: metabolite_identity() frame}.

    Returns `discovery_identity` plus:
        In_<cohort>            bool, per target cohort
        N_replication_cohorts  int
        Replication_cohorts    comma-joined names, or "" when none
        Replication_tier       one of REPLICATION_TIERS

    WHY THIS EXISTS. Only 90 of Yachida's 445 annotated compounds are measurable in
    all three cohorts. Without this tag, a metabolite absent from Sinha's platform is
    indistinguishable in the results from one that was tested there and failed --
    "not replicated" would silently absorb "not measurable". Every result frame that
    makes a replication claim must carry this column.
    """
    out = discovery_identity.copy()
    target_id_sets = {}
    for name, frame in target_identities.items():
        s = set()
        if frame is not None and len(frame):
            for ids in frame["Ids"]:
                s |= set(ids)
        target_id_sets[name] = s

    for name, id_set in target_id_sets.items():
        out[f"In_{name}"] = [bool(ids & id_set) for ids in out["Ids"]]

    hit_cols = [f"In_{n}" for n in target_id_sets]
    out["N_replication_cohorts"] = out[hit_cols].sum(axis=1).astype(int) if hit_cols else 0
    out["Replication_cohorts"] = [
        ",".join(n for n in target_id_sets if row[f"In_{n}"])
        for _, row in out.iterrows()
    ]

    # Tier labels are specific to this study's three cohorts; fall back to a generic
    # count-based label if the cohort names are not the expected two targets.
    kim_col, sinha_col = "In_Kim_adenomas_2020", "In_SINHA_CRC_2016"

    def _tier(row):
        if kim_col in out.columns and sinha_col in out.columns:
            k, s = bool(row[kim_col]), bool(row[sinha_col])
            if k and s:
                return "all_three"
            if k:
                return "yachida_kim"
            if s:
                return "yachida_sinha"
            return "discovery_only"
        return "discovery_only" if row["N_replication_cohorts"] == 0 else "replicable"

    out["Replication_tier"] = [_tier(r) for _, r in out.iterrows()]
    return out


# ============================================================================
# POLYAMINE IDENTIFICATION  (v2.2: superseded -- see W1 of the pivot plan)
# ============================================================================
# RETAINED ONLY so the un-regenerated notebooks still import. The metabolome-wide
# redesign drops the privileged-compound concept entirely; nothing new should call
# these. Removal is gated on the notebooks being regenerated, because the on-disk
# .ipynb files and the edits_nb*.py generators have diverged.

def find_polyamine_columns(mtb_map_df, mtb_df, dataset_name):
    """Find polyamine metabolite columns in a dataset's mtb matrix.

    Returns a dict: {KEGG_ID: column_name_in_mtb} for found polyamines.
    """
    found = {}
    mtb_cols = set(mtb_df.columns) - {"Sample"}

    if mtb_map_df is not None and "KEGG" in mtb_map_df.columns and "Compound" in mtb_map_df.columns:
        for kegg_id in POLYAMINE_KEGG:
            matches = mtb_map_df[mtb_map_df["KEGG"] == kegg_id]
            for _, row in matches.iterrows():
                compound = row["Compound"]
                if compound in mtb_cols:
                    found[kegg_id] = compound

    for kegg_id in POLYAMINE_KEGG:
        if kegg_id in found:
            continue
        for col in mtb_cols:
            if col.startswith(kegg_id + "_"):
                found[kegg_id] = col
                break

    return found


def find_polyamine_columns_with_proxies(mtb_map_df, mtb_df, dataset_name):
    """find_polyamine_columns(), plus pre-specified acetylated proxies for the
    polyamines a cohort does not measure directly (POLYAMINE_PROXY_KEGG).

    Needed because Kim measures N-acetylputrescine but not putrescine, so a
    strict direct-match lookup would silently make that cohort unable to
    replicate the study's headline compound. The proxy is NOT treated as the
    parent compound: it is returned separately, so every table can mark which
    cohort/compound pairs rest on a derivative.

    Returns (direct, proxies) where
      direct  = {kegg_id: column}   the compound itself
      proxies = {kegg_id: (column, proxy_kegg_id)}  for parents measured only
                as a derivative, and only when the parent is absent.
    """
    direct = find_polyamine_columns(mtb_map_df, mtb_df, dataset_name)

    # Map every KEGG id present in this cohort to its column, so proxies can be
    # looked up the same way the direct compounds were.
    available = {}
    mtb_cols = set(mtb_df.columns) - {"Sample"}
    if mtb_map_df is not None and {"KEGG", "Compound"} <= set(mtb_map_df.columns):
        for _, row in mtb_map_df.dropna(subset=["KEGG"]).iterrows():
            if row["Compound"] in mtb_cols:
                available.setdefault(row["KEGG"], row["Compound"])
    for col in mtb_cols:
        head = str(col).split("_", 1)[0]
        if re.fullmatch(r"C\d{5}", head):
            available.setdefault(head, col)

    proxies = {}
    for parent, proxy_ids in POLYAMINE_PROXY_KEGG.items():
        if parent in direct:
            continue
        for pid in proxy_ids:
            if pid in available:
                proxies[parent] = (available[pid], pid)
                break

    # Name fallback, for derivatives the platform left unannotated.
    def _norm(s):
        return re.sub(r"[^a-z0-9]+", " ", str(s).lower()).strip()

    by_name = {_norm(c): c for c in mtb_cols}
    for parent, names in POLYAMINE_PROXY_NAMES.items():
        if parent in direct or parent in proxies:
            continue
        for nm in names:
            if nm in by_name:
                proxies[parent] = (by_name[nm], f"name:{nm}")
                break
    return direct, proxies


# ============================================================================
# TAXONOMY PARSING
# ============================================================================

def parse_taxonomy(full_string):
    """Parse 'd__X;p__Y;...;s__Z' into a dict of ranks."""
    rank_map = {
        "d": "domain", "p": "phylum", "c": "class",
        "o": "order",  "f": "family", "g": "genus",  "s": "species",
    }
    ranks = {}
    for part in full_string.split(";"):
        if "__" in part:
            prefix, name = part.split("__", 1)
            ranks[rank_map.get(prefix, prefix)] = name
    return ranks


def extract_species_name(full_string):
    """Species name; falls back to genus when the species rank is absent or
    empty ("s__"), and returns None when both are empty (unclassifiable)."""
    parsed = parse_taxonomy(full_string)
    if not parsed:
        return full_string
    sp = (parsed.get("species") or "").strip()
    if sp:
        return sp
    gn = (parsed.get("genus") or "").strip()
    return gn or None


def extract_genus(full_string):
    parsed = parse_taxonomy(full_string)
    return parsed.get("genus", "Unknown")


def reduce_species_names(species_df):
    """Reduce full taxonomy column names to species-level, summing duplicates.

    Works on a float64 array of the feature columns only. (The old version
    called species_df.values on a frame that includes the string "Sample"
    column, which materialised an object array of boxed floats -- ~640 MB
    for Yachida's 57,702 columns and the main cause of NB01 running out of
    memory on the 1.5 GB server.) Columns whose species and genus are both
    empty are dropped.
    """
    feat = species_df.drop(columns="Sample")
    X = feat.to_numpy(dtype=np.float64)
    short_names = [extract_species_name(c) for c in feat.columns]

    col_map: dict = {}
    for col, short in zip(feat.columns, short_names):
        if short is not None:
            col_map.setdefault(short, []).append(col)
    n_dropped = sum(s is None for s in short_names)
    if n_dropped:
        warnings.warn(f"reduce_species_names: dropped {n_dropped} unclassifiable columns")

    if all(len(v) == 1 for v in col_map.values()) and not n_dropped:
        out = pd.DataFrame(X, columns=short_names)
    else:
        pos = {c: i for i, c in enumerate(feat.columns)}
        names = list(col_map.keys())
        out_arr = np.empty((X.shape[0], len(names)), dtype=np.float64)
        for j, short in enumerate(names):
            idx = [pos[c] for c in col_map[short]]
            out_arr[:, j] = X[:, idx[0]] if len(idx) == 1 else X[:, idx].sum(axis=1)
        out = pd.DataFrame(out_arr, columns=names)
    out.insert(0, "Sample", species_df["Sample"].values)
    return out, col_map


# Labels that are not genera, compared case-insensitively after stripping. Taxa
# unresolved below family parse to "Unknown"; the cohort tables also ship an
# explicit "Unclassified" catch-all column.
_NON_GENUS_LABELS = {"", "unknown", "none", "nan", "unclassified", "unassigned"}


def normalise_genus_name(g):
    """Strip the rank prefix ("g__Blautia" -> "Blautia"), GTDB placeholder suffixes
    ("Blautia_A" -> "Blautia") and brackets, so the same genus matches across a
    shotgun and a 16S cohort.

    The rank prefix is stripped defensively. Every cohort currently on disk ships
    full semicolon-delimited taxonomy strings, which reach this function already
    parsed by extract_genus(), so the prefix branch is not exercised in production
    -- but a bare "g__Name" column must not be allowed to normalise to itself and
    then silently fail to match a basis built from parsed names.
    """
    g = re.sub(r"^[a-z]__", "", str(g).strip())
    g = re.sub(r"_[A-Z]+$", "", g)
    return g.strip("[]").strip()


def _genus_abundance(species_df, samples=None):
    """Genus-level abundance table: "Sample" + normalised genus columns, duplicate
    genera merged, unclassified columns dropped. No prevalence filter, no CLR.

    Shared by genus_clr_table() and shared_genera() so the two cannot disagree
    about what counts as a genus.
    """
    src = species_df
    if samples is not None:
        src = src[src["Sample"].isin(set(samples))].reset_index(drop=True)

    gen = aggregate_to_genus(src)
    feat = gen.drop(columns="Sample")
    feat.columns = [normalise_genus_name(c) for c in feat.columns]
    feat = feat.T.groupby(level=0).sum().T
    # Non-genus buckets. Sinha and Kim each carry a literal "Unclassified" column,
    # which is not a genus; the previous drop-list missed it, so it survived
    # aggregation and was treated as one. It never reached the cross-cohort basis
    # (Yachida has no such column to intersect with), but it polluted any
    # single-cohort genus analysis of those two.
    feat = feat.loc[:, [c for c in feat.columns
                        if str(c).strip().lower() not in _NON_GENUS_LABELS]]
    feat.insert(0, "Sample", gen["Sample"].values)
    return feat


def genus_clr_table(species_df, prevalence=None, samples=None, genera=None):
    """Genus-level CLR table, built the same way for every cohort.

    Shared by the cross-cohort steps (NB02's replication, NB05's transfer model) so
    the two cannot drift apart. Species-level names are collapsed to genus, GTDB
    suffixes normalised, genera merged by the normalised name, unclassified columns
    dropped, prevalence-filtered, then CLR-transformed.

    `samples` restricts to a sample list (e.g. the primary set) BEFORE the
    prevalence filter, so which genera are retained is not influenced by samples
    that the analysis excludes.

    `genera` FIXES THE FEATURE BASIS BEFORE THE CLR. Pass it for any cross-cohort
    model. CLR divides by the geometric mean of the columns present, so the
    reference depends on which columns those are: CLR-ing each cohort on its own
    genus set and intersecting the result afterwards puts the two cohorts on
    DIFFERENT scales (Yachida's reference spans its own prevalent genera, Kim's
    ~409), which biases every transfer systematically. That is what NB05 stage 3
    used to do. Reindexing to one shared list first makes the reference identical.

    When `genera` is given:
      * the prevalence filter is SKIPPED -- the list already defines the feature
        set, and filtering per cohort on top would hand each cohort different
        columns, which is the thing this argument exists to prevent;
      * genera absent from this cohort are ADDED as zero columns rather than
        dropped, so every cohort returns the same columns in the same order
        (scikit-learn indexes features positionally; an intersection does not
        guarantee a stable order).
    Choose the list with shared_genera(), which keeps the choice leak-free.
    """
    feat = _genus_abundance(species_df, samples=samples)

    if genera is not None:
        genera = list(genera)
        dupes = {g for g in genera if genera.count(g) > 1}
        if dupes:
            raise ValueError(f"genera contains duplicate names: {sorted(dupes)}")
        sample_col = feat["Sample"].values
        # fill_value=0.0 covers genera this cohort never measured: a true zero
        # abundance, which clr_transform's pseudocount then handles.
        sub = feat.drop(columns="Sample").reindex(columns=genera, fill_value=0.0)

        # Reindexing to the shared basis leaves each row summing to an arbitrary
        # fraction of its original mass, and that fraction differs per cohort --
        # Yachida goes from 13,645 genera to ~59, Sinha from 63. clr_transform's
        # multiplicative replacement closes each row before the log, so the two
        # cohorts meet the zero-replacement delta on the same footing and the CLR
        # reference really is identical across them (checks.py A8). Fixing the basis
        # alone was not sufficient; the scale leaked through the old additive
        # pseudocount.
        # Reported here rather than left to clr_transform's generic warning, because
        # at this point the cause is specific and nameable: the basis does not fit
        # this cohort.
        empty = int((sub.to_numpy(dtype=np.float64).sum(axis=1) <= 0).sum())
        if empty:
            warnings.warn(
                f"genus_clr_table: {empty} sample(s) contain none of the "
                f"{len(genera)} basis genera; their CLR row carries no information. "
                "Check the basis against this cohort."
            )

        sub.insert(0, "Sample", sample_col)
        return clr_transform(sub)

    if prevalence is None:
        prevalence = PREVALENCE_SPECIES
    feat = filter_by_prevalence(feat, threshold=prevalence)
    return clr_transform(feat)


def shared_genera(train_species_df, target_species_dfs, prevalence=None,
                  train_samples=None, target_samples=None):
    """Genera that are PREVALENT in the training cohort and merely PRESENT in every
    target cohort. Returns a sorted list, suitable for genus_clr_table(genera=...).

    The asymmetry is deliberate and is what keeps feature selection honest:
    applying a prevalence threshold to a target cohort would let the held-out data
    decide which features the model gets, which is leakage. Presence is a property
    of the measurement platform (did this cohort's assay cover the genus at all),
    not of its outcome, so intersecting on presence is safe.

    `target_samples`, if given, must be a list parallel to `target_species_dfs`.
    """
    if prevalence is None:
        prevalence = PREVALENCE_SPECIES

    train = _genus_abundance(train_species_df, samples=train_samples)
    train = filter_by_prevalence(train, threshold=prevalence)
    shared = {c for c in train.columns if c != "Sample"}

    targets = list(target_species_dfs)
    if target_samples is None:
        target_samples = [None] * len(targets)
    elif len(target_samples) != len(targets):
        raise ValueError(
            f"target_samples has {len(target_samples)} entries for "
            f"{len(targets)} target cohorts"
        )

    for tgt, smp in zip(targets, target_samples):
        present = _genus_abundance(tgt, samples=smp)
        shared &= {c for c in present.columns if c != "Sample"}

    return sorted(shared)


def _zscore_columns(M):
    """Column-wise z-score; a constant column becomes all zeros."""
    mu = M.mean(axis=0, keepdims=True)
    sd = M.std(axis=0, ddof=1, keepdims=True)
    return np.divide(M - mu, sd, out=np.zeros_like(M), where=sd > 0)


def transfer_genus_features(train_species, target_species, genera, train_samples=None,
                            target_samples=None, harmonise_detection=True,
                            standardise=True):
    """Genus features for a train -> target transfer, harmonised for detection depth.

    Pass UNFILTERED species tables for both cohorts (choose `genera` with
    shared_genera on the same unfiltered tables). Prevalence-filtering the target
    first -- as NB05 stage 3 originally did with NB01's filtered tables -- lets the
    held-out cohort choose the model's features (Sinha's shared basis fell from 55 to
    37 genera, losing Fusobacterium among others).

    1. Genus abundances on the shared basis, closed to proportions within each cohort.
    2. Detection-limit harmonisation: training proportions below the target's
       detection limit (its smallest positive proportion on the basis) are set to 0
       and rows re-closed. Shotgun Yachida resolves proportions down to ~2e-7 while 16S
       Sinha and Kim stop near 1e-4 and 5e-5, so without this the training cohort has
       ~2.5% zeros on the basis against the targets' ~38-48% -- and a transferred score
       can track detection depth instead of disease.
    3. CLR with ONE multiplicative-replacement delta shared by both cohorts
       (CLR_DELTA_FRACTION x the smallest positive proportion across both).
    4. Features z-scored WITHIN each cohort. Label-blind: it uses each cohort's own
       feature distribution, never an outcome, and removes per-genus location and scale
       offsets between platforms.

    Returns (train_df, target_df, info): DataFrames with "Sample" + `genera` columns,
    and an info dict with zero fractions before/after harmonisation, the detection
    limit, the shared delta, and presence/absence matrices for a sensitivity model.
    """
    genera = list(genera)
    if len(set(genera)) != len(genera):
        raise ValueError("transfer_genus_features: duplicate names in `genera`")

    def _table(species, samples):
        g = _genus_abundance(species, samples=samples)
        sample_ids = g["Sample"].to_numpy()
        M = g.drop(columns="Sample").reindex(columns=genera, fill_value=0.0).to_numpy(float)
        return sample_ids, M

    def _close(M):
        s = M.sum(axis=1, keepdims=True)
        return np.divide(M, s, out=np.zeros_like(M), where=s > 0)

    tr_ids, A = _table(train_species, train_samples)
    te_ids, B = _table(target_species, target_samples)
    A, B = _close(A), _close(B)
    if not (B > 0).any():
        raise ValueError("transfer_genus_features: the target has no positive value "
                         "on the shared basis")
    info = {"n_genera": len(genera),
            "zero_frac_train_before": float((A == 0).mean()),
            "zero_frac_target": float((B == 0).mean())}
    if harmonise_detection:
        lod = float(B[B > 0].min())
        A = _close(np.where(A < lod, 0.0, A))
        info["target_detection_limit"] = lod
    info["zero_frac_train"] = float((A == 0).mean())
    for name, M in (("train", A), ("target", B)):
        empty = int((M.sum(axis=1) == 0).sum())
        if empty:
            raise ValueError(f"transfer_genus_features: {empty} {name} sample(s) have no "
                             "detected genus on the shared basis")

    delta = CLR_DELTA_FRACTION * float(min(A[A > 0].min(), B[B > 0].min()))
    info["delta"] = delta

    def _clr(M):
        L = np.log(_multiplicative_replacement(M, delta=delta))
        return L - L.mean(axis=1, keepdims=True)

    LA, LB = _clr(A), _clr(B)
    if standardise:
        LA, LB = _zscore_columns(LA), _zscore_columns(LB)
    info["presence_train"] = (A > 0).astype(float)
    info["presence_target"] = (B > 0).astype(float)

    tr_df = pd.DataFrame(LA, columns=genera)
    tr_df.insert(0, "Sample", tr_ids)
    te_df = pd.DataFrame(LB, columns=genera)
    te_df.insert(0, "Sample", te_ids)
    return tr_df, te_df, info


def aggregate_to_genus(species_df):
    genus_map = {}
    for col in species_df.columns:
        if col == "Sample":
            continue
        genus = (
            extract_genus(col) if ";" in col
            else col.split()[0] if " " in col
            else col
        )
        genus_map.setdefault(genus, []).append(col)

    sums = {genus: species_df[cols].sum(axis=1).values for genus, cols in genus_map.items()}
    agg = pd.DataFrame(sums)          # built in one step (column-by-column insert fragments)
    agg.insert(0, "Sample", species_df["Sample"].values)
    return agg


# ============================================================================
# QC FUNCTIONS
# ============================================================================

def _num(df):
    """Feature columns of `df` (all but "Sample") as a float64 array plus the
    column names. Fast path when every feature column is already numeric
    (avoids a per-column apply(pd.to_numeric) copy of very wide tables)."""
    feat = df.drop(columns="Sample", errors="ignore")
    if all(pd.api.types.is_numeric_dtype(t) for t in feat.dtypes):
        return feat.to_numpy(dtype=np.float64), list(feat.columns)
    return (feat.apply(pd.to_numeric, errors="coerce").to_numpy(dtype=np.float64),
            list(feat.columns))


def is_prescaled(df):
    """True if a feature table contains negative values, i.e. it was already
    centred/standardised upstream (e.g. SINHA metabolites are z-scores), so
    half-min imputation, log transforms and ">0 = detected" prevalence
    filtering are not meaningful for it."""
    X, _ = _num(df)
    return bool(np.nanmin(X) < 0)


def compute_sample_qc(df, data_type="metabolite"):
    X, feature_cols = _num(df)

    qc = pd.DataFrame({"Sample": df["Sample"].values})
    qc["total_signal"]   = np.nansum(X, axis=1)
    qc["n_detected"]     = (X > 0).sum(axis=1)
    qc["n_features"]     = len(feature_cols)
    qc["detection_rate"] = qc["n_detected"] / qc["n_features"]

    if data_type == "species":
        # Row by row: full-matrix temporaries cost 3x the table size.
        shannon = np.zeros(X.shape[0])
        for i in range(X.shape[0]):
            x = X[i][X[i] > 0]
            if x.size:
                p = x / x.sum()
                shannon[i] = -(p * np.log(p)).sum()
        qc["shannon_diversity"] = shannon

    return qc


def compute_feature_qc(df):
    X, feature_cols = _num(df)

    fqc = pd.DataFrame({"Feature": feature_cols})
    fqc["detection_rate"] = (X > 0).mean(axis=0)
    fqc["mean"]     = np.nanmean(X, axis=0)
    fqc["std"]      = np.nanstd(X, axis=0, ddof=1)
    fqc["cv"]       = np.where(fqc["mean"] > 0, fqc["std"] / fqc["mean"], np.nan)
    fqc["variance"] = np.nanvar(X, axis=0, ddof=1)
    return fqc


# ============================================================================
# FILTERING
# ============================================================================

def filter_by_prevalence(df, threshold, protected_cols=None):
    """Remove features detected in fewer than `threshold` fraction of samples.

    Callers should pass the module-level PREVALENCE_SPECIES / PREVALENCE_MTB
    constants explicitly (no silent default) so the value used is visible
    at the call site and can't drift out of sync with what a notebook's
    own preview/printed text claims.
    """
    X, feature_cols = _num(df)
    prevalence = (X > 0).mean(axis=0)
    keep = [c for c, p in zip(feature_cols, prevalence) if p >= threshold]

    if protected_cols:
        for pc in protected_cols:
            if pc in feature_cols and pc not in keep:
                keep.append(pc)

    return df[["Sample"] + keep]


def filter_near_zero_variance(df, threshold=1e-10):
    X, feature_cols = _num(df)
    variances = np.nanvar(X, axis=0, ddof=1)
    keep = [c for c, v in zip(feature_cols, variances) if v > threshold]
    return df[["Sample"] + keep]


# ============================================================================
# TRANSFORMATION & NORMALISATION
# ============================================================================

# Multiplicative zero replacement: delta is this fraction of the smallest observed
# positive value, i.e. of the effective detection limit. 0.65 is the conventional
# choice for multiplicative replacement (Martin-Fernandez et al.).
CLR_DELTA_FRACTION = 0.65


def _multiplicative_replacement(X, delta_fraction=None, delta=None):
    """Close rows to sum 1, then replace zeros multiplicatively.

    Returns an array with no zeros, rows summing to 1, and the ratios between
    originally-positive parts unchanged -- which is the property an additive
    constant destroys.

    `delta`, if given, is used as-is instead of delta_fraction * min(positive). Pass
    it when two tables must share one replacement value (cross-cohort transfer):
    a per-table delta gives a shotgun cohort a far smaller floor than a 16S one.
    """
    if delta_fraction is None:
        delta_fraction = CLR_DELTA_FRACTION
    X = np.asarray(X, dtype=np.float64)
    X = np.where(X < 0, 0.0, X)          # negatives are not compositional

    totals = X.sum(axis=1, keepdims=True)
    empty = (totals <= 0).ravel()
    closed = np.divide(X, totals, out=np.zeros_like(X), where=totals > 0)

    positive = closed[closed > 0]
    if positive.size == 0:
        raise ValueError("clr_transform: no positive values anywhere, cannot transform")
    if delta is None:
        delta = delta_fraction * positive.min()

    n_zero = (closed == 0).sum(axis=1, keepdims=True)
    occupied = n_zero * delta            # share handed to the replaced zeros
    bad = (occupied >= 1.0).ravel() & ~empty
    if bad.any():
        raise ValueError(
            f"clr_transform: zero replacement would consume the whole composition "
            f"for {int(bad.sum())} sample(s) -- too many zeros relative to the "
            "detection limit. Filter features by prevalence first."
        )

    out = np.where(closed > 0, closed * (1.0 - occupied), delta)
    if empty.any():
        # No composition at all. A uniform row makes its CLR exactly zero
        # everywhere, which is the honest "carries no information" encoding.
        out[empty, :] = 1.0 / out.shape[1]
        warnings.warn(
            f"clr_transform: {int(empty.sum())} sample(s) have no positive values; "
            "their CLR row is uniformly zero and carries no information."
        )
    return out


def clr_transform(df, pseudocount=None):
    """Centred log-ratio transform with compositional zero replacement.

    `pseudocount=None` (the default) uses **multiplicative replacement**: rows are
    closed to sum 1, zeros become `CLR_DELTA_FRACTION * min(positive)`, and the
    observed parts are rescaled so the row still sums to 1. This preserves the
    ratios between observed parts -- the thing CLR exists to represent -- and is
    scale-invariant, because both the closure and delta are defined relative to
    the data.

    Passing a float restores the previous behaviour, a fixed constant ADDED to
    every cell. It is kept only so earlier results can be reproduced, and is not
    recommended: it is not scale-invariant and it distorts the ratios between
    observed parts. On the real tables the distortion was severe rather than
    academic. Yachida's median positive genus abundance is 1.4e-6 against the old
    constant of 1e-4, so 93% of its positive values sat *below* the constant and
    their real variation was flattened into the floor -- its measured-entry CLR
    minimum (-2.687) was indistinguishable from its true-zero floor (-2.704).
    Sinha, median 3.6e-3, passed through almost untouched. Training on one and
    transferring to the other therefore crossed a systematic scale mismatch
    (checks.py A8).
    """
    feature_cols = [c for c in df.columns if c != "Sample"]
    X = df[feature_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0).to_numpy(float)

    X = X + float(pseudocount) if pseudocount is not None else _multiplicative_replacement(X)

    log_X = np.log(X)
    clr_values = log_X - log_X.mean(axis=1, keepdims=True)

    result = pd.DataFrame(clr_values, columns=feature_cols)
    result.insert(0, "Sample", df["Sample"].values)
    return result


def impute_half_min(df, protected_cols=None):
    """Impute NaN with half the minimum NONZERO value per column (falls back
    to 0.5 only if a column has no positive values at all).

    The previous version took `.min()` over the whole column, which
    includes literal zeros -- for any column with both zeros and NaNs
    that silently collapsed every imputed value to the same disconnected
    constant (0.5) regardless of the feature's real scale.
    """
    feature_cols = [c for c in df.columns if c != "Sample"]
    out = df.copy()
    for col in feature_cols:
        vals = pd.to_numeric(out[col], errors="coerce")
        positive = vals[vals > 0]
        half_min = positive.min() / 2 if len(positive) > 0 else 0.5
        out[col] = vals.fillna(half_min)
    return out


def log_transform(df, pseudocount=1.0, base=2):
    feature_cols = [c for c in df.columns if c != "Sample"]
    X = df[feature_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0).values
    if base == 2:
        log_X = np.log2(X + pseudocount)
    elif base == 10:
        log_X = np.log10(X + pseudocount)
    else:
        log_X = np.log(X + pseudocount)

    result = pd.DataFrame(log_X, columns=feature_cols)
    result.insert(0, "Sample", df["Sample"].values)
    return result


# Metabolites detected in at least this fraction of samples form the PQN reference.
MTB_REFERENCE_DETECTION = 0.95


def transform_metabolites(df, ref_detection=MTB_REFERENCE_DETECTION):
    """log10 -> reference-set PQN -> non-detects tied at a per-metabolite floor.

    For absolute (not pre-scaled) metabolomes. Returns (transformed_df, info) where
    info holds the per-sample offset, the per-metabolite detection rate and the
    reference set, so the normalisation can be audited downstream.

    1. Detected values (x > 0) become log10(x + 1); non-detects are set aside.
    2. Probabilistic quotient normalisation on the reference set: metabolites
       detected in >= `ref_detection` of samples. Each sample's offset is its median
       log deviation from the reference profile over those metabolites only.
    3. Non-detects are placed at one floor per metabolite: half the lowest
       normalised detected value (log10 2 below it). Every non-detect therefore ranks
       below every detected value and all non-detects are tied, so a rank test or a
       median split sees "undetected" as one level rather than a spread.

    Why: the previous transform centred each sample on the median of log10(x + 1)
    over ALL metabolites with non-detects as 0. With 38.8% non-detects in Yachida the
    offset tracked each sample's detection rate (rho = -0.92), not its dilution;
    22 of 307 samples got no centring at all. Advanced CRC samples had slightly more
    non-detects (39.4% vs 38.1%), a 0.136 smaller offset, and so every detected
    value rose by ~0.14: 38 of 40 significant Healthy-vs-Advanced metabolites came
    out "elevated" (median MeanDiff +0.108). Non-detects were also spread by the
    offset (239 distinct values among Spermine's undetected cells), which defeated
    NB03's class-balance guard. checks.py C6 guards the mechanism: non-detects
    outside the reference set must not move any offset.
    """
    cols = [c for c in df.columns if c != "Sample"]
    X = df[cols].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    detected = np.isfinite(X) & (X > 0)
    L = np.full(X.shape, np.nan)
    L[detected] = np.log10(X[detected] + 1.0)
    det_rate = detected.mean(axis=0)
    ref = det_rate >= ref_detection
    if ref.sum() < 10:
        raise ValueError(
            f"transform_metabolites: only {int(ref.sum())} metabolites are detected in "
            f">= {ref_detection:.0%} of samples -- too few for a stable PQN reference")
    ref_profile = np.nanmedian(L[:, ref], axis=0)
    offset = np.nanmedian(L[:, ref] - ref_profile, axis=1)
    if not np.isfinite(offset).all():
        raise ValueError("transform_metabolites: a sample has no detected reference "
                         "metabolite, so it cannot be normalised")
    Ln = L - offset[:, None]
    floor = np.nanmin(np.where(detected, Ln, np.nan), axis=0) - np.log10(2.0)
    if not np.isfinite(floor).all():
        raise ValueError("transform_metabolites: a metabolite is never detected; "
                         "filter by prevalence first")
    out = np.where(detected, Ln, floor[None, :])
    result = pd.DataFrame(out, columns=cols)
    result.insert(0, "Sample", df["Sample"].values)
    info = {"method": f"log10(x+1) -> PQN on metabolites detected in >= "
                      f"{ref_detection:.0%} of samples -> non-detects at floor",
            "offset": pd.Series(offset, index=df["Sample"].values, name="PQN_offset"),
            "detection_rate": pd.Series(det_rate, index=cols, name="Detection_rate"),
            "reference": [c for c, r in zip(cols, ref) if r]}
    return result, info


def partial_spearman_grid(X, Y, Z):
    """Partial Spearman correlation of every column of X with every column of Y,
    adjusting for the columns of Z -- exactly pingouin.partial_corr(method="spearman").

    pingouin ranks x, y and every covariate (average ties) and takes the partial
    Pearson correlation of the ranks. That equals the Pearson correlation of the
    residuals of rank(x) and rank(y) regressed on [1, rank(Z)], so the residuals are
    computed once per column and every pair comes from one matrix product. The
    p-value is the two-sided t-test with n - 2 - k degrees of freedom, as in
    pingouin. Replaces one pingouin call per pair (about 1.8 ms each, wrapped in an
    `except Exception: continue` that silently dropped any pair that raised).
    checks.py C7 asserts the equivalence against pingouin itself.

    This is NOT "Spearman correlation of residuals from regressing raw values", which
    is a different estimator. Covariates are ranked like everything else, as pingouin
    does; a nominal covariate with three or more levels should be one-hot encoded by
    the caller rather than passed as integer codes.

    X: (n, p), Y: (n, q), Z: (n, k) numeric, complete cases.
    Returns (r, p) arrays of shape (p, q); NaN where a residualised column is constant.
    """
    from scipy.stats import t as _student_t

    X, Y, Z = (np.asarray(a, dtype=float) for a in (X, Y, Z))
    if Z.ndim == 1:
        Z = Z[:, None]
    n, k = X.shape[0], Z.shape[1]
    for name, arr in (("X", X), ("Y", Y), ("Z", Z)):
        if arr.shape[0] != n:
            raise ValueError(f"partial_spearman_grid: {name} has {arr.shape[0]} rows, "
                             f"expected {n}")
        if not np.isfinite(arr).all():
            raise ValueError(f"partial_spearman_grid: {name} contains NaN/inf; "
                             "pass complete cases")
    dof = n - 2 - k
    if dof < 1:
        raise ValueError(f"partial_spearman_grid: n={n} with k={k} covariates "
                         "leaves no degrees of freedom")

    ranked_z = rankdata(Z, axis=0) if k else np.empty((n, 0))
    design = np.column_stack([np.ones(n), ranked_z])

    def _residuals(A):
        R = rankdata(A, axis=0)
        beta, *_ = np.linalg.lstsq(design, R, rcond=None)
        return R - design @ beta

    EX, EY = _residuals(X), _residuals(Y)
    nx, ny = np.linalg.norm(EX, axis=0), np.linalg.norm(EY, axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        r = np.clip((EX.T @ EY) / np.outer(nx, ny), -1.0, 1.0)
        t = r * np.sqrt(dof / (1.0 - r ** 2))
    p = 2.0 * _student_t.sf(np.abs(t), dof)
    degenerate = (nx[:, None] <= 1e-9 * np.sqrt(n)) | (ny[None, :] <= 1e-9 * np.sqrt(n))
    r[degenerate] = np.nan
    p[degenerate] = np.nan
    return r, p


def spearman_grid(X, Y):
    """Spearman correlation of every column of X with every column of Y, as one matrix
    product: partial_spearman_grid with no covariates (Pearson of centred ranks, t-test
    with n - 2 df), which is what scipy.stats.spearmanr computes pair by pair."""
    X = np.asarray(X, dtype=float)
    return partial_spearman_grid(X, Y, np.empty((X.shape[0], 0)))


def sequential_permutation_p(observed, null_batch, h=10, b_max=1000, batch_size=64):
    """Besag-Clifford sequential Monte Carlo p-value (Besag & Clifford 1991, Biometrika).

    `null_batch(start, stop)` must return the null statistics of permutation indices
    start..stop-1, deterministically (fixed per-index seeds). Draws are taken in index
    order and the procedure stops at the h-th EXCEEDANCE (null >= observed): p = h / L,
    with L the number drawn by then. Otherwise it runs to b_max and p = (x + 1) /
    (b_max + 1) with x < h exceedances. A non-finite null statistic (a failed fit)
    counts as an exceedance, so failures can only make p larger.

    Why: NB03 ran a fixed 100 permutations and then applied BH, so the smallest q a
    target could reach was about m / 101 for m targets -- with 14 targets, any
    rejection needed at least three of them at the p floor, and the predictive evidence
    source was structurally starved. Null-like targets now stop after about h / p
    draws, and the budget goes to targets whose p is small; b_max >= 20 m lets a single
    target among m pass BH at 0.05. Valid under the null (checks.py C11), and the
    result does not depend on batch_size.

    Returns (p, n_drawn, n_exceed, stopped_early, null_draws).
    """
    if not np.isfinite(observed):
        raise ValueError("sequential_permutation_p: the observed statistic is not finite")
    if h < 1 or b_max < 1 or batch_size < 1:
        raise ValueError("sequential_permutation_p: h, b_max and batch_size must be >= 1")
    draws = np.empty(0, dtype=float)
    while draws.size < b_max:
        start, stop = draws.size, min(draws.size + batch_size, b_max)
        vals = np.asarray(null_batch(start, stop), dtype=float).ravel()
        if vals.size != stop - start:
            raise ValueError(f"sequential_permutation_p: null_batch({start}, {stop}) returned "
                             f"{vals.size} values")
        draws = np.concatenate([draws, vals])
        exceed = np.flatnonzero(~np.isfinite(draws) | (draws >= observed))
        if exceed.size >= h:
            n_drawn = int(exceed[h - 1]) + 1
            return h / n_drawn, n_drawn, h, True, draws[:n_drawn]
    x = int(np.sum(~np.isfinite(draws) | (draws >= observed)))
    return (x + 1) / (b_max + 1), int(b_max), x, False, draws


def residualise_on_confounders(y, Z):
    """OLS residual of the target `y` on an intercept and the confounders `Z` (n x k).

    NB03 asks whether the species predict a metabolite GIVEN the confounders: the
    target becomes this residual (then a median split) and the confounders leave the
    feature matrix. With them IN the features, the label permutation tested "target
    independent of species AND confounders", which the confounders alone reject
    (putrescine: AUC 0.67 from confounders only), so a significant p said nothing
    about the species. Exact when the confounding is linear.

    Residualise the TARGET, not the species: species residuals keep confounder
    information in their zero patterns (a zero-abundance sample's residual is a linear
    function of Z), and a tree model recovers it -- in a calibration study on the real
    species a confounder-only target was "predicted" at AUC 0.73 from residualised
    species. checks.py C10 guards the residual's properties.

    Z must be complete (impute it from the confounders alone first); rank-deficient Z
    (a constant column) is handled by least squares.
    """
    y = np.asarray(y, dtype=float)
    Z = np.asarray(Z, dtype=float)
    if Z.ndim == 1:
        Z = Z[:, None]
    if Z.shape[0] != y.shape[0]:
        raise ValueError(f"residualise_on_confounders: y has {y.shape[0]} rows, Z has "
                         f"{Z.shape[0]}")
    if not (np.isfinite(y).all() and np.isfinite(Z).all()):
        raise ValueError("residualise_on_confounders: y or Z has NaN/inf -- impute Z "
                         "(from the confounders alone) and drop missing targets first")
    design = np.column_stack([np.ones(len(y)), Z])
    beta, *_ = np.linalg.lstsq(design, y, rcond=None)
    return y - design @ beta


def genus_metabolite_replication(discovery, targets, q_threshold=None):
    """Cross-cohort replication of genus x metabolite Spearman associations.

    `discovery` and each element of `targets` are dicts:
        {"name": cohort, "genus": genus CLR table ("Sample" + genera),
         "mtb": metabolite table ("Sample" + compounds),
         "identity": metabolite_identity() frame for that cohort}
    Genus and metabolite tables are aligned on "Sample" (inner join).

    Returns (pairs, summary). `pairs` has one row per discovery (Genus, Compound) pair
    significant at q < q_threshold in the discovery cohort, with the compound's
    Replication_tier and, per target T:
        In_<T>             the pair is MEASURABLE in T: the compound is identity-matched
                           (KEGG or HMDB, by set intersection) and the genus is present
        Compound_<T>       the matched compound column in T
        Rho_<T>, P_<T>     T's Spearman correlation for the pair
        Q_<T>              BH over T's measurable discovery hits (a replication family)
        Replicated_<T>     Q_<T> < threshold with the discovery sign
        Contradicted_<T>   Q_<T> < threshold with the opposite sign
    A pair not measurable in T has In_<T> False and NaN statistics, and is never
    counted as a failed replication -- the distinction replication tiering exists for.
    """
    if q_threshold is None:
        q_threshold = FDR_THRESHOLD

    def _aligned(c):
        g = [x for x in c["genus"].columns if x != "Sample"]
        m = [x for x in c["mtb"].columns if x != "Sample"]
        d = c["genus"].merge(c["mtb"], on="Sample")
        return d, g, m

    d, genera, compounds = _aligned(discovery)
    rho, p = spearman_grid(d[genera].to_numpy(float), d[compounds].to_numpy(float))
    flat = p.ravel()
    q = np.full(flat.shape, np.nan)
    ok = np.isfinite(flat)
    if ok.any():
        q[ok] = multipletests(flat[ok], method="fdr_bh")[1]
    q = q.reshape(p.shape)
    gi, ci = np.where(q < q_threshold)

    disc_ids = {r.Compound: r.Ids for r in discovery["identity"].itertuples()}
    tiers = replication_tier(discovery["identity"],
                             {t["name"]: t["identity"] for t in targets})
    tier_of = dict(zip(tiers["Compound"], tiers["Replication_tier"]))
    canon_of = dict(zip(discovery["identity"]["Compound"], discovery["identity"]["Canonical"]))

    pairs = pd.DataFrame({
        "Genus": [genera[i] for i in gi],
        "Compound": [compounds[j] for j in ci],
        "Rho_discovery": rho[gi, ci], "P_discovery": p[gi, ci], "Q_discovery": q[gi, ci],
    })
    pairs["Canonical"] = pairs["Compound"].map(canon_of)
    pairs["Replication_tier"] = pairs["Compound"].map(tier_of).fillna("discovery_only")

    summary = []
    for t in targets:
        name = t["name"]
        td, t_genera, t_compounds = _aligned(t)
        t_ids = {r.Compound: r.Ids for r in t["identity"].itertuples()}
        by_id = {}
        for comp in t_compounds:
            for i in t_ids.get(comp, frozenset()):
                by_id.setdefault(i, set()).add(comp)

        def _match(comp):
            hits = {}
            for i in disc_ids.get(comp, frozenset()):
                for tc in by_id.get(i, ()):
                    hits[tc] = hits.get(tc, 0) + 1
            if not hits:
                return None
            return sorted(hits, key=lambda tc: (-hits[tc], tc))[0]   # most shared ids

        match = {comp: _match(comp) for comp in pairs["Compound"].unique()}
        t_gset = set(t_genera)
        pairs[f"Compound_{name}"] = pairs["Compound"].map(match)
        measurable = pairs[f"Compound_{name}"].notna() & pairs["Genus"].isin(t_gset)
        pairs[f"In_{name}"] = measurable.to_numpy()

        for col in ("Rho", "P", "Q"):
            pairs[f"{col}_{name}"] = np.nan
        if measurable.any():
            mg = sorted(pairs.loc[measurable, "Genus"].unique())
            mc = sorted(pairs.loc[measurable, f"Compound_{name}"].unique())
            r_t, p_t = spearman_grid(td[mg].to_numpy(float), td[mc].to_numpy(float))
            gpos = {g: k for k, g in enumerate(mg)}
            cpos = {c: k for k, c in enumerate(mc)}
            idx = pairs.index[measurable]
            gk = pairs.loc[idx, "Genus"].map(gpos).to_numpy()
            ck = pairs.loc[idx, f"Compound_{name}"].map(cpos).to_numpy()
            pairs.loc[idx, f"Rho_{name}"] = r_t[gk, ck]
            pairs.loc[idx, f"P_{name}"] = p_t[gk, ck]
            pv = pairs.loc[idx, f"P_{name}"].to_numpy(float)
            fin = np.isfinite(pv)
            qv = np.full(pv.shape, np.nan)
            if fin.any():
                qv[fin] = multipletests(pv[fin], method="fdr_bh")[1]
            pairs.loc[idx, f"Q_{name}"] = qv
        sig_t = pairs[f"Q_{name}"] < q_threshold
        same = np.sign(pairs[f"Rho_{name}"]) == np.sign(pairs["Rho_discovery"])
        pairs[f"Replicated_{name}"] = (sig_t & same).to_numpy()
        pairs[f"Contradicted_{name}"] = (sig_t & ~same).to_numpy()

        no_cmp = pairs[f"Compound_{name}"].isna()
        summary.append({
            "Target": name, "Discovery_hits": len(pairs),
            "Measurable": int(measurable.sum()),
            "Not_measurable_compound": int(no_cmp.sum()),
            "Not_measurable_genus": int((~no_cmp & ~pairs["Genus"].isin(t_gset)).sum()),
            "Replicated": int(pairs[f"Replicated_{name}"].sum()),
            "Contradicted": int(pairs[f"Contradicted_{name}"].sum()),
            "Measurable_not_significant": int((measurable & ~sig_t).sum()),
        })
    return pairs, pd.DataFrame(summary)


def zscore_per_dataset(df):
    feature_cols = [c for c in df.columns if c != "Sample"]
    X = df[feature_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0).values
    means = X.mean(axis=0, keepdims=True)
    stds = X.std(axis=0, keepdims=True)
    stds[stds == 0] = 1.0
    X_z = (X - means) / stds

    result = pd.DataFrame(X_z, columns=feature_cols)
    result.insert(0, "Sample", df["Sample"].values)
    return result


def quantile_normalize(df):
    feature_cols = [c for c in df.columns if c != "Sample"]
    X = df[feature_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0).values

    sorted_means = np.sort(X, axis=0).mean(axis=1)
    ranks = np.apply_along_axis(lambda col: rankdata(col, method="min") - 1, 0, X)
    normalized = sorted_means[ranks.astype(int)]

    result = pd.DataFrame(normalized, columns=feature_cols)
    result.insert(0, "Sample", df["Sample"].values)
    return result


# ============================================================================
# STATISTICAL ANALYSIS
# ============================================================================

def differential_abundance(df, meta, group1, group2, group_col="Study.Group"):
    """Mann-Whitney U test for each feature between two groups with FDR correction.

    Returns DataFrame with Feature, MeanGroup1, MeanGroup2, MeanDiff,
    Statistic, PValue, CliffsD, QValue, Significant columns.

    MeanDiff = MeanGroup2 - MeanGroup1, computed directly in whatever
    (linear or already-log/CLR-transformed) space `df` is in. This
    replaces the old Log2FC = log2((mean2+eps)/(mean1+eps)), which took
    an extra log of a ratio of already-transformed (frequently negative)
    values -- mathematically invalid and NaN whenever a group mean was
    negative. For already-log-transformed input, a plain mean DIFFERENCE
    in that space is the correct fold-change-equivalent (difference of
    logs = log of a ratio); it just isn't itself base-2 unless the input
    already was.
    """
    merged = df.merge(meta[["Sample", group_col]], on="Sample")
    g1 = merged[merged[group_col] == group1]
    g2 = merged[merged[group_col] == group2]

    feature_cols = [c for c in df.columns if c != "Sample"]
    results = []

    for feat in feature_cols:
        vals1 = g1[feat].dropna()
        vals2 = g2[feat].dropna()

        if len(vals1) < 3 or len(vals2) < 3:
            continue

        try:
            stat, pval = mannwhitneyu(vals1, vals2, alternative="two-sided")
        except ValueError:
            continue

        n1, n2 = len(vals1), len(vals2)
        cliffs_d = (2 * stat / (n1 * n2)) - 1
        mean1 = vals1.mean()
        mean2 = vals2.mean()

        results.append({
            "Feature":    feat,
            "MeanGroup1": mean1,
            "MeanGroup2": mean2,
            "MeanDiff":   mean2 - mean1,
            "Statistic":  stat,
            "PValue":     pval,
            "CliffsD":    cliffs_d,
        })

    if not results:
        return pd.DataFrame()

    results_df = pd.DataFrame(results)
    reject, qvals, _, _ = multipletests(results_df["PValue"], method="fdr_bh")
    results_df["QValue"]      = qvals
    results_df["Significant"] = reject

    return results_df.sort_values("QValue")


def _coef_mode(coefs, bw="scott"):
    """Mode of a coefficient distribution, via a Gaussian KDE evaluated on a grid.

    LinDA's bias-correction term. Falls back to the median when the KDE cannot be
    fitted (fewer than 3 finite values, or zero variance).
    """
    from scipy.stats import gaussian_kde

    c = np.asarray(coefs, dtype=float)
    c = c[np.isfinite(c)]
    if c.size < 3 or np.std(c) == 0:
        return float(np.median(c)) if c.size else 0.0
    try:
        kde = gaussian_kde(c, bw_method=bw)
        grid = np.linspace(c.min(), c.max(), 512)
        return float(grid[np.argmax(kde(grid))])
    except Exception:
        return float(np.median(c))


def clr_lm_da(species_df, meta_df, group1, group2, group_col="Stage.Group",
              covars=None, already_clr=True, pseudocount=None,
              bias_correct=True, min_per_group=5):
    """Covariate-adjusted compositional differential abundance (LinDA-style).

    Complements differential_abundance(): that is an unadjusted rank test, this is
    a linear model on CLR values that can adjust for Age/BMI/Gender/smoking/
    alcohol. A feature is only called a robust hit in NB02 when BOTH agree in
    direction at q<0.05, which is a far stricter bar than either alone.

    Method
    ------
    For each feature: OLS of its CLR value on a binary group indicator plus
    `covars`. The group coefficient is then bias-corrected by subtracting the mode
    of all features' group coefficients (Zhou et al. 2022, LinDA): a CLR value is
    a ratio to the sample's geometric mean, so if a subset of taxa genuinely shifts
    the whole composition, EVERY feature's coefficient picks up a shared offset.
    Subtracting the mode removes that offset under the assumption that most
    features are not differentially abundant, which is what makes the result
    interpretable as absolute rather than relative change.

    Parameters
    ----------
    species_df : "Sample" + feature columns. CLR-transformed when
        already_clr=True (NB01's transformed_species); otherwise relative
        abundances, CLR-transformed here with `pseudocount`.
    covars : list of harmonized metadata columns. "Gender" is encoded to 0/1
        automatically; the rest must be numeric. Rows with any missing covariate
        are dropped (complete-case) and the count is reported in the result's
        attrs, never imputed — imputing a covariate to keep a sample in a
        hypothesis test fabricates precision.

    Returns
    -------
    DataFrame: Feature, Coef (group2 - group1, CLR units), Coef_raw, SE, TStat,
    PValue, QValue, Significant, Direction, N_group1, N_group2 — plus
    `.attrs["covariates"]`, `.attrs["n_dropped_missing_covar"]`,
    `.attrs["bias_offset"]`. Empty DataFrame if either arm is too small.
    """
    covars = list(covars or [])
    spc = species_df if already_clr else clr_transform(species_df, pseudocount=pseudocount)

    meta_cols = [group_col] + [c for c in covars if c in meta_df.columns]
    md = meta_df[["Sample"] + meta_cols].copy()
    md = md[md[group_col].isin([group1, group2])]

    if "Gender" in covars and "Gender" in md.columns:
        md["Gender"] = (md["Gender"].astype(str).str.lower() == "male").astype(float)

    use_covars = [c for c in covars if c in md.columns]
    missing_covars = [c for c in covars if c not in md.columns]
    if missing_covars:
        warnings.warn(f"clr_lm_da: covariates absent from metadata, ignored: {missing_covars}")

    n_before = len(md)
    if use_covars:
        md[use_covars] = md[use_covars].apply(pd.to_numeric, errors="coerce")
        md = md.dropna(subset=use_covars)
    n_dropped = n_before - len(md)

    # Drop covariates that are constant within this subset (they make X singular
    # and carry no information here, e.g. Brinkman.Index in a cohort that has none).
    const_covars = [c for c in use_covars if md[c].nunique(dropna=True) < 2]
    if const_covars:
        use_covars = [c for c in use_covars if c not in const_covars]

    merged = md.merge(spc, on="Sample")
    n1 = int((merged[group_col] == group1).sum())
    n2 = int((merged[group_col] == group2).sum())
    if min(n1, n2) < min_per_group:
        warnings.warn(f"clr_lm_da: {group1} vs {group2} has n={n1}/{n2} "
                      f"(< {min_per_group}) -- skipped")
        return pd.DataFrame()

    feature_cols = [c for c in spc.columns if c != "Sample"]
    g = (merged[group_col] == group2).astype(float).to_numpy()
    # Design: intercept + group + covariates, built once and reused per feature.
    X_cols = [np.ones(len(merged)), g] + [merged[c].to_numpy(dtype=float) for c in use_covars]
    X = np.column_stack(X_cols)
    dof = X.shape[0] - X.shape[1]
    if dof <= 0:
        warnings.warn("clr_lm_da: more parameters than samples -- skipped")
        return pd.DataFrame()

    # Solve all features in one least-squares call: Y is (n_samples, n_features),
    # so B is (n_params, n_features). One QR instead of ~500 separate OLS fits.
    Y = merged[feature_cols].to_numpy(dtype=float)
    B, _, _, _ = np.linalg.lstsq(X, Y, rcond=None)
    resid = Y - X @ B
    sigma2 = (resid ** 2).sum(axis=0) / dof
    XtX_inv = np.linalg.pinv(X.T @ X)
    se_group = np.sqrt(np.maximum(XtX_inv[1, 1], 0) * sigma2)

    coef_raw = B[1, :]
    offset = _coef_mode(coef_raw) if bias_correct else 0.0
    coef = coef_raw - offset

    with np.errstate(divide="ignore", invalid="ignore"):
        tstat = np.where(se_group > 0, coef / se_group, np.nan)
    from scipy.stats import t as _t
    pval = 2 * _t.sf(np.abs(tstat), dof)

    out = pd.DataFrame({
        "Feature":  feature_cols,
        "Coef":     coef,
        "Coef_raw": coef_raw,
        "SE":       se_group,
        "TStat":    tstat,
        "PValue":   pval,
        "N_group1": n1,
        "N_group2": n2,
    })
    ok = out["PValue"].notna()
    out["QValue"] = np.nan
    if ok.any():
        out.loc[ok, "QValue"] = multipletests(out.loc[ok, "PValue"], method="fdr_bh")[1]
    out["Significant"] = out["QValue"] < FDR_THRESHOLD
    out["Direction"] = np.where(out["Coef"] > 0, f"Up in {group2}", f"Up in {group1}")

    out.attrs["covariates"] = use_covars
    out.attrs["dropped_constant_covariates"] = const_covars
    out.attrs["n_dropped_missing_covar"] = int(n_dropped)
    out.attrs["bias_offset"] = float(offset)
    return out.sort_values("QValue", na_position="last").reset_index(drop=True)


def progression_correlation_analysis(
    species_df, mtb_df, meta_df,
    stage_col="CRC.Stage.Numeric",
    min_samples=10,
):
    """Correlate every species x metabolite pair with ordinal CRC stage.

    Species-stage and metabolite-stage correlations are computed ONCE per
    unique feature and FDR-corrected on that deduplicated set (the true
    number of tests), then joined back onto the long pairwise table. The
    previous version ran multipletests() on the full long-format table,
    where each unique species-stage p-value was duplicated once per
    metabolite partner (and vice versa) -- inflating the effective test
    count by ~2-3 orders of magnitude and corrupting every reported
    q-value for both marginal tests.

    Returns a long-form DataFrame with QVal_species_stage / QVal_metabolite_stage
    (now correct) and QVal_spc_mtb (was already correct -- each species x
    metabolite pair is a genuinely unique test).
    """
    common = sorted(
        set(species_df["Sample"]) & set(mtb_df["Sample"]) & set(meta_df["Sample"])
    )

    spc = species_df.set_index("Sample").loc[common]
    mtb = mtb_df.set_index("Sample").loc[common]
    meta = meta_df.set_index("Sample").loc[common]

    stage = meta[stage_col]
    valid = stage.notna()
    spc, mtb, stage = spc.loc[valid], mtb.loc[valid], stage.loc[valid]

    if len(stage) < min_samples:
        warnings.warn(
            f"progression_correlation_analysis: only {len(stage)} samples "
            f"with valid stage -- skipping."
        )
        return pd.DataFrame()

    # --- unique per-feature stage correlations, computed once each ---
    spc_stage_rows = []
    for spc_name in spc.columns:
        vals = spc[spc_name].values
        if np.std(vals) == 0:
            continue
        rho, pval = spearmanr(vals, stage.values)
        if np.isnan(rho):
            continue
        spc_stage_rows.append({"Species": spc_name, "Rho_species_stage": rho, "PVal_species_stage": pval})
    spc_stage_df = pd.DataFrame(spc_stage_rows)
    if not spc_stage_df.empty:
        _, q, _, _ = multipletests(spc_stage_df["PVal_species_stage"], method="fdr_bh")
        spc_stage_df["QVal_species_stage"] = q

    mtb_stage_rows = []
    for mtb_name in mtb.columns:
        vals = mtb[mtb_name].values
        if np.std(vals) == 0:
            continue
        rho, pval = spearmanr(vals, stage.values)
        if np.isnan(rho):
            continue
        mtb_stage_rows.append({"Metabolite": mtb_name, "Rho_metabolite_stage": rho, "PVal_metabolite_stage": pval})
    mtb_stage_df = pd.DataFrame(mtb_stage_rows)
    if not mtb_stage_df.empty:
        _, q, _, _ = multipletests(mtb_stage_df["PVal_metabolite_stage"], method="fdr_bh")
        mtb_stage_df["QVal_metabolite_stage"] = q

    if spc_stage_df.empty or mtb_stage_df.empty:
        return pd.DataFrame()

    # --- pairwise species x metabolite correlation (genuinely one test each) ---
    pair_rows = []
    spc_vals_map = {c: spc[c].values for c in spc_stage_df["Species"]}
    mtb_vals_map = {c: mtb[c].values for c in mtb_stage_df["Metabolite"]}
    for spc_name, spc_vals in spc_vals_map.items():
        for mtb_name, mtb_vals in mtb_vals_map.items():
            rho_sm, p_sm = spearmanr(spc_vals, mtb_vals)
            pair_rows.append({"Species": spc_name, "Metabolite": mtb_name,
                               "Rho_spc_mtb": rho_sm, "PVal_spc_mtb": p_sm})
    pair_df = pd.DataFrame(pair_rows)
    _, q, _, _ = multipletests(pair_df["PVal_spc_mtb"].fillna(1.0), method="fdr_bh")
    pair_df["QVal_spc_mtb"] = q

    out = pair_df.merge(spc_stage_df, on="Species").merge(mtb_stage_df, on="Metabolite")
    return out.sort_values("QVal_species_stage")


def compute_correlations(species_df, mtb_df, target_metabolites=None):
    merged = species_df.merge(mtb_df, on="Sample", suffixes=("_spc", "_mtb"))
    spc_cols = [c for c in species_df.columns if c != "Sample"]

    if target_metabolites:
        mtb_cols = [c for c in target_metabolites if c in mtb_df.columns]
    else:
        mtb_cols = [c for c in mtb_df.columns if c != "Sample"]

    results = []
    for spc in spc_cols:
        for mtb in mtb_cols:
            spc_key = spc + "_spc" if spc + "_spc" in merged.columns else spc
            mtb_key = mtb + "_mtb" if mtb + "_mtb" in merged.columns else mtb
            x = merged[spc_key].values
            y = merged[mtb_key].values
            if x.std() == 0 or y.std() == 0:
                continue
            try:
                rho, pval = spearmanr(x, y)
                if not np.isnan(rho):
                    results.append({"Species": spc, "Metabolite": mtb, "Rho": rho, "PValue": pval})
            except Exception:
                continue

    if not results:
        return pd.DataFrame()

    results_df = pd.DataFrame(results)
    _, qvals, _, _ = multipletests(results_df["PValue"], method="fdr_bh")
    results_df["QValue"] = qvals
    return results_df


def find_consistent_correlations(corr_dicts, kegg_map=None, min_datasets=2, q_threshold=0.05):
    """Find species-metabolite pairs significant in >= min_datasets.

    If `kegg_map` (dict: dataset -> {mtb_column_name: kegg_id}) is given,
    matching is done on KEGG ID rather than raw column-name string, so
    cohorts using different metabolomics platforms/naming conventions
    (e.g. "C00134_Putrescine(1,4-Butanediamine)" vs "PUTRESCINE" vs
    "HILp_QI163__putrescine") are correctly recognised as the same
    compound. Without kegg_map, falls back to exact string matching
    (only reliable for cohorts sharing a naming convention).
    """
    all_results = []
    for ds, corr_df in corr_dicts.items():
        if corr_df is None or corr_df.empty:
            continue
        sig = corr_df[corr_df["QValue"] < q_threshold].copy()
        sig["Dataset"] = ds
        if kegg_map is not None:
            ds_map = kegg_map.get(ds, {})
            sig["Metabolite_key"] = sig["Metabolite"].map(ds_map).fillna(sig["Metabolite"])
        else:
            sig["Metabolite_key"] = sig["Metabolite"]
        all_results.append(sig)

    if not all_results:
        return pd.DataFrame()

    combined = pd.concat(all_results, ignore_index=True)
    counts = combined.groupby(["Species", "Metabolite_key"]).agg(
        n_datasets=("Dataset", "nunique"),
        mean_rho=("Rho", "mean"),
        datasets=("Dataset", lambda x: list(x)),
        example_metabolite_names=("Metabolite", lambda x: list(set(x))),
    ).reset_index()

    return counts[counts["n_datasets"] >= min_datasets].sort_values("n_datasets", ascending=False)


# ============================================================================
# AGORA2 SPECIES-LEVEL ENZYME PRESENCE (genomic capability evidence)
# ----------------------------------------------------------------------------
# History: disabled 2026-10-07 (v2.1). Building the EC index scans 1,439 AGORA2
# reconstructions (1.8 GB), which overloaded the 2-core / 1.5 GB analysis server;
# every lookup returned None and NB08 R1/R3, NB09 Modules B/G and synthesis Source B
# were empty. Re-enabled 2026-10-21: that server is gone and the reconstructions are
# still on disk.
#
# The index is an INPUT, built ONCE outside the pipeline:
#     python tools/build_agora2_index.py
# which writes AGORA2_INDEX_CACHE plus a sidecar recording the source fingerprint
# (checks.py D5 verifies it). Notebooks only LOAD it (load_agora2_ec_index); none may
# trigger the scan.
#
# What it is evidence of: a species' draft reconstruction contains a reaction with
# that EC producing the compound. Reconstructions carry gap-filled reactions, so this
# is CAPABILITY, never production; a species without a model is None (no evidence),
# never False. POLYAMINE_EC lists producing reactions only (checks.py D1).
# ============================================================================

AGORA2_ENABLED = True

_EC_RE = re.compile(r"\d+\.\d+\.\d+\.(?:n?\d+|-)")


def _clean_ec_index(index: dict) -> dict:
    """Normalise every EC list in an index to clean, de-duplicated EC strings
    (repairs caches built before the _EC_RE fix in place)."""
    for sp, kmap in index.items():
        for kid, ecs in kmap.items():
            kmap[kid] = sorted(set(_EC_RE.findall(" ".join(map(str, ecs)))))
    return index


def build_agora2_ec_index(agora2_data_dir: Path, cache_path: Path | None = None,
                           force_rebuild: bool = False) -> dict:
    """Build {species_id: {kegg_compound_id: [ec_numbers]}} from AGORA2 species
    reconstructions. `species_id` matches the filename stem (e.g.
    "s__Sutterella_wadsworthensis_A") -- model["id"] itself is a generic
    "model" placeholder in every AGORA2 file, not species-specific.

    Scanning all ~1439 files takes a few minutes; the result is cached as
    JSON at `cache_path` and reused on subsequent calls unless
    force_rebuild=True. A compatible cache already exists at
    RESULTS_DIR/tables/agora2_species_ec_index.json (built during the
    CRC_Gut_Analysis_Git AGORA2 migration) and will be reused automatically
    if `cache_path` points at it.
    """
    import json as _json

    if cache_path is not None and Path(cache_path).exists() and not force_rebuild:
        with open(cache_path, "r", encoding="utf-8") as f:
            return _clean_ec_index(_json.load(f))

    index: dict = {}
    for fpath in sorted(Path(agora2_data_dir).glob("*.json")):
        try:
            with open(fpath, "r", encoding="utf-8") as f:
                model = _json.load(f)
        except Exception:
            continue

        species_id = fpath.stem
        met_to_kegg = {}
        for met in model.get("metabolites", []):
            kid = met.get("annotation", {}).get("kegg.compound")
            if kid and kid != "NaN":
                met_to_kegg[met["id"]] = kid
        if not met_to_kegg:
            continue

        kegg_ecs: dict = {}
        for rxn in model.get("reactions", []):
            ec_raw = rxn.get("annotation", {}).get("ec-code")
            if not ec_raw or ec_raw == "NaN":
                continue
            # Regex, not split(): raw values are often stringified lists like
            # "['2.3.3.1', '4.1.3.12']", which a ;/, split turns into
            # "['2.3.3.1'" / "'4.1.3.12']" that never match a clean EC.
            ecs = _EC_RE.findall(str(ec_raw))
            if not ecs:
                continue
            for met_id, coeff in rxn.get("metabolites", {}).items():
                if coeff <= 0:
                    continue  # only production (product) direction counts as "encodes"
                kid = met_to_kegg.get(met_id)
                if kid:
                    kegg_ecs.setdefault(kid, set()).update(ecs)

        if kegg_ecs:
            index[species_id] = {k: sorted(v) for k, v in kegg_ecs.items()}

    if cache_path is not None:
        cache_path = Path(cache_path)
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        with open(cache_path, "w", encoding="utf-8") as f:
            _json.dump(index, f)

    return index


_GTDB_SUFFIX_RE = re.compile(r"_[A-Z]{1,2}(?=_|$)")


def agora2_species_id(species: str, agora2_ec_index: dict):
    """Map a study species name to its AGORA2 model id: (id, match_type).

    match_type is "exact" when "s__" + the name with non-word characters as "_" is a
    model; otherwise "suffix_normalised" when the two agree after GTDB placeholder
    suffixes (_A, _B, ...) are dropped from both and exactly ONE model matches; else
    (None, None). Ambiguous suffix matches are refused rather than guessed.
    """
    exact = "s__" + re.sub(r"[^\w]", "_", str(species).strip())
    if exact in agora2_ec_index:
        return exact, "exact"
    key = _GTDB_SUFFIX_RE.sub("", exact)
    held = _AGORA2_NORMALISED.get(id(agora2_ec_index))
    if held is None or held[0] is not agora2_ec_index:      # ids are reused after GC
        memo = {}
        for sid in agora2_ec_index:
            memo.setdefault(_GTDB_SUFFIX_RE.sub("", sid), []).append(sid)
        held = (agora2_ec_index, memo)
        _AGORA2_NORMALISED[id(agora2_ec_index)] = held
    cands = held[1].get(key, [])
    return (cands[0], "suffix_normalised") if len(cands) == 1 else (None, None)


_AGORA2_NORMALISED: dict = {}


def species_encodes_enzyme(species: str, kegg_id: str, agora2_ec_index: dict):
    """Does `species` (a plain binomial, e.g. "Sutterella wadsworthensis_A")
    encode a terminal biosynthetic enzyme for `kegg_id`, per its AGORA2
    reconstruction?

    Returns True/False if the species has an AGORA2 model, None if it
    doesn't (no evidence -- not the same as False).
    """
    species_id, _ = agora2_species_id(species, agora2_ec_index)
    if species_id is None:
        return None
    species_kegg_ecs = agora2_ec_index[species_id].get(kegg_id, [])
    required = POLYAMINE_EC.get(kegg_id)
    if required is not None:
        return any(ec in species_kegg_ecs for ec in required)
    return len(species_kegg_ecs) > 0


# An evidence source counts toward an evidence tier only where it is informative: its
# base rate -- how often it fires for an arbitrary modelled species -- must be at most
# this (council review, D7). AGORA2 "encodes the enzyme" fires for 87% of modelled
# study species for putrescine and 97% for arginine, so for those compounds it cannot
# distinguish a producer from any other species.
EVIDENCE_MAX_BASE_RATE = 0.20


def agora2_base_rate(species, kegg_id: str, agora2_ec_index: dict):
    """(base rate, n modelled): the fraction of `species` that HAVE an AGORA2 model and
    encode a producing enzyme for `kegg_id`. Species without a model are left out."""
    calls = [species_encodes_enzyme(s, kegg_id, agora2_ec_index) for s in species]
    known = [c for c in calls if c is not None]
    return (float(np.mean(known)) if known else float("nan")), len(known)


def species_has_any_ec(species: str, ecs, agora2_ec_index: dict):
    """Does `species` encode ANY of `ecs` anywhere in its AGORA2 model,
    regardless of which compound the reaction is indexed under? (Direction-
    agnostic: needed for guild ECs such as butyrate kinase, whose reaction
    is written with butyrate as a substrate.)

    Returns True/False if the species has an AGORA2 model, None if not.
    """
    species_id, _ = agora2_species_id(species, agora2_ec_index)
    if species_id is None:
        return None
    wanted = set(ecs)
    return any(wanted.intersection(v) for v in agora2_ec_index[species_id].values())


# An INPUT, built once by tools/build_agora2_index.py. Not under Results/: a run that
# recomputes NB01 archives Results/tables, which would take the index with it.
AGORA2_INDEX_CACHE = DATA_DIR / "reference" / "agora2_species_ec_index.json"
AGORA2_INDEX_SIDECAR = AGORA2_INDEX_CACHE.with_suffix(".build.json")


_AGORA2_INDEX_MEMO: dict | None = None


def load_agora2_ec_index() -> dict:
    """AGORA2 species->EC index from its JSON cache (an input, see the section header).

    Raises if the cache is absent. It never falls back to building the index from the
    1.8 GB of reconstructions: a scan triggered from inside a notebook is what
    overloaded the old server, and the build belongs to tools/build_agora2_index.py.

    Memoised per process: the parsed index is large, so repeated calls within one
    notebook must not load extra copies.
    """
    global _AGORA2_INDEX_MEMO
    if _AGORA2_INDEX_MEMO is not None:
        return _AGORA2_INDEX_MEMO
    if not AGORA2_ENABLED:
        raise RuntimeError("AGORA2 evidence is disabled (utils.AGORA2_ENABLED is False)")
    if not AGORA2_INDEX_CACHE.exists():
        raise RuntimeError(
            f"the AGORA2 EC index is absent ({AGORA2_INDEX_CACHE}). Build it once, outside "
            "the pipeline:  python tools/build_agora2_index.py")
    _AGORA2_INDEX_MEMO = build_agora2_ec_index(None, cache_path=AGORA2_INDEX_CACHE)
    return _AGORA2_INDEX_MEMO


def find_agora2_data_dir(data_dir=DATA_DIR):
    """Locate the extracted AGORA2 species JSON directory (the manifest.csv's
    parent), wherever the QIIME2-artifact extraction happened to unpack it.
    """
    root = Path(data_dir) / "micom" / "agora2_extracted"
    manifests = list(root.rglob("manifest.csv"))
    if not manifests:
        raise FileNotFoundError(f"No AGORA2 manifest.csv found under {root}")
    return manifests[0].parent


def agora2_source_fingerprint(agora2_data_dir) -> dict:
    """Identity of the reconstructions an index was built from: file count and a
    sha256 over the sorted (name, size) list. Cheap (no file contents read); the
    sidecar records it and checks.py D5 compares it with what is on disk."""
    import hashlib as _hashlib
    files = sorted(Path(agora2_data_dir).glob("*.json"))
    h = _hashlib.sha256()
    for p in files:
        h.update(f"{p.name}\t{p.stat().st_size}\n".encode("utf-8"))
    return {"n_files": len(files), "names_sizes_sha256": h.hexdigest()}


# ============================================================================
# PATSY FORMULA SAFETY
# ============================================================================

def patsy_safe_name(col: str) -> str:
    """Wrap a column name in Q("...") for statsmodels/patsy formulas if it
    contains anything other than letters/digits/underscore (not just
    spaces -- a "." as in "Study.Group" is just as unsafe and was the
    actual cause of every OLS fit silently failing in the old confound
    diagnostic, since that check only looked for spaces).
    """
    if re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", col):
        return col
    return f'Q("{col}")'


def build_formula(outcome: str, predictors: list[str]) -> str:
    """Build a safe OLS-style formula string, quoting every unsafe name."""
    return f"{patsy_safe_name(outcome)} ~ " + " + ".join(patsy_safe_name(p) for p in predictors)


# ============================================================================
# PLOTTING HELPERS
# ============================================================================

def _draw_confidence_ellipse(ax, x, y, color, n_std=1.96, **kwargs):
    from matplotlib.patches import Ellipse
    if len(x) < 3:
        return
    cov = np.cov(x, y)
    eigenvalues, eigenvectors = np.linalg.eigh(cov)
    order = eigenvalues.argsort()[::-1]
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]
    angle = np.degrees(np.arctan2(eigenvectors[1, 0], eigenvectors[0, 0]))
    width, height = 2 * n_std * np.sqrt(np.maximum(eigenvalues, 0))
    ellipse = Ellipse(
        xy=(np.mean(x), np.mean(y)), width=width, height=height,
        angle=angle, facecolor=color, alpha=0.15,
        edgecolor=color, linewidth=1.5, linestyle="--", **kwargs
    )
    ax.add_patch(ellipse)


def _save_figure(fig, save_path):
    """The plotting helpers below save through viz.save_figure: a 300-dpi PNG plus a
    PDF twin (checks.py E12 forbids raw savefig). Imported lazily so loading utils does
    not install the figure theme's machinery."""
    from viz import save_figure
    save_figure(fig, save_path)


def plot_pca(df, meta, title, color_col="Study.Group", n_components=2,
             figsize=(12, 8), save_path=None, show_ellipses=True,
             show_variance_bar=True):
    from matplotlib.gridspec import GridSpec

    feature_cols = [c for c in df.columns if c != "Sample"]
    X = df[feature_cols].apply(pd.to_numeric, errors="coerce").values
    X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)

    n_comp_full = min(10, min(X.shape))
    pca_full = PCA(n_components=n_comp_full)
    pca_full.fit(X)

    pca = PCA(n_components=n_components)
    coords = pca.fit_transform(X)

    plot_df = pd.DataFrame({"PC1": coords[:, 0], "PC2": coords[:, 1], "Sample": df["Sample"].values})
    plot_df = plot_df.merge(meta[["Sample", color_col]], on="Sample", how="left")

    if show_variance_bar:
        fig = plt.figure(figsize=figsize)
        gs = GridSpec(1, 5, figure=fig, width_ratios=[4, 0.1, 0.8, 0.1, 0.1])
        ax = fig.add_subplot(gs[0, 0])
        ax_var = fig.add_subplot(gs[0, 2])
    else:
        fig, ax = plt.subplots(figsize=figsize)
        ax_var = None

    # Colour by IDENTITY, not by rank. This used to be
    #     palette = sns.color_palette("Set2", n_colors=nunique)
    #     color_map = {g: palette[i] for i, g in enumerate(sorted(unique))}
    # so a group's colour depended on which OTHER groups happened to be present:
    # the three plot_pca loops in NB01 run over cohorts with different group sets,
    # and were therefore not colour-comparable with each other. viz.group_color
    # looks each group up on the fixed severity scale instead.
    import viz

    groups = sorted(plot_df[color_col].dropna().unique())
    color_map = viz.group_palette(groups)

    for group in groups:
        subset = plot_df[plot_df[color_col] == group]
        ax.scatter(subset["PC1"], subset["PC2"], c=[color_map[group]],
                   label=f"{group} (n={len(subset)})", alpha=0.7, s=50,
                   edgecolors="white", linewidth=0.5)
        if show_ellipses and len(subset) >= 3:
            _draw_confidence_ellipse(ax, subset["PC1"].values, subset["PC2"].values, color_map[group])

    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.1%} variance)")
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.1%} variance)")
    ax.set_title(title)
    ax.legend(bbox_to_anchor=(1.0, 1), loc="upper left", fontsize=8, framealpha=0.9)

    if ax_var is not None:
        var_ratios = pca_full.explained_variance_ratio_[:n_comp_full] * 100
        ax_var.barh(range(len(var_ratios)), var_ratios, color="steelblue", edgecolor="white")
        ax_var.set_yticks(range(len(var_ratios)))
        ax_var.set_yticklabels([f"PC{i+1}" for i in range(len(var_ratios))], fontsize=7)
        ax_var.set_xlabel("% Var", fontsize=8)
        ax_var.set_title("Scree", fontsize=9)
        ax_var.invert_yaxis()
        for i, v in enumerate(var_ratios):
            ax_var.text(v + 0.3, i, f"{v:.1f}", va="center", fontsize=6)

    plt.tight_layout()
    if save_path:
        _save_figure(fig, save_path)
    return fig, ax


def volcano_plot(results_df, title, fc_thresh=1.0, q_thresh=0.05,
                 fc_col="MeanDiff", highlight_features=None,
                 figsize=(10, 8), save_path=None, ax=None):
    """Volcano plot from differential_abundance() results.

    `fc_col` defaults to "MeanDiff" (see differential_abundance docstring);
    pass a different column name if plotting a differently-computed table.
    """
    _standalone = ax is None
    if _standalone:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()

    neg_log_q = -np.log10(results_df["QValue"].clip(lower=1e-300))
    sig_mask = (results_df["QValue"] < q_thresh) & (results_df[fc_col].abs() > fc_thresh)
    colors = np.where(sig_mask, "red", "grey")

    ax.scatter(results_df[fc_col], neg_log_q, c=colors, alpha=0.5, s=10)
    ax.axhline(-np.log10(q_thresh), ls="--", color="blue", alpha=0.5)
    ax.axvline(fc_thresh, ls="--", color="blue", alpha=0.5)
    ax.axvline(-fc_thresh, ls="--", color="blue", alpha=0.5)

    if highlight_features:
        for feat in highlight_features:
            row = results_df[results_df["Feature"] == feat]
            if not row.empty:
                x = row[fc_col].values[0]
                y = -np.log10(row["QValue"].values[0])
                ax.annotate(feat, (x, y), fontsize=7, alpha=0.8)

    n_sig = int(sig_mask.sum())
    ax.set_xlabel(fc_col)
    ax.set_ylabel("-log10(Q-value)")
    ax.set_title(f"{title}\n(n_sig={n_sig})", fontsize=9 if not _standalone else 11)

    if _standalone:
        plt.tight_layout()
        if save_path:
            _save_figure(fig, save_path)
    return fig, ax


# plot_detection_histogram() was removed in v2.1: it had ZERO call sites anywhere,
# while 01_preprocessing.ipynb re-implemented it inline (a 50-bin hist of
# detection_rate plus a red threshold axvline). Nothing to migrate.


def plot_correlation_heatmap(corr_df, title, top_n=20, figsize=(14, 10), save_path=None):
    if corr_df.empty:
        return None, None

    top = corr_df.nsmallest(top_n * 5, "QValue") if len(corr_df) > top_n * 5 else corr_df
    pivot = top.pivot_table(index="Species", columns="Metabolite", values="Rho", aggfunc="first")

    if pivot.empty:
        return None, None

    if pivot.shape[0] > top_n:
        row_vars = pivot.var(axis=1).nlargest(top_n).index
        pivot = pivot.loc[row_vars]

    import viz

    # Rho is a SIGNED quantity -> diverging scale with a neutral midpoint at 0.
    # (Was a hardcoded "RdBu_r"; viz.diverging_cmap is the single definition.)
    fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(pivot, cmap=viz.diverging_cmap(), center=0, annot=True, fmt=".2f",
                ax=ax, xticklabels=True, yticklabels=True)
    ax.set_title(title)
    plt.tight_layout()
    if save_path:
        _save_figure(fig, save_path)
    return fig, ax


def plot_stage_trend(species_name, metabolite_name,
                     species_df, mtb_df, meta_df,
                     stage_order=CRC_STAGE_ORDER,
                     stage_col="CRC.Stage.Numeric",
                     group_col="Study.Group",
                     figsize=(12, 5), save_path=None):
    common = (
        set(species_df["Sample"]) & set(mtb_df["Sample"]) & set(meta_df["Sample"])
    )
    spc = species_df.set_index("Sample").loc[list(common), species_name]
    mtb = mtb_df.set_index("Sample").loc[list(common), metabolite_name]
    grp = meta_df.set_index("Sample").loc[list(common), group_col]

    plot_df = pd.DataFrame({"Species": spc, "Metabolite": mtb, "Stage": grp}).dropna()
    plot_df = plot_df[plot_df["Stage"].isin(stage_order)]
    plot_df["Stage"] = pd.Categorical(plot_df["Stage"], categories=stage_order, ordered=True)
    plot_df = plot_df.sort_values("Stage")

    import viz

    fig, axes = plt.subplots(1, 2, figsize=figsize)
    # Stage is ORDINAL MAGNITUDE -> sequential ramp. This was
    # sns.color_palette("coolwarm", ...), a DIVERGING map used as a severity
    # palette: it implies a meaningful midpoint and two opposed poles, and its
    # pale middle made mid-stage boxes vanish against the panel.
    palette = viz.severity_colors(order=stage_order)

    for ax, col, label in zip(axes, ["Species", "Metabolite"],
                              [species_name[:40], metabolite_name[:40]]):
        sns.boxplot(data=plot_df, x="Stage", y=col, order=stage_order,
                    palette=palette, ax=ax, linewidth=0.8, fliersize=2)
        sns.stripplot(data=plot_df, x="Stage", y=col, order=stage_order,
                      color="black", alpha=0.3, size=3, ax=ax)
        ax.set_title(label, fontsize=9)
        ax.set_xlabel("")
        ax.set_ylabel("Transformed abundance")
        ax.tick_params(axis="x", rotation=40)

    fig.suptitle(f"{species_name[:35]}  <->  {metabolite_name[:35]}", fontsize=10, y=1.01)
    plt.tight_layout()
    if save_path:
        _save_figure(fig, save_path)
    return fig, axes


# ============================================================================
# VALIDATION HELPERS
# ============================================================================

def validate_sample_alignment(data, dataset_name):
    meta_samples = set(data[dataset_name]["metadata"]["Sample"])
    mtb_samples = set(data[dataset_name]["mtb"]["Sample"])
    spc_samples = set(data[dataset_name]["species"]["Sample"])

    return {
        "dataset": dataset_name,
        "meta_n": len(meta_samples),
        "mtb_n": len(mtb_samples),
        "species_n": len(spc_samples),
        "meta_only": meta_samples - mtb_samples - spc_samples,
        "mtb_only": mtb_samples - meta_samples,
        "spc_only": spc_samples - meta_samples,
        "common": meta_samples & mtb_samples & spc_samples,
    }


def validate_metabolite_alignment(data, dataset_name):
    mtb_cols = set(data[dataset_name]["mtb"].columns) - {"Sample"}
    mtb_map = data[dataset_name]["mtb.map"]
    if mtb_map is None:
        return {
            "dataset": dataset_name, "mtb_features": len(mtb_cols),
            "map_compounds": 0, "in_mtb_not_map": mtb_cols, "in_map_not_mtb": set(),
        }
    map_compounds = set(mtb_map["Compound"])
    return {
        "dataset": dataset_name,
        "mtb_features": len(mtb_cols),
        "map_compounds": len(map_compounds),
        "in_mtb_not_map": mtb_cols - map_compounds,
        "in_map_not_mtb": map_compounds - mtb_cols,
    }


def print_dataset_summary(data):
    rows = []
    for ds in DATASETS:
        if ds not in data:
            continue
        rows.append({
            "Dataset": ds,
            "Samples (meta)": len(data[ds]["metadata"]),
            "Samples (mtb)": len(data[ds]["mtb"]),
            "Samples (species)": len(data[ds]["species"]),
            "Metabolites": len(data[ds]["mtb"].columns) - 1,
            "Species": len(data[ds]["species"].columns) - 1,
            "Map entries": len(data[ds]["mtb.map"]) if data[ds]["mtb.map"] is not None else 0,
        })
    return pd.DataFrame(rows)
