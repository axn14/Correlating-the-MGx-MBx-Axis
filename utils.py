"""
utils.py  —  Generalised Metabolite–Metagenomics Correlation Pipeline
======================================================================
Publication notebooks NB01–NB06.

Design principles
-----------------
* FULLY GENERALISED: all features (species, metabolites) are treated
  equally. No pathway-based pre-filtering is applied at any stage.
* Pathway labels (e.g., polyamine, SCFA, bile acid) are applied
  POST-DISCOVERY in NB02 as optional metadata annotations only.
* Cohort handling: ALL 7 datasets available; primary analysis cohort
  is YACHIDA-CRC-2019 (stage-annotated); SINHA_CRC_2016 and
  Kim_adenomas_2020 used for replication in NB02 only.
* ML (NB03) runs on YACHIDA only to preserve stage stratification.
* GutSMASH benchmark (NB06) uses top-50 SHAP-ranked producers from NB03.
"""

from __future__ import annotations

import warnings
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr, mannwhitneyu, rankdata, kruskal
from statsmodels.stats.multitest import multipletests
from sklearn.decomposition import PCA
from sklearn.preprocessing import LabelEncoder
import re

warnings.filterwarnings("ignore")

# ============================================================================
# PATHS  (Windows paths — WSL notebooks resolve via drive mapping)
# ============================================================================
DATA_DIR    = Path(r"C:\Users\andna\Documents\KI\Data")
RESULTS_DIR = Path(r"C:\Users\andna\Documents\KI\Results")
INTER_DIR   = RESULTS_DIR / "intermediate"
FIG_DIR     = RESULTS_DIR / "figures"
TABLE_DIR   = RESULTS_DIR / "tables"

# ============================================================================
# DATASET REGISTRY
# ============================================================================
DATASETS_ALL = [
    "ERAWIJANTARI-GASTRIC-CANCER-2020",
    "FRANZOSA-IBD-2019",
    "iHMP-IBDMDB-2019",
    "MARS-IBS-2020",
    "WANG_ESRD_2020",
    "YACHIDA-CRC-2019",
    "SINHA_CRC_2016",
    "Kim_adenomas_2020",
]

DATASET_PRIMARY   = "YACHIDA-CRC-2019"
DATASET_SECONDARY = "SINHA_CRC_2016"

# Cohorts with generalised inflammation — excluded from all CRC analyses
IBD_COHORTS = frozenset(["FRANZOSA-IBD-2019", "iHMP-IBDMDB-2019"])

# Cohorts available for Leave-One-Dataset-Out cross-validation
# All must have shotgun metagenomics + paired LC-MS metabolomics at species resolution
DATASETS_LODO = [
    "YACHIDA-CRC-2019",
    "KONG_EOCRC_2023",
    "SINHA_VOGTMANN_SHOTGUN",
    "BORENSTEIN_CRC",
]

# Remove IBD cohorts at source; register new LODO cohorts
DATASETS_ALL = [d for d in DATASETS_ALL if d not in IBD_COHORTS] + [
    "KONG_EOCRC_2023",
    "SINHA_VOGTMANN_SHOTGUN",
    "BORENSTEIN_CRC",
]

FILE_SUFFIXES = {
    "metadata": " metadata.tsv",
    "mtb":      " mtb.tsv",
    "mtb.map":  " mtb.map.tsv",
    "species":  " species.tsv",
}
KIM_SUFFIXES = {
    "metadata": " metadata.tsv",
    "mtb":      " mtb.tsv",
    "mtb.map":  " mtb_map.tsv",
    "species":  " species.tsv",
}

# ============================================================================
# STAGE ORDERINGS  (YACHIDA-specific)
# ============================================================================
CRC_STAGE_ORDER  = ["Healthy", "Stage_0", "HS", "Stage_I_II", "Stage_III_IV", "MP"]
CRC_STAGE_NUMERIC = {s: i for i, s in enumerate(CRC_STAGE_ORDER)}

THREE_GROUP_ORDER = ["Healthy", "Early_CRC", "Advanced_CRC"]
THREE_GROUP_MAP   = {
    "Healthy":      "Healthy",
    "Stage_0":      "Early_CRC",
    "HS":           "Early_CRC",
    "Stage_I_II":   "Advanced_CRC",
    "Stage_III_IV": "Advanced_CRC",
    "MP":           "Advanced_CRC",
}

# ============================================================================
# POST-DISCOVERY PATHWAY REFERENCE KEGG IDs
# ── Applied ONLY to significant hits after discovery. Never used as filters.
# ============================================================================
_PATHWAY_KEGG_REF: dict[str, dict[str, str]] = {
    "Polyamine": {
        "C00134": "Putrescine",    "C00315": "Spermidine",
        "C00750": "Spermine",      "C00179": "Agmatine",
        "C01672": "Cadaverine",    "C00077": "Ornithine",
        "C00062": "Arginine",
    },
    "SCFA": {
        "C00033": "Acetate",       "C00163": "Propionate",
        "C00246": "Butyrate",      "C01134": "Valerate",
    },
    "Bile Acid": {
        "C05122": "Chenodeoxycholic acid", "C04483": "Taurocholic acid",
        "C05466": "Glycocholic acid",      "C00695": "Cholic acid",
        "C04515": "Deoxycholic acid",      "C17231": "Lithocholic acid",
    },
    "Amino Acid": {
        "C00073": "Methionine",    "C00064": "Glutamine",
        "C00025": "Glutamate",     "C00037": "Glycine",
        "C00049": "Aspartate",     "C00041": "Alanine",
        "C00047": "Lysine",        "C00148": "Proline",
        "C00188": "Threonine",     "C00123": "Leucine",
        "C00407": "Isoleucine",
    },
    "Nucleotide": {
        "C00144": "GMP",  "C00020": "AMP",  "C00105": "UMP",
        "C00055": "CMP",  "C00459": "dTMP",
    },
    "Lipid / Fatty Acid": {
        "C00249": "Palmitate",      "C00712": "Oleate",
        "C00350": "Phosphatidylethanolamine", "C00157": "Phosphatidylcholine",
    },
    "Indole / Tryptophan": {
        "C00078": "Tryptophan",     "C00398": "Indole",
        "C02693": "Indole-3-acetate", "C11310": "Indole-3-propionate",
    },
}


def annotate_pathway(kegg_id: str) -> str:
    """Post-discovery pathway annotation. Returns pathway label or 'Other'.
    Handles both bare KEGG IDs ('C00134') and KEGGID_Name format ('C00134_Putrescine').
    """
    kid = kegg_id.split("_")[0] if "_" in kegg_id else kegg_id
    for pathway, kegg_dict in _PATHWAY_KEGG_REF.items():
        if kid in kegg_dict:
            return pathway
    return "Other"


def pathway_kegg_ids(pathway: str) -> dict[str, str]:
    """Return {KEGG_ID: name} for a named pathway (case-insensitive)."""
    for k, v in _PATHWAY_KEGG_REF.items():
        if k.lower() == pathway.lower():
            return v
    return {}


def extract_genus(species_str: str) -> str:
    """
    Extract genus from MetaPhlAn-style or GTDB species string.
    'k__Bacteria|...|s__Fusobacterium_A_mortiferum' → 'Fusobacterium'
    """
    part = species_str.split("|")[-1]
    part = part.split("__")[-1]
    genus = re.split(r"[_\s]", part)[0]
    return genus


# ============================================================================
# QUALITY CONTROL THRESHOLDS
# ============================================================================
PREVALENCE_SPE = 0.10    # species: ≥10% of samples non-zero
PREVALENCE_MTB = 0.15    # metabolites: ≥15% (sparser; stricter reduces FDR inflation)
PREVALENCE_MIN = 0.10    # generic fallback
FDR_THRESHOLD  = 0.05    # BH-corrected q-value threshold
MIN_CORR       = 0.20    # minimum |ρ| pre-filter before BH correction

MAX_SPECIES_NB02 = 500   # top-variance species retained for correlation analysis
MAX_MTB_NB02     = 150   # top-variance metabolites retained for correlation analysis
N_TOP_TARGETS    = 50    # top dysregulated metabolite targets for ML (NB03)

# ============================================================================
# DATA LOADING
# ============================================================================

def _suffixes_for(dataset: str) -> dict:
    return KIM_SUFFIXES if dataset.startswith("Kim") else FILE_SUFFIXES


def load_dataset(dataset: str, data_dir: Path = DATA_DIR) -> dict[str, pd.DataFrame]:
    """Load metadata, mtb, mtb.map and species TSVs for a single dataset."""
    sfx = _suffixes_for(dataset)
    out = {}
    for key, suffix in sfx.items():
        fpath = data_dir / f"{dataset}{suffix}"
        if not fpath.exists():
            warnings.warn(f"Missing file: {fpath}")
            out[key] = pd.DataFrame()
            continue
        idx_col = 1 if key == "metadata" else 0
        out[key] = pd.read_csv(fpath, sep="\t", index_col=idx_col, low_memory=False)
        if key == "mtb":
            out[key] = out[key].apply(pd.to_numeric, errors="coerce")
    return out


def load_all_datasets(datasets: list[str] = DATASETS_ALL,
                      data_dir: Path = DATA_DIR) -> dict[str, dict]:
    """Load all datasets; silently skip those with missing files."""
    data = {}
    for ds in datasets:
        try:
            data[ds] = load_dataset(ds, data_dir)
        except Exception as e:
            warnings.warn(f"Could not load {ds}: {e}")
    print(f"Loaded {len(data)} datasets: {list(data.keys())}")
    return data


# ============================================================================
# METADATA HARMONISATION
# ============================================================================

def harmonize_metadata(meta: pd.DataFrame, dataset: str) -> pd.DataFrame:
    """
    Adds harmonized columns:
      study_group  — original group label
      Stage.6      — 6-level CRC stage (YACHIDA only, else NaN)
      Stage.3Group — Healthy / Early_CRC / Advanced_CRC
      Stage.Num    — ordinal integer (YACHIDA only)
      Dataset      — dataset name
    """
    meta = meta.copy()
    meta["Dataset"] = dataset

    if dataset == "YACHIDA-CRC-2019":
        sg = meta.get("Study.Group", pd.Series(dtype=str))
        meta["study_group"] = sg
        meta["Stage.6"]     = sg.map(lambda x: x if x in CRC_STAGE_ORDER else np.nan)
        meta["Stage.Num"]   = meta["Stage.6"].map(CRC_STAGE_NUMERIC)
        meta["Stage.3Group"] = meta["Stage.6"].map(THREE_GROUP_MAP)

    elif dataset == "SINHA_CRC_2016":
        sg = meta.get("Study.Group", pd.Series(dtype=str))
        meta["study_group"] = sg
        meta["Stage.6"]     = np.nan
        meta["Stage.Num"]   = np.nan
        meta["Stage.3Group"] = sg.map(
            lambda x: "Healthy" if str(x).strip() in ("0", "Healthy", "healthy")
                      else "Advanced_CRC")

    elif dataset == "Kim_adenomas_2020":
        sg = meta.get("Group", meta.get("Study.Group", pd.Series(dtype=str)))
        meta["study_group"] = sg
        meta["Stage.6"]     = np.nan
        meta["Stage.Num"]   = np.nan
        meta["Stage.3Group"] = sg.map(
            lambda x: "Healthy"    if "control" in str(x).lower()
                      else "Early_CRC" if "adenoma" in str(x).lower()
                      else "Advanced_CRC")
    elif dataset == "KONG_EOCRC_2023":
        sg = meta.get("Study.Group", pd.Series(dtype=str))
        meta["study_group"] = sg
        meta["Stage.6"]     = sg.map(lambda x: x if x in CRC_STAGE_ORDER else np.nan)
        meta["Stage.Num"]   = meta["Stage.6"].map(CRC_STAGE_NUMERIC)
        meta["Stage.3Group"] = meta["Stage.6"].map(THREE_GROUP_MAP)

    elif dataset in ("SINHA_VOGTMANN_SHOTGUN", "BORENSTEIN_CRC"):
        sg = meta.get("Study.Group", pd.Series(dtype=str))
        meta["study_group"] = sg
        meta["Stage.6"]     = np.nan
        meta["Stage.Num"]   = np.nan
        meta["Stage.3Group"] = sg.map(
            lambda x: "Healthy" if str(x).strip() in ("0", "Healthy", "healthy", "Control", "control")
                      else "Advanced_CRC")

    else:
        sg = meta.get("Study.Group", pd.Series(dtype=str))
        meta["study_group"] = sg
        meta["Stage.6"]     = np.nan
        meta["Stage.Num"]   = np.nan
        meta["Stage.3Group"] = np.nan

    return meta


# ============================================================================
# SAMPLE ALIGNMENT VALIDATION
# ============================================================================

def validate_sample_alignment(meta: pd.DataFrame,
                               mtb:  pd.DataFrame,
                               spe:  pd.DataFrame) -> dict:
    """Return shared sample sets and counts (all three modalities)."""
    s_meta = set(meta.index)
    s_mtb  = set(mtb.index)
    s_spe  = set(spe.index)
    shared = s_meta & s_mtb & s_spe
    return {
        "n_shared_samples": len(shared),
        "n_meta_only":      len(s_meta - s_mtb - s_spe),
        "n_mtb_only":       len(s_mtb  - s_meta - s_spe),
        "n_species_only":   len(s_spe  - s_meta - s_mtb),
        "shared":           shared,
    }


# ============================================================================
# METABOLITE MAP UTILITY
# ============================================================================

def build_metabolite_name_map(mtb_map: pd.DataFrame) -> dict[str, str]:
    """Return {KEGG_ID: human_readable_name} from the mtb.map table."""
    if mtb_map.empty:
        return {}
    col = mtb_map.columns[0]
    return mtb_map[col].to_dict()


# ============================================================================
# SPECIES NAME REDUCTION
# ============================================================================

def reduce_species_names(species_df: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """
    Collapse full MetaPhlAn taxonomy strings to species-level short names.
    Sums duplicates that arise from collapse.
    Returns (reduced_df, {short_name: [original_cols]}).
    """
    col_map: dict[str, list[str]] = {}
    for col in species_df.columns:
        short = col.split("|")[-1]
        short = short.split("__")[-1] if "__" in short else short
        col_map.setdefault(short, []).append(col)

    data: dict[str, np.ndarray] = {}
    for short, orig_cols in col_map.items():
        arr = species_df[orig_cols].values.astype(float)
        data[short] = arr.sum(axis=1) if len(orig_cols) > 1 else arr[:, 0]

    return pd.DataFrame(data, index=species_df.index), col_map


# ============================================================================
# SAMPLE-LEVEL QC  (includes Shannon diversity)
# ============================================================================

def compute_sample_qc(df: pd.DataFrame, data_type: str = "metabolite") -> pd.DataFrame:
    """
    Compute per-sample QC metrics.
    For species: also computes Shannon alpha diversity and row_sum.
    row_sum for MetaPhlAn data (relative abundances in [0,1]) should be ≈ 1.0.
    """
    from scipy.stats import entropy as _entropy
    x = df.values.astype(float)

    qc = pd.DataFrame(index=df.index)
    qc["total_signal"]   = np.nansum(x, axis=1)
    qc["n_detected"]     = (x > 0).sum(axis=1)
    qc["n_features"]     = x.shape[1]
    qc["detection_rate"] = qc["n_detected"] / qc["n_features"]

    if data_type == "species":
        shannon = np.array([
            _entropy(row[row > 0]) if (row > 0).any() else 0.0
            for row in x
        ])
        qc["shannon_diversity"] = shannon
        # row_sum = total_signal for non-negative data; MetaPhlAn profiles sum to ~1.0
        qc["row_sum"] = qc["total_signal"]

    return qc


def compute_feature_qc(df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute per-feature QC metrics: detection rate, mean, std, CV, variance.
    """
    x = df.values.astype(float)
    fqc = pd.DataFrame(index=df.columns)
    fqc["detection_rate"] = (x > 0).mean(axis=0)
    fqc["mean"]           = np.nanmean(x, axis=0)
    fqc["std"]            = np.nanstd(x, axis=0)
    fqc["cv"]             = np.where(fqc["mean"] > 0, fqc["std"] / fqc["mean"], np.nan)
    fqc["variance"]       = np.nanvar(x, axis=0)
    return fqc


# ============================================================================
# SHANNON DIVERSITY STATISTICAL TEST
# ============================================================================

def shannon_kruskal_test(qc_df: pd.DataFrame,
                          meta: pd.DataFrame,
                          group_col: str = "Stage.3Group",
                          groups: list[str] | None = None) -> pd.DataFrame:
    """
    Kruskal-Wallis test on Shannon diversity across disease groups,
    followed by pairwise Mann-Whitney U with BH correction.

    Parameters
    ----------
    qc_df     : sample-level QC DataFrame with 'shannon_diversity' column
    meta      : metadata DataFrame with group column
    group_col : column in meta defining groups
    groups    : list of groups to test (default: THREE_GROUP_ORDER)

    Returns
    -------
    DataFrame with kruskal_stat, kruskal_pval, and pairwise results.
    """
    from scipy.stats import kruskal, mannwhitneyu
    if "shannon_diversity" not in qc_df.columns:
        return pd.DataFrame()

    if groups is None:
        groups = [g for g in THREE_GROUP_ORDER if g in meta[group_col].values]

    # Align indices
    common = qc_df.index.intersection(meta.index)
    sh = qc_df.loc[common, "shannon_diversity"]
    grp = meta.loc[common, group_col]

    group_vals = [sh[grp == g].dropna().values for g in groups if (grp == g).sum() > 2]
    group_lbls = [g for g in groups if (grp == g).sum() > 2]

    if len(group_vals) < 2:
        return pd.DataFrame()

    kw_stat, kw_pval = kruskal(*group_vals)

    # Pairwise Mann-Whitney U
    pairs = []
    for i in range(len(group_lbls)):
        for j in range(i + 1, len(group_lbls)):
            g1, g2 = group_lbls[i], group_lbls[j]
            u_stat, u_pval = mannwhitneyu(group_vals[i], group_vals[j],
                                           alternative="two-sided")
            pairs.append({"group1": g1, "group2": g2, "U": u_stat, "pval_raw": u_pval,
                           "n1": len(group_vals[i]), "n2": len(group_vals[j])})

    result_df = pd.DataFrame(pairs)
    if not result_df.empty:
        _, qvals, _, _ = multipletests(result_df["pval_raw"], method="fdr_bh")
        result_df["qval"] = qvals

    result_df.attrs["kruskal_stat"] = float(kw_stat)
    result_df.attrs["kruskal_pval"] = float(kw_pval)
    return result_df


# ============================================================================
# QC FILTERING
# ============================================================================

def qc_filter_species(df: pd.DataFrame,
                      prevalence_min: float = PREVALENCE_SPE,
                      min_abundance: float = 1e-4) -> pd.DataFrame:
    """Filter species by prevalence and minimum abundance.

    min_abundance: species whose max relative abundance across all samples
    is below this threshold are removed (default 0.01% — MetaPhlAn detection floor).
    """
    prev = (df > 0).mean(axis=0)
    df   = df.loc[:, prev >= prevalence_min]
    if min_abundance > 0:
        df = df.loc[:, df.max(axis=0) >= min_abundance]
    return df


def qc_filter_metabolites(df: pd.DataFrame,
                           prevalence_min: float = PREVALENCE_MTB,
                           var_quantile: float = 0.10) -> pd.DataFrame:
    """Filter metabolites by prevalence then remove bottom var_quantile by variance."""
    prev = (df > 0).mean(axis=0)
    df   = df.loc[:, prev >= prevalence_min]
    if var_quantile > 0:
        var  = df.var(axis=0)
        df   = df.loc[:, var > var.quantile(var_quantile)]
    return df


def qc_filter(df: pd.DataFrame,
              prevalence_min: float = PREVALENCE_MIN,
              min_abundance: float = 0) -> pd.DataFrame:
    """Generic QC filter. Prefer qc_filter_species/qc_filter_metabolites."""
    prev = (df > 0).mean(axis=0)
    df   = df.loc[:, prev >= prevalence_min]
    if min_abundance > 0:
        df = df.loc[:, df.max(axis=0) >= min_abundance]
    return df


# ============================================================================
# TRANSFORMS
# ============================================================================

def clr_transform(df: pd.DataFrame, pseudocount: float = 1e-4) -> pd.DataFrame:
    """
    Centered log-ratio (CLR) transform for compositional species data.
    pseudocount = 1e-4 keeps absent-species CLR at a biologically plausible level.
    """
    x    = df.values.copy().astype(float)
    x    = np.where(x == 0, pseudocount, x)
    logx = np.log(x)
    clr  = logx - logx.mean(axis=1, keepdims=True)
    return pd.DataFrame(clr, index=df.index, columns=df.columns)


def log10_transform(df: pd.DataFrame) -> pd.DataFrame:
    """
    Three-step metabolite normalisation:
    1. Half-minimum imputation per feature.
    2. log10(x + 1).
    3. Sample-wise median centering.
    """
    x = df.values.copy().astype(float)
    for j in range(x.shape[1]):
        col = x[:, j]
        pos_vals = col[col > 0]
        half_min = pos_vals.min() / 2.0 if len(pos_vals) > 0 else 1e-4
        col[~np.isfinite(col) | (col <= 0)] = half_min
        x[:, j] = col

    x = np.log10(x + 1.0)
    row_medians = np.nanmedian(x, axis=1, keepdims=True)
    x -= row_medians
    return pd.DataFrame(x, index=df.index, columns=df.columns)


# ============================================================================
# DIMENSIONALITY REDUCTION: PCoA (Bray-Curtis + Aitchison)
# ============================================================================

def bray_curtis_pcoa(df: pd.DataFrame,
                     n_components: int = 10) -> dict:
    """
    Principal Coordinates Analysis (PCoA) using Bray-Curtis dissimilarity.
    Input df should be non-negative (raw relative abundances, NOT CLR-transformed).
    Uses classical (metric) MDS via eigendecomposition of the doubly-centred
    distance matrix.

    Returns dict with keys:
      coords        — DataFrame (samples × PCo1..PCok)
      var_explained — array of proportion variance explained per axis
      eigvals       — positive eigenvalues
      dist_mat      — n×n Bray-Curtis distance matrix
    """
    from scipy.spatial.distance import pdist, squareform

    x = df.values.copy().astype(float)
    # Bray-Curtis requires non-negative data; shift if needed (e.g. median-centred data)
    x_shift = x - x.min(axis=0, keepdims=True)
    x_shift = np.where(~np.isfinite(x_shift), 0.0, x_shift)

    dist_mat = squareform(pdist(x_shift, metric="braycurtis"))
    n = len(df)

    # Classical MDS: doubly centre D²
    H = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * H @ (dist_mat ** 2) @ H

    eigvals, eigvecs = np.linalg.eigh(B)
    order   = np.argsort(eigvals)[::-1]
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]

    pos     = eigvals > 1e-9
    k       = min(n_components, int(pos.sum()))
    coords  = eigvecs[:, :k] * np.sqrt(eigvals[:k])
    var_exp = eigvals[:k] / eigvals[pos].sum()

    return {
        "coords":       pd.DataFrame(coords,
                                     index=df.index,
                                     columns=[f"PCo{i+1}" for i in range(k)]),
        "var_explained": var_exp,
        "eigvals":      eigvals[:k],
        "dist_mat":     dist_mat,
    }


def aitchison_pcoa(df_clr: pd.DataFrame,
                   n_components: int = 10) -> dict:
    """
    Aitchison PCoA = PCA on CLR-transformed data (Aitchison distance = CLR Euclidean).
    Input should already be CLR-transformed.
    Returns same dict structure as bray_curtis_pcoa.
    """
    from scipy.spatial.distance import pdist, squareform

    pca = PCA(n_components=min(n_components, df_clr.shape[1], df_clr.shape[0] - 1),
              random_state=42)
    coords = pca.fit_transform(df_clr.values)
    dist_mat = squareform(pdist(df_clr.values, metric="euclidean"))

    return {
        "coords":       pd.DataFrame(coords,
                                     index=df_clr.index,
                                     columns=[f"PCo{i+1}" for i in range(coords.shape[1])]),
        "var_explained": pca.explained_variance_ratio_,
        "eigvals":      pca.explained_variance_,
        "dist_mat":     dist_mat,
    }


# ============================================================================
# MANTEL TEST
# ============================================================================

def mantel_test(D1: np.ndarray,
                D2: np.ndarray,
                n_perm: int = 999,
                random_state: int = 42) -> tuple[float, float]:
    """
    Mantel test: Pearson correlation between upper triangles of two
    distance matrices, tested by row/column permutation.

    Parameters
    ----------
    D1, D2       : symmetric n×n distance matrices (numpy arrays)
    n_perm       : number of permutations (default 999)
    random_state : RNG seed

    Returns
    -------
    (r, p_value)  — Pearson r and permutation-based p-value (two-tailed)
    """
    from scipy.stats import pearsonr
    rng = np.random.default_rng(random_state)
    n   = D1.shape[0]
    idx = np.triu_indices(n, k=1)

    d1_flat = D1[idx].astype(float)
    d2_flat = D2[idx].astype(float)

    r_obs, _ = pearsonr(d1_flat, d2_flat)

    count = 0
    for _ in range(n_perm):
        perm    = rng.permutation(n)
        D1_perm = D1[perm][:, perm]
        r_perm, _ = pearsonr(D1_perm[idx], d2_flat)
        if abs(r_perm) >= abs(r_obs):
            count += 1

    p_val = (count + 1) / (n_perm + 1)
    return float(r_obs), float(p_val)


# ============================================================================
# REDUNDANCY ANALYSIS (RDA)
# ============================================================================

def rda_analysis(Y: pd.DataFrame,
                 X: pd.DataFrame,
                 n_components: int = 2) -> dict:
    """
    Simple Redundancy Analysis (RDA).
    Regresses each column of Y on X (standardised), then runs PCA on
    the fitted (constrained) values.

    Parameters
    ----------
    Y : samples × response variables (species or metabolites), already transformed
    X : samples × explanatory variables (e.g. one-hot encoded Stage.3Group)
    n_components : number of RDA axes to return

    Returns
    -------
    dict with keys:
      rda_scores         — DataFrame (samples × RDA1..RDAk)
      r2_total           — total constrained variance / total variance
      explained_variance — per-axis proportion of constrained variance
      r2_adjusted        — R²adj (Ezekiel 1930 formula; corrects for n and p)
    """
    from sklearn.linear_model import LinearRegression
    from sklearn.preprocessing import StandardScaler

    common = Y.index.intersection(X.index)
    Y_sub  = Y.loc[common]
    X_sub  = X.loc[common]

    n, p = X_sub.shape
    q    = Y_sub.shape[1]

    Y_sc = StandardScaler().fit_transform(Y_sub.values)
    X_sc = StandardScaler().fit_transform(X_sub.values)

    lm = LinearRegression(fit_intercept=True)
    lm.fit(X_sc, Y_sc)
    Y_hat = lm.predict(X_sc)

    ss_total   = np.sum(Y_sc ** 2)
    ss_fitted  = np.sum(Y_hat ** 2)
    r2         = ss_fitted / ss_total if ss_total > 0 else 0.0
    r2_adj     = 1 - (1 - r2) * (n - 1) / (n - p - 1) if n > p + 1 else r2

    k  = min(n_components, p, q, n - 1)
    pca = PCA(n_components=k, random_state=42)
    rda_scores = pca.fit_transform(Y_hat)

    return {
        "rda_scores":         pd.DataFrame(
            rda_scores, index=common,
            columns=[f"RDA{i+1}" for i in range(k)]),
        "r2_total":           float(r2),
        "r2_adjusted":        float(r2_adj),
        "explained_variance": pca.explained_variance_ratio_,
        "n_samples":          n,
        "n_predictors":       p,
    }


# ============================================================================
# DIFFERENTIAL ABUNDANCE
# ============================================================================

def differential_abundance(df:       pd.DataFrame,
                            meta:     pd.DataFrame,
                            group_col: str,
                            g1: str,
                            g2: str,
                            fdr: float = FDR_THRESHOLD) -> pd.DataFrame:
    """
    Mann-Whitney U test for each feature between g1 and g2.
    Returns: feature, U, pval, qval, effect_size (rank-biserial r),
             mean_g1, mean_g2, log2FC, median_g1, median_g2.
    """
    idx1 = meta.index[meta[group_col] == g1].intersection(df.index)
    idx2 = meta.index[meta[group_col] == g2].intersection(df.index)

    results = []
    for feat in df.columns:
        x1 = df.loc[idx1, feat].dropna().values
        x2 = df.loc[idx2, feat].dropna().values
        if len(x1) < 3 or len(x2) < 3:
            continue
        stat, pval = mannwhitneyu(x1, x2, alternative="two-sided")
        n1, n2     = len(x1), len(x2)
        effect     = 1 - (2 * stat) / (n1 * n2)  # rank-biserial correlation
        results.append({
            "feature":    feat,
            "U":          float(stat),
            "pval":       float(pval),
            "effect_size": float(effect),
            "mean_g1":    float(x1.mean()),
            "mean_g2":    float(x2.mean()),
            "median_g1":  float(np.median(x1)),
            "median_g2":  float(np.median(x2)),
            "log2FC":     float(np.log2((x2.mean() + 1e-9) / (x1.mean() + 1e-9))),
        })

    if not results:
        return pd.DataFrame(columns=["feature", "U", "pval", "effect_size",
                                      "mean_g1", "mean_g2", "median_g1", "median_g2",
                                      "log2FC", "qval"])
    res = pd.DataFrame(results)
    _, qvals, _, _ = multipletests(res["pval"], method="fdr_bh")
    res["qval"] = qvals
    return res.sort_values("qval").reset_index(drop=True)


# ============================================================================
# SPEARMAN CORRELATIONS  (vectorised; generalised — no pathway filter)
# ============================================================================

def spearman_correlation_matrix(species_df:    pd.DataFrame,
                                 metabolite_df: pd.DataFrame,
                                 fdr: float   = FDR_THRESHOLD,
                                 min_r: float  = MIN_CORR) -> pd.DataFrame:
    """
    Compute all species × metabolite Spearman correlations using vectorised
    rank-based approach (~100× faster than nested loops).
    Returns tidy DataFrame: species, metabolite, rho, pval, qval.
    Pre-filters by |rho| >= min_r before BH correction.
    """
    from scipy.stats import rankdata, t as t_dist

    common     = species_df.index.intersection(metabolite_df.index)
    sp         = species_df.loc[common].values.astype(float)
    mt         = metabolite_df.loc[common].values.astype(float)
    n          = sp.shape[0]
    spc_names  = list(species_df.columns)
    mtb_names  = list(metabolite_df.columns)

    def rank_cols(mat):
        return np.apply_along_axis(
            lambda col: rankdata(col, method="average"), 0, mat)

    sp_r  = rank_cols(sp)
    mt_r  = rank_cols(mt)
    sp_rc = sp_r - sp_r.mean(axis=0, keepdims=True)
    mt_rc = mt_r - mt_r.mean(axis=0, keepdims=True)

    sp_std  = np.sqrt((sp_rc ** 2).sum(axis=0, keepdims=True))
    mt_std  = np.sqrt((mt_rc ** 2).sum(axis=0, keepdims=True))
    rho_mat = (sp_rc.T @ mt_rc) / (sp_std.T * mt_std)

    rho_safe = np.clip(rho_mat, -0.9999, 0.9999)  # clip before squaring to avoid boundary discontinuity
    t_mat = rho_safe * np.sqrt((n - 2) / (1 - rho_safe ** 2))
    p_mat = 2 * t_dist.sf(np.abs(t_mat), df=n - 2)

    mask = np.abs(rho_mat) >= min_r
    rows_spc, rows_mtb = np.where(mask)

    if len(rows_spc) == 0:
        return pd.DataFrame(columns=["species", "metabolite", "rho", "pval", "qval"])

    df_out = pd.DataFrame({
        "species":    [spc_names[i] for i in rows_spc],
        "metabolite": [mtb_names[j] for j in rows_mtb],
        "rho":        rho_mat[rows_spc, rows_mtb],
        "pval":       p_mat[rows_spc, rows_mtb],
    })
    _, qvals, _, _ = multipletests(df_out["pval"], method="fdr_bh")
    df_out["qval"] = qvals
    df_out = df_out[df_out["qval"] <= fdr]
    return df_out.sort_values("qval").reset_index(drop=True)


# ============================================================================
# PARTIAL CORRELATION (confounder adjustment via OLS residuals)
# ============================================================================

def partial_corr_residuals(x: np.ndarray,
                            y: np.ndarray,
                            covariates: np.ndarray) -> tuple[float, float]:
    """
    Partial Spearman correlation between x and y after regressing out covariates.
    Returns (rho_partial, pval).
    """
    from sklearn.linear_model import LinearRegression
    mask = ~(np.isnan(x) | np.isnan(y) | np.any(np.isnan(covariates), axis=1))
    if mask.sum() < 10:
        return np.nan, np.nan
    lm = LinearRegression()
    lm.fit(covariates[mask], x[mask])
    rx = x[mask] - lm.predict(covariates[mask])
    lm.fit(covariates[mask], y[mask])
    ry = y[mask] - lm.predict(covariates[mask])
    return spearmanr(rx, ry)


# ============================================================================
# BOOTSTRAP MEDIATION (ACME; percentile bootstrap CI)
# ============================================================================

def bootstrap_mediation(
    x: np.ndarray,
    m: np.ndarray,
    y: np.ndarray,
    n_boot: int = 1000,
    ci: float   = 0.95,
    random_state: int = 42,
) -> dict:
    """
    Non-parametric percentile bootstrap for ACME (Average Causal Mediation Effect).

    Model:
        a  path : X → M  (species → metabolite)
        b  path : M → Y controlling for X
        c' path : direct X → Y controlling for M
        ACME    = a × b

    Returns dict with: acme, acme_ci_lo, acme_ci_hi, a, b, c_direct,
                       c_total, prop_mediated, p_value, n_boot, ci.
    """
    from sklearn.linear_model import LinearRegression

    rng  = np.random.default_rng(random_state)
    mask = ~(np.isnan(x) | np.isnan(m) | np.isnan(y))
    x, m, y = x[mask], m[mask], y[mask]
    if len(x) < 15:
        return {k: np.nan for k in
                ["acme", "acme_ci_lo", "acme_ci_hi", "a", "b",
                 "c_direct", "c_total", "prop_mediated", "p_value"]}

    def _acme(xi, mi, yi):
        a_  = LinearRegression().fit(xi.reshape(-1, 1), mi).coef_[0]
        bc  = LinearRegression().fit(np.column_stack([xi, mi]), yi)
        b_  = bc.coef_[1]
        cd  = bc.coef_[0]
        ct  = LinearRegression().fit(xi.reshape(-1, 1), yi).coef_[0]
        return a_ * b_, a_, b_, cd, ct

    acme_obs, a_obs, b_obs, c_dir, c_tot = _acme(x, m, y)

    boot_acme = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, len(x), size=len(x))
        boot_acme[i] = _acme(x[idx], m[idx], y[idx])[0]

    alpha = (1 - ci) / 2
    ci_lo = float(np.nanpercentile(boot_acme, alpha * 100))
    ci_hi = float(np.nanpercentile(boot_acme, (1 - alpha) * 100))

    # Permute the outcome (y) to break the mediation chain under the null.
    # Permuting x would test a different null ("no x→m→y path") which is weaker.
    perm_acme = np.empty(500)
    for i in range(500):
        idx = rng.permutation(len(y))
        perm_acme[i] = _acme(x, m, y[idx])[0]
    p_val = float((np.abs(perm_acme) >= abs(acme_obs)).mean())

    prop_med = acme_obs / c_tot if abs(c_tot) > 1e-9 else np.nan

    return {
        "acme":          float(acme_obs),
        "acme_ci_lo":    ci_lo,
        "acme_ci_hi":    ci_hi,
        "a":             float(a_obs),
        "b":             float(b_obs),
        "c_direct":      float(c_dir),
        "c_total":       float(c_tot),
        "prop_mediated": float(prop_med) if not np.isnan(prop_med) else np.nan,
        "p_value":       p_val,
        "n_boot":        n_boot,
        "ci":            ci,
    }


# ============================================================================
# SPECIES–SPECIES CO-ABUNDANCE CORRELATION
# ============================================================================

def species_coabundance_matrix(
    species_df: pd.DataFrame,
    fdr:   float = FDR_THRESHOLD,
    min_r: float = 0.20,
) -> pd.DataFrame:
    """
    Vectorised species × species Spearman co-abundance (upper triangle).
    Returns tidy DataFrame: species_a, species_b, rho, pval, qval.
    """
    from scipy.stats import rankdata, t as t_dist

    sp    = species_df.values.astype(float)
    n     = sp.shape[0]
    names = list(species_df.columns)

    def rank_cols(mat):
        return np.apply_along_axis(
            lambda col: rankdata(col, method="average"), 0, mat)

    sp_r   = rank_cols(sp)
    sp_rc  = sp_r - sp_r.mean(axis=0, keepdims=True)
    sp_std = np.sqrt((sp_rc ** 2).sum(axis=0, keepdims=True))
    rho_mat = (sp_rc.T @ sp_rc) / (sp_std.T * sp_std)

    rho_safe = np.clip(rho_mat, -0.9999, 0.9999)
    t_mat = rho_safe * np.sqrt((n - 2) / (1 - rho_safe ** 2))
    p_mat = 2 * t_dist.sf(np.abs(t_mat), df=n - 2)

    mask = np.triu(np.abs(rho_mat) >= min_r, k=1)
    rows_a, rows_b = np.where(mask)

    if len(rows_a) == 0:
        return pd.DataFrame(columns=["species_a", "species_b", "rho", "pval", "qval"])

    df_out = pd.DataFrame({
        "species_a": [names[i] for i in rows_a],
        "species_b": [names[j] for j in rows_b],
        "rho":       rho_mat[rows_a, rows_b],
        "pval":      p_mat[rows_a, rows_b],
    })
    _, qvals, _, _ = multipletests(df_out["pval"], method="fdr_bh")
    df_out["qval"] = qvals
    df_out = df_out[df_out["qval"] <= fdr]
    return df_out.sort_values("rho", key=abs, ascending=False).reset_index(drop=True)


# ============================================================================
# PICKLE I/O
# ============================================================================

def save_pickle(obj: object, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "wb") as f:
        pickle.dump(obj, f)
    print(f"Saved: {path}")


def load_pickle(path: Path) -> object:
    with open(path, "rb") as f:
        obj = pickle.load(f)
    print(f"Loaded: {path}")
    return obj


def load_multi_cohort_lodo(
    pkl: dict,
    cohorts: list | None = None,
    top_species: int = 300,
) -> tuple:
    """
    Concatenate CLR species matrices across LODO cohorts for cross-study CV.

    Parameters
    ----------
    pkl        : preprocessed_data.pkl dictionary
    cohorts    : list of dataset names; defaults to DATASETS_LODO ∩ pkl keys
    top_species: number of top-variance species to select per cohort before
                 taking the intersection across cohorts

    Returns
    -------
    X      : pd.DataFrame — concatenated species CLR (samples × common_species)
    y      : pd.Series    — Stage.3Group labels (index aligned to X)
    groups : pd.Series    — cohort label per sample (for LeaveOneGroupOut)

    Notes
    -----
    - IBD_COHORTS are automatically excluded.
    - Species are selected as the intersection of the top-`top_species`
      variance columns across cohorts; ensures all cohorts contribute equally.
    - Samples with NaN Stage.3Group are dropped.
    """
    if cohorts is None:
        cohorts = [c for c in DATASETS_LODO if c in pkl.get("species_clr", {})]

    spe_parts: list[tuple[str, pd.DataFrame, pd.DataFrame]] = []
    grp_parts: list[pd.Series] = []

    for ds in cohorts:
        if ds in IBD_COHORTS:
            continue
        if ds not in pkl.get("species_clr", {}):
            print(f"  WARNING: {ds} not in pkl — skipping (add TSV files to Data/ to activate)")
            continue
        spe  = pkl["species_clr"][ds]
        meta = pkl["harmonized_meta"][ds]
        valid = meta["Stage.3Group"].notna()
        if valid.sum() == 0:
            print(f"  WARNING: {ds} has no valid Stage.3Group labels — skipping")
            continue
        spe_parts.append((ds, spe.loc[valid], meta.loc[valid]))
        grp_parts.append(pd.Series(ds, index=meta.loc[valid].index, name="cohort"))

    if not spe_parts:
        raise ValueError(
            "No LODO cohorts found in pkl. "
            f"Checked: {cohorts}. "
            "Download Kong 2023 / Sinha-Vogtmann / Borenstein TSV files to Data/ to enable LODO."
        )

    # Common species: intersection of per-cohort top-variance columns
    common_sp = None
    for _, spe, _ in spe_parts:
        top = spe.var(axis=0).nlargest(top_species).index
        common_sp = top if common_sp is None else common_sp.intersection(top)

    if len(common_sp) == 0:
        raise ValueError(
            f"No common species found across cohorts after selecting top-{top_species} "
            "by variance per cohort. Consider increasing top_species."
        )

    X      = pd.concat([spe[common_sp]               for _, spe,  _    in spe_parts], axis=0)
    y      = pd.concat([meta["Stage.3Group"]          for _, _,    meta in spe_parts], axis=0)
    groups = pd.concat(grp_parts, axis=0)

    print(f"LODO matrix: {X.shape[0]} samples × {X.shape[1]} species "
          f"from {len(spe_parts)} cohorts")
    return X, y, groups


# ============================================================================
# PLOT HELPERS
# ============================================================================

PALETTE_STAGE6 = {
    "Healthy":      "#4CAF50",
    "Stage_0":      "#8BC34A",
    "HS":           "#FFC107",
    "Stage_I_II":   "#FF9800",
    "Stage_III_IV": "#F44336",
    "MP":           "#9C27B0",
}

PALETTE_3GROUP = {
    "Healthy":      "#4CAF50",
    "Early_CRC":    "#FFC107",
    "Advanced_CRC": "#F44336",
}


def savefig(fig: plt.Figure, subdir: str, filename: str, dpi: int = 150) -> None:
    out = FIG_DIR / subdir
    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / filename, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved figure: {out / filename}")
