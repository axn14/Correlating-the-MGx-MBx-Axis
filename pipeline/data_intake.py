"""
data_intake.py — Format detection and canonical conversion for user-submitted data.

Supported inputs:
  Species/microbiome:
    - MetaPhlAn3/4 profile  (clade names contain k__|p__|s__ patterns)
    - Generic OTU/ASV table (plain taxon names, rows=taxa, cols=samples)
    - Auto-transposes if samples appear to be rows

  Metabolomics:
    - KEGG IDs in first column  (C12345 pattern)
    - HMDB IDs (HMDB00001)  → mapped to KEGG via reference/hmdb_to_kegg.csv
    - Compound names         → mapped via reference/compound_name_map.csv
    - Generic feature names  → pipeline runs without KEGG annotation

  Metadata:
    - Any TSV/CSV with a sample_id / Sample / SampleID column
"""

from __future__ import annotations

import io
import re
from pathlib import Path

import numpy as np
import pandas as pd

_REF_DIR = Path(__file__).parent / "reference"

_METAPHLAN_RE = re.compile(r"k__\w+")
_KEGG_RE      = re.compile(r"^C\d{5}$")
_HMDB_RE      = re.compile(r"^HMDB\d+$", re.IGNORECASE)
_SAMPLE_ID_CANDIDATES = {"sample_id", "sampleid", "sample", "#sampleid", "id"}


# ── Internal reference maps ───────────────────────────────────────────────────

def _load_hmdb_map() -> dict[str, str]:
    p = _REF_DIR / "hmdb_to_kegg.csv"
    if not p.exists():
        return {}
    df = pd.read_csv(p, dtype=str)
    return dict(zip(df.iloc[:, 0].str.upper(), df.iloc[:, 1]))


def _load_name_map() -> dict[str, str]:
    p = _REF_DIR / "compound_name_map.csv"
    if not p.exists():
        return {}
    df = pd.read_csv(p, dtype=str)
    return dict(zip(df.iloc[:, 0].str.lower(), df.iloc[:, 1]))


_HMDB_MAP: dict[str, str] = {}
_NAME_MAP: dict[str, str] = {}


def _get_hmdb_map():
    global _HMDB_MAP
    if not _HMDB_MAP:
        _HMDB_MAP = _load_hmdb_map()
    return _HMDB_MAP


def _get_name_map():
    global _NAME_MAP
    if not _NAME_MAP:
        _NAME_MAP = _load_name_map()
    return _NAME_MAP


# ── File reading ──────────────────────────────────────────────────────────────

def _read_table(file_obj) -> pd.DataFrame:
    """Read TSV or CSV from a file-like object or bytes."""
    if isinstance(file_obj, (bytes, bytearray)):
        file_obj = io.BytesIO(file_obj)
    sample = file_obj.read(1024)
    file_obj.seek(0)
    sep = "\t" if b"\t" in sample else ","
    return pd.read_csv(file_obj, sep=sep, index_col=0)


# ── Species detection ─────────────────────────────────────────────────────────

def _is_metaphlan(index: pd.Index) -> bool:
    return any(_METAPHLAN_RE.search(str(v)) for v in index[:10])


def _filter_metaphlan_species(df: pd.DataFrame) -> pd.DataFrame:
    """Keep only species-level rows (contain |s__) and strip to short name."""
    mask = df.index.str.contains(r"\|s__", na=False)
    df = df[mask].copy()
    # Extract the s__ part as short name
    df.index = df.index.str.extract(r"\|s__([\w]+)$", expand=False).fillna(df.index)
    return df


def detect_and_load_species(file_obj) -> tuple[pd.DataFrame, str]:
    """
    Returns (samples × species DataFrame, format_description).
    Rows are samples, columns are species. Values are relative abundances or counts.
    """
    df = _read_table(file_obj)

    # Auto-transpose: if more rows than columns, species are probably rows already
    if df.shape[0] > df.shape[1]:
        # rows=taxa, cols=samples → transpose to samples×taxa
        df = df.T
    # After transpose: rows=samples, cols=taxa

    # If columns look like MetaPhlAn taxa, filter to species level
    if _is_metaphlan(df.columns):
        # Filter to s__ level rows from the transposed frame
        df = df.T  # back to taxa×samples
        df = _filter_metaphlan_species(df)
        df = df.T  # samples×species
        fmt = "MetaPhlAn3/4 profile (species-level filtered)"
    else:
        fmt = "Generic OTU/ASV table"

    df = df.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    df.columns = df.columns.astype(str)
    df.index   = df.index.astype(str)
    return df, fmt


# ── Metabolomics detection ────────────────────────────────────────────────────

def _map_identifier(identifier: str) -> str | None:
    """Map HMDB ID or compound name to KEGG ID. Returns None if unknown."""
    if _KEGG_RE.match(identifier):
        return identifier
    if _HMDB_RE.match(identifier):
        return _get_hmdb_map().get(identifier.upper())
    return _get_name_map().get(identifier.lower())


def detect_and_load_metabolomics(file_obj) -> tuple[pd.DataFrame, str, dict[str, str]]:
    """
    Returns (samples × metabolites DataFrame, format_description, kegg_map).
    kegg_map: {column_name: KEGG_ID} — empty if no KEGG mapping possible.
    """
    df = _read_table(file_obj)

    if df.shape[0] > df.shape[1]:
        df = df.T  # to samples×metabolites

    df = df.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    df.columns = df.columns.astype(str)
    df.index   = df.index.astype(str)

    # Detect identifier type from column names
    sample_cols = list(df.columns[:5])
    kegg_map: dict[str, str] = {}
    fmt = "Generic feature names (no KEGG annotation)"

    if all(_KEGG_RE.match(c) for c in sample_cols if c):
        fmt = "KEGG IDs (direct)"
        kegg_map = {c: c for c in df.columns}
    else:
        # Attempt mapping
        for col in df.columns:
            k = _map_identifier(col)
            if k:
                kegg_map[col] = k
        if len(kegg_map) > 0.3 * len(df.columns):
            fmt = f"Mapped {len(kegg_map)}/{len(df.columns)} features to KEGG IDs"
        elif kegg_map:
            fmt = f"Partial KEGG mapping ({len(kegg_map)} features)"

    return df, fmt, kegg_map


# ── Metadata detection ────────────────────────────────────────────────────────

def detect_and_load_metadata(file_obj) -> pd.DataFrame:
    """
    Returns DataFrame indexed by sample_id.
    Attempts to auto-detect the sample_id column.
    """
    df = _read_table(file_obj)

    # Check if index already looks like sample IDs
    if df.index.name and df.index.name.lower().replace(" ", "_") in _SAMPLE_ID_CANDIDATES:
        df.index = df.index.astype(str)
        return df

    # Try to find a sample_id column
    for col in df.columns:
        if col.lower().replace(" ", "_").replace("-", "_") in _SAMPLE_ID_CANDIDATES:
            df = df.set_index(col)
            df.index = df.index.astype(str)
            return df

    # Fall back: use the existing index
    df.index = df.index.astype(str)
    return df


# ── Overlap validation ────────────────────────────────────────────────────────

def validate_overlap(
    species_df: pd.DataFrame,
    metabolomics_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
) -> dict:
    """
    Check sample overlap across three DataFrames.
    Returns dict with n_common_samples and list of warnings.
    """
    spe_samples = set(species_df.index)
    mtb_samples = set(metabolomics_df.index)
    meta_samples = set(metadata_df.index)

    common = spe_samples & mtb_samples & meta_samples
    warnings_list = []

    only_spe  = spe_samples  - common
    only_mtb  = mtb_samples  - common
    only_meta = meta_samples - common

    if only_spe:
        warnings_list.append(f"{len(only_spe)} samples only in species table (will be dropped)")
    if only_mtb:
        warnings_list.append(f"{len(only_mtb)} samples only in metabolomics table (will be dropped)")
    if only_meta:
        warnings_list.append(f"{len(only_meta)} samples only in metadata (will be dropped)")
    if len(common) < 10:
        warnings_list.append(f"WARNING: only {len(common)} shared samples — results may be unreliable")

    return {
        "n_common_samples": len(common),
        "common_samples":   sorted(common),
        "warnings":         warnings_list,
    }


def align_to_common_samples(
    species_df: pd.DataFrame,
    metabolomics_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    common_samples: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Subset all three DataFrames to the common sample set, in the same order."""
    return (
        species_df.loc[common_samples],
        metabolomics_df.loc[common_samples],
        metadata_df.loc[common_samples],
    )
