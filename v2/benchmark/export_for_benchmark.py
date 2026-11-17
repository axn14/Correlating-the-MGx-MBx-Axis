"""Benchmark step 1/3 — export Yachida matrices in MIMOSA2 and MelonnPan input formats.

Reads the v2 preprocessed pickle and writes four TSVs under benchmark/inputs/:
  MelonnPan (samples x features, relative abundance):
    mp_metag_samples_x_species.tsv      metab predictor input (microbiome)
    mp_metab_samples_x_metabolites.tsv  metab response input
  MIMOSA2 (features x samples; KEGG compound IDs for mets):
    mimosa2_species_x_samples.tsv       file1 (taxonomic abundances)
    mimosa2_mets_x_samples.tsv          file2 (metabolites, KEGG IDs)

Both tools compare against the SAME cohort the framework used, so the comparison is
like-for-like. Run on the server venv:  python export_for_benchmark.py
"""
import pickle
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # v2/ on path
import numpy as np
import pandas as pd
import utils

RESULTS_DIR = utils.RESULTS_DIR
INTER       = RESULTS_DIR / 'intermediate'
OUT         = Path(__file__).resolve().parent / 'inputs'
OUT.mkdir(parents=True, exist_ok=True)
YACHIDA     = utils.DATASET_DISCOVERY

with open(INTER / 'preprocessed_data_v2.pkl', 'rb') as f:
    pdata = pickle.load(f)

def _get(key):
    d = pdata.get(key, {})
    return d.get(YACHIDA) if isinstance(d, dict) else None

# Prefer pre-transform RELATIVE abundances (both tools want compositions, not CLR/log).
spc = _get('filtered_species')
mtb = _get('filtered_mtb')
if spc is None:  # fall back; note this is the transformed matrix -> VERIFY scale
    spc = _get('transformed_species'); print('WARNING: using transformed_species (no filtered_species)  # >>> VERIFY scale')
if mtb is None:
    mtb = _get('transformed_mtb');     print('WARNING: using transformed_mtb (no filtered_mtb)  # >>> VERIFY scale')
meta = _get('harmonized_meta')
mtb_kegg_map = (pdata.get('mtb_kegg_map', {}) or {}).get(YACHIDA, {})   # {column: kegg}

for df in (spc, mtb, meta):
    if 'Sample' not in df.columns:
        df.reset_index(inplace=True); df.rename(columns={df.columns[0]: 'Sample'}, inplace=True)

# Primary samples only (match the framework's scope).
if hasattr(utils, 'primary_mask'):
    keep = meta.set_index('Sample').index[utils.primary_mask(meta.set_index('Sample')).values]
    meta = meta[meta['Sample'].isin(keep)]
common = sorted(set(spc['Sample']) & set(mtb['Sample']) & set(meta['Sample']))
spc = spc.set_index('Sample').loc[common]
mtb = mtb.set_index('Sample').loc[common]
print(f'Samples: {len(common)} | species: {spc.shape[1]} | metabolites: {mtb.shape[1]}')

# Renormalise species to relative abundance (sum=1 per sample).
spc = spc.div(spc.sum(axis=1), axis=0).fillna(0.0)

# Metabolite columns -> KEGG compound IDs where known (MIMOSA2 metType = KEGG Compound IDs).
def to_kegg(col):
    if col in mtb_kegg_map and str(mtb_kegg_map[col]).startswith('C'):
        return mtb_kegg_map[col]
    if str(col)[:1] == 'C' and str(col)[1:6].isdigit():
        return str(col).split('_')[0]
    return col
mtb_kegg = mtb.rename(columns={c: to_kegg(c) for c in mtb.columns})
# collapse duplicate KEGG columns (name + id mapping to same id) by mean
mtb_kegg = mtb_kegg.T.groupby(level=0).mean().T

# ---- MelonnPan: samples x features, sample id as first column ----
spc.reset_index().rename(columns={'index': 'Sample'}).to_csv(
    OUT / 'mp_metag_samples_x_species.tsv', sep='\t', index=False)
mtb_kegg.reset_index().rename(columns={'index': 'Sample'}).to_csv(
    OUT / 'mp_metab_samples_x_metabolites.tsv', sep='\t', index=False)

# ---- MIMOSA2: features x samples, feature id as first column ----
spc.T.reset_index().rename(columns={'index': 'Species'}).to_csv(
    OUT / 'mimosa2_species_x_samples.tsv', sep='\t', index=False)
mtb_kegg.T.reset_index().rename(columns={'index': 'compound'}).to_csv(
    OUT / 'mimosa2_mets_x_samples.tsv', sep='\t', index=False)

# KEGG map sidecar for the comparator (compound id -> display name).
pd.DataFrame({'KEGG': [to_kegg(c) for c in mtb.columns],
              'Name': list(mtb.columns)}).drop_duplicates('KEGG').to_csv(
    OUT / 'metabolite_kegg_names.tsv', sep='\t', index=False)

print('Wrote:')
for p in sorted(OUT.glob('*.tsv')):
    print(f'  {p.name}')
print('\nPolyamine KEGG IDs to track:', list(utils.POLYAMINE_KEGG.items()))
