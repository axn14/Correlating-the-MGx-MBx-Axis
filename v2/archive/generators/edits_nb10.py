"""NB10 edits: CREATE 10_micom_flux.ipynb from scratch (per CONTRACTS.md section 8).

MICOM community flux balance on AGORA2-GTDB models, producing the single
mechanistic-evidence input NB08 consumes:

    tables/nb10_micom_species_flux_pairs_v2.csv  -> columns: Species, Polyamine, Rho, QValue

NB08 Part B renames Rho->MICOM_flux_rho, QValue->MICOM_flux_q, merges on
(Species, Polyamine), and sets Evidence_mechanistic = (q<0.05) & (rho>0).

This is a *draft staged for a run*, not yet executed. First-run VERIFY points are
marked  # >>> VERIFY  — the name-overlap diagnostic and the smoke gate both run
before the multi-hour batch so a mismatch fails in seconds, not hours.

Run via run_micom.sh in its own tmux session, AFTER the main pipeline. Then re-run
NB08 + NB09 so the mechanistic column fills.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from nbedit import Notebook

cells = []

# ------------------------------------------------------------------ 00 title
cells.append(("markdown", "nb10-00-title", r'''# 10 — MICOM Community Flux (mechanistic evidence)

**Decoupled run (CONTRACTS section 8).** Not part of `run_pipeline.sh`; launched by
`run_micom.sh` in its own tmux session because a community solve is minutes per sample
over hundreds of samples. When it finishes, re-run NB08 + NB09 so the mechanistic
evidence column fills; until then `MICOM_status == 'pending'` and the attainable
maximum is 4 of 5 sources.

## One job

> **`tables/nb10_micom_species_flux_pairs_v2.csv`** — one row per (species, polyamine),
> with the Spearman correlation between that species' MICOM-predicted **export flux** of
> the polyamine (per sample) and the **measured** metabolite (per sample).
> Columns: `Species, Polyamine, Rho, QValue` (+ provenance). NB08 reads exactly these.

The mechanistic claim is deliberately narrow: a positive, significant correlation means
*the species' modelled capacity to secrete the compound tracks the compound's measured
level across patients* — community metabolic modelling and the metabolome agree on who
makes it. It is still observational (abundance-weighted), which is why it is one source
of five, never a proof.

## Fixes vs the earlier MICOM attempts (publication/NB09, Git notebook)

| Earlier problem | Here |
|---|---|
| AGORA103, ~42% abundance coverage | **AGORA2-GTDB** (`agora2_gtdb_species.qza`), GTDB-native names |
| GTDB strains collapsed to NCBI binomials and summed | **no collapsing** — a taxon matches exactly or is reported unmatched |
| tradeoff hardcoded at 0.5 | **sweep** on 10 stratified pilot samples; chosen value + sweep saved |
| `strategy="none"` -> non-unique fluxes | **`strategy="pFBA"`** |
| growth rates discarded; fluxes >=80 deleted | growth rates kept; extreme fluxes **flagged, never deleted** |
| statistics on positive fluxes only (true by construction) | all samples, zeros included |
| HS->Early, MP->Advanced grouping | **v2.1 primary grouping** (MP excluded, HS quarantined) |

## Inputs / outputs

- **Inputs:** `Data/micom/agora2_gtdb_species.qza`, `Data/micom/western_diet_gut.csv`,
  `intermediate/preprocessed_data_v2.pkl` (full filtered species, raw relative abundance),
  and — when present — `Data/micom/models_agora2/` (**347 pre-built AGORA2 community
  models + `manifest.csv`, covering all 277 primary Yachida samples**). That set is used
  in place as a read-only input and the build step is skipped, saving hours of solver
  time and ~36 GB. It is never deleted by the cleanup cell. Move it aside to force a
  rebuild.
- **Checkpoints:** `intermediate/micom_v2/` (manifest, per-batch exchange pickles).
- **Tables:** `nb10_micom_sample_coverage_v2.csv`, `nb10_micom_community_flux_v2.csv`,
  `nb10_micom_species_flux_pairs_v2.csv`.'''))

# ------------------------------------------------------------------ 01 imports / paths
cells.append(("code", "nb10-01-imports", r'''import sys, os, pickle, warnings, re
sys.path.insert(0, '.')
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

import utils

warnings.filterwarnings('ignore')

try:
    import micom
    from micom.workflows import build, grow
    from micom.qiime_formats import load_qiime_model_db
    print(f'micom {micom.__version__} loaded')
except ImportError as e:
    raise ImportError(
        'micom is required for NB10. Install into the pipeline venv:\n'
        '    pip install "micom>=0.33" osqp highspy\n'
        f'(import error: {e})')

RESULTS_DIR = utils.RESULTS_DIR
DATA_DIR    = utils.DATA_DIR
TABLES      = RESULTS_DIR / 'tables'
FIGURES     = RESULTS_DIR / 'figures' / 'micom'
INTER       = RESULTS_DIR / 'intermediate'
MICOM_WORK  = INTER / 'micom_v2'                 # NB10's own checkpoint dir
for d in (TABLES, FIGURES, MICOM_WORK):
    d.mkdir(parents=True, exist_ok=True)

MICOM_DIR  = DATA_DIR / 'micom'
AGORA_DB   = MICOM_DIR / 'agora2_gtdb_species.qza'     # GTDB-native, section 8
DIET_FILE  = MICOM_DIR / 'western_diet_gut.csv'

# Community models. A COMPLETE AGORA2 set already exists at Data/micom/models_agora2:
# 347 pickles + manifest.csv, covering every Yachida sample including all 277 primary,
# built by an earlier run. Rebuilding costs hours of solver time and ~36 GB of disk, so
# when that set is present it is used IN PLACE as a read-only INPUT rather than copied
# or regenerated. MODELS_ARE_PREBUILT then guards the cleanup step at the end of this
# notebook, which must never unlink an input.
PREBUILT_MODELS   = MICOM_DIR / 'models_agora2'
PREBUILT_MANIFEST = PREBUILT_MODELS / 'manifest.csv'
MODELS_ARE_PREBUILT = (PREBUILT_MANIFEST.exists()
                       and any(PREBUILT_MODELS.glob('*.pickle')))

if MODELS_ARE_PREBUILT:
    MODELS_DIR = PREBUILT_MODELS
    _n_prebuilt = len(list(PREBUILT_MODELS.glob('*.pickle')))
    print(f'Using {_n_prebuilt} pre-built AGORA2 community models (READ-ONLY input): {MODELS_DIR}')
else:
    MODELS_DIR = MICOM_WORK / 'models_agora2'   # built here, deleted after fluxes saved
    MODELS_DIR.mkdir(parents=True, exist_ok=True)

YACHIDA = utils.DATASET_DISCOVERY

# Per-sample resource budget (section 8). Threading kept to the 2-core cap.
N_THREADS       = int(os.environ.get('MICOM_THREADS', '2'))
ABUND_CUTOFF    = 1e-3        # drop taxa below 0.1% rel. abundance from the community
MIN_COVERAGE    = 0.50        # samples below this flagged + excluded in the sensitivity pass
EXTREME_FLUX    = 80.0        # flag (never delete) |flux| above this
SEED            = 42

missing = [str(p) for p in (AGORA_DB, DIET_FILE) if not p.exists()]
if missing:
    raise FileNotFoundError(
        'Required MICOM input(s) missing:\n  ' + '\n  '.join(missing) +
        '\nSee CONTRACTS.md section 8 for provenance.')
print(f'AGORA2 DB : {AGORA_DB.name} ({AGORA_DB.stat().st_size/1e6:.0f} MB)')
print(f'Diet      : {DIET_FILE.name}')
print(f'Work dir  : {MICOM_WORK}')
print(f'Threads   : {N_THREADS}')'''))

# ------------------------------------------------------------------ 02 polyamine -> VMH exchange map
cells.append(("code", "nb10-02-vmh-map", r'''# ---------------------------------------------------------------------------
# Polyamine KEGG -> AGORA2/VMH metabolite id -> exchange-reaction match.
# >>> VERIFY on first run: the diagnostic below prints which of these the loaded
#     models actually expose as exchanges. An unmatched compound yields NO
#     mechanistic evidence for its pairs (MICOM_status stays 'not_modelled' for
#     them) rather than a wrong number.
# ---------------------------------------------------------------------------
# VMH abbreviations (https://www.vmh.life). Several synonyms per compound are
# listed; matching is by base id after stripping any compartment suffix.
POLYAMINE_VMH = {
    'C00134': ('Putrescine', ['ptrc']),                 # putrescine
    'C00315': ('Spermidine', ['spmd']),                 # spermidine
    'C00750': ('Spermine',   ['sprm']),                 # spermine
    'C00179': ('Agmatine',   ['agm']),                  # agmatine
    'C01672': ('Cadaverine', ['15dap', 'cadav']),       # 1,5-diaminopentane
    'C00077': ('Ornithine',  ['orn', 'orndp']),         # L-ornithine
    'C00062': ('Arginine',   ['arg__L', 'arg_L']),      # L-arginine
}
# Compound name -> KEGG (NB08 merges on the NAME in the 'Polyamine' column)
NAME_TO_KEGG = {name: k for k, (name, _ids) in POLYAMINE_VMH.items()}

_SUFFIX = re.compile(r'[\[\(_](?:e|m|lu|u|c)\]?\)?$', re.IGNORECASE)   # strip [e]/(e)/_m/_lu etc.

def _base_metabolite_id(met_id: str) -> str:
    """Normalise a model metabolite / exchange id to its bare VMH base id."""
    s = str(met_id)
    s = re.sub(r'^EX_', '', s)
    s = _SUFFIX.sub('', s)
    return s.strip().lower()

# Precompute the set of acceptable base ids per compound (lower-cased)
VMH_BASE = {k: {b.lower() for b in ids} for k, (_n, ids) in POLYAMINE_VMH.items()}

def classify_metabolite(met_id: str):
    """Return (KEGG, name) if this model metabolite is a tracked polyamine, else (None, None)."""
    base = _base_metabolite_id(met_id)
    for k, bases in VMH_BASE.items():
        if base in bases:
            return k, POLYAMINE_VMH[k][0]
    return None, None

print('Tracked polyamines (KEGG -> name -> VMH base ids):')
for k, (name, ids) in POLYAMINE_VMH.items():
    print(f'  {k} {name:11s} <- {ids}')'''))

# ------------------------------------------------------------------ 03 load abundances + measured metabolites
cells.append(("code", "nb10-03-load", r'''# ---------------------------------------------------------------------------
# Raw RELATIVE abundances (not CLR) for the FULL filtered species table, primary
# samples only. MICOM needs compositions, not log-ratios.
# ---------------------------------------------------------------------------
with open(INTER / 'preprocessed_data_v2.pkl', 'rb') as f:
    pdata = pickle.load(f)

# filtered_species = prevalence-filtered relative abundance, pre-CLR (CONTRACTS NB01).
filt = pdata['filtered_species'][YACHIDA].copy()
meta = pdata['harmonized_meta'][YACHIDA].copy()

# Measured metabolites: use the TRANSFORMED (log10) matrix; Spearman is rank-based
# so the monotone transform does not matter, and this matches NB02/NB08 scales.
mtb  = pdata['transformed_mtb'][YACHIDA].copy()
polyamine_columns = pdata.get('polyamine_columns', {}).get(YACHIDA, {})   # {kegg: column}
polyamine_proxies = pdata.get('polyamine_proxies', {}).get(YACHIDA, {})   # {kegg: column}
mtb_kegg_map      = pdata.get('mtb_kegg_map', {}).get(YACHIDA, {})        # {column: kegg}

def _sample_col(df):
    return 'Sample' if 'Sample' in df.columns else df.index.name or 'index'

for df in (filt, meta, mtb):
    if 'Sample' not in df.columns:
        df.reset_index(inplace=True)
        df.rename(columns={df.columns[0]: 'Sample'}, inplace=True)

# PRIMARY samples only (Healthy + Early_CRC + Advanced_CRC): HS/MP never reach MICOM.
prim = utils.primary_mask(meta.set_index('Sample')) if hasattr(utils, 'primary_mask') else None
if prim is not None:
    keep = meta.set_index('Sample').index[prim.values]
    meta = meta[meta['Sample'].isin(keep)].copy()
assert set(meta['Stage.Group'].dropna()) <= set(utils.STAGE_GROUP_ORDER), 'non-primary group reached NB10'

common = sorted(set(filt['Sample']) & set(meta['Sample']) & set(mtb['Sample']))
filt = filt.set_index('Sample').loc[common]
meta = meta.set_index('Sample').loc[common]
mtb  = mtb.set_index('Sample').loc[common]
print(f'Primary samples with species+metabolome: {len(common)}')
print(meta['Stage.Group'].value_counts().to_string())

# Renormalise each sample to a composition (sum=1).
species_cols = [c for c in filt.columns if c != 'Sample']
filt = filt[species_cols].div(filt[species_cols].sum(axis=1), axis=0).fillna(0.0)
print(f'Species in full filtered table: {filt.shape[1]}')

def resolve_metabolite_column(kegg_id):
    """Measured column for a polyamine, falling back to a pre-specified proxy."""
    for src in (polyamine_columns, polyamine_proxies):
        if kegg_id in src and src[kegg_id] in mtb.columns:
            return src[kegg_id], (src is polyamine_proxies)
    inv = {v: k for k, v in (mtb_kegg_map or {}).items()}
    if kegg_id in inv and inv[kegg_id] in mtb.columns:
        return inv[kegg_id], False
    for col in mtb.columns:
        if col.startswith(str(kegg_id) + '_'):
            return col, False
    return None, False'''))

# ------------------------------------------------------------------ 04 tax_long (GTDB-native, no renaming)
cells.append(("code", "nb10-04-taxlong", r'''# ---------------------------------------------------------------------------
# Build MICOM's taxonomy table. GTDB-native: a taxon matches AGORA2 exactly or is
# reported unmatched. We DO NOT strip _A suffixes or collapse strains (section 8).
# ---------------------------------------------------------------------------
def gtdb_binomial(col: str):
    """Species binomial from a column name that may be a full GTDB lineage or a bare name.
    Returns (genus, species_full) keeping GTDB suffixes, or (None, None) if unusable."""
    s = str(col)
    if ';s__' in s or s.startswith('s__'):
        m = re.search(r's__([^;]+)', s)
        name = m.group(1).strip() if m else ''
    else:
        name = s.strip()
    toks = name.split()
    if len(toks) < 2:
        return (None, None)
    genus = toks[0]
    return (genus, name)

long = (filt.reset_index().melt(id_vars='Sample', var_name='_raw', value_name='abundance'))
long = long[long['abundance'] > ABUND_CUTOFF].copy()
gs = long['_raw'].map(gtdb_binomial)
long['genus']   = gs.map(lambda t: t[0])
long['species'] = gs.map(lambda t: t[1])
long = long.dropna(subset=['genus', 'species'])
long = long.rename(columns={'Sample': 'sample_id'})
long['sample_id'] = long['sample_id'].astype(str)
# MICOM expects an 'id' per taxon; keep the exact species string so the taxon label
# that comes back out of grow() is identical to NB02/NB08's Species column.
long['id'] = long['species']
tax_long = long[['sample_id', 'id', 'genus', 'species', 'abundance']].copy()
print(f'tax_long: {tax_long["sample_id"].nunique()} samples, '
      f'{tax_long["species"].nunique()} unique species, {len(tax_long)} rows')'''))

# ------------------------------------------------------------------ 05 name-overlap diagnostic (fail fast)
cells.append(("code", "nb10-05-overlap", r'''# ---------------------------------------------------------------------------
# >>> VERIFY (seconds, before the hours-long build): how many of our species
# actually exist in the AGORA2-GTDB db? Zero overlap = a naming problem to fix
# BEFORE build(), not after.
# ---------------------------------------------------------------------------
from zipfile import ZipFile
with ZipFile(str(AGORA_DB)) as z:
    man_files = [f for f in z.namelist() if f.endswith('manifest.csv')]
    if not man_files:
        raise RuntimeError(f'No manifest.csv inside {AGORA_DB.name}')
    with z.open(man_files[0]) as f:
        agora_mdf = pd.read_csv(f)
print(f'AGORA2 db manifest: {len(agora_mdf)} models, columns={list(agora_mdf.columns)}')

# Match on the species string (GTDB-native). Try the likely columns.
agora_species = set()
for col in ('species', 'id', 'file'):
    if col in agora_mdf.columns:
        agora_species |= set(agora_mdf[col].astype(str).str.replace(r'\.(xml|json|mat)(\.gz)?$', '', regex=True))
ours = set(tax_long['species'])
hit  = ours & agora_species
print(f'Our species: {len(ours)} | AGORA2 species: {len(agora_species)} | matched: {len(hit)} '
      f'({100*len(hit)/max(1,len(ours)):.1f}% of ours)')
if len(hit) == 0:
    print('WARNING: zero overlap. Inspect the two namings below and adjust gtdb_binomial().')
    print('  ours  :', sorted(list(ours))[:5])
    print('  agora :', sorted(list(agora_species))[:5])
else:
    print('  matched examples:', sorted(list(hit))[:8])'''))

# ------------------------------------------------------------------ 06 build community models (checkpointed)
cells.append(("code", "nb10-06-build", r'''# ---------------------------------------------------------------------------
# Build per-sample community models from the pre-extracted AGORA2 dir (unzipping the
# .qza per sample would add ~1 min each). Checkpointed by manifest; existing model
# pickles are reused on a re-run.
# ---------------------------------------------------------------------------
manifest_path = MICOM_WORK / 'manifest_agora2_v2.csv'
AGORA_EXTRACTED = MICOM_DIR / 'agora2_extracted'

_n_needed = tax_long['sample_id'].nunique()

if manifest_path.exists() and pd.read_csv(manifest_path).shape[0] >= _n_needed:
    manifest = pd.read_csv(manifest_path)
    print(f'Reusing manifest ({len(manifest)} samples).')
elif MODELS_ARE_PREBUILT and pd.read_csv(PREBUILT_MANIFEST).shape[0] >= _n_needed:
    # Adopt the pre-built set's manifest instead of rebuilding. Its `file` column holds
    # a bare filename (e.g. "10021.pickle"), which micom resolves against
    # model_folder=MODELS_DIR, so no path rewriting is needed. Its coverage columns
    # (found_fraction / found_abundance_fraction) are the ones the coverage cell expects.
    manifest = pd.read_csv(PREBUILT_MANIFEST)
    _missing = [f for f in manifest['file'] if not (MODELS_DIR / str(f)).exists()]
    if _missing:
        raise FileNotFoundError(
            f'{len(_missing)} model pickle(s) named in {PREBUILT_MANIFEST.name} are absent '
            f'from {MODELS_DIR} (first: {_missing[:3]}). Refusing to run on a partial '
            'model set -- a silently incomplete community would bias every flux.')
    _have = set(manifest['sample_id'].astype(str))
    _want = set(tax_long['sample_id'].astype(str))
    if not _want <= _have:
        raise ValueError(
            f'{len(_want - _have)} sample(s) in this analysis have no pre-built model '
            f'(first: {sorted(_want - _have)[:3]}). Move or delete '
            f'{PREBUILT_MODELS} to force a rebuild.')
    manifest.to_csv(manifest_path, index=False)
    print(f'Adopted pre-built manifest ({len(manifest)} samples); build step skipped.')
else:
    if AGORA_EXTRACTED.exists() and any(AGORA_EXTRACTED.rglob('manifest.csv')):
        model_db = str(next(AGORA_EXTRACTED.rglob('manifest.csv')).parent)
        print(f'Reusing pre-extracted AGORA2 at {model_db}')
    else:
        AGORA_EXTRACTED.mkdir(parents=True, exist_ok=True)
        print('Pre-extracting AGORA2 .qza (one-time) ...')
        _m = load_qiime_model_db(str(AGORA_DB), str(AGORA_EXTRACTED))
        model_db = str(Path(_m['file'].iloc[0]).parent)
    print(f'Building community models (threads={N_THREADS}) ...')
    manifest = build(tax_long, out_folder=str(MODELS_DIR), model_db=model_db,
                     cutoff=ABUND_CUTOFF, threads=N_THREADS)
    manifest.to_csv(manifest_path, index=False)
print(f'Models available: {len(manifest)}')

# Coverage per sample -> sample_coverage table (section 8: flagged, excluded in sensitivity).
cov_col = next((c for c in ('found_fraction', 'found_abundance_fraction', 'coverage')
                if c in manifest.columns), None)
cov = manifest[['sample_id']].copy()
cov['coverage'] = manifest[cov_col] if cov_col else np.nan
cov['below_min_coverage'] = cov['coverage'] < MIN_COVERAGE if cov_col else False
cov = cov.merge(meta.reset_index()[['Sample', 'Stage.Group']],
                left_on='sample_id', right_on='Sample', how='left').drop(columns='Sample')
cov.to_csv(TABLES / 'nb10_micom_sample_coverage_v2.csv', index=False)
if cov_col:
    print(f'Median AGORA2 coverage: {cov["coverage"].median():.1%}; '
          f'{int(cov["below_min_coverage"].sum())} samples < {MIN_COVERAGE:.0%} (flagged).')'''))

# ------------------------------------------------------------------ 07 tradeoff sweep (pilot)
cells.append(("code", "nb10-07-tradeoff", r'''# ---------------------------------------------------------------------------
# Choose the cooperative-tradeoff by a sweep on 10 stratified pilot samples
# (section 8), instead of hardcoding 0.5. Value + full sweep are saved.
# Needs a QP-capable solver (osqp / HiGHS); community.py is stock micom.
# ---------------------------------------------------------------------------
tradeoff_path = MICOM_WORK / 'tradeoff_sweep_v2.csv'
medium = pd.read_csv(DIET_FILE)
if not {'reaction', 'flux'} <= set(medium.columns):
    raise ValueError(f'{DIET_FILE.name} needs reaction,flux columns; got {list(medium.columns)}')
# AGORA community medium reactions are lumen exchanges (_m); convert if the file uses _e.
n_e = medium['reaction'].str.endswith('_e').sum()
medium['reaction'] = medium['reaction'].str.replace(r'_e$', '_m', regex=True)
medium = medium.drop_duplicates('reaction')
print(f'Diet: {len(medium)} medium reactions ({n_e} converted _e->_m).')

if tradeoff_path.exists():
    sweep = pd.read_csv(tradeoff_path)
    print(f'Reusing tradeoff sweep ({len(sweep)} rows).')
else:
    strat = (meta.reset_index().groupby('Stage.Group', group_keys=False)
             .apply(lambda g: g.sample(min(len(g), 4), random_state=SEED)))
    pilot_ids = [s for s in strat['Sample'].astype(str) if s in set(manifest['sample_id'].astype(str))][:10]
    pilot_man = manifest[manifest['sample_id'].astype(str).isin(pilot_ids)]
    print(f'Tradeoff pilot on {len(pilot_man)} samples: {pilot_ids}')
    from micom.workflows import tradeoff as micom_tradeoff
    sweep = micom_tradeoff(pilot_man, model_folder=str(MODELS_DIR), medium=medium,
                           tradeoffs=np.arange(0.1, 1.01, 0.1), threads=N_THREADS)
    sweep.to_csv(tradeoff_path, index=False)

# Pick the tradeoff maximising the fraction of taxa with growth above the micom
# default floor (1e-6) averaged over pilot samples — the standard heuristic.
if 'tradeoff' in sweep.columns and 'growth_rate' in sweep.columns:
    frac = (sweep.assign(alive=sweep['growth_rate'] > 1e-6)
            .groupby('tradeoff')['alive'].mean())
    OPT_TRADEOFF = float(frac.idxmax())
    print(frac.round(3).to_string())
else:
    OPT_TRADEOFF = 0.5
    print('Sweep columns unexpected; falling back to tradeoff=0.5.  # >>> VERIFY')
print(f'Chosen tradeoff: {OPT_TRADEOFF:.2f}')'''))

# ------------------------------------------------------------------ 08 smoke gate (1 sample)
cells.append(("code", "nb10-08-smoke", r'''# ---------------------------------------------------------------------------
# Smoke gate: one sample with pFBA before the batch (section 8 budget: <1.2 GB,
# <15 min/sample). If this errors or blows the budget, stop here — do not batch.
# ---------------------------------------------------------------------------
smoke_path = MICOM_WORK / 'smoke_one_v2.pkl'
if smoke_path.exists():
    print('Smoke checkpoint present — skipping.')
else:
    one = manifest.iloc[[0]]
    print(f'Smoke grow() on sample {one["sample_id"].iloc[0]} (pFBA) ...')
    import time, tracemalloc
    tracemalloc.start(); t0 = time.time()
    gr = grow(one, model_folder=str(MODELS_DIR), medium=medium,
              tradeoff=OPT_TRADEOFF, threads=1, strategy='pFBA', presolve=True)
    dt = time.time() - t0; _, peak = tracemalloc.get_traced_memory(); tracemalloc.stop()
    with open(smoke_path, 'wb') as f:
        pickle.dump(gr, f)
    print(f'Smoke PASSED: {len(gr.exchanges)} exchange rows, {dt:.0f}s, ~{peak/1e6:.0f} MB python heap.')
    # sanity: are any tracked polyamines present as exchanges at all?
    ex = gr.exchanges.copy()
    met_col = 'metabolite' if 'metabolite' in ex.columns else 'reaction'
    found = {classify_metabolite(m)[1] for m in ex[met_col].unique()} - {None}
    print(f'  Polyamine exchanges present in smoke model: {sorted(found) or "NONE  # >>> VERIFY VMH ids"}')'''))

# ------------------------------------------------------------------ 09 full batch grow (checkpointed)
cells.append(("code", "nb10-09-batch", r'''# ---------------------------------------------------------------------------
# Full batch, checkpointed in chunks so a crash resumes instead of restarting.
# Set MICOM_RUN_BATCH=1 to run the (hours-long) batch; unset it to stop after the
# smoke gate. Extreme fluxes are flagged, never deleted. Per-chunk growth rates are
# NOT persisted (v2.1 -- nothing read them); the exchange chunks are resume
# checkpoints and are removed once consolidated.
# ---------------------------------------------------------------------------
exch_path = MICOM_WORK / 'exchanges_all_v2.pkl'
RUN_BATCH = os.environ.get('MICOM_RUN_BATCH') == '1'

if exch_path.exists():
    exchanges = pd.read_pickle(exch_path)
    print(f'Loaded cached exchanges: {len(exchanges):,} rows, '
          f'{exchanges["sample_id"].nunique()} samples.')
elif not RUN_BATCH:
    raise SystemExit(
        'Smoke gate done. Set MICOM_RUN_BATCH=1 (see run_micom.sh) to run the full '
        'batch. Stopping here so the notebook does not launch hours of compute on an '
        '--execute pass by accident.')
else:
    CHUNK = int(os.environ.get('MICOM_CHUNK', '25'))
    done_ids, parts = set(), []
    for ckpt in sorted(MICOM_WORK.glob('exchanges_chunk_*.pkl')):
        p = pd.read_pickle(ckpt); parts.append(p); done_ids |= set(p['sample_id'].astype(str))
    todo = [s for s in manifest['sample_id'].astype(str) if s not in done_ids]
    print(f'Resuming: {len(done_ids)} done, {len(todo)} to go, chunk={CHUNK}.')
    for i in range(0, len(todo), CHUNK):
        ids = todo[i:i+CHUNK]
        sub = manifest[manifest['sample_id'].astype(str).isin(ids)]
        gr = grow(sub, model_folder=str(MODELS_DIR), medium=medium,
                  tradeoff=OPT_TRADEOFF, threads=N_THREADS, strategy='pFBA', presolve=True)
        ex = gr.exchanges.copy()
        # v2.1: growth_chunk_*.pkl is no longer written. It was saved "for the record"
        # but the resume logic below globs only exchanges_chunk_*, and a tree-wide
        # read-vs-write cross-reference found no reader for it anywhere -- ~12 files
        # per batch written and never opened.
        ex.to_pickle(MICOM_WORK / f'exchanges_chunk_{i:04d}.pkl')
        parts.append(ex)
        print(f'  chunk {i//CHUNK+1}: samples {i}-{i+len(ids)-1}, {len(ex):,} rows')
    exchanges = pd.concat(parts, ignore_index=True)
    exchanges.to_pickle(exch_path)
    print(f'Batch complete: {len(exchanges):,} exchange rows total.')

    # The chunks exist only so an interrupted batch can resume. Once the consolidated
    # exchanges_all_v2.pkl is on disk they are a second full copy of it, so drop them.
    # Mirrors the models_agora2/*.pickle cleanup further down.
    _n_ckpt = 0
    for ckpt in sorted(MICOM_WORK.glob('exchanges_chunk_*.pkl')):
        ckpt.unlink()
        _n_ckpt += 1
    if _n_ckpt:
        print(f'  removed {_n_ckpt} resume checkpoint(s) now consolidated into '
              f'{exch_path.name}.')

# Flag (do not delete) extreme fluxes.
exchanges['extreme'] = exchanges['flux'].abs() > EXTREME_FLUX
print(f'Extreme-flux rows flagged (|flux|>{EXTREME_FLUX}): {int(exchanges["extreme"].sum())}')'''))

# ------------------------------------------------------------------ 10 per-species polyamine export flux
cells.append(("code", "nb10-10-species-flux", r'''# ---------------------------------------------------------------------------
# Per (sample, species) EXPORT flux of each polyamine. In micom's exchange table a
# taxon row with direction 'export' (or flux>0 on a taxon exchange) is secretion
# into the shared lumen = production available to the host.
# ---------------------------------------------------------------------------
ex = exchanges.copy()
met_col = 'metabolite' if 'metabolite' in ex.columns else 'reaction'
tax_col = 'taxon' if 'taxon' in ex.columns else ('compartment' if 'compartment' in ex.columns else None)
if tax_col is None:
    raise KeyError(f'No taxon column in exchanges; columns={list(ex.columns)}  # >>> VERIFY')

# Restrict to per-taxon rows (drop the community "medium"/external pseudo-taxon).
ex = ex[ex[tax_col].astype(str).str.lower() != 'medium'].copy()

# Classify each exchange row to a tracked polyamine.
kegg_name = ex[met_col].map(lambda m: classify_metabolite(m))
ex['pa_kegg'] = kegg_name.map(lambda t: t[0])
ex['pa_name'] = kegg_name.map(lambda t: t[1])
ex = ex.dropna(subset=['pa_kegg'])

# Export (production) flux: positive secretion. micom convention: export flux is
# typically negative on the taxon exchange; handle both by taking the secretion side.
if 'direction' in ex.columns:
    prod = ex[ex['direction'].astype(str).str.lower().eq('export')].copy()
    prod['export_flux'] = prod['flux'].abs()
else:
    prod = ex.copy()
    prod['export_flux'] = prod['flux'].clip(lower=0)   # fallback: positive = secretion

# (sample, species, polyamine) export flux (sum duplicate reactions of same compound).
sp_flux = (prod.groupby(['sample_id', tax_col, 'pa_name'])['export_flux']
           .sum().reset_index().rename(columns={tax_col: 'Species'}))
print(f'Per-species polyamine export fluxes: {len(sp_flux)} (sample,species,polyamine) rows, '
      f'{sp_flux["Species"].nunique()} species, {sp_flux["pa_name"].nunique()} polyamines.')

# Community net export per (sample, polyamine) -> community flux table.
comm = (sp_flux.groupby(['sample_id', 'pa_name'])['export_flux'].sum()
        .reset_index().rename(columns={'export_flux': 'community_export_flux'}))
comm_wide = comm.pivot(index='sample_id', columns='pa_name', values='community_export_flux').fillna(0.0)
comm_wide.to_csv(TABLES / 'nb10_micom_community_flux_v2.csv')
print(f'Saved community flux: nb10_micom_community_flux_v2.csv ({comm_wide.shape})')'''))

# ------------------------------------------------------------------ 11 species flux vs measured -> pairs table
cells.append(("code", "nb10-11-pairs", r'''# ---------------------------------------------------------------------------
# The NB08 input: for each (species, polyamine), Spearman(species export flux,
# measured metabolite) across samples. Columns Species, Polyamine, Rho, QValue
# are what NB08 renames and merges. BH across all tested pairs, one family.
# ---------------------------------------------------------------------------
# Measured value per (sample, polyamine name).
meas = {}
proxy_flag = {}
for name, kegg in NAME_TO_KEGG.items():
    col, is_proxy = resolve_metabolite_column(kegg)
    if col is not None:
        meas[name] = mtb[col]
        proxy_flag[name] = is_proxy
print('Measured polyamines available:', sorted(meas.keys()))

rows = []
for (species, pa_name), g in sp_flux.groupby(['Species', 'pa_name']):
    if pa_name not in meas:
        continue
    s = g.set_index('sample_id')['export_flux']
    # align flux and measured on shared samples; samples with no model row => 0 flux
    y = meas[pa_name].reindex(mtb.index)
    x = s.reindex(mtb.index).fillna(0.0)
    ok = y.notna()
    n = int(ok.sum())
    # need variation in flux to correlate at all
    if n < 20 or x[ok].nunique() < 3:
        continue
    rho, p = spearmanr(x[ok], y[ok])
    if np.isnan(rho):
        continue
    rows.append({
        'Species': species, 'Polyamine': pa_name, 'KEGG': NAME_TO_KEGG[pa_name],
        'Rho': rho, 'PValue': p, 'N_samples': n,
        'Mean_export_flux': float(x[ok].mean()),
        'Frac_samples_producing': float((x[ok] > 0).mean()),
        'Measured_is_proxy': bool(proxy_flag.get(pa_name, False)),
        'MICOM_status': 'computed',
    })

pairs = pd.DataFrame(rows)
if len(pairs):
    pairs['QValue'] = multipletests(pairs['PValue'], method='fdr_bh')[1]
else:
    pairs['QValue'] = pd.Series(dtype=float)
pairs = pairs.sort_values(['QValue', 'Rho'], ascending=[True, False]).reset_index(drop=True)
pairs.to_csv(TABLES / 'nb10_micom_species_flux_pairs_v2.csv', index=False)
print(f'Saved: nb10_micom_species_flux_pairs_v2.csv  ({len(pairs)} species-polyamine pairs)')
if len(pairs):
    sig = pairs[(pairs['QValue'] < utils.FDR_THRESHOLD) & (pairs['Rho'] > 0)]
    print(f'Positive & significant (q<{utils.FDR_THRESHOLD}): {len(sig)} pairs — '
          f'these set Evidence_mechanistic=True in NB08.')
    print(pairs.head(15)[['Species', 'Polyamine', 'Rho', 'QValue', 'N_samples',
                          'Frac_samples_producing']].round(4).to_string(index=False))'''))

# ------------------------------------------------------------------ 12 cleanup + handoff
cells.append(("code", "nb10-12-cleanup", r'''# ---------------------------------------------------------------------------
# Models are large; delete them once fluxes are saved (section 8) -- but ONLY if this
# notebook built them. Keep the exchange pickles and the three CSVs. Then re-run
# NB08 + NB09.
# ---------------------------------------------------------------------------
if MODELS_ARE_PREBUILT:
    # MODELS_DIR is Data/micom/models_agora2, a read-only INPUT holding ~36 GB of
    # community models that cost hours to build. Deleting an input is never correct,
    # and this guard is the reason MODELS_ARE_PREBUILT exists.
    print(f'Models left in place: {MODELS_DIR} is a pre-built INPUT, not a work product.')
elif os.environ.get('MICOM_KEEP_MODELS') != '1':
    n = 0
    for p in MODELS_DIR.glob('*.pickle'):
        p.unlink(); n += 1
    print(f'Deleted {n} community model pickles (set MICOM_KEEP_MODELS=1 to retain).')

print('\nNB10 done. Outputs:')
for t in ('nb10_micom_sample_coverage_v2.csv', 'nb10_micom_community_flux_v2.csv',
          'nb10_micom_species_flux_pairs_v2.csv'):
    p = TABLES / t
    print(f'  {"OK " if p.exists() else "MISSING "} {t}')
print('\nNext: re-run NB08 then NB09 so the mechanistic column fills:')
print('  jupyter nbconvert --to notebook --execute --inplace 08_pathway_analysis.ipynb')
print('  jupyter nbconvert --to notebook --execute --inplace 09_pathway_network.ipynb')'''))

nb = Notebook.create(Path(__file__).resolve().parent / "10_micom_flux.ipynb", cells)
nb.save()
