# Retired notebook generators

These `edits_nb*.py` scripts used to patch the v2 notebooks through `../../nbedit.py`
and were nominally the source of truth. **They are retired. The notebooks in `v2/` are
authoritative**, and `checks.py` E4 fails if a generator reappears in `v2/`.

## Why they were retired (code review, 2026-10-20)

The generators and notebooks had drifted apart in both directions:

- 12 of 14 generators differed from the notebooks on disk. NB05, NB06, NB08 and NB09
  had been executing pre-v2.1 code while the generators held the v2.1 design.
- Hand fixes made directly in notebooks (NB03 `c13`/`c19`, NB06, NB08, NB09) were not
  mirrored back, so re-running a generator would have silently reverted them.
- `edits_nb08.py`, `edits_nb08b.py` and `edits_nb08c.py` cannot produce a runnable NB08
  in any order: 08b's inserted cell uses names that only `edits_nb08`'s imports cell
  defines, plus a variable that 08c deletes, so NB08 raises `NameError`.
- `edits_nb09_glasso.py` holds a GLASSO configuration that the statistical reviewer
  measured as *worse* than the one in the notebook: every permutation refit fails to
  converge and the module takes about 21 minutes.

## What was transplanted into the notebooks

Applied to a scratch copy, diffed cell by cell, then transplanted by cell id:

| Generator | Cells taken | Why |
|---|---|---|
| `edits_nb02a.py` | NB02 `cell-002`, `cell-007` | shared theme / setup |
| `edits_nb02b.py` | NB02 `cell-014`, `cell-015` | writes the tidy `differential_abundance_v2.csv` |
| `edits_nb02c.py` | NB02 `cell-054`-`056` | stops writing `analysis_results_v2.pkl` (no readers) |
| `edits_nb03.py` | NB03 `c01-imports`, `c06-top15-select` | reads the tidy DA table |
| `edits_nb06.py` | NB06 `cell-001`, `017`, `030`; inserts `cell-029b`, `029c` | PostSurgery projection; corrected stage-group prose |
| `edits_nb06.py` | NB06 `cell-020`, merged | its labels kept; its hardcoded hex palette replaced by `viz.group_palette` |

## What was not transplanted, and why

- **`edits_nb05.py`** (NB05 rebuild: cross-cohort transfer, severity transfer). Its
  transfer stage has two known defects, so it is transplanted only together with their
  fixes: a held-out-prevalence leak in the genus basis, and a detection-limit mismatch
  between the shotgun training cohort and the 16S targets.
- **`edits_nb08.py`, `08b`, `08c`, `edits_nb09.py`**: reference material for the NB08
  evidence-assembler rebuild only. They are not applicable as written (see above).
- **`edits_nb09_glasso.py`**: superseded by the validated GLASSO configuration applied
  directly to NB09.
- **`edits_nb04.py`**: built `04_crc_stage_tumour_site.ipynb`, which was never created;
  NB04 is discarded. Kept as reference for porting its tumour-site exercise.
- **`edits_nb07.py`**: NB07 is markdown-only and slated to be folded into CONTRACTS.md.
- **`edits_nb10.py`**: NB10 (MICOM) is outside the pipeline and needs `micom`, which is
  not in the pinned environment. Paper 2 work.

To change a notebook now, edit the notebook. For a large mechanical change, write a
one-off script with `nbedit.Notebook` and commit it with the change.
