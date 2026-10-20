# v2.1 Pipeline Contracts

The authoritative description of what each notebook consumes and produces, and of the
conventions every notebook must obey. `utils.py` implements these rules; `checks.py` enforces
the ones that can be tested. If this file and a notebook disagree, the notebook is wrong.

All outputs are suffixed `_v2` and live under the **shared** `RESULTS_DIR` (`<project
root>/Results`), so they never collide with `publication`'s or `CRC_Gut_Analysis_Git`'s
outputs even though the tree is shared. `intermediate/` for pickles, `tables/` for CSVs,
`figures/<topic>/` for figures, `logs/` for run logs.

`utils.py` resolves `PROJECT_ROOT` from its own location — it lives at `<project
root>/v2/utils.py`, so the root is its parent's parent (`D:\KI Research` on the PC, `~/ki`
on the server). Set the `KI_PROJECT_ROOT` environment variable to override. `utils` warns at
import if `DATA_DIR` is missing, which is what a moved tree looks like.

Every notebook imports `utils` via `sys.path.insert(0, '.')` (run from `v2/`) — no absolute
paths. Every notebook that uses stage encoding imports it from `utils` and asserts against
it; none defines a local copy.

Run `python checks.py` before any pipeline run. `run_pipeline.py` does this and aborts on
failure.

---

## 1. Cohorts

| Cohort | Role | n (analysed) | Taxonomy | Metabolome | Staging | Confounders |
|---|---|---|---|---|---|---|
| `YACHIDA-CRC-2019` | **discovery** | 277 primary (+30 sensitivity) | shotgun, species | absolute, log10 | full TNM | Age, BMI, sex, Brinkman, alcohol |
| `SINHA_CRC_2016` | validation | 131 (89 Healthy + 42 CRC) | 16S, genus | pre-scaled | none | Age, BMI, sex |
| `Kim_adenomas_2020` | validation | 138 (102 Healthy + 36 CRC) | 16S, genus | pre-scaled | none | banded age, sex, smoking |

> **Corrected 2026-10-07.** This table previously recorded Sinha as 89 and Kim as 102 —
> the control-only counts produced by the `primary_mask()` include-list bug, which dropped
> every case from both validation cohorts. See `utils.primary_mask()` and `checks.py` A6.

`DATASET_DISCOVERY` / `DATASETS_VALIDATION` name these roles.

**`ERAWIJANTARI-GASTRIC-CANCER-2020` must never be used as a cohort.** 53 of its samples are
the same Healthy subjects as Yachida's (`Shared.w.ERAWIJANTARI_2020 == True`), so it is not
independent and would double-count controls.

**The IBD/IBS/ESRD cohorts were removed from the project (2026-10-07).**
`FRANZOSA-IBD-2019`, `iHMP-IBDMDB-2019`, `MARS-IBS-2020` and `WANG_ESRD_2020` are gone from
disk and from `DATASETS_ALL`, `COMPARISONS`, `POLYAMINE_DATASETS` and `STUDY_GROUP_MAP`.
Nothing downstream of NB01 ever used them. `DATASETS_IBD` survives as an empty list because
NB01 imports it to tag exploratory figures; that test is now always `False`.

**Absolute metabolomes** (Yachida) go through `utils.transform_metabolites`: log10(x+1),
probabilistic quotient normalisation on the metabolites detected in ≥95% of samples (61 in
Yachida), and non-detects tied at one floor per metabolite, below every detected value
(`checks.py` C6). **Errata (2026-10-20):** this replaced a sample-wise median taken over all
metabolites with non-detects as 0. With 38.8% non-detects that offset tracked each sample's
detection rate (ρ −0.92) rather than its dilution, and because Advanced CRC samples had
slightly more non-detects (39.4% vs 38.1%), every detected value in that group was shifted up
by ~0.14: 11 robust Healthy-vs-Advanced metabolites, all "elevated". After the fix (and the
cap removal below) the robust counts are Healthy vs Early 21 (6 up / 15 down), Healthy vs
Advanced 19 (7 up / 12 down), Early vs Advanced 0, and the median difference across all
metabolites is within ±0.05 for every comparison.

**Pre-scaled metabolomes** (Sinha, Kim) are detected by `utils.is_prescaled()` (any negative
value) and skip prevalence filtering and the log10/PQN transform.
Applying `log10(x+1)` to a z-score turns every value below −1 into NaN, which is what
silently voided Sinha in the earlier pipeline.

---

## 2. Group policy — the rule everything else depends on

Yachida's raw `Study.Group` values and their v2.1 treatment:

| Raw group | n | Treatment |
|---|---|---|
| `Healthy` | 127 | primary → `Healthy` |
| `Stage_0` | 27 | primary → `Early_CRC` |
| `Stage_I_II` | 69 | primary → `Early_CRC` |
| `Stage_III_IV` | 54 | primary → `Advanced_CRC` |
| `HS` ("healthy after surgery") | 30 | **sensitivity only** → `PostSurgery` |
| `MP` (multiple polypoid adenomas) | 40 | **excluded from the study** |

- **MP is dropped at load** by `utils.load_all_datasets()` via `EXCLUDED_STUDY_GROUPS`, before
  any QC step, so no downstream notebook can see an adenoma sample. Kim's own `Adenoma` arm
  (102 samples) is dropped for the same reason, leaving Control vs Carcinoma.
- **HS is loaded but quarantined** (`SENSITIVITY_GROUPS`). Those patients had colorectal
  cancer resected and are currently cancer-free: a resected gut is neither a clean control
  nor an active tumour, so pooling it either way biases the comparison. NB01b measures how
  far it sits from Healthy.
- **`PostSurgery` has no position on any ordinal scale.** `STAGE_GROUP_NUMERIC['PostSurgery']`,
  `CRC_STAGE_NUMERIC['HS']` and `STUDY_GROUP_MAP['HS']` are all `NaN`, so it cannot enter a
  trend test even if a caller forgets to mask it.

Counts: **307 loaded, 277 primary, 30 sensitivity, 150 carcinomas.**

Earlier pipelines were wrong in different ways: the `publication/` pipeline put HS in
Early_CRC and MP in Advanced_CRC; the first v2 run dropped HS and called MP "Adenoma".
`publication/` is unchanged and still uses the old grouping.

### Masks (the enforcement mechanism)

| Helper | Selects | Used by |
|---|---|---|
| `primary_mask(meta)` | Healthy + Early_CRC + Advanced_CRC | every primary analysis |
| `sensitivity_mask(meta)` | PostSurgery only | NB01b §4b, NB06 §14 |
| `crc_only_mask(meta)` | Early_CRC + Advanced_CRC (n=150) | NB05 stage 2 |
| `apply_mask(mask, *frames, meta_df=)` | filters frames **by Sample ID**, not row position | all of the above |

**NB02 performs the primary split once**, so `preprocessed_data_reduced_v2.pkl` holds primary
samples under its normal keys and the sensitivity rows under `sensitivity_species` /
`sensitivity_mtb` / `sensitivity_meta`. A downstream notebook that simply loads
`transformed_species` therefore gets the correct sample set by default — the failure mode of
forgetting is absence, not contamination.

### Comparisons

| Constant | Contents |
|---|---|
| `COMPARISONS_GROUP` | Healthy vs Early_CRC, Healthy vs Advanced_CRC, Early_CRC vs Advanced_CRC |
| `COMPARISONS_SENSITIVITY` | Healthy vs PostSurgery — **its own FDR family**, never merged |
| `COMPARISONS` (fine-grained) | Healthy vs Stage_0 / Stage_I_II / Stage_III_IV; validation cohorts: Healthy vs CRC |

### Differential abundance output — one tidy table (v2.1)

`tables/differential_abundance_v2.csv` is the DA table every downstream reader uses:

```
Dataset | Comparison | Scheme | Assay | Feature | MeanDiff | PValue | QValue
        | CliffsD | Coef_adj | QValue_adj | Agreement | Robust
```

It holds the **confounder-adjusted (`grouped`) scheme only**, and replaces the 16-file
`{ds}_{g1}_vs_{g2}_{species|metabolites}_DA_v2.csv` family. That family was two *different*
statistics sharing one filename template — Step 6's unadjusted Mann-Whitney and Step 6B's
adjusted model — so Yachida carried two answers to one question under near-identical names,
and because `Stage_0 ∪ Stage_I_II == Early_CRC` they overlap by construction. 12 of the 16
had no reader at all, and the two live readers disagreed about which was authoritative.
Both now slice this table: NB03 (`c06-top15-select`) and `ml_workbench/engine.py`.

**Source of truth (2026-10-20).** The notebooks in `v2/` are authoritative. The
`edits_nb*.py` generators had drifted from them in both directions and are retired to
`archive/generators/` (see its README for what was transplanted and why the rest was
not); `checks.py` E4 fails if one reappears in `v2/`.

Step 6's fine-grained DA is still *computed* (it feeds the volcano plots and the summary) but
is **no longer written to disk**.

**Two DA tables remain outside it, deliberately:**
`CRC_StageI_II_vs_III_IV_{species,metabolites}_DA_v2.csv` (Step 10B). That contrast excludes
`Stage_0`, whereas the grouped scheme's `Early_CRC_vs_Advanced_CRC` pools it into the early
arm — a genuinely different question, which §10B already states is reported alongside on
purpose. They keep their own names rather than being folded in under a `Scheme` value,
because a reader filtering on `Comparison` would otherwise double-count the same samples.

**There is no `analysis_results_v2.pkl`.** It had zero readers and duplicated these CSVs;
`full_corr` (all species × all metabolites) was its largest key and nothing read it.

---

## 3. Tumour location — admissible only within cancer

`Tumor location` is `"-"` for **every** non-cancer sample (all 127 Healthy, all 30 HS, all 40
MP). On any set containing non-cancer samples, a tumour-site covariate is an exact
disease-status indicator. Measured cost of the old leak: NB05 stage-1 accuracy 0.63 → 0.39,
κ 0.48 → 0.08.

`utils.build_confounder_matrix(meta, include_tumor_location=True)` **raises `ValueError`**
unless every row is `Early_CRC` or `Advanced_CRC`. Enforced in code, not left as a convention.

| Column | Values | Use |
|---|---|---|
| `Tumor.Location.Raw` | original string | provenance |
| `Tumor.Location.Category` | NoTumor / Distal / Proximal / MultiSite / Unknown | descriptive only |
| `Tumor.Location.Site` | Rectum / Left_colon / Right_colon / MultiSite, **NaN** for non-cancer | the analysis column |
| `TumorSite_*` dummies | from `build_confounder_matrix(..., include_tumor_location=True)` | NB05 stage 2, workbench CRC scope |

Within-cancer site counts: Rectum 69, Right_colon 44, Left_colon 32, MultiSite 5.

Permitted uses: **NB05 stage 2** (Early vs Advanced, with and without site on identical
folds), **ML Workbench** in "CRC only" scope. Nowhere else. (NB04, which also used it, is
discarded.)

---

## 4. Differential abundance — two tests, agreement required

| Function | Test | Adjusted | Compositional |
|---|---|---|---|
| `differential_abundance` | Mann-Whitney U + Cliff's δ | no | no |
| `clr_lm_da` | linear model on CLR values, LinDA bias correction | **yes** | **yes** |

- **Species:** `clr_lm_da(..., already_clr=True, bias_correct=True)`.
- **Metabolites:** `clr_lm_da(..., already_clr=True, bias_correct=False)` — the metabolite
  matrix is log10, not a composition, so LinDA's mode subtraction does not apply.
- Covariates: `YACHIDA_COMPLETE_CONFOUNDERS` (Age, BMI, Gender, Brinkman.Index, Alcohol),
  complete-case. Missing covariates are **dropped, never mean-filled**: filling a covariate
  and then partialling it out understates uncertainty for exactly the samples where it was
  unknown.

### Sign conventions (fixed, relied on downstream)

| Column | `> 0` means |
|---|---|
| `MeanDiff` = mean(group2) − mean(group1) | higher in **group2** |
| `Coef` (adjusted model) | higher in **group2** |
| `CliffsD` | higher in **group1** — the **opposite** sign |

There is no `Log2FC` anywhere: the old `log2((mean2+eps)/(mean1+eps))` took a log of a ratio
of already-log/CLR-transformed (frequently negative) means and was NaN whenever a group mean
was negative.

### Robust hits

`Robust == True` ⟺ q < 0.05 in **both** tests **and** matching direction. `Agreement` labels
the alternatives: `both_q_same_dir`, `both_q_opposite_dir`, `mw_only`, `lm_only`, `neither`.
**NB03 selects targets only from `Robust == True`.** `mw_only` features are reported, never
modelled — in this cohort the median Brinkman index is 100 in Healthy vs 196–390 in the CRC
groups, so they are the prime suspects for confounding.

---

## 5. Polyamines and proxies

`POLYAMINE_KEGG`: C00134 Putrescine, C00315 Spermidine, C00750 Spermine, C00179 Agmatine,
C01672 Cadaverine, C00077 Ornithine, C00062 Arginine.

`POLYAMINE_EC` lists **producing** reactions only. The pre-v2 map listed enzymes that
*consume* ornithine and arginine as their producers, which inverted the producer call for
those compounds. `checks.py` D1 asserts both the required producers and the forbidden
consumers; D5 asserts each compound has at least one EC that AGORA2 actually annotates.
Known gap: EC 4.1.1.96 is absent from AGORA2, so spermidine producer calls rest on 2.5.1.16
alone.

`POLYAMINE_PATHWAY_KEGG` excludes ornithine and arginine: their synthesis ECs appear in ~90%
of AGORA2 species and carry no pathway-membership information.

### Proxies

Cohort coverage (`tables/polyamine_cohort_coverage_v2.csv`):

| Compound | Yachida | Sinha | Kim |
|---|---|---|---|
| Putrescine | direct | direct | **proxy** (N-acetylputrescine, C02714) |
| Cadaverine | direct | direct | **proxy** (N-acetyl-cadaverine, by name) |
| Spermidine | direct | — | direct |
| Ornithine / Arginine | direct | direct | direct |
| Spermine / Agmatine | direct | — | — |

Proxies are **pre-specified** (`POLYAMINE_PROXY_KEGG`, `POLYAMINE_PROXY_NAMES`), held in a
separate `polyamine_proxies` dict, and used only where the parent is unmeasured.
N-acetylation is an inactivation step, so a derivative tracks its parent loosely — every
downstream row records `Measured_as` / `Proxy_used_in`, and a replication resting on a proxy
is weaker evidence than one on the parent. Polyamines **and** their proxies are exempt from
every prevalence and variance filter.

---

## 6. Notebook contracts

### NB01 `01_preprocessing.ipynb`
**Reads** `Data/*.tsv`. **Writes** `intermediate/preprocessed_data_v2.pkl` (< 200 MB,
asserted), `tables/polyamine_cohort_coverage_v2.csv`, QC/DA figures.
Pickle keys: `transformed_species`, `transformed_mtb`, `filtered_species`, `filtered_mtb`,
`harmonized_meta`, `polyamine_columns`, `polyamine_proxies`, `mtb_kegg_map`,
`da_results_spc`, `da_results_mtb`, `mtb_prescaled`, `mtb_normalisation`, `excluded_study_groups`,
`sensitivity_groups`, `stage_group_order`, `sensitivity_stage_groups`, `datasets`,
`dataset_discovery`, `datasets_validation`.
All 3 cohorts, 307 Yachida rows (MP already dropped). Asserts every metabolite matrix is
finite and that every polyamine/proxy column survives filtering.

### NB01b `01b_yachida_significant_features.ipynb`
**Reads** NB01's pickle. **Writes** `tables/YACHIDA_v2_significant_*.csv`,
`tables/YACHIDA_v2_sensitivity_healthy_vs_postsurgery_{species,metabolites}.csv`,
`intermediate/yachida_significant_subset_v2.pkl`.
3 primary comparisons plus the Healthy-vs-PostSurgery sensitivity analysis in **its own FDR
family**. The downstream subset is primary-only (asserted). Gate: q<0.05 **and**
|Cliff's δ| > 0.33 — a rank-based effect size, because `MeanDiff` in CLR/log10 space rarely
exceeds 0.7 and the old `|MeanDiff| > 1.0` gate silently zeroed every comparison.

### NB02 `02_analysis.ipynb`
**Reads** NB01's pickle. **Writes** `intermediate/preprocessed_data_reduced_v2.pkl`
and tables:
`differential_abundance_v2.csv` (tidy; `Robust`/`Agreement`/`*_adj`; replaces the
per-comparison `{ds}_{g1}_vs_{g2}_{species|metabolites}_DA_v2.csv` family),
`da_robust_hit_summary_v2.csv`, `{ds}_polyamine_partial_correlations_v2.csv`,
`consistent_polyamine_correlations_genus_v2.csv`, `{ds}_genus_polyamine_correlations_v2.csv`,
`polyamine_DA_summary_v2.csv`, `replication_genus_metabolome_v2.csv`, plus progression and
network tables.
**Performs the primary split.** Feature reduction: top-500 species by variance, computed on
primary samples only. **Metabolites are not capped** — every metabolite surviving NB01's
label-blind prevalence filter is kept (Yachida 279, Sinha 530, Kim 462; asserted).
**Errata (2026-10-20):** a `MAX_MTB_REDUCED=150` variance cap used to apply here. It was the
only filter on the pre-scaled Sinha/Kim metabolomes (discarding ~70%), it left only 15/113 and
24/152 of the compounds identity-matched to Yachida, and it removed most of the robust
early-stage signal. `full_corr` and the binary-stage correlations also lost their 200-species /
500-metabolite caps and two contaminations: species significant against polyamines on the same
data topped up into a "top-variance" set, and forced polyamine inclusion.
Partial correlations adjust for the full confounder set and are **NB09's SHAP-independent
candidate pool**. They are computed by `utils.partial_spearman_grid`, which reproduces
`pingouin.partial_corr(method="spearman")` exactly (`checks.py` C7; on the real data max
|Δr| 1e-15) in one matrix product instead of one call per pair — the old loop's
`except Exception: continue` dropped any pair that raised without a trace. Covariates are
per cohort (`utils.COHORT_CONFOUNDERS`): Yachida age, BMI, sex, Brinkman index, alcohol;
Sinha age, BMI, sex; Kim age, sex, smoking (Kim was adjusted for age and sex only before
2026-10-20, because Yachida's list was applied to every cohort).

Cross-cohort replication is **genus level** (species names cannot match between shotgun and
16S), from `utils.genus_clr_table`, and runs at two scopes:
- **Metabolome-wide** (`replication_genus_metabolome_v2.csv`, `utils.genus_metabolite_replication`,
  `checks.py` C9): every Yachida genus × metabolite Spearman association at q < 0.05 is
  re-tested in each validation cohort where it is **measurable** — the compound matched by
  identity (KEGG ∪ HMDB intersection) and the genus present. Pairs that cannot be measured are
  reported as `Not_measurable_*`, never as failed replications; validation q-values are BH
  over the measurable discovery hits; every row carries `Replication_tier`. A pair significant
  in a validation cohort with the opposite sign is `Contradicted`, not replicated.
- **Polyamine-targeted** (`consistent_polyamine_correlations_genus_v2.csv`): metabolites
  matched by **KEGG parent** compound; `Same_direction` requires sign agreement.

### NB03 `03_ml_benchmarking.ipynb`
**Reads** the reduced pickle + NB02's grouped DA tables. **Writes**
`tables/ml_benchmark_dysregulated_metabolites_v2.csv`,
`tables/polyamine_producer_candidates_v2.csv`, `intermediate/ml_results_v2.pkl`,
`logs/nb03_progress_v2.txt`.
Targets: `Robust == True` metabolites (≤ 8 elevated + 7 depleted per comparison) **plus** the
7 pre-specified polyamines. Protocol: nested CV for `AUC_nested`; fixed-parameter permutation
test (`N_PERM`, default 100, `NB03_N_PERM` overrides) for `Perm_pvalue` → BH `Perm_qvalue`
(the p-value refers to `AUC_fixed`, not `AUC_nested`); signed OOF SHAP with `Direction` and
bootstrap rank CIs. No oversampling (a median split is balanced by construction). Minority
class < 20 → skipped and logged. Verifies the tumour-location guard is live.

**Species given confounders (2026-10-20).** Features are the species only; each target is
the OLS residual of the metabolite on Age, BMI, sex, Brinkman index and alcohol
(`utils.residualise_on_confounders`, `checks.py` C10), split at its median. **Errata:** the
confounders used to be features, so the label permutation tested "target independent of
species *and* confounders" — the confounders alone predict putrescine's split at AUC 0.67 —
and a significant p said nothing about the species. Calibration on the real Yachida species
and confounders, 100 synthetic targets that depend on the confounders only (a valid test
rejects ~5%): the old test rejected **100%** (mean AUC 0.81); the residualised target
rejects **6.0%** (continuous target) and **5.0%** (zero-inflated, 60% at the floor), mean
AUC 0.49, p-values slightly conservative. Residualising the *species* instead is invalid
(99% rejection): a zero-abundance sample's residual is a linear function of the confounders,
and the trees recover it. `R2_confounders`, `AUC_confounders_raw` (confounders alone, raw
split) and `AUC_confounders_adjusted` (confounders alone, adjusted split; ≈ 0.5 unless
confounding is nonlinear) are reported per target. Confounders are imputed once from the
confounders alone (§7).

### NB04 — discarded (2026-10-20)
`04_mp_trajectory.ipynb` was superseded (its subject, the MP arm, is excluded) and its
planned replacement `04_crc_stage_tumour_site.ipynb` was never built; neither exists. Tumour
site as a covariate lives in NB05 stage 2 (CRC-only, with and without site on identical
folds), stage trends in NB02, and the guard demonstration is `checks.py` B2. No notebook
reads an `nb04_*` table.

### NB05 `05_two_stage_model.ipynb`
**Reads** the reduced pickle, `ml_results_v2.pkl`, and the raw TSVs: stage 3 with the
default exclusions, stage 6 with `exclude={}`. Transfer features come from **unfiltered**
species tables through `utils.transfer_genus_features` (`checks.py` C8).
**Writes** `tables/nb05_stage2_early_vs_advanced_v2.csv`,
`tables/nb05_external_transfer_v2.csv`, `tables/nb05_severity_transfer_v2.csv`,
`tables/two_stage_model_results_v2.csv`, `tables/two_stage_shap_<polyamine>_v2.csv`,
`figures/ml/nb05_external_transfer_roc_v2.png`,
`figures/ml/nb05_severity_transfer_v2.png`.
Stage 1 Healthy vs CRC (no site). Stage 2 Early vs Advanced (site allowed; with/without
compared on identical folds). **Stage 3 external transfer**: a genus-level model trained on
Yachida only, applied unchanged to Sinha and Kim, AUC with a 2000-sample bootstrap CI, an
exact one-sided Mann–Whitney p against chance, a presence/absence sensitivity model, and the
within-Yachida CV AUC as the ceiling. Features: genus basis prevalent in Yachida and merely
present in the target (`shared_genera` on unfiltered tables), Yachida proportions below the
target's detection limit set to zero, one shared CLR zero-replacement value, per-cohort
standardisation (label-blind). **Errata (2026-10-20):** the first version chose the basis on
NB01's prevalence-filtered target table — the held-out cohort helped choose the features (Sinha
37 genera instead of 59) — and gave each cohort its own CLR delta although shotgun Yachida had
5–12% zeros on the basis against the targets' 55–84%. Corrected results: Yachida → Sinha AUC
0.655 (p = 0.002; presence-only 0.567); Yachida → Kim AUC 0.452 (p = 0.80, does not transfer). Stage 4 asks whether stage probabilities add anything
to metabolite prediction. Stage 5 is per-polyamine SHAP attribution. Leaf notebook — nothing
downstream reads it.

#### Stage 6 — cross-cohort severity transfer (the one opt-in analysis)

The **only** place in the pipeline where the adenoma groups appear. It re-admits Yachida `MP`
(n=40) and Kim `Adenoma` (n=102) via `exclude={}` passed to **both** `load_all_datasets` and
`harmonize_metadata` — passing it to only the first is silently undone by the second.
`STAGE_GROUP_MAP` is untouched, so the primary set stays n=277 and checks A1/A7 still pass.
Labels come from `SEVERITY_TRAIN_MAP` (4-level, Yachida) and `SEVERITY_EVAL_MAP` (3-level,
cross-cohort), which are **separate constants** from the primary stage encoding.

**What it can and cannot test.** Neither validation cohort has a staging column — Sinha's
`Study.Group` is a bare `0`/`1`, Kim's is `Control`/`Adenoma`/`Carcinoma`. TNM stage is
therefore **unfalsifiable outside Yachida**, and no notebook may claim to predict it there.
3-level severity (Healthy → Adenoma → Carcinoma) is the finest axis Kim can confirm. Sinha
has no adenoma arm and is **not** scored for severity at all.

**Direction and scope.** Yachida → Kim only. Shared genera, measured at runtime (prevalent
in Yachida's labelled training samples, present in the target): stage 3 Yachida ∩ Sinha = 59,
Yachida ∩ Kim = 398; stage 6 (adenomas admitted) Yachida ∩ Kim = 397. The earlier figures of
60 and 408 were not reproduced by any code path and are withdrawn. This is **single-direction external validation, not
LODO**: only Yachida has an adenoma-plus-carcinoma gradient, so holding it out leaves no
staged training data. Binary LODO across all three cohorts is stage 3.

**Reported with negative controls**, because a transfer number is not interpretable without
them: a 999-permutation label null with a p-value (was 20 permutations under an "exceeds
every permutation" rule, whose p floor is 1/21), and a cohort-separability AUC measured on the
harmonised features **before** per-cohort standardisation (after it, cross-validated cohort
membership anti-learns and the AUC is meaningless — the first run reported 0.130). Current
result: within-Yachida ρ +0.190, transferred to Kim ρ −0.059, permutation p = 0.33 — the
severity gradient **does not transfer**; separability AUC 1.000, so shotgun and 16S remain
perfectly distinguishable on these genera even after detection harmonisation. Class cut-points come from
the **training** score distribution; step-wise AUROCs (Control-vs-Adenoma and
Adenoma-vs-Carcinoma) are reported **separately**, as a pooled number hides which step fails.

### NB06 `06_multiomics_factors.ipynb`
**Reads** the reduced pickle. **Writes** `tables/mofa_factor_stage_association_v2.csv`,
`tables/mofa_factor_loadings_v2.csv` (**every** factor × **every** feature — a hard
requirement for NB08/NB09), `tables/mofa_factor_scores_v2.csv`,
`tables/mofa_postsurgery_projection_v2.csv`.
**Single-group** MOFA+: per-stage groups would centre away the between-stage differences the
stage test looks for, which is what the earlier pipeline did. Fitted on primary samples only;
factor sign canonicalised so **higher = more advanced disease**. PostSurgery is **projected
afterwards**, never used in fitting, so comparing it against the axes is not circular.
Undetected metabolites → `NaN` via NB02's `undetected_mtb` mask, for MOFA+'s native
missing-value handling.

### NB07 `07_mediation.ipynb`
Markdown only — documents why mediation is excluded: circular mediator selection, a
hardcoded 0/1 treatment contrast against an unscaled exposure, and cross-sectional data that
cannot order cause and effect.

### NB08 `08_pathway_analysis.ipynb`
> **STATUS (2026-10-20) — the design below is NOT what executes.** The executed NB08 is the
> polyamine-era evidence notebook: it writes the `r1_*`, `r3_*`, `r4_*`, `r5_*` and `r8_*`
> tables but **not** `pair_evidence_nb08_v2.csv`, which only the retired `edits_nb08b.py`
> produced (not applicable as written; see `archive/generators/README.md`). The schema below
> is the design target for the NB08 rebuild as the metabolome-wide evidence assembler. Until
> then, do not cite `pair_evidence_nb08_v2.csv` or its columns.

**Reads** the reduced pickle, `ml_results_v2.pkl`, MOFA tables, NB02 partial correlations and
replication, and optionally NB10's MICOM pairs.
**Primary output: `tables/pair_evidence_nb08_v2.csv`** — one row per (species, polyamine),
**read by NB09**. Also `r1_agora2_producer_enrichment_v2.csv`,
`r1_agora2_enrichment_test_v2.csv`, `r3_scfa_polyamine_enrichment_v2.csv`, and the
supplementary `r4_*`, `r5_*`, `r8_*` tables.

`pair_evidence_nb08_v2.csv` schema (candidate pool = partial correlation q<0.05, both signs):

| Column | Meaning |
|---|---|
| `Species`, `Genus`, `Polyamine`, `KEGG_parent`, `Metabolite`, `Measured_as` | identity |
| `Rho_partial`, `QValue_partial`, `PValue_partial`, `Covariates` | NB02 statistical anchor |
| `Candidate_producer` | `Rho_partial > 0` |
| `Mean_SHAP_OOF`, `Mean_SHAP_signed`, `SHAP_dir_rho`, `Direction`, `Rank`, `Rank_stable`, `Model_AUC_nested`, `Model_perm_q` | NB03 |
| `Evidence_predictive` | positive-direction SHAP **and** `Model_perm_q < 0.05` |
| `AGORA2_encodes_enzyme` (True/False/**None**), `AGORA2_model_available` | AGORA2 |
| `Evidence_genomic` | `AGORA2_encodes_enzyme is True` |
| `MICOM_flux_rho`, `MICOM_flux_q`, `MICOM_status` (`computed`/`pending`/`not_modelled`) | NB10 |
| `Evidence_mechanistic` | MICOM q<0.05 and rho>0 |
| `Replication_n_cohorts`, `Replication_mean_rho`, `Replication_same_direction`, `Replication_cohorts`, `Replication_proxy_used` | NB02, genus level |
| `Evidence_replication` | `Replication_same_direction` |
| `Polyamine_stage_rho/q`, `Species_stage_rho/q` | NB04 |
| `SCFA_producer`, `Is_Fusobacterium` | context |
| `Evidence_count`, `Evidence_sources`, `Tier` | NB08's 4-source count |

R1's enzyme evidence is **tri-state**: `True`/`False` only where AGORA2 has a reconstruction;
`None` means no evidence either way and is excluded from the enrichment test rather than
counted as a negative. The SCFA guild is **EC-defined** (`utils.SCFA_EC`), not a curated genus
list drawn from this cohort's own SHAP output. R4/R5/R8 are **supplementary descriptive
indices** and contribute no evidence column.

### NB09 `09_pathway_network.ipynb`
> **STATUS (2026-10-20):** as executed, NB09 reads NB02's polyamine partial correlations and
> NB08's `r1`/`r8` tables directly and scores **three** sources (NB03 SHAP, AGORA2 — always
> unavailable — and MOFA co-loading). Reading `pair_evidence_nb08_v2.csv` and the five-source
> count below take effect with the NB08 rebuild.

**Reads** `pair_evidence_nb08_v2.csv` (**required** — raises if absent and never recomputes
it), plus the reduced pickle and MOFA loadings.
**Writes** `tables/final_evidence_synthesis_v2.csv`,
`tables/d_glasso_partial_correlation_edges_v2.csv`, Module B/C/G tables,
`figures/pathway_network/*` (including `final_evidence_synthesis_heatmap_v2.png`).
Adds the one source NB08 cannot compute — **MOFA same-sign co-loading** (|weight| > 0.10) —
then counts **five independent sources** and applies pre-declared tiers: `Strong` ≥ 3,
`Moderate` = 2, `Weak` = 1, `Statistical only` = 0. The statistical association is **not**
counted; it is how a pair was nominated, so counting it would give every row a free point.
GLASSO edges use a **Westfall–Young max-T** permutation FWER test (per-edge BH on a precision
matrix was structurally unable to reach significance).

### NB10 `10_micom_flux.ipynb` *(separate run — see §8)*

---

## 7. Conventions that apply everywhere

- **FDR:** `utils.FDR_THRESHOLD = 0.05`, Benjamini-Hochberg. One family per table;
  `bh_family()` corrects several p-columns jointly and emits `Q_<col>`.
- **Missing coverage is never negative evidence.** `None`/`NaN` means "not assessable".
- **No live network calls.** AGORA2 enzyme evidence comes from the local index at
  `Results/tables/agora2_species_ec_index.json`, rebuilt from `Data/micom/agora2_extracted/`
  by `utils.build_agora2_ec_index` (~40 s, 1439 species).
- **Imputation never sees a target or a held-out fold.** Where confounders are model
  features (NB05), their imputation is a pipeline step refit per CV fold. Where they adjust
  a target (NB03), they are imputed once from the confounders alone — no target, no species
  — so no label information can enter the folds. Targets are never imputed — rows with a
  missing target are dropped.
- **Genus matching** uses `utils.normalise_genus_name` (strips GTDB `_A` suffixes) and
  `utils.genus_clr_table`, shared by NB02's replication and NB05's transfer model so the two
  cannot drift apart.
- **Cross-cohort genus features: fix the basis BEFORE the CLR.** Any model spanning cohorts
  must pick its genera with `utils.shared_genera()` and build every cohort's table with
  `utils.genus_clr_table(..., genera=shared)`. CLR divides by the geometric mean of the
  columns present, so CLR-ing each cohort on its own genus set and intersecting the columns
  *afterwards* leaves the cohorts on different scales — Yachida's reference spans its own
  prevalent genera, Kim's ~409 — which biases the transfer systematically. NB05 stage 3 did
  exactly that before v2.1. With `genera=` the prevalence filter is skipped and absent
  genera become zero columns, so every cohort returns identical columns in identical order
  (scikit-learn indexes features positionally).
  `shared_genera()` is deliberately **asymmetric**: prevalence is applied to the *training*
  cohort only, and targets are intersected on mere *presence*. Thresholding the target would
  let held-out data choose the model's features; presence is a property of that cohort's
  assay coverage, not of its outcome.
- **Compute (revised 2026-10-18 — the analysis server is no longer available).** Execution
  moved to the local Windows workstation: **8 physical cores / 16 logical** (Ryzen 7 6800H),
  **31.2 GB RAM**, Python 3.12.10 in `.venv`. The former 1500 MiB / 2-core slice no longer
  applies, and `run_pipeline.sh` is superseded by `run_pipeline.py` because bash is not
  available as the execution host.
  Two of the old constraints are **kept, on their own merits** — they were never purely
  concessions to the small slice:
  - `max_bin=64` for XGBoost. With n=277 primary samples, 256 bins cannot be populated
    meaningfully; 64 is justified by sample size and is ~3× faster per fit.
  - `n_jobs=1` on any estimator fitted *inside* an already-parallel loop
    (`ml_workbench/engine.py`, `run_cv_auc_cached`). This prevents oversubscription at any
    core count, not just two.

  What changes is the **outer** parallelism and the thread budget. BLAS threads multiply with
  joblib workers, so 8 × 8 = 64 threads would thrash 8 cores. `run_pipeline.py` therefore sets
  the thread environment **per notebook**: NB03, the only notebook that parallelises at the
  joblib level, gets `OMP_NUM_THREADS=1` and spends the cores on workers via `KI_N_JOBS=8`;
  every other notebook receives the cores through BLAS (`OMP_NUM_THREADS=8`).
  Resume rule: `run_pipeline.py` reuses a notebook's previous result only if it succeeded
  from the same source, under the same result-affecting profile (`NB03_N_PERM`,
  `NB03_N_ITER`, `NB03_MAX_TARGETS`), **and finished after every upstream notebook's latest
  run** — so a re-run of NB02 always forces NB03 onward to recompute. `--smoke` runs a
  shake-out profile (NB03 capped to 3 targets, `N_ITER=2`, `N_PERM=19`) and writes
  `Results/SMOKE_RUN.txt`; a full run recomputes everything and deletes the marker.
  `GraphicalLassoCV`'s `n_jobs=1` was a 1500 MiB concession and may rise now.

---

## 8. MICOM (NB10) — decoupled run

**Not runnable in the current environment:** `micom`/`cobra` are not in the pinned
environment and `run_micom.sh` (a bash runner for the retired server) refuses to start. Once
NB10 is redesigned and run, **re-run NB08 and NB09** so the mechanistic column fills; until
then `MICOM_status == 'pending'`.

### Attainable evidence ceiling — read this before interpreting any `Tier`

**As executed today, NB09 counts three sources** — NB03 SHAP, AGORA2 genomic and MOFA
co-loading — and genomic is disabled, so the attainable maximum is **2 of 3** and its
`Strong (3/3)` tier is **unreachable**. In the smoke run every supported row was `Weak (1/3)`,
and those came from MOFA co-loading alone, which the integrity review measured firing near
its chance rate (council item C11). The five-source design below applies after the NB08
rebuild.

The five-source design (§7) has two sources currently unavailable:

| Source | State |
|---|---|
| `predictive` (NB03 SHAP) | available |
| `replication` (NB02, genus level) | available |
| MOFA same-sign co-loading (NB09) | available |
| `genomic` (AGORA2) | **DISABLED** — `utils.AGORA2_ENABLED is False`, enforced by `checks.py` D5 |
| `mechanistic` (MICOM, NB10) | **pending** until the decoupled batch runs |

So the attainable maximum is **3 of 5**, not 4 of 5. `Strong` (≥ 3) is therefore reachable
only when *all three* available sources fire on the same pair — a far more fragile tier than
"3 of 5" suggests, and `Moderate`/`Weak` cannot be read as "failed to reach 5".
**Running MICOM restores a fourth source and makes `Strong` 3-of-4 rather than 3-of-3.**
State the ceiling alongside any table or figure that carries a `Tier`.

Inputs: `Data/micom/agora2_gtdb_species.qza` (GTDB-native, so Yachida's names match without
renaming), `Data/micom/western_diet_gut.csv`, and — when present —
`Data/micom/models_agora2/` (**347 pre-built AGORA2 community models + `manifest.csv`,
covering all 277 primary Yachida samples**). NB10 uses that set in place as a **read-only
input**, skips the build, and never deletes it (`MODELS_ARE_PREBUILT`); move it aside to
force a rebuild. It is ~36 GB and hours of solver time already spent — do not prune it.
Outputs: `intermediate/micom_v2/` per-sample checkpoints,
`tables/nb10_micom_sample_coverage_v2.csv`, `tables/nb10_micom_community_flux_v2.csv`,
`tables/nb10_micom_species_flux_pairs_v2.csv` (**read by NB08**).

Fixes relative to the earlier MICOM attempts (publication NB09 and the Git notebook):

| Earlier problem | NB10 |
|---|---|
| AGORA103, ~42% abundance coverage | AGORA2 GTDB (~56% median) |
| GTDB strains collapsed to NCBI binomials and summed | no collapsing; a taxon matches exactly or is reported unmatched |
| tradeoff hardcoded at 0.5 | chosen by a `tradeoff` sweep on 10 stratified pilot samples; value and sweep saved |
| `strategy="none"` → non-unique fluxes reaching 98+ | `strategy="pFBA"` |
| growth rates discarded; fluxes ≥ 80 deleted | growth rates kept; extreme fluxes flagged, never deleted |
| statistics on positive fluxes only (so the check was true by construction) | all samples, zeros included |
| equal-abundance control silently reused the weighted models | separate directory per run |
| HS→Early, MP→Advanced grouping | v2.1 groups |
| diet-file duplicates ignored | duplicates removed and logged; medium checked against model exchanges |
| hand-patched `community.py` | stock `micom` with OSQP/HiGHS from pip |

Per-sample budget: < 1.2 GB, < 15 min, enforced by a smoke gate on one sample before the
batch. Models this notebook *builds* are deleted after their fluxes are saved, so disk use
stays small; a pre-built input set is exempt (see above). Coverage < 0.5 is flagged and
excluded in a sensitivity analysis.

---

## 9. ML Workbench (`ml_workbench/`)

`ENGINE_VERSION = 4` (invalidates the older fold cache). Scope selector: `All (Healthy+CRC)`
or `CRC only`. The `tumour_location` feature group exists **only** in CRC-only scope, with a
rectum/left/right subgroup filter that disables runs where a class has < 20 samples. A
`PostSurgery` toggle (off by default) only ever *scores* those samples, never trains on them.
Candidate targets come from the robust-hit tables, not a nominal-p screen.

> This last rule was **documented but not implemented** until v2.1:
> `candidate_targets()` read the fine-grained *unadjusted* DA family and screened on
> `PValue < 0.05`, while NB03 used the confounder-adjusted `Robust` flag — so the two
> disagreed about which statistics were authoritative on the same cohort. Both now slice
> `differential_abundance_v2.csv` and screen on `Robust`. Keep `ENGINE_VERSION` here and in
> `engine.py` in step; it was 2 in code while this file said 3, so the fold cache was not
> being invalidated as documented.

---

## 10. Pre-flight checks (`checks.py`)

| Group | Asserts |
|---|---|
| **A** group policy | MP excluded, HS quarantined, PostSurgery NaN on every ordinal scale, exact counts (307/277/30/150), primary ∪ sensitivity partitions the cohort, `apply_mask` filters by Sample ID; every validation case survives `primary_mask` (A6); the cross-cohort genus basis is fixed before the CLR, with prevalence applied to the training cohort only (A8/A9); the opt-in severity scale cannot touch the primary encoding (A10) |
| **B** leakage | default confounders carry no tumour column; `include_tumor_location=True` raises on the full cohort / Healthy-only / PostSurgery-only / CRC+1-Healthy; succeeds CRC-only with exact site counts; `Tumor.Location.Site` is NaN for every non-cancer sample |
| **C** statistics | `MeanDiff`/`CliffsD` sign conventions; `clr_lm_da` recovers a planted effect with the right sign, is calibrated on null data, refuses a too-small arm, drops constant covariates; progression FDR corrects each marginal test once; metabolite normalisation is blind to non-detects outside its reference set (C6); `partial_spearman_grid` reproduces pingouin exactly (C7); transfer features harmonise detection depth and share one CLR delta (C8); replication keeps "not measurable" apart from "not replicated" and flags contradictions (C9) |
| **D** biology | `POLYAMINE_EC` lists producing reactions only; the pathway set excludes the universal precursors; proxies are pre-specified and non-overlapping; Kim/Sinha lookups find the expected direct+proxy sets; AGORA2 evidence is disabled and no notebook calls its helpers (D5); the curated EC map survives for re-enabling (D6) |
| **E** hygiene | no notebook uses `SMOTE`, `TumorLoc_`, `Precursor`, or `'Adenoma'`/`'MP'` as a code-level group label; no notebook asserts the old HS-exclusion invariant; this file documents the policy; no active generator (E4); no blanket warning suppression (E5); no notebook-local genus/CLR helper (E6); the cohort counts quoted here match the data (E7); every output named here is declared in `output_contract.py` or listed as not produced (E8); every `utils` function an A–D check guards is reached by an executed notebook (E9) |

**Post-run output contract (`output_contract.py`).** `checks.py` runs before anything
executes; the output contract runs after each notebook. Every output a notebook declares
must exist, be written by that execution, and — for a CSV — carry a header row, the
required columns, enough rows and any required values (e.g. `differential_abundance_v2.csv`
must hold all three grouped comparisons for both assays); a pickle must hold its required
keys. Numbers this document quotes from results (NB05's shared-genus counts and transfer
statistics) are re-read from the run's tables and must match at the quoted precision. A
violation turns the notebook's OK into FAIL, and the same check gates reuse of an earlier
result. An output may be declared empty-able or optional only with a stated reason. Translate
Methods text only from outputs this check has confirmed.

---

## 11. Run order

```bash
python checks.py                       # must pass
python run_pipeline.py                 # 01, 01b, 02, 03, 05, 06, 08, 09; reuses valid results
python run_pipeline.py --smoke         # fast shake-out pass -- NOT publication output
python run_pipeline.py --dry-run       # shows every reuse decision and its reason
# NB10 (MICOM) is not runnable in the pinned environment (section 8).
```

Fatal if they fail: 01, 02, 03, 06, 08, 09. Non-fatal: 01b, 05. A notebook that exits
cleanly but breaks its output contract (§10) has failed. Executed notebooks, HTML,
logs and `provenance.json` land in `Results/runs/<run_id>/`; tracked notebooks are never
modified. Outputs replaced by a run are kept in `Results/_archive/<run_id>/` (newest three).
Status and peak memory per notebook: `Results/logs/summary_v2.txt`.
NB03 live progress: `Results/logs/nb03_progress_v2.txt`.
