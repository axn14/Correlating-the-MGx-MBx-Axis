# `v2/targets` — the target registry (CONTRACTS §5d-ii)

There is one row per trinity target. It carries the gates that decide which evidence
classes a target *can* receive. These are properties of the metabolite, fixed before any
attribution is scored.

```
python v2/targets/build_registry.py --results "D:/KI Research-dardel/<run_id>"
```

Outputs: `out/target_registry.tsv`, `out/target_registry_provenance.json`, and
`out/kegg_cache.json`. The cache holds the KEGG REST answers, so a re-run is offline and
reproducible.

## What it found on run `20261108T133844Z_1720a81` (2026-11-11)

**33 + 16 = 49.**
- 33 targets are `stage_associated`: robust DA in Yachida, Healthy vs Early and/or Advanced CRC.
- 16 are `a_priori_class`: in the pre-specified decarboxylation/polyamine class but not robust DA.
- None of the 16 a-priori targets is a progression target. They are biological controls,
  and every table now says so.

**The export gate removes 23 of 49 targets from MICOM.**
- MICOM's AGORA2 database has 1,440 models. Not one carries an exchange reaction for 23
  of the targets. The registry reads the models' own `kegg.compound` annotations, and
  `grep` on the raw files confirms it. Among them:
  - N-acetylglutamate: 265 of 327 cohort models touch it intracellularly, but none exports it.
  - N-acetylhistamine (C05135).
  - Agmatine, spermine, citrulline and homoserine.
- For comparison, putrescine is exchangeable in 1,171 models.
- For these 23, class M (forced secretion) is unavailable, and the registry flags
  `attribution_requires_lysis`.
- **Two targets hold 7 of the 8 pairs the old tiers called `Strong`: N-acetylglutamate (5) and N-acetylhistamine (2).**

**Origin flags.**

| Flag | Count |
|---|---|
| `host_and_microbial` | 33 |
| `host_dominant` | 5 (hydroxyproline, N,N-dimethylglycine, N1-acetylspermine, trimethyllysine, 1-methylhistamine) |
| `microbial_only` | 4 |
| `unannotated` | 5 (neither KEGG-human nor AGORA2 annotates them, so no call either way) |
| `unassessable` | 2 (no KEGG id) |

The host flag is direction-agnostic: a human enzyme acts on the compound, whether it
makes or consumes it. It is an annotation and does not correct anything.

**Confounder gate.** 8 targets are `confounder_unattributable` in this run: the confounders
alone match or beat the species model (`evidence_tools.confounder_unattributable`).

## Columns
- `label`: `stage_associated` or `a_priori_class`.
- `robust_early` / `robust_advanced`: the direction of NB02's MeanDiff (CRC minus Healthy) where robust.
- `auc_model`, `auc_confounders_adjusted`, `perm_q`, `confounder_unattributable`: from NB03.
- `agora2_models_exporting` / `agora2_models_touching`: cohort models, out of 327.
- `export_gate`, `attribution_requires_lysis`, `human_ec`, `origin_flag`, `note`.

## Sequencing depth (`depth_screen.py`)

Yachida's metadata has no read counts. As a depth proxy, the screen uses the gzipped FASTQ
size per sample from the ENA manifest (`v2/reads/manifest_347.tsv`). The decision rule was
written into the script before it ran: if depth tracks the stage axis at p < 0.05, add it
to NB03's confounders.

**Result:** Spearman ρ = −0.12, p = 0.046 (n = 277). The group test is not significant
(Kruskal p = 0.11). The median is 4.9 GB in Healthy and Stage 0, and 4.5 GB in Stage I–IV.

It is weak, but the rule says add it, so depth joins the confounder set at the Phase G
rerun. It matters most for class G: whether a gene is detected depends on depth, so the
gene layer (Phase D) carries depth as a covariate.
