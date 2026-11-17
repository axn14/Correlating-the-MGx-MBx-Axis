# Microbial sources of CRC-associated faecal metabolites (v2.0)

This repository asks which gut bacteria produce the faecal metabolites that distinguish
colorectal cancer (CRC) from health. It uses public paired metagenome and metabolome data
and an evidence standard that was written down before its results existed.

Version 2.0 replaces the thesis pipeline (tag `v1.0-thesis`). The two versions reach
different conclusions. [CHANGELOG.md](CHANGELOG.md) explains every change, and the
[wiki](../../wiki) page "What changed in v2.0" explains why the earlier results should not
be used.

## What the analysis found

In 277 paired samples from Yachida et al. (2019), 20 metabolites separated healthy people
from early CRC and 18 from advanced CRC, against 0 and 3 species. Nothing separated early
from advanced CRC, or adenoma from health.

We tested 4,293 species–metabolite pairs on 47 target metabolites against five independent
classes of evidence (association, reference genome, gene reads in the patients' own
metagenomes, community metabolic models, and replication in a second cohort). A class counts
once however many tools it contains. 82 pairs met the Supported tier or better, and 1 met
every criterion (*Veillonella atypica* and putrescine).

Shuffling each class within its metabolite produces about 45 such pairs by chance, so about
half of that tier is expected without any real agreement. For the 30 CRC-associated
metabolites, Supported links did not exceed chance. The links the standard does recover,
such as *Ruminococcus gnavus* and tryptamine, are cases where the gene evidence points to a
single species, and they match published biochemistry.

The short version: with current public data, the faecal metabolites that mark CRC cannot be
traced to a microbial producer above chance. The [wiki](../../wiki) explains why, target by
target.

## Repository layout

```
v2/                  the pipeline
  CONTRACTS.md       every rule, dated; amendments were committed before the results they govern
  checks.py          56 pre-flight assertions; the pipeline will not run if one fails
  run_pipeline.py    notebooks NB01 to NB09 in order, with an output contract
  01_ ... 10_*.ipynb analysis notebooks
  targets/           the 47 targets and their export, origin and confounder gates
  genes/             class G: HUMAnN gene families attributed to species
  micom/             class M: community models, certified solver, forced secretion
  replication/       class X: Erawijantari 2020 replication; Nofinska 2026 overlap map
  calibration/       the tier-rule null and replication power analysis
  assemble_classes.py  tiers every pair by evidence class
  dardel/            batch scripts for the PDC Dardel cluster
```

## Running it

Python 3.12. Install `v2/requirements.lock`; the community models need
`v2/requirements-micom.txt` in addition. Data are the public tables of the curated
microbiome–metabolome collection (Muller et al. 2022) and curatedMetagenomicData, placed
under `Data/` as described in `v2/CONTRACTS.md`.

```
cd v2
python checks.py                 # must pass before anything runs
python run_pipeline.py           # NB01 to NB09
python targets/build_registry.py --results <run>
python replication/erawijantari_species.py --results <run>
python assemble_classes.py --results <run>
python calibration/class_null.py --results <run>
python calibration/replication_power.py --results <run>
```

The community-model batch (`v2/micom/`) used about 410 core-hours on a cluster.

## Status

The manuscript is in preparation. Until it is published, treat the results here as
preliminary. Earlier summaries of this work, including the 2026 MSc thesis, reported
"Trinity" producer pairs that do not hold up; see the changelog.

## Citing

Please cite the v2.0 release (DOI to follow) and the original data publications: Yachida et
al. 2019, Erawijantari et al. 2020, Sinha et al. 2016 and Kim et al. 2020.
