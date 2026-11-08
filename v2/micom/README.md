# `v2/micom` — the medium the MICOM batch was blocked on

CONTRACTS.md §8 put the 277-sample MICOM batch on hold on 2026-10-23 behind two
decisions, measured on 12 samples with the batch's exact settings:

1. **Solver.** 5 of 12 samples end infeasible at the pFBA step under the open-source
   HiGHS + OSQP hybrid. Looser tolerances did not help; the same samples solve with
   `strategy="none"`.
2. **Diet.** `western_diet_gut.csv` supplies putrescine and spermidine at
   1.0 mmol/gDW/h, so every community takes them up and none secretes them. Removing
   both made 1 of 4 test samples show secretion — and 3 of 4 become infeasible.

`diet_screen.py` settles the second and bounds the first. It works on **single-species
AGORA2 models**, where a medium question is a linear program that takes 0.03 s, rather
than on community models, where it takes 8 minutes and answers several questions at
once. 327 models, about four minutes.

## Three findings, in order of how much they change

### 1. The medium has been supplying no glucose

`western_diet_gut.csv` lists `EX_glc_e`. Every AGORA2 model spells glucose
`EX_glc_D(e)`. Across the **whole AGORA2 release, 1,326 of 1,439 models carry
`EX_glc_D(e)` and exactly 0 carry `EX_glc(e)`** (308 of our 327 cohort models carry it).
The file does use `_D`/`_L` suffixes elsewhere — it lists both `EX_glu_D_e` and
`EX_glu_L_e` — so glucose specifically lost its stereochemistry somewhere.

A medium in micom and cobra is a dictionary lookup: a reaction id the model does not
have is not an error, not a warning and not a zero-flux entry. It is simply absent, and
the nutrient is silently not supplied. So the single largest carbon source in the file
has been reaching none of the communities, with nothing anywhere saying so. This is the
same failure shape as the `cadA`/cadmium-pump collision in `../reads/README.md` and the
EC-map blind spot in `../genes/README.md`: a wrong number, not an error message.

Three naming dialects are in play and none of them is wrong on its own: micom's own
bundled medium writes glucose `EX_glc__D_m` (BiGG, double underscore), this AGORA2
release writes `EX_glc_D(e)`, and the project's diet file writes `EX_glc_e`. That is why
`diet_screen.py` audits the diet against the model files on every run instead of
hard-coding a correction — swap the model database and the right id changes with it.

Seven further lines reach no cohort model either — `EX_arachd_e`, `EX_hdcea_e`,
`EX_lnlnca_e`, `EX_octa_e`, `EX_adpcbl_e`, `EX_fe3dcit_e`, `EX_amp_e`. Those are long
chain fatty acids, adenosylcobalamin, Fe(III)-dicitrate and AMP, which AGORA2 genuinely
does not exchange. They are left in place and named in the new file's header, because
that is AGORA2's coverage rather than an error in the diet.

### 2. Putrescine and spermidine are oversupplied by 300–5,000×

Withdrawing them outright is not an option, and the screen says why with names:

| withdrawn | cohort species that stop growing |
|---|---|
| putrescine | **64 of 315** |
| spermidine | **172 of 315** |
| both | **174 of 315** |

More than half the modelled community is a spermidine auxotroph. That is the whole
explanation of "3 of 4 samples became infeasible when we removed them" — it is a
property of the medium and the models, not of the solver.

But *how much* they need is the question nobody asked. Spermidine is a biomass precursor
in these models, so uptake is proportional to growth rate, and the curve saturates far
below the diet's value:

```
Bacteroides uniformis     0.0000  0.0014  0.0142  0.1424  0.8591  0.8591  0.8591
Phocaeicola dorei         0.0000  0.0025  0.0252  0.2051  0.2051  0.2051  0.2051
Bacteroides stercoris     0.0000  0.0013  0.0126  0.1259  0.2774  0.2774  0.2774
EX_spmd(e) supply         0       1e-05   1e-04   1e-03   1e-02   0.1     1.0
```

Taking the saturation point per species (99% of the growth the full diet gives) and
weighting by abundance — a taxon's draw on the lumen is its abundance times its own
per-gram flux, the unit rule `micom_tools` encodes and NB10 asserts — gives the
community-level demand the medium actually has to meet:

| | per-gram saturation | worst sample's community demand | diet | new |
|---|---|---|---|---|
| putrescine | up to 0.005974 | 0.0001005 | 1.0 | **0.000201** (4,974× less) |
| spermidine | up to 0.006042 | 0.001568 | 1.0 | **0.003137** (319× less) |

So the diet can be cut by two to three orders of magnitude with **no taxon losing any
growth at all**, which should make secretion visible against a lumen that is no longer
flooded. That is the reformulation, and it is a measurement rather than a preference.

An earlier version of this script bisected for "the smallest supply at which the model
grows at all" and returned 6e-05. A diet built on that number would have been feasible
and meaningless: it holds every auxotroph at about 1% of its growth rate. The quantity
that matters is where the curve flattens, not where it leaves zero.

### 3. Twelve species cannot grow on this medium at all

| species | what `micom.media.complete_medium` says it is missing |
|---|---|
| *Coprococcus* sp000154245 | `EX_nh4(e)` ammonium |
| *Ruminococcus_F champanellensis* | `EX_nh4(e)` ammonium |
| *Porphyromonas_A somerae* | `EX_lac_L(e)` L-lactate |
| *Allisonella pneumosintes* | `EX_mal_L(e)` L-malate |
| *Enterococcus_A avium* | `EX_mal_L(e)` L-malate |
| *Bacteroides_F pectinophilus* | `EX_ins(e)` inosine |
| *Turicibacter* sp001543345 | `EX_ins(e)`, `EX_glytyr(e)` |
| *Phascolarctobacterium_A succinatutens* | `EX_succ(e)`, `EX_anzp(e)` |
| *Akkermansia muciniphila* | `EX_glyasp(e)` Gly-Asp |
| *Acutalibacter timonensis* | `EX_glyasn(e)` Gly-Asn |
| *Alloscardovia omnicolens* | `EX_alaasp(e)` Ala-Asp |
| *Tractidigestivibacter* sp000752675 | `EX_no3(e)` nitrate |

**142 of 277 samples (51.3%) contain at least one.** They are not an infeasibility
cause — micom's `Community.cooperative_tradeoff` takes `min_growth=0.0` and constrains
a fraction of the *community* growth rate, not each taxon's own maximum, so a taxon that
cannot grow simply grows at zero. (An earlier draft of this file asserted the opposite;
micom's source says otherwise and the claim is withdrawn.) They are a validity problem:
a community model whose member contributes no growth and no metabolism is not modelling
that member.

The cost is small — those taxa carry a median 0.11% and at most 4.62% of a sample's
abundance — so `gut_diet_v2.csv` **does not** add these nutrients. Ammonium, L-lactate,
L-malate and succinate are all unarguably present in a colon and adding them is
defensible; the dipeptides and nitrate are more of a judgement call; and every addition
changes the medium for the other 315 species too. That is a scientific decision to take
deliberately, not a side effect of fixing a typo, so it is reported here and left open.
`out/diet_gaps_nongrowing.tsv` has the full list.

## The reformulated medium

`gut_diet_v2.csv` — `western_diet_gut.csv` with exactly three changes, each recorded in
its header:

- `EX_glc_e` → `EX_glc_D_e`, same flux
- putrescine 1.0 → 0.000201
- spermidine 1.0 → 0.003137

Nothing else is touched. Point NB10's `DIET_FILE` at it:

```python
DIET_FILE = V2 / 'micom' / 'gut_diet_v2.csv'      # was DATA_DIR/'micom'/'western_diet_gut.csv'
```

NB10 reads it with `pd.read_csv`, so the `#` header lines need `comment='#'` there, and
`BATCH_PARAMS['diet_sha256']` changes — which is the mechanism that correctly refuses to
mix old and new chunks.

## What this cannot settle

Whether the **remaining** infeasibility is the solver. Single-model LPs are not the
community QP: a sample can fail for reasons no single model shows, because the tradeoff
constraint couples every taxon through a shared lumen. What the screen bounds is the
medium's contribution. `diet_screen.py --community N` measures the rest directly —
it builds N community models once and runs micom's own `grow` on them under both media,
counting solutions returned, solutions that pass `micom_tools.solution_quality`, and
taxa secreting each compound.

The honest prior before that runs: the glucose fix is a real candidate cure for the pFBA
failures, because a carbon-starved community has a near-zero growth maximum and the
tradeoff QP is at its worst there. It is a hypothesis until the numbers come back.

## Files

| | |
|---|---|
| `diet_screen.py` | the whole thing: name audit, auxotrophy screen, dose–response, community demand, the new file, and the optional community check |
| `gut_diet_v2.csv` | the reformulated medium, with its reasoning in the header |
| `out/diet_auxotrophy_screen.tsv` | per species: growth on the diet and without each target compound |
| `out/cohort_species_agora2.tsv` | the 1,562 cohort species, which 327 are in AGORA2, and their abundances |
| `out/diet_gaps_nongrowing.tsv` | the 12 species that cannot grow, and what each is missing |
| `out/diet_screen_provenance.json` | every number above, including the dose–response curves |
| `out/community_check.tsv` | written by `--community`: both media, side by side |

## Reproducing

```bash
python v2/micom/diet_screen.py                 # ~4 min, no HPC
python v2/micom/diet_screen.py --community 8   # adds a real community comparison, hours
```

Inputs are `Results/intermediate/preprocessed_data_v2.pkl` (the abundance table NB10
uses), `Data/micom/agora2_extracted/` (AGORA2 pre-extracted by NB10 cell 5) and
`Data/micom/western_diet_gut.csv`. The species set is derived with NB10's own
`gtdb_binomial` rule and reproduces its reported **64.6% median per-sample abundance
coverage** exactly, which is the check that this screen characterises the models the
batch would actually build.
