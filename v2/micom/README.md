# `v2/micom` — the medium the MICOM batch was blocked on

CONTRACTS.md §8 put the 277-sample MICOM batch on hold on 2026-10-23 behind two
decisions, measured on 12 samples with the batch's exact settings:

1. **Solver.** 5 of 12 samples end infeasible at the pFBA step under the open-source
   HiGHS + OSQP hybrid. Looser tolerances did not help; the same samples solve with
   `strategy="none"`.
2. **Diet.** `western_diet_gut.csv` supplies putrescine and spermidine at
   1.0 mmol/gDW/h, so every community takes them up and none secretes them. Removing
   both made 1 of 4 test samples show secretion — and 3 of 4 become infeasible.

`diet_screen.py` settles both. The diet question it answers on **single-species AGORA2
models**, where a medium question is a linear program that takes 0.03 s rather than the
8 minutes a community model takes to answer several questions at once — 327 models, about
four minutes, no HPC. The solver question needs real communities and gets them under
`--community N`, which is findings 4 and 5 below and takes about an hour for 8 samples.

**The answer, for a reader who wants it before the reasoning.** The diet had three defects
and they are fixed. They were not the cause of the infeasibility: the same single sample
fails under both media, with a lumen imbalance seven orders of magnitude over tolerance
and every one of its 45 taxa at zero growth, which is a numerical failure and not a
biological one. **The solver decision is unchanged and still open.** The diet fix matters
for a different reason than the one it was undertaken for — it changes which taxa are
observed to produce putrescine and spermidine, from 1 sample to 7 and from 0 to 3, while
leaving every compound the diet does not supply untouched.

## Eight findings, in order of how much they change

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

**What fixing it is worth, measured (added 2026-11-08).** Nothing, in growth terms. On
the 8 community models of finding 4, restoring glucose moves community growth by a median
of **0.03%** and at most 2.4% (`out/community_per_sample.tsv`). So the bug is real, the
id is wrong, and the fix belongs in the file — but the communities were not
carbon-limited: this medium has 166 other lines and the pFBA optimum simply routed around
the missing sugar. Any claim that the missing glucose was *costing* these models anything
is unsupported, and the first version of this section implied it by calling glucose "the
single largest carbon source" without saying what its absence did. The honest statement is
narrower: a nutrient the file believed it was supplying was not being supplied, which
would have mattered for any question about carbon source, and did not matter for growth.

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

### 6. Would dropping the failures bias the result? Bounded, and smaller than I claimed

`dropout_risk.py` exists because the alternative to buying a QP is to run the batch and
drop whatever fails, and that is only safe if the dropout is independent of stage. I had
warned that it might not be. The measurement is more nuanced than the warning.

**The premise of the worry is confirmed.** Taxon count — the most obvious driver of how
hard a community QP is — tracks the stage axis across all 277 primary samples, and
survives adjustment for age, BMI, Brinkman index, alcohol and sex:

| driver | ρ (stage) | p | q | adjusted ρ | adjusted q |
|---|---|---|---|---|---|
| **taxa in the model** | **+0.145** | 0.016 | 0.080 | **+0.169** | **0.026** |
| Shannon | +0.111 | 0.065 | 0.161 | +0.136 | 0.062 |
| evenness | +0.035 | 0.57 | 0.94 | +0.057 | 0.58 |
| AGORA2 coverage | −0.003 | 0.96 | 0.96 | −0.031 | 0.61 |
| max abundance | −0.009 | 0.88 | 0.96 | −0.037 | 0.61 |

Median modelled taxa run 47 (Healthy) → 49 (Early_CRC) → 49.5 (Advanced_CRC). So a
dropout driven by community size *would* be stage-linked. One BH family over five
drivers on the three-level axis; the adenoma arm is outside the batch's 277 entirely.

**The consequence is inside the noise.** Take the worst case — failure perfectly monotone
in a driver, so the 115 hardest samples are exactly the ones lost — and compare it to
dropping 115 at random, 2,000 times, which is what "safe to drop" actually means:

| scenario | Healthy | Early_CRC | Advanced_CRC | Adv−Healthy taxa gap |
|---|---|---|---|---|
| all 277 | 45.8% | 34.7% | 19.5% | 2.5 |
| **random 115 dropped, 95% band** | **40.7–50.6%** | **29.6–39.5%** | **15.4–23.5%** | **−1.0 to 6.0** |
| highest taxon count dropped | 42.0% | 35.8% | 22.2% | 3.0 |
| highest Shannon dropped | 41.4% | 35.8% | 22.8% | 4.0 |
| highest evenness dropped | 43.8% | 35.8% | 20.4% | 4.0 |

**Every driver-led scenario, in both directions, falls inside the random band.** The
worst stage-composition shift a richness-linked dropout can produce in this cohort is
smaller than the ordinary sampling variation of losing 115 samples at all. (A single
seeded random draw moved Healthy by 4.8 points on its own, which is why the null here is
2,000 draws and not one — comparing against one draw would have invented a difference.)

**And the one failure we have points away from size.** Sample `10091` is at the *25th*
percentile of the 8 for both reactions (74,676) and taxa (45). The two largest models in
the draw — 113,354 and 110,510 reactions, 61 and 60 taxa — both solved cleanly. Its only
distinguishing feature is evenness, where it is the highest of the 8, and evenness is the
driver that tracks stage *least* (ρ 0.035, q 0.94). n=1, so this is an indication and not
a characterisation, but it points at a stage-independent driver rather than a
stage-linked one.

**So the honest verdict: the bias argument for buying a solver is weak.** What is left is
a power argument, which is a different claim and a smaller one. The worst case retains
68 Healthy / 58 Early / 36 Advanced against 127 / 96 / 54 — and a per-taxon secretion
analysis with BH over hundreds of taxa feels a 40% loss. Against a free academic licence
that is still an easy trade, but it should be made for the right reason: **the dropout
costs power, not validity.** Any write-up that justifies the solver choice by appealing to
confounding would be overstating this evidence.

Two limits worth keeping in view. The point estimates do move consistently (Healthy down,
Advanced up) under both richness and Shannon, so "inside the band" means undetectable at
this cohort size, not zero. And whether the quantity MICOM *reports* scales with taxon
count is untested: on the 7 valid samples, community growth, putrescine and spermidine
production all give p > 0.4, which at n=7 is close to no information.

### 7. The solver, without Gurobi: OSQP, repaired, then certified (2026-11-11)

No Gurobi licence was obtainable, and the project carries no commercial dependency.
`highs_qp.py` replaces micom's QP step in-process (`highs_qp.install("repair")`):

1. OSQP solves the cooperative-tradeoff QP as before. **It stopped at its 20,000-iteration
   limit on all 8 samples, including the 7 that `solution_quality` called valid.** That
   check certifies feasibility, not optimality, and it was the only check we had.
2. One HiGHS LP *repairs* OSQP's point: the exactly feasible flux state whose taxon growth
   rates are nearest (L1) to OSQP's. The pFBA step then always has a feasible point.
3. Wolfe's minimum-norm-point method (HiGHS LPs only) polishes it until the **certified**
   relative optimality gap (Frank-Wolfe bound, an upper bound and not an estimate) is
   <= 1%. Each sample's certificate is recorded.

Rejected on the way, each measured on sample 10723 (74,813 columns x 33,155 rows):
HiGHS's own QP solver ("Unbounded" on a QP bounded below by 0; >30 min), Kelley cutting
planes (stalled on degenerate cuts), Wolfe from scratch (~2.5 min per iteration), and
OSQP with more iterations (>10 min without converging).

`solver_check.py` re-solved the 8 `--community` samples, against criteria fixed before it ran:

| | criterion | result | |
|---|---|---|---|
| P1 | every sample valid | 8/8, residual ~0; sample 10091 from 0 to 45 taxa growing | PASS |
| P2 | no meaningful change where OSQP was valid | growth median 0.6%, community growth within 7e-4, **secretion calls identical in 79% (needed 90%)** | FAIL |
| P3 | certified gap <= 5% | worst 3.1%, 7/8 below 1% | PASS |
| P4 | median <= 5 min per sample | **5.5 min**; the two ~60-taxon samples took 36 and 57 min | FAIL |

Every secretion disagreement is in 10092 and 10334, the two largest communities, where
the repair had to move OSQP's "valid" growth rates by up to 11.3 -- OSQP was far from
optimal there, and the new answers are certified within 3.1% and 0.65%. In 10092 the
putrescine producer changes from *P. vulgatus* (OSQP) to *E. coli* (certified): **producer
identity, the paper's output, depends on solver accuracy in large communities.** P2 still
counts as failed; it is not re-scored.

### 8. The 277-sample batch: class M, and what it does and does not change (2026-11-13)

Dardel jobs 27477263 and 27477264 ran each sample as its own task. The task built the
community model (`micom_batch.py build-one`), then ran the certified cooperative tradeoff
at alpha 0.5 and the forced-secretion FVA (`forced_secretion.py`), in separate processes.
Before that, a 3-sample pilot (27477054/55) matched the laptop's memory predictions to
within 5%. Two earlier build jobs were lost to configuration errors, about 40 core-h:
the `.qza` was unzipped into the memory-backed `/tmp` of every worker, and then the
workers thrashed at the 20 GB ceiling.

| | |
|---|---|
| samples solved | 276 / 277 (10212: the pFBA step did not reach optimal) |
| certified gap <= 1% (the class M requirement, CONTRACTS 5d) | **137** |
| certified gap <= 5% (non-voting sensitivity) | 197 |
| (species, target) pairs assessable | 1,528 (1,741 at 5%) |
| **class M supported** (forced in >= 20% and >= 10 certified samples; base rate <= 0.2) | **8** (9 at 5%) |
| billed | 412 core-h against an estimate of 307 (tasks ran longer: median 21 min) |

The 8 pairs:
- **histidine**: *Anaerobutyricum hallii*, *Eubacterium_G ventriosum*, *Ruminococcus_C callidus*, each forced in 100% of samples;
- **tyrosine**: *F. prausnitzii_D* and *_G*, 99%;
- **ornithine**: *Enterocloster* sp001517625, 95%;
- **spermidine**: *E. bolteae*, 27%;
- **valine**: *Phascolarctobacterium_A succinatutens*, 31%.

The 1% certification rule costs almost nothing: 8 pairs against 9 at 5%.

**Class M moves no tier.** A pair enters the evidence table only through a significant
partial correlation in Yachida (NB08), and none of the 8 has one. Four were tested and
are null (rho -0.12 to +0.13, q 0.14-0.62). The other four involve species outside NB02's
prevalence-filtered set, so they were never tested. **The organisms the community model
says must export a target are not the organisms whose abundance tracks it.** That is a
result in its own right, not a gap to be closed.
- Some reasons are mechanistic. Amino acids are exported by many taxa, and faecal
  concentration integrates host, diet and microbial supply.
- Some are statistical. A forced export at fixed abundance need not co-vary with that
  taxon's abundance across people.

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

### 4. The remaining infeasibility is the solver, and the diet does not touch it

`--community 8` builds 8 real community models (151 distinct taxa, 14 min) and runs
micom's own `grow` on each twice — `tradeoff=0.5`, `strategy="pFBA"`, the batch's
settings — once under each medium, 26 min a side. **The prior stated here before it ran
was wrong, and this is the refutation.**

| | `western_diet_gut` | `gut_diet_v2` |
|---|---|---|
| returned a solution | 8/8 | 8/8 |
| passed `solution_quality` | 7/8 | 7/8 |
| median community growth | 2.0516 /h | 2.0502 /h |

The margin is identical, and a margin cannot be read: 7-of-8 twice is equally consistent
with "the same sample fails both times" and with "a different sample fails each time". The
**paired** join settles it, and being paired it says something at n=8 that the rate never
could — a marginal 7/8 has a Wilson interval from 0.53 to 0.98 and excludes nothing.
**Validity agrees on 8 of 8 samples.** The same sample, `10091`, fails under both media,
and every other sample passes under both.

That failure is not biological:

- **0 of its 45 taxa grow**, at exactly zero — while all 8 other samples have *every*
  taxon growing (61/61, 60/60, 57/57 …). The single-model screen found 315 of 327 species
  grow on this medium alone, so a 45-taxon community on the same medium in which nothing
  grows is not a statement about biology.
- **None of its 45 taxa is one of the 12** that cannot grow alone, so finding 3 does not
  explain it either.
- The returned point **is not mass-balanced**: `max_balance_residual` is 18.3 under the
  new medium and 62.5 under the old, against a tolerance of 1e-6 — seven orders of
  magnitude over. micom returns infeasible points without raising, and this is one.
- Its largest lumen fluxes are through AGORA2's drug-metabolism exchanges
  (`EX_bilr_M12_m`, `EX_chlphncl_glc_m`, `EX_1hibup_S_m`, `EX_Tn_antigen_m`) at 5–18
  units, with 255 compounds being taken up and zero growth produced. That is a solver in
  a garbage corner of the polytope, which is the characteristic numerical failure of
  HiGHS + OSQP on a badly scaled QP.

**So the solver decision stands as it was on 2026-10-23, and the diet was not a
confounder of it.** The glucose fix is not the cure, for the reason finding 1 now records:
the communities were never carbon-limited, so the growth maximum was never near zero, so
the mechanism the prior proposed did not exist. A Gurobi academic licence (or another QP
that is not OSQP) remains the open item for the batch.

What this draw does **not** establish is the *rate*. 1 of 8 here is 12.5%; the 2026-10-23
measurement was 5 of 12, 42%. With n=8 the Wilson interval on 7/8 is [0.53, 0.98], which
contains the old rate comfortably — these two numbers are not in conflict and this run
does not revise the 42%. It identifies the failure mode, not its frequency.

### 5. The diet change does not alter growth. It alters the answer.

The point of cutting putrescine and spermidine was never growth — it was to stop flooding
the lumen so that production becomes observable. Over the **7 samples valid under both
media** (`out/community_secretion_valid.tsv`):

| compound | in the diet? | samples with a secretor | distinct taxa |
|---|---|---|---|
| putrescine | yes, cut 4,974× | **1 → 7** | **1 → 4** |
| spermidine | yes, cut 319× | **0 → 3** | **0 → 2** |
| ornithine | no | 7 → 7 | 27 → 28 |
| GABA (`4abut`) | no | 2 → 2 | 3 → 3 |
| agmatine | no | 0 → 0 | 0 → 0 |
| cadaverine (`15dap`) | no | 0 → 0 | 0 → 0 |

The last four rows are the control, and they are why the first two can be believed: the
only compounds that moved are the two whose supply changed, and nothing was lost anywhere.
Under the old medium this cohort would have been reported as showing **no putrescine
production in 6 of 7 samples and no spermidine production in any of them** — not because
nobody makes them, but because the medium supplied 1.0 mmol/gDW/h of each and every taxon
was a net consumer at that concentration. That is the result the batch exists to produce,
and the old medium would have produced the wrong one.

Two cautions, both load-bearing. First, **the signal is detectable, not large**: in 6 of
the 7 samples putrescine production is a single taxon at 0.0002–0.11 mmol/gDW/h, and only
`10335` carries a big flux (16.5) that was already there under the old medium. Spermidine
is one taxon in each of 3 samples. Second, **n = 7 samples**; the taxa named here
(*E. coli*, *Parasutterella excrementihominis* and *Phocaeicola vulgatus* for putrescine;
*Enterocloster bolteae* and *E. coli* for spermidine) are a demonstration that the
measurement now works, not a finding about the cohort. The full 277-sample batch is what
would make them one, and it is still behind the solver.

A third caution about reading the raw counts: `community_check.tsv` reports 14 → 23
spermidine secretors across all 8 samples, which is larger than the 0 → 2 above. Every
spermidine secretor under the old medium sits in sample `10091` — the invalid one. A
secretor count read off a point with an 18-unit mass imbalance is not evidence of
anything, which is why `community_secretion_valid.tsv` exists and is the table to quote.

## Files

| | |
|---|---|
| `diet_screen.py` | the whole thing: name audit, auxotrophy screen, dose–response, community demand, the new file, and the optional community check |
| `gut_diet_v2.csv` | the reformulated medium, with its reasoning in the header |
| `out/diet_auxotrophy_screen.tsv` | per species: growth on the diet and without each target compound |
| `out/cohort_species_agora2.tsv` | the 1,562 cohort species, which 327 are in AGORA2, and their abundances |
| `out/diet_gaps_nongrowing.tsv` | the 12 species that cannot grow, and what each is missing |
| `out/diet_screen_provenance.json` | every number above, including the dose–response curves |
| `out/community_check.tsv` | written by `--community`: both media, side by side. Counts over **all** samples, so read it with the caution at the end of finding 5 |
| `out/community_per_sample.tsv` | the paired per-sample table: validity, growth and lumen imbalance under each medium. This is what settles finding 4 |
| `out/community_secretion_valid.tsv` | secretion restricted to samples valid under both media, with the taxon names. **The table to quote** |
| `highs_qp.py` | finding 7: the Gurobi-free QP step (OSQP, exact HiGHS-LP repair, Wolfe polish to a certified gap) |
| `solver_check.py` | finding 7: the new solver against OSQP on the 8 `--community` samples, criteria P1-P4 |
| `out/solver_check_per_sample.tsv` | per sample: validity, growth, certified gap, OSQP status, wall time |
| `out/solver_check_growth.tsv` / `out/solver_check_secretion.tsv` | per-taxon growth and per-compound secretor sets, OSQP vs new |
| `dropout_risk.py` | finding 6: whether losing ~40% of samples to the solver would bias the progression analysis |
| `out/dropout_risk.tsv` | per sample: what its community model would contain, plus stage |
| `out/dropout_risk_association.tsv` | each candidate driver of solver difficulty against the stage axis, adjusted and not |
| `out/dropout_risk_scenarios.tsv` | the worst-case dropout against a 2,000-draw random null |

## Reproducing

```bash
python v2/micom/diet_screen.py                 # ~4 min, no HPC
python v2/micom/diet_screen.py --community 8   # adds the community comparison, ~70 min
python v2/micom/dropout_risk.py                # finding 6, ~1 min, no solver needed
```

`--community` caches: `out/community/models/` holds the built community models and
`out/community/grow_*.pkl` the two solves, so a re-run re-does only the analysis in
seconds. Both are gitignored (hundreds of MB) and rebuild from the flag; the three small
tables they produce are committed. Delete `out/community/` to force a genuine rebuild.

Inputs are `Results/intermediate/preprocessed_data_v2.pkl` (the abundance table NB10
uses), `Data/micom/agora2_extracted/` (AGORA2 pre-extracted by NB10 cell 5) and
`Data/micom/western_diet_gut.csv`. The species set is derived with NB10's own
`gtdb_binomial` rule and reproduces its reported **64.6% median per-sample abundance
coverage** exactly, which is the check that this screen characterises the models the
batch would actually build.
