# Running the v2 pipeline on Dardel

Dardel is the PDC/KTH cluster (HPE Cray EX, AMD EPYC). This directory holds the four
pieces needed to run the pipeline there and evaluate the results on the laptop,
without moving the heavy parts of the tree in either direction.

## Why this is cheap

The weight is all on the input side, and the pipeline needs almost none of it:

| | size | needed by the pipeline? |
|---|---|---|
| `Data/micom` (prebuilt AGORA2 models) | 53.5 GB, 2,756 files | **no** — NB10 only, which `run_pipeline.py` does not run |
| `Data/*.tsv` (cohort tables) | 134 MB | yes |
| `Data/reference/` (AGORA2 EC index) | 30 MB | yes |
| code, as a git bundle | ~20 MB | yes |
| **pulled back per run** | **~174 MB** | tables 104 MB, figures 54 MB, run dir 15 MB |

So staging is 164 MB, and a complete evaluable copy of a run is ~174 MB and ~360
files. The rendered notebook HTML in the run directory opens in a browser with no
Python environment on the laptop at all.

## Which shell to run the laptop-side scripts in

Either MobaXterm's **local shell** tab or WSL. MobaXterm is usually the better
choice, because it bundles `git`, `rsync` and `ssh` (under
`Documents\MobaXterm\slash\mx86_64b\bin`) *and* already holds the PDC key, so nothing
has to be copied anywhere. Its drive root differs from WSL's — `/drives/d/...`
against `/mnt/d/...` — which `_paths.sh` resolves, so one `dardel.env` works in both.

Set `PDC_SSH_KEY` in `dardel.env` if your key is not at a default identity name.
MobaXterm's is `~/.ssh/id-ed25519-pdc`, which `ssh` will not find on its own.

In WSL the same key has to be copied into `~/.ssh` with `chmod 600` first: `ssh`
refuses a key that is world-readable, and everything under `/mnt/c` is.

## Sequence

1. **Once, on the laptop.** Copy `dardel.env.example` to `dardel.env` and fill it in.
   `dardel.env` is gitignored: it holds your username, allocation id and absolute
   paths. Confirm `PDC_STORAGE` with `projinfo` on Dardel first — the PDC docs are
   inconsistent between `/cfs/klemming/projects/snic/<project>` and
   `/cfs/klemming/<project>/<user>`, and only your own `projinfo` output settles it.

2. **Stage, from MobaXterm's local shell (or WSL).**
   ```
   cd /drives/d/KI\ Research/v2/dardel && bash 02_stage_inputs.sh
   ```
   Pushes the repository as a git bundle and unpacks it into a working repository on
   Dardel, then the cohort tables and the AGORA2 index. Add `--with-micom` later, only
   when the MICOM work starts.

3. **Build the environment, once, on a Dardel login node.**
   ```
   PDC_STORAGE=/cfs/klemming/projects/snic/<project> bash 01_setup_env.sh
   ```

4. **Submit.**
   ```
   sbatch -A <project> run_pipeline.sbatch --dry-run   # the reuse plan, seconds
   sbatch -A <project> run_pipeline.sbatch --smoke     # shake-out, not publication output
   sbatch -A <project> run_pipeline.sbatch             # full run
   ```

5. **Pull, from WSL on the laptop.**
   ```
   bash pull_results.sh --list
   bash pull_results.sh            # the newest run
   ```

## The things that will bite

**`Results/` is not append-only.** `run_pipeline.py:283` moves every output directory
into `Results/_archive/<run_id>` when a run starts, and `prune_archives` keeps only
the newest three, by `rmtree`. So:

- `02_stage_inputs.sh` never pushes `Results/`. If you uploaded the frozen baseline,
  the next run would relocate it and three runs later delete it.
- `pull_results.sh` is one-way and additive, into a per-run directory, and refuses a
  destination inside the project's own `Results/`. Never sync the two trees.

**Klemming does not keep things for you.** Project storage is not backed up and is
deleted three months after the project ends. Scratch removes anything not modified
for 30 days. Nothing you intend to cite should exist only on Dardel.

**Home is 25 GB and 100,000 files.** The equivalent `.venv` on the laptop is 43,256
files, so an environment in home would consume a large share of the inode quota.
`01_setup_env.sh` puts the venv in project storage instead, which has the room — the
project's file quota is 1,000,000 with roughly 340,000 free.

**The environment is `cray-python/3.12.12` plus a venv, not conda.** Both were
checked on the machine: `miniconda3/26.3.2-2-cpeGNU-26.03` exists under `PDC/26.03`,
and cray-python gives Python 3.12.12 with venv and pip 25.3. venv wins because
`requirements.lock` is a `pip freeze` from Python 3.12.10, so every wheel carries a
`cp312` tag and installs unchanged; because pip wheels bundle OpenBLAS on Linux just
as they did on Windows, keeping numpy and scipy on the same BLAS family as the frozen
baseline where conda-forge would substitute its own; and because it needs no channel
resolution or `.condarc` relocation. The venv is created without
`--system-site-packages` on purpose: the module ships numpy 2.3.5 and the pipeline
pins 2.3.1.

**Two core budgets, not one, because they differ in what they can change.**
`run_pipeline.py` reads both:

- `KI_N_CORES` — worker count (loky/joblib, `KI_N_JOBS`). `PROFILE_KEYS` excludes
  `KI_N_JOBS` because every estimator pins `n_jobs=1` and `random_state`, so a single
  fit does not depend on it. NB03's permutation loop does: it stops early on the
  exceedance count, checked at joblib batch boundaries, and those boundaries move with
  the worker count. Measured against the baseline, `N_perm` for one target differed by
  557 draws. `run_pipeline.sbatch` sets it from `SLURM_CPUS_PER_TASK`, which is 10 for
  a `-c 8 --mem=8G` job, not 8.
- `KI_BLAS_THREADS` — BLAS threads for the non-joblib notebooks. Defaults to 8, the
  value the frozen baseline used, and `run_pipeline.sbatch` deliberately leaves it
  alone. Raising it reorders floating-point reductions inside BLAS, which can move
  the last digits of every number downstream of NB01 — and reuse keys do not track
  it, so a change needs `--fresh` and a re-verification. The runner prints a loud
  warning when it is not at baseline.

**A Dardel run is a new run, not a replacement — and measurably so.** The manuscript's
271 assertions are pinned to the frozen local run at `f14c7f9`. The first full Dardel
run used the same commit, the same `inputs afb1eebccde3`, the same pinned versions
down to the patch and the same bundled OpenBLAS 0.3.29, and still moved every one of
NB03's 36 AUCs, by up to 0.041; one target crossed the `q<0.05` gate and added two
SHAP figures. The cause is the CPU, not the port: OpenBLAS dispatches
microarchitecture-specific kernels (EPYC 7742 is Zen 2, the laptop is Zen 3+), and
nested CV amplifies last-digit differences into percent-level ones because choosing
hyperparameters is an argmax over near-tied fold scores. `R2_confounders`, a plain
ridge fit with no workers involved, differs at 2e-16 — that is the fingerprint.
CONTRACTS.md section 12 has the full account.

So compare, do not substitute:

```
python "/mnt/d/KI Research/v2/dardel/compare_to_baseline.py" "/mnt/d/KI Research-dardel/<run_id>"
```

It reports the gate flips first (the only differences that change what the chapter can
claim), then which quoted AUCs would have to be requoted, then the spread.

`writeup/verify_manuscript_numbers.py` is **not** the tool for this. It reads the
baseline's tables by a hard-coded absolute path, so `KI_PROJECT_ROOT` does not
redirect it: aimed at a pulled tree it checks the laptop's own tables, silently skips
the sections it resolves through `utils`, and still prints "no mismatches".

Failures there are information about the two runs differing, not something to fix by
overwriting files.

## Partition choice

`-p shared`, not `-p main`. `main` allocates whole nodes and bills 128 cores whatever
you use, with a 24-hour ceiling; this job is one sequential chain of notebooks, so a
five-hour run there would cost roughly 640 core-hours against about 160 on `shared`,
which allocates by core and allows 7 days. Memory in `shared` is proportional to the
core count (roughly 0.85 GB per core), so `--mem=24G` sits inside what `-c 32`
already grants. `seff <jobid>` afterwards gives the real peak.

`-p memory` (512 GB / 1 TB / 2 TB nodes) is for the MICOM and GTDB-Tk work later, not
for this pipeline.

## What this does not do

**The NB03 job array.** §4 of `../../writeup/dardel_plan.md` proposes splitting the
benchmark over 36 SLURM tasks, one per target. That is the right next step — it is
near-linear, about fifteen minutes wall-clock, and because each target stays
single-threaded it is *bit-identical* to the frozen run rather than drifting, so the
parallel route is also the reproducible one. It then makes the deferred budgets
affordable: the `b_max` = 1,000 permutation ceiling that currently sets the
attainable *q* floor, and the bootstrap replicates for the rank confidence intervals.

It is not here because it needs a contained refactor of the benchmark loop into a
per-target entry point writing one row and one SHAP table each, plus an extension of
`output_contract.py` to cover the per-target intermediates. That belongs in its own
change, with its own checks, not bundled into a first port.

Expect the port alone to take the run from about five hours to roughly ninety
minutes. Dardel's per-core clock is not better than a modern laptop's — these are
Zen 2 parts. What it buys is core count, memory, a 7-day wall clock, and the ability
to run the read programme at all.

**Apptainer.** It is installed `nosuid` on Dardel, so `.sif` images cannot be run,
only sandbox directories, and recipes cannot be built on the cluster. For a pure
Python stack with a working lock file, cray-python plus a venv is the shorter path.

**A gotcha that cost a round trip.** `ml` is a shell function, so piping it into
`head` or `grep` runs it in a subshell and the environment changes are silently
discarded — the probe reported Python 2.7.18 while claiming the module had loaded.
Both scripts now say so where they call it.

## Sources

PDC documentation, read 2026-11-07:
[Klemming storage](https://support.pdc.kth.se/doc/data_management/klemming/) (paths,
quotas, purge policy),
[job scheduling](https://support.pdc.kth.se/doc/run_jobs/job_scheduling/)
(partitions, node types, limits),
[Python](https://support.pdc.kth.se/doc/applications/python/) (miniconda3,
`.condarc`, when to fall back to `cray-python`),
[Apptainer](https://support.pdc.kth.se/doc/applications/apptainer/),
[SSH login](https://support.pdc.kth.se/doc/login/ssh_login/) (key registration via
the login portal and SUPR),
[data management](https://support.pdc.kth.se/doc/data_management/data_management/)
(transfer node).

Module versions move: confirm `PDC/...` and `miniconda3/...` with `ml spider
miniconda3` before the first run. Both scripts take them as overridable variables.
