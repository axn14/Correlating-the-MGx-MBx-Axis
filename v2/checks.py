"""Pre-flight self-checks for the v2.1 pipeline.

`run_pipeline.py` runs this FIRST and aborts if anything fails, so a run can never
produce results under a broken group definition, a leaked covariate or an inverted
effect sign. Every check is fast (synthetic data or metadata only) -- no notebook
output is required, so this is runnable before the first real run.

Checks are grouped by what they protect:
  A. group policy      -- MP excluded, HS quarantined as PostSurgery, counts exact;
                          every validation case survives primary_mask (A6); the
                          cross-cohort genus basis is fixed before the CLR (A8/A9);
                          the opt-in severity scale cannot touch the primary set (A10)
  B. leakage           -- tumour location refused whenever a non-cancer sample is in scope
  C. statistical signs -- DA / partial-correlation / CLR-LM directions
  D. biology           -- the polyamine EC map lists producing reactions only
  E. notebook hygiene  -- no stale constructs (SMOTE, TumorLoc_, Adenoma, Precursor);
                          CONTRACTS.md agrees with the data and with output_contract.py
                          (E7/E8); every guarded utils function is reachable (E9)

Run directly:  python checks.py        (exit 0 = all passed)
"""

import json
import re
import sys
import traceback
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import utils  # noqa: E402

NB_DIR = Path(__file__).resolve().parent
_results = []


def check(name):
    """Decorator registering a check. A check raises AssertionError to fail."""
    def deco(fn):
        _results.append((name, fn))
        return fn
    return deco


# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------

def _yachida_meta():
    """Harmonized Yachida metadata straight from the raw TSV (no pickle needed)."""
    raw = pd.read_csv(utils.DATA_DIR / "YACHIDA-CRC-2019 metadata.tsv", sep="\t")
    drop = utils.EXCLUDED_STUDY_GROUPS.get("YACHIDA-CRC-2019", [])
    raw = raw[~raw["Study.Group"].isin(drop)].reset_index(drop=True)
    return utils.harmonize_metadata(raw, "YACHIDA-CRC-2019")


def _synthetic_cohort(n=120, p=60, effect_feature=0, effect=6.0, seed=0):
    """Compositional two-group cohort with one planted elevated feature."""
    rng = np.random.default_rng(seed)
    base = rng.lognormal(0, 1, (n, p))
    half = n // 2
    grp = np.array(["Healthy"] * half + ["Early_CRC"] * (n - half))
    base[half:, effect_feature] *= effect
    rel = base / base.sum(axis=1, keepdims=True)
    spc = pd.DataFrame(rel, columns=[f"f{i}" for i in range(p)])
    spc.insert(0, "Sample", range(n))
    meta = pd.DataFrame({
        "Sample": range(n), "Stage.Group": grp, "Study.Group": grp,
        "Age": rng.normal(60, 8, n), "Gender": rng.choice(["Male", "Female"], n),
        "BMI": rng.normal(24, 3, n), "Brinkman.Index": rng.normal(200, 80, n),
        "Alcohol": rng.normal(10, 4, n),
    })
    return spc, meta, grp


# ---------------------------------------------------------------------------
# A. Group policy
# ---------------------------------------------------------------------------

@check("A1 MP is excluded, HS is a sensitivity group (not excluded, not primary)")
def a1():
    assert utils.EXCLUDED_STUDY_GROUPS["YACHIDA-CRC-2019"] == ["MP"], \
        f"expected MP excluded, got {utils.EXCLUDED_STUDY_GROUPS['YACHIDA-CRC-2019']}"
    assert utils.SENSITIVITY_GROUPS["YACHIDA-CRC-2019"] == ["HS"]
    assert "HS" not in utils.EXCLUDED_STUDY_GROUPS["YACHIDA-CRC-2019"], \
        "HS must be loaded, not dropped"
    assert utils.STAGE_GROUP_MAP["HS"] == "PostSurgery"
    assert "MP" not in utils.STAGE_GROUP_MAP, "MP must have no Stage.Group mapping"


@check("A2 PostSurgery is absent from the ordinal severity axis")
def a2():
    assert utils.STAGE_GROUP_ORDER == ["Healthy", "Early_CRC", "Advanced_CRC"], \
        utils.STAGE_GROUP_ORDER
    assert "PostSurgery" not in utils.STAGE_GROUP_ORDER
    assert np.isnan(utils.STAGE_GROUP_NUMERIC["PostSurgery"]), \
        "PostSurgery must be NaN numerically so it cannot enter a trend test"
    assert np.isnan(utils.CRC_STAGE_NUMERIC["HS"])
    assert np.isnan(utils.STUDY_GROUP_MAP["HS"])
    # And the ordinal codes must be strictly increasing over the primary groups.
    codes = [utils.STAGE_GROUP_NUMERIC[g] for g in utils.STAGE_GROUP_ORDER]
    assert codes == sorted(codes) and len(set(codes)) == len(codes), codes


@check("A3 real cohort counts: 307 loaded, 277 primary, 30 PostSurgery, 150 CRC-only")
def a3():
    y = _yachida_meta()
    assert len(y) == 307, f"expected 307 rows after dropping MP, got {len(y)}"
    assert not y["Study.Group"].eq("MP").any(), "MP leaked into the harmonized metadata"
    counts = y["Stage.Group"].value_counts().to_dict()
    assert counts.get("Healthy") == 127, counts
    assert counts.get("Early_CRC") == 96, counts
    assert counts.get("Advanced_CRC") == 54, counts
    assert counts.get("PostSurgery") == 30, counts
    assert int(utils.primary_mask(y).sum()) == 277
    assert int(utils.sensitivity_mask(y).sum()) == 30
    assert int(utils.crc_only_mask(y).sum()) == 150
    # Primary and sensitivity must be disjoint and jointly cover everything.
    assert not (utils.primary_mask(y) & utils.sensitivity_mask(y)).any()
    assert (utils.primary_mask(y) | utils.sensitivity_mask(y)).all()


@check("A4 no comparison mixes PostSurgery into the primary family")
def a4():
    for ds, comps in utils.COMPARISONS_GROUP.items():
        for g1, g2 in comps:
            assert "PostSurgery" not in (g1, g2), f"{ds}: {g1} vs {g2}"
            assert g1 in utils.STAGE_GROUP_ORDER and g2 in utils.STAGE_GROUP_ORDER
    # The sensitivity comparison must exist and be the only place PostSurgery appears.
    sens = utils.COMPARISONS_SENSITIVITY["YACHIDA-CRC-2019"]
    assert ("Healthy", "PostSurgery") in sens, sens
    for ds, comps in utils.COMPARISONS.items():
        for g1, g2 in comps:
            assert "MP" not in (g1, g2), f"{ds}: MP still compared"


@check("A5 apply_mask filters by Sample ID, not row position")
def a5():
    y = _yachida_meta()
    feat = pd.DataFrame({"Sample": y["Sample"].values[::-1],       # reversed order
                         "x": np.arange(len(y), dtype=float)})
    (sub,) = utils.apply_mask(utils.primary_mask(y), feat, meta_df=y)
    expected = set(y.loc[utils.primary_mask(y), "Sample"])
    assert set(sub["Sample"]) == expected, "apply_mask mis-sliced a reordered table"
    assert len(sub) == 277


@check("A6 primary_mask keeps every case in the unstaged validation cohorts")
def a6():
    """Regression test for the 2026-10-07 primary_mask bug.

    primary_mask() used to be an include-list over STAGE_GROUP_ORDER. Sinha and Kim
    harmonize to Healthy/"CRC", "CRC" has no Stage.Group, so their Stage.Group is
    partially NaN and the include-list silently dropped every case:
    Sinha 131 -> 89, Kim 138 -> 102. Nothing downstream raised.
    """
    expected = {"SINHA_CRC_2016": (131, 42), "Kim_adenomas_2020": (138, 36)}
    for ds, (n_total, n_cases) in expected.items():
        raw = pd.read_csv(utils.DATA_DIR / f"{ds} metadata.tsv", sep="\t")
        drop = utils.EXCLUDED_STUDY_GROUPS.get(ds, [])
        raw = raw[~raw["Study.Group"].isin(drop)].reset_index(drop=True)
        meta = utils.harmonize_metadata(raw, ds)
        pm = utils.primary_mask(meta)
        assert int(pm.sum()) == n_total, \
            f"{ds}: primary_mask kept {int(pm.sum())} of {n_total} samples"
        kept_cases = int((meta.loc[pm, "Study.Group"] == "CRC").sum())
        assert kept_cases == n_cases, \
            f"{ds}: primary set holds {kept_cases} cases, expected {n_cases}"
        # The comparison COMPARISONS[ds] asks for must be runnable at all.
        for g1, g2 in utils.COMPARISONS[ds]:
            for g in (g1, g2):
                assert int((meta.loc[pm, "Study.Group"] == g).sum()) >= 3, \
                    f"{ds}: group {g!r} has <3 samples in the primary set"


@check("A7 Yachida's Stage.Group is fully mapped (primary_mask treats NaN as primary)")
def a7():
    """primary_mask() now excludes the sensitivity groups rather than enumerating
    the primary ones, so an UNMAPPED Yachida group would silently become primary.
    That is safe only while every loaded Study.Group has a Stage.Group."""
    y = _yachida_meta()
    unmapped = sorted(set(y.loc[y["Stage.Group"].isna(), "Study.Group"].astype(str)))
    assert not unmapped, \
        f"Yachida Study.Group(s) with no STAGE_GROUP_MAP entry: {unmapped}"
    assert int(utils.primary_mask(y).sum()) == 277


@check("A8 genus_clr_table(genera=...) fixes the feature basis before the CLR")
def a8():
    """Contract for the cross-cohort feature basis.

    Two cohorts with DIFFERENT genus sets must come back with identical columns in
    identical order, absent genera as real zeros, and -- the point of the whole
    argument -- a CLR reference computed over the shared basis, not over each
    cohort's own genus set. Regression guard for the NB05 stage-3 bug, where genera
    were intersected AFTER the CLR and the two cohorts ended up on different scales.
    """
    # Two toy cohorts: A has an extra genus, B has one A lacks.
    a = pd.DataFrame({"Sample": [1, 2, 3],
                      "g__Bacteroides": [10.0, 20.0, 30.0],
                      "g__Blautia_A":   [5.0, 5.0, 5.0],
                      "g__OnlyInA":     [1.0, 2.0, 3.0]})
    b = pd.DataFrame({"Sample": [4, 5],
                      "g__Bacteroides": [10.0, 20.0],
                      "g__Blautia":     [5.0, 5.0],
                      "g__OnlyInB":     [7.0, 8.0]})

    # normalise_genus_name folds Blautia_A -> Blautia, so it is shared.
    basis = ["Bacteroides", "Blautia", "OnlyInA"]
    ta = utils.genus_clr_table(a, genera=basis)
    tb = utils.genus_clr_table(b, genera=basis)

    assert list(ta.columns) == ["Sample"] + basis, list(ta.columns)
    assert list(tb.columns) == list(ta.columns), "column order differs between cohorts"
    assert len(ta) == 3 and len(tb) == 2

    # A genus absent from cohort B must be the ramp's floor, not dropped and not NaN.
    assert tb["OnlyInA"].notna().all(), "absent genus came back NaN instead of zero"
    assert (tb["OnlyInA"] < tb["Bacteroides"]).all(), \
        "absent genus should sit below an observed one after CLR"

    # The CLR reference must come from the SHARED basis. Cohort A's rows 1 and 2 have
    # identical composition over the basis up to a constant factor, so CLR -- which is
    # scale-invariant -- must give them the same Bacteroides value. If the prevalence
    # filter or a per-cohort basis had crept back in, this breaks.
    scaled = a.copy()
    scaled[["g__Bacteroides", "g__Blautia_A", "g__OnlyInA"]] *= 7.0
    ts = utils.genus_clr_table(scaled, genera=basis)
    assert np.allclose(ts[basis].to_numpy(float), ta[basis].to_numpy(float)), \
        "CLR is not scale-invariant on a fixed basis -- the basis is not actually fixed"

    # Duplicates in the basis would silently misalign features.
    try:
        utils.genus_clr_table(a, genera=["Bacteroides", "Bacteroides"])
    except ValueError:
        pass
    else:
        raise AssertionError("duplicate genera in `genera` must raise")


@check("A9 shared_genera applies prevalence to TRAIN only, presence to targets")
def a9():
    """Feature selection must not consult the target's prevalence -- that is the
    held-out cohort choosing the model's features."""
    # "Rare" is in 1 of 10 training samples (10% prevalence, below the 20% default),
    # so it must be excluded by the TRAIN prevalence filter.
    n = 10
    train = pd.DataFrame({
        "Sample": range(n),
        "g__Common": np.linspace(1, 10, n),
        "g__Rare":   [1.0] + [0.0] * (n - 1),
        "g__TrainOnly": np.linspace(1, 5, n),
    })
    # In the target, "Rare" is in every sample -- but that must not rescue it.
    target = pd.DataFrame({
        "Sample": [100, 101, 102],
        "g__Common": [1.0, 2.0, 3.0],
        "g__Rare":   [9.0, 9.0, 9.0],
        "g__TargetOnly": [4.0, 5.0, 6.0],
    })

    got = utils.shared_genera(train, [target], prevalence=0.20)
    assert got == ["Common"], (
        f"expected ['Common'], got {got}. 'Rare' is below the TRAIN prevalence "
        "threshold and must stay out however prevalent it is in the target; "
        "'TrainOnly'/'TargetOnly' are not shared."
    )
    # A target_samples list of the wrong length is a silent mis-pairing otherwise.
    try:
        utils.shared_genera(train, [target], target_samples=[])
    except ValueError:
        pass
    else:
        raise AssertionError("mismatched target_samples length must raise")


@check("A10 severity-transfer scale is separate from the primary stage encoding")
def a10():
    """The opt-in adenoma analysis must not be able to disturb the primary set."""
    # MP must still have NO Stage.Group mapping: that is what keeps adenomas out of
    # primary_mask(), which now EXCLUDES sensitivity groups rather than enumerating.
    assert "MP" not in utils.STAGE_GROUP_MAP, \
        "MP gained a Stage.Group mapping -- adenomas can now enter the primary set"
    assert "Adenoma" not in utils.STAGE_GROUP_ORDER

    # The parallel scales must be ordered and cover both cohorts' vocabularies.
    assert utils.SEVERITY_TRAIN_ORDER == ["Healthy", "Adenoma", "Early_CRC", "Advanced_CRC"]
    assert utils.SEVERITY_EVAL_ORDER == ["Healthy", "Adenoma", "Carcinoma"]
    assert utils.SEVERITY_TRAIN_MAP["MP"] == 1, "Yachida MP must be the adenoma level"
    assert utils.SEVERITY_EVAL_MAP["Adenoma"] == 1, "Kim Adenoma must be the adenoma level"
    assert utils.SEVERITY_EVAL_MAP["MP"] == utils.SEVERITY_EVAL_MAP["Adenoma"], \
        "Yachida MP and Kim Adenoma must land on the same eval level"
    for off in ("HS", "PostSurgery"):
        assert np.isnan(utils.SEVERITY_EVAL_MAP[off]), f"{off} must have no severity position"
        assert off not in utils.SEVERITY_TRAIN_MAP

    # Collapsing train -> eval must be monotone and must merge Early+Advanced.
    t2e = utils.SEVERITY_TRAIN_TO_EVAL
    assert t2e[2] == t2e[3] == 2, "Early_CRC and Advanced_CRC must collapse to Carcinoma"
    assert [t2e[k] for k in sorted(t2e)] == sorted(t2e[k] for k in sorted(t2e))

    # harmonize_metadata must carry Kim's Adenoma through when it is admitted.
    # Omitting it from the remap would turn all 102 adenomas into NaN silently.
    raw = pd.read_csv(utils.DATA_DIR / "Kim_adenomas_2020 metadata.tsv", sep="\t")
    hm = utils.harmonize_metadata(raw, "Kim_adenomas_2020", exclude={})
    counts = hm["Study.Group"].value_counts().to_dict()
    assert counts.get("Adenoma") == 102, (
        f"expected 102 Kim adenomas with exclude={{}}, got {counts}. "
        "harmonize_metadata's Kim remap must include 'Adenoma': 'Adenoma'."
    )
    assert hm["Study.Group"].notna().all(), "a Kim label became NaN in the remap"
    # And the default must still drop them, so the primary pipeline is untouched.
    assert "Adenoma" not in set(
        utils.harmonize_metadata(raw, "Kim_adenomas_2020")["Study.Group"])


# ---------------------------------------------------------------------------
# B. Leakage
# ---------------------------------------------------------------------------

@check("B1 default confounders carry no tumour-location column")
def b1():
    y = _yachida_meta()
    cols = [c for c in utils.build_confounder_matrix(y).columns if c != "Sample"]
    assert cols == ["Age", "BMI", "Brinkman.Index", "Alcohol", "Gender_enc"], cols
    assert not any("Tumor" in c or "TumorLoc" in c or "TumorSite" in c for c in cols)


@check("B2 tumour location is REFUSED whenever a non-cancer sample is in scope")
def b2():
    y = _yachida_meta()
    for label, subset in [
        ("full cohort", y),
        ("Healthy only", y[y["Stage.Group"] == "Healthy"]),
        ("PostSurgery only", y[utils.sensitivity_mask(y)]),
        ("CRC + one Healthy", pd.concat([y[utils.crc_only_mask(y)],
                                         y[y["Stage.Group"] == "Healthy"].head(1)])),
    ]:
        try:
            utils.build_confounder_matrix(subset.reset_index(drop=True),
                                          include_tumor_location=True)
        except ValueError:
            continue
        raise AssertionError(f"tumour location was allowed on: {label}")


@check("B3 tumour location IS allowed on a CRC-only set, with correct site counts")
def b3():
    y = _yachida_meta()
    ycrc = y[utils.crc_only_mask(y)].reset_index(drop=True)
    conf = utils.build_confounder_matrix(ycrc, include_tumor_location=True)
    sites = {c: int(conf[c].sum()) for c in conf.columns if c.startswith("TumorSite_")}
    assert sites == {"TumorSite_Rectum": 69, "TumorSite_Left_colon": 32,
                     "TumorSite_Right_colon": 44, "TumorSite_MultiSite": 5}, sites
    # Dummies must be exhaustive: every CRC sample has exactly one site.
    assert int(sum(sites.values())) == len(ycrc) == 150


@check("B4 tumour site is NaN for every non-cancer sample")
def b4():
    y = _yachida_meta()
    site = y["Tumor.Location.Site"]
    assert site[utils.crc_only_mask(y)].notna().all(), "a CRC sample has no site"
    assert site[~utils.crc_only_mask(y)].isna().all(), \
        "a non-cancer sample carries a tumour site, which would encode disease status"


# ---------------------------------------------------------------------------
# C. Statistical signs
# ---------------------------------------------------------------------------

@check("C1 differential_abundance: MeanDiff > 0 means higher in group2")
def c1():
    n = 60
    rng = np.random.default_rng(1)
    df = pd.DataFrame({"Sample": range(n),
                       "up": np.r_[rng.normal(0, 1, 30), rng.normal(3, 1, 30)],
                       "down": np.r_[rng.normal(3, 1, 30), rng.normal(0, 1, 30)]})
    meta = pd.DataFrame({"Sample": range(n),
                         "Study.Group": ["Healthy"] * 30 + ["Early_CRC"] * 30})
    res = utils.differential_abundance(df, meta, "Healthy", "Early_CRC").set_index("Feature")
    assert res.loc["up", "MeanDiff"] > 0, "feature higher in group2 got MeanDiff <= 0"
    assert res.loc["down", "MeanDiff"] < 0
    # Cliff's delta has the OPPOSITE sign to MeanDiff (documented convention).
    assert res.loc["up", "CliffsD"] < 0, "CliffsD sign convention changed"
    assert res.loc["down", "CliffsD"] > 0


@check("C2 clr_lm_da recovers a planted effect with the right direction")
def c2():
    spc, meta, _ = _synthetic_cohort()
    res = utils.clr_lm_da(spc, meta, "Healthy", "Early_CRC", already_clr=False,
                          covars=["Age", "Gender", "BMI", "Brinkman.Index", "Alcohol"])
    assert not res.empty
    row = res[res["Feature"] == "f0"].iloc[0]
    assert row["Significant"], f"planted effect missed (q={row['QValue']:.3g})"
    assert row["Coef"] > 0, "planted elevated feature got a negative coefficient"
    assert row["Direction"] == "Up in Early_CRC", row["Direction"]
    n_sig = int(res["Significant"].sum())
    assert n_sig <= 5, f"{n_sig} significant features where 1 was planted"


@check("C3 clr_lm_da is calibrated on null data (no planted effect)")
def c3():
    spc, meta, _ = _synthetic_cohort(effect=1.0, seed=7)   # effect=1.0 -> no effect
    res = utils.clr_lm_da(spc, meta, "Healthy", "Early_CRC", already_clr=False,
                          covars=["Age", "Gender"])
    n_sig = int(res["Significant"].sum())
    assert n_sig == 0, f"{n_sig} false positives on null data"


@check("C4 clr_lm_da refuses a too-small arm and drops constant covariates")
def c4():
    spc, meta, _ = _synthetic_cohort(n=20)
    tiny = meta.copy()
    tiny.loc[tiny.index[:18], "Stage.Group"] = "Healthy"
    tiny.loc[tiny.index[18:], "Stage.Group"] = "Early_CRC"   # n=2 in one arm
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = utils.clr_lm_da(spc, tiny, "Healthy", "Early_CRC", already_clr=False)
    assert res.empty, "a 2-sample arm should be refused"

    spc2, meta2, _ = _synthetic_cohort()
    meta2["Brinkman.Index"] = 0.0          # constant -> singular design
    res2 = utils.clr_lm_da(spc2, meta2, "Healthy", "Early_CRC", already_clr=False,
                           covars=["Age", "Brinkman.Index"])
    assert "Brinkman.Index" in res2.attrs["dropped_constant_covariates"]
    assert "Brinkman.Index" not in res2.attrs["covariates"]


@check("C5 progression correlation FDR-corrects each marginal test once")
def c5():
    rng = np.random.default_rng(3)
    n = 80
    stage = np.repeat([0, 1, 2, 3], n // 4).astype(float)
    spc = pd.DataFrame({"Sample": range(n),
                        "s_up": stage + rng.normal(0, 0.5, n),
                        "s_flat": rng.normal(0, 1, n)})
    mtb = pd.DataFrame({"Sample": range(n),
                        "m_up": stage + rng.normal(0, 0.5, n),
                        "m_flat": rng.normal(0, 1, n)})
    meta = pd.DataFrame({"Sample": range(n), "CRC.Stage.Numeric": stage})
    res = utils.progression_correlation_analysis(spc, mtb, meta)
    assert not res.empty
    # 2 species x 2 metabolites = 4 pairs, but only 2 unique species-stage tests.
    assert len(res) == 4, len(res)
    q_up = res.loc[res["Species"] == "s_up", "QVal_species_stage"].unique()
    assert len(q_up) == 1, "the same species-stage test got different q-values"
    assert q_up[0] < 0.05
    assert res.loc[res["Species"] == "s_up", "Rho_species_stage"].iloc[0] > 0


@check("C6 metabolite normalisation is blind to non-detects outside its reference set")
def c6():
    """Regression guard for the zero-inclusive median-centring artefact.

    The old transform centred each sample on the median of log10(x+1) over ALL
    metabolites with non-detects as 0, so a group with slightly more non-detects got
    a smaller offset and every detected value rose: 38 of 40 significant
    Healthy-vs-Advanced metabolites were "elevated". This tests the mechanism, not a
    simulated p-value count: extra non-detects in non-reference metabolites must not
    move any sample's offset or any detected value -- and the test must FAIL on the
    old transform, or it guards nothing.
    """
    rng = np.random.default_rng(6)
    n, p_ref, p_other = 40, 15, 25
    ref = 10 ** rng.normal(3.0, 0.5, (n, p_ref))            # always detected
    other = 10 ** rng.normal(2.0, 0.8, (n, p_other))
    other[rng.random((n, p_other)) < 0.4] = 0.0              # zero-inflated
    cols_ref = [f"r{j}" for j in range(p_ref)]
    cols_other = [f"o{j}" for j in range(p_other)]
    base = pd.DataFrame(np.hstack([ref, other]), columns=cols_ref + cols_other)
    base.insert(0, "Sample", range(n))

    more = base.copy()
    extra = rng.random((n, p_other)) < 0.3
    extra[: n // 2] = False                                  # only the second half
    more[cols_other] = more[cols_other].where(~extra, 0.0)

    t0, i0 = utils.transform_metabolites(base)
    t1, i1 = utils.transform_metabolites(more)
    assert i0["reference"] == cols_ref, "reference set should be the always-detected block"
    assert np.allclose(i0["offset"].to_numpy(), i1["offset"].to_numpy()), \
        "non-detects outside the reference set moved the normalisation offsets"
    still = (more[cols_other].to_numpy() > 0)
    a = t0[cols_other].to_numpy()[still]
    b = t1[cols_other].to_numpy()[still]
    assert np.allclose(a, b), "detected values moved when other metabolites lost detection"

    # The test discriminates: the old zero-inclusive median DOES move.
    old = lambda d: np.median(np.log10(d.drop(columns="Sample").to_numpy(float) + 1.0), axis=1)
    assert not np.allclose(old(base), old(more)), \
        "test would not have caught the old transform -- it guards nothing"

    # Non-detects: one tied floor per metabolite, strictly below every detected value.
    X1 = t1[cols_other].to_numpy()
    for j in range(p_other):
        nd = more[cols_other[j]].to_numpy() == 0
        if nd.any() and (~nd).any():
            assert np.unique(X1[nd, j]).size == 1, "non-detects are not tied"
            assert X1[nd, j][0] < X1[~nd, j].min(), "a non-detect ranks above a detection"


@check("C7 partial_spearman_grid reproduces pingouin.partial_corr(method='spearman')")
def c7():
    import pingouin as pg
    rng = np.random.default_rng(7)
    n = 120
    X = rng.normal(size=(n, 5))
    X[:, 0] = np.round(X[:, 0])                              # heavy ties
    Y = rng.normal(size=(n, 4)) + 0.4 * X[:, :4]
    Z = np.column_stack([rng.normal(size=n), rng.integers(0, 2, n).astype(float)])
    r, p = utils.partial_spearman_grid(X, Y, Z)
    df = pd.DataFrame(np.column_stack([X, Y, Z]),
                      columns=[f"x{i}" for i in range(5)] + [f"y{j}" for j in range(4)]
                      + ["z0", "z1"])
    for i in range(5):
        for j in range(4):
            res = pg.partial_corr(data=df, x=f"x{i}", y=f"y{j}", covar=["z0", "z1"],
                                  method="spearman")
            assert abs(res["r"].iloc[0] - r[i, j]) < 1e-10, (i, j, res["r"].iloc[0], r[i, j])
            pref = res["p-val"].iloc[0]
            assert abs(np.log10(pref) - np.log10(p[i, j])) < 1e-8, (i, j, pref, p[i, j])
    # A constant column yields NaN, never an exception or a silently dropped pair.
    Xc = np.column_stack([X, np.ones(n)])
    rc, pc = utils.partial_spearman_grid(Xc, Y, Z)
    assert np.isnan(rc[-1]).all() and np.isnan(pc[-1]).all()


@check("C8 transfer features harmonise detection depth and share one CLR delta")
def c8():
    """A shotgun training cohort resolves far smaller proportions than a 16S target.
    Unharmonised, the training set has few zeros and the target many, so a transferred
    score can learn detection depth as if it were disease."""
    rng = np.random.default_rng(8)
    genera = [f"G{j}" for j in range(12)]

    def cohort(n, lod):
        P = rng.dirichlet(np.full(len(genera), 0.15), size=n)    # many tiny proportions
        P = np.where(P < lod, 0.0, P)
        df = pd.DataFrame(P, columns=[f"d__B;p__P;c__C;o__O;f__F;g__{g}" for g in genera])
        df.insert(0, "Sample", [f"s{lod}_{i}" for i in range(n)])
        return df

    train = cohort(80, 1e-9)          # shotgun-like: almost nothing censored
    target = cohort(60, 1e-3)         # 16S-like: a coarse detection limit
    tr, te, info = utils.transfer_genus_features(train, target, genera)
    assert list(tr.columns) == list(te.columns) == ["Sample"] + genera
    assert info["zero_frac_train"] > info["zero_frac_train_before"], \
        "detection harmonisation did not censor the training cohort"
    assert abs(info["zero_frac_train"] - info["zero_frac_target"]) < \
        abs(info["zero_frac_train_before"] - info["zero_frac_target"]), \
        "harmonisation did not move the training zero fraction toward the target's"
    assert np.isclose(info["delta"], utils.CLR_DELTA_FRACTION * info["target_detection_limit"]), \
        "the two cohorts do not share one delta set by the target's detection limit"
    for M in (tr[genera].to_numpy(), te[genera].to_numpy()):
        sd = M.std(axis=0, ddof=1)
        assert np.allclose(M.mean(axis=0), 0, atol=1e-9)
        assert np.allclose(sd[sd > 0], 1.0), "features are not standardised within cohort"
    try:
        utils.transfer_genus_features(train, target, genera + [genera[0]])
    except ValueError:
        pass
    else:
        raise AssertionError("duplicate genera must raise")


@check("C9 replication separates 'not measurable' from 'not replicated'; flags contradictions")
def c9():
    from scipy.stats import spearmanr as _sp
    rng = np.random.default_rng(9)
    # spearman_grid == scipy.stats.spearmanr, ties included.
    A = np.round(rng.normal(size=(60, 3)), 1)
    B = rng.normal(size=(60, 2)) + 0.5 * A[:, :2]
    r, p = utils.spearman_grid(A, B)
    for i in range(3):
        for j in range(2):
            ref = _sp(A[:, i], B[:, j])
            assert abs(ref.statistic - r[i, j]) < 1e-10 and abs(ref.pvalue - p[i, j]) < 1e-9

    n = 150
    base = rng.normal(size=n)

    def cohort(name, ids, signs, genera):
        G = pd.DataFrame({"Sample": range(n)})
        for g in genera:
            G[g] = base + rng.normal(0, 0.5, n) if g == "G1" else rng.normal(size=n)
        M = pd.DataFrame({"Sample": range(n)})
        for comp, s in signs.items():
            M[comp] = s * base + rng.normal(0, 0.5, n)
        ident = utils.metabolite_identity(
            pd.DataFrame({"Compound": list(ids), "KEGG": list(ids.values())}), name)
        return {"name": name, "genus": G, "mtb": M, "identity": ident}

    disc = cohort("D", {"dA": "K1", "dB": "K2", "dC": "K3"},
                  {"dA": 1, "dB": 1, "dC": 1}, ["G1", "G2"])
    t1 = cohort("T1", {"x1": "K1", "x2": "K2"}, {"x1": 1, "x2": -1}, ["G1", "G2"])
    t2 = cohort("T2", {"y1": "K1"}, {"y1": 1}, ["G2"])          # G1 absent here
    pairs, summary = utils.genus_metabolite_replication(disc, [t1, t2])
    row = pairs.set_index(["Genus", "Compound"])
    assert {("G1", "dA"), ("G1", "dB"), ("G1", "dC")} <= set(row.index), "discovery hits missing"
    a, b, c = row.loc[("G1", "dA")], row.loc[("G1", "dB")], row.loc[("G1", "dC")]
    assert a["In_T1"] and a["Compound_T1"] == "x1" and a["Replicated_T1"], \
        "an identity-matched, same-sign hit was not replicated"
    assert b["Contradicted_T1"] and not b["Replicated_T1"], "opposite sign not flagged"
    assert not c["In_T1"] and np.isnan(c["Rho_T1"]) and not c["Replicated_T1"] \
        and not c["Contradicted_T1"], "an unmeasured compound was scored as a failure"
    assert not a["In_T2"], "a pair whose genus is absent was scored as measurable"
    s1 = summary.set_index("Target").loc["T1"]
    assert s1["Not_measurable_compound"] >= 1 and s1["Replicated"] >= 1 \
        and s1["Contradicted"] >= 1


# ---------------------------------------------------------------------------
# D. Biology
# ---------------------------------------------------------------------------

@check("D1 POLYAMINE_EC lists only reactions that PRODUCE each compound")
def d1():
    # Curated ground truth: every EC below has the keyed compound as a PRODUCT.
    # The pre-v2 map listed consuming enzymes for ornithine/arginine/putrescine,
    # which inverted the producer call for those three compounds.
    must_have = {
        "C00134": {"4.1.1.17", "3.5.3.11", "3.5.1.53"},   # ODC, agmatinase, NCP amidohydrolase
        "C00315": {"2.5.1.16"},                            # spermidine synthase
        "C00750": {"2.5.1.22"},                            # spermine synthase
        "C00179": {"4.1.1.19"},                            # arginine decarboxylase
        "C01672": {"4.1.1.18"},                            # lysine decarboxylase
        "C00077": {"3.5.3.1"},                             # arginase
        "C00062": {"4.3.2.1"},                             # argininosuccinate lyase
    }
    # Enzymes that CONSUME the keyed compound and must never appear under it.
    must_not = {
        "C00134": {"4.1.1.18", "4.1.1.19", "2.5.1.16"},   # make cadaverine/agmatine/spermidine
        "C00077": {"4.1.1.17"},                            # ODC consumes ornithine
        "C00062": {"4.1.1.19", "3.5.3.11", "3.5.3.1"},     # ADC/agmatinase/arginase consume arginine
        "C00179": {"3.5.3.11"},                            # agmatinase consumes agmatine
    }
    for kid, need in must_have.items():
        have = set(utils.POLYAMINE_EC[kid])
        assert need <= have, f"{kid}: missing producer EC(s) {need - have}"
    for kid, bad in must_not.items():
        have = set(utils.POLYAMINE_EC[kid])
        assert not (bad & have), f"{kid}: lists consuming EC(s) {bad & have}"
    assert set(utils.POLYAMINE_EC) == set(utils.POLYAMINE_KEGG), "EC map / KEGG map mismatch"


@check("D2 pathway-membership KEGG set excludes the near-universal precursors")
def d2():
    assert "C00077" not in utils.POLYAMINE_PATHWAY_KEGG, "ornithine in pathway set"
    assert "C00062" not in utils.POLYAMINE_PATHWAY_KEGG, "arginine in pathway set"
    assert set(utils.POLYAMINE_PATHWAY_KEGG) == {"C00134", "C00315", "C00750",
                                                 "C00179", "C01672"}


@check("D3 every polyamine proxy maps to a distinct parent and is pre-specified")
def d3():
    for parent, ids in utils.POLYAMINE_PROXY_KEGG.items():
        assert parent in utils.POLYAMINE_KEGG, parent
        assert parent not in ids, f"{parent} listed as its own proxy"
    for parent in utils.POLYAMINE_PROXY_NAMES:
        assert parent in utils.POLYAMINE_KEGG, parent
    overlap = set(utils.POLYAMINE_PROXY_KEGG) & set(utils.POLYAMINE_PROXY_NAMES)
    assert not overlap, f"compound has both a KEGG and a name proxy: {overlap}"


@check("D4 Kim/Sinha polyamine lookup finds the expected direct + proxy compounds")
def d4():
    expect = {
        "SINHA_CRC_2016":    ({"C00134", "C01672", "C00077", "C00062"}, set()),
        "Kim_adenomas_2020": ({"C00315", "C00077", "C00062"}, {"C00134", "C01672"}),
    }
    for ds, (exp_direct, exp_proxy) in expect.items():
        mtb = pd.read_csv(utils.DATA_DIR / f"{ds} mtb.tsv", sep="\t")
        mp_path = utils.DATA_DIR / f"{ds} mtb.map.tsv"
        if not mp_path.exists():
            mp_path = utils.DATA_DIR / f"{ds} mtb_map.tsv"
        mmap = pd.read_csv(mp_path, sep="\t")
        direct, proxy = utils.find_polyamine_columns_with_proxies(mmap, mtb, ds)
        assert set(direct) == exp_direct, f"{ds} direct: {set(direct)} != {exp_direct}"
        assert set(proxy) == exp_proxy, f"{ds} proxy: {set(proxy)} != {exp_proxy}"


@check("D5 AGORA2 genome-scale evidence is disabled and nothing silently re-enables it")
def d5():
    # v2.1 removed genome-scale evidence: the EC-index build needs 1.7 GB of
    # reconstructions on disk and MICOM needs ~100 MB of model per sample, which
    # overloaded the analysis server. This check makes the removal verifiable
    # rather than a claim, and it is the guard that stops a notebook from
    # triggering a 1.7 GB scan from inside a pipeline run.
    assert utils.AGORA2_ENABLED is False, \
        "utils.AGORA2_ENABLED is True -- v2.1 does not support genome-scale evidence"

    # The loader must refuse rather than rebuild when the cache is absent.
    if not utils.AGORA2_INDEX_CACHE.exists():
        utils._AGORA2_INDEX_MEMO = None
        try:
            utils.load_agora2_ec_index()
        except RuntimeError:
            pass
        else:
            raise AssertionError(
                "load_agora2_ec_index() did not refuse with the index cache absent -- it "
                "may still fall back to scanning Data/micom, which is what overloaded "
                "the server."
            )

    # No notebook may call the AGORA2 helpers.
    banned = (r"load_agora2_ec_index|species_encodes_enzyme|species_has_any_ec|"
              r"build_agora2_ec_index|find_agora2_data_dir|agora2_evidence")
    hits = []
    for name, idx, ctype, text in _notebook_sources():
        if ctype != "code":
            continue
        for lineno, line in enumerate(text.splitlines(), 1):
            if line.strip().startswith("#"):
                continue
            if re.search(banned, line):
                hits.append(f"{name} cell {idx} line {lineno}")
    assert not hits, ("notebook(s) still call AGORA2 enzyme helpers, which cannot run in "
                      "v2.1: " + ", ".join(hits))


@check("D6 the curated polyamine EC map is preserved for future re-enabling")
def d6():
    # POLYAMINE_EC is unused while AGORA2 is disabled, but the curation is the
    # thing that was previously WRONG (consuming enzymes listed as producers), so
    # it is kept and still guarded by D1. This check only asserts it survives, so
    # a future re-enable starts from the corrected map.
    assert utils.POLYAMINE_EC and set(utils.POLYAMINE_EC) == set(utils.POLYAMINE_KEGG)
    assert utils.SCFA_EC and set(utils.SCFA_EC) == {"butyrate", "propionate"}


# ---------------------------------------------------------------------------
# E. Notebook hygiene
# ---------------------------------------------------------------------------

def _notebook_sources():
    for nb in sorted(NB_DIR.glob("*.ipynb")):
        doc = json.loads(nb.read_text(encoding="utf-8"))
        for i, cell in enumerate(doc.get("cells", [])):
            src = cell.get("source", [])
            text = "".join(src) if isinstance(src, list) else str(src)
            yield nb.name, i, cell.get("cell_type"), text


def _is_guard_line(line):
    """True if this line is defensive code that legitimately NAMES a banned
    construct in order to reject it (an assert, a banned-prefix list, a
    startswith test) or is a comment explaining one.

    Without this, the hygiene checks fire on exactly the code that enforces
    them -- e.g. `assert not any(c.startswith(('TumorSite_', 'TumorLoc_')))`.
    """
    s = line.strip()
    if s.startswith("#"):
        return True
    return any(tok in line for tok in
               ("assert", "startswith", "_banned", "raise ", "LEAK GUARD", "not any("))


@check("E1 no notebook USES a removed construct (SMOTE / TumorLoc_ / Adenoma or MP as a group / Precursor)")
def e1():
    # Banned everywhere: never a legitimate use in v2.1.
    banned_any = {
        "SMOTE":     r"\bSMOTE\b",
        "TumorLoc_": r"TumorLoc_",     # the leaking dummy prefix; TumorSite_ replaced it
        "Precursor": r"\bPrecursor\b",
    }
    # Banned in CODE only: a *quoted* group label, i.e. the stale grouping used as
    # data. Prose may still explain why MP and Kim's Adenoma arm are excluded.
    banned_code = {
        "'Adenoma' as a group label": r"""['"]Adenoma['"]""",
        "'MP' as a group label":      r"""['"]MP['"]""",
    }
    # Banned in MARKDOWN only: phrasing that asserts the stale MP->Adenoma scheme.
    banned_md = {
        "MP described as the Adenoma group":
            r"Adenoma\s*\(MP|MP\s*(?:=|->|→)\s*(?:the\s*)?[`'\"]?Adenoma",
    }
    hits = []
    for name, idx, ctype, text in _notebook_sources():
        patterns = dict(banned_any)
        patterns.update(banned_code if ctype == "code" else banned_md)
        for lineno, line in enumerate(text.splitlines(), 1):
            if ctype == "code" and _is_guard_line(line):
                continue          # defensive code may name what it rejects
            for label, pat in patterns.items():
                if re.search(pat, line):
                    hits.append(f"{name} cell {idx} ({ctype}) line {lineno}: {label}")
    assert not hits, "stale constructs still in use:\n  " + "\n  ".join(hits)


@check("E2 no notebook states the wrong HS definition or the old HS-exclusion claim")
def e2():
    # HS means "healthy after surgery" -- patients whose colorectal cancer was
    # resected and who are currently cancer-free. The earlier pipelines described
    # it as "history of colorectal surgery", which is a different and misleading
    # claim (it omits that these are former cancer patients) and was the basis for
    # excluding the group outright. Both the wrong definition and the stale
    # "excluded from the study" claim are banned, in code and in prose.
    stale = {
        "wrong HS definition ('history of colorectal surgery')":
            r"history of colorectal surgery",
        "stale claim that HS is excluded from the study":
            r"HS[^.\n]{0,60}excluded from (?:v2|every|the study|the analysis)",
        "stale assertion that HS must be absent":
            r"HS\s+(?:samples\s+)?must\s+be\s+excluded",
    }
    hits = []
    for name, idx, ctype, text in _notebook_sources():
        for label, pat in stale.items():
            if re.search(pat, text, flags=re.IGNORECASE):
                hits.append(f"{name} cell {idx} ({ctype}): {label}")
    # utils.py and CONTRACTS.md carry the definition too, so check them as well.
    for extra in ("utils.py", "CONTRACTS.md"):
        p = NB_DIR / extra
        if p.exists():
            txt = p.read_text(encoding="utf-8")
            for label, pat in stale.items():
                if re.search(pat, txt, flags=re.IGNORECASE):
                    hits.append(f"{extra}: {label}")
    assert not hits, "stale HS description(s):\n  " + "\n  ".join(hits)


@check("E3 CONTRACTS.md documents the v2.1 group policy")
def e3():
    txt = (NB_DIR / "CONTRACTS.md").read_text(encoding="utf-8")
    required = ["PostSurgery", "primary_mask", "crc_only_mask", "Kim_adenomas_2020",
                "clr_lm_da", "pair_evidence_nb08_v2.csv", "Robust", "sensitivity_mask"]
    missing = [t for t in required if t not in txt]
    assert not missing, f"CONTRACTS.md does not document: {missing}"
    # The stale scheme listed Adenoma as a live group on the severity axis.
    stale = [p for p in (r"Healthy\s*/\s*Adenoma", r"Adenoma\s*\(MP",
                         r"`Adenoma`\s*group")
             if re.search(p, txt)]
    assert not stale, f"CONTRACTS.md still presents Adenoma as a group: {stale}"


@check("E4 generators are retired: no edits_nb*.py in v2/ (the notebooks are authoritative)")
def e4():
    # The edits_nb*.py generators had drifted from the notebooks in both directions:
    # 12 of 14 differed, the NB08 trio could not produce a runnable NB08 in any order,
    # and the GLASSO one held a configuration measured to make convergence worse.
    # Their reviewed content was transplanted by cell id and they were archived under
    # archive/generators/. Re-running one would silently overwrite hand fixes, so an
    # active generator in v2/ is a failure, not a warning.
    active = sorted(p.name for p in NB_DIR.glob("edits_nb*.py"))
    assert not active, (f"active notebook generators in v2/: {active}. Move them to "
                        "archive/generators/ -- re-running one overwrites hand fixes.")
    assert (NB_DIR / "archive" / "generators" / "README.md").exists(), \
        "archive/generators/README.md is missing (it records what was transplanted)"


@check("E5 no blanket warning suppression in notebooks or shared code")
def e5():
    # NB03 and NB05 used warnings.filterwarnings('ignore'), which muted convergence,
    # numerical and missing-glyph warnings; NB09's GLASSO non-convergence went
    # unreported the same way. Suppression must name a category.
    bare = re.compile(r"""(filterwarnings|simplefilter)\(\s*['"]ignore['"]\s*\)""")
    hits = []
    for name, idx, ctype, text in _notebook_sources():
        if ctype != "code":
            continue
        for lineno, line in enumerate(text.splitlines(), 1):
            if not line.strip().startswith("#") and bare.search(line):
                hits.append(f"{name} cell {idx} line {lineno}")
    for extra in ("utils.py", "viz.py", "run_pipeline.py"):
        for lineno, line in enumerate((NB_DIR / extra).read_text(encoding="utf-8").splitlines(), 1):
            if not line.strip().startswith("#") and bare.search(line):
                hits.append(f"{extra} line {lineno}")
    assert not hits, ("blanket warning suppression (name a category instead): "
                      + ", ".join(hits))


@check("E6 no notebook defines its own genus or CLR helper")
def e6():
    # NB02 had its own _norm_genus / genus_clr, which kept 'Unclassified' as a genus
    # and drifted from utils.genus_clr_table -- the helper the A8/C8 guards test.
    # A private copy is code no check covers.
    rx = re.compile(r"^\s*def\s+(_?norm(alise)?_genus\w*|_?genus_clr\w*|_?clr\w*)\s*\(")
    hits = [f"{name} cell {idx} line {k}"
            for name, idx, ctype, text in _notebook_sources() if ctype == "code"
            for k, line in enumerate(text.splitlines(), 1) if rx.search(line)]
    assert not hits, ("notebook-local genus/CLR helper (use utils.genus_clr_table / "
                      "utils.clr_transform): " + ", ".join(hits))


def _quoted_ints(pattern, txt, what):
    m = re.search(pattern, txt)
    assert m, f"CONTRACTS.md no longer states {what} in the expected form ({pattern!r})"
    return tuple(int(g) for g in m.groups())


@check("E7 CONTRACTS.md quotes the cohort counts the data produce")
def e7():
    # The counts are recomputed from the metadata here and compared with the prose, so
    # the document cannot drift from the data (A3/A6 pin the same numbers in code).
    txt = (NB_DIR / "CONTRACTS.md").read_text(encoding="utf-8")
    y = _yachida_meta()
    yc = (len(y), int(utils.primary_mask(y).sum()), int(utils.sensitivity_mask(y).sum()),
          int(utils.crc_only_mask(y).sum()))
    q = _quoted_ints(r"Counts: \*\*(\d+) loaded, (\d+) primary, (\d+) sensitivity, "
                     r"(\d+) carcinomas\.\*\*", txt, "the Yachida counts line")
    assert q == yc, f"Yachida counts line quotes {q}, the data give {yc}"
    q = _quoted_ints(r"\| `YACHIDA-CRC-2019` \|[^|\n]*\| (\d+) primary \(\+(\d+) sensitivity\)",
                     txt, "Yachida's cohort-table row")
    assert q == yc[1:3], f"Yachida cohort row quotes {q}, the data give {yc[1:3]}"
    q = _quoted_ints(r"`crc_only_mask\(meta\)` \|[^|\n]*\(n=(\d+)\)", txt, "the CRC-only mask row")
    assert q == (yc[3],), f"crc_only_mask row quotes n={q[0]}, the data give {yc[3]}"
    for ds in ("SINHA_CRC_2016", "Kim_adenomas_2020"):
        raw = pd.read_csv(utils.DATA_DIR / f"{ds} metadata.tsv", sep="\t")
        raw = raw[~raw["Study.Group"].isin(utils.EXCLUDED_STUDY_GROUPS.get(ds, []))]
        meta = utils.harmonize_metadata(raw.reset_index(drop=True), ds)
        groups = meta.loc[utils.primary_mask(meta), "Study.Group"]
        got = (len(groups), int((groups == "Healthy").sum()), int((groups == "CRC").sum()))
        q = _quoted_ints(rf"\| `{ds}` \|[^|\n]*\| (\d+) \((\d+) Healthy \+ (\d+) CRC\)",
                         txt, f"{ds}'s cohort-table row")
        assert q == got, f"{ds} cohort row quotes {q}, the data give {got}"


# Files CONTRACTS.md names that are pipeline INPUTS, not outputs (paths under Data/ are
# inputs by location and need no entry).
_DOC_INPUTS = {
    "manifest.csv": "NB10 input: the MICOM model manifest in Data/micom/models_agora2/",
}


@check("E8 every output CONTRACTS.md names is declared in output_contract.py")
def e8():
    # The prose contract once listed outputs no code wrote (pair_evidence_nb08_v2.csv,
    # NB05's transfer tables) while every run reported success. output_contract.py is
    # the executable version run_pipeline enforces; this keeps the two in step.
    import fnmatch
    import output_contract as oc
    import run_pipeline as rp
    assert set(oc.CONTRACT_OUTPUTS) == set(rp.ORDER), (
        "output_contract.py and run_pipeline.py disagree on the notebooks: "
        f"{sorted(set(oc.CONTRACT_OUTPUTS) ^ set(rp.ORDER))}")
    overlap = [f"{o.path} ({nb}) also declared by {oc.owners(o.path.replace('*', 'X'))}"
               for nb, outs in oc.CONTRACT_OUTPUTS.items() for o in outs
               if oc.owners(o.path.replace("*", "X")) != [nb]]
    assert not overlap, "contract patterns overlap between notebooks:\n  " + "\n  ".join(overlap)

    declared = ([o.path.rsplit("/", 1)[-1] for outs in oc.CONTRACT_OUTPUTS.values()
                 for o in outs] + list(oc.NOT_PRODUCED))
    txt = (NB_DIR / "CONTRACTS.md").read_text(encoding="utf-8")
    unknown = []
    for name in sorted(set(re.findall(r"`([^`\s]+\.(?:csv|pkl|png|pdf))`", txt))):
        if name.startswith("Data/") or name in _DOC_INPUTS:
            continue
        pat = re.sub(r"\{[^}]*\}|<[^>]*>", "*", name.rsplit("/", 1)[-1])
        if not any(fnmatch.fnmatchcase(pat.replace("*", "X"), d)
                   or fnmatch.fnmatchcase(d.replace("*", "X"), pat) for d in declared):
            unknown.append(name)
    assert not unknown, ("CONTRACTS.md names output(s) that no notebook's contract declares "
                         "(declare them, or list them in output_contract.NOT_PRODUCED with "
                         "the reason): " + ", ".join(unknown))


# utils functions a behavioural check exercises although no executed notebook may reach
# them -- each a guard on NOT calling the function.
_GUARD_ONLY = {
    "load_agora2_ec_index": "D5 asserts it refuses while AGORA2 is disabled",
}


@check("E9 every utils function a behavioural check guards is reached by an executed notebook")
def e9():
    # The CLR transfer fixes guarded by A8/A9 once protected helpers that no executed
    # notebook called; the pipeline ran a private copy instead. A check on code the
    # pipeline never runs is a guard on dead code.
    import ast
    import run_pipeline as rp

    def names(node, pool):
        found = set()
        for n in ast.walk(node):
            if isinstance(n, ast.Name) and n.id in pool:
                found.add(n.id)
            elif (isinstance(n, ast.Attribute) and n.attr in pool
                  and isinstance(n.value, ast.Name) and n.value.id == "utils"):
                found.add(n.attr)
        return found

    utree = ast.parse((NB_DIR / "utils.py").read_text(encoding="utf-8"))
    funcs = {n.name: n for n in utree.body if isinstance(n, ast.FunctionDef)}
    calls = {f: names(node, funcs) - {f} for f, node in funcs.items()}

    roots = set()
    for nb in rp.ORDER:
        doc = json.loads((NB_DIR / nb).read_text(encoding="utf-8"))
        src = "\n".join("".join(c["source"]) for c in doc["cells"] if c["cell_type"] == "code")
        roots |= {m for m in re.findall(r"\butils\.(\w+)", src) if m in funcs}
        for block in re.findall(r"from utils import \(([^)]*)\)|from utils import ([^\n(]+)", src):
            roots |= {w for w in re.findall(r"\w+", " ".join(block)) if w in funcs}
    roots |= names(ast.parse((NB_DIR / "viz.py").read_text(encoding="utf-8")), funcs)
    reached, todo = set(), list(roots)
    while todo:
        f = todo.pop()
        if f not in reached:
            reached.add(f)
            todo.extend(calls.get(f, ()))

    ctree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    helpers = {n.name: n for n in ctree.body if isinstance(n, ast.FunctionDef)}
    guarded = {}
    for name, fn in _results:
        if name[:1] not in "ABCD":
            continue
        seen, todo, used = set(), [fn.__name__], set()
        while todo:
            h = todo.pop()
            if h in seen or h not in helpers:
                continue
            seen.add(h)
            used |= names(helpers[h], funcs)
            todo.extend(n.id for n in ast.walk(helpers[h])
                        if isinstance(n, ast.Name) and n.id in helpers)
        for f in used:
            guarded.setdefault(f, []).append(name.split()[0])
    dead = {f: chk for f, chk in guarded.items() if f not in reached and f not in _GUARD_ONLY}
    assert not dead, ("checks guard utils functions no executed notebook reaches (guards on "
                      "dead code): " + ", ".join(f"{f} [{'/'.join(sorted(set(c)))}]"
                                                 for f, c in sorted(dead.items())))


# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------

def main(argv=None):
    argv = argv or sys.argv[1:]
    only = [a for a in argv if not a.startswith("-")]
    selected = [(n, f) for n, f in _results
                if not only or any(n.startswith(o) for o in only)]

    passed, failed = 0, []
    print("=" * 78)
    print(f"v2.1 pre-flight checks  ({len(selected)} checks)")
    print("=" * 78)
    for name, fn in selected:
        try:
            fn()
        except Exception as exc:                       # noqa: BLE001 - report, don't raise
            failed.append((name, exc))
            print(f"FAIL  {name}")
            detail = str(exc).strip() or exc.__class__.__name__
            for line in detail.splitlines():
                print(f"      {line}")
            if not isinstance(exc, AssertionError):
                print("      " + traceback.format_exc().splitlines()[-3].strip())
        else:
            passed += 1
            print(f"ok    {name}")

    print("-" * 78)
    print(f"{passed} passed, {len(failed)} failed")
    if failed:
        print("\nPipeline must not run until these are fixed:")
        for name, _ in failed:
            print(f"  - {name}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
