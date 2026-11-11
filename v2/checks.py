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

@check("A1 nothing is dropped at load; HS and the adenoma arms are quarantined")
def a1():
    """v2.2: the adenoma arms are LOADED and quarantined rather than excluded.

    The dangerous way to re-admit them would be to leave them unmapped, because
    primary_mask() excludes rather than enumerates -- an unmapped group counts as
    primary (see A7). So the invariant is no longer "MP has no Stage.Group" but
    "MP has a Stage.Group that is NaN on the primary ordinal scale and is listed
    in NON_PRIMARY_STAGE_GROUPS".
    """
    assert utils.EXCLUDED_STUDY_GROUPS == {}, \
        f"nothing should be dropped at load, got {utils.EXCLUDED_STUDY_GROUPS}"
    assert utils.SENSITIVITY_GROUPS["YACHIDA-CRC-2019"] == ["HS"]
    assert utils.STAGE_GROUP_MAP["HS"] == "PostSurgery"

    # The adenoma arms are declared, mapped, and off the primary axis.
    assert utils.ADENOMA_GROUPS["YACHIDA-CRC-2019"] == ["MP"]
    assert utils.ADENOMA_GROUPS["Kim_adenomas_2020"] == ["Adenoma"]
    assert utils.STAGE_GROUP_MAP["MP"] == "Adenoma"
    assert utils.STAGE_GROUP_MAP["Adenoma"] == "Adenoma"
    assert "Adenoma" not in utils.STAGE_GROUP_ORDER, \
        "Adenoma must not be on the primary ordered axis"
    assert np.isnan(utils.STAGE_GROUP_NUMERIC["Adenoma"]), \
        "Adenoma must be NaN on the PRIMARY scale so it cannot enter a primary trend test"
    assert np.isnan(utils.CRC_STAGE_NUMERIC["MP"]), \
        "MP has no TNM stage; it must be NaN on the carcinoma-staging scale"
    assert "Adenoma" in utils.NON_PRIMARY_STAGE_GROUPS, \
        "Adenoma must be excluded by primary_mask"
    # ...and it is NOT filed as a sensitivity group: it is its own axis.
    assert "Adenoma" not in utils.SENSITIVITY_STAGE_GROUPS


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


@check("A3 real counts: 347 loaded, 277 primary, 317 early, 40 adenoma, 30 PS, 150 CRC")
def a3():
    y = _yachida_meta()
    assert len(y) == 347, f"expected 347 rows loaded (nothing dropped), got {len(y)}"
    counts = y["Stage.Group"].value_counts().to_dict()
    assert counts.get("Healthy") == 127, counts
    assert counts.get("Early_CRC") == 96, counts
    assert counts.get("Advanced_CRC") == 54, counts
    assert counts.get("Adenoma") == 40, counts
    assert counts.get("PostSurgery") == 30, counts

    pm = utils.primary_mask(y)
    sm = utils.sensitivity_mask(y)
    am = utils.adenoma_mask(y)
    em = utils.early_axis_mask(y)
    # THE invariant of the v2.2 change: admitting the adenoma arm must not move
    # the primary set by a single sample.
    assert int(pm.sum()) == 277, f"primary set moved to {int(pm.sum())}"
    assert int(sm.sum()) == 30
    assert int(am.sum()) == 40
    assert int(em.sum()) == 317, f"early axis is {int(em.sum())}, expected 317"
    assert int(utils.crc_only_mask(y).sum()) == 150

    # The three sets partition the cohort: disjoint, and jointly exhaustive.
    assert not (pm & sm).any() and not (pm & am).any() and not (sm & am).any(), \
        "primary / sensitivity / adenoma overlap"
    assert (pm | sm | am).all(), "a loaded row belongs to no analysis set"
    # The early axis is exactly primary + adenoma, and excludes PostSurgery.
    assert (em == (pm | am)).all(), "early axis is not primary + adenoma"
    assert not (em & sm).any(), "PostSurgery leaked onto the early axis"
    # No adenoma may carry a primary ordinal code.
    assert y.loc[am, "Stage.Group"].map(utils.STAGE_GROUP_NUMERIC).isna().all()
    assert y.loc[am, "CRC.Stage.Numeric"].isna().all()


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
    # (n_loaded, n_primary, n_cases_in_primary). v2.2: Kim loads its adenoma arm
    # too, so 240 rows come in and 138 of them are primary.
    expected = {"SINHA_CRC_2016": (131, 131, 42), "Kim_adenomas_2020": (240, 138, 36)}
    for ds, (n_loaded, n_primary, n_cases) in expected.items():
        raw = pd.read_csv(utils.DATA_DIR / f"{ds} metadata.tsv", sep="\t")
        drop = utils.EXCLUDED_STUDY_GROUPS.get(ds, [])
        raw = raw[~raw["Study.Group"].isin(drop)].reset_index(drop=True)
        meta = utils.harmonize_metadata(raw, ds)
        assert len(meta) == n_loaded, f"{ds}: loaded {len(meta)}, expected {n_loaded}"
        pm = utils.primary_mask(meta)
        assert int(pm.sum()) == n_primary, \
            f"{ds}: primary_mask kept {int(pm.sum())} of {n_loaded} samples, expected {n_primary}"
        kept_cases = int((meta.loc[pm, "Study.Group"] == "CRC").sum())
        assert kept_cases == n_cases, \
            f"{ds}: primary set holds {kept_cases} cases, expected {n_cases}"
        # No adenoma may reach a primary analysis through the back door.
        assert not (meta.loc[pm, "Study.Group"] == "Adenoma").any(), \
            f"{ds}: an adenoma survived primary_mask"
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
    """The opt-in adenoma analysis must not be able to disturb the primary set.

    v2.2: MP now HAS a Stage.Group ("Adenoma"), which is what A1 requires -- an
    unmapped group would count as primary. What keeps adenomas out of the primary
    set is NON_PRIMARY_STAGE_GROUPS plus the NaN primary ordinal code, both
    asserted in A1/A3. This check keeps its original subject: the severity scale.
    """
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
    # v2.2: the default no longer DROPS them -- it quarantines them. What must hold
    # is that the default load keeps the primary set untouched.
    default = utils.harmonize_metadata(raw, "Kim_adenomas_2020")
    assert int((default["Study.Group"] == "Adenoma").sum()) == 102, \
        "Kim's adenomas must be loaded under the default policy, not dropped"
    assert int(utils.primary_mask(default).sum()) == 138, \
        "Kim's primary set must stay 138 with the adenomas quarantined"


@check("A11 the early axis is pre-specified, ordered, and its own FDR family")
def a11():
    """The early-detection axis (v2.2). It must be ordered, fully mapped on every
    cohort that carries it, never merged with the primary comparison family, and
    never claimed for a cohort that has no pre-malignant arm."""
    assert utils.EARLY_AXIS_ORDER == ["Healthy", "Adenoma", "Early_CRC", "Advanced_CRC"]
    assert utils.EARLY_AXIS_ORDER_EVAL == ["Healthy", "Adenoma", "CRC"]
    codes = [utils.EARLY_AXIS_NUMERIC[g] for g in utils.EARLY_AXIS_ORDER]
    assert codes == sorted(codes) and len(set(codes)) == len(codes), codes
    assert np.isnan(utils.EARLY_AXIS_NUMERIC["PostSurgery"]), \
        "PostSurgery must have no position on the early axis either"
    # Yachida's two carcinoma arms must collapse onto one cross-cohort level.
    assert (utils.EARLY_AXIS_NUMERIC_EVAL["Early_CRC"]
            == utils.EARLY_AXIS_NUMERIC_EVAL["Advanced_CRC"]
            == utils.EARLY_AXIS_NUMERIC_EVAL["CRC"]), \
        "the eval axis must merge Early/Advanced into CRC -- the finest level Kim can confirm"
    assert utils.EARLY_AXIS_NUMERIC_EVAL["MP"] == utils.EARLY_AXIS_NUMERIC_EVAL["Adenoma"]

    # Sinha has no adenoma arm: claiming it here would yield a silently empty result.
    assert "SINHA_CRC_2016" not in utils.DATASETS_EARLY_AXIS, \
        "Sinha has no pre-malignant arm and must not be scored on the early axis"
    assert set(utils.DATASETS_EARLY_AXIS) == {"YACHIDA-CRC-2019", "Kim_adenomas_2020"}

    # Its own FDR family: no early comparison may appear in the primary family.
    primary_pairs = {p for ps in utils.COMPARISONS_GROUP.values() for p in ps}
    early_pairs = {p for ps in utils.COMPARISONS_EARLY.values() for p in ps}
    assert not (primary_pairs & early_pairs), \
        f"these contrasts sit in both FDR families: {primary_pairs & early_pairs}"
    # Every early comparison must involve the adenoma arm -- that is the point.
    for ds, comps in utils.COMPARISONS_EARLY.items():
        assert ds in utils.DATASETS_EARLY_AXIS, ds
        for g1, g2 in comps:
            assert "Adenoma" in (g1, g2), f"{ds}: {g1} vs {g2} does not test the adenoma arm"

    # Each cohort's labels must be fully mapped in the column the axis reads.
    for ds, col in utils.EARLY_AXIS_GROUP_COL.items():
        raw = pd.read_csv(utils.DATA_DIR / f"{ds} metadata.tsv", sep="\t")
        meta = utils.harmonize_metadata(raw, ds)
        sub = meta.loc[utils.early_axis_mask(meta), col]
        amap = (utils.EARLY_AXIS_NUMERIC if ds == utils.DATASET_DISCOVERY
                else utils.EARLY_AXIS_NUMERIC_EVAL)
        unmapped = sorted(set(sub[sub.map(amap).isna()].astype(str)))
        assert not unmapped, f"{ds}: early-axis labels with no ordinal code: {unmapped}"
        # Each arm of each declared comparison must be large enough to test.
        for g1, g2 in utils.COMPARISONS_EARLY[ds]:
            for g in (g1, g2):
                n = int((sub == g).sum())
                assert n >= 3, f"{ds}: early-axis group {g!r} has {n} samples"


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


@check("C10 NB03's adjusted target carries no confounder signal and keeps a species effect")
def c10():
    """NB03 used to keep the confounders as features and permute the label against
    all of them, testing 'target independent of species AND confounders'. The target
    is now residualised on the confounders (utils.residualise_on_confounders). Here a
    species tracks the confounders and the target depends on the confounders only, so
    the species must lose its association with the adjusted target; a planted direct
    species effect must survive."""
    rng = np.random.default_rng(10)
    n = 400
    Z = np.column_stack([rng.normal(size=n), rng.integers(0, 2, n), rng.gamma(2, 1, n)])
    s_conf = Z @ np.array([0.8, 0.6, 0.5]) + rng.normal(0, 0.6, n)    # confounded species
    s_true = rng.normal(size=n)                                        # causal species
    y_null = Z @ np.array([1.0, 0.8, 0.7]) + rng.normal(0, 1, n)      # confounders only
    y_alt = y_null + 0.6 * s_true

    r = utils.residualise_on_confounders(y_null, Z)
    design = np.column_stack([np.ones(n), Z])
    assert np.abs(design.T @ r).max() < 1e-8 * n, "residual is not orthogonal to [1, Z]"
    raw = abs(np.corrcoef(s_conf, y_null)[0, 1])
    adj = abs(np.corrcoef(s_conf, r)[0, 1])
    assert raw > 0.5 and adj < 0.1, \
        f"confounded species: |rho| raw {raw:.2f} -> adjusted {adj:.2f} (expected < 0.1)"
    kept = np.corrcoef(s_true, utils.residualise_on_confounders(y_alt, Z))[0, 1]
    assert kept > 0.4, f"a direct species effect did not survive adjustment (rho {kept:.2f})"
    # Rank-deficient Z (a duplicated and a constant column) gives the same residual.
    # (C11 below covers the permutation p-values these targets are tested with.)
    Zd = np.column_stack([Z, Z[:, 0], np.ones(n)])
    assert np.allclose(utils.residualise_on_confounders(y_null, Zd), r, atol=1e-8)
    try:
        utils.residualise_on_confounders(y_null, np.where(Z > 2.5, np.nan, Z))
    except ValueError:
        pass
    else:
        raise AssertionError("missing confounder values were accepted")


@check("C11 sequential permutation p-values are valid and spend draws where p is small")
def c11():
    """NB03's Besag-Clifford p-values (utils.sequential_permutation_p), checked with a
    cheap statistic whose null draws are fixed per index, like NB03's seeded permutations."""
    rng = np.random.default_rng(11)
    seq = utils.sequential_permutation_p

    def nulls(seed, n=400):
        v = np.random.default_rng(seed).standard_normal(n)
        return lambda a, b: v[a:b]

    # Valid under the null: P(p <= a) <= a, within Monte Carlo error.
    reps = 600
    ps = np.array([seq(rng.standard_normal(), nulls(1000 + r), h=10, b_max=200,
                       batch_size=16)[0] for r in range(reps)])
    for a in (0.05, 0.1, 0.25, 0.5):
        assert np.mean(ps <= a) <= a + 3 * np.sqrt(a * (1 - a) / reps), \
            f"P(p <= {a}) = {np.mean(ps <= a):.3f} under the null"
    # A strong effect never meets h exceedances: full budget, p at the floor.
    p, n, x, early, _ = seq(10.0, nulls(1), h=10, b_max=200)
    assert (n, x, early) == (200, 0, False) and abs(p - 1 / 201) < 1e-15
    # A null-like effect stops early, at p = h / L exactly.
    p, n, x, early, _ = seq(0.0, nulls(2), h=10, b_max=200)
    assert early and x == 10 and n < 60 and abs(p - 10 / n) < 1e-15
    # The answer does not depend on how the draws are batched.
    assert seq(0.8, nulls(3), h=10, b_max=300, batch_size=7)[:4] == \
        seq(0.8, nulls(3), h=10, b_max=300, batch_size=64)[:4]
    # A failed null fit (non-finite) counts against the observed statistic.
    p, n, x, early, _ = seq(10.0, lambda a, b: np.full(b - a, np.nan), h=3, b_max=50)
    assert early and n == 3 and p == 1.0


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


@check("D5 AGORA2 evidence comes from the prebuilt index, which matches the reconstructions")
def d5():
    # Disabled 2026-10-07 because building the index (a 1.8 GB scan) from inside a run
    # overloaded the old server; re-enabled 2026-10-21 with the build moved to
    # tools/build_agora2_index.py. Notebooks may LOAD the index, never build it.
    assert utils.AGORA2_ENABLED is True, "AGORA2 evidence is disabled in utils"
    cache, side = utils.AGORA2_INDEX_CACHE, utils.AGORA2_INDEX_SIDECAR
    assert cache.exists() and side.exists(), (
        f"AGORA2 index or its build record is missing ({cache}); build it once with "
        "python tools/build_agora2_index.py")
    built = json.loads(side.read_text(encoding="utf-8"))
    try:
        src = utils.find_agora2_data_dir()
    except FileNotFoundError:
        src = None          # reconstructions removed after the build: the index stands alone
    if src is not None:
        now = utils.agora2_source_fingerprint(src)
        assert all(built.get(k) == v for k, v in now.items()), (
            "the AGORA2 index was built from different reconstructions than those on "
            "disk; rebuild: python tools/build_agora2_index.py --force")
    banned = r"build_agora2_ec_index|find_agora2_data_dir"
    hits = [f"{name} cell {idx} line {k}"
            for name, idx, ctype, text in _notebook_sources() if ctype == "code"
            for k, line in enumerate(text.splitlines(), 1)
            if not line.strip().startswith("#") and re.search(banned, line)]
    assert not hits, ("notebook(s) build the AGORA2 index (the scan that overloaded the "
                      "old server); load it with utils.load_agora2_ec_index: "
                      + ", ".join(hits))


@check("D7 AGORA2 lookups are tri-state, match names safely and index producing reactions only")
def d7():
    import tempfile
    idx = {"s__Alpha_beta": {"C00134": ["4.1.1.17"]},
           "s__Gamma_D_delta": {"C00134": ["1.1.1.1"]},
           "s__Eps_A_zeta": {}, "s__Eps_B_zeta": {}}
    sid = utils.agora2_species_id
    assert sid("Alpha beta", idx) == ("s__Alpha_beta", "exact")
    assert sid("Gamma delta", idx) == ("s__Gamma_D_delta", "suffix_normalised")
    assert sid("Eps zeta", idx) == (None, None), "an ambiguous suffix match was guessed"
    enc = utils.species_encodes_enzyme
    assert enc("Alpha beta", "C00134", idx) is True
    assert enc("Gamma delta", "C00134", idx) is False, "a non-producing EC counted"
    assert enc("Unknown sp", "C00134", idx) is None, "no model must be None, not False"
    assert utils.species_has_any_ec("Unknown sp", ["4.1.1.17"], idx) is None
    rate, n = utils.agora2_base_rate(["Alpha beta", "Gamma delta", "Unknown sp"], "C00134", idx)
    assert (rate, n) == (0.5, 2), "base rate must leave out species without a model"
    # The builder indexes an EC under a compound only when the reaction PRODUCES it.
    with tempfile.TemporaryDirectory() as d:
        model = {"metabolites": [{"id": "put_c", "annotation": {"kegg.compound": "C00134"}},
                                 {"id": "orn_c", "annotation": {"kegg.compound": "C00077"}}],
                 "reactions": [{"annotation": {"ec-code": "['4.1.1.17']"},
                                "metabolites": {"orn_c": -1, "put_c": 1}}]}
        Path(d, "s__Toy_species.json").write_text(json.dumps(model), encoding="utf-8")
        built = utils.build_agora2_ec_index(Path(d))
    assert built == {"s__Toy_species": {"C00134": ["4.1.1.17"]}}, built


@check("D8 MICOM helpers: taxon ids map back, fluxes are abundance-weighted, tradeoff rule, failed solves")
def d8():
    # Each assertion is a fact about micom 0.39.1 output that NB10 once got wrong.
    import micom_tools as mt
    # micom rewrites taxon ids (community.py); the names must map back, collisions refused
    assert mt.micom_taxon_id("CAG-110 sp003525905") == "CAG_110_sp003525905"
    assert mt.micom_taxon_id("Prevotella copri_A") == "Prevotella_copri_A"
    assert mt.taxon_name_map(["Dialister hominis", "CAG-110 sp003525905"]) == {
        "Dialister_hominis": "Dialister hominis", "CAG_110_sp003525905": "CAG-110 sp003525905"}
    try:
        mt.taxon_name_map(["Genus sp-1", "Genus sp 1"])
        raise AssertionError("two species sharing one micom id were accepted")
    except ValueError:
        pass
    # metabolite ids in every form micom writes
    for m in ("ptrc[e]", "ptrc_m", "EX_ptrc(e)", "EX_ptrc_m", "PTRC[e]"):
        assert mt.classify_metabolite(m) == ("C00134", "Putrescine"), m
    assert mt.classify_metabolite("arg_L[e]") == ("C00062", "Arginine")
    assert mt.classify_metabolite("co2[e]") == (None, None)
    # tradeoff: low values tie at "every taxon grows"; the rule takes the LARGEST tie
    sweep = pd.DataFrame({"tradeoff": [np.nan] * 2 + [0.1] * 2 + [0.5] * 2 + [0.9] * 2,
                          "growth_rate": [0, 0, .01, .01, .1, .1, .3, 0]})
    assert mt.choose_tradeoff(sweep)[0] == 0.5, "first tie (idxmax) picked, not the largest"
    # a two-taxon lumen: A (abundance .25) secretes putrescine at 2 per gram of A, B (.75)
    # takes up 1 per gram of B; the community imports the difference, .5 - .75 = -.25
    ex = pd.DataFrame({
        "taxon": ["A_x", "B_y", "medium", "A_x"], "sample_id": ["1"] * 4,
        "reaction": ["EX_ptrc(e)", "EX_ptrc(e)", "EX_ptrc_m", "EX_co2(e)"],
        "metabolite": ["ptrc[e]", "ptrc[e]", "ptrc_m", "co2[e]"],
        "flux": [2.0, -1.0, -0.25, 1.0], "abundance": [0.25, 0.75, np.nan, 0.25],
        "direction": ["export", "import", "import", "export"]})
    names = {"A_x": "A x", "B_y": "B y"}
    pa = mt.polyamine_exchanges(ex, names)
    assert set(pa["Species"]) == {"A x", "B y"}, "taxon ids were not mapped back to names"
    sp = mt.species_export(pa)
    assert sp[["Species", "export_community", "export_specific"]].values.tolist() == [["A x", 0.5, 2.0]], \
        "export must be abundance-weighted and secretion-only"
    bal = mt.community_balance(ex, pa).iloc[0]
    assert np.isclose(bal["Net_medium_flux"], bal["Gross_export"] - bal["Gross_uptake"]), \
        "lumen balance: medium flux = sum of abundance-weighted taxon fluxes"
    assert (bal["N_exporting_taxa"], bal["N_uptaking_taxa"]) == (1, 1)
    try:
        mt.polyamine_exchanges(ex, {"A_x": "A x"})
        raise AssertionError("a taxon without a species name was dropped silently")
    except KeyError:
        pass
    assert mt.missing_samples(["1", "2"], ex) == ["2"], "a failed solve went unnoticed"
    # micom returns an infeasible point without an error: zero growth, unbalanced lumen
    ok = pd.DataFrame({"taxon": ["A_x", "B_y", "medium"], "sample_id": "1",
                       "metabolite": ["ptrc[e]", "ptrc[e]", "ptrc_m"], "flux": [2.0, -1.0, -0.25],
                       "abundance": [.25, .75, np.nan], "tolerance": 1e-6})
    bad = ok.assign(sample_id="2", flux=[2.0, -1.0, 7.0])
    growth = pd.DataFrame({"taxon": ["A_x", "B_y"] * 2, "sample_id": ["1", "1", "2", "2"],
                           "abundance": [.25, .75, .25, .75], "growth_rate": [.4, .2, 0., 0.]})
    q = mt.solution_quality(pd.concat([ok, bad]), growth).set_index("sample_id")
    assert q.loc["1", "valid"] and not q.loc["2", "valid"], \
        "an infeasible solution (no growth, unbalanced lumen) passed as solved"


@check("D9 evidence_tools: tri-state calls, base-rate gate, attributability, source counting")
def d9():
    """The assembler's rules. Each assertion is a way the evidence synthesis has been, or
    could silently be, wrong."""
    import evidence_tools as et

    # 1. Genus derivation must use the SHARED normaliser, or a species-level pair will not
    #    join the genus-level replication table for any GTDB-suffixed genus.
    assert et.species_to_genus("Veillonella parvula_A") == "Veillonella"
    assert et.species_to_genus("Blautia_A sp900066165") == "Blautia"
    assert et.species_to_genus("CAG-110 sp003525905") == "CAG-110"
    assert et.species_to_genus("") == utils.normalise_genus_name("")

    # 2. AGORA2 capability is TRI-STATE: no model means None, never False.
    idx = {"s__Fake_species": {"C00134": ["4.1.1.17"], "C00077": []}}
    assert et.species_touches_compound("Fake species", "C00134", idx) is True
    assert et.species_touches_compound("Fake species", "C00077", idx) is False, \
        "a compound present with no reactions must be False, not True"
    assert et.species_touches_compound("Fake species", "C99999", idx) is False
    assert et.species_touches_compound("Absent species", "C00134", idx) is None, \
        "a species with no AGORA2 model must carry NO evidence, not negative evidence"

    # 3. Base rate excludes unmodelled species from the DENOMINATOR. Counting them as
    #    negatives would deflate the rate and make an uninformative source look selective.
    idx2 = {"s__A": {"C1": ["x"]}, "s__B": {"C1": ["x"]}, "s__C": {}}
    rate, n = et.compound_base_rate("C1", ["A", "B", "C", "Unmodelled"], idx2)
    assert (rate, n) == (2 / 3, 3), (rate, n)
    rate_none, n_none = et.compound_base_rate("C1", ["Unmodelled"], idx2)
    assert np.isnan(rate_none) and n_none == 0

    # 4. The base-rate gate: a source nearly every species satisfies earns no point, and
    #    an unassessable one (NaN) can never pass.
    assert et.counts_toward_tier(0.10) and et.counts_toward_tier(0.20)
    assert not et.counts_toward_tier(0.21)
    assert not et.counts_toward_tier(float("nan")), "NaN must not count as informative"
    assert not et.counts_toward_tier(None)
    # The measured polyamine base rates: only spermine is informative.
    measured = {"arginine": 0.974, "ornithine": 0.901, "putrescine": 0.874,
                "agmatine": 0.702, "cadaverine": 0.623, "spermidine": 0.576,
                "spermine": 0.000}
    informative = {k for k, v in measured.items() if et.counts_toward_tier(v)}
    assert informative == {"spermine"}, informative

    # 5. Attributability: the confounders alone must not match the species model.
    #    Measured on the 2026-10-26 run; both failures are polyamines.
    assert et.confounder_unattributable(0.633, 0.823) is True, "spermine must be excluded"
    assert et.confounder_unattributable(0.577, 0.741) is True, "N1-acetylspermine likewise"
    assert et.confounder_unattributable(0.787, 0.484) is False, "Phe must survive"
    assert et.confounder_unattributable(0.735, 0.602) is False, "Ser must survive"
    assert et.confounder_unattributable(0.70, None) is False, "missing data is not a verdict"

    # 6. Predictive support needs direction, rank stability, a model better than chance,
    #    AND attributability. Each omission is a way a non-finding becomes a tier.
    base = {"Direction": "positive", "Rank_stable": True, "Model_perm_q": 0.01,
            "Confounder_unattributable": False}
    assert et.predictive_support(base)
    assert not et.predictive_support({**base, "Direction": "negative"}), \
        "a negative-direction species is not a candidate producer"
    assert not et.predictive_support({**base, "Rank_stable": False})
    assert not et.predictive_support({**base, "Model_perm_q": 0.20}), \
        "a chance-level model's ranking is not evidence"
    assert not et.predictive_support({**base, "Confounder_unattributable": True}), \
        "a model whose advantage belongs to the confounders is not species evidence"
    assert not et.predictive_support({**base, "Model_perm_q": float("nan")})

    # 7. Replication keeps "not measurable" apart from "not replicated" (cf. C9) and
    #    reports a contradiction as its own outcome.
    cohorts = ["SINHA_CRC_2016", "Kim_adenomas_2020"]
    none_row = et.replication_summary(None, cohorts)
    assert not none_row["Replication_available"] and not none_row["Replication_support"]
    unmeasurable = et.replication_summary(
        {"In_SINHA_CRC_2016": False, "In_Kim_adenomas_2020": False}, cohorts)
    assert not unmeasurable["Replication_available"], \
        "an unmeasurable pair must be unavailable, never a failed replication"
    good = et.replication_summary(
        {"In_SINHA_CRC_2016": True, "Rho_SINHA_CRC_2016": 0.3,
         "Replicated_SINHA_CRC_2016": True, "Contradicted_SINHA_CRC_2016": False,
         "In_Kim_adenomas_2020": True, "Rho_Kim_adenomas_2020": 0.4,
         "Replicated_Kim_adenomas_2020": False,
         "Contradicted_Kim_adenomas_2020": False}, cohorts)
    assert good["Replication_available"] and good["Replication_support"]
    assert good["Replication_n_measurable"] == 2 and good["Replication_n_replicated"] == 1
    assert good["Replication_cohorts"] == "SINHA_CRC_2016"
    assert np.isclose(good["Replication_mean_rho"], 0.35)
    contra = et.replication_summary(
        {"In_Kim_adenomas_2020": True, "Rho_Kim_adenomas_2020": -0.4,
         "Replicated_Kim_adenomas_2020": False,
         "Contradicted_Kim_adenomas_2020": True}, cohorts)
    assert contra["Replication_contradicted"] and not contra["Replication_support"]

    # 8. Source counting: an unavailable source can never supply a point, however its
    #    support column reads.
    srcs = [("predictive", "P_avail", "P_sup"), ("genomic", "G_avail", "G_sup"),
            ("replication", "R_avail", "R_sup")]
    n_av, n_sup, names = et.count_sources(
        {"P_avail": True, "P_sup": True, "G_avail": False, "G_sup": True,
         "R_avail": True, "R_sup": False}, srcs)
    assert (n_av, n_sup) == (2, 1), (n_av, n_sup)
    assert names == "predictive", names

    # 9. Tiers state the attainable ceiling, because 1-of-1 is not the claim 1-of-4 is.
    assert et.tier(3, 4).startswith("Strong")
    assert et.tier(2, 4).startswith("Moderate")
    assert et.tier(1, 4).startswith("Weak")
    assert et.tier(0, 4).startswith("Statistical only")
    assert "1/4" in et.tier(1, 4) and "1/1" in et.tier(1, 1)


@check("D12 evidence classes (CONTRACTS 5d): a class counts once, genus replication never votes")
def d12():
    import evidence_tools as et

    def row(**kw):
        r = {}
        for src, (avail, sup) in kw.items():
            r[f"Source_{src}_available"] = avail
            r[f"Source_{src}_support"] = sup
        return r

    # 1. Every member of every class is classifiable; a non-voting source has no class;
    #    an unknown source raises rather than slipping into a tier.
    for cls, members in et.EVIDENCE_CLASSES.items():
        for m in members:
            assert et.source_class(m) == cls
    assert et.source_class("replication") is None
    try:
        et.source_class("brand_new_source")
        raise AssertionError("an unclassified source must raise")
    except KeyError:
        pass

    # 2. The worked example that motivated 5d: NB03 + MOFA + genus replication is ONE class.
    ramosum = row(predictive=(True, True), mofa_coloading=(True, True),
                  replication=(True, True), genomic_capability=(False, False))
    av, sup = et.class_support(ramosum)
    assert sup == ["O"], sup
    assert et.class_tier(ramosum).startswith("Observational only"), et.class_tier(ramosum)

    # 3. Three classes without G or M is not Conclusive.
    no_direct = row(predictive=(True, True), genomic_capability=(True, True),
                    replication_species=(True, True))
    assert et.class_tier(no_direct).startswith("Supported"), et.class_tier(no_direct)

    # 4. Three classes including G is Conclusive -- unless contradicted, unattributable,
    #    or the target failed the export gate.
    full = row(predictive=(True, True), genomic_capability=(True, True),
               gene_level=(True, True))
    assert et.class_tier(full).startswith("Conclusive"), et.class_tier(full)
    assert not et.class_tier({**full, "Replication_species_contradicted": True}).startswith("Conclusive")
    assert not et.class_tier({**full, "Confounder_unattributable": True}).startswith("Conclusive")
    assert not et.class_tier({**full, "Target_export_gate": False}).startswith("Conclusive")

    # 5. An unavailable member cannot supply a point, whatever its support column says.
    stray = row(gene_level=(False, True), predictive=(True, True))
    assert et.class_support(stray)[1] == ["O"]

    # 6. Two members of one class still count once.
    double_r = row(genomic_producer=(True, True), mimosa2=(True, True))
    assert et.class_support(double_r)[1] == ["R"]
    assert et.class_tier(double_r).startswith("Single-class")

    # 7. The label states the classes available.
    assert "classes" in et.class_tier(full) and "OR" in et.class_tier(full).replace("G", "")


@check("D13 the target registry labels every target and states its gates (CONTRACTS 5d-ii)")
def d13():
    reg_path = NB_DIR / "targets" / "out" / "target_registry.tsv"
    if not reg_path.exists():
        print("        (no registry yet: run v2/targets/build_registry.py after NB03)")
        return
    reg = pd.read_csv(reg_path, sep="\t")
    need = {"target", "kegg", "label", "export_gate", "origin_flag",
            "confounder_unattributable", "attribution_requires_lysis"}
    assert need <= set(reg.columns), need - set(reg.columns)
    assert reg["target"].is_unique, "one row per target"
    assert set(reg["label"]) <= {"stage_associated", "a_priori_class"}, set(reg["label"])
    assert set(reg["origin_flag"]) <= {"microbial_only", "host_and_microbial", "host_dominant",
                                       "unannotated", "unassessable"}, set(reg["origin_flag"])
    # A stage-associated target must say in which contrast it is robust.
    sa = reg[reg["label"] == "stage_associated"]
    assert (sa[["robust_early", "robust_advanced"]].fillna("").ne("").any(axis=1)).all(), \
        "a stage_associated target with no robust contrast"
    # The gate and the lysis flag are complements wherever the gate is assessable.
    g = reg.dropna(subset=["export_gate"])
    assert (g["export_gate"].astype(bool) != g["attribution_requires_lysis"].astype(bool)).all()


@check("D10 the decarboxylation target class is pre-specified, paired, and producing-only")
def d10():
    """The v2.2 target class. Its validity rests on being chosen from an enzyme family
    rather than from this cohort's statistics -- E14 enforces that NB03 cannot read the
    correlation table -- and on substrate/product being paired, which is what makes it a
    directional signature rather than a bare association."""
    sp = utils.DECARBOXYLATION_SUBSTRATE_PRODUCT
    assert len(sp) == 8, f"expected 8 substrate->product pairs, got {len(sp)}"

    # Every pair must be substrate -> product with a producing EC, and no compound may be
    # both a substrate and a product of this family (that would double-count a target).
    substrates, products = set(sp), {p for p, _e, _ec in sp.values()}
    assert not (substrates & products), \
        f"compound is both substrate and product: {substrates & products}"
    assert products == set(utils.DECARBOXYLATION_KEGG), \
        "DECARBOXYLATION_KEGG must be exactly the products of the pairs"
    assert set(utils.DECARBOXYLATION_EC) == products

    # Every EC must be a decarboxylase (EC 4.1.1.x). A transferase or a consuming enzyme
    # here would invert the producer call, which is the pre-v2 polyamine-map error.
    for prod, ecs in utils.DECARBOXYLATION_EC.items():
        assert ecs, prod
        for ec in ecs:
            assert ec.startswith("4.1.1."), \
                f"{prod}: EC {ec} is not a decarboxylase (4.1.1.x)"

    # The polyamine products that ARE decarboxylations must agree with POLYAMINE_EC, and
    # the aminopropyl-transfer products must NOT be in this class.
    for shared in ("C00134", "C01672", "C00179"):
        assert shared in utils.DECARBOXYLATION_EC and shared in utils.POLYAMINE_EC
        assert set(utils.DECARBOXYLATION_EC[shared]) <= set(utils.POLYAMINE_EC[shared]), \
            f"{shared}: the decarboxylase must already be a curated producing EC"
    for transfer in ("C00315", "C00750"):
        assert transfer not in utils.DECARBOXYLATION_KEGG, \
            f"{transfer} is an aminopropyl-transfer product, not a decarboxylation"

    # The arginine-route intermediate is kept apart: its enzyme is a transferase.
    assert set(utils.ARGININE_ROUTE_EC) == {"C00624"}
    assert utils.ARGININE_ROUTE_EC["C00624"] == ["2.3.1.1"]
    assert not (set(utils.ARGININE_ROUTE_KEGG) & products)

    # The assembled target class: every substrate, every product, the intermediate; each
    # tagged with its role so a table can say why a target is in the list.
    cls = utils.DECARBOXYLATION_TARGET_KEGG
    assert len(cls) == 8 + 8 + 1 == 17, f"expected 17 class members, got {len(cls)}"
    roles = {r for _n, r in cls.values()}
    assert roles == {"substrate", "product", "arginine_route"}, roles
    for kegg, (name, role) in cls.items():
        assert re.fullmatch(r"C\d{5}", kegg), kegg
        assert name and isinstance(name, str), kegg
    assert {k for k, (_n, r) in cls.items() if r == "substrate"} == substrates
    assert {k for k, (_n, r) in cls.items() if r == "product"} == products


@check("D11 the sensitivity engine refuses unreachable effect sizes and scores the real rule")
def d11():
    """Step 6D turns "0 robust adenoma hits" into a bounded statement, so its two ways of
    being silently wrong are guarded here.

    First, Cliff's delta is a STEP function of a location shift wherever values are tied.
    A metabolite sitting mostly at the detection floor jumps straight past a small target
    to a large effect, and scoring that case reports power for an effect size that was
    never applied -- it showed up as 26% "power" at delta=0.05, which is impossible for a
    near-null effect. evaluate_shift must return None there instead.

    Second, the criterion scored must be the criterion the early axis is scored by. If
    robust_da_call drifted from merge_da's rule, the analysis would bound the wrong test.
    """
    import numpy as _np

    rng = _np.random.default_rng(11)

    # --- the rule: both significant AND same direction, with NaN never counting
    assert utils.robust_da_call(0.01, 0.01, 1.0, 2.0) is True
    assert utils.robust_da_call(0.01, 0.01, 1.0, -2.0) is False, "opposite signs"
    assert utils.robust_da_call(0.20, 0.01, 1.0, 2.0) is False, "rank test not significant"
    assert utils.robust_da_call(0.01, 0.20, 1.0, 2.0) is False, "adjusted not significant"
    assert utils.robust_da_call(_np.nan, 0.01, 1.0, 2.0) is False, "NaN must not count"
    assert utils.robust_da_call(0.01, _np.nan, 1.0, 2.0) is False, "NaN must not count"
    # the threshold is the project's, not a local 0.05
    assert utils.robust_da_call(0.04, 0.04, 1.0, 1.0) is (utils.FDR_THRESHOLD > 0.04)

    # --- Cliff's delta sign convention: positive means group1 exceeds group2, which is
    # the OPPOSITE sign to MeanDiff. Getting this backwards would invert every curve.
    assert utils.cliffs_delta(_np.array([5.0, 6, 7]), _np.array([1.0, 2, 3])) == 1.0
    assert utils.cliffs_delta(_np.array([1.0, 2, 3]), _np.array([5.0, 6, 7])) == -1.0

    # --- a continuous feature: the requested delta is reachable and is actually hit
    x = rng.normal(size=60)
    y = rng.normal(size=40)
    for target in (-0.5, -0.2, 0.3):
        shift, achieved = utils.shift_for_cliffs_delta(x, y, target)
        assert _np.isfinite(shift), f"continuous data: {target} should be reachable"
        assert abs(achieved - target) <= utils.SENSITIVITY_DELTA_TOL, \
            f"continuous data: asked {target}, achieved {achieved}"

    # --- a floor-heavy feature: most values tied at the floor, so delta moves in large
    # steps and a small target cannot be realised by shifting.
    floor = -3.0
    y_tied = _np.concatenate([_np.full(34, floor), rng.normal(2.0, 0.3, 6)])
    x_tied = _np.concatenate([_np.full(100, floor), rng.normal(2.0, 0.3, 27)])
    _, achieved = utils.shift_for_cliffs_delta(x_tied, y_tied, -0.02)
    assert not _np.isfinite(achieved) or abs(achieved - (-0.02)) > utils.SENSITIVITY_DELTA_TOL, \
        ("this tied feature must NOT be able to realise delta=-0.02; if it can, the "
         "fixture no longer exercises the step-function trap")

    # --- the gate that uses it, and the prevalence model's applicability rule
    assert utils.SENSITIVITY_DELTA_TOL > 0
    assert utils.SENSITIVITY_MIN_FLOOR >= 3 and utils.SENSITIVITY_MIN_DETECTED >= 3

    # --- the engine must expose the proof, not just the shortcut
    for name in ("evaluate_shift", "evaluate_prevalence", "prevalence_choice",
                 "verify_against_pipeline", "_verify_baseline"):
        assert hasattr(utils.RobustDASensitivity, name), f"missing {name}"


@check("D6 the curated polyamine and SCFA EC maps cover every compound and guild")
def d6():
    # The curation was previously WRONG (consuming enzymes listed as producers); D1
    # guards its content, this its coverage.
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


@check("E1 no notebook USES a removed construct (SMOTE / TumorLoc_ / Precursor / a pooled adenoma)")
def e1():
    """v2.2 narrowed this check.

    It used to ban a quoted 'Adenoma' or 'MP' anywhere in code, because under v2.1
    those arms were dropped at load and any mention of them as data was stale. They
    are now the early-detection axis (CONTRACTS section 2b), so the blanket ban
    would forbid the correct code.

    What is still banned is the thing the old check was really aiming at: an adenoma
    POOLED INTO a primary disease group, which is what the pre-v2 pipeline did
    (MP -> Advanced_CRC). The stronger guarantees are in code, not in this grep:
    A1 (Adenoma is NaN on the primary ordinal scale and sits in
    NON_PRIMARY_STAGE_GROUPS), A3 (the three analysis sets partition the cohort and
    the primary set is still 277) and A11 (the early axis is its own FDR family).
    """
    # Banned everywhere: never a legitimate use.
    banned_any = {
        "SMOTE":     r"\bSMOTE\b",
        "TumorLoc_": r"TumorLoc_",     # the leaking dummy prefix; TumorSite_ replaced it
        "Precursor": r"\bPrecursor\b",
    }
    # Banned in CODE only: an adenoma arm mapped onto a PRIMARY group label. That is
    # the pre-v2 grouping (MP -> Advanced_CRC, HS -> Early_CRC) and it is what made
    # the old results uninterpretable.
    banned_code = {
        "an adenoma arm pooled into a primary group":
            r"""['"](?:MP|Adenoma)['"]\s*:\s*['"](?:Healthy|Early_CRC|Advanced_CRC)['"]""",
        "HS pooled into a primary group":
            r"""['"]HS['"]\s*:\s*['"](?:Healthy|Early_CRC|Advanced_CRC)['"]""",
        "an off-axis group appended to the primary ordered axis":
            r"""STAGE_GROUP_ORDER\s*\+\s*\[\s*['"](?:Adenoma|MP|PostSurgery)['"]""",
    }
    hits = []
    for name, idx, ctype, text in _notebook_sources():
        patterns = dict(banned_any)
        if ctype == "code":
            patterns.update(banned_code)
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
    "manifest.csv": "the MICOM model manifest of the unused pre-built set in Data/micom/models_agora2/",
    "western_diet_gut.csv": "the source MICOM medium in Data/micom/, read by diet_screen.py and NB10",
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
    "build_agora2_ec_index": "run only by tools/build_agora2_index.py; D5 forbids notebooks",
    "find_agora2_data_dir": "run only by tools/build_agora2_index.py; D5 forbids notebooks",
    "agora2_source_fingerprint": "the build tool and D5 compare it; no notebook needs it",
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
    dead = {f: chk for f, chk in guarded.items()
            if f not in reached and f not in _GUARD_ONLY}
    assert not dead, ("checks guard utils functions no executed notebook reaches (guards on "
                      "dead code): " + ", ".join(f"{f} [{'/'.join(sorted(set(c)))}]"
                                                 for f, c in sorted(dead.items())))


@check("E10 every plotting notebook installs viz.apply_theme and none restyles on its own")
def e10():
    # Five notebooks styled themselves with seaborn (system Arial first on Windows) and
    # three applied no theme at all, so one run rendered two typefaces and subscript
    # labels lost their glyphs.
    plots = re.compile(r"\bplt\.|\bsns\.|savefig\(|save_figure\(")
    restyle = re.compile(r"sns\.set_style\(|sns\.set_theme\(|sns\.set\(|sns\.set_context\(|"
                         r"plt\.style\.use\(|mpl\.rc\(|rcParams\[[^\]]+\]\s*=|"
                         r"rcParams\.update\(")
    hits = []
    for nb in sorted(NB_DIR.glob("*.ipynb")):
        doc = json.loads(nb.read_text(encoding="utf-8"))
        src = "\n".join("".join(c["source"]) for c in doc["cells"] if c["cell_type"] == "code")
        if plots.search(src) and "viz.apply_theme(" not in src:
            hits.append(f"{nb.name}: draws figures without viz.apply_theme()")
        hits += [f"{nb.name} line {k}: {line.strip()[:70]}"
                 for k, line in enumerate(src.splitlines(), 1)
                 if restyle.search(line) and not line.strip().startswith("#")]
    assert not hits, "figure styling outside viz.apply_theme:\n  " + "\n  ".join(hits)


@check("E12 figures are written only through viz.save_figure (300-dpi PNG + PDF twin)")
def e12():
    # 42 call sites typed their own savefig arguments: 150-dpi PNG only, no vector copy,
    # most without the theme's facecolor. viz.save_figure is the one writer; the output
    # contract then requires the PDF twin of every declared PNG.
    rx = re.compile(r"\.savefig\(")
    hits = [f"{name} cell {idx} line {k}"
            for name, idx, ctype, text in _notebook_sources() if ctype == "code"
            for k, line in enumerate(text.splitlines(), 1)
            if rx.search(line) and not line.strip().startswith("#")]
    hits += [f"utils.py line {k}"
             for k, line in enumerate((NB_DIR / "utils.py").read_text(encoding="utf-8")
                                      .splitlines(), 1)
             if rx.search(line) and not line.strip().startswith("#")]
    assert not hits, "raw savefig (use viz.save_figure): " + ", ".join(hits)


@check("E14 NB03 never selects its targets from a within-cohort species-metabolite statistic")
def e14():
    """The anti-circularity guard for NB03's target list (v2.2).

    NB03 asks, per target, whether the SPECIES predict that metabolite, against a label
    permutation null. That question is only meaningful while the target list was chosen
    independently of how strongly species track the metabolite in this same cohort. Choose
    the strongest species-metabolite correlations as targets and the test rejects by
    construction, and the predictive evidence source in NB08/NB09 stops being independent
    of the statistical one -- it becomes a restatement of it.

    So NB03's two admissible criteria are:
      * robust DIFFERENTIAL ABUNDANCE (disease axis -- a different axis from species
        association), via differential_abundance_v2.csv; and
      * a PRE-SPECIFIED class, from an enzyme family in the literature
        (utils.POLYAMINE_KEGG, utils.DECARBOXYLATION_TARGET_KEGG).

    Every table banned below carries a within-cohort species x metabolite association, so
    reading any of them in NB03 would reintroduce the circularity. (NB08 and NB09 read
    them freely -- that is where evidence is combined, not where targets are chosen.)
    """
    banned = (
        "metabolome_partial_correlations",     # the metabolome-wide grid (NB02 7.2b)
        "polyamine_partial_correlations",      # the polyamine grid (NB02 7.2)
        "CRC_binary_stage_correlations",
        "CRC_progression_ordinal_correlations",
        "replication_genus_metabolome",        # rows conditioned on Yachida significance
        "pair_evidence_nb08",
        "CRC_progression_top_pairs",           # ranked species x metabolite pairs (CONTRACTS 5c)
    )
    doc = json.loads((NB_DIR / "03_ml_benchmarking.ipynb").read_text(encoding="utf-8"))
    hits = []
    for i, cell in enumerate(doc.get("cells", [])):
        if cell.get("cell_type") != "code":
            continue
        for lineno, line in enumerate("".join(cell["source"]).splitlines(), 1):
            if line.strip().startswith("#"):
                continue                       # prose may name what it must not read
            for name in banned:
                if name in line:
                    hits.append(f"cell {i} ({cell.get('id')}) line {lineno}: {name}")
    assert not hits, (
        "NB03 reads a within-cohort species-metabolite association table, which would "
        "make its permutation test circular:\n  " + "\n  ".join(hits))


@check("E16 NB03 takes every robust DA metabolite: the trinity target set is uncapped")
def e16():
    """CONTRACTS 5c (v2.3): the MIMOSA2/MICOM/MelonnPan trinity covers every robust DA
    metabolite plus the pre-specified class, and NB03 must cover the same set so each
    target carries the predictive source. Until v2.3, select_robust_balanced kept the
    8 most elevated + 7 most depleted per comparison, leaving 13 of 33 robust metabolites
    unmodelled. This fails if a cap comes back, as a default or at a call site."""
    doc = json.loads((NB_DIR / "03_ml_benchmarking.ipynb").read_text(encoding="utf-8"))
    src = "\n".join("".join(c["source"]) for c in doc.get("cells", [])
                    if c.get("cell_type") == "code")
    assert "def select_robust_balanced(da_df, n_up=None, n_down=None" in src, (
        "select_robust_balanced must default to no cap (n_up=None, n_down=None)")
    calls = re.findall(r"select_robust_balanced\(([^)]*)\)", src)
    calls = [c for c in calls if not c.startswith("da_df, n_up=None")]   # the def itself
    assert calls, "NB03 no longer calls select_robust_balanced"
    capped = [c for c in calls if "n_up" in c or "n_down" in c]
    assert not capped, f"NB03 caps its robust targets at a call site: {capped}"


@check("E15 NB09's synthesis reads NB08's assembler and re-derives none of its sources")
def e15():
    """The final tier must be assigned on ONE table, not a second weaker synthesis.

    Until v2.2 this section scored three sources over the seven polyamines only, and its
    genomic source used utils.species_encodes_enzyme -- which means curated producing ECs
    where POLYAMINE_EC has an entry and "any reaction indexed under this compound"
    everywhere else, two different claims under one name. NB08 separates those into
    genomic_producer and genomic_capability over the whole metabolome. Scoring both would
    leave two different answers to one question under near-identical filenames, which is
    the failure CONTRACTS section 2 records for the retired per-comparison DA family.

    So the synthesis cells must read pair_evidence_nb08_v2.csv and must not recompute a
    source NB08 owns. Module B may still call AGORA2 directly -- it profiles producer
    enzymes rather than scoring a pair -- so this is scoped to the synthesis cells.
    """
    SYNTH_CELLS = {"synthesis-build", "synthesis-score", "synthesis-fig"}
    REDERIVED = ("species_encodes_enzyme", "species_touches_compound",
                 "agora2_base_rate", "load_agora2_ec_index",
                 "genus_metabolite_replication")
    doc = json.loads((NB_DIR / "09_pathway_network.ipynb").read_text(encoding="utf-8"))
    cells = {c.get("id"): "".join(c["source"])
             for c in doc.get("cells", []) if c.get("cell_type") == "code"}

    missing = SYNTH_CELLS - set(cells)
    assert not missing, f"NB09 lost its synthesis cell(s): {sorted(missing)}"

    all_code = "\n".join(ln for ln in "\n".join(cells.values()).splitlines()
                         if not ln.strip().startswith("#"))
    assert "pair_evidence_nb08_v2.csv" in all_code, \
        "NB09 does not read pair_evidence_nb08_v2.csv -- the final tier would be scored " \
        "on a different, weaker candidate pool than the one NB08 assembled"

    hits = []
    for cid in sorted(SYNTH_CELLS):
        for lineno, line in enumerate(cells[cid].splitlines(), 1):
            if line.strip().startswith("#"):
                continue               # prose may name what it must not recompute
            for name in REDERIVED:
                if name in line:
                    hits.append(f"{cid} line {lineno}: {name}")
    assert not hits, (
        "NB09's synthesis re-derives a source NB08 already scored, so the two could "
        "disagree:\n  " + "\n  ".join(hits))

    # The one source NB09 legitimately adds, and the gate it must face like the others.
    score = cells["synthesis-score"] + cells["synthesis-build"]
    assert "mofa_coloading" in score, "Source F (MOFA co-loading) is not scored"
    assert "counts_toward_tier" in score, (
        "the MOFA source does not face the base-rate gate; a source firing for most "
        "pairs cannot discriminate and would raise every tier for free")
    assert "evidence_tools.tier" in score, \
        "the tier must come from evidence_tools.tier, which states the ceiling"


@check("E13 NB10 presolves every micom solve and keeps no local copy of micom_tools' rules")
def e13():
    # Without presolve the tradeoff sweep ran >22 h on one sample; local copies of the
    # flux/id rules are how unweighted fluxes and raw taxon ids reached the outputs.
    doc = json.loads((NB_DIR / "10_micom_flux.ipynb").read_text(encoding="utf-8"))
    src = "\n".join("".join(c["source"]) for c in doc["cells"] if c["cell_type"] == "code")
    code = "\n".join(ln for ln in src.splitlines() if not ln.strip().startswith("#"))
    calls = re.findall(r"\b(?:grow|micom_tradeoff)\((?!\))(?:[^()]|\([^()]*\))*\)", code, re.S)
    assert calls, "no grow()/tradeoff() call found in NB10"
    bad = [c.split("(")[0] for c in calls if "presolve=True" not in c]
    assert not bad, f"micom solve(s) without presolve=True: {bad}"
    local = re.findall(r"^\s*def (classify_metabolite|micom_taxon_id|taxon_name_map|choose_tradeoff|"
                       r"species_export|community_balance|polyamine_exchanges)\b", code, re.M)
    assert not local, f"NB10 redefines micom_tools helpers: {local}"
    assert "idxmax()" not in code, "tradeoff picked with idxmax (first tie) instead of mt.choose_tradeoff"


@check("E11 the theme renders subscript labels in DejaVu Sans with editable PDF fonts")
def e11():
    # Under seaborn's whitegrid alone, Arial resolved and "-log10" with subscript digits
    # rendered as boxes (with only a warning); PDFs embedded Type 3 fonts.
    import io
    import viz
    from matplotlib.figure import Figure
    with warnings.catch_warnings():               # apply_theme's error filter stays local
        viz.apply_theme()
        fig = Figure(figsize=(3, 2))
        ax = fig.subplots()
        ax.set_ylabel("−log₁₀(q)  ρ  δ  ≥  ±  →")
        ax.set_xlabel(r"$\log_{10}$ abundance")
        png, pdf = io.BytesIO(), io.BytesIO()
        fig.savefig(png, format="png")            # raises if a glyph is missing
        fig.savefig(pdf, format="pdf")
    data = pdf.getvalue()
    assert b"/Type3" not in data, "the PDF embeds Type 3 fonts (pdf.fonttype is not 42)"
    assert b"DejaVuSans" in data, "the PDF does not embed DejaVu Sans"


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
