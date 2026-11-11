"""Shared rules for the metabolome-wide evidence assembler (NB08's synthesis section).

Why this module exists rather than more functions in `utils.py`: `run_pipeline.py` hashes
`utils.py` as a shared input (`SHARED_CODE`), so editing it invalidates every cached
notebook result -- including NB03, which costs ~4.5 h to recompute. These helpers are
needed only by the evidence assembler, so they live here, the way the MICOM rules live in
`micom_tools.py`. `checks.py` group D tests them.

Everything here encodes one of three rules the assembler must not restate inline:

1. **Genus derivation** (`species_to_genus`) -- the replication table is genus-level, so a
   species-level pair can only be matched to it through a normalised genus. One
   implementation, shared with `utils.normalise_genus_name`, so a species-level and a
   genus-level table cannot drift apart.

2. **Tri-state AGORA2 calls with a base-rate gate** -- "missing coverage is never negative
   evidence" (CONTRACTS section 7), and a source that nearly every species satisfies
   cannot single anything out. Both rules are mechanical, and both were got wrong before.

3. **Confounder attributability** (`confounder_unattributable`) -- NB03 reports the AUC the
   confounders alone reach on the *adjusted* split, which is ~0.5 when residualisation
   worked. Where it instead meets or exceeds the species model's own AUC, that target's
   SHAP ranking is not evidence about species, and the assembler must withhold predictive
   support rather than report a tier built on it.
"""

from __future__ import annotations

import math

import numpy as np
import pandas as pd

import utils

# A direction-agnostic capability call is weaker evidence than a producer-specific one, so
# the two are kept as separate columns and never merged. These labels name them in output.
GENOMIC_PRODUCER = "producer_ec"      # POLYAMINE_EC: reactions that PRODUCE the compound
GENOMIC_CAPABILITY = "capability"     # any reaction indexed under the compound, either way


def species_to_genus(species: str) -> str:
    """Genus of a species-level name ("Veillonella parvula_A" -> "Veillonella").

    Takes the first whitespace-delimited token and normalises it with
    `utils.normalise_genus_name`, which strips rank prefixes and GTDB placeholder
    suffixes. Using the shared normaliser is the point: the replication table's genera
    come from `utils.genus_clr_table`, which normalises the same way, so a species-level
    pair and a genus-level row that refer to the same organism actually join. Deriving the
    genus by a local `split()[0]` would match "Blautia" but silently miss "Blautia_A".
    """
    token = str(species).strip().split()[0] if str(species).strip() else ""
    return utils.normalise_genus_name(token)


def species_touches_compound(species: str, kegg_id: str, agora2_ec_index: dict):
    """Does `species`' AGORA2 reconstruction carry ANY reaction indexed under `kegg_id`?

    Returns True/False where AGORA2 has a model for the species, else None -- the same
    tri-state contract as `utils.species_encodes_enzyme`, for the same reason: no model
    means no evidence either way, which is not the same as evidence of absence.

    This is **metabolic capability, and it is direction-agnostic**: the reaction may
    consume the compound rather than produce it. It is therefore strictly weaker than
    `utils.species_encodes_enzyme`, which checks the curated producing ECs in
    `utils.POLYAMINE_EC`. The reason it exists at all is coverage: the producing-EC map
    covers 7 compounds, so metabolome-wide it would leave the genomic source empty for
    every other metabolite, whereas the index is keyed by KEGG compound for every compound
    the reconstructions mention. Reconstructions also include gap-filled reactions, so
    even a True is capability, never observed production.
    """
    species_id, _ = utils.agora2_species_id(species, agora2_ec_index)
    if species_id is None:
        return None
    return bool(agora2_ec_index[species_id].get(kegg_id))


def compound_base_rate(kegg_id, species, agora2_ec_index, call=species_touches_compound):
    """(base_rate, n_assessed) for a compound: the fraction of species WITH an AGORA2
    model whose call is True.

    Species without a model return None and are excluded from the denominator rather than
    counted as negatives, which would deflate the rate and make an uninformative source
    look selective.

    `base_rate` is NaN when no species in `species` has a model; a NaN rate can never pass
    `counts_toward_tier`, so an unassessable compound contributes nothing by default.
    """
    calls = [call(s, kegg_id, agora2_ec_index) for s in species]
    assessed = [bool(c) for c in calls if c is not None]
    if not assessed:
        return float("nan"), 0
    return sum(assessed) / len(assessed), len(assessed)


def counts_toward_tier(base_rate, max_base_rate=None) -> bool:
    """May a source with this base rate count toward an evidence tier?

    Pre-declared rule (council D7, CONTRACTS section 8): only where the base rate is
    <= `utils.EVIDENCE_MAX_BASE_RATE` (0.2). Among the polyamines, six of seven modelled
    base rates sit between 0.58 and 0.97 -- nearly every modelled species "encodes the
    enzyme", so the call cannot distinguish a producer from anything else. The raw call is
    still reported; it just does not earn a point.

    A NaN base rate (nothing assessable) is False, never True.
    """
    if max_base_rate is None:
        max_base_rate = utils.EVIDENCE_MAX_BASE_RATE
    try:
        rate = float(base_rate)
    except (TypeError, ValueError):
        return False
    return (not math.isnan(rate)) and rate <= max_base_rate


def confounder_unattributable(auc_model, auc_confounders_adjusted) -> bool:
    """Is this target's species model NOT attributable to the species?

    True when the confounders ALONE reach an AUC on the adjusted split that meets or
    exceeds the species model's own nested AUC. NB03 residualises each target on the
    confounders and reports `AUC_confounders_adjusted`, which is ~0.5 when that worked
    (measured 0.43-0.60 for 25 of 27 targets). Where it instead matches or beats the
    model, OLS residualisation did not remove a nonlinear confounder effect and the
    model's advantage over chance cannot be credited to the microbiome.

    Measured on the 2026-10-26 full run, two targets fail this test -- spermine
    (confounders 0.823 vs model 0.633) and N1-acetylspermine (0.741 vs 0.577) -- and both
    would otherwise have carried a significant permutation q into a tier.

    The comparison is deliberately threshold-free: it asks whether the species add
    anything over the confounders, not whether the confounders clear some fixed bar.
    Either value missing returns False (not assessable, so not a disqualification).
    """
    try:
        model = float(auc_model)
        conf = float(auc_confounders_adjusted)
    except (TypeError, ValueError):
        return False
    if math.isnan(model) or math.isnan(conf):
        return False
    return conf >= model


def predictive_support(row, fdr_threshold=None) -> bool:
    """Does NB03's SHAP ranking support this pair as a candidate producer?

    The definition NB09 has used for Source A, unchanged, plus the attributability gate:

      positive direction  -- more of the species pushes the prediction towards a high
                             metabolite level (a negative-direction species is not a
                             candidate producer)
      stable bootstrap rank
      the model beat its permutation null (q < FDR_THRESHOLD)
      the target is not confounder-unattributable

    A model no better than chance cannot lend its ranking as evidence, and neither can one
    whose advantage belongs to the confounders.
    """
    if fdr_threshold is None:
        fdr_threshold = utils.FDR_THRESHOLD
    if row.get("Confounder_unattributable", False):
        return False
    try:
        q = float(row.get("Model_perm_q", float("nan")))
    except (TypeError, ValueError):
        return False
    if math.isnan(q) or q >= fdr_threshold:
        return False
    return (str(row.get("Direction", "")) == "positive"
            and bool(row.get("Rank_stable", False)))


def replication_summary(repl_row, cohorts) -> dict:
    """Collapse one genus-level replication row into the assembler's columns.

    Keeps the three outcomes the replication table separates and `checks.py` C9 protects:
    measurable-and-replicated, measurable-and-not, and NOT MEASURABLE. A pair the
    validation cohort cannot measure (compound not identity-matched, or genus absent) is
    reported as unavailable, never as a failed replication -- and a pair significant with
    the OPPOSITE sign is a contradiction, which is evidence against, not absence of
    evidence.
    """
    out = {
        "Replication_available": False,
        "Replication_support": False,
        "Replication_contradicted": False,
        "Replication_n_measurable": 0,
        "Replication_n_replicated": 0,
        "Replication_cohorts": "",
        "Replication_mean_rho": float("nan"),
    }
    if repl_row is None:
        return out

    measurable, replicated, contradicted, rhos = [], [], [], []
    for ds in cohorts:
        if not bool(repl_row.get(f"In_{ds}", False)):
            continue
        measurable.append(ds)
        rho = repl_row.get(f"Rho_{ds}", float("nan"))
        if pd.notna(rho):
            rhos.append(float(rho))
        if bool(repl_row.get(f"Replicated_{ds}", False)):
            replicated.append(ds)
        if bool(repl_row.get(f"Contradicted_{ds}", False)):
            contradicted.append(ds)

    out["Replication_available"] = bool(measurable)
    out["Replication_support"] = bool(replicated)
    out["Replication_contradicted"] = bool(contradicted)
    out["Replication_n_measurable"] = len(measurable)
    out["Replication_n_replicated"] = len(replicated)
    out["Replication_cohorts"] = ",".join(replicated)
    if rhos:
        out["Replication_mean_rho"] = float(np.mean(rhos))
    return out


def count_sources(row, sources) -> tuple:
    """(n_available, n_supporting, source_names) over `sources`.

    `sources` is a list of (label, available_col, support_col). A source counts as
    available only where the row says so, and supporting only where available AND
    supported -- so a source that could not be assessed can never be read as a negative,
    and a stray True on a support column cannot smuggle in a point.

    The statistical association is NOT a source here. It is how the pair was nominated, so
    counting it would hand every row in the table a free point.
    """
    available, supporting = [], []
    for label, avail_col, support_col in sources:
        if not bool(row.get(avail_col, False)):
            continue
        available.append(label)
        if bool(row.get(support_col, False)):
            supporting.append(label)
    return len(available), len(supporting), ",".join(supporting)


# ---------------------------------------------------------------------------------------
# Evidence classes (CONTRACTS section 5d, v2.4, 2026-11-11)
# ---------------------------------------------------------------------------------------
#
# Sources are not independent: NB03 prediction and MOFA co-loading are two readings of the
# same Yachida association, and AGORA2 capability, MIMOSA2 and a MICOM exchange's mere
# presence all read the same reference reconstructions. So a source belongs to a class, and
# a class counts once. Each source is named by the `<name>` in its NB08/NB09 columns
# `Source_<name>_available` / `Source_<name>_support`.

EVIDENCE_CLASSES = {
    "O": ("predictive", "mofa_coloading", "melonnpan"),
    "R": ("genomic_producer", "genomic_capability", "mimosa2", "micom_presence"),
    "G": ("gene_level",),
    "M": ("mechanistic",),
    "X": ("replication_species",),
}
CLASS_NAMES = {"O": "observational", "R": "reference potential", "G": "sample-level genetic",
               "M": "contextual mechanism", "X": "external species-level replication"}
# Sources that are reported but can never vote for a species claim (section 5d-i).
NON_VOTING_SOURCES = ("replication",)       # genus-level Sinha/Kim: ecological consistency
DIRECT_CLASSES = ("G", "M")                  # the top tier needs one of these


def source_class(source: str) -> str | None:
    """The class a source belongs to, or None for a non-voting source.

    Raises on a source no class knows, so a new source cannot slip into a tier unclassified.
    """
    if source in NON_VOTING_SOURCES:
        return None
    for cls, members in EVIDENCE_CLASSES.items():
        if source in members:
            return cls
    raise KeyError(f"source {source!r} has no evidence class (CONTRACTS 5d)")


def class_support(row) -> tuple[list, list]:
    """(classes available, classes supporting) for one pair, in O, R, G, M, X order.

    A class is available when any member's `Source_<m>_available` is True, and supports
    when any AVAILABLE member's `Source_<m>_support` is True -- so, as in `count_sources`,
    a stray True on an unavailable source cannot supply a point. Members whose columns are
    absent from the row are simply not available (a tool that has not run is not a vote).
    """
    available, supporting = [], []
    for cls, members in EVIDENCE_CLASSES.items():
        avail = sup = False
        for m in members:
            if bool(row.get(f"Source_{m}_available", False)):
                avail = True
                if bool(row.get(f"Source_{m}_support", False)):
                    sup = True
        if avail:
            available.append(cls)
        if sup:
            supporting.append(cls)
    return available, supporting


def class_tier(row) -> str:
    """The section 5d tier for one pair, with the number of classes available stated.

    Conclusive          >= 3 classes, at least one of G or M, no contradiction, target
                        confounder-attributable and export-gated
    Supported           2 classes, at least one not O
    Single-class        1 class, not O
    Observational only  only O
    Unsupported         none

    Gates read from the row: `Replication_species_contradicted` (X significant with the
    opposite sign), `Confounder_unattributable`, and `Target_export_gate` (False only when
    the registry says the target has no modelled export route; missing counts as passed,
    because the gate removes classes R/M upstream rather than vetoing here).
    """
    available, supporting = class_support(row)
    n_av = len(available)
    sup = set(supporting)
    contradicted = bool(row.get("Replication_species_contradicted", False))
    unattributable = bool(row.get("Confounder_unattributable", False))
    gated = row.get("Target_export_gate", True)
    gated = True if gated is None or (isinstance(gated, float) and math.isnan(gated)) else bool(gated)
    tag = f"({len(sup)}/{n_av} classes: {''.join(supporting) or '-'})"
    if (len(sup) >= 3 and sup & set(DIRECT_CLASSES) and not contradicted
            and not unattributable and gated):
        return f"Conclusive {tag}"
    if len(sup) >= 2 and sup - {"O"}:
        return f"Supported {tag}"
    if len(sup) == 1 and sup != {"O"}:
        return f"Single-class {tag}"
    if sup == {"O"}:
        return f"Observational only {tag}"
    return f"Unsupported {tag}"


def tier(n_supporting, n_available) -> str:
    """Pre-declared tiers (CONTRACTS section 6/8), with the attainable ceiling stated.

    `Strong` >= 3 supporting, `Moderate` == 2, `Weak` == 1, else `Statistical only`. The
    label carries `n_available` because a tier is meaningless without the ceiling: a pair
    with 1 of 1 available source is not the same claim as 1 of 4, and the earlier
    polyamine synthesis reported "Weak (1/3)" for rows whose attainable maximum was 2.
    """
    n_supporting = int(n_supporting)
    n_available = int(n_available)
    if n_supporting >= 3:
        label = "Strong"
    elif n_supporting == 2:
        label = "Moderate"
    elif n_supporting == 1:
        label = "Weak"
    else:
        return f"Statistical only (0/{n_available} available supporting)"
    return f"{label} ({n_supporting}/{n_available} available supporting)"
