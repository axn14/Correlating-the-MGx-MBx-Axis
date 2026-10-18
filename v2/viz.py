"""Centralised figure theme for the v2 pipeline.

WHY THIS MODULE EXISTS
----------------------
Before it, styling lived nowhere: `utils.py` had plotting helpers but no theme, and
each notebook re-declared its own. Measured consequences (audit 2026-10-07):

  * `figure.dpi` was 100 in NB01/01b/02/06, 110 in NB03/05, and unset in NB04/08/09,
    so identical `fontsize=` arguments rendered at different physical sizes.
  * `'Healthy'` was drawn in FIVE different colours across figures
    ('#2ecc71', '#a0a0a0', '#9ec5f4', '#dfe6e9', and Set2 slot 0).
  * '#2ecc71' simultaneously meant 'Healthy', 'Other metabolite' (network node) and
    'species' (SHAP legend). '#3498db' carried five unrelated meanings.
  * Correlation SIGN -- the main message of the network figures -- was encoded
    red/green in five places, the most common colour-vision confusion axis.
  * Severity was drawn with `viridis`, `magma_r` and `Blues_r` as CATEGORICAL fills
    in three adjacent NB08 figures, and with `RdYlBu_r` (a diverging map with a pale
    midpoint) as a SEQUENTIAL map, so mid-stage points vanished against the panel.

The fix is to assign colour by the JOB IT DOES, not per figure:

  | Variable                                   | Job               | Encoding        |
  |--------------------------------------------|-------------------|-----------------|
  | severity (Healthy->Adenoma->Early->Adv)    | ordinal magnitude | SEQUENTIAL      |
  | effect direction (MeanDiff, Cliff's d, rho)| polarity          | DIVERGING       |
  | cohort (Yachida/Sinha/Kim)                 | identity          | CATEGORICAL     |
  | network node type                          | identity          | CATEGORICAL (2) |

Two rules follow and are enforced by the lookups below rather than documented:
colour follows the ENTITY, never its rank (so a group keeps its colour when the set
of present groups changes -- see `group_color`), and categorical hues are assigned
in fixed order, never cycled.

!! PALETTE NOT YET VALIDATED !!
The hexes below are consolidated from values already in the repo, but the six-checks
validator has NOT been run against them (no Node on the authoring machine). Before
relying on these in a published figure, run:

    node <dataviz-skill>/scripts/validate_palette.js "#2a78d6,#eb6834,#2e9e5b" --mode light
    node <dataviz-skill>/scripts/validate_palette.js "#2a78d6,#eb6834,#2e9e5b" --mode dark
    node <dataviz-skill>/scripts/validate_palette.js "#7a4fd4,#1b9aaa,#d68910" --mode light

and re-step any pair that FAILs the CVD / normal-vision separation checks.
"""

from __future__ import annotations

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

# ============================================================================
# TOKENS
# ============================================================================

# Ink / surface tokens. Ported from ml_workbench/app.py:32-40, which was the only
# real token system in the repo but is plotly-only, so the notebooks got none of it.
INK      = "#0b0b0b"   # primary text
INK2     = "#52514e"   # secondary text
MUTED    = "#898781"   # tick labels, de-emphasised text
GRID     = "#e1e0d9"   # gridlines, axis spines
SURFACE  = "#fcfcfb"   # figure/axes background

# --- SEQUENTIAL: severity. One hue, light -> dark. ---------------------------
# Anchors are the existing `seq_blue` ramp (02_analysis.ipynb:2285) extended to four
# steps so the opt-in adenoma level has a slot. This is the ONLY severity encoding;
# do not reach for viridis/magma/Blues_r.
RAMP_SEVERITY = ["#cde2fb", "#6da7ec", "#2a78d6", "#104281"]

# Groups deliberately OFF the severity axis (PostSurgery / HS carry NaN in
# STAGE_GROUP_NUMERIC on purpose). They must not get a ramp step, or the figure
# implies an ordinal position the analysis explicitly refuses to give them.
# Pair with hatching -- see `OFF_AXIS_HATCH`.
COLOR_OFF_AXIS = MUTED
OFF_AXIS_HATCH = "//"

# --- DIVERGING: effect direction. Two hues + a NEUTRAL GREY midpoint. --------
# Replaces every red/green sign encoding. Never put a hue at the midpoint.
DIVERGING = ("#eb6834", GRID, "#2a78d6")   # depleted/negative | zero | enriched/positive

# --- CATEGORICAL: cohort identity. Fixed order, never cycled. ---------------
PALETTE_COHORT = ["#2a78d6", "#eb6834", "#2e9e5b"]
COHORT_COLOR = {
    "YACHIDA-CRC-2019":  PALETTE_COHORT[0],
    "SINHA_CRC_2016":    PALETTE_COHORT[1],
    "Kim_adenomas_2020": PALETTE_COHORT[2],
}

# --- CATEGORICAL: network node type. A SEPARATE slot set. -------------------
# Deliberately not PALETTE_COHORT: reusing those hexes is how '#2ecc71' ended up
# meaning 'Healthy', 'Other metabolite' and 'species' at once. The repo's old triple
# was red/blue/green, which fails red-green colour vision on the figure's own legend.
PALETTE_NODE = ["#7a4fd4", "#1b9aaa", "#d68910"]
NODE_COLOR = {
    "polyamine":  PALETTE_NODE[0],
    "species":    PALETTE_NODE[1],
    "metabolite": PALETTE_NODE[2],
}

# Canonical severity orderings this module knows how to colour. Longest first, so a
# lookup matches the most specific scale a label belongs to.
_SEVERITY_SCALES = [
    ["Healthy", "Adenoma", "Early_CRC", "Advanced_CRC"],   # SEVERITY_TRAIN_ORDER
    ["Healthy", "Stage_0", "Stage_I_II", "Stage_III_IV"],  # CRC_STAGE_ORDER
    ["Healthy", "Adenoma", "Carcinoma"],                   # SEVERITY_EVAL_ORDER
    ["Healthy", "Early_CRC", "Advanced_CRC"],              # STAGE_GROUP_ORDER
    ["Healthy", "CRC"],                                    # validation cohorts
]

# Labels that mean "off the severity axis", in any of the spellings used upstream.
_OFF_AXIS_LABELS = {"PostSurgery", "HS"}

# Aliases, so the same biological group gets the same colour whichever vocabulary a
# notebook happens to be using.
_SEVERITY_ALIASES = {
    "Control": "Healthy",
    "nonIBD": "Healthy",
    "MP": "Adenoma",
    "Carcinoma": "CRC",
}


# ============================================================================
# THEME
# ============================================================================

def apply_theme(context: str = "notebook", dpi: int = 110) -> None:
    """Install the shared theme. Call ONCE per notebook, right after the imports.

    Replaces the per-notebook `sns.set_style(...)` + `plt.rcParams['figure.dpi']`
    blocks, which disagreed on dpi (100/110/unset) and on whether to call
    `set_context` at all -- the reason identical `fontsize=` values rendered at
    different physical sizes in different notebooks.
    """
    sns.set_style("whitegrid")
    sns.set_context(context)
    mpl.rcParams.update({
        "figure.dpi":        dpi,
        "savefig.dpi":       150,
        "figure.facecolor":  SURFACE,
        "axes.facecolor":    SURFACE,
        "savefig.facecolor": SURFACE,
        # Text wears text tokens, never a series colour.
        "text.color":        INK,
        "axes.labelcolor":   INK2,
        "axes.titlecolor":   INK,
        "xtick.color":       MUTED,
        "ytick.color":       MUTED,
        # Recessive grid and axes; thin marks.
        "axes.edgecolor":    GRID,
        "grid.color":        GRID,
        "grid.linewidth":    0.6,
        "axes.linewidth":    0.8,
        "lines.linewidth":   2.0,
        "lines.markersize":  8.0,
        "axes.grid":         True,
        "axes.axisbelow":    True,
        # Never let matplotlib invent a hue: the default cycle is the fixed
        # categorical order, so an unspecified series cannot collide with a
        # semantic colour by accident.
        "axes.prop_cycle":   mpl.cycler(color=PALETTE_COHORT + PALETTE_NODE),
        "legend.frameon":    False,
        "figure.titlesize":  13,
    })


# ============================================================================
# COLOUR LOOKUPS
# ============================================================================

def severity_cmap(name: str = "severity") -> LinearSegmentedColormap:
    """Continuous sequential colormap for severity / any ordinal magnitude."""
    return LinearSegmentedColormap.from_list(name, RAMP_SEVERITY)


def severity_colors(n: int | None = None, order=None) -> list:
    """`n` evenly spaced steps of the severity ramp (or one per label in `order`).

    Use for an ORDERED axis only. A label in `order` that is off the severity axis
    (PostSurgery / HS) gets COLOR_OFF_AXIS instead of a ramp step, because giving it
    one would assert an ordinal position that STAGE_GROUP_NUMERIC deliberately
    leaves as NaN.
    """
    if order is not None:
        order = list(order)
        on_axis = [g for g in order if g not in _OFF_AXIS_LABELS]
        steps = severity_colors(n=len(on_axis))
        it = iter(steps)
        return [COLOR_OFF_AXIS if g in _OFF_AXIS_LABELS else next(it) for g in order]

    if n is None:
        n = len(RAMP_SEVERITY)
    if n <= 1:
        return [RAMP_SEVERITY[-2]]
    if n <= len(RAMP_SEVERITY):
        # Prefer the hand-picked anchors over interpolation when they suffice.
        idx = np.linspace(0, len(RAMP_SEVERITY) - 1, n).round().astype(int)
        return [RAMP_SEVERITY[i] for i in idx]
    cmap = severity_cmap()
    return [mpl.colors.to_hex(cmap(x)) for x in np.linspace(0.0, 1.0, n)]


def group_color(group: str, order=None) -> str:
    """Colour for a Study.Group / Stage.Group label, BY IDENTITY not by rank.

    This is the fix for `utils.plot_pca`, which assigned colours from
    `sorted(unique())` and therefore gave the same group a different colour in a
    different cohort whenever the set of present groups differed -- so NB01's three
    `plot_pca` loops were not colour-comparable with each other.

    Passing `order` pins the ramp to a specific scale; otherwise the longest
    canonical scale containing `group` is used, so the answer does not depend on
    which other groups happen to be present.
    """
    g = _SEVERITY_ALIASES.get(str(group), str(group))
    if g in _OFF_AXIS_LABELS or str(group) in _OFF_AXIS_LABELS:
        return COLOR_OFF_AXIS
    scales = [list(order)] if order is not None else _SEVERITY_SCALES
    for scale in scales:
        canon = [_SEVERITY_ALIASES.get(s, s) for s in scale]
        if g in canon:
            return severity_colors(n=len(canon))[canon.index(g)]
    return MUTED


def group_palette(groups, order=None) -> dict:
    """{group: colour} for seaborn's `palette=`, built by identity."""
    return {g: group_color(g, order=order) for g in groups}


def cohort_color(dataset: str) -> str:
    """Fixed colour per cohort. Identity, so it never shifts with cohort count."""
    return COHORT_COLOR.get(str(dataset), MUTED)


def node_color(kind: str) -> str:
    """Colour for a network node type: 'polyamine' | 'species' | 'metabolite'."""
    return NODE_COLOR.get(str(kind).lower(), MUTED)


def diverging_cmap(name: str = "effect") -> LinearSegmentedColormap:
    """Diverging colormap for signed quantities. Neutral grey midpoint, no hue."""
    return LinearSegmentedColormap.from_list(name, list(DIVERGING))


def diverging_norm(vmin: float, vmax: float, center: float = 0.0) -> TwoSlopeNorm:
    """TwoSlopeNorm anchored at `center`, so zero sits on the neutral midpoint.

    Guards the degenerate cases TwoSlopeNorm rejects (vmin == center, or data
    entirely on one side), which otherwise raise mid-figure.
    """
    vmin, vmax = float(min(vmin, vmax)), float(max(vmin, vmax))
    span = max(abs(vmin - center), abs(vmax - center)) or 1.0
    return TwoSlopeNorm(vmin=center - span, vcenter=center, vmax=center + span)


def sign_color(value: float) -> str:
    """Colour for the SIGN of an effect. Replaces red/green sign encoding."""
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return MUTED
    return DIVERGING[2] if value > 0 else DIVERGING[0] if value < 0 else GRID


# ============================================================================
# LAYOUT
# ============================================================================

def panel_grid(n: int, ncols: int | None = None, panel_size=(4.0, 3.5), **kwargs):
    """`(fig, axes)` for `n` panels, with unused panels hidden.

    Supersedes a pattern copy-pasted eight times across the notebooks in three
    mutually inconsistent forms: `squeeze=False` + `.flatten()`, a bare `axes[0]`,
    and a manual `if len(x) == 1: axes = [axes]` (09_pathway_network.ipynb:1042).
    Only NB05 stage 6 uses this helper so far; the other call sites are NOT yet
    migrated -- see the plan's W4.

    `axes` is always a flat list of length `n`, so callers never branch on the
    single-panel case.
    """
    if n < 1:
        raise ValueError(f"panel_grid needs at least one panel, got {n}")
    if ncols is None:
        ncols = min(4, n)
    ncols = max(1, min(ncols, n))
    nrows = int(np.ceil(n / ncols))
    kwargs.setdefault("figsize", (panel_size[0] * ncols, panel_size[1] * nrows))
    kwargs.setdefault("constrained_layout", True)
    fig, axes = plt.subplots(nrows, ncols, squeeze=False, **kwargs)
    flat = list(axes.flatten())
    for ax in flat[n:]:
        ax.set_visible(False)
    return fig, flat[:n]


def despine(ax, keep=("left", "bottom")) -> None:
    """Hide top/right spines and recede the rest.

    Supersedes a verbatim 3-line spine/tick loop that is still copy-pasted five
    times in NB02 and `edits_nb02c.py` (`:270-275`, `:298-303`). Those call sites
    have NOT been migrated yet -- see the plan's W4.
    """
    for name, spine in ax.spines.items():
        if name in keep:
            spine.set_color(GRID)
        else:
            spine.set_visible(False)
    ax.tick_params(colors=MUTED)


def save_figure(fig, path, also_pdf: bool = True, close: bool = True) -> None:
    """Save a figure consistently, and in vector form alongside the raster.

    Supersedes `dpi=150, bbox_inches='tight'` typed out ~45 times, only two of which
    also passed `facecolor=fig.get_facecolor()` -- so every other figure lost the
    off-white canvas. `also_pdf` exists because the pipeline emitted nothing but
    150 dpi PNG, leaving no figure re-typesettable for a journal.

    Only a couple of call sites use this so far (`edits_nb02a.py`, NB05 stage 6);
    the rest still call `fig.savefig(...)` directly -- see the plan's W4.

    Closes the figure by default: the notebooks previously `savefig`'d AND
    `plt.show()`ed every figure, and `nbconvert --execute --inplace` then wrote a
    second base64 copy of each one back into the .ipynb on disk.
    """
    from pathlib import Path

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor())
    if also_pdf:
        fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight",
                    facecolor=fig.get_facecolor())
    if close:
        plt.close(fig)
