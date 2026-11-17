"""Phase E3, class M (CONTRACTS 5d): is a taxon FORCED to export a target?

A pFBA solution is one point of a degenerate optimal face: a taxon that exports a compound
there may export none at another, equally optimal point. "The reaction exists and pFBA
used it" is reference potential (class R), not mechanism. Class M asks the stronger
question: with every taxon's growth fixed at the certified cooperative-tradeoff optimum,
is the taxon's net export of the target bounded AWAY FROM ZERO (FVA minimum > tolerance)?

The feasible set is reproduced exactly as micom builds it (micom 0.39.1,
problems.cooperative_tradeoff + solution.add_pfba_objective): community growth maximised;
the L2 tradeoff QP solved at fraction alpha (through `highs_qp`, certified); every taxon's
growth then fixed with `_apply_min_growth(rates, atol, rtol)`. Inside that set:

1. pFBA gives the reference fluxes; candidates are taxon->lumen exchanges of a target with
   export > tol there (an export that is zero at a feasible point cannot be forced).
2. Elimination: minimise the SUM of the remaining candidates' exports. Any candidate at
   <= tol in that solution has a feasible point with no export: not forced, dropped.
   Repeat while something drops.
3. Each survivor is minimised on its own; forced iff its minimum > tol.

Every LP goes to HiGHS (the hybrid interface's LP path). Export orientation is read per
reaction from which side the lumen metabolite is on.
"""
from __future__ import annotations

import time

import pandas as pd
from optlang.interface import OPTIMAL
from optlang.symbolics import Zero

EXPORT_TOL = 1e-6


def _lumen_kegg(met) -> set[str]:
    raw = str(getattr(met, "annotation", {}).get("kegg.compound", "NaN"))
    return {x.strip() for x in raw.split(";")} if raw != "NaN" else set()


def target_exchanges(com, kegg_targets: set[str]):
    """[(reaction, taxon, kegg, sign)] for taxon<->lumen exchanges of target compounds.

    sign is +1 when positive flux moves the compound INTO the lumen (export), else -1.
    """
    out = []
    for r in com.internal_exchanges:
        lumen = [m for m in r.metabolites if m.compartment == "m"]
        if not lumen:
            continue
        ks = _lumen_kegg(lumen[0]) & kegg_targets
        if not ks:
            continue
        sign = 1.0 if r.metabolites[lumen[0]] > 0 else -1.0
        for k in ks:
            out.append((r, r.community_id, k, sign))
    return out


def _export_value(r, sign: float) -> float:
    return sign * (r.forward_variable.primal - r.reverse_variable.primal)


def _minimise_exports(com, items) -> bool:
    com.objective = Zero
    com.objective_direction = "min"
    coefs = {}
    for r, _, _, s in items:
        coefs[r.forward_variable] = coefs.get(r.forward_variable, 0.0) + s
        coefs[r.reverse_variable] = coefs.get(r.reverse_variable, 0.0) - s
    com.objective.set_linear_coefficients(coefs)
    com.solver.optimize()
    return com.solver.status == OPTIMAL


def forced_exports(com, alpha: float, kegg_targets: set[str],
                   atol: float | None = None, rtol: float | None = None) -> pd.DataFrame:
    """Per (taxon, target): reference export, and whether the export is forced."""
    from micom.problems import regularize_l2_norm
    from micom.solution import optimize_with_retry
    from micom.util import _apply_min_growth, _format_min_growth

    tol = com.solver.configuration.tolerances.feasibility
    atol = tol if atol is None else atol
    rtol = tol if rtol is None else rtol
    t0 = time.time()
    rows = []
    with com:
        _apply_min_growth(com, _format_min_growth(0.0, com.taxa))
        com.objective = com.scale * com.variables.community_objective
        max_growth = optimize_with_retry(com, "could not get community growth rate.") / com.scale
        regularize_l2_norm(com, 0.0)
        com.variables.community_objective.lb = alpha * max_growth
        com.variables.community_objective.ub = max_growth
        com.solver.optimize()                       # the tradeoff QP (highs_qp)
        qp_status = com.solver.status
        rates = {sp: com.constraints["objective_" + sp].primal for sp in com.taxa}
        _apply_min_growth(com, rates, atol, rtol)

        # 1. reference pFBA point
        com.objective = Zero
        com.objective_direction = "min"
        com.objective.set_linear_coefficients(
            {v: 1.0 for rxn in com.reactions for v in (rxn.forward_variable, rxn.reverse_variable)})
        com.solver.optimize()
        if com.solver.status != OPTIMAL:
            return pd.DataFrame([{"status": f"pfba {com.solver.status}", "qp_status": qp_status}])
        items = target_exchanges(com, kegg_targets)
        ref = {(r.id, k): _export_value(r, s) for r, _, k, s in items}
        cand = [it for it in items if ref[(it[0].id, it[2])] > EXPORT_TOL]
        n_lp, n_candidates = 0, len(cand)

        # 2. elimination by summed minimisation
        alive = list(cand)
        while alive:
            if not _minimise_exports(com, alive):
                break
            n_lp += 1
            keep = [it for it in alive if _export_value(it[0], it[3]) > EXPORT_TOL]
            if len(keep) == len(alive):
                break
            alive = keep
        # 3. individual minima for the survivors
        minima = {}
        for it in alive:
            if _minimise_exports(com, [it]):
                minima[(it[0].id, it[2])] = _export_value(it[0], it[3])
            else:
                minima[(it[0].id, it[2])] = float("nan")
            n_lp += 1

        for r, taxon, k, s in items:
            key = (r.id, k)
            e = ref[key]
            m = minima.get(key, 0.0 if e > EXPORT_TOL else min(e, 0.0))
            rows.append({"taxon": taxon, "kegg": k, "reaction": r.id,
                         "export_pfba": e, "export_min": m,
                         "forced": bool(m > EXPORT_TOL)})
    df = pd.DataFrame(rows)
    # Growth rates at the certified optimum, so a batch needs one solve per sample: the
    # pFBA step above is an LP solved to OPTIMAL inside the fixed-growth set, so the
    # point is mass-balanced by construction (the failure micom_tools.solution_quality
    # exists to catch cannot occur here).
    df.attrs.update({"qp_status": str(qp_status), "n_candidates": n_candidates,
                     "n_lp": n_lp, "seconds": time.time() - t0, "max_growth": max_growth,
                     "growth_rates": rates,
                     "community_growth": float(sum(
                         rates[sp] * com.taxonomy.loc[sp, "abundance"] for sp in com.taxa))})
    return df
