"""Solve micom's cooperative-tradeoff QP exactly with HiGHS LPs -- no OSQP, no Gurobi.

Why this exists (v2/micom/README.md, findings 4-5; CONTRACTS.md section 8):
optlang's "hybrid" interface, which micom uses when no commercial solver is installed,
sends every LP to HiGHS and every QP to OSQP. OSQP is a first-order (ADMM) method, so its
answer is approximate. micom's cooperative tradeoff is a QP whose growth rates are then
fixed as constraints for the pFBA step; on badly scaled communities the approximate point
leaves pFBA infeasible and micom returns a zero-growth, mass-unbalanced solution without
raising (sample 10091: lumen imbalance 18.3 against a 1e-6 tolerance). A Gurobi licence
was not obtainable (2026-11-11), and the project is to carry no commercial dependency.

What it does. The tradeoff QP is  min 0.5 * sum_j P_jj z_j^2  over the polytope F of
feasible community flux states, where z are the k taxon growth variables (k = number of
taxa; P is diagonal). That is the minimum-norm point of the projection of F onto those k
coordinates, so Wolfe's minimum-norm-point algorithm (Wolfe 1976) solves it using only a
linear-minimisation oracle -- one HiGHS LP per major iteration, with the same constraints
every time and only the objective changed, so each LP warm-starts from the last basis.
The returned flux vector is a convex combination of LP-optimal vertices of F, hence
feasible to HiGHS's LP tolerance -- exactly what the following pFBA step needs. The
optimum in z is unique (strictly convex in z), so the result does not depend on the
method; only its accuracy does, and the stopping rule bounds the optimality gap.

Tried first and rejected (2026-11-11, sample 10723, 74,813 columns x 33,155 rows):
- HiGHS 1.15.1's own QP solver (active set): "Unbounded" on a QP bounded below by 0
  after 124 s with presolve on; no answer within 30 min with presolve off.
- Kelley cutting planes on HiGHS LPs: gap fell 0.134 -> 0.0075 in 4 rounds, then the
  5th LP did not finish in 14 min (near-parallel cuts make the LP degenerate).

Use:  import highs_qp; highs_qp.install()   # before micom solves anything
`install()` patches optlang's HybridProblem.solve in the current process only. micom's
`grow(threads>1)` starts worker processes, so call install() at module level of the
script that calls grow (spawned workers re-import that module), or use threads=1.

Conventions checked against optlang 1.9.1 (matrix_interface.MatrixProblem.build): `P`
holds the quadratic coefficients times 2*direction in both triangles and `q` the linear
ones times direction, i.e. the problem is  min 0.5 x'Px + q'x.
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np
import scipy.sparse as sps

import highspy as hs
from optlang import hybrid_interface as hi

_ORIGINAL_SOLVE = hi.HybridProblem.solve

DEBUG = bool(os.environ.get("HIGHS_QP_DEBUG"))
GAP_TOL = 1e-10       # stop when <w,w> - min_p <w,p>  <=  GAP_TOL * max(<w,w>, 1e-12)
MAX_MAJOR = 500
LAST_RUN: dict = {}   # diagnostics from the most recent solve


def _affine_min(B: np.ndarray) -> np.ndarray:
    """Weights a (sum 1) minimising ||B a||^2 over the affine hull of B's columns."""
    m = B.shape[1]
    K = np.zeros((m + 1, m + 1))
    K[:m, :m] = B.T @ B
    K[:m, m] = 1.0
    K[m, :m] = 1.0
    rhs = np.zeros(m + 1)
    rhs[m] = 1.0
    return np.linalg.lstsq(K, rhs, rcond=None)[0][:m]


def min_norm_point(sp, options, x_start=None, cert_tol=None, max_major=MAX_MAJOR):
    """Wolfe's minimum-norm point for  min 0.5 x'Px  (P diagonal PSD, q = 0) over sp's LP.

    x_start: an exactly feasible flux vector to start from (e.g. a repaired OSQP point);
    without it the first corral point is an LP vertex. cert_tol: stop once the CERTIFIED
    relative gap (f - best lower bound) / f is at most this; without it, run to GAP_TOL.
    Every oracle LP at the current w yields the Frank-Wolfe lower bound
    f* >= f(w) - <w, w - p>, and Wolfe's iterates only decrease f, so the best bound seen
    certifies the returned point.
    Returns dict(status, x, major, minor, gap, f, lower_bound, cert_gap_rel,
                 n_square_terms, lp_seconds, n_lp, seconds).
    """
    P = sps.csc_matrix(sp.P)
    diag = P.diagonal()
    if (P - sps.diags(diag)).count_nonzero() > 0:
        raise hi.SolverError("min-norm-point QP needs a diagonal Hessian")
    if np.any(diag < 0):
        raise hi.SolverError("min-norm-point QP needs a convex (PSD) Hessian")
    if np.any(sp.q != 0):
        raise hi.SolverError("min-norm-point QP supports a purely quadratic objective only")
    J = np.flatnonzero(diag > 0)
    root = np.sqrt(diag[J])          # w = root * x_J turns the objective into 0.5 ||w||^2
    nv, nc, k = sp.vbounds.shape[0], sp.bounds.shape[0], len(J)

    env = hs.Highs()
    env.passOptions(options)
    model = hs.HighsModel()
    lp = model.lp_
    lp.sense_ = hs.ObjSense.kMinimize
    lp.num_col_ = nv
    lp.num_row_ = nc
    lp.col_cost_ = np.zeros(nv)
    lp.col_lower_ = sp.vbounds[:, 0]
    lp.col_upper_ = sp.vbounds[:, 1]
    if sp.A is not None:
        A = sps.csc_matrix(sp.A)
        lp.a_matrix_.format_ = hs.MatrixFormat.kColwise
        lp.a_matrix_.start_ = A.indptr
        lp.a_matrix_.index_ = A.indices
        lp.a_matrix_.value_ = A.data
        lp.row_lower_ = sp.bounds[:, 0]
        lp.row_upper_ = sp.bounds[:, 1]
    env.passModel(model)
    Jidx = J.astype(np.int32)
    stats = {"lp_seconds": 0.0, "n_lp": 0}

    def oracle(direction_w):
        """argmin over F of <direction_w, w(x)>; returns (status, w, x)."""
        env.changeColsCost(k, Jidx, direction_w * root)
        t = time.time()
        env.run()
        stats["lp_seconds"] += time.time() - t
        stats["n_lp"] += 1
        # Presolve stays on: on sample 10723 a cold solve with presolve took 2.5-8 s,
        # while a warm start without it took ~25-55 s on the full 74,813-column LP.
        status = env.modelStatusToString(env.getModelStatus())
        if status != "Optimal":
            return status, None, None
        x = np.array(env.getSolution().col_value)
        return status, root * x[J], x

    t0 = time.time()
    if x_start is not None:
        status, x0 = "Optimal", np.asarray(x_start, dtype=float)
        w0 = root * x0[J]
    else:
        status, w0, x0 = oracle(np.ones(k))
    if status != "Optimal":
        out = {"status": status, "x": None, "major": 0, "minor": 0, "gap": np.inf,
               "n_square_terms": k, **stats, "seconds": time.time() - t0}
        LAST_RUN.clear()
        LAST_RUN.update({kk: v for kk, v in out.items() if kk != "x"})
        return out
    W, X = [w0], [x0]                # corral: points in w coordinates and as full fluxes
    lam = np.array([1.0])
    w = w0.copy()
    major = minor = 0
    gap = np.inf
    best_lb = -np.inf
    while major < max_major:
        major += 1
        status, wq, xq = oracle(w)
        if status != "Optimal":
            break
        ww = float(w @ w)
        gap = (ww - float(w @ wq)) / max(ww, 1e-12)
        best_lb = max(best_lb, 0.5 * ww - max(0.0, ww - float(w @ wq)))
        cert = (0.5 * ww - best_lb) / max(0.5 * ww, 1e-12)
        if DEBUG:
            print(f"[highs_qp] major {major}: |w|^2 {ww:.10g}, rel gap {gap:.3g}, "
                  f"certified {cert:.3g}, corral {len(W)}, {time.time() - t0:.1f}s",
                  file=sys.stderr, flush=True)
        if gap <= GAP_TOL or (cert_tol is not None and cert <= cert_tol):
            break
        if any(np.allclose(wq, p, rtol=0, atol=1e-12) for p in W):
            break                     # oracle returned a corral point: no further progress
        W.append(wq)
        X.append(xq)
        lam = np.append(lam, 0.0)
        while True:                   # minor cycle
            minor += 1
            a = _affine_min(np.column_stack(W))
            if np.all(a > 1e-14):
                lam = a
                break
            neg = a <= 1e-14
            theta = np.min(lam[neg] / (lam[neg] - a[neg]))
            lam = theta * a + (1.0 - theta) * lam
            keep = lam > 1e-14
            W = [p for p, kp in zip(W, keep) if kp]
            X = [p for p, kp in zip(X, keep) if kp]
            lam = lam[keep] / lam[keep].sum()
        w = np.column_stack(W) @ lam
    x = np.column_stack(X) @ lam
    f = 0.5 * float(w @ w)
    cert = (f - best_lb) / max(f, 1e-12) if np.isfinite(best_lb) else np.inf
    target = cert_tol if cert_tol is not None else 2e-6
    if status != "Optimal":
        final = status
    elif cert <= target:
        final = "Optimal"
    else:
        final = "Iteration limit reached"           # maps to ITERATION_LIMIT, not OPTIMAL
    out = {"status": final, "x": x, "major": major, "minor": minor, "gap": gap,
           "f": f, "lower_bound": best_lb, "cert_gap_rel": cert,
           "n_square_terms": k, **stats, "seconds": time.time() - t0}
    LAST_RUN.clear()
    LAST_RUN.update({kk: v for kk, v in out.items() if kk != "x"})
    return out


def _solve_min_norm(self) -> None:
    sp = self.build()
    res = min_norm_point(sp, self.highs_settings())
    nv, nc = len(self.variables), len(self.constraints)
    self.status = res["status"]
    if res["x"] is None:
        self.obj_value = None
        return
    x = res["x"]
    cprimals = sps.csc_matrix(sp.A) @ x if sp.A is not None else np.zeros(nc)
    self._solution = {"x": x, "y": np.zeros(nc + nv)}
    self.primals = dict(zip(self.variables, x))
    # Duals are not produced by this method; micom's grow() does not read them.
    self.vduals = dict(zip(self.variables, np.zeros(nv)))
    self.cprimals = dict(zip(self.constraints, cprimals))
    self.duals = dict(zip(self.constraints, np.zeros(nc)))
    P = sps.csc_matrix(sp.P)
    self.obj_value = float(0.5 * x @ (P @ x) + sp.q @ x) * self.direction
    self.info = res


def repair(sp, x_hat, options):
    """Exactly feasible point closest (L1, on the squared variables) to x_hat.

    One HiGHS LP:  min sum_j d_j  s.t.  x in F,  -d_j <= x_j - x_hat_j <= d_j  for the k
    variables that carry the quadratic objective (the taxon growth variables). Every other
    flux is free to move to whatever keeps the point feasible. Returns dict(status, x,
    max_shift, seconds).
    """
    P = sps.csc_matrix(sp.P)
    J = np.flatnonzero(P.diagonal() > 0)
    nv, nc, k = sp.vbounds.shape[0], sp.bounds.shape[0], len(J)
    t0 = time.time()
    env = hs.Highs()
    env.passOptions(options)
    model = hs.HighsModel()
    lp = model.lp_
    lp.sense_ = hs.ObjSense.kMinimize
    lp.num_col_ = nv + k
    lp.num_row_ = nc + 2 * k
    lp.col_cost_ = np.concatenate([np.zeros(nv), np.ones(k)])
    lp.col_lower_ = np.concatenate([sp.vbounds[:, 0], np.zeros(k)])
    lp.col_upper_ = np.concatenate([sp.vbounds[:, 1], np.full(k, np.inf)])
    A = sps.csc_matrix(sp.A) if sp.A is not None else sps.csc_matrix((nc, nv))
    e = sps.csc_matrix((np.ones(k), (np.arange(k), J)), shape=(k, nv))
    I = sps.identity(k, format="csc")
    full = sps.vstack([
        sps.hstack([A, sps.csc_matrix((nc, k))]),
        sps.hstack([e, -I]),          # x_j - d_j <= x_hat_j
        sps.hstack([e, I]),           # x_j + d_j >= x_hat_j
    ]).tocsc()
    zh = x_hat[J]
    lp.a_matrix_.format_ = hs.MatrixFormat.kColwise
    lp.a_matrix_.start_ = full.indptr
    lp.a_matrix_.index_ = full.indices
    lp.a_matrix_.value_ = full.data
    lp.row_lower_ = np.concatenate([sp.bounds[:, 0], np.full(k, -np.inf), zh])
    lp.row_upper_ = np.concatenate([sp.bounds[:, 1], zh, np.full(k, np.inf)])
    env.passModel(model)
    env.run()
    status = env.modelStatusToString(env.getModelStatus())
    if status != "Optimal":
        return {"status": status, "x": None, "max_shift": np.inf,
                "seconds": time.time() - t0}
    x = np.array(env.getSolution().col_value)[:nv]
    return {"status": status, "x": x, "max_shift": float(np.max(np.abs(x[J] - zh))),
            "seconds": time.time() - t0}


def certificate(sp, x, options):
    """Certified optimality gap of a feasible x for  min 0.5 x'Px  over sp's LP.

    With w = sqrt(P_JJ) x_J, f(x) = 0.5 ||w||^2 and, for the LP vertex p minimising
    <w, p>, convexity gives  f(x) - f*  <=  <w, w - p>  (the Frank-Wolfe gap). One LP.
    Returns dict(status, f, gap_abs, gap_rel, lower_bound, seconds).
    """
    P = sps.csc_matrix(sp.P)
    diag = P.diagonal()
    J = np.flatnonzero(diag > 0)
    root = np.sqrt(diag[J])
    nv, nc = sp.vbounds.shape[0], sp.bounds.shape[0]
    t0 = time.time()
    env = hs.Highs()
    env.passOptions(options)
    model = hs.HighsModel()
    lp = model.lp_
    lp.sense_ = hs.ObjSense.kMinimize
    lp.num_col_ = nv
    lp.num_row_ = nc
    w = root * x[J]
    cost = np.zeros(nv)
    cost[J] = w * root
    lp.col_cost_ = cost
    lp.col_lower_ = sp.vbounds[:, 0]
    lp.col_upper_ = sp.vbounds[:, 1]
    if sp.A is not None:
        A = sps.csc_matrix(sp.A)
        lp.a_matrix_.format_ = hs.MatrixFormat.kColwise
        lp.a_matrix_.start_ = A.indptr
        lp.a_matrix_.index_ = A.indices
        lp.a_matrix_.value_ = A.data
        lp.row_lower_ = sp.bounds[:, 0]
        lp.row_upper_ = sp.bounds[:, 1]
    env.passModel(model)
    env.run()
    status = env.modelStatusToString(env.getModelStatus())
    f = 0.5 * float(w @ w)
    if status != "Optimal":
        return {"status": status, "f": f, "gap_abs": np.inf, "gap_rel": np.inf,
                "lower_bound": -np.inf, "seconds": time.time() - t0}
    p = np.array(env.getSolution().col_value)
    gap = max(0.0, float(w @ w) - float(w @ (root * p[J])))
    return {"status": status, "f": f, "gap_abs": gap, "gap_rel": gap / max(f, 1e-12),
            "lower_bound": f - gap, "seconds": time.time() - t0}


CERT_TOL = 0.01    # repair mode: polish until the certified relative gap is <= 1%
POLISH_MAX = 60    # ... or this many Wolfe major iterations (one LP each)


def _solve_osqp_repaired(self) -> None:
    """OSQP for the QP, `repair` to an exactly feasible point, then Wolfe-polish it.

    The polish starts from the repaired point (so it only ever improves on it), stops as
    soon as the certified relative optimality gap is <= CERT_TOL, and records that
    certificate. A sample that hits POLISH_MAX keeps its best point and is reported with
    status ITERATION_LIMIT and its (larger) certified gap -- never silently as optimal.
    """
    t0 = time.time()
    self.solve_osqp()
    osqp_status = self.status
    t_osqp = time.time() - t0
    sp = self.build()
    nv, nc = len(self.variables), len(self.constraints)
    x_hat = np.array([self.primals.get(v, 0.0) for v in self.variables])
    res = repair(sp, x_hat, self.highs_settings())
    pol = None
    if res["x"] is not None:
        pol = min_norm_point(sp, self.highs_settings(), x_start=res["x"],
                             cert_tol=CERT_TOL, max_major=POLISH_MAX)
    LAST_RUN.clear()
    LAST_RUN.update(mode="osqp+repair+polish", osqp_status=str(osqp_status),
                    osqp_seconds=t_osqp, repair_status=res["status"],
                    max_shift=res["max_shift"], repair_seconds=res["seconds"])
    if pol is not None:
        LAST_RUN.update(polish_status=pol["status"], polish_major=pol["major"],
                        polish_seconds=pol["seconds"], f=pol["f"],
                        lower_bound=pol["lower_bound"], gap_rel=pol["cert_gap_rel"])
    LAST_RUN["seconds"] = time.time() - t0
    if DEBUG:
        print(f"[highs_qp] osqp {osqp_status} ({t_osqp:.0f}s), repair {res['status']} "
              f"(shift {res['max_shift']:.3g}), polish "
              f"{pol['status'] if pol else '-'} in {pol['major'] if pol else 0} its, "
              f"certified gap {LAST_RUN.get('gap_rel', float('nan')):.3g}, "
              f"{time.time() - t0:.1f}s", file=sys.stderr, flush=True)
    if pol is not None and pol["x"] is not None:
        res = {"status": pol["status"], "x": pol["x"]}
    self.status = res["status"]
    if res["x"] is None:
        self.obj_value = None
        return
    x = res["x"]
    cprimals = sps.csc_matrix(sp.A) @ x if sp.A is not None else np.zeros(nc)
    self._solution = {"x": x, "y": np.zeros(nc + nv)}
    self.primals = dict(zip(self.variables, x))
    self.vduals = dict(zip(self.variables, np.zeros(nv)))
    self.cprimals = dict(zip(self.constraints, cprimals))
    self.duals = dict(zip(self.constraints, np.zeros(nc)))
    P = sps.csc_matrix(sp.P)
    self.obj_value = float(0.5 * x @ (P @ x) + sp.q @ x) * self.direction


MODE = "exact"   # "exact": Wolfe min-norm point; "repair": OSQP then exact LP repair


def _solve(self) -> None:
    if len(self.obj_quadratic_coefs) == 0:
        t0 = time.time()
        self.solve_highs()
        if DEBUG:
            print(f"[highs_qp] LP: {self.status}, {time.time() - t0:.1f}s",
                  file=sys.stderr, flush=True)
    elif len(self.integer_vars) > 0:
        raise hi.SolverError("MIQPs are not supported here.")
    elif MODE == "repair":
        _solve_osqp_repaired(self)
    else:
        _solve_min_norm(self)


def install(mode: str = "exact") -> None:
    """Patch optlang's hybrid solve in this process. mode: "exact" or "repair"."""
    global MODE
    if mode not in ("exact", "repair"):
        raise ValueError(mode)
    MODE = mode
    hi.HybridProblem.solve = _solve


def uninstall() -> None:
    hi.HybridProblem.solve = _ORIGINAL_SOLVE


def installed() -> bool:
    return hi.HybridProblem.solve is _solve


def describe() -> dict:
    import optlang
    names = {"exact": "Wolfe min-norm point over HiGHS LPs (v2/micom/highs_qp.py)",
             "repair": ("OSQP, then exact HiGHS LP repair and Wolfe polish to a "
                        "certified gap <= CERT_TOL (v2/micom/highs_qp.py)")}
    return {"qp_solver": names[MODE] if installed() else "OSQP",
            "gap_tol": GAP_TOL if installed() else None,
            "highspy": hs.Highs().version(), "optlang": optlang.__version__}
