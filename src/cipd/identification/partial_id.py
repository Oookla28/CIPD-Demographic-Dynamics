"""V3: Goal-Oriented Partial Identification Engine (CPI) with Independent Primal-Dual Certificates."""
from __future__ import annotations
from dataclasses import dataclass
import time
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import csc_matrix, eye, hstack, vstack
from cipd.observation.flow_compiler import FlowPolytope


@dataclass
class BoundCertificate:
    target_name: str
    lower: float
    upper: float
    width: float
    status_low: int
    status_high: int
    primal_residual_low: float
    primal_residual_high: float
    dual_residual_low: float
    dual_residual_high: float
    duality_gap_low: float
    duality_gap_high: float
    top_binding_low: list[tuple[str, float]]
    top_binding_high: list[tuple[str, float]]
    z_low: np.ndarray | None
    z_high: np.ndarray | None
    runtime_sec: float


def _verify_lp_certificate(
    c: np.ndarray,
    A_eq: csc_matrix,
    b_eq: np.ndarray,
    A_ub: csc_matrix,
    b_ub: np.ndarray,
    res,
) -> tuple[float, float, float]:
    """Independently verify primal residual, dual residual, and normalized primal-dual gap.
    For min c^T z s.t. A_eq z = b_eq, A_ub z <= b_ub, z >= 0:
    SciPy HiGHS returns marginals y_eq (unrestricted), y_ub (<= 0) such that:
    dual objective = b_eq^T y_eq + b_ub^T y_ub, and A_eq^T y_eq + A_ub^T y_ub <= c.
    """
    z = res.x
    p_eq = np.max(np.abs(A_eq @ z - b_eq)) if A_eq.shape[0] > 0 else 0.0
    p_ub = np.max(np.maximum(A_ub @ z - b_ub, 0.0)) if A_ub.shape[0] > 0 else 0.0
    p_nn = np.max(np.maximum(-z, 0.0))
    primal_res = float(max(p_eq, p_ub, p_nn))

    y_eq = res.eqlin.marginals if A_eq.shape[0] > 0 else np.zeros(0)
    y_ub = res.ineqlin.marginals if A_ub.shape[0] > 0 else np.zeros(0)
    dual_obj = float((b_eq @ y_eq if len(y_eq) else 0.0) + (b_ub @ y_ub if len(y_ub) else 0.0))
    primal_obj = float(c @ z)

    lhs = np.zeros_like(c)
    if A_eq.shape[0] > 0:
        lhs += A_eq.T @ y_eq
    if A_ub.shape[0] > 0:
        lhs += A_ub.T @ y_ub
    # Dual feasibility: lhs <= c and y_ub <= 0
    d_slack = np.max(np.maximum(lhs - c, 0.0))
    d_sign = np.max(np.maximum(y_ub, 0.0)) if len(y_ub) else 0.0
    dual_res = float(max(d_slack, d_sign))

    denom = max(1.0, abs(primal_obj), abs(dual_obj))
    gap = float(abs(primal_obj - dual_obj) / denom)
    return primal_res, dual_res, gap


def solve_identified_bounds(
    poly: FlowPolytope,
    c: np.ndarray,
    target_name: str = "target",
    top_k_binding: int = 5,
) -> BoundCertificate:
    """Compute certified [lower, upper] bounds for Y = c^T z over FlowPolytope via HiGHS."""
    t0 = time.perf_counter()
    c = np.asarray(c, dtype=float)
    res_lo = linprog(
        c=c, A_ub=poly.A_ub, b_ub=poly.b_ub, A_eq=poly.A_eq, b_eq=poly.b_eq,
        bounds=(0.0, None), method="highs",
    )
    res_hi = linprog(
        c=-c, A_ub=poly.A_ub, b_ub=poly.b_ub, A_eq=poly.A_eq, b_eq=poly.b_eq,
        bounds=(0.0, None), method="highs",
    )
    dt = time.perf_counter() - t0

    if res_lo.status != 0 or res_hi.status != 0:
        return BoundCertificate(
            target_name=target_name,
            lower=float("nan"), upper=float("nan"), width=float("nan"),
            status_low=res_lo.status, status_high=res_hi.status,
            primal_residual_low=float("nan"), primal_residual_high=float("nan"),
            dual_residual_low=float("nan"), dual_residual_high=float("nan"),
            duality_gap_low=float("nan"), duality_gap_high=float("nan"),
            top_binding_low=[], top_binding_high=[],
            z_low=None, z_high=None, runtime_sec=dt,
        )

    pr_lo, dr_lo, gap_lo = _verify_lp_certificate(c, poly.A_eq, poly.b_eq, poly.A_ub, poly.b_ub, res_lo)
    pr_hi, dr_hi, gap_hi = _verify_lp_certificate(-c, poly.A_eq, poly.b_eq, poly.A_ub, poly.b_ub, res_hi)

    mu_lo = np.abs(res_lo.ineqlin.marginals)
    mu_hi = np.abs(res_hi.ineqlin.marginals)
    top_lo_idx = np.argsort(-mu_lo)[:top_k_binding]
    top_hi_idx = np.argsort(-mu_hi)[:top_k_binding]
    top_binding_low = [(poly.ub_labels[i], float(mu_lo[i])) for i in top_lo_idx if mu_lo[i] > 1e-12]
    top_binding_high = [(poly.ub_labels[i], float(mu_hi[i])) for i in top_hi_idx if mu_hi[i] > 1e-12]

    lo_val = float(res_lo.fun)
    hi_val = float(-res_hi.fun)

    return BoundCertificate(
        target_name=target_name,
        lower=lo_val,
        upper=hi_val,
        width=hi_val - lo_val,
        status_low=res_lo.status,
        status_high=res_hi.status,
        primal_residual_low=pr_lo,
        primal_residual_high=pr_hi,
        dual_residual_low=dr_lo,
        dual_residual_high=dr_hi,
        duality_gap_low=gap_lo,
        duality_gap_high=gap_hi,
        top_binding_low=top_binding_low,
        top_binding_high=top_binding_high,
        z_low=res_lo.x.copy(),
        z_high=res_hi.x.copy(),
        runtime_sec=dt,
    )


def diagnose_infeasibility(poly: FlowPolytope) -> dict:
    """L1-slack feasibility diagnostic: min sum(s) s.t. A_eq z = b_eq, A_ub z - s <= b_ub, z,s >= 0."""
    m_ub = poly.A_ub.shape[0]
    n = poly.n_vars
    c = np.concatenate([np.zeros(n), np.ones(m_ub)])
    A_eq_ext = hstack([poly.A_eq, csc_matrix((poly.A_eq.shape[0], m_ub))], format="csc")
    A_ub_ext = hstack([poly.A_ub, -eye(m_ub, format="csc")], format="csc")
    res = linprog(c=c, A_ub=A_ub_ext, b_ub=poly.b_ub, A_eq=A_eq_ext, b_eq=poly.b_eq, bounds=(0.0, None), method="highs")
    if res.status != 0:
        return {"feasible": False, "status": res.status, "total_slack": float("nan"), "conflicts": []}
    slacks = res.x[n:]
    conflicts = [(poly.ub_labels[i], float(slacks[i])) for i in np.argsort(-slacks) if slacks[i] > 1e-8]
    return {
        "feasible": float(res.fun) <= 1e-8,
        "status": res.status,
        "total_slack": float(res.fun),
        "conflicts": conflicts[:10],
    }


def solve_representative_map(poly: FlowPolytope, z_prior: np.ndarray) -> np.ndarray | None:
    """L1-regularized MAP-style representative point inside feasible polytope P: min ||z - z_prior||_1."""
    n = poly.n_vars
    # min sum t_i s.t. z - t <= z_prior, -z - t <= -z_prior, z in P, t >= 0
    c = np.concatenate([np.zeros(n), np.ones(n)])
    I_n = eye(n, format="csc")
    A_eq_ext = hstack([poly.A_eq, csc_matrix((poly.A_eq.shape[0], n))], format="csc")
    A_ub_top = hstack([poly.A_ub, csc_matrix((poly.A_ub.shape[0], n))], format="csc")
    A_ub_p = hstack([I_n, -I_n], format="csc")
    A_ub_m = hstack([-I_n, -I_n], format="csc")
    A_ub_ext = vstack([A_ub_top, A_ub_p, A_ub_m], format="csc")
    b_ub_ext = np.concatenate([poly.b_ub, z_prior, -z_prior])
    res = linprog(c=c, A_ub=A_ub_ext, b_ub=b_ub_ext, A_eq=A_eq_ext, b_eq=poly.b_eq, bounds=(0.0, None), method="highs")
    return res.x[:n].copy() if res.status == 0 else None
