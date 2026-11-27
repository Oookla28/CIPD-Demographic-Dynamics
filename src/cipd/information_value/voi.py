"""V6: Goal-Oriented Data Value & Two-Trajectory Observation Separation (VOI-CPI)."""
from __future__ import annotations
from dataclasses import dataclass
from itertools import combinations
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import csc_matrix, hstack, vstack
from cipd.observation.flow_compiler import DemographicFlowCompiler, FlowPolytope, ObservationSpec
from cipd.identification.partial_id import solve_identified_bounds


@dataclass
class TwoTrajectoryResult:
    max_separation: float
    z1: np.ndarray | None
    z2: np.ndarray | None
    status: int


def solve_two_trajectory_separation(
    base_poly: FlowPolytope,
    c: np.ndarray,
    obs_list: list[ObservationSpec],
    deltas: np.ndarray,
) -> TwoTrajectoryResult:
    """Two-Trajectory Observation Separation LP (Contribution C):
      max_{z, z' in P_base}  c^T (z - z')
      s.t.  -2*delta_k <= h_k^T (z - z') <= 2*delta_k   for all k in obs_list.

    Measures the worst-case target discrepancy between ANY two demographically valid trajectories
    that can produce a common noisy observation vector under box measurement errors |eps_k| <= delta_k.
    Notice this depends ONLY on the observation design (H, delta), NOT on unknown future realized values y!
    """
    n = base_poly.n_vars
    c_ext = np.concatenate([-c, c])  # minimize -c^T(z - z')

    # Block diagonal equality constraints: A_eq z = b_eq and A_eq z' = b_eq
    Z_eq = csc_matrix(base_poly.A_eq.shape)
    A_eq = vstack([
        hstack([base_poly.A_eq, Z_eq], format="csc"),
        hstack([Z_eq, base_poly.A_eq], format="csc"),
    ], format="csc")
    b_eq = np.concatenate([base_poly.b_eq, base_poly.b_eq])

    # Block diagonal base inequality constraints
    Z_ub = csc_matrix(base_poly.A_ub.shape)
    ub_blocks = [
        hstack([base_poly.A_ub, Z_ub], format="csc"),
        hstack([Z_ub, base_poly.A_ub], format="csc"),
    ]
    b_ub_blocks = [base_poly.b_ub, base_poly.b_ub]

    # Coupling observation separation constraints: |h_k^T (z - z')| <= 2 * delta_k
    if len(obs_list) > 0:
        rows, cols, data = [], [], []
        rhs = []
        row_idx = 0
        for obs, d in zip(obs_list, deltas):
            for idx, w in zip(obs.indices, obs.weights):
                # h^T z - h^T z' <= 2d
                rows.extend([row_idx, row_idx])
                cols.extend([int(idx), n + int(idx)])
                data.extend([float(w), -float(w)])
                # -h^T z + h^T z' <= 2d
                rows.extend([row_idx + 1, row_idx + 1])
                cols.extend([int(idx), n + int(idx)])
                data.extend([-float(w), float(w)])
            rhs.extend([2.0 * float(d), 2.0 * float(d)])
            row_idx += 2
        A_sep = csc_matrix((data, (rows, cols)), shape=(row_idx, 2 * n))
        ub_blocks.append(A_sep)
        b_ub_blocks.append(np.array(rhs, dtype=float))

    A_ub = vstack(ub_blocks, format="csc")
    b_ub = np.concatenate(b_ub_blocks)

    res = linprog(c=c_ext, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq, bounds=(0.0, None), method="highs")
    if res.status != 0:
        return TwoTrajectoryResult(max_separation=float("nan"), z1=None, z2=None, status=res.status)
    return TwoTrajectoryResult(
        max_separation=float(-res.fun),
        z1=res.x[:n].copy(),
        z2=res.x[n:].copy(),
        status=res.status,
    )


def select_observations_greedy_voi(
    compiler: DemographicFlowCompiler,
    base_obs: list[ObservationSpec],
    candidate_obs: list[ObservationSpec],
    c: np.ndarray,
    budget: float,
    use_two_trajectory: bool = True,
) -> dict:
    """Goal-oriented greedy observation selection under acquisition cost budget.
    Supports both:
      - Two-Trajectory Worst-Case Separation reduction (prior to observing y_d)
      - Ex-post exact width reduction W(P) - W(P cap C_d)
    """
    selected: list[ObservationSpec] = []
    remaining = list(candidate_obs)
    spent = 0.0

    def eval_width(obs_subset: list[ObservationSpec]) -> float:
        if use_two_trajectory:
            poly_base = compiler.compile(base_obs)
            deltas = np.array([0.5 * (o.upper - o.lower) for o in obs_subset], dtype=float)
            res = solve_two_trajectory_separation(poly_base, c, obs_subset, deltas)
            return res.max_separation
        else:
            poly = compiler.compile(base_obs + obs_subset)
            cert = solve_identified_bounds(poly, c)
            return cert.width

    init_w = eval_width([])
    cur_w = init_w
    history = [{"step": 0, "added": None, "cost_spent": 0.0, "width": cur_w, "reduction": 0.0}]

    step = 0
    while True:
        best_score = -1.0
        best_cand = None
        best_new_w = cur_w
        for cand in remaining:
            if spent + cand.cost > budget + 1e-12:
                continue
            new_w = eval_width(selected + [cand])
            gain = max(0.0, cur_w - new_w)
            score = gain / max(cand.cost, 1e-9)
            if score > best_score:
                best_score = score
                best_cand = cand
                best_new_w = new_w
        if best_cand is None or best_score <= 1e-12:
            break
        selected.append(best_cand)
        remaining.remove(best_cand)
        spent += best_cand.cost
        cur_w = best_new_w
        step += 1
        history.append({
            "step": step,
            "added": best_cand.name,
            "cost_spent": spent,
            "width": cur_w,
            "reduction": 1.0 - cur_w / max(init_w, 1e-12),
        })

    return {
        "selected_names": [o.name for o in selected],
        "cost_spent": spent,
        "initial_width": init_w,
        "final_width": cur_w,
        "reduction": 1.0 - cur_w / max(init_w, 1e-12),
        "history": history,
    }


def select_observations_exhaustive_optimal(
    compiler: DemographicFlowCompiler,
    base_obs: list[ObservationSpec],
    candidate_obs: list[ObservationSpec],
    c: np.ndarray,
    budget: float,
    use_two_trajectory: bool = True,
) -> dict:
    """Exact exhaustive enumeration benchmark for small candidate sets (<= 10 candidates)."""
    poly_base = compiler.compile(base_obs)
    init_w = (
        solve_two_trajectory_separation(poly_base, c, [], np.zeros(0)).max_separation
        if use_two_trajectory
        else solve_identified_bounds(poly_base, c).width
    )
    best_w = init_w
    best_sub: list[ObservationSpec] = []
    best_cost = 0.0

    for k in range(1, len(candidate_obs) + 1):
        for comb in combinations(candidate_obs, k):
            cost = sum(o.cost for o in comb)
            if cost > budget + 1e-12:
                continue
            if use_two_trajectory:
                deltas = np.array([0.5 * (o.upper - o.lower) for o in comb], dtype=float)
                w = solve_two_trajectory_separation(poly_base, c, list(comb), deltas).max_separation
            else:
                w = solve_identified_bounds(compiler.compile(base_obs + list(comb)), c).width
            if w < best_w - 1e-12:
                best_w = w
                best_sub = list(comb)
                best_cost = cost

    return {
        "selected_names": [o.name for o in best_sub],
        "cost_spent": best_cost,
        "initial_width": init_w,
        "final_width": best_w,
        "reduction": 1.0 - best_w / max(init_w, 1e-12),
    }
