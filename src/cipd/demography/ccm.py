"""V0 (Two-Sex CCM) & V1 (Age-Parity Cohort-Component Model, P-CCM) with Adjoint Gradients & Decomposition."""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np


@dataclass
class AccountingResult:
    pop_male: np.ndarray        # (T+1, A)
    pop_female: np.ndarray      # (T+1, A, P) for V1 or (T+1, A) for V0
    births_total: np.ndarray    # (T,)
    births_by_parity: np.ndarray | None  # (T, P_birth)
    deaths_total: np.ndarray    # (T,)
    in_migration: np.ndarray    # (T,)
    out_migration: np.ndarray   # (T,)
    max_abs_residual: float
    max_rel_residual: float


def run_ccm_v0(
    pop_m0: np.ndarray,
    pop_f0: np.ndarray,
    surv_m: np.ndarray,      # (T, A) survival prob from age a to a+1 (last age stays in open group)
    surv_f: np.ndarray,      # (T, A)
    asfr: np.ndarray,        # (T, A) age-specific fertility rate (births per woman)
    srb: float = 1.07,       # male births per female birth
    in_mig_m: np.ndarray | None = None,
    out_mig_m: np.ndarray | None = None,
    in_mig_f: np.ndarray | None = None,
    out_mig_f: np.ndarray | None = None,
) -> AccountingResult:
    """V0: Two-Sex Cohort-Component Model with exact demographic accounting."""
    T, A = asfr.shape
    pop_m = np.zeros((T + 1, A), dtype=float)
    pop_f = np.zeros((T + 1, A), dtype=float)
    pop_m[0] = np.asarray(pop_m0, dtype=float)
    pop_f[0] = np.asarray(pop_f0, dtype=float)

    if np.any(pop_m[0] < 0) or np.any(pop_f[0] < 0):
        raise ValueError("Initial population must be non-negative.")

    in_m = np.zeros((T, A)) if in_mig_m is None else np.asarray(in_mig_m, dtype=float)
    out_m = np.zeros((T, A)) if out_mig_m is None else np.asarray(out_mig_m, dtype=float)
    in_f = np.zeros((T, A)) if in_mig_f is None else np.asarray(in_mig_f, dtype=float)
    out_f = np.zeros((T, A)) if out_mig_f is None else np.asarray(out_mig_f, dtype=float)

    if np.any(in_m < 0) or np.any(out_m < 0) or np.any(in_f < 0) or np.any(out_f < 0):
        raise ValueError("In-migration and out-migration flows must be non-negative.")

    births_t = np.zeros(T)
    deaths_t = np.zeros(T)
    in_t = np.zeros(T)
    out_t = np.zeros(T)
    p_male_birth = srb / (1.0 + srb)

    for t in range(T):
        # Enforce out-migration <= available surviving population
        surv_m_pop = pop_m[t] * surv_m[t]
        surv_f_pop = pop_f[t] * surv_f[t]
        d_m = pop_m[t] * (1.0 - surv_m[t])
        d_f = pop_f[t] * (1.0 - surv_f[t])

        om = np.minimum(out_m[t], surv_m_pop)
        of = np.minimum(out_f[t], surv_f_pop)

        # Births from mid-year female population approximation or start-of-year exposure
        b_total = float(np.sum(pop_f[t] * asfr[t]))
        b_m = b_total * p_male_birth
        b_f = b_total * (1.0 - p_male_birth)

        # Age progression
        next_m = np.zeros(A)
        next_f = np.zeros(A)
        next_m[0] = b_m + in_m[t, 0]
        next_f[0] = b_f + in_f[t, 0]

        # Ages 0..A-3 advance to 1..A-2; Ages A-2 and A-1 merge into open age group A-1
        next_m[1 : A - 1] = surv_m_pop[0 : A - 2] - om[0 : A - 2] + in_m[t, 1 : A - 1]
        next_f[1 : A - 1] = surv_f_pop[0 : A - 2] - of[0 : A - 2] + in_f[t, 1 : A - 1]

        next_m[A - 1] = (
            surv_m_pop[A - 2] - om[A - 2]
            + surv_m_pop[A - 1] - om[A - 1]
            + in_m[t, A - 1]
        )
        next_f[A - 1] = (
            surv_f_pop[A - 2] - of[A - 2]
            + surv_f_pop[A - 1] - of[A - 1]
            + in_f[t, A - 1]
        )

        # Accounting: note om[A-1] is already subtracted above; verify exact conservation
        pop_m[t + 1] = next_m
        pop_f[t + 1] = next_f
        births_t[t] = b_total
        deaths_t[t] = float(np.sum(d_m) + np.sum(d_f))
        in_t[t] = float(np.sum(in_m[t]) + np.sum(in_f[t]))
        out_t[t] = float(np.sum(om) + np.sum(of))

    N_total = pop_m.sum(axis=1) + pop_f.sum(axis=1)
    residuals = (N_total[1:] - N_total[:-1]) - (births_t - deaths_t + in_t - out_t)
    max_abs = float(np.max(np.abs(residuals)))
    max_rel = float(np.max(np.abs(residuals) / np.maximum(N_total[:-1], 1e-12)))

    return AccountingResult(
        pop_male=pop_m,
        pop_female=pop_f,
        births_total=births_t,
        births_by_parity=None,
        deaths_total=deaths_t,
        in_migration=in_t,
        out_migration=out_t,
        max_abs_residual=max_abs,
        max_rel_residual=max_rel,
    )


def run_pccm_v1(
    pop_m0: np.ndarray,      # (A,)
    pop_f0: np.ndarray,      # (A, P) female stock by age a and parity r in {0, 1, ..., P-1}
    surv_m: np.ndarray,      # (T, A)
    surv_f: np.ndarray,      # (T, A)
    p_birth: np.ndarray,     # (T, A, P) probability of giving birth in year t given (a, r); terminal P-1 allows higher-order births
    srb: float = 1.07,
    twin_factor: float = 1.0,  # >= 1.0 multiplier for multiple live births per delivery event
) -> AccountingResult:
    """V1: Age-Parity Cohort-Component Model (P-CCM).
    - Preserves exact adult female conservation across parity transitions (r -> min(r+1, P-1)).
    - Terminal parity r = P-1 remains open to higher-order births (self-loop with birth generation).
    - Births of order 1 come from r=0, order 2 from r=1, order 3+ from r>=2.
    """
    T, A, P = p_birth.shape
    if np.any(p_birth < 0) or np.any(p_birth > 1.0 + 1e-12):
        raise ValueError("Parity birth probabilities must lie in [0, 1].")

    pop_m = np.zeros((T + 1, A), dtype=float)
    pop_f = np.zeros((T + 1, A, P), dtype=float)
    pop_m[0] = np.asarray(pop_m0, dtype=float)
    pop_f[0] = np.asarray(pop_f0, dtype=float)

    births_t = np.zeros(T)
    births_par = np.zeros((T, P))  # order 1..P (where last is P+)
    deaths_t = np.zeros(T)
    p_male_birth = srb / (1.0 + srb)

    for t in range(T):
        # Surviving adults
        sm = pop_m[t] * surv_m[t]
        sf = pop_f[t] * surv_f[t][:, None]  # (A, P)
        deaths_t[t] = float(np.sum(pop_m[t] * (1.0 - surv_m[t])) + np.sum(pop_f[t] * (1.0 - surv_f[t][:, None])))

        # Delivery events among surviving women (conserving adult headcount)
        # Note: births from women who die mid-year are folded into effective rate or excluded cleanly;
        # here delivery probability applies to surviving women so adult parity transitions sum identically to sf.
        deliveries = sf * p_birth[t]          # (A, P)
        stay = sf * (1.0 - p_birth[t])        # (A, P)

        # Birth counts by order (r=0 -> 1st birth, r=1 -> 2nd birth, ..., r=P-1 -> P+ birth)
        b_by_r = deliveries.sum(axis=0) * twin_factor  # (P,)
        births_par[t] = b_by_r
        b_total = float(b_by_r.sum())
        births_t[t] = b_total

        # Parity transition for surviving women: r -> r+1 for r < P-1; r=P-1 stays at P-1
        sf_next_parity = np.zeros((A, P), dtype=float)
        for r in range(P):
            sf_next_parity[:, r] += stay[:, r]
            r_next = min(r + 1, P - 1)
            sf_next_parity[:, r_next] += deliveries[:, r]

        # Age progression
        next_m = np.zeros(A)
        next_f = np.zeros((A, P))
        next_m[0] = b_total * p_male_birth
        next_f[0, 0] = b_total * (1.0 - p_male_birth)  # All newborn girls start at parity 0

        next_m[1 : A - 1] = sm[0 : A - 2]
        next_m[A - 1] = sm[A - 2] + sm[A - 1]

        next_f[1 : A - 1, :] = sf_next_parity[0 : A - 2, :]
        next_f[A - 1, :] = sf_next_parity[A - 2, :] + sf_next_parity[A - 1, :]

        pop_m[t + 1] = next_m
        pop_f[t + 1] = next_f

    N_total = pop_m.sum(axis=1) + pop_f.sum(axis=(1, 2))
    residuals = (N_total[1:] - N_total[:-1]) - (births_t - deaths_t)
    max_abs = float(np.max(np.abs(residuals)))
    max_rel = float(np.max(np.abs(residuals) / np.maximum(N_total[:-1], 1e-12)))

    return AccountingResult(
        pop_male=pop_m,
        pop_female=pop_f,
        births_total=births_t,
        births_by_parity=births_par,
        deaths_total=deaths_t,
        in_migration=np.zeros(T),
        out_migration=np.zeros(T),
        max_abs_residual=max_abs,
        max_rel_residual=max_rel,
    )


def adjoint_policy_gradient(
    pop_f0: np.ndarray,      # (A, P)
    surv_f: np.ndarray,      # (T, A)
    p_base: np.ndarray,      # (T, A, P)
    u: np.ndarray,           # (T, P) additive policy shift on fertility hazard by parity
    srb: float = 1.07,
    weights_t: np.ndarray | None = None,
) -> tuple[float, np.ndarray]:
    """Exact reverse-mode discrete adjoint gradient dJ/du for cumulative weighted births J = sum_t w_t B_t.
    Demonstrates O(1) backward pass independent of number of policy parameters.
    """
    T, A, P = p_base.shape
    w = np.ones(T) if weights_t is None else np.asarray(weights_t, dtype=float)
    p_girl = 1.0 / (1.0 + srb)

    # Forward pass storing states
    x = np.zeros((T + 1, A, P))
    x[0] = pop_f0
    p_eff = np.clip(p_base + u[:, None, :], 0.0, 1.0)
    mask_active = ((p_base + u[:, None, :]) > 0.0) & ((p_base + u[:, None, :]) < 1.0)
    B = np.zeros(T)

    for t in range(T):
        sf = x[t] * surv_f[t][:, None]
        deliv = sf * p_eff[t]
        stay = sf * (1.0 - p_eff[t])
        B[t] = deliv.sum()
        sf_next = np.zeros((A, P))
        for r in range(P):
            sf_next[:, r] += stay[:, r]
            sf_next[:, min(r + 1, P - 1)] += deliv[:, r]
        x[t + 1, 0, 0] = B[t] * p_girl
        x[t + 1, 1 : A - 1, :] = sf_next[0 : A - 2, :]
        x[t + 1, A - 1, :] = sf_next[A - 2, :] + sf_next[A - 1, :]

    J = float(np.dot(w, B))

    # Backward adjoint pass: lam[t] = dJ / dx[t]
    lam = np.zeros((T + 1, A, P))
    grad_u = np.zeros((T, P))

    for t in range(T - 1, -1, -1):
        # Direct + indirect sensitivity w.r.t. B[t]
        dJ_dBt = w[t] + p_girl * lam[t + 1, 0, 0]
        # Sensitivity w.r.t. sf_next[a, r]
        dJ_dsf_next = np.zeros((A, P))
        dJ_dsf_next[0 : A - 2, :] = lam[t + 1, 1 : A - 1, :]
        dJ_dsf_next[A - 2, :] += lam[t + 1, A - 1, :]
        dJ_dsf_next[A - 1, :] += lam[t + 1, A - 1, :]

        sf = x[t] * surv_f[t][:, None]
        # dJ/dp_eff[t, a, r] = sf[a,r] * (dJ_dBt + dJ_dsf_next[a, min(r+1,P-1)] - dJ_dsf_next[a, r])
        dJ_dp = np.zeros((A, P))
        dJ_dsf = np.zeros((A, P))
        for r in range(P):
            rn = min(r + 1, P - 1)
            dJ_dp[:, r] = sf[:, r] * (dJ_dBt + dJ_dsf_next[:, rn] - dJ_dsf_next[:, r])
            dJ_dsf[:, r] = (
                p_eff[t, :, r] * (dJ_dBt + dJ_dsf_next[:, rn])
                + (1.0 - p_eff[t, :, r]) * dJ_dsf_next[:, r]
            )
        grad_u[t] = np.sum(dJ_dp * mask_active[t], axis=0)
        lam[t] = dJ_dsf * surv_f[t][:, None]

    return J, grad_u


def decompose_births_direct_structure(
    f0: np.ndarray, x0: np.ndarray, f1: np.ndarray, x1: np.ndarray
) -> dict[str, float]:
    """Exact 3-way identity: dB = Direct + Structure + Interaction."""
    B0 = float(np.sum(f0 * x0))
    B1 = float(np.sum(f1 * x1))
    D = float(np.sum((f1 - f0) * x0))
    S = float(np.sum(f0 * (x1 - x0)))
    I = float(np.sum((f1 - f0) * (x1 - x0)))
    return {
        "B0": B0,
        "B1": B1,
        "delta_B": B1 - B0,
        "direct": D,
        "structure": S,
        "interaction": I,
        "identity_residual": abs((B1 - B0) - (D + S + I)),
    }


def shapley_quantum_tempo_structure(
    Q0: float, g0: np.ndarray, W0: np.ndarray,
    Q1: float, g1: np.ndarray, W1: np.ndarray,
) -> dict[str, float]:
    """Order-invariant Shapley decomposition of Delta B into Quantum (Q), Tempo (g), and Structure (W).
    Model: B(Q, g, W) = Q * sum_a g_a * W_a, where sum_a g_a = 1.
    """
    def B_fn(q, g, w):
        return float(q * np.dot(g, w))

    # 3! = 6 permutations of {Q, g, W}
    from itertools import permutations
    factors = ["Q", "g", "W"]
    phi = {k: 0.0 for k in factors}
    for perm in permutations(factors):
        state = {"Q": Q0, "g": g0, "W": W0}
        prev_val = B_fn(state["Q"], state["g"], state["W"])
        for f in perm:
            state[f] = {"Q": Q1, "g": g1, "W": W1}[f]
            cur_val = B_fn(state["Q"], state["g"], state["W"])
            phi[f] += (cur_val - prev_val) / 6.0
            prev_val = cur_val

    delta_B = B_fn(Q1, g1, W1) - B_fn(Q0, g0, W0)
    return {
        "B0": B_fn(Q0, g0, W0),
        "B1": B_fn(Q1, g1, W1),
        "delta_B": delta_B,
        "quantum": phi["Q"],
        "tempo": phi["g"],
        "structure": phi["W"],
        "identity_residual": abs(delta_B - (phi["Q"] + phi["g"] + phi["W"])),
    }
