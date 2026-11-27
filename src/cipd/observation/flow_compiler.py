"""V2: Multi-Source Observation Flow Compiler (OFR) — Time-Expanded Demographic Flow Polytope."""
from __future__ import annotations
from dataclasses import dataclass, field
import numpy as np
from scipy.sparse import csc_matrix, lil_matrix


@dataclass
class ObservationSpec:
    name: str
    indices: np.ndarray        # variable indices in z
    weights: np.ndarray        # coefficients in row of H
    lower: float               # lower bound y_low
    upper: float               # upper bound y_high
    cost: float = 1.0          # acquisition cost C_d for VOI experiments
    group: str = "default"     # observation group / statistical table ID


@dataclass
class FlowPolytope:
    """Polytope P = {z >= 0 : A_eq z = b_eq, A_ub z <= b_ub}."""
    n_vars: int
    T: int
    A: int
    P: int
    A_eq: csc_matrix
    b_eq: np.ndarray
    A_ub: csc_matrix
    b_ub: np.ndarray
    ub_labels: list[str]
    observations: list[ObservationSpec]
    # Variable indexing helpers
    x_offset: int
    z_stay_offset: int
    z_birth_offset: int
    n_x: int
    n_z_stay: int
    n_z_birth: int

    def idx_x(self, t: int, a: int, r: int) -> int:
        return self.x_offset + (t * self.A + a) * self.P + r

    def idx_stay(self, t: int, a: int, r: int) -> int:
        return self.z_stay_offset + (t * self.A + a) * self.P + r

    def idx_birth(self, t: int, a: int, r: int) -> int:
        return self.z_birth_offset + (t * self.A + a) * self.P + r


class DemographicFlowCompiler:
    """Compiles an age-parity female demographic flow network into sparse LP constraints.
    Variables z = [x (T+1, A, P), z_stay (T, A, P), z_birth (T, A, P)] >= 0.
    Conservation equations (Az = b):
      1. Source conservation: z_stay[t,a,r] + z_birth[t,a,r] = s[t,a] * x[t,a,r]
      2. Target conservation:
         - Age 0, parity 0: x[t+1, 0, 0] = p_girl * sum_{a,r} z_birth[t,a,r]
         - Age 0, parity > 0: x[t+1, 0, r] = 0
         - Age a in 1..A-2: x[t+1, a, r] = z_stay[t, a-1, r] + sum_{r': next(r')=r} z_birth[t, a-1, r']
         - Open age A-1:    x[t+1, A-1, r] = (flows from A-2 and A-1 into parity r)
    Rate box constraints (optional linear bounds p_lo * s * x <= z_birth <= p_hi * s * x).
    """

    def __init__(
        self,
        T: int,
        A: int,
        P: int,
        surv_f: np.ndarray,  # (T, A)
        srb: float = 1.07,
        fertile_age_idx: list[int] | None = None,
        p_lower: np.ndarray | None = None,  # (T, A, P)
        p_upper: np.ndarray | None = None,  # (T, A, P)
    ):
        self.T = T
        self.A = A
        self.P = P
        self.surv_f = np.asarray(surv_f, dtype=float)
        self.srb = float(srb)
        self.p_girl = 1.0 / (1.0 + self.srb)
        self.fertile_age_idx = set(range(A) if fertile_age_idx is None else fertile_age_idx)
        self.p_lower = p_lower
        self.p_upper = p_upper

        self.n_x = (T + 1) * A * P
        self.n_z_stay = T * A * P
        self.n_z_birth = T * A * P
        self.x_offset = 0
        self.z_stay_offset = self.n_x
        self.z_birth_offset = self.n_x + self.n_z_stay
        self.n_vars = self.n_x + self.n_z_stay + self.n_z_birth

    def idx_x(self, t: int, a: int, r: int) -> int:
        return self.x_offset + (t * self.A + a) * self.P + r

    def idx_stay(self, t: int, a: int, r: int) -> int:
        return self.z_stay_offset + (t * self.A + a) * self.P + r

    def idx_birth(self, t: int, a: int, r: int) -> int:
        return self.z_birth_offset + (t * self.A + a) * self.P + r

    def compile(self, observations: list[ObservationSpec]) -> FlowPolytope:
        T, A, P = self.T, self.A, self.P
        eq_rows = []
        eq_rhs = []

        # 1. Source conservation: z_stay + z_birth - s * x = 0
        for t in range(T):
            for a in range(A):
                s = self.surv_f[t, a]
                for r in range(P):
                    row = {
                        self.idx_stay(t, a, r): 1.0,
                        self.idx_birth(t, a, r): 1.0,
                        self.idx_x(t, a, r): -s,
                    }
                    eq_rows.append(row)
                    eq_rhs.append(0.0)
                    # Non-fertile ages cannot give birth
                    if a not in self.fertile_age_idx:
                        eq_rows.append({self.idx_birth(t, a, r): 1.0})
                        eq_rhs.append(0.0)

        # 2. Target conservation for t+1
        for t in range(T):
            # Age 0, parity 0: x[t+1, 0, 0] - p_girl * sum_{a,r} z_birth[t,a,r] = 0
            row0 = {self.idx_x(t + 1, 0, 0): 1.0}
            for a in self.fertile_age_idx:
                for r in range(P):
                    row0[self.idx_birth(t, a, r)] = -self.p_girl
            eq_rows.append(row0)
            eq_rhs.append(0.0)

            for r in range(1, P):
                eq_rows.append({self.idx_x(t + 1, 0, r): 1.0})
                eq_rhs.append(0.0)

            # Ages 1..A-2
            for a in range(1, A - 1):
                for r in range(P):
                    row = {self.idx_x(t + 1, a, r): 1.0, self.idx_stay(t, a - 1, r): -1.0}
                    if r > 0:
                        row[self.idx_birth(t, a - 1, r - 1)] = row.get(self.idx_birth(t, a - 1, r - 1), 0.0) - 1.0
                    if r == P - 1:
                        row[self.idx_birth(t, a - 1, P - 1)] = row.get(self.idx_birth(t, a - 1, P - 1), 0.0) - 1.0
                    eq_rows.append(row)
                    eq_rhs.append(0.0)

            # Open age A-1 receives from A-2 and A-1
            for r in range(P):
                row = {self.idx_x(t + 1, A - 1, r): 1.0}
                for a_src in (A - 2, A - 1):
                    row[self.idx_stay(t, a_src, r)] = row.get(self.idx_stay(t, a_src, r), 0.0) - 1.0
                    if r > 0:
                        row[self.idx_birth(t, a_src, r - 1)] = row.get(self.idx_birth(t, a_src, r - 1), 0.0) - 1.0
                    if r == P - 1:
                        row[self.idx_birth(t, a_src, P - 1)] = row.get(self.idx_birth(t, a_src, P - 1), 0.0) - 1.0
                eq_rows.append(row)
                eq_rhs.append(0.0)

        A_eq = lil_matrix((len(eq_rows), self.n_vars), dtype=float)
        for i, row in enumerate(eq_rows):
            for col, val in row.items():
                A_eq[i, col] = val
        b_eq = np.array(eq_rhs, dtype=float)

        # Inequalities: rate bounds + observation box constraints
        ub_rows = []
        ub_rhs = []
        ub_labels = []

        if self.p_lower is not None or self.p_upper is not None:
            for t in range(T):
                for a in self.fertile_age_idx:
                    s = self.surv_f[t, a]
                    for r in range(P):
                        if self.p_lower is not None and self.p_lower[t, a, r] > 0:
                            # p_lo * s * x - z_birth <= 0
                            ub_rows.append({
                                self.idx_x(t, a, r): self.p_lower[t, a, r] * s,
                                self.idx_birth(t, a, r): -1.0,
                            })
                            ub_rhs.append(0.0)
                            ub_labels.append(f"rate_lo_t{t}_a{a}_r{r}")
                        if self.p_upper is not None and self.p_upper[t, a, r] < 1.0:
                            # z_birth - p_hi * s * x <= 0
                            ub_rows.append({
                                self.idx_birth(t, a, r): 1.0,
                                self.idx_x(t, a, r): -self.p_upper[t, a, r] * s,
                            })
                            ub_rhs.append(0.0)
                            ub_labels.append(f"rate_hi_t{t}_a{a}_r{r}")

        for obs in observations:
            # H z <= upper
            row_hi = {int(idx): float(w) for idx, w in zip(obs.indices, obs.weights)}
            ub_rows.append(row_hi)
            ub_rhs.append(float(obs.upper))
            ub_labels.append(f"obs_hi:{obs.name}")
            # -H z <= -lower
            row_lo = {int(idx): -float(w) for idx, w in zip(obs.indices, obs.weights)}
            ub_rows.append(row_lo)
            ub_rhs.append(-float(obs.lower))
            ub_labels.append(f"obs_lo:{obs.name}")

        A_ub = lil_matrix((len(ub_rows), self.n_vars), dtype=float)
        for i, row in enumerate(ub_rows):
            for col, val in row.items():
                A_ub[i, col] = val
        b_ub = np.array(ub_rhs, dtype=float)

        return FlowPolytope(
            n_vars=self.n_vars,
            T=T, A=A, P=P,
            A_eq=A_eq.tocsc(),
            b_eq=b_eq,
            A_ub=A_ub.tocsc(),
            b_ub=b_ub,
            ub_labels=ub_labels,
            observations=list(observations),
            x_offset=self.x_offset,
            z_stay_offset=self.z_stay_offset,
            z_birth_offset=self.z_birth_offset,
            n_x=self.n_x,
            n_z_stay=self.n_z_stay,
            n_z_birth=self.n_z_birth,
        )

    # Convenience builders for standard demographic observations
    def obs_total_births(self, t: int, lower: float, upper: float, name: str | None = None, cost: float = 1.0, group: str = "total_births") -> ObservationSpec:
        idxs = [self.idx_birth(t, a, r) for a in self.fertile_age_idx for r in range(self.P)]
        return ObservationSpec(
            name=name or f"total_births_t{t}",
            indices=np.array(idxs, dtype=int),
            weights=np.ones(len(idxs), dtype=float),
            lower=lower, upper=upper, cost=cost, group=group,
        )

    def obs_parity_births(self, t: int, parities: list[int], lower: float, upper: float, name: str | None = None, cost: float = 1.5, group: str = "parity_births") -> ObservationSpec:
        idxs = [self.idx_birth(t, a, r) for a in self.fertile_age_idx for r in parities]
        return ObservationSpec(
            name=name or f"parity_births_t{t}_r{''.join(map(str, parities))}",
            indices=np.array(idxs, dtype=int),
            weights=np.ones(len(idxs), dtype=float),
            lower=lower, upper=upper, cost=cost, group=group,
        )

    def obs_age_births(self, t: int, ages: list[int], lower: float, upper: float, name: str | None = None, cost: float = 1.5, group: str = "age_births") -> ObservationSpec:
        idxs = [self.idx_birth(t, a, r) for a in ages if a in self.fertile_age_idx for r in range(self.P)]
        return ObservationSpec(
            name=name or f"age_births_t{t}_a{ages[0]}_{ages[-1]}",
            indices=np.array(idxs, dtype=int),
            weights=np.ones(len(idxs), dtype=float),
            lower=lower, upper=upper, cost=cost, group=group,
        )

    def obs_age_stock(self, t: int, ages: list[int], lower: float, upper: float, name: str | None = None, cost: float = 1.2, group: str = "age_stock") -> ObservationSpec:
        idxs = [self.idx_x(t, a, r) for a in ages for r in range(self.P)]
        return ObservationSpec(
            name=name or f"age_stock_t{t}_a{ages[0]}_{ages[-1]}",
            indices=np.array(idxs, dtype=int),
            weights=np.ones(len(idxs), dtype=float),
            lower=lower, upper=upper, cost=cost, group=group,
        )

    def obs_parity_stock(self, t: int, ages: list[int], parities: list[int], lower: float, upper: float, name: str | None = None, cost: float = 2.0, group: str = "parity_stock") -> ObservationSpec:
        idxs = [self.idx_x(t, a, r) for a in ages for r in parities]
        return ObservationSpec(
            name=name or f"parity_stock_t{t}_a{ages[0]}_{ages[-1]}_r{''.join(map(str, parities))}",
            indices=np.array(idxs, dtype=int),
            weights=np.ones(len(idxs), dtype=float),
            lower=lower, upper=upper, cost=cost, group=group,
        )
