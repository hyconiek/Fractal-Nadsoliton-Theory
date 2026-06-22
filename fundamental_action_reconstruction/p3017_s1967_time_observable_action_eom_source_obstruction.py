#!/usr/bin/env python3
"""P3017/S1967: time-observable action/EOM source obstruction.

After P3016 closes the U(12)-quotient successor route, this attacks exactly one
remaining atom named by the guardrail: a unit-bearing action/EOM source for one
already typed observable.  The chosen typed observable is the P3014 strict-kernel
finite-difference scalar T_K(d).

The constructed candidate is a discrete quadratic variational source

    S_formal[q] = (lambda/2) * sum_d (q_{d+1} - q_d)^2,

with cyclic labels.  This is a real formal EOM object: its Hessian is the cycle
Laplacian with rank 11 and nullity 1.  The finite obstruction is twofold: the
P3014 observable does not solve the formal EOM, and the action remains non-unit-
bearing because no strict source fixes lambda, an action unit, a clock unit, or a
Hamiltonian normalization.  Therefore no physical EOM/Hamiltonian, L_total, or
ToE closure is exported.
"""
from __future__ import annotations

import hashlib, json, math
from fractions import Fraction
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3016_s1966_quotient_clock_successor_semigroup_exhaustion import OUT as P3016

OUT = GEN / "p3017_s1967_time_observable_action_eom_source_obstruction.json"
MD = GEN / "p3017_s1967_time_observable_action_eom_source_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

N = 12
OMEGA = 0.18575
PHI = 0.16250
BETA = 1.0
ETA = 1.8


def k_strict(d: int) -> float:
    return math.cos(OMEGA * d + PHI) / (1.0 + BETA * (float(d) ** ETA))


def t_obs(d: int) -> float:
    return k_strict(d + 1) - k_strict(d)


def cycle_laplacian(n: int) -> list[list[Fraction]]:
    mat = [[Fraction(0) for _ in range(n)] for _ in range(n)]
    for i in range(n):
        mat[i][i] = Fraction(2)
        mat[i][(i - 1) % n] = Fraction(-1)
        mat[i][(i + 1) % n] = Fraction(-1)
    return mat


def rank_fraction(matrix: list[list[Fraction]]) -> int:
    a = [row[:] for row in matrix]
    rows, cols = len(a), len(a[0]) if a else 0
    rank = 0
    pivot_col = 0
    while rank < rows and pivot_col < cols:
        pivot = next((r for r in range(rank, rows) if a[r][pivot_col] != 0), None)
        if pivot is None:
            pivot_col += 1
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        pv = a[rank][pivot_col]
        a[rank] = [x / pv for x in a[rank]]
        for r in range(rows):
            if r != rank and a[r][pivot_col] != 0:
                factor = a[r][pivot_col]
                a[r] = [x - factor * y for x, y in zip(a[r], a[rank])]
        rank += 1
        pivot_col += 1
    return rank


def mat_vec_float(mat: list[list[Fraction]], vec: list[float]) -> list[float]:
    return [sum(float(c) * v for c, v in zip(row, vec)) for row in mat]


def build_action_eom_matrix() -> dict[str, Any]:
    lap = cycle_laplacian(N)
    rank = rank_fraction(lap)
    nullity = N - rank
    q = [t_obs(d) for d in range(1, N + 1)]
    residual = mat_vec_float(lap, q)
    residual_rows = [
        {"d": i + 1, "T_K(d)": round(q[i], 12), "laplacian_residual": round(residual[i], 12), "zero_within_1e-12": abs(residual[i]) <= 1e-12}
        for i in range(N)
    ]
    nonzero_residual_count = sum(1 for r in residual_rows if not r["zero_within_1e-12"])
    unit_rows = [
        {"unit_obligation": "lambda_action_per_observable_square", "satisfied": False, "detail": "requires strict source for [lambda]=[action]/[T_K]^2"},
        {"unit_obligation": "clock_unit_for_Hamiltonian", "satisfied": False, "detail": "requires strict clock unit to convert action to energy/Hamiltonian"},
        {"unit_obligation": "observable_unit_normalization", "satisfied": False, "detail": "T_K is dimensionless kernel-difference unless a readout unit is exported"},
        {"unit_obligation": "physical_boundary_conditions", "satisfied": False, "detail": "cyclic formal boundary is a computation choice, not a strict physical time boundary theorem"},
    ]
    obligations = [
        {"obligation": "typed_observable_input", "satisfied": True, "detail": "uses the P3014 scalar T_K(d)"},
        {"obligation": "formal_variational_eom", "satisfied": rank == N - 1 and nullity == 1, "detail": f"cycle Laplacian rank/nullity = {rank}/{nullity}"},
        {"obligation": "observable_solves_formal_eom", "satisfied": nonzero_residual_count == 0, "detail": f"nonzero residual rows = {nonzero_residual_count}/{N}"},
        {"obligation": "unit_bearing_action_source", "satisfied": all(row["satisfied"] for row in unit_rows), "detail": "lambda/action/clock/observable units are not strict-sourced"},
        {"obligation": "hamiltonian_export", "satisfied": False, "detail": "no clock unit or canonical momentum normalization is exported"},
    ]
    return {
        "object": "TimeObservableQuadraticActionEOM_UnitSourceObstructionMatrix",
        "formal_action": "S_formal[q]=(lambda/2)*sum_d (q_{d+1}-q_d)^2 on cyclic Z/12 labels",
        "typed_observable": "q_d=T_K(d)=K_strict_gate(d+1)-K_strict_gate(d)",
        "laplacian_rank": rank,
        "laplacian_nullity": nullity,
        "residual_rows": residual_rows,
        "nonzero_residual_count": nonzero_residual_count,
        "unit_obligation_rows": unit_rows,
        "proof_obligations": obligations,
        "accepted_as_unit_bearing_action_eom_source": all(row["satisfied"] for row in obligations),
    }


def build_payload(p3016_path: Any) -> dict[str, Any]:
    matrix = build_action_eom_matrix()
    return {
        "status": "P3017_TIME_OBSERVABLE_ACTION_EOM_SOURCE_OBSTRUCTION_NO_CLOSURE",
        "input_hashes": {"P3016": hashlib.sha256(p3016_path.read_bytes()).hexdigest() if p3016_path.exists() else None},
        "constructed_theoretical_objects": matrix,
        "finite_certificate": {
            "laplacian_rank": matrix["laplacian_rank"],
            "laplacian_nullity": matrix["laplacian_nullity"],
            "residual_row_count": len(matrix["residual_rows"]),
            "nonzero_residual_count": matrix["nonzero_residual_count"],
            "unit_obligation_count": len(matrix["unit_obligation_rows"]),
            "satisfied_unit_obligation_count": sum(1 for row in matrix["unit_obligation_rows"] if row["satisfied"]),
            "accepted_as_unit_bearing_action_eom_source": matrix["accepted_as_unit_bearing_action_eom_source"],
        },
        "decision": {
            "breakthrough": "A formal quadratic action/EOM object for the typed time observable was constructed: the cyclic Laplacian has rank 11 and nullity 1.  The bounded obstruction is that T_K does not solve the formal EOM and the action is not unit-bearing because lambda, clock, observable, and Hamiltonian normalizations remain unsourced.",
            "negative_export_flags": {k: False for k in ["unit_bearing_action_eom_source_exported", "hamiltonian_exported", "time_observable_generator_exported", "time_arrow_exported", "ltotal_exported", "observed_physics_exported", "selector_source_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not promote formal quadratic actions to physical EOM/Hamiltonian closure.  The next proof-grade move should attack exactly one unit source atom, preferably a strict source for lambda/action-unit normalization for the typed observable, or introduce a new strict time-order object with its own unit theorem; keep selector, bridge, role-transfer, L_total, observed-physics, and ToE closure blocked.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = [
        "# P3017/S1967 time-observable action/EOM source obstruction", "",
        f"Status: `{payload['status']}`", "", "## Finite certificate",
        f"- Laplacian rank/nullity: `{c['laplacian_rank']}/{c['laplacian_nullity']}`",
        f"- residual rows/nonzero: `{c['residual_row_count']}/{c['nonzero_residual_count']}`",
        f"- unit obligations satisfied/total: `{c['satisfied_unit_obligation_count']}/{c['unit_obligation_count']}`",
        f"- accepted as unit-bearing action/EOM source: `{c['accepted_as_unit_bearing_action_eom_source']}`", "",
        "## Decision", payload["decision"]["breakthrough"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    read_json(P3016)
    payload = build_payload(P3016)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3017/S1967 time-observable action/EOM source obstruction", "## P3017/S1967 time-observable action/EOM source obstruction\n\n`P3017/S1967` pivots from quotient time-order replay to exactly one unit-bearing action/EOM source atom for an already typed observable.  It constructs the formal quadratic action `S_formal[q]=(lambda/2) sum_d (q_{d+1}-q_d)^2` for `q_d=T_K(d)`.  The resulting cyclic Laplacian is a real formal EOM object with rank/nullity `11/1`, but the audited `T_K` vector has nonzero EOM residuals and the action remains non-unit-bearing: no strict source fixes `lambda`, an action unit, a clock unit, observable-unit normalization, or Hamiltonian normalization.  Therefore no unit-bearing action/EOM source, Hamiltonian, time arrow, `L_total`, observed-physics export, bridge/role transfer, selector closure, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3017/S1967 formal time-observable action `L_total` guard", "## P3017/S1967 formal time-observable action `L_total` guard\n\n`P3017/S1967` adds no physical `L_total` term.  Although the formal quadratic action for `T_K` has a valid cyclic Laplacian EOM with rank/nullity `11/1`, the observable does not solve that EOM and the construction lacks strict `lambda`, action-unit, clock-unit, observable-unit, Hamiltonian, selector-source, bridge, and role-transfer theorems.\n")
    append_once(AGENTS, "Current formal time-observable action/EOM guardrail (P3017/S1967, 2026-06-22)", "## Current formal time-observable action/EOM guardrail (P3017/S1967, 2026-06-22)\n\n- P3017 constructs a formal quadratic action/EOM candidate for the typed P3014 time observable `T_K(d)`.\n- The finite formal EOM object is real: the cyclic Laplacian has rank/nullity `11/1`; however `T_K` has nonzero residuals and no strict source fixes `lambda`, action unit, clock unit, observable-unit normalization, or Hamiltonian normalization.\n- Do not promote formal quadratic time-observable actions to unit-bearing EOM/Hamiltonian, `L_total`, observed-physics, selector, bridge/role-transfer, or ToE closure.\n- The next honest move should attack exactly one unit source atom, preferably strict `lambda`/action-unit normalization for a typed observable, or introduce a genuinely new strict time-order object with its own unit theorem.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
