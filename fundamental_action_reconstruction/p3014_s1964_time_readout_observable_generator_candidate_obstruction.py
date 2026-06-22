#!/usr/bin/env python3
"""P3014/S1964: time-readout observable-generator candidate obstruction.

This is the next bounded, proof-grade step after P3013.  It attacks exactly one
readout atom: an observer-independent observable generator for the `time` row.
The new candidate is deliberately narrow: use the strict kernel values on the
Z/12Z label cycle to form a signed finite-difference clock observable

    T_K(d) = K_strict_gate(d + 1) - K_strict_gate(d).

The formula has explicit input/output types and needs no observer premise, but a
ToE-grade time observable must also be compatible with relabeling/unit action,
select a directed arrow without importing selector closure, and be installable
into a unit-bearing EOM/Hamiltonian readout.  The finite matrix below checks
those obligations and records the bounded no-go rather than promoting a closure.
"""
from __future__ import annotations

import hashlib, json, math
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3013_s1963_observer_physics_readout_strict_kernel_selector_obligation_matrix import OUT as P3013

OUT = GEN / "p3014_s1964_time_readout_observable_generator_candidate_obstruction.json"
MD = GEN / "p3014_s1964_time_readout_observable_generator_candidate_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

N = 12
UNITS = [1, 5, 7, 11]
OMEGA = 0.18575
PHI = 0.16250
BETA = 1.0
ETA = 1.8


def k_strict(d: int) -> float:
    return math.cos(OMEGA * d + PHI) / (1.0 + BETA * (float(d) ** ETA))


def t_obs(d: int) -> float:
    return k_strict(d + 1) - k_strict(d)


def relabel(d: int, unit: int) -> int:
    r = (unit * d) % N
    return N if r == 0 else r


def build_witness() -> dict[str, Any]:
    values = {d: round(k_strict(d), 12) for d in range(1, N + 1)}
    obs = {d: round(t_obs(d), 12) for d in range(1, N + 1)}
    signs = {d: (1 if t_obs(d) > 0 else -1 if t_obs(d) < 0 else 0) for d in range(1, N + 1)}
    monotone_nonincreasing = all(t_obs(d) <= 0 for d in range(1, N + 1))
    monotone_nondecreasing = all(t_obs(d) >= 0 for d in range(1, N + 1))

    unit_rows = []
    failures = 0
    for u in UNITS:
        for d in range(1, N + 1):
            lhs = t_obs(relabel(d, u))
            rhs = t_obs(d)
            invariant = math.isclose(lhs, rhs, rel_tol=0.0, abs_tol=1e-12)
            if not invariant:
                failures += 1
            unit_rows.append({
                "unit": u,
                "d": d,
                "u_dot_d": relabel(d, u),
                "T_K(u*d)": round(lhs, 12),
                "T_K(d)": round(rhs, 12),
                "unit_invariant": invariant,
            })

    arrow_rows = [
        {"test": "nonzero_signed_values", "passed": any(v != 0 for v in signs.values()), "evidence": signs},
        {"test": "single_global_arrow_sign", "passed": len(set(signs.values()) - {0}) == 1, "evidence": signs},
        {"test": "monotone_clock_on_1_to_12", "passed": monotone_nonincreasing or monotone_nondecreasing, "evidence": {"nonincreasing": monotone_nonincreasing, "nondecreasing": monotone_nondecreasing}},
        {"test": "Aut_Z12_unit_invariant_readout", "passed": failures == 0, "evidence": {"unit_rows": len(unit_rows), "failures": failures}},
    ]
    obligations = [
        {"obligation": "explicit_input_output_types", "satisfied": True, "detail": "input d in {1,...,12}; output real finite-difference scalar T_K(d)"},
        {"obligation": "observer_independent_formula", "satisfied": True, "detail": "formula uses only K_strict_gate values and finite difference"},
        {"obligation": "unit_action_compatible", "satisfied": failures == 0, "detail": f"{failures} of {len(unit_rows)} U(12) relabeling checks fail"},
        {"obligation": "directed_time_arrow_without_selector_import", "satisfied": len(set(signs.values()) - {0}) == 1, "detail": "requires one consistent arrow sign on the audited cycle"},
        {"obligation": "eom_hamiltonian_installation", "satisfied": False, "detail": "no unit-bearing action/EOM/Hamiltonian source is exported by this observable formula"},
    ]
    return {
        "candidate": "StrictKernelFiniteDifference_TimeObservableGenerator",
        "formula": "T_K(d)=K_strict_gate(d+1)-K_strict_gate(d), K=cos(omega*d+phi)/(1+beta*d^eta)",
        "parameters": {"omega": OMEGA, "phi": PHI, "beta": BETA, "eta": ETA, "cycle": N},
        "kernel_values": values,
        "observable_values": obs,
        "signs": signs,
        "unit_equivariance_rows": unit_rows,
        "unit_equivariance_failure_count": failures,
        "arrow_rows": arrow_rows,
        "proof_obligations": obligations,
        "accepted_as_time_observable_generator": all(o["satisfied"] for o in obligations),
    }


def build_payload(p3013_path: Any) -> dict[str, Any]:
    w = build_witness()
    return {
        "status": "P3014_TIME_READOUT_OBSERVABLE_GENERATOR_CANDIDATE_BOUNDED_NO_GO",
        "input_hashes": {"P3013": hashlib.sha256(p3013_path.read_bytes()).hexdigest() if p3013_path.exists() else None},
        "constructed_theoretical_objects": {
            "object": "TimeReadoutObservableGenerator_FiniteDifferenceObstructionMatrix",
            "targeted_readout_atom": "observer_independent_observable_generator_for_time_row",
            "witness": w,
        },
        "finite_certificate": {
            "explicit_io": True,
            "observer_independent_formula": True,
            "unit_equivariance_rows": len(w["unit_equivariance_rows"]),
            "unit_equivariance_failure_count": w["unit_equivariance_failure_count"],
            "arrow_passed_rows": sum(1 for r in w["arrow_rows"] if r["passed"]),
            "arrow_row_count": len(w["arrow_rows"]),
            "accepted_as_time_observable_generator": w["accepted_as_time_observable_generator"],
        },
        "decision": {
            "breakthrough": "A concrete observer-independent time-readout formula was constructed from the strict kernel finite difference, but the finite test blocks it as a ToE-grade time observable: it is not U(12)-unit compatible, does not install an EOM/Hamiltonian, and cannot be promoted to a directed physical time arrow without the missing selector/unit/action atoms.",
            "negative_export_flags": {k: False for k in ["time_observable_generator_exported", "time_arrow_exported", "strict_selector_source_exported", "unit_action_exported", "eom_hamiltonian_exported", "observed_physics_exported", "ltotal_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not replay finite-difference kernel-clock observables.  The next proof-grade move should attack one different atom: either construct a unit-action-compatible observable functor for a single readout row, or provide a genuine unit-bearing action/EOM source for the already explicit time observable; keep selector, bridge, role-transfer, L_total, and ToE closure blocked unless that atom passes a finite acceptance test.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = [
        "# P3014/S1964 time-readout observable-generator candidate obstruction", "",
        f"Status: `{payload['status']}`", "", "## Finite certificate",
        f"- explicit input/output types: `{c['explicit_io']}`",
        f"- observer-independent formula: `{c['observer_independent_formula']}`",
        f"- unit equivariance rows/failures: `{c['unit_equivariance_rows']}/{c['unit_equivariance_failure_count']}`",
        f"- arrow rows passed/total: `{c['arrow_passed_rows']}/{c['arrow_row_count']}`",
        f"- accepted as time observable generator: `{c['accepted_as_time_observable_generator']}`", "",
        "## Decision", payload["decision"]["breakthrough"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    read_json(P3013)
    payload = build_payload(P3013)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3014/S1964 time-readout observable-generator candidate obstruction", "## P3014/S1964 time-readout observable-generator candidate obstruction\n\n`P3014/S1964` attacks exactly one P3013 readout atom: an observer-independent observable generator for the time row.  It constructs the explicit finite-difference candidate `T_K(d)=K_strict_gate(d+1)-K_strict_gate(d)` on the twelve-label strict-kernel cycle.  The formula has explicit input/output types and is observer-independent, but the finite acceptance matrix is negative: U(12) relabeling/unit compatibility fails, no unit-bearing action/EOM/Hamiltonian source is installed, and no selector-free directed time-arrow export follows.  Therefore no time observable generator, time arrow, `L_total`, bridge/role transfer, observed-physics export, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3014/S1964 finite-difference time observable `L_total` guard", "## P3014/S1964 finite-difference time observable `L_total` guard\n\n`P3014/S1964` does not add an `L_total` term.  The candidate `T_K(d)=K_strict_gate(d+1)-K_strict_gate(d)` is a typed time-readout observable formula only; it fails U(12) unit compatibility and carries no unit-bearing variational action, EOM, Hamiltonian, stress-energy, selector source, bridge closure, or role-transfer theorem.\n")
    append_once(AGENTS, "Current finite-difference time observable guardrail (P3014/S1964, 2026-06-22)", "## Current finite-difference time observable guardrail (P3014/S1964, 2026-06-22)\n\n- P3014 attacks one P3013 atom: an observer-independent observable generator for the time row.\n- The constructed candidate `T_K(d)=K_strict_gate(d+1)-K_strict_gate(d)` has explicit input/output types and uses only strict-kernel values, but its finite acceptance matrix is bounded no-go: U(12) unit compatibility fails and no unit-bearing EOM/Hamiltonian or selector-free directed arrow is exported.\n- Do not replay finite-difference kernel-clock observables as time-arrow, observed-physics, `L_total`, bridge/role-transfer, or ToE closure evidence.\n- The next honest move should target exactly one different atom: either a unit-action-compatible observable functor for one readout row, or a genuine unit-bearing action/EOM source for an already typed observable.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
