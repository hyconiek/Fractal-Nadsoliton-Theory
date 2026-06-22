#!/usr/bin/env python3
"""P3015/S1965: unit-action-compatible time observable functor obstruction.

This step follows P3014 without replaying the finite-difference kernel-clock lane.
It attacks the different atom named by P3014: construct a U(12)-unit-action
compatible observable functor for one readout row (`time`).

The constructed object is the strict-kernel orbit-quotient observable functor:
labels d in Z/12Z are mapped to U(12)-orbits, and a real observable is assigned
by averaging K_strict_gate over each orbit.  This succeeds as a unit-compatible,
observer-independent scalar readout.  The finite obstruction is that the quotient
does not carry a well-defined successor/clock map: members of the same U(12)
orbit can have successors in different orbits.  Thus unit compatibility is gained
only by forgetting directed clock structure, so no time-arrow, EOM/Hamiltonian,
L_total, or ToE closure is exported.
"""
from __future__ import annotations

import hashlib, json, math
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3014_s1964_time_readout_observable_generator_candidate_obstruction import OUT as P3014

OUT = GEN / "p3015_s1965_unit_action_compatible_time_observable_functor_obstruction.json"
MD = GEN / "p3015_s1965_unit_action_compatible_time_observable_functor_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

N = 12
UNITS = [1, 5, 7, 11]
OMEGA = 0.18575
PHI = 0.16250
BETA = 1.0
ETA = 1.8


def normalize_label(x: int) -> int:
    r = x % N
    return N if r == 0 else r


def act(unit: int, d: int) -> int:
    return normalize_label(unit * d)


def k_strict(d: int) -> float:
    return math.cos(OMEGA * d + PHI) / (1.0 + BETA * (float(d) ** ETA))


def orbit_of(d: int) -> tuple[int, ...]:
    return tuple(sorted({act(u, d) for u in UNITS}))


def orbit_name(orbit: tuple[int, ...]) -> str:
    return "{" + ",".join(str(x) for x in orbit) + "}"


def build_orbit_functor() -> dict[str, Any]:
    orbits = sorted({orbit_of(d) for d in range(1, N + 1)}, key=lambda o: (len(o), o))
    label_to_orbit = {d: orbit_of(d) for d in range(1, N + 1)}
    orbit_values = {}
    for orbit in orbits:
        vals = [k_strict(d) for d in orbit]
        orbit_values[orbit_name(orbit)] = {
            "members": list(orbit),
            "mean_K": round(sum(vals) / len(vals), 12),
            "min_K": round(min(vals), 12),
            "max_K": round(max(vals), 12),
        }

    unit_rows = []
    failures = 0
    for u in UNITS:
        for d in range(1, N + 1):
            invariant = label_to_orbit[act(u, d)] == label_to_orbit[d]
            failures += 0 if invariant else 1
            unit_rows.append({
                "unit": u,
                "d": d,
                "u_dot_d": act(u, d),
                "orbit_d": orbit_name(label_to_orbit[d]),
                "orbit_u_dot_d": orbit_name(label_to_orbit[act(u, d)]),
                "orbit_functor_invariant": invariant,
            })

    successor_rows = []
    bad_successor_orbits = 0
    for orbit in orbits:
        successor_orbits = sorted({orbit_of(normalize_label(d + 1)) for d in orbit}, key=lambda o: (len(o), o))
        well_defined = len(successor_orbits) == 1
        bad_successor_orbits += 0 if well_defined else 1
        successor_rows.append({
            "orbit": orbit_name(orbit),
            "member_successors": {str(d): normalize_label(d + 1) for d in orbit},
            "successor_orbits": [orbit_name(o) for o in successor_orbits],
            "successor_well_defined_on_quotient": well_defined,
        })

    obligations = [
        {"obligation": "explicit_input_output_types", "satisfied": True, "detail": "input U(12)-orbit [d]; output real scalar mean_{x in [d]} K_strict_gate(x)"},
        {"obligation": "observer_independent_formula", "satisfied": True, "detail": "orbit mean uses strict-kernel values and finite group action only"},
        {"obligation": "unit_action_compatible", "satisfied": failures == 0, "detail": f"{failures} of {len(unit_rows)} orbit-invariance checks fail"},
        {"obligation": "well_defined_clock_successor", "satisfied": bad_successor_orbits == 0, "detail": f"{bad_successor_orbits} of {len(orbits)} quotient orbits have nonunique successor orbit"},
        {"obligation": "directed_time_arrow_without_selector_import", "satisfied": False, "detail": "orbit quotient is invariant but orientation-blind; it does not select a directed time branch"},
        {"obligation": "eom_hamiltonian_installation", "satisfied": False, "detail": "no unit-bearing variational action/EOM/Hamiltonian source is supplied by orbit averaging"},
    ]
    return {
        "candidate": "StrictKernelOrbitQuotient_TimeObservableFunctor",
        "formula": "O_K([d]) = |[d]|^{-1} * sum_{x in U(12).d} K_strict_gate(x)",
        "parameters": {"omega": OMEGA, "phi": PHI, "beta": BETA, "eta": ETA, "cycle": N, "units": UNITS},
        "orbits": [orbit_name(o) for o in orbits],
        "orbit_values": orbit_values,
        "unit_compatibility_rows": unit_rows,
        "unit_compatibility_failure_count": failures,
        "successor_rows": successor_rows,
        "bad_successor_orbit_count": bad_successor_orbits,
        "proof_obligations": obligations,
        "accepted_as_time_observable_generator": all(o["satisfied"] for o in obligations),
    }


def build_payload(p3014_path: Any) -> dict[str, Any]:
    functor = build_orbit_functor()
    return {
        "status": "P3015_UNIT_ACTION_COMPATIBLE_TIME_OBSERVABLE_FUNCTOR_OBSTRUCTION_NO_CLOSURE",
        "input_hashes": {"P3014": hashlib.sha256(p3014_path.read_bytes()).hexdigest() if p3014_path.exists() else None},
        "constructed_theoretical_objects": {
            "object": "UnitActionCompatibleTimeObservableFunctor_QuotientSuccessorObstructionMatrix",
            "targeted_atom": "unit_action_compatible_observable_functor_for_time_row",
            "functor": functor,
        },
        "finite_certificate": {
            "orbit_count": len(functor["orbits"]),
            "unit_compatibility_rows": len(functor["unit_compatibility_rows"]),
            "unit_compatibility_failure_count": functor["unit_compatibility_failure_count"],
            "successor_row_count": len(functor["successor_rows"]),
            "bad_successor_orbit_count": functor["bad_successor_orbit_count"],
            "accepted_as_time_observable_generator": functor["accepted_as_time_observable_generator"],
        },
        "decision": {
            "breakthrough": "A genuine U(12)-unit-compatible, observer-independent strict-kernel observable functor was constructed by orbit-averaging.  The bounded obstruction is that the orbit quotient destroys the clock successor: same-orbit labels can advance to different successor orbits, so the functor is not a directed time observable generator and cannot install EOM/Hamiltonian or ToE closure.",
            "negative_export_flags": {k: False for k in ["time_observable_generator_exported", "time_arrow_exported", "strict_selector_source_exported", "unit_bearing_eom_hamiltonian_exported", "ltotal_exported", "observed_physics_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not replay orbit-average quotient observables as time-arrow closure.  The next proof-grade move should attack exactly one missing successor/evolution atom: construct a strict unit-compatible directed successor/semigroup on the quotient, or pivot to a unit-bearing action/EOM source for one typed observable; keep selector, bridge, role-transfer, L_total, observed-physics, and ToE closure blocked until that single atom passes finite acceptance.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = [
        "# P3015/S1965 unit-action-compatible time observable functor obstruction", "",
        f"Status: `{payload['status']}`", "", "## Finite certificate",
        f"- orbit count: `{c['orbit_count']}`",
        f"- unit compatibility rows/failures: `{c['unit_compatibility_rows']}/{c['unit_compatibility_failure_count']}`",
        f"- successor rows / bad successor orbits: `{c['successor_row_count']}/{c['bad_successor_orbit_count']}`",
        f"- accepted as time observable generator: `{c['accepted_as_time_observable_generator']}`", "",
        "## Decision", payload["decision"]["breakthrough"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    read_json(P3014)
    payload = build_payload(P3014)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3015/S1965 unit-action-compatible time observable functor obstruction", "## P3015/S1965 unit-action-compatible time observable functor obstruction\n\n`P3015/S1965` attacks the P3014 recommendation without replaying finite-difference kernel-clock observables.  It constructs a U(12)-unit-compatible observer-independent observable functor for the time row by passing to U(12)-orbits and assigning `O_K([d]) = |[d]|^{-1} sum_{x in [d]} K_strict_gate(x)`.  The finite matrix verifies zero unit-compatibility failures, but the quotient has no well-defined clock successor: same-orbit labels can advance to different successor orbits.  Thus unit compatibility is achieved only by losing directed clock/evolution structure.  No time-observable-generator export, time arrow, EOM/Hamiltonian, `L_total`, observed-physics, bridge/role transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3015/S1965 orbit-quotient time observable `L_total` guard", "## P3015/S1965 orbit-quotient time observable `L_total` guard\n\n`P3015/S1965` adds no `L_total` term.  The orbit-average functor is unit-action-compatible as a scalar quotient readout, but it lacks a well-defined clock successor on the quotient and supplies no unit-bearing action, variational density, EOM, Hamiltonian, selector source, bridge closure, or role-transfer theorem.\n")
    append_once(AGENTS, "Current orbit-quotient time observable functor guardrail (P3015/S1965, 2026-06-22)", "## Current orbit-quotient time observable functor guardrail (P3015/S1965, 2026-06-22)\n\n- P3015 constructs a genuine U(12)-unit-compatible observer-independent scalar observable functor for the time row by orbit-averaging strict-kernel values.\n- The finite obstruction is that the U(12) quotient has no well-defined clock successor: same-orbit labels can advance to different successor orbits, so directed time/evolution structure is lost.\n- Do not promote orbit-average quotient observables to time-arrow, EOM/Hamiltonian, `L_total`, observed-physics, bridge/role-transfer, selector closure, or ToE closure.\n- The next honest move should target exactly one missing successor/evolution atom, such as a strict unit-compatible directed successor/semigroup on the quotient, or pivot to one unit-bearing action/EOM source for a typed observable.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
