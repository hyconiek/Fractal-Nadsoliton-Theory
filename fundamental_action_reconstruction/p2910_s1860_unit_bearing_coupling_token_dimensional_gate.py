#!/usr/bin/env python3
"""P2910/S1860: unit-bearing coupling-token dimensional gate.

P2909 says the next admissible move must either supply a new strict construction
for J_{0,+} or pivot to a different typed object with a bounded acceptance test.
P2910 takes the pivot route: it constructs a new typed object, a dimensioned
coupling token Gamma_9_5, and checks the purely dimensional part of the missing
unit-bearing L_total coupling problem.

This is deliberately not an Xi/J/defect-placement replay.  The token is a typed
unit carrier for a generic local 9/5 density slot.  The finite vector algebra
shows what dimension Gamma_9_5 would need to have for a symbolic q_9_5 slot to
enter an action-density term.  The gate passes dimensional consistency but fails
strict provenance, localization/pullback, variational-chain, and nonproxy L_total
exports.
"""
from __future__ import annotations

import hashlib
import json
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN

P2909 = GEN / "p2909_s1859_post_joint_source_state_map_no_new_live_frontier.json"
OUT = GEN / "p2910_s1860_unit_bearing_coupling_token_dimensional_gate.json"
MD = GEN / "p2910_s1860_unit_bearing_coupling_token_dimensional_gate.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

BASES = ("Action", "Length", "Time")
ZERO = (0, 0, 0)
ACTION_DENSITY = (1, -1, 0)  # action per one strict length/site unit
DIMENSIONLESS = ZERO
EDGE_DELTA = (0, -1, 0)      # localized carrier contributes one inverse site-length
Q_9_5 = DIMENSIONLESS


def add(a: tuple[int, int, int], b: tuple[int, int, int]) -> tuple[int, int, int]:
    return tuple(x + y for x, y in zip(a, b))  # type: ignore[return-value]


def sub(a: tuple[int, int, int], b: tuple[int, int, int]) -> tuple[int, int, int]:
    return tuple(x - y for x, y in zip(a, b))  # type: ignore[return-value]


def fmt(dim: tuple[int, int, int]) -> dict[str, int]:
    return dict(zip(BASES, dim))


def required_gamma_dimension() -> tuple[int, int, int]:
    # Gamma + delta_edge + q_9_5 = action_density.
    return sub(ACTION_DENSITY, add(EDGE_DELTA, Q_9_5))


def build_payload(p2909: dict[str, Any]) -> dict[str, Any]:
    gamma_dim = required_gamma_dimension()
    reconstructed = add(gamma_dim, add(EDGE_DELTA, Q_9_5))
    candidate = {
        "name": "Gamma_9_5",
        "typed_object_class": "unit-bearing coupling token candidate, not a strict source theorem",
        "dimension_vector_basis": BASES,
        "required_dimension": fmt(gamma_dim),
        "formal_density_slot": "Gamma_9_5 * delta_edge * q_9_5",
        "reconstructed_density_dimension": fmt(reconstructed),
    }
    obligations = [
        {"obligation": "dimension algebra", "passed": reconstructed == ACTION_DENSITY, "exported_strictly": False},
        {"obligation": "strict source/provenance for Gamma_9_5", "passed": False, "exported_strictly": False},
        {"obligation": "localization/pullback map to continuum or strict site measure", "passed": False, "exported_strictly": False},
        {"obligation": "variational chain rule into nonproxy L_total", "passed": False, "exported_strictly": False},
        {"obligation": "compatibility with J provenance", "passed": False, "exported_strictly": False},
    ]
    return {
        "status": "P2910_UNIT_BEARING_COUPLING_TOKEN_DIMENSIONAL_GATE_READINESS_NO_EXPORT",
        "input_hashes": {"P2909": hashlib.sha256(P2909.read_bytes()).hexdigest() if P2909.exists() else None},
        "constructed_theoretical_objects": {
            "dimension_basis": BASES,
            "dimension_vectors": {
                "target_action_density": fmt(ACTION_DENSITY),
                "localized_edge_delta": fmt(EDGE_DELTA),
                "q_9_5_slot": fmt(Q_9_5),
                "required_gamma_9_5": fmt(gamma_dim),
            },
            "coupling_token_candidate": candidate,
            "acceptance_obligation_rows": obligations,
        },
        "acceptance_matrix": {
            "p2909_rechecked_no_live_frontier": p2909.get("acceptance_matrix", {}).get("no_new_live_frontier_certificate") is True,
            "new_typed_object_constructed": True,
            "outside_xi_j_defect_placement_lane": True,
            "dimension_algebra_passed": reconstructed == ACTION_DENSITY,
            "required_gamma_dimension": fmt(gamma_dim),
            "strict_gamma_source_exported": False,
            "localization_pullback_exported": False,
            "variational_chain_rule_exported": False,
            "compatible_with_j_provenance": False,
            "accepted_as_unit_bearing_ltotal_coupling": False,
        },
        "decision": {
            "positive_witnesses": {
                "typed_coupling_token_constructed": True,
                "required_unit_dimension_computed": True,
                "dimension_consistency_witnessed": reconstructed == ACTION_DENSITY,
            },
            "negative_export_flags": {
                "strict_gamma_9_5_source_exported": False,
                "unit_bearing_nonproxy_ltotal_exported": False,
                "localization_pullback_exported": False,
                "variational_chain_rule_exported": False,
                "eom_closure_exported": False,
                "hamiltonian_closure_exported": False,
                "role_transfer_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2910 pivots outside the Xi/J/defect-placement replay by constructing a new typed coupling-token candidate Gamma_9_5.  Finite dimension algebra computes the required unit vector and verifies that Gamma_9_5 * delta_edge * q_9_5 can match an action-density dimension.  This is only dimensional readiness: no strict source for Gamma_9_5, no localization/pullback, no variational chain rule, no J provenance compatibility, and no nonproxy L_total export are present.",
            "next_honest_step": "The next proof-grade move, if staying in this new typed-object lane, must provide exactly one strict source/provenance theorem for Gamma_9_5 or a localization/pullback map with a variational chain rule.  Do not promote dimension matching alone to L_total/EOM/Hamiltonian/ToE closure.  If no such theorem is supplied, preserve no-new-live-frontier or pivot to another genuinely new typed object.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    acc = payload["acceptance_matrix"]
    lines = [
        "# P2910/S1860 unit-bearing coupling-token dimensional gate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Dimensional gate",
        f"- new typed object constructed: `{acc['new_typed_object_constructed']}`",
        f"- outside Xi/J/defect-placement lane: `{acc['outside_xi_j_defect_placement_lane']}`",
        f"- dimension algebra passed: `{acc['dimension_algebra_passed']}`",
        f"- required Gamma_9_5 dimension: `{acc['required_gamma_dimension']}`",
        f"- accepted as unit-bearing L_total coupling: `{acc['accepted_as_unit_bearing_ltotal_coupling']}`",
        "",
        "## Boundary",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload(read_json(P2909))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2910/S1860 unit-bearing coupling-token dimensional gate", "## P2910/S1860 unit-bearing coupling-token dimensional gate\n\n`P2910/S1860` pivots outside Xi/J/defect-placement by constructing a new typed coupling-token candidate `Gamma_9_5`.  Finite dimension algebra computes the required unit vector so that `Gamma_9_5 * delta_edge * q_9_5` has action-density dimension; this passes only the dimensional readiness gate.  No strict source/provenance for `Gamma_9_5`, no localization/pullback map, no variational chain rule, no compatibility with `J` provenance, and no nonproxy `L_total`, EOM, Hamiltonian, role transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2910/S1860 coupling-token dimensional `L_total` guard", "## P2910/S1860 coupling-token dimensional `L_total` guard\n\n`P2910/S1860` identifies the dimension a candidate `Gamma_9_5` would need for a local `9/5` density slot to have action-density units.  Dimension matching is necessary but not sufficient: without strict source/provenance, localization/pullback, and a variational chain rule, the token cannot be inserted into nonproxy `L_total`, EOM, Hamiltonian, role transfer, or ToE closure.\n")
    append_once(AGENTS, "Current unit-bearing coupling-token dimensional gate guardrail (P2910/S1860, 2026-06-19)", "## Current unit-bearing coupling-token dimensional gate guardrail (P2910/S1860, 2026-06-19)\n\n- P2910 introduces a new typed object outside Xi/J/defect-placement: the coupling-token candidate `Gamma_9_5`.\n- Finite dimension algebra computes the required action-density unit vector and passes dimensional consistency, but this is readiness only.\n- Do not promote `Gamma_9_5`, dimension matching, symbolic `q_9_5`, `rho_9/5`, or `U_9_5` to strict sourcehood, nonproxy `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE closure without strict source/provenance, localization/pullback, and variational chain-rule theorems.\n- A next admissible move in this lane must prove exactly one of those missing theorems; otherwise preserve no-new-live-frontier or pivot to a different genuinely new typed object.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
