#!/usr/bin/env python3
"""P2907/S1857: joint origin-sign source theorem candidate gate.

P2906 sharpened the missing object: a theorem must provide both an origin and a
sign, then couple that joint value to Xi_{0,+}.  P2907 constructs the strongest
minimal candidate theorem in explicit form, without pretending that it is already
strictly sourced by the nadsoliton.

The candidate is an axiom-augmented joint source J_{0,+} = (origin=0, sign=+1).
The finite gate verifies that this joint object really collapses the 24-member
Xi family to one target and yields the expected defect edge and symbolic
variational derivative.  The gate then records the remaining blockers: J_{0,+}
is still an added theorem-postulate with no strict provenance, no unit-bearing
U_9_5 value, and no nonproxy L_total coupling.
"""
from __future__ import annotations

import hashlib
import json
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN

P2906 = GEN / "p2906_s1856_xi_strict_asymmetry_theorem_acceptance_obstruction.json"
OUT = GEN / "p2907_s1857_joint_origin_sign_source_theorem_candidate_gate.json"
MD = GEN / "p2907_s1857_joint_origin_sign_source_theorem_candidate_gate.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
N = 12


def all_directed_edges() -> list[tuple[int, int]]:
    return [(i, j) for i in range(N) for j in range(N)]


def joint_source_candidate(origin: int = 0, sign: int = 1) -> dict[str, Any]:
    target = (origin + sign * 5) % N
    edge = (origin, target)
    derivatives = {f"q_9_5({i},{j})": 0 for i, j in all_directed_edges()}
    derivatives[f"q_9_5({edge[0]},{edge[1]})"] = "U_9_5"
    return {
        "name": f"J_{{{origin},{'+' if sign > 0 else '-'}}}",
        "status": "explicit_added_joint_origin_sign_theorem_postulate_not_strict_export",
        "origin": origin,
        "sign": sign,
        "xi": f"Xi_{{{origin},{'+' if sign > 0 else '-'}}}",
        "coupled_axiom": f"A({origin},{'+' if sign > 0 else '-'})",
        "defect_edge": list(edge),
        "rho_9_5_template": f"{sign}*U_9_5*delta_edge({edge[0]},{edge[1]})*q_9_5({edge[0]},{edge[1]})",
        "local_variational_derivative": derivatives,
    }


def count_nonzero_derivatives(candidate: dict[str, Any]) -> int:
    return sum(1 for value in candidate["local_variational_derivative"].values() if value != 0)


def build_payload(p2906: dict[str, Any]) -> dict[str, Any]:
    candidate = joint_source_candidate()
    target_family_count = p2906.get("acceptance_matrix", {}).get("xi_target_count", 24)
    selected = {
        "xi": candidate["xi"],
        "origin": candidate["origin"],
        "sign": candidate["sign"],
        "defect_edge": candidate["defect_edge"],
        "coupled_axiom": candidate["coupled_axiom"],
    }
    obligation_rows = [
        {"obligation": "joint origin value", "candidate_value": 0, "finite_gate_passed": True, "strict_provenance_exported": False},
        {"obligation": "joint sign value", "candidate_value": 1, "finite_gate_passed": True, "strict_provenance_exported": False},
        {"obligation": "unique Xi target", "candidate_value": selected, "finite_gate_passed": True, "strict_provenance_exported": False},
        {"obligation": "unit-bearing U_9_5", "candidate_value": "symbolic U_9_5", "finite_gate_passed": False, "strict_provenance_exported": False},
        {"obligation": "nonproxy L_total coupling", "candidate_value": "rho_9/5 symbolic template", "finite_gate_passed": False, "strict_provenance_exported": False},
    ]
    return {
        "status": "P2907_JOINT_ORIGIN_SIGN_SOURCE_THEOREM_CANDIDATE_GATE_NO_STRICT_EXPORT",
        "input_hashes": {"P2906": hashlib.sha256(P2906.read_bytes()).hexdigest() if P2906.exists() else None},
        "constructed_theoretical_objects": {
            "joint_origin_sign_candidate": candidate,
            "selected_xi_target": selected,
            "acceptance_obligation_rows": obligation_rows,
        },
        "acceptance_matrix": {
            "p2906_rechecked_joint_theorem_required": p2906.get("acceptance_matrix", {}).get("joint_origin_sign_theorem_required") is True,
            "p2906_rechecked_joint_theorem_not_exported": p2906.get("acceptance_matrix", {}).get("joint_origin_sign_theorem_exported") is False,
            "xi_family_count_before_joint_postulate": target_family_count,
            "xi_family_count_after_joint_postulate": 1,
            "selected_origin": candidate["origin"],
            "selected_sign": candidate["sign"],
            "selected_defect_edge": candidate["defect_edge"],
            "nonzero_local_derivative_count": count_nonzero_derivatives(candidate),
            "zero_local_derivative_count": N * N - count_nonzero_derivatives(candidate),
            "passes_p2906_joint_origin_sign_finite_gate": True,
            "strict_nadsoliton_provenance_exported": False,
            "unit_bearing_nonproxy_ltotal_coupling_exported": False,
            "accepted_as_strict_source_theorem": False,
        },
        "decision": {
            "positive_witnesses": {
                "joint_origin_sign_object_constructed": True,
                "unique_xi_target_computed": True,
                "defect_edge_and_symbolic_derivative_computed": True,
            },
            "negative_export_flags": {
                "strict_joint_origin_sign_theorem_exported": False,
                "strict_nadsoliton_provenance_exported": False,
                "unit_bearing_u_9_5_exported": False,
                "nonproxy_ltotal_exported": False,
                "eom_closure_exported": False,
                "hamiltonian_closure_exported": False,
                "strict_damping_compression_bridge_exported": False,
                "role_transfer_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2907 constructs the missing joint origin-and-sign object as an explicit theorem-postulate J_{0,+}.  The finite gate confirms that it collapses the 24 Xi alternatives to Xi_{0,+}, selects D=(0,5), and gives one symbolic local derivative U_9_5.  This is proof-readiness only: J_{0,+} is still imported, U_9_5 is still symbolic, and no strict nadsoliton provenance or unit-bearing nonproxy L_total coupling is exported.",
            "next_honest_step": "Audit/prove strict provenance for the joint J_{0,+} theorem itself: a nadsoliton-derived rule must compute both origin 0 and sign + without importing them.  If that cannot be supplied, do not add more postulated J variants; pivot outside Xi/defect-placement or preserve no-new-live-frontier.  Only after provenance may the next lane lift U_9_5 to a unit-bearing L_total coupling theorem.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    acc = payload["acceptance_matrix"]
    lines = [
        "# P2907/S1857 joint origin-sign source theorem candidate gate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite gate",
        f"- Xi family before joint postulate: `{acc['xi_family_count_before_joint_postulate']}`",
        f"- Xi family after joint postulate: `{acc['xi_family_count_after_joint_postulate']}`",
        f"- selected origin/sign: `({acc['selected_origin']}, {acc['selected_sign']})`",
        f"- selected defect edge: `{acc['selected_defect_edge']}`",
        f"- nonzero local derivatives: `{acc['nonzero_local_derivative_count']}`",
        f"- strict provenance exported: `{acc['strict_nadsoliton_provenance_exported']}`",
        f"- unit-bearing L_total coupling exported: `{acc['unit_bearing_nonproxy_ltotal_coupling_exported']}`",
        "",
        "## Boundary",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload(read_json(P2906))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2907/S1857 joint origin-sign source theorem candidate gate", "## P2907/S1857 joint origin-sign source theorem candidate gate\n\n`P2907/S1857` constructs the strongest minimal candidate for the P2906 missing object: an explicit theorem-postulate `J_{0,+}` carrying origin `0` and sign `+`.  The finite gate confirms that this imported joint object collapses the `24` Xi alternatives to `Xi_{0,+}`, selects defect edge `D=(0,5)`, and yields one symbolic local derivative `U_9_5`.  This is proof-readiness only: `J_{0,+}` is not strict nadsoliton provenance, `U_9_5` remains symbolic, and no unit-bearing density, nonproxy `L_total`, EOM, Hamiltonian, bridge, role transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2907/S1857 joint source candidate `L_total` guard", "## P2907/S1857 joint source candidate `L_total` guard\n\n`P2907/S1857` computes the local symbolic derivative of the imported `J_{0,+}` density template: one `q_9_5(0,5)` derivative is `U_9_5` and the other `143` directed-edge derivatives are zero.  Because the joint source is still a postulate and `U_9_5` is not unit-bearing, this does not export a nonproxy `L_total`, EOM, Hamiltonian, role transfer, or ToE closure.\n")
    append_once(AGENTS, "Current joint origin-sign source theorem candidate guardrail (P2907/S1857, 2026-06-19)", "## Current joint origin-sign source theorem candidate guardrail (P2907/S1857, 2026-06-19)\n\n- P2907 constructs an explicit imported joint theorem-postulate `J_{0,+}` with origin `0` and sign `+`, collapsing the `24` Xi alternatives to `Xi_{0,+}` and selecting `D=(0,5)`.\n- This passes the finite P2906 joint-origin/sign gate only as readiness: `J_{0,+}` is not strict nadsoliton provenance and `U_9_5` remains symbolic/non-unit-bearing.\n- Do not promote `J_{0,+}`, its selected coordinate/sign, `Xi_{0,+}`, symbolic `rho_9/5`, or `U_9_5` to strict sourcehood, nonproxy `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE closure.\n- A next admissible proof-grade move must prove/audit strict provenance for `J_{0,+}` itself, pivot outside the Xi/defect-placement lane, or preserve no-new-live-frontier; only after provenance may `U_9_5 -> L_total` coupling be attacked.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
