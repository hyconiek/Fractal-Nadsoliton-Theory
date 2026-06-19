#!/usr/bin/env python3
"""P2925/S1875: damping delta-source linear-system frontier certificate.

P2924 showed that the strict target slope delta=4/5 is compatible with the
prime-log character shape but not selected by it.  P2925 makes that obstruction
more algebraic: it builds the finite linear system with variables
(y_2,...,y_11, delta).  Existing character-shape equations have rank 10 and
nullity 1, leaving exactly the global slope line free.  The single row that
would select delta=4/5 is therefore a new anchor row, not a consequence of the
current finite data.
"""
from __future__ import annotations

import hashlib
import json
import math
from fractions import Fraction
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN

P2923 = GEN / "p2923_s1873_prime_log_proportionality_source_matrix.json"
P2924 = GEN / "p2924_s1874_slope_prime_anchor_source_obstruction_matrix.json"
OUT = GEN / "p2925_s1875_damping_delta_source_linear_system_frontier_certificate.json"
MD = GEN / "p2925_s1875_damping_delta_source_linear_system_frontier_certificate.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

NODES = list(range(2, 12))
VARIABLES = [f"y_{d}" for d in NODES] + ["delta"]
STRICT_DELTA = Fraction(4, 5)
STRICT_ETA = Fraction(9, 5)


def float_rank(matrix: list[list[float]], tol: float = 1e-10) -> int:
    mat = [row[:] for row in matrix]
    rank = 0
    col_count = len(mat[0]) if mat else 0
    for col in range(col_count):
        pivot = next((r for r in range(rank, len(mat)) if abs(mat[r][col]) > tol), None)
        if pivot is None:
            continue
        mat[rank], mat[pivot] = mat[pivot], mat[rank]
        pivot_value = mat[rank][col]
        mat[rank] = [value / pivot_value for value in mat[rank]]
        for r, row in enumerate(mat):
            if r != rank and abs(row[col]) > tol:
                factor = row[col]
                mat[r] = [x - factor * y for x, y in zip(row, mat[rank])]
        rank += 1
    return rank


def shape_rows() -> list[dict[str, Any]]:
    rows = []
    for d in NODES:
        coeffs = [0.0 for _ in VARIABLES]
        coeffs[VARIABLES.index(f"y_{d}")] = 1.0
        coeffs[VARIABLES.index("delta")] = -math.log(d)
        rows.append({
            "row": f"shape_y_{d}_minus_delta_log_{d}",
            "equation": f"y_{d} - delta*log({d}) = 0",
            "coefficients_over_variables_y2_to_y11_delta": coeffs,
            "rhs": 0.0,
            "sources_delta_value": False,
        })
    return rows


def target_anchor_row() -> dict[str, Any]:
    coeffs = [0.0 for _ in VARIABLES]
    coeffs[VARIABLES.index("delta")] = 1.0
    return {
        "row": "external_target_anchor_delta_4_5",
        "equation": "delta = 4/5",
        "coefficients_over_variables_y2_to_y11_delta": coeffs,
        "rhs": float(STRICT_DELTA),
        "sources_delta_value": False,
        "strict_source_status": "missing_new_anchor_source_law",
    }


def source_packet_rows() -> list[dict[str, Any]]:
    return [
        {
            "packet_row": "P2601_unital_multiplicative_source",
            "positive_content": "sources y_1=0 / unital multiplicative identity-action boundary",
            "missing_for_beta_eta": "does not provide prime-log atom values L_p or delta=4/5 anchor",
            "unlocks_strict_beta_eta": False,
        },
        {
            "packet_row": "P2923_prime_log_carrier_readiness",
            "positive_content": "finite factor-exponent carrier and product additivity over nodes 1..11",
            "missing_for_beta_eta": "does not source analytic L_p values or slope multiplier",
            "unlocks_strict_beta_eta": False,
        },
        {
            "packet_row": "P2924_slope_anchor_compatibility",
            "positive_content": "delta=4/5 is compatible with y_d=delta*log(d)",
            "missing_for_beta_eta": "same equations admit non-target deltas, so the target is not selected",
            "unlocks_strict_beta_eta": False,
        },
    ]


def candidate_next_objects() -> list[dict[str, Any]]:
    return [
        {
            "candidate": "Strict_Damping_Delta_Source_Law",
            "required_formula_shape": "strict_nadsoliton_data -> delta=4/5",
            "must_add_independent_rank": True,
            "currently_exported": False,
            "accepted_now": False,
            "reason": "the exact missing row; no current artifact supplies it",
        },
        {
            "candidate": "Strict_Prime_Log_Value_Source_Law",
            "required_formula_shape": "strict_nadsoliton_data -> (L_2,L_3,L_5,L_7,L_11)",
            "must_add_independent_rank": True,
            "currently_exported": False,
            "accepted_now": False,
            "reason": "the other P2923 residual atom; still absent after P2924",
        },
        {
            "candidate": "Damping_Beta_Eta_Source_Packet",
            "required_formula_shape": "prime-log values + delta source + coupling to beta/eta",
            "must_add_independent_rank": True,
            "currently_exported": False,
            "accepted_now": False,
            "reason": "a complete packet would be admissible, but current rows only define readiness/obstruction",
        },
    ]


def build_payload(p2923: dict[str, Any], p2924: dict[str, Any]) -> dict[str, Any]:
    rows = shape_rows()
    matrix = [row["coefficients_over_variables_y2_to_y11_delta"] for row in rows]
    rank_without_anchor = float_rank(matrix)
    augmented_matrix = matrix + [target_anchor_row()["coefficients_over_variables_y2_to_y11_delta"]]
    rank_with_external_anchor = float_rank(augmented_matrix)
    nullity_without_anchor = len(VARIABLES) - rank_without_anchor
    nullity_with_external_anchor = len(VARIABLES) - rank_with_external_anchor
    candidates = candidate_next_objects()
    return {
        "status": "P2925_DAMPING_DELTA_SOURCE_LINEAR_SYSTEM_FRONTIER_CERTIFICATE_NO_NEW_UNLOCK",
        "input_hashes": {
            "P2923": hashlib.sha256(P2923.read_bytes()).hexdigest() if P2923.exists() else None,
            "P2924": hashlib.sha256(P2924.read_bytes()).hexdigest() if P2924.exists() else None,
        },
        "constructed_theoretical_objects": {
            "linear_system_name": "Damping_Delta_Source_Linear_System",
            "variables": VARIABLES,
            "shape_rows": rows,
            "missing_anchor_row": target_anchor_row(),
            "minimal_unlock_packet_name": "Strict_Damping_Beta_Eta_Source_Packet",
            "minimal_unlock_packet_obligations": [
                "strict source of prime-log atom values L_p or an equivalent intrinsic log scale",
                "strict source law fixing delta=4/5 / eta=9/5 rather than importing a fitted anchor",
                "coupling theorem from the sourced log character and slope anchor to damping beta/eta",
                "explicit nonpromotion boundary for L_total, bridge, role transfer, and ToE until all obligations pass",
            ],
            "source_packet_rows": source_packet_rows(),
            "candidate_next_objects": candidates,
        },
        "linear_algebra_certificate": {
            "variable_count": len(VARIABLES),
            "shape_equation_count": len(rows),
            "rank_without_target_anchor": rank_without_anchor,
            "nullity_without_target_anchor": nullity_without_anchor,
            "rank_with_external_target_anchor": rank_with_external_anchor,
            "nullity_with_external_target_anchor": nullity_with_external_anchor,
            "target_anchor_adds_independent_row": rank_with_external_anchor == rank_without_anchor + 1,
            "target_anchor_is_sourced_by_current_rows": False,
        },
        "acceptance_matrix": {
            "p2923_prime_log_readiness_inherited": p2923.get("acceptance_matrix", {}).get("formal_log_character_readiness_exported") is True,
            "p2924_no_anchor_inherited": p2924.get("acceptance_matrix", {}).get("accepted_anchor_source_count") == 0,
            "source_packet_row_count": len(source_packet_rows()),
            "candidate_next_object_count": len(candidates),
            "accepted_current_unlock_object_count": sum(1 for row in candidates if row["accepted_now"]),
            "strict_prime_log_value_source_exported": False,
            "strict_delta_source_law_exported": False,
            "strict_damping_beta_eta_source_exported": False,
            "no_new_live_frontier_certificate_exported": True,
        },
        "decision": {
            "positive_witnesses": {
                "linear_free_slope_line_computed": nullity_without_anchor == 1,
                "external_target_anchor_would_close_nullity": nullity_with_external_anchor == 0,
                "minimal_unlock_packet_formulated": True,
            },
            "negative_export_flags": {
                "strict_prime_log_value_source_exported": False,
                "strict_delta_source_law_exported": False,
                "strict_damping_beta_eta_source_exported": False,
                "nonproxy_ltotal_exported": False,
                "eom_hamiltonian_exported": False,
                "bridge_closure_exported": False,
                "role_transfer_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2925 turns P2924 into a finite linear-algebra certificate.  The current character-shape rows have rank 10 in 11 variables (y_2..y_11, delta), so they leave exactly one free slope line.  Adding delta=4/5 would add the missing independent row and close the nullity, but that row is not sourced by current strict artifacts.",
            "next_honest_step": "Do not replay P2601/P2923/P2924 readiness.  The next admissible proof-grade move must introduce one genuinely new source object: either a strict prime-log value source law, a strict delta=4/5 source law, or a combined Strict_Damping_Beta_Eta_Source_Packet that passes the listed obligations.  If no such object is supplied, preserve this no-new-live-frontier certificate.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    la = payload["linear_algebra_certificate"]
    acc = payload["acceptance_matrix"]
    lines = [
        "# P2925/S1875 damping delta-source linear-system frontier certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Linear-algebra certificate",
        f"- variable count: `{la['variable_count']}`",
        f"- shape equation count: `{la['shape_equation_count']}`",
        f"- rank without target anchor: `{la['rank_without_target_anchor']}`",
        f"- nullity without target anchor: `{la['nullity_without_target_anchor']}`",
        f"- rank with external target anchor: `{la['rank_with_external_target_anchor']}`",
        f"- nullity with external target anchor: `{la['nullity_with_external_target_anchor']}`",
        f"- target anchor sourced by current rows: `{la['target_anchor_is_sourced_by_current_rows']}`",
        "",
        "## Acceptance",
        f"- P2923 readiness inherited: `{acc['p2923_prime_log_readiness_inherited']}`",
        f"- P2924 no-anchor inherited: `{acc['p2924_no_anchor_inherited']}`",
        f"- accepted current unlock objects: `{acc['accepted_current_unlock_object_count']}`",
        f"- no-new-live-frontier certificate exported: `{acc['no_new_live_frontier_certificate_exported']}`",
        "",
        "## Boundary",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload(read_json(P2923), read_json(P2924))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2925/S1875 damping delta-source linear-system frontier certificate", "## P2925/S1875 damping delta-source linear-system frontier certificate\n\n`P2925/S1875` strengthens P2924 with a finite linear system in variables `(y_2,...,y_11, delta)`.  The current prime-log character-shape equations `y_d - delta*log(d)=0` have rank `10` in `11` variables, leaving nullity `1`: exactly the free global slope line.  Adding the target row `delta=4/5` would raise the rank to `11` and close the nullity, but that row is an unsourced anchor, not a consequence of P2601/P2923/P2924 readiness.  The minimal unlock packet is a new strict prime-log value source, a new strict `delta=4/5` source law, or a combined strict damping `beta/eta` source packet.  No nonproxy `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2925/S1875 damping delta-source linear-system `L_total` guard", "## P2925/S1875 damping delta-source linear-system `L_total` guard\n\n`P2925/S1875` records a no-new-live-frontier certificate for the current damping `delta` source lane: the finite shape equations leave one free slope dimension, and the row selecting `delta=4/5` remains a missing strict source law.  Therefore the damping contribution cannot be inserted as role-bearing nonproxy `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE without a genuinely new source object.\n")
    append_once(AGENTS, "Current damping delta-source linear-system frontier guardrail (P2925/S1875, 2026-06-19)", "## Current damping delta-source linear-system frontier guardrail (P2925/S1875, 2026-06-19)\n\n- P2925 converts the P2924 slope/prime-anchor obstruction into a finite linear-algebra certificate on variables `(y_2,...,y_11, delta)`.\n- Existing character-shape rows have rank `10` and nullity `1`, leaving the global slope line free; adding `delta=4/5` would be an independent anchor row, but current artifacts do not source that row.\n- No strict prime-log value source, no strict `delta=4/5` source law, no strict damping `beta/eta`, no nonproxy `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE is exported.\n- A next admissible move must introduce one genuinely new source object: a strict prime-log value source law, a strict `delta=4/5` source law, or a combined `Strict_Damping_Beta_Eta_Source_Packet`; otherwise preserve the no-new-live-frontier certificate.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
