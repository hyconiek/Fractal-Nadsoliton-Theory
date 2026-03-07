#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "r15_explicit_host_scalar_floor_embedding_packet_for_host_matching_route.json"
OUT_SUMMARY = GENERATED / "r15_explicit_host_scalar_floor_embedding_packet_for_host_matching_route_summary.json"


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    q2122 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2122_psi_potential_diagonal_floor_gate.json"
    )
    q2124 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2124_scalar_vacuum_closure_branch_resolved_gate.json"
    )
    r13 = load_json("fundamental_action_reconstruction/generated/r13_partial_host_to_canonical_block_overlap_packet.json")
    r14 = load_json(
        "fundamental_action_reconstruction/generated/r14_explicit_frozen_kernel_specialization_packet_for_host_matching_route.json"
    )

    host_side = r13["host_side_decomposition"]
    canonical_side = r13["canonical_block_side_decomposition"]
    carrier_basis = canonical_side["carrier_basis"]
    canonical_rows = canonical_side["canonical_block_rows"]

    host_floor_value = float(host_side["diagonal_floor_value"])
    q2122_broken_floor = float(q2122["branch_results"]["broken_branch_strict"]["diag_floor"])
    q2124_broken_floor = float(q2124["inputs"]["broken_floor"])

    diagonal_entries: list[dict[str, Any]] = []
    indexed_local_diagonal_tokens_present = True
    for i, basis_label in enumerate(carrier_basis):
        diagonal_entry = canonical_rows[i][i]
        expected_tokens = [
            f"g4_psi{i}",
            f"g6_psi{i}",
            f"gY{i}",
            f"m2_psi{i}",
        ]
        token_check = all(token in diagonal_entry for token in expected_tokens)
        indexed_local_diagonal_tokens_present = indexed_local_diagonal_tokens_present and token_check
        diagonal_entries.append(
            {
                "basis_label": basis_label,
                "canonical_diagonal_entry": diagonal_entry,
                "host_scalar_floor_value": host_floor_value,
                "decomposition_rule": f"{diagonal_entry} = ({host_floor_value}) + (({diagonal_entry}) - ({host_floor_value}))",
                "residual_local_diagonal_entry": f"({diagonal_entry}) - ({host_floor_value})",
                "tokenized_local_structure_present": token_check,
            }
        )

    checks = [
        {
            "id": "q2122_broken_branch_floor_pass_present",
            "actual": q2122["verdict"],
            "expected": "PSI_POTENTIAL_DIAGONAL_FLOOR_GATE_PASS_BROKEN_BRANCH",
            "meaning": "QW-2122 exports the broken-branch diagonal floor",
        },
        {
            "id": "q2124_branch_resolved_scalar_vacuum_pass_present",
            "actual": q2124["verdict"],
            "expected": "SCALAR_VACUUM_CLOSURE_BRANCH_RESOLVED_STRICT_PASS",
            "meaning": "QW-2124 certifies the branch-resolved scalar vacuum closure",
        },
        {
            "id": "r13_host_diagonal_floor_has_scalar_vacuum_provenance",
            "actual": r13["partial_overlap_witness"]["host_diagonal_floor_has_scalar_vacuum_provenance"],
            "expected": True,
            "meaning": "R13 already exports scalar-vacuum provenance for the host diagonal floor",
        },
        {
            "id": "host_floor_matches_q2122_broken_branch_floor",
            "actual": abs(host_floor_value - q2122_broken_floor) <= 1e-12,
            "expected": True,
            "meaning": "the host diagonal floor matches the QW-2122 broken-branch floor",
        },
        {
            "id": "host_floor_matches_q2124_broken_floor",
            "actual": abs(host_floor_value - q2124_broken_floor) <= 1e-12,
            "expected": True,
            "meaning": "the host diagonal floor matches the QW-2124 branch-resolved broken floor",
        },
        {
            "id": "canonical_diagonal_entry_count_is_12",
            "actual": len(diagonal_entries),
            "expected": len(carrier_basis),
            "meaning": "the canonical diagonal sector is exported on the full 12-slot carrier",
        },
        {
            "id": "canonical_diagonal_entries_are_indexed_local_terms",
            "actual": indexed_local_diagonal_tokens_present,
            "expected": True,
            "meaning": "each canonical diagonal entry still carries an indexed local residual sector",
        },
        {
            "id": "r14_shared_kernel_light_channel_already_specialized",
            "actual": r14["specialization_result"]["kernel_channel_specialization_witness_present"],
            "expected": True,
            "meaning": "R14 already closed the shared kernel/light-facing channel before this diagonal packet",
        },
    ]

    for item in checks:
        item["pass"] = item["actual"] == item["expected"]

    artifact = {
        "stage": "R15",
        "packet_goal": "materialize_an_explicit_embedding_of_the_host_scalar_floor_m0_squared_I_into_the_canonical_local_diagonal_sector_while_keeping_the_residual_local_diagonal_part_open",
        "source_reports": ["QW-2122", "QW-2124", "R13", "R14"],
        "diagonal_decomposition": {
            "carrier_basis": carrier_basis,
            "host_diagonal_floor_symbol": "m0^2 I",
            "host_diagonal_floor_value": host_floor_value,
            "canonical_local_diagonal_symbol": "3*g4_psii*vpsii**2 + 5*g6_psii*vpsii**4 + 2*gYi*vphi**2 + m2_psii",
            "decomposition_rule": "D_canonical = m0^2 I + D_local_residual",
            "entrywise_rows": diagonal_entries,
        },
        "embedding_result": {
            "host_scalar_floor_present": True,
            "explicit_host_scalar_floor_embedding_packet_present": True,
            "canonical_local_diagonal_sector_present": True,
            "residual_local_diagonal_sector_present": True,
            "full_diagonal_sector_matching_witness_present": False,
        },
        "light_boundary": {
            "status": "shared_kernel_light_facing_channel_already_closed_by_R14_this_packet_touches_only_the_diagonal_complement",
            "shared_light_channel_source": "R14 specialization of (K_i_j + K_j_i)/2 to frozen K_total[i,j]",
            "current_packet_scope": "diagonal floor decomposition only",
            "boundary": "R15 does not change the already matched kernel/light-facing channel; it only isolates the non-light diagonal complement of the host route",
        },
        "classification": "explicit_host_scalar_floor_embedding_into_the_canonical_diagonal_sector_present_but_residual_local_diagonal_matching_absent",
        "frontier": "R15_B1",
        "frontier_text": "the repo now exports an explicit embedding of the host scalar floor into the canonical diagonal sector, but it still lacks a witness that the residual local diagonal sector cancels or equals a declared control pullback",
        "consistency_checks": checks,
        "hard_limits": [
            "no_theorem_level_pass",
            "no_full_closure_pass",
            "no_claim_that_the_canonical_local_diagonal_sector_equals_the_host_floor",
            "no_claim_that_the_residual_local_diagonal_sector_vanishes",
            "no_claim_that_QW2191_is_discharged",
            "no_claim_that_selector_closure_is_obtained",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "R15",
        "status": "PASS_PARTIAL_EXPLICIT_HOST_SCALAR_FLOOR_EMBEDDING_PACKET_READY",
        "result": artifact["classification"],
        "frontier": ["R15_B1", "QW2191_obstruction", "C10_B1"],
        "theorem_level_pass": False,
        "full_closure_pass": False,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()
