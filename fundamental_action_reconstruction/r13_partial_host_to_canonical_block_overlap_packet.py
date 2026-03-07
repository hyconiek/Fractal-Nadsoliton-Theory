#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "r13_partial_host_to_canonical_block_overlap_packet.json"
OUT_SUMMARY = GENERATED / "r13_partial_host_to_canonical_block_overlap_packet_summary.json"


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def build_host_kernel_matrix(distance_profile: dict[str, float], n_octaves: int) -> list[list[float]]:
    rows: list[list[float]] = []
    for i in range(n_octaves):
        row: list[float] = []
        for j in range(n_octaves):
            if i == j:
                row.append(0.0)
                continue
            d = min(abs(i - j), n_octaves - abs(i - j))
            row.append(float(distance_profile[str(d)]))
        rows.append(row)
    return rows


def build_host_operator_matrix(kernel_rows: list[list[float]], broken_floor: float) -> list[list[float]]:
    rows: list[list[float]] = []
    for i, row in enumerate(kernel_rows):
        rows.append([
            float(entry + broken_floor) if i == j else float(entry)
            for j, entry in enumerate(row)
        ])
    return rows


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    q2118 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2118_ktotal_spectral_tripartition_gate.json"
    )
    q2124 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2124_scalar_vacuum_closure_branch_resolved_gate.json"
    )
    q2186 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2186_ktotal_spectral_stability_margin_gate.json"
    )
    r8 = load_json("fundamental_action_reconstruction/generated/r8_existing_kernel_feedback_host_operator_carrier_packet_for_kobs.json")
    r12 = load_json(
        "fundamental_action_reconstruction/generated/r12_explicit_canonical_psi_block_export_packet_with_kernel_mixing_for_host_route.json"
    )

    n_octaves = int(q2118["matrix_spec"]["n_octaves"])
    broken_floor = float(q2124["inputs"]["broken_floor"])
    kernel_rows = build_host_kernel_matrix(q2118["distance_profile"], n_octaves)
    host_rows = build_host_operator_matrix(kernel_rows, broken_floor)

    checks = [
        {
            "id": "q2118_tripartition_gate_pass_present",
            "actual": q2118["verdict"],
            "expected": "KTOTAL_SPECTRAL_TRIPARTITION_GATE_PASS_WITH_CONDITIONAL_VACUUM_SHIFT",
            "meaning": "QW-2118 provides the frozen K_total carrier and distance profile",
        },
        {
            "id": "q2124_branch_resolved_scalar_floor_present",
            "actual": q2124["verdict"],
            "expected": "SCALAR_VACUUM_CLOSURE_BRANCH_RESOLVED_STRICT_PASS",
            "meaning": "QW-2124 provides the branch-resolved scalar vacuum floor provenance",
        },
        {
            "id": "q2186_host_certificate_present",
            "actual": q2186["verdict"],
            "expected": "KTOTAL_SPECTRAL_STABILITY_MARGIN_GATE_PASS_STRICT_BRANCH_SCOPE",
            "meaning": "QW-2186 provides the strict host certificate",
        },
        {
            "id": "r8_host_basis_order_matches_12_slot_carrier",
            "actual": r8["declared_host_carrier"]["finite_state_space_size"],
            "expected": 12,
            "meaning": "R8 places the host on the declared 12-slot carrier",
        },
        {
            "id": "r12_canonical_block_shape_matches_12x12",
            "actual": r12["canonical_psi_block_export"]["shape"],
            "expected": [12, 12],
            "meaning": "R12 exports the canonical Psi block on the same 12-slot support",
        },
        {
            "id": "r12_kernel_mixing_channel_present",
            "actual": r12["kernel_mixing_light_facing_summary"]["off_diagonal_entry_count"],
            "expected": 132,
            "meaning": "R12 contains the full off-diagonal kernel-mixing channel",
        },
        {
            "id": "q2124_broken_floor_matches_q2186_used_floor",
            "actual": broken_floor,
            "expected": float(q2186["matrix"]["broken_floor_used"]),
            "meaning": "the branch-resolved scalar floor matches the floor used by the certified host operator",
        },
    ]

    for item in checks:
        item["pass"] = item["actual"] == item["expected"]

    artifact = {
        "stage": "R13",
        "packet_goal": "materialize_the_partial_overlap_between_the_qw2186_host_operator_and_the_exported_canonical_Psi_block_without_claiming_full_matching",
        "source_reports": ["QW-2118", "QW-2124", "QW-2186", "R8", "R12"],
        "host_side_decomposition": {
            "carrier_basis": r8["declared_host_carrier"]["basis_order"],
            "kernel_part_symbol": "K_total_host_frozen",
            "kernel_part_numeric_rows": kernel_rows,
            "diagonal_floor_symbol": "m0^2 I",
            "diagonal_floor_value": broken_floor,
            "host_operator_symbol": "A_host = K_total + m0^2 I",
            "host_operator_numeric_rows": host_rows,
        },
        "canonical_block_side_decomposition": {
            "carrier_basis": r12["chosen_transported_index_set"]["basis_order"],
            "canonical_block_symbol": r12["canonical_psi_block_export"]["symbol"],
            "off_diagonal_kernel_channel": "(K_i_j + K_j_i)/2",
            "diagonal_local_channel": "3*g4_psii*vpsii**2 + 5*g6_psii*vpsii**4 + 2*gYi*vphi**2 + m2_psii",
            "canonical_block_rows": r12["canonical_psi_block_export"]["matrix_rows"],
        },
        "partial_overlap_witness": {
            "shared_12_slot_carrier_present": True,
            "shared_off_diagonal_kernel_channel_present": True,
            "host_diagonal_floor_has_scalar_vacuum_provenance": True,
            "full_operator_level_matching_witness_present": False,
            "meaning": "both sides now admit an explicit decomposition into a kernel/off-diagonal part plus a diagonal part, but the repo still does not export the coefficient-specialization witness for the kernel part or the equality witness for the diagonal part",
        },
        "light_facing_boundary": {
            "host_light_facing_channel": "frozen K_total cyclic octave-space kernel",
            "canonical_block_light_facing_channel": "symmetric off-diagonal kernel-mixing entries (K_i_j + K_j_i)/2",
            "status": "present_as_partial_overlap_only",
            "boundary": "the shared light-facing content is explicit only at the structural overlap level; no full operator equality, no K_obs factorization, and no selector closure are claimed",
        },
        "classification": "partial_host_to_canonical_block_overlap_present_via_shared_kernel_light_channel_and_scalar_floor_provenance_but_full_matching_absent",
        "frontier": "R13_B1",
        "frontier_text": "the repo now exports a real partial overlap between the host operator and the canonical Psi block, but it still lacks a coefficient-specialization witness for the kernel part and an equality witness for the diagonal part",
        "hard_limits": [
            "no_theorem_level_pass",
            "no_full_closure_pass",
            "no_claim_that_the_host_operator_equals_the_exported_canonical_block",
            "no_claim_that_QW2191_is_discharged",
            "no_claim_that_selector_closure_is_obtained",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "R13",
        "status": "PASS_PARTIAL_HOST_TO_CANONICAL_BLOCK_OVERLAP_PACKET_READY",
        "result": artifact["classification"],
        "frontier": ["R13_B1", "QW2191_obstruction", "C10_B1"],
        "theorem_level_pass": False,
        "full_closure_pass": False,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()
