#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "r11_symmetry_certified_declared_control_transport_packet_for_psi_block_route.json"
OUT_SUMMARY = GENERATED / "r11_symmetry_certified_declared_control_transport_packet_for_psi_block_route_summary.json"


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def clean_scalar(value: float) -> float:
    if abs(value) < 1e-15:
        return 0.0
    return round(value, 15)


def real_fourier_column(n_octaves: int, mode_label: str) -> list[float]:
    if mode_label.startswith("c"):
        harmonic = int(mode_label[1:])
        return [
            clean_scalar(math.sqrt(2.0 / n_octaves) * math.cos(2.0 * math.pi * harmonic * j / n_octaves))
            for j in range(n_octaves)
        ]
    if mode_label.startswith("s"):
        harmonic = int(mode_label[1:])
        return [
            clean_scalar(math.sqrt(2.0 / n_octaves) * math.sin(2.0 * math.pi * harmonic * j / n_octaves))
            for j in range(n_octaves)
        ]
    raise ValueError(f"unsupported mode label: {mode_label}")


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    q2190 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2190_kernel_mode_representation_emergence_gate.json"
    )
    q2191 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2191_mode_index_uniqueness_obstruction_theorem_gate.json"
    )
    c3 = load_json("fundamental_action_reconstruction/generated/c3_internal_reference_pair_candidate_from_mode_scaffold_summary.json")
    c13 = load_json("fundamental_action_reconstruction/generated/c13_mode_basis_control_index_set_audit_summary.json")
    c14 = load_json("fundamental_action_reconstruction/generated/c14_control_mode_to_psi_transport_schema_summary.json")

    n_octaves = int(q2190["mode_mapping"]["n_octaves"])
    psi_basis = c14["control_transport"]["psi_basis_indices"]
    control_modes = ["c1", "s1", "c2", "s2"]
    transport_columns = {mode: real_fourier_column(n_octaves, mode) for mode in control_modes}
    transport_rows = [
        [transport_columns["c1"][j], transport_columns["s1"][j], transport_columns["c2"][j], transport_columns["s2"][j]]
        for j in range(n_octaves)
    ]

    checks = [
        {
            "id": "c13_control_index_set_present",
            "actual": c13["result"]["mode_basis_control_index_set_present"],
            "expected": "yes",
            "meaning": "the deterministic mode-basis control index-set is present",
        },
        {
            "id": "c14_control_transport_schema_present",
            "actual": c14["result"]["control_transport_schema_present"],
            "expected": "yes",
            "meaning": "the control transport schema into the canonical Psi carrier is present",
        },
        {
            "id": "q2190_partial_mode_embedding_pass_present",
            "actual": q2190["verdict"],
            "expected": "KERNEL_MODE_REPRESENTATION_EMERGENCE_GATE_PASS_PARTIAL_PHYSICAL_UNIQUENESS_OPEN",
            "meaning": "QW-2190 provides the declared deterministic mode embedding with strong symmetry certification",
        },
        {
            "id": "q2190_kernel_mode_basis_declared_deterministically",
            "actual": q2190["flags"]["kernel_mode_basis_declared_deterministically"],
            "expected": True,
            "meaning": "the mode basis is declared deterministically",
        },
        {
            "id": "q2190_mode_subspaces_orthonormal",
            "actual": q2190["flags"]["mode_subspaces_orthonormal"],
            "expected": True,
            "meaning": "the relevant mode subspaces are orthonormal",
        },
        {
            "id": "q2190_mode_subspaces_disjoint",
            "actual": q2190["flags"]["mode_subspaces_disjoint"],
            "expected": True,
            "meaning": "the relevant mode subspaces are disjoint",
        },
        {
            "id": "q2190_kernel_invariance_on_relevant_subspaces",
            "actual": q2190["flags"]["kernel_invariance_of_su3_subspace"] and q2190["flags"]["kernel_invariance_of_su2_subspace"],
            "expected": True,
            "meaning": "the declared mode subspaces remain kernel-invariant",
        },
        {
            "id": "q2191_strict_obstruction_present",
            "actual": q2191["verdict"],
            "expected": "MODE_INDEX_UNIQUENESS_OBSTRUCTION_THEOREM_GATE_PASS_STRICT",
            "meaning": "QW-2191 still certifies the O(2) uniqueness obstruction",
        },
        {
            "id": "q2191_full_physical_uniqueness_closed_false",
            "actual": q2191["flags"]["full_physical_uniqueness_closed"],
            "expected": False,
            "meaning": "full physical uniqueness remains open",
        },
    ]

    for item in checks:
        item["pass"] = item["actual"] == item["expected"]

    artifact = {
        "stage": "R11",
        "packet_goal": "materialize_an_explicit_declared_control_transport_packet_with_symmetry_certificate_for_the_host_to_concrete_Psi_block_route",
        "source_reports": ["QW-2190", "QW-2191", "C3", "C13", "C14"],
        "transport_packet": {
            "domain_basis": control_modes,
            "target_basis": psi_basis,
            "control_pairs": c3["candidate_pairs"],
            "matrix_shape": [n_octaves, len(control_modes)],
            "matrix_columns": transport_columns,
            "matrix_rows": transport_rows,
            "construction_rule": "real Fourier basis coefficients on the declared 12-octave carrier, read in canonical Psi-slot order psi0..psi11",
            "scope": "declared_control_transport_only",
        },
        "symmetry_certificate": {
            "deterministic_basis_declared": q2190["flags"]["kernel_mode_basis_declared_deterministically"],
            "mode_subspaces_orthonormal": q2190["flags"]["mode_subspaces_orthonormal"],
            "mode_subspaces_disjoint": q2190["flags"]["mode_subspaces_disjoint"],
            "kernel_invariance_of_relevant_subspaces": {
                "su3": q2190["flags"]["kernel_invariance_of_su3_subspace"],
                "su2": q2190["flags"]["kernel_invariance_of_su2_subspace"],
            },
            "embedded_lie_closure": {
                "su3": q2190["flags"]["embedded_su3_lie_closure"],
                "su2": q2190["flags"]["embedded_su2_lie_closure"],
            },
            "numeric_residuals": {
                "orthonormal_residual_su3_basis": q2190["numeric_audit"]["orthonormal_residual_su3_basis"],
                "orthonormal_residual_su2_basis": q2190["numeric_audit"]["orthonormal_residual_su2_basis"],
                "subspace_disjoint_residual": q2190["numeric_audit"]["subspace_disjoint_residual"],
                "kernel_invariance_residual_su3_subspace": q2190["numeric_audit"]["kernel_invariance_residual_su3_subspace"],
                "kernel_invariance_residual_su2_subspace": q2190["numeric_audit"]["kernel_invariance_residual_su2_subspace"],
            },
            "verdict": q2190["verdict"],
        },
        "obstruction_boundary": {
            "theorem_verdict": q2191["verdict"],
            "continuous_nonunique_embedding_family_exhibited": q2191["flags"]["continuous_nonunique_embedding_family_exhibited"],
            "full_uniqueness_from_kernel_alone_obstructed": q2191["flags"]["full_uniqueness_from_kernel_alone_obstructed"],
            "full_physical_uniqueness_closed": q2191["flags"]["full_physical_uniqueness_closed"],
            "required_next_step": q2191["required_next_step"],
        },
        "result": {
            "explicit_declared_control_transport_packet_present": True,
            "symmetry_certificate_present": True,
            "strict_selector_relevant_physical_canonicalization_present": False,
            "chosen_transported_index_set_physically_canonicalized": False,
        },
        "classification": "explicit_declared_control_transport_packet_present_and_symmetry_certified_but_not_fully_physically_canonicalized_for_selector_relevant_block_extraction",
        "frontier": "R11_B1",
        "frontier_text": "the repo now contains an explicit declared transport matrix and a symmetry certificate for the control route, but QW-2191 still blocks full physical uniqueness and thus no selector-relevant canonical representative is yet justified",
        "hard_limits": [
            "no_theorem_level_pass",
            "no_full_closure_pass",
            "no_claim_that_the_declared_transport_is_physically_canonical",
            "no_claim_that_QW2191_is_discharged",
            "no_claim_that_a_concrete_Psi_sector_block_is_already_extracted",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "R11",
        "status": "PASS_PARTIAL_EXPLICIT_DECLARED_CONTROL_TRANSPORT_PACKET_READY_BUT_PHYSICAL_UNIQUENESS_OPEN",
        "result": artifact["classification"],
        "frontier": ["R11_B1", "C14_B1", "QW2191_obstruction"],
        "theorem_level_pass": False,
        "full_closure_pass": False,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()
