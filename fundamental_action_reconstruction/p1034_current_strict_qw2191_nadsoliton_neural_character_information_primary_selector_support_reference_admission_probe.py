#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-23"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_F952 = GENERATED / "f952_current_strict_qw2191_lambda_branch_info_primary_scpc_like_same_lane_stagnation_stop_packet_summary.json"
IN_P1011 = GENERATED / "p1011_current_strict_qw2191_lambda_branch_info_primary_scpc_like_selector_provider_shift_candidate_reference_lane_admission_probe.json"
IN_P1012 = GENERATED / "p1012_current_strict_qw2191_lambda_branch_info_primary_scpc_like_selector_interface_nonexport_audit_probe_summary.json"
IN_NADSOLITON_NEURAL = REPO / "NADSOLITON_NEURAL_ANALYSIS.md"
IN_INTERNAL_NEURAL = REPO / "INTERNAL_NEURAL_REPORT.md"
IN_MASTER_NEURAL = REPO / "MASTER_NEURAL_REPORT.md"
IN_QW540_REPORT = REPO / "raport_qw540_544_neural.md"

OUT_JSON = GENERATED / "p1034_current_strict_qw2191_nadsoliton_neural_character_information_primary_selector_support_reference_admission_probe.json"
OUT_SUMMARY = GENERATED / "p1034_current_strict_qw2191_nadsoliton_neural_character_information_primary_selector_support_reference_admission_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [
        IN_F952,
        IN_P1011,
        IN_P1012,
        IN_NADSOLITON_NEURAL,
        IN_INTERNAL_NEURAL,
        IN_MASTER_NEURAL,
        IN_QW540_REPORT,
    ]
    missing = [rel(path) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1034",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    f952 = load_json(IN_F952)
    p1011 = load_json(IN_P1011)
    p1012 = load_json(IN_P1012)
    nadsoliton_neural = load_text(IN_NADSOLITON_NEURAL)
    internal_neural = load_text(IN_INTERNAL_NEURAL)
    master_neural = load_text(IN_MASTER_NEURAL)
    qw540_report = load_text(IN_QW540_REPORT)

    theorem_result = p1011.get("theorem_result") or {}
    candidate_lane = p1011.get("candidate_reference_lane") or {}

    checks: list[dict[str, Any]] = []
    blocking: list[str] = []

    def add_check(check_id: str, actual: Any, expected: Any, meaning: str) -> None:
        passed = actual == expected
        checks.append(
            {
                "id": check_id,
                "actual": actual,
                "expected": expected,
                "pass": passed,
                "meaning": meaning,
            }
        )
        if not passed:
            blocking.append(check_id)

    same_lane_stop_packet_exported = (
        f952.get("status")
        == "PASS_CURRENT_STRICT_QW2191_LAMBDA_BRANCH_INFO_PRIMARY_SCPC_LIKE_SAME_LANE_STAGNATION_STOP_PACKET_EXPORTED"
        and f952.get("packet_exported_on_current_repo_state") is True
        and f952.get("same_lane_stagnation_boundary_reached") is True
    )

    current_candidate_lane_reference_only = (
        p1011.get("status")
        == "P1011_CURRENT_STRICT_QW2191_LAMBDA_BRANCH_INFO_PRIMARY_SCPC_LIKE_SELECTOR_PROVIDER_SHIFT_CANDIDATE_REFERENCE_LANE_ADMITTED_SELECTOR_INTERFACE_BLOCKED"
        and theorem_result.get("info_primary_scpc_like_selector_provider_shift_candidate_reference_lane_admitted")
        is True
        and theorem_result.get("candidate_reference_lane_grade") == "reference_context_candidate_only"
        and theorem_result.get("strict_selector_interface_exported") is False
        and theorem_result.get("strict_selector_source_exported") is False
        and theorem_result.get("neural_network_identity_claim_admitted") is False
        and theorem_result.get("energy_based_identity_claim_admitted") is False
    )

    exact_selector_interface_still_blocked = (
        p1012.get("status")
        == "P1012_CURRENT_STRICT_QW2191_LAMBDA_BRANCH_INFO_PRIMARY_SCPC_LIKE_SELECTOR_INTERFACE_NONEXPORT_AUDITED"
        and p1012.get("exact_selector_interface_exported_on_current_repo_state") is False
    )

    nadsoliton_neural_shape_exported = all(
        needle in nadsoliton_neural
        for needle in [
            "behaves fundamentally as a Neural Network",
            "Modified Oja's Rule",
            "Hopfield Network / Attractor Network",
            "spatially distributed, self-organizing neural network",
        ]
    )

    internal_neural_shape_exported = all(
        needle in internal_neural
        for needle in [
            "Lateral Inhibition",
            "Sub-critical",
            "Associative Memory",
        ]
    )

    master_neural_shape_exported = all(
        needle in master_neural
        for needle in [
            "Hebbian Emergence",
            "Lateral Inhibition",
            "E/I Balance",
        ]
    )

    qw540_neural_dynamics_exported = all(
        needle in qw540_report
        for needle in [
            "Grawitacja Hebbowska",
            "Ciemna Energia (Neural Forgetting)",
            "Pamiec Czastek",
        ]
    ) or all(
        needle in qw540_report
        for needle in [
            "Grawitacja Hebbowska",
            "Ciemna Energia (Neural Forgetting)",
            "Pami\u0119\u0107 Cz\u0105stek",
        ]
    )

    neural_character_support_reference_admitted = (
        same_lane_stop_packet_exported
        and current_candidate_lane_reference_only
        and exact_selector_interface_still_blocked
        and nadsoliton_neural_shape_exported
        and internal_neural_shape_exported
        and master_neural_shape_exported
        and qw540_neural_dynamics_exported
    )

    add_check(
        "same_lane_stop_packet_exported",
        same_lane_stop_packet_exported,
        True,
        "F952 already stops further same-lane SCPC-like descent as a primary move.",
    )
    add_check(
        "current_candidate_lane_reference_only",
        current_candidate_lane_reference_only,
        True,
        "The active Lambda lane remains admitted only as an information-primary reference-context candidate lane.",
    )
    add_check(
        "exact_selector_interface_still_blocked",
        exact_selector_interface_still_blocked,
        True,
        "The exact strict selector interface remains unexported on the current repo state.",
    )
    add_check(
        "nadsoliton_neural_shape_exported",
        nadsoliton_neural_shape_exported,
        True,
        "The dedicated nadsoliton neural analysis exports Hebbian, Oja-type, and attractor/neural-network language.",
    )
    add_check(
        "internal_neural_shape_exported",
        internal_neural_shape_exported,
        True,
        "The internal neural report exports lateral inhibition, subcritical inhibition, and associative-memory language.",
    )
    add_check(
        "master_neural_shape_exported",
        master_neural_shape_exported,
        True,
        "The master neural report exports Hebbian emergence, inhibition, and E/I-balance framing.",
    )
    add_check(
        "qw540_neural_dynamics_exported",
        qw540_neural_dynamics_exported,
        True,
        "The QW-540..544 neural report exports Hebbian gravity, neural forgetting, and particle-memory language.",
    )
    add_check(
        "neural_character_support_reference_admitted",
        neural_character_support_reference_admitted,
        True,
        "Therefore the neural-character corpus may lawfully support the current information-primary selector hypothesis at support-reference grade only.",
    )

    status = (
        "P1034_CURRENT_STRICT_QW2191_NADSOLITON_NEURAL_CHARACTER_INFORMATION_PRIMARY_SELECTOR_SUPPORT_REFERENCE_ADMITTED_INTERFACE_STILL_BLOCKED"
        if not blocking and neural_character_support_reference_admitted
        else "FAIL_CURRENT_STRICT_QW2191_NADSOLITON_NEURAL_CHARACTER_INFORMATION_PRIMARY_SELECTOR_SUPPORT_REFERENCE_ADMISSION"
    )

    artifact = {
        "stage": "P1034",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "f952_same_lane_stagnation_stop_packet_summary": rel(IN_F952),
            "p1011_candidate_reference_lane_admission_probe": rel(IN_P1011),
            "p1012_selector_interface_nonexport_summary": rel(IN_P1012),
            "nadsoliton_neural_analysis": rel(IN_NADSOLITON_NEURAL),
            "internal_neural_report": rel(IN_INTERNAL_NEURAL),
            "master_neural_report": rel(IN_MASTER_NEURAL),
            "qw540_544_neural_report": rel(IN_QW540_REPORT),
        },
        "theorem_result": {
            "nadsoliton_neural_character_support_reference_admitted": neural_character_support_reference_admitted,
            "support_reference_grade": "cross_repo_support_reference_only",
            "support_reference_id": "nadsoliton_neural_character_information_primary_selector_support_reference_v1",
            "supported_candidate_lane_ref": theorem_result.get("candidate_reference_lane_id"),
            "same_lane_stop_constraint_active": same_lane_stop_packet_exported,
            "strict_selector_interface_exported": False,
            "strict_selector_source_exported": False,
            "theorem_level_neural_identity_exported": False,
            "energy_based_selector_identity_exported": False,
            "admitted_character_features": [
                "hebbian_learning",
                "attractor_memory",
                "lateral_inhibition",
                "subcritical_inhibition",
                "information_primary_processing",
            ],
            "next_honest_move_requires_exact_bridge_from_support_reference_to_selector_interface": True,
            "no_false_pass": True,
        },
        "support_reference": {
            "support_reference_id": "nadsoliton_neural_character_information_primary_selector_support_reference_v1",
            "supported_candidate_lane_ref": candidate_lane.get("candidate_id"),
            "support_reference_grade": "cross_repo_support_reference_only",
            "reading_contract": "use_neural_character_as_information_primary_selector_support_reference_only",
        },
        "checks": checks,
        "blocking_checks": blocking,
        "current_honest_reading": [
            "The repo does contain explicit neural-character evidence relevant to selector work.",
            "That evidence may now lawfully support the active information-primary selector hypothesis only at cross-repo support-reference grade.",
            "The exact strict selector interface and strict selector source remain blocked.",
        ],
        "recommended_next_packet": {
            "id": "F953_CURRENT_STRICT_QW2191_NADSOLITON_NEURAL_CHARACTER_INFORMATION_PRIMARY_SELECTOR_SUPPORT_REFERENCE_PACKET",
            "goal": "Package the neural-character corpus only as a support-reference packet for the active information-primary selector hypothesis.",
            "export_object_id": "nadsoliton_neural_character_information_primary_selector_support_reference_v1",
        },
        "no_false_pass": True,
    }
    write_json(OUT_JSON, artifact)

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "support_reference_id": artifact["theorem_result"]["support_reference_id"],
        "supported_candidate_lane_ref": artifact["theorem_result"]["supported_candidate_lane_ref"],
        "support_reference_grade": artifact["theorem_result"]["support_reference_grade"],
        "nadsoliton_neural_character_support_reference_admitted": artifact["theorem_result"][
            "nadsoliton_neural_character_support_reference_admitted"
        ],
        "strict_selector_interface_exported": artifact["theorem_result"]["strict_selector_interface_exported"],
        "strict_selector_source_exported": artifact["theorem_result"]["strict_selector_source_exported"],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
