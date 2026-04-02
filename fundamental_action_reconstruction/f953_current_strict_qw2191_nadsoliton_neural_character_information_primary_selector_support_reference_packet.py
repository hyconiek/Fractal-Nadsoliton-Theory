#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-23"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1034 = GENERATED / "p1034_current_strict_qw2191_nadsoliton_neural_character_information_primary_selector_support_reference_admission_probe.json"
IN_F952 = GENERATED / "f952_current_strict_qw2191_lambda_branch_info_primary_scpc_like_same_lane_stagnation_stop_packet_summary.json"

OUT_JSON = GENERATED / "f953_current_strict_qw2191_nadsoliton_neural_character_information_primary_selector_support_reference_packet.json"
OUT_SUMMARY = GENERATED / "f953_current_strict_qw2191_nadsoliton_neural_character_information_primary_selector_support_reference_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P1034, IN_F952]
    missing = [rel(path) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "F953",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1034 = load_json(IN_P1034)
    f952 = load_json(IN_F952)

    theorem_result = p1034.get("theorem_result") or {}
    support_reference = p1034.get("support_reference") or {}

    status = (
        "F953_EXECUTED_CURRENT_STRICT_QW2191_NADSOLITON_NEURAL_CHARACTER_INFORMATION_PRIMARY_SELECTOR_SUPPORT_REFERENCE_PACKET_NO_FALSE_PASS"
        if p1034.get("status")
        == "P1034_CURRENT_STRICT_QW2191_NADSOLITON_NEURAL_CHARACTER_INFORMATION_PRIMARY_SELECTOR_SUPPORT_REFERENCE_ADMITTED_INTERFACE_STILL_BLOCKED"
        and theorem_result.get("nadsoliton_neural_character_support_reference_admitted") is True
        and theorem_result.get("support_reference_grade") == "cross_repo_support_reference_only"
        and theorem_result.get("strict_selector_interface_exported") is False
        and theorem_result.get("strict_selector_source_exported") is False
        and f952.get("status")
        == "PASS_CURRENT_STRICT_QW2191_LAMBDA_BRANCH_INFO_PRIMARY_SCPC_LIKE_SAME_LANE_STAGNATION_STOP_PACKET_EXPORTED"
        else "F953_REQUIRES_REVIEW"
    )

    artifact = {
        "stage": "F953",
        "packet_name": "CurrentStrictQW2191NadsolitonNeuralCharacterInformationPrimarySelectorSupportReference_v1",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p1034_support_reference_admission_probe": rel(IN_P1034),
            "f952_same_lane_stagnation_stop_packet_summary": rel(IN_F952),
        },
        "exported_object": {
            "object_id": theorem_result.get("support_reference_id"),
            "goal": "Package the nadsoliton neural-character corpus only as a support-reference for the current information-primary selector hypothesis.",
            "support_reference_grade": theorem_result.get("support_reference_grade"),
            "supported_candidate_lane_ref": theorem_result.get("supported_candidate_lane_ref"),
            "same_lane_stop_constraint_active": theorem_result.get("same_lane_stop_constraint_active"),
            "admitted_character_features": theorem_result.get("admitted_character_features"),
            "strict_selector_interface_status": "blocked_nonexport",
            "strict_selector_source_status": "blocked_nonexport",
            "identity_contract": {
                "theorem_level_neural_identity_exported": False,
                "energy_based_selector_identity_exported": False,
                "support_reference_only": True,
            },
            "reading_contract": support_reference.get("reading_contract"),
        },
        "current_honest_reading": [
            "The repo now exports one lawful neural-character support-reference for the active information-primary selector hypothesis.",
            "That support-reference remains cross-repo and reference-only; it does not realize selector interface or selector source.",
            "The current work contract is to use neural character as support-reference, not as fake strict closure.",
        ],
        "recommended_next_move": "Audit whether any exact bridge from this support-reference into an exact strict selector-interface question already exists before freezing a new bridge target.",
        "hard_limits": [
            "Does not claim strict selector interface export.",
            "Does not claim strict selector source export.",
            "Does not claim T176 export.",
            "Does not claim QW-2191 discharge.",
            "Does not claim theorem-level neural identity in the strict front.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F953",
        "status": status,
        "as_of": AS_OF,
        "exported_object_id": artifact["exported_object"]["object_id"],
        "support_reference_grade": artifact["exported_object"]["support_reference_grade"],
        "supported_candidate_lane_ref": artifact["exported_object"]["supported_candidate_lane_ref"],
        "strict_selector_interface_status": artifact["exported_object"]["strict_selector_interface_status"],
        "strict_selector_source_status": artifact["exported_object"]["strict_selector_source_status"],
        "recommended_next_move": artifact["recommended_next_move"],
        "no_false_pass": True,
    }
    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
