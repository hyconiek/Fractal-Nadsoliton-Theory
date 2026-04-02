#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-23"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1035 = GENERATED / "p1035_current_strict_qw2191_nadsoliton_neural_character_support_reference_to_exact_selector_interface_bridge_nonexport_audit_probe.json"
IN_F953 = GENERATED / "f953_current_strict_qw2191_nadsoliton_neural_character_information_primary_selector_support_reference_packet_summary.json"

OUT_JSON = GENERATED / "f954_current_strict_qw2191_nadsoliton_neural_character_support_reference_to_exact_selector_interface_bridge_target_packet.json"
OUT_SUMMARY = GENERATED / "f954_current_strict_qw2191_nadsoliton_neural_character_support_reference_to_exact_selector_interface_bridge_target_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P1035, IN_F953]
    missing = [rel(path) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "F954",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1035 = load_json(IN_P1035)
    f953 = load_json(IN_F953)

    theorem_result = p1035.get("theorem_result") or {}
    bridge_target = p1035.get("exact_missing_bridge_target_candidate") or {}

    status = (
        "F954_EXECUTED_CURRENT_STRICT_QW2191_NADSOLITON_NEURAL_CHARACTER_SUPPORT_REFERENCE_TO_EXACT_SELECTOR_INTERFACE_BRIDGE_TARGET_PACKET_NO_FALSE_PASS"
        if p1035.get("status")
        == "P1035_CURRENT_STRICT_QW2191_NADSOLITON_NEURAL_CHARACTER_SUPPORT_REFERENCE_TO_EXACT_SELECTOR_INTERFACE_BRIDGE_NONEXPORT_AUDITED"
        and theorem_result.get("next_honest_move_requires_freezing_exact_bridge_target") is True
        and f953.get("status")
        == "F953_EXECUTED_CURRENT_STRICT_QW2191_NADSOLITON_NEURAL_CHARACTER_INFORMATION_PRIMARY_SELECTOR_SUPPORT_REFERENCE_PACKET_NO_FALSE_PASS"
        else "F954_REQUIRES_REVIEW"
    )

    artifact = {
        "stage": "F954",
        "packet_name": "CurrentStrictQW2191NadsolitonNeuralCharacterSupportReferenceToExactSelectorInterfaceBridgeTarget_v1",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p1035_bridge_nonexport_audit_probe": rel(IN_P1035),
            "f953_support_reference_packet_summary": rel(IN_F953),
        },
        "exported_object": {
            "object_id": bridge_target.get("exact_bridge_target_id"),
            "goal": "Freeze the exact missing bridge target from the neural-character support-reference into the exact strict selector-interface question for the active candidate lane.",
            "bridge_target_grade": "future_only_exact_bridge_target",
            "support_reference_ref": bridge_target.get("support_reference_ref"),
            "supported_candidate_lane_ref": bridge_target.get("supported_candidate_lane_ref"),
            "target_contract": "exact_bridge_from_support_reference_to_selector_interface_question_only",
            "strict_selector_interface_status": "blocked_nonexport",
            "strict_selector_source_status": "blocked_nonexport",
            "reading_contract": "freeze_exact_bridge_target_not_fake_interface_realization",
        },
        "current_honest_reading": [
            "The repo now names one exact missing bridge target between the neural-character support-reference and the strict selector-interface question.",
            "That bridge remains future-only and does not yet realize selector interface or selector source.",
            "The current work contract is to freeze the bridge target, not to fake interface realization.",
        ],
        "recommended_next_move": "Audit whether the frozen bridge target already has any actual realization on the current repo state before attempting a new bridge realization packet.",
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
        "stage": "F954",
        "status": status,
        "as_of": AS_OF,
        "exported_object_id": artifact["exported_object"]["object_id"],
        "bridge_target_grade": artifact["exported_object"]["bridge_target_grade"],
        "support_reference_ref": artifact["exported_object"]["support_reference_ref"],
        "supported_candidate_lane_ref": artifact["exported_object"]["supported_candidate_lane_ref"],
        "strict_selector_interface_status": artifact["exported_object"]["strict_selector_interface_status"],
        "recommended_next_move": artifact["recommended_next_move"],
        "no_false_pass": True,
    }
    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
