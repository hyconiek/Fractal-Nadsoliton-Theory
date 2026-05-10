#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-23"
ACTUAL_REALIZATION_ID = (
    "nadsoliton_neural_character_support_reference_to_exact_selector_interface_bridge_target_actual_realization_v1"
)

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_F954 = GENERATED / "f954_current_strict_qw2191_nadsoliton_neural_character_support_reference_to_exact_selector_interface_bridge_target_packet_summary.json"
IN_P1035 = GENERATED / "p1035_current_strict_qw2191_nadsoliton_neural_character_support_reference_to_exact_selector_interface_bridge_nonexport_audit_probe.json"
IN_P1012 = GENERATED / "p1012_current_strict_qw2191_lambda_branch_info_primary_scpc_like_selector_interface_nonexport_audit_probe_summary.json"

OUT_JSON = GENERATED / "p1036_current_strict_qw2191_nadsoliton_neural_character_support_reference_to_exact_selector_interface_bridge_target_actual_realization_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p1036_current_strict_qw2191_nadsoliton_neural_character_support_reference_to_exact_selector_interface_bridge_target_actual_realization_nonexport_audit_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def repo_scan_hits_for_actual_realization() -> list[str]:
    hits: list[str] = []
    ignored_markers = (
        "P1036_CURRENT_STRICT_QW2191_NADSOLITON_NEURAL_CHARACTER_SUPPORT_REFERENCE_TO_EXACT_SELECTOR_INTERFACE_BRIDGE_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDIT_PROBE",
        "N869_CURRENT_STRICT_QW2191_NADSOLITON_NEURAL_CHARACTER_SUPPORT_REFERENCE_TO_EXACT_SELECTOR_INTERFACE_BRIDGE_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDIT_THEOREM",
        "T297_CURRENT_STRICT_QW2191_NADSOLITON_NEURAL_CHARACTER_SUPPORT_REFERENCE_TO_EXACT_SELECTOR_INTERFACE_BRIDGE_TARGET_ACTUAL_REALIZATION_ATTEMPT_SPEC",
        "N870_CURRENT_STRICT_QW2191_NADSOLITON_NEURAL_CHARACTER_SUPPORT_REFERENCE_TO_EXACT_SELECTOR_INTERFACE_BRIDGE_TARGET_ACTUAL_REALIZATION_ATTEMPT_THEOREM",
        "p1036_current_strict_qw2191_nadsoliton_neural_character_support_reference_to_exact_selector_interface_bridge_target_actual_realization_nonexport_audit_probe",
        "p1037_current_strict_qw2191_nadsoliton_neural_character_support_reference_to_exact_selector_interface_bridge_target_actual_realization_attempt_probe",
    )
    for path in ROOT.rglob("*"):
        if not path.is_file():
            continue
        if path.suffix.lower() not in {".md", ".py", ".json", ".txt"}:
            continue
        rel_path = rel(path)
        if any(marker in rel_path for marker in ignored_markers):
            continue
        try:
            text = path.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            continue
        if ACTUAL_REALIZATION_ID in text:
            hits.append(rel_path)
    return sorted(hits)


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_F954, IN_P1035, IN_P1012]
    missing = [rel(path) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1036",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    f954 = load_json(IN_F954)
    p1035 = load_json(IN_P1035)
    p1012 = load_json(IN_P1012)

    theorem_result = p1035.get("theorem_result") or {}
    bridge_target = p1035.get("exact_missing_bridge_target_candidate") or {}
    scan_hits = repo_scan_hits_for_actual_realization()

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

    bridge_target_packet_exported = (
        f954.get("status")
        == "F954_EXECUTED_CURRENT_STRICT_QW2191_NADSOLITON_NEURAL_CHARACTER_SUPPORT_REFERENCE_TO_EXACT_SELECTOR_INTERFACE_BRIDGE_TARGET_PACKET_NO_FALSE_PASS"
        and f954.get("exported_object_id")
        == "nadsoliton_neural_character_support_reference_to_exact_selector_interface_bridge_target_v1"
    )

    bridge_still_unexported = (
        p1035.get("status")
        == "P1035_CURRENT_STRICT_QW2191_NADSOLITON_NEURAL_CHARACTER_SUPPORT_REFERENCE_TO_EXACT_SELECTOR_INTERFACE_BRIDGE_NONEXPORT_AUDITED"
        and theorem_result.get("support_reference_present_but_unbridged") is True
        and theorem_result.get("exact_support_reference_to_selector_interface_bridge_exported") is False
    )

    exact_selector_interface_still_unexported = (
        p1012.get("status")
        == "P1012_CURRENT_STRICT_QW2191_LAMBDA_BRANCH_INFO_PRIMARY_SCPC_LIKE_SELECTOR_INTERFACE_NONEXPORT_AUDITED"
        and p1012.get("exact_selector_interface_exported_on_current_repo_state") is False
    )

    exact_bridge_target_actual_realization_exported = len(scan_hits) > 0

    next_honest_move_requires_first_actual_realization_attempt = (
        bridge_target_packet_exported
        and bridge_still_unexported
        and exact_selector_interface_still_unexported
        and not exact_bridge_target_actual_realization_exported
    )

    add_check(
        "bridge_target_packet_exported",
        bridge_target_packet_exported,
        True,
        "F954 already exports the exact bridge target packet.",
    )
    add_check(
        "bridge_still_unexported",
        bridge_still_unexported,
        True,
        "P1035 still keeps the bridge unexported.",
    )
    add_check(
        "exact_selector_interface_still_unexported",
        exact_selector_interface_still_unexported,
        True,
        "P1012 still keeps the exact strict selector interface unexported.",
    )
    add_check(
        "exact_bridge_target_actual_realization_exported",
        exact_bridge_target_actual_realization_exported,
        False,
        "The repo still does not export any exact actual realization object for the frozen bridge target.",
    )
    add_check(
        "next_honest_move_requires_first_actual_realization_attempt",
        next_honest_move_requires_first_actual_realization_attempt,
        True,
        "Therefore the honest next move is to freeze one first exact actual-realization attempt rather than claim realized interface.",
    )

    status = (
        "P1036_CURRENT_STRICT_QW2191_NADSOLITON_NEURAL_CHARACTER_SUPPORT_REFERENCE_TO_EXACT_SELECTOR_INTERFACE_BRIDGE_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
        if not blocking and next_honest_move_requires_first_actual_realization_attempt
        else "FAIL_CURRENT_STRICT_QW2191_NADSOLITON_NEURAL_CHARACTER_SUPPORT_REFERENCE_TO_EXACT_SELECTOR_INTERFACE_BRIDGE_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDIT"
    )

    artifact = {
        "stage": "P1036",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "f954_bridge_target_packet_summary": rel(IN_F954),
            "p1035_bridge_nonexport_audit_probe": rel(IN_P1035),
            "p1012_selector_interface_nonexport_summary": rel(IN_P1012),
        },
        "exact_missing_actual_realization_candidate": {
            "candidate_id": "nadsoliton_neural_character_support_reference_to_exact_selector_interface_bridge_target_actual_realization_missing_v1",
            "bridge_target_ref": bridge_target.get("exact_bridge_target_id"),
            "support_reference_ref": bridge_target.get("support_reference_ref"),
            "supported_candidate_lane_ref": bridge_target.get("supported_candidate_lane_ref"),
            "exact_actual_realization_id": ACTUAL_REALIZATION_ID,
            "repo_scan_hits_for_exact_actual_realization": scan_hits,
        },
        "theorem_result": {
            "exact_bridge_target_frozen": bridge_target_packet_exported,
            "bridge_target_present_but_future_only": True,
            "exact_selector_interface_still_unexported": exact_selector_interface_still_unexported,
            "exact_bridge_target_actual_realization_exported": exact_bridge_target_actual_realization_exported,
            "next_honest_move_requires_first_actual_realization_attempt": next_honest_move_requires_first_actual_realization_attempt,
            "no_false_pass": True,
        },
        "checks": checks,
        "blocking_checks": blocking,
        "current_honest_reading": [
            "The repo now freezes one exact bridge target between the neural-character support-reference and the exact selector-interface question.",
            "But no exact actual realization yet exists for that bridge target.",
            "So the next honest move is to freeze the first exact actual-realization attempt, not to claim realized interface.",
        ],
        "recommended_next_packet": {
            "id": "T297_CURRENT_STRICT_QW2191_NADSOLITON_NEURAL_CHARACTER_SUPPORT_REFERENCE_TO_EXACT_SELECTOR_INTERFACE_BRIDGE_TARGET_ACTUAL_REALIZATION_ATTEMPT_SPEC",
            "goal": "Freeze the first exact actual-realization attempt for the bridge target from the neural-character support-reference to the exact strict selector-interface question.",
            "export_object_id": "nadsoliton_neural_character_support_reference_to_exact_selector_interface_bridge_target_actual_realization_attempt_v1",
        },
        "no_false_pass": True,
    }
    write_json(OUT_JSON, artifact)

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "exact_actual_realization_id": ACTUAL_REALIZATION_ID,
        "bridge_target_present_but_future_only": artifact["theorem_result"]["bridge_target_present_but_future_only"],
        "exact_selector_interface_still_unexported": artifact["theorem_result"]["exact_selector_interface_still_unexported"],
        "exact_bridge_target_actual_realization_exported": artifact["theorem_result"][
            "exact_bridge_target_actual_realization_exported"
        ],
        "next_honest_move_requires_first_actual_realization_attempt": artifact["theorem_result"][
            "next_honest_move_requires_first_actual_realization_attempt"
        ],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
