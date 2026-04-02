#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-23"
BRIDGE_TARGET_ID = (
    "nadsoliton_neural_character_support_reference_to_exact_selector_interface_bridge_target_v1"
)

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_F953 = GENERATED / "f953_current_strict_qw2191_nadsoliton_neural_character_information_primary_selector_support_reference_packet_summary.json"
IN_P1034 = GENERATED / "p1034_current_strict_qw2191_nadsoliton_neural_character_information_primary_selector_support_reference_admission_probe.json"
IN_P1012 = GENERATED / "p1012_current_strict_qw2191_lambda_branch_info_primary_scpc_like_selector_interface_nonexport_audit_probe_summary.json"

OUT_JSON = GENERATED / "p1035_current_strict_qw2191_nadsoliton_neural_character_support_reference_to_exact_selector_interface_bridge_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p1035_current_strict_qw2191_nadsoliton_neural_character_support_reference_to_exact_selector_interface_bridge_nonexport_audit_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def repo_scan_hits_for_exact_bridge_target() -> list[str]:
    hits: list[str] = []
    ignored_markers = (
        "P1035_CURRENT_STRICT_QW2191_NADSOLITON_NEURAL_CHARACTER_SUPPORT_REFERENCE_TO_EXACT_SELECTOR_INTERFACE_BRIDGE_NONEXPORT_AUDIT_PROBE",
        "N868_CURRENT_STRICT_QW2191_NADSOLITON_NEURAL_CHARACTER_SUPPORT_REFERENCE_TO_EXACT_SELECTOR_INTERFACE_BRIDGE_NONEXPORT_AUDIT_THEOREM",
        "F954_CURRENT_STRICT_QW2191_NADSOLITON_NEURAL_CHARACTER_SUPPORT_REFERENCE_TO_EXACT_SELECTOR_INTERFACE_BRIDGE_TARGET_PACKET",
        "p1035_current_strict_qw2191_nadsoliton_neural_character_support_reference_to_exact_selector_interface_bridge_nonexport_audit_probe",
        "f954_current_strict_qw2191_nadsoliton_neural_character_support_reference_to_exact_selector_interface_bridge_target_packet",
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
        if BRIDGE_TARGET_ID in text:
            hits.append(rel_path)
    return sorted(hits)


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_F953, IN_P1034, IN_P1012]
    missing = [rel(path) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1035",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    f953 = load_json(IN_F953)
    p1034 = load_json(IN_P1034)
    p1012 = load_json(IN_P1012)

    theorem_result = p1034.get("theorem_result") or {}
    scan_hits = repo_scan_hits_for_exact_bridge_target()

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

    neural_character_support_reference_exported = (
        f953.get("status")
        == "F953_EXECUTED_CURRENT_STRICT_QW2191_NADSOLITON_NEURAL_CHARACTER_INFORMATION_PRIMARY_SELECTOR_SUPPORT_REFERENCE_PACKET_NO_FALSE_PASS"
        and f953.get("exported_object_id")
        == "nadsoliton_neural_character_information_primary_selector_support_reference_v1"
    )

    support_reference_is_reference_only = (
        p1034.get("status")
        == "P1034_CURRENT_STRICT_QW2191_NADSOLITON_NEURAL_CHARACTER_INFORMATION_PRIMARY_SELECTOR_SUPPORT_REFERENCE_ADMITTED_INTERFACE_STILL_BLOCKED"
        and theorem_result.get("nadsoliton_neural_character_support_reference_admitted") is True
        and theorem_result.get("support_reference_grade") == "cross_repo_support_reference_only"
        and theorem_result.get("strict_selector_interface_exported") is False
        and theorem_result.get("strict_selector_source_exported") is False
    )

    exact_selector_interface_still_unexported = (
        p1012.get("status")
        == "P1012_CURRENT_STRICT_QW2191_LAMBDA_BRANCH_INFO_PRIMARY_SCPC_LIKE_SELECTOR_INTERFACE_NONEXPORT_AUDITED"
        and p1012.get("exact_selector_interface_exported_on_current_repo_state") is False
    )

    exact_bridge_target_already_exists_on_current_repo_state = len(scan_hits) > 0
    exact_bridge_to_selector_interface_exported = False

    next_honest_move_requires_freezing_exact_bridge_target = (
        neural_character_support_reference_exported
        and support_reference_is_reference_only
        and exact_selector_interface_still_unexported
        and not exact_bridge_target_already_exists_on_current_repo_state
        and not exact_bridge_to_selector_interface_exported
    )

    add_check(
        "neural_character_support_reference_exported",
        neural_character_support_reference_exported,
        True,
        "F953 already exports the neural-character support-reference packet.",
    )
    add_check(
        "support_reference_is_reference_only",
        support_reference_is_reference_only,
        True,
        "P1034 keeps that support-reference only at cross-repo support-reference grade.",
    )
    add_check(
        "exact_selector_interface_still_unexported",
        exact_selector_interface_still_unexported,
        True,
        "P1012 still keeps the exact strict selector interface unexported.",
    )
    add_check(
        "exact_bridge_target_already_exists_on_current_repo_state",
        exact_bridge_target_already_exists_on_current_repo_state,
        False,
        "The repo still does not export any exact support-reference-to-selector-interface bridge target under the exact bridge identifier.",
    )
    add_check(
        "next_honest_move_requires_freezing_exact_bridge_target",
        next_honest_move_requires_freezing_exact_bridge_target,
        True,
        "Therefore the honest next move is to freeze one exact missing bridge target rather than claim interface realization.",
    )

    status = (
        "P1035_CURRENT_STRICT_QW2191_NADSOLITON_NEURAL_CHARACTER_SUPPORT_REFERENCE_TO_EXACT_SELECTOR_INTERFACE_BRIDGE_NONEXPORT_AUDITED"
        if not blocking and next_honest_move_requires_freezing_exact_bridge_target
        else "FAIL_CURRENT_STRICT_QW2191_NADSOLITON_NEURAL_CHARACTER_SUPPORT_REFERENCE_TO_EXACT_SELECTOR_INTERFACE_BRIDGE_NONEXPORT_AUDIT"
    )

    artifact = {
        "stage": "P1035",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "f953_support_reference_packet_summary": rel(IN_F953),
            "p1034_support_reference_admission_probe": rel(IN_P1034),
            "p1012_selector_interface_nonexport_summary": rel(IN_P1012),
        },
        "exact_missing_bridge_target_candidate": {
            "candidate_id": "nadsoliton_neural_character_support_reference_to_exact_selector_interface_bridge_target_missing_v1",
            "support_reference_ref": theorem_result.get("support_reference_id"),
            "supported_candidate_lane_ref": theorem_result.get("supported_candidate_lane_ref"),
            "exact_bridge_target_id": BRIDGE_TARGET_ID,
            "repo_scan_hits_for_exact_bridge_target": scan_hits,
        },
        "theorem_result": {
            "exact_support_reference_frozen": neural_character_support_reference_exported,
            "support_reference_present_but_unbridged": True,
            "exact_selector_interface_still_unexported": exact_selector_interface_still_unexported,
            "exact_support_reference_to_selector_interface_bridge_exported": exact_bridge_to_selector_interface_exported,
            "next_honest_move_requires_freezing_exact_bridge_target": next_honest_move_requires_freezing_exact_bridge_target,
            "no_false_pass": True,
        },
        "checks": checks,
        "blocking_checks": blocking,
        "current_honest_reading": [
            "The repo now freezes one lawful neural-character support-reference and one active information-primary candidate lane.",
            "But no exact bridge yet connects that support-reference to one exact strict selector-interface question.",
            "So the next honest move is to freeze the exact missing bridge target, not to claim interface realization.",
        ],
        "recommended_next_packet": {
            "id": "F954_CURRENT_STRICT_QW2191_NADSOLITON_NEURAL_CHARACTER_SUPPORT_REFERENCE_TO_EXACT_SELECTOR_INTERFACE_BRIDGE_TARGET_PACKET",
            "goal": "Freeze the exact missing bridge target from the neural-character support-reference into the exact strict selector-interface question for the active candidate lane.",
            "export_object_id": BRIDGE_TARGET_ID,
        },
        "no_false_pass": True,
    }
    write_json(OUT_JSON, artifact)

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "exact_bridge_target_id": BRIDGE_TARGET_ID,
        "support_reference_present_but_unbridged": artifact["theorem_result"]["support_reference_present_but_unbridged"],
        "exact_selector_interface_still_unexported": artifact["theorem_result"]["exact_selector_interface_still_unexported"],
        "exact_support_reference_to_selector_interface_bridge_exported": artifact["theorem_result"][
            "exact_support_reference_to_selector_interface_bridge_exported"
        ],
        "next_honest_move_requires_freezing_exact_bridge_target": artifact["theorem_result"][
            "next_honest_move_requires_freezing_exact_bridge_target"
        ],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
