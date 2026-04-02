#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_F828 = GENERATED / "f828_current_strict_alpha_s_different_required_form_statement_provider_class_shift_candidate_reference_packet.json"
IN_F825 = GENERATED / "f825_current_strict_alpha_s_exact_required_form_statement_target_packet.json"
IN_P759 = GENERATED / "p759_current_strict_t213_pair12_source_side_branch_selection_provider_target_probe_summary.json"
IN_P762 = GENERATED / "p762_current_strict_t216_pair12_source_side_branch_selection_provider_actual_realization_attempt_probe_summary.json"
IN_P763 = GENERATED / "p763_current_strict_t217_pair12_source_side_branch_selection_provider_actual_realization_attempt_immediate_missing_interface_nonexport_audit_probe_summary.json"
IN_P764 = GENERATED / "p764_current_strict_t218_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_target_probe_summary.json"

OUT = GENERATED / "p829_current_strict_alpha_s_no_exact_pair12_source_side_branch_selection_provider_shift_to_exact_required_form_statement_interface_target_exported_target_freeze_required.json"
OUT_SUMMARY = GENERATED / "p829_current_strict_alpha_s_no_exact_pair12_source_side_branch_selection_provider_shift_to_exact_required_form_statement_interface_target_exported_target_freeze_required_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def normalize(text: str) -> str:
    return " ".join(
        text.lower()
        .replace("“", '"')
        .replace("”", '"')
        .replace("’", "'")
        .replace("‑", "-")
        .replace("–", "-")
        .replace("—", "-")
        .replace("->", " ")
        .replace("→", " ")
        .replace("/", " ")
        .replace("-", " ")
        .replace("_", " ")
        .split()
    )


def find_repo_hits(candidate_id: str, downstream_id: str, needles: list[str]) -> list[str]:
    hits: list[str] = []
    excluded_suffixes = [
        "generated/p828_current_strict_alpha_s_different_required_form_statement_provider_class_shift_candidate_reference_lane_admission_probe.json",
        "generated/p828_current_strict_alpha_s_different_required_form_statement_provider_class_shift_candidate_reference_lane_admission_probe_summary.json",
        "generated/f828_current_strict_alpha_s_different_required_form_statement_provider_class_shift_candidate_reference_packet.json",
        "generated/f828_current_strict_alpha_s_different_required_form_statement_provider_class_shift_candidate_reference_packet_summary.json",
        "P828_CURRENT_STRICT_ALPHA_S_DIFFERENT_REQUIRED_FORM_STATEMENT_PROVIDER_CLASS_SHIFT_CANDIDATE_REFERENCE_LANE_ADMISSION_PROBE.md",
        "F828_CURRENT_STRICT_ALPHA_S_DIFFERENT_REQUIRED_FORM_STATEMENT_PROVIDER_CLASS_SHIFT_CANDIDATE_REFERENCE_PACKET.md",
        "p828_current_strict_alpha_s_different_required_form_statement_provider_class_shift_candidate_reference_lane_admission_probe.py",
        "f828_current_strict_alpha_s_different_required_form_statement_provider_class_shift_candidate_reference_packet.py",
        "generated/p829_current_strict_alpha_s_no_exact_pair12_source_side_branch_selection_provider_shift_to_exact_required_form_statement_interface_target_exported_target_freeze_required.json",
        "generated/p829_current_strict_alpha_s_no_exact_pair12_source_side_branch_selection_provider_shift_to_exact_required_form_statement_interface_target_exported_target_freeze_required_summary.json",
        "generated/f829_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_exact_required_form_statement_interface_target_packet.json",
        "generated/f829_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_exact_required_form_statement_interface_target_packet_summary.json",
        "P829_CURRENT_STRICT_ALPHA_S_NO_EXACT_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_SHIFT_TO_EXACT_REQUIRED_FORM_STATEMENT_INTERFACE_TARGET_EXPORTED_TARGET_FREEZE_REQUIRED.md",
        "F829_CURRENT_STRICT_ALPHA_S_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_SHIFT_TO_EXACT_REQUIRED_FORM_STATEMENT_INTERFACE_TARGET_PACKET.md",
        "p829_current_strict_alpha_s_no_exact_pair12_source_side_branch_selection_provider_shift_to_exact_required_form_statement_interface_target_exported_target_freeze_required.py",
        "f829_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_exact_required_form_statement_interface_target_packet.py",
    ]
    norm_candidate = normalize(candidate_id)
    norm_downstream = normalize(downstream_id)
    norm_needles = [normalize(needle) for needle in needles]
    for path in ROOT.rglob("*"):
        if not path.is_file():
            continue
        if path.suffix not in {".md", ".json", ".py"}:
            continue
        relpath = rel(path)
        if any(relpath.endswith(excluded) for excluded in excluded_suffixes):
            continue
        try:
            text = path.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            continue
        hay = normalize(text)
        if norm_candidate in hay and norm_downstream in hay and all(needle in hay for needle in norm_needles):
            hits.append(relpath)
    return sorted(hits)


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_F828, IN_F825, IN_P759, IN_P762, IN_P763, IN_P764]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P829",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    f828 = load_json(IN_F828)
    f825 = load_json(IN_F825)
    p759 = load_json(IN_P759)
    p762 = load_json(IN_P762)
    p763 = load_json(IN_P763)
    p764 = load_json(IN_P764)

    candidate_lane = f828.get("exported_object") or {}
    exact_target = f825.get("target_object") or {}

    interface_target_hits = find_repo_hits(
        "alpha_s_pair12_source_side_branch_selection_provider_exact_required_form_statement_shift_candidate_reference_lane_v1",
        "alpha_s_pair12_source_side_branch_selection_provider_lawful_schema_domain_admission_exact_required_form_statement_target_v1",
        ["interface target"],
    )
    adapter_or_rule_hits = find_repo_hits(
        "alpha_s_pair12_source_side_branch_selection_provider_exact_required_form_statement_shift_candidate_reference_lane_v1",
        "alpha_s_pair12_source_side_branch_selection_provider_lawful_schema_domain_admission_exact_required_form_statement_target_v1",
        ["adapter"],
    ) + find_repo_hits(
        "alpha_s_pair12_source_side_branch_selection_provider_exact_required_form_statement_shift_candidate_reference_lane_v1",
        "alpha_s_pair12_source_side_branch_selection_provider_lawful_schema_domain_admission_exact_required_form_statement_target_v1",
        ["carrier identification"],
    )
    adapter_or_rule_hits = sorted(set(adapter_or_rule_hits))

    checks = [
        {
            "id": "f828_exports_only_candidate_reference_lane_with_blocked_alpha_s_shift_interface",
            "pass": (
                f828.get("status")
                == "F828_EXECUTED_CURRENT_STRICT_ALPHA_S_DIFFERENT_REQUIRED_FORM_STATEMENT_PROVIDER_CLASS_SHIFT_CANDIDATE_REFERENCE_PACKET_NO_FALSE_PASS"
                and candidate_lane.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_exact_required_form_statement_shift_candidate_reference_lane_v1"
                and candidate_lane.get("candidate_reference_lane_grade") == "reference_context_candidate_only"
                and candidate_lane.get("alpha_s_shift_interface_status") == "blocked_nonexport"
                and candidate_lane.get("exact_required_form_statement_status") == "blocked_nonexport"
                and candidate_lane.get("provider_class_shift_realization_status") == "not_realized"
            ),
            "details": "F828 exports the T213/T216 lane only as a candidate reference lane and keeps the alpha_s-side shift interface blocked.",
        },
        {
            "id": "f825_already_freezes_the_exact_downstream_exact_required_form_statement_target",
            "pass": (
                f825.get("status")
                == "F825_EXECUTED_CURRENT_STRICT_ALPHA_S_EXACT_REQUIRED_FORM_STATEMENT_TARGET_PACKET_NO_FALSE_PASS"
                and exact_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_schema_domain_admission_exact_required_form_statement_target_v1"
                and any(
                    isinstance(item, dict) and item.get("name") == "exact_required_form_statement_ref"
                    for item in (exact_target.get("required_fields") or [])
                )
            ),
            "details": "F825 already freezes the exact downstream exact-required-form-statement problem that any lawful shift interface would have to reach.",
        },
        {
            "id": "t213_t216_lane_still_exports_only_own_lane_target_attempt_and_missing_interface",
            "pass": (
                p759.get("status") == "PASS_STRICT_T213_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_TARGET_EXPORTED"
                and p759.get("current_t213_target_is_future_route_only") is True
                and p762.get("status")
                == "PASS_STRICT_T216_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_EXPORTED"
                and p762.get("first_actual_t213_realization_attempt_keeps_success_failure_open") is True
                and p763.get("status")
                == "PASS_STRICT_T217_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_IMMEDIATE_MISSING_INTERFACE_NONEXPORT_AUDITED"
                and p763.get("current_t216_attempt_immediate_missing_interface_is_still_unexported") is True
                and p764.get("status")
                == "PASS_STRICT_T218_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_TARGET_EXPORTED"
                and p764.get("current_t218_target_is_future_route_only") is True
            ),
            "details": "The T213/T216 lane remains self-contained: own future-only target, own attempt, and own missing interface only.",
        },
        {
            "id": "no_exact_shift_to_exact_required_form_statement_interface_target_exported_for_t213_t216_lane",
            "pass": interface_target_hits == [],
            "details": "No current export names one exact T213/T216 -> F825 shift-to-exact-required-form-statement interface target.",
        },
        {
            "id": "no_exact_adapter_or_carrier_rule_exported_for_that_missing_interface",
            "pass": adapter_or_rule_hits == [],
            "details": "No current export names one exact adapter or carrier-identification rule for the still-missing T213/T216 -> F825 interface.",
        },
    ]

    failed_checks = [item["id"] for item in checks if not item["pass"]]
    all_pass = all(item["pass"] for item in checks)

    theorem_result = {
        "exact_pair12_source_side_branch_selection_provider_shift_to_exact_required_form_statement_interface_target_exported": False
        if all_pass
        else None,
        "exact_pair12_source_side_branch_selection_provider_shift_to_exact_required_form_statement_adapter_or_carrier_rule_exported": False
        if all_pass
        else None,
        "next_honest_move_requires_freezing_exact_shift_to_exact_required_form_statement_interface_target": all_pass,
        "no_false_pass": True,
    }

    status = (
        "P829_CURRENT_STRICT_ALPHA_S_NO_EXACT_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_SHIFT_TO_EXACT_REQUIRED_FORM_STATEMENT_INTERFACE_TARGET_EXPORTED_TARGET_FREEZE_REQUIRED"
        if all_pass
        else "P829_REQUIRES_REVIEW"
    )

    artifact = {
        "stage": "P829",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "f828_different_required_form_statement_provider_class_shift_candidate_reference_packet": rel(IN_F828),
            "f825_downstream_exact_required_form_statement_target_packet": rel(IN_F825),
            "p759_t213_target_probe_summary": rel(IN_P759),
            "p762_t216_attempt_probe_summary": rel(IN_P762),
            "p763_t217_missing_interface_boundary_summary": rel(IN_P763),
            "p764_t218_missing_interface_target_summary": rel(IN_P764),
        },
        "exact_missing_interface_target_candidate": {
            "candidate_id": "alpha_s_pair12_source_side_branch_selection_provider_shift_to_exact_required_form_statement_interface_target_missing_v1",
            "different_provider_class_shift_candidate_reference_lane_ref": candidate_lane.get("object_id"),
            "downstream_exact_required_form_statement_target_ref": exact_target.get("object_id"),
            "repo_scan_hits_for_exact_interface_target": interface_target_hits,
            "repo_scan_hits_for_exact_adapter_or_carrier_rule": adapter_or_rule_hits,
        },
        "theorem_result": theorem_result,
        "checks": checks,
        "failed_checks": failed_checks,
        "current_honest_reading": [
            "The repo now exports a real different-provider-class candidate lane for the alpha_s exact-required-form-statement problem.",
            "But that lane still exports only its own-lane target, attempt, and missing interface, not one exact interface object into the current F825 downstream problem.",
            "So the next honest move is to freeze the exact shift-to-exact-required-form-statement interface target itself before claiming any adapter or carrier-safe entry rule.",
        ],
        "recommended_next_packet": {
            "id": "F829_CURRENT_STRICT_ALPHA_S_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_SHIFT_TO_EXACT_REQUIRED_FORM_STATEMENT_INTERFACE_TARGET_PACKET",
            "goal": "Freeze the exact missing T213/T216 -> F825 shift-to-exact-required-form-statement interface target without claiming realization.",
            "export_object_id": "alpha_s_pair12_source_side_branch_selection_provider_shift_to_exact_required_form_statement_interface_target_v1",
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "P829",
        "status": status,
        "as_of": AS_OF,
        "exact_pair12_source_side_branch_selection_provider_shift_to_exact_required_form_statement_interface_target_exported": theorem_result[
            "exact_pair12_source_side_branch_selection_provider_shift_to_exact_required_form_statement_interface_target_exported"
        ],
        "exact_pair12_source_side_branch_selection_provider_shift_to_exact_required_form_statement_adapter_or_carrier_rule_exported": theorem_result[
            "exact_pair12_source_side_branch_selection_provider_shift_to_exact_required_form_statement_adapter_or_carrier_rule_exported"
        ],
        "recommended_next_packet": artifact["recommended_next_packet"]["id"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
