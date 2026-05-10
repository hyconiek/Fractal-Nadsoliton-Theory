#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P825 = GENERATED / "p825_current_strict_alpha_s_required_form_statement_class_candidate_supported_exact_required_form_statement_blocked.json"
IN_P826 = GENERATED / "p826_current_strict_alpha_s_exact_required_form_statement_same_lane_exhaustion_boundary_audit_probe.json"
IN_F825 = GENERATED / "f825_current_strict_alpha_s_exact_required_form_statement_target_packet.json"
IN_F826 = GENERATED / "f826_current_strict_alpha_s_exact_required_form_statement_continuation_boundary_packet.json"

OUT = GENERATED / "p827_current_strict_alpha_s_neighboring_statement_form_stack_reuse_as_same_lane_exact_required_form_statement_source_negative_provider_shift_required.json"
OUT_SUMMARY = GENERATED / "p827_current_strict_alpha_s_neighboring_statement_form_stack_reuse_as_same_lane_exact_required_form_statement_source_negative_provider_shift_required_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def has_required_field(obj: dict[str, Any], name: str) -> bool:
    return any(
        isinstance(item, dict) and item.get("name") == name
        for item in (obj.get("required_fields") or [])
    )


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P825, IN_P826, IN_F825, IN_F826]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P827",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p825 = load_json(IN_P825)
    p826 = load_json(IN_P826)
    f825 = load_json(IN_F825)
    f826 = load_json(IN_F826)

    p825_theorem = p825.get("theorem_result") or {}
    p826_theorem = p826.get("theorem_result") or {}
    f825_target = f825.get("target_object") or {}
    f826_export = f826.get("exported_object") or {}

    support_refs = (f825.get("target_refs") or {}).get("required_form_statement_class_candidate_support_refs") or []
    slot_refs = (f825.get("target_refs") or {}).get("neighboring_statement_or_form_slot_refs") or []
    hard_limits = f825.get("hard_limits") or []

    checks = [
        {
            "id": "f826_already_requires_new_same_lane_source_or_provider_shift",
            "pass": (
                f826.get("status")
                == "F826_EXECUTED_CURRENT_STRICT_ALPHA_S_EXACT_REQUIRED_FORM_STATEMENT_CONTINUATION_BOUNDARY_PACKET_NO_FALSE_PASS"
                and f826_export.get("object_id") == "alpha_s_exact_required_form_statement_continuation_boundary_v1"
                and f826_export.get("admitted_next_move_classes")
                == [
                    "export_one_genuinely_new_same_lane_exact_required_form_statement_source",
                    "shift_to_a_different_required_form_statement_provider_class_lane",
                ]
            ),
            "details": "F826 already narrows continuation to one genuinely new same-lane exact required-form statement source or a provider-class shift.",
        },
        {
            "id": "f825_exact_required_form_statement_target_is_real",
            "pass": (
                f825.get("status")
                == "F825_EXECUTED_CURRENT_STRICT_ALPHA_S_EXACT_REQUIRED_FORM_STATEMENT_TARGET_PACKET_NO_FALSE_PASS"
                and f825_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_schema_domain_admission_exact_required_form_statement_target_v1"
                and has_required_field(f825_target, "exact_required_form_statement_ref")
            ),
            "details": "F825 exports a real exact required-form statement target object for the current lane.",
        },
        {
            "id": "p825_keeps_exact_statement_unexported_and_support_only_candidate_grade",
            "pass": (
                p825.get("status")
                == "P825_CURRENT_STRICT_ALPHA_S_REQUIRED_FORM_STATEMENT_CLASS_CANDIDATE_SUPPORTED_EXACT_REQUIRED_FORM_STATEMENT_BLOCKED"
                and p825_theorem.get("required_form_statement_class_candidate_supported_now") is True
                and p825_theorem.get("exact_required_form_statement_exported_now") is False
            ),
            "details": "P825 already keeps the exact required-form statement unexported while preserving only candidate-grade neighboring support.",
        },
        {
            "id": "current_stack_is_explicitly_support_or_slot_scaffolding",
            "pass": (
                has_required_field(f825_target, "required_form_statement_class_candidate_support_refs")
                and has_required_field(f825_target, "neighboring_statement_or_form_slot_refs")
                and len(support_refs) >= 9
                and len(slot_refs) >= 8
                and any("silently discharges the new lane" in item for item in hard_limits)
            ),
            "details": "The current stack is explicitly exported only as neighboring support and slot scaffolding, not as an exact source object.",
        },
        {
            "id": "same_lane_repetition_is_already_closed_below_current_blocker",
            "pass": (
                p826.get("status")
                == "P826_CURRENT_STRICT_ALPHA_S_EXACT_REQUIRED_FORM_STATEMENT_SAME_LANE_EXHAUSTION_BOUNDARY_AUDIT_PROBE"
                and p826_theorem.get("same_lane_passive_groundwork_exhausted") is True
                and p826_theorem.get("next_honest_move_requires_continuation_boundary_export") is True
            ),
            "details": "P826 already closes further same-lane passive repetition below the current blocker.",
        },
        {
            "id": "no_current_neighboring_stack_object_counts_as_genuinely_new_same_lane_source",
            "pass": (
                f826_export.get("forbidden_repetition_clause")
                == "no_further_same_level_passive_split_under_unchanged_exact_required_form_statement_ref_as_primary_move"
                and p825_theorem.get("exact_required_form_statement_exported_now") is False
                and has_required_field(f825_target, "exact_required_form_statement_ref")
                and not has_required_field(f825_target, "same_lane_exact_required_form_statement_source_ref")
            ),
            "details": "Nothing currently exported upgrades the neighboring statement/form stack into a genuinely new same-lane exact required-form statement source for alpha_s.",
        },
    ]

    failed_checks = [item["id"] for item in checks if not item["pass"]]
    all_pass = all(item["pass"] for item in checks)

    theorem_result = {
        "neighboring_statement_form_stack_exists": checks[3]["pass"],
        "neighboring_statement_form_stack_is_same_lane_exact_required_form_statement_source_for_alpha_s": False
        if all_pass
        else None,
        "current_repo_exports_no_genuinely_new_same_lane_exact_required_form_statement_source_for_alpha_s": all_pass,
        "next_honest_move_requires_required_form_statement_provider_class_shift": all_pass,
        "no_false_pass": True,
    }

    status = (
        "P827_CURRENT_STRICT_ALPHA_S_NEIGHBORING_STATEMENT_FORM_STACK_REUSE_AS_SAME_LANE_EXACT_REQUIRED_FORM_STATEMENT_SOURCE_NEGATIVE_PROVIDER_SHIFT_REQUIRED"
        if all_pass
        else "P827_REQUIRES_REVIEW"
    )

    artifact = {
        "stage": "P827",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p825_required_form_statement_probe": rel(IN_P825),
            "p826_same_lane_exhaustion_boundary_probe": rel(IN_P826),
            "f825_exact_required_form_statement_target_packet": rel(IN_F825),
            "f826_exact_required_form_statement_continuation_boundary_packet": rel(IN_F826),
        },
        "rejected_same_lane_source_candidate_class": {
            "candidate_id": "neighboring_statement_form_stack_reuse_same_lane_exact_required_form_statement_source_candidate_v1",
            "current_exact_required_form_statement_target_ref": f825_target.get("object_id"),
            "continuation_boundary_ref": f826_export.get("object_id"),
            "rejection_basis_refs": [
                rel(IN_P825),
                rel(IN_P826),
                rel(IN_F825),
                rel(IN_F826),
            ],
        },
        "theorem_result": theorem_result,
        "checks": checks,
        "failed_checks": failed_checks,
        "current_honest_reading": [
            "Neighboring statement/form supports and slot scaffolding do exist on the current repo state.",
            "But those artifacts remain support-only and slot-only; they are not currently exported as a genuinely new same-lane exact required-form statement source for the alpha_s lane.",
            "So the remaining honest continuation class is provider-class shift, not same-lane verbal promotion of neighboring scaffolding.",
        ],
        "recommended_next_packet": {
            "id": "F827_CURRENT_STRICT_ALPHA_S_EXACT_REQUIRED_FORM_STATEMENT_PROVIDER_CLASS_SHIFT_REQUIREMENT_PACKET",
            "goal": "Freeze the current-repo-state requirement that continuation now proceed by provider-class shift rather than same-lane verbal promotion of neighboring statement/form scaffolding.",
            "export_object_id": "alpha_s_exact_required_form_statement_provider_class_shift_requirement_v1",
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "P827",
        "status": status,
        "as_of": AS_OF,
        "current_repo_exports_no_genuinely_new_same_lane_exact_required_form_statement_source_for_alpha_s": theorem_result[
            "current_repo_exports_no_genuinely_new_same_lane_exact_required_form_statement_source_for_alpha_s"
        ],
        "next_honest_move_requires_required_form_statement_provider_class_shift": theorem_result[
            "next_honest_move_requires_required_form_statement_provider_class_shift"
        ],
        "recommended_next_packet": artifact["recommended_next_packet"]["id"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
