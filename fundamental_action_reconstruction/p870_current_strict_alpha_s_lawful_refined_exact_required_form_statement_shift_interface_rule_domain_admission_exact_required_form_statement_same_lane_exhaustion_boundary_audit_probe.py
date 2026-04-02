#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-20"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P867 = GENERATED / "p867_current_strict_alpha_s_output_schema_statement_class_candidate_supported_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_exact_output_schema_statement_blocked.json"
IN_P868 = GENERATED / "p868_current_strict_alpha_s_statement_form_class_candidate_supported_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_exact_statement_required_form_blocked.json"
IN_P869 = GENERATED / "p869_current_strict_alpha_s_required_form_statement_class_candidate_supported_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_exact_required_form_statement_blocked.json"
IN_F867 = GENERATED / "f867_current_strict_alpha_s_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_exact_output_schema_statement_target_packet.json"
IN_F868 = GENERATED / "f868_current_strict_alpha_s_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_exact_statement_required_form_target_packet.json"
IN_F869 = GENERATED / "f869_current_strict_alpha_s_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_exact_required_form_statement_target_packet.json"
IN_S2 = ROOT / "S2_CURRENT_FAR_STRATEGIC_PRIORITY_REORIENTATION_PACKET.md"

OUT = GENERATED / "p870_current_strict_alpha_s_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_exact_required_form_statement_same_lane_exhaustion_boundary_audit_probe.json"
OUT_SUMMARY = GENERATED / "p870_current_strict_alpha_s_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_exact_required_form_statement_same_lane_exhaustion_boundary_audit_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def has_required_field(obj: dict[str, Any], name: str) -> bool:
    return any(
        isinstance(item, dict) and item.get("name") == name
        for item in (obj.get("required_fields") or [])
    )


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P867, IN_P868, IN_P869, IN_F867, IN_F868, IN_F869, IN_S2]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P870",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p867 = load_json(IN_P867)
    p868 = load_json(IN_P868)
    p869 = load_json(IN_P869)
    f867 = load_json(IN_F867)
    f868 = load_json(IN_F868)
    f869 = load_json(IN_F869)
    s2_text = load_text(IN_S2)

    p867_theorem = p867.get("theorem_result") or {}
    p868_theorem = p868.get("theorem_result") or {}
    p869_theorem = p869.get("theorem_result") or {}
    f867_target = f867.get("target_object") or {}
    f868_target = f868.get("target_object") or {}
    f869_target = f869.get("target_object") or {}

    support_refs = (f869.get("target_refs") or {}).get("required_form_statement_class_candidate_support_refs") or []
    slot_refs = (f869.get("target_refs") or {}).get("neighboring_statement_or_form_slot_refs") or []

    lane_token_hits: list[str] = []
    lane_export_token = "lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_exact_required_form_statement"
    for path in sorted(GENERATED.glob("*.json")):
        name = path.name
        if name.startswith(("p867_", "f867_", "p868_", "f868_", "p869_", "f869_", "p870_", "f870_")):
            continue
        text = path.read_text(encoding="utf-8")
        if lane_export_token in text:
            lane_token_hits.append(rel(path))

    checks = [
        {
            "id": "f867_refined_exact_output_schema_statement_target_already_frozen",
            "pass": (
                f867.get("status")
                == "F867_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_EXACT_REQUIRED_FORM_STATEMENT_SHIFT_INTERFACE_RULE_DOMAIN_ADMISSION_EXACT_OUTPUT_SCHEMA_STATEMENT_TARGET_PACKET_NO_FALSE_PASS"
                and f867_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_exact_output_schema_statement_target_v1"
            ),
            "details": "F867 already freezes the lawful refined shift-interface-rule exact output-schema statement target for the current lane.",
        },
        {
            "id": "f868_refined_exact_statement_required_form_target_already_frozen",
            "pass": (
                f868.get("status")
                == "F868_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_EXACT_REQUIRED_FORM_STATEMENT_SHIFT_INTERFACE_RULE_DOMAIN_ADMISSION_EXACT_STATEMENT_REQUIRED_FORM_TARGET_PACKET_NO_FALSE_PASS"
                and f868_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_exact_statement_required_form_target_v1"
                and has_required_field(f868_target, "exact_required_form_statement_ref")
            ),
            "details": "F868 already freezes the lawful refined shift-interface-rule exact statement-required-form target and names exact_required_form_statement_ref as one exact missing field.",
        },
        {
            "id": "f869_refined_exact_required_form_statement_target_already_frozen",
            "pass": (
                f869.get("status")
                == "F869_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_EXACT_REQUIRED_FORM_STATEMENT_SHIFT_INTERFACE_RULE_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_TARGET_PACKET_NO_FALSE_PASS"
                and f869_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_exact_required_form_statement_target_v1"
                and has_required_field(f869_target, "exact_required_form_statement_ref")
            ),
            "details": "F869 already freezes the lawful refined shift-interface-rule exact required-form statement target itself.",
        },
        {
            "id": "p869_keeps_refined_exact_required_form_statement_unexported",
            "pass": (
                p869.get("status")
                == "P869_CURRENT_STRICT_ALPHA_S_REQUIRED_FORM_STATEMENT_CLASS_CANDIDATE_SUPPORTED_LAWFUL_REFINED_EXACT_REQUIRED_FORM_STATEMENT_SHIFT_INTERFACE_RULE_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_BLOCKED"
                and p869_theorem.get("exact_required_form_statement_exported_now") is False
                and p869_theorem.get("sharp_blocker_field") == "exact_required_form_statement_ref"
            ),
            "details": "P869 already keeps the lawful refined shift-interface-rule exact required-form statement unexported on the current repo state.",
        },
        {
            "id": "f869_already_packs_neighboring_support_and_slot_context",
            "pass": (
                len(support_refs) >= 7
                and len(slot_refs) >= 7
                and has_required_field(f869_target, "required_form_statement_class_candidate_support_refs")
                and has_required_field(f869_target, "neighboring_statement_or_form_slot_refs")
            ),
            "details": "F869 already packs the neighboring statement/form-class support refs and neighboring slot refs into one exact lawful refined target context.",
        },
        {
            "id": "repo_scan_finds_no_current_lane_exact_required_form_statement_export",
            "pass": lane_token_hits == [],
            "details": "Repo scan finds no generated artifact exporting the current lawful refined shift-interface-rule exact required-form statement outside the new frozen F869 lineage itself.",
        },
        {
            "id": "blocker_chain_has_already_descended_to_current_floor",
            "pass": (
                p867.get("status")
                == "P867_CURRENT_STRICT_ALPHA_S_OUTPUT_SCHEMA_STATEMENT_CLASS_CANDIDATE_SUPPORTED_LAWFUL_REFINED_EXACT_REQUIRED_FORM_STATEMENT_SHIFT_INTERFACE_RULE_DOMAIN_ADMISSION_EXACT_OUTPUT_SCHEMA_STATEMENT_BLOCKED"
                and p867_theorem.get("sharp_blocker_field") == "exact_output_schema_statement"
                and p868.get("status")
                == "P868_CURRENT_STRICT_ALPHA_S_STATEMENT_FORM_CLASS_CANDIDATE_SUPPORTED_LAWFUL_REFINED_EXACT_REQUIRED_FORM_STATEMENT_SHIFT_INTERFACE_RULE_DOMAIN_ADMISSION_EXACT_STATEMENT_REQUIRED_FORM_BLOCKED"
                and p868_theorem.get("sharp_blocker_field") == "exact_statement_required_form_ref"
                and p869.get("status")
                == "P869_CURRENT_STRICT_ALPHA_S_REQUIRED_FORM_STATEMENT_CLASS_CANDIDATE_SUPPORTED_LAWFUL_REFINED_EXACT_REQUIRED_FORM_STATEMENT_SHIFT_INTERFACE_RULE_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_BLOCKED"
                and p869_theorem.get("sharp_blocker_field") == "exact_required_form_statement_ref"
            ),
            "details": "The lawful refined shift-interface-rule blocker chain has already been pushed down from exact output-schema statement to exact statement-required form to exact required-form statement, so the current lane is already at its present passive floor.",
        },
        {
            "id": "s2_noncyclic_continuation_discipline_applies",
            "pass": (
                "strict-core ToE closure using only strict-side sources" in s2_text
                and "new provider class and noncyclic anchor, not a repetition of L5/L12." in s2_text
            ),
            "details": "S2 still enforces noncyclic continuation discipline against repeated same-blocker passive splitting.",
        },
        {
            "id": "no_residual_passive_same_lane_loophole_below_current_blocker",
            "pass": (
                p869_theorem.get("required_form_statement_class_candidate_supported_now") is True
                and p869_theorem.get("exact_required_form_statement_exported_now") is False
                and len(support_refs) >= 7
                and len(slot_refs) >= 7
                and lane_token_hits == []
            ),
            "details": "Given that the lawful refined shift-interface-rule support stack is already packed and the exact statement still remains unexported, no residual passive same-lane loophole remains below the current blocker.",
        },
    ]

    failed_checks = [item["id"] for item in checks if not item["pass"]]
    all_pass = all(item["pass"] for item in checks)

    theorem_result = {
        "boundary_name": "CurrentStrictAlphaSLawfulRefinedExactRequiredFormStatementShiftInterfaceRuleDomainAdmissionExactRequiredFormStatementSameLaneExhaustionBoundary_strict_v1",
        "exact_required_form_statement_exported_on_current_repo_state": False if all_pass else None,
        "same_lane_passive_groundwork_exhausted": True if all_pass else None,
        "next_honest_move_requires_continuation_boundary_export": all_pass,
        "no_false_pass": True,
    }

    status = (
        "P870_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_EXACT_REQUIRED_FORM_STATEMENT_SHIFT_INTERFACE_RULE_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_SAME_LANE_EXHAUSTION_BOUNDARY_AUDIT_PROBE"
        if all_pass
        else "P870_REQUIRES_REVIEW"
    )

    artifact = {
        "stage": "P870",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p867_refined_exact_output_schema_statement_probe": rel(IN_P867),
            "p868_refined_exact_statement_required_form_probe": rel(IN_P868),
            "p869_refined_exact_required_form_statement_probe": rel(IN_P869),
            "f867_refined_exact_output_schema_statement_target_packet": rel(IN_F867),
            "f868_refined_exact_statement_required_form_target_packet": rel(IN_F868),
            "f869_refined_exact_required_form_statement_target_packet": rel(IN_F869),
            "s2_priority_packet": rel(IN_S2),
        },
        "repo_scan_token_hits_for_current_lane_exact_required_form_statement": lane_token_hits,
        "theorem_result": theorem_result,
        "checks": checks,
        "failed_checks": failed_checks,
        "current_honest_reading": [
            "The exact target chain is already frozen down to the exact required-form statement object itself for the current lawful refined shift-interface-rule domain-admission lane.",
            "The current repo state still keeps that exact required-form statement unexported, while the neighboring support and slot context is already packed inside F869.",
            "So no residual passive same-lane loophole remains below the current blocker on this lawful refined shift-interface-rule lane.",
        ],
        "recommended_next_packet": {
            "id": "F870_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_EXACT_REQUIRED_FORM_STATEMENT_SHIFT_INTERFACE_RULE_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_CONTINUATION_BOUNDARY_PACKET",
            "goal": "Export the continuation boundary after the current same-lane passive groundwork is exhausted but the exact required-form statement still remains missing on the lawful refined shift-interface-rule exact-required-form-statement lane.",
            "export_object_id": "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_exact_required_form_statement_continuation_boundary_v1",
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "P870",
        "status": status,
        "as_of": AS_OF,
        "exact_required_form_statement_exported_on_current_repo_state": theorem_result[
            "exact_required_form_statement_exported_on_current_repo_state"
        ],
        "same_lane_passive_groundwork_exhausted": theorem_result["same_lane_passive_groundwork_exhausted"],
        "recommended_next_packet": artifact["recommended_next_packet"]["id"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
