#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P845 = GENERATED / "p845_current_strict_alpha_s_output_schema_statement_class_candidate_supported_lawful_exact_required_form_statement_domain_admission_exact_output_schema_statement_blocked.json"
IN_P846 = GENERATED / "p846_current_strict_alpha_s_statement_form_class_candidate_supported_lawful_exact_required_form_statement_domain_admission_exact_statement_required_form_blocked.json"
IN_P847 = GENERATED / "p847_current_strict_alpha_s_required_form_statement_class_candidate_supported_lawful_exact_required_form_statement_domain_admission_exact_required_form_statement_blocked.json"
IN_F845 = GENERATED / "f845_current_strict_alpha_s_lawful_exact_required_form_statement_domain_admission_exact_output_schema_statement_target_packet.json"
IN_F846 = GENERATED / "f846_current_strict_alpha_s_lawful_exact_required_form_statement_domain_admission_exact_statement_required_form_target_packet.json"
IN_F847 = GENERATED / "f847_current_strict_alpha_s_lawful_exact_required_form_statement_domain_admission_exact_required_form_statement_target_packet.json"
IN_S2 = ROOT / "S2_CURRENT_FAR_STRATEGIC_PRIORITY_REORIENTATION_PACKET.md"

OUT = GENERATED / "p848_current_strict_alpha_s_refined_exact_required_form_statement_same_lane_exhaustion_boundary_audit_probe.json"
OUT_SUMMARY = GENERATED / "p848_current_strict_alpha_s_refined_exact_required_form_statement_same_lane_exhaustion_boundary_audit_probe_summary.json"


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

    prereq = [IN_P845, IN_P846, IN_P847, IN_F845, IN_F846, IN_F847, IN_S2]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P848",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p845 = load_json(IN_P845)
    p846 = load_json(IN_P846)
    p847 = load_json(IN_P847)
    f845 = load_json(IN_F845)
    f846 = load_json(IN_F846)
    f847 = load_json(IN_F847)
    s2_text = load_text(IN_S2)

    p845_theorem = p845.get("theorem_result") or {}
    p846_theorem = p846.get("theorem_result") or {}
    p847_theorem = p847.get("theorem_result") or {}
    f845_target = f845.get("target_object") or {}
    f846_target = f846.get("target_object") or {}
    f847_target = f847.get("target_object") or {}

    support_refs = (f847.get("target_refs") or {}).get("required_form_statement_class_candidate_support_refs") or []
    slot_refs = (f847.get("target_refs") or {}).get("neighboring_statement_or_form_slot_refs") or []

    lane_token_hits: list[str] = []
    lane_export_token = "lawful_exact_required_form_statement_domain_admission_refined_exact_required_form_statement"
    for path in sorted(GENERATED.glob("*.json")):
        name = path.name
        if name.startswith(("p845_", "f845_", "p846_", "f846_", "p847_", "f847_", "p848_", "f848_")):
            continue
        text = path.read_text(encoding="utf-8")
        if lane_export_token in text:
            lane_token_hits.append(rel(path))

    checks = [
        {
            "id": "f845_refined_exact_output_schema_statement_target_already_frozen",
            "pass": (
                f845.get("status")
                == "F845_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_OUTPUT_SCHEMA_STATEMENT_TARGET_PACKET_NO_FALSE_PASS"
                and f845_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_exact_required_form_statement_domain_admission_refined_exact_output_schema_statement_target_v1"
            ),
            "details": "F845 already freezes the refined exact output-schema statement target for the current lane.",
        },
        {
            "id": "f846_refined_exact_statement_required_form_target_already_frozen",
            "pass": (
                f846.get("status")
                == "F846_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_STATEMENT_REQUIRED_FORM_TARGET_PACKET_NO_FALSE_PASS"
                and f846_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_exact_required_form_statement_domain_admission_refined_exact_statement_required_form_target_v1"
                and has_required_field(f846_target, "exact_required_form_statement_ref")
            ),
            "details": "F846 already freezes the refined exact statement-required-form target and names exact_required_form_statement_ref as one exact missing field.",
        },
        {
            "id": "f847_refined_exact_required_form_statement_target_already_frozen",
            "pass": (
                f847.get("status")
                == "F847_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_TARGET_PACKET_NO_FALSE_PASS"
                and f847_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_exact_required_form_statement_domain_admission_refined_exact_required_form_statement_target_v1"
                and has_required_field(f847_target, "exact_required_form_statement_ref")
            ),
            "details": "F847 already freezes the refined exact required-form statement target itself.",
        },
        {
            "id": "p847_keeps_refined_exact_required_form_statement_unexported",
            "pass": (
                p847.get("status")
                == "P847_CURRENT_STRICT_ALPHA_S_REQUIRED_FORM_STATEMENT_CLASS_CANDIDATE_SUPPORTED_LAWFUL_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_BLOCKED"
                and p847_theorem.get("exact_required_form_statement_exported_now") is False
                and p847_theorem.get("sharp_blocker_field") == "exact_required_form_statement_ref"
            ),
            "details": "P847 already keeps the refined exact required-form statement unexported on the current repo state.",
        },
        {
            "id": "f847_already_packs_neighboring_support_and_slot_context",
            "pass": (
                len(support_refs) >= 7
                and len(slot_refs) >= 7
                and has_required_field(f847_target, "required_form_statement_class_candidate_support_refs")
                and has_required_field(f847_target, "neighboring_statement_or_form_slot_refs")
            ),
            "details": "F847 already packs the neighboring statement/form-class support refs and neighboring slot refs into one exact refined target context.",
        },
        {
            "id": "repo_scan_finds_no_current_lane_exact_required_form_statement_export",
            "pass": lane_token_hits == [],
            "details": "Repo scan finds no generated artifact exporting the refined current-lane lawful exact-required-form-statement domain-admission exact required-form statement outside the new frozen F847 lineage itself.",
        },
        {
            "id": "blocker_chain_has_already_descended_to_current_floor",
            "pass": (
                p845.get("status")
                == "P845_CURRENT_STRICT_ALPHA_S_OUTPUT_SCHEMA_STATEMENT_CLASS_CANDIDATE_SUPPORTED_LAWFUL_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_OUTPUT_SCHEMA_STATEMENT_BLOCKED"
                and p845_theorem.get("sharp_blocker_field") == "exact_output_schema_statement"
                and p846.get("status")
                == "P846_CURRENT_STRICT_ALPHA_S_STATEMENT_FORM_CLASS_CANDIDATE_SUPPORTED_LAWFUL_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_STATEMENT_REQUIRED_FORM_BLOCKED"
                and p846_theorem.get("sharp_blocker_field") == "exact_statement_required_form_ref"
                and p847.get("status")
                == "P847_CURRENT_STRICT_ALPHA_S_REQUIRED_FORM_STATEMENT_CLASS_CANDIDATE_SUPPORTED_LAWFUL_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_BLOCKED"
                and p847_theorem.get("sharp_blocker_field") == "exact_required_form_statement_ref"
            ),
            "details": "The refined blocker chain has already been pushed down from exact output-schema statement to exact statement-required form to exact required-form statement, so the current refined lane is already at its present passive floor.",
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
                p847_theorem.get("required_form_statement_class_candidate_supported_now") is True
                and p847_theorem.get("exact_required_form_statement_exported_now") is False
                and len(support_refs) >= 7
                and len(slot_refs) >= 7
                and lane_token_hits == []
            ),
            "details": "Given that the refined support stack is already packed and the exact statement still remains unexported, no residual passive same-lane loophole remains below the current blocker.",
        },
    ]

    failed_checks = [item["id"] for item in checks if not item["pass"]]
    all_pass = all(item["pass"] for item in checks)

    theorem_result = {
        "boundary_name": "CurrentStrictAlphaSLawfulExactRequiredFormStatementDomainAdmissionRefinedExactRequiredFormStatementSameLaneExhaustionBoundary_strict_v1",
        "exact_required_form_statement_exported_on_current_repo_state": False if all_pass else None,
        "same_lane_passive_groundwork_exhausted": True if all_pass else None,
        "next_honest_move_requires_continuation_boundary_export": all_pass,
        "no_false_pass": True,
    }

    status = (
        "P848_CURRENT_STRICT_ALPHA_S_REFINED_EXACT_REQUIRED_FORM_STATEMENT_SAME_LANE_EXHAUSTION_BOUNDARY_AUDIT_PROBE"
        if all_pass
        else "P848_REQUIRES_REVIEW"
    )

    artifact = {
        "stage": "P848",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p845_refined_exact_output_schema_statement_probe": rel(IN_P845),
            "p846_refined_exact_statement_required_form_probe": rel(IN_P846),
            "p847_refined_exact_required_form_statement_probe": rel(IN_P847),
            "f845_refined_exact_output_schema_statement_target_packet": rel(IN_F845),
            "f846_refined_exact_statement_required_form_target_packet": rel(IN_F846),
            "f847_refined_exact_required_form_statement_target_packet": rel(IN_F847),
            "s2_priority_packet": rel(IN_S2),
        },
        "repo_scan_token_hits_for_current_lane_exact_required_form_statement": lane_token_hits,
        "theorem_result": theorem_result,
        "checks": checks,
        "failed_checks": failed_checks,
        "current_honest_reading": [
            "The exact target chain is already frozen down to the exact required-form statement object itself for the current refined lawful domain-admission lane.",
            "The current repo state still keeps that exact required-form statement unexported, while the neighboring support and slot context is already packed inside F847.",
            "So no residual passive same-lane loophole remains below the current blocker on this refined lane.",
        ],
        "recommended_next_packet": {
            "id": "F848_CURRENT_STRICT_ALPHA_S_REFINED_EXACT_REQUIRED_FORM_STATEMENT_CONTINUATION_BOUNDARY_PACKET",
            "goal": "Export the refined continuation boundary after the current same-lane passive groundwork is exhausted but the exact required-form statement still remains missing on the lawful exact-required-form-statement domain-admission lane.",
            "export_object_id": "alpha_s_pair12_source_side_branch_selection_provider_lawful_exact_required_form_statement_domain_admission_refined_exact_required_form_statement_continuation_boundary_v1",
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "P848",
        "status": status,
        "as_of": AS_OF,
        "exact_required_form_statement_exported_on_current_repo_state": theorem_result[
            "exact_required_form_statement_exported_on_current_repo_state"
        ],
        "same_lane_passive_groundwork_exhausted": theorem_result["same_lane_passive_groundwork_exhausted"],
        "next_honest_move_requires_continuation_boundary_export": theorem_result[
            "next_honest_move_requires_continuation_boundary_export"
        ],
        "recommended_next_packet": artifact["recommended_next_packet"]["id"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
