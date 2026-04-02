#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-21"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P911 = GENERATED / "p911_current_strict_alpha_s_output_schema_statement_class_candidate_supported_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_output_schema_statement_blocked.json"
IN_P912 = GENERATED / "p912_current_strict_alpha_s_statement_form_class_candidate_supported_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_statement_required_form_blocked.json"
IN_P913 = GENERATED / "p913_current_strict_alpha_s_required_form_statement_class_candidate_supported_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_blocked.json"
IN_F911 = GENERATED / "f911_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_output_schema_statement_target_packet.json"
IN_F912 = GENERATED / "f912_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_statement_required_form_target_packet.json"
IN_F913 = GENERATED / "f913_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_target_packet.json"
IN_S2 = ROOT / "S2_CURRENT_FAR_STRATEGIC_PRIORITY_REORIENTATION_PACKET.md"

OUT = GENERATED / "p914_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_same_lane_exhaustion_boundary_audit_probe.json"
OUT_SUMMARY = GENERATED / "p914_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_same_lane_exhaustion_boundary_audit_probe_summary.json"


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

    prereq = [IN_P911, IN_P912, IN_P913, IN_F911, IN_F912, IN_F913, IN_S2]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P914",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p911 = load_json(IN_P911)
    p912 = load_json(IN_P912)
    p913 = load_json(IN_P913)
    f911 = load_json(IN_F911)
    f912 = load_json(IN_F912)
    f913 = load_json(IN_F913)
    s2_text = load_text(IN_S2)

    p911_theorem = p911.get("theorem_result") or {}
    p912_theorem = p912.get("theorem_result") or {}
    p913_theorem = p913.get("theorem_result") or {}
    f911_target = f911.get("target_object") or {}
    f912_target = f912.get("target_object") or {}
    f913_target = f913.get("target_object") or {}

    support_refs = (f913.get("target_refs") or {}).get("required_form_statement_class_candidate_support_refs") or []
    slot_refs = (f913.get("target_refs") or {}).get("neighboring_statement_or_form_slot_refs") or []

    lane_token_hits: list[str] = []
    lane_export_token = (
        "lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_"
        "domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_"
        "domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement"
    )
    for path in sorted(GENERATED.glob("*.json")):
        name = path.name
        if name.startswith(("p911_", "f911_", "p912_", "f912_", "p913_", "f913_", "p914_", "f914_")):
            continue
        text = path.read_text(encoding="utf-8")
        if lane_export_token in text:
            lane_token_hits.append(rel(path))

    checks = [
        {
            "id": "f911_exact_output_schema_statement_target_already_frozen",
            "pass": (
                f911.get("status")
                == "F911_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_OUTPUT_SCHEMA_STATEMENT_TARGET_PACKET_NO_FALSE_PASS"
                and f911_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_output_schema_statement_target_v1"
            ),
            "details": "F911 already freezes the deeper lawful refined domain-admission exact output-schema statement target for the current lane.",
        },
        {
            "id": "f912_exact_statement_required_form_target_already_frozen",
            "pass": (
                f912.get("status")
                == "F912_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_STATEMENT_REQUIRED_FORM_TARGET_PACKET_NO_FALSE_PASS"
                and f912_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_statement_required_form_target_v1"
                and has_required_field(f912_target, "exact_required_form_statement_ref")
            ),
            "details": "F912 already freezes the deeper lawful refined exact statement-required-form target and names exact_required_form_statement_ref as one exact missing field.",
        },
        {
            "id": "f913_exact_required_form_statement_target_already_frozen",
            "pass": (
                f913.get("status")
                == "F913_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_TARGET_PACKET_NO_FALSE_PASS"
                and f913_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_target_v1"
                and has_required_field(f913_target, "exact_required_form_statement_ref")
            ),
            "details": "F913 already freezes the deeper lawful refined domain-admission exact required-form statement target itself.",
        },
        {
            "id": "p913_keeps_exact_required_form_statement_unexported",
            "pass": (
                p913.get("status")
                == "P913_CURRENT_STRICT_ALPHA_S_REQUIRED_FORM_STATEMENT_CLASS_CANDIDATE_SUPPORTED_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_BLOCKED"
                and p913_theorem.get("exact_required_form_statement_exported_now") is False
                and p913_theorem.get("sharp_blocker_field") == "exact_required_form_statement_ref"
            ),
            "details": "P913 already keeps the deeper lawful refined domain-admission exact required-form statement unexported on the current repo state.",
        },
        {
            "id": "f913_already_packs_neighboring_support_and_slot_context",
            "pass": (
                len(support_refs) >= 7
                and len(slot_refs) >= 7
                and has_required_field(f913_target, "required_form_statement_class_candidate_support_refs")
                and has_required_field(f913_target, "neighboring_statement_or_form_slot_refs")
            ),
            "details": "F913 already packs the neighboring statement/form-class support refs and neighboring slot refs into one exact deeper lawful refined target context.",
        },
        {
            "id": "repo_scan_finds_no_current_lane_exact_required_form_statement_export",
            "pass": lane_token_hits == [],
            "details": "Repo scan finds no generated artifact exporting the current deeper lawful refined domain-admission exact required-form statement outside the new frozen F913 lineage itself.",
        },
        {
            "id": "blocker_chain_has_already_descended_to_current_floor",
            "pass": (
                p911.get("status")
                == "P911_CURRENT_STRICT_ALPHA_S_OUTPUT_SCHEMA_STATEMENT_CLASS_CANDIDATE_SUPPORTED_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_OUTPUT_SCHEMA_STATEMENT_BLOCKED"
                and p911_theorem.get("sharp_blocker_field") == "exact_output_schema_statement"
                and p912.get("status")
                == "P912_CURRENT_STRICT_ALPHA_S_STATEMENT_FORM_CLASS_CANDIDATE_SUPPORTED_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_STATEMENT_REQUIRED_FORM_BLOCKED"
                and p912_theorem.get("sharp_blocker_field") == "exact_statement_required_form_ref"
                and p913.get("status")
                == "P913_CURRENT_STRICT_ALPHA_S_REQUIRED_FORM_STATEMENT_CLASS_CANDIDATE_SUPPORTED_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_BLOCKED"
                and p913_theorem.get("sharp_blocker_field") == "exact_required_form_statement_ref"
            ),
            "details": "The deeper lawful refined domain-admission blocker chain has already been pushed down from exact output-schema statement to exact statement-required form to exact required-form statement, so the current lane is already at its present passive floor.",
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
                p913_theorem.get("required_form_statement_class_candidate_supported_now") is True
                and p913_theorem.get("exact_required_form_statement_exported_now") is False
                and len(support_refs) >= 7
                and len(slot_refs) >= 7
                and lane_token_hits == []
            ),
            "details": "Given that the deeper lawful refined support stack is already packed and the exact required-form statement still remains unexported, no residual passive same-lane loophole remains below the current blocker.",
        },
    ]

    failed_checks = [item["id"] for item in checks if not item["pass"]]
    all_pass = all(item["pass"] for item in checks)

    theorem_result = {
        "boundary_name": "CurrentStrictAlphaSLawfulRefinedSirDaErfsDaErfsDaExactRequiredFormStatementDomainAdmissionExactRequiredFormStatementDomainAdmissionExactRequiredFormStatementSameLaneExhaustionBoundary_strict_v1",
        "exact_required_form_statement_exported_on_current_repo_state": False if all_pass else None,
        "same_lane_passive_groundwork_exhausted": True if all_pass else None,
        "next_honest_move_requires_continuation_boundary_export": all_pass,
        "no_false_pass": True,
    }

    status = (
        "P914_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_SAME_LANE_EXHAUSTION_BOUNDARY_AUDIT_PROBE"
        if all_pass
        else "P914_REQUIRES_REVIEW"
    )

    artifact = {
        "stage": "P914",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p911_exact_output_schema_statement_probe": rel(IN_P911),
            "p912_exact_statement_required_form_probe": rel(IN_P912),
            "p913_exact_required_form_statement_probe": rel(IN_P913),
            "f911_exact_output_schema_statement_target_packet": rel(IN_F911),
            "f912_exact_statement_required_form_target_packet": rel(IN_F912),
            "f913_exact_required_form_statement_target_packet": rel(IN_F913),
            "s2_priority_packet": rel(IN_S2),
        },
        "repo_scan_token_hits_for_current_lane_exact_required_form_statement": lane_token_hits,
        "theorem_result": theorem_result,
        "checks": checks,
        "failed_checks": failed_checks,
        "current_honest_reading": [
            "The exact target chain is already frozen down to the exact required-form statement object itself for the current deeper lawful refined domain-admission lane.",
            "The current repo state still keeps that exact required-form statement unexported, while the neighboring support and slot context is already packed inside F913.",
            "So no residual passive same-lane loophole remains below the current blocker on this deeper lawful refined domain-admission lane.",
        ],
        "recommended_next_packet": {
            "id": "F914_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_CONTINUATION_BOUNDARY_PACKET",
            "goal": "Export the continuation boundary after the current same-lane passive groundwork is exhausted but the exact required-form statement still remains missing on the deeper lawful refined domain-admission lane.",
            "export_object_id": "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_continuation_boundary_v1",
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "P914",
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
