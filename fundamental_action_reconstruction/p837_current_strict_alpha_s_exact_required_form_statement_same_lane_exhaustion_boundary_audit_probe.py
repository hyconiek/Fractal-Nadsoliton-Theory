#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P834 = GENERATED / "p834_current_strict_alpha_s_output_schema_statement_class_candidate_supported_lawful_exact_required_form_statement_domain_admission_exact_output_schema_statement_blocked.json"
IN_P835 = GENERATED / "p835_current_strict_alpha_s_statement_form_class_candidate_supported_lawful_exact_required_form_statement_domain_admission_exact_statement_required_form_blocked.json"
IN_P836 = GENERATED / "p836_current_strict_alpha_s_required_form_statement_class_candidate_supported_lawful_exact_required_form_statement_domain_admission_exact_required_form_statement_blocked.json"
IN_F834 = GENERATED / "f834_current_strict_alpha_s_lawful_exact_required_form_statement_domain_admission_exact_output_schema_statement_target_packet.json"
IN_F835 = GENERATED / "f835_current_strict_alpha_s_lawful_exact_required_form_statement_domain_admission_exact_statement_required_form_target_packet.json"
IN_F836 = GENERATED / "f836_current_strict_alpha_s_lawful_exact_required_form_statement_domain_admission_exact_required_form_statement_target_packet.json"
IN_S2 = ROOT / "S2_CURRENT_FAR_STRATEGIC_PRIORITY_REORIENTATION_PACKET.md"

OUT = GENERATED / "p837_current_strict_alpha_s_exact_required_form_statement_same_lane_exhaustion_boundary_audit_probe.json"
OUT_SUMMARY = GENERATED / "p837_current_strict_alpha_s_exact_required_form_statement_same_lane_exhaustion_boundary_audit_probe_summary.json"


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

    prereq = [IN_P834, IN_P835, IN_P836, IN_F834, IN_F835, IN_F836, IN_S2]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P837",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p834 = load_json(IN_P834)
    p835 = load_json(IN_P835)
    p836 = load_json(IN_P836)
    f834 = load_json(IN_F834)
    f835 = load_json(IN_F835)
    f836 = load_json(IN_F836)
    s2_text = load_text(IN_S2)

    p834_theorem = p834.get("theorem_result") or {}
    p835_theorem = p835.get("theorem_result") or {}
    p836_theorem = p836.get("theorem_result") or {}
    f834_target = f834.get("target_object") or {}
    f835_target = f835.get("target_object") or {}
    f836_target = f836.get("target_object") or {}

    support_refs = (f836.get("target_refs") or {}).get("required_form_statement_class_candidate_support_refs") or []
    slot_refs = (f836.get("target_refs") or {}).get("neighboring_statement_or_form_slot_refs") or []

    lane_token_hits: list[str] = []
    lane_export_token = "lawful_exact_required_form_statement_domain_admission_exact_required_form_statement"
    for path in sorted(GENERATED.glob("*.json")):
        name = path.name
        if name.startswith(("p834_", "f834_", "p835_", "f835_", "p836_", "f836_", "p837_", "f837_")):
            continue
        text = path.read_text(encoding="utf-8")
        if lane_export_token in text:
            lane_token_hits.append(rel(path))

    checks = [
        {
            "id": "f834_exact_output_schema_statement_target_already_frozen",
            "pass": (
                f834.get("status")
                == "F834_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_OUTPUT_SCHEMA_STATEMENT_TARGET_PACKET_NO_FALSE_PASS"
                and f834_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_exact_required_form_statement_domain_admission_exact_output_schema_statement_target_v1"
            ),
            "details": "F834 already freezes the exact output-schema statement target for the current lane.",
        },
        {
            "id": "f835_exact_statement_required_form_target_already_frozen",
            "pass": (
                f835.get("status")
                == "F835_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_STATEMENT_REQUIRED_FORM_TARGET_PACKET_NO_FALSE_PASS"
                and f835_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_exact_required_form_statement_domain_admission_exact_statement_required_form_target_v1"
                and has_required_field(f835_target, "exact_required_form_statement_ref")
            ),
            "details": "F835 already freezes the exact statement-required-form target and names exact_required_form_statement_ref as one exact missing field.",
        },
        {
            "id": "f836_exact_required_form_statement_target_already_frozen",
            "pass": (
                f836.get("status")
                == "F836_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_TARGET_PACKET_NO_FALSE_PASS"
                and f836_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_exact_required_form_statement_domain_admission_exact_required_form_statement_target_v1"
                and has_required_field(f836_target, "exact_required_form_statement_ref")
            ),
            "details": "F836 already freezes the exact required-form statement target itself.",
        },
        {
            "id": "p836_keeps_exact_required_form_statement_unexported",
            "pass": (
                p836.get("status")
                == "P836_CURRENT_STRICT_ALPHA_S_REQUIRED_FORM_STATEMENT_CLASS_CANDIDATE_SUPPORTED_LAWFUL_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_BLOCKED"
                and p836_theorem.get("exact_required_form_statement_exported_now") is False
                and p836_theorem.get("sharp_blocker_field") == "exact_required_form_statement_ref"
            ),
            "details": "P836 already keeps the exact required-form statement unexported on the current repo state.",
        },
        {
            "id": "f836_already_packs_neighboring_support_and_slot_context",
            "pass": (
                len(support_refs) >= 8
                and len(slot_refs) >= 8
                and has_required_field(f836_target, "required_form_statement_class_candidate_support_refs")
                and has_required_field(f836_target, "neighboring_statement_or_form_slot_refs")
            ),
            "details": "F836 already packs the neighboring statement/form-class support refs and neighboring slot refs into one exact target context.",
        },
        {
            "id": "repo_scan_finds_no_current_lane_exact_required_form_statement_export",
            "pass": lane_token_hits == [],
            "details": "Repo scan finds no generated artifact exporting the current-lane lawful exact-required-form-statement domain-admission exact required-form statement outside the new frozen F836 lineage itself.",
        },
        {
            "id": "blocker_chain_has_already_descended_to_current_floor",
            "pass": (
                p834.get("status")
                == "P834_CURRENT_STRICT_ALPHA_S_OUTPUT_SCHEMA_STATEMENT_CLASS_CANDIDATE_SUPPORTED_LAWFUL_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_OUTPUT_SCHEMA_STATEMENT_BLOCKED"
                and p834_theorem.get("sharp_blocker_field") == "exact_output_schema_statement"
                and p835.get("status")
                == "P835_CURRENT_STRICT_ALPHA_S_STATEMENT_FORM_CLASS_CANDIDATE_SUPPORTED_LAWFUL_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_STATEMENT_REQUIRED_FORM_BLOCKED"
                and p835_theorem.get("sharp_blocker_field") == "exact_statement_required_form_ref"
                and p836.get("status")
                == "P836_CURRENT_STRICT_ALPHA_S_REQUIRED_FORM_STATEMENT_CLASS_CANDIDATE_SUPPORTED_LAWFUL_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_BLOCKED"
                and p836_theorem.get("sharp_blocker_field") == "exact_required_form_statement_ref"
            ),
            "details": "The blocker chain has already been pushed down from exact output-schema statement to exact statement-required form to exact required-form statement, so the current lane is already at its present passive floor.",
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
                p836_theorem.get("required_form_statement_class_candidate_supported_now") is True
                and p836_theorem.get("exact_required_form_statement_exported_now") is False
                and len(support_refs) >= 8
                and len(slot_refs) >= 8
                and lane_token_hits == []
            ),
            "details": "Given that the support stack is already packed and the exact statement still remains unexported, no residual passive same-lane loophole remains below the current blocker.",
        },
    ]

    failed_checks = [item["id"] for item in checks if not item["pass"]]
    all_pass = all(item["pass"] for item in checks)

    theorem_result = {
        "boundary_name": "CurrentStrictAlphaSLawfulExactRequiredFormStatementDomainAdmissionExactRequiredFormStatementSameLaneExhaustionBoundary_strict_v1",
        "exact_required_form_statement_exported_on_current_repo_state": False if all_pass else None,
        "same_lane_passive_groundwork_exhausted": True if all_pass else None,
        "next_honest_move_requires_continuation_boundary_export": all_pass,
        "no_false_pass": True,
    }

    status = (
        "P837_CURRENT_STRICT_ALPHA_S_EXACT_REQUIRED_FORM_STATEMENT_SAME_LANE_EXHAUSTION_BOUNDARY_AUDIT_PROBE"
        if all_pass
        else "P837_REQUIRES_REVIEW"
    )

    artifact = {
        "stage": "P837",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p834_exact_output_schema_statement_probe": rel(IN_P834),
            "p835_exact_statement_required_form_probe": rel(IN_P835),
            "p836_exact_required_form_statement_probe": rel(IN_P836),
            "f834_exact_output_schema_statement_target_packet": rel(IN_F834),
            "f835_exact_statement_required_form_target_packet": rel(IN_F835),
            "f836_exact_required_form_statement_target_packet": rel(IN_F836),
            "s2_priority_packet": rel(IN_S2),
        },
        "repo_scan_token_hits_for_current_lane_exact_required_form_statement": lane_token_hits,
        "theorem_result": theorem_result,
        "checks": checks,
        "failed_checks": failed_checks,
        "current_honest_reading": [
            "The exact target chain is already frozen down to the exact required-form statement object itself for the current lawful domain-admission lane.",
            "The current repo state still keeps that exact required-form statement unexported, while the neighboring support and slot context is already packed inside F836.",
            "So no residual passive same-lane loophole remains below the current blocker on this lane.",
        ],
        "recommended_next_packet": {
            "id": "F837_CURRENT_STRICT_ALPHA_S_EXACT_REQUIRED_FORM_STATEMENT_CONTINUATION_BOUNDARY_PACKET",
            "goal": "Export the continuation boundary after the current same-lane passive groundwork is exhausted but the exact required-form statement still remains missing on the lawful exact-required-form-statement domain-admission lane.",
            "export_object_id": "alpha_s_pair12_source_side_branch_selection_provider_lawful_exact_required_form_statement_domain_admission_exact_required_form_statement_continuation_boundary_v1",
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "P837",
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
