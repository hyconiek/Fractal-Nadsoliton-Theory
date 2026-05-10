#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-20"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_F857 = GENERATED / "f857_current_strict_alpha_s_lawful_refined_exact_required_form_statement_domain_admission_exact_statement_required_form_target_packet.json"
IN_F858 = GENERATED / "f858_current_strict_alpha_s_lawful_refined_exact_required_form_statement_domain_admission_exact_required_form_statement_target_packet.json"
IN_P858 = GENERATED / "p858_current_strict_alpha_s_required_form_statement_class_candidate_supported_lawful_refined_exact_required_form_statement_domain_admission_exact_required_form_statement_blocked.json"
IN_F859 = GENERATED / "f859_current_strict_alpha_s_lawful_refined_exact_required_form_statement_continuation_boundary_packet.json"

OUT = GENERATED / "p860_current_strict_alpha_s_lawful_refined_required_form_statement_class_support_reuse_as_same_lane_exact_required_form_statement_source_audit_probe.json"
OUT_SUMMARY = GENERATED / "p860_current_strict_alpha_s_lawful_refined_required_form_statement_class_support_reuse_as_same_lane_exact_required_form_statement_source_audit_probe_summary.json"


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

    prereq = [IN_F857, IN_F858, IN_P858, IN_F859]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P860",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    f857 = load_json(IN_F857)
    f858 = load_json(IN_F858)
    p858 = load_json(IN_P858)
    f859 = load_json(IN_F859)

    f857_target = f857.get("target_object") or {}
    f858_target = f858.get("target_object") or {}
    p858_theorem = p858.get("theorem_result") or {}
    p858_support_stack = p858.get("required_form_statement_class_support_stack") or {}
    f859_export = f859.get("exported_object") or {}

    support_refs = (f858.get("target_refs") or {}).get("required_form_statement_class_candidate_support_refs") or []
    slot_refs = (f858.get("target_refs") or {}).get("neighboring_statement_or_form_slot_refs") or []

    token_hits: list[str] = []
    current_lane_candidate_id = (
        "required_form_statement_class_support_stack_reuse_same_lane_"
        "lawful_refined_exact_required_form_statement_source_candidate_v1"
    )
    for path in sorted(GENERATED.glob("*.json")):
        name = path.name
        if name.startswith(("p860_", "f860_")):
            continue
        text = path.read_text(encoding="utf-8")
        if current_lane_candidate_id in text:
            token_hits.append(rel(path))

    checks = [
        {
            "id": "f859_already_requires_new_same_lane_source_or_provider_shift",
            "pass": (
                f859.get("status")
                == "F859_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_EXACT_REQUIRED_FORM_STATEMENT_CONTINUATION_BOUNDARY_PACKET_NO_FALSE_PASS"
                and f859_export.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_exact_required_form_statement_continuation_boundary_v1"
                and f859_export.get("admitted_next_move_classes")
                == [
                    "export_one_genuinely_new_same_lane_exact_required_form_statement_source",
                    "shift_to_a_different_required_form_statement_provider_class_lane",
                ]
            ),
            "details": "F859 already narrows continuation to one genuinely new same-lane exact required-form statement source or a provider-class shift.",
        },
        {
            "id": "f858_current_exact_required_form_statement_target_and_support_stack_exist",
            "pass": (
                f858.get("status")
                == "F858_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_TARGET_PACKET_NO_FALSE_PASS"
                and f858_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_exact_required_form_statement_domain_admission_exact_required_form_statement_target_v1"
                and has_required_field(f858_target, "required_form_statement_class_candidate_support_refs")
                and has_required_field(f858_target, "neighboring_statement_or_form_slot_refs")
                and len(support_refs) >= 7
                and len(slot_refs) >= 7
            ),
            "details": "F858 already freezes the current lawful refined exact required-form statement target and packs the neighboring support/slot context.",
        },
        {
            "id": "p858_keeps_support_stack_below_exact_export",
            "pass": (
                p858.get("status")
                == "P858_CURRENT_STRICT_ALPHA_S_REQUIRED_FORM_STATEMENT_CLASS_CANDIDATE_SUPPORTED_LAWFUL_REFINED_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_BLOCKED"
                and p858_theorem.get("required_form_statement_class_candidate_supported_now") is True
                and p858_theorem.get("exact_required_form_statement_exported_now") is False
                and p858_support_stack.get("support_grade") == "candidate_supported_not_yet_exported"
                and "none exports the exact refined required-form statement needed by F857."
                in (p858_support_stack.get("scope_limit") or "")
            ),
            "details": "P858 already keeps the lawful refined support stack explicitly below exact export and denies that it exports the exact required-form statement.",
        },
        {
            "id": "f857_upstream_target_still_names_exact_required_form_statement_as_missing",
            "pass": (
                f857.get("status")
                == "F857_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_STATEMENT_REQUIRED_FORM_TARGET_PACKET_NO_FALSE_PASS"
                and f857_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_exact_required_form_statement_domain_admission_exact_statement_required_form_target_v1"
                and has_required_field(f857_target, "exact_required_form_statement_ref")
            ),
            "details": "The upstream lawful refined exact statement-required-form target still names exact_required_form_statement_ref as missing, so it is not itself the exact statement source.",
        },
        {
            "id": "same_lane_support_stack_remains_nonidentical_support_only",
            "pass": (
                "only through neighboring target fields and neighboring old-lane targets"
                in " ".join(f858.get("current_honest_reading") or [])
                and "No current export yet names the exact required-form statement required by F857."
                in " ".join(f858.get("current_honest_reading") or [])
                and "Does not claim that any neighboring statement or form slot silently discharges the refined new lane."
                in (f858.get("hard_limits") or [])
            ),
            "details": "The same-lane lawful refined support stack remains neighboring nonidentical support only and does not discharge the new lane.",
        },
        {
            "id": "no_current_same_lane_exact_required_form_statement_source_is_exported",
            "pass": (
                p858_theorem.get("exact_required_form_statement_exported_now") is False
                and p858_support_stack.get("support_grade") == "candidate_supported_not_yet_exported"
                and len(support_refs) >= 7
            ),
            "details": "Nothing currently exported upgrades the neighboring lawful refined support stack into a genuinely new same-lane exact required-form statement source for the current alpha_s lane.",
        },
    ]

    failed_checks = [item["id"] for item in checks if not item["pass"]]
    all_pass = all(item["pass"] for item in checks)

    theorem_result = {
        "required_form_statement_class_support_stack_exists": checks[1]["pass"],
        "required_form_statement_class_support_stack_is_same_lane_exact_required_form_statement_source_for_alpha_s": False
        if all_pass
        else None,
        "current_repo_exports_no_genuinely_new_same_lane_exact_required_form_statement_source_for_alpha_s": all_pass,
        "next_honest_move_requires_required_form_statement_provider_class_shift": all_pass,
        "no_false_pass": True,
    }

    status = (
        "P860_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_REQUIRED_FORM_STATEMENT_CLASS_SUPPORT_REUSE_AS_SAME_LANE_EXACT_REQUIRED_FORM_STATEMENT_SOURCE_NEGATIVE_PROVIDER_SHIFT_REQUIRED"
        if all_pass
        else "P860_REQUIRES_REVIEW"
    )

    artifact = {
        "stage": "P860",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "f857_refined_exact_statement_required_form_target_packet": rel(IN_F857),
            "f858_refined_exact_required_form_statement_target_packet": rel(IN_F858),
            "p858_refined_required_form_statement_probe": rel(IN_P858),
            "f859_lawful_refined_exact_required_form_statement_continuation_boundary_packet": rel(IN_F859),
        },
        "repo_scan_token_hits_for_same_lane_exact_required_form_statement_source_candidate": token_hits,
        "rejected_same_lane_source_candidate_class": {
            "candidate_id": "required_form_statement_class_support_stack_reuse_same_lane_lawful_refined_exact_required_form_statement_source_candidate_v1",
            "exact_required_form_statement_target_ref": f858_target.get("object_id"),
            "continuation_boundary_ref": f859_export.get("object_id"),
            "rejection_basis_refs": [
                rel(IN_F857),
                rel(IN_F858),
                rel(IN_P858),
                rel(IN_F859),
            ],
        },
        "theorem_result": theorem_result,
        "checks": checks,
        "failed_checks": failed_checks,
        "current_honest_reading": [
            "Required-form statement-class support does exist on the lawful refined current repo state.",
            "But that support stack remains candidate-only, neighboring, and explicitly nonidentical to the missing exact required-form statement.",
            "So it is not currently exported as a genuinely new same-lane exact required-form statement source for the current alpha_s lawful refined lane.",
        ],
        "recommended_next_packet": {
            "id": "F860_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_EXACT_REQUIRED_FORM_STATEMENT_PROVIDER_CLASS_SHIFT_REQUIREMENT_PACKET",
            "goal": "Freeze the current-repo-state requirement that lawful refined alpha_s continuation now proceed by provider-class shift rather than same-lane promotion of neighboring required-form statement support.",
            "export_object_id": "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_exact_required_form_statement_provider_class_shift_requirement_v1",
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "P860",
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
