#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-21"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_F923 = GENERATED / "f923_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_domain_admission_exact_statement_required_form_target_packet.json"
IN_F924 = GENERATED / "f924_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_domain_admission_exact_required_form_statement_target_packet.json"
IN_P924 = GENERATED / "p924_current_strict_alpha_s_required_form_statement_class_candidate_supported_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_domain_admission_exact_required_form_statement_blocked.json"
IN_F925 = GENERATED / "f925_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_exact_required_form_statement_continuation_boundary_packet.json"

OUT = GENERATED / "p926_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_rfs_class_support_reuse_as_same_lane_exact_rfs_source_audit_probe.json"
OUT_SUMMARY = GENERATED / "p926_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_rfs_class_support_reuse_as_same_lane_exact_rfs_source_audit_probe_summary.json"


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

    prereq = [IN_F923, IN_F924, IN_P924, IN_F925]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P926",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    f923 = load_json(IN_F923)
    f924 = load_json(IN_F924)
    p924 = load_json(IN_P924)
    f925 = load_json(IN_F925)

    f923_target = f923.get("target_object") or {}
    f924_target = f924.get("target_object") or {}
    p924_theorem = p924.get("theorem_result") or {}
    p924_support_stack = p924.get("required_form_statement_class_support_stack") or {}
    f925_export = f925.get("exported_object") or {}

    support_refs = (f924.get("target_refs") or {}).get("required_form_statement_class_candidate_support_refs") or []
    slot_refs = (f924.get("target_refs") or {}).get("neighboring_statement_or_form_slot_refs") or []

    token_hits: list[str] = []
    current_lane_candidate_id = (
        "required_form_statement_class_support_stack_reuse_same_lane_"
        "lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_"
        "domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_"
        "domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_"
        "domain_admission_exact_required_form_statement_domain_admission_source_candidate_v1"
    )
    for path in sorted(GENERATED.glob("*.json")):
        name = path.name
        if name.startswith(("p926_", "f926_")):
            continue
        text = path.read_text(encoding="utf-8")
        if current_lane_candidate_id in text:
            token_hits.append(rel(path))

    checks = [
        {
            "id": "f925_already_requires_new_same_lane_source_or_provider_shift",
            "pass": (
                f925.get("status")
                == "F925_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_CONTINUATION_BOUNDARY_PACKET_NO_FALSE_PASS"
                and f925_export.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_continuation_boundary_v1"
                and f925_export.get("admitted_next_move_classes")
                == [
                    "export_one_genuinely_new_same_lane_exact_required_form_statement_source",
                    "shift_to_a_different_required_form_statement_provider_class_lane",
                ]
            ),
            "details": "F925 already narrows continuation to one genuinely new same-lane exact required-form statement source or a provider-class shift.",
        },
        {
            "id": "f924_current_exact_required_form_statement_target_and_support_stack_exist",
            "pass": (
                f924.get("status")
                == "F924_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_TARGET_PACKET_NO_FALSE_PASS"
                and f924_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_target_v1"
                and has_required_field(f924_target, "required_form_statement_class_candidate_support_refs")
                and has_required_field(f924_target, "neighboring_statement_or_form_slot_refs")
                and len(support_refs) >= 7
                and len(slot_refs) >= 7
            ),
            "details": "F924 already freezes the current deeper lawful refined domain-admission exact required-form statement target and packs the neighboring support/slot context.",
        },
        {
            "id": "p924_keeps_support_stack_below_exact_export",
            "pass": (
                p924.get("status")
                == "P924_CURRENT_STRICT_ALPHA_S_REQUIRED_FORM_STATEMENT_CLASS_CANDIDATE_SUPPORTED_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_BLOCKED"
                and p924_theorem.get("required_form_statement_class_candidate_supported_now") is True
                and p924_theorem.get("exact_required_form_statement_exported_now") is False
                and p924_support_stack.get("support_grade") == "candidate_supported_not_yet_exported"
                and "none exports the exact deeper lawful refined domain-admission required-form statement needed by F923."
                in (p924_support_stack.get("scope_limit") or "")
            ),
            "details": "P924 already keeps the deeper lawful refined domain-admission support stack explicitly below exact export and denies that it exports the exact required-form statement.",
        },
        {
            "id": "f923_upstream_target_still_names_exact_required_form_statement_as_missing",
            "pass": (
                f923.get("status")
                == "F923_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_STATEMENT_REQUIRED_FORM_TARGET_PACKET_NO_FALSE_PASS"
                and f923_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement.domain_admission_exact_statement_required_form_target_v1".replace(".", "_")
                and has_required_field(f923_target, "exact_required_form_statement_ref")
            ),
            "details": "The upstream deeper lawful refined domain-admission exact statement-required-form target still names exact_required_form_statement_ref as missing, so it is not itself the exact statement source.",
        },
        {
            "id": "same_lane_support_stack_remains_nonidentical_support_only",
            "pass": (
                "only through neighboring target fields and neighboring old-lane targets."
                in " ".join(f924.get("current_honest_reading") or [])
                and "No current export yet names the exact required-form statement required by F923"
                in " ".join(f924.get("current_honest_reading") or [])
                and "Does not claim that any neighboring statement or form slot silently discharges the lawful refined new lane."
                in (f924.get("hard_limits") or [])
            ),
            "details": "The same-lane deeper lawful refined domain-admission support stack remains neighboring nonidentical support only and does not discharge the new lane.",
        },
        {
            "id": "no_current_same_lane_exact_required_form_statement_source_is_exported",
            "pass": (
                p924_theorem.get("exact_required_form_statement_exported_now") is False
                and p924_support_stack.get("support_grade") == "candidate_supported_not_yet_exported"
                and len(support_refs) >= 7
                and token_hits == []
            ),
            "details": "Nothing currently exported upgrades the neighboring deeper lawful refined domain-admission support stack into a genuinely new same-lane exact required-form statement source for the current alpha_s lane.",
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
        "P926_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_REQUIRED_FORM_STATEMENT_CLASS_SUPPORT_REUSE_AS_SAME_LANE_EXACT_REQUIRED_FORM_STATEMENT_SOURCE_NEGATIVE_PROVIDER_SHIFT_REQUIRED"
        if all_pass
        else "P926_REQUIRES_REVIEW"
    )

    artifact = {
        "stage": "P926",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "f923_exact_statement_required_form_target_packet": rel(IN_F923),
            "f924_exact_required_form_statement_target_packet": rel(IN_F924),
            "p924_required_form_statement_probe": rel(IN_P924),
            "f925_continuation_boundary_packet": rel(IN_F925),
        },
        "repo_scan_token_hits_for_same_lane_exact_required_form_statement_source_candidate": token_hits,
        "rejected_same_lane_source_candidate_class": {
            "candidate_id": current_lane_candidate_id,
            "exact_required_form_statement_target_ref": f924_target.get("object_id"),
            "continuation_boundary_ref": f925_export.get("object_id"),
            "rejection_basis_refs": [
                rel(IN_F923),
                rel(IN_F924),
                rel(IN_P924),
                rel(IN_F925),
            ],
        },
        "theorem_result": theorem_result,
        "checks": checks,
        "failed_checks": failed_checks,
        "current_honest_reading": [
            "Required-form statement-class support does exist on the deeper lawful refined domain-admission current repo state.",
            "But that support stack remains candidate-only, neighboring, and explicitly nonidentical to the missing exact required-form statement.",
            "So it is not currently exported as a genuinely new same-lane exact required-form statement source for the current alpha_s deeper lawful refined domain-admission lane.",
        ],
        "recommended_next_packet": {
            "id": "F926_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_PROVIDER_CLASS_SHIFT_REQUIREMENT_PACKET",
            "goal": "Freeze the current-repo-state requirement that deeper lawful refined alpha_s continuation now proceed by provider-class shift rather than same-lane promotion of neighboring required-form statement support.",
            "export_object_id": "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement.domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_provider_class_shift_requirement_v1".replace(".", "_"),
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "P926",
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
