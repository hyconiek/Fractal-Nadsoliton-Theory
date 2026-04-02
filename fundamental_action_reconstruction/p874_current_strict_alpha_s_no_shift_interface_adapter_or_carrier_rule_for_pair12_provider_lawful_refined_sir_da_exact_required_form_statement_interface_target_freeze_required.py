#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-20"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_F873 = GENERATED / "f873_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_interface_target_packet.json"
IN_F872 = GENERATED / "f872_current_strict_alpha_s_lawful_refined_shift_interface_rule_domain_admission_different_required_form_statement_provider_class_shift_candidate_reference_packet.json"
IN_F863 = GENERATED / "f863_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_lawful_refined_exact_required_form_statement_adapter_or_carrier_identification_rule_target_packet.json"
IN_F852 = GENERATED / "f852_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_refined_exact_required_form_statement_adapter_or_carrier_identification_rule_target_packet.json"
IN_F841 = GENERATED / "f841_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_exact_required_form_statement_adapter_or_carrier_identification_rule_target_packet.json"
IN_F830 = GENERATED / "f830_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_exact_required_form_statement_adapter_or_carrier_identification_rule_target_packet.json"
IN_F819 = GENERATED / "f819_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_schema_adapter_or_carrier_identification_rule_target_packet.json"
IN_F810 = GENERATED / "f810_current_strict_alpha_s_strict_source_shannon_shift_to_domain_carrier_identification_or_adapter_rule_target_packet.json"
IN_P788 = GENERATED / "p788_current_alpha_s_dimensionless_or_normalized_replacement_route_probe.json"
IN_P764 = GENERATED / "p764_current_strict_t218_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_target_probe.json"

OUT = GENERATED / "p874_current_strict_alpha_s_no_shift_interface_adapter_or_carrier_rule_for_pair12_provider_lawful_refined_sir_da_exact_required_form_statement_interface_target_freeze_required.json"
OUT_SUMMARY = GENERATED / "p874_current_strict_alpha_s_no_shift_interface_adapter_or_carrier_rule_for_pair12_provider_lawful_refined_sir_da_exact_required_form_statement_interface_target_freeze_required_summary.json"


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


def find_exact_rule_hits(interface_target_id: str, rule_terms: list[str]) -> list[str]:
    hits: list[str] = []
    excluded_suffixes = [
        "generated/p873_current_strict_alpha_s_no_exact_pair12_source_side_branch_selection_provider_shift_to_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_interface_target_exported_target_freeze_required.json",
        "generated/p873_current_strict_alpha_s_no_exact_pair12_source_side_branch_selection_provider_shift_to_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_interface_target_exported_target_freeze_required_summary.json",
        "generated/f873_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_interface_target_packet.json",
        "generated/f873_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_interface_target_packet_summary.json",
        "P873_CURRENT_STRICT_ALPHA_S_NO_EXACT_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_SHIFT_TO_LAWFUL_REFINED_SHIFT_INTERFACE_RULE_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_INTERFACE_TARGET_EXPORTED_TARGET_FREEZE_REQUIRED.md",
        "F873_CURRENT_STRICT_ALPHA_S_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_SHIFT_TO_LAWFUL_REFINED_SHIFT_INTERFACE_RULE_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_INTERFACE_TARGET_PACKET.md",
        "p873_current_strict_alpha_s_no_exact_pair12_source_side_branch_selection_provider_shift_to_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_interface_target_exported_target_freeze_required.py",
        "f873_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_interface_target_packet.py",
        "generated/p874_current_strict_alpha_s_no_shift_interface_adapter_or_carrier_rule_for_pair12_provider_lawful_refined_sir_da_exact_required_form_statement_interface_target_freeze_required.json",
        "generated/p874_current_strict_alpha_s_no_shift_interface_adapter_or_carrier_rule_for_pair12_provider_lawful_refined_sir_da_exact_required_form_statement_interface_target_freeze_required_summary.json",
        "generated/f874_current_strict_alpha_s_pair12_provider_shift_to_lawful_refined_sir_da_exact_required_form_statement_adapter_or_carrier_rule_target_packet.json",
        "generated/f874_current_strict_alpha_s_pair12_provider_shift_to_lawful_refined_sir_da_exact_required_form_statement_adapter_or_carrier_rule_target_packet_summary.json",
        "P874_CURRENT_STRICT_ALPHA_S_NO_SHIFT_INTERFACE_ADAPTER_OR_CARRIER_RULE_FOR_PAIR12_PROVIDER_LAWFUL_REFINED_SIR_DA_EXACT_REQUIRED_FORM_STATEMENT_INTERFACE_TARGET_FREEZE_REQUIRED.md",
        "F874_CURRENT_STRICT_ALPHA_S_PAIR12_PROVIDER_SHIFT_TO_LAWFUL_REFINED_SIR_DA_EXACT_REQUIRED_FORM_STATEMENT_ADAPTER_OR_CARRIER_RULE_TARGET_PACKET.md",
        "p874_current_strict_alpha_s_no_shift_interface_adapter_or_carrier_rule_for_pair12_provider_lawful_refined_sir_da_exact_required_form_statement_interface_target_freeze_required.py",
        "f874_current_strict_alpha_s_pair12_provider_shift_to_lawful_refined_sir_da_exact_required_form_statement_adapter_or_carrier_rule_target_packet.py",
    ]
    norm_interface = normalize(interface_target_id)
    norm_terms = [normalize(term) for term in rule_terms]
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
        if norm_interface in hay and any(term in hay for term in norm_terms):
            hits.append(relpath)
    return sorted(hits)


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_F873, IN_F872, IN_F863, IN_F852, IN_F841, IN_F830, IN_F819, IN_F810, IN_P788, IN_P764]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P874",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    f873 = load_json(IN_F873)
    f872 = load_json(IN_F872)
    f863 = load_json(IN_F863)
    f852 = load_json(IN_F852)
    f841 = load_json(IN_F841)
    f830 = load_json(IN_F830)
    f819 = load_json(IN_F819)
    f810 = load_json(IN_F810)
    p788 = load_json(IN_P788)
    p764 = load_json(IN_P764)

    interface_target = f873.get("target_object") or {}
    interface_refs = f873.get("target_refs") or {}
    candidate_lane = f872.get("exported_object") or {}
    old_current_rule_target = f863.get("target_object") or {}
    old_current_rule_refs = f863.get("target_refs") or {}
    old_refined_rule_target = f852.get("target_object") or {}
    old_refined_rule_refs = f852.get("target_refs") or {}
    old_lawful_rule_target = f841.get("target_object") or {}
    old_lawful_rule_refs = f841.get("target_refs") or {}
    old_exact_rule_target = f830.get("target_object") or {}
    old_exact_rule_refs = f830.get("target_refs") or {}
    old_schema_rule_target = f819.get("target_object") or {}
    old_schema_rule_refs = f819.get("target_refs") or {}
    old_shannon_rule_target = f810.get("target_object") or {}
    old_shannon_rule_refs = f810.get("target_refs") or {}
    p764_theorem = p764.get("theorem_result") or {}

    current_interface_target_id = (
        "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_"
        "domain_admission_exact_required_form_statement_shift_interface_target_v1"
    )

    exact_rule_scan_hits = find_exact_rule_hits(
        current_interface_target_id,
        [
            "adapter",
            "carrier identification",
            "adapter action schema",
            "carrier-safe",
            "domain admission rule",
        ],
    )

    checks = [
        {
            "id": "f873_freezes_lawful_refined_shift_interface_rule_interface_target_and_requires_adapter_or_carrier_rule",
            "pass": (
                f873.get("status")
                == "F873_EXECUTED_CURRENT_STRICT_ALPHA_S_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_SHIFT_TO_LAWFUL_REFINED_SHIFT_INTERFACE_RULE_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_INTERFACE_TARGET_PACKET_NO_FALSE_PASS"
                and interface_target.get("object_id") == current_interface_target_id
                and any(
                    isinstance(item, dict)
                    and item.get("name") == "shift_interface_adapter_or_carrier_identification_rule_ref"
                    for item in (interface_target.get("required_fields") or [])
                )
            ),
            "details": "F873 already freezes the lawful refined shift-interface-rule interface target and names the still-missing adapter/carrier-identification rule slot.",
        },
        {
            "id": "f872_still_keeps_t213_t216_lane_at_shift_interface_rule_candidate_reference_only_grade",
            "pass": (
                f872.get("status")
                == "F872_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_SHIFT_INTERFACE_RULE_DOMAIN_ADMISSION_DIFFERENT_REQUIRED_FORM_STATEMENT_PROVIDER_CLASS_SHIFT_CANDIDATE_REFERENCE_PACKET_NO_FALSE_PASS"
                and candidate_lane.get("candidate_reference_lane_grade") == "reference_context_candidate_only"
                and candidate_lane.get("alpha_s_shift_interface_status") == "blocked_nonexport"
                and candidate_lane.get("provider_class_shift_realization_status") == "not_realized"
                and candidate_lane.get("exact_required_form_statement_status") == "blocked_nonexport"
            ),
            "details": "The T213/T216 lane still remains shift-interface-rule candidate-reference only and does not realize shift.",
        },
        {
            "id": "p764_keeps_the_lane_own_missing_interface_below_actual_export",
            "pass": (
                p764.get("status")
                == "PASS_STRICT_T218_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_TARGET_EXPORTED"
                and p764_theorem.get("current_t218_target_is_future_route_only") is True
                and p764_theorem.get("current_t218_target_remains_below_actual_interface_export_and_below_t176_discharge")
                is True
            ),
            "details": "The candidate lane still has its own unresolved future-only interface target below actual export.",
        },
        {
            "id": "p788_still_reports_no_generic_exported_alpha_s_adapter",
            "pass": (
                p788.get("status")
                == "P788_CURRENT_ALPHA_S_DIMENSIONLESS_OR_NORMALIZED_REPLACEMENT_ROUTE_BLOCKED_ON_CURRENT_REPO_STATE"
                and any(
                    isinstance(item, dict)
                    and item.get("id") == "no_exported_dimensionless_alpha_s_adapter_detected"
                    and item.get("pass") is True
                    for item in (p788.get("checks") or [])
                )
                and (p788.get("adapter_token_hits") or {}) == {}
            ),
            "details": "The repo still exports no generic alpha_s adapter even on the older canonical lane.",
        },
        {
            "id": "older_f863_f852_f841_f830_f819_and_f810_rule_targets_are_lane_specific_and_not_silently_reusable",
            "pass": (
                f863.get("status")
                == "F863_EXECUTED_CURRENT_STRICT_ALPHA_S_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_SHIFT_TO_LAWFUL_REFINED_EXACT_REQUIRED_FORM_STATEMENT_ADAPTER_OR_CARRIER_IDENTIFICATION_RULE_TARGET_PACKET_NO_FALSE_PASS"
                and old_current_rule_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_exact_required_form_statement_shift_interface_adapter_or_carrier_identification_rule_target_v1"
                and old_current_rule_refs.get("shift_to_lawful_refined_exact_required_form_statement_interface_target_ref")
                != interface_target.get("object_id")
                and old_current_rule_refs.get("different_required_form_statement_provider_class_shift_candidate_reference_lane_ref")
                != interface_refs.get("different_required_form_statement_provider_class_shift_candidate_reference_lane_ref")
                and f852.get("status")
                == "F852_EXECUTED_CURRENT_STRICT_ALPHA_S_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_SHIFT_TO_REFINED_EXACT_REQUIRED_FORM_STATEMENT_ADAPTER_OR_CARRIER_IDENTIFICATION_RULE_TARGET_PACKET_NO_FALSE_PASS"
                and old_refined_rule_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_exact_required_form_statement_domain_admission_refined_exact_required_form_statement_shift_interface_adapter_or_carrier_identification_rule_target_v1"
                and old_refined_rule_refs.get("shift_to_refined_exact_required_form_statement_interface_target_ref")
                != interface_target.get("object_id")
                and f841.get("status")
                == "F841_EXECUTED_CURRENT_STRICT_ALPHA_S_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_SHIFT_TO_EXACT_REQUIRED_FORM_STATEMENT_ADAPTER_OR_CARRIER_IDENTIFICATION_RULE_TARGET_PACKET_NO_FALSE_PASS"
                and old_lawful_rule_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_exact_required_form_statement_domain_admission_exact_required_form_statement_shift_interface_adapter_or_carrier_identification_rule_target_v1"
                and old_lawful_rule_refs.get("shift_to_exact_required_form_statement_interface_target_ref")
                != interface_target.get("object_id")
                and f830.get("status")
                == "F830_EXECUTED_CURRENT_STRICT_ALPHA_S_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_SHIFT_TO_EXACT_REQUIRED_FORM_STATEMENT_ADAPTER_OR_CARRIER_IDENTIFICATION_RULE_TARGET_PACKET_NO_FALSE_PASS"
                and old_exact_rule_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_shift_to_exact_required_form_statement_adapter_or_carrier_identification_rule_target_v1"
                and old_exact_rule_refs.get("shift_to_exact_required_form_statement_interface_target_ref")
                != interface_target.get("object_id")
                and f810.get("status")
                == "F810_EXECUTED_CURRENT_STRICT_ALPHA_S_STRICT_SOURCE_SHANNON_SHIFT_TO_DOMAIN_CARRIER_IDENTIFICATION_OR_ADAPTER_RULE_TARGET_PACKET_NO_FALSE_PASS"
                and old_shannon_rule_target.get("object_id")
                == "alpha_s_strict_source_shannon_shift_to_domain_carrier_identification_or_adapter_rule_target_v1"
                and old_shannon_rule_refs.get("provider_shift_candidate_reference_lane_ref")
                != interface_refs.get("different_required_form_statement_provider_class_shift_candidate_reference_lane_ref")
                and f819.get("status")
                == "F819_EXECUTED_CURRENT_STRICT_ALPHA_S_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_SHIFT_TO_SCHEMA_ADAPTER_OR_CARRIER_IDENTIFICATION_RULE_TARGET_PACKET_NO_FALSE_PASS"
                and old_schema_rule_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_shift_to_schema_adapter_or_carrier_identification_rule_target_v1"
                and old_schema_rule_refs.get("shift_to_schema_interface_target_ref")
                != interface_target.get("object_id")
            ),
            "details": "The older F863, F852, F841, F830, F819, and F810 rule targets are frozen for different lanes and cannot silently fill the new lawful refined F873 slot.",
        },
        {
            "id": "no_current_export_names_exact_rule_for_f873_target",
            "pass": exact_rule_scan_hits == [],
            "details": "No current export names one exact adapter or carrier-identification rule for the frozen F873 lawful refined interface target.",
        },
    ]

    failed_checks = [item["id"] for item in checks if not item["pass"]]
    all_pass = all(item["pass"] for item in checks)

    theorem_result = {
        "exact_shift_interface_adapter_or_carrier_identification_rule_exported_for_f873_target": False
        if all_pass
        else None,
        "next_honest_move_requires_freezing_exact_shift_interface_adapter_or_carrier_identification_rule_target": all_pass,
        "no_false_pass": True,
    }

    status = (
        "P874_CURRENT_STRICT_ALPHA_S_NO_SHIFT_INTERFACE_ADAPTER_OR_CARRIER_RULE_FOR_PAIR12_PROVIDER_LAWFUL_REFINED_SIR_DA_EXACT_REQUIRED_FORM_STATEMENT_INTERFACE_TARGET_FREEZE_REQUIRED"
        if all_pass
        else "P874_REQUIRES_REVIEW"
    )

    artifact = {
        "stage": "P874",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "f873_shift_to_lawful_refined_shift_interface_rule_exact_required_form_statement_interface_target_packet": rel(IN_F873),
            "f872_lawful_refined_shift_interface_rule_different_required_form_statement_provider_class_shift_candidate_reference_packet": rel(IN_F872),
            "f863_older_current_rule_target_packet": rel(IN_F863),
            "f852_older_refined_rule_target_packet": rel(IN_F852),
            "f841_older_lawful_rule_target_packet": rel(IN_F841),
            "f830_older_exact_rule_target_packet": rel(IN_F830),
            "f819_older_schema_rule_target_packet": rel(IN_F819),
            "f810_older_shannon_rule_target_packet": rel(IN_F810),
            "p788_generic_alpha_s_adapter_probe": rel(IN_P788),
            "p764_own_lane_missing_interface_target_probe": rel(IN_P764),
        },
        "exact_missing_rule_target_candidate": {
            "candidate_id": "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_shift_interface_adapter_or_carrier_identification_rule_target_missing_v1",
            "shift_to_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_interface_target_ref": interface_target.get("object_id"),
            "different_required_form_statement_provider_class_shift_candidate_reference_lane_ref": interface_refs.get("different_required_form_statement_provider_class_shift_candidate_reference_lane_ref"),
            "repo_scan_hits_for_exact_rule": exact_rule_scan_hits,
        },
        "theorem_result": theorem_result,
        "checks": checks,
        "failed_checks": failed_checks,
        "current_honest_reading": [
            "F873 already freezes the exact lawful refined T213/T216 -> lawful refined alpha_s shift-interface-rule domain-admission exact-required-form-statement interface target.",
            "P874 shows that no current export names the adapter or carrier-identification rule required to instantiate that target.",
            "Therefore the next honest object is the exact missing lawful refined rule target itself, not a claim that the rule already exists.",
        ],
        "recommended_next_packet": {
            "id": "F874_CURRENT_STRICT_ALPHA_S_PAIR12_PROVIDER_SHIFT_TO_LAWFUL_REFINED_SIR_DA_EXACT_REQUIRED_FORM_STATEMENT_ADAPTER_OR_CARRIER_RULE_TARGET_PACKET",
            "goal": "Freeze the exact missing adapter-or-carrier-identification rule target required to instantiate the frozen F873 lawful refined interface target without silent domain identification.",
            "export_object_id": "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_shift_interface_adapter_or_carrier_identification_rule_target_v1",
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "P874",
        "status": status,
        "as_of": AS_OF,
        "exact_shift_interface_adapter_or_carrier_identification_rule_exported_for_f873_target": theorem_result[
            "exact_shift_interface_adapter_or_carrier_identification_rule_exported_for_f873_target"
        ],
        "recommended_next_packet": artifact["recommended_next_packet"]["id"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
