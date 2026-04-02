#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-21"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_F939 = GENERATED / "f939_current_strict_alpha_s_pair12_provider_shift_to_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_interface_target_packet.json"
IN_F938 = GENERATED / "f938_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_diff_rfs_provider_shift_candidate_reference_packet.json"
IN_F929 = GENERATED / "f929_current_strict_alpha_s_pair12_provider_shift_to_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_adapter_or_carrier_rule_target_packet.json"
IN_F918 = GENERATED / "f918_current_strict_alpha_s_pair12_provider_shift_to_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_adapter_or_carrier_rule_target_packet.json"
IN_F907 = GENERATED / "f907_current_strict_alpha_s_pair12_provider_shift_to_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_domain_admission_adapter_or_carrier_rule_target_packet.json"
IN_F896 = GENERATED / "f896_current_strict_alpha_s_pair12_provider_shift_to_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_adapter_or_carrier_rule_target_packet.json"
IN_F885 = GENERATED / "f885_current_strict_alpha_s_pair12_provider_shift_to_lawful_refined_sir_da_exact_required_form_statement_da_exact_required_form_statement_adapter_or_carrier_rule_target_packet.json"
IN_F874 = GENERATED / "f874_current_strict_alpha_s_pair12_provider_shift_to_lawful_refined_sir_da_exact_required_form_statement_adapter_or_carrier_rule_target_packet.json"
IN_F863 = GENERATED / "f863_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_lawful_refined_exact_required_form_statement_adapter_or_carrier_identification_rule_target_packet.json"
IN_F852 = GENERATED / "f852_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_refined_exact_required_form_statement_adapter_or_carrier_identification_rule_target_packet.json"
IN_F841 = GENERATED / "f841_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_exact_required_form_statement_adapter_or_carrier_identification_rule_target_packet.json"
IN_F830 = GENERATED / "f830_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_exact_required_form_statement_adapter_or_carrier_identification_rule_target_packet.json"
IN_F819 = GENERATED / "f819_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_schema_adapter_or_carrier_identification_rule_target_packet.json"
IN_F810 = GENERATED / "f810_current_strict_alpha_s_strict_source_shannon_shift_to_domain_carrier_identification_or_adapter_rule_target_packet.json"
IN_P788 = GENERATED / "p788_current_alpha_s_dimensionless_or_normalized_replacement_route_probe.json"
IN_P764 = GENERATED / "p764_current_strict_t218_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_target_probe.json"

OUT = GENERATED / "p940_current_strict_alpha_s_no_shift_interface_adapter_or_carrier_rule_for_pair12_provider_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_interface_target_freeze_required.json"
OUT_SUMMARY = GENERATED / "p940_current_strict_alpha_s_no_shift_interface_adapter_or_carrier_rule_for_pair12_provider_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_interface_target_freeze_required_summary.json"


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


def has_required_field(obj: dict[str, Any], name: str) -> bool:
    return any(
        isinstance(item, dict) and item.get("name") == name
        for item in (obj.get("required_fields") or [])
    )


def find_exact_rule_hits(interface_target_id: str, rule_terms: list[str]) -> list[str]:
    hits: list[str] = []
    excluded_suffixes = [
        "generated/p939_current_strict_alpha_s_no_exact_pair12_provider_shift_to_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_interface_target_freeze_required.json",
        "generated/p939_current_strict_alpha_s_no_exact_pair12_provider_shift_to_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_interface_target_freeze_required_summary.json",
        "generated/f939_current_strict_alpha_s_pair12_provider_shift_to_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_interface_target_packet.json",
        "generated/f939_current_strict_alpha_s_pair12_provider_shift_to_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_interface_target_packet_summary.json",
        "P939_CURRENT_STRICT_ALPHA_S_NO_EXACT_PAIR12_PROVIDER_SHIFT_TO_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_DA_ERFS_INTERFACE_TARGET_FREEZE_REQUIRED.md",
        "F939_CURRENT_STRICT_ALPHA_S_PAIR12_PROVIDER_SHIFT_TO_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_DA_ERFS_INTERFACE_TARGET_PACKET.md",
        "p939_current_strict_alpha_s_no_exact_pair12_provider_shift_to_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_interface_target_freeze_required.py",
        "f939_current_strict_alpha_s_pair12_provider_shift_to_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_interface_target_packet.py",
        "generated/p940_current_strict_alpha_s_no_shift_interface_adapter_or_carrier_rule_for_pair12_provider_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_interface_target_freeze_required.json",
        "generated/p940_current_strict_alpha_s_no_shift_interface_adapter_or_carrier_rule_for_pair12_provider_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_interface_target_freeze_required_summary.json",
        "generated/f940_current_strict_alpha_s_pair12_provider_shift_to_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_adapter_or_carrier_rule_target_packet.json",
        "generated/f940_current_strict_alpha_s_pair12_provider_shift_to_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_adapter_or_carrier_rule_target_packet_summary.json",
        "P940_CURRENT_STRICT_ALPHA_S_NO_SHIFT_INTERFACE_ADAPTER_OR_CARRIER_RULE_FOR_PAIR12_PROVIDER_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_DA_ERFS_INTERFACE_TARGET_FREEZE_REQUIRED.md",
        "F940_CURRENT_STRICT_ALPHA_S_PAIR12_PROVIDER_SHIFT_TO_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_DA_ERFS_ADAPTER_OR_CARRIER_RULE_TARGET_PACKET.md",
        "p940_current_strict_alpha_s_no_shift_interface_adapter_or_carrier_rule_for_pair12_provider_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_interface_target_freeze_required.py",
        "f940_current_strict_alpha_s_pair12_provider_shift_to_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_adapter_or_carrier_rule_target_packet.py",
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

    prereq = [
        IN_F939,
        IN_F938,
        IN_F929,
        IN_F918,
        IN_F907,
        IN_F896,
        IN_F885,
        IN_F874,
        IN_F863,
        IN_F852,
        IN_F841,
        IN_F830,
        IN_F819,
        IN_F810,
        IN_P788,
        IN_P764,
    ]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P940",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    f939 = load_json(IN_F939)
    f938 = load_json(IN_F938)
    f929 = load_json(IN_F929)
    f918 = load_json(IN_F918)
    f907 = load_json(IN_F907)
    f896 = load_json(IN_F896)
    f885 = load_json(IN_F885)
    f874 = load_json(IN_F874)
    f863 = load_json(IN_F863)
    f852 = load_json(IN_F852)
    f841 = load_json(IN_F841)
    f830 = load_json(IN_F830)
    f819 = load_json(IN_F819)
    f810 = load_json(IN_F810)
    p788 = load_json(IN_P788)
    p764 = load_json(IN_P764)

    interface_target = f939.get("target_object") or {}
    interface_refs = f939.get("target_refs") or {}
    candidate_lane = f938.get("exported_object") or {}
    p764_theorem = p764.get("theorem_result") or {}

    previous_packets = [
        ("f929", f929, "target_object", "target_refs"),
        ("f918", f918, "target_object", "target_refs"),
        ("f907", f907, "target_object", "target_refs"),
        ("f896", f896, "target_object", "target_refs"),
        ("f885", f885, "target_object", "target_refs"),
        ("f874", f874, "target_object", "target_refs"),
        ("f863", f863, "target_object", "target_refs"),
        ("f852", f852, "target_object", "target_refs"),
        ("f841", f841, "target_object", "target_refs"),
        ("f830", f830, "target_object", "target_refs"),
        ("f819", f819, "target_object", "target_refs"),
        ("f810", f810, "target_object", "target_refs"),
    ]

    current_interface_target_id = (
        "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_"
        "domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_"
        "domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_"
        "domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_"
        "domain_admission_exact_required_form_statement_shift_interface_target_v1"
    )

    exact_rule_scan_hits = find_exact_rule_hits(
        current_interface_target_id,
        [
            "adapter",
            "carrier identification",
            "adapter action schema",
            "carrier safe",
            "domain admission rule",
        ],
    )

    current_candidate_lane_ref = interface_refs.get(
        "different_required_form_statement_provider_class_shift_candidate_reference_lane_ref"
    )
    current_downstream_ref = interface_refs.get(
        "downstream_lawful_refined_deeper_exact_required_form_statement_domain_admission_target_ref"
    )

    previous_targets_nonreusable = True
    for _, packet, object_key, refs_key in previous_packets:
        obj = packet.get(object_key) or {}
        refs = packet.get(refs_key) or {}
        status = packet.get("status")
        if not (isinstance(status, str) and status.endswith("_NO_FALSE_PASS")):
            previous_targets_nonreusable = False
            break
        if obj.get("object_id") == current_interface_target_id:
            previous_targets_nonreusable = False
            break
        for key, value in refs.items():
            if key.endswith("interface_target_ref") and value == current_interface_target_id:
                previous_targets_nonreusable = False
                break
            if key.endswith("domain_admission_target_ref") and value == current_downstream_ref:
                previous_targets_nonreusable = False
                break
            if key == "provider_shift_candidate_reference_lane_ref" and value == current_candidate_lane_ref:
                previous_targets_nonreusable = False
                break
            if key == "different_required_form_statement_provider_class_shift_candidate_reference_lane_ref" and value == current_candidate_lane_ref:
                previous_targets_nonreusable = False
                break
        if not previous_targets_nonreusable:
            break

    checks = [
        {
            "id": "f939_freezes_current_deeper_interface_target_and_requires_adapter_or_carrier_rule",
            "pass": (
                f939.get("status")
                == "F939_EXECUTED_CURRENT_STRICT_ALPHA_S_PAIR12_PROVIDER_SHIFT_TO_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_INTERFACE_TARGET_PACKET_NO_FALSE_PASS"
                and interface_target.get("object_id") == current_interface_target_id
                and has_required_field(interface_target, "shift_interface_adapter_or_carrier_identification_rule_ref")
            ),
            "details": "F939 already freezes the current deeper interface target and names the still-missing adapter-or-carrier-identification rule slot.",
        },
        {
            "id": "f938_still_keeps_t213_t216_lane_at_candidate_reference_only_grade",
            "pass": (
                f938.get("status")
                == "F938_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_DIFFERENT_REQUIRED_FORM_STATEMENT_PROVIDER_CLASS_SHIFT_CANDIDATE_REFERENCE_PACKET_NO_FALSE_PASS"
                and candidate_lane.get("candidate_reference_lane_grade") == "reference_context_candidate_only"
                and candidate_lane.get("alpha_s_shift_interface_status") == "blocked_nonexport"
                and candidate_lane.get("provider_class_shift_realization_status") == "not_realized"
                and candidate_lane.get("exact_required_form_statement_status") == "blocked_nonexport"
            ),
            "details": "The T213/T216 lane still remains candidate-reference only and does not realize shift.",
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
            "id": "older_rule_targets_remain_lane_specific_and_not_silently_reusable",
            "pass": previous_targets_nonreusable,
            "details": "All older rule targets remain lane-specific and cannot silently instantiate the new deeper lawful refined interface target.",
        },
        {
            "id": "repo_scan_finds_no_exact_rule_for_the_f939_interface_target",
            "pass": exact_rule_scan_hits == [],
            "details": "Repo scan finds no exact rule export tied to the frozen F939 interface target.",
        },
    ]

    failed_checks = [item["id"] for item in checks if not item["pass"]]
    all_pass = all(item["pass"] for item in checks)

    theorem_result = {
        "exact_shift_interface_adapter_or_carrier_identification_rule_exported_for_f939_target": False
        if all_pass
        else None,
        "next_honest_move_requires_freezing_exact_shift_interface_adapter_or_carrier_identification_rule_target": all_pass,
        "no_false_pass": True,
    }

    status = (
        "P940_CURRENT_STRICT_ALPHA_S_NO_SHIFT_INTERFACE_ADAPTER_OR_CARRIER_RULE_FOR_PAIR12_PROVIDER_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_INTERFACE_TARGET_FREEZE_REQUIRED"
        if all_pass
        else "P940_REQUIRES_REVIEW"
    )

    artifact = {
        "stage": "P940",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "f939_current_deeper_interface_target_packet": rel(IN_F939),
            "f938_current_deeper_candidate_reference_packet": rel(IN_F938),
            "f929_old_adjacent_rule_target_packet": rel(IN_F929),
            "f918_old_shallower_rule_target_packet": rel(IN_F918),
            "f907_old_current_rule_target_packet": rel(IN_F907),
            "f896_old_refined_rule_target_packet": rel(IN_F896),
            "f885_old_lawful_rule_target_packet": rel(IN_F885),
            "f874_old_exact_rule_target_packet": rel(IN_F874),
            "f863_old_schema_rule_target_packet": rel(IN_F863),
            "f852_old_refined_exact_rule_target_packet": rel(IN_F852),
            "f841_old_shannon_rule_target_packet": rel(IN_F841),
            "f830_oldest_shannon_rule_target_packet": rel(IN_F830),
            "f819_oldest_domain_rule_target_packet": rel(IN_F819),
            "f810_oldest_gate_rule_target_packet": rel(IN_F810),
            "p788_alpha_s_adapter_probe": rel(IN_P788),
            "p764_own_lane_missing_interface_target_probe": rel(IN_P764),
        },
        "exact_missing_rule_target_candidate": {
            "candidate_id": "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_shift_interface_adapter_or_carrier_identification_rule_target_missing_v1",
            "shift_interface_target_ref": interface_target.get("object_id"),
            "different_provider_class_shift_candidate_reference_lane_ref": candidate_lane.get("object_id"),
            "repo_scan_hits_for_exact_rule": exact_rule_scan_hits,
        },
        "theorem_result": theorem_result,
        "checks": checks,
        "failed_checks": failed_checks,
        "current_honest_reading": [
            "F939 already freezes the exact deeper lawful refined interface target between the admitted T213/T216 candidate lane and the current alpha_s problem.",
            "P940 shows that no current export names the adapter or carrier-identification rule required to instantiate that target.",
            "P940 therefore freezes the exact missing rule target without claiming that the rule already exists.",
        ],
        "recommended_next_packet": {
            "id": "F940_CURRENT_STRICT_ALPHA_S_PAIR12_PROVIDER_SHIFT_TO_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_ADAPTER_OR_CARRIER_RULE_TARGET_PACKET",
            "goal": "Freeze the exact missing adapter-or-carrier-identification rule target required to instantiate the frozen F939 interface target without silent domain identification.",
            "export_object_id": "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_shift_interface_adapter_or_carrier_identification_rule_target_v1",
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "P940",
        "status": status,
        "as_of": AS_OF,
        "exact_shift_interface_adapter_or_carrier_identification_rule_exported_for_f939_target": theorem_result[
            "exact_shift_interface_adapter_or_carrier_identification_rule_exported_for_f939_target"
        ],
        "recommended_next_packet": artifact["recommended_next_packet"]["id"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
