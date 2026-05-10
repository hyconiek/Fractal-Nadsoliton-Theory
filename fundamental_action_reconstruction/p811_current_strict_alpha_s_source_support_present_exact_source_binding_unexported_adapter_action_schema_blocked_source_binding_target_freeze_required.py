#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_F810 = GENERATED / "f810_current_strict_alpha_s_strict_source_shannon_shift_to_domain_carrier_identification_or_adapter_rule_target_packet.json"
IN_P755 = GENERATED / "p755_current_strict_t209_t26_component2_minimal_designated_pair12_noncyclic_entry_object_target_probe.json"
IN_P756 = GENERATED / "p756_current_strict_t210_t26_component2_minimal_designated_pair12_noncyclic_entry_object_actual_realization_nonexport_audit_probe.json"
IN_F321 = ROOT / "F321_FIRST_ACTUAL_STRICT_NAD12_SIGMA_STRICT_SOURCE_SHANNON_WEIGHTED_PAIR_POPULATION_REFINEMENT_CANDIDATE_PACKET.md"

OUT = GENERATED / "p811_current_strict_alpha_s_source_support_present_exact_source_binding_unexported_adapter_action_schema_blocked_source_binding_target_freeze_required.json"
OUT_SUMMARY = GENERATED / "p811_current_strict_alpha_s_source_support_present_exact_source_binding_unexported_adapter_action_schema_blocked_source_binding_target_freeze_required_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


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


def find_exact_source_binding_hits(target_id: str, source_terms: list[str]) -> list[str]:
    hits: list[str] = []
    excluded_suffixes = [
        "generated/f810_current_strict_alpha_s_strict_source_shannon_shift_to_domain_carrier_identification_or_adapter_rule_target_packet.json",
        "generated/f810_current_strict_alpha_s_strict_source_shannon_shift_to_domain_carrier_identification_or_adapter_rule_target_packet_summary.json",
        "generated/p810_current_strict_alpha_s_strict_source_shannon_shift_to_domain_carrier_identification_or_adapter_rule_audit_probe.json",
        "generated/p810_current_strict_alpha_s_strict_source_shannon_shift_to_domain_carrier_identification_or_adapter_rule_audit_probe_summary.json",
        "F810_CURRENT_STRICT_ALPHA_S_STRICT_SOURCE_SHANNON_SHIFT_TO_DOMAIN_CARRIER_IDENTIFICATION_OR_ADAPTER_RULE_TARGET_PACKET.md",
        "P810_CURRENT_STRICT_ALPHA_S_STRICT_SOURCE_SHANNON_SHIFT_TO_DOMAIN_CARRIER_IDENTIFICATION_OR_ADAPTER_RULE_AUDIT_PROBE.md",
        "f810_current_strict_alpha_s_strict_source_shannon_shift_to_domain_carrier_identification_or_adapter_rule_target_packet.py",
        "p810_current_strict_alpha_s_strict_source_shannon_shift_to_domain_carrier_identification_or_adapter_rule_audit_probe.py",
        "generated/p811_current_strict_alpha_s_source_support_present_exact_source_binding_unexported_adapter_action_schema_blocked_source_binding_target_freeze_required.json",
        "generated/p811_current_strict_alpha_s_source_support_present_exact_source_binding_unexported_adapter_action_schema_blocked_source_binding_target_freeze_required_summary.json",
        "generated/f811_current_strict_alpha_s_strict_source_shannon_shift_to_domain_source_binding_target_packet.json",
        "generated/f811_current_strict_alpha_s_strict_source_shannon_shift_to_domain_source_binding_target_packet_summary.json",
        "P811_CURRENT_STRICT_ALPHA_S_SOURCE_SUPPORT_PRESENT_EXACT_SOURCE_BINDING_UNEXPORTED_ADAPTER_ACTION_SCHEMA_BLOCKED_SOURCE_BINDING_TARGET_FREEZE_REQUIRED.md",
        "F811_CURRENT_STRICT_ALPHA_S_STRICT_SOURCE_SHANNON_SHIFT_TO_DOMAIN_SOURCE_BINDING_TARGET_PACKET.md",
        "p811_current_strict_alpha_s_source_support_present_exact_source_binding_unexported_adapter_action_schema_blocked_source_binding_target_freeze_required.py",
        "f811_current_strict_alpha_s_strict_source_shannon_shift_to_domain_source_binding_target_packet.py",
    ]
    norm_target = normalize(target_id)
    norm_terms = [normalize(term) for term in source_terms]
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
        if norm_target in hay and any(term in hay for term in norm_terms):
            hits.append(relpath)
    return sorted(hits)


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_F810, IN_P755, IN_P756, IN_F321]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P811",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    f810 = load_json(IN_F810)
    p755 = load_json(IN_P755)
    p756 = load_json(IN_P756)
    f321_text = load_text(IN_F321)

    target_object = f810.get("target_object") or {}
    target_refs = f810.get("target_refs") or {}
    p755_theorem = p755.get("theorem_result") or {}

    exact_source_binding_hits = find_exact_source_binding_hits(
        "alpha_s_strict_source_shannon_shift_to_domain_carrier_identification_or_adapter_rule_target_v1",
        [
            "source candidate lane or entry ref",
            "source binding",
            "entry object",
            "pair population refinement candidate",
        ],
    )

    f321_norm = normalize(f321_text)
    p756_theorem = p756.get("theorem_result") or {}

    checks = [
        {
            "id": "f810_requires_exact_source_field_and_adapter_action_schema",
            "pass": (
                f810.get("status")
                == "F810_EXECUTED_CURRENT_STRICT_ALPHA_S_STRICT_SOURCE_SHANNON_SHIFT_TO_DOMAIN_CARRIER_IDENTIFICATION_OR_ADAPTER_RULE_TARGET_PACKET_NO_FALSE_PASS"
                and target_object.get("object_id")
                == "alpha_s_strict_source_shannon_shift_to_domain_carrier_identification_or_adapter_rule_target_v1"
                and any(
                    isinstance(item, dict) and item.get("name") == "source_candidate_lane_or_entry_ref"
                    for item in (target_object.get("required_fields") or [])
                )
                and any(
                    isinstance(item, dict) and item.get("name") == "adapter_action_schema"
                    for item in (target_object.get("required_fields") or [])
                )
            ),
            "details": "F810 already freezes the exact missing rule target and requires both source_candidate_lane_or_entry_ref and adapter_action_schema.",
        },
        {
            "id": "f321_exports_real_current_strict_source_candidate_support",
            "pass": (
                "f321 current actual strict nad12 sigma strict source shannon weighted pair population refinement candidate packet no false pass"
                in f321_norm
                and "basispair strict nad12 sigma residual nonequality shannon weighted population refinement candidate v1"
                in f321_norm
                and "still below actual pair population" in f321_norm
            ),
            "details": "F321 exports one real current strict-source Shannon pair-population refinement candidate, but only at candidate-only grade.",
        },
        {
            "id": "p755_exports_lawful_future_entry_object_support",
            "pass": (
                p755.get("status") == "PASS_STRICT_T26_COMPONENT2_MINIMAL_DESIGNATED_PAIR12_NONCYCLIC_ENTRY_OBJECT_TARGET_EXPORTED"
                and p755_theorem.get("t209_target_exported_on_current_repo_state") is True
                and p755_theorem.get("current_t209_target_is_future_only") is True
                and p755_theorem.get("current_t209_target_is_source_side_observer_free") is True
                and p755_theorem.get("current_t209_target_is_kobs_independent_and_kernel_split_safe") is True
            ),
            "details": "P755/T209 export one lawful future-only source-side entry-object target that stays observer-free and kernel-split-safe.",
        },
        {
            "id": "p756_keeps_future_entry_object_below_actual_realization",
            "pass": (
                p756.get("status")
                == "PASS_STRICT_T26_COMPONENT2_MINIMAL_DESIGNATED_PAIR12_NONCYCLIC_ENTRY_OBJECT_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
                and p756_theorem.get("t210_target_exported_on_current_repo_state") is False
                and p756_theorem.get("current_repo_still_does_not_export_actual_realization_of_t209_target")
                is True
            ),
            "details": "P756 keeps the T209 future entry object below actual realization on current repo state.",
        },
        {
            "id": "no_current_export_binds_one_exact_source_object_to_f810_target",
            "pass": exact_source_binding_hits == [],
            "details": "No current export binds one exact source-side candidate or future entry object to the F810 rule target.",
        },
        {
            "id": "adapter_action_schema_still_blocked_downstream",
            "pass": (
                target_refs.get("shift_to_domain_interface_target_ref")
                == "alpha_s_strict_source_shannon_shift_to_domain_interface_target_v1"
                and any(
                    isinstance(item, dict) and item.get("name") == "adapter_action_schema"
                    for item in (target_object.get("required_fields") or [])
                )
            ),
            "details": "Adapter action schema remains a downstream blocked field and is not supplied by the current source-side supports.",
        },
    ]

    failed_checks = [item["id"] for item in checks if not item["pass"]]
    all_pass = all(item["pass"] for item in checks)

    theorem_result = {
        "current_source_candidate_support_present": checks[1]["pass"],
        "lawful_future_entry_object_support_present": checks[2]["pass"],
        "exact_source_candidate_lane_or_entry_ref_exported_for_f810_target": False if all_pass else None,
        "exact_adapter_action_schema_exported_for_f810_target": False if all_pass else None,
        "next_honest_move_requires_freezing_exact_source_binding_target": all_pass,
        "no_false_pass": True,
    }

    status = (
        "P811_CURRENT_STRICT_ALPHA_S_SOURCE_SUPPORT_PRESENT_EXACT_SOURCE_BINDING_UNEXPORTED_ADAPTER_ACTION_SCHEMA_BLOCKED_SOURCE_BINDING_TARGET_FREEZE_REQUIRED"
        if all_pass
        else "P811_REQUIRES_REVIEW"
    )

    artifact = {
        "stage": "P811",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "f810_rule_target_packet": rel(IN_F810),
            "f321_current_source_candidate_packet": rel(IN_F321),
            "p755_future_entry_object_target_probe": rel(IN_P755),
            "p756_future_entry_object_actual_realization_nonexport_audit": rel(IN_P756),
        },
        "support_objects": {
            "current_source_candidate_support_ref": "BasisPair_strict_nad12_sigma_residual_nonequality_shannon_weighted_population_refinement_candidate_v1",
            "lawful_future_entry_object_support_ref": p755_theorem.get("t209_target_name"),
            "repo_scan_hits_for_exact_source_binding": exact_source_binding_hits,
        },
        "theorem_result": theorem_result,
        "checks": checks,
        "failed_checks": failed_checks,
        "current_honest_reading": [
            "The current repo does have real strict-source Shannon support objects on the source side.",
            "But it still does not bind one exact source object to the frozen F810 rule target.",
            "So the next honest move is to freeze the exact source-binding target before pretending that any adapter action schema can act on a selected source.",
        ],
        "recommended_next_packet": {
            "id": "F811_CURRENT_STRICT_ALPHA_S_STRICT_SOURCE_SHANNON_SHIFT_TO_DOMAIN_SOURCE_BINDING_TARGET_PACKET",
            "goal": "Freeze the exact source-binding target required before any later adapter action schema could lawfully act on the F810 rule target.",
            "export_object_id": "alpha_s_strict_source_shannon_shift_to_domain_source_binding_target_v1",
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "P811",
        "status": status,
        "as_of": AS_OF,
        "current_source_candidate_support_present": theorem_result["current_source_candidate_support_present"],
        "lawful_future_entry_object_support_present": theorem_result["lawful_future_entry_object_support_present"],
        "exact_source_candidate_lane_or_entry_ref_exported_for_f810_target": theorem_result[
            "exact_source_candidate_lane_or_entry_ref_exported_for_f810_target"
        ],
        "exact_adapter_action_schema_exported_for_f810_target": theorem_result[
            "exact_adapter_action_schema_exported_for_f810_target"
        ],
        "next_honest_move_requires_freezing_exact_source_binding_target": theorem_result[
            "next_honest_move_requires_freezing_exact_source_binding_target"
        ],
        "recommended_next_packet": artifact["recommended_next_packet"]["id"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
