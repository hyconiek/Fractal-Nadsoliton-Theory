#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P811 = GENERATED / "p811_current_strict_alpha_s_source_support_present_exact_source_binding_unexported_adapter_action_schema_blocked_source_binding_target_freeze_required.json"
IN_F810 = GENERATED / "f810_current_strict_alpha_s_strict_source_shannon_shift_to_domain_carrier_identification_or_adapter_rule_target_packet.json"
IN_P755 = GENERATED / "p755_current_strict_t209_t26_component2_minimal_designated_pair12_noncyclic_entry_object_target_probe.json"

OUT = GENERATED / "f811_current_strict_alpha_s_strict_source_shannon_shift_to_domain_source_binding_target_packet.json"
OUT_SUMMARY = GENERATED / "f811_current_strict_alpha_s_strict_source_shannon_shift_to_domain_source_binding_target_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P811, IN_F810, IN_P755]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "F811",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p811 = load_json(IN_P811)
    f810 = load_json(IN_F810)
    p755 = load_json(IN_P755)

    if (
        p811.get("status")
        == "P811_CURRENT_STRICT_ALPHA_S_SOURCE_SUPPORT_PRESENT_EXACT_SOURCE_BINDING_UNEXPORTED_ADAPTER_ACTION_SCHEMA_BLOCKED_SOURCE_BINDING_TARGET_FREEZE_REQUIRED"
        and (p811.get("theorem_result") or {}).get("exact_source_candidate_lane_or_entry_ref_exported_for_f810_target")
        is False
        and (p811.get("theorem_result") or {}).get("next_honest_move_requires_freezing_exact_source_binding_target")
        is True
        and f810.get("status")
        == "F810_EXECUTED_CURRENT_STRICT_ALPHA_S_STRICT_SOURCE_SHANNON_SHIFT_TO_DOMAIN_CARRIER_IDENTIFICATION_OR_ADAPTER_RULE_TARGET_PACKET_NO_FALSE_PASS"
    ):
        status = "F811_EXECUTED_CURRENT_STRICT_ALPHA_S_STRICT_SOURCE_SHANNON_SHIFT_TO_DOMAIN_SOURCE_BINDING_TARGET_PACKET_NO_FALSE_PASS"
    else:
        status = "F811_REQUIRES_REVIEW"

    p811_support = p811.get("support_objects") or {}
    f810_target = f810.get("target_object") or {}
    f810_refs = f810.get("target_refs") or {}
    p755_theorem = p755.get("theorem_result") or {}

    artifact = {
        "stage": "F811",
        "packet_name": "CurrentStrictAlphaSStrictSourceShannonShiftToDomainSourceBindingTarget_v1",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p811_source_binding_audit_probe": rel(IN_P811),
            "f810_rule_target_packet": rel(IN_F810),
            "p755_future_entry_object_target_probe": rel(IN_P755),
        },
        "target_object": {
            "object_id": "alpha_s_strict_source_shannon_shift_to_domain_source_binding_target_v1",
            "goal": "Freeze the exact source-binding target required before any later adapter action schema could lawfully act on the F810 rule target.",
            "required_fields": [
                {
                    "name": "carrier_identification_or_adapter_rule_target_ref",
                    "required": True,
                    "hard_limit": "Must point to the exact F810 rule target and not silently replace it."
                },
                {
                    "name": "current_source_candidate_support_ref",
                    "required": True,
                    "hard_limit": "Must identify one real current strict-source Shannon candidate support object without promoting it into alpha_s entry."
                },
                {
                    "name": "lawful_future_entry_object_support_ref",
                    "required": True,
                    "hard_limit": "Must identify one lawful future-only source-side entry-object support target without promoting it into actual entry."
                },
                {
                    "name": "source_binding_selection_or_preference_rule_ref",
                    "required": True,
                    "hard_limit": "Must state how one exact source support is selected or bound for the F810 route."
                },
                {
                    "name": "selected_source_binding_schema",
                    "required": True,
                    "hard_limit": "Must state what exact binding relation is exported between the chosen source object and the F810 route."
                },
                {
                    "name": "downstream_adapter_action_schema_target_ref",
                    "required": True,
                    "hard_limit": "Must keep adapter action schema downstream until exact source binding is frozen."
                },
                {
                    "name": "same_domain_admission_or_nonidentification_boundary_ref",
                    "required": True,
                    "hard_limit": "Must explicitly keep silent alpha_s-domain identification forbidden."
                },
                {
                    "name": "future_route_grade_ref",
                    "required": True,
                    "hard_limit": "Must keep this object at future-route target grade until a real source binding is exported."
                },
                {
                    "name": "hard_limits",
                    "required": True,
                    "hard_limit": "Must explicitly deny provider-shift success, QCD closure, and ToE closure."
                },
            ],
        },
        "target_refs": {
            "carrier_identification_or_adapter_rule_target_ref": f810_target.get("object_id"),
            "current_source_candidate_support_ref": p811_support.get("current_source_candidate_support_ref"),
            "lawful_future_entry_object_support_ref": p811_support.get("lawful_future_entry_object_support_ref") or p755_theorem.get("t209_target_name"),
            "shift_to_domain_interface_target_ref": f810_refs.get("shift_to_domain_interface_target_ref"),
        },
        "current_honest_reading": [
            "Source-side Shannon support is real on the current repo state.",
            "But no current export binds one exact source object to the F810 rule target.",
            "F811 therefore freezes the exact missing source-binding target and keeps adapter action schema downstream.",
        ],
        "recommended_next_move": "Build one narrow probe testing whether any current export can lawfully supply source_binding_selection_or_preference_rule_ref for alpha_s_strict_source_shannon_shift_to_domain_source_binding_target_v1 without silent domain identification.",
        "hard_limits": [
            "Does not claim that the source binding already exists.",
            "Does not claim that the F810 rule target is already realized.",
            "Does not claim that T209 already gives actual entry.",
            "Does not claim that any adapter action schema is already exported.",
            "Does not claim that provider shift has already succeeded.",
            "Does not claim alpha_s boundary export readiness.",
            "Does not claim QCD closure.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F811",
        "status": status,
        "as_of": AS_OF,
        "target_object_id": artifact["target_object"]["object_id"],
        "carrier_identification_or_adapter_rule_target_ref": artifact["target_refs"][
            "carrier_identification_or_adapter_rule_target_ref"
        ],
        "recommended_next_move": artifact["recommended_next_move"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
