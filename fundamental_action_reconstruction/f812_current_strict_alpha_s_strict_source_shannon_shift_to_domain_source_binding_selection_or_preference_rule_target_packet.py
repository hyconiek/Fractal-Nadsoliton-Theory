#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P812 = GENERATED / "p812_current_strict_alpha_s_no_source_binding_selection_or_preference_rule_for_f811_target_exported_target_freeze_required.json"
IN_F811 = GENERATED / "f811_current_strict_alpha_s_strict_source_shannon_shift_to_domain_source_binding_target_packet.json"
IN_P811 = GENERATED / "p811_current_strict_alpha_s_source_support_present_exact_source_binding_unexported_adapter_action_schema_blocked_source_binding_target_freeze_required.json"

OUT = GENERATED / "f812_current_strict_alpha_s_strict_source_shannon_shift_to_domain_source_binding_selection_or_preference_rule_target_packet.json"
OUT_SUMMARY = GENERATED / "f812_current_strict_alpha_s_strict_source_shannon_shift_to_domain_source_binding_selection_or_preference_rule_target_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P812, IN_F811, IN_P811]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "F812",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p812 = load_json(IN_P812)
    f811 = load_json(IN_F811)
    p811 = load_json(IN_P811)

    if (
        p812.get("status")
        == "P812_CURRENT_STRICT_ALPHA_S_NO_SOURCE_BINDING_SELECTION_OR_PREFERENCE_RULE_FOR_F811_TARGET_EXPORTED_TARGET_FREEZE_REQUIRED"
        and (p812.get("theorem_result") or {}).get(
            "exact_source_binding_selection_or_preference_rule_exported_for_f811_target"
        )
        is False
        and (p812.get("theorem_result") or {}).get(
            "next_honest_move_requires_freezing_exact_source_binding_selection_or_preference_rule_target"
        )
        is True
        and f811.get("status")
        == "F811_EXECUTED_CURRENT_STRICT_ALPHA_S_STRICT_SOURCE_SHANNON_SHIFT_TO_DOMAIN_SOURCE_BINDING_TARGET_PACKET_NO_FALSE_PASS"
    ):
        status = "F812_EXECUTED_CURRENT_STRICT_ALPHA_S_STRICT_SOURCE_SHANNON_SHIFT_TO_DOMAIN_SOURCE_BINDING_SELECTION_OR_PREFERENCE_RULE_TARGET_PACKET_NO_FALSE_PASS"
    else:
        status = "F812_REQUIRES_REVIEW"

    p812_missing = p812.get("exact_missing_rule_target_candidate") or {}
    f811_target = f811.get("target_object") or {}
    p811_support = p811.get("support_objects") or {}

    artifact = {
        "stage": "F812",
        "packet_name": "CurrentStrictAlphaSStrictSourceShannonShiftToDomainSourceBindingSelectionOrPreferenceRuleTarget_v1",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p812_selection_or_preference_rule_audit_probe": rel(IN_P812),
            "f811_source_binding_target_packet": rel(IN_F811),
            "p811_source_support_probe": rel(IN_P811),
        },
        "target_object": {
            "object_id": "alpha_s_strict_source_shannon_shift_to_domain_source_binding_selection_or_preference_rule_target_v1",
            "goal": "Freeze the exact missing source-binding selection/preference-rule target required before one source can be lawfully chosen for the F811 route.",
            "required_fields": [
                {
                    "name": "source_binding_target_ref",
                    "required": True,
                    "hard_limit": "Must point to the exact F811 source-binding target and not silently replace it."
                },
                {
                    "name": "current_source_candidate_support_ref",
                    "required": True,
                    "hard_limit": "Must identify one real current strict-source Shannon candidate support object without promoting it into alpha_s entry."
                },
                {
                    "name": "lawful_future_entry_object_support_ref",
                    "required": True,
                    "hard_limit": "Must identify one lawful future-only entry-object support target without promoting it into actual entry."
                },
                {
                    "name": "candidate_source_support_domain_ref",
                    "required": True,
                    "hard_limit": "Must define the exact finite support domain on which the source-binding selection/preference rule acts."
                },
                {
                    "name": "selection_or_preference_rule_schema",
                    "required": True,
                    "hard_limit": "Must state how one exact source support is selected, preferred, or proven uniquely bindable for the F811 route."
                },
                {
                    "name": "nontransfer_boundary_ref",
                    "required": True,
                    "hard_limit": "Must block silent reuse of foreign-domain selection templates and probe-local family-order rules."
                },
                {
                    "name": "selected_source_binding_output_schema",
                    "required": True,
                    "hard_limit": "Must state what exact selected source binding object or relation the future rule would output."
                },
                {
                    "name": "downstream_selected_source_binding_schema_target_ref",
                    "required": True,
                    "hard_limit": "Must keep exact source-binding schema downstream until the rule is exported."
                },
                {
                    "name": "future_route_grade_ref",
                    "required": True,
                    "hard_limit": "Must keep this object at future-route target grade until a real rule is exported."
                },
                {
                    "name": "hard_limits",
                    "required": True,
                    "hard_limit": "Must explicitly deny provider-shift success, QCD closure, and ToE closure."
                },
            ],
        },
        "target_refs": {
            "source_binding_target_ref": f811_target.get("object_id"),
            "current_source_candidate_support_ref": p811_support.get("current_source_candidate_support_ref"),
            "lawful_future_entry_object_support_ref": p811_support.get("lawful_future_entry_object_support_ref"),
            "missing_target_candidate_id": p812_missing.get("candidate_id"),
        },
        "current_honest_reading": [
            "The repo now freezes the exact source-binding target and already has the source-side supports on hand.",
            "But no current export lawfully chooses between those supports for the F811 route.",
            "F812 therefore freezes the exact missing source-binding selection/preference-rule target and keeps actual source binding downstream.",
        ],
        "recommended_next_move": "Build one narrow probe testing whether any current export can lawfully supply candidate_source_support_domain_ref or selection_or_preference_rule_schema for alpha_s_strict_source_shannon_shift_to_domain_source_binding_selection_or_preference_rule_target_v1 without silent domain identification.",
        "hard_limits": [
            "Does not claim that the source-binding selection/preference rule already exists.",
            "Does not claim that the F811 source binding already exists.",
            "Does not claim that any exact source has already been chosen.",
            "Does not claim that any adapter action schema is already exported.",
            "Does not claim that provider shift has already succeeded.",
            "Does not claim alpha_s boundary export readiness.",
            "Does not claim QCD closure.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F812",
        "status": status,
        "as_of": AS_OF,
        "target_object_id": artifact["target_object"]["object_id"],
        "source_binding_target_ref": artifact["target_refs"]["source_binding_target_ref"],
        "recommended_next_move": artifact["recommended_next_move"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
