#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P814 = GENERATED / "p814_current_strict_alpha_s_no_selection_or_preference_rule_schema_for_f812_target_on_exported_f813_support_domain_narrow_schema_target_freeze_required.json"
IN_F812 = GENERATED / "f812_current_strict_alpha_s_strict_source_shannon_shift_to_domain_source_binding_selection_or_preference_rule_target_packet.json"
IN_F813 = GENERATED / "f813_current_strict_alpha_s_strict_source_shannon_source_binding_candidate_support_domain_packet.json"

OUT = GENERATED / "f814_current_strict_alpha_s_strict_source_shannon_source_binding_selection_or_preference_rule_schema_target_packet.json"
OUT_SUMMARY = GENERATED / "f814_current_strict_alpha_s_strict_source_shannon_source_binding_selection_or_preference_rule_schema_target_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P814, IN_F812, IN_F813]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "F814",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p814 = load_json(IN_P814)
    f812 = load_json(IN_F812)
    f813 = load_json(IN_F813)

    p814_theorem = p814.get("theorem_result") or {}
    f812_target = f812.get("target_object") or {}
    f813_export = f813.get("exported_object") or {}

    if (
        p814.get("status")
        == "P814_CURRENT_STRICT_ALPHA_S_NO_SELECTION_OR_PREFERENCE_RULE_SCHEMA_FOR_F812_TARGET_ON_EXPORTED_F813_SUPPORT_DOMAIN_NARROW_SCHEMA_TARGET_FREEZE_REQUIRED"
        and p814_theorem.get("candidate_source_support_domain_now_exactly_exported") is True
        and p814_theorem.get("exact_selection_or_preference_rule_schema_exported_on_f813_domain") is False
        and p814_theorem.get("next_honest_move_requires_freezing_exact_schema_only_target") is True
        and f812.get("status")
        == "F812_EXECUTED_CURRENT_STRICT_ALPHA_S_STRICT_SOURCE_SHANNON_SHIFT_TO_DOMAIN_SOURCE_BINDING_SELECTION_OR_PREFERENCE_RULE_TARGET_PACKET_NO_FALSE_PASS"
        and f813.get("status")
        == "F813_EXECUTED_CURRENT_STRICT_ALPHA_S_STRICT_SOURCE_SHANNON_SOURCE_BINDING_CANDIDATE_SUPPORT_DOMAIN_PACKET_NO_FALSE_PASS"
    ):
        status = "F814_EXECUTED_CURRENT_STRICT_ALPHA_S_STRICT_SOURCE_SHANNON_SOURCE_BINDING_SELECTION_OR_PREFERENCE_RULE_SCHEMA_TARGET_PACKET_NO_FALSE_PASS"
    else:
        status = "F814_REQUIRES_REVIEW"

    artifact = {
        "stage": "F814",
        "packet_name": "CurrentStrictAlphaSStrictSourceShannonSourceBindingSelectionOrPreferenceRuleSchemaTarget_v1",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p814_schema_absence_probe": rel(IN_P814),
            "f812_selection_or_preference_rule_target_packet": rel(IN_F812),
            "f813_candidate_support_domain_packet": rel(IN_F813),
        },
        "target_object": {
            "object_id": "alpha_s_strict_source_shannon_source_binding_selection_or_preference_rule_schema_target_v1",
            "goal": "Freeze the exact missing schema-only target now that the F813 candidate support domain is explicit but no lawful source-binding selection/preference schema is exported.",
            "required_fields": [
                {
                    "name": "source_binding_selection_or_preference_rule_target_ref",
                    "required": True,
                    "hard_limit": "Must point to the exact F812 source-binding selection/preference-rule target and not silently replace it.",
                },
                {
                    "name": "candidate_source_support_domain_ref",
                    "required": True,
                    "hard_limit": "Must point to the exact F813 finite two-support domain and not silently replace either member.",
                },
                {
                    "name": "selection_or_preference_rule_schema",
                    "required": True,
                    "hard_limit": "Must state how one exact source support is selected, preferred, or proven uniquely bindable on the exported F813 domain.",
                },
                {
                    "name": "domain_member_grade_handling_clause",
                    "required": True,
                    "hard_limit": "Must keep explicit that one support is current-candidate-only and the other is future-only entry target; no silent grade promotion is allowed.",
                },
                {
                    "name": "nontransfer_boundary_ref",
                    "required": True,
                    "hard_limit": "Must block silent reuse of foreign-domain selection templates and probe-local family-order rules.",
                },
                {
                    "name": "selected_source_binding_output_schema",
                    "required": True,
                    "hard_limit": "Must state what exact selected source-binding object or relation a future schema would output.",
                },
                {
                    "name": "future_route_grade_ref",
                    "required": True,
                    "hard_limit": "Must keep this schema at future-route target grade until a real schema is exported.",
                },
                {
                    "name": "hard_limits",
                    "required": True,
                    "hard_limit": "Must explicitly deny provider-shift success, QCD closure, and ToE closure.",
                },
            ],
        },
        "target_refs": {
            "source_binding_selection_or_preference_rule_target_ref": f812_target.get("object_id"),
            "candidate_source_support_domain_ref": f813_export.get("object_id"),
            "missing_target_candidate_id": "alpha_s_strict_source_shannon_source_binding_selection_or_preference_rule_schema_target_missing_v1",
        },
        "current_honest_reading": [
            "The repo now exports the exact candidate support domain for the strict-source Shannon -> alpha_s source-binding route.",
            "But no current export supplies the remaining selection_or_preference_rule_schema on that domain.",
            "F814 therefore freezes the exact narrower schema-only target and keeps actual source binding downstream.",
        ],
        "recommended_next_move": "Build one narrow probe testing whether any current export can lawfully supply selection_or_preference_rule_schema or domain_member_grade_handling_clause for alpha_s_strict_source_shannon_source_binding_selection_or_preference_rule_schema_target_v1 on the exported F813 support domain without silent domain identification.",
        "hard_limits": [
            "Does not claim that the selection/preference rule schema already exists.",
            "Does not claim that any source has already been selected.",
            "Does not claim that the F811 source binding already exists.",
            "Does not claim that any adapter action schema is already exported.",
            "Does not claim that foreign-domain selection templates are reusable here.",
            "Does not claim that provider shift has already succeeded.",
            "Does not claim alpha_s boundary export readiness.",
            "Does not claim QCD closure.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F814",
        "status": status,
        "as_of": AS_OF,
        "target_object_id": artifact["target_object"]["object_id"],
        "candidate_source_support_domain_ref": artifact["target_refs"]["candidate_source_support_domain_ref"],
        "recommended_next_move": artifact["recommended_next_move"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
