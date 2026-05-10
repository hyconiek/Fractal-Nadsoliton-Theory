#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P816 = GENERATED / "p816_current_strict_alpha_s_selection_or_preference_rule_schema_same_lane_exhaustion_boundary_audit_probe.json"
IN_F814 = GENERATED / "f814_current_strict_alpha_s_strict_source_shannon_source_binding_selection_or_preference_rule_schema_target_packet.json"
IN_F815 = GENERATED / "f815_current_strict_alpha_s_strict_source_shannon_source_binding_domain_member_grade_handling_clause_packet.json"

OUT = GENERATED / "f816_current_strict_alpha_s_selection_or_preference_rule_schema_continuation_boundary_packet.json"
OUT_SUMMARY = GENERATED / "f816_current_strict_alpha_s_selection_or_preference_rule_schema_continuation_boundary_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P816, IN_F814, IN_F815]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "F816",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p816 = load_json(IN_P816)
    f814 = load_json(IN_F814)
    f815 = load_json(IN_F815)

    p816_theorem = p816.get("theorem_result") or {}
    f814_target = f814.get("target_object") or {}
    f815_export = f815.get("exported_object") or {}

    if (
        p816.get("status")
        == "P816_CURRENT_STRICT_ALPHA_S_SELECTION_OR_PREFERENCE_RULE_SCHEMA_SAME_LANE_EXHAUSTION_BOUNDARY_AUDIT_PROBE"
        and p816_theorem.get("selection_or_preference_rule_schema_exported_on_current_repo_state") is False
        and p816_theorem.get("same_lane_passive_groundwork_exhausted") is True
        and p816_theorem.get("next_honest_move_requires_continuation_boundary_export") is True
        and f814.get("status")
        == "F814_EXECUTED_CURRENT_STRICT_ALPHA_S_STRICT_SOURCE_SHANNON_SOURCE_BINDING_SELECTION_OR_PREFERENCE_RULE_SCHEMA_TARGET_PACKET_NO_FALSE_PASS"
        and f815.get("status")
        == "F815_EXECUTED_CURRENT_STRICT_ALPHA_S_STRICT_SOURCE_SHANNON_SOURCE_BINDING_DOMAIN_MEMBER_GRADE_HANDLING_CLAUSE_PACKET_NO_FALSE_PASS"
    ):
        status = "F816_EXECUTED_CURRENT_STRICT_ALPHA_S_SELECTION_OR_PREFERENCE_RULE_SCHEMA_CONTINUATION_BOUNDARY_PACKET_NO_FALSE_PASS"
    else:
        status = "F816_REQUIRES_REVIEW"

    artifact = {
        "stage": "F816",
        "packet_name": "CurrentStrictAlphaSSelectionOrPreferenceRuleSchemaContinuationBoundary_v1",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p816_schema_same_lane_exhaustion_audit": rel(IN_P816),
            "f814_schema_only_target_packet": rel(IN_F814),
            "f815_grade_clause_packet": rel(IN_F815),
        },
        "exported_object": {
            "object_id": "alpha_s_selection_or_preference_rule_schema_continuation_boundary_v1",
            "goal": "Export the continuation boundary after the current same-lane passive groundwork under the missing selection_or_preference_rule_schema is exhausted.",
            "source_binding_selection_or_preference_rule_schema_target_ref": f814_target.get("object_id"),
            "candidate_source_support_domain_ref": f815_export.get("candidate_source_support_domain_ref"),
            "grade_clause_ref": f815_export.get("object_id"),
            "allowed_next_move_classes": [
                "export_one_genuinely_new_same_domain_selection_or_preference_schema_source",
                "shift_to_different_selection_provider_class_lane",
            ],
        },
        "current_honest_reading": [
            "The repo now exports an explicit continuation boundary for the source-binding selection/preference-schema lane.",
            "The exact schema-only target is still missing, but the current same-lane passive groundwork beneath it is exhausted.",
            "So the next honest move must be a genuinely new schema source or a different selection-provider class, not another fabricated passive split.",
        ],
        "recommended_next_move": "Attack one genuinely new same-domain selection_or_preference_schema source candidate, or justify a shift to a different selection-provider class lane without silent domain identification or grade promotion.",
        "hard_limits": [
            "Does not claim that the selection/preference rule schema already exists.",
            "Does not claim that any source has already been selected.",
            "Does not claim that the F811 source binding already exists.",
            "Does not claim that provider shift has already succeeded.",
            "Does not claim alpha_s boundary export readiness.",
            "Does not claim QCD closure.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F816",
        "status": status,
        "as_of": AS_OF,
        "exported_object_id": artifact["exported_object"]["object_id"],
        "recommended_next_move": artifact["recommended_next_move"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
