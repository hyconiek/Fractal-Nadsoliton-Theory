#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P815 = GENERATED / "p815_current_strict_alpha_s_domain_member_grade_handling_clause_export_admitted_selection_or_preference_rule_schema_still_blocked.json"
IN_F814 = GENERATED / "f814_current_strict_alpha_s_strict_source_shannon_source_binding_selection_or_preference_rule_schema_target_packet.json"
IN_F813 = GENERATED / "f813_current_strict_alpha_s_strict_source_shannon_source_binding_candidate_support_domain_packet.json"

OUT = GENERATED / "f815_current_strict_alpha_s_strict_source_shannon_source_binding_domain_member_grade_handling_clause_packet.json"
OUT_SUMMARY = GENERATED / "f815_current_strict_alpha_s_strict_source_shannon_source_binding_domain_member_grade_handling_clause_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P815, IN_F814, IN_F813]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "F815",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p815 = load_json(IN_P815)
    f814 = load_json(IN_F814)
    f813 = load_json(IN_F813)

    p815_theorem = p815.get("theorem_result") or {}
    f814_target = f814.get("target_object") or {}
    f813_export = f813.get("exported_object") or {}

    if (
        p815.get("status")
        == "P815_CURRENT_STRICT_ALPHA_S_DOMAIN_MEMBER_GRADE_HANDLING_CLAUSE_EXPORT_ADMITTED_SELECTION_OR_PREFERENCE_RULE_SCHEMA_STILL_BLOCKED"
        and p815_theorem.get("domain_member_grade_handling_clause_exportable_now") is True
        and p815_theorem.get("selection_or_preference_rule_schema_exported_now") is False
        and p815_theorem.get("next_honest_move_is_export_grade_clause_leave_selection_schema_blocked") is True
        and f814.get("status")
        == "F814_EXECUTED_CURRENT_STRICT_ALPHA_S_STRICT_SOURCE_SHANNON_SOURCE_BINDING_SELECTION_OR_PREFERENCE_RULE_SCHEMA_TARGET_PACKET_NO_FALSE_PASS"
        and f813.get("status")
        == "F813_EXECUTED_CURRENT_STRICT_ALPHA_S_STRICT_SOURCE_SHANNON_SOURCE_BINDING_CANDIDATE_SUPPORT_DOMAIN_PACKET_NO_FALSE_PASS"
    ):
        status = "F815_EXECUTED_CURRENT_STRICT_ALPHA_S_STRICT_SOURCE_SHANNON_SOURCE_BINDING_DOMAIN_MEMBER_GRADE_HANDLING_CLAUSE_PACKET_NO_FALSE_PASS"
    else:
        status = "F815_REQUIRES_REVIEW"

    grade_clause = (p815.get("grade_clause_candidate") or {}).get("fields") or {}

    artifact = {
        "stage": "F815",
        "packet_name": "CurrentStrictAlphaSStrictSourceShannonSourceBindingDomainMemberGradeHandlingClause_v1",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p815_grade_clause_probe": rel(IN_P815),
            "f814_schema_only_target_packet": rel(IN_F814),
            "f813_candidate_support_domain_packet": rel(IN_F813),
        },
        "exported_object": {
            "object_id": "alpha_s_strict_source_shannon_source_binding_domain_member_grade_handling_clause_v1",
            "goal": "Export the exact no-grade-promotion clause that governs the two members of the F813 support domain while leaving the source-binding selection/preference schema unresolved.",
            "source_binding_selection_or_preference_rule_schema_target_ref": f814_target.get("object_id"),
            "candidate_source_support_domain_ref": f813_export.get("object_id"),
            "clause_fields": grade_clause,
            "scope": "strict_source_side_domain_grade_discipline_only",
        },
        "current_honest_reading": [
            "The repo now exports one exact grade-handling clause for the two members of the F813 support domain.",
            "That clause keeps the F321 support at current-candidate-only grade and the T209/P755 support at future-only target grade below actual realization.",
            "It still does not choose one source and still does not export the selection/preference rule schema.",
        ],
        "recommended_next_move": "Build one narrow probe testing whether any current export can lawfully supply selection_or_preference_rule_schema for alpha_s_strict_source_shannon_source_binding_selection_or_preference_rule_schema_target_v1 on the already grade-disciplined F813 support domain without silent domain identification or grade promotion.",
        "hard_limits": [
            "Does not claim that any source has already been selected.",
            "Does not claim that the selection/preference rule schema already exists.",
            "Does not claim that the F811 source binding already exists.",
            "Does not claim that any adapter action schema is already exported.",
            "Does not claim that provider shift has already succeeded.",
            "Does not claim alpha_s boundary export readiness.",
            "Does not claim QCD closure.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F815",
        "status": status,
        "as_of": AS_OF,
        "exported_object_id": artifact["exported_object"]["object_id"],
        "candidate_source_support_domain_ref": artifact["exported_object"]["candidate_source_support_domain_ref"],
        "recommended_next_move": artifact["recommended_next_move"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
