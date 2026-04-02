#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P813 = GENERATED / "p813_current_strict_alpha_s_source_binding_candidate_support_domain_export_admitted_selection_or_preference_rule_schema_still_blocked.json"
IN_F812 = GENERATED / "f812_current_strict_alpha_s_strict_source_shannon_shift_to_domain_source_binding_selection_or_preference_rule_target_packet.json"

OUT = GENERATED / "f813_current_strict_alpha_s_strict_source_shannon_source_binding_candidate_support_domain_packet.json"
OUT_SUMMARY = GENERATED / "f813_current_strict_alpha_s_strict_source_shannon_source_binding_candidate_support_domain_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P813, IN_F812]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "F813",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p813 = load_json(IN_P813)
    f812 = load_json(IN_F812)

    if (
        p813.get("status")
        == "P813_CURRENT_STRICT_ALPHA_S_SOURCE_BINDING_CANDIDATE_SUPPORT_DOMAIN_EXPORT_ADMITTED_SELECTION_OR_PREFERENCE_RULE_SCHEMA_STILL_BLOCKED"
        and (p813.get("theorem_result") or {}).get("candidate_source_support_domain_exportable_now")
        is True
        and (p813.get("theorem_result") or {}).get("selection_or_preference_rule_schema_exported_now")
        is False
        and (p813.get("theorem_result") or {}).get("next_honest_move_is_export_support_domain_leave_rule_schema_blocked")
        is True
        and f812.get("status")
        == "F812_EXECUTED_CURRENT_STRICT_ALPHA_S_STRICT_SOURCE_SHANNON_SHIFT_TO_DOMAIN_SOURCE_BINDING_SELECTION_OR_PREFERENCE_RULE_TARGET_PACKET_NO_FALSE_PASS"
    ):
        status = "F813_EXECUTED_CURRENT_STRICT_ALPHA_S_STRICT_SOURCE_SHANNON_SOURCE_BINDING_CANDIDATE_SUPPORT_DOMAIN_PACKET_NO_FALSE_PASS"
    else:
        status = "F813_REQUIRES_REVIEW"

    support_domain = p813.get("support_domain_candidate") or {}
    f812_target = f812.get("target_object") or {}
    f812_refs = f812.get("target_refs") or {}

    artifact = {
        "stage": "F813",
        "packet_name": "CurrentStrictAlphaSStrictSourceShannonSourceBindingCandidateSupportDomain_v1",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p813_support_domain_probe": rel(IN_P813),
            "f812_selection_or_preference_rule_target_packet": rel(IN_F812),
        },
        "exported_object": {
            "object_id": "alpha_s_strict_source_shannon_source_binding_candidate_support_domain_v1",
            "goal": "Export the exact finite support domain on which a future source-binding selection/preference rule would act, without selecting one source.",
            "source_binding_selection_or_preference_rule_target_ref": f812_target.get("object_id"),
            "support_domain_kind": "finite_two_support_domain",
            "member_count": support_domain.get("member_count"),
            "members": support_domain.get("members"),
            "unresolved_selection_or_preference_rule_schema_ref": "alpha_s_strict_source_shannon_shift_to_domain_source_binding_selection_or_preference_rule_target_v1",
            "scope": "strict_source_side_candidate_support_domain_only",
        },
        "current_honest_reading": [
            "The repo now exports one explicit finite support domain for source binding on the strict-source Shannon -> alpha_s route.",
            "That domain contains the current candidate-only support and the lawful future-only entry support, with their grades kept explicit.",
            "It still does not choose one source and still does not export the selection/preference rule schema.",
        ],
        "recommended_next_move": "Build one narrow probe testing whether any current export can lawfully supply selection_or_preference_rule_schema for alpha_s_strict_source_shannon_shift_to_domain_source_binding_selection_or_preference_rule_target_v1 on the exported candidate support domain without silent domain identification.",
        "hard_limits": [
            "Does not claim that any source has already been selected.",
            "Does not claim that the selection/preference rule already exists.",
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
        "stage": "F813",
        "status": status,
        "as_of": AS_OF,
        "exported_object_id": artifact["exported_object"]["object_id"],
        "source_binding_selection_or_preference_rule_target_ref": artifact["exported_object"][
            "source_binding_selection_or_preference_rule_target_ref"
        ],
        "recommended_next_move": artifact["recommended_next_move"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
