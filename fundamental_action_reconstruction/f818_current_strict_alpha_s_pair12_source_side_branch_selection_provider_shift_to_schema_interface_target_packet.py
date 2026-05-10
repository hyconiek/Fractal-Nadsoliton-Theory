#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P818 = GENERATED / "p818_current_strict_alpha_s_no_exact_pair12_source_side_branch_selection_provider_shift_to_schema_interface_target_exported_target_freeze_required.json"
IN_F817 = GENERATED / "f817_current_strict_alpha_s_different_selection_provider_class_shift_candidate_reference_packet.json"
IN_F814 = GENERATED / "f814_current_strict_alpha_s_strict_source_shannon_source_binding_selection_or_preference_rule_schema_target_packet.json"
IN_F815 = GENERATED / "f815_current_strict_alpha_s_strict_source_shannon_source_binding_domain_member_grade_handling_clause_packet.json"
IN_P764 = GENERATED / "p764_current_strict_t218_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_target_probe_summary.json"

OUT = GENERATED / "f818_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_schema_interface_target_packet.json"
OUT_SUMMARY = GENERATED / "f818_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_schema_interface_target_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P818, IN_F817, IN_F814, IN_F815, IN_P764]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "F818",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p818 = load_json(IN_P818)
    f817 = load_json(IN_F817)
    f814 = load_json(IN_F814)
    f815 = load_json(IN_F815)
    p764 = load_json(IN_P764)

    p818_theorem = p818.get("theorem_result") or {}
    missing_target = p818.get("exact_missing_interface_target_candidate") or {}
    f817_export = f817.get("exported_object") or {}
    f814_target = f814.get("target_object") or {}
    f815_export = f815.get("exported_object") or {}

    if (
        p818.get("status")
        == "P818_CURRENT_STRICT_ALPHA_S_NO_EXACT_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_SHIFT_TO_SCHEMA_INTERFACE_TARGET_EXPORTED_TARGET_FREEZE_REQUIRED"
        and p818_theorem.get("exact_pair12_source_side_branch_selection_provider_shift_to_schema_interface_target_exported")
        is False
        and p818_theorem.get("next_honest_move_requires_freezing_exact_shift_to_schema_interface_target")
        is True
        and f817.get("status")
        == "F817_EXECUTED_CURRENT_STRICT_ALPHA_S_DIFFERENT_SELECTION_PROVIDER_CLASS_SHIFT_CANDIDATE_REFERENCE_PACKET_NO_FALSE_PASS"
        and f814.get("status")
        == "F814_EXECUTED_CURRENT_STRICT_ALPHA_S_STRICT_SOURCE_SHANNON_SOURCE_BINDING_SELECTION_OR_PREFERENCE_RULE_SCHEMA_TARGET_PACKET_NO_FALSE_PASS"
        and f815.get("status")
        == "F815_EXECUTED_CURRENT_STRICT_ALPHA_S_STRICT_SOURCE_SHANNON_SOURCE_BINDING_DOMAIN_MEMBER_GRADE_HANDLING_CLAUSE_PACKET_NO_FALSE_PASS"
    ):
        status = "F818_EXECUTED_CURRENT_STRICT_ALPHA_S_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_SHIFT_TO_SCHEMA_INTERFACE_TARGET_PACKET_NO_FALSE_PASS"
    else:
        status = "F818_REQUIRES_REVIEW"

    artifact = {
        "stage": "F818",
        "packet_name": "CurrentStrictAlphaSPair12SourceSideBranchSelectionProviderShiftToSchemaInterfaceTarget_v1",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p818_shift_to_schema_interface_audit_probe": rel(IN_P818),
            "f817_different_provider_class_shift_candidate_reference_packet": rel(IN_F817),
            "f814_downstream_schema_target_packet": rel(IN_F814),
            "f815_downstream_grade_clause_packet": rel(IN_F815),
            "p764_own_lane_missing_interface_target_summary": rel(IN_P764),
        },
        "target_object": {
            "object_id": "alpha_s_pair12_source_side_branch_selection_provider_shift_to_schema_interface_target_v1",
            "goal": "Freeze the exact missing interface target from the T213/T216 provider-class shift candidate lane into the current alpha_s schema problem without claiming interface realization.",
            "required_fields": [
                {
                    "name": "different_selection_provider_class_shift_candidate_reference_lane_ref",
                    "required": True,
                    "hard_limit": "Must point to the exact F817 candidate-reference lane and not silently replace it."
                },
                {
                    "name": "downstream_selection_or_preference_rule_schema_target_ref",
                    "required": True,
                    "hard_limit": "Must point to the exact F814 downstream schema-only target and not silently replace the problem."
                },
                {
                    "name": "downstream_domain_member_grade_handling_clause_ref",
                    "required": True,
                    "hard_limit": "Must preserve the exact F815 grade discipline on the downstream domain; no silent grade promotion is allowed."
                },
                {
                    "name": "self_lane_missing_interface_target_ref",
                    "required": True,
                    "hard_limit": "Must keep explicit that the candidate lane still has its own unresolved T218 interface target and cannot silently bypass it."
                },
                {
                    "name": "shift_interface_adapter_or_carrier_identification_rule_ref",
                    "required": True,
                    "hard_limit": "Must export the exact adapter or carrier-safe identification rule; silent domain identification is forbidden."
                },
                {
                    "name": "nontransfer_boundary_ref",
                    "required": True,
                    "hard_limit": "Must block silent reuse of unrelated alpha_s or foreign-domain interface artifacts."
                },
                {
                    "name": "future_route_grade_ref",
                    "required": True,
                    "hard_limit": "Must keep this interface target at future-route target grade until a real interface is exported."
                },
                {
                    "name": "exact_interface_output_schema",
                    "required": True,
                    "hard_limit": "Must state what exact lawful output would enter the current alpha_s schema problem if the interface were exported."
                },
                {
                    "name": "hard_limits",
                    "required": True,
                    "hard_limit": "Must explicitly deny interface realization, provider-class shift success, QCD closure, and ToE closure."
                },
            ],
        },
        "target_refs": {
            "different_selection_provider_class_shift_candidate_reference_lane_ref": f817_export.get("object_id"),
            "downstream_selection_or_preference_rule_schema_target_ref": f814_target.get("object_id"),
            "downstream_domain_member_grade_handling_clause_ref": f815_export.get("object_id"),
            "self_lane_missing_interface_target_ref": p764.get("t218_target_name"),
            "missing_target_candidate_id": missing_target.get("candidate_id"),
        },
        "current_honest_reading": [
            "The repo now exports the exact missing interface target between the admitted T213/T216 candidate lane and the current alpha_s schema problem.",
            "This sits above the admitted different-provider-class candidate lane and below any adapter/rule that would instantiate it.",
            "It does not claim that the interface exists; it only localizes the missing object sharply.",
        ],
        "recommended_next_move": "Build one narrow probe testing whether any current export can lawfully supply shift_interface_adapter_or_carrier_identification_rule_ref for alpha_s_pair12_source_side_branch_selection_provider_shift_to_schema_interface_target_v1 without silent domain identification.",
        "hard_limits": [
            "Does not claim that the alpha_s selection/preference schema already exists.",
            "Does not claim that the T213/T216 lane already enters the alpha_s schema domain.",
            "Does not claim that any adapter or carrier-identification rule is already exported.",
            "Does not claim that provider-class shift has already succeeded.",
            "Does not claim alpha_s boundary export readiness.",
            "Does not claim QCD closure.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F818",
        "status": status,
        "as_of": AS_OF,
        "target_object_id": artifact["target_object"]["object_id"],
        "different_selection_provider_class_shift_candidate_reference_lane_ref": artifact["target_refs"][
            "different_selection_provider_class_shift_candidate_reference_lane_ref"
        ],
        "recommended_next_move": artifact["recommended_next_move"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
