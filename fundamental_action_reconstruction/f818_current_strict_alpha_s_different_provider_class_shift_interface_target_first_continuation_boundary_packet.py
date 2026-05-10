#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P818 = GENERATED / "p818_current_strict_alpha_s_different_provider_class_shift_interface_target_first_continuation_class_audit_probe.json"
IN_F817 = GENERATED / "f817_current_strict_alpha_s_different_selection_provider_class_shift_candidate_reference_packet.json"
IN_F816 = GENERATED / "f816_current_strict_alpha_s_selection_or_preference_rule_schema_continuation_boundary_packet.json"
IN_F807 = GENERATED / "f807_current_strict_alpha_s_provider_class_shift_requirement_packet.json"

OUT = GENERATED / "f818_current_strict_alpha_s_different_provider_class_shift_interface_target_first_continuation_boundary_packet.json"
OUT_SUMMARY = GENERATED / "f818_current_strict_alpha_s_different_provider_class_shift_interface_target_first_continuation_boundary_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P818, IN_F817, IN_F816, IN_F807]
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
    f816 = load_json(IN_F816)
    f807 = load_json(IN_F807)

    p818_theorem = p818.get("theorem_result") or {}
    f817_export = f817.get("exported_object") or {}
    f816_export = f816.get("exported_object") or {}
    f807_export = f807.get("exported_object") or {}

    if (
        p818.get("status")
        == "P818_CURRENT_STRICT_ALPHA_S_DIFFERENT_PROVIDER_CLASS_SHIFT_INTERFACE_TARGET_FIRST_CONTINUATION_CLASS_AUDITED"
        and p818_theorem.get("interface_target_first_continuation_class_admitted") is True
        and p818_theorem.get("immediate_next_move_class")
        == "freeze_exact_alpha_s_side_shift_interface_target"
        and f817.get("status")
        == "F817_EXECUTED_CURRENT_STRICT_ALPHA_S_DIFFERENT_SELECTION_PROVIDER_CLASS_SHIFT_CANDIDATE_REFERENCE_PACKET_NO_FALSE_PASS"
        and f816.get("status")
        == "F816_EXECUTED_CURRENT_STRICT_ALPHA_S_SELECTION_OR_PREFERENCE_RULE_SCHEMA_CONTINUATION_BOUNDARY_PACKET_NO_FALSE_PASS"
        and f807.get("status")
        == "F807_EXECUTED_CURRENT_STRICT_ALPHA_S_PROVIDER_CLASS_SHIFT_REQUIREMENT_PACKET_NO_FALSE_PASS"
    ):
        status = "F818_EXECUTED_CURRENT_STRICT_ALPHA_S_DIFFERENT_PROVIDER_CLASS_SHIFT_INTERFACE_TARGET_FIRST_CONTINUATION_BOUNDARY_PACKET_NO_FALSE_PASS"
    else:
        status = "F818_REQUIRES_REVIEW"

    artifact = {
        "stage": "F818",
        "packet_name": "CurrentStrictAlphaSDifferentProviderClassShiftInterfaceTargetFirstContinuationBoundary_v1",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p818_continuation_class_audit_probe": rel(IN_P818),
            "f817_different_provider_class_candidate_reference_packet": rel(IN_F817),
            "f816_schema_lane_continuation_boundary_packet": rel(IN_F816),
            "f807_provider_class_shift_requirement_packet": rel(IN_F807),
        },
        "exported_object": {
            "object_id": "alpha_s_different_provider_class_shift_interface_target_first_continuation_boundary_v1",
            "goal": "Export the continuation boundary that fixes interface-target-first ordering after F817 for the different provider-class alpha_s shift route.",
            "shift_candidate_reference_lane_ref": f817_export.get("object_id"),
            "upstream_continuation_boundary_ref": f816_export.get("object_id"),
            "provider_class_shift_requirement_ref": f807_export.get("object_id"),
            "immediate_admitted_next_move_class": p818_theorem.get("immediate_next_move_class"),
            "deferred_move_classes": [
                "freeze_carrier_identification_or_adapter_rule_target_after_interface_target",
                "freeze_adapter_action_schema_after_interface_target_if_still_missing",
            ],
            "forbidden_immediate_move_classes": [
                "freeze_carrier_identification_target_first",
                "freeze_adapter_rule_target_first",
                "verbal_promotion_of_t213_t216_candidate_lane_into_alpha_s_interface",
            ],
        },
        "current_honest_reading": [
            "The repo now exports the continuation boundary that fixes interface-target-first ordering for the different provider-class alpha_s shift route.",
            "This does not claim that the exact alpha_s-side shift interface already exists.",
            "It only says that interface-target work comes before adapter-rule or carrier-identification targeting on this route.",
        ],
        "recommended_next_move": "Run one narrow interface-target-oriented probe for the T213/T216 lane into alpha_s; do not jump directly to adapter-rule or carrier-identification targeting.",
        "hard_limits": [
            "Does not claim that the exact alpha_s-side shift interface already exists.",
            "Does not claim that any adapter rule already exists.",
            "Does not claim that carrier identification is already solved.",
            "Does not claim that the T213/T216 lane already enters the alpha_s domain.",
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
        "exported_object_id": artifact["exported_object"]["object_id"],
        "immediate_admitted_next_move_class": artifact["exported_object"][
            "immediate_admitted_next_move_class"
        ],
        "recommended_next_move": artifact["recommended_next_move"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
