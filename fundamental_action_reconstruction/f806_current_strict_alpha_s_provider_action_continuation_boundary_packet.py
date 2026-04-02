#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P806 = GENERATED / "p806_current_strict_alpha_s_provider_action_same_lane_exhaustion_boundary_audit_probe.json"
IN_F802 = GENERATED / "f802_current_strict_alpha_s_provider_action_rule_target_packet.json"
IN_F805 = GENERATED / "f805_current_strict_alpha_s_acting_input_bundle_packet.json"

OUT = GENERATED / "f806_current_strict_alpha_s_provider_action_continuation_boundary_packet.json"
OUT_SUMMARY = GENERATED / "f806_current_strict_alpha_s_provider_action_continuation_boundary_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P806, IN_F802, IN_F805]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "F806",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p806 = load_json(IN_P806)
    f802 = load_json(IN_F802)
    f805 = load_json(IN_F805)

    if (
        p806.get("status")
        == "P806_CURRENT_STRICT_ALPHA_S_PROVIDER_ACTION_SAME_LANE_EXHAUSTION_BOUNDARY_AUDITED"
        and (p806.get("theorem_result") or {}).get("boundary_exported_on_current_repo_state") is True
        and (p806.get("theorem_result") or {}).get(
            "same_level_alpha_s_provider_action_lane_continuation_no_longer_admitted_primary_move"
        )
        is True
        and (p806.get("theorem_result") or {}).get(
            "next_honest_move_requires_genuinely_new_provider_action_source_or_provider_shift"
        )
        is True
        and f802.get("status")
        == "F802_EXECUTED_CURRENT_STRICT_ALPHA_S_PROVIDER_ACTION_RULE_TARGET_PACKET_NO_FALSE_PASS"
        and f805.get("status")
        == "F805_EXECUTED_CURRENT_STRICT_ALPHA_S_ACTING_INPUT_BUNDLE_PACKET_NO_FALSE_PASS"
    ):
        status = "F806_EXECUTED_CURRENT_STRICT_ALPHA_S_PROVIDER_ACTION_CONTINUATION_BOUNDARY_PACKET_NO_FALSE_PASS"
    else:
        status = "F806_REQUIRES_REVIEW"

    theorem_result = p806.get("theorem_result") or {}
    acting_input = f805.get("exported_object") or {}
    action_target = f802.get("target_object") or {}

    artifact = {
        "stage": "F806",
        "packet_name": "CurrentStrictAlphaSProviderActionContinuationBoundary_v1",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p806_same_lane_exhaustion_boundary_probe": rel(IN_P806),
            "f802_provider_action_rule_target_packet": rel(IN_F802),
            "f805_acting_input_bundle_packet": rel(IN_F805),
        },
        "exported_object": {
            "object_id": "alpha_s_provider_action_continuation_boundary_v1",
            "goal": "Export the current continuation boundary after the alpha_s same-domain passive lane is exhausted but the provider action rule remains missing.",
            "same_lane_exhaustion_boundary_ref": theorem_result.get("boundary_name"),
            "current_acting_input_bundle_ref": acting_input.get("object_id"),
            "current_provider_action_target_ref": action_target.get("object_id"),
            "admitted_next_move_classes": [
                "export_one_genuinely_new_same_domain_provider_action_source",
                "shift_to_a_different_provider_class_lane",
            ],
            "forbidden_repetition_clause": "no_further_same_level_passive_split_under_unchanged_provider_action_rule_ref_as_primary_move",
            "scope": "strict_alpha_s_provider_action_continuation_boundary_only",
        },
        "current_honest_reading": [
            "The repo now exports an explicit continuation boundary for the alpha_s provider-action lane.",
            "The current passive same-domain chain has been exported as far as the present honest local split goes.",
            "The next honest move is now constrained to a genuinely new provider-action source or a provider-class shift.",
        ],
        "recommended_next_move": "Attack one genuinely new same-domain provider-action source candidate before attempting any further alpha_s provider-action promotion on the current lane.",
        "hard_limits": [
            "Does not claim that the provider action rule already exists.",
            "Does not claim that a genuinely new provider-action source already exists.",
            "Does not claim that provider-class shift has already succeeded.",
            "Does not claim alpha_s boundary export readiness.",
            "Does not claim QCD closure.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F806",
        "status": status,
        "as_of": AS_OF,
        "exported_object_id": artifact["exported_object"]["object_id"],
        "current_provider_action_target_ref": artifact["exported_object"]["current_provider_action_target_ref"],
        "recommended_next_move": artifact["recommended_next_move"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
