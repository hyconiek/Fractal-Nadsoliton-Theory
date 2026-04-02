#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P802 = GENERATED / "p802_current_strict_alpha_s_provider_class_reorganization_audit_probe.json"
IN_F801 = GENERATED / "f801_current_strict_alpha_s_same_domain_provider_skeleton_packet.json"

OUT = GENERATED / "f802_current_strict_alpha_s_provider_action_rule_target_packet.json"
OUT_SUMMARY = GENERATED / "f802_current_strict_alpha_s_provider_action_rule_target_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P802, IN_F801]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "F802",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p802 = load_json(IN_P802)
    f801 = load_json(IN_F801)

    if (
        p802.get("status")
        == "P802_CURRENT_STRICT_ALPHA_S_PASSIVE_PROVIDER_SKELETON_SUPPORTED_ACTIVE_PROVIDER_ACTION_RULE_BLOCKED"
        and (p802.get("clause_split") or {}).get("sharp_blocker_clause") == "provider_action_rule_ref"
        and f801.get("status")
        == "F801_EXECUTED_CURRENT_STRICT_ALPHA_S_SAME_DOMAIN_PROVIDER_SKELETON_PACKET_NO_FALSE_PASS"
    ):
        status = "F802_EXECUTED_CURRENT_STRICT_ALPHA_S_PROVIDER_ACTION_RULE_TARGET_PACKET_NO_FALSE_PASS"
    else:
        status = "F802_REQUIRES_REVIEW"

    skeleton = f801.get("exported_object") or {}

    artifact = {
        "stage": "F802",
        "packet_name": "CurrentStrictAlphaSProviderActionRuleTarget_v1",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p802_provider_class_reorganization_probe": rel(IN_P802),
            "f801_provider_skeleton_packet": rel(IN_F801),
        },
        "why_this_packet_exists": [
            "F801 exports the passive same-domain provider skeleton for the alpha_s reference-scale lane.",
            "P802 shows that no active provider action rule is currently exported on that same domain.",
        ],
        "target_object": {
            "object_id": "alpha_s_reference_scale_provider_action_rule_target_v1",
            "goal": "Freeze the exact active-rule object still missing before the passive same-domain skeleton can count as a real provider class.",
            "required_fields": [
                {
                    "name": "provider_skeleton_ref",
                    "required": True,
                    "hard_limit": "Must point to one explicit same-domain provider skeleton.",
                },
                {
                    "name": "support_bundle_ref",
                    "required": True,
                    "hard_limit": "Must point to the already-exported support bundle carried by that skeleton.",
                },
                {
                    "name": "acting_same_domain_input_ref",
                    "required": True,
                    "hard_limit": "Must identify the same-domain input on which the active rule acts; foreign-domain inputs are forbidden.",
                },
                {
                    "name": "provider_action_rule_ref",
                    "required": True,
                    "hard_limit": "Must export the active rule that makes the passive skeleton act as a provider class.",
                },
                {
                    "name": "semantic_principle_supply_ref",
                    "required": True,
                    "hard_limit": "Must state how the active rule supplies the missing semantic principle on the same domain.",
                },
                {
                    "name": "passive_to_active_upgrade_block_ref",
                    "required": True,
                    "hard_limit": "Must explicitly record why passive support does not auto-upgrade into an active provider class.",
                },
                {
                    "name": "foreign_domain_exclusion_ref",
                    "required": True,
                    "hard_limit": "Must explicitly fence off foreign-domain reference structures and host-semantics imports.",
                },
                {
                    "name": "selected_provider_class_output_schema",
                    "required": True,
                    "hard_limit": "Must output the admitted provider-class action and its scope on the alpha_s lane.",
                },
                {
                    "name": "hard_limits",
                    "required": True,
                    "hard_limit": "Must explicitly deny Standard Model identification, QCD closure, and ToE closure.",
                },
            ],
        },
        "provider_skeleton_ref": skeleton.get("object_id"),
        "current_honest_reading": [
            "The current sharp blocker is no longer whether the alpha_s lane has a passive same-domain provider base.",
            "It is the missing active rule that would make that passive base act as a provider class.",
            "F802 freezes that exact missing object without promoting the current skeleton into alpha_s closure.",
        ],
        "recommended_next_move": "Build one narrow probe testing whether any current same-domain relation inside the passive provider skeleton can be strengthened into alpha_s_reference_scale_provider_action_rule_target_v1 without importing foreign-domain reference structures or new host semantics.",
        "hard_limits": [
            "Does not claim that the active rule already exists.",
            "Does not claim that the provider class already exists.",
            "Does not claim that the semantic principle already exists.",
            "Does not claim that the F704 maximum is already a lawful reference-scale point.",
            "Does not claim alpha_s boundary export readiness.",
            "Does not claim QCD closure.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F802",
        "status": status,
        "as_of": AS_OF,
        "target_object_id": artifact["target_object"]["object_id"],
        "provider_skeleton_ref": artifact["provider_skeleton_ref"],
        "recommended_next_move": artifact["recommended_next_move"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
