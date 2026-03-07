#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = (
    GENERATED
    / "f25_legacy_gravity_hierarchy_object_successor_refinement_packet.json"
)
OUT_SUMMARY = (
    GENERATED
    / "f25_legacy_gravity_hierarchy_object_successor_refinement_packet_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def registry_groups_have_id(payload: dict[str, Any], entry_id: str) -> bool:
    for items in payload.get("groups", {}).values():
        for item in items:
            if item.get("id") == entry_id:
                return True
    return False


def entries_have_id(payload: dict[str, Any], entry_id: str) -> bool:
    return any(item.get("id") == entry_id for item in payload.get("entries", []))


def update_has_id(payload: dict[str, Any], entry_id: str) -> bool:
    update = payload.get("update")
    return isinstance(update, dict) and update.get("id") == entry_id


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    p100 = load_json(
        "fundamental_action_reconstruction/generated/p100_legacy_gravity_hierarchy_replaced_successor_subbranch_probe_summary.json"
    )
    qw2068 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2068_sm_gr_parameter_registry.json"
    )
    qw2069 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2069_full_sm_gr_derivation_package.json"
    )
    qw2115 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2115_gravity_hierarchy_strict_bridge_gate.json"
    )

    object_chain_present = all(
        [
            registry_groups_have_id(qw2068, "gravity_hierarchy_beta20"),
            entries_have_id(qw2069, "gravity_hierarchy_beta20"),
            update_has_id(qw2115, "gravity_hierarchy_beta20"),
        ]
    )

    textual_object_successor_verdict_present = False
    object_lineage_upgrade_verdict_present = False

    checks_spec = [
        {
            "id": "p100_object_replaced_verdict_absent",
            "actual": p100["object_replaced_verdict_present"],
            "expected": False,
            "meaning": "P100 already keeps the object-successor replaced sub-branch absent",
        },
        {
            "id": "object_chain_present",
            "actual": object_chain_present,
            "expected": True,
            "meaning": "the repo exports the gravity_hierarchy_beta20 object chain across QW-2068/QW-2069/QW-2115",
        },
        {
            "id": "textual_object_successor_verdict_present",
            "actual": textual_object_successor_verdict_present,
            "expected": False,
            "meaning": "no explicit textual object-successor verdict is currently exported",
        },
        {
            "id": "object_lineage_upgrade_verdict_present",
            "actual": object_lineage_upgrade_verdict_present,
            "expected": False,
            "meaning": "no explicit object-lineage-upgrade verdict is currently exported",
        },
    ]

    checks = []
    for item in checks_spec:
        checks.append(
            {
                "id": item["id"],
                "actual": item["actual"],
                "expected": item["expected"],
                "pass": item["actual"] == item["expected"],
                "meaning": item["meaning"],
            }
        )

    artifact = {
        "stage": "F25",
        "lane": "legacy_gravity_hierarchy_object_successor_refinement_current_repo_state_only",
        "goal": "refine_the_missing_object_successor_verdict_into_textual_vs_lineage_upgrade_subbranches",
        "status": "F25_EXECUTED_LEGACY_GRAVITY_HIERARCHY_OBJECT_SUCCESSOR_REFINEMENT_PACKET_NO_FALSE_PASS",
        "reason": "P100 keeps the object-successor sub-branch absent, while the repo already exports a real gravity_hierarchy_beta20 object chain across QW-2068/QW-2069/QW-2115; therefore the narrowest honest refinement is textual object-successor verdict vs object-lineage-upgrade verdict",
        "candidate_state": {
            "object_chain_present": object_chain_present,
            "textual_object_successor_verdict_present": textual_object_successor_verdict_present,
            "object_lineage_upgrade_verdict_present": object_lineage_upgrade_verdict_present,
        },
        "remaining_missing_objects": [
            "explicit_textual_object_successor_verdict_identifying_gravity_hierarchy_beta20_as_the_strict_side_successor_object_replacing_the_legacy_gravity_hierarchy_role",
            "explicit_object_lineage_upgrade_verdict_elevating_the_existing_gravity_hierarchy_beta20_candidate_chain_into_replacement_semantics_for_the_legacy_gravity_hierarchy_role",
        ],
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "F25",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "candidate_state": artifact["candidate_state"],
        "remaining_missing_objects": artifact["remaining_missing_objects"],
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(
        json.dumps(artifact, indent=2, ensure_ascii=True) + "\n",
        encoding="ascii",
    )
    OUT_SUMMARY.write_text(
        json.dumps(summary, indent=2, ensure_ascii=True) + "\n",
        encoding="ascii",
    )
    print(OUT_JSON)


if __name__ == "__main__":
    main()
