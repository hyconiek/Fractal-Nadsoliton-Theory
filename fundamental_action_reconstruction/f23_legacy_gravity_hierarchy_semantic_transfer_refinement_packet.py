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
    / "f23_legacy_gravity_hierarchy_semantic_transfer_refinement_packet.json"
)
OUT_SUMMARY = (
    GENERATED
    / "f23_legacy_gravity_hierarchy_semantic_transfer_refinement_packet_summary.json"
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


def any_entry_with_id(obj: Any, target_id: str) -> bool:
    if isinstance(obj, dict):
        if obj.get("id") == target_id:
            return True
        return any(any_entry_with_id(value, target_id) for value in obj.values())
    if isinstance(obj, list):
        return any(any_entry_with_id(value, target_id) for value in obj)
    return False


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    p96 = load_json(
        "fundamental_action_reconstruction/generated/p96_strict_side_role_equivalence_probe_for_legacy_gravity_hierarchy_role_summary.json"
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
            any_entry_with_id(qw2115, "gravity_hierarchy_beta20"),
        ]
    )

    checks_spec = [
        {
            "id": "p96_candidate_present",
            "actual": p96["strict_side_candidate_object"]["present"],
            "expected": True,
            "meaning": "P96 already confirms that a real strict-side gravity-hierarchy candidate object is exported",
        },
        {
            "id": "object_chain_present",
            "actual": object_chain_present,
            "expected": True,
            "meaning": "the repo exports the gravity_hierarchy_beta20 candidate chain across QW-2068/QW-2069/QW-2115",
        },
    ]

    checks: list[dict[str, Any]] = []
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
        "stage": "F23",
        "lane": "legacy_gravity_hierarchy_semantic_transfer_refinement_current_repo_state_only",
        "goal": "refine_the_missing_retained_semantic_transfer_verdict_into_textual_vs_object_lineage_upgrade_subbranches",
        "status": "F23_EXECUTED_LEGACY_GRAVITY_HIERARCHY_SEMANTIC_TRANSFER_REFINEMENT_PACKET_NO_FALSE_PASS",
        "reason": "P96 keeps the retained semantic-transfer verdict absent, while the repo already exports a real gravity_hierarchy_beta20 candidate object and a real strict-side candidate chain; therefore the narrowest honest refinement is textual retained-successor verdict vs object-lineage-upgrade verdict",
        "candidate_state": {
            "strict_side_candidate_object_present": p96["strict_side_candidate_object"]["present"],
            "strict_side_candidate_object": "gravity_hierarchy_beta20",
            "object_chain_present": object_chain_present,
        },
        "remaining_missing_objects": [
            "explicit_textual_retained_successor_verdict_identifying_gravity_hierarchy_beta20_as_the_retained_strict_side_successor_of_the_legacy_gravity_hierarchy_role",
            "explicit_object_lineage_upgrade_verdict_elevating_the_existing_gravity_hierarchy_beta20_candidate_chain_into_retained_strict_side_gravity_hierarchy_role_transfer",
        ],
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "F23",
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
