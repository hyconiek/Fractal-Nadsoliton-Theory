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
    / "f24_legacy_gravity_hierarchy_replaced_branch_refinement_packet.json"
)
OUT_SUMMARY = (
    GENERATED
    / "f24_legacy_gravity_hierarchy_replaced_branch_refinement_packet_summary.json"
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


def method_present_in_entries(payload: dict[str, Any], target_method: str) -> bool:
    return any(item.get("method") == target_method for item in payload.get("entries", []))


def update_has_method(payload: dict[str, Any], target_method: str) -> bool:
    update = payload.get("update")
    return isinstance(update, dict) and update.get("method") == target_method


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    n106 = load_json(
        "fundamental_action_reconstruction/generated/n106_current_legacy_gravity_hierarchy_retained_branch_full_negative_closure_theorem_summary.json"
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

    method_id = "qw2115_micro_supported_beta_hierarchy_bridge"

    strict_object_candidate_present = all(
        [
            registry_groups_have_id(qw2068, "gravity_hierarchy_beta20"),
            entries_have_id(qw2069, "gravity_hierarchy_beta20"),
            update_has_id(qw2115, "gravity_hierarchy_beta20"),
        ]
    )
    strict_method_candidate_present = all(
        [
            method_present_in_entries(qw2069, method_id),
            update_has_method(qw2115, method_id),
            method_id in json.dumps(qw2115, ensure_ascii=True),
        ]
    )

    checks_spec = [
        {
            "id": "n106_retained_branch_closed",
            "actual": n106["theorem_result"]["retained_branch_closed_negatively_on_current_repo_state"],
            "expected": True,
            "meaning": "N106 already closes the retained branch negatively on the current repo state",
        },
        {
            "id": "strict_object_candidate_present",
            "actual": strict_object_candidate_present,
            "expected": True,
            "meaning": "the strict side exports gravity_hierarchy_beta20 as a real successor-candidate object",
        },
        {
            "id": "strict_method_candidate_present",
            "actual": strict_method_candidate_present,
            "expected": True,
            "meaning": "the strict side exports qw2115_micro_supported_beta_hierarchy_bridge as a real successor-candidate method",
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
        "stage": "F24",
        "lane": "legacy_gravity_hierarchy_replaced_branch_refinement_current_repo_state_only",
        "goal": "refine_the_missing_replaced_branch_verdict_into_object_successor_vs_method_successor_subbranches",
        "status": "F24_EXECUTED_LEGACY_GRAVITY_HIERARCHY_REPLACED_BRANCH_REFINEMENT_PACKET_NO_FALSE_PASS",
        "reason": "N106 closes the retained branch negatively, while the repo already exports gravity_hierarchy_beta20 as a real successor-candidate object and qw2115_micro_supported_beta_hierarchy_bridge as a real successor-candidate method lineage; therefore the narrowest honest refinement is object-successor verdict vs method-successor-semantics verdict",
        "candidate_state": {
            "strict_object_candidate_present": strict_object_candidate_present,
            "strict_method_candidate_present": strict_method_candidate_present,
            "strict_object_candidate": "gravity_hierarchy_beta20",
            "strict_method_candidate": method_id,
        },
        "remaining_missing_objects": [
            "explicit_object_successor_verdict_identifying_gravity_hierarchy_beta20_as_the_strict_side_successor_object_replacing_the_legacy_gravity_hierarchy_role",
            "explicit_method_successor_semantics_verdict_identifying_qw2115_micro_supported_beta_hierarchy_bridge_as_the_strict_side_successor_semantics_replacing_the_legacy_gravity_hierarchy_role",
        ],
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "F24",
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
