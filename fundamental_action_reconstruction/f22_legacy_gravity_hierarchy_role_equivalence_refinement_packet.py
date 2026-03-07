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
    / "f22_legacy_gravity_hierarchy_role_equivalence_refinement_packet.json"
)
OUT_SUMMARY = (
    GENERATED
    / "f22_legacy_gravity_hierarchy_role_equivalence_refinement_packet_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


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

    p95 = load_json(
        "fundamental_action_reconstruction/generated/p95_strict_side_literal_retention_probe_for_legacy_gravity_hierarchy_claim_summary.json"
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

    candidate_chain = {
        "q2068_registry_has_gravity_hierarchy_beta20": any_entry_with_id(
            qw2068, "gravity_hierarchy_beta20"
        ),
        "q2069_package_has_gravity_hierarchy_beta20": any_entry_with_id(
            qw2069, "gravity_hierarchy_beta20"
        ),
        "q2115_gate_has_gravity_hierarchy_beta20": any_entry_with_id(
            qw2115, "gravity_hierarchy_beta20"
        ),
    }
    candidate_object_present = all(candidate_chain.values())

    checks_spec = [
        {
            "id": "p95_literal_retention_path_closed",
            "actual": p95["literal_retention_present"],
            "expected": False,
            "meaning": "P95 already closes the literal-retention path negatively on the current repo state",
        },
        {
            "id": "strict_side_candidate_object_present",
            "actual": candidate_object_present,
            "expected": True,
            "meaning": "the current strict-side sources export a real gravity-hierarchy candidate object via gravity_hierarchy_beta20",
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
        "stage": "F22",
        "lane": "legacy_gravity_hierarchy_role_equivalence_refinement_current_repo_state_only",
        "goal": "refine_the_remaining_gravity_hierarchy_role_equivalence_frontier_into_candidate_object_presence_vs_explicit_semantic_transfer_verdict",
        "status": "F22_EXECUTED_LEGACY_GRAVITY_HIERARCHY_ROLE_EQUIVALENCE_REFINEMENT_PACKET_NO_FALSE_PASS",
        "reason": "the strict side exports gravity_hierarchy_beta20 as a real candidate object through QW-2068/QW-2069/QW-2115, but that candidate is still weaker than an explicit retained-side semantic-transfer verdict for the old gravity-hierarchy role",
        "role_equivalence_state": {
            "strict_side_candidate_object_present": candidate_object_present,
            "strict_side_candidate_object": "gravity_hierarchy_beta20",
            "candidate_chain": candidate_chain,
        },
        "remaining_missing_objects": [
            "explicit_legacy_to_strict_semantic_transfer_verdict_identifying_gravity_hierarchy_beta20_as_the_retained_strict_side_successor_of_the_legacy_gravity_hierarchy_role"
        ],
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "F22",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "role_equivalence_state": artifact["role_equivalence_state"],
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
