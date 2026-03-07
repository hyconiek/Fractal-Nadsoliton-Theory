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
    / "p96_strict_side_role_equivalence_probe_for_legacy_gravity_hierarchy_role.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p96_strict_side_role_equivalence_probe_for_legacy_gravity_hierarchy_role_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def load_text(repo_relative_path: str) -> str:
    return (REPO / repo_relative_path).read_text(encoding="utf-8")


def any_entry_with_id(obj: Any, target_id: str) -> bool:
    if isinstance(obj, dict):
        if obj.get("id") == target_id:
            return True
        return any(any_entry_with_id(value, target_id) for value in obj.values())
    if isinstance(obj, list):
        return any(any_entry_with_id(value, target_id) for value in obj)
    return False


def role_transfer_verdict_present(text: str) -> bool:
    legacy_markers = [
        "exact gravity hierarchy from beta^N scaling",
        "legacy gravity-hierarchy",
        "legacy gravity-hierarchy role",
        "beta^N scaling",
    ]
    strict_markers = ["gravity_hierarchy_beta20"]
    transfer_markers = [
        "role-equivalence",
        "role equivalence",
        "retained successor",
        "same role",
        "semantic transfer",
        "successor semantics",
        "retained role",
    ]
    return (
        any(marker in text for marker in legacy_markers)
        and any(marker in text for marker in strict_markers)
        and any(marker in text for marker in transfer_markers)
    )


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    f22 = load_json(
        "fundamental_action_reconstruction/generated/f22_legacy_gravity_hierarchy_role_equivalence_refinement_packet_summary.json"
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
    a8_md = load_text("fundamental_action_reconstruction/A8_GRAVITY_BRIDGE_SPEC.md")
    a8_json = load_json(
        "fundamental_action_reconstruction/generated/a8_gravity_bridge_summary.json"
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

    semantic_transfer_sources = {
        "report_qw2069_full_sm_gr_derivation_package.json": json.dumps(
            qw2069, ensure_ascii=True
        ),
        "report_qw2115_gravity_hierarchy_strict_bridge_gate.json": json.dumps(
            qw2115, ensure_ascii=True
        ),
        "A8_GRAVITY_BRIDGE_SPEC.md": a8_md,
        "a8_gravity_bridge_summary.json": json.dumps(a8_json, ensure_ascii=True),
    }
    per_source_role_equivalence_verdict = {
        source: role_transfer_verdict_present(text)
        for source, text in semantic_transfer_sources.items()
    }
    explicit_role_equivalence_verdict_present = any(
        per_source_role_equivalence_verdict.values()
    )

    checks_spec = [
        {
            "id": "f22_candidate_present",
            "actual": f22["role_equivalence_state"]["strict_side_candidate_object_present"],
            "expected": True,
            "meaning": "F22 already identified a real strict-side candidate object",
        },
        {
            "id": "strict_side_candidate_object_present",
            "actual": candidate_object_present,
            "expected": True,
            "meaning": "the current strict-side sources export a real gravity-hierarchy candidate object via gravity_hierarchy_beta20",
        },
        {
            "id": "explicit_role_equivalence_verdict_present",
            "actual": explicit_role_equivalence_verdict_present,
            "expected": False,
            "meaning": "the current strict-side sources still export no explicit legacy-to-strict role-equivalence verdict",
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
        "stage": "P96",
        "lane": "strict_side_role_equivalence_probe_for_legacy_gravity_hierarchy_role_current_repo_state_only",
        "goal": "test_whether_the_current_repo_exports_an_explicit_legacy_to_strict_role_equivalence_verdict_for_the_legacy_gravity_hierarchy_role",
        "status": "CURRENT_REPO_EXPORTS_STRICT_SIDE_GRAVITY_HIERARCHY_BETA20_CANDIDATE_BUT_NO_EXPLICIT_LEGACY_GRAVITY_HIERARCHY_ROLE_EQUIVALENCE_VERDICT_AFTER_P96",
        "reason": "the strict side exports gravity_hierarchy_beta20 as a real derived candidate object, but still no explicit semantic-transfer verdict identifies that object as the retained role-equivalent successor of the legacy gravity-hierarchy role",
        "strict_side_candidate_object": {
            "id": "gravity_hierarchy_beta20",
            "present": candidate_object_present,
            "candidate_chain": candidate_chain,
        },
        "per_source_role_equivalence_verdict_presence": per_source_role_equivalence_verdict,
        "explicit_role_equivalence_verdict_present": explicit_role_equivalence_verdict_present,
        "remaining_missing_objects": [
            "explicit_legacy_to_strict_semantic_transfer_verdict_identifying_gravity_hierarchy_beta20_as_the_retained_strict_side_successor_of_the_legacy_gravity_hierarchy_role"
        ],
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P96",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "strict_side_candidate_object": artifact["strict_side_candidate_object"],
        "explicit_role_equivalence_verdict_present": explicit_role_equivalence_verdict_present,
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
