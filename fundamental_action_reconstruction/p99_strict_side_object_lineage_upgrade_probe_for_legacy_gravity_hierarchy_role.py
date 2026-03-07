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
    / "p99_strict_side_object_lineage_upgrade_probe_for_legacy_gravity_hierarchy_role.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p99_strict_side_object_lineage_upgrade_probe_for_legacy_gravity_hierarchy_role_summary.json"
)

STRICT_SIDE_OBJECT_SOURCES = [
    "RELEASE_4_9_TEXTBOOK_EN_PL.md",
    "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2068_sm_gr_parameter_registry.json",
    "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2069_full_sm_gr_derivation_package.json",
    "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2115_gravity_hierarchy_strict_bridge_gate.json",
    "fundamental_action_reconstruction/A8_GRAVITY_BRIDGE_SPEC.md",
    "fundamental_action_reconstruction/generated/a8_gravity_bridge_summary.json",
]

LEGACY_MARKERS = [
    "legacy gravity-hierarchy",
    "legacy gravity-hierarchy role",
    "exact gravity hierarchy from beta^n scaling",
    "beta^n scaling",
]
OBJECT_MARKERS = ["gravity_hierarchy_beta20"]
LINEAGE_MARKERS = [
    "strict-derived",
    "strict derived",
    "target parameter",
    "gravity_hierarchy_beta20",
]
UPGRADE_MARKERS = [
    "upgraded to retained",
    "semantic transfer",
    "retained successor",
    "role-equivalence",
    "role equivalence",
    "legacy-to-strict",
]


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def load_text(repo_relative_path: str) -> str:
    return (REPO / repo_relative_path).read_text(encoding="utf-8")


def object_lineage_upgrade_verdict_present(text: str) -> bool:
    lowered = text.lower()
    return (
        any(marker in lowered for marker in LEGACY_MARKERS)
        and any(marker in lowered for marker in OBJECT_MARKERS)
        and any(marker in lowered for marker in LINEAGE_MARKERS)
        and any(marker in lowered for marker in UPGRADE_MARKERS)
    )


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    p98 = load_json(
        "fundamental_action_reconstruction/generated/p98_strict_side_textual_retained_successor_probe_for_legacy_gravity_hierarchy_role_summary.json"
    )

    per_source = {}
    any_lineage_upgrade = False
    for source in STRICT_SIDE_OBJECT_SOURCES:
        present = object_lineage_upgrade_verdict_present(load_text(source))
        per_source[source] = present
        any_lineage_upgrade = any_lineage_upgrade or present

    checks_spec = [
        {
            "id": "p98_textual_successor_subbranch_closed",
            "actual": p98["textual_successor_verdict_present"],
            "expected": False,
            "meaning": "P98 already closes the textual retained-successor sub-branch negatively on the current repo state",
        },
        {
            "id": "strict_side_object_lineage_upgrade_verdict_present",
            "actual": any_lineage_upgrade,
            "expected": False,
            "meaning": "the current strict-side object source set exports no explicit object-lineage-upgrade verdict for the legacy gravity-hierarchy role",
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
        "stage": "P99",
        "lane": "strict_side_object_lineage_upgrade_probe_for_legacy_gravity_hierarchy_role_current_repo_state_only",
        "goal": "test_whether_the_current_strict_side_object_source_set_exports_an_explicit_object_lineage_upgrade_verdict_for_the_legacy_gravity_hierarchy_role",
        "status": "CURRENT_STRICT_SIDE_SOURCES_DO_NOT_EXPORT_AN_EXPLICIT_OBJECT_LINEAGE_UPGRADE_VERDICT_FOR_THE_LEGACY_GRAVITY_HIERARCHY_ROLE_AFTER_P99",
        "reason": "the strict-side object sources do export the gravity_hierarchy_beta20 chain, but none of those sources upgrades that chain into explicit retained semantics for the legacy gravity-hierarchy role",
        "strict_side_sources_checked": STRICT_SIDE_OBJECT_SOURCES,
        "per_source_object_lineage_upgrade_verdict_presence": per_source,
        "object_lineage_upgrade_verdict_present": any_lineage_upgrade,
        "retained_branch_closed_on_current_repo_state": True,
        "remaining_missing_objects": [
            "explicit_strict_side_replaced_verdict_for_the_legacy_gravity_hierarchy_role_by_an_explicit_strict_successor_semantics"
        ],
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P99",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "per_source_object_lineage_upgrade_verdict_presence": per_source,
        "object_lineage_upgrade_verdict_present": any_lineage_upgrade,
        "retained_branch_closed_on_current_repo_state": True,
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
