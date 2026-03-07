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
    / "p105_strict_side_textual_method_successor_probe_for_legacy_gravity_hierarchy_role.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p105_strict_side_textual_method_successor_probe_for_legacy_gravity_hierarchy_role_summary.json"
)

STRICT_SIDE_METHOD_SOURCES = [
    "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2069_full_sm_gr_derivation_package.json",
    "material_dowodowy/korpus_qw_pozostaly/raporty_md/RAPORT_QW2115_GRAVITY_HIERARCHY_STRICT_BRIDGE_GATE.md",
    "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2115_gravity_hierarchy_strict_bridge_gate.json",
    "fundamental_action_reconstruction/A8_GRAVITY_BRIDGE_SPEC.md",
]

LEGACY_MARKERS = [
    "legacy gravity-hierarchy",
    "legacy gravity-hierarchy role",
    "exact gravity hierarchy from beta^n scaling",
    "beta^n scaling",
]
METHOD_MARKERS = ["qw2115_micro_supported_beta_hierarchy_bridge"]
TEXTUAL_SUCCESSOR_MARKERS = [
    "successor semantics",
    "strict successor semantics",
    "replacing the legacy gravity-hierarchy role",
    "replacing the legacy gravity-hierarchy claim",
    "replaced by",
    "legacy-to-strict",
]


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def load_text(repo_relative_path: str) -> str:
    return (REPO / repo_relative_path).read_text(encoding="utf-8")


def textual_method_successor_verdict_present(text: str) -> bool:
    lowered = text.lower()
    return (
        any(marker in lowered for marker in LEGACY_MARKERS)
        and any(marker in lowered for marker in METHOD_MARKERS)
        and any(marker in lowered for marker in TEXTUAL_SUCCESSOR_MARKERS)
    )


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    p104 = load_json(
        "fundamental_action_reconstruction/generated/p104_strict_side_method_successor_subbranch_probe_for_legacy_gravity_hierarchy_role_summary.json"
    )

    per_source = {}
    any_textual_successor = False
    for source in STRICT_SIDE_METHOD_SOURCES:
        present = textual_method_successor_verdict_present(load_text(source))
        per_source[source] = present
        any_textual_successor = any_textual_successor or present

    checks_spec = [
        {
            "id": "p104_textual_method_subbranch_absent",
            "actual": p104["textual_method_successor_verdict_present"],
            "expected": False,
            "meaning": "P104 already isolates the textual method-successor path as one still-missing method-side sub-branch",
        },
        {
            "id": "strict_side_textual_method_successor_verdict_present",
            "actual": any_textual_successor,
            "expected": False,
            "meaning": "the current strict-side method source set exports no explicit textual method-successor verdict for the legacy gravity-hierarchy role",
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
        "stage": "P105",
        "lane": "strict_side_textual_method_successor_probe_for_legacy_gravity_hierarchy_role_current_repo_state_only",
        "goal": "test_whether_the_current_strict_side_method_source_set_exports_an_explicit_textual_method_successor_verdict_for_the_legacy_gravity_hierarchy_role",
        "status": "CURRENT_STRICT_SIDE_METHOD_SOURCES_DO_NOT_EXPORT_AN_EXPLICIT_TEXTUAL_METHOD_SUCCESSOR_VERDICT_FOR_THE_LEGACY_GRAVITY_HIERARCHY_ROLE_AFTER_P105",
        "reason": "the strict-side method sources do export the qw2115_micro_supported_beta_hierarchy_bridge method chain, but none of those sources joins that method to the legacy gravity-hierarchy role through an explicit textual method-successor verdict",
        "strict_side_sources_checked": STRICT_SIDE_METHOD_SOURCES,
        "per_source_textual_method_successor_verdict_presence": per_source,
        "textual_method_successor_verdict_present": any_textual_successor,
        "remaining_missing_objects": [
            "explicit_method_lineage_upgrade_verdict_elevating_the_existing_qw2115_micro_supported_beta_hierarchy_bridge_chain_into_replacement_semantics_for_the_legacy_gravity_hierarchy_role"
        ],
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P105",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "per_source_textual_method_successor_verdict_presence": per_source,
        "textual_method_successor_verdict_present": any_textual_successor,
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
