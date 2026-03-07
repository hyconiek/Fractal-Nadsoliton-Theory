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
    / "p98_strict_side_textual_retained_successor_probe_for_legacy_gravity_hierarchy_role.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p98_strict_side_textual_retained_successor_probe_for_legacy_gravity_hierarchy_role_summary.json"
)

STRICT_SIDE_SOURCES = [
    "RELEASE_4_9_TEXTBOOK_EN_PL.md",
    "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2068_sm_gr_parameter_registry.json",
    "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2069_full_sm_gr_derivation_package.json",
    "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2115_gravity_hierarchy_strict_bridge_gate.json",
    "fundamental_action_reconstruction/A8_GRAVITY_BRIDGE_SPEC.md",
    "fundamental_action_reconstruction/generated/a8_gravity_bridge_summary.json",
]

STRICT_MARKERS = ["gravity_hierarchy_beta20"]
LEGACY_MARKERS = [
    "legacy gravity-hierarchy",
    "legacy gravity-hierarchy role",
    "exact gravity hierarchy from beta^N scaling",
    "beta^N scaling",
]
SUCCESSOR_MARKERS = [
    "retained successor",
    "retained strict-side successor",
    "same legacy gravity-hierarchy role",
    "successor semantics",
    "retained role",
    "legacy-to-strict",
    "role-equivalence",
    "role equivalence",
]


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def load_text(repo_relative_path: str) -> str:
    return (REPO / repo_relative_path).read_text(encoding="utf-8")


def textual_successor_verdict_present(text: str) -> bool:
    lowered = text.lower()
    return (
        any(marker in lowered for marker in STRICT_MARKERS)
        and any(marker in lowered for marker in LEGACY_MARKERS)
        and any(marker in lowered for marker in SUCCESSOR_MARKERS)
    )


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    p97 = load_json(
        "fundamental_action_reconstruction/generated/p97_legacy_gravity_hierarchy_semantic_transfer_subbranch_probe_summary.json"
    )

    per_source = {}
    any_textual_successor = False
    for source in STRICT_SIDE_SOURCES:
        present = textual_successor_verdict_present(load_text(source))
        per_source[source] = present
        any_textual_successor = any_textual_successor or present

    checks_spec = [
        {
            "id": "p97_textual_subbranch_absent",
            "actual": p97["textual_retained_successor_verdict_present"],
            "expected": False,
            "meaning": "P97 already isolates the textual successor path as one still-missing retained semantic-transfer sub-branch",
        },
        {
            "id": "strict_side_textual_successor_verdict_present",
            "actual": any_textual_successor,
            "expected": False,
            "meaning": "the current strict-side source set exports no explicit textual retained-successor verdict for the legacy gravity-hierarchy role",
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
        "stage": "P98",
        "lane": "strict_side_textual_retained_successor_probe_for_legacy_gravity_hierarchy_role_current_repo_state_only",
        "goal": "test_whether_the_current_strict_side_source_set_exports_an_explicit_textual_retained_successor_verdict_for_the_legacy_gravity_hierarchy_role",
        "status": "CURRENT_STRICT_SIDE_SOURCES_DO_NOT_EXPORT_AN_EXPLICIT_TEXTUAL_RETAINED_SUCCESSOR_VERDICT_FOR_THE_LEGACY_GRAVITY_HIERARCHY_ROLE_AFTER_P98",
        "reason": "the strict-side source set does promote gravity_hierarchy_beta20 as a strict-derived observable, but none of those sources joins that object to the legacy gravity-hierarchy role through an explicit retained-successor textual verdict",
        "strict_side_sources_checked": STRICT_SIDE_SOURCES,
        "per_source_textual_successor_verdict_presence": per_source,
        "textual_successor_verdict_present": any_textual_successor,
        "remaining_missing_objects": [
            "explicit_object_lineage_upgrade_verdict_elevating_the_existing_gravity_hierarchy_beta20_candidate_chain_into_retained_strict_side_gravity_hierarchy_role_transfer"
        ],
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P98",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "per_source_textual_successor_verdict_presence": per_source,
        "textual_successor_verdict_present": any_textual_successor,
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
