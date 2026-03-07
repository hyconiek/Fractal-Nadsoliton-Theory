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
    / "p104_strict_side_method_successor_subbranch_probe_for_legacy_gravity_hierarchy_role.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p104_strict_side_method_successor_subbranch_probe_for_legacy_gravity_hierarchy_role_summary.json"
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
LINEAGE_MARKERS = [
    "derived",
    "strict bridge gate",
    "gravity_hierarchy_beta20",
]
UPGRADE_MARKERS = [
    "upgraded to replacement",
    "replacement semantics",
    "replaced by",
    "strict successor semantics",
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


def method_lineage_upgrade_verdict_present(text: str) -> bool:
    lowered = text.lower()
    return (
        any(marker in lowered for marker in LEGACY_MARKERS)
        and any(marker in lowered for marker in METHOD_MARKERS)
        and any(marker in lowered for marker in LINEAGE_MARKERS)
        and any(marker in lowered for marker in UPGRADE_MARKERS)
    )


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    f26 = load_json(
        "fundamental_action_reconstruction/generated/f26_legacy_gravity_hierarchy_method_successor_semantics_refinement_packet_summary.json"
    )

    per_source_textual = {}
    per_source_lineage = {}
    any_textual = False
    any_lineage = False
    for source in STRICT_SIDE_METHOD_SOURCES:
        text = load_text(source)
        textual = textual_method_successor_verdict_present(text)
        lineage = method_lineage_upgrade_verdict_present(text)
        per_source_textual[source] = textual
        per_source_lineage[source] = lineage
        any_textual = any_textual or textual
        any_lineage = any_lineage or lineage

    checks_spec = [
        {
            "id": "f26_method_chain_present",
            "actual": f26["candidate_state"]["method_chain_present"],
            "expected": True,
            "meaning": "F26 already confirms that the repo exports the qw2115 gravity method chain",
        },
        {
            "id": "textual_method_successor_verdict_present",
            "actual": any_textual,
            "expected": False,
            "meaning": "the current repo exports no explicit textual method-successor verdict for the legacy gravity-hierarchy role",
        },
        {
            "id": "method_lineage_upgrade_verdict_present",
            "actual": any_lineage,
            "expected": False,
            "meaning": "the current repo exports no explicit method-lineage-upgrade verdict for the legacy gravity-hierarchy role",
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
        "stage": "P104",
        "lane": "strict_side_method_successor_subbranch_probe_for_legacy_gravity_hierarchy_role_current_repo_state_only",
        "goal": "test_whether_the_current_repo_exports_either_the_textual_method_successor_or_method_lineage_upgrade_subbranch_for_the_legacy_gravity_hierarchy_role",
        "status": "CURRENT_REPO_EXPORTS_THE_STRICT_QW2115_GRAVITY_METHOD_CHAIN_BUT_NEITHER_TEXTUAL_METHOD_SUCCESSOR_NOR_METHOD_LINEAGE_UPGRADE_VERDICT_FOR_THE_LEGACY_GRAVITY_HIERARCHY_ROLE_AFTER_P104",
        "reason": "the repo exports the qw2115_micro_supported_beta_hierarchy_bridge method chain, but no current source upgrades that chain into replacement semantics for the legacy gravity-hierarchy role",
        "strict_side_sources_checked": STRICT_SIDE_METHOD_SOURCES,
        "per_source_textual_method_successor_verdict_presence": per_source_textual,
        "per_source_method_lineage_upgrade_verdict_presence": per_source_lineage,
        "textual_method_successor_verdict_present": any_textual,
        "method_lineage_upgrade_verdict_present": any_lineage,
        "remaining_missing_objects": [
            "explicit_textual_method_successor_semantics_verdict_identifying_qw2115_micro_supported_beta_hierarchy_bridge_as_the_strict_side_successor_semantics_replacing_the_legacy_gravity_hierarchy_role",
            "explicit_method_lineage_upgrade_verdict_elevating_the_existing_qw2115_micro_supported_beta_hierarchy_bridge_chain_into_replacement_semantics_for_the_legacy_gravity_hierarchy_role",
        ],
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P104",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "per_source_textual_method_successor_verdict_presence": per_source_textual,
        "per_source_method_lineage_upgrade_verdict_presence": per_source_lineage,
        "textual_method_successor_verdict_present": any_textual,
        "method_lineage_upgrade_verdict_present": any_lineage,
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
