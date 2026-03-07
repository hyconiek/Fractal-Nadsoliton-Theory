#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = (
    GENERATED
    / "p97_legacy_gravity_hierarchy_semantic_transfer_subbranch_probe.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p97_legacy_gravity_hierarchy_semantic_transfer_subbranch_probe_summary.json"
)

STRICT_SIDE_GRAVITY_SOURCES = [
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
    "exact gravity hierarchy from beta^N scaling",
    "beta^N scaling",
]
OBJECT_MARKERS = ["gravity_hierarchy_beta20"]
TEXTUAL_SUCCESSOR_MARKERS = [
    "retained successor",
    "retained strict-side successor",
    "role-equivalence",
    "role equivalence",
    "same role",
    "semantic transfer",
    "retained role",
]
LINEAGE_MARKERS = [
    "strict_internal_gate",
    "strict bridge gate",
    "strict-derived",
    "target parameter",
]
UPGRADE_MARKERS = [
    "upgraded to retained",
    "semantic transfer",
    "retained successor",
    "role-equivalence",
    "legacy-to-strict",
]


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def load_text(repo_relative_path: str) -> str:
    return (REPO / repo_relative_path).read_text(encoding="utf-8")


def textual_retained_successor_verdict_present(text: str) -> bool:
    lowered = text.lower()
    return (
        any(marker in lowered for marker in LEGACY_MARKERS)
        and any(marker in lowered for marker in OBJECT_MARKERS)
        and any(marker in lowered for marker in TEXTUAL_SUCCESSOR_MARKERS)
    )


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

    f23 = load_json(
        "fundamental_action_reconstruction/generated/f23_legacy_gravity_hierarchy_semantic_transfer_refinement_packet_summary.json"
    )

    per_source_textual = {}
    per_source_lineage = {}
    any_textual = False
    any_lineage = False
    for source in STRICT_SIDE_GRAVITY_SOURCES:
        text = load_text(source)
        textual = textual_retained_successor_verdict_present(text)
        lineage = object_lineage_upgrade_verdict_present(text)
        per_source_textual[source] = textual
        per_source_lineage[source] = lineage
        any_textual = any_textual or textual
        any_lineage = any_lineage or lineage

    checks_spec = [
        {
            "id": "f23_candidate_present",
            "actual": f23["candidate_state"]["strict_side_candidate_object_present"],
            "expected": True,
            "meaning": "F23 already confirms that a real strict-side gravity-hierarchy candidate object is present",
        },
        {
            "id": "f23_object_chain_present",
            "actual": f23["candidate_state"]["object_chain_present"],
            "expected": True,
            "meaning": "F23 already confirms that a real strict-side gravity_hierarchy_beta20 chain is present",
        },
        {
            "id": "textual_retained_successor_verdict_present",
            "actual": any_textual,
            "expected": False,
            "meaning": "the current repo exports no explicit textual retained-successor verdict for the legacy gravity-hierarchy role",
        },
        {
            "id": "object_lineage_upgrade_verdict_present",
            "actual": any_lineage,
            "expected": False,
            "meaning": "the current repo exports no explicit object-lineage-upgrade verdict for retained gravity-hierarchy transfer",
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
        "stage": "P97",
        "lane": "legacy_gravity_hierarchy_semantic_transfer_subbranch_probe_current_repo_state_only",
        "goal": "test_whether_the_current_repo_exports_either_the_textual_retained_successor_or_object_lineage_upgrade_subbranch_for_the_legacy_gravity_hierarchy_role",
        "status": "CURRENT_REPO_EXPORTS_NEITHER_TEXTUAL_SUCCESSOR_NOR_OBJECT_LINEAGE_UPGRADE_TRANSFER_VERDICT_FOR_THE_LEGACY_GRAVITY_HIERARCHY_ROLE_AFTER_P97",
        "reason": "the repo exports the gravity_hierarchy_beta20 candidate object and chain, but no current source upgrades that chain into retained semantic transfer for the legacy gravity-hierarchy role",
        "strict_side_sources_checked": STRICT_SIDE_GRAVITY_SOURCES,
        "per_source_textual_retained_successor_verdict_presence": per_source_textual,
        "per_source_object_lineage_upgrade_verdict_presence": per_source_lineage,
        "textual_retained_successor_verdict_present": any_textual,
        "object_lineage_upgrade_verdict_present": any_lineage,
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
        "stage": "P97",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "per_source_textual_retained_successor_verdict_presence": per_source_textual,
        "per_source_object_lineage_upgrade_verdict_presence": per_source_lineage,
        "textual_retained_successor_verdict_present": any_textual,
        "object_lineage_upgrade_verdict_present": any_lineage,
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
