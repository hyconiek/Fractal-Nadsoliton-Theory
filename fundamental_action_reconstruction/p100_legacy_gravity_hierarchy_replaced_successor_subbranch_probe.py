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
    / "p100_legacy_gravity_hierarchy_replaced_successor_subbranch_probe.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p100_legacy_gravity_hierarchy_replaced_successor_subbranch_probe_summary.json"
)

STRICT_SIDE_SOURCES = [
    "RELEASE_4_9_TEXTBOOK_EN_PL.md",
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
OBJECT_MARKERS = ["gravity_hierarchy_beta20"]
METHOD_MARKERS = ["qw2115_micro_supported_beta_hierarchy_bridge"]
REPLACEMENT_MARKERS = [
    "replaced",
    "replacement",
    "replaced by",
    "superseded",
    "supersession",
    "strict successor semantics",
    "successor semantics",
]


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def load_text(repo_relative_path: str) -> str:
    return (REPO / repo_relative_path).read_text(encoding="utf-8")


def object_replaced_verdict_present(text: str) -> bool:
    lowered = text.lower()
    return (
        any(marker in lowered for marker in LEGACY_MARKERS)
        and any(marker in lowered for marker in OBJECT_MARKERS)
        and any(marker in lowered for marker in REPLACEMENT_MARKERS)
    )


def method_replaced_verdict_present(text: str) -> bool:
    lowered = text.lower()
    return (
        any(marker in lowered for marker in LEGACY_MARKERS)
        and any(marker in lowered for marker in METHOD_MARKERS)
        and any(marker in lowered for marker in REPLACEMENT_MARKERS)
    )


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    f24 = load_json(
        "fundamental_action_reconstruction/generated/f24_legacy_gravity_hierarchy_replaced_branch_refinement_packet_summary.json"
    )

    per_source_object = {}
    per_source_method = {}
    any_object_verdict = False
    any_method_verdict = False
    for source in STRICT_SIDE_SOURCES:
        text = load_text(source)
        object_present = object_replaced_verdict_present(text)
        method_present = method_replaced_verdict_present(text)
        per_source_object[source] = object_present
        per_source_method[source] = method_present
        any_object_verdict = any_object_verdict or object_present
        any_method_verdict = any_method_verdict or method_present

    checks_spec = [
        {
            "id": "f24_object_candidate_present",
            "actual": f24["candidate_state"]["strict_object_candidate_present"],
            "expected": True,
            "meaning": "F24 already confirms that the strict side exports gravity_hierarchy_beta20",
        },
        {
            "id": "f24_method_candidate_present",
            "actual": f24["candidate_state"]["strict_method_candidate_present"],
            "expected": True,
            "meaning": "F24 already confirms that the strict side exports qw2115_micro_supported_beta_hierarchy_bridge",
        },
        {
            "id": "object_replaced_verdict_present",
            "actual": any_object_verdict,
            "expected": False,
            "meaning": "the current repo exports no explicit object-successor replaced verdict for the legacy gravity-hierarchy role",
        },
        {
            "id": "method_replaced_verdict_present",
            "actual": any_method_verdict,
            "expected": False,
            "meaning": "the current repo exports no explicit method-successor-semantics replaced verdict for the legacy gravity-hierarchy role",
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
        "stage": "P100",
        "lane": "legacy_gravity_hierarchy_replaced_successor_subbranch_probe_current_repo_state_only",
        "goal": "test_whether_the_current_repo_exports_either_the_object_successor_or_method_successor_replaced_subbranch_for_the_legacy_gravity_hierarchy_role",
        "status": "CURRENT_REPO_EXPORTS_STRICT_GRAVITY_HIERARCHY_SUCCESSOR_CANDIDATES_BUT_NEITHER_OBJECT_NOR_METHOD_REPLACED_SUCCESSOR_VERDICT_FOR_THE_LEGACY_GRAVITY_HIERARCHY_ROLE_AFTER_P100",
        "reason": "the strict side exports both gravity_hierarchy_beta20 and qw2115_micro_supported_beta_hierarchy_bridge as real successor-candidate structures, but no current source upgrades either one into an explicit replaced verdict for the legacy gravity-hierarchy role",
        "strict_side_sources_checked": STRICT_SIDE_SOURCES,
        "per_source_object_replaced_verdict_presence": per_source_object,
        "per_source_method_replaced_verdict_presence": per_source_method,
        "object_replaced_verdict_present": any_object_verdict,
        "method_replaced_verdict_present": any_method_verdict,
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
        "stage": "P100",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "per_source_object_replaced_verdict_presence": per_source_object,
        "per_source_method_replaced_verdict_presence": per_source_method,
        "object_replaced_verdict_present": any_object_verdict,
        "method_replaced_verdict_present": any_method_verdict,
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
