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
    / "p74_strict_side_textual_object_successor_probe_for_legacy_weinberg_role.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p74_strict_side_textual_object_successor_probe_for_legacy_weinberg_role_summary.json"
)

STRICT_SIDE_OBJECT_SOURCES = [
    "RELEASE_4_9_TEXTBOOK_EN_PL.md",
    "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2068_sm_gr_parameter_registry.json",
    "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2069_full_sm_gr_derivation_package.json",
    "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2094_strict_rigor_defect_sweep.json",
    "material_dowodowy/korpus_qw_pozostaly/raporty_md/RAPORT_QW2098_EW_SECONDARY_NONANCHOR_CLOSURE_GATE.md",
    "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2098_ew_secondary_nonanchor_closure_gate.json",
]

LEGACY_MARKERS = [
    "legacy weinberg",
    "legacy weinberg-angle role",
    "alpha_geo/12",
    "sin^2(theta_W)=alpha_geo/12",
]
OBJECT_MARKERS = ["sin2_theta_w_mz"]
TEXTUAL_SUCCESSOR_MARKERS = [
    "successor object",
    "strict-side successor object",
    "replacing the legacy weinberg-angle role",
    "replacing the legacy weinberg role",
    "replaced by",
    "replacement semantics",
    "strict successor semantics",
]


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def load_text(repo_relative_path: str) -> str:
    return (REPO / repo_relative_path).read_text(encoding="utf-8")


def textual_object_successor_verdict_present(text: str) -> bool:
    return (
        any(marker in text for marker in LEGACY_MARKERS)
        and any(marker in text for marker in OBJECT_MARKERS)
        and any(marker in text for marker in TEXTUAL_SUCCESSOR_MARKERS)
    )


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    p73 = load_json(
        "fundamental_action_reconstruction/generated/p73_strict_side_object_successor_subbranch_probe_for_legacy_weinberg_role_summary.json"
    )

    per_source = {}
    any_textual_successor = False
    for source in STRICT_SIDE_OBJECT_SOURCES:
        present = textual_object_successor_verdict_present(load_text(source))
        per_source[source] = present
        any_textual_successor = any_textual_successor or present

    checks_spec = [
        {
            "id": "p73_textual_object_subbranch_absent",
            "actual": p73["textual_object_successor_verdict_present"],
            "expected": False,
            "meaning": "P73 already isolates the textual object-successor path as one still-missing object-side sub-branch",
        },
        {
            "id": "strict_side_textual_object_successor_verdict_present",
            "actual": any_textual_successor,
            "expected": False,
            "meaning": "the current strict-side object source set exports no explicit textual object-successor verdict for the legacy Weinberg role",
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
        "stage": "P74",
        "lane": "strict_side_textual_object_successor_probe_for_legacy_weinberg_role_current_repo_state_only",
        "goal": "test_whether_the_current_strict_side_object_source_set_exports_an_explicit_textual_object_successor_verdict_for_the_legacy_weinberg_role",
        "status": "CURRENT_STRICT_SIDE_OBJECT_SOURCES_DO_NOT_EXPORT_AN_EXPLICIT_TEXTUAL_OBJECT_SUCCESSOR_VERDICT_FOR_THE_LEGACY_WEINBERG_ROLE_AFTER_P74",
        "reason": "the strict-side object sources do export the sin2_theta_w_mz object chain, but none of those sources joins that object to the legacy Weinberg role through an explicit textual object-successor verdict",
        "strict_side_sources_checked": STRICT_SIDE_OBJECT_SOURCES,
        "per_source_textual_object_successor_verdict_presence": per_source,
        "textual_object_successor_verdict_present": any_textual_successor,
        "remaining_missing_objects": [
            "explicit_object_lineage_upgrade_verdict_elevating_the_existing_sin2_theta_w_mz_candidate_chain_into_replacement_semantics_for_the_legacy_weinberg_angle_role"
        ],
        "method_branch_still_open": [
            "explicit_method_successor_semantics_verdict_identifying_qw2098_sin2_from_nonanchor_ew_pole_chain_as_the_strict_side_successor_semantics_replacing_the_legacy_weinberg_angle_role"
        ],
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P74",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "per_source_textual_object_successor_verdict_presence": per_source,
        "textual_object_successor_verdict_present": any_textual_successor,
        "remaining_missing_objects": artifact["remaining_missing_objects"],
        "method_branch_still_open": artifact["method_branch_still_open"],
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
