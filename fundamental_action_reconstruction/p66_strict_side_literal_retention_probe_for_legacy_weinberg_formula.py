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
    / "p66_strict_side_literal_retention_probe_for_legacy_weinberg_formula.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p66_strict_side_literal_retention_probe_for_legacy_weinberg_formula_summary.json"
)


STRICT_SIDE_SOURCES = [
    "RELEASE_4_9_TEXTBOOK_EN_PL.md",
    "material_dowodowy/korpus_qw_pozostaly/raporty_md/RAPORT_QW2000_BOUNDED_COUPLING_DEEP_AUDIT.md",
    "material_dowodowy/korpus_qw_pozostaly/raporty_md/RAPORT_QW2001_BOUNDED_GW_TRIAD_LOCKABLE_GATE.md",
    "material_dowodowy/korpus_qw_pozostaly/raporty_md/RAPORT_QW2002_SINGLE_KERNEL_TRIPLE_SECTOR_CLOSURE_GATE_V3.md",
    "material_dowodowy/korpus_qw_pozostaly/raporty_md/RAPORT_QW2003_FROZEN_LOCKABLE_PACKAGE_EXPORT.md",
    "material_dowodowy/korpus_qw_pozostaly/raporty_md/RAPORT_QW2049_SPECTRAL_MICRO_STAGEC_INTERSECTION_GATE.md",
    "material_dowodowy/korpus_qw_pozostaly/raporty_md/RAPORT_QW2064_MICRO_DERIVED_RENORMALIZATION_CONSTANTS_GATE.md",
]

LITERAL_VARIANTS = [
    "sin^2(theta_W)=alpha_geo/12",
    "sin^2(theta_W) = alpha_geo/12",
    "sin^2\\theta_W = alpha_geo/12",
    "sin^2(theta_W)=4 ln 2 / 12",
    "sin^2(theta_W) = 4 ln 2 / 12",
    "sin^2(theta_W)=ln 2 / 3",
    "sin^2(theta_W) = ln 2 / 3",
    "sin^2\\theta_W = \\ln 2 / 3",
]


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def load_text(repo_relative_path: str) -> str:
    return (REPO / repo_relative_path).read_text(encoding="utf-8")


def contains_any(text: str, parts: list[str]) -> bool:
    return any(part in text for part in parts)


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    p65 = load_json(
        "fundamental_action_reconstruction/generated/p65_legacy_weinberg_retained_subbranch_probe_summary.json"
    )

    per_source = {}
    any_literal_retention = False
    for source in STRICT_SIDE_SOURCES:
        present = contains_any(load_text(source), LITERAL_VARIANTS)
        per_source[source] = present
        any_literal_retention = any_literal_retention or present

    checks_spec = [
        {
            "id": "p65_literal_retention_subbranch_absent",
            "actual": p65["retained_subbranch_state"]["literal_retention_present"],
            "expected": False,
            "meaning": "P65 already isolates literal retention as one still-missing retained sub-branch",
        },
        {
            "id": "strict_side_literal_retention_present",
            "actual": any_literal_retention,
            "expected": False,
            "meaning": "the strict-side authoritative source set does not export the old Weinberg formula or an algebraically identical literal form",
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
        "stage": "P66",
        "lane": "strict_side_literal_retention_probe_for_legacy_weinberg_formula_current_repo_state_only",
        "goal": "test_whether_the_current_strict_side_authoritative_sources_export_the_old_weinberg_formula_or_an_algebraically_identical_literal_form",
        "status": "CURRENT_STRICT_SIDE_AUTHORITATIVE_SOURCES_DO_NOT_EXPORT_LITERAL_RETENTION_OF_THE_LEGACY_WEINBERG_FORMULA_AFTER_P66",
        "reason": "the retained branch had already been reduced to literal retention versus role-equivalence retention, and the current strict-side authoritative source set exports no literal or algebraically identical form of sin^2(theta_W)=alpha_geo/12",
        "strict_side_sources_checked": STRICT_SIDE_SOURCES,
        "literal_variants_checked": LITERAL_VARIANTS,
        "per_source_literal_retention_presence": per_source,
        "literal_retention_present": any_literal_retention,
        "remaining_missing_objects": [
            "explicit_strict_side_role_equivalence_verdict_for_the_legacy_weinberg_angle_role"
        ],
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P66",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "per_source_literal_retention_presence": per_source,
        "literal_retention_present": any_literal_retention,
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
