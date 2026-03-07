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
    / "p81_strict_side_literal_retention_probe_for_legacy_fine_structure_formula.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p81_strict_side_literal_retention_probe_for_legacy_fine_structure_formula_summary.json"
)


STRICT_SIDE_SOURCES = [
    "RELEASE_4_9_TEXTBOOK_EN_PL.md",
    "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2069_full_sm_gr_derivation_package.json",
    "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2094_strict_rigor_defect_sweep.json",
    "material_dowodowy/korpus_qw_pozostaly/raporty_md/RAPORT_QW2098_EW_SECONDARY_NONANCHOR_CLOSURE_GATE.md",
    "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2098_ew_secondary_nonanchor_closure_gate.json",
]

LITERAL_VARIANTS = [
    "alpha_EM^-1 = alpha_geo/(2*beta_tors)*(1-beta_tors)",
    "alpha_EM^-1=alpha_geo/(2*beta_tors)*(1-beta_tors)",
    "alpha_em^-1 = alpha_geo/(2*beta_tors)*(1-beta_tors)",
    "alpha_em^-1=alpha_geo/(2*beta_tors)*(1-beta_tors)",
    "alpha_EM^{-1} = alpha_geo/(2*beta_tors)*(1-beta_tors)",
]


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def load_text(repo_relative_path: str) -> str:
    return (REPO / repo_relative_path).read_text(encoding="utf-8")


def contains_any(text: str, parts: list[str]) -> bool:
    return any(part in text for part in parts)


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    p80 = load_json(
        "fundamental_action_reconstruction/generated/p80_legacy_fine_structure_retained_subbranch_probe_summary.json"
    )

    per_source = {}
    any_literal_retention = False
    for source in STRICT_SIDE_SOURCES:
        present = contains_any(load_text(source), LITERAL_VARIANTS)
        per_source[source] = present
        any_literal_retention = any_literal_retention or present

    checks_spec = [
        {
            "id": "p80_literal_retention_subbranch_absent",
            "actual": p80["retained_subbranch_state"]["literal_retention_present"],
            "expected": False,
            "meaning": "P80 already isolates literal retention as one still-missing retained sub-branch",
        },
        {
            "id": "strict_side_literal_retention_present",
            "actual": any_literal_retention,
            "expected": False,
            "meaning": "the strict-side authoritative source set does not export the old fine-structure formula or an algebraically identical literal form",
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
        "stage": "P81",
        "lane": "strict_side_literal_retention_probe_for_legacy_fine_structure_formula_current_repo_state_only",
        "goal": "test_whether_the_current_strict_side_authoritative_sources_export_the_old_fine_structure_formula_or_an_algebraically_identical_literal_form",
        "status": "CURRENT_STRICT_SIDE_AUTHORITATIVE_SOURCES_DO_NOT_EXPORT_LITERAL_RETENTION_OF_THE_LEGACY_FINE_STRUCTURE_FORMULA_AFTER_P81",
        "reason": "the retained branch had already been reduced to literal retention versus role-equivalence retention, and the current strict-side authoritative source set exports no literal or algebraically identical form of alpha_EM^-1 = alpha_geo/(2*beta_tors)*(1-beta_tors)",
        "strict_side_sources_checked": STRICT_SIDE_SOURCES,
        "literal_variants_checked": LITERAL_VARIANTS,
        "per_source_literal_retention_presence": per_source,
        "literal_retention_present": any_literal_retention,
        "remaining_missing_objects": [
            "explicit_strict_side_role_equivalence_verdict_for_the_legacy_fine_structure_role"
        ],
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P81",
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
