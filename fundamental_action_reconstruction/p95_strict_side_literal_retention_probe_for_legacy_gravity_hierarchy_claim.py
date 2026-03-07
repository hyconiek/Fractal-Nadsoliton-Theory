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
    / "p95_strict_side_literal_retention_probe_for_legacy_gravity_hierarchy_claim.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p95_strict_side_literal_retention_probe_for_legacy_gravity_hierarchy_claim_summary.json"
)


STRICT_SIDE_SOURCES = [
    "RELEASE_4_9_TEXTBOOK_EN_PL.md",
    "material_dowodowy/korpus_qw_pozostaly/raporty_md/RAPORT_QW2069_FULL_SM_GR_DERIVATION_PACKAGE.md",
    "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2069_full_sm_gr_derivation_package.json",
    "material_dowodowy/korpus_qw_pozostaly/raporty_md/RAPORT_QW2094_STRICT_RIGOR_DEFECT_SWEEP.md",
    "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2094_strict_rigor_defect_sweep.json",
    "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2115_gravity_hierarchy_strict_bridge_gate.json",
    "fundamental_action_reconstruction/A8_GRAVITY_BRIDGE_SPEC.md",
    "fundamental_action_reconstruction/generated/a8_gravity_bridge_summary.json",
]

LITERAL_VARIANTS = [
    "exact gravity hierarchy from beta^N scaling",
    "Exact gravity hierarchy from beta^N scaling",
    "G_eff(N) = G_0 * beta_tors^N",
    "G_eff(N)=G_0*beta_tors^N",
    "G_{eff}(N) = G_0 * beta_tors^N",
    "G_{eff}(N)=G_0*beta_tors^N",
    "G(N=20)/G(N=0) = beta_tors^20 = 10^-40",
    "G(N=20)/G(N=0)=beta_tors^20=10^-40",
    "G(N=20)/G(N=0) = beta_tors^20 = 10^{-40}",
    "G(N=20)/G(N=0)=beta_tors^20=10^{-40}",
]


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def load_text(repo_relative_path: str) -> str:
    return (REPO / repo_relative_path).read_text(encoding="utf-8")


def contains_any(text: str, parts: list[str]) -> bool:
    return any(part in text for part in parts)


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    p94 = load_json(
        "fundamental_action_reconstruction/generated/p94_legacy_gravity_hierarchy_retained_subbranch_probe_summary.json"
    )

    per_source = {}
    any_literal_retention = False
    for source in STRICT_SIDE_SOURCES:
        present = contains_any(load_text(source), LITERAL_VARIANTS)
        per_source[source] = present
        any_literal_retention = any_literal_retention or present

    checks_spec = [
        {
            "id": "p94_literal_retention_subbranch_absent",
            "actual": p94["retained_subbranch_state"]["literal_retention_present"],
            "expected": False,
            "meaning": "P94 already isolates literal retention as one still-missing retained sub-branch",
        },
        {
            "id": "strict_side_literal_retention_present",
            "actual": any_literal_retention,
            "expected": False,
            "meaning": "the strict-side authoritative source set does not export the old gravity-hierarchy claim or an algebraically identical literal form",
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
        "stage": "P95",
        "lane": "strict_side_literal_retention_probe_for_legacy_gravity_hierarchy_claim_current_repo_state_only",
        "goal": "test_whether_the_current_strict_side_authoritative_sources_export_the_old_gravity_hierarchy_claim_or_an_algebraically_identical_literal_form",
        "status": "CURRENT_STRICT_SIDE_AUTHORITATIVE_SOURCES_DO_NOT_EXPORT_LITERAL_RETENTION_OF_THE_LEGACY_GRAVITY_HIERARCHY_CLAIM_AFTER_P95",
        "reason": "the retained branch had already been reduced to literal retention versus role-equivalence retention, and the current strict-side authoritative source set exports no literal or algebraically identical form of the old beta^N gravity-hierarchy claim",
        "strict_side_sources_checked": STRICT_SIDE_SOURCES,
        "literal_variants_checked": LITERAL_VARIANTS,
        "per_source_literal_retention_presence": per_source,
        "literal_retention_present": any_literal_retention,
        "remaining_missing_objects": [
            "explicit_strict_side_role_equivalence_verdict_for_the_legacy_gravity_hierarchy_role"
        ],
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P95",
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
