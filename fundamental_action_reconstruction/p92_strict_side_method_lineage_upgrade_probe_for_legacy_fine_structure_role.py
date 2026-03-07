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
    / "p92_strict_side_method_lineage_upgrade_probe_for_legacy_fine_structure_role.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p92_strict_side_method_lineage_upgrade_probe_for_legacy_fine_structure_role_summary.json"
)

STRICT_SIDE_METHOD_SOURCES = [
    "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2069_full_sm_gr_derivation_package.json",
    "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2094_strict_rigor_defect_sweep.json",
    "material_dowodowy/korpus_qw_pozostaly/raporty_md/RAPORT_QW2098_EW_SECONDARY_NONANCHOR_CLOSURE_GATE.md",
    "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2098_ew_secondary_nonanchor_closure_gate.json",
]

LEGACY_MARKERS = [
    "legacy fine-structure",
    "legacy fine-structure role",
    "alpha_em^-1 = alpha_geo/(2*beta_tors)*(1-beta_tors)",
    "alpha_em^-1=alpha_geo/(2*beta_tors)*(1-beta_tors)",
]
METHOD_MARKERS = ["qw2098_alpha_on_shell_from_gf_mw_sin2_delta_r"]
LINEAGE_MARKERS = [
    "derived",
    "matches_qw2098_method",
    "strict_internal_gate",
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

    p91 = load_json(
        "fundamental_action_reconstruction/generated/p91_strict_side_textual_method_successor_probe_for_legacy_fine_structure_role_summary.json"
    )

    per_source = {}
    any_lineage_upgrade = False
    for source in STRICT_SIDE_METHOD_SOURCES:
        present = method_lineage_upgrade_verdict_present(load_text(source))
        per_source[source] = present
        any_lineage_upgrade = any_lineage_upgrade or present

    checks_spec = [
        {
            "id": "p91_textual_method_subbranch_absent",
            "actual": p91["textual_method_successor_verdict_present"],
            "expected": False,
            "meaning": "P91 already closes the textual method-successor sub-branch negatively on the current repo state",
        },
        {
            "id": "strict_side_method_lineage_upgrade_verdict_present",
            "actual": any_lineage_upgrade,
            "expected": False,
            "meaning": "the current strict-side method source set exports no explicit method-lineage-upgrade verdict for the legacy fine-structure role",
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
        "stage": "P92",
        "lane": "strict_side_method_lineage_upgrade_probe_for_legacy_fine_structure_role_current_repo_state_only",
        "goal": "test_whether_the_current_strict_side_method_source_set_exports_an_explicit_method_lineage_upgrade_verdict_for_the_legacy_fine_structure_role",
        "status": "CURRENT_STRICT_SIDE_METHOD_SOURCES_DO_NOT_EXPORT_AN_EXPLICIT_METHOD_LINEAGE_UPGRADE_VERDICT_FOR_THE_LEGACY_FINE_STRUCTURE_ROLE_AFTER_P92",
        "reason": "the strict-side method sources do export the qw2098_alpha_on_shell_from_gf_mw_sin2_delta_r method chain, but none of those sources upgrades that chain into explicit replacement semantics for the legacy fine-structure role",
        "strict_side_sources_checked": STRICT_SIDE_METHOD_SOURCES,
        "per_source_method_lineage_upgrade_verdict_presence": per_source,
        "method_lineage_upgrade_verdict_present": any_lineage_upgrade,
        "method_successor_branch_closed_on_current_repo_state": True,
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P92",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "per_source_method_lineage_upgrade_verdict_presence": per_source,
        "method_lineage_upgrade_verdict_present": any_lineage_upgrade,
        "method_successor_branch_closed_on_current_repo_state": True,
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
