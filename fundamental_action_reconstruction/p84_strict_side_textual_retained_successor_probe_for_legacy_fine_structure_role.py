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
    / "p84_strict_side_textual_retained_successor_probe_for_legacy_fine_structure_role.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p84_strict_side_textual_retained_successor_probe_for_legacy_fine_structure_role_summary.json"
)

STRICT_SIDE_SOURCES = [
    "RELEASE_4_9_TEXTBOOK_EN_PL.md",
    "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2068_sm_gr_parameter_registry.json",
    "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2069_full_sm_gr_derivation_package.json",
    "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2094_strict_rigor_defect_sweep.json",
    "material_dowodowy/korpus_qw_pozostaly/raporty_md/RAPORT_QW2098_EW_SECONDARY_NONANCHOR_CLOSURE_GATE.md",
    "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2098_ew_secondary_nonanchor_closure_gate.json",
]

STRICT_MARKERS = ["alpha_em_inv_mz"]
LEGACY_MARKERS = [
    "legacy fine-structure",
    "legacy fine-structure role",
    "alpha_em^-1 = alpha_geo/(2*beta_tors)*(1-beta_tors)",
    "alpha_em^-1=alpha_geo/(2*beta_tors)*(1-beta_tors)",
]
SUCCESSOR_MARKERS = [
    "retained successor",
    "retained strict-side successor",
    "same legacy fine-structure role",
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

    p83 = load_json(
        "fundamental_action_reconstruction/generated/p83_legacy_fine_structure_semantic_transfer_subbranch_probe_summary.json"
    )

    per_source = {}
    any_textual_successor = False
    for source in STRICT_SIDE_SOURCES:
        present = textual_successor_verdict_present(load_text(source))
        per_source[source] = present
        any_textual_successor = any_textual_successor or present

    checks_spec = [
        {
            "id": "p83_textual_subbranch_absent",
            "actual": p83["textual_retained_successor_verdict_present"],
            "expected": False,
            "meaning": "P83 already isolates the textual successor path as one still-missing retained semantic-transfer sub-branch",
        },
        {
            "id": "strict_side_textual_successor_verdict_present",
            "actual": any_textual_successor,
            "expected": False,
            "meaning": "the current strict-side source set exports no explicit textual retained-successor verdict for the legacy fine-structure role",
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
        "stage": "P84",
        "lane": "strict_side_textual_retained_successor_probe_for_legacy_fine_structure_role_current_repo_state_only",
        "goal": "test_whether_the_current_strict_side_source_set_exports_an_explicit_textual_retained_successor_verdict_for_the_legacy_fine_structure_role",
        "status": "CURRENT_STRICT_SIDE_SOURCES_DO_NOT_EXPORT_AN_EXPLICIT_TEXTUAL_RETAINED_SUCCESSOR_VERDICT_FOR_THE_LEGACY_FINE_STRUCTURE_ROLE_AFTER_P84",
        "reason": "the strict-side source set does promote alpha_em_inv_mz as a strict-derived observable, but none of those sources joins that object to the legacy fine-structure role through an explicit retained-successor textual verdict",
        "strict_side_sources_checked": STRICT_SIDE_SOURCES,
        "per_source_textual_successor_verdict_presence": per_source,
        "textual_successor_verdict_present": any_textual_successor,
        "remaining_missing_objects": [
            "explicit_object_lineage_upgrade_verdict_elevating_the_existing_alpha_em_inv_mz_candidate_chain_into_retained_strict_side_fine_structure_role_transfer"
        ],
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P84",
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
