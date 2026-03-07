#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = (
    GENERATED / "p71_strict_side_replaced_branch_probe_for_legacy_weinberg_role.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p71_strict_side_replaced_branch_probe_for_legacy_weinberg_role_summary.json"
)

STRICT_SIDE_SOURCES = [
    "RELEASE_4_9_TEXTBOOK_EN_PL.md",
    "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2069_full_sm_gr_derivation_package.json",
    "material_dowodowy/korpus_qw_pozostaly/raporty_md/RAPORT_QW2098_EW_SECONDARY_NONANCHOR_CLOSURE_GATE.md",
    "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2098_ew_secondary_nonanchor_closure_gate.json",
    "t1_nonanchor_observables_input_qw2085_2086.json",
]

LEGACY_MARKERS = [
    "legacy weinberg",
    "legacy weinberg-angle role",
    "alpha_geo/12",
    "sin^2(theta_W)=alpha_geo/12",
]
STRICT_SUCCESSOR_MARKERS = [
    "sin2_theta_w_mz",
    "sin2_theta_w_eff",
]
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


def replaced_verdict_present(text: str) -> bool:
    return (
        any(marker in text for marker in LEGACY_MARKERS)
        and any(marker in text for marker in STRICT_SUCCESSOR_MARKERS)
        and any(marker in text for marker in REPLACEMENT_MARKERS)
    )


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    n73 = load_json(
        "fundamental_action_reconstruction/generated/n73_current_legacy_weinberg_retained_branch_full_negative_closure_theorem_summary.json"
    )

    per_source = {}
    any_replaced_verdict = False
    for source in STRICT_SIDE_SOURCES:
        present = replaced_verdict_present(load_text(source))
        per_source[source] = present
        any_replaced_verdict = any_replaced_verdict or present

    checks_spec = [
        {
            "id": "n73_retained_branch_closed",
            "actual": n73["theorem_result"]["retained_branch_closed_negatively_on_current_repo_state"],
            "expected": True,
            "meaning": "N73 already closes the retained branch negatively on the current repo state",
        },
        {
            "id": "strict_side_replaced_verdict_present",
            "actual": any_replaced_verdict,
            "expected": False,
            "meaning": "the current strict-side source set exports no explicit replaced-branch verdict for the legacy Weinberg role",
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
        "stage": "P71",
        "lane": "strict_side_replaced_branch_probe_for_legacy_weinberg_role_current_repo_state_only",
        "goal": "test_whether_the_current_strict_side_source_set_exports_an_explicit_replaced_branch_verdict_for_the_legacy_weinberg_angle_role",
        "status": "CURRENT_STRICT_SIDE_SOURCES_DO_NOT_EXPORT_AN_EXPLICIT_REPLACED_BRANCH_VERDICT_FOR_THE_LEGACY_WEINBERG_ROLE_AFTER_P71",
        "reason": "the retained branch is already closed negatively, but the current strict-side source set still exports no explicit replacement verdict linking the legacy Weinberg role to a named strict successor semantics",
        "strict_side_sources_checked": STRICT_SIDE_SOURCES,
        "per_source_replaced_verdict_presence": per_source,
        "replaced_verdict_present": any_replaced_verdict,
        "remaining_missing_objects": [
            "explicit_strict_side_replaced_verdict_for_the_legacy_weinberg_angle_role_by_an_explicit_strict_successor_semantics"
        ],
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P71",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "per_source_replaced_verdict_presence": per_source,
        "replaced_verdict_present": any_replaced_verdict,
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
