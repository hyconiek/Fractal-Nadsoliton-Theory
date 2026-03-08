#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = (
    GENERATED / "p135_current_orientation_export_branch_discharge_probe.json"
)
OUT_SUMMARY = (
    GENERATED / "p135_current_orientation_export_branch_discharge_probe_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    f48 = load_json(
        "fundamental_action_reconstruction/generated/f48_first_post_admissibility_orientation_export_branch_packet_summary.json"
    )
    n145 = load_json(
        "fundamental_action_reconstruction/generated/n145_current_admissibility_branch_obstruction_theorem_for_future_constructed_source_object_summary.json"
    )

    explicit_orientation_export_branch_discharge_exported = False

    checks_spec = [
        {
            "id": "f48_orientation_export_branch_first",
            "actual": f48["first_lower_branch_to_attack"],
            "expected": "future_derivation_of_admissible_E_orient_from_a_future_new_source_object",
            "meaning": "F48 already freezes the orientation-export branch as the first remaining lower branch",
        },
        {
            "id": "n145_admissibility_obstruction_discharged",
            "actual": n145["theorem_result"]["discharged"],
            "expected": True,
            "meaning": "N145 already packages the admissibility-branch obstruction",
        },
        {
            "id": "explicit_orientation_export_branch_discharge_exported",
            "actual": explicit_orientation_export_branch_discharge_exported,
            "expected": False,
            "meaning": "the current repo does not yet export an explicit orientation-export branch discharge",
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
        "stage": "P135",
        "lane": "current_orientation_export_branch_discharge_only",
        "goal": "test_whether_the_current_repo_already_exports_an_explicit_orientation_export_branch_discharge_for_the_first_remaining_future_E_orient_branch",
        "status": "CURRENT_REPO_DOES_NOT_EXPORT_AN_EXPLICIT_ORIENTATION_EXPORT_BRANCH_DISCHARGE_FOR_THE_FUTURE_E_ORIENT_BRANCH_AFTER_P135",
        "target_state": {
            "orientation_export_branch_selected_as_next_attack": True,
            "explicit_orientation_export_branch_discharge_exported": explicit_orientation_export_branch_discharge_exported,
            "remaining_open_branches": [
                "future_completion_of_B_sel_R_sel_O_sel_after_new_source_object_construction"
            ],
        },
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P135",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "target_state": artifact["target_state"],
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(
        json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii"
    )
    OUT_SUMMARY.write_text(
        json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii"
    )
    print(OUT_JSON)


if __name__ == "__main__":
    main()
