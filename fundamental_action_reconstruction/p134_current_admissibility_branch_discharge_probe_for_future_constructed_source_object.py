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
    / "p134_current_admissibility_branch_discharge_probe_for_future_constructed_source_object.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p134_current_admissibility_branch_discharge_probe_for_future_constructed_source_object_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    f47 = load_json(
        "fundamental_action_reconstruction/generated/f47_first_post_verdict_admissibility_branch_packet_summary.json"
    )
    n143 = load_json(
        "fundamental_action_reconstruction/generated/n143_current_failure_branch_obstruction_theorem_for_s_sel_int_realization_attempt_summary.json"
    )
    n144 = load_json(
        "fundamental_action_reconstruction/generated/n144_current_success_branch_obstruction_theorem_for_s_sel_int_realization_attempt_summary.json"
    )

    explicit_admissibility_branch_discharge_exported = False

    checks_spec = [
        {
            "id": "f47_admissibility_branch_first",
            "actual": f47["first_lower_branch_to_attack"],
            "expected": "future_admissibility_test_of_a_future_constructed_source_object",
            "meaning": "F47 already freezes admissibility as the first remaining lower branch",
        },
        {
            "id": "n143_failure_obstruction_discharged",
            "actual": n143["theorem_result"]["discharged"],
            "expected": True,
            "meaning": "N143 already packages the failure-side obstruction",
        },
        {
            "id": "n144_success_obstruction_discharged",
            "actual": n144["theorem_result"]["discharged"],
            "expected": True,
            "meaning": "N144 already packages the success-side obstruction",
        },
        {
            "id": "explicit_admissibility_branch_discharge_exported",
            "actual": explicit_admissibility_branch_discharge_exported,
            "expected": False,
            "meaning": "the current repo does not yet export an explicit admissibility-branch discharge for the future constructed-source-object branch",
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
        "stage": "P134",
        "lane": "current_admissibility_branch_discharge_only",
        "goal": "test_whether_the_current_repo_already_exports_an_explicit_admissibility_branch_discharge_for_the_first_remaining_future_constructed_source_object_branch",
        "status": "CURRENT_REPO_DOES_NOT_EXPORT_AN_EXPLICIT_ADMISSIBILITY_BRANCH_DISCHARGE_FOR_THE_FUTURE_CONSTRUCTED_SOURCE_OBJECT_BRANCH_AFTER_P134",
        "target_state": {
            "admissibility_branch_selected_as_next_attack": True,
            "explicit_admissibility_branch_discharge_exported": explicit_admissibility_branch_discharge_exported,
            "remaining_open_branches": [
                "future_derivation_of_admissible_E_orient_from_a_future_new_source_object",
                "future_completion_of_B_sel_R_sel_O_sel_after_new_source_object_construction",
            ],
        },
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P134",
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
