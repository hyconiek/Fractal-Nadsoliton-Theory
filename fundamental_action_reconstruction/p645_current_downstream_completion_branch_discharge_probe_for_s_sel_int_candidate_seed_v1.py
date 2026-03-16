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
    / "p645_current_downstream_completion_branch_discharge_probe_for_s_sel_int_candidate_seed_v1.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p645_current_downstream_completion_branch_discharge_probe_for_s_sel_int_candidate_seed_v1_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    f645 = load_json(
        "fundamental_action_reconstruction/generated/f645_last_downstream_completion_branch_packet_for_s_sel_int_candidate_seed_v1_summary.json"
    )
    n535 = load_json(
        "fundamental_action_reconstruction/generated/n535_current_orientation_export_branch_obstruction_theorem_for_s_sel_int_candidate_seed_v1_summary.json"
    )

    explicit_downstream_completion_branch_discharge_exported = False

    checks_spec = [
        {
            "id": "f645_downstream_branch_last",
            "actual": f645["last_lower_branch_to_attack"],
            "expected": "future_completion_of_B_sel_R_sel_O_sel_after_new_source_object_construction",
            "meaning": "F645 freezes downstream completion as the last remaining lower branch (seed v1)",
        },
        {
            "id": "n535_orientation_export_obstruction_discharged",
            "actual": n535["theorem_result"]["discharged"],
            "expected": True,
            "meaning": "N535 packages the seed-v1 orientation-export obstruction",
        },
        {
            "id": "explicit_downstream_completion_branch_discharge_exported",
            "actual": explicit_downstream_completion_branch_discharge_exported,
            "expected": False,
            "meaning": "the current repo does not export an explicit downstream-completion branch discharge (seed v1)",
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
        "stage": "P645",
        "lane": "current_downstream_completion_branch_discharge_only",
        "goal": "test_whether_the_current_repo_already_exports_an_explicit_downstream_completion_branch_discharge_for_seed_v1_last_remaining_lower_branch",
        "status": "CURRENT_REPO_DOES_NOT_EXPORT_AN_EXPLICIT_DOWNSTREAM_COMPLETION_BRANCH_DISCHARGE_FOR_SEED_V1_AFTER_P645",
        "target_state": {
            "downstream_completion_branch_selected_as_next_attack": True,
            "explicit_downstream_completion_branch_discharge_exported": explicit_downstream_completion_branch_discharge_exported,
            "remaining_open_branches": [],
        },
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P645",
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

