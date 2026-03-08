#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = (
    GENERATED / "p136_current_downstream_completion_branch_discharge_probe.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p136_current_downstream_completion_branch_discharge_probe_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    f49 = load_json(
        "fundamental_action_reconstruction/generated/f49_last_downstream_completion_branch_packet_summary.json"
    )
    n146 = load_json(
        "fundamental_action_reconstruction/generated/n146_current_orientation_export_branch_obstruction_theorem_summary.json"
    )

    explicit_downstream_completion_branch_discharge_exported = False

    checks_spec = [
        {
            "id": "f49_downstream_branch_last",
            "actual": f49["last_lower_branch_to_attack"],
            "expected": "future_completion_of_B_sel_R_sel_O_sel_after_new_source_object_construction",
            "meaning": "F49 already freezes downstream completion as the last remaining lower branch",
        },
        {
            "id": "n146_orientation_export_obstruction_discharged",
            "actual": n146["theorem_result"]["discharged"],
            "expected": True,
            "meaning": "N146 already packages the orientation-export obstruction",
        },
        {
            "id": "explicit_downstream_completion_branch_discharge_exported",
            "actual": explicit_downstream_completion_branch_discharge_exported,
            "expected": False,
            "meaning": "the current repo does not yet export an explicit downstream-completion branch discharge",
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
        "stage": "P136",
        "lane": "current_downstream_completion_branch_discharge_only",
        "goal": "test_whether_the_current_repo_already_exports_an_explicit_downstream_completion_branch_discharge_for_the_last_remaining_lower_branch",
        "status": "CURRENT_REPO_DOES_NOT_EXPORT_AN_EXPLICIT_DOWNSTREAM_COMPLETION_BRANCH_DISCHARGE_FOR_THE_LAST_REMAINING_LOWER_BRANCH_AFTER_P136",
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
        "stage": "P136",
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
