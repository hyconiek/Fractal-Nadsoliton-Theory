#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "f49_last_downstream_completion_branch_packet.json"
OUT_SUMMARY = (
    GENERATED / "f49_last_downstream_completion_branch_packet_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    n145 = load_json(
        "fundamental_action_reconstruction/generated/n145_current_admissibility_branch_obstruction_theorem_for_future_constructed_source_object_summary.json"
    )
    n146 = load_json(
        "fundamental_action_reconstruction/generated/n146_current_orientation_export_branch_obstruction_theorem_summary.json"
    )

    checks_spec = [
        {
            "id": "n145_admissibility_obstruction_discharged",
            "actual": n145["theorem_result"]["discharged"],
            "expected": True,
            "meaning": "N145 already packages the admissibility obstruction",
        },
        {
            "id": "n146_orientation_export_obstruction_discharged",
            "actual": n146["theorem_result"]["discharged"],
            "expected": True,
            "meaning": "N146 already packages the orientation-export obstruction",
        },
        {
            "id": "n146_last_remaining_open_branch",
            "actual": n146["remaining_open_branches"][0],
            "expected": "future_completion_of_B_sel_R_sel_O_sel_after_new_source_object_construction",
            "meaning": "after orientation-export, the only remaining lower branch is downstream completion",
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
        "stage": "F49",
        "lane": "last_downstream_completion_branch_only",
        "goal": "freeze_the_downstream_completion_branch_as_the_only_remaining_lower_branch_after_current_orientation_export_obstruction",
        "status": "F49_EXECUTED_LAST_DOWNSTREAM_COMPLETION_BRANCH_PACKET_NO_FALSE_PASS",
        "last_lower_branch_to_attack": "future_completion_of_B_sel_R_sel_O_sel_after_new_source_object_construction",
        "branch_ordering_basis": [
            "current_admissibility_branch_obstruction_already_packaged_by_N145",
            "current_orientation_export_branch_obstruction_already_packaged_by_N146",
            "downstream_requires_source_object_and_orientation_export_stages",
        ],
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "F49",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "last_lower_branch_to_attack": artifact["last_lower_branch_to_attack"],
        "branch_ordering_basis": artifact["branch_ordering_basis"],
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
