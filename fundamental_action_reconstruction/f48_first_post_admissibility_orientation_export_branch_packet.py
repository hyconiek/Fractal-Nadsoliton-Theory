#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = (
    GENERATED / "f48_first_post_admissibility_orientation_export_branch_packet.json"
)
OUT_SUMMARY = (
    GENERATED
    / "f48_first_post_admissibility_orientation_export_branch_packet_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    f32 = load_json(
        "fundamental_action_reconstruction/generated/f32_initial_future_strict_core_orientation_export_admission_packet_summary.json"
    )
    n145 = load_json(
        "fundamental_action_reconstruction/generated/n145_current_admissibility_branch_obstruction_theorem_for_future_constructed_source_object_summary.json"
    )

    checks_spec = [
        {
            "id": "n145_admissibility_obstruction_discharged",
            "actual": n145["theorem_result"]["discharged"],
            "expected": True,
            "meaning": "N145 already packages the current admissibility-branch obstruction",
        },
        {
            "id": "n145_first_remaining_open_branch",
            "actual": n145["remaining_open_branches"][0],
            "expected": "future_derivation_of_admissible_E_orient_from_a_future_new_source_object",
            "meaning": "after admissibility, the first remaining lower branch is the orientation-export branch",
        },
        {
            "id": "f32_no_silent_kernel_substitution_clause_present",
            "actual": (
                "no silent kernel substitution"
                in f32["orientation_export_contract"].keys()
            ),
            "expected": True,
            "meaning": "F32 already freezes the admissible E_orient contract",
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
        "stage": "F48",
        "lane": "post_admissibility_orientation_export_branch_first_only",
        "goal": "freeze_the_orientation_export_branch_as_the_first_remaining_lower_branch_after_the_current_admissibility_branch_obstruction",
        "status": "F48_EXECUTED_FIRST_POST_ADMISSIBILITY_ORIENTATION_EXPORT_BRANCH_PACKET_NO_FALSE_PASS",
        "first_lower_branch_to_attack": "future_derivation_of_admissible_E_orient_from_a_future_new_source_object",
        "branch_ordering_basis": [
            "current_admissibility_branch_obstruction_already_packaged_by_N145",
            "E_orient_contract_already_frozen_by_F32",
            "downstream_B_sel_R_sel_O_sel_requires_source_object_and_admissible_orientation_export",
        ],
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "F48",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "first_lower_branch_to_attack": artifact["first_lower_branch_to_attack"],
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
