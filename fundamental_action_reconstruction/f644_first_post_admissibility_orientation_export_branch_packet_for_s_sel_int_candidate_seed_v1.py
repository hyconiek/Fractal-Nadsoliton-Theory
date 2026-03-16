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
    / "f644_first_post_admissibility_orientation_export_branch_packet_for_s_sel_int_candidate_seed_v1.json"
)
OUT_SUMMARY = (
    GENERATED
    / "f644_first_post_admissibility_orientation_export_branch_packet_for_s_sel_int_candidate_seed_v1_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    f32 = load_json(
        "fundamental_action_reconstruction/generated/f32_initial_future_strict_core_orientation_export_admission_packet_summary.json"
    )
    n534 = load_json(
        "fundamental_action_reconstruction/generated/n534_current_admissibility_branch_obstruction_theorem_for_future_constructed_source_object_for_s_sel_int_candidate_seed_v1_summary.json"
    )

    checks_spec = [
        {
            "id": "n534_admissibility_obstruction_discharged",
            "actual": n534["theorem_result"]["discharged"],
            "expected": True,
            "meaning": "N534 packages the current seed-v1 admissibility-branch obstruction",
        },
        {
            "id": "n534_first_remaining_open_branch",
            "actual": n534["remaining_open_branches"][0],
            "expected": "future_derivation_of_admissible_E_orient_from_a_future_new_source_object",
            "meaning": "after admissibility, the first remaining lower branch is the orientation-export branch (seed v1)",
        },
        {
            "id": "f32_no_silent_kernel_substitution_clause_present",
            "actual": ("no silent kernel substitution" in f32["orientation_export_contract"].keys()),
            "expected": True,
            "meaning": "F32 freezes the admissible E_orient contract (kernel split discipline present)",
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
        "stage": "F644",
        "lane": "post_admissibility_orientation_export_branch_first_only",
        "goal": "freeze_the_orientation_export_branch_as_the_first_remaining_lower_branch_after_the_seed_v1_admissibility_branch_obstruction",
        "status": "F644_EXECUTED_FIRST_POST_ADMISSIBILITY_ORIENTATION_EXPORT_BRANCH_PACKET_FOR_SEED_V1_NO_FALSE_PASS",
        "first_lower_branch_to_attack": "future_derivation_of_admissible_E_orient_from_a_future_new_source_object",
        "branch_ordering_basis": [
            "current_admissibility_branch_obstruction_already_packaged_by_N534",
            "E_orient_contract_already_frozen_by_F32",
            "downstream_B_sel_R_sel_O_sel_requires_source_object_and_admissible_orientation_export",
        ],
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "F644",
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

