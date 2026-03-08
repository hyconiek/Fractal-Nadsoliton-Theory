#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "n129_last_positive_branch_reduced_to_one_initial_source_plus_orientation_package_theorem_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    p118 = load_json(
        "fundamental_action_reconstruction/generated/p118_initial_source_admission_and_orientation_export_contract_probe_summary.json"
    )

    checks_spec = [
        {
            "id": "p118_initial_package_reduced",
            "actual": p118["target_state"][
                "last_positive_branch_reduced_to_one_initial_package"
            ],
            "expected": True,
            "meaning": "P118 already proves the last branch is reduced to one initial source-plus-orientation package",
        },
        {
            "id": "package_source_name",
            "actual": p118["target_state"]["initial_package"]["strict_core_source_object"],
            "expected": "S_sel_int",
            "meaning": "the initial package begins with the future source object",
        },
        {
            "id": "package_orientation_name",
            "actual": p118["target_state"]["initial_package"]["internal_orientation_export"],
            "expected": "E_orient",
            "meaning": "the initial package includes the future orientation export",
        },
    ]

    checks = []
    mismatches = []
    for item in checks_spec:
        ok = item["actual"] == item["expected"]
        checks.append(
            {
                "id": item["id"],
                "actual": item["actual"],
                "expected": item["expected"],
                "pass": ok,
                "meaning": item["meaning"],
            }
        )
        if not ok:
            mismatches.append(item["id"])

    if mismatches:
        summary = {
            "step": "N129",
            "status": "N129_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_INITIAL_SOURCE_PLUS_ORIENTATION_PACKAGE_STATE",
            "scope": "initial_source_plus_orientation_package_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N129",
            "status": "N129_DISCHARGED_LAST_POSITIVE_BRANCH_REDUCED_TO_ONE_INITIAL_SOURCE_PLUS_ORIENTATION_PACKAGE_THEOREM_NO_FALSE_PASS",
            "scope": "initial_source_plus_orientation_package_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "last_positive_branch_reduced_to_one_initial_package": True,
                "initial_package": p118["target_state"]["initial_package"],
                "downstream_chain_left_open": p118["target_state"][
                    "downstream_chain_left_open"
                ],
                "full_closure_pass": False,
            },
            "remaining_open_branches": [
                "future_construction_of_admissible_S_sel_int",
                "future_derivation_of_admissible_E_orient_from_S_sel_int",
                "future_completion_of_B_sel_R_sel_O_sel_after_seed_package",
            ],
            "hard_limits": [
                "initial_package_not_yet_constructed",
                "downstream_chain_not_yet_constructed",
                "no_strict_core_selector_closure",
                "no_QW2191_discharge",
                "no_ToE_closure",
            ],
        }

    OUT.write_text(
        json.dumps(summary, indent=2, ensure_ascii=True) + "\n",
        encoding="ascii",
    )
    print(OUT)


if __name__ == "__main__":
    main()
