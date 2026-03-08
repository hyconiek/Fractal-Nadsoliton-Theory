#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "n150_only_remaining_positive_move_class_reduced_to_one_additive_construction_contract_theorem_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    n149 = load_json(
        "fundamental_action_reconstruction/generated/n149_current_repo_constructive_selector_frontier_exhaustion_theorem_summary.json"
    )
    p137 = load_json(
        "fundamental_action_reconstruction/generated/p137_next_positive_move_class_reduction_probe_summary.json"
    )

    checks_spec = [
        {
            "id": "n149_constructive_frontier_exhausted",
            "actual": n149["theorem_result"]["constructive_selector_frontier_exhausted_on_current_repo_state"],
            "expected": True,
            "meaning": "N149 already exhausts the current constructive selector frontier",
        },
        {
            "id": "n149_only_remaining_positive_move_class",
            "actual": n149["theorem_result"]["only_remaining_positive_move_class"],
            "expected": "future_genuinely_additive_new_strict_core_source_object_construction",
            "meaning": "N149 already fixes the only remaining positive move class",
        },
        {
            "id": "p137_reduced_to_one_additive_contract",
            "actual": p137["target_state"]["reduced_to_one_explicit_minimal_additive_construction_contract"],
            "expected": True,
            "meaning": "P137 already reduces that positive move class to one explicit additive-construction contract",
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
            "step": "N150",
            "status": "N150_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_POSITIVE_MOVE_CLASS_REDUCTION_STATE",
            "scope": "only_remaining_positive_move_class_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N150",
            "status": "N150_DISCHARGED_ONLY_REMAINING_POSITIVE_MOVE_CLASS_REDUCED_TO_ONE_ADDITIVE_CONSTRUCTION_CONTRACT_THEOREM_NO_FALSE_PASS",
            "scope": "only_remaining_positive_move_class_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "only_remaining_positive_move_class": "future_genuinely_additive_new_strict_core_source_object_construction",
                "reduced_to_one_explicit_minimal_additive_construction_contract": True,
                "full_closure_pass": False,
            },
            "remaining_open_branches": [
                "future_attempted_genuinely_additive_new_strict_core_source_object_construction"
            ],
            "hard_limits": [
                "constructed_source_object_not_yet_exported",
                "admissible_S_sel_int_not_yet_constructed",
                "admissible_E_orient_not_yet_constructed",
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
