#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "n134_next_constructive_move_reduced_to_one_first_source_seed_admissibility_upgrade_target_theorem_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    p123 = load_json(
        "fundamental_action_reconstruction/generated/p123_first_source_seed_admissibility_upgrade_target_probe_summary.json"
    )

    checks_spec = [
        {
            "id": "p123_target_reduced",
            "actual": p123["target_state"]["next_constructive_move_reduced_to_one_first_admissibility_upgrade_target"],
            "expected": True,
            "meaning": "P123 already proves the next constructive move is reduced to one first admissibility-upgrade target",
        },
        {
            "id": "candidate_seed_name",
            "actual": p123["target_state"]["admissibility_upgrade_target"]["candidate_seed_name"],
            "expected": "S_sel_int_candidate_seed_v0",
            "meaning": "the admissibility-upgrade target is anchored on the one candidate seed instance",
        },
        {
            "id": "upgrade_target_name",
            "actual": p123["target_state"]["admissibility_upgrade_target"]["target_admissible_source_name"],
            "expected": "S_sel_int",
            "meaning": "the admissibility-upgrade target is anchored on future S_sel_int",
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
            "step": "N134",
            "status": "N134_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_FIRST_SOURCE_SEED_ADMISSIBILITY_UPGRADE_TARGET_STATE",
            "scope": "first_source_seed_admissibility_upgrade_target_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N134",
            "status": "N134_DISCHARGED_NEXT_CONSTRUCTIVE_MOVE_REDUCED_TO_ONE_FIRST_SOURCE_SEED_ADMISSIBILITY_UPGRADE_TARGET_THEOREM_NO_FALSE_PASS",
            "scope": "first_source_seed_admissibility_upgrade_target_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "next_constructive_move_reduced_to_one_first_admissibility_upgrade_target": True,
                "admissibility_upgrade_target": p123["target_state"]["admissibility_upgrade_target"],
                "later_open_branches": p123["target_state"]["later_open_branches"],
                "full_closure_pass": False,
            },
            "remaining_open_branches": [
                "future_clause_by_clause_attempted_admissibility_upgrade_of_S_sel_int_candidate_seed_v0",
                "future_derivation_of_admissible_E_orient_from_an_upgraded_source_seed",
                "future_completion_of_B_sel_R_sel_O_sel_after_seed_upgrade",
            ],
            "hard_limits": [
                "admissibility_upgrade_not_yet_successful",
                "admissible_S_sel_int_not_yet_constructed",
                "admissible_E_orient_not_yet_constructed",
                "downstream_chain_not_yet_constructed",
                "no_strict_core_selector_closure",
                "no_QW2191_discharge",
                "no_ToE_closure",
            ],
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
