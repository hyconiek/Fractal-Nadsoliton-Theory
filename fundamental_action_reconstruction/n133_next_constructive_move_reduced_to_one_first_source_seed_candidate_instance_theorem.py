#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "n133_next_constructive_move_reduced_to_one_first_source_seed_candidate_instance_theorem_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    p122 = load_json(
        "fundamental_action_reconstruction/generated/p122_first_source_seed_candidate_construction_instance_probe_summary.json"
    )

    checks_spec = [
        {
            "id": "p122_candidate_instance_reduced",
            "actual": p122["target_state"]["next_constructive_move_reduced_to_one_first_candidate_instance"],
            "expected": True,
            "meaning": "P122 already proves the next constructive move is reduced to one first candidate instance",
        },
        {
            "id": "candidate_seed_name",
            "actual": p122["target_state"]["candidate_seed_instance"]["candidate_seed_name"],
            "expected": "S_sel_int_candidate_seed_v0",
            "meaning": "the first candidate instance is explicitly named",
        },
        {
            "id": "candidate_not_counted_as_admissible_source",
            "actual": "admissible_S_sel_int"
            in p122["target_state"]["candidate_seed_instance"]["does_not_count_as"],
            "expected": True,
            "meaning": "the candidate instance is explicitly blocked from being counted as admissible S_sel_int",
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
            "step": "N133",
            "status": "N133_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_FIRST_SOURCE_SEED_CANDIDATE_INSTANCE_STATE",
            "scope": "first_source_seed_candidate_instance_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N133",
            "status": "N133_DISCHARGED_NEXT_CONSTRUCTIVE_MOVE_REDUCED_TO_ONE_FIRST_SOURCE_SEED_CANDIDATE_INSTANCE_THEOREM_NO_FALSE_PASS",
            "scope": "first_source_seed_candidate_instance_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "next_constructive_move_reduced_to_one_first_candidate_instance": True,
                "candidate_seed_instance": p122["target_state"]["candidate_seed_instance"],
                "later_open_branches": p122["target_state"]["later_open_branches"],
                "full_closure_pass": False,
            },
            "remaining_open_branches": [
                "future_admissibility_upgrade_of_S_sel_int_candidate_seed_v0_to_admissible_S_sel_int",
                "future_derivation_of_admissible_E_orient_from_an_upgraded_source_seed",
                "future_completion_of_B_sel_R_sel_O_sel_after_seed_upgrade",
            ],
            "hard_limits": [
                "candidate_seed_not_yet_admissible_source_object",
                "candidate_seed_not_yet_E_orient",
                "candidate_seed_not_yet_B_sel_R_sel_O_sel",
                "no_strict_core_selector_closure",
                "no_QW2191_discharge",
                "no_ToE_closure",
            ],
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
