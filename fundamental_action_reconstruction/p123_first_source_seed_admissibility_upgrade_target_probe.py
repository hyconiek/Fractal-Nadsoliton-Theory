#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "p123_first_source_seed_admissibility_upgrade_target_probe.json"
OUT_SUMMARY = (
    GENERATED / "p123_first_source_seed_admissibility_upgrade_target_probe_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    n133 = load_json(
        "fundamental_action_reconstruction/generated/n133_next_constructive_move_reduced_to_one_first_source_seed_candidate_instance_theorem_summary.json"
    )
    f34 = load_json(
        "fundamental_action_reconstruction/generated/f34_minimal_admissible_strict_core_source_seed_construction_contract_packet_summary.json"
    )
    f37 = load_json(
        "fundamental_action_reconstruction/generated/f37_first_source_seed_admissibility_upgrade_target_packet_summary.json"
    )

    reduced_to_one_admissibility_upgrade_target = (
        n133["theorem_result"]["next_constructive_move_reduced_to_one_first_candidate_instance"] is True
        and f34["minimal_source_seed_construction_contract"]["strict_core_source_object"] == "S_sel_int"
        and f37["admissibility_upgrade_target"]["candidate_seed_name"] == "S_sel_int_candidate_seed_v0"
        and f37["admissibility_upgrade_target"]["target_admissible_source_name"] == "S_sel_int"
        and f37["admissibility_upgrade_target"]["counts_only_as"]
        == "first_source_seed_admissibility_upgrade_target"
    )

    checks_spec = [
        {
            "id": "n133_first_candidate_instance_fixed",
            "actual": n133["theorem_result"]["candidate_seed_instance"]["candidate_seed_name"],
            "expected": "S_sel_int_candidate_seed_v0",
            "meaning": "N133 already fixes the first candidate source-seed instance",
        },
        {
            "id": "f34_contract_target_is_S_sel_int",
            "actual": f34["minimal_source_seed_construction_contract"]["strict_core_source_object"],
            "expected": "S_sel_int",
            "meaning": "F34 already fixes the minimal admissibility contract target",
        },
        {
            "id": "f37_target_is_not_successful_upgrade",
            "actual": "successful_admissibility_upgrade"
            in f37["admissibility_upgrade_target"]["does_not_count_as"],
            "expected": True,
            "meaning": "F37 keeps the target below a successful admissibility upgrade claim",
        },
        {
            "id": "reduced_to_one_admissibility_upgrade_target",
            "actual": reduced_to_one_admissibility_upgrade_target,
            "expected": True,
            "meaning": "the next constructive move is reduced to one first source-seed admissibility-upgrade target",
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
        "stage": "P123",
        "lane": "first_source_seed_admissibility_upgrade_target_only",
        "goal": "test_whether_the_next_constructive_move_is_now_reduced_to_one_first_source_seed_admissibility_upgrade_target",
        "status": "CURRENT_REPO_REDUCES_THE_NEXT_CONSTRUCTIVE_MOVE_TO_ONE_FIRST_SOURCE_SEED_ADMISSIBILITY_UPGRADE_TARGET_AFTER_P123",
        "target_state": {
            "next_constructive_move_reduced_to_one_first_admissibility_upgrade_target": reduced_to_one_admissibility_upgrade_target,
            "admissibility_upgrade_target": f37["admissibility_upgrade_target"],
            "later_open_branches": [
                "future_clause_by_clause_attempted_admissibility_upgrade_of_S_sel_int_candidate_seed_v0",
                "future_derivation_of_admissible_E_orient_from_an_upgraded_source_seed",
                "future_completion_of_B_sel_R_sel_O_sel_after_seed_upgrade",
            ],
        },
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P123",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "target_state": artifact["target_state"],
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()
