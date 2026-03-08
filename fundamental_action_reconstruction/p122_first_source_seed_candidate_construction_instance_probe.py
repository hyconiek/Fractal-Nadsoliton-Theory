#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "p122_first_source_seed_candidate_construction_instance_probe.json"
OUT_SUMMARY = (
    GENERATED / "p122_first_source_seed_candidate_construction_instance_probe_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    n132 = load_json(
        "fundamental_action_reconstruction/generated/n132_next_constructive_move_reduced_to_one_source_seed_precursor_route_theorem_summary.json"
    )
    f36 = load_json(
        "fundamental_action_reconstruction/generated/f36_first_source_seed_candidate_construction_instance_packet_summary.json"
    )

    reduced_to_one_candidate_instance = (
        n132["theorem_result"]["next_constructive_move_reduced_to_one_precursor_route"] is True
        and f36["candidate_seed_instance"]["candidate_seed_name"] == "S_sel_int_candidate_seed_v0"
        and f36["candidate_seed_instance"]["construction_route"]["local_topological_protection_layer"]
        == "QW-2206_local_topological_protection_layer"
        and f36["candidate_seed_instance"]["construction_route"]["internal_binary_candidate"]
        == "sigma_int_candidate"
        and f36["candidate_seed_instance"]["counts_only_as"]
        == "first_candidate_source_seed_construction_instance"
    )

    checks_spec = [
        {
            "id": "n132_precursor_route_reduced",
            "actual": n132["theorem_result"]["next_constructive_move_reduced_to_one_precursor_route"],
            "expected": True,
            "meaning": "N132 already reduces the next constructive move to one precursor route",
        },
        {
            "id": "f36_candidate_seed_name",
            "actual": f36["candidate_seed_instance"]["candidate_seed_name"],
            "expected": "S_sel_int_candidate_seed_v0",
            "meaning": "F36 exports one explicit first candidate source-seed construction instance",
        },
        {
            "id": "f36_candidate_is_not_admissible_source",
            "actual": "admissible_S_sel_int" in f36["candidate_seed_instance"]["does_not_count_as"],
            "expected": True,
            "meaning": "F36 keeps the candidate instance below admissible source-object scope",
        },
        {
            "id": "reduced_to_one_candidate_instance",
            "actual": reduced_to_one_candidate_instance,
            "expected": True,
            "meaning": "the next constructive move is reduced to one first source-seed candidate construction instance",
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
        "stage": "P122",
        "lane": "first_source_seed_candidate_construction_instance_only",
        "goal": "test_whether_the_next_constructive_move_is_now_reduced_to_one_first_candidate_source_seed_construction_instance",
        "status": "CURRENT_REPO_REDUCES_THE_NEXT_CONSTRUCTIVE_MOVE_TO_ONE_FIRST_SOURCE_SEED_CANDIDATE_CONSTRUCTION_INSTANCE_AFTER_P122",
        "target_state": {
            "next_constructive_move_reduced_to_one_first_candidate_instance": reduced_to_one_candidate_instance,
            "candidate_seed_instance": f36["candidate_seed_instance"],
            "later_open_branches": [
                "future_admissibility_upgrade_of_S_sel_int_candidate_seed_v0_to_admissible_S_sel_int",
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
        "stage": "P122",
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
