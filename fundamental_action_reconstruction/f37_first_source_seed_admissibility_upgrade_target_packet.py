#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "f37_first_source_seed_admissibility_upgrade_target_packet.json"
OUT_SUMMARY = (
    GENERATED / "f37_first_source_seed_admissibility_upgrade_target_packet_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    f34 = load_json(
        "fundamental_action_reconstruction/generated/f34_minimal_admissible_strict_core_source_seed_construction_contract_packet_summary.json"
    )
    n133 = load_json(
        "fundamental_action_reconstruction/generated/n133_next_constructive_move_reduced_to_one_first_source_seed_candidate_instance_theorem_summary.json"
    )

    checks_spec = [
        {
            "id": "n133_candidate_seed_name",
            "actual": n133["theorem_result"]["candidate_seed_instance"]["candidate_seed_name"],
            "expected": "S_sel_int_candidate_seed_v0",
            "meaning": "N133 fixes one first candidate source-seed instance",
        },
        {
            "id": "f34_contract_target_name",
            "actual": f34["minimal_source_seed_construction_contract"]["strict_core_source_object"],
            "expected": "S_sel_int",
            "meaning": "F34 fixes one minimal admissibility contract targeted at S_sel_int",
        },
        {
            "id": "f34_contract_strict_core_only",
            "actual": f34["minimal_source_seed_construction_contract"]["strict_core_only_required"],
            "expected": True,
            "meaning": "the admissibility contract remains strict-core only",
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

    admissibility_upgrade_target = {
        "candidate_seed_name": n133["theorem_result"]["candidate_seed_instance"]["candidate_seed_name"],
        "target_admissible_source_name": f34["minimal_source_seed_construction_contract"]["strict_core_source_object"],
        "candidate_seed_instance": n133["theorem_result"]["candidate_seed_instance"],
        "minimal_admissibility_contract": f34["minimal_source_seed_construction_contract"],
        "counts_only_as": "first_source_seed_admissibility_upgrade_target",
        "does_not_count_as": [
            "successful_admissibility_upgrade",
            "admissible_S_sel_int",
            "admissible_E_orient",
            "B_sel",
            "R_sel",
            "O_sel",
            "strict_core_selector_closure",
        ],
    }

    artifact = {
        "stage": "F37",
        "lane": "first_source_seed_admissibility_upgrade_target_only",
        "goal": "freeze_the_narrowest_first_attempted_admissibility_upgrade_target_for_S_sel_int_candidate_seed_v0",
        "status": "F37_EXECUTED_FIRST_SOURCE_SEED_ADMISSIBILITY_UPGRADE_TARGET_PACKET_NO_FALSE_PASS",
        "admissibility_upgrade_target": admissibility_upgrade_target,
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "F37",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "admissibility_upgrade_target": artifact["admissibility_upgrade_target"],
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()
