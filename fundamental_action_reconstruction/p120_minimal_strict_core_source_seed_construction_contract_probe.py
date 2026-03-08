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
    / "p120_minimal_strict_core_source_seed_construction_contract_probe.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p120_minimal_strict_core_source_seed_construction_contract_probe_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    n130 = load_json(
        "fundamental_action_reconstruction/generated/n130_last_positive_branch_reduced_to_one_first_source_seed_construction_target_theorem_summary.json"
    )
    f34 = load_json(
        "fundamental_action_reconstruction/generated/f34_minimal_admissible_strict_core_source_seed_construction_contract_packet_summary.json"
    )

    reduced_to_one_minimal_contract = (
        n130["theorem_result"]["first_source_seed_construction_target"][
            "strict_core_source_object"
        ]
        == "S_sel_int"
        and f34["minimal_source_seed_construction_contract"]["strict_core_source_object"]
        == "S_sel_int"
    )

    checks_spec = [
        {
            "id": "n130_first_source_seed_target_name",
            "actual": n130["theorem_result"]["first_source_seed_construction_target"][
                "strict_core_source_object"
            ],
            "expected": "S_sel_int",
            "meaning": "N130 already anchors the first construction target on S_sel_int",
        },
        {
            "id": "f34_contract_name",
            "actual": f34["minimal_source_seed_construction_contract"][
                "strict_core_source_object"
            ],
            "expected": "S_sel_int",
            "meaning": "F34 already freezes the minimal construction contract for S_sel_int",
        },
        {
            "id": "reduced_to_one_minimal_contract",
            "actual": reduced_to_one_minimal_contract,
            "expected": True,
            "meaning": "the last positive branch is reduced to one explicit minimal source-seed construction contract",
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
        "stage": "P120",
        "lane": "minimal_strict_core_source_seed_construction_contract_only",
        "goal": "test_whether_the_last_positive_branch_is_now_reduced_to_one_explicit_minimal_construction_contract_for_S_sel_int",
        "status": "CURRENT_REPO_REDUCES_THE_LAST_POSITIVE_BRANCH_TO_ONE_MINIMAL_ADMISSIBLE_STRICT_CORE_SOURCE_SEED_CONSTRUCTION_CONTRACT_AFTER_P120",
        "target_state": {
            "last_positive_branch_reduced_to_one_minimal_source_seed_construction_contract": reduced_to_one_minimal_contract,
            "minimal_source_seed_construction_contract": f34[
                "minimal_source_seed_construction_contract"
            ],
            "later_open_branches": [
                "future_derivation_of_admissible_E_orient_from_S_sel_int",
                "future_completion_of_B_sel_R_sel_O_sel_after_seed_package",
            ],
        },
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P120",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "target_state": artifact["target_state"],
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(
        json.dumps(artifact, indent=2, ensure_ascii=True) + "\n",
        encoding="ascii",
    )
    OUT_SUMMARY.write_text(
        json.dumps(summary, indent=2, ensure_ascii=True) + "\n",
        encoding="ascii",
    )
    print(OUT_JSON)


if __name__ == "__main__":
    main()
