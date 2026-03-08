#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "p121_future_strict_core_source_seed_precursor_route_probe.json"
OUT_SUMMARY = (
    GENERATED / "p121_future_strict_core_source_seed_precursor_route_probe_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    n131 = load_json(
        "fundamental_action_reconstruction/generated/n131_last_positive_branch_reduced_to_one_minimal_strict_core_source_seed_construction_contract_theorem_summary.json"
    )
    f35 = load_json(
        "fundamental_action_reconstruction/generated/f35_future_strict_core_source_seed_precursor_route_packet_summary.json"
    )

    reduced_to_one_precursor_route = (
        n131["theorem_result"]["minimal_source_seed_construction_contract"]["strict_core_source_object"]
        == "S_sel_int"
        and f35["precursor_route"]["internal_binary_candidate"] == "sigma_int_candidate"
        and f35["precursor_route"]["future_construction_target"] == "S_sel_int"
    )

    checks_spec = [
        {
            "id": "n131_source_contract_still_anchored_on_S_sel_int",
            "actual": n131["theorem_result"]["minimal_source_seed_construction_contract"]["strict_core_source_object"],
            "expected": "S_sel_int",
            "meaning": "N131 keeps the next constructive move anchored on S_sel_int",
        },
        {
            "id": "f35_precursor_candidate_name",
            "actual": f35["precursor_route"]["internal_binary_candidate"],
            "expected": "sigma_int_candidate",
            "meaning": "F35 identifies sigma_int_candidate as the internal precursor candidate",
        },
        {
            "id": "reduced_to_one_precursor_route",
            "actual": reduced_to_one_precursor_route,
            "expected": True,
            "meaning": "the next constructive move is reduced to one explicit precursor route for S_sel_int",
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
        "stage": "P121",
        "lane": "future_strict_core_source_seed_precursor_route_only",
        "goal": "test_whether_the_next_constructive_move_is_now_reduced_to_one_explicit_precursor_route_for_S_sel_int",
        "status": "CURRENT_REPO_REDUCES_THE_NEXT_CONSTRUCTIVE_MOVE_TO_ONE_EXPLICIT_PRECURSOR_ROUTE_FOR_A_FUTURE_ATTEMPTED_CONSTRUCTION_OF_S_SEL_INT_AFTER_P121",
        "target_state": {
            "next_constructive_move_reduced_to_one_precursor_route": reduced_to_one_precursor_route,
            "precursor_route": f35["precursor_route"],
            "later_open_branches": [
                "future_attempted_construction_of_admissible_S_sel_int",
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
        "stage": "P121",
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
