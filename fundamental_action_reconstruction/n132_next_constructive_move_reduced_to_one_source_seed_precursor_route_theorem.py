#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED / "n132_next_constructive_move_reduced_to_one_source_seed_precursor_route_theorem_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    p121 = load_json(
        "fundamental_action_reconstruction/generated/p121_future_strict_core_source_seed_precursor_route_probe_summary.json"
    )

    checks_spec = [
        {
            "id": "p121_precursor_route_reduced",
            "actual": p121["target_state"][
                "next_constructive_move_reduced_to_one_precursor_route"
            ],
            "expected": True,
            "meaning": "P121 already proves the next constructive move is reduced to one precursor route",
        },
        {
            "id": "precursor_target_name",
            "actual": p121["target_state"]["precursor_route"]["future_construction_target"],
            "expected": "S_sel_int",
            "meaning": "the precursor route is explicitly targeted at future S_sel_int",
        },
        {
            "id": "precursor_internal_candidate_name",
            "actual": p121["target_state"]["precursor_route"]["internal_binary_candidate"],
            "expected": "sigma_int_candidate",
            "meaning": "the precursor route is explicitly anchored on sigma_int_candidate",
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
            "step": "N132",
            "status": "N132_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_SOURCE_SEED_PRECURSOR_ROUTE_STATE",
            "scope": "source_seed_precursor_route_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N132",
            "status": "N132_DISCHARGED_NEXT_CONSTRUCTIVE_MOVE_REDUCED_TO_ONE_SOURCE_SEED_PRECURSOR_ROUTE_THEOREM_NO_FALSE_PASS",
            "scope": "source_seed_precursor_route_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "next_constructive_move_reduced_to_one_precursor_route": True,
                "precursor_route": p121["target_state"]["precursor_route"],
                "later_open_branches": p121["target_state"]["later_open_branches"],
                "full_closure_pass": False,
            },
            "remaining_open_branches": [
                "future_attempted_construction_of_admissible_S_sel_int_from_the_precursor_route",
                "future_derivation_of_admissible_E_orient_from_S_sel_int",
                "future_completion_of_B_sel_R_sel_O_sel_after_seed_package",
            ],
            "hard_limits": [
                "source_seed_not_yet_constructed",
                "precursor_route_not_yet_upgraded_to_source_object",
                "orientation_export_not_yet_constructed",
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
