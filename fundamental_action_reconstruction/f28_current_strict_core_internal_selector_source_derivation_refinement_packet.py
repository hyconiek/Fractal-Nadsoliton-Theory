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
    / "f28_current_strict_core_internal_selector_source_derivation_refinement_packet.json"
)
OUT_SUMMARY = (
    GENERATED
    / "f28_current_strict_core_internal_selector_source_derivation_refinement_packet_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    b2 = load_json(
        "fundamental_action_reconstruction/generated/b2_internal_orientation_datum_source_audit_summary.json"
    )
    n4 = load_json(
        "fundamental_action_reconstruction/generated/n4_current_repo_psi0_strict_core_nonderivation_theorem_summary.json"
    )
    n5 = load_json(
        "fundamental_action_reconstruction/generated/n5_current_strict_core_psi0_route_obstruction_theorem_summary.json"
    )
    n6 = load_json(
        "fundamental_action_reconstruction/generated/n6_current_strict_core_fr_route_nonderivation_theorem_summary.json"
    )
    n7 = load_json(
        "fundamental_action_reconstruction/generated/n7_current_strict_core_sigma_int_residual_datum_nonderivation_theorem_summary.json"
    )
    n8 = load_json(
        "fundamental_action_reconstruction/generated/n8_current_strict_core_sigma_int_residual_datum_obstruction_after_target_slot_export_theorem_summary.json"
    )

    generic_hidden_source_branch_closed = (
        b2["b2"]["strict_internal_selector_derivations_found"] == 0
    )
    psi0_branch_closed = (
        not n4["current_repo_theorem_result"][
            "current_repo_strict_core_psi0_selector_source_present"
        ]
        and not n5["theorem_result"]["current_strict_core_psi0_route_closes_selector"]
    )
    fr_branch_closed = not n6["theorem_result"][
        "strict_core_fr_route_serves_as_internal_selector_source"
    ]
    sigma_int_bridge_branch_closed = (
        not n7["theorem_result"][
            "strict_core_sigma_int_route_serves_as_internal_residual_datum_bridge"
        ]
        and not n8["theorem_result"][
            "strict_core_sigma_int_route_derives_residual_datum_bridge"
        ]
    )

    checks_spec = [
        {
            "id": "b2_generic_hidden_source_branch_closed",
            "actual": generic_hidden_source_branch_closed,
            "expected": True,
            "meaning": "B2 already closes the generic hidden internal-orientation branch negatively on the current repo state",
        },
        {
            "id": "n4_n5_psi0_branch_closed",
            "actual": psi0_branch_closed,
            "expected": True,
            "meaning": "N4 and N5 already close the current strict-core psi0 branch negatively on the current repo state",
        },
        {
            "id": "n6_fr_branch_closed",
            "actual": fr_branch_closed,
            "expected": True,
            "meaning": "N6 already closes the current strict-core FR route negatively on the current repo state",
        },
        {
            "id": "n7_n8_sigma_int_bridge_branch_closed",
            "actual": sigma_int_bridge_branch_closed,
            "expected": True,
            "meaning": "N7 and N8 already close the current strict-core sigma-int bridge route negatively on the current repo state",
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
        "stage": "F28",
        "lane": "current_strict_core_internal_selector_source_refinement_only",
        "goal": "refine_the_last_remaining_higher_order_strict_core_frontier_into_the_narrowest_current_repo_state_source_branches",
        "status": "F28_EXECUTED_CURRENT_STRICT_CORE_INTERNAL_SELECTOR_SOURCE_DERIVATION_REFINEMENT_PACKET_NO_FALSE_PASS",
        "branch_state": {
            "generic_hidden_source_branch_closed_negatively_on_current_repo_state": generic_hidden_source_branch_closed,
            "psi0_branch_closed_negatively_on_current_repo_state": psi0_branch_closed,
            "fr_branch_closed_negatively_on_current_repo_state": fr_branch_closed,
            "sigma_int_bridge_branch_closed_negatively_on_current_repo_state": sigma_int_bridge_branch_closed,
        },
        "remaining_missing_objects": [
            "explicit_strict_core_internal_selector_source_derivation_discharge",
        ],
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "F28",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "branch_state": artifact["branch_state"],
        "remaining_missing_objects": artifact["remaining_missing_objects"],
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
