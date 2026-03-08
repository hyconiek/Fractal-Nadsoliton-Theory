#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "f35_future_strict_core_source_seed_precursor_route_packet.json"
OUT_SUMMARY = (
    GENERATED / "f35_future_strict_core_source_seed_precursor_route_packet_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    b4 = load_json("fundamental_action_reconstruction/generated/b4_minimal_sigma_int_candidate_summary.json")
    b5 = load_json("fundamental_action_reconstruction/generated/b5_sigma_int_local_stability_audit_summary.json")
    n126 = load_json(
        "fundamental_action_reconstruction/generated/n126_current_repo_exports_no_admissible_strict_core_internal_selector_source_object_theorem_summary.json"
    )
    n131 = load_json(
        "fundamental_action_reconstruction/generated/n131_last_positive_branch_reduced_to_one_minimal_strict_core_source_seed_construction_contract_theorem_summary.json"
    )

    checks_spec = [
        {
            "id": "b4_sigma_candidate_exists",
            "actual": b4["b4"]["candidate_name"],
            "expected": "sigma_int_candidate",
            "meaning": "B4 exports one canonical internal binary candidate",
        },
        {
            "id": "b5_local_support_present",
            "actual": b5["b5"]["findings"][0]["status"],
            "expected": "supported_partial",
            "meaning": "B5 exports local topological support for sigma_int_candidate",
        },
        {
            "id": "n126_no_current_admissible_source_object",
            "actual": n126["theorem_result"]["admissible_strict_core_internal_selector_source_object_present"],
            "expected": False,
            "meaning": "N126 blocks treating any current object as already admissible source",
        },
        {
            "id": "n131_next_move_anchored_on_S_sel_int",
            "actual": n131["theorem_result"]["minimal_source_seed_construction_contract"]["strict_core_source_object"],
            "expected": "S_sel_int",
            "meaning": "N131 keeps the next constructive move anchored on future S_sel_int",
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

    precursor_route = {
        "local_topological_protection_layer": "QW-2206_local_topological_protection_layer",
        "internal_binary_candidate": "sigma_int_candidate",
        "future_construction_target": "S_sel_int",
        "counts_only_as_precursor_route": True,
        "does_not_count_as_constructed_source_object": True,
    }

    artifact = {
        "stage": "F35",
        "lane": "future_strict_core_source_seed_precursor_route_only",
        "goal": "freeze_the_narrowest_existing_internal_precursor_route_for_a_future_attempted_construction_of_S_sel_int",
        "status": "F35_EXECUTED_FUTURE_STRICT_CORE_SOURCE_SEED_PRECURSOR_ROUTE_PACKET_NO_FALSE_PASS",
        "precursor_route": precursor_route,
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "F35",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "precursor_route": artifact["precursor_route"],
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()
