#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "f36_first_source_seed_candidate_construction_instance_packet.json"
OUT_SUMMARY = (
    GENERATED / "f36_first_source_seed_candidate_construction_instance_packet_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    b4 = load_json("fundamental_action_reconstruction/generated/b4_minimal_sigma_int_candidate_summary.json")
    b5 = load_json("fundamental_action_reconstruction/generated/b5_sigma_int_local_stability_audit_summary.json")
    n132 = load_json(
        "fundamental_action_reconstruction/generated/n132_next_constructive_move_reduced_to_one_source_seed_precursor_route_theorem_summary.json"
    )

    checks_spec = [
        {
            "id": "n132_precursor_target_is_S_sel_int",
            "actual": n132["theorem_result"]["precursor_route"]["future_construction_target"],
            "expected": "S_sel_int",
            "meaning": "N132 keeps the future source target anchored on S_sel_int",
        },
        {
            "id": "n132_precursor_candidate_is_sigma_int_candidate",
            "actual": n132["theorem_result"]["precursor_route"]["internal_binary_candidate"],
            "expected": "sigma_int_candidate",
            "meaning": "N132 keeps the precursor route anchored on sigma_int_candidate",
        },
        {
            "id": "b4_candidate_name",
            "actual": b4["b4"]["candidate_name"],
            "expected": "sigma_int_candidate",
            "meaning": "B4 exports sigma_int_candidate as the canonical internal binary candidate",
        },
        {
            "id": "b4_local_value",
            "actual": b4["b4"]["local_unit_sector_value"],
            "expected": "-1",
            "meaning": "B4 fixes the local unit-sector seed value for sigma_int_candidate",
        },
        {
            "id": "b5_local_stability_support",
            "actual": b5["b5"]["findings"][0]["status"],
            "expected": "supported_partial",
            "meaning": "B5 keeps only local topological support and does not overclaim gauge-complete safety",
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

    candidate_seed_instance = {
        "candidate_seed_name": "S_sel_int_candidate_seed_v0",
        "construction_route": {
            "local_topological_protection_layer": "QW-2206_local_topological_protection_layer",
            "internal_binary_candidate": "sigma_int_candidate",
            "future_construction_target": "S_sel_int",
        },
        "assembled_seed_data": {
            "protected_local_sector": "local_B_tilde_1_topological_sector",
            "topological_protection_status": "PASS_PARTIAL_LOCAL_PROTECTION_ONLY",
            "binary_internal_datum_name": "sigma_int_candidate",
            "binary_internal_datum_codomain": b4["b4"]["candidate_codomain"],
            "binary_internal_datum_local_value": b4["b4"]["local_unit_sector_value"],
            "local_stability_support_status": b5["b5"]["findings"][0]["status"],
        },
        "counts_only_as": "first_candidate_source_seed_construction_instance",
        "does_not_count_as": [
            "admissible_S_sel_int",
            "strict_core_internal_selector_source_object",
            "admissible_E_orient",
            "B_sel",
            "R_sel",
            "O_sel",
            "strict_core_selector_closure",
        ],
    }

    artifact = {
        "stage": "F36",
        "lane": "first_source_seed_candidate_construction_instance_only",
        "goal": "freeze_the_narrowest_first_candidate_construction_instance_on_the_fixed_precursor_route_to_future_S_sel_int",
        "status": "F36_EXECUTED_FIRST_SOURCE_SEED_CANDIDATE_CONSTRUCTION_INSTANCE_PACKET_NO_FALSE_PASS",
        "candidate_seed_instance": candidate_seed_instance,
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "F36",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "candidate_seed_instance": artifact["candidate_seed_instance"],
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()
