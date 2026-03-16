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
    / "f318_first_current_strict_side_source_seed_strict_sigma_int_upgrade_candidate_instance_packet.json"
)
OUT_SUMMARY = (
    GENERATED
    / "f318_first_current_strict_side_source_seed_strict_sigma_int_upgrade_candidate_instance_packet_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    sigma_int_strict = load_json(
        "fundamental_action_reconstruction/generated/sigma_int_strict_derived_v1.json"
    )
    f36 = load_json(
        "fundamental_action_reconstruction/generated/f36_first_source_seed_candidate_construction_instance_packet_summary.json"
    )

    sigma_value = sigma_int_strict.get("value")
    sigma_in_codomain = sigma_value in (-1, 1, "+1", "-1")

    checks_spec = [
        {
            "id": "sigma_int_strict_derived_v1_object_name",
            "actual": sigma_int_strict.get("object"),
            "expected": "sigma_int_strict_derived_v1",
            "meaning": "strict sigma-int source-upgrade value object is present",
        },
        {
            "id": "sigma_int_strict_derived_v1_value_in_codomain",
            "actual": bool(sigma_in_codomain),
            "expected": True,
            "meaning": "sigma_int_strict_derived_v1 lies in {+1,-1}",
        },
        {
            "id": "f36_seed_v0_present",
            "actual": f36.get("candidate_seed_instance", {}).get("candidate_seed_name"),
            "expected": "S_sel_int_candidate_seed_v0",
            "meaning": "historical seed-v0 instance is exported (baseline continuity check)",
        },
        {
            "id": "f36_local_topological_protection_layer_label",
            "actual": f36.get("candidate_seed_instance", {})
            .get("construction_route", {})
            .get("local_topological_protection_layer"),
            "expected": "QW-2206_local_topological_protection_layer",
            "meaning": "the local topological protection layer label stays anchored on QW-2206",
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
        "candidate_seed_name": "S_sel_int_candidate_seed_v1",
        "construction_route": {
            "local_topological_protection_layer": "QW-2206_local_topological_protection_layer",
            "internal_binary_candidate": "sigma_int_strict_derived_v1",
            "future_construction_target": "S_sel_int",
        },
        "assembled_seed_data": {
            "protected_local_sector": "local_B_tilde_1_topological_sector",
            "local_topological_label_only": True,
            "binary_internal_datum_name": "sigma_int_strict_derived_v1",
            "binary_internal_datum_codomain": ["+1", "-1"],
            "binary_internal_datum_value": sigma_value,
            "binary_internal_datum_status": sigma_int_strict.get("status"),
            "binary_internal_datum_definition": sigma_int_strict.get("definition"),
            "binary_internal_datum_declared_domain": sigma_int_strict.get("declared_domain"),
            "binary_internal_datum_provenance": sigma_int_strict.get("provenance"),
        },
        "counts_only_as": "strict_sigma_int_upgraded_strict_side_source_seed_candidate_instance",
        "does_not_count_as": [
            "sigma_int_candidate_identification",
            "admissible_S_sel_int",
            "strict_core_internal_selector_source_object",
            "E_orient",
            "B_sel",
            "R_sel",
            "O_sel",
            "strict_core_selector_closure",
            "QW-2191_discharge",
            "ToE_closure",
        ],
    }

    artifact = {
        "stage": "F318",
        "lane": "strict_sigma_int_upgraded_seed_candidate_instance_only",
        "goal": "export_one_strict_sigma_int_upgraded_strict_side_seed_candidate_instance_for_future_S_sel_int_work",
        "status": "F318_EXECUTED_FIRST_CURRENT_STRICT_SIDE_SOURCE_SEED_STRICT_SIGMA_INT_UPGRADE_CANDIDATE_INSTANCE_PACKET_NO_FALSE_PASS",
        "candidate_seed_instance": candidate_seed_instance,
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "lane": artifact["lane"],
        "candidate_seed_instance": artifact["candidate_seed_instance"],
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(
        json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii"
    )
    OUT_SUMMARY.write_text(
        json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii"
    )
    print(OUT_JSON)


if __name__ == "__main__":
    main()

