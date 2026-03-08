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
    / "p118_initial_source_admission_and_orientation_export_contract_probe.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p118_initial_source_admission_and_orientation_export_contract_probe_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    f29 = load_json(
        "fundamental_action_reconstruction/generated/f29_genuine_strict_core_internal_selector_source_admission_packet_summary.json"
    )
    f31 = load_json(
        "fundamental_action_reconstruction/generated/f31_initial_future_strict_core_selector_source_seed_target_packet_summary.json"
    )
    f32 = load_json(
        "fundamental_action_reconstruction/generated/f32_initial_future_strict_core_orientation_export_admission_packet_summary.json"
    )

    reduced_to_one_initial_package = (
        f29["admission_contract"]["strict_core_source_export_required"]
        and f31["initial_seed_target"]["strict_core_source_object"] == "S_sel_int"
        and f31["initial_seed_target"]["internal_orientation_export"] == "E_orient"
        and f32["orientation_export_name"] == "E_orient"
        and f32["upstream_source_name"] == "S_sel_int"
    )

    checks_spec = [
        {
            "id": "f29_source_admission_present",
            "actual": f29["admission_contract"]["strict_core_source_export_required"],
            "expected": True,
            "meaning": "the future source admission package is present",
        },
        {
            "id": "f31_seed_source_name",
            "actual": f31["initial_seed_target"]["strict_core_source_object"],
            "expected": "S_sel_int",
            "meaning": "the frozen seed package names the future source object",
        },
        {
            "id": "f31_seed_orientation_name",
            "actual": f31["initial_seed_target"]["internal_orientation_export"],
            "expected": "E_orient",
            "meaning": "the frozen seed package names the future orientation export",
        },
        {
            "id": "f32_orientation_contract_bound_to_seed",
            "actual": [
                f32["upstream_source_name"],
                f32["orientation_export_name"],
            ],
            "expected": ["S_sel_int", "E_orient"],
            "meaning": "the orientation-export contract is bound to the current seed pair",
        },
        {
            "id": "reduced_to_one_initial_package",
            "actual": reduced_to_one_initial_package,
            "expected": True,
            "meaning": "the last positive branch is reduced to one explicit initial source-plus-orientation package",
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
        "stage": "P118",
        "lane": "initial_source_admission_and_orientation_export_contract_only",
        "goal": "test_whether_the_last_positive_branch_is_now_reduced_to_one_initial_source_plus_orientation_contract_package",
        "status": "CURRENT_REPO_REDUCES_THE_LAST_POSITIVE_BRANCH_TO_ONE_INITIAL_SOURCE_ADMISSION_AND_ORIENTATION_EXPORT_CONTRACT_PACKAGE_AFTER_P118",
        "target_state": {
            "last_positive_branch_reduced_to_one_initial_package": reduced_to_one_initial_package,
            "initial_package": {
                "strict_core_source_object": "S_sel_int",
                "source_admission_contract": f29["admission_contract"],
                "internal_orientation_export": "E_orient",
                "orientation_export_contract": f32["orientation_export_contract"],
            },
            "downstream_chain_left_open": f31["downstream_chain_left_open"],
        },
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P118",
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
