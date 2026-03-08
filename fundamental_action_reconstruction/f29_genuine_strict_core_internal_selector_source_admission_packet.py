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
    / "f29_genuine_strict_core_internal_selector_source_admission_packet.json"
)
OUT_SUMMARY = (
    GENERATED
    / "f29_genuine_strict_core_internal_selector_source_admission_packet_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    n123 = load_json(
        "fundamental_action_reconstruction/generated/n123_current_legacy_to_strict_kernel_package_level_nonbridge_theorem_summary.json"
    )
    n124 = load_json(
        "fundamental_action_reconstruction/generated/n124_current_strict_core_internal_selector_source_derivation_full_negative_closure_theorem_summary.json"
    )
    n125 = load_json(
        "fundamental_action_reconstruction/generated/n125_current_selector_requirement_theory_level_acceptance_theorem_summary.json"
    )
    p2 = load_json(
        "fundamental_action_reconstruction/generated/p2_strict_core_sigma_int_to_a1_pair1_probe_summary.json"
    )

    checks_spec = [
        {
            "id": "n124_no_current_source_discharge",
            "actual": n124["theorem_result"][
                "strict_core_internal_selector_source_derivation_discharge_present"
            ],
            "expected": False,
            "meaning": "N124 already shows that no current strict-core internal selector source discharge is exported",
        },
        {
            "id": "p2_no_downstream_operator_reachability",
            "actual": p2["status"],
            "expected": "NOT_COMPUTABLE_FROM_CURRENT_STRICT_CORE_ROUTE",
            "meaning": "P2 still shows no downstream strict-core operator reachability to A1(pair1)",
        },
        {
            "id": "n123_nonbridge_explicit",
            "actual": n123["theorem_result"][
                "package_level_nonbridge_on_current_repo_state"
            ],
            "expected": True,
            "meaning": "N123 keeps hidden legacy-to-strict substitution forbidden",
        },
        {
            "id": "n125_selector_acceptance_outside_strict_core",
            "actual": n125["theorem_result"]["accepted_scope"],
            "expected": "axiom_augmented_only",
            "meaning": "N125 keeps selector acceptance outside current strict core",
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
        "stage": "F29",
        "lane": "genuine_strict_core_internal_selector_source_admission_only",
        "goal": "define_the_minimal_admission_contract_for_any_future_genuine_strict_core_internal_selector_source_object",
        "status": "F29_EXECUTED_GENUINE_STRICT_CORE_INTERNAL_SELECTOR_SOURCE_ADMISSION_PACKET_NO_FALSE_PASS",
        "admission_contract": {
            "strict_core_source_export_required": True,
            "internal_orientation_discharge_required": True,
            "strict_core_bridge_discharge_required": True,
            "selector_reduction_discharge_required": True,
            "downstream_operator_reachability_required": True,
            "silent_legacy_to_strict_substitution_forbidden": True,
            "selector_acceptance_outside_strict_core_may_not_count_as_source_derivation": True,
        },
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "F29",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "admission_contract": artifact["admission_contract"],
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
