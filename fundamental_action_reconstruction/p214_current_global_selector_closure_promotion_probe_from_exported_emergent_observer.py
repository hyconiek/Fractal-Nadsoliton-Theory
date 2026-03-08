#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "p214_current_global_selector_closure_promotion_probe_from_exported_emergent_observer_summary.json"
)


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def check_map(checks: list[dict]) -> dict:
    return {item["id"]: item for item in checks}


def main() -> None:
    p213 = load_json(
        GENERATED / "p213_current_next_actual_emergent_observer_closure_readout_operator_probe_summary.json"
    )
    n118 = load_json(
        GENERATED
        / "n118_current_selector_or_symmetry_breaking_requirement_theorem_for_qw2191_summary.json"
    )
    n163 = load_json(
        GENERATED
        / "n163_current_observer_information_deficit_downstream_symptom_theorem_summary.json"
    )
    n233 = load_json(
        GENERATED
        / "n233_current_next_actual_emergent_observer_closure_readout_operator_theorem_summary.json"
    )

    p213_checks = check_map(p213["checks"])

    checks_spec = [
        {
            "id": "admissible_local_positive_composite_exported",
            "actual": p213["actual_closure_readout_exported"],
            "expected": True,
        },
        {
            "id": "local_asymmetry_stable_commit_positive_and_gap_zero",
            "actual": p213_checks["positive_commit_amplitude"]["pass"]
            and p213_checks["zero_gap_channel"]["pass"],
            "expected": True,
        },
        {
            "id": "theorem_level_actual_emergent_observer_closure_exported",
            "actual": False,
            "expected": True,
        },
        {
            "id": "observer_not_downstream_only",
            "actual": not n163["theorem_result"][
                "observer_information_deficit_downstream_symptom_on_current_repo_state"
            ],
            "expected": True,
        },
        {
            "id": "theorem_level_strict_core_selector_closure_exported",
            "actual": p213["strict_core_selector_closure"],
            "expected": True,
        },
        {
            "id": "basis_independent_global_promotion_witness_exported",
            "actual": False,
            "expected": True,
        },
        {
            "id": "theorem_level_global_qw2191_discharge_exported",
            "actual": "no_QW2191_discharge" not in n118["hard_limits"]
            and "no_QW2191_discharge" not in n233["hard_limits"],
            "expected": True,
        },
        {
            "id": "kernel_split_safe_promotion_semantics",
            "actual": p213_checks["kernel_split_safe"]["pass"],
            "expected": True,
        },
    ]

    checks = []
    failed = []
    for item in checks_spec:
        ok = item["actual"] == item["expected"]
        checks.append(
            {
                "id": item["id"],
                "actual": item["actual"],
                "expected": item["expected"],
                "pass": ok,
            }
        )
        if not ok:
            failed.append(item["id"])

    if failed:
        status = "CURRENT_REPO_DOES_NOT_JUSTIFY_GLOBAL_PROMOTION_OF_THE_EXPORTED_EMERGENT_OBSERVER_COMPOSITE_TO_SELECTOR_CLOSURE_OR_QW2191_DISCHARGE_AFTER_P214"
    else:
        status = "CURRENT_REPO_JUSTIFIES_GLOBAL_PROMOTION_OF_THE_EXPORTED_EMERGENT_OBSERVER_COMPOSITE_TO_SELECTOR_CLOSURE_AND_QW2191_DISCHARGE_AFTER_P214"

    summary = {
        "stage": "P214",
        "lane": "current_global_promotion_from_exported_emergent_observer_only",
        "status": status,
        "checks": checks,
        "failed_global_promotion_requirements": failed,
        "local_positive_composite_exists": p213["actual_closure_readout_exported"],
        "local_asymmetry_stable": p213_checks["positive_commit_amplitude"]["pass"]
        and p213_checks["zero_gap_channel"]["pass"],
        "observer_downstream_only": n163["theorem_result"][
            "observer_information_deficit_downstream_symptom_on_current_repo_state"
        ],
        "global_selector_closure_promotable_on_current_repo_state": not failed,
        "global_qw2191_discharge_promotable_on_current_repo_state": not failed,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
