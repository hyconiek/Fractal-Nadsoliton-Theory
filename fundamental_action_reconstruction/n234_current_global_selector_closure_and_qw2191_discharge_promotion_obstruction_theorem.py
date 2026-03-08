#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "n234_current_global_selector_closure_and_qw2191_discharge_promotion_obstruction_theorem_summary.json"
)


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    p214 = load_json(
        GENERATED
        / "p214_current_global_selector_closure_promotion_probe_from_exported_emergent_observer_summary.json"
    )

    expected_status = (
        "CURRENT_REPO_DOES_NOT_JUSTIFY_GLOBAL_PROMOTION_OF_THE_EXPORTED_EMERGENT_OBSERVER_COMPOSITE_TO_SELECTOR_CLOSURE_OR_QW2191_DISCHARGE_AFTER_P214"
    )
    status_ok = p214["status"] == expected_status

    summary = {
        "step": "N234",
        "status": "N234_DISCHARGED_CURRENT_GLOBAL_SELECTOR_CLOSURE_AND_QW2191_DISCHARGE_PROMOTION_OBSTRUCTION_THEOREM_NO_FALSE_PASS",
        "scope": "current_global_selector_closure_and_qw2191_discharge_promotion_question_only",
        "checks": [
            {
                "id": "p214_obstruction_status",
                "actual": p214["status"],
                "expected": expected_status,
                "pass": status_ok,
            }
        ],
        "theorem_result": {
            "discharged": status_ok,
            "local_exported_emergent_observer_composite_exists": p214[
                "local_positive_composite_exists"
            ],
            "local_asymmetry_stable": p214["local_asymmetry_stable"],
            "observer_downstream_only": p214["observer_downstream_only"],
            "global_selector_closure_justified_on_current_repo_state": False,
            "global_qw2191_discharge_justified_on_current_repo_state": False,
            "full_closure_pass": False,
        },
        "hard_limits": [
            "actual_emergent_observer_closure_not_yet_constructed",
            "observer_chain_remains_downstream_only",
            "no_basis_independent_global_selector_promotion_witness",
            "no_strict_core_selector_closure",
            "no_QW2191_discharge",
            "no_ToE_closure",
        ],
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
