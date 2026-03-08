#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "n259_current_declared_scope_source_topology_selector_theorem_promotion_obstruction_theorem_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    p239 = load_json(
        GENERATED / "p239_current_declared_scope_source_topology_selector_theorem_promotion_probe_summary.json"
    )

    expected_status = (
        "CURRENT_REPO_DOES_NOT_JUSTIFY_PROMOTION_OF_THE_DECLARED_SCOPE_SOURCE_TOPOLOGY_SELECTOR_THEOREM_TO_STRICT_CORE_SELECTOR_CLOSURE_OR_GLOBAL_QW2191_DISCHARGE_AFTER_P239"
    )
    status_ok = p239["status"] == expected_status

    summary = {
        "step": "N259",
        "status": "N259_DISCHARGED_CURRENT_DECLARED_SCOPE_SOURCE_TOPOLOGY_SELECTOR_THEOREM_PROMOTION_OBSTRUCTION_THEOREM_NO_FALSE_PASS",
        "scope": "current_declared_scope_source_topology_selector_theorem_promotion_question_only",
        "checks": [
            {
                "id": "p239_promotion_obstruction_status",
                "actual": p239["status"],
                "expected": expected_status,
                "pass": status_ok,
            }
        ],
        "theorem_result": {
            "discharged": status_ok,
            "declared_scope_theorem_exported": p239["declared_scope_theorem_exported"],
            "declared_scope_only": p239["declared_scope_only"],
            "admissible_strict_core_internal_selector_source_object_present": p239["admissible_strict_core_internal_selector_source_object_present"],
            "selector_requirement_still_active_at_closure_frontier": p239["selector_requirement_still_active_at_closure_frontier"],
            "observer_downstream_only": p239["observer_downstream_only"],
            "strict_core_selector_closure_justified_on_current_repo_state": False,
            "global_selector_closure_justified_on_current_repo_state": False,
            "global_qw2191_discharge_justified_on_current_repo_state": False,
            "full_closure_pass": False,
        },
        "hard_limits": [
            "declared_scope_theorem_is_real_but_not_yet_closure_level",
            "no_future_stronger_source_topology_result_is_ruled_out_forever",
            "no_admissible_strict_core_internal_selector_source_object_is_currently_exported",
            "no_strict_core_selector_closure",
            "no_global_selector_closure",
            "no_global_QW2191_discharge",
            "no_ToE_closure",
        ],
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
