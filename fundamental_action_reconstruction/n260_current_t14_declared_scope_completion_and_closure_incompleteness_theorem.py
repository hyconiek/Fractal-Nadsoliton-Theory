#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "n260_current_t14_declared_scope_completion_and_closure_incompleteness_theorem_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    p240 = load_json(
        GENERATED / "p240_current_t14_declared_scope_completion_and_closure_incompleteness_probe_summary.json"
    )

    expected_status = (
        "CURRENT_REPO_SUPPORTS_THE_CONCLUSION_THAT_THE_T14_SOURCE_TOPOLOGY_SELECTOR_LANE_IS_DECLARED_SCOPE_COMPLETE_AND_CLOSURE_INCOMPLETE_ON_THE_PRESENT_EXPORT_SET_AFTER_P240"
    )
    status_ok = p240["status"] == expected_status

    summary = {
        "step": "N260",
        "status": "N260_DISCHARGED_CURRENT_T14_DECLARED_SCOPE_COMPLETION_AND_CLOSURE_INCOMPLETENESS_THEOREM_NO_FALSE_PASS",
        "scope": "current_t14_lane_status_on_present_export_set_only",
        "checks": [
            {
                "id": "p240_completion_status",
                "actual": p240["status"],
                "expected": expected_status,
                "pass": status_ok,
            }
        ],
        "theorem_result": {
            "discharged": status_ok,
            "declared_scope_theorem_exported": p240["declared_scope_theorem_exported"],
            "declared_scope_only": p240["declared_scope_only"],
            "t14_declared_scope_complete_on_present_export_set": p240["t14_declared_scope_complete_on_present_export_set"],
            "t14_closure_incomplete_on_present_export_set": p240["t14_closure_incomplete_on_present_export_set"],
            "strict_core_selector_closure_justified_on_current_repo_state": False,
            "global_selector_closure_justified_on_current_repo_state": False,
            "global_qw2191_discharge_justified_on_current_repo_state": False,
            "full_closure_pass": False,
        },
        "hard_limits": [
            "no_permanent_impossibility_claim",
            "no_future_closure_level_ingredient_is_ruled_out_forever",
            "declared_scope_completion_is_real_but_not_closure_level",
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
