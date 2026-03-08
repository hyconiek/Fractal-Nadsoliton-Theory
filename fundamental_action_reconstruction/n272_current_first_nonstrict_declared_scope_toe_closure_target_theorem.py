#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "n272_current_first_nonstrict_declared_scope_toe_closure_target_theorem_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    p252 = load_json(
        GENERATED / "p252_current_nonstrict_declared_scope_toe_closure_target_probe_summary.json"
    )

    expected_status = (
        "CURRENT_REPO_EXPORTS_ONE_FUTURE_ONLY_NONSTRICT_DECLARED_SCOPE_TOE_CLOSURE_TARGET_BELOW_ACTUAL_TOE_CLOSURE_AFTER_P252"
    )
    status_ok = p252["status"] == expected_status

    summary = {
        "step": "N272",
        "status": "N272_DISCHARGED_CURRENT_FIRST_NONSTRICT_DECLARED_SCOPE_TOE_CLOSURE_TARGET_THEOREM_NO_FALSE_PASS",
        "scope": "axiom_augmented_declared_scope_toe_future_target_only",
        "checks": [
            {
                "id": "p252_nonstrict_declared_scope_toe_target_status",
                "actual": p252["status"],
                "expected": expected_status,
                "pass": status_ok,
            }
        ],
        "theorem_result": {
            "discharged": status_ok and len(p252.get("blocking_mismatches", [])) == 0,
            "future_only_target_exported": p252["future_only_target_exported"],
            "accepted_scope": p252["accepted_scope"],
            "strict_core_changed": p252["strict_core_changed"],
            "observer_role": p252["observer_role"],
            "actual_nonstrict_declared_scope_toe_closure_discharged": p252["actual_nonstrict_declared_scope_toe_closure_discharged"],
            "actual_strict_core_toe_closure_discharged": p252["actual_strict_core_toe_closure_discharged"],
            "actual_global_toe_closure_discharged": p252["actual_global_toe_closure_discharged"],
            "toe_closure_claimed": p252["toe_closure_claimed"],
        },
        "hard_limits": [
            "no_actual_nonstrict_declared_scope_toe_closure",
            "no_actual_strict_core_toe_closure",
            "no_actual_global_toe_closure",
        ],
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
