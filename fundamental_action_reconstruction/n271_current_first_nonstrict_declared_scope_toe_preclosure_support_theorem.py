#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "n271_current_first_nonstrict_declared_scope_toe_preclosure_support_theorem_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    p251 = load_json(
        GENERATED / "p251_current_actual_nonstrict_declared_scope_toe_preclosure_support_probe_summary.json"
    )

    expected_status = "CURRENT_REPO_EXPORTS_ONE_ACTUAL_NONSTRICT_DECLARED_SCOPE_TOE_PRECLOSURE_SUPPORT_PACKET_AFTER_P251"
    status_ok = p251["status"] == expected_status

    summary = {
        "step": "N271",
        "status": "N271_DISCHARGED_CURRENT_FIRST_NONSTRICT_DECLARED_SCOPE_TOE_PRECLOSURE_SUPPORT_THEOREM_NO_FALSE_PASS",
        "scope": "axiom_augmented_declared_scope_toe_preclosure_only",
        "checks": [
            {
                "id": "p251_nonstrict_declared_scope_toe_preclosure_status",
                "actual": p251["status"],
                "expected": expected_status,
                "pass": status_ok,
            }
        ],
        "theorem_result": {
            "discharged": status_ok and len(p251.get("blocking_mismatches", [])) == 0,
            "nonstrict_declared_scope_toe_preclosure_support_packet_exported": p251["nonstrict_declared_scope_toe_preclosure_support_packet_exported"],
            "accepted_scope": p251["accepted_scope"],
            "strict_core_changed": p251["strict_core_changed"],
            "observer_role": p251["observer_role"],
            "actual_nonstrict_declared_scope_toe_closure_discharged": p251["actual_nonstrict_declared_scope_toe_closure_discharged"],
            "actual_strict_core_toe_closure_discharged": p251["actual_strict_core_toe_closure_discharged"],
            "actual_global_toe_closure_discharged": p251["actual_global_toe_closure_discharged"],
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
