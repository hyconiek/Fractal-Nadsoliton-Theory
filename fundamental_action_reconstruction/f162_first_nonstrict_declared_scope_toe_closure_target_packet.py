#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f162_first_nonstrict_declared_scope_toe_closure_target_packet_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    n271 = load_json(
        GENERATED / "n271_current_first_nonstrict_declared_scope_toe_preclosure_support_theorem_summary.json"
    )

    theorem = n271["theorem_result"]

    summary = {
        "packet_id": "F162",
        "status": "F162_EXECUTED_FIRST_NONSTRICT_DECLARED_SCOPE_TOE_CLOSURE_TARGET_PACKET_FUTURE_ROUTE_ONLY_NO_FALSE_PASS",
        "as_of": "2026-03-09",
        "support_packet_name": "Gamma_nonstrict_declared_scope_toe_closure_target_support_v1",
        "future_target_name": "C_toe_declared_scope_nonstrict_future_target_v1",
        "future_target_codomain": "axiom_augmented_declared_scope_toe_closure_target_v1",
        "nonstrict_declared_scope_toe_preclosure_support_packet_exported": theorem[
            "nonstrict_declared_scope_toe_preclosure_support_packet_exported"
        ],
        "accepted_scope": theorem["accepted_scope"],
        "strict_core_changed": theorem["strict_core_changed"],
        "observer_role": theorem["observer_role"],
        "future_theorem_route_spec_present": True,
        "future_only_target_exported": True,
        "actual_nonstrict_declared_scope_toe_closure_discharged": False,
        "actual_strict_core_toe_closure_discharged": False,
        "actual_global_toe_closure_discharged": False,
        "toe_closure_claimed": False,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
