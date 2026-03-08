#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f161_first_actual_nonstrict_declared_scope_toe_preclosure_support_packet_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    n258 = load_json(
        GENERATED / "n258_current_first_declared_scope_source_topology_selector_theorem_summary.json"
    )
    n269 = load_json(
        GENERATED / "n269_current_first_nadsoliton_macroscopic_identification_role_separation_theorem_summary.json"
    )
    n270 = load_json(
        GENERATED / "n270_current_first_nonstrict_declared_scope_selector_closure_theorem_summary.json"
    )

    summary = {
        "packet_id": "F161",
        "status": "F161_EXECUTED_FIRST_ACTUAL_NONSTRICT_DECLARED_SCOPE_TOE_PRECLOSURE_SUPPORT_PACKET_NO_FALSE_PASS",
        "as_of": "2026-03-09",
        "support_packet_name": "Lambda_nonstrict_declared_scope_toe_preclosure_support_v1",
        "nonstrict_declared_scope_selector_closure_exported": n270["theorem_result"]["nonstrict_declared_scope_selector_closure_exported"],
        "accepted_scope": n270["theorem_result"]["accepted_scope"],
        "strict_core_changed": n270["theorem_result"]["strict_core_changed"],
        "bridge_nonbridge_nonmandatory_for_t14_after_n269": n269["theorem_result"]["post_n260_bridge_nonbridge_impasse_withdrawn_as_mandatory_t14_closure_gate"],
        "declared_scope_source_topology_selector_theorem_exported": n258["declared_scope_source_topology_selector_theorem_exported"],
        "observer_role": n258["observer_role"],
        "nonstrict_declared_scope_toe_future_route_spec_present": True,
        "nonstrict_declared_scope_toe_preclosure_support_packet_exported": True,
        "actual_nonstrict_declared_scope_toe_closure_discharged": False,
        "actual_strict_core_toe_closure_discharged": False,
        "actual_global_toe_closure_discharged": False,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
