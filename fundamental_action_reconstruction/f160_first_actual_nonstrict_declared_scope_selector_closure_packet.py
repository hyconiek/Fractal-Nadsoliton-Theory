#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f160_first_actual_nonstrict_declared_scope_selector_closure_packet_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    n125 = load_json(
        GENERATED / "n125_current_selector_requirement_theory_level_acceptance_theorem_summary.json"
    )
    n258 = load_json(
        GENERATED / "n258_current_first_declared_scope_source_topology_selector_theorem_summary.json"
    )
    n269 = load_json(
        GENERATED / "n269_current_first_nadsoliton_macroscopic_identification_role_separation_theorem_summary.json"
    )

    summary = {
        "packet_id": "F160",
        "status": "F160_EXECUTED_FIRST_ACTUAL_NONSTRICT_DECLARED_SCOPE_SELECTOR_CLOSURE_PACKET_NO_FALSE_PASS",
        "as_of": "2026-03-09",
        "support_packet_name": "Lambda_nonstrict_declared_scope_selector_closure_support_v1",
        "closure_witness_name": "C_sel_declared_scope_nonstrict_actual_witness_v1",
        "closure_target_name": "axiom_augmented_declared_scope_selector_closure_target_v1",
        "declared_scope_source_topology_selector_theorem_exported": n258["declared_scope_source_topology_selector_theorem_exported"],
        "selector_requirement_accepted_at_theory_level": n125["theorem_result"]["selector_requirement_accepted_at_theory_level"],
        "accepted_scope": n125["theorem_result"]["accepted_scope"],
        "strict_core_changed": n125["theorem_result"]["strict_core_changed"],
        "bridge_nonbridge_nonmandatory_for_t14_after_n269": n269["theorem_result"]["post_n260_bridge_nonbridge_impasse_withdrawn_as_mandatory_t14_closure_gate"],
        "nonstrict_declared_scope_selector_closure_exported": True,
        "strict_core_selector_closure_claimed": False,
        "global_selector_closure_claimed": False,
        "global_qw2191_discharge_claimed": False,
        "toe_closure_claimed": False,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
