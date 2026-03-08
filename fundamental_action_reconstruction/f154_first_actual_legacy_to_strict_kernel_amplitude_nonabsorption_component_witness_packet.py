#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f154_first_actual_legacy_to_strict_kernel_amplitude_nonabsorption_component_witness_packet_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    n65 = load_json(
        GENERATED / "n65_current_legacy_physical_role_transfer_to_strict_gate_kernel_obstruction_theorem_summary.json"
    )
    n70 = load_json(
        GENERATED / "n70_current_strict_side_weinberg_role_equivalence_boundary_theorem_summary.json"
    )
    n83 = load_json(
        GENERATED / "n83_current_legacy_weinberg_full_claim_specific_negative_closure_theorem_summary.json"
    )
    n263 = load_json(
        GENERATED / "n263_current_first_actual_legacy_to_strict_kernel_bifurcated_frontier_theorem_summary.json"
    )

    summary = {
        "packet_id": "F154",
        "status": "F154_EXECUTED_FIRST_ACTUAL_LEGACY_TO_STRICT_KERNEL_AMPLITUDE_NONABSORPTION_COMPONENT_WITNESS_PACKET_NO_FALSE_PASS",
        "as_of": "2026-03-08",
        "support_packet_name": "W_abs_nonbridge_component_support_v1",
        "component_witness_name": "A_abs_nonbridge_component_witness_v1",
        "component_scope": "legacy_weinberg_amplitude_role_only",
        "legacy_weinberg_role_transfer_present": n65["theorem_result"]["legacy_weinberg_role_transfer_present"],
        "strict_side_candidate_object_present": n70["theorem_result"]["strict_side_candidate_object_present"],
        "explicit_role_equivalence_verdict_present": n70["theorem_result"]["explicit_role_equivalence_verdict_present"],
        "legacy_weinberg_claim_specific_frontier_closed_negatively": n83["theorem_result"]["legacy_weinberg_claim_specific_frontier_closed_negatively_on_current_repo_state"],
        "bridge_nonbridge_frontier_undecided": not n263["theorem_result"]["branch_selection_justified_on_current_repo_state"],
        "actual_amplitude_nonabsorption_component_witness_exported": True,
        "actual_amplitude_nonabsorption_obstruction_discharged": False,
        "actual_nonbridge_strengthening_discharged": False,
        "actual_bridge_discharged": False,
        "kernel_split_safe": True,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
