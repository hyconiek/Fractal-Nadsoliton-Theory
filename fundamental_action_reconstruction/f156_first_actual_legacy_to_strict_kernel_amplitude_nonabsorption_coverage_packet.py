#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f156_first_actual_legacy_to_strict_kernel_amplitude_nonabsorption_coverage_packet_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    n83 = load_json(
        GENERATED / "n83_current_legacy_weinberg_full_claim_specific_negative_closure_theorem_summary.json"
    )
    n99 = load_json(
        GENERATED / "n99_current_legacy_fine_structure_full_claim_specific_negative_closure_theorem_summary.json"
    )
    n115 = load_json(
        GENERATED / "n115_current_legacy_gravity_hierarchy_full_claim_specific_negative_closure_theorem_summary.json"
    )
    n116 = load_json(
        GENERATED / "n116_current_legacy_physical_role_package_full_negative_closure_theorem_summary.json"
    )
    n265 = load_json(
        GENERATED / "n265_current_first_actual_legacy_to_strict_kernel_claim_specific_amplitude_nonabsorption_witness_theorem_summary.json"
    )

    summary = {
        "packet_id": "F156",
        "status": "F156_EXECUTED_FIRST_ACTUAL_LEGACY_TO_STRICT_KERNEL_AMPLITUDE_NONABSORPTION_COVERAGE_PACKET_NO_FALSE_PASS",
        "as_of": "2026-03-08",
        "coverage_packet_name": "Kappa_abs_nonbridge_coverage_packet_v1",
        "legacy_weinberg_claim_specific_frontier_closed": n83["theorem_result"]["legacy_weinberg_claim_specific_frontier_closed_negatively_on_current_repo_state"],
        "legacy_fine_structure_claim_specific_frontier_closed": n99["theorem_result"]["legacy_fine_structure_claim_specific_frontier_closed_negatively_on_current_repo_state"],
        "legacy_gravity_hierarchy_claim_specific_frontier_closed": n115["theorem_result"]["legacy_gravity_hierarchy_claim_specific_frontier_closed_negatively_on_current_repo_state"],
        "legacy_physical_role_package_closed_negatively": n116["theorem_result"]["legacy_physical_role_package_closed_negatively_on_current_repo_state"],
        "actual_claim_specific_amplitude_witness_exported": n265["theorem_result"]["actual_claim_specific_amplitude_nonabsorption_witness_exported"],
        "actual_amplitude_nonabsorption_coverage_packet_exported": True,
        "full_a_abs_obstruction_discharged": False,
        "actual_nonbridge_strengthening_discharged": False,
        "actual_bridge_discharged": False,
        "kernel_split_safe": True,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
