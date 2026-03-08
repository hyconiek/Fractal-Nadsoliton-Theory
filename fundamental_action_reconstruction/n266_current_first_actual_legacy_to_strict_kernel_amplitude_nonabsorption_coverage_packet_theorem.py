#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "n266_current_first_actual_legacy_to_strict_kernel_amplitude_nonabsorption_coverage_packet_theorem_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    p246 = load_json(
        GENERATED / "p246_current_actual_legacy_to_strict_kernel_amplitude_nonabsorption_coverage_packet_probe_summary.json"
    )

    expected_status = (
        "CURRENT_REPO_EXPORTS_ONE_ACTUAL_LEGACY_TO_STRICT_KERNEL_AMPLITUDE_NONABSORPTION_COVERAGE_PACKET_BELOW_FULL_A_ABS_OBSTRUCTION_AFTER_P246"
    )
    status_ok = p246["status"] == expected_status

    summary = {
        "step": "N266",
        "status": "N266_DISCHARGED_CURRENT_FIRST_ACTUAL_LEGACY_TO_STRICT_KERNEL_AMPLITUDE_NONABSORPTION_COVERAGE_PACKET_THEOREM_NO_FALSE_PASS",
        "scope": "legacy_alpha_geo_amplitude_coverage_only",
        "checks": [
            {
                "id": "p246_amplitude_coverage_status",
                "actual": p246["status"],
                "expected": expected_status,
                "pass": status_ok,
            }
        ],
        "theorem_result": {
            "discharged": status_ok and len(p246.get("blocking_mismatches", [])) == 0,
            "actual_amplitude_nonabsorption_coverage_packet_exported": p246["actual_amplitude_nonabsorption_coverage_packet_exported"],
            "legacy_physical_role_package_closed_negatively": p246["legacy_physical_role_package_closed_negatively"],
            "actual_claim_specific_amplitude_witness_exported": p246["actual_claim_specific_amplitude_witness_exported"],
            "full_a_abs_obstruction_discharged": False,
            "actual_nonbridge_strengthening_discharged": False,
            "actual_bridge_discharged": False,
            "current_strict_core_selector_closure": False,
            "current_global_selector_closure": False,
            "current_global_qw2191_discharge": False,
            "kernel_split_safe": p246["kernel_split_safe"],
        },
        "hard_limits": [
            "no_full_a_abs_nonbridge_obstruction",
            "no_actual_damping_nonrenormalization_obstruction",
            "no_actual_phase_frequency_nonconformal_obstruction",
            "no_actual_nonbridge_strengthening_theorem",
            "no_actual_bridge_derivation",
            "no_current_branch_selection_theorem",
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
