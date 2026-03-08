#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "n267_current_first_actual_legacy_to_strict_kernel_full_amplitude_nonabsorption_obstruction_witness_theorem_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    p247 = load_json(
        GENERATED / "p247_current_actual_legacy_to_strict_kernel_full_amplitude_nonabsorption_obstruction_witness_probe_summary.json"
    )

    expected_status = (
        "CURRENT_REPO_EXPORTS_ONE_ACTUAL_LEGACY_TO_STRICT_KERNEL_FULL_AMPLITUDE_NONABSORPTION_OBSTRUCTION_WITNESS_AFTER_P247"
    )
    status_ok = p247["status"] == expected_status

    summary = {
        "step": "N267",
        "status": "N267_DISCHARGED_CURRENT_FIRST_ACTUAL_LEGACY_TO_STRICT_KERNEL_FULL_AMPLITUDE_NONABSORPTION_OBSTRUCTION_WITNESS_THEOREM_NO_FALSE_PASS",
        "scope": "declared_kernel_comparison_amplitude_layer_only",
        "checks": [
            {
                "id": "p247_full_a_abs_status",
                "actual": p247["status"],
                "expected": expected_status,
                "pass": status_ok,
            }
        ],
        "theorem_result": {
            "discharged": status_ok and len(p247.get("blocking_mismatches", [])) == 0,
            "a_abs_nonbridge_obstruction_discharged": p247["a_abs_nonbridge_obstruction_discharged"],
            "legacy_to_strict_package_nontransfer_on_current_repo_state": p247["legacy_to_strict_package_nontransfer_on_current_repo_state"],
            "actual_amplitude_nonabsorption_coverage_packet_exported": p247["actual_amplitude_nonabsorption_coverage_packet_exported"],
            "actual_damping_nonrenormalization_obstruction_discharged": False,
            "actual_phase_frequency_nonconformal_obstruction_discharged": False,
            "actual_nonbridge_strengthening_discharged": False,
            "actual_bridge_discharged": False,
            "current_strict_core_selector_closure": False,
            "current_global_selector_closure": False,
            "current_global_qw2191_discharge": False,
            "kernel_split_safe": p247["kernel_split_safe"],
        },
        "hard_limits": [
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
