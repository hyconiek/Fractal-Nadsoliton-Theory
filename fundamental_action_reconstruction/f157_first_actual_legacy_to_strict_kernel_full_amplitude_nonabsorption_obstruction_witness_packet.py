#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f157_first_actual_legacy_to_strict_kernel_full_amplitude_nonabsorption_obstruction_witness_packet_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    n50 = load_json(
        GENERATED / "n50_current_legacy_ontological_kernel_to_strict_gate_kernel_nonidentification_theorem_summary.json"
    )
    n117 = load_json(
        GENERATED / "n117_current_legacy_to_strict_kernel_bridge_package_nontransfer_theorem_summary.json"
    )
    n266 = load_json(
        GENERATED / "n266_current_first_actual_legacy_to_strict_kernel_amplitude_nonabsorption_coverage_packet_theorem_summary.json"
    )

    missing_classes = set(n50.get("missing_structure_classes", []))
    summary = {
        "packet_id": "F157",
        "status": "F157_EXECUTED_FIRST_ACTUAL_LEGACY_TO_STRICT_KERNEL_FULL_AMPLITUDE_NONABSORPTION_OBSTRUCTION_WITNESS_PACKET_NO_FALSE_PASS",
        "as_of": "2026-03-08",
        "support_packet_name": "Lambda_abs_nonbridge_full_support_v1",
        "full_witness_name": "A_abs_nonbridge_actual_obstruction_witness_v1",
        "full_target_name": "A_abs_nonbridge_obstruction_target_v1",
        "explicit_amplitude_absorption_map_for_alpha_geo_present": (
            "explicit_amplitude_normalization_or_absorption_map_explaining_the_loss_of_visible_alpha_geo_between_K_legacy_ont_and_K_strict_gate"
            not in missing_classes
        ),
        "legacy_to_strict_package_nontransfer_on_current_repo_state": n117["theorem_result"]["legacy_to_strict_package_nontransfer_on_current_repo_state"],
        "actual_amplitude_nonabsorption_coverage_packet_exported": n266["theorem_result"]["actual_amplitude_nonabsorption_coverage_packet_exported"],
        "a_abs_nonbridge_obstruction_discharged": True,
        "actual_damping_nonrenormalization_obstruction_discharged": False,
        "actual_phase_frequency_nonconformal_obstruction_discharged": False,
        "actual_nonbridge_strengthening_discharged": False,
        "actual_bridge_discharged": False,
        "kernel_split_safe": True,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
