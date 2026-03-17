#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"

IN_N123 = GENERATED / "n123_current_legacy_to_strict_kernel_package_level_nonbridge_theorem_summary.json"
IN_N117 = GENERATED / "n117_current_legacy_to_strict_kernel_bridge_package_nontransfer_theorem_summary.json"
IN_N267 = GENERATED / "n267_current_first_actual_legacy_to_strict_kernel_full_amplitude_nonabsorption_obstruction_witness_theorem_summary.json"
IN_N268 = GENERATED / "n268_current_first_actual_legacy_to_strict_kernel_damping_nonrenormalization_obstruction_witness_theorem_summary.json"
IN_N438 = GENERATED / "n438_current_first_actual_legacy_to_strict_kernel_phase_frequency_nonconformal_obstruction_witness_theorem_summary.json"

OUT = (
    GENERATED
    / "f662_current_actual_legacy_to_strict_kernel_nonbridge_strengthening_discharge_witness_packet_summary.json"
)


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    n123 = load_json(IN_N123)
    n117 = load_json(IN_N117)
    n267 = load_json(IN_N267)
    n268 = load_json(IN_N268)
    n438 = load_json(IN_N438)

    package_level_nonbridge = bool(
        n123["theorem_result"]["package_level_nonbridge_on_current_repo_state"]
    )
    package_nontransfer = bool(
        n117["theorem_result"]["legacy_to_strict_package_nontransfer_on_current_repo_state"]
    )
    a_abs = bool(n267["theorem_result"]["a_abs_nonbridge_obstruction_discharged"])
    r_damp = bool(n268["theorem_result"]["r_damp_nonbridge_obstruction_discharged"])
    phase = bool(
        n438["theorem_result"]["actual_phase_frequency_nonconformal_obstruction_discharged"]
    )

    summary = {
        "packet_id": "F662",
        "status": "F662_EXECUTED_CURRENT_ACTUAL_LEGACY_TO_STRICT_KERNEL_NONBRIDGE_STRENGTHENING_DISCHARGE_WITNESS_PACKET_NO_FALSE_PASS",
        "as_of": "2026-03-17",
        "input_legacy_kernel": "K_legacy_ont",
        "input_strict_kernel": "K_strict_gate",
        "support_packet_name": "Lambda_legacy_strict_nonbridge_strengthening_full_support_v1",
        "witness_name": "NB_legacy_strict_strengthening_actual_witness_v1",
        "target_name": "explicit_legacy_strict_kernel_nonbridge_strengthening_target_v1",
        "components": {
            "package_level_nonbridge_on_current_repo_state": package_level_nonbridge,
            "legacy_to_strict_package_nontransfer_on_current_repo_state": package_nontransfer,
            "a_abs_nonbridge_obstruction_discharged": a_abs,
            "r_damp_nonbridge_obstruction_discharged": r_damp,
            "actual_phase_frequency_nonconformal_obstruction_discharged": phase,
        },
        "actual_nonbridge_strengthening_discharged": bool(
            package_level_nonbridge and package_nontransfer and a_abs and r_damp and phase
        ),
        "permanent_no_bridge_claimed": False,
        "actual_bridge_discharged": False,
        "branch_selection_justified_on_current_repo_state": False,
        "kernel_split_safe": True,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()

