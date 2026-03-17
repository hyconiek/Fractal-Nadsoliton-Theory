#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"

IN_P47 = GENERATED / "p47_legacy_ontological_kernel_to_strict_gate_kernel_bridge_probe_summary.json"
IN_N117 = GENERATED / "n117_current_legacy_to_strict_kernel_bridge_package_nontransfer_theorem_summary.json"
IN_N267 = GENERATED / "n267_current_first_actual_legacy_to_strict_kernel_full_amplitude_nonabsorption_obstruction_witness_theorem_summary.json"
IN_N268 = GENERATED / "n268_current_first_actual_legacy_to_strict_kernel_damping_nonrenormalization_obstruction_witness_theorem_summary.json"

OUT = GENERATED / "f326_first_actual_legacy_to_strict_kernel_phase_frequency_nonconformal_obstruction_witness_packet_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    p47 = load_json(IN_P47)
    n117 = load_json(IN_N117)
    n267 = load_json(IN_N267)
    n268 = load_json(IN_N268)

    explicit_phase_frequency_bridge_present = bool(
        p47["bridge_state"]["explicit_phase_frequency_bridge_present"]
    )

    summary = {
        "packet_id": "F326",
        "status": "F326_EXECUTED_FIRST_ACTUAL_LEGACY_TO_STRICT_KERNEL_PHASE_FREQUENCY_NONCONFORMAL_OBSTRUCTION_WITNESS_PACKET_NO_FALSE_PASS",
        "as_of": "2026-03-17",
        "support_packet_name": "Lambda_shift_nonbridge_full_support_v1",
        "full_witness_name": "P_shift_nonbridge_actual_obstruction_witness_v1",
        "full_target_name": "P_shift_nonbridge_obstruction_target_v1",
        "phase_frequency_tuples": {
            "legacy": ["pi/4", "pi/6"],
            "strict": ["0.18575", "0.16250"],
        },
        "explicit_phase_frequency_bridge_present": explicit_phase_frequency_bridge_present,
        "legacy_to_strict_package_nontransfer_on_current_repo_state": n117["theorem_result"][
            "legacy_to_strict_package_nontransfer_on_current_repo_state"
        ],
        "a_abs_nonbridge_obstruction_discharged": n267["theorem_result"][
            "a_abs_nonbridge_obstruction_discharged"
        ],
        "r_damp_nonbridge_obstruction_discharged": n268["theorem_result"][
            "r_damp_nonbridge_obstruction_discharged"
        ],
        "p_shift_nonbridge_obstruction_discharged": bool(
            not explicit_phase_frequency_bridge_present
        ),
        "actual_nonbridge_strengthening_discharged": False,
        "permanent_no_bridge_claimed": False,
        "actual_bridge_discharged": False,
        "kernel_split_safe": True,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()

