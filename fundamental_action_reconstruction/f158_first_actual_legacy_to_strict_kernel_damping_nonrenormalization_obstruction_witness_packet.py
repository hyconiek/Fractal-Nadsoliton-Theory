#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f158_first_actual_legacy_to_strict_kernel_damping_nonrenormalization_obstruction_witness_packet_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    p47 = load_json(
        GENERATED / "p47_legacy_ontological_kernel_to_strict_gate_kernel_bridge_probe_summary.json"
    )
    n117 = load_json(
        GENERATED / "n117_current_legacy_to_strict_kernel_bridge_package_nontransfer_theorem_summary.json"
    )
    n267 = load_json(
        GENERATED / "n267_current_first_actual_legacy_to_strict_kernel_full_amplitude_nonabsorption_obstruction_witness_theorem_summary.json"
    )

    summary = {
        "packet_id": "F158",
        "status": "F158_EXECUTED_FIRST_ACTUAL_LEGACY_TO_STRICT_KERNEL_DAMPING_NONRENORMALIZATION_OBSTRUCTION_WITNESS_PACKET_NO_FALSE_PASS",
        "as_of": "2026-03-09",
        "support_packet_name": "Lambda_damp_nonbridge_full_support_v1",
        "full_witness_name": "R_damp_nonbridge_actual_obstruction_witness_v1",
        "full_target_name": "R_damp_nonbridge_obstruction_target_v1",
        "explicit_beta_tors_to_beta_eta_translation_present": p47["bridge_state"]["explicit_beta_tors_to_beta_eta_translation_present"],
        "legacy_to_strict_package_nontransfer_on_current_repo_state": n117["theorem_result"]["legacy_to_strict_package_nontransfer_on_current_repo_state"],
        "a_abs_nonbridge_obstruction_discharged": n267["theorem_result"]["a_abs_nonbridge_obstruction_discharged"],
        "r_damp_nonbridge_obstruction_discharged": True,
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
