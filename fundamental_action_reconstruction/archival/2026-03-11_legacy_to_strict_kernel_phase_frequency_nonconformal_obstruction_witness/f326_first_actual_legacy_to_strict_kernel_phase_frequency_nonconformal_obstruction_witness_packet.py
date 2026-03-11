#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[2]
ACTIVE_GENERATED = REPO / "fundamental_action_reconstruction" / "generated"
GENERATED = ROOT / "generated"

IN_P47 = ACTIVE_GENERATED / "p47_legacy_ontological_kernel_to_strict_gate_kernel_bridge_probe_summary.json"
IN_N117 = ACTIVE_GENERATED / "n117_current_legacy_to_strict_kernel_bridge_package_nontransfer_theorem_summary.json"
IN_N267 = ACTIVE_GENERATED / "n267_current_first_actual_legacy_to_strict_kernel_full_amplitude_nonabsorption_obstruction_witness_theorem_summary.json"
IN_N268 = ACTIVE_GENERATED / "n268_current_first_actual_legacy_to_strict_kernel_damping_nonrenormalization_obstruction_witness_theorem_summary.json"

OUT_JSON = GENERATED / "f326_phase_frequency_nonconformal_obstruction_witness_packet.json"
OUT_SUMMARY = GENERATED / "f326_phase_frequency_nonconformal_obstruction_witness_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    p47 = load_json(IN_P47)
    n117 = load_json(IN_N117)
    n267 = load_json(IN_N267)
    n268 = load_json(IN_N268)

    explicit_phase_frequency_bridge_present = bool(
        p47["bridge_state"]["explicit_phase_frequency_bridge_present"]
    )
    legacy_to_strict_package_nontransfer = bool(
        n117["theorem_result"]["legacy_to_strict_package_nontransfer_on_current_repo_state"]
    )
    a_abs_discharged = bool(n267["theorem_result"]["a_abs_nonbridge_obstruction_discharged"])
    r_damp_discharged = bool(
        n268["theorem_result"]["r_damp_nonbridge_obstruction_discharged"]
    )

    p_shift_obstruction_on_current_exports = bool(not explicit_phase_frequency_bridge_present)

    packet = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "stage": "F326_ARCHIVAL",
        "archival": True,
        "status": "ARCHIVAL_SCRATCH_DRAFT_DO_NOT_CITE_AS_ACTIVE_EXECUTION_LANE",
        "goal": "Preserve a draft phase/frequency nonconformal obstruction witness (P_shift) logic for legacy→strict kernel comparison, frozen by S2.",
        "inputs": {
            "p47_bridge_probe_summary": str(IN_P47.relative_to(REPO)),
            "n117_package_nontransfer_summary": str(IN_N117.relative_to(REPO)),
            "n267_amplitude_obstruction_summary": str(IN_N267.relative_to(REPO)),
            "n268_damping_obstruction_summary": str(IN_N268.relative_to(REPO)),
        },
        "computed": {
            "explicit_phase_frequency_bridge_present": explicit_phase_frequency_bridge_present,
            "legacy_to_strict_package_nontransfer_on_current_repo_state": legacy_to_strict_package_nontransfer,
            "a_abs_nonbridge_obstruction_discharged": a_abs_discharged,
            "r_damp_nonbridge_obstruction_discharged": r_damp_discharged,
        },
        "draft_object": {
            "name": "P_shift_nonbridge_obstruction_witness_draft_v1",
            "typed_shape": "(K_legacy_ont, K_strict_gate) -> P_shift_nonbridge_obstruction_target_v1",
            "meaning": "On the current export set, no explicit phase/frequency bridge is exported between the legacy and strict kernel tuples; therefore the phase/frequency layer is obstructed on that export set.",
        },
        "verdict": {
            "phase_frequency_layer_obstructed_on_current_exports": p_shift_obstruction_on_current_exports,
            "note": "This is an archival scratch computation and does not constitute an active-lane discharge.",
        },
        "hard_limits": [
            "ARCHIVAL ONLY: route frozen by S2 (legacy kernel retirement decree).",
            "Does not revive legacy→strict bridge/non-bridge frontier as active.",
            "Does not claim strict-core selector closure or QW-2191 discharge.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": packet["stage"],
        "archival": True,
        "phase_frequency_layer_obstructed_on_current_exports": bool(
            packet["verdict"]["phase_frequency_layer_obstructed_on_current_exports"]
        ),
        "explicit_phase_frequency_bridge_present": explicit_phase_frequency_bridge_present,
        "strict_core_promotion": False,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(
        json.dumps(packet, indent=2, ensure_ascii=True) + "\n", encoding="ascii"
    )
    OUT_SUMMARY.write_text(
        json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii"
    )
    print(json.dumps(summary, ensure_ascii=True))


if __name__ == "__main__":
    main()

