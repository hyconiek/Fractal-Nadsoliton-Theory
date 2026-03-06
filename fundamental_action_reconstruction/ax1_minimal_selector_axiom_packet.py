#!/usr/bin/env python3
import json
from pathlib import Path

root = Path(__file__).resolve().parent
out = root / "generated" / "ax1_minimal_selector_axiom_packet_summary.json"
out.parent.mkdir(parents=True, exist_ok=True)

payload = {
    "step": "AX1",
    "status": "AX1_PACKET_READY_MINIMAL_SELECTOR_AXIOM_LANE_NO_FALSE_PASS",
    "date": "2026-03-06",
    "result_kind": "axiom_augmented_positive_lane_opened",
    "axiom": "minimum_harmonic_alignment_with_orientation_convention",
    "lane": "axiom_augmented_only",
    "consequences": {
        "theta_1": "0_mod_2pi",
        "theta_2": "0_mod_2pi",
        "u_1": "c_1",
        "u_2": "c_2",
        "S_orient_axiom": "span{c_1,c_2}"
    },
    "closes_in_lane_only": [
        "actual_theta_values",
        "actual_basis_pair",
        "actual_orientation_slice_carrier"
    ],
    "does_not_close": [
        "strict_core_selector_closure",
        "T12_B1",
        "N2",
        "QW_2191",
        "axiom_free_uniqueness",
        "full_ToE_closure"
    ],
    "theorem_level_pass": False,
    "full_closure_pass": False
}

out.write_text(json.dumps(payload, indent=2) + "\n", encoding="ascii")
print(out)
