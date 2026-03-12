from __future__ import annotations

import json
from pathlib import Path

root = Path(__file__).resolve().parent
generated = root / "generated"

in_ax19 = generated / "strict_extension_lane_sigma_int_residual_orientation_closure_packet.json"
in_ax21 = generated / "strict_extension_lane_sigma_int_to_theta_selector_ingredient_closure_packet.json"

out_packet = generated / "strict_extension_lane_publication_ready_summary_packet.json"
out_summary = generated / "ax22_strict_extension_lane_publication_ready_summary_packet_summary.json"

ax19 = json.loads(in_ax19.read_text(encoding="utf-8"))
ax21 = json.loads(in_ax21.read_text(encoding="utf-8"))

theta = ax21.get("theta_pair_extension_lane_representative", {})

packet = {
    "lane": "strict_extension_only",
    "step": "AX22",
    "status": "AX22_EXECUTED_STRICT_EXTENSION_LANE_PUBLICATION_READY_SUMMARY_PACKET_NO_FALSE_PASS",
    "as_of": "2026-03-12",
    "scope": {
        "strict_core_changed": False,
        "strict_extension_only": True,
        "axiom_augmented_only": False,
    },
    "summary": {
        "slot_closure_premises": {
            "eps": ax21.get("inputs", {}).get("eps_sigma_int_E_pair_amplitude"),
            "delta_d": ax21.get("inputs", {}).get("delta_d_sigma_int_positive_window_step"),
            "notes": [
                "Both eps and delta_d are premise-based strict provenance objects on this lane (not strict-derived).",
                "No additional tuning permitted in this scope (AX21 packet).",
            ],
        },
        "theta_pair_representative": {
            "theta_1": theta.get("theta_1"),
            "theta_2": theta.get("theta_2"),
            "u_1": theta.get("u_1"),
            "u_2": theta.get("u_2"),
            "artifact": theta.get("artifact"),
        },
        "residual_orientation_role": {
            "closure_packet": "generated/strict_extension_lane_sigma_int_residual_orientation_closure_packet.json (AX19)",
            "bridge_instance_artifact": ax19.get("r1_target_slot_population", {}).get("artifact"),
        },
    },
    "assembled_from": {
        "AX19": "generated/strict_extension_lane_sigma_int_residual_orientation_closure_packet.json",
        "AX21": "generated/strict_extension_lane_sigma_int_to_theta_selector_ingredient_closure_packet.json",
        "notes": "AX22 is a publication-ready assembly only; it adds no new strict-core claims.",
    },
    "strict_core_frontier": {
        "T159_discharged": False,
        "T160_T161_discharged": False,
        "strict_core_theta_export_present": False,
        "QW_2191_discharged": False,
    },
    "forbidden_overclaim_set": sorted(
        set(ax19.get("forbidden_claims", [])) | set(ax21.get("forbidden_overclaim_set", []))
    ),
    "no_false_pass": True,
}

out_packet.write_text(json.dumps(packet, indent=2) + "\n", encoding="ascii")

summary = {
    "step": "AX22",
    "status": packet["status"],
    "goal": "Assemble strict-extension lane (AX19 + AX21) into one publication-ready summary carrier without promoting to strict core.",
    "created_files": [
        "generated/strict_extension_lane_publication_ready_summary_packet.json",
        "generated/ax22_strict_extension_lane_publication_ready_summary_packet_summary.json",
    ],
    "strict_core_changed": False,
    "no_false_pass": True,
}

out_summary.write_text(json.dumps(summary, indent=2) + "\n", encoding="ascii")
print(out_summary)

