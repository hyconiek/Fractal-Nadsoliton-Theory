from __future__ import annotations

import json
from pathlib import Path

root = Path(__file__).resolve().parent
generated = root / "generated"

in_theta_candidate = generated / "theta_pair_sigma_int_strict_selector_ingredient_o2_cut_candidate_v1.json"
in_eps = generated / "eps_sigma_int_E_pair_amplitude_strict_provenance_v1.json"
in_delta_d = generated / "delta_d_sigma_int_positive_window_step_strict_provenance_v1.json"
in_sigma_int = generated / "sigma_int_strict_derived_v1.json"

out_packet = generated / "strict_extension_lane_sigma_int_to_theta_selector_ingredient_closure_packet.json"
out_summary = (
    generated / "ax21_strict_extension_lane_sigma_int_to_theta_selector_ingredient_closure_packet_summary.json"
)

theta_cand = json.loads(in_theta_candidate.read_text(encoding="utf-8"))
eps_obj = json.loads(in_eps.read_text(encoding="utf-8"))
delta_obj = json.loads(in_delta_d.read_text(encoding="utf-8"))
sigma_int_obj = json.loads(in_sigma_int.read_text(encoding="utf-8"))

theta_1 = float(theta_cand["outputs"]["pair1"]["theta_1_cand"])
theta_2 = float(theta_cand["outputs"]["pair2"]["theta_2_cand"])

packet = {
    "lane": "strict_extension_only",
    "step": "AX21",
    "status": "strict_extension_lane_sigma_int_to_theta_slot_choices_frozen__closure_packet_constructed__no_false_pass",
    "as_of": "2026-03-12",
    "assembled_from": {
        "extension_scope_acceptance": "AX16 (strict_extension_only)",
        "slot_closure_selector_premise_packet": "AX21_STRICT_EXTENSION_LANE_SIGMA_INT_TO_THETA_SLOT_CLOSURE_SELECTOR_PREMISE_PACKET.md",
        "delta_d_maximal_step_convention": "AX17 + generated/delta_d_sigma_int_positive_window_step_strict_provenance_v1.json (F328/N440)",
        "eps_value_object": "generated/eps_sigma_int_E_pair_amplitude_strict_provenance_v1.json (F317/N428)",
        "theta_candidate_selector_ingredient": "generated/theta_pair_sigma_int_strict_selector_ingredient_o2_cut_candidate_v1.json (F325/N436)",
        "q2191_compatibility_witness": "P401",
        "strict_admissibility_boundaries": ["P408", "N446", "N447", "N448", "N442", "N443"],
    },
    "inputs": {
        "sigma_int_strict_derived_v1": sigma_int_obj.get("value"),
        "eps_sigma_int_E_pair_amplitude": eps_obj.get("value"),
        "delta_d_sigma_int_positive_window_step": delta_obj.get("value"),
    },
    "extension_scope_premises": {
        "lane": "strict_extension_only",
        "slot_fixes": {
            "eps": "eps_sigma_int_E_pair_amplitude_strict_provenance_v1 := 1/2 (premise-based; not strict-derived)",
            "delta_d": "delta_d_sigma_int_positive_window_step_strict_provenance_v1 := delta_max (premise-based; not strict-derived)",
        },
        "no_additional_tuning_permitted": True,
    },
    "theta_pair_extension_lane_representative": {
        "artifact": "generated/theta_pair_sigma_int_strict_selector_ingredient_o2_cut_candidate_v1.json",
        "theta_1": theta_1,
        "theta_2": theta_2,
        "u_1": theta_cand["outputs"]["pair1"]["u_1_cand"],
        "u_2": theta_cand["outputs"]["pair2"]["u_2_cand"],
        "o2_cut_argument": theta_cand.get("o2_cut_argument", {}),
        "scope_note": "This is an extension-lane representative point on the QW-2191 O(2) family; strict-core canonical upgrade (T159) remains open.",
    },
    "strict_core_status": {
        "T159_strict_core_upgrade_discharged": False,
        "QW_2191_discharged": False,
        "strict_core_theta_export_present": False,
    },
    "forbidden_overclaim_set": [
        "strict-derived eps selection (T160)",
        "strict-derived delta_d selection (T161)",
        "strict-core theta export",
        "strict-core selector closure",
        "strict-core QW-2191 discharge",
        "ToE closure",
    ],
    "notes": [
        "AX21 packages a strict-extension-only slot-closure premise: eps=1/2 and delta_d=delta_max.",
        "No strict-core promotion is performed; P408/N446/N447/N448 remain in force.",
    ],
}

out_packet.write_text(json.dumps(packet, indent=2) + "\n", encoding="ascii")

summary = {
    "step": "AX21",
    "status": "AX21_EXECUTED_STRICT_EXTENSION_LANE_SIGMA_INT_TO_THETA_SLOT_CLOSURE_SELECTOR_PREMISE_PACKET_NO_FALSE_PASS",
    "goal": "Freeze eps and delta_d selector slots in strict_extension_only scope and assemble one reproducible theta-pair representative without promoting to strict core.",
    "created_files": [
        "generated/strict_extension_lane_sigma_int_to_theta_selector_ingredient_closure_packet.json",
        "generated/ax21_strict_extension_lane_sigma_int_to_theta_selector_ingredient_closure_packet_summary.json",
    ],
    "theta_pair": {"theta_1": theta_1, "theta_2": theta_2},
    "strict_core_changed": False,
    "T159_discharged": False,
    "QW_2191_discharged": False,
    "ToE_closed": False,
    "no_false_pass": True,
}

out_summary.write_text(json.dumps(summary, indent=2) + "\n", encoding="ascii")
print(out_summary)

