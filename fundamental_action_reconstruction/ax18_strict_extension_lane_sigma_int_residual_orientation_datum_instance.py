from __future__ import annotations

import json
from pathlib import Path

root = Path(__file__).resolve().parent
generated = root / "generated"

in_sigma_int = generated / "sigma_int_strict_derived_v1.json"
in_delta_d = generated / "delta_d_sigma_int_positive_window_step_strict_provenance_v1.json"
in_theta_pair = generated / "theta_pair_sigma_int_strict_selector_ingredient_o2_cut_candidate_v1.json"
in_r1_population = generated / "r1_residual_orientation_datum_target_slot_population_candidate_from_sigma_int_theta_pair_v1.json"

out_instance = generated / "strict_extension_lane_sigma_int_residual_orientation_datum_instance.json"
out_summary = generated / "ax18_strict_extension_lane_sigma_int_residual_orientation_datum_instance_summary.json"

sigma_int_obj = json.loads(in_sigma_int.read_text(encoding="utf-8"))
delta_d_obj = json.loads(in_delta_d.read_text(encoding="utf-8"))
theta_pair_obj = json.loads(in_theta_pair.read_text(encoding="utf-8"))
r1_obj = json.loads(in_r1_population.read_text(encoding="utf-8"))

instance = {
    "lane": "strict_extension_only",
    "step": "AX18",
    "status": "strict_extension_lane_instance_constructed__no_false_pass",
    "as_of": "2026-03-11",
    "intent": (
        "Materialize one explicit strict-extension-lane sigma_int -> residual orientation datum instance "
        "by citing the exported strict sigma_int datum and the exported sigma_int->theta candidate selector ingredient, "
        "citing the exported strict-provenance delta_d value object (premise-based), and attaching the exported R1 target-slot candidate inhabitant "
        "constructed from that theta pair (P402)."
    ),
    "inputs": {
        "sigma_int_strict_derived_v1": {
            "value": sigma_int_obj.get("value"),
            "artifact": "generated/sigma_int_strict_derived_v1.json",
            "export": "F307/N418",
            "note": "Strict-side source upgrade by explicit premise; not a legacy->strict bridge.",
        },
        "delta_d_sigma_int_positive_window_step_strict_provenance_v1": {
            "value": delta_d_obj.get("value"),
            "artifact": "generated/delta_d_sigma_int_positive_window_step_strict_provenance_v1.json",
            "export": "F328/N440",
            "note": "Strict-side delta_d value object (strict-source-upgraded by explicit premise); not a uniqueness claim.",
        }
    },
    "extension_scope_premises": {
        "strict_side_admissibility_principle": {"accepted": True, "packet": "AX16", "scope": "strict_extension_only"},
        "delta_d_value_object": {
            "accepted": True,
            "packet": "F328/N440",
            "scope": "strict_side__strict_source_upgraded",
            "artifact": "generated/delta_d_sigma_int_positive_window_step_strict_provenance_v1.json",
            "value": delta_d_obj.get("value"),
            "rule": "delta_d := delta_max = d_local/11 (corridor saturation premise)",
            "note": "Dedicated strict-side delta_d value object; premise-based; not a strict-derived uniqueness claim.",
        },
        "theta_supply_selector_ingredient": {
            "accepted": True,
            "packet": "AX18",
            "scope": "strict_extension_only",
            "artifact": "generated/theta_pair_sigma_int_strict_selector_ingredient_o2_cut_candidate_v1.json",
            "export": "F325/N436",
            "note": (
                "A strict-side candidate selector ingredient is cited as the chosen representative theta supply in "
                "this extension scope. It does not discharge QW-2191 in strict core."
            ),
        },
    },
    "theta_pair": {
        "theta_1_cand": float(theta_pair_obj["outputs"]["pair1"]["theta_1_cand"]),
        "theta_2_cand": float(theta_pair_obj["outputs"]["pair2"]["theta_2_cand"]),
        "provenance": theta_pair_obj.get("provenance", {}),
        "delta_d_recorded": theta_pair_obj.get("positive_window_corridor", {}).get("numeric", {}).get("delta_d"),
        "delta_d_value_object": {"artifact": "generated/delta_d_sigma_int_positive_window_step_strict_provenance_v1.json", "export": "F328/N440"},
    },
    "r1_target_slot_population": {
        "artifact": "generated/r1_residual_orientation_datum_target_slot_population_candidate_from_sigma_int_theta_pair_v1.json",
        "export": "P402",
        "u_1_cand": r1_obj["outputs"]["u_1_cand"],
        "u_2_cand": r1_obj["outputs"]["u_2_cand"],
        "orientation_slice_candidate": r1_obj["outputs"]["orientation_slice_candidate"],
        "audits": r1_obj.get("audits", {}),
    },
    "strict_core_status": {
        "strict_core_theta_export_present": False,
        "strict_core_target_slot_population_present": False,
        "admissible_S_sel_int_present": False,
        "strict_core_selector_closure_present": False,
        "QW_2191_discharged": False,
        "object_support_above_export_map_object_discharged": False,
        "ToE_closed": False,
    },
    "forbidden_overclaim_set": [
        "strict_extension_is_not_strict_core",
        "no_strict_core_theta_export",
        "no_strict_core_target_slot_population_by_export_map_object",
        "no_object_support_discharge_above_export_map_object",
        "no_admissible_S_sel_int",
        "no_strict_core_selector_closure",
        "no_QW_2191_discharge",
        "no_ToE_closure",
    ],
    "notes": [
        "This instance exists only to keep one explicit reproducible strict-extension sigma_int->R1 population record.",
        "It does not claim equivalence with the axiom-augmented lane (AX1..AX6).",
    ],
}

out_instance.write_text(json.dumps(instance, indent=2) + "\n", encoding="ascii")

summary = {
    "step": "AX18",
    "status": "AX18_EXECUTED_STRICT_EXTENSION_LANE_SIGMA_INT_RESIDUAL_ORIENTATION_DATUM_INSTANCE_PACKET_NO_FALSE_PASS",
    "goal": "Create one strict-extension-only sigma_int -> residual orientation datum instance from the exported sigma_int->theta candidate ingredient, citing the exported delta_d value object.",
    "created_files": [
        "generated/strict_extension_lane_sigma_int_residual_orientation_datum_instance.json",
        "generated/ax18_strict_extension_lane_sigma_int_residual_orientation_datum_instance_summary.json",
    ],
    "result": {
        "extension_scope": "strict_extension_only",
        "theta_supply_present_extension_only": True,
        "r1_population_present_extension_only": True,
        "strict_core_changed": False,
    },
    "hard_limits": instance["forbidden_overclaim_set"],
    "next_honest_step": "AX19_OR_DECLARE_NEW_STRICT_CORE_SELECTOR_SOURCE",
}

out_summary.write_text(json.dumps(summary, indent=2) + "\n", encoding="ascii")
print(out_summary)
