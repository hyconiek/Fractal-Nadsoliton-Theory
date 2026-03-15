from __future__ import annotations

import json
import math
from pathlib import Path

root = Path(__file__).resolve().parent
generated = root / "generated"

in_alpha_geo = generated / "alpha_geo_strict_derived_v1.json"
in_p439 = generated / "p439_current_strict_qw2191_weighted_kl_reference_objective_o2_cut_audit_probe.json"

out_packet = generated / "strict_extension_lane_z12_generator_fix_and_shannon_kl_o2_cut_closure_packet.json"
out_summary = (
    generated / "ax23_strict_extension_lane_z12_generator_fix_and_shannon_kl_o2_cut_closure_packet_summary.json"
)

alpha_geo_obj = json.loads(in_alpha_geo.read_text(encoding="utf-8"))
p439 = json.loads(in_p439.read_text(encoding="utf-8"))

alpha_geo_value = alpha_geo_obj.get("value")
alpha_geo_numeric = float(4.0 * math.log(2.0))

ref_id = "r_exp_minus_alpha_geo_directed_d"
obj_id = "KL_u1_to_r"

ref_block = p439["results_by_reference"][ref_id]
row = None
for r in ref_block["objective_rows_compact"]:
    if r["objective_id"] == obj_id:
        row = r
        break
if row is None:
    raise RuntimeError(f"missing objective row {obj_id} for reference {ref_id} in P439 output")

centers = [float(x) for x in row["cluster_centers_theta"]]
centers_sorted = sorted((c % (2.0 * math.pi)) for c in centers)
theta_fix = float(centers_sorted[0] % math.pi)  # representative in [0,pi)

packet = {
    "lane": "strict_extension_only",
    "step": "AX23",
    "status": "strict_extension_lane_z12_generator_fix_and_shannon_kl_pair1_o2_cut_premise_packet_constructed__no_false_pass",
    "as_of": "2026-03-14",
    "assembled_from": {
        "extension_scope_acceptance": "AX16 (strict_extension_only)",
        "premise_packet": "AX23_STRICT_EXTENSION_LANE_Z12_GENERATOR_FIX_AND_SHANNON_KL_O2_CUT_PREMISE_PACKET.md",
        "alpha_geo_source_upgrade": "generated/alpha_geo_strict_derived_v1.json (F309/N420)",
        "p439_objective_audit": "generated/p439_current_strict_qw2191_weighted_kl_reference_objective_o2_cut_audit_probe.json",
        "strict_core_noncanonicity_boundary": "N462 (no Aut(Z12)-invariant canonical generator choice); T164 remains open",
        "strict_core_selector_target": "T165 remains open",
    },
    "inputs": {
        "alpha_geo_strict_derived_v1": alpha_geo_value,
        "alpha_geo_numeric_used": alpha_geo_numeric,
        "p439_reference_id": ref_id,
        "p439_objective_id": obj_id,
    },
    "extension_scope_premises": {
        "lane": "strict_extension_only",
        "z12_generator_orientation_fix": {
            "generator_fixed": 1,
            "successor_map": "suc_fix(k) := (k+1) mod 12",
            "note": "This fixes a marked direction on Z_12 and is not strict-core canonical under N462; it is accepted only in extension scope.",
        },
        "reference_distribution": {
            "r_dir(x)": "exp(-alpha_geo * x) / Σ_{t=0..11} exp(-alpha_geo * t)",
            "x_domain": "I_12_v1 = {0,1,...,11}",
            "note": "Marked-site + marked-direction reference distribution (non-translation-invariant).",
        },
    },
    "pair1_o2_cut_extension_lane_representative": {
        "objective": "J_dir(theta) := KL( u1(theta)^2 || r_dir )",
        "theta_pair1_fix": theta_fix,
        "residual_symmetry": "Z2 (theta -> theta + π) unavoidable for squared-amplitude objectives",
        "p439_grid_audit": {
            "cluster_count": int(row["cluster_count"]),
            "cluster_centers_theta": [float(x) for x in row["cluster_centers_theta"]],
            "z2_unique_on_grid": bool(row.get("z2_unique_on_grid")),
        },
        "scope_note": "This is an extension-lane representative point on the QW-2191 O(2) family on pair1; strict-core upgrade requires T164+T165 discharge.",
    },
    "strict_core_status": {
        "T164_discharged": False,
        "T165_discharged": False,
        "QW_2191_discharged": False,
        "strict_core_theta_export_present": False,
    },
    "forbidden_overclaim_set": [
        "strict-core canonical generator/orientation fixing datum",
        "strict-core Shannon symmetry-breaking selector ingredient",
        "strict-core theta export",
        "strict-core selector closure",
        "strict-core QW-2191 discharge",
        "ToE closure",
    ],
    "no_false_pass": True,
}

generated.mkdir(exist_ok=True)
out_packet.write_text(json.dumps(packet, indent=2) + "\n", encoding="ascii")

summary = {
    "step": "AX23",
    "status": "AX23_EXECUTED_STRICT_EXTENSION_LANE_Z12_GENERATOR_FIX_AND_SHANNON_KL_O2_CUT_PREMISE_PACKET_NO_FALSE_PASS",
    "goal": "Freeze a Z_12 generator/orientation and one Shannon/KL pair1 representative theta in strict_extension_only scope, without promoting to strict core.",
    "created_files": [
        "generated/strict_extension_lane_z12_generator_fix_and_shannon_kl_o2_cut_closure_packet.json",
        "generated/ax23_strict_extension_lane_z12_generator_fix_and_shannon_kl_o2_cut_closure_packet_summary.json",
    ],
    "theta_pair1_fix": theta_fix,
    "strict_core_changed": False,
    "T164_discharged": False,
    "T165_discharged": False,
    "QW_2191_discharged": False,
    "ToE_closed": False,
    "no_false_pass": True,
}

out_summary.write_text(json.dumps(summary, indent=2) + "\n", encoding="ascii")
print(out_summary)

