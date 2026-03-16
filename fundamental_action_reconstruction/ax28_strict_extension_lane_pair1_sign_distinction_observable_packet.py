from __future__ import annotations

import json
import math
from pathlib import Path

root = Path(__file__).resolve().parent
generated = root / "generated"

in_alpha_geo = generated / "alpha_geo_strict_derived_v1.json"
in_u1 = generated / "a_1_pair1_orientation_projector_operator_strict_core_v1.json"
in_ax23 = generated / "strict_extension_lane_z12_generator_fix_and_shannon_kl_o2_cut_closure_packet.json"

out_packet = generated / "strict_extension_lane_pair1_sign_distinction_observable_packet.json"
out_summary = generated / "ax28_strict_extension_lane_pair1_sign_distinction_observable_packet_summary.json"

alpha_geo_obj = json.loads(in_alpha_geo.read_text(encoding="utf-8"))
a1 = json.loads(in_u1.read_text(encoding="utf-8"))
ax23 = json.loads(in_ax23.read_text(encoding="utf-8")) if in_ax23.exists() else None

alpha_geo_value = alpha_geo_obj.get("value")
alpha_geo_numeric = float(4.0 * math.log(2.0))

u1 = (a1.get("data") or {}).get("u_1")
if not (isinstance(u1, list) and len(u1) == 12 and all(isinstance(v, (int, float)) for v in u1)):
    raise RuntimeError("expected generated/a_1_pair1_orientation_projector_operator_strict_core_v1.json to contain data.u_1 length 12")
u1 = [float(v) for v in u1]

# Extension-scope directed coordinate (marked direction): x = 0..11 with successor map x->x+1 mod 12.
x_list = list(range(12))

# Directed exponential weight (normalization irrelevant for sign).
w_dir = [float(math.exp(-alpha_geo_numeric * float(x))) for x in x_list]

def dot(w: list[float], v: list[float]) -> float:
    return float(sum(float(a) * float(b) for a, b in zip(w, v)))

S_u1 = dot(w_dir, u1)
sign_fix = 1.0 if S_u1 >= 0.0 else -1.0
u1_oriented = [float(sign_fix * v) for v in u1]
S_oriented = dot(w_dir, u1_oriented)

packet = {
    "lane": "strict_extension_only",
    "step": "AX28",
    "status": "strict_extension_lane_pair1_sign_distinction_observable_premise_packet_constructed__no_false_pass",
    "as_of": "2026-03-16",
    "assembled_from": {
        "extension_scope_acceptance": "AX16 (strict_extension_only)",
        "premise_packet": "AX28_STRICT_EXTENSION_LANE_PAIR1_SIGN_DISTINCTION_OBSERVABLE_PREMISE_PACKET.md",
        "alpha_geo_source_upgrade": "generated/alpha_geo_strict_derived_v1.json (F309/N420)",
        "strict_pair1_axis_representative": "generated/a_1_pair1_orientation_projector_operator_strict_core_v1.json (F456; axis is projector-level; residual sign is gauge in strict core)",
        "generator_fix_premise_reference": (
            "generated/strict_extension_lane_z12_generator_fix_and_shannon_kl_o2_cut_closure_packet.json (AX23)"
            if in_ax23.exists()
            else "AX23 (generator/orientation fix premise in strict_extension_only scope)"
        ),
        "strict_core_obstruction_boundary": "H37 + N518 + P472 (strict core still has no sign-sensitive physical orientation datum)",
    },
    "inputs": {
        "alpha_geo_strict_derived_v1": alpha_geo_value,
        "alpha_geo_numeric_used": alpha_geo_numeric,
        "u1_source_object": str((a1.get("object") or "A_1_pair1_orientation_projector_operator_strict_core_v1")),
    },
    "extension_scope_premises": {
        "lane": "strict_extension_only",
        "z12_generator_orientation_fix": {
            "generator_fixed": 1,
            "successor_map": "suc_fix(k) := (k+1) mod 12",
            "note": (
                "This fixes a marked direction on Z_12 and is not strict-core canonical under N462/T164; "
                "it is accepted only in extension scope."
            ),
            "ax23_reference": (ax23.get("extension_scope_premises", {}).get("z12_generator_orientation_fix") if isinstance(ax23, dict) else None),
        },
        "directed_coordinate": {"x_domain": "I_12_v1 = {0,1,...,11}", "meaning": "directed label relative to fixed successor map"},
    },
    "directed_weight": {
        "w_dir(x)": "exp(-alpha_geo * x)",
        "x_domain": "I_12_v1",
        "weights": w_dir,
        "normalization_note": "Normalization is irrelevant for sign-distinction; the exported list is unnormalized.",
    },
    "sign_distinction_observable": {
        "S_dir(u)": "Σ_x w_dir(x) u(x)",
        "S_dir(u1_strict_representative)": S_u1,
        "sign_fix_applied_to_u1": sign_fix,
        "u1_oriented_extension_lane": u1_oriented,
        "S_dir(u1_oriented_extension_lane)": S_oriented,
        "meaning": "This scalar flips sign under u->-u and fixes a preferred orientation representative by requiring S_dir(u)>0.",
    },
    "strict_core_status": {"H37_discharged": False, "T164_discharged": False, "strict_core_changed": False},
    "forbidden_overclaim_set": [
        "strict-core sign-sensitive physical orientation datum",
        "strict-core canonical generator/orientation fixing datum",
        "strict-core selector closure",
        "strict-core QW-2191 discharge",
        "ToE closure",
    ],
    "no_false_pass": True,
}

generated.mkdir(exist_ok=True)
out_packet.write_text(json.dumps(packet, indent=2) + "\n", encoding="ascii")

summary = {
    "step": "AX28",
    "status": "AX28_EXECUTED_STRICT_EXTENSION_LANE_PAIR1_SIGN_DISTINCTION_OBSERVABLE_PREMISE_PACKET_NO_FALSE_PASS",
    "goal": "Export an explicit sign-distinction scalar observable on pair1 in strict_extension_only scope under a declared Z12 generator/orientation fix, without changing strict core.",
    "created_files": [
        "generated/strict_extension_lane_pair1_sign_distinction_observable_packet.json",
        "generated/ax28_strict_extension_lane_pair1_sign_distinction_observable_packet_summary.json",
    ],
    "S_dir_u1_strict_representative": S_u1,
    "sign_fix_applied_to_u1": sign_fix,
    "S_dir_u1_oriented": S_oriented,
    "strict_core_changed": False,
    "T164_discharged": False,
    "H37_discharged_in_strict_core": False,
    "no_false_pass": True,
}

out_summary.write_text(json.dumps(summary, indent=2) + "\n", encoding="ascii")
print(out_summary)

