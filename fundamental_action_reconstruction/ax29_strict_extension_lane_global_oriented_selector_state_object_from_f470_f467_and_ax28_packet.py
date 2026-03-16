from __future__ import annotations

import json
import math
from pathlib import Path

root = Path(__file__).resolve().parent
generated = root / "generated"

in_f470 = generated / "selector_state_global_c_v1_projective_strict_v1.json"
in_f467 = generated / "u_12345_pair12345_chart_glued_orientation_vector_section_oriented_mod_2pi_strict_convention_v1.json"
in_ax28 = generated / "strict_extension_lane_pair1_sign_distinction_observable_packet.json"

out_packet = generated / "strict_extension_lane_selector_state_global_c_v1_oriented_vector_v1.json"
out_summary = (
    generated
    / "ax29_strict_extension_lane_global_oriented_selector_state_object_from_f470_f467_and_ax28_packet_summary.json"
)

f470 = json.loads(in_f470.read_text(encoding="utf-8"))
f467 = json.loads(in_f467.read_text(encoding="utf-8"))
ax28 = json.loads(in_ax28.read_text(encoding="utf-8"))

u_vecs = (f467.get("outputs") or {}).get("u_vectors") or {}
u1 = u_vecs.get("u_1")
if not (isinstance(u1, list) and len(u1) == 12 and all(isinstance(v, (int, float)) for v in u1)):
    raise RuntimeError("expected F467 outputs.u_vectors.u_1 length 12")

sign_fix = float((ax28.get("sign_distinction_observable") or {}).get("sign_fix_applied_to_u1") or 1.0)
if sign_fix not in (-1.0, 1.0):
    raise RuntimeError("expected AX28 sign_fix_applied_to_u1 to be +/-1")

w_dir = (ax28.get("directed_weight") or {}).get("weights")
if not (isinstance(w_dir, list) and len(w_dir) == 12 and all(isinstance(v, (int, float)) for v in w_dir)):
    raise RuntimeError("expected AX28 directed_weight.weights length 12")
w_dir = [float(v) for v in w_dir]

def dot(w: list[float], v: list[float]) -> float:
    return float(sum(float(a) * float(b) for a, b in zip(w, v)))

u_oriented = {}
for key, vec in u_vecs.items():
    if not (isinstance(vec, list) and len(vec) == 12 and all(isinstance(v, (int, float)) for v in vec)):
        continue
    u_oriented[key] = [float(sign_fix * float(v)) for v in vec]

S_u1_oriented = dot(w_dir, u_oriented["u_1"])
if not (S_u1_oriented > 0.0):
    raise RuntimeError("AX29 construction failed: expected S_dir(u1_oriented)>0 by AX28 sign-fix")

packet = {
    "lane": "strict_extension_only",
    "step": "AX29",
    "status": "strict_extension_lane_global_oriented_selector_state_object_constructed_from_strict_projective_state_plus_convention_vector_section_plus_ax28_sign_fix__no_false_pass",
    "as_of": "2026-03-16",
    "assembled_from": {
        "strict_global_projective_state": "generated/selector_state_global_c_v1_projective_strict_v1.json (F470/N516)",
        "strict_convention_vector_section": "generated/u_12345_pair12345_chart_glued_orientation_vector_section_oriented_mod_2pi_strict_convention_v1.json (F467/N511)",
        "extension_sign_fix": "generated/strict_extension_lane_pair1_sign_distinction_observable_packet.json (AX28)",
        "notes": "AX29 adds an oriented vector representative in extension scope by applying the AX28 sign-fix globally to the existing F467 vector section; the underlying strict state remains projective.",
    },
    "strict_inputs": {
        "C_v1_projective_selector_state_object": f470.get("object"),
        "F467_vector_section_object": f467.get("object"),
    },
    "extension_scope_inputs": {
        "sign_fix_applied_to_u1": sign_fix,
        "directed_weight_w_dir": w_dir,
    },
    "outputs": {
        "u_vectors_oriented_extension_lane": u_oriented,
        "S_dir_u1_oriented": S_u1_oriented,
    },
    "meaning": [
        "This is a global oriented vector representative of the already exported strict projective selector state, constructed only in strict_extension_only scope.",
        "Strict core remains unchanged: the physical strict state is still projective/ray-level; this packet does not export a strict sign-sensitive physical orientation datum.",
    ],
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
    "step": "AX29",
    "status": "AX29_EXECUTED_STRICT_EXTENSION_LANE_GLOBAL_ORIENTED_SELECTOR_STATE_OBJECT_FROM_F470_F467_AND_AX28_PACKET_NO_FALSE_PASS",
    "goal": "Export a global oriented vector representative of the strict projective selector state on C_v1 in strict_extension_only scope by applying AX28 sign-fix to the F467 vector section.",
    "created_files": [
        "generated/strict_extension_lane_selector_state_global_c_v1_oriented_vector_v1.json",
        "generated/ax29_strict_extension_lane_global_oriented_selector_state_object_from_f470_f467_and_ax28_packet_summary.json",
    ],
    "sign_fix_applied_to_u1": sign_fix,
    "S_dir_u1_oriented": S_u1_oriented,
    "strict_core_changed": False,
    "T164_discharged": False,
    "H37_discharged_in_strict_core": False,
    "no_false_pass": True,
}

out_summary.write_text(json.dumps(summary, indent=2) + "\n", encoding="ascii")
print(out_summary)

