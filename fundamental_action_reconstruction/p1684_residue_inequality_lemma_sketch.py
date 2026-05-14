#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "generated" / "p1684_s634_residue_inequality_lemma_sketch.json"

R_max = 1.0
params = {
    "Z2_0": 1.0,
    "Z0_0": 1.0,
    "a2": 0.10,
    "a0": 0.06,
    "c1": 0.08,
    "c2": 0.04,
    "xi": 0.03,
}

z2_min = params["Z2_0"] - abs(params["a2"]) * R_max
z0_min = params["Z0_0"] - abs(params["a0"]) * R_max
k2_min = z2_min * (1.0 - 0.04 * R_max)
k0_min = z0_min * (1.0 - 0.02 * R_max)

payload = {
    "checkpoint": "P1684_S634_RESIDUE_INEQUALITY_LEMMA_SKETCH",
    "strict_only": True,
    "legacy_bridge_used": False,
    "pipeline": "K_strict -> coeff -> full_L_total -> linearized_EOM -> residue inequality sketch",
    "domain": {"R_max": R_max, "curvature_domain": "|R|<=R_max"},
    "parameter_point": params,
    "inequality_sketch": {
        "Z2_lower_bound": f"Z2(R)>= {z2_min:.3f}",
        "Z0_lower_bound": f"Z0(R)>= {z0_min:.3f}",
        "detK2_lower_bound": f"det(K2)(R)>= {k2_min:.3f}",
        "detK0_lower_bound": f"det(K0)(R)>= {k0_min:.3f}"
    },
    "status": "OPEN_OBLIGATION",
    "limitation": "Lemma sketch only; full theorem-level variational proof remains missing.",
    "next_honest_step": "S635: symbolic bound certificate for (c1,c2,xi) with boundary counterexample map."
}

OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
print(f"Wrote {OUT}")
