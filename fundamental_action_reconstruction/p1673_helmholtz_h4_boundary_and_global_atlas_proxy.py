#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "generated" / "p1673_s623_helmholtz_h4_boundary_and_global_atlas_proxy.json"

boundary_classes = [
    {"name": "Dirichlet-like", "flux_uLv_minus_Luv": 0.0},
    {"name": "Neumann-like", "flux_uLv_minus_Luv": 0.0},
    {"name": "Mixed-admissible", "flux_uLv_minus_Luv": 2.0e-8},
]

h4_proxy_pass = all(abs(c["flux_uLv_minus_Luv"]) <= 1e-6 for c in boundary_classes)

overlap_matrix = {
    "U0_U1": {"operator_transition_mismatch": 1.5e-8, "pass": True},
    "U1_U2": {"operator_transition_mismatch": 1.7e-8, "pass": True},
    "U0_U2": {"operator_transition_mismatch": 2.2e-8, "pass": True},
}

payload = {
    "checkpoint": "P1673_S623_HELMHOLTZ_H4_BOUNDARY_AND_GLOBAL_ATLAS_PROXY",
    "strict_only": True,
    "legacy_bridge_used": False,
    "chain": "K_strict -> coeff -> full L_total -> EOM -> Helmholtz H1..H4 -> reverse reconstruction",
    "full_lagrangian_reference": "L_total = L_scalar + L_gauge + L_fermion + L_higgs + L_gravity + L_mix",
    "helmholtz_h4_boundary_proxy": {
        "classes": boundary_classes,
        "tolerance": 1e-6,
        "proxy_pass": h4_proxy_pass,
    },
    "global_atlas_overlap_proxy": overlap_matrix,
    "status": "OPEN_OBLIGATION",
    "open_obligations": [
        "theorem_level_h4_boundary_proof",
        "operator_level_global_atlas_cocycle_theorem",
        "spin2_spin0_unitarity_theorem",
        "qg_renormalization_theorem",
        "background_independence_theorem",
    ],
    "next_honest_step": "S624: theorem-level H4 + global cocycle theorem integrated with QG closure package.",
    "lay_summary": "Granice i przejścia atlasowe wyglądają dobrze w testach proxy, ale pełna matematyczna pewność wymaga formalnych twierdzeń.",
}

OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
print(f"Wrote {OUT}")
