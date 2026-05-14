#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "generated" / "p1677_s627_combined_unitarity_renormalization_theorem_object.json"

restricted_domain = {
    "beta_range": [0.60, 0.90],
    "eta_range": [1.10, 1.50],
    "curvature_R_range": [0.0, 0.15],
}

ur_conditions = [
    {"id": "U1", "statement": "Z_spin2>0 and Z_spin0>0", "status": "PROXY_READY"},
    {"id": "U2", "statement": "m2_spin2>0 and m2_spin0>0", "status": "PROXY_READY"},
    {"id": "R1", "statement": "counterterm basis closed on restricted domain", "status": "PARTIAL"},
    {"id": "R2", "statement": "running consistency without sign-flip instabilities", "status": "PARTIAL"},
    {"id": "UR_link", "statement": "unitarity-safe poles preserved under admissible running", "status": "MISSING_THEOREM"},
]

payload = {
    "checkpoint": "P1677_S627_COMBINED_UNITARITY_RENORMALIZATION_THEOREM_OBJECT",
    "strict_only": True,
    "legacy_bridge_used": False,
    "chain": "K_strict -> coeff -> full L_total -> EOM -> {unitarity, renormalization} -> combined gate",
    "full_lagrangian_reference": "L_scalar + L_gauge + L_fermion + L_higgs + L_gravity + L_mix",
    "restricted_domain": restricted_domain,
    "conditions": ur_conditions,
    "status": "OPEN_OBLIGATION",
    "missing_theorems": [
        "UR_link_theorem",
        "global_extension_beyond_restricted_domain",
        "background_independence_coupling_theorem",
    ],
    "next_honest_step": "S628: couple combined UR object with background-independence theorem and atlas globalization proof.",
    "lay_summary": "Połączyliśmy warunki unitarności i renormalizacji w jeden formalny obiekt roboczy, ale brakuje końcowego dowodu spinającego je globalnie.",
}

OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False)+"\n", encoding="utf-8")
print(f"Wrote {OUT}")
