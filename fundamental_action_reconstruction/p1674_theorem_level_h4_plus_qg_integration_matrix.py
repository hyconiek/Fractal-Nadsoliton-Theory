#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "generated" / "p1674_s624_theorem_level_h4_plus_qg_integration_matrix.json"

requirements = [
    {"gate": "H4_boundary_theorem", "status": "THEOREM_MISSING", "upstream_proxy": "P1673_PASS"},
    {"gate": "global_atlas_cocycle_theorem", "status": "THEOREM_MISSING", "upstream_proxy": "P1673_PASS"},
    {"gate": "spin2_spin0_unitarity_theorem", "status": "THEOREM_MISSING", "upstream_proxy": "P1665_P1666_PROXY"},
    {"gate": "qg_renormalization_theorem", "status": "THEOREM_MISSING", "upstream_proxy": "P1663_OBLIGATION_MAP"},
    {"gate": "background_independence_theorem", "status": "THEOREM_MISSING", "upstream_proxy": "P1663_OBLIGATION_MAP"},
]

payload = {
    "checkpoint": "P1674_S624_THEOREM_LEVEL_H4_PLUS_QG_INTEGRATION_MATRIX",
    "strict_only": True,
    "legacy_bridge_used": False,
    "forward_chain": "K_strict -> coeff -> full L_total -> EOM",
    "reverse_chain": "EOM -> Helmholtz H1..H4 -> L_total -> coeff -> K_strict",
    "full_lagrangian_reference": "L_scalar + L_gauge + L_fermion + L_higgs + L_gravity + L_mix",
    "requirements": requirements,
    "status": "OPEN_OBLIGATION",
    "critical_path": [
        "spin2_spin0_unitarity_theorem",
        "qg_renormalization_theorem",
        "background_independence_theorem"
    ],
    "next_honest_step": "S625: theorem-level spin2/spin0 unitarity on curved background with ghost-free domain proof.",
    "lay_summary": "Mamy kompletną mapę brakujących twierdzeń: model jest spójny roboczo, ale finalna pewność wymaga ścisłych dowodów QG i globalnej odwracalności.",
}

OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False)+"\n", encoding="utf-8")
print(f"Wrote {OUT}")
