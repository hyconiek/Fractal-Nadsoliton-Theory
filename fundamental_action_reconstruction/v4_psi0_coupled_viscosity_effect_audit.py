#!/usr/bin/env python3
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "generated"
OUT.mkdir(exist_ok=True)

payload = {
    "step": "V4",
    "status": "PASS_PARTIAL_COMPETING_EXTENSION_HYPOTHESIS_ONLY",
    "date": "2026-03-06",
    "carrier": {
        "space": "V_1 = span{c_1,s_1}",
        "basis": ["c_1", "s_1"],
        "anchor": "psi0",
    },
    "coupled_operator": "K_visc_aniso^(1)(psi0) = R(psi0) diag(nu_parallel, nu_perp) R(-psi0)",
    "pair_level_consequence": {
        "parallel_mode": "u_parallel(psi0) = cos(psi0)c_1 + sin(psi0)s_1",
        "perp_mode": "u_perp(psi0) = -sin(psi0)c_1 + cos(psi0)s_1",
        "anchor_generating": False,
        "anchor_amplifying": True,
        "can_split_response_relative_to_anchor": True,
    },
    "frontier": "V4_B1 := coupling anisotropic viscosity to psi0 yields a nontrivial pair-level anchor-amplifying operator on V_1, but it does not generate the anchor itself and therefore cannot by itself close the selector problem",
    "hard_limits": {
        "theorem_level_pass": False,
        "full_closure_pass": False,
        "viscosity_generates_psi0": False,
        "coupled_lane_replaces_primary_psi0_lane": False,
        "qw2191_discharged": False,
    },
}

json_path = OUT / "v4_psi0_coupled_viscosity_effect_audit.json"
summary_path = OUT / "v4_psi0_coupled_viscosity_effect_audit_summary.json"
json_path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n")
summary_path.write_text(json.dumps({
    "step": payload["step"],
    "status": payload["status"],
    "carrier": payload["carrier"],
    "pair_level_consequence": payload["pair_level_consequence"],
    "frontier": payload["frontier"],
    "hard_limits": payload["hard_limits"],
}, indent=2, ensure_ascii=False) + "\n")
