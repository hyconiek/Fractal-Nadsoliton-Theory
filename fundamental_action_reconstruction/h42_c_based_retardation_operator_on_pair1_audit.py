#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "generated"
OUT.mkdir(exist_ok=True)

omega = 1.0
L = 1.0
c = 1.0
psi0 = 0.3
L_parallel = 0.8
L_perp = 1.2

phi_iso = omega * L / c
iso_value = math.cos(phi_iso)

delta_ret = math.cos(omega * L_parallel / c) - math.cos(omega * L_perp / c)

payload = {
    "step": "H42",
    "status": "PASS_PARTIAL_C_BASED_RETARDATION_SPLIT_ESTABLISHED",
    "date": "2026-03-07",
    "pair": "pair1=(c_1,s_1)",
    "minimal_phase": "phi_ret(L;c)=omega*L/c",
    "case_without_psi0": {
        "operator": "K_ret,iso^(1)(c;L)=cos(omega*L/c)*I_2",
        "omega": omega,
        "L": L,
        "c": c,
        "phi_ret": phi_iso,
        "common_scalar": iso_value,
        "selector_effect": "trivial",
        "breaks_O2": False,
    },
    "case_with_psi0": {
        "operator": "K_ret,psi0^(1)=R(psi0) diag(cos(omega*L_parallel/c), cos(omega*L_perp/c)) R(-psi0)",
        "omega": omega,
        "c": c,
        "psi0": psi0,
        "L_parallel": L_parallel,
        "L_perp": L_perp,
        "delta_ret": delta_ret,
        "new_effect_possible": abs(delta_ret) > 0.0,
        "anchor_origin": "imported_psi0",
    },
    "frontier": "H42_B1 := a minimal c-based retardation operator on pair1 is selector-trivial without an imported anchor and can yield a genuine new spectral/response split only through psi0-dependent anisotropic path data, so c controls coupling speed but does not by itself generate selector orientation",
    "hard_limits": {
        "theorem_level_pass": False,
        "full_closure_pass": False,
        "c_alone_selects_orientation": False,
        "psi0_plus_c_is_strict_core": False,
        "qw_2191_discharged": False,
    },
}

summary = {
    "step": "H42",
    "status": payload["status"],
    "frontier": payload["frontier"],
    "without_psi0_breaks_O2": payload["case_without_psi0"]["breaks_O2"],
    "with_psi0_new_effect_possible": payload["case_with_psi0"]["new_effect_possible"],
    "hard_limits": payload["hard_limits"],
}

(OUT / "h42_c_based_retardation_operator_on_pair1_audit.json").write_text(
    json.dumps(payload, indent=2, ensure_ascii=True) + "\n",
    encoding="utf-8",
)
(OUT / "h42_c_based_retardation_operator_on_pair1_audit_summary.json").write_text(
    json.dumps(summary, indent=2, ensure_ascii=True) + "\n",
    encoding="utf-8",
)

print(json.dumps(summary, indent=2, ensure_ascii=True))
