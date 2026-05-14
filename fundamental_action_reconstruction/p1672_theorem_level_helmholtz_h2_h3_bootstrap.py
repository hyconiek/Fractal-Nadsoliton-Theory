#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "generated" / "p1672_s622_theorem_level_helmholtz_h2_h3_bootstrap.json"

# Closed-form proxy for coupled scalar-gauge channel (strict-side local model)
# L = 1/2 h'^2 + 1/2 zeta a'^2 -1/2 omega^2 h^2 -1/2 a^2 -xi h a + lam h^2 a^2
# lam := beta/(1+eta)

E_h = "h'' + omega^2*h + xi*a - 2*lam*h*a^2"
E_a = "zeta*a'' + a + xi*h - 2*lam*h^2*a"

# H2 proxy: dE_h/da - dE_a/dh = (xi - 4*lam*h*a) - (xi - 4*lam*h*a) = 0
h2_expression = "(xi - 4*lam*h*a) - (xi - 4*lam*h*a)"
h2_pass = True

# H3 proxy for algebraic forcing part after removing principal second-derivative terms
# F_h = omega^2*h + xi*a - 2*lam*h*a^2
# F_a = a + xi*h - 2*lam*h^2*a
# dF_h/da - dF_a/dh = (xi - 4*lam*h*a) - (xi - 4*lam*h*a) = 0
h3_expression = "(xi - 4*lam*h*a) - (xi - 4*lam*h*a)"
h3_pass = True

payload = {
    "checkpoint": "P1672_S622_THEOREM_LEVEL_HELMHOLTZ_H2_H3_BOOTSTRAP",
    "strict_only": True,
    "legacy_bridge_used": False,
    "chain": "K_strict -> coeff -> full L_total -> EOM; reverse bootstrap: EOM -> L_total",
    "full_lagrangian_density_proxy": "1/2 h'^2 + 1/2 zeta a'^2 -1/2 omega^2 h^2 -1/2 a^2 -xi h a + lam h^2 a^2",
    "eom": {"E_h": E_h, "E_a": E_a},
    "helmholtz_h2": {"expression": h2_expression, "proxy_pass": h2_pass},
    "helmholtz_h3": {"expression": h3_expression, "proxy_pass": h3_pass},
    "status": "OPEN_OBLIGATION",
    "open_obligations": [
        "helmholtz_h4_global_atlas_overlap_consistency",
        "full_operator_class_extension_beyond_proxy",
        "spin2_spin0_unitarity_theorem",
        "qg_renormalization_theorem",
        "background_independence_theorem"
    ],
    "next_honest_step": "S623: H4 global atlas theorem + integration with strict QG theorem package.",
    "lay_summary": "Lokalny bootstrap Helmholtza H2/H3 jest zgodny dla sprzężonego kanału, ale pełne globalne domknięcie i bramki QG pozostają otwarte."
}

OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
print(f"Wrote {OUT}")
