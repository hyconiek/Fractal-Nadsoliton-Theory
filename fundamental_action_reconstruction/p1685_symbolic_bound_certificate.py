#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "generated" / "p1685_s635_symbolic_bound_certificate.json"

R_max = 1.0
admissible_domain = {
    "c1": [0.02, 0.20],
    "c2": [0.01, 0.12],
    "xi": [0.00, 0.10],
}

# sample boundary probes (certificate scaffold)
probes = [
    {"c1": 0.08, "c2": 0.04, "xi": 0.03},
    {"c1": 0.20, "c2": 0.12, "xi": 0.10},
    {"c1": 0.01, "c2": 0.04, "xi": 0.03},  # counterexample: c1 below bound
    {"c1": 0.08, "c2": 0.13, "xi": 0.03},  # counterexample: c2 above bound
]

def in_domain(p: dict[str, float]) -> bool:
    return (
        admissible_domain["c1"][0] <= p["c1"] <= admissible_domain["c1"][1]
        and admissible_domain["c2"][0] <= p["c2"] <= admissible_domain["c2"][1]
        and admissible_domain["xi"][0] <= p["xi"] <= admissible_domain["xi"][1]
    )

rows = []
for p in probes:
    ok = in_domain(p)
    rows.append({**p, "admissible": ok})

counterexamples = [r for r in rows if not r["admissible"]]

payload = {
    "checkpoint": "P1685_S635_SYMBOLIC_BOUND_CERTIFICATE",
    "strict_only": True,
    "legacy_bridge_used": False,
    "pipeline": "K_strict -> coeff -> full_L_total -> linearized_EOM -> symbolic_bound_certificate",
    "curvature_domain": f"|R|<={R_max}",
    "admissible_domain": admissible_domain,
    "probe_rows": rows,
    "boundary_counterexamples": counterexamples,
    "status": "OPEN_OBLIGATION",
    "limitation": "Certificate scaffold only; theorem-level global completeness proof is still missing.",
    "next_honest_step": "S636: bridge certificate to global cocycle background-independence and full unitarity theorem statement."
}

OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
print(f"Wrote {OUT}")
