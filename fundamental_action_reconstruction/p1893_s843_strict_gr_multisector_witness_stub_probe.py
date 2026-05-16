#!/usr/bin/env python3
"""P1893 S843 strict renormalization-gate multisector witness stub probe."""
from __future__ import annotations
import json
from pathlib import Path
import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"


def load(name: str) -> dict:
    p = GEN / name
    if not p.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(p.read_text(encoding="utf-8"))


def main() -> None:
    p1892 = load("p1892_s842_strict_gr_first_witness_stub_probe.json")

    beta, eta, omega, phi = sp.symbols("beta eta omega phi", positive=True, real=True)
    lam, y, xi, kappa2 = sp.symbols("lam y xi kappa2", positive=True, real=True)
    eps, mu2, s = sp.symbols("eps mu2 s", positive=True, real=True)

    c0 = sp.cos(omega + phi) / (beta + 1)
    c1 = (-beta * eta * sp.cos(omega + phi) - omega * (beta + 1) * sp.sin(omega + phi)) / (beta + 1) ** 2
    lam_eff = lam * (1 + c1**2)
    y_eff = y * (1 + c0 / 2)

    A_div, A_fin = sp.symbols("A_div A_fin", real=True)
    B_div, B_fin = sp.symbols("B_div B_fin", real=True)
    C_div, C_fin = sp.symbols("C_div C_fin", real=True)

    I = 1 / (16 * sp.pi**2)

    # Sector stubs: scalar quartic, Yukawa, nonminimal/gravity-mixed slot.
    loop_scalar4 = lam_eff**2 * I * (A_div / eps + A_fin + sp.log(mu2 / s))
    ct_scalar4 = -lam_eff**2 * I * (A_div / eps + A_fin)

    loop_yukawa = y_eff**2 * I * (B_div / eps + B_fin + sp.log(mu2 / s))
    ct_yukawa = -y_eff**2 * I * (B_div / eps + B_fin)

    loop_nonminimal = xi * (y_eff**2 + kappa2) * I * (C_div / eps + C_fin + sp.log(mu2 / s))
    ct_nonminimal = -xi * (y_eff**2 + kappa2) * I * (C_div / eps + C_fin)

    ren_scalar4 = sp.expand(loop_scalar4 + ct_scalar4)
    ren_yukawa = sp.expand(loop_yukawa + ct_yukawa)
    ren_nonminimal = sp.expand(loop_nonminimal + ct_nonminimal)

    gate_stub_sum = sp.expand(ren_scalar4 + ren_yukawa + ren_nonminimal)
    pass_metric = sp.simplify(
        gate_stub_sum
        - I * sp.log(mu2 / s) * (lam_eff**2 + y_eff**2 + xi * (y_eff**2 + kappa2))
    )

    qw2049 = {beta: 1.0, eta: 1.8, omega: 0.18575, phi: 0.16250}

    out = {
        "packet_id": "P1893",
        "stage_id": "S843",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "depends_on": {"p1892_present": "first_gate_witness_stub" in p1892},
        "strict_chain_step": "K_strict -> coefficients -> full-L multisector one-loop stubs -> renormalization gate promotion candidate",
        "effective_coefficients": {
            "lam_eff": str(lam_eff),
            "y_eff": str(y_eff),
        },
        "multisector_witness_stub": {
            "scalar4": {"loop": str(loop_scalar4), "ct": str(ct_scalar4), "ren": str(ren_scalar4)},
            "yukawa": {"loop": str(loop_yukawa), "ct": str(ct_yukawa), "ren": str(ren_yukawa)},
            "nonminimal": {"loop": str(loop_nonminimal), "ct": str(ct_nonminimal), "ren": str(ren_nonminimal)},
            "gate_stub_sum": str(gate_stub_sum),
            "pass_metric_symbol": str(pass_metric),
        },
        "promotion_gate_target": "G_R",
        "promotion_readiness_condition": "pass_metric_symbol == 0 and sector-coupled scheme coherence",
        "qw2049_trace": {
            "lam_eff_over_lam": str(sp.N((lam_eff / lam).subs(qw2049), 12)),
            "y_eff_over_y": str(sp.N((y_eff / y).subs(qw2049), 12)),
        },
        "strict_core_closure_missing_items": {
            "renormalization": "Need diagram-fixed (A_*, B_*, C_*) and extension to remaining gauge/gravity channels.",
            "unitarity": "G_U still open.",
            "background_independence": "G_BI still open.",
            "selector_obstruction": "QW-2191 remains open.",
        },
        "false_pass_guard": "Multisector stub is still pre-witness and not full G_R closure theorem.",
        "next_honest_step": "Attach explicit diagram coefficients and check consistency of this multisector G_R stub with Cutkosky equations on same scheme.",
        "lay_explanation": "To rozszerzenie pierwszego testu renormalizacji na kilka sektorów naraz, żeby zbliżyć się do realnej promocji bramy G_R.",
    }

    path = GEN / "p1893_s843_strict_gr_multisector_witness_stub_probe.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
