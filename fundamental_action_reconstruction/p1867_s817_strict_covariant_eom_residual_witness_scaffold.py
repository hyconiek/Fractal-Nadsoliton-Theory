#!/usr/bin/env python3
"""P1867 S817 strict covariant EOM residual witness scaffold (strict-only, no false pass)."""
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
    p1866 = load("p1866_s816_strict_kernel_to_full_lagrangian_and_eom_symbolic_export.json")

    beta, eta, omega, phi = sp.symbols("beta eta omega phi", positive=True, real=True)
    m2, lam, y, g = sp.symbols("m2 lam y g", positive=True, real=True)
    x = sp.symbols("x", real=True)
    phix, Ax, psix, psibx = (sp.Function("phi")(x), sp.Function("A")(x), sp.Function("psi")(x), sp.Function("psib")(x))

    c0 = sp.cos(omega + phi) / (beta + 1)
    c1 = (-beta * eta * sp.cos(omega + phi) - omega * (beta + 1) * sp.sin(omega + phi)) / (beta + 1) ** 2
    c2 = (2 * beta * eta * omega * (beta + 1) * sp.sin(omega + phi)
          + beta * eta * (2 * beta * eta + (1 - eta) * (beta + 1)) * sp.cos(omega + phi)
          - omega ** 2 * (beta + 1) ** 2 * sp.cos(omega + phi)) / (2 * (beta + 1) ** 3)

    m2_eff = m2 * (1 + c0)
    lam_eff = lam * (1 + c1**2)
    y_eff = y * (1 + c0 / 2)
    g_eff = g * (1 + c2)

    # Covariant-shape local equations (flat proxy with explicit source structure).
    eom_phi_shape = (
        sp.Symbol("Box_g_phi")
        + m2_eff * sp.Symbol("phi")
        + lam_eff * sp.Symbol("phi")**3
        + y_eff * sp.Symbol("psibpsi")
        - g_eff * sp.Symbol("A2") * sp.Symbol("phi")
    )
    eom_A_shape = sp.Symbol("Nabla_mu_F_mu_nu") - g_eff * sp.Symbol("A_nu") * sp.Symbol("phi2")

    # 1D executable residual from same chain for regression visibility.
    dphi, dA = sp.diff(phix, x), sp.diff(Ax, x)
    L = (sp.Rational(1, 2) * dphi**2 - sp.Rational(1, 2) * m2_eff * phix**2 - lam_eff * phix**4 / 4
         + sp.I * psibx * sp.diff(psix, x) - y_eff * phix * psibx * psix
         - sp.Rational(1, 2) * dA**2 + g_eff * Ax**2 * phix**2 / 2)
    eom_phi_1d = sp.simplify(sp.diff(sp.diff(L, dphi), x) - sp.diff(L, phix))
    eom_A_1d = sp.simplify(sp.diff(sp.diff(L, dA), x) - sp.diff(L, Ax))

    tuple_qw2049 = {beta: 1.0, eta: 1.8, omega: 0.18575, phi: 0.16250}

    out = {
        "packet_id": "P1867",
        "stage_id": "S817",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "depends_on": {"p1866_present": "coefficient_export" in p1866},
        "strict_chain_scope": "K_strict -> coefficients -> L_SM+L_GR -> EOM residual witnesses",
        "covariant_eom_shape_export": {
            "E_phi_shape": str(sp.expand(eom_phi_shape)),
            "E_A_shape": str(sp.expand(eom_A_shape)),
            "note": "Shape-level covariant operators exported; theorem-grade tensor evaluation remains open.",
        },
        "residual_regression_proxy_1d": {
            "E_phi_1d_symbolic": str(sp.expand(eom_phi_1d)),
            "E_A_1d_symbolic": str(sp.expand(eom_A_1d)),
            "qw2049_parameter_injection": {
                "E_phi_1d_numeric_coeff_form": str(sp.N(eom_phi_1d.subs(tuple_qw2049), 12)),
                "E_A_1d_numeric_coeff_form": str(sp.N(eom_A_1d.subs(tuple_qw2049), 12)),
            },
        },
        "bidirectional_rule": {
            "forward": "Kernel tuple fixes effective couplings and changes every EOM coefficient.",
            "reverse": "Any claimed EOM closure must be traceable back to admissible kernel tuple constraints.",
        },
        "strict_core_closure_blockers": [
            "full 4D componentwise covariant E_phi/E_A/E_g residual-zero tables not yet exported",
            "counterterm cancellation theorem on full operator basis not yet exported",
            "Cutkosky discontinuity theorem witness not yet exported",
            "background-independence lift theorem not yet exported",
        ],
        "false_pass_guard": "Residual scaffold and numeric coefficient forms are not theorem discharge.",
        "next_honest_step": "Export 4D componentwise residual tables for E_phi, E_A, and Einstein equations under one shared renormalization scheme.",
        "lay_explanation": "To jest krok z planu: mamy już wzory mówiące, jak z kernela strict dostaje się konkretne równania ruchu; teraz trzeba policzyć pełne 4D tabele zer reszt, żeby zrobić prawdziwy dowód.",
    }

    path = GEN / "p1867_s817_strict_covariant_eom_residual_witness_scaffold.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
