#!/usr/bin/env python3
"""P1871 S821 strict FRW constitutive closure probe from strict coefficient lane."""
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
    p1870 = load("p1870_s820_strict_frw_background_metric_residual_probe.json")

    beta, eta, omega, phi = sp.symbols("beta eta omega phi", positive=True, real=True)
    m2, lam = sp.symbols("m2 lam", positive=True, real=True)
    H, kappa2, Lambda_cc = sp.symbols("H kappa2 Lambda_cc", real=True)

    c0 = sp.cos(omega + phi) / (beta + 1)
    c1 = (-beta * eta * sp.cos(omega + phi) - omega * (beta + 1) * sp.sin(omega + phi)) / (beta + 1) ** 2
    m2_eff = sp.simplify(m2 * (1 + c0))
    lam_eff = sp.simplify(lam * (1 + c1**2))

    # Minimal strict constitutive ansatz on FRW slice: homogeneous scalar condensate phi0.
    phi0 = sp.symbols("phi0", real=True)
    Veff = sp.simplify(sp.Rational(1, 2) * m2_eff * phi0**2 + lam_eff * phi0**4 / 4)
    rho_eff = Veff
    p_eff = -Veff  # slow-roll/static condensate limit

    R00 = sp.simplify(3 * H**2 + Lambda_cc - kappa2 * rho_eff)
    Rii = sp.simplify(-3 * H**2 + Lambda_cc - kappa2 * p_eff)

    closure_branch = sp.solve([sp.Eq(R00, 0), sp.Eq(Rii, 0)], [H**2, Lambda_cc], dict=True)

    qw2049 = {beta: 1.0, eta: 1.8, omega: 0.18575, phi: 0.16250}

    out = {
        "packet_id": "P1871",
        "stage_id": "S821",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "depends_on": {
            "p1866_present": "coefficient_export" in p1866,
            "p1870_present": "einstein_component_residual_probe" in p1870,
        },
        "strict_chain_step": "kernel strict -> effective couplings -> FRW constitutive map (rho,p) -> Einstein residual closure branch",
        "constitutive_ansatz": {
            "rho_eff": str(rho_eff),
            "p_eff": str(p_eff),
            "assumption": "homogeneous scalar condensate, static/slow-roll slice",
        },
        "frw_residuals_after_constitutive_injection": {
            "R_g_00": str(R00),
            "R_g_ii": str(Rii),
            "symbolic_closure_branches": [str(b) for b in closure_branch],
        },
        "qw2049_trace": {
            "rho_eff_factor_expr": str(sp.N((rho_eff / (m2 * phi0**2 / 2)).subs(qw2049), 12)),
            "p_eff_relation": "p_eff = -rho_eff",
        },
        "bidirectional_note": {
            "forward": "Strict kernel tuple fixes m2_eff/lambda_eff and thus the FRW energy-pressure source terms.",
            "reverse": "Residual-zero FRW branch constrains admissible cosmological sector parameters and feeds back to closure viability.",
        },
        "qg_obstruction_register": {
            "renormalization": "Need loop-stable constitutive map and counterterm closure beyond tree-level ansatz.",
            "unitarity": "Need Cutkosky/BRST witness compatibility for condensate background sector.",
            "background_independence": "Need atlas-level proof that closure branch is not gauge/background artifact.",
        },
        "false_pass_guard": "Existence of algebraic FRW closure branch under ansatz is not theorem-grade ToE closure.",
        "next_honest_step": "Replace ansatz with explicit stress-energy export from full L_SM+L_GR sector terms and recompute FRW component residual table.",
        "lay_explanation": "To pierwszy krok łączący kernel strict z kosmologią: z wyliczonych współczynników budujemy energię i ciśnienie, a potem sprawdzamy warunki równań grawitacji FRW.",
    }

    path = GEN / "p1871_s821_strict_frw_constitutive_closure_probe.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
