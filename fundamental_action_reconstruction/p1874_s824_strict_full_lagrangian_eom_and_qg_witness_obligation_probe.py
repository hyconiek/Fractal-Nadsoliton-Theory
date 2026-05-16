#!/usr/bin/env python3
"""P1874 S824 strict full-Lagrangian EOM + QG witness obligation export probe."""
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
    p1873 = load("p1873_s823_strict_frw_full_component_stress_energy_table_probe.json")

    beta, eta, omega, phi = sp.symbols("beta eta omega phi", positive=True, real=True)
    m2, lam, y, xi = sp.symbols("m2 lam y xi", real=True)
    kappa2, Lambda_cc = sp.symbols("kappa2 Lambda_cc", real=True)

    # Strict-kernel driven effective couplings.
    c0 = sp.cos(omega + phi) / (beta + 1)
    c1 = (-beta * eta * sp.cos(omega + phi) - omega * (beta + 1) * sp.sin(omega + phi)) / (beta + 1) ** 2
    m2_eff = sp.simplify(m2 * (1 + c0))
    lam_eff = sp.simplify(lam * (1 + c1**2))
    y_eff = sp.simplify(y * (1 + c0 / 2))

    L_total = {
        "L_gravity": "sqrt(-g) * ((R - 2*Lambda_cc)/(2*kappa2))",
        "L_scalar": "sqrt(-g) * (1/2*g^{mu nu}*d_mu(phi)*d_nu(phi) - 1/2*m2_eff*phi^2 - lam_eff/4*phi^4)",
        "L_gauge": "sqrt(-g) * (-1/4*F^a_{mu nu}*F_a^{mu nu})",
        "L_fermion": "sqrt(-g) * (i*psi_bar*gamma^mu*D_mu*psi - y_eff*phi*psi_bar*psi)",
        "L_nonminimal": "sqrt(-g) * (-1/2*xi*R*phi^2)",
    }

    eom = {
        "phi": "Box(phi) + m2_eff*phi + lam_eff*phi^3 + y_eff*psi_bar*psi + xi*R*phi = 0",
        "gauge": "D_mu F^{a mu nu} - J^{a nu}_fermion = 0",
        "fermion": "i*gamma^mu*D_mu*psi - y_eff*phi*psi = 0",
        "einstein": "G_{mu nu} + Lambda_cc*g_{mu nu} - kappa2*(T_scalar + T_gauge + T_fermion + T_nonminimal)_{mu nu} = 0",
    }

    qg_obligations = {
        "renormalization_export": {
            "required_object": "One-loop counterterm dictionary DeltaZ_{phi,F,psi,xi,kappa,Lambda} on the same FRW branch",
            "status": "OPEN_MISSING_EXPORT",
        },
        "unitarity_export": {
            "required_object": "Cutkosky/optical witness with the same effective couplings and field content",
            "status": "OPEN_MISSING_WITNESS",
        },
        "background_independence_export": {
            "required_object": "Atlas-level theorem transporting closure beyond one FRW slicing",
            "status": "OPEN_MISSING_THEOREM",
        },
    }

    qw2049 = {beta: 1.0, eta: 1.8, omega: 0.18575, phi: 0.16250}

    out = {
        "packet_id": "P1874",
        "stage_id": "S824",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "depends_on": {
            "p1866_present": "full_lagrangian_non_skeleton" in p1866,
            "p1873_present": "stress_energy_components" in p1873,
        },
        "strict_chain": {
            "kernel": "K_strict(d)=cos(omega*d+phi)/(1+beta*d^eta)",
            "coefficients": {
                "m2_eff": str(m2_eff),
                "lam_eff": str(lam_eff),
                "y_eff": str(y_eff),
            },
            "lagrangian_full_non_skeleton": L_total,
            "equations_of_motion": eom,
        },
        "bidirectional_theory_chain": {
            "forward": "K_strict fixes effective couplings (m2_eff, lam_eff, y_eff), then full L_total, then EOM and Einstein residuals.",
            "reverse": "Consistency of EOM/residuals constrains admissible strict-kernel parameter windows and effective-coupling sectors.",
        },
        "qw2049_coefficient_trace": {
            "m2_eff_qw2049": str(sp.N(m2_eff.subs(qw2049), 12)),
            "lam_eff_qw2049_over_lam": str(sp.N((lam_eff / lam).subs(qw2049), 12)),
            "y_eff_qw2049_over_y": str(sp.N((y_eff / y).subs(qw2049), 12)),
        },
        "strict_core_closure_obligations": qg_obligations,
        "selector_obstruction": "QW-2191 remains open; no strict selector closure implied.",
        "false_pass_guard": "Full strict chain export is still not ToE closure without missing renormalization/unitarity/background-independence witnesses/theorems.",
        "next_honest_step": "Export explicit one-loop counterterm equations (per sector) and bind them to a same-branch Cutkosky witness table.",
        "lay_explanation": "Mamy teraz pełniejszą mapę: od kernela strict przechodzimy przez współczynniki do pełnego lagranżianu i równań ruchu, ale do finalnej teorii nadal trzeba dowodów kwantowych (renormalizacja, unitarność, niezależność od tła).",
    }

    path = GEN / "p1874_s824_strict_full_lagrangian_eom_and_qg_witness_obligation_probe.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
