#!/usr/bin/env python3
"""P1880 S830 strict full-Lagrangian term-registry -> EOM -> transport contract probe."""
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
    p1879 = load("p1879_s829_strict_frw_to_bianchiI_transport_contract_probe.json")

    beta, eta, omega, phi = sp.symbols("beta eta omega phi", positive=True, real=True)
    m2, lam, y, xi, kappa2, Lambda_cc = sp.symbols("m2 lam y xi kappa2 Lambda_cc", real=True)

    c0 = sp.cos(omega + phi) / (beta + 1)
    c1 = (-beta * eta * sp.cos(omega + phi) - omega * (beta + 1) * sp.sin(omega + phi)) / (beta + 1) ** 2

    m2_eff = m2 * (1 + c0)
    lam_eff = lam * (1 + c1**2)
    y_eff = y * (1 + c0 / 2)

    term_registry = {
        "L_gravity": "sqrt(-g) * (R - 2*Lambda_cc)/(2*kappa2)",
        "L_scalar_kinetic": "sqrt(-g) * (1/2*g^{mu nu}*d_mu(phi)*d_nu(phi))",
        "L_scalar_mass": "sqrt(-g) * (-1/2*m2_eff*phi^2)",
        "L_scalar_quartic": "sqrt(-g) * (-lam_eff/4*phi^4)",
        "L_yukawa": "sqrt(-g) * (-y_eff*phi*psi_bar*psi)",
        "L_fermion_kinetic": "sqrt(-g) * (i*psi_bar*gamma^mu*D_mu*psi)",
        "L_gauge": "sqrt(-g) * (-1/4*F^a_{mu nu}F_a^{mu nu})",
        "L_nonminimal": "sqrt(-g) * (-1/2*xi*R*phi^2)",
    }

    eom_registry = {
        "EOM_phi": "Box(phi) + m2_eff*phi + lam_eff*phi^3 + y_eff*psi_bar*psi + xi*R*phi = 0",
        "EOM_psi": "i*gamma^mu*D_mu*psi - y_eff*phi*psi = 0",
        "EOM_A": "D_mu F^{a mu nu} - J^{a nu}_matter = 0",
        "EOM_g": "G_{mu nu} + Lambda_cc*g_{mu nu} - kappa2*T_{mu nu}^{total} = 0",
    }

    qw2049 = {beta: 1.0, eta: 1.8, omega: 0.18575, phi: 0.16250}

    out = {
        "packet_id": "P1880",
        "stage_id": "S830",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "depends_on": {
            "p1866_present": "full_lagrangian_non_skeleton" in p1866,
            "p1879_present": "bianchiI_transport_ansatz" in p1879,
        },
        "strict_chain_forward": "K_strict -> (m2_eff,lam_eff,y_eff) -> full L_total term registry -> covariant EOM registry -> FRW/BianchiI transport contracts",
        "strict_chain_reverse": "transport+unitarity+renormalization witnesses -> EOM consistency -> full L_total integrity -> admissible strict-kernel corridor",
        "effective_coefficients": {
            "m2_eff": str(m2_eff),
            "lam_eff": str(lam_eff),
            "y_eff": str(y_eff),
        },
        "full_lagrangian_term_registry": term_registry,
        "covariant_eom_registry": eom_registry,
        "qw2049_trace": {
            "m2_eff_over_m2": str(sp.N((m2_eff / m2).subs(qw2049), 12)),
            "lam_eff_over_lam": str(sp.N((lam_eff / lam).subs(qw2049), 12)),
            "y_eff_over_y": str(sp.N((y_eff / y).subs(qw2049), 12)),
        },
        "strict_core_closure_missing_items": {
            "renormalization": "Need explicit 1-loop finite-part-fixed counterterm solution attached to every registry term.",
            "unitarity": "Need channel ImM computations closing optical defects with this exact EOM/term registry.",
            "background_independence": "Need theorem-grade FRW<->BianchiI transport witness for the same registry.",
            "selector_obstruction": "QW-2191 remains open; selector source/theorem still missing.",
        },
        "false_pass_guard": "Term registry + EOM registry are explicit but do not constitute strict-core closure theorem.",
        "next_honest_step": "Attach per-term one-loop integrals and prove defect transport consistency between FRW and BianchiI on identical renormalization data.",
        "lay_explanation": "To katalog pełnego lagranżianu i równań ruchu w wersji strict. Dzięki temu wiadomo dokładnie, które składniki muszą przejść testy kwantowe i transport między tłami.",
    }

    path = GEN / "p1880_s830_strict_full_lagrangian_term_registry_to_eom_transport_probe.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
