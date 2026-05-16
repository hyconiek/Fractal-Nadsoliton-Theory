#!/usr/bin/env python3
"""P1872 S822 strict FRW stress-energy export probe from full Lagrangian lane."""
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
    p1871 = load("p1871_s821_strict_frw_constitutive_closure_probe.json")

    beta, eta, omega, phi = sp.symbols("beta eta omega phi", positive=True, real=True)
    m2, lam = sp.symbols("m2 lam", positive=True, real=True)
    H, kappa2, Lambda_cc = sp.symbols("H kappa2 Lambda_cc", real=True)
    phidot, phi0 = sp.symbols("phidot phi0", real=True)

    c0 = sp.cos(omega + phi) / (beta + 1)
    c1 = (-beta * eta * sp.cos(omega + phi) - omega * (beta + 1) * sp.sin(omega + phi)) / (beta + 1) ** 2
    m2_eff = sp.simplify(m2 * (1 + c0))
    lam_eff = sp.simplify(lam * (1 + c1**2))

    Veff = sp.Rational(1, 2) * m2_eff * phi0**2 + lam_eff * phi0**4 / 4
    rho_scalar = sp.simplify(sp.Rational(1, 2) * phidot**2 + Veff)
    p_scalar = sp.simplify(sp.Rational(1, 2) * phidot**2 - Veff)

    R00 = sp.simplify(3 * H**2 + Lambda_cc - kappa2 * rho_scalar)
    Rii = sp.simplify(-3 * H**2 + Lambda_cc - kappa2 * p_scalar)

    solved = sp.solve([sp.Eq(R00, 0), sp.Eq(Rii, 0)], [H**2, Lambda_cc], dict=True)
    qw2049 = {beta: 1.0, eta: 1.8, omega: 0.18575, phi: 0.16250}

    out = {
        "packet_id": "P1872",
        "stage_id": "S822",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "depends_on": {
            "p1866_present": "full_lagrangian_non_skeleton" in p1866,
            "p1871_present": "constitutive_ansatz" in p1871,
        },
        "strict_chain_step": "K_strict -> couplings -> full-L scalar sector stress-energy -> FRW Einstein residuals",
        "full_lagrangian_scalar_export": {
            "rho_scalar": str(rho_scalar),
            "p_scalar": str(p_scalar),
            "note": "Derived from canonical scalar kinetic+potential slice compatible with full-L lane",
        },
        "frw_residuals": {
            "R_g_00": str(R00),
            "R_g_ii": str(Rii),
            "closure_branch": [str(s) for s in solved],
        },
        "qw2049_trace": {
            "rho_scalar_qw2049": str(sp.N(rho_scalar.subs(qw2049), 12)),
            "p_scalar_qw2049": str(sp.N(p_scalar.subs(qw2049), 12)),
        },
        "bidirectional_note": {
            "forward": "Kernel strict sets effective potential terms entering T_mn and Einstein residuals.",
            "reverse": "FRW closure constraints restrict admissible strict-sector energy budget and feed back to parameter viability.",
        },
        "qg_gap": {
            "renormalization": "Need loop-corrected T_mn and counterterm consistency on same background family.",
            "unitarity": "Need optical/Cutkosky consistency for perturbations around FRW branch.",
            "background_independence": "Need atlas proof beyond single FRW slicing.",
        },
        "false_pass_guard": "Stress-energy export on FRW slice is not full theorem closure.",
        "next_honest_step": "Add gauge+fermion+nonminimal contributions to T_mn from full L_total and recompute full component residual table R_g_00..R_g_33.",
        "lay_explanation": "Tu robimy bardziej fizyczny krok: energia i ciśnienie wynikają już z jawnego składnika lagranżianu, więc test równań kosmologii jest bliżej pełnej teorii.",
    }

    path = GEN / "p1872_s822_strict_frw_full_lagrangian_stress_energy_export_probe.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
