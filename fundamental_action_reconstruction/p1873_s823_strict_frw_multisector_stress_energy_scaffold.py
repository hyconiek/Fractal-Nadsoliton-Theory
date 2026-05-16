#!/usr/bin/env python3
"""P1873 S823 strict FRW multisector stress-energy scaffold from full L_total."""
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
    p1872 = load("p1872_s822_strict_frw_full_lagrangian_stress_energy_export_probe.json")

    beta, eta, omega, phi = sp.symbols("beta eta omega phi", positive=True, real=True)
    m2, lam, g = sp.symbols("m2 lam g", positive=True, real=True)
    H, kappa2, Lambda_cc = sp.symbols("H kappa2 Lambda_cc", real=True)
    phidot, phi0 = sp.symbols("phidot phi0", real=True)
    E2, B2 = sp.symbols("E2 B2", nonnegative=True, real=True)
    rho_f, p_f = sp.symbols("rho_f p_f", real=True)

    c0 = sp.cos(omega + phi) / (beta + 1)
    c1 = (-beta * eta * sp.cos(omega + phi) - omega * (beta + 1) * sp.sin(omega + phi)) / (beta + 1) ** 2
    c2 = (2 * beta * eta * omega * (beta + 1) * sp.sin(omega + phi)
          + beta * eta * (2 * beta * eta + (1 - eta) * (beta + 1)) * sp.cos(omega + phi)
          - omega ** 2 * (beta + 1) ** 2 * sp.cos(omega + phi)) / (2 * (beta + 1) ** 3)

    m2_eff = sp.simplify(m2 * (1 + c0))
    lam_eff = sp.simplify(lam * (1 + c1**2))
    g_eff = sp.simplify(g * (1 + c2))

    Veff = sp.Rational(1, 2) * m2_eff * phi0**2 + lam_eff * phi0**4 / 4
    rho_scalar = sp.simplify(sp.Rational(1, 2) * phidot**2 + Veff)
    p_scalar = sp.simplify(sp.Rational(1, 2) * phidot**2 - Veff)

    # FRW-averaged gauge radiation-like sector (seed form).
    rho_gauge = sp.simplify(g_eff * (E2 + B2) / 2)
    p_gauge = sp.simplify(rho_gauge / 3)

    rho_total = rho_scalar + rho_gauge + rho_f
    p_total = p_scalar + p_gauge + p_f

    R00 = 3 * H**2 + Lambda_cc - kappa2 * rho_total
    Rii = -3 * H**2 + Lambda_cc - kappa2 * p_total

    out = {
        "packet_id": "P1873",
        "stage_id": "S823",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "depends_on": {
            "p1866_present": "full_lagrangian_non_skeleton" in p1866,
            "p1872_present": "full_lagrangian_scalar_export" in p1872,
        },
        "strict_chain_step": "K_strict -> effective couplings -> multisector T_mn seed -> FRW residual equations",
        "multisector_stress_energy_seed": {
            "rho_scalar": str(rho_scalar),
            "p_scalar": str(p_scalar),
            "rho_gauge_seed": str(rho_gauge),
            "p_gauge_seed": str(p_gauge),
            "rho_fermion_seed": str(rho_f),
            "p_fermion_seed": str(p_f),
            "rho_total": str(rho_total),
            "p_total": str(p_total),
        },
        "frw_component_residuals": {
            "R_g_00": str(R00),
            "R_g_ii": str(Rii),
            "closure_requirements": {
                "friedmann": "3*H**2 + Lambda_cc = kappa2*rho_total",
                "pressure": "-3*H**2 + Lambda_cc = kappa2*p_total",
            },
        },
        "bidirectional_rule": {
            "forward": "Kernel tuple changes m2_eff/lambda_eff/g_eff and propagates into rho_total,p_total and FRW residuals.",
            "reverse": "Any accepted residual-zero branch constrains admissible multisector energy budget and therefore strict parameter viability.",
        },
        "qg_gap": {
            "renormalization": "Need loop-corrected multisector T_mn closure in one renormalization scheme.",
            "unitarity": "Need BRST/Cutkosky consistency for FRW perturbative channels.",
            "background_independence": "Need atlas continuation beyond single FRW class.",
        },
        "false_pass_guard": "Multisector stress-energy seed scaffold is not theorem-level closure.",
        "next_honest_step": "Instantiate rho_f,p_f and gauge E2/B2 from explicit mode sums/operator expectation exports and recompute R_g_00..R_g_33 table.",
        "lay_explanation": "To kolejny krok realizmu: do energii i ciśnienia dokładamy sektory cechowania i fermionów, aby test grawitacyjny obejmował pełniejszy obraz fizyki.",
    }

    path = GEN / "p1873_s823_strict_frw_multisector_stress_energy_scaffold.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
