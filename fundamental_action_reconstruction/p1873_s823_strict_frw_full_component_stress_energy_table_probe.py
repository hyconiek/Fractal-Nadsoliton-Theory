#!/usr/bin/env python3
"""P1873 S823 strict FRW full-component stress-energy table from full L_total lane."""
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
    m2, lam, xi, R, kappa2, Lambda_cc = sp.symbols("m2 lam xi R kappa2 Lambda_cc", real=True)
    H, phidot, phi0 = sp.symbols("H phidot phi0", real=True)
    E2, B2 = sp.symbols("E2 B2", nonnegative=True, real=True)
    psi_bar_psi, i_psi_dpsi = sp.symbols("psi_bar_psi i_psi_dpsi", real=True)

    c0 = sp.cos(omega + phi) / (beta + 1)
    c1 = (-beta * eta * sp.cos(omega + phi) - omega * (beta + 1) * sp.sin(omega + phi)) / (beta + 1) ** 2
    m2_eff = m2 * (1 + c0)
    lam_eff = lam * (1 + c1**2)

    Veff = sp.Rational(1, 2) * m2_eff * phi0**2 + lam_eff * phi0**4 / 4
    rho_scalar = sp.Rational(1, 2) * phidot**2 + Veff
    p_scalar = sp.Rational(1, 2) * phidot**2 - Veff

    rho_gauge = sp.Rational(1, 2) * (E2 + B2)
    p_gauge = sp.Rational(1, 6) * (E2 + B2)

    rho_fermion = i_psi_dpsi
    p_fermion = sp.Integer(0)

    rho_nonminimal = 3 * xi * H**2 * phi0**2
    p_nonminimal = -xi * (3 * H**2 + R) * phi0**2

    rho_total = rho_scalar + rho_gauge + rho_fermion + rho_nonminimal
    p_total = p_scalar + p_gauge + p_fermion + p_nonminimal

    residuals = {
        "R_g_00": 3 * H**2 + Lambda_cc - kappa2 * rho_total,
        "R_g_11": -3 * H**2 + Lambda_cc - kappa2 * p_total,
        "R_g_22": -3 * H**2 + Lambda_cc - kappa2 * p_total,
        "R_g_33": -3 * H**2 + Lambda_cc - kappa2 * p_total,
    }

    symbolic_closure = [{
        "H2": kappa2 * (rho_total + p_total) / 6,
        "Lambda_cc": kappa2 * (rho_total - p_total) / 2,
    }]

    qw2049 = {beta: 1.0, eta: 1.8, omega: 0.18575, phi: 0.16250}

    out = {
        "packet_id": "P1873",
        "stage_id": "S823",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "depends_on": {
            "p1866_present": "full_lagrangian_non_skeleton" in p1866,
            "p1872_present": "full_lagrangian_scalar_export" in p1872,
        },
        "strict_chain_step": "K_strict -> coefficients -> full L_total sector split -> T_mn components -> FRW Einstein residual table",
        "full_lagrangian_density_non_skeleton": {
            "L_scalar": "1/2*(d phi)^2 - V_eff(phi)",
            "L_gauge": "-1/4 * F_{mu nu} F^{mu nu}",
            "L_fermion": "i*psi_bar*gamma^mu*D_mu*psi - y*phi*psi_bar*psi",
            "L_nonminimal": "-1/2*xi*R*phi^2",
            "L_gravity": "(1/(2*kappa2))*(R - 2*Lambda_cc)",
            "note": "Sector decomposition exported as full-L_total working form (not skeleton-only).",
        },
        "stress_energy_components": {
            "rho_scalar": str(rho_scalar),
            "p_scalar": str(p_scalar),
            "rho_gauge": str(rho_gauge),
            "p_gauge": str(p_gauge),
            "rho_fermion": str(rho_fermion),
            "p_fermion": str(p_fermion),
            "rho_nonminimal": str(rho_nonminimal),
            "p_nonminimal": str(p_nonminimal),
            "rho_total": str(rho_total),
            "p_total": str(p_total),
        },
        "frw_component_residual_table": {k: str(v) for k, v in residuals.items()},
        "symbolic_closure_branch": [{k: str(v) for k, v in b.items()} for b in symbolic_closure],
        "qw2049_trace": {
            "rho_total_qw2049": str(sp.N(rho_total.subs(qw2049), 12)),
            "p_total_qw2049": str(sp.N(p_total.subs(qw2049), 12)),
        },
        "bidirectional_chain_note": {
            "forward": "Strict kernel fixes effective couplings, which feed all L_total sectors and then T_mn and Einstein residuals.",
            "reverse": "Component-level residual constraints feed back onto admissible strict-sector coefficient ranges and sector-weight consistency.",
        },
        "strict_core_closure_missing_exports": {
            "renormalization": "Missing loop-level counterterm map on same FRW branch for scalar+gauge+fermion+nonminimal sectors.",
            "unitarity": "Missing Cutkosky/optical witness with the same sector content around the FRW background.",
            "background_independence": "Missing atlas-level closure witness beyond one FRW foliation.",
            "selector_obstruction": "QW-2191 remains open; no strict selector closure is implied.",
        },
        "false_pass_guard": "Full component table export is not strict-core closure theorem.",
        "next_honest_step": "Export one-loop renormalized T_mn counterterm dictionary for these sectors and attach matching Cutkosky witness on same background branch.",
        "lay_explanation": "To kolejny krok: nie tylko pole skalarne, ale też pola cechowania i fermiony dokładają energię i ciśnienie. Liczymy ich wspólny wpływ na równania grawitacji i jawnie pokazujemy, czego jeszcze brakuje do pełnego domknięcia teorii.",
    }

    path = GEN / "p1873_s823_strict_frw_full_component_stress_energy_table_probe.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
