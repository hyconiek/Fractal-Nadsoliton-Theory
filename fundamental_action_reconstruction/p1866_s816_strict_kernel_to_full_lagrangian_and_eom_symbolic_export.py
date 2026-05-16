#!/usr/bin/env python3
"""P1866 S816 strict kernel->coefficients->full Lagrangian->EOM symbolic export."""
from __future__ import annotations
import json
from pathlib import Path
import sympy as sp
ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"

def main() -> None:
    d, beta, eta, omega, phi = sp.symbols("d beta eta omega phi", positive=True, real=True)
    m2, lam, y, g, kappa = sp.symbols("m2 lam y g kappa", positive=True, real=True)
    xi, Lambda_cc = sp.symbols("xi Lambda_cc", real=True)
    x = sp.symbols("x", real=True)
    k_strict = sp.cos(omega * d + phi) / (1 + beta * d**eta)

    # Moments evaluated at finite reference scale d=1 to avoid non-analytic d->0 for generic eta.
    d_ref = sp.Integer(1)
    c0 = sp.simplify(k_strict.subs(d, d_ref))
    c1 = sp.simplify(sp.diff(k_strict, d).subs(d, d_ref))
    c2 = sp.simplify(sp.diff(k_strict, d, 2).subs(d, d_ref) / 2)

    m2_eff = sp.simplify(m2 * (1 + c0))
    lam_eff = sp.simplify(lam * (1 + c1**2))
    y_eff = sp.simplify(y * (1 + c0 / 2))
    g_eff = sp.simplify(g * (1 + c2))
    xi_eff = sp.simplify(xi * (1 + c0))

    phi_f, psi, psib, A = (sp.Function("phi")(x), sp.Function("psi")(x), sp.Function("psib")(x), sp.Function("A")(x))
    dphi, dA = sp.diff(phi_f, x), sp.diff(A, x)
    L_scalar = sp.Rational(1, 2) * dphi**2 - sp.Rational(1, 2) * m2_eff * phi_f**2 - lam_eff * phi_f**4 / 4
    L_fermion = sp.I * psib * sp.diff(psi, x) - y_eff * phi_f * psib * psi
    L_gauge = -sp.Rational(1, 2) * dA**2 + g_eff * A**2 * phi_f**2 / 2
    R = sp.Symbol("R", real=True)
    sqrt_minus_g = sp.Symbol("sqrt_minus_g", positive=True)
    L_gravity = sqrt_minus_g * ((R - 2 * Lambda_cc) / (2 * kappa**2) + xi_eff * phi_f**2 * R)

    L_matter = L_scalar + L_fermion + L_gauge
    eom_phi = sp.simplify(sp.diff(sp.diff(L_matter, dphi), x) - sp.diff(L_matter, phi_f))
    eom_A = sp.simplify(sp.diff(sp.diff(L_matter, dA), x) - sp.diff(L_matter, A))

    out = {"packet_id":"P1866","stage_id":"S816","status":"OPEN_OBSTRUCTION_WITH_TRACE","route":"strict_only",
           "strict_kernel":{"definition":"K_strict(d)=cos(omega*d+phi)/(1+beta*d^eta)","symbolic":str(k_strict),"reference_scale":"d_ref=1"},
           "coefficient_export":{"kernel_moments_at_d_ref":{"c0":str(c0),"c1":str(c1),"c2":str(c2)},"strict_effective_couplings":{"m2_eff":str(m2_eff),"lam_eff":str(lam_eff),"y_eff":str(y_eff),"g_eff":str(g_eff),"xi_eff":str(xi_eff)}},
           "full_lagrangian_non_skeleton":{"L_total_decomposition":{"L_scalar":str(sp.expand(L_scalar)),"L_fermion":str(sp.expand(L_fermion)),"L_gauge":str(sp.expand(L_gauge)),"L_gravity_density":str(L_gravity),"L_SM_plus_GR_statement":"L_total = L_scalar + L_fermion + L_gauge + sqrt(-g)*[(R-2*Lambda)/(2*kappa^2) + xi_eff*phi^2*R]"},"directionality":{"forward":"K_strict -> (c0,c1,c2) -> effective couplings -> explicit L_SM+L_GR","reverse_requirement":"EOM/QG witness residuals must constrain admissible kernel tuple"}},
           "eom_export":{"eom_phi_proxy_1d":str(sp.expand(eom_phi)),"eom_A_proxy_1d":str(sp.expand(eom_A))},
           "strict_core_closure_blockers":["4D covariant Euler-Lagrange+Einstein witness tables still missing","full counterterm cancellation proof on operator basis still missing","exact Cutkosky discontinuity integrals still missing","background-independence lift theorem still missing"],
           "false_pass_guard":"This symbolic chain export is not TG2/TG3/ToE closure.",
           "next_honest_step":"Compute 4D covariant EOM tensors from this L_total and export residual-zero certificates term by term.",
           "lay_explanation":"To jest pełny, policzalny tor: kernel strict daje współczynniki, one budują pełny lagranżian SM+GR, a z niego liczymy równania ruchu; pełne domknięcie ToE wciąż wymaga twardych dowodów QG."}
    path = GEN / "p1866_s816_strict_kernel_to_full_lagrangian_and_eom_symbolic_export.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False)+"\n", encoding="utf-8")
    print(path)
if __name__ == "__main__":
    main()
