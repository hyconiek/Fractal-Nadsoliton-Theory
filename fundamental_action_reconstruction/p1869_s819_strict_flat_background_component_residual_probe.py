#!/usr/bin/env python3
"""P1869 S819 strict flat-background component residual probe (first concrete fill)."""
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
    p1868 = load("p1868_s818_strict_4d_componentwise_residual_table_scaffold.json")

    beta, eta, omega, phi = sp.symbols("beta eta omega phi", positive=True, real=True)
    m2, lam, y, g = sp.symbols("m2 lam y g", positive=True, real=True)

    c0 = sp.cos(omega + phi) / (beta + 1)
    c1 = (-beta * eta * sp.cos(omega + phi) - omega * (beta + 1) * sp.sin(omega + phi)) / (beta + 1) ** 2
    c2 = (2 * beta * eta * omega * (beta + 1) * sp.sin(omega + phi)
          + beta * eta * (2 * beta * eta + (1 - eta) * (beta + 1)) * sp.cos(omega + phi)
          - omega ** 2 * (beta + 1) ** 2 * sp.cos(omega + phi)) / (2 * (beta + 1) ** 3)

    m2_eff = sp.simplify(m2 * (1 + c0))
    lam_eff = sp.simplify(lam * (1 + c1**2))
    y_eff = sp.simplify(y * (1 + c0 / 2))
    g_eff = sp.simplify(g * (1 + c2))

    # Flat/vacuum probe assumptions: eta_{mu nu}, A_mu=0, psi=0, phi=0, no matter stress.
    # Keep cosmological term explicit as obstruction channel.
    Lambda_cc, kappa2 = sp.symbols("Lambda_cc kappa2", real=True)

    E_phi = sp.Symbol("Box_g_phi") + m2_eff * sp.Symbol("phi") + lam_eff * sp.Symbol("phi")**3 + y_eff * sp.Symbol("psibpsi") - g_eff * sp.Symbol("A2") * sp.Symbol("phi")
    E_A_mu = sp.Symbol("Nabla_nu_F_nu_mu") - g_eff * sp.Symbol("A_mu") * sp.Symbol("phi2")

    subs_vac = {
        sp.Symbol("Box_g_phi"): 0,
        sp.Symbol("phi"): 0,
        sp.Symbol("psibpsi"): 0,
        sp.Symbol("A2"): 0,
        sp.Symbol("Nabla_nu_F_nu_mu"): 0,
        sp.Symbol("A_mu"): 0,
        sp.Symbol("phi2"): 0,
    }

    E_phi_vac = sp.simplify(E_phi.subs(subs_vac))
    E_A_vac = sp.simplify(E_A_mu.subs(subs_vac))

    # Einstein residual under flat background: G_mn=0 and T_mn=0 leaves cosmological obstruction.
    R_g_mn_flat = -Lambda_cc * sp.Symbol("eta_mn")

    qw2049 = {beta: 1.0, eta: 1.8, omega: 0.18575, phi: 0.16250}

    out = {
        "packet_id": "P1869",
        "stage_id": "S819",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "depends_on": {"p1868_present": "componentwise_residual_table_scaffold" in p1868},
        "declared_background_family": "flat_minkowski_vacuum_probe",
        "assumptions": [
            "g_mn = eta_mn",
            "A_mu = 0",
            "psi = 0",
            "phi = 0",
            "T_mn_scalar = T_mn_fermion = T_mn_gauge = T_mn_nonminimal = 0",
        ],
        "residual_fill_result": {
            "E_phi_under_assumptions": str(E_phi_vac),
            "E_A_mu_under_assumptions": str(E_A_vac),
            "R_g_mn_under_assumptions": str(R_g_mn_flat),
            "component_status": {
                "scalar": "PASS_ZERO_LOCAL_BACKGROUND",
                "gauge": "PASS_ZERO_LOCAL_BACKGROUND",
                "metric": "OPEN_OBSTRUCTION_WITH_TRACE",
            },
            "metric_obstruction_reason": "Flat vacuum requires Lambda_cc=0; otherwise Einstein residual remains nonzero.",
        },
        "coefficient_trace_qw2049": {
            "m2_eff_factor": str(sp.N((1 + c0).subs(qw2049), 12)),
            "lam_eff_factor": str(sp.N((1 + c1**2).subs(qw2049), 12)),
            "y_eff_factor": str(sp.N((1 + c0 / 2).subs(qw2049), 12)),
            "g_eff_factor": str(sp.N((1 + c2).subs(qw2049), 12)),
        },
        "bidirectional_note": {
            "forward": "Kernel tuple fixes effective-coupling factors used in residual equations.",
            "reverse": "Background residual conditions (e.g., Lambda_cc constraint) restrict admissible closure claims.",
        },
        "false_pass_guard": "Local flat-vacuum partial zeros do not imply global TG2/TG3/ToE closure.",
        "next_honest_step": "Repeat component fill for nontrivial curved background atlas and export explicit R_g_00..R_g_33 numeric/symbolic tables with obstruction/zero verdicts.",
        "lay_explanation": "W najprostszym tle (płaska próżnia) równania dla pola skalarnego i cechowania się zerują, ale część grawitacyjna nadal blokuje domknięcie, jeśli stała kosmologiczna nie spełni warunku.",
    }

    path = GEN / "p1869_s819_strict_flat_background_component_residual_probe.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
