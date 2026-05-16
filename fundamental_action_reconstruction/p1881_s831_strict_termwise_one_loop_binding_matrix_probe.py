#!/usr/bin/env python3
"""P1881 S831 strict termwise one-loop binding matrix from full L_total registry."""
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
    p1880 = load("p1880_s830_strict_full_lagrangian_term_registry_to_eom_transport_probe.json")

    beta, eta, omega, phi = sp.symbols("beta eta omega phi", positive=True, real=True)
    m2, lam, y, xi, kappa2 = sp.symbols("m2 lam y xi kappa2", real=True)
    eps = sp.symbols("eps", positive=True, real=True)

    c0 = sp.cos(omega + phi) / (beta + 1)
    c1 = (-beta * eta * sp.cos(omega + phi) - omega * (beta + 1) * sp.sin(omega + phi)) / (beta + 1) ** 2
    m2_eff = m2 * (1 + c0)
    lam_eff = lam * (1 + c1**2)
    y_eff = y * (1 + c0 / 2)

    # termwise one-loop coefficients (symbolic placeholders, explicitly open)
    z_kin, z_m2, z_lam, z_y, z_F, z_xi, z_kappa = sp.symbols("z_kin z_m2 z_lam z_y z_F z_xi z_kappa", real=True)

    binding = {
        "L_scalar_kinetic": z_kin * lam_eff / (16 * sp.pi**2 * eps),
        "L_scalar_mass": z_m2 * m2_eff / (16 * sp.pi**2 * eps),
        "L_scalar_quartic": z_lam * lam_eff**2 / (16 * sp.pi**2 * eps),
        "L_yukawa": z_y * y_eff**2 / (16 * sp.pi**2 * eps),
        "L_gauge": z_F * y_eff**2 / (16 * sp.pi**2 * eps),
        "L_nonminimal": z_xi * xi * lam_eff / (16 * sp.pi**2 * eps),
        "L_gravity": z_kappa * m2_eff / (16 * sp.pi**2 * eps * kappa2),
    }

    qw2049 = {beta: 1.0, eta: 1.8, omega: 0.18575, phi: 0.16250}

    out = {
        "packet_id": "P1881",
        "stage_id": "S831",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "depends_on": {"p1880_present": "full_lagrangian_term_registry" in p1880},
        "strict_chain_step": "K_strict -> effective coefficients -> full L_total term registry -> termwise one-loop binding matrix",
        "effective_coefficients": {
            "m2_eff": str(m2_eff),
            "lam_eff": str(lam_eff),
            "y_eff": str(y_eff),
        },
        "termwise_one_loop_binding_matrix": {k: str(v) for k, v in binding.items()},
        "binding_interpretation": {
            "forward": "Each full-L term receives an explicit 1/eps binding slot in the same strict scheme.",
            "reverse": "Channel/transport witness failures map back to concrete term slots and strict-kernel corridors.",
        },
        "qw2049_trace": {
            "lam_eff_over_lam": str(sp.N((lam_eff / lam).subs(qw2049), 12)),
            "y_eff_over_y": str(sp.N((y_eff / y).subs(qw2049), 12)),
            "scalar_quartic_prefactor_qw2049": str(sp.N((binding["L_scalar_quartic"] * eps / z_lam / lam**2).subs(qw2049), 12)),
        },
        "strict_core_closure_missing_items": {
            "renormalization": "Need explicit loop integrals fixing z_* finite parts and consistency conditions.",
            "unitarity": "Need Cutkosky channel data to validate this termwise binding matrix.",
            "background_independence": "Need FRW/BianchiI transport of bound termwise counterterms.",
            "selector_obstruction": "QW-2191 remains open.",
        },
        "false_pass_guard": "Termwise binding matrix is a contract artifact, not theorem-level strict-core closure.",
        "next_honest_step": "Attach computed loop integrals to z_* and verify channel defects plus FRW<->BianchiI transport in one shared scheme.",
        "lay_explanation": "To rozpiska, który składnik pełnego lagranżianu dostaje jaką poprawkę kwantową 1-pętlową. Dzięki temu łatwiej sprawdzić, gdzie dokładnie teoria może się nie domykać.",
    }

    path = GEN / "p1881_s831_strict_termwise_one_loop_binding_matrix_probe.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
