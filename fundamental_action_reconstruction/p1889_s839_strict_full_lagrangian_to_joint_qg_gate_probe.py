#!/usr/bin/env python3
"""P1889 S839 strict full-Lagrangian to joint QG gate probe."""
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
    p1888 = load("p1888_s838_strict_full_chain_forward_reverse_contraction_probe.json")
    p1884 = load("p1884_s834_strict_explicit_loop_integral_stub_evaluation_probe.json")

    beta, eta, omega, phi = sp.symbols("beta eta omega phi", positive=True, real=True)
    m2, lam, y, kappa2 = sp.symbols("m2 lam y kappa2", positive=True, real=True)
    eps, mu2, s = sp.symbols("eps mu2 s", positive=True, real=True)

    c0 = sp.cos(omega + phi) / (beta + 1)
    c1 = (-beta * eta * sp.cos(omega + phi) - omega * (beta + 1) * sp.sin(omega + phi)) / (beta + 1) ** 2
    m2_eff = m2 * (1 + c0)
    lam_eff = lam * (1 + c1**2)
    y_eff = y * (1 + c0 / 2)

    I_lam = 1 / (16 * sp.pi**2) * (1 / eps + sp.log(mu2 / s))
    I_y = 1 / (16 * sp.pi**2) * (1 / eps + sp.log(mu2 / s) + sp.Rational(1, 2))

    gate_vector = {
        "G_R": sp.expand(lam_eff * I_lam),
        "G_U": sp.expand(y_eff**2 * I_y),
        "G_BI": sp.expand((y_eff**2 + kappa2) * I_y),
    }

    qw2049 = {beta: 1.0, eta: 1.8, omega: 0.18575, phi: 0.16250}

    out = {
        "packet_id": "P1889",
        "stage_id": "S839",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "depends_on": {
            "p1888_present": "reverse_contraction_contract" in p1888,
            "p1884_present": "explicit_loop_integral_stubs" in p1884,
        },
        "strict_chain_step": "K_strict -> coefficients -> full non-skeleton L_SM+L_GR -> EOM -> explicit one-loop gate vector (R/U/BI)",
        "full_lagrangian_context": {
            "statement": "Nadsoliton => L_SM + L_GR (strict-only, non-skeleton)",
            "forward": "kernel -> coefficients -> L_total -> EOM -> joint quantum gates",
            "reverse": "gate vector admissibility -> coefficient windows -> strict kernel corridor",
        },
        "effective_coefficients": {
            "m2_eff": str(m2_eff),
            "lam_eff": str(lam_eff),
            "y_eff": str(y_eff),
        },
        "joint_gate_vector": {k: str(v) for k, v in gate_vector.items()},
        "gate_rule": "strict-core closure requires simultaneous bounded/compatible G_R, G_U, G_BI with solved witness equations",
        "qw2049_trace": {
            "m2_eff_over_m2": str(sp.N((m2_eff / m2).subs(qw2049), 12)),
            "lam_eff_over_lam": str(sp.N((lam_eff / lam).subs(qw2049), 12)),
            "y_eff_over_y": str(sp.N((y_eff / y).subs(qw2049), 12)),
        },
        "strict_core_closure_missing_items": {
            "renormalization": "Need full diagrammatic values for G_R and finite parts.",
            "unitarity": "Need solved Cutkosky ImM equations linked to G_U.",
            "background_independence": "Need transport theorem/witness linked to G_BI.",
            "selector_obstruction": "QW-2191 remains open.",
        },
        "false_pass_guard": "Gate vector export is not theorem-level strict-core closure.",
        "next_honest_step": "Solve coupled gate equations numerically/symbolically and test reverse contraction stability under branch perturbations.",
        "lay_explanation": "To kompresja całego toru do trzech bram kwantowych: renormalizacja, unitarność i niezależność od tła muszą zadziałać razem.",
    }

    path = GEN / "p1889_s839_strict_full_lagrangian_to_joint_qg_gate_probe.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
