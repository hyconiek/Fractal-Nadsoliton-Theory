#!/usr/bin/env python3
"""P1885 S835 strict joint witness solver contract for R/U/BI blocks."""
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
    p1884 = load("p1884_s834_strict_explicit_loop_integral_stub_evaluation_probe.json")

    beta, eta, omega, phi = sp.symbols("beta eta omega phi", positive=True, real=True)
    lam, y, kappa2 = sp.symbols("lam y kappa2", positive=True, real=True)
    eps, mu2, s = sp.symbols("eps mu2 s", positive=True, real=True)
    sigma2 = sp.symbols("sigma2", nonnegative=True, real=True)

    c0 = sp.cos(omega + phi) / (beta + 1)
    c1 = (-beta * eta * sp.cos(omega + phi) - omega * (beta + 1) * sp.sin(omega + phi)) / (beta + 1) ** 2
    lam_eff = lam * (1 + c1**2)
    y_eff = y * (1 + c0 / 2)

    I_lam = 1 / (16 * sp.pi**2) * (1 / eps + sp.log(mu2 / s))
    I_y = 1 / (16 * sp.pi**2) * (1 / eps + sp.log(mu2 / s) + sp.Rational(1, 2))

    Mss = 6 * lam_eff * (1 + lam_eff * I_lam)
    Msf = y_eff**2 * (1 + y_eff**2 * I_y)
    Mff = (y_eff**2 + kappa2) * (1 + y_eff**2 * I_y)

    rho2_ss, rho2_sf, rho2_ff = sp.symbols("rho2_ss rho2_sf rho2_ff", nonnegative=True, real=True)
    ImM_ss, ImM_sf, ImM_ff = sp.symbols("ImM_ss ImM_sf ImM_ff", real=True)
    a_ss, a_sf, a_ff = sp.symbols("a_ss a_sf a_ff", real=True)

    eq_cut = {
        "ss": sp.expand(ImM_ss - rho2_ss * Mss**2),
        "sf": sp.expand(ImM_sf - rho2_sf * Msf**2),
        "ff": sp.expand(ImM_ff - rho2_ff * Mff**2),
    }
    eq_transport = {
        "ss": sp.expand(Mss * (1 + a_ss * sigma2) - Mss),
        "sf": sp.expand(Msf * (1 + a_sf * sigma2) - Msf),
        "ff": sp.expand(Mff * (1 + a_ff * sigma2) - Mff),
    }

    out = {
        "packet_id": "P1885",
        "stage_id": "S835",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "depends_on": {"p1884_present": "renormalized_amplitudes" in p1884},
        "strict_chain_step": "K_strict -> coefficients -> full-L one-loop data -> joint R/U/BI witness equation system",
        "joint_witness_solver_contract": {
            "inputs_required": [
                "diagram_resolved_integrals_for_I_lam_I_y",
                "ImM_ss_ImM_sf_ImM_ff",
                "a_ss_a_sf_a_ff",
            ],
            "equation_blocks": {
                "cutkosky": {k: str(v) for k, v in eq_cut.items()},
                "transport": {k: str(v) for k, v in eq_transport.items()},
            },
            "target_conditions": [
                "eq_cut_ss=0, eq_cut_sf=0, eq_cut_ff=0",
                "eq_transport_ss=eq_transport_sf=eq_transport_ff=0 in FRW limit and controlled in sigma2 corridor",
                "same renormalization data across all blocks",
            ],
        },
        "full_lagrangian_context": {
            "statement": "Nadsoliton => L_SM + L_GR with strict-kernel effective coefficients and one-loop corrections tested jointly.",
            "no_bridge_note": "strict-only route; no legacy bridge assumption used.",
        },
        "strict_core_closure_missing_items": {
            "renormalization": "Need computed finite parts from full loop diagrams.",
            "unitarity": "Need explicit ImM values satisfying Cutkosky block.",
            "background_independence": "Need Bianchi-I transport coefficients from same loop dataset.",
            "selector_obstruction": "QW-2191 remains open.",
        },
        "false_pass_guard": "Solver contract is not solved witness; strict-core closure remains open.",
        "next_honest_step": "Numerically/symbolically solve the joint block with computed loop data, then run reverse admissibility contraction to strict kernel corridor.",
        "lay_explanation": "To jest plan rozwiązania całości na raz: te same dane kwantowe mają jednocześnie spełnić warunki unitarności, renormalizacji i zgodności między tłami.",
    }

    path = GEN / "p1885_s835_strict_joint_witness_solver_contract_probe.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
