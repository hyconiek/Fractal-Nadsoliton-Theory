#!/usr/bin/env python3
"""P1882 S832 strict kernel->QG closure obligation ledger probe."""
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
    p1881 = load("p1881_s831_strict_termwise_one_loop_binding_matrix_probe.json")
    p1879 = load("p1879_s829_strict_frw_to_bianchiI_transport_contract_probe.json")

    beta, eta, omega, phi = sp.symbols("beta eta omega phi", positive=True, real=True)
    m2, lam, y = sp.symbols("m2 lam y", real=True)

    c0 = sp.cos(omega + phi) / (beta + 1)
    c1 = (-beta * eta * sp.cos(omega + phi) - omega * (beta + 1) * sp.sin(omega + phi)) / (beta + 1) ** 2

    m2_eff = m2 * (1 + c0)
    lam_eff = lam * (1 + c1**2)
    y_eff = y * (1 + c0 / 2)

    qw2049 = {beta: 1.0, eta: 1.8, omega: 0.18575, phi: 0.16250}

    out = {
        "packet_id": "P1882",
        "stage_id": "S832",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "depends_on": {
            "p1881_present": "termwise_one_loop_binding_matrix" in p1881,
            "p1879_present": "bianchiI_transport_ansatz" in p1879,
        },
        "strict_chain": {
            "forward": [
                "K_strict",
                "effective_coefficients",
                "full_non_skeleton_L_total",
                "covariant_EOM",
                "channel_ImM_and_optical_defects",
                "FRW_to_BianchiI_transport",
                "QG_closure_ledger",
            ],
            "reverse": [
                "QG_witnesses",
                "transport_stability",
                "optical_defect_zero",
                "renormalized_EOM_consistency",
                "admissible_strict_kernel_window",
            ],
        },
        "effective_coefficients": {
            "m2_eff": str(m2_eff),
            "lam_eff": str(lam_eff),
            "y_eff": str(y_eff),
        },
        "qg_closure_obligation_ledger": {
            "R1_renormalization": {
                "required_export": "termwise finite-part-fixed one-loop integrals attached to z_* matrix",
                "status": "OPEN_MISSING_WITNESS",
            },
            "U1_unitarity": {
                "required_export": "computed ImM channel tables with defect_ss/sf/ff=0",
                "status": "OPEN_MISSING_WITNESS",
            },
            "B1_background_independence": {
                "required_export": "theorem-grade FRW<->BianchiI transport closure on identical renormalization data",
                "status": "OPEN_MISSING_THEOREM",
            },
            "S1_selector": {
                "required_export": "strict selector source/theorem discharging QW-2191",
                "status": "OPEN_QW2191",
            },
        },
        "qw2049_trace": {
            "m2_eff_over_m2": str(sp.N((m2_eff / m2).subs(qw2049), 12)),
            "lam_eff_over_lam": str(sp.N((lam_eff / lam).subs(qw2049), 12)),
            "y_eff_over_y": str(sp.N((y_eff / y).subs(qw2049), 12)),
        },
        "false_pass_guard": "Ledger export is a strict-core closure contract, not ToE closure proof.",
        "next_honest_step": "Fill R1/U1/B1 with explicit computed witnesses and then run reverse-kernel admissibility contraction.",
        "lay_explanation": "To checklista końcowa: co dokładnie musi być jeszcze policzone i udowodnione, żeby teoria strict mogła być uczciwie uznana za domkniętą.",
    }

    path = GEN / "p1882_s832_strict_kernel_to_qg_closure_obligation_ledger_probe.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
