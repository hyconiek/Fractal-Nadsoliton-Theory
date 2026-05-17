#!/usr/bin/env python3
"""P1878 S828 strict one-loop imaginary amplitude binding probe for channel defects."""
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
    p1876 = load("p1876_s826_strict_channel_amplitude_discontinuity_table_probe.json")
    p1877 = load("p1877_s827_strict_bidirectional_closure_matrix_probe.json")

    beta, eta, omega, phi = sp.symbols("beta eta omega phi", positive=True, real=True)
    lam, y, kappa2 = sp.symbols("lam y kappa2", positive=True, real=True)
    rho2_ss, rho2_sf, rho2_ff = sp.symbols("rho2_ss rho2_sf rho2_ff", positive=True, real=True)

    c0 = sp.cos(omega + phi) / (beta + 1)
    c1 = (-beta * eta * sp.cos(omega + phi) - omega * (beta + 1) * sp.sin(omega + phi)) / (beta + 1) ** 2
    lam_eff = lam * (1 + c1**2)
    y_eff = y * (1 + c0 / 2)

    M_tree_ss = 6 * lam_eff
    M_tree_sf = y_eff**2
    M_tree_ff = y_eff**2 + kappa2

    # Explicit one-loop ImM exports (symbolic strict-side binding ansatz).
    b_ss, b_sf, b_ff = sp.symbols("b_ss b_sf b_ff", real=True)
    ImM_ss_1loop = b_ss * rho2_ss * M_tree_ss**2
    ImM_sf_1loop = b_sf * rho2_sf * M_tree_sf**2
    ImM_ff_1loop = b_ff * rho2_ff * M_tree_ff**2

    defect_ss = sp.expand(ImM_ss_1loop - rho2_ss * M_tree_ss**2)
    defect_sf = sp.expand(ImM_sf_1loop - rho2_sf * M_tree_sf**2)
    defect_ff = sp.expand(ImM_ff_1loop - rho2_ff * M_tree_ff**2)

    closure_solution = {
        "b_ss": "1",
        "b_sf": "1",
        "b_ff": "1",
    }

    qw2049 = {beta: 1.0, eta: 1.8, omega: 0.18575, phi: 0.16250}

    out = {
        "packet_id": "P1878",
        "stage_id": "S828",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "depends_on": {
            "p1876_present": "cutkosky_discontinuity_table" in p1876,
            "p1877_present": "strict_core_closure_matrix" in p1877,
        },
        "strict_chain_step": "K_strict -> effective couplings -> channel tree amplitudes -> explicit ImM_1loop bindings -> channel defect equations",
        "effective_couplings": {
            "lam_eff": str(lam_eff),
            "y_eff": str(y_eff),
        },
        "channel_tree_amplitudes": {
            "M_tree_ss": str(M_tree_ss),
            "M_tree_sf": str(M_tree_sf),
            "M_tree_ff": str(M_tree_ff),
        },
        "imaginary_amplitude_export_ansatz": {
            "ImM_ss_1loop": str(ImM_ss_1loop),
            "ImM_sf_1loop": str(ImM_sf_1loop),
            "ImM_ff_1loop": str(ImM_ff_1loop),
            "binding_note": "b_* encode missing explicit loop-integral/channel integration coefficients in same scheme.",
        },
        "channel_defect_equations": {
            "defect_ss": str(defect_ss),
            "defect_sf": str(defect_sf),
            "defect_ff": str(defect_ff),
            "formal_zero_solution": closure_solution,
        },
        "qw2049_trace": {
            "M_tree_ss_over_lam": str(sp.N((M_tree_ss / lam).subs(qw2049), 12)),
            "M_tree_sf_over_y2": str(sp.N((M_tree_sf / y**2).subs(qw2049), 12)),
        },
        "strict_core_closure_missing_items": {
            "renormalization": "Need explicit loop integrals fixing finite parts consistently with b_*.",
            "unitarity": "Need computed ImM exports proving b_ss=b_sf=b_ff=1 physically, not by ansatz.",
            "background_independence": "Need FRW->BianchiI transport of same ImM/defect system.",
            "selector_obstruction": "QW-2191 remains open.",
        },
        "false_pass_guard": "Binding ansatz is not a computed loop result; strict closure remains open.",
        "next_honest_step": "Replace b_* ansatz by explicit one-loop integrals and run first FRW->BianchiI defect transport check.",
        "lay_explanation": "Tu precyzujemy, jakie dokładnie części liczb zespolonych amplitud trzeba policzyć, żeby sprawdzić unitarność. Na razie to kontrakt obliczeniowy, nie gotowy dowód.",
    }

    path = GEN / "p1878_s828_strict_one_loop_imaginary_amplitude_binding_probe.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
