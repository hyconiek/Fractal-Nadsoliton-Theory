#!/usr/bin/env python3
"""P1877 S827 strict bidirectional closure matrix probe (kernel->coeff->L->EOM and reverse)."""
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
    p1874 = load("p1874_s824_strict_full_lagrangian_eom_and_qg_witness_obligation_probe.json")
    p1876 = load("p1876_s826_strict_channel_amplitude_discontinuity_table_probe.json")

    beta, eta, omega, phi = sp.symbols("beta eta omega phi", positive=True, real=True)
    m2, lam, y = sp.symbols("m2 lam y", real=True)

    c0 = sp.cos(omega + phi) / (beta + 1)
    c1 = (-beta * eta * sp.cos(omega + phi) - omega * (beta + 1) * sp.sin(omega + phi)) / (beta + 1) ** 2

    m2_eff = m2 * (1 + c0)
    lam_eff = lam * (1 + c1**2)
    y_eff = y * (1 + c0 / 2)

    qw2049 = {beta: 1.0, eta: 1.8, omega: 0.18575, phi: 0.16250}

    out = {
        "packet_id": "P1877",
        "stage_id": "S827",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "depends_on": {
            "p1874_present": "strict_chain" in p1874,
            "p1876_present": "channel_amplitude_table" in p1876,
        },
        "strict_forward_chain": [
            "K_strict(d)=cos(omega*d+phi)/(1+beta*d^eta)",
            "effective_coefficients(m2_eff, lam_eff, y_eff)",
            "full_non_skeleton_L_total(L_SM + L_GR)",
            "covariant_EOM_and_Einstein_residuals",
            "channel_amplitudes_and_optical_defects",
        ],
        "strict_reverse_chain": [
            "optical_defect_zero_per_channel",
            "counterterm_consistency_in_same_scheme",
            "residual_zero_stability_under_background_transport",
            "admissible_strict_kernel_parameter_window",
        ],
        "effective_coefficients": {
            "m2_eff": str(m2_eff),
            "lam_eff": str(lam_eff),
            "y_eff": str(y_eff),
        },
        "qw2049_trace": {
            "m2_eff_over_m2": str(sp.N((m2_eff / m2).subs(qw2049), 12)),
            "lam_eff_over_lam": str(sp.N((lam_eff / lam).subs(qw2049), 12)),
            "y_eff_over_y": str(sp.N((y_eff / y).subs(qw2049), 12)),
        },
        "strict_core_closure_matrix": {
            "renormalization": {
                "forward_requirement": "explicit one-loop amplitudes + finite-part-fixed counterterms",
                "reverse_requirement": "renormalized residual invariance back-constrains coefficient admissibility",
                "status": "OPEN_MISSING_WITNESS",
            },
            "unitarity_cutkosky": {
                "forward_requirement": "ImM channel exports on same branch",
                "reverse_requirement": "defect_ss/sf/ff=0 constrains strict coupling corridors",
                "status": "OPEN_MISSING_WITNESS",
            },
            "background_independence": {
                "forward_requirement": "FRW->(BianchiI,static-spherical) transport theorem",
                "reverse_requirement": "transport-stable closure returns same strict parameter lane",
                "status": "OPEN_MISSING_THEOREM",
            },
            "selector": {
                "forward_requirement": "explicit strict selector source or symmetry-breaking theorem",
                "reverse_requirement": "QW-2191 obstruction discharged",
                "status": "OPEN_QW2191",
            },
        },
        "false_pass_guard": "Bidirectional matrix is a closure contract; not a proof of strict-core ToE closure.",
        "next_honest_step": "Bind explicit one-loop ImM channel exports to this matrix and run first FRW->BianchiI transport check under identical renormalization scheme.",
        "lay_explanation": "To mapa domknięcia w obie strony: od kernela do równań i z powrotem z testów kwantowych do ograniczeń kernela. Pokazuje dokładnie, czego jeszcze brakuje do pełnej teorii.",
    }

    path = GEN / "p1877_s827_strict_bidirectional_closure_matrix_probe.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
