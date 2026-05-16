#!/usr/bin/env python3
"""P1875 S825 strict one-loop counterterm and Cutkosky witness obligation probe."""
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

    beta, eta, omega, phi = sp.symbols("beta eta omega phi", positive=True, real=True)
    m2, lam, y, xi, kappa2, Lambda_cc = sp.symbols("m2 lam y xi kappa2 Lambda_cc", real=True)
    eps, mu = sp.symbols("eps mu", positive=True, real=True)
    a_phi, a_F, a_psi, a_xi, a_kappa, a_Lambda = sp.symbols(
        "a_phi a_F a_psi a_xi a_kappa a_Lambda", real=True
    )

    c0 = sp.cos(omega + phi) / (beta + 1)
    c1 = (-beta * eta * sp.cos(omega + phi) - omega * (beta + 1) * sp.sin(omega + phi)) / (beta + 1) ** 2
    m2_eff = sp.simplify(m2 * (1 + c0))
    lam_eff = sp.simplify(lam * (1 + c1**2))
    y_eff = sp.simplify(y * (1 + c0 / 2))

    # Minimal strict-side one-loop dictionary (symbolic placeholder, no false closure).
    DeltaZ_phi = a_phi * lam_eff / (16 * sp.pi**2 * eps)
    DeltaZ_F = a_F * y_eff**2 / (16 * sp.pi**2 * eps)
    DeltaZ_psi = a_psi * y_eff**2 / (16 * sp.pi**2 * eps)
    Delta_xi = a_xi * xi * lam_eff / (16 * sp.pi**2 * eps)
    Delta_kappa_inv = a_kappa * m2_eff / (16 * sp.pi**2 * eps)
    Delta_Lambda = a_Lambda * m2_eff**2 / (16 * sp.pi**2 * eps)

    # Cutkosky/optical consistency placeholders on same branch.
    Im_M_2to2 = sp.symbols("Im_M_2to2", real=True)
    phase_space_2 = sp.symbols("phase_space_2", nonnegative=True, real=True)
    M_tree_sq = sp.symbols("M_tree_sq", nonnegative=True, real=True)
    optical_defect = sp.simplify(Im_M_2to2 - phase_space_2 * M_tree_sq)

    qw2049 = {beta: 1.0, eta: 1.8, omega: 0.18575, phi: 0.16250}

    out = {
        "packet_id": "P1875",
        "stage_id": "S825",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "depends_on": {
            "p1874_present": "strict_chain" in p1874,
        },
        "strict_chain_step": "K_strict -> effective couplings -> one-loop counterterm dictionary candidate -> Cutkosky witness obligations",
        "effective_couplings": {
            "m2_eff": str(m2_eff),
            "lam_eff": str(lam_eff),
            "y_eff": str(y_eff),
        },
        "one_loop_counterterm_dictionary_candidate": {
            "DeltaZ_phi": str(DeltaZ_phi),
            "DeltaZ_F": str(DeltaZ_F),
            "DeltaZ_psi": str(DeltaZ_psi),
            "Delta_xi": str(Delta_xi),
            "Delta_kappa_inverse": str(Delta_kappa_inv),
            "Delta_Lambda_cc": str(Delta_Lambda),
            "scheme_note": "Symbolic MS-like 1/eps skeleton awaiting strict-sector amplitude matching coefficients a_*.",
        },
        "cutkosky_witness_table": {
            "channel": "2to2 scalar-sector proxy around same FRW branch",
            "optical_identity_form": "Im M(2->2) = Integral dPi_2 |M_tree|^2",
            "optical_defect_symbol": str(optical_defect),
            "pass_condition": "optical_defect_symbol == 0 after strict one-loop amplitude export",
        },
        "qw2049_trace": {
            "lam_eff_over_lam_qw2049": str(sp.N((lam_eff / lam).subs(qw2049), 12)),
            "y_eff_over_y_qw2049": str(sp.N((y_eff / y).subs(qw2049), 12)),
            "DeltaZ_phi_prefactor_qw2049": str(sp.N((DeltaZ_phi * eps / a_phi).subs(qw2049), 12)),
        },
        "strict_core_closure_missing_items": {
            "renormalization": "Need explicit loop integrals to fix a_* and verify cancellation in renormalized FRW residuals.",
            "unitarity": "Need explicit discontinuity computation proving optical_defect_symbol=0 channel-by-channel.",
            "background_independence": "Need atlas transport theorem showing branch-independent witness persistence.",
            "selector_obstruction": "QW-2191 remains open and unaffected by this packet.",
        },
        "false_pass_guard": "Counterterm dictionary candidate + optical form are not final witnesses without loop integral export and channel closure proofs.",
        "next_honest_step": "Export strict one-loop amplitude/discontinuity tables that determine a_* and evaluate optical_defect_symbol per channel.",
        "lay_explanation": "To krok kwantowy: zapisujemy, jakie poprawki 1-pętlowe muszą pojawić się w teorii i jak formalnie sprawdzić unitarność przez tożsamość optyczną. Nadal brakuje pełnych obliczeń całek i kanałów.",
    }

    path = GEN / "p1875_s825_strict_one_loop_counterterm_and_cutkosky_witness_probe.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
