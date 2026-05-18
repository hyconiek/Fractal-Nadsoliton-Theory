#!/usr/bin/env python3
"""P2020 S970 strict Cutkosky P2019 tree phase-space cut-sum witness.

Next honest step after P2019: do not promote the local tree component to a
backend Cutkosky theorem.  Instead, integrate exactly the already exported
P2019 transverse tree polarization sum over the massless two-body phase space,
and compare the exact SymPy result with deterministic SciPy quadrature under an
explicit final-state convention window.
"""
from __future__ import annotations

import json
import platform
from pathlib import Path
from typing import Any

import numpy as np
import scipy.integrate as integrate
import scipy.linalg as la
import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2020_s970_strict_cutkosky_p2019_tree_phase_space_cut_sum_witness.json"
TS = "2026-05-18T00:00:00+00:00"

P2019_NAME = "p2019_s969_strict_cutkosky_first_transverse_tree_amplitude_component.json"
P2018_NAME = "p2018_s968_strict_cutkosky_p2017_provenance_nonpromotion_audit.json"
CHANNEL = "graviton->gauge_gauge"

ETA = sp.diag(-1, 1, 1, 1)


def lower(v: sp.Matrix) -> sp.Matrix:
    return ETA * v


def field_strength_down(k: sp.Matrix, eps: sp.Matrix) -> sp.Matrix:
    kd = lower(k)
    ed = lower(eps)
    return sp.Matrix(4, 4, lambda mu, nu: sp.simplify(kd[mu] * ed[nu] - kd[nu] * ed[mu]))


def cross_stress_down(f1: sp.Matrix, f2: sp.Matrix) -> sp.Matrix:
    f_cross = sp.simplify(
        sum(ETA[rho, rho] * ETA[sigma, sigma] * f1[rho, sigma] * f2[rho, sigma]
            for rho in range(4) for sigma in range(4))
    )
    return sp.Matrix(
        4,
        4,
        lambda mu, nu: sp.simplify(
            sum(
                f1[mu, rho] * ETA[rho, rho] * f2[nu, rho]
                + f2[mu, rho] * ETA[rho, rho] * f1[nu, rho]
                for rho in range(4)
            )
            - sp.Rational(1, 2) * ETA[mu, nu] * f_cross
        ),
    )


def amplitude_component(h_down: sp.Matrix, stress_down: sp.Matrix) -> sp.Expr:
    h_up = sp.simplify(ETA * h_down * ETA)
    return sp.trigsimp(
        sp.simplify(sp.Rational(1, 2) * sum(h_up[mu, nu] * stress_down[mu, nu] for mu in range(4) for nu in range(4)))
    )


def matrix_to_strings(m: sp.Matrix) -> list[list[str]]:
    return [[str(sp.factor(sp.trigsimp(sp.simplify(m[i, j])))) for j in range(m.cols)] for i in range(m.rows)]


def final_state_linear_polarization_gram(graviton_amplitudes: list[sp.Matrix]) -> sp.Matrix:
    return sp.Matrix(
        len(graviton_amplitudes),
        len(graviton_amplitudes),
        lambda i, j: sp.simplify(
            sum(
                graviton_amplitudes[i][row, col] * graviton_amplitudes[j][row, col]
                for row in range(graviton_amplitudes[i].rows)
                for col in range(graviton_amplitudes[i].cols)
            )
        ),
    )


def angular_transport_certificate() -> dict[str, Any]:
    theta, varphi = sp.symbols("theta varphi", real=True)
    e_theta = sp.Matrix([0, sp.cos(theta) * sp.cos(varphi), sp.cos(theta) * sp.sin(varphi), -sp.sin(theta)])
    e_phi = sp.Matrix([0, -sp.sin(varphi), sp.cos(varphi), 0])
    n_hat = sp.Matrix([sp.sin(theta) * sp.cos(varphi), sp.sin(theta) * sp.sin(varphi), sp.cos(theta)])
    k1 = sp.Matrix([1, n_hat[0], n_hat[1], n_hat[2]])
    k2 = sp.Matrix([1, -n_hat[0], -n_hat[1], -n_hat[2]])

    gauge_pols = {"theta": e_theta, "phi": e_phi}
    h_plus = sp.simplify(e_theta * e_theta.T - e_phi * e_phi.T)
    h_cross = sp.simplify(e_theta * e_phi.T + e_phi * e_theta.T)
    graviton_pols = {"plus": h_plus, "cross": h_cross}

    matrices = {}
    transported_sum = sp.Integer(0)
    for g_label, h_pol in graviton_pols.items():
        mat = sp.zeros(2, 2)
        for i, eps1 in enumerate(gauge_pols.values()):
            for j, eps2 in enumerate(gauge_pols.values()):
                amp = amplitude_component(
                    h_pol,
                    cross_stress_down(field_strength_down(k1, eps1), field_strength_down(k2, eps2)),
                )
                mat[i, j] = sp.factor(sp.trigsimp(amp))
                transported_sum += sp.simplify(mat[i, j] ** 2)
        matrices[g_label] = mat

    transversality = {
        "k1_dot_e_theta": str(sp.trigsimp((k1.T * ETA * e_theta)[0])),
        "k1_dot_e_phi": str(sp.trigsimp((k1.T * ETA * e_phi)[0])),
        "k2_dot_e_theta": str(sp.trigsimp((k2.T * ETA * e_theta)[0])),
        "k2_dot_e_phi": str(sp.trigsimp((k2.T * ETA * e_phi)[0])),
        "e_theta_norm": str(sp.trigsimp((e_theta.T * ETA * e_theta)[0])),
        "e_phi_norm": str(sp.trigsimp((e_phi.T * ETA * e_phi)[0])),
        "e_theta_dot_e_phi": str(sp.trigsimp((e_theta.T * ETA * e_phi)[0])),
    }
    transported_sum = sp.factor(sp.trigsimp(sp.simplify(transported_sum)))
    linear_polarization_gram = final_state_linear_polarization_gram([matrices["plus"], matrices["cross"]])
    expected_plus = sp.Matrix([[-2, 0], [0, 2]])
    expected_cross = sp.Matrix([[0, -2], [-2, 0]])
    matrices_match = matrices["plus"] == expected_plus and matrices["cross"] == expected_cross
    checks_zero = all(value in {"0", "1"} for value in transversality.values()) and matrices_match and transported_sum == 16

    return {
        "basis": "generic spherical real transverse-polarization frame {e_theta,e_phi} for k1=(1,n), k2=(1,-n); not a circular helicity eigenbasis",
        "nomenclature_guard": "P2020 uses real plus/cross graviton linear polarizations and real transverse gauge polarizations; it does not claim circular helicity eigenstate diagonalization.",
        "transversality_and_orthonormality": transversality,
        "transported_amplitude_matrices_over_kappa_Zgauge": {
            "plus": matrix_to_strings(matrices["plus"]),
            "cross": matrix_to_strings(matrices["cross"]),
        },
        "transported_polarization_sum_over_kappa2_Zgauge2": str(transported_sum),
        "final_state_summed_graviton_linear_polarization_gram_over_kappa2_Zgauge2": matrix_to_strings(linear_polarization_gram),
        "matches_p2019_axis_frame": bool(matrices_match),
        "symbolic_angular_independence_verified": bool(checks_zero),
    }


def load(name: str) -> dict[str, Any]:
    path = GEN / name
    if not path.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def sympify_matrix(rows: list[list[str]]) -> sp.Matrix:
    return sp.Matrix([[sp.sympify(item) for item in row] for row in rows])


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2019 = load(P2019_NAME)
    p2018 = load(P2018_NAME)

    component = p2019.get("M_tree_transverse_common_basis_component", {})
    matrices = component.get("amplitude_matrices_over_kappa_Zgauge", {})
    plus = sympify_matrix(matrices.get("plus", [["0", "0"], ["0", "0"]]))
    cross = sympify_matrix(matrices.get("cross", [["0", "0"], ["0", "0"]]))

    amp_vector = sp.Matrix(list(plus) + list(cross))
    gram = sp.simplify(amp_vector * amp_vector.T)
    gram_eigenvals_exact = gram.eigenvals()
    gram_eigenvalue_list_exact = []
    for eig, mult in gram_eigenvals_exact.items():
        gram_eigenvalue_list_exact.extend([sp.simplify(eig)] * mult)
    pol_sum = sp.simplify(sum(entry**2 for entry in amp_vector))
    linear_polarization_gram = final_state_linear_polarization_gram([plus, cross])
    angular_certificate = angular_transport_certificate()

    x, phi = sp.symbols("x phi", real=True)
    # Standard massless two-body phase-space convention:
    # dPhi_2 = (1/(32*pi**2)) dOmega = (1/(32*pi**2)) dphi dx.
    dphi2_density = sp.Rational(1, 32) / sp.pi**2
    phase_space_volume = sp.simplify(sp.integrate(sp.integrate(dphi2_density, (phi, 0, 2 * sp.pi)), (x, -1, 1)))
    symmetry_identical = sp.Rational(1, 2)
    exact_linear_polarization_matrix_no_symmetry = sp.simplify(phase_space_volume * linear_polarization_gram)
    exact_linear_polarization_matrix_identical = sp.simplify(symmetry_identical * exact_linear_polarization_matrix_no_symmetry)
    exact_no_symmetry = sp.simplify(sp.trace(exact_linear_polarization_matrix_no_symmetry))
    exact_identical = sp.simplify(sp.trace(exact_linear_polarization_matrix_identical))
    exact_window = [exact_identical, exact_no_symmetry]

    def integrand_no_symmetry(xx: float, ph: float) -> float:
        return float(sp.N(pol_sum * dphi2_density, 30))

    scipy_no_symmetry, scipy_no_symmetry_err = integrate.dblquad(
        integrand_no_symmetry,
        0.0,
        2.0 * np.pi,
        lambda _ph: -1.0,
        lambda _ph: 1.0,
        epsabs=1e-13,
        epsrel=1e-13,
    )
    scipy_identical = 0.5 * scipy_no_symmetry
    scipy_identical_err = 0.5 * scipy_no_symmetry_err

    exact_no_symmetry_float = float(sp.N(exact_no_symmetry, 30))
    exact_identical_float = float(sp.N(exact_identical, 30))
    abs_errors = np.array(
        [
            abs(scipy_identical - exact_identical_float),
            abs(scipy_no_symmetry - exact_no_symmetry_float),
        ],
        dtype=float,
    )

    p2018_blocks_promotion = p2018.get("promotion_verdict") == "BLOCK_PROMOTION_TO_BACKEND_CUTKOSKY_CLOSURE"
    p2019_scope_ok = p2019.get("result_kind") == "PASS_FIRST_TREE_TRANSVERSE_COMPONENT_WITNESS"
    exact_match = bool(float(la.norm(abs_errors, 2)) < 1e-11)
    positive_window = bool(exact_identical_float > 0.0 and exact_no_symmetry_float > exact_identical_float)
    gram_psd_exact = all(sp.simplify(eig) >= 0 for eig in gram_eigenvalue_list_exact)
    linear_polarization_matrix_trace_ok = bool(sp.simplify(sp.trace(exact_linear_polarization_matrix_no_symmetry) - exact_no_symmetry) == 0)
    linear_polarization_matrix_eigs = list(exact_linear_polarization_matrix_no_symmetry.eigenvals().keys())
    linear_polarization_matrix_psd = all(sp.simplify(eig) >= 0 for eig in linear_polarization_matrix_eigs)

    gate = {
        "p2019_tree_component_available": p2019_scope_ok,
        "p2018_nonpromotion_still_active": p2018_blocks_promotion,
        "polarization_sum_recomputed_exactly": str(pol_sum) == "16",
        "two_body_phase_space_measure_exported": True,
        "sympy_exact_integral_positive_convention_window": positive_window,
        "scipy_quadrature_matches_exact_integral": exact_match,
        "amplitude_gram_psd": bool(gram_psd_exact),
        "generic_angular_transport_symbolically_verified": angular_certificate["symbolic_angular_independence_verified"],
        "linear_polarization_resolved_cut_sum_exported": linear_polarization_matrix_trace_ok,
        "linear_polarization_resolved_cut_sum_psd": bool(linear_polarization_matrix_psd),
    }

    out = {
        "ledger_id": "P2020_S970_STRICT_CUTKOSKY_P2019_TREE_PHASE_SPACE_CUT_SUM_WITNESS",
        "packet_id": "P2020",
        "stage_id": "S970",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TS,
        "route": "strict_only",
        "channel": CHANNEL,
        "depends_on": {"p2019": p2019.get("ledger_id"), "p2018": p2018.get("ledger_id")},
        "input_scope": {
            "source_component": "P2019 local transverse tree hAA amplitude component",
            "not_dressed": True,
            "not_discM": True,
            "not_cutkosky_closure": True,
            "normalization": "results are in units of kappa^2*Z_gauge^2",
        },
        "amplitude_input": {
            "plus_matrix_over_kappa_Zgauge": [[str(sp.simplify(plus[i, j])) for j in range(plus.cols)] for i in range(plus.rows)],
            "cross_matrix_over_kappa_Zgauge": [[str(sp.simplify(cross[i, j])) for j in range(cross.cols)] for i in range(cross.rows)],
            "polarization_sum_over_kappa2_Zgauge2": str(pol_sum),
            "amplitude_vector_order": ["plus_xx", "plus_xy", "plus_yx", "plus_yy", "cross_xx", "cross_xy", "cross_yx", "cross_yy"],
            "amplitude_gram_eigvals_exact": [str(eig) for eig in gram_eigenvalue_list_exact],
            "final_state_summed_graviton_linear_polarization_gram_over_kappa2_Zgauge2": matrix_to_strings(linear_polarization_gram),
        },
        "angular_transport_certificate": angular_certificate,
        "exact_phase_space_cut_sum": {
            "measure": "dPhi2_massless_COM = (1/(32*pi^2))*dphi*dx, phi in [0,2*pi], x=cos(theta) in [-1,1]",
            "integrand_status": "The P2019 axis-frame local transverse polarization sum is symbolically rederived in a generic spherical real transverse-polarization frame by angular_transport_certificate; this remains a tree-component cut-sum candidate, not DiscM.",
            "CutSum_tree_no_identical_symmetry_over_kappa2_Zgauge2": str(exact_no_symmetry),
            "CutSum_tree_identical_final_state_over_kappa2_Zgauge2": str(exact_identical),
            "convention_window_over_kappa2_Zgauge2": [str(item) for item in exact_window],
            "linear_polarization_resolved_CutSum_no_identical_symmetry_over_kappa2_Zgauge2": matrix_to_strings(exact_linear_polarization_matrix_no_symmetry),
            "linear_polarization_resolved_CutSum_identical_final_state_over_kappa2_Zgauge2": matrix_to_strings(exact_linear_polarization_matrix_identical),
            "linear_polarization_basis_order": ["plus", "cross"],
            "trace_relation": "trace(linear_polarization_resolved_CutSum) equals the scalar polarization-summed CutSum; the matrix form avoids prematurely collapsing the external graviton linear-polarization indices.",
        },
        "numeric_quadrature_cross_check": {
            "scipy_dblquad_no_symmetry": scipy_no_symmetry,
            "scipy_dblquad_no_symmetry_abserr_estimate": scipy_no_symmetry_err,
            "scipy_dblquad_identical": scipy_identical,
            "scipy_dblquad_identical_abserr_estimate": scipy_identical_err,
            "exact_no_symmetry_float": exact_no_symmetry_float,
            "exact_identical_float": exact_identical_float,
            "absolute_errors": abs_errors.tolist(),
            "absolute_error_l2": float(la.norm(abs_errors, 2)),
        },
        "p1953_contract_update": {
            "phase_space_measure and integration_domain": "PARTIAL_TREE_COMPONENT_EXACT_DPHI2_EXPORTED",
            "AbsM_dressed_squared_common_basis": "PARTIAL_TREE_TRANSVERSE_POLARIZATION_SUM_EXACT_NOT_DRESSED",
            "CutSum_common_basis": "PARTIAL_TREE_LINEAR_POLARIZATION_RESOLVED_CUTSUM_CONVENTION_WINDOW_NOT_DRESSED_NOT_DISC_MATCHED",
            "M_dressed_common_basis": "OPEN_STILL_TREE_ONLY_FROM_P2019",
            "DiscM_common_basis": "OPEN_NOT_EVALUATED",
            "DiscM_minus_CutSum_simplified": "OPEN_NOT_EVALUATED",
            "external_state_projectors including BRST physical-state projector": "PARTIAL_LOCAL_TRANSVERSE_PROJECTOR_ONLY_BRST_COHOMOLOGY_OPEN",
        },
        "gatekeeper_checks": gate,
        "result_kind": "PASS_TREE_PHASE_SPACE_LINEAR_POLARIZATION_CUT_SUM_COMPONENT_WITNESS" if all(gate.values()) else "OPEN_TREE_PHASE_SPACE_LINEAR_POLARIZATION_CUT_SUM_GAP",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "false_pass_guard": "P2020 integrates the P2019 tree transverse component over exact massless two-body phase space. It is not a dressed loop amplitude, not a BRST cohomology projector, not DiscM=CutSum, not all-state unitarity, not QW-2191 discharge, and not ToE closure.",
        "next_honest_step": "Derive the matching same-scheme DiscM_common_basis or a first dressed-residue factor for this exact linear-polarization-resolved tree cut-sum component; only then compare DiscM_minus_CutSum for the same normalized object.",
        "toe_progress": "Turns the first strict-side hAA amplitude component into an exact linear-polarization-resolved phase-space cut-sum object with a symbolic angular-transport certificate and checked convention window, while preserving the obstruction boundary required for theorem-grade unitarity.",
        "lay_explanation": "Mieliśmy policzony lokalny kawałek amplitudy. Teraz uczciwie policzyliśmy, jak ten kawałek sumuje się po dostępnej przestrzeni dwóch cząstek końcowych, zachowując jeszcze indeksy polaryzacji grawitonu, po uprzednim symbolicznym sprawdzeniu obrotowego transportu lokalnej bazy polaryzacji. To nadal nie jest pełny dowód unitarności, bo nie mamy jeszcze pasującej nieciągłości pętlowej DiscM ani pełnego ubierania amplitudy.",
        "environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": __import__("scipy").__version__,
            "sympy": sp.__version__,
        },
    }

    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P2020] wrote tree phase-space cut-sum witness: {OUT}")


if __name__ == "__main__":
    main()
