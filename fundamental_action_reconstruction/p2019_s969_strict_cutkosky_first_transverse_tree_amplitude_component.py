#!/usr/bin/env python3
"""P2019 S969 strict Cutkosky first transverse tree-amplitude component.

Next honest step after P2018: stop adding quadrature ansatz layers and compute
one actual, local, same-channel amplitude component from already exported strict
objects: the P1955 minimal hAA vertex and the P1956 local transverse two-gauge
projector.  This is intentionally tree-level/local-transverse only; it is not a
dressed amplitude, not a BRST cohomology theorem, and not Cutkosky closure.
"""
from __future__ import annotations

import json
import platform
from pathlib import Path
from typing import Any

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2019_s969_strict_cutkosky_first_transverse_tree_amplitude_component.json"
TS = "2026-05-18T00:00:00+00:00"

ETA = sp.diag(-1, 1, 1, 1)
CHANNEL = "graviton->gauge_gauge"


def load(name: str) -> dict[str, Any]:
    path = GEN / name
    if not path.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def vec(entries: tuple[int, int, int, int]) -> sp.Matrix:
    return sp.Matrix([sp.Integer(item) for item in entries])


def lower(v: sp.Matrix) -> sp.Matrix:
    return ETA * v


def field_strength_down(k: sp.Matrix, eps: sp.Matrix) -> sp.Matrix:
    kd = lower(k)
    ed = lower(eps)
    return sp.Matrix(4, 4, lambda mu, nu: sp.simplify(kd[mu] * ed[nu] - kd[nu] * ed[mu]))


def cross_stress_down(f1: sp.Matrix, f2: sp.Matrix) -> sp.Matrix:
    """Cross term in Maxwell stress tensor for two gauge waves.

    T_cross_{mu nu} = F1_{mu rho} F2_nu^rho + F2_{mu rho} F1_nu^rho
                      - 1/2 eta_{mu nu} F1_{rho sigma} F2^{rho sigma}.
    The overall Z_gauge factor is kept outside.
    """
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


def raise2(h_down: sp.Matrix) -> sp.Matrix:
    return sp.simplify(ETA * h_down * ETA)


def amplitude_component(h_down: sp.Matrix, stress_down: sp.Matrix) -> sp.Expr:
    # From P1955: L_hAA = (kappa/2) h^{mu nu} T_mu_nu.
    h_up = raise2(h_down)
    return sp.simplify(sp.Rational(1, 2) * sum(h_up[mu, nu] * stress_down[mu, nu] for mu in range(4) for nu in range(4)))


def matrix_to_strings(m: sp.Matrix) -> list[list[str]]:
    return [[str(sp.simplify(m[i, j])) for j in range(m.cols)] for i in range(m.rows)]


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p1955 = load("p1955_s905_strict_minimal_hAA_vertex_export.json")
    p1956 = load("p1956_s906_strict_gauge_gauge_physical_projector_certificate.json")
    p2018 = load("p2018_s968_strict_cutkosky_p2017_provenance_nonpromotion_audit.json")

    # Center-of-mass massless two-body cut frame matching P1956, E=1.
    k1 = vec((1, 0, 0, 1))
    k2 = vec((1, 0, 0, -1))
    eps_x = vec((0, 1, 0, 0))
    eps_y = vec((0, 0, 1, 0))
    gauge_pols = {"x": eps_x, "y": eps_y}

    h_plus = sp.zeros(4, 4)
    h_plus[1, 1] = 1
    h_plus[2, 2] = -1
    h_cross = sp.zeros(4, 4)
    h_cross[1, 2] = h_cross[2, 1] = 1
    graviton_pols = {"plus": h_plus, "cross": h_cross}

    amplitude_table: list[dict[str, Any]] = []
    amp_matrices: dict[str, list[list[str]]] = {}
    abs_sq_sum = sp.Integer(0)

    for g_label, h_pol in graviton_pols.items():
        mat = sp.zeros(2, 2)
        for i, (a_label, eps1) in enumerate(gauge_pols.items()):
            for j, (b_label, eps2) in enumerate(gauge_pols.items()):
                f1 = field_strength_down(k1, eps1)
                f2 = field_strength_down(k2, eps2)
                t_cross = cross_stress_down(f1, f2)
                amp = amplitude_component(h_pol, t_cross)
                mat[i, j] = amp
                abs_sq_sum += sp.simplify(amp**2)
                amplitude_table.append({
                    "graviton_polarization": g_label,
                    "gauge_polarizations": [a_label, b_label],
                    "M_tree_transverse_component_over_kappa_Zgauge": str(amp),
                })
        amp_matrices[g_label] = matrix_to_strings(mat)

    # Gauge transversality and traceless graviton checks for the local component.
    transversality_checks = []
    for label, eps in gauge_pols.items():
        transversality_checks.append({
            "polarization": label,
            "k1_dot_eps": str((k1.T * ETA * eps)[0]),
            "k2_dot_eps": str((k2.T * ETA * eps)[0]),
        })
    graviton_checks = {
        name: {
            "trace_eta_h": str(sp.trace(ETA * h)),
            "h_00": str(h[0, 0]),
            "h_03": str(h[0, 3]),
            "h_33": str(h[3, 3]),
        }
        for name, h in graviton_pols.items()
    }

    numeric_matrix_plus = np.array([[float(sp.N(sp.sympify(x))) for x in row] for row in amp_matrices["plus"]])
    numeric_matrix_cross = np.array([[float(sp.N(sp.sympify(x))) for x in row] for row in amp_matrices["cross"]])

    p1955_vertex_ok = (p1955.get("machine_check", {}) or {}).get("all_local_checks_zero") is True
    p1956_projector_ok = (p1956.get("machine_check_summary", {}) or {}).get("all_local_projector_checks_zero") is True
    nonpromotion_still_required = p2018.get("promotion_verdict") == "BLOCK_PROMOTION_TO_BACKEND_CUTKOSKY_CLOSURE"

    gate = {
        "p1955_minimal_vertex_available": p1955_vertex_ok,
        "p1956_local_transverse_projector_available": p1956_projector_ok,
        "p2018_nonpromotion_still_active": nonpromotion_still_required,
        "amplitude_table_complete": len(amplitude_table) == 8,
        "nonzero_component_exported": any(row["M_tree_transverse_component_over_kappa_Zgauge"] != "0" for row in amplitude_table),
        "gauge_transversality_reference_frame": all(item["k1_dot_eps"] == "0" and item["k2_dot_eps"] == "0" for item in transversality_checks),
        "graviton_polarizations_traceless_reference_frame": all(item["trace_eta_h"] == "0" for item in graviton_checks.values()),
        "abs_sq_sum_positive": bool(sp.simplify(abs_sq_sum) > 0),
    }

    out = {
        "ledger_id": "P2019_S969_STRICT_CUTKOSKY_FIRST_TRANSVERSE_TREE_AMPLITUDE_COMPONENT",
        "packet_id": "P2019",
        "stage_id": "S969",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TS,
        "route": "strict_only",
        "channel": CHANNEL,
        "depends_on": {
            "p1955": p1955.get("packet_id"),
            "p1956": p1956.get("packet_id"),
            "p2018": p2018.get("ledger_id"),
        },
        "component_scope": {
            "object_name": "TreeTransverse_hAA_Component_graviton_to_gauge_gauge_strict_local_v1",
            "scheme_scope": "minimal tree-level hAA vertex from P1955 plus local transverse projector from P1956; not MSbar_B1_seed dressed loop scheme",
            "normalization": "amplitudes are reported in units of kappa*Z_gauge",
            "frame": "massless two-body cut center-of-mass frame, E=1",
            "not_dressed": True,
            "not_brst_cohomology": True,
        },
        "kinematics": {
            "metric_signature": "diag(-1,1,1,1)",
            "k1": [str(x) for x in k1],
            "k2": [str(x) for x in k2],
            "gauge_polarizations": {name: [str(x) for x in eps] for name, eps in gauge_pols.items()},
            "graviton_polarizations": {name: matrix_to_strings(h) for name, h in graviton_pols.items()},
        },
        "M_tree_transverse_common_basis_component": {
            "basis": "{h_plus,h_cross} x {eps_x,eps_y}^{⊗2}",
            "amplitude_matrices_over_kappa_Zgauge": amp_matrices,
            "table": amplitude_table,
            "AbsM_tree_transverse_sum_over_local_polarizations_over_kappa2_Zgauge2": str(sp.simplify(abs_sq_sum)),
            "numeric_matrices": {
                "plus": numeric_matrix_plus.tolist(),
                "cross": numeric_matrix_cross.tolist(),
            },
        },
        "machine_checks": {
            "gauge_transversality": transversality_checks,
            "graviton_polarization_checks": graviton_checks,
            "gatekeeper_checks": gate,
        },
        "p1953_contract_update": {
            "M_dressed_common_basis": "PARTIAL_TREE_TRANSVERSE_COMPONENT_ONLY_NOT_DRESSED",
            "AbsM_dressed_squared_common_basis": "PARTIAL_TREE_TRANSVERSE_POLARIZATION_SUM_ONLY_NOT_DRESSED",
            "external_state_projectors including BRST physical-state projector": "PARTIAL_LOCAL_TRANSVERSE_PROJECTOR_ONLY_BRST_COHOMOLOGY_OPEN",
            "scheme = MSbar_B1_seed": "OPEN_NOT_LOCKED_TREE_LOCAL_SCOPE",
            "DiscM_common_basis": "OPEN_NOT_EXPORTED",
            "CutSum_common_basis": "OPEN_NOT_EXPORTED",
            "DiscM_minus_CutSum_simplified": "OPEN_NOT_EVALUATED",
        },
        "gatekeeper_checks": gate,
        "result_kind": "PASS_FIRST_TREE_TRANSVERSE_COMPONENT_WITNESS" if all(gate.values()) else "OPEN_COMPONENT_WITNESS_GAP",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "false_pass_guard": "P2019 exports only a local tree-level transverse hAA amplitude component. It is not a dressed loop amplitude, not a BRST cohomology projection, not DiscM=CutSum, not all-state unitarity, and not ToE closure.",
        "next_honest_step": "Promote this local component toward the P1953 contract by adding the same-scheme dressing/residue layer or by deriving a genuine BRST cohomology projector; then rerun P2018.",
        "toe_progress": "Provides the first explicit non-quadrature local amplitude component for the Cutkosky channel, replacing pure governance with a concrete strict-side tree-level calculation while keeping theorem closure blocked.",
        "lay_explanation": "Po raz pierwszy liczymy konkretny lokalny element amplitudy z wierzchołka grawiton–dwa bozony cechowania i poprzecznych polaryzacji. To realny kawałek rachunku, ale jeszcze nie pełna ubrana amplituda pętlowa ani dowód unitarności.",
        "environment": {"python": platform.python_version(), "numpy": np.__version__, "sympy": sp.__version__},
    }

    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P2019] wrote component witness: {OUT}")


if __name__ == "__main__":
    main()
