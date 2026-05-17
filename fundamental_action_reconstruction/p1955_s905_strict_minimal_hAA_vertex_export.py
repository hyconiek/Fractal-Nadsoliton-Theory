#!/usr/bin/env python3
"""P1955 S905 strict minimal hAA vertex export.

The required bilingual grep found P1657, which already exports the covariant
minimal gauge-metric density

    sqrt(-g) * [-1/4 g^{mu alpha} g^{nu beta} F_{mu nu} F_{alpha beta}].

This script derives the tree-level hAA vertex density from that term by the
metric perturbation g_{mu nu}=eta_{mu nu}+kappa*h_{mu nu}.  It exports a
machine-checked linearization identity:

    L_gauge[g] = L_gauge[eta] + kappa/2 h^{mu nu} T^gauge_{mu nu} + O(kappa^2).

This is not a dressed amplitude and not a Cutkosky theorem.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p1955_s905_strict_minimal_hAA_vertex_export.json"


def load(name: str) -> dict[str, Any]:
    path = GEN / name
    if not path.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def digest(obj: object) -> str:
    blob = json.dumps(obj, sort_keys=True, ensure_ascii=True).encode("utf-8")
    return hashlib.sha256(blob).hexdigest()


def minkowski_diag() -> list[int]:
    return [-1, 1, 1, 1]


def main() -> None:
    p1657 = load("p1657_s607_strict_helmholtz_h1_gauge_metric_covariant_summary.json")
    p1907 = load("p1907_s857_strict_full_lagrangian_to_eom_witness_matrix_probe.json")
    p1954 = load("p1954_s904_strict_dressed_amplitude_nonavailability_theorem.json")

    kappa, z_gauge = sp.symbols("kappa Z_gauge")
    eta = minkowski_diag()

    h_down = [[sp.symbols(f"h{mu}{nu}") for nu in range(4)] for mu in range(4)]
    f_down: list[list[sp.Expr]] = [[sp.Integer(0) for _ in range(4)] for _ in range(4)]
    for mu in range(4):
        for nu in range(4):
            if mu < nu:
                f_down[mu][nu] = sp.symbols(f"F{mu}{nu}")
            elif mu > nu:
                f_down[mu][nu] = -f_down[nu][mu]

    h_up = [[eta[mu] * eta[nu] * h_down[mu][nu] for nu in range(4)] for mu in range(4)]
    g_inv_linear = [
        [
            (sp.Integer(eta[mu]) if mu == nu else sp.Integer(0)) - kappa * h_up[mu][nu]
            for nu in range(4)
        ]
        for mu in range(4)
    ]
    h_trace = sum(eta[mu] * h_down[mu][mu] for mu in range(4))
    sqrt_minus_g_linear = 1 + kappa * sp.Rational(1, 2) * h_trace

    contraction = sp.Integer(0)
    for mu in range(4):
        for alpha in range(4):
            for nu in range(4):
                for beta in range(4):
                    contraction += (
                        g_inv_linear[mu][alpha]
                        * g_inv_linear[nu][beta]
                        * f_down[mu][nu]
                        * f_down[alpha][beta]
                    )

    lagrangian_linearized = -z_gauge * sp.Rational(1, 4) * sqrt_minus_g_linear * contraction
    flat_lagrangian = sp.simplify(lagrangian_linearized.subs(kappa, 0))
    linear_coefficient = sp.simplify(sp.diff(lagrangian_linearized, kappa).subs(kappa, 0))

    f_squared = sum(
        eta[rho] * eta[sigma] * f_down[rho][sigma] * f_down[rho][sigma]
        for rho in range(4)
        for sigma in range(4)
    )
    t_gauge = [[sp.Integer(0) for _ in range(4)] for _ in range(4)]
    for mu in range(4):
        for nu in range(4):
            f_mu_rho_f_nu_up_rho = sum(
                f_down[mu][rho] * eta[rho] * f_down[nu][rho]
                for rho in range(4)
            )
            eta_down = sp.Integer(eta[mu]) if mu == nu else sp.Integer(0)
            t_gauge[mu][nu] = sp.simplify(
                z_gauge * (f_mu_rho_f_nu_up_rho - sp.Rational(1, 4) * eta_down * f_squared)
            )

    hAA_density_from_stress = sp.Rational(1, 2) * sum(
        h_up[mu][nu] * t_gauge[mu][nu] for mu in range(4) for nu in range(4)
    )
    residual = sp.simplify(sp.expand(linear_coefficient - hAA_density_from_stress))

    # Gauge stress tensor should be symmetric and traceless in 4D Maxwell theory.
    symmetry_residuals = []
    for mu in range(4):
        for nu in range(4):
            symmetry_residuals.append(sp.simplify(t_gauge[mu][nu] - t_gauge[nu][mu]))
    trace_residual = sp.simplify(sum(eta[mu] * t_gauge[mu][mu] for mu in range(4)))

    all_zero = residual == 0 and trace_residual == 0 and all(item == 0 for item in symmetry_residuals)

    out = {
        "packet_id": "P1955",
        "stage_id": "S905",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "local_verdict": "MINIMAL_TREE_LEVEL_hAA_VERTEX_EXPORTED__DRESSED_CUTKOSKY_STILL_OPEN",
        "route": "strict_only",
        "legacy_bridge_used": False,
        "higher_reasoning_used": True,
        "pre_execution_grep_summary": {
            "english_terms": [
                "hAA",
                "graviton gauge",
                "metric perturbation",
                "sqrt(-g)",
                "F_{mu nu}",
                "physical-state projector",
            ],
            "polish_terms": [
                "wierzcholek/wierzcholk",
                "amplituda ubrana",
                "metryka/metryczny",
                "zaburzenie metryczne",
                "projektor stanu fizycznego",
            ],
            "key_existing_sources_found": [
                "P1657: covariant minimal gauge-metric density",
                "P1646: gauge-metric cross-variation scaffold",
                "P1714: stress-tensor split candidate",
                "P1954: dressed amplitude nonavailability theorem",
            ],
        },
        "depends_on": {
            "p1657_present": "sector" in p1657,
            "p1907_present": "full_lagrangian_term_registry_non_skeleton" in p1907,
            "p1954_present": "formal_nonavailability_theorem" in p1954,
        },
        "input_hashes": {
            "p1657_sha256": digest(p1657),
            "p1907_sha256": digest(p1907),
            "p1954_sha256": digest(p1954),
        },
        "metric_perturbation_convention": {
            "background_metric": "eta_mu_nu = diag(-1,1,1,1)",
            "metric_expansion": "g_mu_nu = eta_mu_nu + kappa*h_mu_nu",
            "inverse_metric_expansion": "g^mu_nu = eta^mu_nu - kappa*h^mu_nu + O(kappa^2)",
            "density_expansion": "sqrt(-g) = 1 + (kappa/2)*h + O(kappa^2)",
            "h_trace": "h = eta^{mu nu} h_mu_nu",
        },
        "source_lagrangian": {
            "from_p1657": "sqrt(-g)*[-1/4 g^{mu alpha}g^{nu beta}F_{mu nu}F_{alpha beta}]",
            "flat_lagrangian_symbolic": str(flat_lagrangian),
        },
        "V_hAA_tensor_strict_B1_v1": {
            "object_type": "minimal_tree_level_field_strength_vertex_density",
            "formula": "L_hAA_minimal_tree = (kappa/2) h^{mu nu} T^gauge_{mu nu}",
            "stress_tensor": "T^gauge_{mu nu}=Z_gauge*(F_{mu rho}F_nu^rho - (1/4) eta_{mu nu}F_{rho sigma}F^{rho sigma})",
            "linearity": "first order in kappa and quadratic in gauge field strength",
            "scheme_scope": "tree_level_minimal_gauge_metric_sector; no one-loop dressing",
        },
        "machine_check": {
            "linear_coefficient_minus_stress_tensor_vertex": str(residual),
            "stress_tensor_trace_residual": str(trace_residual),
            "stress_tensor_symmetry_residuals_all_zero": all(item == 0 for item in symmetry_residuals),
            "all_local_checks_zero": bool(all_zero),
        },
        "p1954_requirement_update": {
            "M1_4D_hAA_vertex": "DISCHARGED_FOR_MINIMAL_TREE_LEVEL_FIELD_STRENGTH_VERTEX",
            "remaining_M1_limitations": [
                "no R*F^2 or higher-curvature contact hAA terms",
                "no gauge-group normalization beyond a single Z_gauge delta_ab factor",
                "no momentum-space polarization tensor contraction",
            ],
            "M2_BRST_projector": "OPEN",
            "M3_dressed_propagator_residues": "OPEN",
            "M4_DiscM_CutSum_same_scheme": "OPEN",
            "M5_integral_reductions": "OPEN",
        },
        "false_pass_guard": "This exports only the minimal tree-level hAA vertex density; it is not a dressed amplitude, not a BRST projection, and not a Cutkosky closure theorem.",
        "next_honest_step": "Build P1956: attach a same-scheme BRST physical-state projector and polarization-sum certificate for the gauge_gauge final state, or export the exact missing ghost/projector data obstruction.",
        "lay_explanation": "Z istniejacego kowariantnego czlonu pola cechowania wyprowadzilismy pierwszy prawdziwy wierzcholek grawiton-dwa-bozony-cechowania. To usuwa jeden brak, ale pelny test unitarności nadal potrzebuje projektora BRST, ubranego propagatora i rownosci DiscM=CutSum.",
    }

    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(OUT)


if __name__ == "__main__":
    main()
