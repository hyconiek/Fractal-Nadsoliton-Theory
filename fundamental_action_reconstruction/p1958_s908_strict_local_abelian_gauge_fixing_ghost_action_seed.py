#!/usr/bin/env python3
"""P1958 S908 strict local abelian gauge-fixing and ghost-action seed.

This executor takes the explicit minimal gauge sector already present in the
strict repository state and exports the smallest honest gauge-fixing/ghost
object that can be derived without adding unexported non-Abelian structure:

  G[A] = partial_mu A^mu
  delta A_mu = partial_mu alpha
  M_FP = delta G[A^alpha] / delta alpha = Box
  L_GF = -(1/(2 xi)) G[A]^2
  L_ghost = - cbar Box c  (equivalently partial_mu cbar partial^mu c up to boundary)

It is not a full BV/BRST operator map, not a Q^2=0 certificate, and not a
ghost-cancellation theorem for the dressed gauge_gauge Cutkosky channel.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p1958_s908_strict_local_abelian_gauge_fixing_ghost_action_seed.json"


def load(name: str) -> dict[str, Any]:
    path = GEN / name
    if not path.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def digest(obj: object) -> str:
    blob = json.dumps(obj, sort_keys=True, ensure_ascii=True).encode("utf-8")
    return hashlib.sha256(blob).hexdigest()


def d(expr: sp.Expr, coord: sp.Symbol) -> sp.Expr:
    return sp.diff(expr, coord)


def box(expr: sp.Expr, coords: list[sp.Symbol], eta_diag: list[int]) -> sp.Expr:
    return sp.simplify(sum(eta_diag[i] * d(d(expr, coords[i]), coords[i]) for i in range(4)))


def matrix_to_strings(items: list[list[sp.Expr]]) -> list[list[str]]:
    return [[str(sp.simplify(item)) for item in row] for row in items]


def all_zero(items: list[sp.Expr]) -> bool:
    return all(sp.simplify(item) == 0 for item in items)


def main() -> None:
    p1657 = load("p1657_s607_strict_helmholtz_h1_gauge_metric_covariant_summary.json")
    p1764 = load("p1764_s714_strict_covariant_nonproxy_ea_eh_explicit_export_checkpoint.json")
    p1907 = load("p1907_s857_strict_full_lagrangian_to_eom_witness_matrix_probe.json")
    p1957 = load("p1957_s907_strict_bv_brst_ghost_sector_nonavailability_theorem.json")

    t, x, y, z, eps, xi = sp.symbols("t x y z eps xi", nonzero=True)
    coords = [t, x, y, z]
    eta_diag = [-1, 1, 1, 1]

    a = [sp.Function(f"A{mu}")(*coords) for mu in range(4)]  # lower A_mu
    alpha = sp.Function("alpha")(*coords)
    c = sp.Function("c")(*coords)
    cbar = sp.Function("cbar")(*coords)

    # Lower-index Abelian gauge transformation: delta A_mu = partial_mu alpha.
    delta_a = [d(alpha, coord) for coord in coords]
    a_transformed = [a[mu] + eps * delta_a[mu] for mu in range(4)]

    def lorenz_gauge(field: list[sp.Expr]) -> sp.Expr:
        # G[A] = partial_mu A^mu = partial_mu(eta^{mu mu} A_mu).
        return sp.simplify(sum(eta_diag[mu] * d(field[mu], coords[mu]) for mu in range(4)))

    def field_strength(field: list[sp.Expr]) -> list[list[sp.Expr]]:
        return [
            [sp.simplify(d(field[nu], coords[mu]) - d(field[mu], coords[nu])) for nu in range(4)]
            for mu in range(4)
        ]

    gauge_functional = lorenz_gauge(a)
    transformed_gauge_functional = lorenz_gauge(a_transformed)
    fp_operator_on_alpha = sp.simplify(sp.diff(transformed_gauge_functional, eps).subs(eps, 0))
    box_alpha = box(alpha, coords, eta_diag)
    fp_residual = sp.simplify(fp_operator_on_alpha - box_alpha)

    f_original = field_strength(a)
    f_transformed = field_strength(a_transformed)
    f_variation_residuals = [
        sp.simplify(sp.diff(f_transformed[mu][nu], eps).subs(eps, 0))
        for mu in range(4)
        for nu in range(4)
    ]

    l_gf = sp.simplify(-sp.Rational(1, 2) * gauge_functional**2 / xi)
    l_ghost_operator = sp.simplify(-cbar * box(c, coords, eta_diag))
    l_ghost_by_parts = sp.simplify(
        sum(eta_diag[mu] * d(cbar, coords[mu]) * d(c, coords[mu]) for mu in range(4))
    )
    # The by-parts form differs from -cbar Box c by a total derivative:
    # -cbar Box c = partial_mu(-cbar partial^mu c) + partial_mu cbar partial^mu c.
    boundary_current = [-cbar * eta_diag[mu] * d(c, coords[mu]) for mu in range(4)]
    boundary_divergence = sp.simplify(sum(d(boundary_current[mu], coords[mu]) for mu in range(4)))
    ghost_by_parts_residual = sp.simplify(l_ghost_operator - (l_ghost_by_parts + boundary_divergence))

    checks_zero = [
        fp_residual,
        ghost_by_parts_residual,
        *f_variation_residuals,
    ]

    out = {
        "packet_id": "P1958",
        "stage_id": "S908",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "local_verdict": "LOCAL_ABELIAN_LORENZ_GAUGE_FIXING_AND_GHOST_ACTION_SEED_EXPORTED__FULL_BV_BRST_STILL_OPEN",
        "route": "strict_only",
        "legacy_bridge_used": False,
        "higher_reasoning_used": True,
        "ignored_files": ["TOE_FINAL_DOCUMENTATION_RELEASE_8_STRICT_FULL.pdf"],
        "pre_execution_grep_summary": {
            "english_terms": [
                "gauge sector",
                "L_gauge",
                "F_{mu nu}",
                "Yang/Mills",
                "structure constants",
                "gauge fixing",
                "Faddeev/Popov",
                "ghost action",
                "antighost",
                "BRST differential",
                "gauge transformation",
            ],
            "polish_terms": [
                "sektor cechowania",
                "cechowanie",
                "ustalenie cechowania",
                "duchy",
                "antyduch",
                "akcja duchow",
                "transformacja cechowania",
                "pochodna kowariantna",
            ],
            "key_existing_sources_found": [
                "P1657 exports the minimal covariant gauge-metric density sqrt(-g)*[-1/4 g^{mu alpha}g^{nu beta}F_{mu nu}F_{alpha beta}] and notes gauge fixing is not discharged.",
                "P1907 exports a sector-level SM gauge registry: -1/4 F^2_SU3, -1/4 F^2_SU2, -1/4 F^2_U1, with ghost/BRST constraints still OPEN_SYMBOLIC.",
                "P1764 exports covariant operator-level E_A^mu with ghost/BRST still open.",
                "P1957 proves that full BV/BRST ghost-sector data are not available in the current export state.",
            ],
            "negative_search_result": (
                "No theorem-grade structure-constant table, BV_BRST_operator_map, "
                "or nonproxy ghost/antighost action was found."
            ),
        },
        "depends_on": {
            "p1657_present": "sector" in p1657,
            "p1764_present": "d1_d2_upgrade" in p1764,
            "p1907_present": "full_lagrangian_term_registry_non_skeleton" in p1907,
            "p1957_present": "formal_nonavailability_theorem" in p1957,
        },
        "input_hashes": {
            "p1657_sha256": digest(p1657),
            "p1764_sha256": digest(p1764),
            "p1907_sha256": digest(p1907),
            "p1957_sha256": digest(p1957),
        },
        "scope": {
            "field_content": "single local Abelianized gauge generator over flat B1 tangent patch",
            "metric_signature": "eta_mu_nu = diag(-1,1,1,1)",
            "gauge_potential_convention": "A_mu lower-index components",
            "gauge_transformation": "delta A_mu = partial_mu alpha",
            "not_in_scope": [
                "SU(3)/SU(2) non-Abelian structure constants",
                "diffeomorphism ghosts",
                "BV antifield sector",
                "BRST charge Q",
                "Q^2=0 theorem",
                "ghost-cancelled Cutkosky equality",
            ],
        },
        "GaugeFixingFunctional_strict_B1_v1": {
            "status": "EXPORTED_LOCAL_FREE_ABELIANIZED_ONLY",
            "formula": "G[A] = partial_mu A^mu = -partial_t A_0 + partial_x A_1 + partial_y A_2 + partial_z A_3",
            "gauge_fixing_lagrangian": "L_GF = -(1/(2*xi))*G[A]^2",
            "symbolic_expression": str(gauge_functional),
            "gauge_parameter": "xi",
        },
        "FaddeevPopovOperator_strict_B1_v1": {
            "status": "EXPORTED_LOCAL_FREE_ABELIANIZED_ONLY",
            "derivation": "delta G[A + eps*partial alpha]/delta eps at eps=0",
            "operator_on_alpha": str(fp_operator_on_alpha),
            "box_alpha": str(box_alpha),
            "residual_operator_minus_box": str(fp_residual),
        },
        "GhostAntighostAction_strict_B1_v1": {
            "status": "EXPORTED_LOCAL_FREE_ABELIANIZED_ONLY",
            "operator_form": "L_ghost = - cbar * Box c",
            "by_parts_form": "L_ghost = partial_mu cbar partial^mu c modulo boundary",
            "operator_expression": str(l_ghost_operator),
            "by_parts_expression": str(l_ghost_by_parts),
            "boundary_current_J_mu": [str(sp.simplify(item)) for item in boundary_current],
            "operator_minus_byparts_minus_boundary_residual": str(ghost_by_parts_residual),
        },
        "machine_checks": {
            "abelian_field_strength_variation_zero": all_zero(f_variation_residuals),
            "fp_operator_equals_box": fp_residual == 0,
            "ghost_operator_equals_byparts_plus_boundary": ghost_by_parts_residual == 0,
            "all_local_checks_zero": all_zero(checks_zero),
            "field_strength_variation_residuals": matrix_to_strings(
                [
                    [sp.simplify(sp.diff(f_transformed[mu][nu], eps).subs(eps, 0)) for nu in range(4)]
                    for mu in range(4)
                ]
            ),
        },
        "p1957_requirement_update": {
            "B1_ghost_sector_nonproxy_export": "PARTIAL_SEED_ONLY_FOR_LOCAL_FREE_ABELIANIZED_GAUGE_FIELD",
            "B2_BV_BRST_operator_map": "OPEN",
            "B3_explicit_BRST_charge_Q": "OPEN",
            "B4_Q_squared_simplified_zero": "OPEN",
            "B5_physical_state_cohomology": "OPEN",
            "B6_gauge_gauge_ghost_cancellation_trace": "OPEN",
        },
        "gatekeeper_status": {
            "TG2_BRST_GLOBAL_NILPOTENCY": "NOT_PROMOTED",
            "TG3_CUTKOSKY_GLOBAL_UNITARITY": "NOT_PROMOTED",
            "reason": (
                "P1958 exports a local Abelian gauge-fixing/ghost seed only. "
                "It does not export the non-Abelian SM ghost sector, BV map, Q, Q^2, "
                "cohomology map, or ghost-cancelled Cutkosky trace."
            ),
        },
        "false_pass_guard": (
            "This packet discharges only the minimal local Abelian FP operator and ghost "
            "action seed derivable from the exported F^2 gauge sector. It is not a full "
            "strict BV/BRST witness pack and cannot close TG2/TG3."
        ),
        "next_honest_step": (
            "Build P1959 with high reasoning: decide whether the strict lane has enough "
            "exported gauge-group data to extend this seed to SU(3)xSU(2)xU(1) structure "
            "constants and BRST differential, or prove a non-Abelian structure-data obstruction."
        ),
        "lay_explanation": (
            "Zbudowalismy pierwszy maly, sprawdzalny fragment sektora duchow: dla prostego "
            "pola cechowania warunek Lorentza prowadzi do operatora Box i odpowiadajacej "
            "akcji duchow. To nie jest jeszcze pelny BRST Standard Modelu; to pierwszy "
            "uczciwy klocek, ktory pokazuje, gdzie zaczyna sie rachunek duchow."
        ),
    }

    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(OUT)


if __name__ == "__main__":
    main()
