#!/usr/bin/env python3
"""P1970 S920 strict Higgs/Yukawa/curvature variational normal-form attempt.

This is the next honest step after P1969.  It does not claim full PO2 closure.
It freezes a narrow scalar-Higgs proxy convention for the Higgs/Yukawa/nonminimal
curvature subsector, derives an Euler-Lagrange row symbolically, and tests
whether the P1965 leftover Omega_unexported is forced to zero.

Result: a useful variational row is exported, but Omega is not forced to zero
without an additional branch-lock/invertibility condition for the Higgs
background.  This sharpens P1965 rather than closing it falsely.
"""

from __future__ import annotations

import hashlib
import json
import platform
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p1970_s920_strict_higgs_yukawa_curvature_variational_normal_form_attempt.json"


def load_generated(name: str) -> dict[str, Any]:
    path = GEN / name
    if not path.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def digest(obj: object) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True).encode("utf-8")).hexdigest()


def variational_row() -> dict[str, Any]:
    x = sp.symbols("x", real=True)
    mu2, lam, xi, y = sp.symbols("mu_H_sq lambda_H xi_H y_f", real=True)
    h = sp.Function("h")(x)
    R = sp.Function("R")(x)
    J = sp.Function("J_yukawa")(x)

    # Frozen scalar-Higgs proxy for the local subsector.  The 1D coordinate is
    # only a symbolic carrier for the integration-by-parts Euler operator; the
    # exported dictionary records the covariant replacements.
    lagrangian = (
        sp.Rational(1, 2) * sp.diff(h, x) ** 2
        - sp.Rational(1, 2) * mu2 * h**2
        - sp.Rational(1, 4) * lam * h**4
        - sp.Rational(1, 2) * xi * R * h**2
        - y * J * h
    )
    d_l_d_h = sp.diff(lagrangian, h)
    d_l_d_hx = sp.diff(lagrangian, sp.diff(h, x))
    euler_operator = sp.simplify(sp.diff(d_l_d_hx, x) - d_l_d_h)

    covariant_eom = "Box H + mu_H_sq*H + lambda_H*(H^dagger H)H + xi_H*R*H + y_f*J_yukawa = 0"

    return {
        "coordinate": "x (symbolic integration-by-parts carrier only)",
        "lagrangian_proxy": sp.sstr(lagrangian),
        "dL_dh": sp.sstr(d_l_d_h),
        "dL_dh_x": sp.sstr(d_l_d_hx),
        "euler_lagrange_operator_d_dx_dL_dhx_minus_dL_dh": sp.sstr(euler_operator),
        "covariant_dictionary": {
            "Derivative(h(x), (x, 2))": "Box H",
            "h(x)": "H radial/background component or local H representative",
            "R(x)": "R",
            "J_yukawa(x)": "Yukawa source bilinear qbar/f source in frozen channel",
        },
        "covariant_eom_row_exported": covariant_eom,
        "scope_limit": "single scalar-Higgs/Yukawa/nonminimal-curvature proxy row; not full SM spinor/gauge/metric EOM",
    }


def branch_difference_and_leftover() -> dict[str, Any]:
    hA, hB, boxA, boxB, RA, RB, JA, JB = sp.symbols("H_A H_B BoxH_A BoxH_B R_A R_B J_A J_B")
    mu2, lam, xi, y = sp.symbols("mu_H_sq lambda_H xi_H y_f")
    delta_H, delta_BoxH, delta_R, delta_J = sp.symbols("delta_H delta_BoxH delta_R delta_J")

    eom_A = boxA + mu2 * hA + lam * hA**3 + xi * RA * hA + y * JA
    eom_B = boxB + mu2 * hB + lam * hB**3 + xi * RB * hB + y * JB
    diff_raw = sp.expand(eom_A - eom_B)
    substitutions = {
        hA: hB + delta_H,
        boxA: boxB + delta_BoxH,
        RA: RB + delta_R,
        JA: JB + delta_J,
    }
    diff_delta = sp.expand(diff_raw.subs(substitutions))
    # C1-C3-like background equality removes curvature/source/laplacian
    # differences, but not the Higgs branch value difference itself.
    c_constraints = {delta_R: 0, delta_J: 0, delta_BoxH: 0}
    diff_after_constraints = sp.factor(diff_delta.subs(c_constraints))

    delta_yukawa_background = y * delta_H

    # Concrete underdetermination witness: in the massless flat no-source local
    # proxy row, every constant H is stationary, so two branches can satisfy the
    # row while their Yukawa background differs.  This is not a physical theorem;
    # it is a proof-obligation witness showing that a branch-lock/invertibility
    # premise is required before full PO2 promotion.
    a, b = sp.symbols("a b")
    witness_subs = {
        mu2: 0,
        lam: 0,
        xi: 0,
        y: sp.symbols("y_nonzero"),
        hA: a,
        hB: b,
        boxA: 0,
        boxB: 0,
        RA: 0,
        RB: 0,
        JA: 0,
        JB: 0,
    }
    eom_A_witness = sp.simplify(eom_A.subs(witness_subs))
    eom_B_witness = sp.simplify(eom_B.subs(witness_subs))
    delta_yuk_witness = sp.simplify((sp.symbols("y_nonzero") * (a - b)))

    return {
        "eom_A_minus_eom_B_raw": sp.sstr(diff_raw),
        "delta_substitution": {str(k): sp.sstr(v) for k, v in substitutions.items()},
        "eom_difference_in_delta_variables": sp.sstr(diff_delta),
        "constraints_applied": {
            "delta_R": "0",
            "delta_J": "0",
            "delta_BoxH": "0",
            "note": "These are narrower scalar-row analogues of equal-curvature/equal-source/equal-kinetic branch constraints; they do not include delta_H=0.",
        },
        "eom_difference_after_constraints": sp.sstr(diff_after_constraints),
        "delta_yukawa_background": sp.sstr(delta_yukawa_background),
        "omega_leftover_candidate": "delta_H (or equivalently y_f*delta_H in the Yukawa background)",
        "omega_forced_zero_by_this_row": bool(diff_after_constraints == 0 and delta_yukawa_background == 0),
        "sufficient_extra_conditions_for_zero": [
            "delta_H = 0 same-Higgs-background branch lock",
            "or an exported invertibility theorem proving the only same-source/same-curvature solution has delta_H=0 on the admissible branch",
        ],
        "underdetermination_witness": {
            "parameter_choice": "mu_H_sq=lambda_H=xi_H=0, R_A=R_B=J_A=J_B=BoxH_A=BoxH_B=0, y_f=y_nonzero",
            "H_A": "a",
            "H_B": "b",
            "EOM_A": sp.sstr(eom_A_witness),
            "EOM_B": sp.sstr(eom_B_witness),
            "Delta_Yukawa_background": sp.sstr(delta_yuk_witness),
            "nonzero_if": "y_nonzero != 0 and a != b",
            "meaning": "Current scalar-row variational data alone does not force the Yukawa background difference to vanish.",
        },
    }


def main() -> None:
    GEN.mkdir(parents=True, exist_ok=True)
    p1907 = load_generated("p1907_s857_strict_full_lagrangian_to_eom_witness_matrix_probe.json")
    p1930 = load_generated("p1930_s880_strict_b1_invariant_triplet_branch_evaluation_probe.json")
    p1964 = load_generated("p1964_s914_strict_po2_conditional_sufficiency_and_eom_gap_certificate.json")
    p1965 = load_generated("p1965_s915_strict_delta_bg_yf_eom_normal_form_extraction_nonavailability.json")
    p1969 = load_generated("p1969_s919_strict_toe_progress_dashboard_after_selector_axis_chain.json")

    row = variational_row()
    leftover = branch_difference_and_leftover()
    omega_forced_zero = leftover["omega_forced_zero_by_this_row"]

    out = {
        "packet_id": "P1970",
        "stage_id": "S920",
        "status": "PARTIAL_VARIATIONAL_ROW_EXPORTED__HIGGS_BRANCH_LEFTOVER_NOT_FORCED_ZERO__FULL_PO2_STILL_OPEN",
        "route": "strict_only_higgs_yukawa_curvature_subsector",
        "legacy_bridge_used": False,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "python": platform.python_version(),
        "sympy": sp.__version__,
        "depends_on": {
            "p1907_present": "full_lagrangian_term_registry_non_skeleton" in p1907,
            "p1930_present": "tensorial_b1_witness_form" in p1930,
            "p1964_conditional_po2_pass": p1964.get("symbolic_po2_sufficiency_check", {}).get("sympy_zero_check") is True,
            "p1965_nonavailability_status": p1965.get("status"),
            "p1969_recommended_next": p1969.get("professorial_decision", {}).get("recommended_next_packet"),
        },
        "input_hashes": {
            "p1907_sha256": digest(p1907),
            "p1930_sha256": digest(p1930),
            "p1964_sha256": digest(p1964),
            "p1965_sha256": digest(p1965),
            "p1969_sha256": digest(p1969),
        },
        "frozen_conventions": {
            "density": "sqrt(-g) factored out locally; scalar proxy uses integration-by-parts carrier x",
            "higgs_proxy": "single real H representative for H^dagger H channel; full doublet/gauge/spinor indices not claimed",
            "euler_operator": "d/dx(dL/dh_x) - dL/dh = 0, mapped to covariant Box H row",
            "yukawa_source": "J_yukawa is an external bilinear source placeholder for this narrow row",
            "nonminimal_curvature": "-1/2 xi_H R h^2 proxy for -xi_H H^dagger H R up to representation normalization",
        },
        "variational_row_export": row,
        "branch_difference_and_omega_test": leftover,
        "p1965_requirement_update": {
            "E1_full_non_skeleton_L_total_density": "PARTIAL_SUBSECTOR_ONLY_FROM_P1907_REGISTRY",
            "E2_convention_fixed_covariant_field_calculus": "PARTIAL_FROZEN_FOR_SCALAR_HIGGS_PROXY",
            "E3_termwise_euler_lagrange_derivatives": "PARTIAL_PASS_HIGGS_YUKAWA_CURVATURE_ROW_ONLY",
            "E4_projection_to_DELTA_BG_Yf_basis": "OPEN_DELTA_H_BRANCH_LEFTOVER_NOT_IN_P1964_BASIS",
            "E5_normal_form_extraction_theorem": "OPEN_REQUIRES_BRANCH_LOCK_OR_INVERTIBILITY_THEOREM",
            "E6_same_branch_quantifier_transport": "OPEN",
        },
        "machine_checks": {
            "variational_row_exported": True,
            "omega_unexported_forced_zero": omega_forced_zero,
            "full_po2_closed": False,
            "sharp_new_blocker_identified": "Higgs background branch-lock/invertibility needed to eliminate delta_H leftover",
        },
        "theorem_export": {
            "positive_statement": "The Higgs/Yukawa/nonminimal-curvature scalar proxy row has an explicit Euler-Lagrange derivation with frozen conventions.",
            "negative_statement": "That row alone does not force the Yukawa-background leftover y_f*delta_H to vanish under curvature/source/kinetic equality constraints that omit delta_H=0.",
            "not_promoted_to": [
                "full PO2 sufficiency from L_total",
                "global background-independence closure",
                "full strict-core ToE closure",
                "QW-2191 discharge",
            ],
        },
        "false_pass_guard": "P1970 is a constructive partial variational extraction and a sharper obstruction witness, not PO2 closure.",
        "next_honest_step": "Build P1971: prove or refute a strict branch-lock/invertibility lemma for the Higgs background on the admissible PO3 domain; if unavailable, export the exact delta_H leftover as the PO2 normal-form obstruction.",
        "lay_explanation": "Udało się wykonać kawałek rachunku wariacyjnego dla sektora Higgsa/Yukawy/krzywizny. Ten rachunek pokazuje jednak, że same równania w tej wąskiej postaci nie zmuszają dwóch gałęzi do posiadania tego samego tła Higgsa. Bez tego może pozostać różnica Yukawy, więc pełne domknięcie PO2 nadal byłoby fałszywym zaliczeniem.",
    }
    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(OUT)


if __name__ == "__main__":
    main()
