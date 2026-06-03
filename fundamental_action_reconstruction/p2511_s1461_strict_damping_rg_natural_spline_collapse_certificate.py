#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from typing import Any

import mpmath as mp
import sympy as sp

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import (
    DOC_FILES,
    REPO,
    ROOT,
    load_json,
    rel,
)
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2511_s1461_strict_damping_rg_natural_spline_collapse_certificate.json"
MD = GEN / "p2511_s1461_strict_damping_rg_natural_spline_collapse_certificate.md"

SOURCE_FILES = {
    "P2510_KKT_STATIONARITY": GEN / "p2510_s1460_strict_damping_rg_roughness_kkt_stationarity_certificate.json",
}

mp.mp.dps = 100
DOMAIN = list(range(1, 12))


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")).hexdigest()


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run(
        [
            "rg", "-n", pattern, "fundamental_action_reconstruction", "material_dowodowy",
            "-g", "*.py", "-g", "*.md", "-g", "*.tex", "-g", "!fundamental_action_reconstruction/generated/**",
        ],
        cwd=REPO,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:40]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2511|S1461|natural spline collapse|minimum-roughness spline|biharmonic spline collapse|tridiagonal spline certificate",
        "precursor_packets": "P2510|S1460|roughness KKT stationarity|P2509|minimum-roughness variational well-posedness|P2508|Sobolev node coercivity",
        "spline_language": "natural cubic spline|second-derivative knots|spline tridiagonal|piecewise cubic|interpolation slopes",
        "guardrails": "legacy -> strict completion bridge|role-transfer audit|K_legacy_ont|K_strict_gate|QW-2191|ToE closure",
        "closure_blockers": "source theorem|bridge theorem|role-transfer theorem|physical-value generator|role-bearing L_total|selector closure",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def text(value: mp.mpf, digits: int = 70) -> str:
    return mp.nstr(value, digits)


def symbolic_collinearity_audit() -> dict[str, Any]:
    delta = sp.Rational(14, 5) - 4 * sp.log(2)
    rows = []
    all_zero = True
    for d in range(1, 11):
        left = sp.log(d)
        right = sp.log(d + 1)
        y_left = delta * left
        y_right = delta * right
        slope_residual = sp.simplify((y_right - y_left) / (right - left) - delta)
        rows.append({
            "interval": f"[{d},{d + 1}]",
            "slope_residual": sp.sstr(slope_residual),
            "residual_zero": slope_residual == 0,
        })
        all_zero = all_zero and slope_residual == 0
    return {
        "symbolic_backend": "sympy",
        "sympy_version": sp.__version__,
        "delta": sp.sstr(delta),
        "theorem_step": "All node data y_d=delta*log(d) are exactly collinear in ell=log(d); therefore all divided-slope residuals equal zero.",
        "rows": rows,
        "all_interval_slope_residuals_zero": all_zero,
        "symbolic_fingerprint_sha256": sha256_json({"delta": sp.srepr(delta), "rows": rows}),
    }


def natural_spline_matrix_audit() -> dict[str, Any]:
    xs = [mp.log(d) for d in DOMAIN]
    delta = mp.mpf(14) / 5 - 4 * mp.log(2)
    ys = [delta * x for x in xs]
    hs = [xs[i + 1] - xs[i] for i in range(len(xs) - 1)]
    interior_count = len(xs) - 2
    matrix = mp.matrix(interior_count, interior_count)
    rhs = mp.matrix(interior_count, 1)
    slope_rows = []
    for i in range(1, len(xs) - 1):
        left_h = hs[i - 1]
        right_h = hs[i]
        row = i - 1
        if row - 1 >= 0:
            matrix[row, row - 1] = left_h
        matrix[row, row] = 2 * (left_h + right_h)
        if row + 1 < interior_count:
            matrix[row, row + 1] = right_h
        left_slope = (ys[i] - ys[i - 1]) / left_h
        right_slope = (ys[i + 1] - ys[i]) / right_h
        rhs[row] = 6 * (right_slope - left_slope)
        slope_rows.append({
            "interior_node_d": DOMAIN[i],
            "left_h": text(left_h, 50),
            "right_h": text(right_h, 50),
            "left_slope": text(left_slope, 50),
            "right_slope": text(right_slope, 50),
            "rhs": text(rhs[row], 50),
        })
    solution = mp.lu_solve(matrix, rhs)
    residual = matrix * solution - rhs
    pivots = []
    leading_minors = []
    for k in range(1, interior_count + 1):
        sub = mp.matrix(k, k)
        for i in range(k):
            for j in range(k):
                sub[i, j] = matrix[i, j]
        det = mp.det(sub)
        leading_minors.append(det)
        if k == 1:
            pivots.append(det)
        else:
            pivots.append(det / leading_minors[k - 2])
    full_second_derivatives = [mp.mpf("0")] + [solution[i] for i in range(interior_count)] + [mp.mpf("0")]
    piecewise_rows = []
    max_interp_residual = mp.mpf("0")
    max_quadratic_coeff = mp.mpf("0")
    max_cubic_coeff = mp.mpf("0")
    max_slope_error = mp.mpf("0")
    for i in range(len(hs)):
        h = hs[i]
        a = ys[i]
        b = (ys[i + 1] - ys[i]) / h - h * (2 * full_second_derivatives[i] + full_second_derivatives[i + 1]) / 6
        c = full_second_derivatives[i] / 2
        cubic = (full_second_derivatives[i + 1] - full_second_derivatives[i]) / (6 * h)
        left_value = a
        right_value = a + b * h + c * h ** 2 + cubic * h ** 3
        max_interp_residual = max(max_interp_residual, abs(left_value - ys[i]), abs(right_value - ys[i + 1]))
        max_quadratic_coeff = max(max_quadratic_coeff, abs(c))
        max_cubic_coeff = max(max_cubic_coeff, abs(cubic))
        max_slope_error = max(max_slope_error, abs(b - delta))
        piecewise_rows.append({
            "interval_d_to_d_plus_1": f"{DOMAIN[i]}->{DOMAIN[i + 1]}",
            "linear_coefficient_b": text(b, 50),
            "quadratic_coefficient_c": text(c, 50),
            "cubic_coefficient": text(cubic, 50),
            "right_endpoint_residual": text(right_value - ys[i + 1], 50),
        })
    energy = mp.fsum(
        (full_second_derivatives[i] ** 2 + full_second_derivatives[i] * full_second_derivatives[i + 1] + full_second_derivatives[i + 1] ** 2) * hs[i] / 3
        for i in range(len(hs))
    )
    max_rhs = max(abs(rhs[i]) for i in range(interior_count))
    max_residual = max(abs(residual[i]) for i in range(interior_count))
    max_second_derivative = max(abs(value) for value in full_second_derivatives)
    return {
        "node_count": len(xs),
        "interior_second_derivative_unknown_count": interior_count,
        "matrix_kind": "natural cubic spline second-derivative tridiagonal system",
        "slope_rows": slope_rows,
        "max_abs_rhs": text(max_rhs, 60),
        "tridiagonal_solution_second_derivatives": [text(value, 40) for value in full_second_derivatives],
        "max_abs_second_derivative_knot": text(max_second_derivative, 60),
        "max_abs_linear_solve_residual": text(max_residual, 60),
        "leading_principal_minors": [text(value, 50) for value in leading_minors],
        "cholesky_equivalent_pivots": [text(value, 50) for value in pivots],
        "all_leading_principal_minors_positive": all(value > 0 for value in leading_minors),
        "all_cholesky_equivalent_pivots_positive": all(value > 0 for value in pivots),
        "piecewise_coefficients": piecewise_rows,
        "max_abs_interpolation_residual": text(max_interp_residual, 60),
        "max_abs_slope_error_vs_delta": text(max_slope_error, 60),
        "max_abs_quadratic_coefficient": text(max_quadratic_coeff, 60),
        "max_abs_cubic_coefficient": text(max_cubic_coeff, 60),
        "closed_form_roughness_energy": text(energy, 60),
        "spline_collapses_to_affine_delta_ell": max_rhs < mp.mpf("1e-90") and max_second_derivative < mp.mpf("1e-90") and max_quadratic_coeff < mp.mpf("1e-90") and max_cubic_coeff < mp.mpf("1e-90"),
    }


def build_spline_collapse_certificate(p2510: dict[str, Any]) -> dict[str, Any]:
    symbolic = symbolic_collinearity_audit()
    spline = natural_spline_matrix_audit()
    return {
        "frontier_atom_under_attack": "strict_damping_beta_eta_source",
        "p2510_kkt_stationarity_inherited": p2510.get("kkt_stationarity_confirmed_for_postulated_functional") is True,
        "selector_type": "conditional natural-cubic-spline collapse theorem for the postulated minimum-roughness selector",
        "symbolic_collinearity_audit": symbolic,
        "natural_spline_tridiagonal_audit": spline,
        "natural_spline_collapse_theorem_for_postulated_functional": symbolic["all_interval_slope_residuals_zero"] and spline["spline_collapses_to_affine_delta_ell"],
        "roughness_action_still_postulated_not_derived": True,
        "strict_damping_beta_eta_source_exported": False,
        "strict_dynamical_source_for_A_P_D_exported": False,
        "strict_phase_frequency_source_exported": False,
        "bridge_theorem_exported": False,
        "role_transfer_theorem_exported": False,
        "selector_closure_exported": False,
        "qw2191_discharged_by_this_certificate": False,
        "toe_closure_claimed": False,
    }


def append_doc_sections() -> None:
    eq_section = """
## P2511/S1461 strict damping RG natural spline collapse certificate

`P2511/S1461` recasts the P2509/P2510 postulated roughness selector as the classical natural-cubic-spline interpolation problem on nodes `ell_d=log(d)` with data `y_d=delta log(d)`.  Since every divided slope `(y_{d+1}-y_d)/(ell_{d+1}-ell_d)` is exactly `delta`, the natural-spline second-derivative tridiagonal system has zero right-hand side.  The tridiagonal matrix has positive leading principal minors and positive Cholesky-equivalent pivots, hence the only natural second-derivative knot vector is zero.  The resulting piecewise cubic has zero quadratic/cubic coefficients and collapses to `y0(ell)=delta ell`.

This is an independent finite spline-form theorem witness for the conditional selector, not a source theorem.  It does not derive the roughness action from strict nadsoliton dynamics and does not export `strict_damping_beta_eta_source`, a bridge theorem, role-transfer theorem, QW-2191 discharge, physical-value generator, or ToE closure.
"""
    lag_section = """
## P2511/S1461 natural spline collapse guard

`P2511/S1461` verifies that the postulated strict damping roughness minimizer can also be obtained through the natural-cubic-spline normal form: linear-in-`ell` node data give zero spline second-derivative knots and a collapse to `y0(ell)=delta ell`.  This strengthens the conditional variational bookkeeping but remains dependent on the roughness functional being supplied; it is not a nonlinear compression-flow source theorem and does not license a role-bearing `L_total` term.
"""
    append_once(DOC_FILES["equation_sheet"], "## P2511/S1461 strict damping RG natural spline collapse certificate", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2511/S1461 natural spline collapse guard", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    p2510 = theorem(sources["P2510_KKT_STATIONARITY"], "strict_damping_rg_roughness_kkt_stationarity_certificate")
    cert = build_spline_collapse_certificate(p2510)
    spline = cert["natural_spline_tridiagonal_audit"]
    theorem_export = {
        "theorem_name": "P2511_T1_strict_damping_rg_natural_spline_collapse_certificate",
        "audited_chain": ["P2506/S1456", "P2507/S1457", "P2508/S1458", "P2509/S1459", "P2510/S1460"],
        "strict_damping_rg_natural_spline_collapse_certificate": cert,
        "frontier_atom_under_attack": cert["frontier_atom_under_attack"],
        "p2510_kkt_stationarity_inherited": cert["p2510_kkt_stationarity_inherited"],
        "symbolic_slope_residuals_zero": cert["symbolic_collinearity_audit"]["all_interval_slope_residuals_zero"],
        "tridiagonal_leading_minors_positive": spline["all_leading_principal_minors_positive"],
        "tridiagonal_pivots_positive": spline["all_cholesky_equivalent_pivots_positive"],
        "max_abs_rhs": spline["max_abs_rhs"],
        "max_abs_second_derivative_knot": spline["max_abs_second_derivative_knot"],
        "max_abs_quadratic_coefficient": spline["max_abs_quadratic_coefficient"],
        "max_abs_cubic_coefficient": spline["max_abs_cubic_coefficient"],
        "closed_form_roughness_energy": spline["closed_form_roughness_energy"],
        "natural_spline_collapse_theorem_for_postulated_functional": cert["natural_spline_collapse_theorem_for_postulated_functional"],
        "roughness_action_still_postulated_not_derived": cert["roughness_action_still_postulated_not_derived"],
        "strict_damping_beta_eta_source_exported": False,
        "strict_dynamical_source_for_A_P_D_exported": False,
        "strict_phase_frequency_source_exported": False,
        "bridge_theorem_exported": False,
        "role_transfer_theorem_exported": False,
        "selector_closure_exported": False,
        "qw2191_discharged_by_this_certificate": False,
        "toe_closure_claimed": False,
        "not_licensed": [
            "P2511 proves natural-spline collapse only for the postulated P2506/P2509 roughness minimization problem.",
            "It does not derive the roughness action, delta, beta0, or strict_damping_beta_eta_source from strict nadsoliton dynamics.",
            "No bridge theorem, role-transfer theorem, QW-2191 discharge, role-bearing L_total, physical-value generator, or ToE closure is exported.",
        ],
        "next_honest_step": "The selector route now has well-posedness, KKT stationarity, and natural-spline collapse witnesses; the remaining open step is still a strict source derivation for the roughness action or equivalent damping-flow source.",
    }
    gatekeepers = {
        "rg_audit_performed": grep["tool"] == "rg",
        "p2510_kkt_inherited": theorem_export["p2510_kkt_stationarity_inherited"],
        "symbolic_slope_residuals_zero": theorem_export["symbolic_slope_residuals_zero"],
        "tridiagonal_spd_witnesses_positive": theorem_export["tridiagonal_leading_minors_positive"] and theorem_export["tridiagonal_pivots_positive"],
        "spline_collapse_confirmed": theorem_export["natural_spline_collapse_theorem_for_postulated_functional"],
        "source_not_exported": not theorem_export["strict_damping_beta_eta_source_exported"],
        "negative_controls_preserved": not any(theorem_export[key] for key in [
            "strict_damping_beta_eta_source_exported",
            "strict_dynamical_source_for_A_P_D_exported",
            "strict_phase_frequency_source_exported",
            "bridge_theorem_exported",
            "role_transfer_theorem_exported",
            "selector_closure_exported",
            "qw2191_discharged_by_this_certificate",
            "toe_closure_claimed",
        ]),
    }
    return {
        "packet_id": "P2511",
        "stage_id": "S1461",
        "status": "STRICT_DAMPING_RG_NATURAL_SPLINE_COLLAPSE_FOR_POSTULATED_FUNCTIONAL_NO_SOURCE_EXPORT_NO_BRIDGE_THEOREM_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_damping_rg_natural_spline_collapse_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_damping_rg_natural_spline_collapse_certificate"]["theorem_export"]
    lines = [
        "# P2511/S1461 strict damping RG natural spline collapse certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        f"- Frontier atom under attack: `{t['frontier_atom_under_attack']}`.",
        f"- P2510 KKT stationarity inherited: `{t['p2510_kkt_stationarity_inherited']}`.",
        f"- Symbolic slope residuals zero: `{t['symbolic_slope_residuals_zero']}`.",
        f"- Tridiagonal leading minors positive: `{t['tridiagonal_leading_minors_positive']}`.",
        f"- Tridiagonal pivots positive: `{t['tridiagonal_pivots_positive']}`.",
        f"- Max spline RHS: `{t['max_abs_rhs']}`.",
        f"- Max second-derivative knot: `{t['max_abs_second_derivative_knot']}`.",
        f"- Max quadratic coefficient: `{t['max_abs_quadratic_coefficient']}`.",
        f"- Max cubic coefficient: `{t['max_abs_cubic_coefficient']}`.",
        f"- Closed-form roughness energy: `{t['closed_form_roughness_energy']}`.",
        f"- Natural spline collapse theorem for postulated functional: `{t['natural_spline_collapse_theorem_for_postulated_functional']}`.",
        f"- Source theorem exported: `{t['strict_damping_beta_eta_source_exported']}`.",
        "",
        "## Negative controls",
        "",
        "This packet proves natural-cubic-spline collapse only for the postulated roughness selector. It does not derive the roughness action from strict dynamics and exports no strict source atom, bridge theorem, role-transfer theorem, selector/QW-2191 closure, role-bearing L_total term, physical-value generator, or ToE closure.",
        "",
        "## Fingerprint",
        "",
        f"`{payload['strict_damping_rg_natural_spline_collapse_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload["strict_damping_rg_natural_spline_collapse_certificate"]["theorem_export"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
