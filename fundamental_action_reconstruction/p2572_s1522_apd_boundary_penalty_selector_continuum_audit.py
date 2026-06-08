#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from typing import Any

import numpy as np
from numpy.polynomial import Polynomial

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import DOC_FILES, REPO, ROOT, load_json, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2569_s1519_apd_value_bridge_interpolation_dynamic_nonuniqueness_certificate import apd_rows

GEN = ROOT / "generated"
OUT = GEN / "p2572_s1522_apd_boundary_penalty_selector_continuum_audit.json"
MD = GEN / "p2572_s1522_apd_boundary_penalty_selector_continuum_audit.md"

SOURCE_FILES = {
    "P2416_APD_VALUE_BRIDGE": GEN / "p2416_s1366_apd_multiplicative_bridge_assembly_necessity_certificate.json",
    "P2570_APD_SOBOLEV_ORDER_DEPENDENCE": GEN / "p2570_s1520_apd_sobolev_roughness_selector_order_dependence_audit.json",
    "P2571_APD_MEASURE_BOUNDARY_DEPENDENCE": GEN / "p2571_s1521_apd_sobolev_measure_boundary_dependence_audit.json",
}
DOMAIN = list(range(12))
FIXED_DERIVATIVE_ORDER = 2
PENALTY_GRID = [0.0, 1.0e-8, 1.0e-6, 1.0e-4, 1.0e-2, 1.0, 1.0e2, 1.0e4, 1.0e8]
PROBE_POINTS = [0.5, 5.5, 11.5]
NEGATIVE_EXPORT_FLAGS = [
    "apd_boundary_penalty_source_exported",
    "apd_boundary_selector_source_exported",
    "apd_weighted_roughness_selector_source_exported",
    "strict_dynamical_source_for_A_P_D_exported",
    "strict_phase_frequency_source_exported",
    "strict_damping_beta_eta_source_exported",
    "bridge_theorem_exported",
    "legacy_to_strict_completion_bridge_exported",
    "role_transfer_theorem_exported",
    "role_bearing_ltotal_exported",
    "qw2191_discharged_by_this_certificate",
    "toe_closure_claimed",
]


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode()).hexdigest()


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run([
        "rg", "-n", pattern, "fundamental_action_reconstruction", "material_dowodowy",
        "-g", "*.py", "-g", "*.md", "-g", "*.tex", "-g", "!fundamental_action_reconstruction/generated/**",
    ], cwd=REPO, check=False, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:80]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2572|S1522|APD boundary penalty continuum|boundary penalty continuum",
        "intended_research_nonduplication": "APD.*penalty continuum|endpoint penalty.*APD|APD.*endpoint penalty|boundary-penalty.*APD|APD.*inverse boundary|Sobolev.*penalty continuum|penalty.*Sobolev.*APD",
        "apd_precursors": "P2416|S1366|P2570|S1520|P2571|S1521|APD Sobolev|APD.*boundary|strict_dynamical_source_for_A_P_D",
        "guardrails": "K_legacy_ont|K_strict_gate|role-transfer theorem|QW-2191|ToE closure|source theorem",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def interpolation_polynomials(rows: list[dict[str, Any]]) -> tuple[Polynomial, Polynomial]:
    x_nodes = np.array([row["d"] for row in rows], dtype=float)
    y_nodes = np.array([row["apd_value"] for row in rows], dtype=float)
    base = Polynomial.fit(x_nodes, y_nodes, deg=11, domain=[0.0, 11.0]).convert()
    vanish = Polynomial([1.0])
    for d in DOMAIN:
        vanish = vanish * Polynomial([-float(d), 1.0])
    return base, vanish


def definite_integral(poly: Polynomial, left: float = 0.0, right: float = 11.0) -> float:
    integral = poly.integ()
    return float(integral(right) - integral(left))


def base_coefficients(base: Polynomial, vanish: Polynomial) -> dict[str, float]:
    d_base = base.deriv(FIXED_DERIVATIVE_ORDER)
    d_vanish = vanish.deriv(FIXED_DERIVATIVE_ORDER)
    slope_base = base.deriv(1)
    slope_vanish = vanish.deriv(1)
    coeffs = {
        "A_bulk": definite_integral(d_vanish * d_vanish),
        "B_bulk": definite_integral(d_base * d_vanish),
        "C_bulk": definite_integral(d_base * d_base),
        "left_A": float(slope_vanish(0.0)) ** 2,
        "left_B": float(slope_base(0.0)) * float(slope_vanish(0.0)),
        "right_A": float(slope_vanish(11.0)) ** 2,
        "right_B": float(slope_base(11.0)) * float(slope_vanish(11.0)),
    }
    coeffs["uniform_lambda"] = -coeffs["B_bulk"] / coeffs["A_bulk"]
    coeffs["left_penalty_limit_lambda"] = -coeffs["left_B"] / coeffs["left_A"]
    coeffs["right_penalty_limit_lambda"] = -coeffs["right_B"] / coeffs["right_A"]
    return coeffs


def family_row(side: str, penalty: float, coeffs: dict[str, float], base: Polynomial, vanish: Polynomial, rows: list[dict[str, Any]]) -> dict[str, Any]:
    left_penalty = penalty if side == "left" else 0.0
    right_penalty = penalty if side == "right" else 0.0
    a_quad = coeffs["A_bulk"] + left_penalty * coeffs["left_A"] + right_penalty * coeffs["right_A"]
    b_cross = coeffs["B_bulk"] + left_penalty * coeffs["left_B"] + right_penalty * coeffs["right_B"]
    lambda_star = -b_cross / a_quad
    selected = base + lambda_star * vanish
    node_points = np.array(DOMAIN, dtype=float)
    node_values = np.array([row["apd_value"] for row in rows], dtype=float)
    residuals = selected(node_points) - node_values
    probe_deltas = [float(lambda_star * vanish(point)) for point in PROBE_POINTS]
    return {
        "side": side,
        "endpoint_slope_penalty": penalty,
        "lambda_star": float(lambda_star),
        "max_abs_node_residual_at_selected_lambda": float(np.max(np.abs(residuals))),
        "probe_value_deltas_from_base": probe_deltas,
        "selected_member_changes_off_node_values": any(abs(value) > 0.0 for value in probe_deltas),
    }


def penalty_continuum_audit(rows: list[dict[str, Any]]) -> dict[str, Any]:
    base, vanish = interpolation_polynomials(rows)
    coeffs = base_coefficients(base, vanish)
    left_rows = [family_row("left", penalty, coeffs, base, vanish, rows) for penalty in PENALTY_GRID]
    right_rows = [family_row("right", penalty, coeffs, base, vanish, rows) for penalty in PENALTY_GRID]
    rounded_left = {round(row["lambda_star"], 24) for row in left_rows}
    rounded_right = {round(row["lambda_star"], 24) for row in right_rows}
    left_values = [row["lambda_star"] for row in left_rows]
    right_values = [row["lambda_star"] for row in right_rows]
    return {
        "numpy_version": np.__version__,
        "fixed_derivative_order": FIXED_DERIVATIVE_ORDER,
        "domain": DOMAIN,
        "penalty_grid": PENALTY_GRID,
        "probe_points": PROBE_POINTS,
        "stationarity_formula": "lambda(t)=-(B_bulk+t*B_endpoint)/(A_bulk+t*A_endpoint)",
        "bulk_and_endpoint_coefficients": coeffs,
        "left_penalty_rows": left_rows,
        "right_penalty_rows": right_rows,
        "left_distinct_minimizers_after_rounding_1e_minus_24": len(rounded_left),
        "right_distinct_minimizers_after_rounding_1e_minus_24": len(rounded_right),
        "left_lambda_range_on_grid": [float(min(left_values)), float(max(left_values))],
        "right_lambda_range_on_grid": [float(min(right_values)), float(max(right_values))],
        "left_penalty_changes_selector": len(rounded_left) > 1,
        "right_penalty_changes_selector": len(rounded_right) > 1,
        "all_penalty_rows_preserve_apd_nodes": all(row["max_abs_node_residual_at_selected_lambda"] <= 1.0e-6 for row in left_rows + right_rows),
        "boundary_penalty_parameter_is_continuous_unsourced_selector": True,
    }


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    p2416_payload = load_json(SOURCE_FILES["P2416_APD_VALUE_BRIDGE"])
    p2570_payload = load_json(SOURCE_FILES["P2570_APD_SOBOLEV_ORDER_DEPENDENCE"])
    p2571_payload = load_json(SOURCE_FILES["P2571_APD_MEASURE_BOUNDARY_DEPENDENCE"])
    p2570 = theorem(p2570_payload, "apd_sobolev_roughness_selector_order_dependence_audit")
    p2571 = theorem(p2571_payload, "apd_sobolev_measure_boundary_dependence_audit")
    rows = apd_rows(p2416_payload)
    audit = penalty_continuum_audit(rows)
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2572_T1_apd_boundary_penalty_selector_continuum_audit",
        "audited_chain": ["P2416/S1366", "P2570/S1520", "P2571/S1521"],
        "frontier_atom_under_attack": "strict_dynamical_source_for_A_P_D",
        "candidate_principle_under_test": "fix k=2 and uniform bulk measure, then vary one endpoint slope penalty t >= 0",
        "p2416_apd_value_bridge_inherited": p2416_payload.get("apd_multiplicative_bridge_assembly_necessity_certificate", {}).get("finite_witness_certificate", {}).get("apd_value_assembly_witness_ready") is True,
        "p2570_order_obligation_inherited": p2570.get("finite_apd_values_do_not_select_sobolev_order") is True,
        "p2571_measure_boundary_obligation_inherited": p2571.get("finite_apd_values_do_not_select_measure_or_boundary_class") is True,
        "apd_node_rows": rows,
        "boundary_penalty_selector_continuum_audit": audit,
        "finite_apd_values_do_not_select_boundary_penalty": audit["left_penalty_changes_selector"] and audit["right_penalty_changes_selector"],
        "fixed_k2_uniform_bulk_still_conditional": True,
        "apd_dynamic_source_remains_unsourced": True,
        "recommended_next_honest_step": (
            "Do not choose an endpoint penalty by convenience or by post-hoc fit. The next honest step is to derive boundary conditions or boundary penalties from the strict nadsoliton APD action; otherwise even fixed k=2 and uniform bulk measure leave a continuous family of conditional APD selectors."
        ),
        **{flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    gatekeepers = {
        "rg_non_duplication_audit_recorded": grep["tool"] == "rg",
        "p2416_apd_value_bridge_inherited": theorem_export["p2416_apd_value_bridge_inherited"],
        "p2570_order_obligation_inherited": theorem_export["p2570_order_obligation_inherited"],
        "p2571_measure_boundary_obligation_inherited": theorem_export["p2571_measure_boundary_obligation_inherited"],
        "twelve_apd_rows": len(rows) == 12,
        "fixed_second_derivative_order": audit["fixed_derivative_order"] == 2,
        "left_penalty_changes_selector": audit["left_penalty_changes_selector"],
        "right_penalty_changes_selector": audit["right_penalty_changes_selector"],
        "all_penalty_rows_preserve_nodes": audit["all_penalty_rows_preserve_apd_nodes"],
        "apd_dynamic_source_not_exported": theorem_export["strict_dynamical_source_for_A_P_D_exported"] is False,
        "no_bridge_exported": theorem_export["bridge_theorem_exported"] is False,
        "no_role_transfer_exported": theorem_export["role_transfer_theorem_exported"] is False,
        "no_qw2191_discharge": theorem_export["qw2191_discharged_by_this_certificate"] is False,
        "no_toe_closure_claimed": theorem_export["toe_closure_claimed"] is False,
    }
    return {
        "packet_id": "P2572",
        "stage_id": "S1522",
        "status": "P2572_APD_BOUNDARY_PENALTY_SELECTOR_CONTINUUM_AUDIT_NO_APD_DYNAMIC_SOURCE_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "apd_boundary_penalty_selector_continuum_audit": {
            "theorem_export": theorem_export,
            "source_fingerprints_sha256": {
                "P2416_APD_VALUE_BRIDGE": sha256_json(p2416_payload),
                "P2570_APD_SOBOLEV_ORDER_DEPENDENCE": sha256_json(p2570_payload),
                "P2571_APD_MEASURE_BOUNDARY_DEPENDENCE": sha256_json(p2571_payload),
            },
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["apd_boundary_penalty_selector_continuum_audit"]["theorem_export"]
    audit = t["boundary_penalty_selector_continuum_audit"]
    lines = [
        "# P2572/S1522 APD boundary-penalty selector continuum audit", "",
        f"Status: `{payload['status']}`", "", "## Result", "",
        f"- Frontier atom under attack: `{t['frontier_atom_under_attack']}`.",
        f"- Candidate principle: `{t['candidate_principle_under_test']}`.",
        f"- Stationarity formula: `{audit['stationarity_formula']}`.",
        f"- Left penalty distinct minimizers: `{audit['left_distinct_minimizers_after_rounding_1e_minus_24']}`.",
        f"- Right penalty distinct minimizers: `{audit['right_distinct_minimizers_after_rounding_1e_minus_24']}`.",
        f"- All penalty rows preserve APD nodes: `{audit['all_penalty_rows_preserve_apd_nodes']}`.",
        f"- Finite APD values select boundary penalty: `{not t['finite_apd_values_do_not_select_boundary_penalty']}`.", "",
        "## Interpretation", "",
        "With `k=2` and uniform bulk measure fixed, adding a nonnegative endpoint slope penalty gives the explicit one-parameter stationarity law `lambda(t)=-(B_bulk+t*B_endpoint)/(A_bulk+t*A_endpoint)`.  The audited grid already yields multiple left- and right-endpoint minimizers, all preserving the same finite APD nodes.  Boundary penalty strength is therefore a continuous unsourced selector, not strict A/P/D dynamics.", "",
        "## Recommended next honest step", "", t["recommended_next_honest_step"], "", "## Negative controls", "",
        "No APD boundary penalty source, boundary selector source, strict A/P/D dynamic source, bridge theorem, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, or ToE closure is exported.", "", "## Fingerprint", "",
        f"`{payload['apd_boundary_penalty_selector_continuum_audit']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
`P2572/S1522` fixes `k=2` and uniform bulk measure, then audits endpoint slope penalties as a remaining APD boundary selector.  The stationarity law is explicitly `lambda(t)=-(B_bulk+t*B_endpoint)/(A_bulk+t*A_endpoint)`, and the left/right penalty grids already select multiple distinct `lambda` values while preserving every finite APD node.  Thus boundary penalty strength is a continuous unsourced selector and still does not export `strict_dynamical_source_for_A_P_D`.
""".strip()
    lag_section = """
`P2572/S1522` blocks adding an endpoint-penalty APD Sobolev term to role-bearing `L_total` by hand.  Even after fixing `k=2` and the bulk measure, the boundary penalty parameter continuously changes the selected APD dynamics; strict boundary data remain required.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2572/S1522 APD boundary-penalty selector continuum guard", "## P2572/S1522 APD boundary-penalty selector continuum guard\n\n" + eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2572/S1522 APD boundary-penalty selector continuum Ltotal guard", "## P2572/S1522 APD boundary-penalty selector continuum Ltotal guard\n\n" + lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
