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
OUT = GEN / "p2571_s1521_apd_sobolev_measure_boundary_dependence_audit.json"
MD = GEN / "p2571_s1521_apd_sobolev_measure_boundary_dependence_audit.md"

SOURCE_FILES = {
    "P2416_APD_VALUE_BRIDGE": GEN / "p2416_s1366_apd_multiplicative_bridge_assembly_necessity_certificate.json",
    "P2569_APD_INTERPOLATION_NONUNIQUENESS": GEN / "p2569_s1519_apd_value_bridge_interpolation_dynamic_nonuniqueness_certificate.json",
    "P2570_APD_SOBOLEV_ORDER_DEPENDENCE": GEN / "p2570_s1520_apd_sobolev_roughness_selector_order_dependence_audit.json",
}
DOMAIN = list(range(12))
FIXED_DERIVATIVE_ORDER = 2
PROBE_POINTS = [0.5, 5.5, 11.5]
VARIANTS = [
    {"name": "uniform_free", "weight_coefficients": [1.0], "slope_penalty_left": 0.0, "slope_penalty_right": 0.0},
    {"name": "left_tilt_free", "weight_coefficients": [2.0, -1.0 / 11.0], "slope_penalty_left": 0.0, "slope_penalty_right": 0.0},
    {"name": "right_tilt_free", "weight_coefficients": [1.0, 1.0 / 11.0], "slope_penalty_left": 0.0, "slope_penalty_right": 0.0},
    {"name": "center_quadratic_free", "weight_coefficients": [2.0, -2.0 / 5.5, 1.0 / (5.5 * 5.5)], "slope_penalty_left": 0.0, "slope_penalty_right": 0.0},
    {"name": "uniform_both_endpoint_slope_penalty", "weight_coefficients": [1.0], "slope_penalty_left": 1.0, "slope_penalty_right": 1.0},
    {"name": "uniform_left_endpoint_slope_penalty", "weight_coefficients": [1.0], "slope_penalty_left": 1.0, "slope_penalty_right": 0.0},
    {"name": "uniform_right_endpoint_slope_penalty", "weight_coefficients": [1.0], "slope_penalty_left": 0.0, "slope_penalty_right": 1.0},
]
NEGATIVE_EXPORT_FLAGS = [
    "apd_measure_source_exported",
    "apd_boundary_class_source_exported",
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
        "new_packet": "P2571|S1521|APD Sobolev measure|APD boundary selector|APD measure boundary",
        "intended_research_nonduplication": "APD.*measure.*boundary|APD.*Sobolev.*measure|measure.*APD.*roughness|boundary.*APD.*roughness|APD.*boundary class|Sobolev.*boundary.*APD|weighted roughness.*APD|APD.*weighted roughness",
        "apd_precursors": "P2416|S1366|P2569|S1519|P2570|S1520|APD interpolation|APD Sobolev|strict_dynamical_source_for_A_P_D",
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


def variant_row(variant: dict[str, Any], base: Polynomial, vanish: Polynomial, rows: list[dict[str, Any]]) -> dict[str, Any]:
    weight = Polynomial(variant["weight_coefficients"])
    d_base = base.deriv(FIXED_DERIVATIVE_ORDER)
    d_vanish = vanish.deriv(FIXED_DERIVATIVE_ORDER)
    a_quad = definite_integral(weight * d_vanish * d_vanish)
    b_cross = definite_integral(weight * d_base * d_vanish)
    c_base = definite_integral(weight * d_base * d_base)
    slope_base = base.deriv(1)
    slope_vanish = vanish.deriv(1)
    boundary_terms = []
    for point, penalty, label in [(0.0, variant["slope_penalty_left"], "left"), (11.0, variant["slope_penalty_right"], "right")]:
        b_value = float(slope_base(point))
        v_value = float(slope_vanish(point))
        a_quad += penalty * v_value * v_value
        b_cross += penalty * b_value * v_value
        c_base += penalty * b_value * b_value
        boundary_terms.append({"endpoint": label, "point": point, "slope_penalty": penalty, "base_slope": b_value, "vanishing_slope": v_value})
    lambda_star = -b_cross / a_quad if abs(a_quad) > 0.0 else 0.0
    selected = base + lambda_star * vanish
    node_points = np.array(DOMAIN, dtype=float)
    node_values = np.array([row["apd_value"] for row in rows], dtype=float)
    residuals = selected(node_points) - node_values
    selected_cost = c_base + 2.0 * b_cross * lambda_star + a_quad * lambda_star * lambda_star
    probe_deltas = [float(lambda_star * vanish(point)) for point in PROBE_POINTS]
    weight_samples = [float(weight(point)) for point in [0.0, 5.5, 11.0]]
    return {
        "variant": variant["name"],
        "fixed_derivative_order": FIXED_DERIVATIVE_ORDER,
        "weight_coefficients_in_power_basis": variant["weight_coefficients"],
        "weight_samples_at_0_5p5_11": weight_samples,
        "boundary_terms": boundary_terms,
        "quadratic_A": float(a_quad),
        "linear_cross_B": float(b_cross),
        "base_cost_C": float(c_base),
        "lambda_star": float(lambda_star),
        "selected_cost": float(selected_cost),
        "base_cost": float(c_base),
        "cost_improvement_over_base": float(c_base - selected_cost),
        "max_abs_node_residual_at_selected_lambda": float(np.max(np.abs(residuals))),
        "probe_value_deltas_from_base": probe_deltas,
        "selected_member_changes_off_node_values": any(abs(value) > 0.0 for value in probe_deltas),
    }


def measure_boundary_audit(rows: list[dict[str, Any]]) -> dict[str, Any]:
    base, vanish = interpolation_polynomials(rows)
    audit_rows = [variant_row(variant, base, vanish, rows) for variant in VARIANTS]
    rounded_lambdas = {round(row["lambda_star"], 20) for row in audit_rows}
    return {
        "numpy_version": np.__version__,
        "fixed_derivative_order": FIXED_DERIVATIVE_ORDER,
        "domain": DOMAIN,
        "probe_points": PROBE_POINTS,
        "variant_count": len(audit_rows),
        "variant_rows": audit_rows,
        "distinct_minimizers_after_rounding_1e_minus_20": len(rounded_lambdas),
        "all_variants_preserve_apd_nodes": all(row["max_abs_node_residual_at_selected_lambda"] <= 1.0e-6 for row in audit_rows),
        "all_variants_change_off_node_values": all(row["selected_member_changes_off_node_values"] for row in audit_rows),
        "measure_or_boundary_changes_selector": len(rounded_lambdas) > 1,
        "measure_and_boundary_not_sourced_by_finite_apd_values": True,
    }


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    p2416_payload = load_json(SOURCE_FILES["P2416_APD_VALUE_BRIDGE"])
    p2569_payload = load_json(SOURCE_FILES["P2569_APD_INTERPOLATION_NONUNIQUENESS"])
    p2570_payload = load_json(SOURCE_FILES["P2570_APD_SOBOLEV_ORDER_DEPENDENCE"])
    p2569 = theorem(p2569_payload, "apd_value_bridge_interpolation_dynamic_nonuniqueness_certificate")
    p2570 = theorem(p2570_payload, "apd_sobolev_roughness_selector_order_dependence_audit")
    rows = apd_rows(p2416_payload)
    audit = measure_boundary_audit(rows)
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2571_T1_apd_sobolev_measure_boundary_dependence_audit",
        "audited_chain": ["P2416/S1366", "P2569/S1519", "P2570/S1520"],
        "frontier_atom_under_attack": "strict_dynamical_source_for_A_P_D",
        "candidate_principle_under_test": "fix Sobolev derivative order k=2, then vary positive measure weights and endpoint slope penalties on the P2569 interpolation family",
        "p2416_apd_value_bridge_inherited": p2416_payload.get("apd_multiplicative_bridge_assembly_necessity_certificate", {}).get("finite_witness_certificate", {}).get("apd_value_assembly_witness_ready") is True,
        "p2569_interpolation_family_inherited": p2569.get("finite_apd_values_do_not_select_dynamic_law") is True,
        "p2570_order_obligation_inherited": p2570.get("finite_apd_values_do_not_select_sobolev_order") is True,
        "apd_node_rows": rows,
        "sobolev_measure_boundary_dependence_audit": audit,
        "finite_apd_values_do_not_select_measure_or_boundary_class": audit["measure_or_boundary_changes_selector"],
        "fixed_order_two_still_conditional": True,
        "apd_dynamic_source_remains_unsourced": True,
        "recommended_next_honest_step": (
            "Do not treat k=2 hydrodynamic/Sobolev intuition as an APD source theorem. The next honest step is to derive the measure and boundary class of the APD action from strict nadsoliton dynamics; without that derivation, even a fixed second-order roughness law is a family of conditional selectors rather than strict A/P/D dynamics."
        ),
        **{flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    gatekeepers = {
        "rg_non_duplication_audit_recorded": grep["tool"] == "rg",
        "p2416_apd_value_bridge_inherited": theorem_export["p2416_apd_value_bridge_inherited"],
        "p2569_interpolation_family_inherited": theorem_export["p2569_interpolation_family_inherited"],
        "p2570_order_obligation_inherited": theorem_export["p2570_order_obligation_inherited"],
        "twelve_apd_rows": len(rows) == 12,
        "fixed_second_derivative_order": audit["fixed_derivative_order"] == 2,
        "seven_variants_audited": audit["variant_count"] == 7,
        "all_variants_preserve_nodes": audit["all_variants_preserve_apd_nodes"],
        "measure_or_boundary_changes_selector": audit["measure_or_boundary_changes_selector"],
        "apd_dynamic_source_not_exported": theorem_export["strict_dynamical_source_for_A_P_D_exported"] is False,
        "no_bridge_exported": theorem_export["bridge_theorem_exported"] is False,
        "no_role_transfer_exported": theorem_export["role_transfer_theorem_exported"] is False,
        "no_qw2191_discharge": theorem_export["qw2191_discharged_by_this_certificate"] is False,
        "no_toe_closure_claimed": theorem_export["toe_closure_claimed"] is False,
    }
    return {
        "packet_id": "P2571",
        "stage_id": "S1521",
        "status": "P2571_APD_SOBOLEV_MEASURE_BOUNDARY_DEPENDENCE_AUDIT_NO_APD_DYNAMIC_SOURCE_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "apd_sobolev_measure_boundary_dependence_audit": {
            "theorem_export": theorem_export,
            "source_fingerprints_sha256": {
                "P2416_APD_VALUE_BRIDGE": sha256_json(p2416_payload),
                "P2569_APD_INTERPOLATION_NONUNIQUENESS": sha256_json(p2569_payload),
                "P2570_APD_SOBOLEV_ORDER_DEPENDENCE": sha256_json(p2570_payload),
            },
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["apd_sobolev_measure_boundary_dependence_audit"]["theorem_export"]
    audit = t["sobolev_measure_boundary_dependence_audit"]
    lines = [
        "# P2571/S1521 APD Sobolev measure/boundary dependence audit", "",
        f"Status: `{payload['status']}`", "", "## Result", "",
        f"- Frontier atom under attack: `{t['frontier_atom_under_attack']}`.",
        f"- Candidate principle: `{t['candidate_principle_under_test']}`.",
        f"- Fixed derivative order: `{audit['fixed_derivative_order']}`.",
        f"- Variants audited: `{audit['variant_count']}`.",
        f"- Distinct minimizers after rounding: `{audit['distinct_minimizers_after_rounding_1e_minus_20']}`.",
        f"- All variants preserve APD nodes: `{audit['all_variants_preserve_apd_nodes']}`.",
        f"- Measure or boundary changes selector: `{audit['measure_or_boundary_changes_selector']}`.",
        f"- Finite APD values select measure/boundary class: `{not t['finite_apd_values_do_not_select_measure_or_boundary_class']}`.", "",
        "## Interpretation", "",
        "Even after fixing the Sobolev derivative order to `k=2`, changing the positive measure weight or endpoint slope-penalty class changes the minimizing `lambda` on the same APD interpolation family.  The finite APD nodes remain exact in every audited variant, so measure and boundary data are additional strict source obligations.", "",
        "## Recommended next honest step", "", t["recommended_next_honest_step"], "", "## Negative controls", "",
        "No APD measure source, boundary-class source, weighted roughness selector source, strict A/P/D dynamic source, bridge theorem, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, or ToE closure is exported.", "", "## Fingerprint", "",
        f"`{payload['apd_sobolev_measure_boundary_dependence_audit']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
`P2571/S1521` fixes the APD Sobolev derivative order at `k=2` and audits the next hidden choices: positive measure weights and endpoint slope-penalty boundary classes.  The minimizing `lambda` changes across `7` audited variants while every selected member preserves the finite APD nodes.  Thus even second-order roughness does not export `strict_dynamical_source_for_A_P_D` unless the measure and boundary class are strict-sourced.
""".strip()
    lag_section = """
`P2571/S1521` blocks importing a hydrodynamic `k=2` APD roughness term into role-bearing `L_total` without its measure and boundary data.  The same APD values support distinct weighted/boundary Sobolev minimizers, so a sourced APD action remains required.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2571/S1521 APD Sobolev measure/boundary dependence guard", "## P2571/S1521 APD Sobolev measure/boundary dependence guard\n\n" + eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2571/S1521 APD Sobolev measure/boundary dependence Ltotal guard", "## P2571/S1521 APD Sobolev measure/boundary dependence Ltotal guard\n\n" + lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
