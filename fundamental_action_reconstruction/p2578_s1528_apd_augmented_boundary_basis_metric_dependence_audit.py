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
OUT = GEN / "p2578_s1528_apd_augmented_boundary_basis_metric_dependence_audit.json"
MD = GEN / "p2578_s1528_apd_augmented_boundary_basis_metric_dependence_audit.md"

SOURCE_FILES = {
    "P2416_APD_VALUE_BRIDGE": GEN / "p2416_s1366_apd_multiplicative_bridge_assembly_necessity_certificate.json",
    "P2575_APD_AUGMENTED_BOUNDARY_NULLSPACE": GEN / "p2575_s1525_apd_augmented_boundary_nullspace_nonuniqueness_audit.json",
    "P2577_APD_GRID_MEASURE_DEPENDENCE": GEN / "p2577_s1527_apd_boundary_nullspace_grid_measure_dependence_audit.json",
}
DOMAIN = list(range(12))
PROBE_POINTS = [0.5, 5.5, 11.5]
BOUNDARY_TARGETS = [
    {"name": "zero_zero_neumann", "left_slope_target": 0.0, "right_slope_target": 0.0},
    {"name": "equal_positive_slopes", "left_slope_target": 1.0e-6, "right_slope_target": 1.0e-6},
    {"name": "equal_negative_slopes", "left_slope_target": -1.0e-6, "right_slope_target": -1.0e-6},
]
NEGATIVE_EXPORT_FLAGS = [
    "apd_vanishing_basis_source_exported",
    "apd_coordinate_metric_source_exported",
    "apd_minimum_coefficient_norm_source_exported",
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
        "new_packet": "P2578|S1528|APD basis metric|basis metric.*APD",
        "intended_research_nonduplication": "APD.*coordinate norm|coordinate norm.*APD|APD.*minimum coefficient|minimum coefficient.*APD|APD.*basis dependence|basis dependence.*APD|vanishing-mode basis|basis.*nullspace.*APD|APD.*coordinate metric",
        "apd_precursors": "P2416|S1366|P2575|S1525|P2577|S1527|APD.*nullspace|strict_dynamical_source_for_A_P_D",
        "guardrails": "K_legacy_ont|K_strict_gate|role-transfer theorem|QW-2191|ToE closure|source theorem",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def base_and_vanish(rows: list[dict[str, Any]]) -> tuple[Polynomial, Polynomial]:
    x_nodes = np.array([row["d"] for row in rows], dtype=float)
    y_nodes = np.array([row["apd_value"] for row in rows], dtype=float)
    base = Polynomial.fit(x_nodes, y_nodes, deg=11, domain=[0.0, 11.0]).convert()
    vanish = Polynomial([1.0])
    for d in DOMAIN:
        vanish = vanish * Polynomial([-float(d), 1.0])
    return base, vanish


def basis_variants(vanish: Polynomial) -> list[dict[str, Any]]:
    x_poly = Polynomial([0.0, 1.0])
    return [
        {"name": "monomial_V_xV_x2V", "basis": [vanish, x_poly * vanish, x_poly * x_poly * vanish]},
        {"name": "centered_V_xminus5p5V_xminus5p5sqV", "basis": [vanish, (x_poly - 5.5) * vanish, (x_poly - 5.5) * (x_poly - 5.5) * vanish]},
        {"name": "scaled_V_xover11V_xover11sqV", "basis": [vanish, (x_poly / 11.0) * vanish, (x_poly / 11.0) * (x_poly / 11.0) * vanish]},
        {"name": "shifted_mixed_V_xminus3V_xminus3_xminus7V", "basis": [vanish, (x_poly - 3.0) * vanish, (x_poly - 3.0) * (x_poly - 7.0) * vanish]},
    ]


def polynomial_combination(coefficients: np.ndarray, basis: list[Polynomial]) -> Polynomial:
    total = Polynomial([0.0])
    for coefficient, poly in zip(coefficients, basis):
        total = total + float(coefficient) * poly
    return total


def variant_row(variant: dict[str, Any], base: Polynomial, rows: list[dict[str, Any]], target: dict[str, Any]) -> dict[str, Any]:
    basis = variant["basis"]
    matrix = np.array([[poly.deriv(1)(0.0) for poly in basis], [poly.deriv(1)(11.0) for poly in basis]], dtype=float)
    rhs = np.array([target["left_slope_target"] - base.deriv(1)(0.0), target["right_slope_target"] - base.deriv(1)(11.0)], dtype=float)
    coeffs = np.linalg.lstsq(matrix, rhs, rcond=None)[0]
    correction = polynomial_combination(coeffs, basis)
    selected = base + correction
    node_points = np.array(DOMAIN, dtype=float)
    node_values = np.array([row["apd_value"] for row in rows], dtype=float)
    node_residuals = selected(node_points) - node_values
    left_residual = selected.deriv(1)(0.0) - target["left_slope_target"]
    right_residual = selected.deriv(1)(11.0) - target["right_slope_target"]
    probe_values = [float(selected(point)) for point in PROBE_POINTS]
    return {
        "basis_name": variant["name"],
        "boundary_matrix_rank": int(np.linalg.matrix_rank(matrix)),
        "minimum_coordinate_norm_coefficients": [float(value) for value in coeffs],
        "coefficient_l2_norm": float(np.linalg.norm(coeffs)),
        "probe_values": probe_values,
        "max_abs_node_residual": float(np.max(np.abs(node_residuals))),
        "max_abs_boundary_slope_residual": float(max(abs(left_residual), abs(right_residual))),
    }


def target_row(target: dict[str, Any], base: Polynomial, vanish: Polynomial, rows: list[dict[str, Any]]) -> dict[str, Any]:
    variants = [variant_row(variant, base, rows, target) for variant in basis_variants(vanish)]
    probe_at_middle = {round(row["probe_values"][1], 18) for row in variants}
    return {
        "target_name": target["name"],
        "variant_rows": variants,
        "variant_count": len(variants),
        "distinct_middle_probe_values_after_rounding_1e_minus_18": len(probe_at_middle),
        "all_basis_variants_preserve_apd_nodes": all(row["max_abs_node_residual"] <= 1.0e-6 for row in variants),
        "all_basis_variants_satisfy_boundary_targets": all(row["max_abs_boundary_slope_residual"] <= 1.0e-6 for row in variants),
        "minimum_coordinate_norm_depends_on_basis": len(probe_at_middle) > 1,
    }


def basis_metric_audit(rows: list[dict[str, Any]]) -> dict[str, Any]:
    base, vanish = base_and_vanish(rows)
    target_rows = [target_row(target, base, vanish, rows) for target in BOUNDARY_TARGETS]
    return {
        "numpy_version": np.__version__,
        "domain": DOMAIN,
        "probe_points": PROBE_POINTS,
        "basis_variants": [variant["name"] for variant in basis_variants(vanish)],
        "target_rows": target_rows,
        "target_count": len(target_rows),
        "basis_variant_count": len(basis_variants(vanish)),
        "all_targets_basis_metric_dependent": all(row["minimum_coordinate_norm_depends_on_basis"] for row in target_rows),
        "all_basis_metric_solutions_preserve_nodes_and_boundaries": all(row["all_basis_variants_preserve_apd_nodes"] and row["all_basis_variants_satisfy_boundary_targets"] for row in target_rows),
        "coordinate_metric_is_unsourced_selector": True,
    }


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    p2416_payload = load_json(SOURCE_FILES["P2416_APD_VALUE_BRIDGE"])
    p2575_payload = load_json(SOURCE_FILES["P2575_APD_AUGMENTED_BOUNDARY_NULLSPACE"])
    p2577_payload = load_json(SOURCE_FILES["P2577_APD_GRID_MEASURE_DEPENDENCE"])
    p2575 = theorem(p2575_payload, "apd_augmented_boundary_nullspace_nonuniqueness_audit")
    p2577 = theorem(p2577_payload, "apd_boundary_nullspace_grid_measure_dependence_audit")
    rows = apd_rows(p2416_payload)
    audit = basis_metric_audit(rows)
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2578_T1_apd_augmented_boundary_basis_metric_dependence_audit",
        "audited_chain": ["P2416/S1366", "P2575/S1525", "P2577/S1527"],
        "frontier_atom_under_attack": "strict_dynamical_source_for_A_P_D",
        "candidate_principle_under_test": "select the augmented APD boundary solution by minimum Euclidean coefficient norm in a chosen vanishing-mode basis",
        "p2416_apd_value_bridge_inherited": p2416_payload.get("apd_multiplicative_bridge_assembly_necessity_certificate", {}).get("finite_witness_certificate", {}).get("apd_value_assembly_witness_ready") is True,
        "p2575_augmented_nullspace_inherited": p2575.get("but_augmented_boundary_solution_is_nonunique") is True,
        "p2577_grid_measure_dependence_inherited": p2577.get("finite_apd_values_and_boundary_targets_do_not_select_grid_measure") is True,
        "apd_node_rows": rows,
        "augmented_boundary_basis_metric_dependence_audit": audit,
        "finite_apd_values_and_boundary_targets_do_not_select_coordinate_metric": audit["all_targets_basis_metric_dependent"],
        "minimum_coefficient_norm_selector_is_conditional_not_strict_source": True,
        "apd_dynamic_source_remains_unsourced": True,
        "recommended_next_honest_step": (
            "Do not select the APD augmented-boundary solution by a minimum coefficient norm until the vanishing-mode basis and coordinate metric are strict-sourced. The next honest step is to derive the APD function-space inner product from the strict action; otherwise the same nodes and boundary targets produce basis-dependent APD dynamics."
        ),
        **{flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    gatekeepers = {
        "rg_non_duplication_audit_recorded": grep["tool"] == "rg",
        "p2416_apd_value_bridge_inherited": theorem_export["p2416_apd_value_bridge_inherited"],
        "p2575_augmented_nullspace_inherited": theorem_export["p2575_augmented_nullspace_inherited"],
        "p2577_grid_measure_dependence_inherited": theorem_export["p2577_grid_measure_dependence_inherited"],
        "three_targets_audited": audit["target_count"] == 3,
        "four_basis_variants": audit["basis_variant_count"] == 4,
        "all_targets_basis_metric_dependent": audit["all_targets_basis_metric_dependent"],
        "all_basis_metric_solutions_preserve_nodes_and_boundaries": audit["all_basis_metric_solutions_preserve_nodes_and_boundaries"],
        "apd_dynamic_source_not_exported": theorem_export["strict_dynamical_source_for_A_P_D_exported"] is False,
        "no_bridge_exported": theorem_export["bridge_theorem_exported"] is False,
        "no_role_transfer_exported": theorem_export["role_transfer_theorem_exported"] is False,
        "no_qw2191_discharge": theorem_export["qw2191_discharged_by_this_certificate"] is False,
        "no_toe_closure_claimed": theorem_export["toe_closure_claimed"] is False,
    }
    return {
        "packet_id": "P2578",
        "stage_id": "S1528",
        "status": "P2578_APD_AUGMENTED_BOUNDARY_BASIS_METRIC_DEPENDENCE_AUDIT_NO_APD_DYNAMIC_SOURCE_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "apd_augmented_boundary_basis_metric_dependence_audit": {
            "theorem_export": theorem_export,
            "source_fingerprints_sha256": {
                "P2416_APD_VALUE_BRIDGE": sha256_json(p2416_payload),
                "P2575_APD_AUGMENTED_BOUNDARY_NULLSPACE": sha256_json(p2575_payload),
                "P2577_APD_GRID_MEASURE_DEPENDENCE": sha256_json(p2577_payload),
            },
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["apd_augmented_boundary_basis_metric_dependence_audit"]["theorem_export"]
    audit = t["augmented_boundary_basis_metric_dependence_audit"]
    lines = [
        "# P2578/S1528 APD augmented-boundary basis-metric dependence audit", "",
        f"Status: `{payload['status']}`", "", "## Result", "",
        f"- Frontier atom under attack: `{t['frontier_atom_under_attack']}`.",
        f"- Basis variants audited: `{audit['basis_variant_count']}`.",
        f"- Targets audited: `{audit['target_count']}`.",
        f"- All targets basis/metric dependent: `{audit['all_targets_basis_metric_dependent']}`.",
        f"- Solutions preserve nodes and boundaries: `{audit['all_basis_metric_solutions_preserve_nodes_and_boundaries']}`.",
        f"- Strict APD dynamic source exported: `{t['strict_dynamical_source_for_A_P_D_exported']}`.", "",
        "## Interpretation", "",
        "Minimum Euclidean coefficient norm is not coordinate-free on the augmented vanishing-mode space.  Changing from monomial to centered, scaled, or shifted bases preserves the finite APD nodes and endpoint slopes but changes the selected off-node APD values.  Therefore a coefficient-norm selector needs a sourced basis/inner product.", "",
        "## Recommended next honest step", "", t["recommended_next_honest_step"], "", "## Negative controls", "",
        "No APD vanishing-basis source, coordinate-metric source, minimum-coefficient-norm source, strict A/P/D dynamic source, bridge theorem, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, or ToE closure is exported.", "", "## Fingerprint", "",
        f"`{payload['apd_augmented_boundary_basis_metric_dependence_audit']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
`P2578/S1528` audits a different APD nullspace selector: choose the augmented-boundary solution with minimum Euclidean coefficient norm in a vanishing-mode basis.  Changing the basis among monomial, centered, scaled, and shifted variants preserves finite APD nodes and endpoint slopes but changes the selected off-node APD values.  Thus coefficient-norm minimization is basis/metric dependent and still does not export `strict_dynamical_source_for_A_P_D`.
""".strip()
    lag_section = """
`P2578/S1528` blocks adding a minimum-coefficient-norm APD selector to role-bearing `L_total` without sourcing the vanishing-mode basis and coordinate metric.  A strict APD action must provide the function-space inner product; otherwise the same constraints yield basis-dependent dynamics.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2578/S1528 APD augmented-boundary basis-metric dependence guard", "## P2578/S1528 APD augmented-boundary basis-metric dependence guard\n\n" + eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2578/S1528 APD augmented-boundary basis-metric dependence Ltotal guard", "## P2578/S1528 APD augmented-boundary basis-metric dependence Ltotal guard\n\n" + lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
