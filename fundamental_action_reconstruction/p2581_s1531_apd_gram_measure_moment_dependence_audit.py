#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from typing import Any

import numpy as np

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import DOC_FILES, REPO, ROOT, load_json, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2569_s1519_apd_value_bridge_interpolation_dynamic_nonuniqueness_certificate import apd_rows
from p2578_s1528_apd_augmented_boundary_basis_metric_dependence_audit import BOUNDARY_TARGETS, DOMAIN, PROBE_POINTS, base_and_vanish, polynomial_combination
from p2580_s1530_apd_inner_product_basis_covariance_requirement_certificate import boundary_matrix, constrained_minimizer

GEN = ROOT / "generated"
OUT = GEN / "p2581_s1531_apd_gram_measure_moment_dependence_audit.json"
MD = GEN / "p2581_s1531_apd_gram_measure_moment_dependence_audit.md"

SOURCE_FILES = {
    "P2416_APD_VALUE_BRIDGE": GEN / "p2416_s1366_apd_multiplicative_bridge_assembly_necessity_certificate.json",
    "P2579_APD_INVERSE_METRIC_TUNABILITY": GEN / "p2579_s1529_apd_inner_product_inverse_metric_tunability_audit.json",
    "P2580_APD_BASIS_COVARIANCE_REQUIREMENT": GEN / "p2580_s1530_apd_inner_product_basis_covariance_requirement_certificate.json",
}
MEASURE_VARIANTS = [
    {"name": "uniform_three_midpoints", "points": [0.25, 5.25, 10.75], "weights": [1.0, 1.0, 1.0]},
    {"name": "left_heavy_three_midpoints", "points": [0.25, 5.25, 10.75], "weights": [9.0, 1.0, 1.0]},
    {"name": "center_heavy_three_midpoints", "points": [0.25, 5.25, 10.75], "weights": [1.0, 9.0, 1.0]},
    {"name": "right_heavy_three_midpoints", "points": [0.25, 5.25, 10.75], "weights": [1.0, 1.0, 9.0]},
    {"name": "staggered_four_point_measure", "points": [0.5, 3.5, 7.5, 10.5], "weights": [2.0, 5.0, 3.0, 7.0]},
]
NEGATIVE_EXPORT_FLAGS = [
    "apd_positive_measure_source_exported",
    "apd_gram_moment_source_exported",
    "apd_l2_inner_product_source_exported",
    "apd_measure_selector_source_exported",
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
        "new_packet": "P2581|S1531|APD Gram measure|Gram measure.*APD",
        "intended_research_nonduplication": "APD.*measure moment|measure moment.*APD|APD.*inner product measure|inner product measure.*APD|APD.*positive measure.*metric|positive measure.*APD.*metric|APD.*L2 Gram|L2 Gram.*APD|APD.*quadrature Gram|quadrature Gram.*APD|APD.*moment metric|moment metric.*APD",
        "apd_precursors": "P2416|S1366|P2579|S1529|P2580|S1530|APD.*inner product|APD.*metric covariance|strict_dynamical_source_for_A_P_D",
        "guardrails": "K_legacy_ont|K_strict_gate|role-transfer theorem|QW-2191|ToE closure|source theorem",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def reference_basis(vanish: Any) -> list[Any]:
    x_poly = type(vanish)([0.0, 1.0])
    return [vanish, x_poly * vanish, x_poly * x_poly * vanish]


def gram_metric_for_measure(basis: list[Any], measure: dict[str, Any]) -> np.ndarray:
    points = np.array(measure["points"], dtype=float)
    weights = np.array(measure["weights"], dtype=float)
    samples = np.array([[poly(point) for poly in basis] for point in points], dtype=float)
    return samples.T @ np.diag(weights) @ samples


def measure_row(measure: dict[str, Any], base: Any, basis: list[Any], matrix: np.ndarray, target: dict[str, Any], rows: list[dict[str, Any]]) -> dict[str, Any]:
    rhs = np.array([target["left_slope_target"] - base.deriv(1)(0.0), target["right_slope_target"] - base.deriv(1)(11.0)], dtype=float)
    gram = gram_metric_for_measure(basis, measure)
    coefficients = constrained_minimizer(matrix, rhs, gram)
    selected = base + polynomial_combination(coefficients, basis)
    node_points = np.array(DOMAIN, dtype=float)
    node_values = np.array([row["apd_value"] for row in rows], dtype=float)
    node_residuals = selected(node_points) - node_values
    left_residual = selected.deriv(1)(0.0) - target["left_slope_target"]
    right_residual = selected.deriv(1)(11.0) - target["right_slope_target"]
    eigenvalues = np.linalg.eigvalsh(gram)
    return {
        "measure_name": measure["name"],
        "points": [float(value) for value in measure["points"]],
        "weights": [float(value) for value in measure["weights"]],
        "gram_metric_3x3": [[float(value) for value in row] for row in gram],
        "gram_metric_eigenvalues": [float(value) for value in eigenvalues],
        "gram_metric_is_positive_definite": bool(np.min(eigenvalues) > 0.0),
        "selected_coefficients": [float(value) for value in coefficients],
        "probe_values": [float(selected(point)) for point in PROBE_POINTS],
        "max_abs_node_residual": float(np.max(np.abs(node_residuals))),
        "max_abs_boundary_slope_residual": float(max(abs(left_residual), abs(right_residual))),
    }


def target_row(target: dict[str, Any], base: Any, basis: list[Any], matrix: np.ndarray, rows: list[dict[str, Any]]) -> dict[str, Any]:
    measure_rows = [measure_row(measure, base, basis, matrix, target, rows) for measure in MEASURE_VARIANTS]
    middle_probe_values = {round(row["probe_values"][1], 12) for row in measure_rows}
    return {
        "target_name": target["name"],
        "measure_rows": measure_rows,
        "measure_count": len(measure_rows),
        "distinct_middle_probe_values_after_rounding_1e_minus_12": len(middle_probe_values),
        "all_gram_metrics_positive_definite": all(row["gram_metric_is_positive_definite"] for row in measure_rows),
        "all_measure_minimizers_preserve_apd_nodes": all(row["max_abs_node_residual"] <= 1.0e-6 for row in measure_rows),
        "all_measure_minimizers_satisfy_boundary_targets": all(row["max_abs_boundary_slope_residual"] <= 1.0e-6 for row in measure_rows),
        "positive_measure_gram_metric_is_measure_dependent": len(middle_probe_values) > 1,
    }


def gram_measure_audit(rows: list[dict[str, Any]]) -> dict[str, Any]:
    base, vanish = base_and_vanish(rows)
    basis = reference_basis(vanish)
    matrix = boundary_matrix(basis)
    target_rows = [target_row(target, base, basis, matrix, rows) for target in BOUNDARY_TARGETS]
    return {
        "numpy_version": np.__version__,
        "domain": DOMAIN,
        "probe_points": PROBE_POINTS,
        "coefficient_basis": ["V", "xV", "x^2V"],
        "measure_variants": MEASURE_VARIANTS,
        "target_rows": target_rows,
        "target_count": len(target_rows),
        "measure_variant_count": len(MEASURE_VARIANTS),
        "all_targets_measure_dependent": all(row["positive_measure_gram_metric_is_measure_dependent"] for row in target_rows),
        "all_gram_metrics_positive_definite": all(row["all_gram_metrics_positive_definite"] for row in target_rows),
        "all_measure_minimizers_preserve_nodes_and_boundaries": all(row["all_measure_minimizers_preserve_apd_nodes"] and row["all_measure_minimizers_satisfy_boundary_targets"] for row in target_rows),
        "gram_measure_is_unsourced_inner_product_data": True,
    }


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    p2416_payload = load_json(SOURCE_FILES["P2416_APD_VALUE_BRIDGE"])
    p2579_payload = load_json(SOURCE_FILES["P2579_APD_INVERSE_METRIC_TUNABILITY"])
    p2580_payload = load_json(SOURCE_FILES["P2580_APD_BASIS_COVARIANCE_REQUIREMENT"])
    p2579 = theorem(p2579_payload, "apd_inner_product_inverse_metric_tunability_audit")
    p2580 = theorem(p2580_payload, "apd_inner_product_basis_covariance_requirement_certificate")
    rows = apd_rows(p2416_payload)
    audit = gram_measure_audit(rows)
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2581_T1_apd_gram_measure_moment_dependence_audit",
        "audited_chain": ["P2416/S1366", "P2579/S1529", "P2580/S1530"],
        "frontier_atom_under_attack": "strict_dynamical_source_for_A_P_D",
        "candidate_principle_under_test": "derive APD inner product from a positive L2 Gram measure on the augmented vanishing-mode coefficient space",
        "p2416_apd_value_bridge_inherited": p2416_payload.get("apd_multiplicative_bridge_assembly_necessity_certificate", {}).get("finite_witness_certificate", {}).get("apd_value_assembly_witness_ready") is True,
        "p2579_inner_product_tunability_inherited": p2579.get("finite_apd_values_and_boundary_targets_do_not_select_function_space_inner_product") is True,
        "p2580_metric_covariance_requirement_inherited": p2580.get("coordinate_covariance_removes_basis_artifact_but_not_metric_choice") is True,
        "apd_node_rows": rows,
        "apd_gram_measure_moment_dependence_audit": audit,
        "finite_apd_values_and_boundary_targets_do_not_select_positive_measure": audit["all_targets_measure_dependent"],
        "positive_measure_gram_selector_is_conditional_not_strict_source": True,
        "apd_dynamic_source_remains_unsourced": True,
        "recommended_next_honest_step": (
            "Do not treat an L2/Gram APD inner product as strict unless the positive measure and its moments are derived from strict nadsoliton dynamics. P2581 shows that multiple positive measures give positive-definite, covariance-compatible Gram metrics while selecting different APD off-node dynamics under the same finite values and boundary targets. The next honest step is to derive the APD measure density or kinetic weight, not merely choose a quadrature rule."
        ),
        **{flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    gatekeepers = {
        "rg_non_duplication_audit_recorded": grep["tool"] == "rg",
        "p2416_apd_value_bridge_inherited": theorem_export["p2416_apd_value_bridge_inherited"],
        "p2579_inner_product_tunability_inherited": theorem_export["p2579_inner_product_tunability_inherited"],
        "p2580_metric_covariance_requirement_inherited": theorem_export["p2580_metric_covariance_requirement_inherited"],
        "three_targets_audited": audit["target_count"] == 3,
        "five_positive_measures": audit["measure_variant_count"] == 5,
        "all_targets_measure_dependent": audit["all_targets_measure_dependent"],
        "all_gram_metrics_positive_definite": audit["all_gram_metrics_positive_definite"],
        "all_measure_minimizers_preserve_nodes_and_boundaries": audit["all_measure_minimizers_preserve_nodes_and_boundaries"],
        "apd_dynamic_source_not_exported": theorem_export["strict_dynamical_source_for_A_P_D_exported"] is False,
        "no_bridge_exported": theorem_export["bridge_theorem_exported"] is False,
        "no_role_transfer_exported": theorem_export["role_transfer_theorem_exported"] is False,
        "no_qw2191_discharge": theorem_export["qw2191_discharged_by_this_certificate"] is False,
        "no_toe_closure_claimed": theorem_export["toe_closure_claimed"] is False,
    }
    return {
        "packet_id": "P2581",
        "stage_id": "S1531",
        "status": "P2581_APD_GRAM_MEASURE_MOMENT_DEPENDENCE_AUDIT_NO_APD_DYNAMIC_SOURCE_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "apd_gram_measure_moment_dependence_audit": {
            "theorem_export": theorem_export,
            "source_fingerprints_sha256": {
                "P2416_APD_VALUE_BRIDGE": sha256_json(p2416_payload),
                "P2579_APD_INVERSE_METRIC_TUNABILITY": sha256_json(p2579_payload),
                "P2580_APD_BASIS_COVARIANCE_REQUIREMENT": sha256_json(p2580_payload),
            },
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["apd_gram_measure_moment_dependence_audit"]["theorem_export"]
    audit = t["apd_gram_measure_moment_dependence_audit"]
    lines = [
        "# P2581/S1531 APD Gram-measure moment dependence audit", "",
        f"Status: `{payload['status']}`", "", "## Result", "",
        f"- Frontier atom under attack: `{t['frontier_atom_under_attack']}`.",
        f"- Boundary targets audited: `{audit['target_count']}`.",
        f"- Positive measures audited: `{audit['measure_variant_count']}`.",
        f"- All targets measure dependent: `{audit['all_targets_measure_dependent']}`.",
        f"- All Gram metrics positive definite: `{audit['all_gram_metrics_positive_definite']}`.",
        f"- Strict APD dynamic source exported: `{t['strict_dynamical_source_for_A_P_D_exported']}`.", "",
        "## Interpretation", "",
        "Positive L2/Gram metrics are covariance-compatible candidates, but the positive measure is extra data.  P2581 audits five positive measures and finds that they preserve the same finite APD nodes and boundary slopes while selecting different off-node APD dynamics.  Thus the missing source is pushed from coordinates to the measure/moment law.", "",
        "## Recommended next honest step", "", t["recommended_next_honest_step"], "", "## Negative controls", "",
        "No APD positive-measure source, Gram-moment source, L2-inner-product source, measure-selector source, strict A/P/D dynamic source, bridge theorem, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, or ToE closure is exported.", "", "## Fingerprint", "",
        f"`{payload['apd_gram_measure_moment_dependence_audit']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2581/S1531 APD Gram-measure moment dependence guard

`P2581/S1531` tests the next source candidate after P2580: choose an APD inner product as a positive L2/Gram metric on the augmented vanishing-mode basis.  Across five positive measures, all Gram metrics are positive definite and preserve the same finite APD nodes plus endpoint slopes, but the selected off-node APD dynamics changes with the measure.  Thus covariance-compatible Gram metrics still require a strict source for the measure/moment law before they can export `strict_dynamical_source_for_A_P_D`.
""".strip()
    lag_section = """
## P2581/S1531 APD Gram-measure moment dependence Ltotal guard

`P2581/S1531` blocks inserting a role-bearing APD L2/Gram kinetic term into `L_total` by choosing a convenient positive measure.  The measure and its moments must be derived from strict nadsoliton dynamics; otherwise the same APD value/boundary constraints generate measure-dependent off-node dynamics.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2581/S1531 APD Gram-measure moment dependence guard", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2581/S1531 APD Gram-measure moment dependence Ltotal guard", lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
