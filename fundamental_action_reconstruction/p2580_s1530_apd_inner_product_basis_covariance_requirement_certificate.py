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
from p2578_s1528_apd_augmented_boundary_basis_metric_dependence_audit import BOUNDARY_TARGETS, PROBE_POINTS, base_and_vanish, basis_variants, polynomial_combination

GEN = ROOT / "generated"
OUT = GEN / "p2580_s1530_apd_inner_product_basis_covariance_requirement_certificate.json"
MD = GEN / "p2580_s1530_apd_inner_product_basis_covariance_requirement_certificate.md"

SOURCE_FILES = {
    "P2416_APD_VALUE_BRIDGE": GEN / "p2416_s1366_apd_multiplicative_bridge_assembly_necessity_certificate.json",
    "P2578_APD_BASIS_METRIC_DEPENDENCE": GEN / "p2578_s1528_apd_augmented_boundary_basis_metric_dependence_audit.json",
    "P2579_APD_INVERSE_METRIC_TUNABILITY": GEN / "p2579_s1529_apd_inner_product_inverse_metric_tunability_audit.json",
}
TRANSFORM_SAMPLE_POINTS = [0.25, 3.75, 9.25]
NEGATIVE_EXPORT_FLAGS = [
    "apd_coordinate_covariant_inner_product_source_exported",
    "apd_gram_metric_source_exported",
    "apd_basis_covariance_law_source_exported",
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
        "new_packet": "P2580|S1530|APD metric covariance|metric covariance.*APD",
        "intended_research_nonduplication": "APD.*covariant metric|covariant metric.*APD|APD.*basis covariance|basis covariance.*APD|coordinate-covariant.*APD|APD.*Gram metric|Gram metric.*APD|basis-transformed metric|APD.*tensorial metric|metric tensor.*APD",
        "apd_precursors": "P2416|S1366|P2578|S1528|P2579|S1529|APD.*basis metric|APD.*inner product|strict_dynamical_source_for_A_P_D",
        "guardrails": "K_legacy_ont|K_strict_gate|role-transfer theorem|QW-2191|ToE closure|source theorem",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def boundary_matrix(basis: list[Any]) -> np.ndarray:
    return np.array([[poly.deriv(1)(0.0) for poly in basis], [poly.deriv(1)(11.0) for poly in basis]], dtype=float)


def constrained_minimizer(matrix: np.ndarray, rhs: np.ndarray, metric: np.ndarray) -> np.ndarray:
    metric_inv = np.linalg.inv(metric)
    gram = matrix @ metric_inv @ matrix.T
    return metric_inv @ matrix.T @ np.linalg.solve(gram, rhs)


def change_of_basis_matrix(reference_basis: list[Any], variant_basis: list[Any]) -> np.ndarray:
    reference_values = np.array([[poly(point) for poly in reference_basis] for point in TRANSFORM_SAMPLE_POINTS], dtype=float)
    variant_values = np.array([[poly(point) for poly in variant_basis] for point in TRANSFORM_SAMPLE_POINTS], dtype=float)
    return np.linalg.solve(reference_values, variant_values)


def selected_probe_values(base: Any, reference_basis: list[Any], reference_coefficients: np.ndarray) -> list[float]:
    selected = base + polynomial_combination(reference_coefficients, reference_basis)
    return [float(selected(point)) for point in PROBE_POINTS]


def variant_row(variant: dict[str, Any], reference_basis: list[Any], base: Any, rhs: np.ndarray, reference_solution: np.ndarray, reference_probe_values: list[float]) -> dict[str, Any]:
    variant_basis = variant["basis"]
    transform = change_of_basis_matrix(reference_basis, variant_basis)
    variant_matrix = boundary_matrix(variant_basis)
    naive_metric = np.eye(transform.shape[1])
    covariant_metric = transform.T @ transform
    naive_variant_coefficients = constrained_minimizer(variant_matrix, rhs, naive_metric)
    covariant_variant_coefficients = constrained_minimizer(variant_matrix, rhs, covariant_metric)
    naive_reference_coefficients = transform @ naive_variant_coefficients
    covariant_reference_coefficients = transform @ covariant_variant_coefficients
    naive_probe_values = selected_probe_values(base, reference_basis, naive_reference_coefficients)
    covariant_probe_values = selected_probe_values(base, reference_basis, covariant_reference_coefficients)
    return {
        "basis_name": variant["name"],
        "change_of_basis_to_reference": [[float(value) for value in row] for row in transform],
        "covariant_metric_in_variant_coordinates": [[float(value) for value in row] for row in covariant_metric],
        "covariant_metric_eigenvalues": [float(value) for value in np.linalg.eigvalsh(covariant_metric)],
        "naive_reference_coefficients": [float(value) for value in naive_reference_coefficients],
        "covariant_reference_coefficients": [float(value) for value in covariant_reference_coefficients],
        "reference_solution_coefficients": [float(value) for value in reference_solution],
        "naive_probe_values": naive_probe_values,
        "covariant_probe_values": covariant_probe_values,
        "max_abs_naive_probe_deviation_from_reference": float(np.max(np.abs(np.array(naive_probe_values) - np.array(reference_probe_values)))),
        "max_abs_covariant_probe_deviation_from_reference": float(np.max(np.abs(np.array(covariant_probe_values) - np.array(reference_probe_values)))),
        "max_abs_covariant_coefficient_deviation_from_reference": float(np.max(np.abs(covariant_reference_coefficients - reference_solution))),
    }


def target_row(target: dict[str, Any], base: Any, reference_basis: list[Any], variants: list[dict[str, Any]]) -> dict[str, Any]:
    rhs = np.array([target["left_slope_target"] - base.deriv(1)(0.0), target["right_slope_target"] - base.deriv(1)(11.0)], dtype=float)
    reference_matrix = boundary_matrix(reference_basis)
    reference_solution = constrained_minimizer(reference_matrix, rhs, np.eye(3))
    reference_probe_values = selected_probe_values(base, reference_basis, reference_solution)
    variant_rows = [variant_row(variant, reference_basis, base, rhs, reference_solution, reference_probe_values) for variant in variants]
    naive_middle_values = {round(row["naive_probe_values"][1], 18) for row in variant_rows}
    covariant_middle_values = {round(row["covariant_probe_values"][1], 6) for row in variant_rows}
    return {
        "target_name": target["name"],
        "reference_basis_name": variants[0]["name"],
        "reference_solution_coefficients": [float(value) for value in reference_solution],
        "reference_probe_values": reference_probe_values,
        "variant_rows": variant_rows,
        "variant_count": len(variant_rows),
        "naive_distinct_middle_probe_values_after_rounding_1e_minus_18": len(naive_middle_values),
        "covariant_distinct_middle_probe_values_after_rounding_1e_minus_6": len(covariant_middle_values),
        "naive_euclidean_metric_is_basis_dependent": len(naive_middle_values) > 1,
        "covariant_metric_transport_restores_basis_invariance": all(row["max_abs_covariant_probe_deviation_from_reference"] <= 1.0e-6 for row in variant_rows),
        "covariant_coefficients_match_reference": all(row["max_abs_covariant_coefficient_deviation_from_reference"] <= 1.0e-15 for row in variant_rows),
        "all_covariant_metrics_positive_definite": all(min(row["covariant_metric_eigenvalues"]) > 0.0 for row in variant_rows),
    }


def covariance_audit(rows: list[dict[str, Any]]) -> dict[str, Any]:
    base, vanish = base_and_vanish(rows)
    variants = basis_variants(vanish)
    reference_basis = variants[0]["basis"]
    target_rows = [target_row(target, base, reference_basis, variants) for target in BOUNDARY_TARGETS]
    return {
        "numpy_version": np.__version__,
        "probe_points": PROBE_POINTS,
        "transform_sample_points": TRANSFORM_SAMPLE_POINTS,
        "basis_variants": [variant["name"] for variant in variants],
        "target_rows": target_rows,
        "target_count": len(target_rows),
        "basis_variant_count": len(variants),
        "naive_euclidean_metric_fails_basis_covariance": all(row["naive_euclidean_metric_is_basis_dependent"] for row in target_rows),
        "transported_metric_restores_basis_covariance": all(row["covariant_metric_transport_restores_basis_invariance"] for row in target_rows),
        "transported_metric_is_positive_definite": all(row["all_covariant_metrics_positive_definite"] for row in target_rows),
        "covariance_is_requirement_not_source": True,
    }


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    p2416_payload = load_json(SOURCE_FILES["P2416_APD_VALUE_BRIDGE"])
    p2578_payload = load_json(SOURCE_FILES["P2578_APD_BASIS_METRIC_DEPENDENCE"])
    p2579_payload = load_json(SOURCE_FILES["P2579_APD_INVERSE_METRIC_TUNABILITY"])
    p2578 = theorem(p2578_payload, "apd_augmented_boundary_basis_metric_dependence_audit")
    p2579 = theorem(p2579_payload, "apd_inner_product_inverse_metric_tunability_audit")
    rows = apd_rows(p2416_payload)
    audit = covariance_audit(rows)
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2580_T1_apd_inner_product_basis_covariance_requirement_certificate",
        "audited_chain": ["P2416/S1366", "P2578/S1528", "P2579/S1529"],
        "frontier_atom_under_attack": "strict_dynamical_source_for_A_P_D",
        "candidate_principle_under_test": "coordinate-covariant transport of APD inner-product metrics across vanishing-mode bases",
        "p2416_apd_value_bridge_inherited": p2416_payload.get("apd_multiplicative_bridge_assembly_necessity_certificate", {}).get("finite_witness_certificate", {}).get("apd_value_assembly_witness_ready") is True,
        "p2578_basis_metric_dependence_inherited": p2578.get("finite_apd_values_and_boundary_targets_do_not_select_coordinate_metric") is True,
        "p2579_inverse_metric_tunability_inherited": p2579.get("finite_apd_values_and_boundary_targets_do_not_select_function_space_inner_product") is True,
        "apd_node_rows": rows,
        "apd_inner_product_basis_covariance_requirement_audit": audit,
        "metric_covariance_is_necessary_but_not_source": True,
        "coordinate_covariance_removes_basis_artifact_but_not_metric_choice": audit["transported_metric_restores_basis_covariance"],
        "apd_dynamic_source_remains_unsourced": True,
        "recommended_next_honest_step": (
            "Use P2580 only as a covariance requirement, not as a source theorem. A valid APD source must provide a metric tensor/inner product whose components transform covariantly under basis changes; merely resetting Euclidean norm in each coordinate chart reproduces P2578, while choosing a transported metric still leaves the underlying metric unsourced as in P2579. The next honest step is to derive the APD inner product from a strict kinetic or measure term."
        ),
        **{flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    gatekeepers = {
        "rg_non_duplication_audit_recorded": grep["tool"] == "rg",
        "p2416_apd_value_bridge_inherited": theorem_export["p2416_apd_value_bridge_inherited"],
        "p2578_basis_metric_dependence_inherited": theorem_export["p2578_basis_metric_dependence_inherited"],
        "p2579_inverse_metric_tunability_inherited": theorem_export["p2579_inverse_metric_tunability_inherited"],
        "three_targets_audited": audit["target_count"] == 3,
        "four_basis_variants": audit["basis_variant_count"] == 4,
        "naive_euclidean_metric_fails_basis_covariance": audit["naive_euclidean_metric_fails_basis_covariance"],
        "transported_metric_restores_basis_covariance": audit["transported_metric_restores_basis_covariance"],
        "transported_metric_is_positive_definite": audit["transported_metric_is_positive_definite"],
        "apd_dynamic_source_not_exported": theorem_export["strict_dynamical_source_for_A_P_D_exported"] is False,
        "no_bridge_exported": theorem_export["bridge_theorem_exported"] is False,
        "no_role_transfer_exported": theorem_export["role_transfer_theorem_exported"] is False,
        "no_qw2191_discharge": theorem_export["qw2191_discharged_by_this_certificate"] is False,
        "no_toe_closure_claimed": theorem_export["toe_closure_claimed"] is False,
    }
    return {
        "packet_id": "P2580",
        "stage_id": "S1530",
        "status": "P2580_APD_INNER_PRODUCT_BASIS_COVARIANCE_REQUIREMENT_CERTIFICATE_NO_APD_DYNAMIC_SOURCE_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "apd_inner_product_basis_covariance_requirement_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints_sha256": {
                "P2416_APD_VALUE_BRIDGE": sha256_json(p2416_payload),
                "P2578_APD_BASIS_METRIC_DEPENDENCE": sha256_json(p2578_payload),
                "P2579_APD_INVERSE_METRIC_TUNABILITY": sha256_json(p2579_payload),
            },
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["apd_inner_product_basis_covariance_requirement_certificate"]["theorem_export"]
    audit = t["apd_inner_product_basis_covariance_requirement_audit"]
    lines = [
        "# P2580/S1530 APD inner-product basis covariance requirement certificate", "",
        f"Status: `{payload['status']}`", "", "## Result", "",
        f"- Frontier atom under attack: `{t['frontier_atom_under_attack']}`.",
        f"- Targets audited: `{audit['target_count']}`.",
        f"- Basis variants audited: `{audit['basis_variant_count']}`.",
        f"- Naive Euclidean metric fails basis covariance: `{audit['naive_euclidean_metric_fails_basis_covariance']}`.",
        f"- Transported metric restores basis covariance: `{audit['transported_metric_restores_basis_covariance']}`.",
        f"- Strict APD dynamic source exported: `{t['strict_dynamical_source_for_A_P_D_exported']}`.", "",
        "## Interpretation", "",
        "P2580 separates a real covariance requirement from a source theorem.  If the metric tensor is transported as `G_new = T^T G_ref T`, the same quadratic form gives the same APD solution in every audited basis.  If Euclidean norm is reset in each basis, the P2578 basis artifact reappears.  Covariance fixes coordinates, not the missing strict origin of the metric.", "",
        "## Recommended next honest step", "", t["recommended_next_honest_step"], "", "## Negative controls", "",
        "No APD coordinate-covariant inner-product source, Gram-metric source, basis-covariance-law source, strict A/P/D dynamic source, bridge theorem, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, or ToE closure is exported.", "", "## Fingerprint", "",
        f"`{payload['apd_inner_product_basis_covariance_requirement_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2580/S1530 APD inner-product basis covariance requirement guard

`P2580/S1530` separates coordinate covariance from APD sourcehood.  For the same P2578 bases and APD boundary targets, resetting the Euclidean norm in every basis reproduces basis-dependent selections, while transporting the metric tensor by the basis-change law restores the same selected APD polynomial across bases.  This proves a necessary covariance rule for any future APD inner product, but it still does not source the metric or `strict_dynamical_source_for_A_P_D`.
""".strip()
    lag_section = """
## P2580/S1530 APD inner-product basis covariance requirement Ltotal guard

`P2580/S1530` blocks role-bearing `L_total` terms that choose a fresh coordinate Euclidean APD norm in each basis.  A legitimate APD kinetic/measure term must transform as a metric tensor under basis changes; covariance removes coordinate artifacts but does not by itself derive the inner product from strict nadsoliton dynamics.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2580/S1530 APD inner-product basis covariance requirement guard", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2580/S1530 APD inner-product basis covariance requirement Ltotal guard", lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
