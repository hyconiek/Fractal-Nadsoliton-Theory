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
OUT = GEN / "p2574_s1524_apd_two_endpoint_boundary_compatibility_obstruction_audit.json"
MD = GEN / "p2574_s1524_apd_two_endpoint_boundary_compatibility_obstruction_audit.md"

SOURCE_FILES = {
    "P2416_APD_VALUE_BRIDGE": GEN / "p2416_s1366_apd_multiplicative_bridge_assembly_necessity_certificate.json",
    "P2572_APD_BOUNDARY_PENALTY_CONTINUUM": GEN / "p2572_s1522_apd_boundary_penalty_selector_continuum_audit.json",
    "P2573_APD_INVERSE_BOUNDARY_TUNABILITY": GEN / "p2573_s1523_apd_boundary_penalty_inverse_target_tunability_audit.json",
}
DOMAIN = list(range(12))
PROBE_POINTS = [0.5, 5.5, 11.5]
BOUNDARY_TARGETS = [
    {"name": "zero_zero_neumann", "left_slope_target": 0.0, "right_slope_target": 0.0},
    {"name": "equal_positive_slopes", "left_slope_target": 1.0e-6, "right_slope_target": 1.0e-6},
    {"name": "equal_negative_slopes", "left_slope_target": -1.0e-6, "right_slope_target": -1.0e-6},
]
NEGATIVE_EXPORT_FLAGS = [
    "apd_two_endpoint_boundary_source_exported",
    "apd_neumann_boundary_source_exported",
    "apd_boundary_compatibility_source_exported",
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
        "new_packet": "P2574|S1524|APD boundary compatibility|two-boundary compatibility",
        "intended_research_nonduplication": "APD.*Neumann|Neumann.*APD|APD.*slope compatibility|slope compatibility.*APD|APD.*two endpoint|two endpoint.*APD|APD.*boundary obstruction|boundary.*compatibility.*APD",
        "apd_precursors": "P2416|S1366|P2572|S1522|P2573|S1523|APD.*boundary|APD.*penalty|strict_dynamical_source_for_A_P_D",
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


def boundary_coefficients(base: Polynomial, vanish: Polynomial) -> dict[str, float]:
    slope_base = base.deriv(1)
    slope_vanish = vanish.deriv(1)
    return {
        "base_left_slope": float(slope_base(0.0)),
        "base_right_slope": float(slope_base(11.0)),
        "vanish_left_slope": float(slope_vanish(0.0)),
        "vanish_right_slope": float(slope_vanish(11.0)),
    }


def compatibility_row(target: dict[str, Any], coeffs: dict[str, float], base: Polynomial, vanish: Polynomial, rows: list[dict[str, Any]]) -> dict[str, Any]:
    left_lambda = (target["left_slope_target"] - coeffs["base_left_slope"]) / coeffs["vanish_left_slope"]
    right_lambda = (target["right_slope_target"] - coeffs["base_right_slope"]) / coeffs["vanish_right_slope"]
    design = np.array([[coeffs["vanish_left_slope"]], [coeffs["vanish_right_slope"]]], dtype=float)
    rhs = np.array([
        target["left_slope_target"] - coeffs["base_left_slope"],
        target["right_slope_target"] - coeffs["base_right_slope"],
    ], dtype=float)
    least_squares_lambda = float(np.linalg.lstsq(design, rhs, rcond=None)[0][0])
    selected = base + least_squares_lambda * vanish
    node_points = np.array(DOMAIN, dtype=float)
    node_values = np.array([row["apd_value"] for row in rows], dtype=float)
    residuals = selected(node_points) - node_values
    achieved_left = coeffs["base_left_slope"] + least_squares_lambda * coeffs["vanish_left_slope"]
    achieved_right = coeffs["base_right_slope"] + least_squares_lambda * coeffs["vanish_right_slope"]
    probe_deltas = [float(least_squares_lambda * vanish(point)) for point in PROBE_POINTS]
    return {
        "target_name": target["name"],
        "left_slope_target": target["left_slope_target"],
        "right_slope_target": target["right_slope_target"],
        "lambda_required_by_left_endpoint": float(left_lambda),
        "lambda_required_by_right_endpoint": float(right_lambda),
        "endpoint_lambda_gap": float(abs(left_lambda - right_lambda)),
        "two_endpoint_exact_compatibility_holds": abs(left_lambda - right_lambda) <= 1.0e-18,
        "least_squares_lambda": least_squares_lambda,
        "least_squares_left_slope": float(achieved_left),
        "least_squares_right_slope": float(achieved_right),
        "least_squares_max_abs_boundary_error": float(max(abs(achieved_left - target["left_slope_target"]), abs(achieved_right - target["right_slope_target"]))),
        "max_abs_node_residual_at_least_squares_lambda": float(np.max(np.abs(residuals))),
        "probe_value_deltas_from_base": probe_deltas,
    }


def compatibility_audit(rows: list[dict[str, Any]]) -> dict[str, Any]:
    base, vanish = interpolation_polynomials(rows)
    coeffs = boundary_coefficients(base, vanish)
    audit_rows = [compatibility_row(target, coeffs, base, vanish, rows) for target in BOUNDARY_TARGETS]
    return {
        "numpy_version": np.__version__,
        "domain": DOMAIN,
        "probe_points": PROBE_POINTS,
        "boundary_slope_formula": "q_lambda'(endpoint)=q_interp'(endpoint)+lambda*V'(endpoint)",
        "boundary_coefficients": coeffs,
        "target_rows": audit_rows,
        "target_count": len(audit_rows),
        "all_targets_fail_exact_two_endpoint_compatibility": all(not row["two_endpoint_exact_compatibility_holds"] for row in audit_rows),
        "all_least_squares_members_preserve_apd_nodes": all(row["max_abs_node_residual_at_least_squares_lambda"] <= 1.0e-6 for row in audit_rows),
        "zero_zero_neumann_incompatible": next(row for row in audit_rows if row["target_name"] == "zero_zero_neumann")["two_endpoint_exact_compatibility_holds"] is False,
        "two_endpoint_boundary_conditions_are_extra_source_obligation": True,
    }


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    p2416_payload = load_json(SOURCE_FILES["P2416_APD_VALUE_BRIDGE"])
    p2572_payload = load_json(SOURCE_FILES["P2572_APD_BOUNDARY_PENALTY_CONTINUUM"])
    p2573_payload = load_json(SOURCE_FILES["P2573_APD_INVERSE_BOUNDARY_TUNABILITY"])
    p2572 = theorem(p2572_payload, "apd_boundary_penalty_selector_continuum_audit")
    p2573 = theorem(p2573_payload, "apd_boundary_penalty_inverse_target_tunability_audit")
    rows = apd_rows(p2416_payload)
    audit = compatibility_audit(rows)
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2574_T1_apd_two_endpoint_boundary_compatibility_obstruction_audit",
        "audited_chain": ["P2416/S1366", "P2572/S1522", "P2573/S1523"],
        "frontier_atom_under_attack": "strict_dynamical_source_for_A_P_D",
        "candidate_principle_under_test": "impose two endpoint slope boundary targets on the one-parameter APD interpolation family",
        "p2416_apd_value_bridge_inherited": p2416_payload.get("apd_multiplicative_bridge_assembly_necessity_certificate", {}).get("finite_witness_certificate", {}).get("apd_value_assembly_witness_ready") is True,
        "p2572_boundary_penalty_continuum_inherited": p2572.get("finite_apd_values_do_not_select_boundary_penalty") is True,
        "p2573_inverse_boundary_tunability_inherited": p2573.get("finite_apd_values_do_not_select_target_or_penalty_law") is True,
        "apd_node_rows": rows,
        "two_endpoint_boundary_compatibility_audit": audit,
        "finite_apd_values_do_not_supply_compatible_two_endpoint_boundary_law": audit["all_targets_fail_exact_two_endpoint_compatibility"],
        "zero_neumann_pair_is_not_available_on_this_family": audit["zero_zero_neumann_incompatible"],
        "apd_dynamic_source_remains_unsourced": True,
        "recommended_next_honest_step": (
            "Do not impose attractive two-endpoint APD boundary conditions by analogy. The next honest step is to derive admissible APD boundary data from the strict action and then test compatibility with the interpolation family; common Neumann-style targets are not automatically available on the one-parameter bridge family."
        ),
        **{flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    gatekeepers = {
        "rg_non_duplication_audit_recorded": grep["tool"] == "rg",
        "p2416_apd_value_bridge_inherited": theorem_export["p2416_apd_value_bridge_inherited"],
        "p2572_boundary_penalty_continuum_inherited": theorem_export["p2572_boundary_penalty_continuum_inherited"],
        "p2573_inverse_boundary_tunability_inherited": theorem_export["p2573_inverse_boundary_tunability_inherited"],
        "twelve_apd_rows": len(rows) == 12,
        "three_boundary_targets": audit["target_count"] == 3,
        "all_targets_fail_exact_two_endpoint_compatibility": audit["all_targets_fail_exact_two_endpoint_compatibility"],
        "zero_neumann_pair_incompatible": audit["zero_zero_neumann_incompatible"],
        "apd_dynamic_source_not_exported": theorem_export["strict_dynamical_source_for_A_P_D_exported"] is False,
        "no_bridge_exported": theorem_export["bridge_theorem_exported"] is False,
        "no_role_transfer_exported": theorem_export["role_transfer_theorem_exported"] is False,
        "no_qw2191_discharge": theorem_export["qw2191_discharged_by_this_certificate"] is False,
        "no_toe_closure_claimed": theorem_export["toe_closure_claimed"] is False,
    }
    return {
        "packet_id": "P2574",
        "stage_id": "S1524",
        "status": "P2574_APD_TWO_ENDPOINT_BOUNDARY_COMPATIBILITY_OBSTRUCTION_AUDIT_NO_APD_DYNAMIC_SOURCE_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "apd_two_endpoint_boundary_compatibility_obstruction_audit": {
            "theorem_export": theorem_export,
            "source_fingerprints_sha256": {
                "P2416_APD_VALUE_BRIDGE": sha256_json(p2416_payload),
                "P2572_APD_BOUNDARY_PENALTY_CONTINUUM": sha256_json(p2572_payload),
                "P2573_APD_INVERSE_BOUNDARY_TUNABILITY": sha256_json(p2573_payload),
            },
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["apd_two_endpoint_boundary_compatibility_obstruction_audit"]["theorem_export"]
    audit = t["two_endpoint_boundary_compatibility_audit"]
    lines = [
        "# P2574/S1524 APD two-endpoint boundary compatibility obstruction audit", "",
        f"Status: `{payload['status']}`", "", "## Result", "",
        f"- Frontier atom under attack: `{t['frontier_atom_under_attack']}`.",
        f"- Boundary slope formula: `{audit['boundary_slope_formula']}`.",
        f"- Boundary targets audited: `{audit['target_count']}`.",
        f"- All targets fail exact two-endpoint compatibility: `{audit['all_targets_fail_exact_two_endpoint_compatibility']}`.",
        f"- Zero/zero Neumann target incompatible: `{audit['zero_zero_neumann_incompatible']}`.",
        f"- Finite APD values supply compatible two-endpoint law: `{not t['finite_apd_values_do_not_supply_compatible_two_endpoint_boundary_law']}`.", "",
        "## Interpretation", "",
        "The one-parameter APD family can satisfy at most one independent endpoint slope target unless the left-required and right-required `lambda` values agree.  In the audited targets, including zero/zero Neumann slopes, the required endpoint lambdas differ.  Thus attractive two-endpoint boundary slogans are compatibility constraints, not strict APD source theorems.", "",
        "## Recommended next honest step", "", t["recommended_next_honest_step"], "", "## Negative controls", "",
        "No APD two-endpoint boundary source, Neumann boundary source, strict A/P/D dynamic source, bridge theorem, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, or ToE closure is exported.", "", "## Fingerprint", "",
        f"`{payload['apd_two_endpoint_boundary_compatibility_obstruction_audit']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
`P2574/S1524` audits whether simple two-endpoint APD slope boundary laws can replace the missing source.  On `q_lambda=q_interp+lambda*V`, each endpoint slope target demands its own `lambda`; for the audited target pairs, including zero/zero Neumann data, the left-required and right-required `lambda` values disagree.  Thus two-endpoint boundary slogans are compatibility constraints, not `strict_dynamical_source_for_A_P_D`.
""".strip()
    lag_section = """
`P2574/S1524` blocks inserting two-endpoint Neumann-style APD boundary conditions into role-bearing `L_total` by analogy.  Such boundary data must be derived and compatibility-tested; the audited one-parameter bridge family does not automatically satisfy common two-endpoint slope targets.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2574/S1524 APD two-endpoint boundary compatibility guard", "## P2574/S1524 APD two-endpoint boundary compatibility guard\n\n" + eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2574/S1524 APD two-endpoint boundary compatibility Ltotal guard", "## P2574/S1524 APD two-endpoint boundary compatibility Ltotal guard\n\n" + lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
