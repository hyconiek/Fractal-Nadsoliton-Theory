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
OUT = GEN / "p2573_s1523_apd_boundary_penalty_inverse_target_tunability_audit.json"
MD = GEN / "p2573_s1523_apd_boundary_penalty_inverse_target_tunability_audit.md"

SOURCE_FILES = {
    "P2416_APD_VALUE_BRIDGE": GEN / "p2416_s1366_apd_multiplicative_bridge_assembly_necessity_certificate.json",
    "P2571_APD_MEASURE_BOUNDARY_DEPENDENCE": GEN / "p2571_s1521_apd_sobolev_measure_boundary_dependence_audit.json",
    "P2572_APD_BOUNDARY_PENALTY_CONTINUUM": GEN / "p2572_s1522_apd_boundary_penalty_selector_continuum_audit.json",
}
DOMAIN = list(range(12))
FIXED_DERIVATIVE_ORDER = 2
TARGET_FRACTIONS = [0.0, 0.25, 0.5, 0.75, 0.9]
PROBE_POINTS = [0.5, 5.5, 11.5]
NEGATIVE_EXPORT_FLAGS = [
    "apd_boundary_penalty_law_source_exported",
    "apd_inverse_boundary_selector_source_exported",
    "apd_target_lambda_source_exported",
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
        "new_packet": "P2573|S1523|APD boundary penalty inverse|boundary penalty inverse",
        "intended_research_nonduplication": "APD.*inverse.*penalty|inverse.*boundary penalty.*APD|APD.*target lambda|target lambda.*APD|post-hoc boundary penalty|boundary penalty.*target|APD.*tunability|selector tunability.*APD",
        "apd_precursors": "P2416|S1366|P2571|S1521|P2572|S1522|APD.*boundary|APD.*penalty|strict_dynamical_source_for_A_P_D",
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


def stationarity_coefficients(base: Polynomial, vanish: Polynomial) -> dict[str, float]:
    d_base = base.deriv(FIXED_DERIVATIVE_ORDER)
    d_vanish = vanish.deriv(FIXED_DERIVATIVE_ORDER)
    slope_base = base.deriv(1)
    slope_vanish = vanish.deriv(1)
    coeffs = {
        "A_bulk": definite_integral(d_vanish * d_vanish),
        "B_bulk": definite_integral(d_base * d_vanish),
        "left_A_endpoint": float(slope_vanish(0.0)) ** 2,
        "left_B_endpoint": float(slope_base(0.0)) * float(slope_vanish(0.0)),
        "right_A_endpoint": float(slope_vanish(11.0)) ** 2,
        "right_B_endpoint": float(slope_base(11.0)) * float(slope_vanish(11.0)),
    }
    coeffs["uniform_bulk_lambda"] = -coeffs["B_bulk"] / coeffs["A_bulk"]
    coeffs["left_infinite_penalty_limit_lambda"] = -coeffs["left_B_endpoint"] / coeffs["left_A_endpoint"]
    coeffs["right_infinite_penalty_limit_lambda"] = -coeffs["right_B_endpoint"] / coeffs["right_A_endpoint"]
    return coeffs


def lambda_from_penalty(coeffs: dict[str, float], side: str, penalty: float) -> float:
    endpoint_a = coeffs[f"{side}_A_endpoint"]
    endpoint_b = coeffs[f"{side}_B_endpoint"]
    return -(coeffs["B_bulk"] + penalty * endpoint_b) / (coeffs["A_bulk"] + penalty * endpoint_a)


def penalty_for_target(coeffs: dict[str, float], side: str, target_lambda: float) -> float:
    endpoint_a = coeffs[f"{side}_A_endpoint"]
    endpoint_b = coeffs[f"{side}_B_endpoint"]
    return -(target_lambda * coeffs["A_bulk"] + coeffs["B_bulk"]) / (target_lambda * endpoint_a + endpoint_b)


def inverse_row(side: str, fraction: float, coeffs: dict[str, float], base: Polynomial, vanish: Polynomial, rows: list[dict[str, Any]]) -> dict[str, Any]:
    start = coeffs["uniform_bulk_lambda"]
    limit = coeffs[f"{side}_infinite_penalty_limit_lambda"]
    target_lambda = start + fraction * (limit - start)
    penalty = penalty_for_target(coeffs, side, target_lambda)
    recovered_lambda = lambda_from_penalty(coeffs, side, penalty)
    selected = base + recovered_lambda * vanish
    node_points = np.array(DOMAIN, dtype=float)
    node_values = np.array([row["apd_value"] for row in rows], dtype=float)
    residuals = selected(node_points) - node_values
    probe_deltas = [float(recovered_lambda * vanish(point)) for point in PROBE_POINTS]
    return {
        "side": side,
        "fraction_from_uniform_to_infinite_penalty_limit": fraction,
        "target_lambda": float(target_lambda),
        "solved_endpoint_penalty": float(penalty),
        "recovered_lambda": float(recovered_lambda),
        "lambda_abs_error": float(abs(recovered_lambda - target_lambda)),
        "penalty_is_nonnegative": penalty >= -1.0e-12,
        "max_abs_node_residual_at_recovered_lambda": float(np.max(np.abs(residuals))),
        "probe_value_deltas_from_base": probe_deltas,
        "selected_member_changes_off_node_values": any(abs(value) > 0.0 for value in probe_deltas),
    }


def inverse_tunability_audit(rows: list[dict[str, Any]]) -> dict[str, Any]:
    base, vanish = interpolation_polynomials(rows)
    coeffs = stationarity_coefficients(base, vanish)
    left_rows = [inverse_row("left", fraction, coeffs, base, vanish, rows) for fraction in TARGET_FRACTIONS]
    right_rows = [inverse_row("right", fraction, coeffs, base, vanish, rows) for fraction in TARGET_FRACTIONS]
    all_rows = left_rows + right_rows
    return {
        "numpy_version": np.__version__,
        "fixed_derivative_order": FIXED_DERIVATIVE_ORDER,
        "domain": DOMAIN,
        "target_fractions": TARGET_FRACTIONS,
        "probe_points": PROBE_POINTS,
        "inverse_stationarity_formula": "t(lambda_target)=-(lambda_target*A_bulk+B_bulk)/(lambda_target*A_endpoint+B_endpoint)",
        "stationarity_coefficients": coeffs,
        "left_inverse_rows": left_rows,
        "right_inverse_rows": right_rows,
        "inverse_rows_count": len(all_rows),
        "all_targets_recovered_with_small_error": all(row["lambda_abs_error"] <= 1.0e-20 for row in all_rows),
        "all_solved_penalties_nonnegative": all(row["penalty_is_nonnegative"] for row in all_rows),
        "all_inverse_rows_preserve_apd_nodes": all(row["max_abs_node_residual_at_recovered_lambda"] <= 1.0e-6 for row in all_rows),
        "inverse_boundary_penalty_tunes_selector_post_hoc": True,
    }


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    p2416_payload = load_json(SOURCE_FILES["P2416_APD_VALUE_BRIDGE"])
    p2571_payload = load_json(SOURCE_FILES["P2571_APD_MEASURE_BOUNDARY_DEPENDENCE"])
    p2572_payload = load_json(SOURCE_FILES["P2572_APD_BOUNDARY_PENALTY_CONTINUUM"])
    p2571 = theorem(p2571_payload, "apd_sobolev_measure_boundary_dependence_audit")
    p2572 = theorem(p2572_payload, "apd_boundary_penalty_selector_continuum_audit")
    rows = apd_rows(p2416_payload)
    audit = inverse_tunability_audit(rows)
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2573_T1_apd_boundary_penalty_inverse_target_tunability_audit",
        "audited_chain": ["P2416/S1366", "P2571/S1521", "P2572/S1522"],
        "frontier_atom_under_attack": "strict_dynamical_source_for_A_P_D",
        "candidate_principle_under_test": "given a desired APD lambda in the endpoint-penalty interval, solve the boundary penalty t that realizes it",
        "p2416_apd_value_bridge_inherited": p2416_payload.get("apd_multiplicative_bridge_assembly_necessity_certificate", {}).get("finite_witness_certificate", {}).get("apd_value_assembly_witness_ready") is True,
        "p2571_measure_boundary_obligation_inherited": p2571.get("finite_apd_values_do_not_select_measure_or_boundary_class") is True,
        "p2572_boundary_penalty_continuum_inherited": p2572.get("finite_apd_values_do_not_select_boundary_penalty") is True,
        "apd_node_rows": rows,
        "boundary_penalty_inverse_target_tunability_audit": audit,
        "finite_apd_values_do_not_select_target_or_penalty_law": audit["inverse_boundary_penalty_tunes_selector_post_hoc"],
        "inverse_penalty_fit_is_not_source_theorem": True,
        "apd_dynamic_source_remains_unsourced": True,
        "recommended_next_honest_step": (
            "Do not fit boundary penalties inversely to obtain a desired APD interpolation. The next honest step is to derive admissible boundary conditions and penalty strengths before choosing the APD member; otherwise the boundary term is a post-hoc selector and not strict A/P/D dynamics."
        ),
        **{flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    gatekeepers = {
        "rg_non_duplication_audit_recorded": grep["tool"] == "rg",
        "p2416_apd_value_bridge_inherited": theorem_export["p2416_apd_value_bridge_inherited"],
        "p2571_measure_boundary_obligation_inherited": theorem_export["p2571_measure_boundary_obligation_inherited"],
        "p2572_boundary_penalty_continuum_inherited": theorem_export["p2572_boundary_penalty_continuum_inherited"],
        "twelve_apd_rows": len(rows) == 12,
        "ten_inverse_rows": audit["inverse_rows_count"] == 10,
        "all_targets_recovered": audit["all_targets_recovered_with_small_error"],
        "all_solved_penalties_nonnegative": audit["all_solved_penalties_nonnegative"],
        "all_inverse_rows_preserve_nodes": audit["all_inverse_rows_preserve_apd_nodes"],
        "apd_dynamic_source_not_exported": theorem_export["strict_dynamical_source_for_A_P_D_exported"] is False,
        "no_bridge_exported": theorem_export["bridge_theorem_exported"] is False,
        "no_role_transfer_exported": theorem_export["role_transfer_theorem_exported"] is False,
        "no_qw2191_discharge": theorem_export["qw2191_discharged_by_this_certificate"] is False,
        "no_toe_closure_claimed": theorem_export["toe_closure_claimed"] is False,
    }
    return {
        "packet_id": "P2573",
        "stage_id": "S1523",
        "status": "P2573_APD_BOUNDARY_PENALTY_INVERSE_TARGET_TUNABILITY_AUDIT_NO_APD_DYNAMIC_SOURCE_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "apd_boundary_penalty_inverse_target_tunability_audit": {
            "theorem_export": theorem_export,
            "source_fingerprints_sha256": {
                "P2416_APD_VALUE_BRIDGE": sha256_json(p2416_payload),
                "P2571_APD_MEASURE_BOUNDARY_DEPENDENCE": sha256_json(p2571_payload),
                "P2572_APD_BOUNDARY_PENALTY_CONTINUUM": sha256_json(p2572_payload),
            },
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["apd_boundary_penalty_inverse_target_tunability_audit"]["theorem_export"]
    audit = t["boundary_penalty_inverse_target_tunability_audit"]
    lines = [
        "# P2573/S1523 APD boundary-penalty inverse-target tunability audit", "",
        f"Status: `{payload['status']}`", "", "## Result", "",
        f"- Frontier atom under attack: `{t['frontier_atom_under_attack']}`.",
        f"- Candidate principle: `{t['candidate_principle_under_test']}`.",
        f"- Inverse formula: `{audit['inverse_stationarity_formula']}`.",
        f"- Inverse rows: `{audit['inverse_rows_count']}`.",
        f"- All targets recovered: `{audit['all_targets_recovered_with_small_error']}`.",
        f"- All solved penalties nonnegative: `{audit['all_solved_penalties_nonnegative']}`.",
        f"- All inverse rows preserve APD nodes: `{audit['all_inverse_rows_preserve_apd_nodes']}`.",
        f"- Finite APD values select target/penalty law: `{not t['finite_apd_values_do_not_select_target_or_penalty_law']}`.", "",
        "## Interpretation", "",
        "The P2572 stationarity law can be inverted: for any audited target `lambda` between the uniform-bulk minimizer and the endpoint-penalty limit, the boundary penalty `t` is solved explicitly and recovers the target while preserving all finite APD nodes.  This makes endpoint penalties post-hoc tunable unless their law is derived before selecting the APD dynamics.", "",
        "## Recommended next honest step", "", t["recommended_next_honest_step"], "", "## Negative controls", "",
        "No APD boundary penalty law source, inverse boundary selector source, target-lambda source, strict A/P/D dynamic source, bridge theorem, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, or ToE closure is exported.", "", "## Fingerprint", "",
        f"`{payload['apd_boundary_penalty_inverse_target_tunability_audit']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
`P2573/S1523` inverts the P2572 endpoint-penalty stationarity law.  For targets between the uniform-bulk minimizer and an endpoint-penalty limit, `t(lambda_target)=-(lambda_target*A_bulk+B_bulk)/(lambda_target*A_endpoint+B_endpoint)` recovers the chosen APD `lambda` while preserving all finite APD nodes.  Thus endpoint penalties are post-hoc tunable selectors unless the boundary law is strict-sourced before choosing the APD dynamics.
""".strip()
    lag_section = """
`P2573/S1523` blocks inverse-fitting endpoint penalties into role-bearing `L_total`.  A boundary term that can be solved after choosing a target APD member is not a strict APD action source; admissible boundary conditions and penalty strengths must be derived independently.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2573/S1523 APD boundary-penalty inverse-target tunability guard", "## P2573/S1523 APD boundary-penalty inverse-target tunability guard\n\n" + eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2573/S1523 APD boundary-penalty inverse-target tunability Ltotal guard", "## P2573/S1523 APD boundary-penalty inverse-target tunability Ltotal guard\n\n" + lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
