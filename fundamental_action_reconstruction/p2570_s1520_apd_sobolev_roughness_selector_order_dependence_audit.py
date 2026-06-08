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
OUT = GEN / "p2570_s1520_apd_sobolev_roughness_selector_order_dependence_audit.json"
MD = GEN / "p2570_s1520_apd_sobolev_roughness_selector_order_dependence_audit.md"

SOURCE_FILES = {
    "P2416_APD_VALUE_BRIDGE": GEN / "p2416_s1366_apd_multiplicative_bridge_assembly_necessity_certificate.json",
    "P2561_POST_DAMPING_RESIDUAL_BRIDGE": GEN / "p2561_s1511_strict_completion_post_damping_residual_bridge_two_key_certificate.json",
    "P2569_APD_INTERPOLATION_NONUNIQUENESS": GEN / "p2569_s1519_apd_value_bridge_interpolation_dynamic_nonuniqueness_certificate.json",
}
DOMAIN = list(range(12))
DERIVATIVE_ORDERS = [0, 1, 2, 3, 4, 12]
PROBE_POINTS = [0.5, 5.5, 11.5]
NEGATIVE_EXPORT_FLAGS = [
    "apd_sobolev_order_source_exported",
    "apd_roughness_selector_source_exported",
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
        "new_packet": "P2570|S1520|APD Sobolev roughness|APD roughness selector order|A/P/D roughness selector",
        "intended_research_nonduplication": "APD.*roughness|roughness.*APD|Sobolev.*APD|APD.*regularity|minimal-action.*APD|vanishing polynomial.*roughness|regularity.*vanishing polynomial",
        "apd_precursors": "P2416|S1366|P2569|S1519|APD value bridge|APD interpolation|strict_dynamical_source_for_A_P_D",
        "frontier_precursors": "P2561|S1511|post-damping residual bridge|strict_phase_frequency_source|strict_damping_beta_eta_source",
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


def roughness_rows(rows: list[dict[str, Any]]) -> dict[str, Any]:
    base, vanish = interpolation_polynomials(rows)
    node_points = np.array(DOMAIN, dtype=float)
    node_values = np.array([row["apd_value"] for row in rows], dtype=float)
    audit_rows = []
    for order in DERIVATIVE_ORDERS:
        d_base = base.deriv(order)
        d_vanish = vanish.deriv(order)
        a_quad = definite_integral(d_vanish * d_vanish)
        b_cross = definite_integral(d_base * d_vanish)
        c_base = definite_integral(d_base * d_base)
        lambda_star = -b_cross / a_quad if abs(a_quad) > 0.0 else 0.0
        selected = base + lambda_star * vanish
        residuals = selected(node_points) - node_values
        base_objective = c_base
        selected_objective = c_base + 2.0 * b_cross * lambda_star + a_quad * lambda_star * lambda_star
        probe_delta = [float(lambda_star * vanish(point)) for point in PROBE_POINTS]
        audit_rows.append({
            "derivative_order": order,
            "objective": f"J_{order}(lambda)=int_0^11 (d^{order}/dx^{order} q_lambda(x))^2 dx",
            "quadratic_A": a_quad,
            "linear_cross_B": b_cross,
            "base_cost_C": c_base,
            "lambda_star": float(lambda_star),
            "lambda_star_abs": abs(float(lambda_star)),
            "selects_base_lambda_zero": abs(float(lambda_star)) <= 1.0e-18,
            "selected_cost": float(selected_objective),
            "base_cost": float(base_objective),
            "cost_improvement_over_base": float(base_objective - selected_objective),
            "max_abs_node_residual_at_selected_lambda": float(np.max(np.abs(residuals))),
            "probe_value_deltas_from_base": probe_delta,
            "selected_member_changes_off_node_values": any(abs(value) > 0.0 for value in probe_delta),
        })
    nonzero_orders = [row["derivative_order"] for row in audit_rows if not row["selects_base_lambda_zero"]]
    zero_orders = [row["derivative_order"] for row in audit_rows if row["selects_base_lambda_zero"]]
    rounded_lambdas = {round(row["lambda_star"], 24) for row in audit_rows}
    return {
        "numpy_version": np.__version__,
        "domain": DOMAIN,
        "derivative_orders": DERIVATIVE_ORDERS,
        "probe_points": PROBE_POINTS,
        "roughness_audit_rows": audit_rows,
        "orders_selecting_nonzero_lambda": nonzero_orders,
        "orders_selecting_base_lambda_zero": zero_orders,
        "distinct_minimizers_after_rounding_1e_minus_24": len(rounded_lambdas),
        "all_selected_members_preserve_apd_nodes": all(row["max_abs_node_residual_at_selected_lambda"] <= 1.0e-6 for row in audit_rows),
        "roughness_order_changes_selector": len(rounded_lambdas) > 1 and len(nonzero_orders) > 0 and len(zero_orders) > 0,
        "regularity_principle_not_sourced_by_finite_apd_values": True,
    }


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    p2416_payload = load_json(SOURCE_FILES["P2416_APD_VALUE_BRIDGE"])
    p2561_payload = load_json(SOURCE_FILES["P2561_POST_DAMPING_RESIDUAL_BRIDGE"])
    p2569_payload = load_json(SOURCE_FILES["P2569_APD_INTERPOLATION_NONUNIQUENESS"])
    p2561 = theorem(p2561_payload, "strict_completion_post_damping_residual_bridge_two_key_certificate")
    p2569 = theorem(p2569_payload, "apd_value_bridge_interpolation_dynamic_nonuniqueness_certificate")
    rows = apd_rows(p2416_payload)
    roughness = roughness_rows(rows)
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2570_T1_apd_sobolev_roughness_selector_order_dependence_audit",
        "audited_chain": ["P2416/S1366", "P2561/S1511", "P2569/S1519"],
        "frontier_atom_under_attack": "strict_dynamical_source_for_A_P_D",
        "candidate_principle_under_test": "choose q_lambda minimizing a Sobolev roughness J_k(lambda)=int_0^11 |d^k q_lambda/dx^k|^2 dx on the P2569 interpolation family",
        "p2416_apd_value_bridge_inherited": p2416_payload.get("apd_multiplicative_bridge_assembly_necessity_certificate", {}).get("finite_witness_certificate", {}).get("apd_value_assembly_witness_ready") is True,
        "p2561_apd_residual_atom_inherited": "strict_dynamical_source_for_A_P_D" in p2561.get("residual_atoms_after_hypothetical_damping_source", []),
        "p2569_interpolation_family_inherited": p2569.get("finite_apd_values_do_not_select_dynamic_law") is True,
        "apd_node_rows": rows,
        "sobolev_roughness_order_dependence_audit": roughness,
        "finite_apd_values_do_not_select_sobolev_order": roughness["roughness_order_changes_selector"],
        "sobolev_selector_is_conditional_not_strict_source": True,
        "apd_dynamic_source_remains_unsourced": True,
        "recommended_next_honest_step": (
            "Do not promote a minimum-roughness APD interpolation to strict dynamics unless the derivative order, measure, and boundary class are derived from nadsoliton dynamics. The next honest step is to derive a sourced APD variational functional; if no such functional is available, the A/P/D residual atom remains open even when finite APD values are exact."
        ),
        **{flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    gatekeepers = {
        "rg_non_duplication_audit_recorded": grep["tool"] == "rg",
        "p2416_apd_value_bridge_inherited": theorem_export["p2416_apd_value_bridge_inherited"],
        "p2561_apd_atom_inherited": theorem_export["p2561_apd_residual_atom_inherited"],
        "p2569_interpolation_family_inherited": theorem_export["p2569_interpolation_family_inherited"],
        "twelve_apd_rows": len(rows) == 12,
        "six_derivative_orders_audited": len(roughness["roughness_audit_rows"]) == 6,
        "all_selected_members_preserve_nodes": roughness["all_selected_members_preserve_apd_nodes"],
        "roughness_order_changes_selector": roughness["roughness_order_changes_selector"],
        "apd_dynamic_source_not_exported": theorem_export["strict_dynamical_source_for_A_P_D_exported"] is False,
        "no_bridge_exported": theorem_export["bridge_theorem_exported"] is False,
        "no_role_transfer_exported": theorem_export["role_transfer_theorem_exported"] is False,
        "no_qw2191_discharge": theorem_export["qw2191_discharged_by_this_certificate"] is False,
        "no_toe_closure_claimed": theorem_export["toe_closure_claimed"] is False,
    }
    return {
        "packet_id": "P2570",
        "stage_id": "S1520",
        "status": "P2570_APD_SOBOLEV_ROUGHNESS_SELECTOR_ORDER_DEPENDENCE_AUDIT_NO_APD_DYNAMIC_SOURCE_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "apd_sobolev_roughness_selector_order_dependence_audit": {
            "theorem_export": theorem_export,
            "source_fingerprints_sha256": {
                "P2416_APD_VALUE_BRIDGE": sha256_json(p2416_payload),
                "P2561_POST_DAMPING_RESIDUAL_BRIDGE": sha256_json(p2561_payload),
                "P2569_APD_INTERPOLATION_NONUNIQUENESS": sha256_json(p2569_payload),
            },
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["apd_sobolev_roughness_selector_order_dependence_audit"]["theorem_export"]
    roughness = t["sobolev_roughness_order_dependence_audit"]
    lines = [
        "# P2570/S1520 APD Sobolev roughness selector order-dependence audit", "",
        f"Status: `{payload['status']}`", "", "## Result", "",
        f"- Frontier atom under attack: `{t['frontier_atom_under_attack']}`.",
        f"- Candidate principle: `{t['candidate_principle_under_test']}`.",
        f"- P2569 interpolation family inherited: `{t['p2569_interpolation_family_inherited']}`.",
        f"- Derivative orders audited: `{roughness['derivative_orders']}`.",
        f"- Orders selecting nonzero lambda: `{roughness['orders_selecting_nonzero_lambda']}`.",
        f"- Orders selecting base lambda zero: `{roughness['orders_selecting_base_lambda_zero']}`.",
        f"- Roughness order changes selector: `{roughness['roughness_order_changes_selector']}`.",
        f"- Finite APD values select Sobolev order: `{not t['finite_apd_values_do_not_select_sobolev_order']}`.", "",
        "## Interpretation", "",
        "On the P2569 family `q_lambda=q_interp+lambda*prod_{d=0}^{11}(x-d)`, the audited Sobolev roughness objective has a closed-form quadratic in `lambda`.  Different derivative orders choose different minimizers: low-order roughness selects nonzero lambda values while the twelfth-derivative roughness selects the base interpolant.  Because every selected member still preserves the finite APD nodes, the derivative order/measure is an extra source obligation.", "",
        "## Recommended next honest step", "", t["recommended_next_honest_step"], "", "## Negative controls", "",
        "No APD roughness selector source, strict A/P/D dynamic source, bridge theorem, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, or ToE closure is exported.", "", "## Fingerprint", "",
        f"`{payload['apd_sobolev_roughness_selector_order_dependence_audit']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
`P2570/S1520` follows P2569 by testing a concrete APD regularity selector: minimize Sobolev roughness `J_k(lambda)=int_0^11 |d^k q_lambda/dx^k|^2 dx` on the same vanishing-polynomial interpolation family.  The audit finds that the selected `lambda` depends on the derivative order: lower orders choose nonzero family members while order `12` chooses the base interpolant.  Thus APD finite exactness plus an unsourced roughness slogan still does not export `strict_dynamical_source_for_A_P_D`.
""".strip()
    lag_section = """
`P2570/S1520` blocks promotion of APD minimum-roughness interpolation into role-bearing `L_total` dynamics.  A legitimate APD variational term must source the derivative order, measure, and boundary class; otherwise the Sobolev selector is metric-dependent bookkeeping rather than strict nadsoliton dynamics.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2570/S1520 APD Sobolev roughness selector order-dependence guard", "## P2570/S1520 APD Sobolev roughness selector order-dependence guard\n\n" + eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2570/S1520 APD Sobolev roughness selector order-dependence Ltotal guard", "## P2570/S1520 APD Sobolev roughness selector order-dependence Ltotal guard\n\n" + lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
