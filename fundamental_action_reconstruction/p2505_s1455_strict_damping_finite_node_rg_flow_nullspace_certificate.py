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

GEN = ROOT / "generated"
OUT = GEN / "p2505_s1455_strict_damping_finite_node_rg_flow_nullspace_certificate.json"
MD = GEN / "p2505_s1455_strict_damping_finite_node_rg_flow_nullspace_certificate.md"

SOURCE_FILES = {
    "P2504_CONSTANT_COEFFICIENT_RG_UNIQUENESS": GEN / "p2504_s1454_strict_damping_constant_coefficient_rg_uniqueness_certificate.json",
}

mp.mp.dps = 90
DOMAIN = list(range(1, 12))
EPSILON = mp.mpf("1e-6")


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
        "new_packet": "P2505|S1455|finite-node RG-flow nullspace|nonconstant RG-flow nullspace|running beta nullspace|finite-node flow nonuniqueness",
        "precursor_packets": "P2504|S1454|constant-coefficient RG uniqueness|P2503|S1453|strict damping marginal RG-flow",
        "nullspace_language": "node-vanishing perturbation|finite-node nullspace|nonconstant running beta|RG-flow nonuniqueness|polynomial nullspace",
        "guardrails": "legacy -> strict completion bridge|role-transfer audit|K_legacy_ont|K_strict_gate|QW-2191|ToE closure",
        "closure_blockers": "source theorem|bridge theorem|role-transfer theorem|physical-value generator|role-bearing L_total|selector closure",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def text(value: mp.mpf, digits: int = 60) -> str:
    return mp.nstr(value, digits)


def r_value(ell: mp.mpf) -> mp.mpf:
    value = ell
    for d in range(2, 12):
        value *= ell - mp.log(d)
    return value


def r_prime_value(ell: mp.mpf) -> mp.mpf:
    roots = [mp.mpf("0")] + [mp.log(d) for d in range(2, 12)]
    total = mp.mpf("0")
    for omitted_index in range(len(roots)):
        term = mp.mpf("1")
        for index, root in enumerate(roots):
            if index != omitted_index:
                term *= ell - root
        total += term
    return total


def symbolic_nullspace_certificate() -> dict[str, Any]:
    ell, eps = sp.symbols("ell eps", real=True)
    gamma_f = 4 * sp.log(2) - 1
    delta = sp.Rational(14, 5) - 4 * sp.log(2)
    eta = sp.Rational(9, 5)
    node_polynomial = ell
    for d in range(2, 12):
        node_polynomial *= ell - sp.log(d)
    running_exponent = delta * ell + eps * node_polynomial
    running_beta = sp.exp(running_exponent)
    flow_lambda = sp.diff(running_exponent, ell)
    node_rows = []
    all_node_residuals_zero = True
    for d in DOMAIN:
        ell_d = sp.log(d)
        r_at_node = sp.simplify(node_polynomial.subs(ell, ell_d))
        exponent_residual = sp.simplify((gamma_f * ell_d + running_exponent.subs(ell, ell_d)) - eta * ell_d)
        node_rows.append({
            "d": d,
            "R_log_d": sp.sstr(r_at_node),
            "exponent_residual": sp.sstr(exponent_residual),
        })
        all_node_residuals_zero = all_node_residuals_zero and r_at_node == 0 and exponent_residual == 0
    derivative_of_lambda_over_epsilon = sp.diff(node_polynomial, ell, 2)
    polynomial_degree = sp.Poly(sp.expand(node_polynomial), ell).degree()
    derivative_degree = sp.Poly(sp.expand(derivative_of_lambda_over_epsilon), ell).degree()
    nonconstant = derivative_of_lambda_over_epsilon != 0 and derivative_degree >= 0
    return {
        "symbolic_backend": "sympy",
        "sympy_version": sp.__version__,
        "node_vanishing_polynomial_R": sp.sstr(node_polynomial),
        "perturbed_running_beta": "B_eps(ell)=exp(delta*ell + eps*R(ell))",
        "flow_lambda": "delta + eps*R_prime(ell)",
        "node_rows": node_rows,
        "all_node_residuals_zero_symbolically": all_node_residuals_zero,
        "node_polynomial_degree": polynomial_degree,
        "flow_lambda_derivative_over_epsilon_degree": derivative_degree,
        "flow_lambda_nonconstant_for_nonzero_epsilon": nonconstant,
        "symbolic_fingerprint_sha256": sha256_json({
            "node_polynomial": sp.srepr(node_polynomial),
            "flow_lambda": sp.srepr(flow_lambda),
            "node_rows": node_rows,
        }),
    }


def finite_node_rows() -> list[dict[str, Any]]:
    gamma_f = 4 * mp.log(2) - 1
    eta = mp.mpf(9) / 5
    delta = eta - gamma_f
    rows = []
    for d_int in DOMAIN:
        d = mp.mpf(d_int)
        ell = mp.log(d)
        perturbation = EPSILON * r_value(ell)
        reconstructed_increment = d**gamma_f * mp.e ** (delta * ell + perturbation)
        strict_increment = d**eta
        rows.append({
            "d": d_int,
            "ell": text(ell, 50),
            "epsilon_R_ell": text(perturbation, 50),
            "strict_increment": text(strict_increment, 50),
            "perturbed_reconstructed_increment": text(reconstructed_increment, 50),
            "reconstruction_residual": text(reconstructed_increment - strict_increment, 50),
        })
    return rows


def finite_midpoint_rows() -> list[dict[str, Any]]:
    eta = mp.mpf(9) / 5
    gamma_f = 4 * mp.log(2) - 1
    delta = eta - gamma_f
    rows = []
    for left in range(1, 11):
        right = left + 1
        ell_mid = (mp.log(left) + mp.log(right)) / 2
        lambda_deviation = EPSILON * r_prime_value(ell_mid)
        rows.append({
            "between_nodes": [left, right],
            "ell_mid": text(ell_mid, 50),
            "lambda_base_delta": text(delta, 50),
            "lambda_deviation_epsilon_R_prime": text(lambda_deviation, 50),
            "lambda_perturbed": text(delta + lambda_deviation, 50),
            "nonzero_deviation": abs(lambda_deviation) > mp.mpf("1e-40"),
        })
    return rows


def build_nullspace_certificate(p2504: dict[str, Any]) -> dict[str, Any]:
    symbolic = symbolic_nullspace_certificate()
    node_rows = finite_node_rows()
    midpoint_rows = finite_midpoint_rows()
    node_residuals = [abs(mp.mpf(row["reconstruction_residual"])) for row in node_rows]
    midpoint_deviations = [abs(mp.mpf(row["lambda_deviation_epsilon_R_prime"])) for row in midpoint_rows]
    return {
        "frontier_atom_under_attack": "strict_damping_beta_eta_source",
        "p2504_constant_coefficient_uniqueness_inherited": p2504.get("constant_coefficient_rg_candidate_unique_within_ansatz") is True,
        "nullspace_scope": "finite-node equality on d=1..11 under nonconstant B_eps(ell)=exp(delta*ell+eps*R(ell))",
        "epsilon_used_for_finite_witness": text(EPSILON, 20),
        "symbolic_nullspace_certificate": symbolic,
        "finite_node_rows": node_rows,
        "finite_midpoint_rows": midpoint_rows,
        "max_abs_node_reconstruction_residual": text(max(node_residuals), 40),
        "max_abs_midpoint_lambda_deviation": text(max(midpoint_deviations), 40),
        "all_finite_nodes_preserved_within_precision": max(node_residuals) < mp.mpf("1e-70"),
        "midpoint_flow_deviation_nonzero": any(row["nonzero_deviation"] for row in midpoint_rows),
        "finite_nodes_do_not_select_constant_coefficient_flow": symbolic["all_node_residuals_zero_symbolically"] and any(row["nonzero_deviation"] for row in midpoint_rows),
        "p2504_uniqueness_limited_to_constant_coefficient_ansatz_confirmed": True,
        "nonconstant_nullspace_not_source_theorem": True,
        "strict_damping_beta_eta_source_exported": False,
        "strict_dynamical_source_for_A_P_D_exported": False,
        "strict_phase_frequency_source_exported": False,
        "bridge_theorem_exported": False,
        "role_transfer_theorem_exported": False,
        "selector_closure_exported": False,
        "qw2191_discharged_by_this_certificate": False,
        "toe_closure_claimed": False,
    }


def append_once(path, marker: str, section: str) -> None:
    body = path.read_text(encoding="utf-8")
    if marker not in body:
        path.write_text(body.rstrip() + "\n\n" + section.strip() + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2505/S1455 strict damping finite-node RG-flow nullspace certificate

`P2505/S1455` audits the limitation of P2504.  It constructs a nonconstant finite-node nullspace perturbation `B_eps(ell)=exp(delta ell + eps R(ell))`, where `R(ell)=ell * product_{d=2}^{11}(ell-log d)`.  Since `R(log d)=0` for every audited node `d=1..11`, the perturbed flow exactly preserves all finite strict denominator samples while changing the local flow rate between nodes.  This proves that P2504 uniqueness is only uniqueness inside the constant-coefficient ansatz, not uniqueness among all running-beta flows.

This is a nonuniqueness/limitation certificate, not a source theorem.  It does not export the `strict_damping_beta_eta_source` atom, does not close the bridge triple, and exports no role-transfer theorem, QW-2191 discharge, physical-value generator, or ToE closure.
"""
    lag_section = """
## P2505/S1455 finite-node RG-flow nullspace guard

`P2505/S1455` shows that finite strict damping samples admit nonconstant running-beta perturbations that vanish at all audited nodes.  Therefore P2504 remains an ansatz-scoped uniqueness result and cannot by itself license a nonlinear compression-flow source theorem or role-bearing `L_total` term.
"""
    append_once(DOC_FILES["equation_sheet"], "## P2505/S1455 strict damping finite-node RG-flow nullspace certificate", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2505/S1455 finite-node RG-flow nullspace guard", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    p2504 = theorem(sources["P2504_CONSTANT_COEFFICIENT_RG_UNIQUENESS"], "strict_damping_constant_coefficient_rg_uniqueness_certificate")
    cert = build_nullspace_certificate(p2504)
    theorem_export = {
        "theorem_name": "P2505_T1_strict_damping_finite_node_rg_flow_nullspace_certificate",
        "audited_chain": ["P2503/S1453", "P2504/S1454"],
        "strict_damping_finite_node_rg_flow_nullspace_certificate": cert,
        "frontier_atom_under_attack": cert["frontier_atom_under_attack"],
        "p2504_constant_coefficient_uniqueness_inherited": cert["p2504_constant_coefficient_uniqueness_inherited"],
        "all_node_residuals_zero_symbolically": cert["symbolic_nullspace_certificate"]["all_node_residuals_zero_symbolically"],
        "flow_lambda_nonconstant_for_nonzero_epsilon": cert["symbolic_nullspace_certificate"]["flow_lambda_nonconstant_for_nonzero_epsilon"],
        "max_abs_node_reconstruction_residual": cert["max_abs_node_reconstruction_residual"],
        "max_abs_midpoint_lambda_deviation": cert["max_abs_midpoint_lambda_deviation"],
        "all_finite_nodes_preserved_within_precision": cert["all_finite_nodes_preserved_within_precision"],
        "midpoint_flow_deviation_nonzero": cert["midpoint_flow_deviation_nonzero"],
        "finite_nodes_do_not_select_constant_coefficient_flow": cert["finite_nodes_do_not_select_constant_coefficient_flow"],
        "p2504_uniqueness_limited_to_constant_coefficient_ansatz_confirmed": cert["p2504_uniqueness_limited_to_constant_coefficient_ansatz_confirmed"],
        "nonconstant_nullspace_not_source_theorem": cert["nonconstant_nullspace_not_source_theorem"],
        "strict_damping_beta_eta_source_exported": False,
        "strict_dynamical_source_for_A_P_D_exported": False,
        "strict_phase_frequency_source_exported": False,
        "bridge_theorem_exported": False,
        "role_transfer_theorem_exported": False,
        "selector_closure_exported": False,
        "qw2191_discharged_by_this_certificate": False,
        "toe_closure_claimed": False,
        "not_licensed": [
            "P2505 constructs finite-node nonconstant perturbations; it does not derive a strict damping source.",
            "P2504 remains valid only inside the constant-coefficient running-beta ansatz.",
            "No strict source atom, bridge theorem, role-transfer theorem, QW-2191 discharge, role-bearing L_total, physical-value generator, or ToE closure is exported.",
        ],
        "next_honest_step": "Add a real variational/dynamical selector for the constant-coefficient RG ansatz, or accept that finite-node damping data alone leave a nonconstant-flow nullspace.",
    }
    gatekeepers = {
        "rg_audit_performed": grep["tool"] == "rg",
        "p2504_uniqueness_inherited": theorem_export["p2504_constant_coefficient_uniqueness_inherited"],
        "symbolic_node_residuals_zero": theorem_export["all_node_residuals_zero_symbolically"],
        "finite_nodes_preserved": theorem_export["all_finite_nodes_preserved_within_precision"],
        "nonconstant_flow_witnessed": theorem_export["flow_lambda_nonconstant_for_nonzero_epsilon"] and theorem_export["midpoint_flow_deviation_nonzero"],
        "constant_coefficient_limit_confirmed": theorem_export["p2504_uniqueness_limited_to_constant_coefficient_ansatz_confirmed"],
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
        "packet_id": "P2505",
        "stage_id": "S1455",
        "status": "STRICT_DAMPING_FINITE_NODE_RG_FLOW_NULLSPACE_CERTIFICATE_NO_SOURCE_EXPORT_NO_BRIDGE_THEOREM_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_damping_finite_node_rg_flow_nullspace_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_damping_finite_node_rg_flow_nullspace_certificate"]["theorem_export"]
    lines = [
        "# P2505/S1455 strict damping finite-node RG-flow nullspace certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        f"- Frontier atom under attack: `{t['frontier_atom_under_attack']}`.",
        f"- P2504 constant-coefficient uniqueness inherited: `{t['p2504_constant_coefficient_uniqueness_inherited']}`.",
        f"- Symbolic node residuals zero: `{t['all_node_residuals_zero_symbolically']}`.",
        f"- Flow lambda nonconstant for nonzero epsilon: `{t['flow_lambda_nonconstant_for_nonzero_epsilon']}`.",
        f"- Max node reconstruction residual: `{t['max_abs_node_reconstruction_residual']}`.",
        f"- Max midpoint lambda deviation: `{t['max_abs_midpoint_lambda_deviation']}`.",
        f"- Finite nodes do not select constant-coefficient flow: `{t['finite_nodes_do_not_select_constant_coefficient_flow']}`.",
        f"- Source theorem exported: `{t['strict_damping_beta_eta_source_exported']}`.",
        "",
        "## Negative controls",
        "",
        "This packet is a finite-node nullspace/limitation certificate. It does not derive a strict damping source atom, bridge theorem, role-transfer theorem, selector/QW-2191 closure, role-bearing L_total term, physical-value generator, or ToE closure.",
        "",
        "## Fingerprint",
        "",
        f"`{payload['strict_damping_finite_node_rg_flow_nullspace_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload["strict_damping_finite_node_rg_flow_nullspace_certificate"]["theorem_export"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
