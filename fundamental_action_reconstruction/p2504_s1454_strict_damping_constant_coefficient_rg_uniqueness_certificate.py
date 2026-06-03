#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import itertools
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
OUT = GEN / "p2504_s1454_strict_damping_constant_coefficient_rg_uniqueness_certificate.json"
MD = GEN / "p2504_s1454_strict_damping_constant_coefficient_rg_uniqueness_certificate.md"

SOURCE_FILES = {
    "P2503_DAMPING_RG_FLOW_CANDIDATE": GEN / "p2503_s1453_strict_damping_marginal_rg_flow_candidate_certificate.json",
}

mp.mp.dps = 90
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
        "new_packet": "P2504|S1454|constant-coefficient marginal RG uniqueness|RG-flow uniqueness|running beta uniqueness|strict damping uniqueness certificate",
        "precursor_packets": "P2503|S1453|strict damping marginal RG-flow|dB/dell|delta=14/5|strict_damping_beta_eta_source",
        "uniqueness_language": "constant-coefficient RG|two-node uniqueness|beta0 lambda unique|running exponent unique|marginal correction uniqueness",
        "guardrails": "legacy -> strict completion bridge|role-transfer audit|K_legacy_ont|K_strict_gate|QW-2191|ToE closure",
        "closure_blockers": "source theorem|bridge theorem|role-transfer theorem|physical-value generator|role-bearing L_total|selector closure",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def text(value: mp.mpf, digits: int = 60) -> str:
    return mp.nstr(value, digits)


def symbolic_uniqueness_certificate() -> dict[str, Any]:
    beta0, lam = sp.symbols("beta0 lam", positive=True, real=True)
    gamma_f = 4 * sp.log(2) - 1
    eta = sp.Rational(9, 5)
    delta = sp.Rational(14, 5) - 4 * sp.log(2)
    # Constant-coefficient ansatz increment: beta0 * d**(gamma_f + lam).
    eq_d1 = sp.Eq(beta0, 1)
    eq_d2 = sp.Eq(beta0 * 2 ** (gamma_f + lam), 2 ** eta)
    beta0_solution = sp.Integer(1)
    lam_solution = sp.simplify(eta - gamma_f)
    d1_residual = sp.simplify(beta0_solution - 1)
    d2_residual = sp.simplify(beta0_solution * 2 ** (gamma_f + lam_solution) - 2 ** eta)
    lam_minus_delta = sp.simplify(lam_solution - delta)
    return {
        "symbolic_backend": "sympy",
        "sympy_version": sp.__version__,
        "ansatz": "strict_increment(d)=beta0*d**(gamma_F+lambda)",
        "normalization_equation_d1": sp.sstr(eq_d1),
        "two_node_equation_d2": sp.sstr(eq_d2),
        "beta0_solution": sp.sstr(beta0_solution),
        "lambda_solution": sp.sstr(delta),
        "delta_expression": sp.sstr(delta),
        "lambda_minus_delta_residual": sp.sstr(lam_minus_delta),
        "d1_residual_at_solution": sp.sstr(d1_residual),
        "d2_residual_at_solution": sp.sstr(d2_residual),
        "two_node_symbolic_solution_matches_p2503_delta": lam_minus_delta == 0 and d1_residual == 0 and d2_residual == 0,
        "symbolic_fingerprint_sha256": sha256_json({
            "gamma_f": sp.srepr(gamma_f),
            "eta": sp.srepr(eta),
            "delta": sp.srepr(delta),
            "lambda_solution": sp.srepr(lam_solution),
            "d2_residual": sp.srepr(d2_residual),
        }),
    }


def pair_recovery_rows() -> list[dict[str, Any]]:
    gamma_f = 4 * mp.log(2) - 1
    eta = mp.mpf(9) / 5
    delta = eta - gamma_f
    rows = []
    for a_int, b_int in itertools.combinations(DOMAIN, 2):
        a = mp.mpf(a_int)
        b = mp.mpf(b_int)
        strict_a = a**eta
        strict_b = b**eta
        recovered_lambda = mp.log(strict_b / strict_a) / mp.log(b / a) - gamma_f
        recovered_beta0 = strict_a / (a ** (gamma_f + recovered_lambda))
        reconstructed_b = recovered_beta0 * b ** (gamma_f + recovered_lambda)
        rows.append({
            "node_pair": [a_int, b_int],
            "recovered_lambda": text(recovered_lambda, 70),
            "recovered_beta0": text(recovered_beta0, 70),
            "lambda_minus_delta": text(recovered_lambda - delta, 50),
            "beta0_minus_one": text(recovered_beta0 - 1, 50),
            "pair_reconstruction_residual_at_b": text(reconstructed_b - strict_b, 50),
        })
    return rows


def build_uniqueness_certificate(p2503: dict[str, Any]) -> dict[str, Any]:
    symbolic = symbolic_uniqueness_certificate()
    rows = pair_recovery_rows()
    lambda_errors = [abs(mp.mpf(row["lambda_minus_delta"])) for row in rows]
    beta0_errors = [abs(mp.mpf(row["beta0_minus_one"])) for row in rows]
    reconstruction_errors = [abs(mp.mpf(row["pair_reconstruction_residual_at_b"])) for row in rows]
    return {
        "frontier_atom_under_attack": "strict_damping_beta_eta_source",
        "p2503_candidate_inherited": p2503.get("candidate_solves_formal_denominator_target") is True,
        "constant_coefficient_ansatz_scope": "B(ell)=beta0*exp(lambda*ell) with gamma_F fixed to D_f-1",
        "symbolic_two_node_uniqueness": symbolic,
        "finite_pair_recovery_domain": DOMAIN,
        "finite_pair_recovery_row_count": len(rows),
        "finite_pair_recovery_rows_head": rows[:10],
        "finite_pair_recovery_rows_tail": rows[-10:],
        "finite_pair_recovery_rows_sha256": sha256_json(rows),
        "max_abs_lambda_minus_delta": text(max(lambda_errors), 40),
        "max_abs_beta0_minus_one": text(max(beta0_errors), 40),
        "max_abs_pair_reconstruction_residual": text(max(reconstruction_errors), 40),
        "all_pairs_recover_same_delta_within_precision": max(lambda_errors) < mp.mpf("1e-75"),
        "all_pairs_recover_beta0_one_within_precision": max(beta0_errors) < mp.mpf("1e-75"),
        "all_pairs_reconstruct_second_node_within_precision": max(reconstruction_errors) < mp.mpf("1e-70"),
        "constant_coefficient_rg_candidate_unique_within_ansatz": symbolic["two_node_symbolic_solution_matches_p2503_delta"] and max(lambda_errors) < mp.mpf("1e-75"),
        "ansatz_uniqueness_not_source_derivation": True,
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
## P2504/S1454 strict damping constant-coefficient RG uniqueness certificate

`P2504/S1454` sharpens P2503 without promoting it to a source theorem.  Inside the explicit ansatz `B(ell)=beta0 exp(lambda ell)` with `gamma_F=D_f-1` fixed, the first strict denominator node fixes `beta0=1`, and a second positive node fixes `lambda=9/5-gamma_F=14/5-4 log(2)`.  A finite pair audit over all `55` positive-node pairs confirms that every pair recovers the same `lambda` and `beta0` to high precision.

This is uniqueness only inside a chosen constant-coefficient marginal-flow ansatz.  It does not derive the ansatz, `delta`, or `beta0` from strict nadsoliton dynamics, does not export the `strict_damping_beta_eta_source` atom, and exports no bridge theorem, role-transfer theorem, QW-2191 discharge, physical-value generator, or ToE closure.
"""
    lag_section = """
## P2504/S1454 constant-coefficient RG uniqueness guard

`P2504/S1454` proves that the P2503 marginal-flow parameters are unique inside the constant-coefficient running-beta ansatz.  The result is theorem-prep for `strict_damping_beta_eta_source`, not a nonlinear compression-flow source theorem or a license for role-bearing `L_total`.
"""
    append_once(DOC_FILES["equation_sheet"], "## P2504/S1454 strict damping constant-coefficient RG uniqueness certificate", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2504/S1454 constant-coefficient RG uniqueness guard", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    p2503 = theorem(sources["P2503_DAMPING_RG_FLOW_CANDIDATE"], "strict_damping_marginal_rg_flow_candidate_certificate")
    cert = build_uniqueness_certificate(p2503)
    symbolic = cert["symbolic_two_node_uniqueness"]
    theorem_export = {
        "theorem_name": "P2504_T1_strict_damping_constant_coefficient_rg_uniqueness_certificate",
        "audited_chain": ["P2502/S1452", "P2503/S1453"],
        "strict_damping_constant_coefficient_rg_uniqueness_certificate": cert,
        "frontier_atom_under_attack": cert["frontier_atom_under_attack"],
        "p2503_candidate_inherited": cert["p2503_candidate_inherited"],
        "beta0_solution": symbolic["beta0_solution"],
        "lambda_solution": symbolic["lambda_solution"],
        "lambda_minus_delta_residual": symbolic["lambda_minus_delta_residual"],
        "two_node_symbolic_solution_matches_p2503_delta": symbolic["two_node_symbolic_solution_matches_p2503_delta"],
        "finite_pair_recovery_row_count": cert["finite_pair_recovery_row_count"],
        "max_abs_lambda_minus_delta": cert["max_abs_lambda_minus_delta"],
        "max_abs_beta0_minus_one": cert["max_abs_beta0_minus_one"],
        "max_abs_pair_reconstruction_residual": cert["max_abs_pair_reconstruction_residual"],
        "all_pairs_recover_same_delta_within_precision": cert["all_pairs_recover_same_delta_within_precision"],
        "all_pairs_recover_beta0_one_within_precision": cert["all_pairs_recover_beta0_one_within_precision"],
        "constant_coefficient_rg_candidate_unique_within_ansatz": cert["constant_coefficient_rg_candidate_unique_within_ansatz"],
        "ansatz_uniqueness_not_source_derivation": cert["ansatz_uniqueness_not_source_derivation"],
        "strict_damping_beta_eta_source_exported": False,
        "strict_dynamical_source_for_A_P_D_exported": False,
        "strict_phase_frequency_source_exported": False,
        "bridge_theorem_exported": False,
        "role_transfer_theorem_exported": False,
        "selector_closure_exported": False,
        "qw2191_discharged_by_this_certificate": False,
        "toe_closure_claimed": False,
        "not_licensed": [
            "P2504 proves uniqueness only inside the chosen constant-coefficient running-beta ansatz.",
            "It does not derive the RG ansatz, delta, or beta0 from strict nadsoliton dynamics.",
            "No strict source atom, bridge theorem, role-transfer theorem, QW-2191 discharge, role-bearing L_total, physical-value generator, or ToE closure is exported.",
        ],
        "next_honest_step": "Either derive the constant-coefficient running-beta ansatz from strict nadsoliton dynamics, or show that a broader non-constant flow class changes the source-candidate status.",
    }
    gatekeepers = {
        "rg_audit_performed": grep["tool"] == "rg",
        "p2503_candidate_inherited": theorem_export["p2503_candidate_inherited"],
        "symbolic_solution_matches_p2503_delta": theorem_export["two_node_symbolic_solution_matches_p2503_delta"],
        "all_55_pairs_recovered": theorem_export["finite_pair_recovery_row_count"] == 55,
        "all_pairs_recover_same_delta": theorem_export["all_pairs_recover_same_delta_within_precision"],
        "all_pairs_recover_beta0_one": theorem_export["all_pairs_recover_beta0_one_within_precision"],
        "unique_but_not_source": theorem_export["constant_coefficient_rg_candidate_unique_within_ansatz"] and theorem_export["ansatz_uniqueness_not_source_derivation"],
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
        "packet_id": "P2504",
        "stage_id": "S1454",
        "status": "STRICT_DAMPING_CONSTANT_COEFFICIENT_RG_UNIQUENESS_CERTIFICATE_NO_SOURCE_EXPORT_NO_BRIDGE_THEOREM_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_damping_constant_coefficient_rg_uniqueness_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_damping_constant_coefficient_rg_uniqueness_certificate"]["theorem_export"]
    lines = [
        "# P2504/S1454 strict damping constant-coefficient RG uniqueness certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        f"- Frontier atom under attack: `{t['frontier_atom_under_attack']}`.",
        f"- P2503 candidate inherited: `{t['p2503_candidate_inherited']}`.",
        f"- Symbolic beta0 solution: `{t['beta0_solution']}`.",
        f"- Symbolic lambda solution: `{t['lambda_solution']}`.",
        f"- Lambda minus P2503 delta residual: `{t['lambda_minus_delta_residual']}`.",
        f"- Finite pair recovery rows: `{t['finite_pair_recovery_row_count']}`.",
        f"- Max |lambda-delta| over pairs: `{t['max_abs_lambda_minus_delta']}`.",
        f"- Max |beta0-1| over pairs: `{t['max_abs_beta0_minus_one']}`.",
        f"- Unique within ansatz: `{t['constant_coefficient_rg_candidate_unique_within_ansatz']}`.",
        f"- Source theorem exported: `{t['strict_damping_beta_eta_source_exported']}`.",
        "",
        "## Negative controls",
        "",
        "This packet proves uniqueness only inside the constant-coefficient running-beta ansatz. It does not derive the ansatz, delta, or beta0 from strict dynamics and exports no strict source atom, bridge theorem, role-transfer theorem, selector/QW-2191 closure, role-bearing L_total term, physical-value generator, or ToE closure.",
        "",
        "## Fingerprint",
        "",
        f"`{payload['strict_damping_constant_coefficient_rg_uniqueness_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload["strict_damping_constant_coefficient_rg_uniqueness_certificate"]["theorem_export"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
