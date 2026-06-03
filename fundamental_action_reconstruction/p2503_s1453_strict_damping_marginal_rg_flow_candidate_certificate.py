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
SCRATCH = ROOT / "scratch"
OUT = GEN / "p2503_s1453_strict_damping_marginal_rg_flow_candidate_certificate.json"
MD = GEN / "p2503_s1453_strict_damping_marginal_rg_flow_candidate_certificate.md"

SOURCE_FILES = {
    "P2502_BRIDGE_TRIPLE_FRONTIER": GEN / "p2502_s1452_strict_completion_bridge_minimal_triple_frontier_certificate.json",
    "SCRATCH_DF_LOG_CORRECTION": SCRATCH / "bridge_df_log_correction_report.json",
    "SCRATCH_DF_MARGINAL_RG_FLOW": SCRATCH / "bridge_df_marginal_rg_flow_report.json",
    "SCRATCH_DAMPING_PARAMETER_IDENTIFIABILITY": SCRATCH / "bridge_strict_completion_legacy_to_strict_damping_parameter_identifiability_certificate_report.json",
}

mp.mp.dps = 80
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
        "new_packet": "P2503|S1453|strict damping marginal RG-flow|RG-flow candidate certificate|running beta source candidate|damping source candidate",
        "precursor_packets": "P2502|S1452|strict-source triple|strict_damping_beta_eta_source|bridge minimal triple frontier|P2501",
        "rg_flow_language": "marginal RG-flow|running beta|D_f marginal|delta = 14/5|14/5 - 4 log|dB/dell|RG law",
        "damping_identifiability": "beta eta identifiability|fifth-power cover|strict denominator|damping parameter identifiability|eta=9/5",
        "guardrails": "legacy -> strict completion bridge|role-transfer audit|K_legacy_ont|K_strict_gate|QW-2191|ToE closure",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def text(value: mp.mpf, digits: int = 50) -> str:
    return mp.nstr(value, digits)


def symbolic_flow_certificate() -> dict[str, Any]:
    ell, beta0 = sp.symbols("ell beta0", positive=True)
    gamma_f = 4 * sp.log(2) - 1
    delta = sp.Rational(14, 5) - 4 * sp.log(2)
    eta = sp.Rational(9, 5)
    beta_running = beta0 * sp.exp(delta * ell)
    denominator_after_flow = 1 + beta_running * sp.exp(gamma_f * ell)
    strict_denominator = 1 + beta0 * sp.exp(eta * ell)
    ode_residual = sp.simplify(sp.diff(beta_running, ell) - delta * beta_running)
    denominator_residual = sp.simplify(denominator_after_flow - strict_denominator)
    exponent_residual = sp.simplify(gamma_f + delta - eta)
    return {
        "symbolic_backend": "sympy",
        "sympy_version": sp.__version__,
        "gamma_f_expression": sp.sstr(gamma_f),
        "delta_expression": sp.sstr(delta),
        "eta_expression": sp.sstr(eta),
        "gamma_f_plus_delta_minus_eta_residual": sp.sstr(exponent_residual),
        "running_beta_solution": sp.sstr(beta_running),
        "rg_ode_residual": sp.sstr(ode_residual),
        "denominator_after_flow": sp.sstr(denominator_after_flow),
        "strict_denominator": sp.sstr(strict_denominator),
        "denominator_residual": sp.sstr(denominator_residual),
        "all_symbolic_residuals_zero": exponent_residual == 0 and ode_residual == 0 and denominator_residual == 0,
        "symbolic_fingerprint_sha256": sha256_json({
            "gamma_f": sp.srepr(gamma_f),
            "delta": sp.srepr(delta),
            "eta": sp.srepr(eta),
            "ode_residual": sp.srepr(ode_residual),
            "denominator_residual": sp.srepr(denominator_residual),
        }),
    }


def finite_delta_omission_rows() -> list[dict[str, Any]]:
    gamma_f = 4 * mp.log(2) - 1
    eta = mp.mpf(9) / 5
    delta = eta - gamma_f
    rows = []
    for d_int in DOMAIN:
        d = mp.mpf(d_int)
        strict_increment = d**eta
        no_flow_increment = d**gamma_f
        flow_factor = mp.e ** (delta * mp.log(d))
        reconstructed_increment = no_flow_increment * flow_factor
        rows.append({
            "d": d_int,
            "strict_increment_d_eta": text(strict_increment, 40),
            "no_flow_increment_d_gamma_f": text(no_flow_increment, 40),
            "running_factor_exp_delta_log_d": text(flow_factor, 40),
            "reconstructed_increment": text(reconstructed_increment, 40),
            "delta_omitted_residual_strict_minus_no_flow": text(strict_increment - no_flow_increment, 40),
            "flow_reconstruction_residual": text(reconstructed_increment - strict_increment, 40),
            "delta_omitted_residual_positive": strict_increment - no_flow_increment >= 0,
        })
    return rows


def build_rg_flow_candidate(sources: dict[str, dict[str, Any]]) -> dict[str, Any]:
    p2502 = theorem(sources["P2502_BRIDGE_TRIPLE_FRONTIER"], "strict_completion_bridge_minimal_triple_frontier_certificate")
    rg_flow = sources["SCRATCH_DF_MARGINAL_RG_FLOW"]
    ident = sources["SCRATCH_DAMPING_PARAMETER_IDENTIFIABILITY"]
    symbolic = symbolic_flow_certificate()
    rows = finite_delta_omission_rows()
    omitted_residuals = [mp.mpf(row["delta_omitted_residual_strict_minus_no_flow"]) for row in rows]
    flow_residuals = [abs(mp.mpf(row["flow_reconstruction_residual"])) for row in rows]
    strict_triple = p2502.get("exact_strict_source_triple", [])
    return {
        "frontier_atom_under_attack": "strict_damping_beta_eta_source",
        "p2502_strict_source_triple_inherited": strict_triple,
        "frontier_atom_is_in_p2502_bridge_triple": "strict_damping_beta_eta_source" in strict_triple,
        "symbolic_rg_flow_candidate": symbolic,
        "finite_delta_omission_rows": rows,
        "finite_domain": DOMAIN,
        "max_delta_omitted_residual_on_domain": text(max(omitted_residuals), 40),
        "all_delta_omitted_residuals_nonnegative": all(row["delta_omitted_residual_positive"] for row in rows),
        "max_flow_reconstruction_abs_residual_on_domain": text(max(flow_residuals), 30),
        "scratch_rg_flow_exact_residual_inherited": rg_flow.get("aggregate_summary", {}).get("exact_rg_matches_strict_to_numeric_precision") is True,
        "scratch_damping_identifiability_inherited": ident.get("damping_parameter_identifiability_summary", {}).get("candidate_grid_unique_match") is True,
        "beta_eta_source_candidate_factorization": {
            "base_exponent_visible": "gamma_F = D_f - 1 = 4*log(2)-1",
            "required_marginal_correction": "delta = 14/5 - 4*log(2)",
            "target_exponent": "eta = gamma_F + delta = 9/5",
            "running_beta_law": "dB/dell = delta*B, B(ell)=beta0*exp(delta*ell)",
            "strict_denominator_after_candidate_flow": "1 + beta0*exp((gamma_F+delta)*ell)",
            "remaining_source_obligation": "derive delta and beta0 from nadsoliton/strict dynamics rather than postulating the RG-flow candidate",
        },
        "candidate_solves_formal_denominator_target": symbolic["all_symbolic_residuals_zero"] and max(flow_residuals) < mp.mpf("1e-70"),
        "candidate_is_not_derived_source_theorem": True,
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
## P2503/S1453 strict damping marginal RG-flow candidate certificate

`P2503/S1453` attacks the `strict_damping_beta_eta_source` atom from the P2502 bridge triple by auditing an exact candidate marginal flow.  With `gamma_F=D_f-1=4 log(2)-1` and `delta=14/5-4 log(2)`, the symbolic law `dB/dell=delta*B` gives `B(ell)=beta0 exp(delta ell)` and exactly reconstructs the strict denominator exponent: `gamma_F+delta=9/5`.  The symbolic ODE residual and denominator residual are zero, and finite rows show that omitting `delta` leaves a positive strict-minus-base damping residual on the audited positive nodes.

This is a source-candidate/factorization certificate only.  It does not derive `delta` or `beta0` from nadsoliton dynamics, does not export the `strict_damping_beta_eta_source` atom, does not close the P2502 bridge triple, and exports no role-transfer theorem, QW-2191 discharge, physical-value generator, or ToE closure.
"""
    lag_section = """
## P2503/S1453 strict damping marginal RG-flow candidate guard

`P2503/S1453` records a symbolic candidate flow that would supply the strict damping exponent if a future theorem derives `delta=14/5-4 log(2)`.  It is useful theorem-prep for the P2502 strict-source triple, but it is not yet a nonlinear compression-flow source theorem and does not license a role-bearing `L_total` term.
"""
    append_once(DOC_FILES["equation_sheet"], "## P2503/S1453 strict damping marginal RG-flow candidate certificate", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2503/S1453 strict damping marginal RG-flow candidate guard", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    cert = build_rg_flow_candidate(sources)
    symbolic = cert["symbolic_rg_flow_candidate"]
    theorem_export = {
        "theorem_name": "P2503_T1_strict_damping_marginal_rg_flow_candidate_certificate",
        "audited_chain": ["P2502/S1452", "scratch D_f marginal RG flow", "damping beta/eta identifiability"],
        "strict_damping_marginal_rg_flow_candidate_certificate": cert,
        "frontier_atom_under_attack": cert["frontier_atom_under_attack"],
        "frontier_atom_is_in_p2502_bridge_triple": cert["frontier_atom_is_in_p2502_bridge_triple"],
        "gamma_f_plus_delta_minus_eta_residual": symbolic["gamma_f_plus_delta_minus_eta_residual"],
        "rg_ode_residual": symbolic["rg_ode_residual"],
        "denominator_residual": symbolic["denominator_residual"],
        "all_symbolic_residuals_zero": symbolic["all_symbolic_residuals_zero"],
        "max_delta_omitted_residual_on_domain": cert["max_delta_omitted_residual_on_domain"],
        "all_delta_omitted_residuals_nonnegative": cert["all_delta_omitted_residuals_nonnegative"],
        "max_flow_reconstruction_abs_residual_on_domain": cert["max_flow_reconstruction_abs_residual_on_domain"],
        "candidate_solves_formal_denominator_target": cert["candidate_solves_formal_denominator_target"],
        "candidate_is_not_derived_source_theorem": cert["candidate_is_not_derived_source_theorem"],
        "strict_damping_beta_eta_source_exported": False,
        "strict_dynamical_source_for_A_P_D_exported": False,
        "strict_phase_frequency_source_exported": False,
        "bridge_theorem_exported": False,
        "role_transfer_theorem_exported": False,
        "selector_closure_exported": False,
        "qw2191_discharged_by_this_certificate": False,
        "toe_closure_claimed": False,
        "not_licensed": [
            "P2503 records an exact candidate marginal RG flow; it does not derive the flow from strict nadsoliton dynamics.",
            "The strict_damping_beta_eta_source atom remains open until delta and beta0 have a real source theorem.",
            "No bridge theorem, role-transfer theorem, QW-2191 discharge, role-bearing L_total, physical-value generator, or ToE closure is exported.",
        ],
        "next_honest_step": "Try to derive delta=14/5-4*log(2) and beta0=1 from a strict nadsoliton source, or prove that this marginal RG-flow candidate is only a bookkeeping ansatz.",
    }
    gatekeepers = {
        "rg_audit_performed": grep["tool"] == "rg",
        "frontier_atom_is_in_p2502_bridge_triple": theorem_export["frontier_atom_is_in_p2502_bridge_triple"],
        "symbolic_residuals_zero": theorem_export["all_symbolic_residuals_zero"],
        "delta_omission_audited_nonnegative": theorem_export["all_delta_omitted_residuals_nonnegative"],
        "candidate_solves_formal_denominator_target": theorem_export["candidate_solves_formal_denominator_target"],
        "candidate_marked_not_derived_source": theorem_export["candidate_is_not_derived_source_theorem"],
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
        "packet_id": "P2503",
        "stage_id": "S1453",
        "status": "STRICT_DAMPING_MARGINAL_RG_FLOW_CANDIDATE_CERTIFICATE_NO_SOURCE_EXPORT_NO_BRIDGE_THEOREM_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_damping_marginal_rg_flow_candidate_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_damping_marginal_rg_flow_candidate_certificate"]["theorem_export"]
    lines = [
        "# P2503/S1453 strict damping marginal RG-flow candidate certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        f"- Frontier atom under attack: `{t['frontier_atom_under_attack']}`.",
        f"- Atom is in P2502 bridge triple: `{t['frontier_atom_is_in_p2502_bridge_triple']}`.",
        f"- `gamma_F + delta - eta` residual: `{t['gamma_f_plus_delta_minus_eta_residual']}`.",
        f"- RG ODE residual: `{t['rg_ode_residual']}`.",
        f"- Denominator residual: `{t['denominator_residual']}`.",
        f"- Candidate solves formal denominator target: `{t['candidate_solves_formal_denominator_target']}`.",
        f"- Max residual if delta is omitted on audited domain: `{t['max_delta_omitted_residual_on_domain']}`.",
        f"- Source theorem exported: `{t['strict_damping_beta_eta_source_exported']}`.",
        "",
        "## Negative controls",
        "",
        "This packet records an exact candidate flow but does not derive the strict damping source atom, bridge theorem, role-transfer theorem, selector/QW-2191 closure, role-bearing L_total term, physical-value generator, or ToE closure.",
        "",
        "## Fingerprint",
        "",
        f"`{payload['strict_damping_marginal_rg_flow_candidate_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload["strict_damping_marginal_rg_flow_candidate_certificate"]["theorem_export"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
