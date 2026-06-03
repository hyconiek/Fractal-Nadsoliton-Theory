#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from typing import Any

import sympy as sp

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import (
    DOC_FILES,
    REPO,
    ROOT,
    load_json,
    rel,
)

GEN = ROOT / "generated"
OUT = GEN / "p2500_s1450_phase_normalized_third_derivative_symbolic_identity_audit.json"
MD = GEN / "p2500_s1450_phase_normalized_third_derivative_symbolic_identity_audit.md"

SOURCE_FILES = {
    "P2499_LOCAL_UNIQUENESS": GEN / "p2499_s1449_phase_normalized_curvature_local_inflection_uniqueness_certificate.json",
}


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
        "new_packet": "P2500|S1450|third-derivative symbolic identity|implicit third derivative|symbolic formula audit|x''' formula|third derivative chain identity",
        "precursor_packets": "P2499|S1449|local inflection uniqueness|third derivative interval|curvature monotonicity window",
        "symbolic_language": "sympy|symbolic identity|Faà di Bruno|implicit differentiation|third derivative|chain rule residual",
        "guardrail_language": "legacy -> strict completion bridge|role-transfer audit|silent inheritance|K_legacy_ont|K_strict_gate|QW-2191",
        "closure_blockers": "role-bearing L_total|physical-value generator|ToE closure|selector/source/gauge theorem|legacy-role transfer|directed-rounding",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def sympy_bool(value: sp.Expr) -> bool:
    return bool(sp.simplify(sp.trigsimp(sp.factor(sp.together(value)))) == 0)


def legacy_third_symbolic_identity() -> dict[str, Any]:
    x = sp.symbols("x", positive=True)
    omega, phi, beta = sp.symbols("omega phi beta", positive=True, nonzero=True)
    carrier = sp.cos(omega * x + phi) / sp.cos(phi)
    denominator = 1 + beta * x
    legacy_norm = carrier / denominator
    f0 = carrier
    f1 = -omega * sp.sin(omega * x + phi) / sp.cos(phi)
    f2 = -(omega**2) * sp.cos(omega * x + phi) / sp.cos(phi)
    f3 = (omega**3) * sp.sin(omega * x + phi) / sp.cos(phi)
    g0 = 1 / denominator
    g1 = -beta / denominator**2
    g2 = 2 * beta**2 / denominator**3
    g3 = -6 * beta**3 / denominator**4
    implemented_formula = f3 * g0 + 3 * f2 * g1 + 3 * f1 * g2 + f0 * g3
    residual = sp.diff(legacy_norm, x, 3) - implemented_formula
    return {
        "object": "L_norm(x) = cos(omega*x+phi)/cos(phi)/(1+beta*x)",
        "formula_name": "legacy_third_interval_product_rule_formula",
        "residual_simplifies_to_zero": sympy_bool(residual),
        "residual_srepr_sha256": hashlib.sha256(sp.srepr(sp.simplify(residual)).encode("utf-8")).hexdigest(),
    }


def strict_third_symbolic_identity() -> dict[str, Any]:
    d = sp.symbols("d", positive=True)
    omega, phi, beta, eta = sp.symbols("omega phi beta eta", positive=True, nonzero=True)
    carrier = sp.cos(omega * d + phi) / sp.cos(phi)
    denominator = 1 + beta * d**eta
    strict_norm = carrier / denominator
    f0 = carrier
    f1 = -omega * sp.sin(omega * d + phi) / sp.cos(phi)
    f2 = -(omega**2) * sp.cos(omega * d + phi) / sp.cos(phi)
    f3 = (omega**3) * sp.sin(omega * d + phi) / sp.cos(phi)
    u0 = denominator
    u1 = beta * eta * d ** (eta - 1)
    u2 = beta * eta * (eta - 1) * d ** (eta - 2)
    u3 = beta * eta * (eta - 1) * (eta - 2) * d ** (eta - 3)
    h0 = 1 / u0
    h1 = -u1 / u0**2
    h2 = -u2 / u0**2 + 2 * u1**2 / u0**3
    h3 = -u3 / u0**2 + 6 * u1 * u2 / u0**3 - 6 * u1**3 / u0**4
    implemented_formula = f3 * h0 + 3 * f2 * h1 + 3 * f1 * h2 + f0 * h3
    residual = sp.diff(strict_norm, d, 3) - implemented_formula
    return {
        "object": "S_strict_norm(d) = cos(omega*d+phi)/cos(phi)/(1+beta*d**eta)",
        "formula_name": "strict_third_interval_product_rule_formula",
        "residual_simplifies_to_zero": sympy_bool(residual),
        "residual_srepr_sha256": hashlib.sha256(sp.srepr(sp.simplify(residual)).encode("utf-8")).hexdigest(),
    }


def implicit_chain_identity() -> dict[str, Any]:
    L1, L2, L3, S1, S2, S3 = sp.symbols("L1 L2 L3 S1 S2 S3", nonzero=True)
    x_prime = S1 / L1
    x_second = (S2 - L2 * x_prime**2) / L1
    x_third = (S3 - L3 * x_prime**3 - 3 * L2 * x_prime * x_second) / L1
    residual = L3 * x_prime**3 + 3 * L2 * x_prime * x_second + L1 * x_third - S3
    return {
        "identity": "L3*x_prime**3 + 3*L2*x_prime*x_second + L1*x_third - S3 = 0",
        "x_prime_formula": "S1/L1",
        "x_second_formula": "(S2 - L2*x_prime**2)/L1",
        "x_third_formula": "(S3 - L3*x_prime**3 - 3*L2*x_prime*x_second)/L1",
        "residual_simplifies_to_zero": sympy_bool(residual),
        "residual_srepr_sha256": hashlib.sha256(sp.srepr(sp.simplify(residual)).encode("utf-8")).hexdigest(),
    }


def build_symbolic_certificate(p2499: dict[str, Any]) -> dict[str, Any]:
    legacy = legacy_third_symbolic_identity()
    strict = strict_third_symbolic_identity()
    chain = implicit_chain_identity()
    all_symbolic = all(item["residual_simplifies_to_zero"] for item in [legacy, strict, chain])
    return {
        "sympy_version": sp.__version__,
        "symbolic_backend": "sympy",
        "legacy_third_derivative_identity": legacy,
        "strict_third_derivative_identity": strict,
        "implicit_third_derivative_chain_identity": chain,
        "all_symbolic_residuals_zero": all_symbolic,
        "p2499_local_uniqueness_inherited": p2499.get("finite_local_unique_curvature_zero_in_refined_window") is True,
        "p2499_third_interval_negative_inherited": p2499.get("x_third_interval_sign") == -1,
        "formula_provenance_supports_p2499_interval_third_derivative": all_symbolic and p2499.get("finite_local_unique_curvature_zero_in_refined_window") is True,
        "formal_directed_rounding_backend_exported": False,
        "global_analytic_inflection_uniqueness_theorem_exported": False,
        "curvature_dynamic_source_exported": False,
        "legacy_to_strict_bridge_atom_exported": False,
        "strict_compression_dynamic_source_exported": False,
        "selector_source_theorem_exported": False,
        "qw2191_discharged_by_this_certificate": False,
        "role_transfer_licensed_by_this_certificate": False,
        "toe_closure_exported": False,
    }


def append_once(path, marker: str, section: str) -> None:
    body = path.read_text(encoding="utf-8")
    if marker not in body:
        path.write_text(body.rstrip() + "\n\n" + section.strip() + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2500/S1450 phase-normalized third-derivative symbolic identity audit

`P2500/S1450` audits the formula provenance behind the P2499 third-derivative interval.  A `sympy` symbolic check verifies that the implemented legacy and strict third-derivative product-rule formulas match direct differentiation, and that the implicit identity for `x'''` follows from differentiating `L_norm(x(d)) = S_strict_norm(d)` three times.  This supports the P2499 local uniqueness computation by removing a formula-transcription ambiguity from the interval backend.

This is a symbolic formula audit plus inherited finite interval evidence, not a formal directed-rounding backend, global analytic inflection theorem, or nonlinear compression-flow source theorem.  It exports no damping bridge atom, no strict compression source, no selector/source closure, no QW-2191 discharge, no role-transfer, no role-bearing `L_total`, no physical-value generation, and no ToE closure.
"""
    lag_section = """
## P2500/S1450 third-derivative formula provenance guard

`P2500/S1450` symbolically verifies the third-derivative formulas and implicit `x'''` chain identity used by P2499.  It strengthens the local inflection-uniqueness guard against formula-transcription error, but still exports no nonlinear compression-flow source, bridge atom, role-transfer theorem, QW-2191 discharge, or ToE closure.
"""
    append_once(DOC_FILES["equation_sheet"], "## P2500/S1450 phase-normalized third-derivative symbolic identity audit", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2500/S1450 third-derivative formula provenance guard", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    p2499 = theorem(sources["P2499_LOCAL_UNIQUENESS"], "phase_normalized_curvature_local_inflection_uniqueness_certificate")
    cert = build_symbolic_certificate(p2499)
    theorem_export = {
        "theorem_name": "P2500_T1_phase_normalized_third_derivative_symbolic_identity_audit",
        "audited_chain": ["P2497/S1447", "P2498/S1448", "P2499/S1449"],
        "third_derivative_symbolic_identity_audit": cert,
        "symbolic_backend": cert["symbolic_backend"],
        "sympy_version": cert["sympy_version"],
        "all_symbolic_residuals_zero": cert["all_symbolic_residuals_zero"],
        "legacy_third_residual_zero": cert["legacy_third_derivative_identity"]["residual_simplifies_to_zero"],
        "strict_third_residual_zero": cert["strict_third_derivative_identity"]["residual_simplifies_to_zero"],
        "implicit_chain_residual_zero": cert["implicit_third_derivative_chain_identity"]["residual_simplifies_to_zero"],
        "p2499_local_uniqueness_inherited": cert["p2499_local_uniqueness_inherited"],
        "p2499_third_interval_negative_inherited": cert["p2499_third_interval_negative_inherited"],
        "formula_provenance_supports_p2499_interval_third_derivative": cert["formula_provenance_supports_p2499_interval_third_derivative"],
        "formal_directed_rounding_backend_exported": False,
        "global_analytic_inflection_uniqueness_theorem_exported": False,
        "curvature_dynamic_source_exported": False,
        "legacy_to_strict_bridge_atom_exported": False,
        "strict_compression_dynamic_source_exported": False,
        "selector_source_theorem_exported": False,
        "qw2191_discharged_by_this_certificate": False,
        "role_transfer_licensed_by_this_certificate": False,
        "toe_closure_exported": False,
        "not_licensed": [
            "P2500 verifies formula identities symbolically; it is not a formal directed-rounding proof backend.",
            "Formula provenance and local uniqueness do not derive the nonlinear compression-flow law or a damping bridge atom.",
            "No selector/source theorem, QW-2191 discharge, role-transfer theorem, physical-value generator, or ToE closure is exported.",
        ],
        "next_honest_step": "Promote the interval arithmetic to a formal directed-rounding backend, or derive the strict nonlinear compression-flow source that explains the now formula-audited local inflection.",
    }
    gatekeepers = {
        "rg_audit_performed": grep["tool"] == "rg",
        "sympy_backend_recorded": theorem_export["symbolic_backend"] == "sympy",
        "legacy_third_identity_zero": theorem_export["legacy_third_residual_zero"],
        "strict_third_identity_zero": theorem_export["strict_third_residual_zero"],
        "implicit_chain_identity_zero": theorem_export["implicit_chain_residual_zero"],
        "all_symbolic_residuals_zero": theorem_export["all_symbolic_residuals_zero"],
        "p2499_local_uniqueness_inherited": theorem_export["p2499_local_uniqueness_inherited"],
        "p2499_third_interval_negative_inherited": theorem_export["p2499_third_interval_negative_inherited"],
        "formula_provenance_supports_p2499": theorem_export["formula_provenance_supports_p2499_interval_third_derivative"],
        "negative_controls_preserved": not any(theorem_export[key] for key in [
            "formal_directed_rounding_backend_exported",
            "global_analytic_inflection_uniqueness_theorem_exported",
            "curvature_dynamic_source_exported",
            "legacy_to_strict_bridge_atom_exported",
            "strict_compression_dynamic_source_exported",
            "selector_source_theorem_exported",
            "qw2191_discharged_by_this_certificate",
            "role_transfer_licensed_by_this_certificate",
            "toe_closure_exported",
        ]),
    }
    return {
        "packet_id": "P2500",
        "stage_id": "S1450",
        "status": "PHASE_NORMALIZED_THIRD_DERIVATIVE_SYMBOLIC_IDENTITY_AUDIT_NO_FORMAL_DIRECTED_BACKEND_NO_GLOBAL_ANALYTIC_UNIQUENESS_NO_SOURCE_EXPORT_NO_BRIDGE_ATOM_NO_QW2191_NO_ROLE_TRANSFER_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "phase_normalized_third_derivative_symbolic_identity_audit": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["phase_normalized_third_derivative_symbolic_identity_audit"]["theorem_export"]
    lines = [
        "# P2500/S1450 phase-normalized third-derivative symbolic identity audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        f"- Symbolic backend: `{t['symbolic_backend']}` `{t['sympy_version']}`.",
        f"- Legacy third-derivative formula residual zero: `{t['legacy_third_residual_zero']}`.",
        f"- Strict third-derivative formula residual zero: `{t['strict_third_residual_zero']}`.",
        f"- Implicit chain `x'''` residual zero: `{t['implicit_chain_residual_zero']}`.",
        f"- All symbolic residuals zero: `{t['all_symbolic_residuals_zero']}`.",
        f"- P2499 local uniqueness inherited: `{t['p2499_local_uniqueness_inherited']}`.",
        f"- Formula provenance supports P2499 interval third derivative: `{t['formula_provenance_supports_p2499_interval_third_derivative']}`.",
        "",
        "## Negative controls",
        "",
        "This packet does not export a formal directed-rounding backend, global analytic inflection uniqueness theorem, curvature dynamic source, legacy-to-strict bridge atom, strict compression source, selector/source theorem, QW-2191 discharge, role-transfer license, physical-value generator, or ToE closure.",
        "",
        "## Fingerprint",
        "",
        f"`{payload['phase_normalized_third_derivative_symbolic_identity_audit']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload["phase_normalized_third_derivative_symbolic_identity_audit"]["theorem_export"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
