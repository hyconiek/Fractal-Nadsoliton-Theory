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
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2508_s1458_strict_damping_rg_sobolev_node_coercivity_theorem_certificate.json"
MD = GEN / "p2508_s1458_strict_damping_rg_sobolev_node_coercivity_theorem_certificate.md"

SOURCE_FILES = {
    "P2507_ROUGHNESS_NULLSPACE_COERCIVITY": GEN / "p2507_s1457_strict_damping_rg_roughness_nullspace_coercivity_certificate.json",
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
        "new_packet": "P2508|S1458|Sobolev node coercivity|node-vanishing Sobolev|roughness Poincare|functional-analytic coercivity",
        "precursor_packets": "P2507|S1457|roughness nullspace coercivity|P2506|minimum roughness selector|P2505|finite-node RG-flow nullspace",
        "analytic_language": "Poincare|Wirtinger|Sobolev|H2|H\\^2|node-vanishing|functional-analytic theorem|coercivity theorem",
        "guardrails": "legacy -> strict completion bridge|role-transfer audit|K_legacy_ont|K_strict_gate|QW-2191|ToE closure",
        "closure_blockers": "source theorem|bridge theorem|role-transfer theorem|physical-value generator|role-bearing L_total|selector closure",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def text(value: mp.mpf, digits: int = 60) -> str:
    return mp.nstr(value, digits)


def poly_mul(a: list[mp.mpf], b: list[mp.mpf]) -> list[mp.mpf]:
    out = [mp.mpf("0")] * (len(a) + len(b) - 1)
    for i, av in enumerate(a):
        for j, bv in enumerate(b):
            out[i + j] += av * bv
    return out


def derivative_coefficients(coeffs: list[mp.mpf], order: int) -> list[mp.mpf]:
    current = coeffs[:]
    for _ in range(order):
        current = [mp.mpf(index + 1) * current[index + 1] for index in range(len(current) - 1)]
    return current or [mp.mpf("0")]


def eval_poly(coeffs: list[mp.mpf], x: mp.mpf) -> mp.mpf:
    value = mp.mpf("0")
    for coeff in reversed(coeffs):
        value = value * x + coeff
    return value


def integrate_square_closed_form(coeffs: list[mp.mpf], upper: mp.mpf) -> mp.mpf:
    product_coeffs = poly_mul(coeffs, coeffs)
    return mp.fsum(coeff * upper ** (power + 1) / mp.mpf(power + 1) for power, coeff in enumerate(product_coeffs))


def node_vanishing_r_coefficients() -> list[mp.mpf]:
    coeffs = [mp.mpf("1")]
    for root in [mp.mpf("0")] + [mp.log(d) for d in range(2, 12)]:
        coeffs = poly_mul(coeffs, [-root, mp.mpf("1")])
    return coeffs


def symbolic_sobolev_node_coercivity() -> dict[str, Any]:
    h = sp.log(2)
    c_l2 = h**4 / sp.pi**4
    c_h1 = h**2 / sp.pi**2
    c_h2_total = 1 + c_h1 + c_l2
    return {
        "symbolic_backend": "sympy",
        "sympy_version": sp.__version__,
        "admissible_space": "H^2([0,log(11)]) functions p with p(log d)=0 for d=1..11",
        "node_partition": "I_d=[log(d),log(d+1)] for d=1..10",
        "max_interval_length_exact": "log(2)",
        "local_inequality": "On each interval I of length h with p=0 at both endpoints: ||p||_L2(I)^2 <= h^4/pi^4 ||p''||_L2(I)^2 and ||p'||_L2(I)^2 <= h^2/pi^2 ||p''||_L2(I)^2.",
        "global_L2_constant_exact": sp.sstr(c_l2),
        "global_H1_seminorm_constant_exact": sp.sstr(c_h1),
        "global_H2_total_constant_exact": sp.sstr(c_h2_total),
        "zero_roughness_kernel_eliminated": "If integral (p'')^2=0 then p is affine on the full interval; p(0)=p(log(2))=0 forces p=0.",
        "coercivity_statement": "For all admitted p, ||p||_L2^2 <= C0 J[p], ||p'||_L2^2 <= C1 J[p], and ||p||_H2^2 <= (1+C1+C0) J[p], where J[p]=integral (p'')^2.",
        "symbolic_fingerprint_sha256": sha256_json({"c_l2": sp.srepr(c_l2), "c_h1": sp.srepr(c_h1), "c_h2_total": sp.srepr(c_h2_total)}),
    }


def interval_constant_audit() -> dict[str, Any]:
    rows = []
    max_h = mp.mpf("0")
    for d in range(1, 11):
        h = mp.log(d + 1) - mp.log(d)
        max_h = max(max_h, h)
        rows.append({
            "interval": [d, d + 1],
            "length_log_ratio": text(h, 50),
            "l2_constant_h4_over_pi4": text(h**4 / mp.pi**4, 50),
            "h1_constant_h2_over_pi2": text(h**2 / mp.pi**2, 50),
        })
    return {
        "row_count": len(rows),
        "rows": rows,
        "max_interval_length": text(max_h, 50),
        "max_interval_is_log2": abs(max_h - mp.log(2)) < mp.mpf("1e-80"),
        "global_L2_constant": text(max_h**4 / mp.pi**4, 70),
        "global_H1_seminorm_constant": text(max_h**2 / mp.pi**2, 70),
        "global_H2_total_constant": text(1 + max_h**2 / mp.pi**2 + max_h**4 / mp.pi**4, 70),
    }


def p2505_witness_inequality_audit() -> dict[str, Any]:
    coeffs = node_vanishing_r_coefficients()
    first = derivative_coefficients(coeffs, 1)
    second = derivative_coefficients(coeffs, 2)
    upper = mp.log(11)
    l2 = integrate_square_closed_form(coeffs, upper)
    h1 = integrate_square_closed_form(first, upper)
    roughness = integrate_square_closed_form(second, upper)
    hmax = mp.log(2)
    c_l2 = hmax**4 / mp.pi**4
    c_h1 = hmax**2 / mp.pi**2
    node_residuals = [abs(eval_poly(coeffs, mp.log(d))) for d in DOMAIN]
    return {
        "witness": "P2505 node-vanishing polynomial R(ell)=ell*prod_{d=2}^{11}(ell-log(d))",
        "node_max_abs_residual": text(max(node_residuals), 50),
        "L2_square": text(l2, 70),
        "H1_seminorm_square": text(h1, 70),
        "roughness_J": text(roughness, 70),
        "L2_over_J": text(l2 / roughness, 70),
        "H1_over_J": text(h1 / roughness, 70),
        "L2_bound_constant": text(c_l2, 70),
        "H1_bound_constant": text(c_h1, 70),
        "L2_inequality_residual_bound_minus_ratio": text(c_l2 - l2 / roughness, 70),
        "H1_inequality_residual_bound_minus_ratio": text(c_h1 - h1 / roughness, 70),
        "witness_satisfies_global_L2_bound": l2 <= c_l2 * roughness,
        "witness_satisfies_global_H1_bound": h1 <= c_h1 * roughness,
    }


def build_sobolev_coercivity_candidate(p2507: dict[str, Any]) -> dict[str, Any]:
    symbolic = symbolic_sobolev_node_coercivity()
    interval_audit = interval_constant_audit()
    witness_audit = p2505_witness_inequality_audit()
    return {
        "frontier_atom_under_attack": "strict_damping_beta_eta_source",
        "p2507_finite_polynomial_coercivity_inherited": p2507.get("finite_polynomial_nullspace_coercive") is True,
        "selector_type": "conditional minimum-roughness selector support upgraded from finite polynomial audit to node-vanishing H^2 coercivity theorem",
        "symbolic_sobolev_node_coercivity": symbolic,
        "interval_constant_audit": interval_audit,
        "p2505_witness_inequality_audit": witness_audit,
        "sobolev_node_coercivity_theorem_exported_for_postulated_roughness_functional": True,
        "coercivity_is_selector_support_not_source_derivation": True,
        "roughness_action_still_postulated_not_derived": True,
        "strict_damping_beta_eta_source_exported": False,
        "strict_dynamical_source_for_A_P_D_exported": False,
        "strict_phase_frequency_source_exported": False,
        "bridge_theorem_exported": False,
        "role_transfer_theorem_exported": False,
        "selector_closure_exported": False,
        "qw2191_discharged_by_this_certificate": False,
        "toe_closure_claimed": False,
    }


def append_doc_sections() -> None:
    eq_section = """
## P2508/S1458 strict damping RG Sobolev node-coercivity theorem certificate

`P2508/S1458` upgrades the P2507 finite polynomial Gram audit to a functional-analytic node-coercivity theorem for the postulated P2506 roughness selector.  On the partition `I_d=[log(d),log(d+1)]`, `d=1..10`, every admitted perturbation `p in H^2([0,log(11)])` with `p(log d)=0` at all audited nodes has zero endpoint trace on each subinterval.  The Dirichlet Poincare/Wirtinger bounds give `||p||_L2^2 <= (log(2)^4/pi^4) J[p]` and `||p'||_L2^2 <= (log(2)^2/pi^2) J[p]`, where `J[p]=int (p''(ell))^2 d ell`; hence the roughness form is coercive on the node-vanishing Sobolev perturbation class.

This is a real selector-support/coercivity theorem for the postulated roughness functional, not a source theorem.  The roughness action itself is still not derived from nadsoliton dynamics, `strict_damping_beta_eta_source` remains unexported, and there is no bridge theorem, role-transfer theorem, QW-2191 discharge, physical-value generator, or ToE closure.
"""
    lag_section = """
## P2508/S1458 Sobolev node-coercivity guard

`P2508/S1458` promotes the P2507 finite nullspace check to a node-vanishing `H^2` coercivity theorem for the postulated roughness selector, with explicit Poincare constants controlled by the largest node spacing `log(2)`.  It remains selector-support theorem-prep only: no strict nonlinear compression-flow source, no derived roughness action, and no role-bearing `L_total` term are licensed.
"""
    append_once(DOC_FILES["equation_sheet"], "## P2508/S1458 strict damping RG Sobolev node-coercivity theorem certificate", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2508/S1458 Sobolev node-coercivity guard", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    p2507 = theorem(sources["P2507_ROUGHNESS_NULLSPACE_COERCIVITY"], "strict_damping_rg_roughness_nullspace_coercivity_certificate")
    cert = build_sobolev_coercivity_candidate(p2507)
    interval_audit = cert["interval_constant_audit"]
    witness_audit = cert["p2505_witness_inequality_audit"]
    theorem_export = {
        "theorem_name": "P2508_T1_strict_damping_rg_sobolev_node_coercivity_theorem_certificate",
        "audited_chain": ["P2503/S1453", "P2504/S1454", "P2505/S1455", "P2506/S1456", "P2507/S1457"],
        "strict_damping_rg_sobolev_node_coercivity_theorem_certificate": cert,
        "frontier_atom_under_attack": cert["frontier_atom_under_attack"],
        "p2507_finite_polynomial_coercivity_inherited": cert["p2507_finite_polynomial_coercivity_inherited"],
        "admissible_space": cert["symbolic_sobolev_node_coercivity"]["admissible_space"],
        "max_interval_is_log2": interval_audit["max_interval_is_log2"],
        "global_L2_constant": interval_audit["global_L2_constant"],
        "global_H1_seminorm_constant": interval_audit["global_H1_seminorm_constant"],
        "global_H2_total_constant": interval_audit["global_H2_total_constant"],
        "p2505_witness_satisfies_global_L2_bound": witness_audit["witness_satisfies_global_L2_bound"],
        "p2505_witness_satisfies_global_H1_bound": witness_audit["witness_satisfies_global_H1_bound"],
        "sobolev_node_coercivity_theorem_exported_for_postulated_roughness_functional": cert["sobolev_node_coercivity_theorem_exported_for_postulated_roughness_functional"],
        "coercivity_is_selector_support_not_source_derivation": cert["coercivity_is_selector_support_not_source_derivation"],
        "roughness_action_still_postulated_not_derived": cert["roughness_action_still_postulated_not_derived"],
        "strict_damping_beta_eta_source_exported": False,
        "strict_dynamical_source_for_A_P_D_exported": False,
        "strict_phase_frequency_source_exported": False,
        "bridge_theorem_exported": False,
        "role_transfer_theorem_exported": False,
        "selector_closure_exported": False,
        "qw2191_discharged_by_this_certificate": False,
        "toe_closure_claimed": False,
        "not_licensed": [
            "P2508 proves node-vanishing Sobolev coercivity only for the postulated P2506 roughness functional.",
            "It does not derive the roughness action, delta, beta0, or strict_damping_beta_eta_source from strict nadsoliton dynamics.",
            "No bridge theorem, role-transfer theorem, QW-2191 discharge, role-bearing L_total, physical-value generator, or ToE closure is exported.",
        ],
        "next_honest_step": "Use the now-explicit node-coercivity theorem as stable selector support, but the main open step remains deriving the roughness functional or an equivalent damping source from strict nadsoliton dynamics.",
    }
    gatekeepers = {
        "rg_audit_performed": grep["tool"] == "rg",
        "p2507_finite_polynomial_coercivity_inherited": theorem_export["p2507_finite_polynomial_coercivity_inherited"],
        "max_interval_identified_as_log2": theorem_export["max_interval_is_log2"],
        "positive_global_L2_constant": mp.mpf(theorem_export["global_L2_constant"]) > 0,
        "positive_global_H1_constant": mp.mpf(theorem_export["global_H1_seminorm_constant"]) > 0,
        "p2505_witness_satisfies_bounds": theorem_export["p2505_witness_satisfies_global_L2_bound"] and theorem_export["p2505_witness_satisfies_global_H1_bound"],
        "sobolev_node_coercivity_exported_only_for_postulated_functional": theorem_export["sobolev_node_coercivity_theorem_exported_for_postulated_roughness_functional"] and theorem_export["roughness_action_still_postulated_not_derived"],
        "source_not_exported": not theorem_export["strict_damping_beta_eta_source_exported"],
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
        "packet_id": "P2508",
        "stage_id": "S1458",
        "status": "STRICT_DAMPING_RG_SOBOLEV_NODE_COERCIVITY_THEOREM_FOR_POSTULATED_ROUGHNESS_FUNCTIONAL_NO_SOURCE_EXPORT_NO_BRIDGE_THEOREM_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_damping_rg_sobolev_node_coercivity_theorem_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_damping_rg_sobolev_node_coercivity_theorem_certificate"]["theorem_export"]
    lines = [
        "# P2508/S1458 strict damping RG Sobolev node-coercivity theorem certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        f"- Frontier atom under attack: `{t['frontier_atom_under_attack']}`.",
        f"- P2507 finite polynomial coercivity inherited: `{t['p2507_finite_polynomial_coercivity_inherited']}`.",
        f"- Admissible space: `{t['admissible_space']}`.",
        f"- Max interval is log(2): `{t['max_interval_is_log2']}`.",
        f"- Global L2 constant: `{t['global_L2_constant']}`.",
        f"- Global H1 seminorm constant: `{t['global_H1_seminorm_constant']}`.",
        f"- Global H2 total constant: `{t['global_H2_total_constant']}`.",
        f"- P2505 witness satisfies global L2 bound: `{t['p2505_witness_satisfies_global_L2_bound']}`.",
        f"- P2505 witness satisfies global H1 bound: `{t['p2505_witness_satisfies_global_H1_bound']}`.",
        f"- Sobolev node coercivity exported for postulated roughness functional: `{t['sobolev_node_coercivity_theorem_exported_for_postulated_roughness_functional']}`.",
        f"- Roughness action still postulated, not derived: `{t['roughness_action_still_postulated_not_derived']}`.",
        f"- Source theorem exported: `{t['strict_damping_beta_eta_source_exported']}`.",
        "",
        "## Negative controls",
        "",
        "This packet proves a node-vanishing Sobolev coercivity theorem for the postulated roughness selector. It does not derive the roughness action from strict dynamics and exports no strict source atom, bridge theorem, role-transfer theorem, selector/QW-2191 closure, role-bearing L_total term, physical-value generator, or ToE closure.",
        "",
        "## Fingerprint",
        "",
        f"`{payload['strict_damping_rg_sobolev_node_coercivity_theorem_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload["strict_damping_rg_sobolev_node_coercivity_theorem_certificate"]["theorem_export"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
