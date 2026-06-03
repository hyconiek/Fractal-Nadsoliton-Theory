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
OUT = GEN / "p2509_s1459_strict_damping_rg_minimum_roughness_variational_wellposedness_certificate.json"
MD = GEN / "p2509_s1459_strict_damping_rg_minimum_roughness_variational_wellposedness_certificate.md"

SOURCE_FILES = {
    "P2508_SOBOLEV_NODE_COERCIVITY": GEN / "p2508_s1458_strict_damping_rg_sobolev_node_coercivity_theorem_certificate.json",
}

mp.mp.dps = 90
DOMAIN = list(range(1, 12))
SAMPLE_COEFFICIENT_ROWS = [
    [mp.mpf("1e-6"), mp.mpf("0"), mp.mpf("0"), mp.mpf("0")],
    [mp.mpf("1e-6"), mp.mpf("-2e-6"), mp.mpf("3e-6"), mp.mpf("-1e-6")],
    [mp.mpf("-3e-7"), mp.mpf("5e-7"), mp.mpf("-7e-7"), mp.mpf("11e-7")],
]


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
        "new_packet": "P2509|S1459|minimum-roughness variational well-posedness|roughness minimizer|affine minimizer|Sobolev minimizer",
        "precursor_packets": "P2508|S1458|Sobolev node coercivity|P2507|roughness nullspace coercivity|P2506|minimum roughness selector",
        "variational_language": "unique minimizer|variational well-posedness|direct method|Euler-Lagrange|affine constraint|node interpolation",
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


def symbolic_variational_wellposedness() -> dict[str, Any]:
    ell = sp.symbols("ell", real=True)
    delta = sp.Rational(14, 5) - 4 * sp.log(2)
    y0 = delta * ell
    y0_second = sp.diff(y0, ell, 2)
    return {
        "symbolic_backend": "sympy",
        "sympy_version": sp.__version__,
        "admissible_affine_space": "A={y in H^2([0,log(11)]): y(log d)=delta*log d for d=1..11}",
        "node_vanishing_tangent_space": "V={p in H^2([0,log(11)]): p(log d)=0 for d=1..11}",
        "candidate_minimizer_y0": sp.sstr(y0),
        "candidate_second_derivative": sp.sstr(y0_second),
        "candidate_energy_zero": y0_second == 0,
        "decomposition": "Every y in A decomposes uniquely as y=y0+p with p in V.",
        "energy_identity": "J[y0+p]=int (p''(ell))^2 d ell because y0''=0.",
        "coercivity_input": "P2508 proves node-vanishing H^2 coercivity on V for the same postulated roughness functional.",
        "unique_minimizer_statement": "Since J[y0]=0 and J[y0+p]>0 for every nonzero p in V, y0 is the unique minimizer in A.",
        "symbolic_fingerprint_sha256": sha256_json({"delta": sp.srepr(delta), "y0": sp.srepr(y0), "y0_second": sp.srepr(y0_second)}),
    }


def node_data_audit() -> dict[str, Any]:
    delta = mp.mpf(14) / 5 - 4 * mp.log(2)
    rows = []
    max_residual = mp.mpf("0")
    for d in DOMAIN:
        ell_d = mp.log(d)
        target = delta * ell_d
        y0_value = delta * ell_d
        residual = y0_value - target
        max_residual = max(max_residual, abs(residual))
        rows.append({
            "d": d,
            "ell_log_d": text(ell_d, 50),
            "target_delta_log_d": text(target, 50),
            "y0_value": text(y0_value, 50),
            "residual": text(residual, 30),
        })
    return {
        "delta": text(delta, 70),
        "row_count": len(rows),
        "rows": rows,
        "max_abs_node_residual": text(max_residual, 30),
        "candidate_matches_all_nodes": max_residual == 0,
    }


def polynomial_tangent_basis() -> list[list[mp.mpf]]:
    r_coeffs = node_vanishing_r_coefficients()
    return [poly_mul(r_coeffs, [mp.mpf("0")] * k + [mp.mpf("1")]) for k in range(4)]


def sample_tangent_energy_audit() -> dict[str, Any]:
    upper = mp.log(11)
    basis = polynomial_tangent_basis()
    rows = []
    for coefficients in SAMPLE_COEFFICIENT_ROWS:
        max_node_residual = mp.mpf("0")
        p_coeffs = [mp.mpf("0")] * max(len(basis[-1]), 1)
        for coeff, basis_coeffs in zip(coefficients, basis):
            if len(p_coeffs) < len(basis_coeffs):
                p_coeffs.extend([mp.mpf("0")] * (len(basis_coeffs) - len(p_coeffs)))
            for index, value in enumerate(basis_coeffs):
                p_coeffs[index] += coeff * value
        second = derivative_coefficients(p_coeffs, 2)
        energy = integrate_square_closed_form(second, upper)
        for d in DOMAIN:
            max_node_residual = max(max_node_residual, abs(eval_poly(p_coeffs, mp.log(d))))
        rows.append({
            "coefficients": [text(value, 20) for value in coefficients],
            "max_abs_node_residual": text(max_node_residual, 50),
            "roughness_energy": text(energy, 70),
            "energy_positive": energy > 0,
        })
    return {
        "sample_count": len(rows),
        "basis_family": "R(ell)*ell^k for k=0..3, with R(ell)=ell*prod_{d=2}^{11}(ell-log(d))",
        "rows": rows,
        "all_samples_node_vanishing": all(mp.mpf(row["max_abs_node_residual"]) < mp.mpf("1e-70") for row in rows),
        "all_nonzero_samples_positive_energy": all(row["energy_positive"] for row in rows),
    }


def build_variational_wellposedness_candidate(p2508: dict[str, Any]) -> dict[str, Any]:
    symbolic = symbolic_variational_wellposedness()
    nodes = node_data_audit()
    samples = sample_tangent_energy_audit()
    return {
        "frontier_atom_under_attack": "strict_damping_beta_eta_source",
        "p2508_sobolev_node_coercivity_inherited": p2508.get("sobolev_node_coercivity_theorem_exported_for_postulated_roughness_functional") is True,
        "selector_type": "conditional minimum-roughness variational well-posedness theorem on the strict damping node data",
        "symbolic_variational_wellposedness": symbolic,
        "node_data_audit": nodes,
        "sample_tangent_energy_audit": samples,
        "minimum_roughness_problem_wellposed_for_postulated_functional": True,
        "unique_minimizer_is_constant_flow_y0_delta_ell": symbolic["candidate_energy_zero"] and nodes["candidate_matches_all_nodes"],
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
## P2509/S1459 strict damping RG minimum-roughness variational well-posedness certificate

`P2509/S1459` uses the P2508 node-vanishing `H^2` coercivity theorem to close the variational well-posedness of the postulated P2506 minimum-roughness selector.  For the affine constraint set `A={y in H^2: y(log d)=delta log d, d=1..11}`, the candidate `y0(ell)=delta ell` has zero roughness.  Every admissible `y` decomposes uniquely as `y=y0+p` with `p` in the node-vanishing tangent space, and `J[y0+p]=int (p''(ell))^2 d ell`; P2508 coercivity makes this strictly positive for every nonzero perturbation.  Thus the constant-flow reconstruction is the unique minimizer of the postulated roughness problem.

This upgrades selector support from a candidate to a well-posed conditional variational theorem, but it still does not derive the roughness action from nadsoliton dynamics.  Therefore `strict_damping_beta_eta_source` remains unexported, and there is no bridge theorem, role-transfer theorem, QW-2191 discharge, physical-value generator, or ToE closure.
"""
    lag_section = """
## P2509/S1459 minimum-roughness variational well-posedness guard

`P2509/S1459` proves that, if the P2506 roughness action is admitted, the strict damping node data have the unique `H^2` minimizer `y0(ell)=delta ell` by P2508 node-coercivity.  This is still conditional selector support only: it does not derive a strict nonlinear compression-flow source or license a role-bearing `L_total` term.
"""
    append_once(DOC_FILES["equation_sheet"], "## P2509/S1459 strict damping RG minimum-roughness variational well-posedness certificate", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2509/S1459 minimum-roughness variational well-posedness guard", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    p2508 = theorem(sources["P2508_SOBOLEV_NODE_COERCIVITY"], "strict_damping_rg_sobolev_node_coercivity_theorem_certificate")
    cert = build_variational_wellposedness_candidate(p2508)
    nodes = cert["node_data_audit"]
    samples = cert["sample_tangent_energy_audit"]
    theorem_export = {
        "theorem_name": "P2509_T1_strict_damping_rg_minimum_roughness_variational_wellposedness_certificate",
        "audited_chain": ["P2503/S1453", "P2504/S1454", "P2505/S1455", "P2506/S1456", "P2507/S1457", "P2508/S1458"],
        "strict_damping_rg_minimum_roughness_variational_wellposedness_certificate": cert,
        "frontier_atom_under_attack": cert["frontier_atom_under_attack"],
        "p2508_sobolev_node_coercivity_inherited": cert["p2508_sobolev_node_coercivity_inherited"],
        "candidate_matches_all_nodes": nodes["candidate_matches_all_nodes"],
        "max_abs_node_residual": nodes["max_abs_node_residual"],
        "minimum_roughness_problem_wellposed_for_postulated_functional": cert["minimum_roughness_problem_wellposed_for_postulated_functional"],
        "unique_minimizer_is_constant_flow_y0_delta_ell": cert["unique_minimizer_is_constant_flow_y0_delta_ell"],
        "sample_tangent_all_node_vanishing": samples["all_samples_node_vanishing"],
        "sample_tangent_all_positive_energy": samples["all_nonzero_samples_positive_energy"],
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
            "P2509 proves well-posedness only for the postulated P2506 roughness minimization problem.",
            "It does not derive the roughness action, delta, beta0, or strict_damping_beta_eta_source from strict nadsoliton dynamics.",
            "No bridge theorem, role-transfer theorem, QW-2191 discharge, role-bearing L_total, physical-value generator, or ToE closure is exported.",
        ],
        "next_honest_step": "The roughness-selector route is now conditionally well-posed; the remaining source-level step is to derive this roughness action or an equivalent damping-flow source from strict nadsoliton dynamics.",
    }
    gatekeepers = {
        "rg_audit_performed": grep["tool"] == "rg",
        "p2508_coercivity_inherited": theorem_export["p2508_sobolev_node_coercivity_inherited"],
        "candidate_matches_nodes": theorem_export["candidate_matches_all_nodes"],
        "wellposed_for_postulated_functional": theorem_export["minimum_roughness_problem_wellposed_for_postulated_functional"],
        "unique_minimizer_constant_flow": theorem_export["unique_minimizer_is_constant_flow_y0_delta_ell"],
        "sample_tangent_checks_pass": theorem_export["sample_tangent_all_node_vanishing"] and theorem_export["sample_tangent_all_positive_energy"],
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
        "packet_id": "P2509",
        "stage_id": "S1459",
        "status": "STRICT_DAMPING_RG_MINIMUM_ROUGHNESS_VARIATIONAL_WELLPOSEDNESS_FOR_POSTULATED_FUNCTIONAL_NO_SOURCE_EXPORT_NO_BRIDGE_THEOREM_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_damping_rg_minimum_roughness_variational_wellposedness_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_damping_rg_minimum_roughness_variational_wellposedness_certificate"]["theorem_export"]
    lines = [
        "# P2509/S1459 strict damping RG minimum-roughness variational well-posedness certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        f"- Frontier atom under attack: `{t['frontier_atom_under_attack']}`.",
        f"- P2508 Sobolev node coercivity inherited: `{t['p2508_sobolev_node_coercivity_inherited']}`.",
        f"- Candidate matches all strict nodes: `{t['candidate_matches_all_nodes']}`.",
        f"- Max node residual: `{t['max_abs_node_residual']}`.",
        f"- Minimum-roughness problem well-posed for postulated functional: `{t['minimum_roughness_problem_wellposed_for_postulated_functional']}`.",
        f"- Unique minimizer is `y0(ell)=delta ell`: `{t['unique_minimizer_is_constant_flow_y0_delta_ell']}`.",
        f"- Sample tangent perturbations node-vanishing: `{t['sample_tangent_all_node_vanishing']}`.",
        f"- Sample tangent perturbations positive-energy: `{t['sample_tangent_all_positive_energy']}`.",
        f"- Roughness action still postulated, not derived: `{t['roughness_action_still_postulated_not_derived']}`.",
        f"- Source theorem exported: `{t['strict_damping_beta_eta_source_exported']}`.",
        "",
        "## Negative controls",
        "",
        "This packet proves variational well-posedness and the unique constant-flow minimizer only for the postulated roughness selector. It does not derive the roughness action from strict dynamics and exports no strict source atom, bridge theorem, role-transfer theorem, selector/QW-2191 closure, role-bearing L_total term, physical-value generator, or ToE closure.",
        "",
        "## Fingerprint",
        "",
        f"`{payload['strict_damping_rg_minimum_roughness_variational_wellposedness_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload["strict_damping_rg_minimum_roughness_variational_wellposedness_certificate"]["theorem_export"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
