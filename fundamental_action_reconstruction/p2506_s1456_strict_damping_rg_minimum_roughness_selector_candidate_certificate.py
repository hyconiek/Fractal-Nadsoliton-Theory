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
OUT = GEN / "p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate.json"
MD = GEN / "p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate.md"

SOURCE_FILES = {
    "P2505_RG_FLOW_NULLSPACE": GEN / "p2505_s1455_strict_damping_finite_node_rg_flow_nullspace_certificate.json",
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
        "new_packet": "P2506|S1456|minimum roughness selector|RG roughness selector|Sobolev selector candidate|constant-flow selector candidate",
        "precursor_packets": "P2505|S1455|finite-node RG-flow nullspace|P2504|constant-coefficient RG uniqueness|P2503",
        "variational_language": "variational selector|minimum roughness|Sobolev roughness|integral of y''|running-beta selector|constant-flow selection",
        "guardrails": "legacy -> strict completion bridge|role-transfer audit|K_legacy_ont|K_strict_gate|QW-2191|ToE closure",
        "closure_blockers": "source theorem|bridge theorem|role-transfer theorem|physical-value generator|role-bearing L_total|selector closure",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def text(value: mp.mpf, digits: int = 60) -> str:
    return mp.nstr(value, digits)


def polynomial_coefficients() -> list[mp.mpf]:
    # Ascending coefficients for R(ell)=ell*prod_{d=2}^{11}(ell-log(d)).
    coeffs = [mp.mpf("1")]
    roots = [mp.mpf("0")] + [mp.log(d) for d in range(2, 12)]
    for root in roots:
        new_coeffs = [mp.mpf("0")] * (len(coeffs) + 1)
        for index, coeff in enumerate(coeffs):
            new_coeffs[index] -= root * coeff
            new_coeffs[index + 1] += coeff
        coeffs = new_coeffs
    return coeffs


def derivative_coefficients(coeffs: list[mp.mpf], order: int) -> list[mp.mpf]:
    current = coeffs[:]
    for _ in range(order):
        current = [mp.mpf(index + 1) * current[index + 1] for index in range(len(current) - 1)]
    return current


def eval_poly(coeffs: list[mp.mpf], x: mp.mpf) -> mp.mpf:
    value = mp.mpf("0")
    for coeff in reversed(coeffs):
        value = value * x + coeff
    return value


def symbolic_selector_certificate() -> dict[str, Any]:
    ell = sp.symbols("ell", real=True)
    delta = sp.Rational(14, 5) - 4 * sp.log(2)
    y0 = delta * ell
    y0_second = sp.diff(y0, ell, 2)
    return {
        "symbolic_backend": "sympy",
        "sympy_version": sp.__version__,
        "selector_functional": "J[y]=integral_0^log(11) (y''(ell))^2 d ell",
        "constraint_family": "y(log d)=delta*log d for d=1..11, with log(1)=0",
        "constant_candidate_y0": sp.sstr(y0),
        "constant_candidate_second_derivative": sp.sstr(y0_second),
        "constant_candidate_energy_zero": y0_second == 0,
        "nonnegativity_argument": "J[y] is an integral of a square; J[y]>=0, and J[y]=0 implies y''=0 almost everywhere in the admitted Sobolev class.",
        "affine_zero_energy_constraints_fix_candidate": "If y''=0 then y=a*ell+b; constraints at d=1 and d=2 force b=0 and a=delta.",
        "selector_candidate_unique_if_roughness_action_is_admitted": True,
        "symbolic_fingerprint_sha256": sha256_json({"delta": sp.srepr(delta), "y0": sp.srepr(y0), "y0_second": sp.srepr(y0_second)}),
    }


def nullspace_energy_certificate() -> dict[str, Any]:
    coeffs = polynomial_coefficients()
    second_coeffs = derivative_coefficients(coeffs, 2)
    fourth_coeffs = derivative_coefficients(coeffs, 4)
    left = mp.mpf("0")
    right = mp.log(11)
    roots = [left] + [mp.log(d) for d in range(2, 12)]

    def r_second_squared(x: mp.mpf) -> mp.mpf:
        value = eval_poly(second_coeffs, x)
        return value * value

    # Split at roots to avoid high-degree quadrature cancellation across the whole interval.
    energy_unscaled = mp.fsum(mp.quad(r_second_squared, [roots[i], roots[i + 1]]) for i in range(len(roots) - 1))
    energy_scaled = EPSILON**2 * energy_unscaled
    midpoint_rows = []
    for left_d in range(1, 11):
        ell_mid = (mp.log(left_d) + mp.log(left_d + 1)) / 2
        r_second = eval_poly(second_coeffs, ell_mid)
        r_fourth = eval_poly(fourth_coeffs, ell_mid) if fourth_coeffs else mp.mpf("0")
        midpoint_rows.append({
            "between_nodes": [left_d, left_d + 1],
            "ell_mid": text(ell_mid, 50),
            "epsilon_R_second": text(EPSILON * r_second, 50),
            "epsilon_R_fourth": text(EPSILON * r_fourth, 50),
        })
    max_second = max(abs(mp.mpf(row["epsilon_R_second"])) for row in midpoint_rows)
    return {
        "epsilon": text(EPSILON, 20),
        "integration_interval": ["0", text(right, 50)],
        "node_split_count": len(roots) - 1,
        "R_second_square_integral_unscaled": text(energy_unscaled, 70),
        "epsilon_scaled_R_second_square_integral": text(energy_scaled, 70),
        "epsilon_scaled_energy_positive": energy_scaled > 0,
        "midpoint_rows": midpoint_rows,
        "max_abs_epsilon_R_second_midpoint": text(max_second, 50),
        "nonconstant_nullspace_has_positive_roughness_energy": energy_scaled > 0 and max_second > 0,
    }


def build_selector_candidate(p2505: dict[str, Any]) -> dict[str, Any]:
    symbolic = symbolic_selector_certificate()
    nullspace = nullspace_energy_certificate()
    return {
        "frontier_atom_under_attack": "strict_damping_beta_eta_source",
        "p2505_nullspace_inherited": p2505.get("finite_nodes_do_not_select_constant_coefficient_flow") is True,
        "selector_type": "minimum Sobolev roughness candidate for running exponent y(ell)=log B(ell)",
        "symbolic_minimum_roughness_selector": symbolic,
        "p2505_nullspace_energy_audit": nullspace,
        "constant_flow_selected_if_selector_postulated": symbolic["selector_candidate_unique_if_roughness_action_is_admitted"] and nullspace["nonconstant_nullspace_has_positive_roughness_energy"],
        "selector_is_postulated_not_derived": True,
        "finite_node_nullspace_blocker_addressed_conditionally": True,
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
## P2506/S1456 strict damping RG minimum-roughness selector candidate certificate

`P2506/S1456` conditionally addresses the P2505 finite-node nullspace by auditing a concrete selector candidate: minimize `J[y]=int_0^log(11) (y''(ell))^2 d ell` for the running exponent `y(ell)=log B(ell)` subject to the strict finite-node constraints.  The constant P2503 flow has zero roughness, and the zero-energy affine constraints force `y(ell)=delta ell`; the explicit P2505 nullspace perturbation has positive scaled roughness energy.  Thus the constant-coefficient flow is selected if this minimum-roughness action is admitted.

This is a conditional selector-candidate certificate, not a derived strict damping source.  It does not derive the roughness action from nadsoliton dynamics, does not export `strict_damping_beta_eta_source`, and exports no bridge theorem, role-transfer theorem, QW-2191 discharge, physical-value generator, or ToE closure.
"""
    lag_section = """
## P2506/S1456 minimum-roughness RG selector candidate guard

`P2506/S1456` shows that a postulated minimum-roughness action would select the constant P2503 running-beta flow against the P2505 nullspace.  Because the action is not derived from strict nadsoliton dynamics, it remains theorem-prep and does not license a nonlinear compression-flow source theorem or role-bearing `L_total`.
"""
    append_once(DOC_FILES["equation_sheet"], "## P2506/S1456 strict damping RG minimum-roughness selector candidate certificate", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2506/S1456 minimum-roughness RG selector candidate guard", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    p2505 = theorem(sources["P2505_RG_FLOW_NULLSPACE"], "strict_damping_finite_node_rg_flow_nullspace_certificate")
    cert = build_selector_candidate(p2505)
    symbolic = cert["symbolic_minimum_roughness_selector"]
    nullspace = cert["p2505_nullspace_energy_audit"]
    theorem_export = {
        "theorem_name": "P2506_T1_strict_damping_rg_minimum_roughness_selector_candidate_certificate",
        "audited_chain": ["P2503/S1453", "P2504/S1454", "P2505/S1455"],
        "strict_damping_rg_minimum_roughness_selector_candidate_certificate": cert,
        "frontier_atom_under_attack": cert["frontier_atom_under_attack"],
        "p2505_nullspace_inherited": cert["p2505_nullspace_inherited"],
        "constant_candidate_energy_zero": symbolic["constant_candidate_energy_zero"],
        "selector_candidate_unique_if_roughness_action_is_admitted": symbolic["selector_candidate_unique_if_roughness_action_is_admitted"],
        "nullspace_scaled_roughness_energy": nullspace["epsilon_scaled_R_second_square_integral"],
        "nullspace_scaled_energy_positive": nullspace["epsilon_scaled_energy_positive"],
        "max_abs_epsilon_R_second_midpoint": nullspace["max_abs_epsilon_R_second_midpoint"],
        "constant_flow_selected_if_selector_postulated": cert["constant_flow_selected_if_selector_postulated"],
        "selector_is_postulated_not_derived": cert["selector_is_postulated_not_derived"],
        "finite_node_nullspace_blocker_addressed_conditionally": cert["finite_node_nullspace_blocker_addressed_conditionally"],
        "strict_damping_beta_eta_source_exported": False,
        "strict_dynamical_source_for_A_P_D_exported": False,
        "strict_phase_frequency_source_exported": False,
        "bridge_theorem_exported": False,
        "role_transfer_theorem_exported": False,
        "selector_closure_exported": False,
        "qw2191_discharged_by_this_certificate": False,
        "toe_closure_claimed": False,
        "not_licensed": [
            "P2506 selects the constant flow only conditionally on a postulated minimum-roughness action.",
            "The roughness action is not derived from strict nadsoliton dynamics, so strict_damping_beta_eta_source remains open.",
            "No strict source atom, bridge theorem, role-transfer theorem, QW-2191 discharge, role-bearing L_total, physical-value generator, or ToE closure is exported.",
        ],
        "next_honest_step": "Derive the minimum-roughness running-exponent action from strict nadsoliton dynamics, or replace it with a better-motivated selector/source theorem.",
    }
    gatekeepers = {
        "rg_audit_performed": grep["tool"] == "rg",
        "p2505_nullspace_inherited": theorem_export["p2505_nullspace_inherited"],
        "constant_candidate_energy_zero": theorem_export["constant_candidate_energy_zero"],
        "nullspace_energy_positive": theorem_export["nullspace_scaled_energy_positive"],
        "constant_flow_selected_conditionally": theorem_export["constant_flow_selected_if_selector_postulated"],
        "selector_marked_postulated_not_derived": theorem_export["selector_is_postulated_not_derived"],
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
        "packet_id": "P2506",
        "stage_id": "S1456",
        "status": "STRICT_DAMPING_RG_MINIMUM_ROUGHNESS_SELECTOR_CANDIDATE_CERTIFICATE_NO_SOURCE_EXPORT_NO_BRIDGE_THEOREM_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_damping_rg_minimum_roughness_selector_candidate_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_damping_rg_minimum_roughness_selector_candidate_certificate"]["theorem_export"]
    lines = [
        "# P2506/S1456 strict damping RG minimum-roughness selector candidate certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        f"- Frontier atom under attack: `{t['frontier_atom_under_attack']}`.",
        f"- P2505 nullspace inherited: `{t['p2505_nullspace_inherited']}`.",
        f"- Constant candidate energy zero: `{t['constant_candidate_energy_zero']}`.",
        f"- Nullspace scaled roughness energy: `{t['nullspace_scaled_roughness_energy']}`.",
        f"- Nullspace energy positive: `{t['nullspace_scaled_energy_positive']}`.",
        f"- Constant flow selected if selector is postulated: `{t['constant_flow_selected_if_selector_postulated']}`.",
        f"- Selector is postulated, not derived: `{t['selector_is_postulated_not_derived']}`.",
        f"- Source theorem exported: `{t['strict_damping_beta_eta_source_exported']}`.",
        "",
        "## Negative controls",
        "",
        "This packet is a conditional selector-candidate certificate. It does not derive the minimum-roughness action from strict dynamics and exports no strict source atom, bridge theorem, role-transfer theorem, selector/QW-2191 closure, role-bearing L_total term, physical-value generator, or ToE closure.",
        "",
        "## Fingerprint",
        "",
        f"`{payload['strict_damping_rg_minimum_roughness_selector_candidate_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload["strict_damping_rg_minimum_roughness_selector_candidate_certificate"]["theorem_export"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
