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
from p2509_s1459_strict_damping_rg_minimum_roughness_variational_wellposedness_certificate import (
    derivative_coefficients,
    eval_poly,
    integrate_square_closed_form,
    node_vanishing_r_coefficients,
    poly_mul,
)

GEN = ROOT / "generated"
OUT = GEN / "p2512_s1462_strict_damping_rg_quadratic_source_admissibility_audit.json"
MD = GEN / "p2512_s1462_strict_damping_rg_quadratic_source_admissibility_audit.md"

SOURCE_FILES = {
    "P2511_NATURAL_SPLINE_COLLAPSE": GEN / "p2511_s1461_strict_damping_rg_natural_spline_collapse_certificate.json",
}

mp.mp.dps = 100
DOMAIN = list(range(1, 12))
BASIS_DEGREES = list(range(4))


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
        "new_packet": "P2512|S1462|source operator admissibility|derivative-only quadratic|mass term obstruction|roughness source ambiguity",
        "precursor_packets": "P2511|S1461|natural spline collapse|P2510|roughness KKT stationarity|P2509|minimum-roughness variational well-posedness",
        "source_operator_language": "quadratic source|local quadratic operator|first variation|stationarity residual|mass term|derivative-only",
        "guardrails": "legacy -> strict completion bridge|role-transfer audit|K_legacy_ont|K_strict_gate|QW-2191|ToE closure",
        "closure_blockers": "source theorem|bridge theorem|role-transfer theorem|physical-value generator|role-bearing L_total|selector closure",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def text(value: mp.mpf, digits: int = 70) -> str:
    return mp.nstr(value, digits)


def integrate_product_closed_form(a: list[mp.mpf], b: list[mp.mpf], upper: mp.mpf) -> mp.mpf:
    product = poly_mul(a, b)
    return mp.fsum(coeff * upper ** (power + 1) / mp.mpf(power + 1) for power, coeff in enumerate(product))


def monomial_times(coeffs: list[mp.mpf], degree: int) -> list[mp.mpf]:
    return [mp.mpf("0")] * degree + coeffs[:]


def basis_polynomials() -> list[list[mp.mpf]]:
    r = node_vanishing_r_coefficients()
    return [monomial_times(r, degree) for degree in BASIS_DEGREES]


def symbolic_first_variation_audit() -> dict[str, Any]:
    ell, c0, c1, c2 = sp.symbols("ell c0 c1 c2", real=True)
    delta = sp.Rational(14, 5) - 4 * sp.log(2)
    y0 = delta * ell
    formal_residual = c0 * sp.Integral(y0 * sp.Function("p")(ell), ell) + c1 * delta * sp.Integral(sp.diff(sp.Function("p")(ell), ell), ell)
    return {
        "symbolic_backend": "sympy",
        "sympy_version": sp.__version__,
        "candidate_y0": sp.sstr(y0),
        "quadratic_source_family": "S[y]=1/2 int (c0*y^2 + c1*(y')^2 + c2*(y'')^2) d ell under fixed node constraints",
        "first_variation_on_node_vanishing_p": "delta S[y0](p)=c0*int y0*p + c1*delta*(p(L)-p(0)) + c2*int y0''*p''; node vanishing gives p(0)=p(L)=0 and y0''=0.",
        "derivative_only_stationarity_condition": "If c0=0, every derivative-only quadratic c1*(y')^2+c2*(y'')^2 has zero first variation at y0 on the node-vanishing tangent space.",
        "mass_term_obstruction_condition": "If c0 != 0, stationarity would require int ell*p(ell) d ell=0 for every node-vanishing p, which is false by the polynomial witness audit.",
        "formal_residual_template": sp.sstr(formal_residual),
        "symbolic_fingerprint_sha256": sha256_json({"delta": sp.srepr(delta), "template": sp.srepr(formal_residual)}),
    }


def stationarity_witness_audit() -> dict[str, Any]:
    upper = mp.log(11)
    delta = mp.mpf(14) / 5 - 4 * mp.log(2)
    ell_coeffs = [mp.mpf("0"), mp.mpf("1")]
    y0_coeffs = [mp.mpf("0"), delta]
    y0_prime_coeffs = [delta]
    y0_second_coeffs = [mp.mpf("0")]
    rows = []
    max_node_residual = mp.mpf("0")
    max_abs_mass_moment = mp.mpf("0")
    max_abs_derivative_residual = mp.mpf("0")
    max_abs_roughness_residual = mp.mpf("0")
    nonzero_mass_rows = 0
    for degree, p_coeffs in zip(BASIS_DEGREES, basis_polynomials()):
        p_prime = derivative_coefficients(p_coeffs, 1)
        p_second = derivative_coefficients(p_coeffs, 2)
        mass_moment = integrate_product_closed_form(y0_coeffs, p_coeffs, upper)
        ell_moment = integrate_product_closed_form(ell_coeffs, p_coeffs, upper)
        derivative_residual = integrate_product_closed_form(y0_prime_coeffs, p_prime, upper)
        roughness_residual = integrate_product_closed_form(y0_second_coeffs, p_second, upper)
        node_residual = max(abs(eval_poly(p_coeffs, mp.log(d))) for d in DOMAIN)
        energy = integrate_square_closed_form(p_second, upper)
        max_node_residual = max(max_node_residual, node_residual)
        max_abs_mass_moment = max(max_abs_mass_moment, abs(mass_moment))
        max_abs_derivative_residual = max(max_abs_derivative_residual, abs(derivative_residual))
        max_abs_roughness_residual = max(max_abs_roughness_residual, abs(roughness_residual))
        if abs(mass_moment) > mp.mpf("1e-40"):
            nonzero_mass_rows += 1
        rows.append({
            "basis_element": f"R(ell)*ell^{degree}",
            "max_abs_node_residual": text(node_residual, 60),
            "mass_variation_moment_int_y0_p": text(mass_moment, 70),
            "ell_moment_int_ell_p": text(ell_moment, 70),
            "first_derivative_variation_int_y0_prime_p_prime": text(derivative_residual, 70),
            "roughness_variation_int_y0_second_p_second": text(roughness_residual, 70),
            "roughness_energy_of_witness": text(energy, 70),
            "mass_term_obstructs_stationarity_on_this_witness": abs(mass_moment) > mp.mpf("1e-40"),
            "derivative_only_terms_stationary_on_this_witness": abs(derivative_residual) < mp.mpf("1e-80") and abs(roughness_residual) < mp.mpf("1e-80"),
        })
    return {
        "basis_family": "R(ell)*ell^k, k=0..3, with R(ell)=ell*prod_{d=2}^{11}(ell-log(d))",
        "row_count": len(rows),
        "rows": rows,
        "max_abs_node_residual": text(max_node_residual, 60),
        "max_abs_mass_variation_moment": text(max_abs_mass_moment, 70),
        "max_abs_first_derivative_variation": text(max_abs_derivative_residual, 70),
        "max_abs_roughness_variation": text(max_abs_roughness_residual, 70),
        "nonzero_mass_obstruction_row_count": nonzero_mass_rows,
        "all_basis_node_vanishing": max_node_residual < mp.mpf("1e-80"),
        "mass_term_obstruction_witness_exists": nonzero_mass_rows > 0,
        "derivative_only_quadratics_stationary_on_audited_tangent_basis": max_abs_derivative_residual < mp.mpf("1e-80") and max_abs_roughness_residual < mp.mpf("1e-80"),
    }


def coefficient_acceptance_grid() -> dict[str, Any]:
    witness = stationarity_witness_audit()
    m = mp.mpf(witness["max_abs_mass_variation_moment"])
    d1 = mp.mpf(witness["max_abs_first_derivative_variation"])
    d2 = mp.mpf(witness["max_abs_roughness_variation"])
    rows = []
    for c0, c1, c2 in [(0, 0, 1), (0, 1, 0), (0, 1, 1), (1, 0, 0), (1, 1, 1)]:
        upper_bound = abs(c0) * m + abs(c1) * d1 + abs(c2) * d2
        rows.append({
            "c0_mass": c0,
            "c1_first_derivative": c1,
            "c2_roughness": c2,
            "audited_stationarity_residual_bound": text(upper_bound, 70),
            "admitted_by_stationarity_audit": upper_bound < mp.mpf("1e-80"),
            "interpretation": "derivative-only ambiguity remains" if c0 == 0 else "mass term obstructed without extra forcing/counterterm",
        })
    return {
        "rows": rows,
        "derivative_only_rows_admitted": all(row["admitted_by_stationarity_audit"] for row in rows if row["c0_mass"] == 0),
        "mass_term_rows_rejected_without_extra_forcing": all(not row["admitted_by_stationarity_audit"] for row in rows if row["c0_mass"] != 0),
        "source_ambiguity_found": True,
        "source_acceptance_boundary": "A future strict source theorem must either export a derivative-only quadratic operator (or a compensating forcing/counterterm for any mass component) and must still choose the roughness order/coefficient from strict dynamics rather than from node fitting.",
    }


def build_source_admissibility_audit(p2511: dict[str, Any]) -> dict[str, Any]:
    symbolic = symbolic_first_variation_audit()
    witness = stationarity_witness_audit()
    grid = coefficient_acceptance_grid()
    return {
        "frontier_atom_under_attack": "strict_damping_beta_eta_source",
        "p2511_spline_collapse_inherited": p2511.get("natural_spline_collapse_theorem_for_postulated_functional") is True,
        "audit_type": "quadratic local source operator admissibility and obstruction audit for the conditional selector",
        "symbolic_first_variation_audit": symbolic,
        "stationarity_witness_audit": witness,
        "coefficient_acceptance_grid": grid,
        "mass_term_obstruction_theorem_exported_for_unforced_quadratic_sources": witness["mass_term_obstruction_witness_exists"] and grid["mass_term_rows_rejected_without_extra_forcing"],
        "derivative_only_source_ambiguity_identified": grid["derivative_only_rows_admitted"],
        "roughness_order_not_uniquely_sourced_by_stationarity": True,
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
## P2512/S1462 strict damping RG quadratic source admissibility audit

`P2512/S1462` audits what a future strict source theorem must still supply after the P2509--P2511 conditional selector chain.  For the local quadratic family `S[y]=1/2 int (c0*y^2+c1*(y')^2+c2*(y'')^2)d ell`, the first variation at `y0(ell)=delta ell` on node-vanishing perturbations is `c0 int y0*p + c1 delta (p(L)-p(0)) + c2 int y0''*p''`.  Because node-vanishing gives `p(0)=p(L)=0` and `y0''=0`, derivative-only terms are stationary, but an unforced mass term requires `int ell*p=0` for all node-vanishing `p`, which is false.  The finite polynomial witness family `R(ell)*ell^k` supplies explicit nonzero mass moments while preserving zero derivative-only residuals.

This narrows the source acceptance target but does not close it.  The audit identifies a real obstruction for unforced mass-like quadratic sources and a real ambiguity among derivative-only quadratic sources; it still does not derive the roughness action/order/coefficient from strict nadsoliton dynamics and exports no `strict_damping_beta_eta_source`, bridge theorem, role-transfer theorem, QW-2191 discharge, physical-value generator, or ToE closure.
"""
    lag_section = """
## P2512/S1462 quadratic source admissibility guard

`P2512/S1462` tests candidate local quadratic source operators for the postulated strict damping selector.  It shows that unforced `y^2` mass terms violate stationarity on explicit node-vanishing witnesses, while derivative-only quadratic terms remain stationary at `y0(ell)=delta ell`; hence the source problem is narrowed but not closed, because the strict dynamics still must choose the derivative order/coefficient and supply the action rather than assume it.  No nonlinear compression-flow source theorem or role-bearing `L_total` term is licensed.
"""
    append_once(DOC_FILES["equation_sheet"], "## P2512/S1462 strict damping RG quadratic source admissibility audit", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2512/S1462 quadratic source admissibility guard", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    p2511 = theorem(sources["P2511_NATURAL_SPLINE_COLLAPSE"], "strict_damping_rg_natural_spline_collapse_certificate")
    cert = build_source_admissibility_audit(p2511)
    witness = cert["stationarity_witness_audit"]
    grid = cert["coefficient_acceptance_grid"]
    theorem_export = {
        "theorem_name": "P2512_T1_strict_damping_rg_quadratic_source_admissibility_audit",
        "audited_chain": ["P2509/S1459", "P2510/S1460", "P2511/S1461"],
        "strict_damping_rg_quadratic_source_admissibility_audit": cert,
        "frontier_atom_under_attack": cert["frontier_atom_under_attack"],
        "p2511_spline_collapse_inherited": cert["p2511_spline_collapse_inherited"],
        "mass_term_obstruction_witness_exists": witness["mass_term_obstruction_witness_exists"],
        "derivative_only_quadratics_stationary_on_audited_tangent_basis": witness["derivative_only_quadratics_stationary_on_audited_tangent_basis"],
        "derivative_only_rows_admitted": grid["derivative_only_rows_admitted"],
        "mass_term_rows_rejected_without_extra_forcing": grid["mass_term_rows_rejected_without_extra_forcing"],
        "mass_term_obstruction_theorem_exported_for_unforced_quadratic_sources": cert["mass_term_obstruction_theorem_exported_for_unforced_quadratic_sources"],
        "derivative_only_source_ambiguity_identified": cert["derivative_only_source_ambiguity_identified"],
        "roughness_order_not_uniquely_sourced_by_stationarity": cert["roughness_order_not_uniquely_sourced_by_stationarity"],
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
            "P2512 supplies a source acceptance/obstruction audit, not a strict source derivation.",
            "Derivative-only quadratic stationarity is still ambiguous and does not choose the P2506 roughness order/coefficient from nadsoliton dynamics.",
            "No bridge theorem, role-transfer theorem, QW-2191 discharge, role-bearing L_total, physical-value generator, or ToE closure is exported.",
        ],
        "next_honest_step": grid["source_acceptance_boundary"],
    }
    gatekeepers = {
        "rg_audit_performed": grep["tool"] == "rg",
        "p2511_spline_collapse_inherited": theorem_export["p2511_spline_collapse_inherited"],
        "mass_obstruction_witness_exists": theorem_export["mass_term_obstruction_witness_exists"],
        "derivative_only_stationarity_confirmed": theorem_export["derivative_only_quadratics_stationary_on_audited_tangent_basis"],
        "source_ambiguity_not_hidden": theorem_export["derivative_only_source_ambiguity_identified"] and theorem_export["roughness_order_not_uniquely_sourced_by_stationarity"],
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
        "packet_id": "P2512",
        "stage_id": "S1462",
        "status": "STRICT_DAMPING_RG_QUADRATIC_SOURCE_ADMISSIBILITY_AUDIT_OBSTRUCTION_AND_AMBIGUITY_NO_SOURCE_EXPORT_NO_BRIDGE_THEOREM_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_damping_rg_quadratic_source_admissibility_audit": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_damping_rg_quadratic_source_admissibility_audit"]["theorem_export"]
    cert = t["strict_damping_rg_quadratic_source_admissibility_audit"]
    witness = cert["stationarity_witness_audit"]
    lines = [
        "# P2512/S1462 strict damping RG quadratic source admissibility audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        f"- Frontier atom under attack: `{t['frontier_atom_under_attack']}`.",
        f"- P2511 spline collapse inherited: `{t['p2511_spline_collapse_inherited']}`.",
        f"- Mass-term obstruction witness exists: `{t['mass_term_obstruction_witness_exists']}`.",
        f"- Nonzero mass obstruction rows: `{witness['nonzero_mass_obstruction_row_count']}`.",
        f"- Max mass variation moment: `{witness['max_abs_mass_variation_moment']}`.",
        f"- Max first-derivative variation: `{witness['max_abs_first_derivative_variation']}`.",
        f"- Max roughness variation: `{witness['max_abs_roughness_variation']}`.",
        f"- Derivative-only rows admitted: `{t['derivative_only_rows_admitted']}`.",
        f"- Mass-term rows rejected without extra forcing: `{t['mass_term_rows_rejected_without_extra_forcing']}`.",
        f"- Roughness order not uniquely sourced by stationarity: `{t['roughness_order_not_uniquely_sourced_by_stationarity']}`.",
        f"- Source theorem exported: `{t['strict_damping_beta_eta_source_exported']}`.",
        "",
        "## Negative controls",
        "",
        "This packet narrows the local quadratic source acceptance target but does not derive the strict damping source. It identifies an obstruction for unforced mass terms and an ambiguity among derivative-only quadratic sources; no strict source atom, bridge theorem, role-transfer theorem, selector/QW-2191 closure, role-bearing L_total term, physical-value generator, or ToE closure is exported.",
        "",
        "## Fingerprint",
        "",
        f"`{payload['strict_damping_rg_quadratic_source_admissibility_audit']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload["strict_damping_rg_quadratic_source_admissibility_audit"]["theorem_export"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
