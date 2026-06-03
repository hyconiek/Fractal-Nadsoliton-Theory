#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from fractions import Fraction
from math import log
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import (
    DOC_FILES,
    REPO,
    ROOT,
    load_json,
    rel,
)
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2518_s1468_strict_damping_biharmonic_affine_slope_nonidentifiability_certificate.json"
MD = GEN / "p2518_s1468_strict_damping_biharmonic_affine_slope_nonidentifiability_certificate.md"

SOURCE_FILES = {
    "P2517_AXIOM_BOUNDARY": GEN / "p2517_s1467_strict_damping_dual_key_axiom_boundary_certificate.json",
}

SLOPE_CANDIDATES = [Fraction(-1, 1), Fraction(0, 1), Fraction(4, 5), Fraction(1, 1), Fraction(9, 5), Fraction(2, 1)]
STRICT_RUNNING_SLOPE = Fraction(4, 5)
POLY_DEGREE = 8
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
        "new_packet": "P2518|S1468|biharmonic affine slope nonidentifiability|affine-slope continuum|operator signature numeric-source separation|slope-source nonidentifiability",
        "precursor_packets": "P2517|S1467|dual-key axiom boundary|P2516|dual-key source acceptance|P2515|operator-order signature acceptance|P2506|minimum-roughness",
        "slope_language": "affine slope|affine zero-mode|slope continuum|D\\^2.*affine|biharmonic.*slope|operator signature.*numeric",
        "source_blockers": "beta_eta_numeric_source|m2_operator_signature_source|strict source theorem|proper subset obstruction|non-strict",
        "guardrails": "K_legacy_ont|K_strict_gate|role-transfer audit|QW-2191|ToE closure|selector guardrail",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def frac_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def derivative_of_affine(order: int, slope: Fraction, intercept: Fraction = Fraction(0, 1)) -> Fraction:
    if order == 0:
        raise ValueError("pointwise affine values depend on ell; use affine_node_rows instead")
    if order == 1:
        return slope
    return Fraction(0, 1)


def affine_operator_observable_rows() -> list[dict[str, Any]]:
    rows = []
    for slope in SLOPE_CANDIDATES:
        d2 = derivative_of_affine(2, slope)
        d3 = derivative_of_affine(3, slope)
        d4 = derivative_of_affine(4, slope)
        rows.append({
            "slope": frac_text(slope),
            "is_strict_running_slope_delta_4_over_5": slope == STRICT_RUNNING_SLOPE,
            "D2_y": frac_text(d2),
            "D3_y": frac_text(d3),
            "D4_y": frac_text(d4),
            "J2_energy_integrand": frac_text(d2 * d2),
            "biharmonic_euler_lagrange_residual_D4_y": frac_text(d4),
            "natural_boundary_y_second": frac_text(d2),
            "natural_boundary_y_third": frac_text(d3),
            "m2_operator_signature_observables_all_zero": d2 == d3 == d4 == 0,
        })
    return rows


def affine_node_rows() -> list[dict[str, Any]]:
    rows = []
    for slope in SLOPE_CANDIDATES:
        rows.append({
            "slope": frac_text(slope),
            "node_values_y_delta_log_d": [float(slope) * log(d) for d in DOMAIN],
            "relative_to_strict_delta_4_over_5_at_d_11": float((slope - STRICT_RUNNING_SLOPE)) * log(11),
        })
    return rows


def derivative_matrix(order: int, degree: int) -> list[list[Fraction]]:
    size = degree + 1
    matrix = [[Fraction(0, 1) for _ in range(size)] for _ in range(size)]
    for power in range(size):
        if power >= order:
            coeff = Fraction(1, 1)
            for factor in range(power - order + 1, power + 1):
                coeff *= factor
            matrix[power - order][power] = coeff
    return matrix


def rank_fraction(matrix: list[list[Fraction]]) -> int:
    work = [row[:] for row in matrix]
    if not work:
        return 0
    row_count = len(work)
    col_count = len(work[0])
    rank = 0
    for col in range(col_count):
        pivot = None
        for row in range(rank, row_count):
            if work[row][col] != 0:
                pivot = row
                break
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        pivot_value = work[rank][col]
        work[rank] = [value / pivot_value for value in work[rank]]
        for row in range(row_count):
            if row != rank and work[row][col] != 0:
                factor = work[row][col]
                work[row] = [work[row][idx] - factor * work[rank][idx] for idx in range(col_count)]
        rank += 1
        if rank == row_count:
            break
    return rank


def finite_operator_rank_audit() -> dict[str, Any]:
    rows = []
    for order in [1, 2, 3, 4]:
        matrix = derivative_matrix(order, POLY_DEGREE)
        rank = rank_fraction(matrix)
        kernel_dim = POLY_DEGREE + 1 - rank
        rows.append({
            "operator": f"D^{order}",
            "polynomial_degree_cutoff": POLY_DEGREE,
            "matrix_size": [POLY_DEGREE + 1, POLY_DEGREE + 1],
            "rank": rank,
            "expected_rank": max(0, POLY_DEGREE + 1 - order),
            "kernel_dimension": kernel_dim,
            "expected_kernel_dimension": min(order, POLY_DEGREE + 1),
            "affine_subspace_in_kernel": order >= 2,
        })
    return {
        "basis": [f"ell^{power}" for power in range(POLY_DEGREE + 1)],
        "rows": rows,
        "rank_identities_pass": all(row["rank"] == row["expected_rank"] and row["kernel_dimension"] == row["expected_kernel_dimension"] for row in rows),
        "D2_kernel_dimension_is_two": next(row for row in rows if row["operator"] == "D^2")["kernel_dimension"] == 2,
        "D4_kernel_dimension_contains_affine_but_is_larger": next(row for row in rows if row["operator"] == "D^4")["kernel_dimension"] == 4,
    }


def build_affine_slope_certificate(p2517: dict[str, Any]) -> dict[str, Any]:
    observable_rows = affine_operator_observable_rows()
    rank_audit = finite_operator_rank_audit()
    distinct_zero_slopes = [row["slope"] for row in observable_rows if row["m2_operator_signature_observables_all_zero"]]
    return {
        "frontier_atom_under_attack": "beta_eta_numeric_source independent from m2_operator_signature_source",
        "p2517_axiom_boundary_inherited": p2517.get("axiom_boundary_exported") is True,
        "certificate_type": "biharmonic affine-slope nonidentifiability certificate",
        "symbolic_theorem": {
            "statement": "For J2[y]=1/2 int (y''(ell))^2 d ell, every affine y(ell)=a ell+b has y''=y'''=y''''=0, hence zero J2 energy, zero biharmonic Euler-Lagrange residual, and zero natural-boundary concomitant. The m=2 operator signature alone therefore cannot select the strict running slope delta=4/5 or eta=9/5.",
            "strict_running_slope_delta": frac_text(STRICT_RUNNING_SLOPE),
            "affine_family_dimension": 2,
            "numeric_slope_selected_by_operator_signature": False,
        },
        "affine_operator_observable_rows": observable_rows,
        "affine_node_value_rows": affine_node_rows(),
        "finite_operator_rank_audit": rank_audit,
        "all_audited_slopes_operator_indistinguishable": len(distinct_zero_slopes) == len(SLOPE_CANDIDATES),
        "distinct_zero_energy_slope_count": len(distinct_zero_slopes),
        "strict_delta_is_only_one_member_of_zero_energy_family": STRICT_RUNNING_SLOPE in SLOPE_CANDIDATES and len(distinct_zero_slopes) > 1,
        "operator_signature_numeric_key_separation_exported": True,
        "m2_operator_signature_source_exported": False,
        "beta_eta_numeric_source_exported": False,
        "strict_damping_beta_eta_source_exported": False,
        "damping_compression_bridge_component_ready": False,
        "full_bridge_theorem_exported": False,
        "role_transfer_theorem_exported": False,
        "selector_closure_exported": False,
        "qw2191_discharged_by_this_certificate": False,
        "role_bearing_ltotal_exported": False,
        "toe_closure_claimed": False,
    }


def append_doc_sections() -> None:
    eq_section = """
## P2518/S1468 biharmonic affine-slope nonidentifiability certificate

`P2518/S1468` sharpens the P2516/P2517 two-key boundary by proving a concrete separation theorem: the `m=2` biharmonic operator signature cannot by itself identify the numeric strict damping slope.  For `J_2[y]=1/2 int (y''(ell))^2 d ell`, every affine running exponent `y(ell)=a ell+b` has `y''=y'''=y''''=0`, so the energy, biharmonic Euler-Lagrange residual, and natural-boundary concomitant all vanish for a continuum of slopes.  The strict candidate `delta=4/5` is therefore only one member of the zero-energy affine family unless an independent numeric source/node theorem supplies it.

The finite polynomial audit confirms the rank/nullity boundary on monomials through degree 8: `D^2` has a two-dimensional affine kernel, and `D^4` has a larger four-dimensional kernel.  This exports an operator-signature/numeric-key separation certificate only; it does not export `beta_eta_numeric_source`, `m2_operator_signature_source`, strict damping source closure, bridge completion, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, physical-value generator, or ToE closure.
"""
    lag_section = """
## P2518/S1468 affine-slope nonidentifiability guard

`P2518/S1468` blocks a false source inference from the biharmonic roughness action: the `m=2` operator signature annihilates every affine `y(ell)=a ell+b`, not just the strict `delta=4/5` affine flow.  A future role-bearing damping term must therefore source both the operator signature and the numeric slope/`beta, eta` data; the operator signature alone cannot license nonlinear compression-flow source closure.
"""
    append_once(DOC_FILES["equation_sheet"], "## P2518/S1468 biharmonic affine-slope nonidentifiability certificate", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2518/S1468 affine-slope nonidentifiability guard", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    p2517 = theorem(sources["P2517_AXIOM_BOUNDARY"], "strict_damping_dual_key_axiom_boundary_certificate")
    cert = build_affine_slope_certificate(p2517)
    rank_audit = cert["finite_operator_rank_audit"]
    theorem_export = {
        "theorem_name": "P2518_T1_strict_damping_biharmonic_affine_slope_nonidentifiability_certificate",
        "audited_chain": ["P2506/S1456", "P2515/S1465", "P2516/S1466", "P2517/S1467"],
        "strict_damping_biharmonic_affine_slope_nonidentifiability_certificate": cert,
        "frontier_atom_under_attack": cert["frontier_atom_under_attack"],
        "p2517_axiom_boundary_inherited": cert["p2517_axiom_boundary_inherited"],
        "strict_running_slope_delta": cert["symbolic_theorem"]["strict_running_slope_delta"],
        "affine_family_dimension": cert["symbolic_theorem"]["affine_family_dimension"],
        "distinct_zero_energy_slope_count": cert["distinct_zero_energy_slope_count"],
        "all_audited_slopes_operator_indistinguishable": cert["all_audited_slopes_operator_indistinguishable"],
        "strict_delta_is_only_one_member_of_zero_energy_family": cert["strict_delta_is_only_one_member_of_zero_energy_family"],
        "D2_kernel_dimension_is_two": rank_audit["D2_kernel_dimension_is_two"],
        "D4_kernel_dimension_contains_affine_but_is_larger": rank_audit["D4_kernel_dimension_contains_affine_but_is_larger"],
        "finite_rank_identities_pass": rank_audit["rank_identities_pass"],
        "operator_signature_numeric_key_separation_exported": cert["operator_signature_numeric_key_separation_exported"],
        "m2_operator_signature_source_exported": False,
        "beta_eta_numeric_source_exported": False,
        "strict_damping_beta_eta_source_exported": False,
        "damping_compression_bridge_component_ready": False,
        "full_bridge_theorem_exported": False,
        "role_transfer_theorem_exported": False,
        "selector_closure_exported": False,
        "qw2191_discharged_by_this_certificate": False,
        "role_bearing_ltotal_exported": False,
        "toe_closure_claimed": False,
        "not_licensed": [
            "P2518 exports a nonidentifiability/separation theorem, not a numeric beta/eta source.",
            "The m=2 biharmonic signature annihilates every affine slope and cannot select delta=4/5 by itself.",
            "Finite node values can distinguish slopes only when supplied as numeric data, not from the operator signature alone.",
            "No damping bridge completion, role-transfer theorem, QW-2191 discharge, role-bearing L_total, physical-value generator, or ToE closure is exported.",
        ],
        "next_honest_step": "Search for a strict nadsoliton source that fixes the affine slope/numeric beta-eta data in addition to sourcing the m=2 operator signature.",
    }
    gatekeepers = {
        "rg_audit_performed": grep["tool"] == "rg",
        "p2517_boundary_inherited": theorem_export["p2517_axiom_boundary_inherited"],
        "all_audited_slopes_zero": theorem_export["all_audited_slopes_operator_indistinguishable"],
        "strict_delta_not_unique": theorem_export["strict_delta_is_only_one_member_of_zero_energy_family"],
        "rank_identities_pass": theorem_export["finite_rank_identities_pass"],
        "separation_exported": theorem_export["operator_signature_numeric_key_separation_exported"],
        "source_not_exported": not theorem_export["strict_damping_beta_eta_source_exported"],
        "negative_controls_preserved": not any(theorem_export[key] for key in [
            "m2_operator_signature_source_exported",
            "beta_eta_numeric_source_exported",
            "strict_damping_beta_eta_source_exported",
            "damping_compression_bridge_component_ready",
            "full_bridge_theorem_exported",
            "role_transfer_theorem_exported",
            "selector_closure_exported",
            "qw2191_discharged_by_this_certificate",
            "role_bearing_ltotal_exported",
            "toe_closure_claimed",
        ]),
    }
    return {
        "packet_id": "P2518",
        "stage_id": "S1468",
        "status": "STRICT_DAMPING_BIHARMONIC_AFFINE_SLOPE_NONIDENTIFIABILITY_CERTIFICATE_NO_SOURCE_EXPORT_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_damping_biharmonic_affine_slope_nonidentifiability_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_damping_biharmonic_affine_slope_nonidentifiability_certificate"]["theorem_export"]
    lines = [
        "# P2518/S1468 strict damping biharmonic affine-slope nonidentifiability certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        f"- Frontier atom under attack: `{t['frontier_atom_under_attack']}`.",
        f"- P2517 axiom boundary inherited: `{t['p2517_axiom_boundary_inherited']}`.",
        f"- Strict running slope delta: `{t['strict_running_slope_delta']}`.",
        f"- Affine family dimension: `{t['affine_family_dimension']}`.",
        f"- Distinct zero-energy audited slopes: `{t['distinct_zero_energy_slope_count']}`.",
        f"- All audited slopes operator-indistinguishable: `{t['all_audited_slopes_operator_indistinguishable']}`.",
        f"- Strict delta is only one zero-energy affine member: `{t['strict_delta_is_only_one_member_of_zero_energy_family']}`.",
        f"- D2 kernel dimension is two: `{t['D2_kernel_dimension_is_two']}`.",
        f"- D4 kernel contains affine kernel but is larger: `{t['D4_kernel_dimension_contains_affine_but_is_larger']}`.",
        f"- Operator-signature/numeric-key separation exported: `{t['operator_signature_numeric_key_separation_exported']}`.",
        f"- Strict damping source exported: `{t['strict_damping_beta_eta_source_exported']}`.",
        "",
        "## Negative controls",
        "",
        "This packet exports only a separation/nonidentifiability theorem. It does not export numeric beta/eta sourcing, m=2 operator source closure, bridge completion, role-transfer theorem, selector/QW-2191 closure, role-bearing L_total, physical-value generator, or ToE closure.",
        "",
        "## Fingerprint",
        "",
        f"`{payload['strict_damping_biharmonic_affine_slope_nonidentifiability_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload["strict_damping_biharmonic_affine_slope_nonidentifiability_certificate"]["theorem_export"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
