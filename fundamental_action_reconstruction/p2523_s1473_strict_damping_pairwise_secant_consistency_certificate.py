#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
import subprocess
from itertools import combinations
from statistics import mean, pstdev
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
OUT = GEN / "p2523_s1473_strict_damping_pairwise_secant_consistency_certificate.json"
MD = GEN / "p2523_s1473_strict_damping_pairwise_secant_consistency_certificate.md"

SOURCE_FILES = {
    "P2522_TWO_NODE_ANCHOR_BASIS": GEN / "p2522_s1472_strict_damping_two_node_anchor_basis_equivalence_certificate.json",
}

STRICT_DELTA = 4.0 / 5.0
STRICT_ETA = 9.0 / 5.0
NODE_DOMAIN = list(range(1, 12))
TOL = 1e-14


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
        "new_packet": "P2523|S1473|pairwise secant consistency|secant slope consistency|triangle secant cocycle|basis-independent affine invariant|affine data consistency",
        "precursor_packets": "P2522|S1472|two-node anchor basis equivalence|P2521|single-node anchor equivalence|P2520|endpoint anchor subkey lattice",
        "secant_language": "pairwise secant|secant slope|triangle cocycle|affine consistency|basis-independent|overdetermined node consistency",
        "source_blockers": "beta_eta_numeric_source|node value source|anchor basis source|strict source theorem|proper subset obstruction",
        "guardrails": "K_legacy_ont|K_strict_gate|role-transfer audit|QW-2191|ToE closure|selector guardrail",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def ell(d: int) -> float:
    return math.log(d)


def y_value(d: int) -> float:
    return STRICT_DELTA * ell(d)


def pairwise_secant_rows() -> list[dict[str, Any]]:
    rows = []
    for d_i, d_j in combinations(NODE_DOMAIN, 2):
        ell_i = ell(d_i)
        ell_j = ell(d_j)
        y_i = y_value(d_i)
        y_j = y_value(d_j)
        secant = (y_j - y_i) / (ell_j - ell_i)
        rows.append({
            "node_pair": [d_i, d_j],
            "ell_gap": ell_j - ell_i,
            "y_gap": y_j - y_i,
            "secant_delta": secant,
            "secant_delta_matches_4_over_5": abs(secant - STRICT_DELTA) < TOL,
            "eta_from_secant": 1.0 + secant,
            "eta_matches_9_over_5": abs((1.0 + secant) - STRICT_ETA) < TOL,
        })
    return rows


def triangle_cocycle_rows() -> list[dict[str, Any]]:
    rows = []
    for d_i, d_j, d_k in combinations(NODE_DOMAIN, 3):
        y_i = y_value(d_i)
        y_j = y_value(d_j)
        y_k = y_value(d_k)
        additive_cocycle = (y_j - y_i) + (y_k - y_j) - (y_k - y_i)
        secant_ij = (y_j - y_i) / (ell(d_j) - ell(d_i))
        secant_jk = (y_k - y_j) / (ell(d_k) - ell(d_j))
        secant_ik = (y_k - y_i) / (ell(d_k) - ell(d_i))
        rows.append({
            "node_triple": [d_i, d_j, d_k],
            "additive_y_cocycle": additive_cocycle,
            "additive_y_cocycle_zero": abs(additive_cocycle) < TOL,
            "secant_delta_ij": secant_ij,
            "secant_delta_jk": secant_jk,
            "secant_delta_ik": secant_ik,
            "secant_spread": max(secant_ij, secant_jk, secant_ik) - min(secant_ij, secant_jk, secant_ik),
            "all_three_secants_match_delta": all(abs(value - STRICT_DELTA) < TOL for value in [secant_ij, secant_jk, secant_ik]),
        })
    return rows


def affine_design_projection_audit() -> dict[str, Any]:
    xs = [ell(d) for d in NODE_DOMAIN]
    ys = [y_value(d) for d in NODE_DOMAIN]
    n = len(xs)
    sum_x = sum(xs)
    sum_x2 = sum(x * x for x in xs)
    sum_y = sum(ys)
    sum_xy = sum(x * y for x, y in zip(xs, ys))
    determinant = n * sum_x2 - sum_x * sum_x
    intercept = (sum_x2 * sum_y - sum_x * sum_xy) / determinant
    slope = (n * sum_xy - sum_x * sum_y) / determinant
    residuals = [y - (intercept + slope * x) for x, y in zip(xs, ys)]
    return {
        "design": "columns [1, log(d)] over d=1..11",
        "normal_matrix": [[n, sum_x], [sum_x, sum_x2]],
        "normal_matrix_determinant": determinant,
        "normal_matrix_determinant_positive": determinant > 0,
        "projected_intercept_log_beta": intercept,
        "projected_intercept_matches_zero": abs(intercept) < TOL,
        "projected_delta": slope,
        "projected_delta_matches_4_over_5": abs(slope - STRICT_DELTA) < TOL,
        "projected_eta": 1.0 + slope,
        "projected_eta_matches_9_over_5": abs((1.0 + slope) - STRICT_ETA) < TOL,
        "max_abs_projection_residual": max(abs(value) for value in residuals),
        "projection_residual_zero": max(abs(value) for value in residuals) < TOL,
    }


def build_secant_consistency_certificate(p2522: dict[str, Any]) -> dict[str, Any]:
    pair_rows = pairwise_secant_rows()
    triangle_rows = triangle_cocycle_rows()
    design = affine_design_projection_audit()
    secants = [row["secant_delta"] for row in pair_rows]
    return {
        "frontier_atom_under_attack": "basis-independent numeric node-data consistency for beta_eta target",
        "p2522_anchor_basis_equivalence_inherited": p2522.get("anchor_basis_equivalence_exported") is True,
        "certificate_type": "pairwise secant consistency and affine-data cocycle certificate",
        "symbolic_statement": "If strict node data are y_d=(4/5)log d on d=1..11, then every pairwise secant slope equals 4/5 and every triangle additive y-cocycle vanishes. This provides a basis-independent consistency check for the numeric target, but it does not source the node data.",
        "pairwise_secant_rows": pair_rows,
        "triangle_cocycle_rows": triangle_rows,
        "affine_design_projection_audit": design,
        "pair_count": len(pair_rows),
        "triangle_count": len(triangle_rows),
        "all_pairwise_secants_match_delta": all(row["secant_delta_matches_4_over_5"] for row in pair_rows),
        "all_pairwise_etas_match_9_over_5": all(row["eta_matches_9_over_5"] for row in pair_rows),
        "secant_delta_min": min(secants),
        "secant_delta_max": max(secants),
        "secant_delta_mean": mean(secants),
        "secant_delta_population_std": pstdev(secants),
        "all_triangle_additive_cocycles_zero": all(row["additive_y_cocycle_zero"] for row in triangle_rows),
        "all_triangle_secants_match_delta": all(row["all_three_secants_match_delta"] for row in triangle_rows),
        "basis_independent_affine_consistency_exported": True,
        "node_data_source_exported": False,
        "anchor_basis_source_exported": False,
        "beta_eta_numeric_source_exported": False,
        "m2_operator_signature_source_exported": False,
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
## P2523/S1473 pairwise secant consistency certificate

`P2523/S1473` turns the P2522 node-pair basis equivalence into a basis-independent finite consistency audit.  For the strict target node data `y_d=(4/5)log d` on `d=1..11`, all 55 pairwise secants have slope `4/5`, all 165 triangle additive `y`-cocycles vanish, and the affine projection on columns `[1, log(d)]` recovers intercept `0`, slope `4/5`, hence `beta=1`, `eta=9/5`, with zero residual within arithmetic tolerance.

This is a consistency certificate for already supplied node data, not a source theorem for those data.  It does not export `beta_eta_numeric_source`, does not source the `m=2` operator signature, and exports no strict damping source closure, bridge completion, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, physical-value generator, or ToE closure.
"""
    lag_section = """
## P2523/S1473 pairwise secant consistency guard

`P2523/S1473` gives a basis-independent check that all strict node values are mutually affine with slope `4/5`, but it assumes the node data being checked.  Therefore the calculation still does not source the numeric damping data and does not license a nonlinear compression-flow source or role-bearing `L_total` term.
"""
    append_once(DOC_FILES["equation_sheet"], "## P2523/S1473 pairwise secant consistency certificate", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2523/S1473 pairwise secant consistency guard", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    p2522 = theorem(sources["P2522_TWO_NODE_ANCHOR_BASIS"], "strict_damping_two_node_anchor_basis_equivalence_certificate")
    cert = build_secant_consistency_certificate(p2522)
    design = cert["affine_design_projection_audit"]
    theorem_export = {
        "theorem_name": "P2523_T1_strict_damping_pairwise_secant_consistency_certificate",
        "audited_chain": ["P2522/S1472", "P2521/S1471", "P2520/S1470"],
        "strict_damping_pairwise_secant_consistency_certificate": cert,
        "frontier_atom_under_attack": cert["frontier_atom_under_attack"],
        "p2522_anchor_basis_equivalence_inherited": cert["p2522_anchor_basis_equivalence_inherited"],
        "pair_count": cert["pair_count"],
        "triangle_count": cert["triangle_count"],
        "all_pairwise_secants_match_delta": cert["all_pairwise_secants_match_delta"],
        "all_pairwise_etas_match_9_over_5": cert["all_pairwise_etas_match_9_over_5"],
        "all_triangle_additive_cocycles_zero": cert["all_triangle_additive_cocycles_zero"],
        "all_triangle_secants_match_delta": cert["all_triangle_secants_match_delta"],
        "projected_intercept_matches_zero": design["projected_intercept_matches_zero"],
        "projected_delta_matches_4_over_5": design["projected_delta_matches_4_over_5"],
        "projected_eta_matches_9_over_5": design["projected_eta_matches_9_over_5"],
        "projection_residual_zero": design["projection_residual_zero"],
        "basis_independent_affine_consistency_exported": cert["basis_independent_affine_consistency_exported"],
        "node_data_source_exported": False,
        "anchor_basis_source_exported": False,
        "beta_eta_numeric_source_exported": False,
        "m2_operator_signature_source_exported": False,
        "strict_damping_beta_eta_source_exported": False,
        "damping_compression_bridge_component_ready": False,
        "full_bridge_theorem_exported": False,
        "role_transfer_theorem_exported": False,
        "selector_closure_exported": False,
        "qw2191_discharged_by_this_certificate": False,
        "role_bearing_ltotal_exported": False,
        "toe_closure_claimed": False,
        "not_licensed": [
            "P2523 exports basis-independent consistency of supplied node data, not a source for those data.",
            "The pairwise secant and triangle cocycle checks cannot replace a strict node-data source theorem.",
            "The m=2 operator signature source remains separate and unsourced by this packet.",
            "No damping bridge completion, role-transfer theorem, QW-2191 discharge, role-bearing L_total, physical-value generator, or ToE closure is exported.",
        ],
        "next_honest_step": "Search for a strict mechanism that generates the affine node-data vector itself, rather than only checking its pairwise consistency after it is supplied.",
    }
    gatekeepers = {
        "rg_audit_performed": grep["tool"] == "rg",
        "p2522_inherited": theorem_export["p2522_anchor_basis_equivalence_inherited"],
        "all_55_pairs_checked": theorem_export["pair_count"] == 55,
        "all_165_triangles_checked": theorem_export["triangle_count"] == 165,
        "secants_and_cocycles_pass": theorem_export["all_pairwise_secants_match_delta"] and theorem_export["all_triangle_additive_cocycles_zero"],
        "projection_passes": theorem_export["projected_intercept_matches_zero"] and theorem_export["projected_delta_matches_4_over_5"] and theorem_export["projection_residual_zero"],
        "negative_controls_preserved": not any(theorem_export[key] for key in [
            "node_data_source_exported",
            "anchor_basis_source_exported",
            "beta_eta_numeric_source_exported",
            "m2_operator_signature_source_exported",
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
        "packet_id": "P2523",
        "stage_id": "S1473",
        "status": "STRICT_DAMPING_PAIRWISE_SECANT_CONSISTENCY_CERTIFICATE_BASIS_INDEPENDENT_CHECK_NO_SOURCE_EXPORT_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_damping_pairwise_secant_consistency_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_damping_pairwise_secant_consistency_certificate"]["theorem_export"]
    lines = [
        "# P2523/S1473 strict damping pairwise secant consistency certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        f"- Frontier atom under attack: `{t['frontier_atom_under_attack']}`.",
        f"- P2522 anchor-basis equivalence inherited: `{t['p2522_anchor_basis_equivalence_inherited']}`.",
        f"- Pair count: `{t['pair_count']}`.",
        f"- Triangle count: `{t['triangle_count']}`.",
        f"- All pairwise secants match delta: `{t['all_pairwise_secants_match_delta']}`.",
        f"- All triangle additive cocycles zero: `{t['all_triangle_additive_cocycles_zero']}`.",
        f"- Projected delta matches 4/5: `{t['projected_delta_matches_4_over_5']}`.",
        f"- Projection residual zero: `{t['projection_residual_zero']}`.",
        f"- Basis-independent affine consistency exported: `{t['basis_independent_affine_consistency_exported']}`.",
        f"- Strict damping source exported: `{t['strict_damping_beta_eta_source_exported']}`.",
        "",
        "## Negative controls",
        "",
        "This packet exports only consistency of supplied node data. It does not source node data, anchor basis, m=2 operator signature, bridge completion, role-transfer theorem, selector/QW-2191 closure, role-bearing L_total, physical-value generator, or ToE closure.",
        "",
        "## Fingerprint",
        "",
        f"`{payload['strict_damping_pairwise_secant_consistency_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload["strict_damping_pairwise_secant_consistency_certificate"]["theorem_export"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
