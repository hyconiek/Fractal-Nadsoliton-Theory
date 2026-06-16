#!/usr/bin/env python3
"""P2785/S1735: exact characteristic-polynomial certificate.

P2784 separated the seven P2781/P2783 local quotient representatives by three
floating-point spectral probes.  This follow-up removes the numerical-eigenvalue
weakness: it recomputes the same seven representatives and uses exact integer
characteristic polynomials for the adjacency, Laplacian, and signless-Laplacian
matrices.  It also computes Kirchhoff spanning-tree counts by an exact integer
cofactor determinant.

The result is a local exact-algebra certificate only.  It does not expand the
graph class, certify a full 16-node 4-regular generator, or export any strict
spectral action/source law coupled to K/L_total.
"""
from __future__ import annotations

import hashlib
import itertools
import json
from pathlib import Path
from typing import Any

import sympy as sp

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2778_s1728_max_symmetry_16node_geometry_source_audit import N, candidate_edge_sets
from p2779_s1729_16node_circulant_full_spectrum_quotient_audit import quotient_classes
from p2781_s1731_enumerated_two_layer_c8_spectrum_collision_audit import all_two_layer_c8_candidates

GEN = ROOT / "generated"
P2784 = GEN / "p2784_s1734_seven_class_multispectral_robustness_audit.json"
OUT = GEN / "p2785_s1735_exact_characteristic_polynomial_certificate.json"
MD = GEN / "p2785_s1735_exact_characteristic_polynomial_certificate.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

NEGATIVE_EXPORT_FLAGS = [
    "canonical_geometry_source_exported",
    "strict_spectral_source_law_exported",
    "global_full_spectrum_geometry_theorem_exported",
    "kernel_geometry_closure_exported",
    "kernel_fully_expresses_nadsoliton_characteristics",
    "role_bearing_ltotal_promoted",
    "bridge_closure_exported",
    "selector_closure_exported",
    "toe_closure_exported",
]


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"missing": True, "path": rel(path)}


def adjacency_matrix(edges: list[tuple[int, int]]) -> sp.Matrix:
    mat = sp.zeros(N, N)
    for a, b in edges:
        mat[a, b] = 1
        mat[b, a] = 1
    return mat


def charpoly_coefficients(mat: sp.Matrix) -> list[int]:
    return [int(coeff) for coeff in mat.charpoly().all_coeffs()]


def exact_row(rep: dict[str, Any], by_name: dict[str, dict[str, Any]]) -> dict[str, Any]:
    edges = by_name[rep["representative"]]["edges"]
    adj = adjacency_matrix(edges)
    deg = sp.diag(*[sum(adj[row, col] for col in range(N)) for row in range(N)])
    lap = deg - adj
    signless = deg + adj
    cofactor = lap[1:, 1:]
    return {
        "representative": rep["representative"],
        "member_count": len(rep["members"]),
        "members": rep["members"],
        "adjacency_charpoly_coefficients": charpoly_coefficients(adj),
        "laplacian_charpoly_coefficients": charpoly_coefficients(lap),
        "signless_laplacian_charpoly_coefficients": charpoly_coefficients(signless),
        "kirchhoff_spanning_tree_count_exact": int(cofactor.det()),
    }


def exact_polynomial_witness() -> dict[str, Any]:
    candidates = candidate_edge_sets() + all_two_layer_c8_candidates()
    by_name = {row["geometry"]: row for row in candidates}
    classes = quotient_classes(candidates)
    rows = [exact_row(rep, by_name) for rep in classes]
    invariant_keys = [
        "adjacency_charpoly_coefficients",
        "laplacian_charpoly_coefficients",
        "signless_laplacian_charpoly_coefficients",
        "kirchhoff_spanning_tree_count_exact",
    ]
    pair_rows: list[dict[str, Any]] = []
    for left, right in itertools.combinations(rows, 2):
        collisions = {key: left[key] == right[key] for key in invariant_keys}
        pair_rows.append({
            "left": left["representative"],
            "right": right["representative"],
            "exact_invariant_collisions": collisions,
            "all_three_charpolys_distinct": not collisions["adjacency_charpoly_coefficients"] and not collisions["laplacian_charpoly_coefficients"] and not collisions["signless_laplacian_charpoly_coefficients"],
        })
    collision_counts = {key: sum(1 for row in pair_rows if row["exact_invariant_collisions"][key]) for key in invariant_keys}
    return {
        "source_class": "Same seven P2781/P2783/P2784 local quotient representatives; no graph-class expansion is introduced.",
        "representative_count": len(rows),
        "pair_count": len(pair_rows),
        "exact_rows": rows,
        "pair_rows": pair_rows,
        "collision_counts_by_exact_invariant": collision_counts,
        "exact_charpoly_collision_counts": {key: collision_counts[key] for key in invariant_keys[:3]},
        "all_pairs_separated_by_all_three_exact_charpolys": all(row["all_three_charpolys_distinct"] for row in pair_rows),
        "finite_certificate_statement": "All 21 representative pairs are separated by exact integer adjacency, Laplacian, and signless-Laplacian characteristic polynomials.",
    }


def acceptance_matrix(witness: dict[str, Any], p2784: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "p2784_multispectral_certificate_present": p2784.get("status") == "P2784_SEVEN_CLASS_MULTISPECTRAL_ROBUSTNESS_AUDIT_NO_CLOSURE",
        "seven_representatives_reused_without_class_expansion": witness["representative_count"] == 7,
        "all_21_pairs_checked_exactly": witness["pair_count"] == 21,
        "zero_exact_adjacency_charpoly_collisions": witness["collision_counts_by_exact_invariant"]["adjacency_charpoly_coefficients"] == 0,
        "zero_exact_laplacian_charpoly_collisions": witness["collision_counts_by_exact_invariant"]["laplacian_charpoly_coefficients"] == 0,
        "zero_exact_signless_charpoly_collisions": witness["collision_counts_by_exact_invariant"]["signless_laplacian_charpoly_coefficients"] == 0,
        "strict_nadsoliton_spectral_source_law_exported": False,
        "canonical_full_16node_4regular_generator_supplied": False,
        "kernel_or_ltotal_variational_coupling_exported": False,
    }
    return {
        "facts": facts,
        "accepted_as_exact_local_algebra_certificate": all(facts[key] for key in [
            "p2784_multispectral_certificate_present",
            "seven_representatives_reused_without_class_expansion",
            "all_21_pairs_checked_exactly",
            "zero_exact_adjacency_charpoly_collisions",
            "zero_exact_laplacian_charpoly_collisions",
            "zero_exact_signless_charpoly_collisions",
        ]),
        "accepted_as_strict_spectral_source_law": False,
        "accepted_as_canonical_nadsoliton_geometry_source": False,
        "accepted_as_ltotal_or_toe_promotion": False,
        "missing_criteria": [key for key, value in facts.items() if not value],
        "blocker": "Exact integer characteristic polynomials remove numerical spectral ambiguity inside the seven-class local quotient, but no canonical full graph generator, strict source law, or K/L_total variational coupling is exported.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    w = payload["exact_polynomial_witness"]
    lines = [
        "# P2785/S1735 exact characteristic-polynomial certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Exact algebra result",
        f"- representative_count={w['representative_count']}",
        f"- pair_count={w['pair_count']}",
        f"- exact_charpoly_collision_counts={w['exact_charpoly_collision_counts']}",
        f"- all_pairs_separated_by_all_three_exact_charpolys={w['all_pairs_separated_by_all_three_exact_charpolys']}",
        "",
        "## Decision",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    p2784 = read_json(P2784)
    witness = exact_polynomial_witness()
    acceptance = acceptance_matrix(witness, p2784)
    payload = {
        "status": "P2785_EXACT_CHARACTERISTIC_POLYNOMIAL_CERTIFICATE_NO_CLOSURE",
        "input_hashes": {"P2784": sha(P2784)},
        "input_statuses": {"P2784": p2784.get("status")},
        "audited_question": "Does the P2784 multispectral separation survive exact integer characteristic-polynomial computation?",
        "exact_polynomial_witness": witness,
        "acceptance_matrix": acceptance,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": acceptance["blocker"],
            "next_honest_step": "Use P2785 as an exact local algebra certificate only.  The next honest move is exactly one of: supply a canonical-generation theorem/tool certificate for the full connected 16-node 4-regular graph class with reproducible graph6/hash provenance and then rerun exact-polynomial collision auditing; or export a strict nadsoliton spectral action/source law fixing the admissible class, target spectrum, and K/L_total coupling before testing.  Otherwise preserve the P2697-P2785 no-canonical-geometry/no-closure certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2785/S1735 exact characteristic-polynomial certificate", "## P2785/S1735 exact characteristic-polynomial certificate\n\n`P2785/S1735` removes the floating-eigenvalue ambiguity in P2784 by computing exact integer characteristic polynomials for adjacency, Laplacian, and signless-Laplacian matrices on the same seven local representatives.  All 21 pairs are separated by all three exact characteristic polynomial families, and exact Kirchhoff cofactor determinants are recorded.  This is still only a local exact-algebra certificate: it does not expand the graph class, certify a canonical full 16-node 4-regular generator, source the target spectrum from `K`/`L_total`, or export canonical geometry, bridge closure, role transfer, selector closure, role-bearing `L_total`, or ToE closure.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2785/S1735 exact-polynomial Ltotal guard", "## P2785/S1735 exact-polynomial Ltotal guard\n\n`P2785/S1735` adds no variational source term.  Exact characteristic polynomials strengthen the seven-class spectral witness, but the computation remains a local algebraic audit rather than a sourced nonproxy `K`/`L_total` spectral action.\n")
    append_once(AGENTS, "Current exact characteristic-polynomial certificate guardrail (P2785/S1735, 2026-06-16)", "## Current exact characteristic-polynomial certificate guardrail (P2785/S1735, 2026-06-16)\n\n- P2785 removes P2784's numerical-eigenvalue dependence by computing exact integer adjacency, Laplacian, and signless-Laplacian characteristic polynomials on the same seven local representatives.\n- All 21 representative pairs are separated by all three exact characteristic polynomial families; exact Kirchhoff spanning-tree counts are also recorded.\n- Do not promote this local exact-algebra witness to canonical geometry, strict spectral source law, selector closure, kernel full-expression, role-bearing `L_total`, bridge closure, role transfer, or ToE closure.  A next admissible move must supply a canonical-generation theorem/tool certificate for the full graph class, or export a strict spectral action/source law before testing.\n")
    return payload


if __name__ == "__main__":
    main()
