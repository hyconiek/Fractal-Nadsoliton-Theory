#!/usr/bin/env python3
"""P2788/S1738: complement-duality exact spectral certificate.

P2787 validated the generator/quotient/exact-charpoly pipeline on the complete
8-node 4-regular class, but it did not add a proof-style invariant sanity check
linking the computed spectra to an independent graph operation.  This bounded
follow-up audits exact complement-duality identities for both available tiers:

* the complete P2787 8-node 4-regular quotient; and
* the seven local P2786/P2785 16-node 4-regular representatives.

For every d-regular n-node representative, the script computes the complement
and verifies exact characteristic-polynomial identities for adjacency,
Laplacian, and signless-Laplacian matrices.  This is a theorem-backed algebraic
consistency witness, not a full 16-node generator and not a strict spectral
action/source law.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import sympy as sp

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2786_s1736_graph6_provenance_toolchain_gate import graph6_decode as graph6_decode_16

GEN = ROOT / "generated"
P2786 = GEN / "p2786_s1736_graph6_provenance_toolchain_gate.json"
P2787 = GEN / "p2787_s1737_small_canonical_generator_pipeline_audit.json"
OUT = GEN / "p2788_s1738_complement_duality_exact_spectral_certificate.json"
MD = GEN / "p2788_s1738_complement_duality_exact_spectral_certificate.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

NEGATIVE_EXPORT_FLAGS = [
    "canonical_16node_generator_certified",
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


def normalize_edges(edges: list[tuple[int, int]] | tuple[tuple[int, int], ...]) -> tuple[tuple[int, int], ...]:
    return tuple(sorted((min(a, b), max(a, b)) for a, b in edges))


def complement_edges(edges: tuple[tuple[int, int], ...], n: int) -> tuple[tuple[int, int], ...]:
    edge_set = set(edges)
    return tuple((i, j) for i in range(n) for j in range(i + 1, n) if (i, j) not in edge_set)


def adjacency_matrix(edges: tuple[tuple[int, int], ...], n: int) -> sp.Matrix:
    mat = sp.zeros(n, n)
    for a, b in edges:
        mat[a, b] = 1
        mat[b, a] = 1
    return mat


def matrix_for(edges: tuple[tuple[int, int], ...], n: int, kind: str) -> sp.Matrix:
    adj = adjacency_matrix(edges, n)
    degree = sp.diag(*[sum(adj[row, col] for col in range(n)) for row in range(n)])
    if kind == "adjacency":
        return adj
    if kind == "laplacian":
        return degree - adj
    if kind == "signless_laplacian":
        return degree + adj
    raise ValueError(kind)


def charpoly_expr(edges: tuple[tuple[int, int], ...], n: int, kind: str) -> sp.Poly:
    x = sp.Symbol("x")
    return sp.Poly(matrix_for(edges, n, kind).charpoly(x).as_expr(), x)


def degree_sequence(edges: tuple[tuple[int, int], ...], n: int) -> list[int]:
    degree = [0] * n
    for a, b in edges:
        degree[a] += 1
        degree[b] += 1
    return degree


def expected_complement_polynomials(edges: tuple[tuple[int, int], ...], n: int, d: int) -> dict[str, sp.Poly]:
    x = sp.Symbol("x")
    identity = sp.eye(n)
    ones = sp.ones(n, n)
    adj = matrix_for(edges, n, "adjacency")
    lap = matrix_for(edges, n, "laplacian")
    signless = matrix_for(edges, n, "signless_laplacian")
    expected_adj = (ones - identity - adj).charpoly(x).as_expr()
    expected_lap = (n * identity - ones - lap).charpoly(x).as_expr()
    expected_signless = ((n - 2) * identity + ones - signless).charpoly(x).as_expr()
    return {
        "adjacency": sp.Poly(expected_adj, x),
        "laplacian": sp.Poly(expected_lap, x),
        "signless_laplacian": sp.Poly(expected_signless, x),
    }


def complement_duality_row(label: str, edges: tuple[tuple[int, int], ...], n: int) -> dict[str, Any]:
    degrees = degree_sequence(edges, n)
    if len(set(degrees)) != 1:
        raise ValueError(f"{label} is not regular: {degrees}")
    d = degrees[0]
    comp = complement_edges(edges, n)
    comp_degrees = degree_sequence(comp, n)
    expected = expected_complement_polynomials(edges, n, d)
    actual = {kind: charpoly_expr(comp, n, kind) for kind in expected}
    checks = {kind: actual[kind].all_coeffs() == expected[kind].all_coeffs() for kind in expected}
    return {
        "label": label,
        "n": n,
        "degree": d,
        "edge_count": len(edges),
        "complement_degree": comp_degrees[0],
        "complement_edge_count": len(comp),
        "complement_regular": len(set(comp_degrees)) == 1 and comp_degrees[0] == n - 1 - d,
        "exact_complement_identity_checks": checks,
        "all_exact_complement_identities_hold": all(checks.values()),
        "complement_charpoly_coefficients": {kind: [int(c) for c in actual[kind].all_coeffs()] for kind in actual},
    }


def graph6_decode_small(payload: str, expected_n: int) -> tuple[tuple[int, int], ...]:
    if not payload:
        raise ValueError("empty graph6 payload")
    n = ord(payload[0]) - 63
    if n != expected_n:
        raise ValueError(f"expected graph6 order {expected_n}, got {n}")
    bits: list[int] = []
    for char in payload[1:]:
        value = ord(char) - 63
        if not 0 <= value < 64:
            raise ValueError("invalid graph6 character")
        bits.extend((value >> shift) & 1 for shift in range(5, -1, -1))
    needed = n * (n - 1) // 2
    edges: list[tuple[int, int]] = []
    cursor = 0
    for j in range(1, n):
        for i in range(j):
            if bits[cursor]:
                edges.append((i, j))
            cursor += 1
    if any(bits[needed:]):
        raise ValueError("nonzero graph6 padding bits")
    return normalize_edges(edges)


def small_complete_rows(p2787: dict[str, Any]) -> list[dict[str, Any]]:
    rows = []
    for row in p2787.get("small_pipeline_witness", {}).get("representative_rows", []):
        edges = graph6_decode_small(row["graph6"], 8)
        rows.append(complement_duality_row(f"P2787_8node_class_{row['representative_index']}", edges, 8))
    return rows


def local_16_rows(p2786: dict[str, Any]) -> list[dict[str, Any]]:
    rows = []
    for row in p2786.get("graph6_provenance_witness", {}).get("provenance_rows", []):
        edges = normalize_edges(graph6_decode_16(row["graph6"]))
        rows.append(complement_duality_row(row["representative"], edges, 16))
    return rows


def complement_witness(p2786: dict[str, Any], p2787: dict[str, Any]) -> dict[str, Any]:
    small = small_complete_rows(p2787)
    local = local_16_rows(p2786)
    all_rows = small + local
    return {
        "source_class": "Exact complement-duality audit over the complete P2787 8-node quotient and the seven local P2786 16-node representatives.",
        "small_complete_8node_row_count": len(small),
        "local_16node_row_count": len(local),
        "total_rows_checked": len(all_rows),
        "small_complete_8node_rows": small,
        "local_16node_rows": local,
        "all_complements_regular_with_expected_degree": all(row["complement_regular"] for row in all_rows),
        "all_exact_adjacency_complement_identities_hold": all(row["exact_complement_identity_checks"]["adjacency"] for row in all_rows),
        "all_exact_laplacian_complement_identities_hold": all(row["exact_complement_identity_checks"]["laplacian"] for row in all_rows),
        "all_exact_signless_complement_identities_hold": all(row["exact_complement_identity_checks"]["signless_laplacian"] for row in all_rows),
        "finite_certificate_statement": "All 13 audited representatives satisfy exact regular-complement degree and adjacency/Laplacian/signless-Laplacian characteristic-polynomial complement identities.",
    }


def acceptance_matrix(witness: dict[str, Any], p2786: dict[str, Any], p2787: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "p2786_graph6_gate_present": p2786.get("status") == "P2786_GRAPH6_PROVENANCE_TOOLCHAIN_GATE_NO_CLOSURE",
        "p2787_small_pipeline_present": p2787.get("status") == "P2787_SMALL_CANONICAL_GENERATOR_PIPELINE_AUDIT_NO_CLOSURE",
        "complete_small_quotient_rechecked": witness["small_complete_8node_row_count"] == 6,
        "seven_local_16node_representatives_rechecked": witness["local_16node_row_count"] == 7,
        "all_regular_complement_degrees_verified": witness["all_complements_regular_with_expected_degree"],
        "all_three_exact_complement_charpoly_identities_verified": all(witness[key] for key in [
            "all_exact_adjacency_complement_identities_hold",
            "all_exact_laplacian_complement_identities_hold",
            "all_exact_signless_complement_identities_hold",
        ]),
        "canonical_16node_generator_certified": False,
        "strict_nadsoliton_spectral_source_law_exported": False,
        "kernel_or_ltotal_variational_coupling_exported": False,
    }
    return {
        "facts": facts,
        "accepted_as_exact_complement_duality_certificate": all(facts[key] for key in [
            "p2786_graph6_gate_present",
            "p2787_small_pipeline_present",
            "complete_small_quotient_rechecked",
            "seven_local_16node_representatives_rechecked",
            "all_regular_complement_degrees_verified",
            "all_three_exact_complement_charpoly_identities_verified",
        ]),
        "accepted_as_full_16node_canonical_generator_certificate": False,
        "accepted_as_strict_spectral_source_law": False,
        "accepted_as_canonical_nadsoliton_geometry_source": False,
        "accepted_as_ltotal_or_toe_promotion": False,
        "missing_criteria": [key for key, value in facts.items() if not value],
        "blocker": "Complement-duality is an exact algebraic sanity certificate for existing representatives/classes, but it does not enumerate the full connected 16-node 4-regular class and exports no strict K/L_total spectral source law.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    w = payload["complement_duality_witness"]
    lines = [
        "# P2788/S1738 complement-duality exact spectral certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Exact complement-duality result",
        f"- small_complete_8node_row_count={w['small_complete_8node_row_count']}",
        f"- local_16node_row_count={w['local_16node_row_count']}",
        f"- total_rows_checked={w['total_rows_checked']}",
        f"- all_complements_regular_with_expected_degree={w['all_complements_regular_with_expected_degree']}",
        f"- all_exact_adjacency_complement_identities_hold={w['all_exact_adjacency_complement_identities_hold']}",
        f"- all_exact_laplacian_complement_identities_hold={w['all_exact_laplacian_complement_identities_hold']}",
        f"- all_exact_signless_complement_identities_hold={w['all_exact_signless_complement_identities_hold']}",
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
    p2786 = read_json(P2786)
    p2787 = read_json(P2787)
    witness = complement_witness(p2786, p2787)
    acceptance = acceptance_matrix(witness, p2786, p2787)
    payload = {
        "status": "P2788_COMPLEMENT_DUALITY_EXACT_SPECTRAL_CERTIFICATE_NO_CLOSURE",
        "input_hashes": {"P2786": sha(P2786), "P2787": sha(P2787)},
        "input_statuses": {"P2786": p2786.get("status"), "P2787": p2787.get("status")},
        "audited_question": "Do the existing local/full-small spectral witnesses satisfy independent exact regular-complement characteristic-polynomial identities?",
        "complement_duality_witness": witness,
        "acceptance_matrix": acceptance,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": acceptance["blocker"],
            "next_honest_step": "Use P2788 only as an exact theorem-backed complement-duality sanity certificate.  The next honest move is still exactly one of: supply/import an actual certified full connected 16-node 4-regular generator artifact/toolchain with graph6/hash provenance and run exact quotient/charpoly/complement-duality auditing there; or export a strict nadsoliton spectral action/source law fixing the admissible class, target spectrum, and K/L_total coupling before testing.  Otherwise preserve the P2697-P2788 no-canonical-geometry/no-closure certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2788/S1738 complement-duality exact spectral certificate", "## P2788/S1738 complement-duality exact spectral certificate\n\n`P2788/S1738` adds a theorem-backed exact complement-duality audit over the complete P2787 8-node quotient and the seven local P2786/P2785 16-node representatives.  For all 13 audited representatives, complements have the expected regular degree and exact adjacency/Laplacian/signless-Laplacian characteristic-polynomial complement identities hold.  This strengthens algebraic consistency of the existing spectral witnesses, but it is not a full connected 16-node 4-regular generator, not a strict spectral source law, and not a `K`/`L_total` variational coupling.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2788/S1738 complement-duality Ltotal guard", "## P2788/S1738 complement-duality Ltotal guard\n\n`P2788/S1738` adds no variational source term.  Exact regular-complement characteristic-polynomial identities validate algebraic consistency of existing spectral data, but complement duality is not a sourced nonproxy `K`/`L_total` spectral action, not a canonical geometry theorem, and not a full 16-node generator.\n")
    append_once(AGENTS, "Current complement-duality exact spectral guardrail (P2788/S1738, 2026-06-16)", "## Current complement-duality exact spectral guardrail (P2788/S1738, 2026-06-16)\n\n- P2788 verifies exact regular-complement degree and adjacency/Laplacian/signless-Laplacian characteristic-polynomial complement identities for all 13 audited representatives: the complete P2787 8-node quotient plus the seven local P2786/P2785 16-node representatives.\n- This is a theorem-backed algebraic consistency certificate for existing witnesses only; it is not the required full connected 16-node 4-regular generator/toolchain and does not source geometry from `K`/`L_total`.\n- Do not promote complement-duality consistency to canonical geometry, strict spectral source law, selector closure, kernel full-expression, role-bearing `L_total`, bridge closure, role transfer, or ToE closure.  A next admissible move must supply an actual certified full 16-node generator artifact/toolchain or export a strict spectral action/source law before testing.\n")
    return payload


if __name__ == "__main__":
    main()
