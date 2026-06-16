#!/usr/bin/env python3
"""P2787/S1737: small canonical generator pipeline audit.

P2786 showed that the current environment lacks a certified full 16-node
4-regular generator/toolchain.  This bounded follow-up does not jump to the
blocked 16-node problem.  Instead, it proves the in-repo canonical-generation
pipeline on the largest small regular class that is still cheap enough to audit
self-containedly here: all connected 8-node 4-regular simple graphs.

The script enumerates every labeled 8-node 4-regular graph by degree-constrained
backtracking, verifies connectivity, quotients by direct graph isomorphism, and
then runs exact adjacency/Laplacian/signless characteristic-polynomial collision
audits on the quotient representatives.  It is a pipeline validation and scale
bridge only, not a canonical nadsoliton geometry theorem.
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
from p2786_s1736_graph6_provenance_toolchain_gate import graph6_encode

GEN = ROOT / "generated"
P2786 = GEN / "p2786_s1736_graph6_provenance_toolchain_gate.json"
OUT = GEN / "p2787_s1737_small_canonical_generator_pipeline_audit.json"
MD = GEN / "p2787_s1737_small_canonical_generator_pipeline_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

N_SMALL = 8
D_SMALL = 4
EDGE_TARGET = N_SMALL * D_SMALL // 2

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


def connected(edges: tuple[tuple[int, int], ...], n: int = N_SMALL) -> bool:
    adj = [set() for _ in range(n)]
    for a, b in edges:
        adj[a].add(b)
        adj[b].add(a)
    seen = {0}
    stack = [0]
    while stack:
        node = stack.pop()
        for nxt in adj[node] - seen:
            seen.add(nxt)
            stack.append(nxt)
    return len(seen) == n


def enumerate_labeled_regular_graphs(n: int = N_SMALL, degree: int = D_SMALL) -> list[tuple[tuple[int, int], ...]]:
    possible_edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    out: list[tuple[tuple[int, int], ...]] = []

    def rec(index: int, degrees: tuple[int, ...], chosen: tuple[tuple[int, int], ...]) -> None:
        if any(value > degree for value in degrees):
            return
        if len(chosen) > n * degree // 2:
            return
        for vertex in range(n):
            remaining = sum(1 for cursor in range(index, len(possible_edges)) if vertex in possible_edges[cursor])
            if degrees[vertex] + remaining < degree:
                return
        if index == len(possible_edges):
            if all(value == degree for value in degrees):
                out.append(normalize_edges(chosen))
            return
        a, b = possible_edges[index]
        if degrees[a] < degree and degrees[b] < degree:
            updated = list(degrees)
            updated[a] += 1
            updated[b] += 1
            rec(index + 1, tuple(updated), chosen + ((a, b),))
        rec(index + 1, degrees, chosen)

    rec(0, tuple([0] * n), tuple())
    return out


def adjacency_sets(edges: tuple[tuple[int, int], ...], n: int = N_SMALL) -> list[set[int]]:
    adj = [set() for _ in range(n)]
    for a, b in edges:
        adj[a].add(b)
        adj[b].add(a)
    return adj


def isomorphic(left: tuple[tuple[int, int], ...], right: tuple[tuple[int, int], ...], n: int = N_SMALL) -> bool:
    adj_left = adjacency_sets(left, n)
    adj_right = adjacency_sets(right, n)
    order = sorted(range(n), key=lambda node: (-len(adj_left[node]), node))
    mapping: dict[int, int] = {}
    used: set[int] = set()

    def compatible(source: int, target: int) -> bool:
        return all((mapped_source in adj_left[source]) == (mapped_target in adj_right[target]) for mapped_source, mapped_target in mapping.items())

    def rec(index: int) -> bool:
        if index == n:
            return True
        source = order[index]
        for target in range(n):
            if target in used or not compatible(source, target):
                continue
            mapping[source] = target
            used.add(target)
            if rec(index + 1):
                return True
            used.remove(target)
            del mapping[source]
        return False

    return rec(0)


def quotient_classes(graphs: list[tuple[tuple[int, int], ...]]) -> list[dict[str, Any]]:
    classes: list[dict[str, Any]] = []
    for graph in graphs:
        for row in classes:
            if isomorphic(graph, row["representative_edges"]):
                row["labeled_member_count"] += 1
                break
        else:
            classes.append({
                "representative_index": len(classes),
                "representative_edges": graph,
                "labeled_member_count": 1,
            })
    return classes


def charpoly_coefficients(edges: tuple[tuple[int, int], ...], kind: str) -> list[int]:
    adj = sp.zeros(N_SMALL, N_SMALL)
    for a, b in edges:
        adj[a, b] = 1
        adj[b, a] = 1
    degree = sp.diag(*[sum(adj[row, col] for col in range(N_SMALL)) for row in range(N_SMALL)])
    if kind == "adjacency":
        mat = adj
    elif kind == "laplacian":
        mat = degree - adj
    elif kind == "signless_laplacian":
        mat = degree + adj
    else:
        raise ValueError(kind)
    return [int(coeff) for coeff in mat.charpoly().all_coeffs()]


def exact_pipeline_witness() -> dict[str, Any]:
    labeled = enumerate_labeled_regular_graphs()
    connected_graphs = [graph for graph in labeled if connected(graph)]
    classes = quotient_classes(connected_graphs)
    rows = []
    for row in classes:
        edges = row["representative_edges"]
        exact_payload = {
            "adjacency_charpoly_coefficients": charpoly_coefficients(edges, "adjacency"),
            "laplacian_charpoly_coefficients": charpoly_coefficients(edges, "laplacian"),
            "signless_laplacian_charpoly_coefficients": charpoly_coefficients(edges, "signless_laplacian"),
        }
        rows.append({
            "representative_index": row["representative_index"],
            "labeled_member_count": row["labeled_member_count"],
            "edge_count": len(edges),
            "graph6": graph6_encode(list(edges), n=N_SMALL),
            **exact_payload,
        })
    invariant_keys = [
        "adjacency_charpoly_coefficients",
        "laplacian_charpoly_coefficients",
        "signless_laplacian_charpoly_coefficients",
    ]
    pair_rows = []
    for left, right in itertools.combinations(rows, 2):
        collisions = {key: left[key] == right[key] for key in invariant_keys}
        pair_rows.append({
            "left": left["representative_index"],
            "right": right["representative_index"],
            "exact_charpoly_collisions": collisions,
            "all_three_charpolys_distinct": all(not value for value in collisions.values()),
        })
    collision_counts = {key: sum(1 for row in pair_rows if row["exact_charpoly_collisions"][key]) for key in invariant_keys}
    return {
        "source_class": "All connected simple 8-node 4-regular graphs generated self-containedly by degree-constrained backtracking.",
        "labeled_candidate_count": len(labeled),
        "connected_labeled_candidate_count": len(connected_graphs),
        "isomorphism_class_count": len(classes),
        "pair_count_after_quotient": len(pair_rows),
        "representative_rows": rows,
        "pair_rows": pair_rows,
        "exact_charpoly_collision_counts": collision_counts,
        "all_pairs_separated_by_all_three_exact_charpolys": all(row["all_three_charpolys_distinct"] for row in pair_rows),
        "finite_certificate_statement": "The self-contained 8-node 4-regular generator pipeline yields 19,355 connected labeled graphs, 6 isomorphism classes, and zero exact charpoly collisions across all 15 quotient pairs.",
    }


def acceptance_matrix(witness: dict[str, Any], p2786: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "p2786_graph6_gate_present": p2786.get("status") == "P2786_GRAPH6_PROVENANCE_TOOLCHAIN_GATE_NO_CLOSURE",
        "self_contained_small_generator_executed": witness["labeled_candidate_count"] == 19355,
        "small_connected_class_quotiented": witness["isomorphism_class_count"] == 6,
        "all_15_small_quotient_pairs_checked": witness["pair_count_after_quotient"] == 15,
        "zero_exact_charpoly_collisions_on_small_quotient": all(value == 0 for value in witness["exact_charpoly_collision_counts"].values()),
        "canonical_16node_generator_certified": False,
        "strict_nadsoliton_spectral_source_law_exported": False,
        "kernel_or_ltotal_variational_coupling_exported": False,
    }
    return {
        "facts": facts,
        "accepted_as_small_self_contained_generator_pipeline_certificate": all(facts[key] for key in [
            "p2786_graph6_gate_present",
            "self_contained_small_generator_executed",
            "small_connected_class_quotiented",
            "all_15_small_quotient_pairs_checked",
            "zero_exact_charpoly_collisions_on_small_quotient",
        ]),
        "accepted_as_full_16node_canonical_generator_certificate": False,
        "accepted_as_strict_spectral_source_law": False,
        "accepted_as_canonical_nadsoliton_geometry_source": False,
        "accepted_as_ltotal_or_toe_promotion": False,
        "missing_criteria": [key for key, value in facts.items() if not value],
        "blocker": "The exact generator/quotient/charpoly pipeline is validated on the complete 8-node 4-regular class, but it is not the required full 16-node class and exports no strict K/L_total source law.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    w = payload["small_pipeline_witness"]
    lines = [
        "# P2787/S1737 small canonical generator pipeline audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Exact small-pipeline result",
        f"- labeled_candidate_count={w['labeled_candidate_count']}",
        f"- connected_labeled_candidate_count={w['connected_labeled_candidate_count']}",
        f"- isomorphism_class_count={w['isomorphism_class_count']}",
        f"- pair_count_after_quotient={w['pair_count_after_quotient']}",
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
    p2786 = read_json(P2786)
    witness = exact_pipeline_witness()
    acceptance = acceptance_matrix(witness, p2786)
    payload = {
        "status": "P2787_SMALL_CANONICAL_GENERATOR_PIPELINE_AUDIT_NO_CLOSURE",
        "input_hashes": {"P2786": sha(P2786)},
        "input_statuses": {"P2786": p2786.get("status")},
        "audited_question": "Can the in-repo generator/quotient/exact-charpoly pipeline be validated on a complete smaller regular graph class before the blocked 16-node class?",
        "small_pipeline_witness": witness,
        "acceptance_matrix": acceptance,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": acceptance["blocker"],
            "next_honest_step": "Use P2787 only as a complete small-class pipeline validation.  The next honest move is exactly one of: supply/import an actual certified full connected 16-node 4-regular generator artifact/toolchain with graph6/hash provenance and run the same exact quotient/charpoly audit there; or export a strict nadsoliton spectral action/source law fixing the admissible class, target spectrum, and K/L_total coupling before testing.  Otherwise preserve the P2697-P2787 no-canonical-geometry/no-closure certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2787/S1737 small canonical generator pipeline audit", "## P2787/S1737 small canonical generator pipeline audit\n\n`P2787/S1737` validates the exact generator/quotient/characteristic-polynomial pipeline on the complete smaller class of connected 8-node 4-regular simple graphs.  Degree-constrained backtracking gives 19,355 connected labeled candidates, direct isomorphism quotienting gives 6 classes, and exact adjacency/Laplacian/signless-Laplacian characteristic polynomials have zero collisions across all 15 quotient pairs.  This is only a small-class pipeline certificate: it is not the blocked full connected 16-node 4-regular generator, not a strict spectral source law, and not a `K`/`L_total` variational coupling.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2787/S1737 small-pipeline Ltotal guard", "## P2787/S1737 small-pipeline Ltotal guard\n\n`P2787/S1737` adds no variational source term.  The exact generator/quotient/charpoly pipeline works on the complete 8-node 4-regular class, but this is a computational pipeline validation rather than a sourced nonproxy `K`/`L_total` spectral action or a full 16-node canonical generator.\n")
    append_once(AGENTS, "Current small canonical generator pipeline guardrail (P2787/S1737, 2026-06-16)", "## Current small canonical generator pipeline guardrail (P2787/S1737, 2026-06-16)\n\n- P2787 validates the self-contained generator/quotient/exact-charpoly pipeline on all connected 8-node 4-regular simple graphs: 19,355 labeled connected candidates, 6 isomorphism classes, and zero exact charpoly collisions across 15 quotient pairs.\n- This is a complete small-class certificate only; it is not the required full connected 16-node 4-regular generator/toolchain and does not source geometry from `K`/`L_total`.\n- Do not promote this small-pipeline witness to canonical geometry, strict spectral source law, selector closure, kernel full-expression, role-bearing `L_total`, bridge closure, role transfer, or ToE closure.  A next admissible move must supply an actual certified full 16-node generator artifact/toolchain or export a strict spectral action/source law before testing.\n")
    return payload


if __name__ == "__main__":
    main()
