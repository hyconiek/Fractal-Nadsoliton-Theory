#!/usr/bin/env python3
"""P2790/S1740: eighth 16-node witness no-exhaustion certificate.

P2789 validated quotient arithmetic for the seven local 16-node representatives,
but the honest blocker remains that those representatives are not a full
connected 16-node 4-regular class.  This bounded follow-up supplies one explicit
additional connected 16-node 4-regular graph and verifies it is not isomorphic
to any of the seven local representatives.

The result is deliberately negative for closure: it proves the seven-local set
is not exhaustive, while still not providing a full generator or a strict
spectral source law.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2786_s1736_graph6_provenance_toolchain_gate import graph6_decode, graph6_encode
from p2787_s1737_small_canonical_generator_pipeline_audit import isomorphic
from p2788_s1738_complement_duality_exact_spectral_certificate import charpoly_expr, normalize_edges
from p2789_s1739_orbit_stabilizer_exact_quotient_certificate import automorphism_count

GEN = ROOT / "generated"
P2786 = GEN / "p2786_s1736_graph6_provenance_toolchain_gate.json"
P2789 = GEN / "p2789_s1739_orbit_stabilizer_exact_quotient_certificate.json"
OUT = GEN / "p2790_s1740_eighth_16node_witness_no_exhaustion_certificate.json"
MD = GEN / "p2790_s1740_eighth_16node_witness_no_exhaustion_certificate.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

N = 16
D = 4
EIGHTH_WITNESS_EDGES = normalize_edges([
    (0, 1), (0, 5), (0, 8), (0, 10),
    (1, 2), (1, 6), (1, 7),
    (2, 3), (2, 11), (2, 13),
    (3, 6), (3, 14), (3, 15),
    (4, 8), (4, 12), (4, 14), (4, 15),
    (5, 7), (5, 9), (5, 11),
    (6, 9), (6, 13),
    (7, 9), (7, 10),
    (8, 12), (8, 15),
    (9, 13),
    (10, 11), (10, 14),
    (11, 13),
    (12, 14), (12, 15),
])

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


def degree_sequence(edges: tuple[tuple[int, int], ...], n: int = N) -> list[int]:
    degree = [0] * n
    for a, b in edges:
        degree[a] += 1
        degree[b] += 1
    return degree


def connected(edges: tuple[tuple[int, int], ...], n: int = N) -> bool:
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


def local_rows(p2786: dict[str, Any]) -> list[dict[str, Any]]:
    rows = []
    for row in p2786.get("graph6_provenance_witness", {}).get("provenance_rows", []):
        rows.append({
            "representative": row["representative"],
            "edges": normalize_edges(graph6_decode(row["graph6"])),
        })
    return rows


def eighth_witness(p2786: dict[str, Any]) -> dict[str, Any]:
    rows = local_rows(p2786)
    pair_checks = []
    for row in rows:
        same_iso_class = isomorphic(EIGHTH_WITNESS_EDGES, row["edges"], N)
        pair_checks.append({
            "against_representative": row["representative"],
            "isomorphic_to_eighth_witness": same_iso_class,
        })
    aut_size = automorphism_count(EIGHTH_WITNESS_EDGES, N)
    graph6 = graph6_encode(list(EIGHTH_WITNESS_EDGES), n=N)
    charpolys = {
        "adjacency": [int(c) for c in charpoly_expr(EIGHTH_WITNESS_EDGES, N, "adjacency").all_coeffs()],
        "laplacian": [int(c) for c in charpoly_expr(EIGHTH_WITNESS_EDGES, N, "laplacian").all_coeffs()],
        "signless_laplacian": [int(c) for c in charpoly_expr(EIGHTH_WITNESS_EDGES, N, "signless_laplacian").all_coeffs()],
    }
    return {
        "source_class": "One explicit connected 16-node 4-regular witness outside the seven local P2786/P2785 representatives.",
        "edge_count": len(EIGHTH_WITNESS_EDGES),
        "degree_sequence": degree_sequence(EIGHTH_WITNESS_EDGES),
        "is_connected": connected(EIGHTH_WITNESS_EDGES),
        "graph6": graph6,
        "graph6_sha256": hashlib.sha256(graph6.encode("ascii")).hexdigest(),
        "automorphism_group_size": aut_size,
        "orbit_size_by_orbit_stabilizer": math.factorial(N) // aut_size,
        "exact_charpoly_coefficients": charpolys,
        "local_pair_isomorphism_checks": pair_checks,
        "is_nonisomorphic_to_all_seven_local_representatives": all(not row["isomorphic_to_eighth_witness"] for row in pair_checks),
        "finite_certificate_statement": "The explicit eighth witness is connected 16-node 4-regular and non-isomorphic to all seven local representatives, so the seven-local set is not exhaustive.",
    }


def acceptance_matrix(witness: dict[str, Any], p2786: dict[str, Any], p2789: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "p2786_graph6_gate_present": p2786.get("status") == "P2786_GRAPH6_PROVENANCE_TOOLCHAIN_GATE_NO_CLOSURE",
        "p2789_orbit_stabilizer_present": p2789.get("status") == "P2789_ORBIT_STABILIZER_EXACT_QUOTIENT_CERTIFICATE_NO_CLOSURE",
        "witness_has_16_vertices_degree_4_edge_count_32": witness["edge_count"] == 32 and witness["degree_sequence"] == [4] * 16,
        "witness_connected": witness["is_connected"],
        "witness_nonisomorphic_to_all_seven": witness["is_nonisomorphic_to_all_seven_local_representatives"],
        "witness_stabilizer_positive": witness["automorphism_group_size"] > 0,
        "canonical_16node_generator_certified": False,
        "strict_nadsoliton_spectral_source_law_exported": False,
        "kernel_or_ltotal_variational_coupling_exported": False,
    }
    return {
        "facts": facts,
        "accepted_as_eighth_witness_no_exhaustion_certificate": all(facts[key] for key in [
            "p2786_graph6_gate_present",
            "p2789_orbit_stabilizer_present",
            "witness_has_16_vertices_degree_4_edge_count_32",
            "witness_connected",
            "witness_nonisomorphic_to_all_seven",
            "witness_stabilizer_positive",
        ]),
        "accepted_as_full_16node_canonical_generator_certificate": False,
        "accepted_as_strict_spectral_source_law": False,
        "accepted_as_canonical_nadsoliton_geometry_source": False,
        "accepted_as_ltotal_or_toe_promotion": False,
        "missing_criteria": [key for key, value in facts.items() if not value],
        "blocker": "An eighth non-isomorphic connected 16-node 4-regular witness proves the seven-local set is not exhaustive, but one extra witness is still not a full generator and exports no strict K/L_total spectral source law.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    w = payload["eighth_witness"]
    lines = [
        "# P2790/S1740 eighth 16-node witness no-exhaustion certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Exact eighth-witness result",
        f"- edge_count={w['edge_count']}",
        f"- degree_sequence={w['degree_sequence']}",
        f"- is_connected={w['is_connected']}",
        f"- automorphism_group_size={w['automorphism_group_size']}",
        f"- is_nonisomorphic_to_all_seven_local_representatives={w['is_nonisomorphic_to_all_seven_local_representatives']}",
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
    p2789 = read_json(P2789)
    witness = eighth_witness(p2786)
    acceptance = acceptance_matrix(witness, p2786, p2789)
    payload = {
        "status": "P2790_EIGHTH_16NODE_WITNESS_NO_EXHAUSTION_CERTIFICATE_NO_CLOSURE",
        "input_hashes": {"P2786": sha(P2786), "P2789": sha(P2789)},
        "input_statuses": {"P2786": p2786.get("status"), "P2789": p2789.get("status")},
        "audited_question": "Can the seven local 16-node representatives be honestly treated as exhaustive?",
        "eighth_witness": witness,
        "acceptance_matrix": acceptance,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": acceptance["blocker"],
            "next_honest_step": "Use P2790 as a no-exhaustion certificate: the seven-local 16-node set is now explicitly known to be incomplete.  The next honest move is exactly one of: supply/import an actual certified full connected 16-node 4-regular generator artifact/toolchain with graph6/hash provenance and run full exact quotient/charpoly/complement/orbit auditing; or export a strict nadsoliton spectral action/source law fixing the admissible class, target spectrum, and K/L_total coupling before testing.  Otherwise preserve the P2697-P2790 no-canonical-geometry/no-closure certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2790/S1740 eighth 16-node witness no-exhaustion certificate", "## P2790/S1740 eighth 16-node witness no-exhaustion certificate\n\n`P2790/S1740` adds one explicit connected 16-node 4-regular graph outside the seven local P2786/P2785 representatives.  Exact pairwise isomorphism checks show it is non-isomorphic to all seven local classes; its automorphism group has size 2 and exact spectral data are recorded.  This is a no-exhaustion certificate for the seven-local set, not a full connected 16-node 4-regular generator, not a strict spectral source law, and not a `K`/`L_total` variational coupling.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2790/S1740 eighth-witness Ltotal guard", "## P2790/S1740 eighth-witness Ltotal guard\n\n`P2790/S1740` adds no variational source term.  The eighth non-isomorphic 16-node graph proves the seven-local graph set is incomplete, but incompleteness evidence is not a sourced nonproxy `K`/`L_total` spectral action, not a canonical geometry theorem, and not a full 16-node generator.\n")
    append_once(AGENTS, "Current eighth 16-node witness no-exhaustion guardrail (P2790/S1740, 2026-06-16)", "## Current eighth 16-node witness no-exhaustion guardrail (P2790/S1740, 2026-06-16)\n\n- P2790 constructs one explicit connected 16-node 4-regular graph with 32 edges, automorphism group size 2, and exact spectral data; pairwise isomorphism checks show it is non-isomorphic to all seven local P2786/P2785 representatives.\n- This proves the seven-local representative set is not exhaustive, but it is not the required full connected 16-node 4-regular generator/toolchain and does not source geometry from `K`/`L_total`.\n- Do not promote the eighth-witness no-exhaustion result to canonical geometry, strict spectral source law, selector closure, kernel full-expression, role-bearing `L_total`, bridge closure, role transfer, or ToE closure.  A next admissible move must supply an actual certified full 16-node generator artifact/toolchain or export a strict spectral action/source law before testing.\n")
    return payload


if __name__ == "__main__":
    main()
