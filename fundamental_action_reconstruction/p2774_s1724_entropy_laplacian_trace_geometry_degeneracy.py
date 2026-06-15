#!/usr/bin/env python3
"""P2774/S1724: entropy-plus-Laplacian-trace geometry degeneracy.

P2773 showed that Shannon entropy alone does not force geometry.  This bounded
follow-up tests a natural next strengthening: add a graph-Laplacian/degree-energy
constraint.  On a 16-point support with H=ln(16)=4 ln 2, we compare two connected
4-regular geometries.  They have the same entropy, edge count, degree sequence,
and Laplacian trace, but different distance histograms.  Therefore entropy plus
this simple Laplacian trace/regularity principle still does not force a unique
nadsoliton geometry.
"""
from __future__ import annotations

import hashlib
import json
import math
from collections import deque
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
P2773 = GEN / "p2773_s1723_shannon_entropy_geometry_forcing_obstruction.json"
OUT = GEN / "p2774_s1724_entropy_laplacian_trace_geometry_degeneracy.json"
MD = GEN / "p2774_s1724_entropy_laplacian_trace_geometry_degeneracy.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

N = 16
ALPHA_GEO = 4.0 * math.log(2.0)
NEGATIVE_EXPORT_FLAGS = [
    "entropy_laplacian_trace_forces_geometry_exported",
    "canonical_geometry_from_entropy_laplacian_exported",
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


def shannon_entropy_uniform(n: int = N) -> float:
    return math.log(n)


def normalized_edges(edges: set[tuple[int, int]]) -> list[tuple[int, int]]:
    return sorted({(min(a, b), max(a, b)) for a, b in edges if a != b})


def graph_edges(kind: str) -> list[tuple[int, int]]:
    edges: set[tuple[int, int]] = set()
    if kind == "torus_4x4":
        for r in range(4):
            for c in range(4):
                node = 4 * r + c
                right = 4 * r + ((c + 1) % 4)
                down = 4 * ((r + 1) % 4) + c
                edges.add((node, right))
                edges.add((node, down))
        return normalized_edges(edges)
    if kind == "circulant_pm1_pm2":
        for i in range(N):
            for step in (1, 2):
                edges.add((i, (i + step) % N))
                edges.add((i, (i - step) % N))
        return normalized_edges(edges)
    raise ValueError(kind)


def adjacency(edges: list[tuple[int, int]], n: int = N) -> list[list[int]]:
    adj = [[] for _ in range(n)]
    for a, b in edges:
        adj[a].append(b)
        adj[b].append(a)
    return adj


def all_pairs_distances(edges: list[tuple[int, int]], n: int = N) -> list[list[int]]:
    adj = adjacency(edges, n)
    distances: list[list[int]] = []
    for source in range(n):
        dist = [-1] * n
        dist[source] = 0
        queue: deque[int] = deque([source])
        while queue:
            node = queue.popleft()
            for nxt in adj[node]:
                if dist[nxt] == -1:
                    dist[nxt] = dist[node] + 1
                    queue.append(nxt)
        if any(value < 0 for value in dist):
            raise ValueError("disconnected graph")
        distances.append(dist)
    return distances


def distance_histogram(distances: list[list[int]]) -> dict[int, int]:
    hist: dict[int, int] = {}
    for i in range(len(distances)):
        for j in range(i + 1, len(distances)):
            d = distances[i][j]
            hist[d] = hist.get(d, 0) + 1
    return dict(sorted(hist.items()))


def graph_summary(kind: str) -> dict[str, Any]:
    edges = graph_edges(kind)
    adj = adjacency(edges)
    degrees = sorted(len(row) for row in adj)
    distances = all_pairs_distances(edges)
    hist = distance_histogram(distances)
    pair_count = sum(hist.values())
    return {
        "geometry": kind,
        "node_count": N,
        "edge_count": len(edges),
        "degree_sequence": degrees,
        "is_4_regular": all(degree == 4 for degree in degrees),
        "laplacian_trace": sum(degrees),
        "diameter": max(hist),
        "average_pair_distance": sum(distance * count for distance, count in hist.items()) / pair_count,
        "distance_histogram": {str(key): value for key, value in hist.items()},
    }


def degeneracy_witness() -> dict[str, Any]:
    entropy = shannon_entropy_uniform()
    rows = [graph_summary("torus_4x4"), graph_summary("circulant_pm1_pm2")]
    distance_signatures = {json.dumps(row["distance_histogram"], sort_keys=True) for row in rows}
    degree_signatures = {json.dumps(row["degree_sequence"]) for row in rows}
    trace_values = {row["laplacian_trace"] for row in rows}
    return {
        "entropy_object": "uniform distribution on 16 support points",
        "shannon_entropy_nats": entropy,
        "alpha_geo_4_ln_2": ALPHA_GEO,
        "entropy_matches_alpha_geo": abs(entropy - ALPHA_GEO) < 1e-12,
        "candidate_extra_principle": "same 4-regular graph-Laplacian trace / degree-energy constraint",
        "geometry_rows": rows,
        "geometry_row_count": len(rows),
        "same_degree_sequence": len(degree_signatures) == 1,
        "same_laplacian_trace": len(trace_values) == 1,
        "distinct_distance_histogram_count": len(distance_signatures),
        "entropy_plus_laplacian_trace_forces_unique_geometry": len(distance_signatures) == 1,
        "finite_obstruction_statement": "H=4 ln 2 plus equal 4-regular Laplacian trace still permits inequivalent finite metric geometries.",
    }


def acceptance_matrix(witness: dict[str, Any], p2773: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "p2773_entropy_alone_no_closure_boundary_present": p2773.get("status") == "P2773_SHANNON_ENTROPY_GEOMETRY_FORCING_OBSTRUCTION_NO_CLOSURE",
        "entropy_matches_alpha_geo": witness["entropy_matches_alpha_geo"],
        "extra_laplacian_trace_principle_supplied": True,
        "same_degree_sequence_and_laplacian_trace": witness["same_degree_sequence"] and witness["same_laplacian_trace"],
        "entropy_plus_laplacian_trace_forces_unique_geometry": witness["entropy_plus_laplacian_trace_forces_unique_geometry"],
        "stronger_metric_source_supplied": False,
    }
    return {
        "facts": facts,
        "accepted_as_entropy_laplacian_geometry_forcing_theorem": False,
        "accepted_as_canonical_nadsoliton_geometry_source": False,
        "accepted_as_kernel_full_expression_theorem": False,
        "accepted_as_ltotal_or_toe_promotion": False,
        "missing_criteria": [key for key, value in facts.items() if not value],
        "blocker": "Adding 4-regular Laplacian trace/degree-energy data to H=4 ln 2 still leaves at least two inequivalent metric geometries.  A stronger metric, transport, curvature, spectrum, or adjacency source is required.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    witness = payload["entropy_laplacian_degeneracy_witness"]
    lines = [
        "# P2774/S1724 entropy-plus-Laplacian-trace geometry degeneracy",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Shared constraints",
        f"- H_uniform_16={witness['shannon_entropy_nats']}",
        f"- alpha_geo_4_ln_2={witness['alpha_geo_4_ln_2']}",
        f"- same_degree_sequence={witness['same_degree_sequence']}",
        f"- same_laplacian_trace={witness['same_laplacian_trace']}",
        "",
        "## Degeneracy result",
        f"- geometry_row_count={witness['geometry_row_count']}",
        f"- distinct_distance_histogram_count={witness['distinct_distance_histogram_count']}",
        f"- entropy_plus_laplacian_trace_forces_unique_geometry={witness['entropy_plus_laplacian_trace_forces_unique_geometry']}",
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
    p2773 = read_json(P2773)
    witness = degeneracy_witness()
    acceptance = acceptance_matrix(witness, p2773)
    payload = {
        "status": "P2774_ENTROPY_LAPLACIAN_TRACE_GEOMETRY_DEGENERACY_NO_CLOSURE",
        "input_hashes": {"P2773": sha(P2773)},
        "input_statuses": {"P2773": p2773.get("status")},
        "audited_question": "Does Shannon entropy plus a simple graph-Laplacian trace/regularity principle force nadsoliton geometry?",
        "entropy_laplacian_degeneracy_witness": witness,
        "acceptance_matrix": acceptance,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": acceptance["blocker"],
            "next_honest_step": "Do not claim geometry forcing from H=4 ln 2 plus Laplacian trace/regularity.  The next honest move must add a stronger sourced metric principle, such as a full Laplacian-spectrum selector, Ollivier/Forman curvature functional, transport-cost minimizer, or an explicit strict adjacency theorem, and rerun degeneracy-breaking against the torus_4x4 versus circulant_pm1_pm2 pair; otherwise preserve the P2697-P2774 no-geometry-forcing/no-closure certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2774/S1724 entropy-plus-Laplacian-trace geometry degeneracy", "## P2774/S1724 entropy-plus-Laplacian-trace geometry degeneracy\n\n`P2774/S1724` strengthens P2773 by adding a simple graph-Laplacian trace/regularity constraint to `H=ln(16)=4 ln 2=alpha_geo`.  Two connected 4-regular 16-node geometries, `torus_4x4` and `circulant_pm1_pm2`, have the same entropy, degree sequence, edge count, and Laplacian trace, but different graph-distance histograms.  Thus entropy plus this Laplacian trace/degree-energy datum still does not select a canonical nadsoliton metric/adjacency.  No kernel full-expression theorem, role-bearing `L_total`, bridge closure, role transfer, selector closure, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2774/S1724 entropy-Laplacian Ltotal guard", "## P2774/S1724 entropy-Laplacian Ltotal guard\n\n`P2774/S1724` adds no variational source term.  It shows that `H=4 ln 2` plus equal 4-regular Laplacian trace remains geometrically degenerate, so this strengthened entropy principle cannot promote the kernels to role-bearing `L_total` or canonical nadsoliton geometry.\n")
    append_once(AGENTS, "Current entropy-plus-Laplacian-trace geometry degeneracy guardrail (P2774/S1724, 2026-06-15)", "## Current entropy-plus-Laplacian-trace geometry degeneracy guardrail (P2774/S1724, 2026-06-15)\n\n- P2774 strengthens P2773 by adding a graph-Laplacian trace/regularity constraint to `H=ln(16)=4 ln 2=alpha_geo`.\n- The connected 4-regular geometries `torus_4x4` and `circulant_pm1_pm2` share entropy, degree sequence, edge count, and Laplacian trace, but have different distance histograms; entropy plus this simple Laplacian datum still does not force a canonical geometry.\n- Do not promote entropy-plus-Laplacian-trace data to canonical nadsoliton geometry, kernel full-expression, role-bearing `L_total`, bridge closure, role transfer, selector closure, or ToE closure.  A next admissible move must add a stronger sourced metric principle and pass degeneracy-breaking against the P2774 pair.\n")
    return payload


if __name__ == "__main__":
    main()
