#!/usr/bin/env python3
"""P2773/S1723: Shannon-entropy geometry-forcing obstruction.

The user asked whether the nadsoliton geometry might be forced by Shannon
entropy.  This bounded audit tests the strongest simple version suggested by
`alpha_geo = 4 ln 2`: if Shannon entropy alone forced the geometry, then a
fixed entropy value should not allow mutually different finite geometries on the
same support.  On a 16-point support, the uniform distribution has
H = ln(16) = 4 ln 2, but many inequivalent graph geometries carry exactly that
same probability entropy.  We compute graph-distance witnesses for four such
geometries.  Different distance histograms/diameters at identical entropy give
a finite obstruction to entropy-alone geometry forcing.
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
P2772 = GEN / "p2772_s1722_self_learning_kernel_update_law_stationarity_witness.json"
OUT = GEN / "p2773_s1723_shannon_entropy_geometry_forcing_obstruction.json"
MD = GEN / "p2773_s1723_shannon_entropy_geometry_forcing_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

N = 16
ALPHA_GEO = 4.0 * math.log(2.0)
NEGATIVE_EXPORT_FLAGS = [
    "shannon_entropy_forces_geometry_exported",
    "canonical_geometry_from_entropy_exported",
    "entropy_source_closes_kernel_geometry",
    "kernel_fully_expresses_nadsoliton_characteristics",
    "self_learning_kernel_update_theorem_exported",
    "role_bearing_ltotal_promoted",
    "bridge_closure_exported",
    "toe_closure_exported",
]


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"missing": True, "path": rel(path)}


def shannon_entropy(probabilities: list[float]) -> float:
    return -sum(p * math.log(p) for p in probabilities if p > 0.0)


def graph_edges(kind: str, n: int = N) -> list[tuple[int, int]]:
    if kind == "complete":
        return [(i, j) for i in range(n) for j in range(i + 1, n)]
    if kind == "cycle":
        return [(i, (i + 1) % n) for i in range(n)]
    if kind == "path":
        return [(i, i + 1) for i in range(n - 1)]
    if kind == "star":
        return [(0, i) for i in range(1, n)]
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


def graph_metric_summary(kind: str) -> dict[str, Any]:
    edges = graph_edges(kind)
    distances = all_pairs_distances(edges)
    hist = distance_histogram(distances)
    pair_count = sum(hist.values())
    diameter = max(hist)
    average_distance = sum(distance * count for distance, count in hist.items()) / pair_count
    return {
        "geometry": kind,
        "node_count": N,
        "edge_count": len(edges),
        "diameter": diameter,
        "average_pair_distance": average_distance,
        "distance_histogram": {str(key): value for key, value in hist.items()},
    }


def entropy_geometry_witness() -> dict[str, Any]:
    probabilities = [1.0 / N] * N
    entropy = shannon_entropy(probabilities)
    rows = [graph_metric_summary(kind) for kind in ["complete", "cycle", "path", "star"]]
    signatures = {json.dumps(row["distance_histogram"], sort_keys=True) for row in rows}
    return {
        "entropy_object": "uniform distribution on 16 support points",
        "support_size": N,
        "shannon_entropy_nats": entropy,
        "alpha_geo_4_ln_2": ALPHA_GEO,
        "entropy_matches_alpha_geo": abs(entropy - ALPHA_GEO) < 1e-12,
        "geometry_rows": rows,
        "geometry_row_count": len(rows),
        "distinct_distance_histogram_count": len(signatures),
        "all_rows_same_entropy": all(abs(entropy - ALPHA_GEO) < 1e-12 for _ in rows),
        "entropy_forces_unique_geometry_on_this_class": len(signatures) == 1,
        "finite_obstruction_statement": "The same H=4 ln 2 probability entropy supports multiple inequivalent finite graph geometries, so Shannon entropy alone does not force a unique nadsoliton geometry without extra geometric/source data.",
    }


def acceptance_matrix(witness: dict[str, Any], p2772: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "p2772_no_closure_boundary_present": p2772.get("status") == "P2772_SELF_LEARNING_KERNEL_UPDATE_LAW_STATIONARITY_WITNESS_BOUNDED_NO_GO_NO_CLOSURE",
        "entropy_matches_alpha_geo": witness["entropy_matches_alpha_geo"],
        "multiple_inequivalent_geometries_at_same_entropy": witness["distinct_distance_histogram_count"] > 1,
        "entropy_forces_unique_geometry_on_test_class": witness["entropy_forces_unique_geometry_on_this_class"],
        "additional_geometric_source_supplied": False,
        "kernel_geometry_closure_exported": False,
    }
    return {
        "facts": facts,
        "accepted_as_shannon_entropy_geometry_forcing_theorem": False,
        "accepted_as_canonical_nadsoliton_geometry_source": False,
        "accepted_as_kernel_full_expression_theorem": False,
        "accepted_as_ltotal_or_toe_promotion": False,
        "missing_criteria": [key for key, value in facts.items() if not value],
        "blocker": "Although H=ln(16)=4 ln 2 matches alpha_geo, identical Shannon entropy is carried by several inequivalent finite geometries.  Entropy can normalize/count information, but it does not by itself select the metric, adjacency, or geometric self-coupling law.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    witness = payload["entropy_geometry_witness"]
    lines = [
        "# P2773/S1723 Shannon-entropy geometry-forcing obstruction",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Entropy witness",
        f"- H_uniform_16={witness['shannon_entropy_nats']}",
        f"- alpha_geo_4_ln_2={witness['alpha_geo_4_ln_2']}",
        f"- entropy_matches_alpha_geo={witness['entropy_matches_alpha_geo']}",
        "",
        "## Geometry obstruction",
        f"- geometry_row_count={witness['geometry_row_count']}",
        f"- distinct_distance_histogram_count={witness['distinct_distance_histogram_count']}",
        f"- entropy_forces_unique_geometry_on_this_class={witness['entropy_forces_unique_geometry_on_this_class']}",
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
    p2772 = read_json(P2772)
    witness = entropy_geometry_witness()
    acceptance = acceptance_matrix(witness, p2772)
    payload = {
        "status": "P2773_SHANNON_ENTROPY_GEOMETRY_FORCING_OBSTRUCTION_NO_CLOSURE",
        "input_hashes": {"P2772": sha(P2772)},
        "input_statuses": {"P2772": p2772.get("status")},
        "audited_question": "Does Shannon entropy alone force the nadsoliton geometry?",
        "entropy_geometry_witness": witness,
        "acceptance_matrix": acceptance,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": acceptance["blocker"],
            "next_honest_step": "Do not claim that Shannon entropy alone forces the nadsoliton geometry.  The next honest move must add one extra geometric/source principle beyond H=4 ln 2: for example an entropy-plus-transport variational principle, an entropy-plus-graph-Laplacian minimizer, or a strict nadsoliton adjacency/metric source theorem.  Then rerun a degeneracy-breaking test against the complete/cycle/path/star witness class; otherwise preserve the P2697-P2773 no-geometry-forcing/no-closure certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2773/S1723 Shannon-entropy geometry-forcing obstruction", "## P2773/S1723 Shannon-entropy geometry-forcing obstruction\n\n`P2773/S1723` tests whether Shannon entropy alone can force nadsoliton geometry.  The uniform distribution on 16 points has `H=ln(16)=4 ln 2=alpha_geo`, but complete, cycle, path, and star graph geometries on the same support have inequivalent distance histograms and diameters at the same entropy.  Thus Shannon entropy can support the information-count normalization, but it does not by itself select the metric/adjacency/geometric self-coupling law.  No canonical geometry source, kernel full-expression theorem, role-bearing `L_total`, bridge closure, role transfer, selector closure, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2773/S1723 entropy-geometry Ltotal guard", "## P2773/S1723 entropy-geometry Ltotal guard\n\n`P2773/S1723` adds no variational source term.  It shows that `H=4 ln 2` is compatible with multiple inequivalent finite geometries, so Shannon entropy alone cannot promote the kernels to role-bearing `L_total` or a canonical nadsoliton geometry.\n")
    append_once(AGENTS, "Current Shannon-entropy geometry-forcing obstruction guardrail (P2773/S1723, 2026-06-15)", "## Current Shannon-entropy geometry-forcing obstruction guardrail (P2773/S1723, 2026-06-15)\n\n- P2773 tests whether Shannon entropy alone forces nadsoliton geometry: the uniform 16-point distribution has `H=ln(16)=4 ln 2=alpha_geo`.\n- Complete, cycle, path, and star geometries on the same 16-point support have the same Shannon entropy but inequivalent distance histograms/diameters, so entropy alone does not select metric/adjacency/geometric self-coupling.\n- Do not promote `H=4 ln 2` to canonical nadsoliton geometry, kernel full-expression, self-learning closure, role-bearing `L_total`, bridge closure, role transfer, selector closure, or ToE closure.  A next admissible move must add one extra geometric/source principle and pass a degeneracy-breaking test against this witness class.\n")
    return payload


if __name__ == "__main__":
    main()
