#!/usr/bin/env python3
"""P2778/S1728: max-symmetry 16-node geometry source audit.

P2777 found that automorphism-group size distinguishes the P2774 pair, but also
warned that no strict law says to maximize automorphism count.  This follow-up
supplies exactly that bounded candidate law on a declared 16-node class: connected
4-regular circulant/Cayley graphs on Z16, plus the P2774 torus_4x4 reference.

The candidate law does not support canonical nadsoliton geometry.  It selects
highly symmetric non-torus circulants instead of torus_4x4, and the selected
connection labels are not unique before quotienting.  Thus max-symmetry is not a
safe geometry source on this declared finite class.
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
from p2774_s1724_entropy_laplacian_trace_geometry_degeneracy import graph_edges

GEN = ROOT / "generated"
P2777 = GEN / "p2777_s1727_symmetry_source_selector_geometry_audit.json"
OUT = GEN / "p2778_s1728_max_symmetry_16node_geometry_source_audit.json"
MD = GEN / "p2778_s1728_max_symmetry_16node_geometry_source_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

N = 16
NEGATIVE_EXPORT_FLAGS = [
    "max_symmetry_geometry_source_exported",
    "canonical_geometry_source_exported",
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


def normalized_edges(edges: set[tuple[int, int]]) -> list[tuple[int, int]]:
    return sorted({(min(a, b), max(a, b)) for a, b in edges if a != b})


def circulant_edges(a: int, b: int) -> list[tuple[int, int]]:
    edges: set[tuple[int, int]] = set()
    for node in range(N):
        for step in (a, b):
            edges.add((node, (node + step) % N))
            edges.add((node, (node - step) % N))
    return normalized_edges(edges)


def candidate_edge_sets() -> list[dict[str, Any]]:
    rows = [{"geometry": "torus_4x4", "family": "p2774_reference", "steps": None, "edges": graph_edges("torus_4x4")}]
    for a in range(1, 8):
        for b in range(a + 1, 8):
            if math.gcd(math.gcd(N, a), b) != 1:
                continue
            rows.append({"geometry": f"circulant_pm{a}_pm{b}", "family": "connected_4_regular_circulant_z16", "steps": [a, b], "edges": circulant_edges(a, b)})
    return rows


def adjacency(edges: list[tuple[int, int]]) -> list[set[int]]:
    adj = [set() for _ in range(N)]
    for a, b in edges:
        adj[a].add(b)
        adj[b].add(a)
    return adj


def is_connected(edges: list[tuple[int, int]]) -> bool:
    adj = adjacency(edges)
    seen = {0}
    queue: deque[int] = deque([0])
    while queue:
        node = queue.popleft()
        for nxt in adj[node]:
            if nxt not in seen:
                seen.add(nxt)
                queue.append(nxt)
    return len(seen) == N


def automorphisms_from_edges(edges: list[tuple[int, int]]) -> list[dict[int, int]]:
    adj = adjacency(edges)
    degrees = [len(row) for row in adj]
    order = sorted(range(N), key=lambda node: (-degrees[node], node))
    candidates_by_degree: dict[int, list[int]] = {}
    for node, degree in enumerate(degrees):
        candidates_by_degree.setdefault(degree, []).append(node)
    mapping: dict[int, int] = {}
    used: set[int] = set()
    out: list[dict[int, int]] = []

    def compatible(source: int, target: int) -> bool:
        for mapped_source, mapped_target in mapping.items():
            if (mapped_source in adj[source]) != (mapped_target in adj[target]):
                return False
        return True

    def rec(index: int) -> None:
        if index == N:
            out.append(dict(mapping))
            return
        source = order[index]
        for target in candidates_by_degree[degrees[source]]:
            if target in used:
                continue
            if not compatible(source, target):
                continue
            mapping[source] = target
            used.add(target)
            rec(index + 1)
            used.remove(target)
            del mapping[source]

    rec(0)
    return out


def row_summary(candidate: dict[str, Any]) -> dict[str, Any]:
    edges = candidate["edges"]
    adj = adjacency(edges)
    autos = automorphisms_from_edges(edges)
    orbit_sizes = [len({auto[node] for auto in autos}) for node in range(N)]
    degree_sequence = sorted(len(row) for row in adj)
    return {
        "geometry": candidate["geometry"],
        "family": candidate["family"],
        "steps": candidate["steps"],
        "node_count": N,
        "edge_count": len(edges),
        "connected": is_connected(edges),
        "degree_sequence": degree_sequence,
        "is_4_regular": all(degree == 4 for degree in degree_sequence),
        "automorphism_count": len(autos),
        "vertex_transitive": all(size == N for size in orbit_sizes),
    }


def max_symmetry_witness() -> dict[str, Any]:
    rows = [row_summary(candidate) for candidate in candidate_edge_sets()]
    max_count = max(row["automorphism_count"] for row in rows)
    max_rows = [row for row in rows if row["automorphism_count"] == max_count]
    torus_row = next(row for row in rows if row["geometry"] == "torus_4x4")
    return {
        "candidate_law": "choose the connected 16-node 4-regular candidate with maximal automorphism-group size",
        "candidate_class": "torus_4x4 plus all connected 4-regular circulant Cayley graphs C16({±a,±b}) with 1<=a<b<=7 and gcd(16,a,b)=1",
        "candidate_count": len(rows),
        "all_candidates_connected_4_regular": all(row["connected"] and row["is_4_regular"] for row in rows),
        "all_candidates_vertex_transitive": all(row["vertex_transitive"] for row in rows),
        "rows": sorted(rows, key=lambda row: (-row["automorphism_count"], row["geometry"])),
        "max_automorphism_count": max_count,
        "max_geometry_names": [row["geometry"] for row in max_rows],
        "max_geometry_count": len(max_rows),
        "torus_4x4_automorphism_count": torus_row["automorphism_count"],
        "max_symmetry_selects_torus_4x4": any(row["geometry"] == "torus_4x4" for row in max_rows),
        "max_symmetry_selects_unique_labeled_geometry": len(max_rows) == 1,
        "finite_obstruction_statement": "On this declared 16-node class, maximal automorphism count selects non-torus circulant labels rather than torus_4x4, so max-symmetry does not source the intended geometry.",
    }


def acceptance_matrix(witness: dict[str, Any], p2777: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "p2777_symmetry_boundary_present": p2777.get("status") == "P2777_SYMMETRY_SOURCE_SELECTOR_GEOMETRY_AUDIT_NO_CLOSURE",
        "strict_candidate_law_supplied": True,
        "declared_16node_class_audited": True,
        "all_candidates_connected_4_regular": witness["all_candidates_connected_4_regular"],
        "max_symmetry_selects_unique_labeled_geometry": witness["max_symmetry_selects_unique_labeled_geometry"],
        "max_symmetry_selects_torus_4x4": witness["max_symmetry_selects_torus_4x4"],
        "kernel_or_ltotal_variational_coupling_exported": False,
    }
    return {
        "facts": facts,
        "accepted_as_max_symmetry_geometry_source": False,
        "accepted_as_canonical_nadsoliton_geometry_source": False,
        "accepted_as_kernel_full_expression_theorem": False,
        "accepted_as_ltotal_or_toe_promotion": False,
        "missing_criteria": [key for key, value in facts.items() if not value],
        "blocker": "The supplied max-automorphism candidate law fails on the declared 16-node class: it does not select torus_4x4 and is not unique on labeled circulant candidates.  No K/L_total variational coupling is exported.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    witness = payload["max_symmetry_16node_witness"]
    lines = [
        "# P2778/S1728 max-symmetry 16-node geometry source audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite class result",
        f"- candidate_count={witness['candidate_count']}",
        f"- max_automorphism_count={witness['max_automorphism_count']}",
        f"- max_geometry_names={witness['max_geometry_names']}",
        f"- torus_4x4_automorphism_count={witness['torus_4x4_automorphism_count']}",
        f"- max_symmetry_selects_torus_4x4={witness['max_symmetry_selects_torus_4x4']}",
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
    p2777 = read_json(P2777)
    witness = max_symmetry_witness()
    acceptance = acceptance_matrix(witness, p2777)
    payload = {
        "status": "P2778_MAX_SYMMETRY_16NODE_GEOMETRY_SOURCE_AUDIT_NO_CLOSURE",
        "input_hashes": {"P2777": sha(P2777)},
        "input_statuses": {"P2777": p2777.get("status")},
        "audited_question": "Does a concrete maximal-automorphism symmetry law source the intended 16-node geometry?",
        "max_symmetry_16node_witness": witness,
        "acceptance_matrix": acceptance,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": acceptance["blocker"],
            "next_honest_step": "Do not use maximal automorphism count as the geometry source.  The next honest move is either to test a different sourced symmetry functional with a declared target and quotient rule, or pivot from symmetry to a genuine strict metric/variational source such as a spectral action with fixed target spectrum and K/L_total coupling.  Otherwise preserve the P2697-P2778 no-canonical-geometry/no-closure certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2778/S1728 max-symmetry 16-node geometry source audit", "## P2778/S1728 max-symmetry 16-node geometry source audit\n\n`P2778/S1728` supplies the concrete symmetry-source law left open by P2777: maximize automorphism-group size on a declared 16-node connected 4-regular class consisting of `torus_4x4` plus connected circulant Cayley graphs `C16({±a,±b})`.  The law fails as a source for the intended geometry: the maximum automorphism count is attained by non-torus circulant labels, not by `torus_4x4`, and no `K`/`L_total` variational coupling is exported.  No canonical nadsoliton geometry, kernel full-expression theorem, role-bearing `L_total`, bridge closure, role transfer, selector closure, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2778/S1728 max-symmetry Ltotal guard", "## P2778/S1728 max-symmetry Ltotal guard\n\n`P2778/S1728` adds no variational source term.  A concrete maximal-automorphism law on a declared 16-node class selects non-torus circulant labels rather than a sourced `K`/`L_total` geometry, so it cannot promote role-bearing `L_total` or canonical nadsoliton geometry.\n")
    append_once(AGENTS, "Current max-symmetry 16-node geometry source audit guardrail (P2778/S1728, 2026-06-15)", "## Current max-symmetry 16-node geometry source audit guardrail (P2778/S1728, 2026-06-15)\n\n- P2778 tests the concrete symmetry-source law left open by P2777: maximize automorphism-group size on `torus_4x4` plus connected 16-node 4-regular circulant Cayley candidates.\n- The law fails to source the intended geometry: maximum automorphism count is attained by non-torus circulant labels, not by `torus_4x4`, and no `K`/`L_total` variational coupling is exported.\n- Do not promote maximal symmetry to canonical nadsoliton geometry, selector closure, kernel full-expression, role-bearing `L_total`, bridge closure, role transfer, or ToE closure.  A next admissible move must supply a different sourced symmetry functional with target/quotient rule, or pivot to a genuine strict metric/variational source.\n")
    return payload


if __name__ == "__main__":
    main()
