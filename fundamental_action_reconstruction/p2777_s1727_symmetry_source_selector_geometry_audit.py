#!/usr/bin/env python3
"""P2777/S1727: symmetry-source selector/geometry audit.

This bounded audit tests the user's proposed next source: symmetry.  It checks
both live meanings separately:

1. selector source: whether Aut(Z12)-symmetry can choose +1 over -1;
2. geometry source: whether graph symmetry can select the P2774/P2775 geometry.

The result is mixed and guarded.  Aut(Z12) symmetry cannot select orientation
because inversion puts +1 and -1 in one orbit.  On the P2774 geometry pair,
automorphism-group size distinguishes torus_4x4 from circulant_pm1_pm2, but
mere pair-local maximal symmetry is not a strict sourced geometry theorem.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2774_s1724_entropy_laplacian_trace_geometry_degeneracy import N, graph_edges

GEN = ROOT / "generated"
P2776 = GEN / "p2776_s1726_small_graph_full_spectrum_uniqueness_audit.json"
OUT = GEN / "p2777_s1727_symmetry_source_selector_geometry_audit.json"
MD = GEN / "p2777_s1727_symmetry_source_selector_geometry_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

NEGATIVE_EXPORT_FLAGS = [
    "aut_z12_selector_source_exported",
    "canonical_geometry_source_exported",
    "max_symmetry_geometry_theorem_exported",
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


def aut_z12_units() -> list[int]:
    return [unit for unit in range(12) if gcd(unit, 12) == 1]


def gcd(a: int, b: int) -> int:
    while b:
        a, b = b, a % b
    return abs(a)


def z12_selector_witness() -> dict[str, Any]:
    units = aut_z12_units()
    directed_units = [1, 11]
    orbit = sorted({(unit * value) % 12 for unit in units for value in directed_units})
    singleton_selectors = [{1}, {11}]
    invariant_singletons = [sorted(sel) for sel in singleton_selectors if all(((unit * next(iter(sel))) % 12) in sel for unit in units)]
    return {
        "group": "Aut(Z12)=U(12)",
        "aut_units": units,
        "directed_selector_candidates": {"+1": 1, "-1": 11},
        "orbit_of_directed_candidates": orbit,
        "inversion_unit_present": 11 in units,
        "plus_maps_to_minus_under_inversion": (11 * 1) % 12 == 11,
        "minus_maps_to_plus_under_inversion": (11 * 11) % 12 == 1,
        "aut_invariant_singleton_selector_count": len(invariant_singletons),
        "aut_invariant_singleton_selectors": invariant_singletons,
        "symmetry_selects_orientation": len(invariant_singletons) == 1,
        "finite_obstruction_statement": "Aut(Z12) symmetry places +1 and -1 in the same directed-unit orbit, so symmetry alone cannot provide a non-premise selector orientation.",
    }


def adjacency_sets(kind: str) -> list[set[int]]:
    adj = [set() for _ in range(N)]
    for a, b in graph_edges(kind):
        adj[a].add(b)
        adj[b].add(a)
    return adj


def automorphisms(kind: str) -> list[dict[int, int]]:
    """Enumerate automorphisms of the 16-node P2774 graph by backtracking."""
    adj = adjacency_sets(kind)
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


def geometry_symmetry_row(kind: str) -> dict[str, Any]:
    autos = automorphisms(kind)
    orbit_sizes = []
    for node in range(N):
        orbit_sizes.append(len({auto[node] for auto in autos}))
    return {
        "geometry": kind,
        "node_count": N,
        "edge_count": len(graph_edges(kind)),
        "automorphism_count": len(autos),
        "orbit_sizes": orbit_sizes,
        "vertex_transitive": all(size == N for size in orbit_sizes),
    }


def geometry_symmetry_witness() -> dict[str, Any]:
    rows = [geometry_symmetry_row("torus_4x4"), geometry_symmetry_row("circulant_pm1_pm2")]
    max_count = max(row["automorphism_count"] for row in rows)
    max_rows = [row["geometry"] for row in rows if row["automorphism_count"] == max_count]
    vertex_transitive_rows = [row["geometry"] for row in rows if row["vertex_transitive"]]
    return {
        "audited_pair": [row["geometry"] for row in rows],
        "geometry_rows": rows,
        "all_pair_geometries_vertex_transitive": len(vertex_transitive_rows) == len(rows),
        "vertex_transitivity_selects_unique_geometry": len(vertex_transitive_rows) == 1,
        "max_automorphism_count": max_count,
        "max_automorphism_count_geometries": max_rows,
        "max_symmetry_selects_unique_geometry_on_pair": len(max_rows) == 1,
        "finite_positive_statement": "Automorphism-group size selects torus_4x4 over circulant_pm1_pm2 on this pair, while vertex-transitivity alone does not select because both graphs are vertex-transitive.",
    }


def combined_witness() -> dict[str, Any]:
    selector = z12_selector_witness()
    geometry = geometry_symmetry_witness()
    return {
        "selector_symmetry_witness": selector,
        "geometry_symmetry_witness": geometry,
        "symmetry_as_selector_source_passes": selector["symmetry_selects_orientation"],
        "symmetry_as_pair_local_geometry_discriminator_passes": geometry["max_symmetry_selects_unique_geometry_on_pair"],
    }


def acceptance_matrix(witness: dict[str, Any], p2776: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "p2776_small_class_spectral_boundary_present": p2776.get("status") == "P2776_SMALL_GRAPH_FULL_SPECTRUM_UNIQUENESS_AUDIT_NO_CLOSURE",
        "aut_z12_symmetry_selects_orientation": witness["symmetry_as_selector_source_passes"],
        "vertex_transitivity_selects_p2774_geometry": witness["geometry_symmetry_witness"]["vertex_transitivity_selects_unique_geometry"],
        "max_automorphism_count_selects_p2774_pair_geometry": witness["symmetry_as_pair_local_geometry_discriminator_passes"],
        "strict_source_law_for_max_symmetry_exported": False,
        "global_strict_graph_class_symmetry_uniqueness_audited": False,
        "kernel_or_ltotal_variational_coupling_exported": False,
    }
    return {
        "facts": facts,
        "accepted_as_selector_source_theorem": False,
        "accepted_as_pair_local_max_symmetry_geometry_discriminator": facts["max_automorphism_count_selects_p2774_pair_geometry"],
        "accepted_as_canonical_nadsoliton_geometry_source": False,
        "accepted_as_kernel_full_expression_theorem": False,
        "accepted_as_ltotal_or_toe_promotion": False,
        "missing_criteria": [key for key, value in facts.items() if not value],
        "blocker": "Aut(Z12) symmetry cannot select +1 over -1.  Maximal automorphism count distinguishes the P2774 pair, but no strict law says to maximize automorphism count, no global strict graph class is audited, and no K/L_total variational coupling is exported.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    witness = payload["symmetry_source_witness"]
    selector = witness["selector_symmetry_witness"]
    geometry = witness["geometry_symmetry_witness"]
    lines = [
        "# P2777/S1727 symmetry-source selector/geometry audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Selector result",
        f"- orbit_of_directed_candidates={selector['orbit_of_directed_candidates']}",
        f"- aut_invariant_singleton_selector_count={selector['aut_invariant_singleton_selector_count']}",
        f"- symmetry_selects_orientation={selector['symmetry_selects_orientation']}",
        "",
        "## Geometry result",
        f"- all_pair_geometries_vertex_transitive={geometry['all_pair_geometries_vertex_transitive']}",
        f"- max_automorphism_count_geometries={geometry['max_automorphism_count_geometries']}",
        f"- max_symmetry_selects_unique_geometry_on_pair={geometry['max_symmetry_selects_unique_geometry_on_pair']}",
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
    p2776 = read_json(P2776)
    witness = combined_witness()
    acceptance = acceptance_matrix(witness, p2776)
    payload = {
        "status": "P2777_SYMMETRY_SOURCE_SELECTOR_GEOMETRY_AUDIT_NO_CLOSURE",
        "input_hashes": {"P2776": sha(P2776)},
        "input_statuses": {"P2776": p2776.get("status")},
        "audited_question": "Can symmetry source either the selector orientation or canonical nadsoliton geometry?",
        "symmetry_source_witness": witness,
        "acceptance_matrix": acceptance,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": acceptance["blocker"],
            "next_honest_step": "Do not claim selector or geometry closure from symmetry alone.  The next honest move is exactly one of: supply a strict source law that justifies a concrete symmetry functional such as maximal automorphism count and audit it over a declared 16-point/strict graph class, or supply a symmetry-breaking/chiral/pseudoscalar source that couples to the Aut(Z12) +1/-1 selector orbit.  Without one of those, preserve the P2697-P2777 no-selector/no-canonical-geometry/no-closure certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2777/S1727 symmetry-source selector/geometry audit", "## P2777/S1727 symmetry-source selector/geometry audit\n\n`P2777/S1727` tests symmetry as a source for both selector orientation and geometry.  On `Aut(Z12)`, inversion sends `+1` to `-1` and the full action sends the directed candidates across the unit orbit, so symmetry alone leaves zero invariant singleton selectors.  On the P2774 pair, both geometries are vertex-transitive, while automorphism-group size distinguishes `torus_4x4` from `circulant_pm1_pm2`; this is only a pair-local maximal-symmetry discriminator.  No strict law choosing maximal automorphism count, global strict graph-class theorem, kernel full-expression theorem, role-bearing `L_total`, bridge closure, role transfer, selector closure, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2777/S1727 symmetry-source Ltotal guard", "## P2777/S1727 symmetry-source Ltotal guard\n\n`P2777/S1727` adds no variational source term.  Symmetry does not fix the `Aut(Z12)` selector orientation; maximal automorphism count distinguishes only the P2774 pair and is not sourced as a variational principle coupled to `K`/`L_total`.  Therefore it cannot promote role-bearing `L_total`, selector closure, or canonical nadsoliton geometry.\n")
    append_once(AGENTS, "Current symmetry-source selector/geometry audit guardrail (P2777/S1727, 2026-06-15)", "## Current symmetry-source selector/geometry audit guardrail (P2777/S1727, 2026-06-15)\n\n- P2777 tests symmetry as a source for selector orientation and geometry.  `Aut(Z12)` inversion maps `+1` to `-1` and the full action spreads directed candidates across the unit orbit, leaving no invariant singleton selector.\n- On the P2774 pair, both geometries are vertex-transitive; automorphism-group size selects `torus_4x4` over `circulant_pm1_pm2` only as a pair-local maximal-symmetry discriminator.\n- Do not promote symmetry alone to `QW-2191` discharge, selector closure, canonical nadsoliton geometry, kernel full-expression, role-bearing `L_total`, bridge closure, role transfer, or ToE closure.  A next admissible move must supply a strict source law for a concrete symmetry functional over a declared graph class, or a genuine symmetry-breaking/chiral source coupled to the `Aut(Z12)` selector orbit.\n")
    return payload


if __name__ == "__main__":
    main()
