#!/usr/bin/env python3
"""P2789/S1739: orbit-stabilizer exact quotient certificate.

P2788 checked complement-duality identities for the already available spectral
witnesses.  This follow-up adds a more proof-style finite group certificate:
compute exact graph automorphism group sizes and apply orbit-stabilizer checks.

For the complete P2787 8-node 4-regular quotient, the script verifies that each
stored labeled-member count equals ``8! / |Aut(G)|`` and that the six orbit sizes
sum to the known 19,355 labeled connected candidates.  For the seven local
P2786 16-node representatives, it computes exact stabilizer sizes and implied
labeled orbit sizes under S_16 without pretending that these seven orbits exhaust
the full 16-node class.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2786_s1736_graph6_provenance_toolchain_gate import graph6_decode as graph6_decode_16
from p2788_s1738_complement_duality_exact_spectral_certificate import graph6_decode_small, normalize_edges

GEN = ROOT / "generated"
P2786 = GEN / "p2786_s1736_graph6_provenance_toolchain_gate.json"
P2787 = GEN / "p2787_s1737_small_canonical_generator_pipeline_audit.json"
P2788 = GEN / "p2788_s1738_complement_duality_exact_spectral_certificate.json"
OUT = GEN / "p2789_s1739_orbit_stabilizer_exact_quotient_certificate.json"
MD = GEN / "p2789_s1739_orbit_stabilizer_exact_quotient_certificate.md"
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


def adjacency_sets(edges: tuple[tuple[int, int], ...], n: int) -> list[set[int]]:
    adj = [set() for _ in range(n)]
    for a, b in edges:
        adj[a].add(b)
        adj[b].add(a)
    return adj


def vertex_invariants(adj: list[set[int]]) -> list[tuple[int, int, int, tuple[int, ...]]]:
    """Cheap isomorphism-safe invariants used only for pruning automorphism search."""
    out = []
    for vertex, neighbours in enumerate(adj):
        triangle_count = sum(1 for a in neighbours for b in neighbours if a < b and b in adj[a])
        distance_two = set()
        for neighbour in neighbours:
            distance_two.update(adj[neighbour])
        distance_two.discard(vertex)
        distance_two.difference_update(neighbours)
        neighbour_triangle_profile = tuple(sorted(sum(1 for a in adj[n] for b in adj[n] if a < b and b in adj[a]) for n in neighbours))
        out.append((len(neighbours), triangle_count, len(distance_two), neighbour_triangle_profile))
    return out


def automorphism_count(edges: tuple[tuple[int, int], ...], n: int) -> int:
    """Count automorphisms by exact backtracking with adjacency/non-adjacency checks."""
    adj = adjacency_sets(edges, n)
    invariants = vertex_invariants(adj)
    candidate_targets = {source: [target for target in range(n) if invariants[target] == invariants[source]] for source in range(n)}
    preferred_order = sorted(range(n), key=lambda source: (len(candidate_targets[source]), invariants[source], source))
    mapping: dict[int, int] = {}
    used: set[int] = set()
    count = 0

    def compatible(source: int, target: int) -> bool:
        return all((mapped_source in adj[source]) == (mapped_target in adj[target]) for mapped_source, mapped_target in mapping.items())

    def recurse() -> None:
        nonlocal count
        if len(mapping) == n:
            count += 1
            return
        best_source = None
        best_options: list[int] | None = None
        for source in preferred_order:
            if source in mapping:
                continue
            options = [target for target in candidate_targets[source] if target not in used and compatible(source, target)]
            if best_options is None or len(options) < len(best_options):
                best_source = source
                best_options = options
            if best_options == []:
                break
        if best_source is None or not best_options:
            return
        for target in best_options:
            mapping[best_source] = target
            used.add(target)
            recurse()
            used.remove(target)
            del mapping[best_source]

    recurse()
    return count


def small_orbit_rows(p2787: dict[str, Any]) -> list[dict[str, Any]]:
    rows = []
    for row in p2787.get("small_pipeline_witness", {}).get("representative_rows", []):
        edges = graph6_decode_small(row["graph6"], 8)
        aut_size = automorphism_count(edges, 8)
        orbit_size = math.factorial(8) // aut_size
        rows.append({
            "representative_index": row["representative_index"],
            "n": 8,
            "edge_count": len(edges),
            "automorphism_group_size": aut_size,
            "orbit_size_by_orbit_stabilizer": orbit_size,
            "stored_labeled_member_count": row["labeled_member_count"],
            "orbit_stabilizer_matches_stored_member_count": orbit_size == row["labeled_member_count"],
        })
    return rows


def local_16_orbit_rows(p2786: dict[str, Any]) -> list[dict[str, Any]]:
    rows = []
    for row in p2786.get("graph6_provenance_witness", {}).get("provenance_rows", []):
        edges = normalize_edges(graph6_decode_16(row["graph6"]))
        aut_size = automorphism_count(edges, 16)
        orbit_size = math.factorial(16) // aut_size
        rows.append({
            "representative": row["representative"],
            "n": 16,
            "edge_count": len(edges),
            "automorphism_group_size": aut_size,
            "orbit_size_by_orbit_stabilizer": orbit_size,
            "member_count_in_local_quotient": row["member_count"],
            "scope_note": "This is the exact S_16 orbit size of one local representative, not a full-class enumeration count.",
        })
    return rows


def orbit_stabilizer_witness(p2786: dict[str, Any], p2787: dict[str, Any]) -> dict[str, Any]:
    small = small_orbit_rows(p2787)
    local = local_16_orbit_rows(p2786)
    small_orbit_sum = sum(row["orbit_size_by_orbit_stabilizer"] for row in small)
    stored_small_total = p2787.get("small_pipeline_witness", {}).get("connected_labeled_candidate_count")
    return {
        "source_class": "Exact automorphism-group/orbit-stabilizer audit over the complete P2787 8-node quotient and seven local P2786 16-node representatives.",
        "small_complete_8node_row_count": len(small),
        "local_16node_row_count": len(local),
        "small_8node_rows": small,
        "local_16node_rows": local,
        "small_orbit_size_sum": small_orbit_sum,
        "stored_small_connected_labeled_candidate_count": stored_small_total,
        "all_small_orbit_stabilizer_counts_match_stored_members": all(row["orbit_stabilizer_matches_stored_member_count"] for row in small),
        "small_orbit_sum_matches_stored_connected_labeled_total": small_orbit_sum == stored_small_total == 19355,
        "all_local_16node_stabilizers_positive": all(row["automorphism_group_size"] > 0 for row in local),
        "finite_certificate_statement": "The complete P2787 8-node quotient satisfies exact orbit-stabilizer member counts summing to 19,355; the seven local 16-node representatives have exact stabilizer/orbit sizes but do not exhaust the full class.",
    }


def acceptance_matrix(witness: dict[str, Any], p2786: dict[str, Any], p2787: dict[str, Any], p2788: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "p2786_graph6_gate_present": p2786.get("status") == "P2786_GRAPH6_PROVENANCE_TOOLCHAIN_GATE_NO_CLOSURE",
        "p2787_small_pipeline_present": p2787.get("status") == "P2787_SMALL_CANONICAL_GENERATOR_PIPELINE_AUDIT_NO_CLOSURE",
        "p2788_complement_certificate_present": p2788.get("status") == "P2788_COMPLEMENT_DUALITY_EXACT_SPECTRAL_CERTIFICATE_NO_CLOSURE",
        "complete_small_quotient_orbits_checked": witness["small_complete_8node_row_count"] == 6,
        "seven_local_16node_stabilizers_checked": witness["local_16node_row_count"] == 7,
        "small_orbit_stabilizer_counts_match": witness["all_small_orbit_stabilizer_counts_match_stored_members"],
        "small_orbit_sum_matches_19355": witness["small_orbit_sum_matches_stored_connected_labeled_total"],
        "local_16node_stabilizers_positive": witness["all_local_16node_stabilizers_positive"],
        "canonical_16node_generator_certified": False,
        "strict_nadsoliton_spectral_source_law_exported": False,
        "kernel_or_ltotal_variational_coupling_exported": False,
    }
    return {
        "facts": facts,
        "accepted_as_exact_orbit_stabilizer_quotient_certificate": all(facts[key] for key in [
            "p2786_graph6_gate_present",
            "p2787_small_pipeline_present",
            "p2788_complement_certificate_present",
            "complete_small_quotient_orbits_checked",
            "seven_local_16node_stabilizers_checked",
            "small_orbit_stabilizer_counts_match",
            "small_orbit_sum_matches_19355",
            "local_16node_stabilizers_positive",
        ]),
        "accepted_as_full_16node_canonical_generator_certificate": False,
        "accepted_as_strict_spectral_source_law": False,
        "accepted_as_canonical_nadsoliton_geometry_source": False,
        "accepted_as_ltotal_or_toe_promotion": False,
        "missing_criteria": [key for key, value in facts.items() if not value],
        "blocker": "Orbit-stabilizer validates quotient/member arithmetic for the complete 8-node class and computes exact stabilizers for seven local 16-node representatives, but it still does not enumerate the full connected 16-node 4-regular class or export a strict K/L_total spectral source law.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    w = payload["orbit_stabilizer_witness"]
    lines = [
        "# P2789/S1739 orbit-stabilizer exact quotient certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Exact orbit-stabilizer result",
        f"- small_complete_8node_row_count={w['small_complete_8node_row_count']}",
        f"- local_16node_row_count={w['local_16node_row_count']}",
        f"- small_orbit_size_sum={w['small_orbit_size_sum']}",
        f"- stored_small_connected_labeled_candidate_count={w['stored_small_connected_labeled_candidate_count']}",
        f"- all_small_orbit_stabilizer_counts_match_stored_members={w['all_small_orbit_stabilizer_counts_match_stored_members']}",
        f"- small_orbit_sum_matches_stored_connected_labeled_total={w['small_orbit_sum_matches_stored_connected_labeled_total']}",
        f"- all_local_16node_stabilizers_positive={w['all_local_16node_stabilizers_positive']}",
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
    p2788 = read_json(P2788)
    witness = orbit_stabilizer_witness(p2786, p2787)
    acceptance = acceptance_matrix(witness, p2786, p2787, p2788)
    payload = {
        "status": "P2789_ORBIT_STABILIZER_EXACT_QUOTIENT_CERTIFICATE_NO_CLOSURE",
        "input_hashes": {"P2786": sha(P2786), "P2787": sha(P2787), "P2788": sha(P2788)},
        "input_statuses": {"P2786": p2786.get("status"), "P2787": p2787.get("status"), "P2788": p2788.get("status")},
        "audited_question": "Do the stored quotient/member counts satisfy exact finite group orbit-stabilizer checks, and what are the exact stabilizers of the seven local 16-node representatives?",
        "orbit_stabilizer_witness": witness,
        "acceptance_matrix": acceptance,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": acceptance["blocker"],
            "next_honest_step": "Use P2789 only as an exact orbit-stabilizer quotient-arithmetic certificate.  The next honest move is still exactly one of: supply/import an actual certified full connected 16-node 4-regular generator artifact/toolchain with graph6/hash provenance and run exact quotient/charpoly/complement/orbit-stabilizer auditing there; or export a strict nadsoliton spectral action/source law fixing the admissible class, target spectrum, and K/L_total coupling before testing.  Otherwise preserve the P2697-P2789 no-canonical-geometry/no-closure certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2789/S1739 orbit-stabilizer exact quotient certificate", "## P2789/S1739 orbit-stabilizer exact quotient certificate\n\n`P2789/S1739` adds an exact finite-group quotient check.  For the complete P2787 8-node 4-regular quotient, automorphism-group sizes give orbit-stabilizer labeled-member counts that match every stored class and sum to 19,355 connected labeled candidates.  For the seven local P2786/P2785 16-node representatives, exact stabilizer sizes and implied `S_16` orbit sizes are recorded.  This validates quotient arithmetic and local stabilizers, but it is not a full connected 16-node 4-regular generator, not a strict spectral source law, and not a `K`/`L_total` variational coupling.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2789/S1739 orbit-stabilizer Ltotal guard", "## P2789/S1739 orbit-stabilizer Ltotal guard\n\n`P2789/S1739` adds no variational source term.  Exact automorphism and orbit-stabilizer arithmetic validates finite quotient bookkeeping for existing graph witnesses, but it is not a sourced nonproxy `K`/`L_total` spectral action, not a canonical geometry theorem, and not a full 16-node generator.\n")
    append_once(AGENTS, "Current orbit-stabilizer exact quotient guardrail (P2789/S1739, 2026-06-16)", "## Current orbit-stabilizer exact quotient guardrail (P2789/S1739, 2026-06-16)\n\n- P2789 verifies exact automorphism/orbit-stabilizer quotient arithmetic for the complete P2787 8-node quotient: all six stored labeled-member counts match `8! / |Aut(G)|` and sum to 19,355.  It also records exact stabilizer and `S_16` orbit sizes for the seven local P2786/P2785 16-node representatives.\n- This is a finite-group bookkeeping certificate for existing witnesses only; it is not the required full connected 16-node 4-regular generator/toolchain and does not source geometry from `K`/`L_total`.\n- Do not promote orbit-stabilizer consistency to canonical geometry, strict spectral source law, selector closure, kernel full-expression, role-bearing `L_total`, bridge closure, role transfer, or ToE closure.  A next admissible move must supply an actual certified full 16-node generator artifact/toolchain or export a strict spectral action/source law before testing.\n")
    return payload


if __name__ == "__main__":
    main()
