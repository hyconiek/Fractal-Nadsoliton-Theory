#!/usr/bin/env python3
"""P2798/S1748: external girth>=4 subtarget and local girth-filter certificate.

P2797 sized the full connected 16-node 4-regular gap using an external target
count.  P2798 takes a more computational but still honest subtarget step: use the
linked Meringer girth>=4 count (16,828) only as a target count, then compute the
exact girth of every current P2791 representative available through P2786/P2790
graph6 provenance.

The result is not a graph-list import.  It is a local exact girth filter plus an
external subtarget gap certificate: current witnesses cover only 6 of the 16,828
reported connected 16-node 4-regular girth>=4 classes.
"""
from __future__ import annotations

import hashlib
import json
from collections import deque
from fractions import Fraction
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2786_s1736_graph6_provenance_toolchain_gate import graph6_decode
from p2788_s1738_complement_duality_exact_spectral_certificate import normalize_edges

GEN = ROOT / "generated"
P2786 = GEN / "p2786_s1736_graph6_provenance_toolchain_gate.json"
P2790 = GEN / "p2790_s1740_eighth_16node_witness_no_exhaustion_certificate.json"
P2791 = GEN / "p2791_s1741_eight_class_orbit_lower_bound_certificate.json"
P2797 = GEN / "p2797_s1747_external_full_count_gap_certificate.json"
OUT = GEN / "p2798_s1748_external_girth4_subtarget_local_girth_filter_certificate.json"
MD = GEN / "p2798_s1748_external_girth4_subtarget_local_girth_filter_certificate.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

N = 16
EXTERNAL_GIRTH4_COUNT = 16_828
EXTERNAL_SOURCE = {
    "name": "Markus Meringer Regular Graphs page",
    "url": "https://www.mathe2.uni-bayreuth.de/markus/reggraphs.html",
    "linked_detail_url": "https://www.mathe2.uni-bayreuth.de/markus/REGGRAPHS/16_4_4.html",
    "retrieved_utc_date": "2026-06-16",
    "table": "Connected regular graphs with girth at least 4",
    "row": "vertices=16",
    "column": "degree=4",
    "reported_count": EXTERNAL_GIRTH4_COUNT,
    "detail_context": "The count is a link to a detail page offering a shortcode file for the chosen graphs; P2798 records only the target count and does not import the shortcode file.",
}
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


def graph_girth(edges: tuple[tuple[int, int], ...], n: int = N) -> int | None:
    adj = [set() for _ in range(n)]
    for left, right in edges:
        adj[left].add(right)
        adj[right].add(left)
    best = n + 1
    for source in range(n):
        distance = [-1] * n
        parent = [-1] * n
        distance[source] = 0
        queue: deque[int] = deque([source])
        while queue:
            node = queue.popleft()
            for nxt in adj[node]:
                if distance[nxt] == -1:
                    distance[nxt] = distance[node] + 1
                    parent[nxt] = node
                    queue.append(nxt)
                elif parent[node] != nxt and parent[nxt] != node:
                    best = min(best, distance[node] + distance[nxt] + 1)
    return None if best == n + 1 else best


def current_representative_edges(p2786: dict[str, Any], p2790: dict[str, Any]) -> list[dict[str, Any]]:
    rows = []
    for row in p2786.get("graph6_provenance_witness", {}).get("provenance_rows", []):
        rows.append({
            "label": row["representative"],
            "source": "P2786_local_graph6",
            "graph6": row["graph6"],
            "edges": normalize_edges(graph6_decode(row["graph6"])),
        })
    eighth = p2790.get("eighth_witness", {})
    if "graph6" in eighth:
        rows.append({
            "label": "p2790_eighth_witness",
            "source": "P2790_explicit_witness_graph6",
            "graph6": eighth["graph6"],
            "edges": normalize_edges(graph6_decode(eighth["graph6"])),
        })
    return rows


def girth4_witness(p2786: dict[str, Any], p2790: dict[str, Any], p2791: dict[str, Any], p2797: dict[str, Any]) -> dict[str, Any]:
    representatives = current_representative_edges(p2786, p2790)
    rows = []
    for row in representatives:
        girth = graph_girth(row["edges"])
        rows.append({
            "label": row["label"],
            "source": row["source"],
            "edge_count": len(row["edges"]),
            "girth": girth,
            "girth_at_least_4": girth is not None and girth >= 4,
            "graph6_sha256": hashlib.sha256(row["graph6"].encode("ascii")).hexdigest(),
        })
    girth4_labels = [row["label"] for row in rows if row["girth_at_least_4"]]
    triangle_labels = [row["label"] for row in rows if row["girth"] == 3]
    fraction = Fraction(len(girth4_labels), EXTERNAL_GIRTH4_COUNT)
    return {
        "external_source": EXTERNAL_SOURCE,
        "external_connected_16node_4regular_girth4_class_count": EXTERNAL_GIRTH4_COUNT,
        "p2791_representative_count": p2791.get("orbit_lower_bound_witness", {}).get("representative_count"),
        "p2797_full_count_gap_status": p2797.get("status"),
        "local_representative_girth_rows": rows,
        "local_representative_count": len(rows),
        "local_girth_at_least_4_count": len(girth4_labels),
        "local_girth_at_least_4_labels": girth4_labels,
        "local_triangle_containing_count": len(triangle_labels),
        "local_triangle_containing_labels": triangle_labels,
        "external_girth4_gap_after_current_witnesses": EXTERNAL_GIRTH4_COUNT - len(girth4_labels),
        "local_girth4_coverage_fraction_exact": f"{fraction.numerator}/{fraction.denominator}",
        "local_girth4_coverage_fraction_decimal": f"{len(girth4_labels) / EXTERNAL_GIRTH4_COUNT:.12f}",
        "finite_certificate_statement": "Exact local girth filtering leaves 6 current representatives in the girth>=4 subtarget, compared with the external target count 16,828; 16,822 girth>=4 classes remain outside current local witnesses.",
    }


def acceptance_matrix(witness: dict[str, Any], p2786: dict[str, Any], p2790: dict[str, Any], p2791: dict[str, Any], p2797: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "p2786_graph6_gate_present": p2786.get("status") == "P2786_GRAPH6_PROVENANCE_TOOLCHAIN_GATE_NO_CLOSURE",
        "p2790_eighth_witness_present": p2790.get("status") == "P2790_EIGHTH_16NODE_WITNESS_NO_EXHAUSTION_CERTIFICATE_NO_CLOSURE",
        "p2791_lower_bound_present": p2791.get("status") == "P2791_EIGHT_CLASS_ORBIT_LOWER_BOUND_CERTIFICATE_NO_CLOSURE",
        "p2797_full_gap_present": p2797.get("status") == "P2797_EXTERNAL_FULL_COUNT_GAP_CERTIFICATE_NO_CLOSURE",
        "external_girth4_count_recorded": witness["external_connected_16node_4regular_girth4_class_count"] == EXTERNAL_GIRTH4_COUNT,
        "all_eight_current_representatives_girth_filtered": witness["local_representative_count"] == 8,
        "six_current_representatives_in_girth4_subtarget": witness["local_girth_at_least_4_count"] == 6,
        "girth4_gap_is_positive": witness["external_girth4_gap_after_current_witnesses"] == 16_822,
        "full_girth4_shortcode_graph_list_imported": False,
        "strict_nadsoliton_spectral_source_law_exported": False,
        "kernel_or_ltotal_variational_coupling_exported": False,
    }
    return {
        "facts": facts,
        "accepted_as_external_girth4_subtarget_local_girth_filter_certificate": all(facts[key] for key in [
            "p2786_graph6_gate_present",
            "p2790_eighth_witness_present",
            "p2791_lower_bound_present",
            "p2797_full_gap_present",
            "external_girth4_count_recorded",
            "all_eight_current_representatives_girth_filtered",
            "six_current_representatives_in_girth4_subtarget",
            "girth4_gap_is_positive",
        ]),
        "accepted_as_full_16node_canonical_generator_certificate": False,
        "accepted_as_strict_spectral_source_law": False,
        "accepted_as_canonical_nadsoliton_geometry_source": False,
        "accepted_as_ltotal_or_toe_promotion": False,
        "missing_criteria": [key for key, value in facts.items() if not value],
        "blocker": "P2798 computes exact local girths and sizes a reachable external girth>=4 subtarget gap, but it does not import the shortcode graph list and therefore remains a target/gap certificate rather than generator provenance or a strict source law.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    w = payload["external_girth4_subtarget_witness"]
    lines = [
        "# P2798/S1748 external girth>=4 subtarget local girth-filter certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## External girth>=4 subtarget and local filter",
        f"- external_connected_16node_4regular_girth4_class_count={w['external_connected_16node_4regular_girth4_class_count']}",
        f"- local_representative_count={w['local_representative_count']}",
        f"- local_girth_at_least_4_count={w['local_girth_at_least_4_count']}",
        f"- local_girth_at_least_4_labels={w['local_girth_at_least_4_labels']}",
        f"- local_triangle_containing_labels={w['local_triangle_containing_labels']}",
        f"- external_girth4_gap_after_current_witnesses={w['external_girth4_gap_after_current_witnesses']}",
        f"- local_girth4_coverage_fraction_exact={w['local_girth4_coverage_fraction_exact']}",
        "",
        "## External source",
        f"- source={w['external_source']['name']}",
        f"- url={w['external_source']['url']}",
        f"- linked_detail_url={w['external_source']['linked_detail_url']}",
        f"- retrieved_utc_date={w['external_source']['retrieved_utc_date']}",
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
    p2790 = read_json(P2790)
    p2791 = read_json(P2791)
    p2797 = read_json(P2797)
    witness = girth4_witness(p2786, p2790, p2791, p2797)
    acceptance = acceptance_matrix(witness, p2786, p2790, p2791, p2797)
    payload = {
        "status": "P2798_EXTERNAL_GIRTH4_SUBTARGET_LOCAL_GIRTH_FILTER_CERTIFICATE_NO_CLOSURE",
        "input_hashes": {"P2786": sha(P2786), "P2790": sha(P2790), "P2791": sha(P2791), "P2797": sha(P2797)},
        "input_statuses": {"P2786": p2786.get("status"), "P2790": p2790.get("status"), "P2791": p2791.get("status"), "P2797": p2797.get("status")},
        "audited_question": "How much of the externally linked girth>=4 16-node 4-regular subtarget is covered by exact local girth filtering of the current P2791 representatives?",
        "external_girth4_subtarget_witness": witness,
        "acceptance_matrix": acceptance,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": acceptance["blocker"],
            "next_honest_step": "Use P2798 only as a local exact-girth filter plus external girth>=4 target-gap certificate.  The next proof-grade move is to import the actual linked shortcode/graph-list artifact for the 16,828 girth>=4 graphs with hash/provenance and run exact quotient/charpoly/complement/orbit auditing, or else import the full 8,037,418-class toolchain/list, or export a strict spectral source law fixing the admissible class and K/L_total coupling.  Otherwise preserve the P2697-P2798 no-canonical-geometry/no-closure certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2798/S1748 external girth4 subtarget local filter certificate", "## P2798/S1748 external girth>=4 subtarget local filter certificate\n\n`P2798/S1748` uses the linked Meringer girth>=4 subtarget only as a count target: the connected regular graph with girth at least 4 table reports `16,828` connected 4-regular graphs on 16 vertices and links to a shortcode detail page.  The script computes exact girth for all eight current P2791 representatives available through P2786/P2790 graph6 provenance: six have girth at least 4 and two contain triangles, leaving a subtarget gap of `16,822`.  This is a local girth-filter and target-gap certificate only; it does not import the shortcode graph list, does not provide full graph-list provenance, and does not export a strict spectral source law or `K`/`L_total` coupling.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2798/S1748 girth4 subtarget Ltotal guard", "## P2798/S1748 girth>=4 subtarget Ltotal guard\n\n`P2798/S1748` adds no variational source term.  Exact local girth filtering and an external girth>=4 target count do not import the graph list, source a spectral action, or license any `K`/`L_total` promotion.\n")
    append_once(AGENTS, "Current external girth4 subtarget guardrail (P2798/S1748, 2026-06-16)", "## Current external girth>=4 subtarget guardrail (P2798/S1748, 2026-06-16)\n\n- P2798 records the external linked girth>=4 target count `16,828` and computes exact girth for all eight current P2791 representatives: six satisfy girth>=4 and two contain triangles, leaving a girth>=4 subtarget gap of `16,822`.\n- This is a local girth-filter/target-gap certificate only; it is not the linked shortcode graph-list import, not full `graph6` provenance, and not a strict spectral source law.\n- Do not promote the girth>=4 target count or local filter to canonical geometry, strict spectral source law, selector closure, kernel full-expression, role-bearing `L_total`, bridge closure, role transfer, or ToE closure.  The next admissible move is importing the actual shortcode/list provenance or exporting a strict spectral source law.\n")
    return payload


if __name__ == "__main__":
    main()
