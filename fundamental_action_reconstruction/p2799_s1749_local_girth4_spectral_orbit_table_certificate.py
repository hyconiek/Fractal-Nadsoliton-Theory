#!/usr/bin/env python3
"""P2799/S1749: local girth>=4 spectral/orbit table certificate.

P2798 filtered the current eight P2791 representatives into a six-class local
intersection with the external girth>=4 subtarget.  P2799 makes that local
intersection more proof-grade by computing exact spectral and orbit-stabilizer
bookkeeping for precisely those six graph6-provenanced witnesses.

This still does not import the external 16,828 graph shortcode/list.  It is a
local table for the already-present girth>=4 witnesses and a guard against
mistaking that local table for subtarget or full-class closure.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2786_s1736_graph6_provenance_toolchain_gate import graph6_decode
from p2788_s1738_complement_duality_exact_spectral_certificate import charpoly_expr, complement_edges, degree_sequence, normalize_edges
from p2789_s1739_orbit_stabilizer_exact_quotient_certificate import automorphism_count

GEN = ROOT / "generated"
P2786 = GEN / "p2786_s1736_graph6_provenance_toolchain_gate.json"
P2790 = GEN / "p2790_s1740_eighth_16node_witness_no_exhaustion_certificate.json"
P2798 = GEN / "p2798_s1748_external_girth4_subtarget_local_girth_filter_certificate.json"
OUT = GEN / "p2799_s1749_local_girth4_spectral_orbit_table_certificate.json"
MD = GEN / "p2799_s1749_local_girth4_spectral_orbit_table_certificate.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

N = 16
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


def graph6_rows(p2786: dict[str, Any], p2790: dict[str, Any]) -> dict[str, dict[str, Any]]:
    rows: dict[str, dict[str, Any]] = {}
    for row in p2786.get("graph6_provenance_witness", {}).get("provenance_rows", []):
        rows[row["representative"]] = {"source": "P2786_local_graph6", "graph6": row["graph6"]}
    eighth = p2790.get("eighth_witness", {})
    if "graph6" in eighth:
        rows["p2790_eighth_witness"] = {"source": "P2790_explicit_witness_graph6", "graph6": eighth["graph6"]}
    return rows


def coeffs(edges: tuple[tuple[int, int], ...], kind: str) -> list[int]:
    return [int(value) for value in charpoly_expr(edges, N, kind).all_coeffs()]


def spectral_orbit_witness(p2786: dict[str, Any], p2790: dict[str, Any], p2798: dict[str, Any]) -> dict[str, Any]:
    rows_by_label = graph6_rows(p2786, p2790)
    girth4_labels = p2798.get("external_girth4_subtarget_witness", {}).get("local_girth_at_least_4_labels", [])
    table_rows = []
    for label in girth4_labels:
        graph6 = rows_by_label[label]["graph6"]
        edges = normalize_edges(graph6_decode(graph6))
        comp = complement_edges(edges, N)
        aut_size = automorphism_count(edges, N)
        comp_aut_size = automorphism_count(comp, N)
        table_rows.append({
            "label": label,
            "source": rows_by_label[label]["source"],
            "graph6": graph6,
            "graph6_sha256": hashlib.sha256(graph6.encode("ascii")).hexdigest(),
            "edge_count": len(edges),
            "degree_sequence": degree_sequence(edges, N),
            "complement_edge_count": len(comp),
            "complement_degree_sequence": degree_sequence(comp, N),
            "automorphism_group_size": aut_size,
            "complement_automorphism_group_size": comp_aut_size,
            "automorphism_matches_complement": aut_size == comp_aut_size,
            "orbit_size_by_orbit_stabilizer": math.factorial(N) // aut_size,
            "adjacency_charpoly_coefficients": coeffs(edges, "adjacency"),
            "laplacian_charpoly_coefficients": coeffs(edges, "laplacian"),
            "signless_laplacian_charpoly_coefficients": coeffs(edges, "signless_laplacian"),
            "complement_adjacency_charpoly_coefficients": coeffs(comp, "adjacency"),
            "complement_laplacian_charpoly_coefficients": coeffs(comp, "laplacian"),
            "complement_signless_laplacian_charpoly_coefficients": coeffs(comp, "signless_laplacian"),
        })
    return {
        "input_girth4_label_count": len(girth4_labels),
        "table_row_count": len(table_rows),
        "table_rows": table_rows,
        "all_rows_16node_4regular": all(row["edge_count"] == 32 and row["degree_sequence"] == [4] * N for row in table_rows),
        "all_complements_16node_11regular": all(row["complement_edge_count"] == 88 and row["complement_degree_sequence"] == [11] * N for row in table_rows),
        "all_automorphism_counts_match_complements": all(row["automorphism_matches_complement"] for row in table_rows),
        "distinct_adjacency_charpoly_count": len({tuple(row["adjacency_charpoly_coefficients"]) for row in table_rows}),
        "total_labeled_orbit_size_for_local_girth4_rows": sum(row["orbit_size_by_orbit_stabilizer"] for row in table_rows),
        "finite_certificate_statement": "The six current girth>=4 witnesses have exact spectral/orbit/complement bookkeeping, but this table covers only the local filtered witnesses and does not import the external 16,828-class shortcode list.",
    }


def acceptance_matrix(witness: dict[str, Any], p2786: dict[str, Any], p2790: dict[str, Any], p2798: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "p2786_graph6_gate_present": p2786.get("status") == "P2786_GRAPH6_PROVENANCE_TOOLCHAIN_GATE_NO_CLOSURE",
        "p2790_eighth_witness_present": p2790.get("status") == "P2790_EIGHTH_16NODE_WITNESS_NO_EXHAUSTION_CERTIFICATE_NO_CLOSURE",
        "p2798_girth4_filter_present": p2798.get("status") == "P2798_EXTERNAL_GIRTH4_SUBTARGET_LOCAL_GIRTH_FILTER_CERTIFICATE_NO_CLOSURE",
        "six_local_girth4_rows_tabulated": witness["table_row_count"] == 6,
        "all_rows_16node_4regular": witness["all_rows_16node_4regular"],
        "all_complements_16node_11regular": witness["all_complements_16node_11regular"],
        "all_automorphism_counts_match_complements": witness["all_automorphism_counts_match_complements"],
        "external_girth4_shortcode_graph_list_imported": False,
        "full_generator_artifact_imported": False,
        "strict_nadsoliton_spectral_source_law_exported": False,
        "kernel_or_ltotal_variational_coupling_exported": False,
    }
    return {
        "facts": facts,
        "accepted_as_local_girth4_spectral_orbit_table_certificate": all(facts[key] for key in [
            "p2786_graph6_gate_present",
            "p2790_eighth_witness_present",
            "p2798_girth4_filter_present",
            "six_local_girth4_rows_tabulated",
            "all_rows_16node_4regular",
            "all_complements_16node_11regular",
            "all_automorphism_counts_match_complements",
        ]),
        "accepted_as_full_16node_canonical_generator_certificate": False,
        "accepted_as_strict_spectral_source_law": False,
        "accepted_as_canonical_nadsoliton_geometry_source": False,
        "accepted_as_ltotal_or_toe_promotion": False,
        "missing_criteria": [key for key, value in facts.items() if not value],
        "blocker": "P2799 strengthens the local girth>=4 witness table, but it still covers only six current graph6-provenanced rows and imports no external shortcode graph list or strict spectral source law.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    w = payload["local_girth4_spectral_orbit_witness"]
    lines = [
        "# P2799/S1749 local girth>=4 spectral/orbit table certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Local girth>=4 exact table summary",
        f"- input_girth4_label_count={w['input_girth4_label_count']}",
        f"- table_row_count={w['table_row_count']}",
        f"- all_rows_16node_4regular={w['all_rows_16node_4regular']}",
        f"- all_complements_16node_11regular={w['all_complements_16node_11regular']}",
        f"- all_automorphism_counts_match_complements={w['all_automorphism_counts_match_complements']}",
        f"- distinct_adjacency_charpoly_count={w['distinct_adjacency_charpoly_count']}",
        f"- total_labeled_orbit_size_for_local_girth4_rows={w['total_labeled_orbit_size_for_local_girth4_rows']}",
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
    p2798 = read_json(P2798)
    witness = spectral_orbit_witness(p2786, p2790, p2798)
    acceptance = acceptance_matrix(witness, p2786, p2790, p2798)
    payload = {
        "status": "P2799_LOCAL_GIRTH4_SPECTRAL_ORBIT_TABLE_CERTIFICATE_NO_CLOSURE",
        "input_hashes": {"P2786": sha(P2786), "P2790": sha(P2790), "P2798": sha(P2798)},
        "input_statuses": {"P2786": p2786.get("status"), "P2790": p2790.get("status"), "P2798": p2798.get("status")},
        "audited_question": "Can the current local girth>=4 intersection be upgraded with exact spectral/orbit/complement bookkeeping without importing or claiming the external shortcode list?",
        "local_girth4_spectral_orbit_witness": witness,
        "acceptance_matrix": acceptance,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": acceptance["blocker"],
            "next_honest_step": "Use P2799 only as the exact local table for the six current girth>=4 witnesses.  The next proof-grade move remains importing the actual linked shortcode/graph-list artifact for the 16,828 girth>=4 graphs with hash/provenance and running exact quotient/charpoly/complement/orbit auditing, or importing the full 8,037,418-class toolchain/list, or exporting a strict spectral source law fixing the admissible class and K/L_total coupling.  Otherwise preserve the P2697-P2799 no-canonical-geometry/no-closure certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2799/S1749 local girth4 spectral-orbit table certificate", "## P2799/S1749 local girth>=4 spectral-orbit table certificate\n\n`P2799/S1749` computes exact spectral/orbit/complement bookkeeping for precisely the six current P2798 girth>=4 witnesses.  Each row is a 16-node 4-regular graph with an 11-regular complement; automorphism counts match the complement automorphism counts; adjacency, Laplacian, signless-Laplacian, and complement characteristic polynomials are recorded exactly.  This strengthens the local filtered witness table, but it still does not import the external 16,828-class shortcode graph list, does not provide full graph-list provenance, and does not export a strict spectral source law or `K`/`L_total` coupling.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2799/S1749 local girth4 spectral-orbit Ltotal guard", "## P2799/S1749 local girth>=4 spectral-orbit Ltotal guard\n\n`P2799/S1749` adds no variational source term.  Exact spectral/orbit bookkeeping for six local girth>=4 witnesses does not import the external graph list, source a spectral action, or license any `K`/`L_total` promotion.\n")
    append_once(AGENTS, "Current local girth4 spectral-orbit guardrail (P2799/S1749, 2026-06-16)", "## Current local girth>=4 spectral-orbit guardrail (P2799/S1749, 2026-06-16)\n\n- P2799 computes exact spectral/orbit/complement bookkeeping for the six current P2798 girth>=4 witnesses, including exact characteristic polynomials, automorphism counts, orbit sizes, and 11-regular complement checks.\n- This is a local filtered witness-table certificate only; it is not the linked 16,828-class shortcode graph-list import, not full `graph6` provenance, and not a strict spectral source law.\n- Do not promote the local girth>=4 spectral/orbit table to canonical geometry, strict spectral source law, selector closure, kernel full-expression, role-bearing `L_total`, bridge closure, role transfer, or ToE closure.\n")
    return payload


if __name__ == "__main__":
    main()
