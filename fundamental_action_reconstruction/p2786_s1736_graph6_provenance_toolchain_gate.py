#!/usr/bin/env python3
"""P2786/S1736: graph6 provenance and toolchain gate.

P2785 gave an exact local algebra certificate for the seven P2781/P2783/P2784
representatives.  The next recommended global move asks for a canonical
full-class generator with reproducible graph6/hash provenance.  This bounded
step does not pretend to provide that full generator.  Instead, it builds the
smallest honest provenance layer on the already-certified seven representatives:

* self-contained graph6 encoding/decoding for the seven 16-node simple graphs;
* round-trip edge-set equality checks;
* SHA-256 hashes of the graph6 strings and exact P2785 characteristic data;
* an explicit local toolchain gate for nauty-style commands and common Python
  graph packages.

The output is therefore a reproducibility/provenance certificate for the local
seven-class witness and a blocker certificate for full-class canonical
generation in the current environment.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import shutil
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2778_s1728_max_symmetry_16node_geometry_source_audit import N, candidate_edge_sets
from p2779_s1729_16node_circulant_full_spectrum_quotient_audit import quotient_classes
from p2781_s1731_enumerated_two_layer_c8_spectrum_collision_audit import all_two_layer_c8_candidates
from p2785_s1735_exact_characteristic_polynomial_certificate import exact_polynomial_witness

GEN = ROOT / "generated"
P2785 = GEN / "p2785_s1735_exact_characteristic_polynomial_certificate.json"
OUT = GEN / "p2786_s1736_graph6_provenance_toolchain_gate.json"
MD = GEN / "p2786_s1736_graph6_provenance_toolchain_gate.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

NEGATIVE_EXPORT_FLAGS = [
    "canonical_full_class_generator_certified",
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


def normalize_edges(edges: list[tuple[int, int]]) -> list[tuple[int, int]]:
    return sorted((min(a, b), max(a, b)) for a, b in edges)


def graph6_encode(edges: list[tuple[int, int]], n: int = N) -> str:
    if not 0 <= n <= 62:
        raise ValueError("P2786 graph6 encoder only covers the small-n graph6 header used here")
    edge_set = set(normalize_edges(edges))
    bits = []
    for j in range(1, n):
        for i in range(j):
            bits.append(1 if (i, j) in edge_set else 0)
    while len(bits) % 6:
        bits.append(0)
    chars = [chr(n + 63)]
    for index in range(0, len(bits), 6):
        value = 0
        for bit in bits[index:index + 6]:
            value = (value << 1) | bit
        chars.append(chr(value + 63))
    return "".join(chars)


def graph6_decode(payload: str) -> list[tuple[int, int]]:
    if not payload:
        raise ValueError("empty graph6 payload")
    n = ord(payload[0]) - 63
    if n != N:
        raise ValueError(f"expected graph6 order {N}, got {n}")
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


def toolchain_status() -> dict[str, Any]:
    commands = ["geng", "labelg", "showg", "countg"]
    modules = ["networkx", "igraph", "graph_tool"]
    return {
        "commands": {command: shutil.which(command) for command in commands},
        "python_modules": {module: importlib.util.find_spec(module) is not None for module in modules},
    }


def provenance_witness() -> dict[str, Any]:
    candidates = candidate_edge_sets() + all_two_layer_c8_candidates()
    by_name = {row["geometry"]: row for row in candidates}
    classes = quotient_classes(candidates)
    exact = exact_polynomial_witness()
    exact_by_name = {row["representative"]: row for row in exact["exact_rows"]}
    rows = []
    for rep in classes:
        name = rep["representative"]
        edges = normalize_edges(by_name[name]["edges"])
        graph6 = graph6_encode(edges)
        decoded = graph6_decode(graph6)
        exact_row = exact_by_name[name]
        exact_payload = json.dumps({
            "adjacency": exact_row["adjacency_charpoly_coefficients"],
            "laplacian": exact_row["laplacian_charpoly_coefficients"],
            "signless": exact_row["signless_laplacian_charpoly_coefficients"],
            "trees": exact_row["kirchhoff_spanning_tree_count_exact"],
        }, sort_keys=True)
        rows.append({
            "representative": name,
            "member_count": len(rep["members"]),
            "graph6": graph6,
            "graph6_sha256": hashlib.sha256(graph6.encode("ascii")).hexdigest(),
            "roundtrip_edges_equal": decoded == edges,
            "edge_count": len(edges),
            "exact_payload_sha256": hashlib.sha256(exact_payload.encode("utf-8")).hexdigest(),
        })
    toolchain = toolchain_status()
    command_available = {key: value is not None for key, value in toolchain["commands"].items()}
    module_available = toolchain["python_modules"]
    return {
        "source_class": "Same seven local quotient representatives; graph6 provenance is local and not a full-class generator.",
        "representative_count": len(rows),
        "provenance_rows": rows,
        "all_graph6_roundtrips_exact": all(row["roundtrip_edges_equal"] for row in rows),
        "all_rows_have_16node_4regular_edge_count": all(row["edge_count"] == 32 for row in rows),
        "distinct_graph6_hash_count": len({row["graph6_sha256"] for row in rows}),
        "toolchain_status": toolchain,
        "canonical_generation_tool_available": bool(command_available.get("geng") and command_available.get("labelg")),
        "common_python_graph_generator_available": bool(module_available.get("networkx") or module_available.get("igraph") or module_available.get("graph_tool")),
        "finite_certificate_statement": "The seven local representatives have exact graph6 round-trips and reproducible hashes, but the current environment does not certify a full canonical 16-node 4-regular generator.",
    }


def acceptance_matrix(witness: dict[str, Any], p2785: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "p2785_exact_certificate_present": p2785.get("status") == "P2785_EXACT_CHARACTERISTIC_POLYNOMIAL_CERTIFICATE_NO_CLOSURE",
        "seven_representatives_reused_without_class_expansion": witness["representative_count"] == 7,
        "all_graph6_roundtrips_exact": witness["all_graph6_roundtrips_exact"],
        "all_rows_have_16node_4regular_edge_count": witness["all_rows_have_16node_4regular_edge_count"],
        "distinct_graph6_hash_count_is_seven": witness["distinct_graph6_hash_count"] == 7,
        "canonical_generation_tool_available": witness["canonical_generation_tool_available"],
        "common_python_graph_generator_available": witness["common_python_graph_generator_available"],
        "strict_nadsoliton_spectral_source_law_exported": False,
        "kernel_or_ltotal_variational_coupling_exported": False,
    }
    return {
        "facts": facts,
        "accepted_as_local_graph6_provenance_certificate": all(facts[key] for key in [
            "p2785_exact_certificate_present",
            "seven_representatives_reused_without_class_expansion",
            "all_graph6_roundtrips_exact",
            "all_rows_have_16node_4regular_edge_count",
            "distinct_graph6_hash_count_is_seven",
        ]),
        "accepted_as_full_canonical_generation_tool_certificate": False,
        "accepted_as_strict_spectral_source_law": False,
        "accepted_as_canonical_nadsoliton_geometry_source": False,
        "accepted_as_ltotal_or_toe_promotion": False,
        "missing_criteria": [key for key, value in facts.items() if not value],
        "blocker": "Local graph6 provenance is reproducible, but current artifacts/environment still do not provide a certified full-class canonical generator or strict K/L_total spectral source law.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    w = payload["graph6_provenance_witness"]
    lines = [
        "# P2786/S1736 graph6 provenance and toolchain gate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Local graph6 provenance result",
        f"- representative_count={w['representative_count']}",
        f"- all_graph6_roundtrips_exact={w['all_graph6_roundtrips_exact']}",
        f"- all_rows_have_16node_4regular_edge_count={w['all_rows_have_16node_4regular_edge_count']}",
        f"- distinct_graph6_hash_count={w['distinct_graph6_hash_count']}",
        f"- canonical_generation_tool_available={w['canonical_generation_tool_available']}",
        f"- common_python_graph_generator_available={w['common_python_graph_generator_available']}",
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
    p2785 = read_json(P2785)
    witness = provenance_witness()
    acceptance = acceptance_matrix(witness, p2785)
    payload = {
        "status": "P2786_GRAPH6_PROVENANCE_TOOLCHAIN_GATE_NO_CLOSURE",
        "input_hashes": {"P2785": sha(P2785)},
        "input_statuses": {"P2785": p2785.get("status")},
        "audited_question": "Can the P2785 exact local certificate be upgraded to reproducible graph6 provenance or a full canonical-generation tool certificate?",
        "graph6_provenance_witness": witness,
        "acceptance_matrix": acceptance,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": acceptance["blocker"],
            "next_honest_step": "Use P2786 as local graph6/hash provenance only.  The next honest move is exactly one of: add an actual certified full-class generator artifact/toolchain (for example nauty/geng plus labelg output with reproducible graph6 hashes) and then run exact-polynomial collision auditing on that imported class; or export a strict nadsoliton spectral action/source law fixing the admissible class, target spectrum, and K/L_total coupling before testing.  Otherwise preserve the P2697-P2786 no-canonical-geometry/no-closure certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2786/S1736 graph6 provenance and toolchain gate", "## P2786/S1736 graph6 provenance and toolchain gate\n\n`P2786/S1736` builds a local reproducibility layer for the P2785 seven-class exact certificate: each representative is encoded to graph6, decoded back to the same edge set, and assigned graph6/exact-payload SHA-256 hashes.  The current environment does not expose a certified full-class canonical generator (`geng`/`labelg`) or common Python graph-generation package, so this remains local graph6 provenance rather than a full 16-node 4-regular generator certificate.  No strict spectral source law, canonical geometry, kernel full-expression theorem, role-bearing `L_total`, bridge closure, role transfer, selector closure, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2786/S1736 graph6 provenance Ltotal guard", "## P2786/S1736 graph6 provenance Ltotal guard\n\n`P2786/S1736` adds no variational source term.  graph6/hash provenance improves reproducibility of the local seven-class exact certificate, but it is not a sourced nonproxy `K`/`L_total` spectral action and not a full canonical graph-class generator.\n")
    append_once(AGENTS, "Current graph6 provenance and toolchain gate guardrail (P2786/S1736, 2026-06-16)", "## Current graph6 provenance and toolchain gate guardrail (P2786/S1736, 2026-06-16)\n\n- P2786 encodes the seven P2785 local representatives to graph6, verifies exact graph6 decode round-trips, and records graph6/exact-payload hashes.\n- The current environment does not provide a certified full-class canonical generator/toolchain such as `geng` plus `labelg`; local graph6 provenance is not a full graph-class certificate.\n- Do not promote this provenance gate to canonical geometry, strict spectral source law, selector closure, kernel full-expression, role-bearing `L_total`, bridge closure, role transfer, or ToE closure.  A next admissible move must supply an actual certified full-class generator artifact/toolchain or export a strict spectral action/source law before testing.\n")
    return payload


if __name__ == "__main__":
    main()
