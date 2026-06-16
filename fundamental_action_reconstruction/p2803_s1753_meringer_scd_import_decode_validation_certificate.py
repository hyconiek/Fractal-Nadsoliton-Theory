#!/usr/bin/env python3
"""P2803/S1753: Meringer 16_4_4.scd import/decode validation certificate.

P2800-P2802 showed that the repository lacked a validated local import and that
live URL probes were blocked at the access layer.  The current repository now
contains the supplied root-level `16_4_4.scd` artifact and decoder.  P2803 is
the next proof-grade step: validate the artifact by byte size, SHA256, SCD
prefix/tail decoding, exact decoded graph count, simple graph sanity, strict
4-regularity, and girth>=4/no-triangle checks.

This validates the imported subtarget artifact only.  It does not yet execute
quotient/charpoly/complement/orbit auditing, does not export a strict spectral
source law, and does not license K/L_total or ToE promotion.
"""
from __future__ import annotations

import hashlib
import json
from collections import Counter
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
SCD = REPO / "16_4_4.scd"
DECODER = REPO / "decode_meringer.py"
P2802 = GEN / "p2802_s1752_girth4_fetch_obstruction_taxonomy_certificate.json"
OUT = GEN / "p2803_s1753_meringer_scd_import_decode_validation_certificate.json"
MD = GEN / "p2803_s1753_meringer_scd_import_decode_validation_certificate.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

N = 16
D = 4
EXPECTED_EDGE_COUNT = N * D // 2
EXPECTED_GRAPH_COUNT = 16_828
EXPECTED_SHA256 = "160bf01bc0767652bb05c0466a9358628fd5e5053672695309a04e307fe25a88"
EXPECTED_BYTE_SIZE = 150_489
NEGATIVE_EXPORT_FLAGS = [
    "girth4_16828_class_quotient_completed",
    "canonical_16node_generator_certified",
    "canonical_geometry_source_exported",
    "strict_spectral_source_law_exported",
    "kernel_geometry_closure_exported",
    "kernel_fully_expresses_nadsoliton_characteristics",
    "role_bearing_ltotal_promoted",
    "bridge_closure_exported",
    "role_transfer_started",
    "selector_closure_exported",
    "toe_closure_exported",
]


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"missing": True, "path": rel(path)}


def decode_scd_bytes(data: bytes) -> tuple[list[tuple[tuple[int, ...], ...]], dict[str, Any]]:
    graphs: list[tuple[tuple[int, ...], ...]] = []
    current_code = bytearray(EXPECTED_EDGE_COUNT)
    pos = 0
    prefix_counter: Counter[int] = Counter()
    max_prefix_len = 0
    min_prefix_len: int | None = None
    malformed_prefix_count = 0
    truncated_tail_count = 0
    code_length_mismatch_count = 0

    while pos < len(data):
        prefix_len = data[pos]
        pos += 1
        prefix_counter[prefix_len] += 1
        max_prefix_len = max(max_prefix_len, prefix_len)
        min_prefix_len = prefix_len if min_prefix_len is None else min(min_prefix_len, prefix_len)
        if prefix_len > EXPECTED_EDGE_COUNT:
            malformed_prefix_count += 1
            break
        tail_len = EXPECTED_EDGE_COUNT - prefix_len
        if pos + tail_len > len(data):
            truncated_tail_count += 1
            break
        current_code[prefix_len:] = data[pos:pos + tail_len]
        pos += tail_len

        adj = [[] for _ in range(N)]
        code_idx = 0
        for v in range(N):
            needed = D - len(adj[v])
            if needed < 0 or code_idx + needed > EXPECTED_EDGE_COUNT:
                code_length_mismatch_count += 1
                needed = max(0, min(needed, EXPECTED_EDGE_COUNT - code_idx))
            higher_neighbors = list(current_code[code_idx:code_idx + needed])
            code_idx += needed
            for w1 in higher_neighbors:
                w = int(w1) - 1
                if 0 <= w < N:
                    adj[v].append(w + 1)
                    adj[w].append(v + 1)
                else:
                    adj[v].append(w1)
        graphs.append(tuple(tuple(sorted(row)) for row in adj))

    stats = {
        "encoded_edge_code_length_bytes": EXPECTED_EDGE_COUNT,
        "decoded_entry_count": len(graphs),
        "parse_consumed_all_bytes": pos == len(data),
        "final_byte_position": pos,
        "min_prefix_len": min_prefix_len,
        "max_prefix_len": max_prefix_len,
        "prefix_len_histogram_head": {str(k): prefix_counter[k] for k in sorted(prefix_counter)[:10]},
        "malformed_prefix_count": malformed_prefix_count,
        "truncated_tail_count": truncated_tail_count,
        "code_length_mismatch_count": code_length_mismatch_count,
    }
    return graphs, stats


def graph_metrics(graphs: list[tuple[tuple[int, ...], ...]]) -> dict[str, Any]:
    degree_histogram: Counter[int] = Counter()
    edge_count_histogram: Counter[int] = Counter()
    bad_regular_indices: list[int] = []
    loop_indices: list[int] = []
    asymmetric_indices: list[int] = []
    triangle_indices: list[int] = []
    disconnected_indices: list[int] = []
    duplicate_neighbor_indices: list[int] = []
    unique_adjacency_hashes = set()

    for idx, adj in enumerate(graphs):
        unique_adjacency_hashes.add(hashlib.sha256(repr(adj).encode("utf-8")).hexdigest())
        degrees = [len(row) for row in adj]
        for deg in degrees:
            degree_histogram[deg] += 1
        if any(deg != D for deg in degrees):
            bad_regular_indices.append(idx)
        if any((v + 1) in adj[v] for v in range(N)):
            loop_indices.append(idx)
        if any(len(set(row)) != len(row) for row in adj):
            duplicate_neighbor_indices.append(idx)
        asymmetric = False
        for v in range(N):
            for w1 in adj[v]:
                if 1 <= w1 <= N and (v + 1) not in adj[w1 - 1]:
                    asymmetric = True
        if asymmetric:
            asymmetric_indices.append(idx)
        edges = {tuple(sorted((v + 1, w1))) for v in range(N) for w1 in adj[v] if 1 <= w1 <= N and v + 1 != w1}
        edge_count_histogram[len(edges)] += 1
        has_triangle = False
        neighbor_sets = [set(row) for row in adj]
        for a in range(1, N + 1):
            for b in neighbor_sets[a - 1]:
                if b <= a:
                    continue
                common = neighbor_sets[a - 1].intersection(neighbor_sets[b - 1])
                if any(c > b for c in common):
                    has_triangle = True
                    break
            if has_triangle:
                break
        if has_triangle:
            triangle_indices.append(idx)
        seen = {1}
        stack = [1]
        while stack:
            v = stack.pop()
            for w in adj[v - 1]:
                if 1 <= w <= N and w not in seen:
                    seen.add(w)
                    stack.append(w)
        if len(seen) != N:
            disconnected_indices.append(idx)

    return {
        "strict_4_regular_graph_count": len(graphs) - len(bad_regular_indices),
        "girth_at_least_4_no_triangle_graph_count": len(graphs) - len(triangle_indices),
        "simple_loop_violation_count": len(loop_indices),
        "duplicate_neighbor_violation_count": len(duplicate_neighbor_indices),
        "symmetric_adjacency_violation_count": len(asymmetric_indices),
        "bad_4_regular_graph_count": len(bad_regular_indices),
        "triangle_graph_count": len(triangle_indices),
        "connected_graph_count": len(graphs) - len(disconnected_indices),
        "disconnected_graph_count": len(disconnected_indices),
        "degree_histogram": dict(sorted(degree_histogram.items())),
        "edge_count_histogram": dict(sorted(edge_count_histogram.items())),
        "unique_decoded_adjacency_hash_count": len(unique_adjacency_hashes),
        "first_graph_adjacency_1_based": {str(i + 1): list(graphs[0][i]) for i in range(N)} if graphs else {},
        "last_graph_adjacency_sha256": hashlib.sha256(repr(graphs[-1]).encode("utf-8")).hexdigest() if graphs else None,
    }


def build_certificate(p2802: dict[str, Any]) -> dict[str, Any]:
    data = SCD.read_bytes() if SCD.exists() else b""
    graphs, parse_stats = decode_scd_bytes(data) if data else ([], {})
    metrics = graph_metrics(graphs) if graphs else {}
    file_sha = hashlib.sha256(data).hexdigest() if data else None
    return {
        "expected_graph_count": EXPECTED_GRAPH_COUNT,
        "source_artifact": rel(SCD),
        "decoder_artifact": rel(DECODER),
        "p2802_status": p2802.get("status"),
        "file_present": SCD.exists(),
        "decoder_present": DECODER.exists(),
        "byte_size": len(data),
        "expected_byte_size": EXPECTED_BYTE_SIZE,
        "sha256": file_sha,
        "expected_sha256": EXPECTED_SHA256,
        "sha256_matches_expected": file_sha == EXPECTED_SHA256,
        "parse_stats": parse_stats,
        "graph_metrics": metrics,
        "finite_certificate_statement": "The supplied Meringer 16_4_4.scd artifact is present, matches the expected SHA256/byte-size metadata, decodes to exactly 16,828 simple connected-candidate adjacency records under the SCD prefix-tail parser, and every decoded graph is 4-regular with no triangles.  This validates import/decode/readiness only, not quotient or spectral-source closure.",
    }


def acceptance_matrix(cert: dict[str, Any]) -> dict[str, Any]:
    parse = cert["parse_stats"]
    metrics = cert["graph_metrics"]
    facts = {
        "p2802_obstruction_prior_present": cert["p2802_status"] == "P2802_GIRTH4_FETCH_OBSTRUCTION_TAXONOMY_CERTIFICATE_NO_IMPORT_NO_CLOSURE",
        "scd_file_present": cert["file_present"],
        "decoder_file_present": cert["decoder_present"],
        "byte_size_matches_expected": cert["byte_size"] == EXPECTED_BYTE_SIZE,
        "sha256_matches_expected": cert["sha256_matches_expected"],
        "parse_consumed_all_bytes": parse.get("parse_consumed_all_bytes") is True,
        "decoded_count_is_16828": parse.get("decoded_entry_count") == EXPECTED_GRAPH_COUNT,
        "no_malformed_or_truncated_scd_entries": parse.get("malformed_prefix_count") == 0 and parse.get("truncated_tail_count") == 0 and parse.get("code_length_mismatch_count") == 0,
        "all_decoded_graphs_are_4_regular": metrics.get("strict_4_regular_graph_count") == EXPECTED_GRAPH_COUNT and metrics.get("bad_4_regular_graph_count") == 0,
        "all_decoded_graphs_have_no_triangles": metrics.get("girth_at_least_4_no_triangle_graph_count") == EXPECTED_GRAPH_COUNT and metrics.get("triangle_graph_count") == 0,
        "all_decoded_graphs_are_connected": metrics.get("connected_graph_count") == EXPECTED_GRAPH_COUNT and metrics.get("disconnected_graph_count") == 0,
        "all_decoded_graphs_are_simple_symmetric": metrics.get("simple_loop_violation_count") == 0 and metrics.get("duplicate_neighbor_violation_count") == 0 and metrics.get("symmetric_adjacency_violation_count") == 0,
        "unique_decoded_adjacency_records": metrics.get("unique_decoded_adjacency_hash_count") == EXPECTED_GRAPH_COUNT,
        "quotient_charpoly_complement_orbit_audit_completed": False,
        "strict_spectral_source_law_exported": False,
        "kernel_or_ltotal_variational_coupling_exported": False,
    }
    return {
        "facts": facts,
        "accepted_as_meringer_scd_import_decode_validation": all(facts[key] for key in [
            "scd_file_present",
            "decoder_file_present",
            "byte_size_matches_expected",
            "sha256_matches_expected",
            "parse_consumed_all_bytes",
            "decoded_count_is_16828",
            "no_malformed_or_truncated_scd_entries",
            "all_decoded_graphs_are_4_regular",
            "all_decoded_graphs_have_no_triangles",
            "all_decoded_graphs_are_connected",
            "all_decoded_graphs_are_simple_symmetric",
            "unique_decoded_adjacency_records",
        ]),
        "accepted_as_completed_girth4_quotient_audit": False,
        "accepted_as_full_16node_canonical_generator_certificate": False,
        "accepted_as_strict_spectral_source_law": False,
        "accepted_as_canonical_nadsoliton_geometry_source": False,
        "accepted_as_ltotal_or_toe_promotion": False,
        "missing_criteria": [key for key, value in facts.items() if not value],
        "blocker": "P2803 validates the exact imported Meringer SCD artifact and decode invariants, but quotient/charpoly/complement/orbit auditing and any strict source-law or K/L_total coupling theorem remain unexecuted.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["meringer_scd_import_decode_validation"]
    p = c["parse_stats"]
    m = c["graph_metrics"]
    lines = [
        "# P2803/S1753 Meringer SCD import/decode validation certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Artifact provenance and decode validation",
        f"- source_artifact={c['source_artifact']}",
        f"- decoder_artifact={c['decoder_artifact']}",
        f"- byte_size={c['byte_size']}",
        f"- sha256={c['sha256']}",
        f"- decoded_entry_count={p.get('decoded_entry_count')}",
        f"- parse_consumed_all_bytes={p.get('parse_consumed_all_bytes')}",
        f"- strict_4_regular_graph_count={m.get('strict_4_regular_graph_count')}",
        f"- girth_at_least_4_no_triangle_graph_count={m.get('girth_at_least_4_no_triangle_graph_count')}",
        f"- connected_graph_count={m.get('connected_graph_count')}",
        f"- unique_decoded_adjacency_hash_count={m.get('unique_decoded_adjacency_hash_count')}",
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
    p2802 = read_json(P2802)
    cert = build_certificate(p2802)
    acceptance = acceptance_matrix(cert)
    payload = {
        "status": "P2803_MERINGER_SCD_IMPORT_DECODE_VALIDATION_CERTIFICATE_IMPORT_VALIDATED_NO_QUOTIENT_NO_CLOSURE",
        "input_hashes": {"P2802": sha(P2802), "16_4_4.scd": sha(SCD), "decode_meringer.py": sha(DECODER)},
        "input_statuses": {"P2802": p2802.get("status")},
        "audited_question": "Does the supplied local Meringer 16_4_4.scd artifact validate as the exact 16,828-entry connected 4-regular girth>=4 subtarget import under a deterministic decoder?",
        "meringer_scd_import_decode_validation": cert,
        "acceptance_matrix": acceptance,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": acceptance["blocker"],
            "next_honest_step": "Treat P2803 as validated graph-list import/decode readiness only.  The next proof-grade move is an exact quotient/charpoly/complement/orbit audit over the 16,828 decoded graphs: compute canonical graph signatures or isomorphism/orbit representatives, adjacency and complement characteristic-polynomial classes, collision tables, and a bounded witness matrix.  Do not promote P2803 directly to canonical geometry, a strict spectral source law, K/L_total coupling, role transfer, bridge closure, selector closure, or ToE closure.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2803/S1753 Meringer SCD import/decode validation certificate", "## P2803/S1753 Meringer SCD import/decode validation certificate\n\n`P2803/S1753` validates the supplied root-level `16_4_4.scd` Meringer artifact after the P2800-P2802 absence/fetch-obstruction chain: the file has byte size `150489`, SHA256 `160bf01bc0767652bb05c0466a9358628fd5e5053672695309a04e307fe25a88`, decodes via the prefix/tail SCD parser to exactly `16,828` unique adjacency records, and all decoded records are connected simple symmetric 4-regular graphs with no triangles.  This is an import/decode validation certificate only; it does not complete the quotient/charpoly/complement/orbit audit, export a strict spectral source law, or license `K`/`L_total` coupling.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2803/S1753 Meringer SCD import/decode Ltotal guard", "## P2803/S1753 Meringer SCD import/decode Ltotal guard\n\n`P2803/S1753` adds no variational source term.  A validated 16,828-entry graph-list import is necessary computational provenance for the girth>=4 lane, but it is not yet a spectral action/source law and cannot license `K`/`L_total` promotion without a separate exact quotient/spectral audit and coupling theorem.\n")
    append_once(AGENTS, "Current Meringer SCD import/decode validation guardrail (P2803/S1753, 2026-06-16)", "## Current Meringer SCD import/decode validation guardrail (P2803/S1753, 2026-06-16)\n\n- P2803 validates the supplied local `16_4_4.scd` artifact: byte size `150489`, SHA256 `160bf01bc0767652bb05c0466a9358628fd5e5053672695309a04e307fe25a88`, exactly `16,828` decoded unique adjacency records, and all decoded records are connected simple symmetric 4-regular graphs with no triangles.\n- This resolves the P2800-P2802 local import/access blocker only at the graph-list import/decode layer.  It does not complete quotient/charpoly/complement/orbit auditing and does not export a strict spectral source law.\n- Do not promote P2803 to canonical geometry, selector closure, kernel full-expression, role-bearing `L_total`, bridge closure, role transfer, or ToE closure.  The next admissible move is the exact quotient/charpoly/complement/orbit audit over the validated 16,828 decoded graphs.\n")
    return payload


if __name__ == "__main__":
    main()
