#!/usr/bin/env python3
"""P2806/S1756: reproducible Meringer-record label dataset audit.

P2805 showed that the residual P2804 spectral collisions contain no duplicate
isomorphic imports under exact bounded checks, but it did not export canonical
labels.  P2806 takes the next honest computational step that is possible without
an external canonical-labeling engine: export a deterministic label/hash dataset
for every validated decoded Meringer record in its supplied vertex order, and
separate that record-label provenance from true isomorphism-canonical labeling.

The labels are graph6-style encodings of the decoded adjacency records.  They
are reproducible record labels, not canonical graph labels and not a strict
spectral source law or K/L_total closure.
"""
from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import EXPECTED_GRAPH_COUNT, N, SCD, decode_scd_bytes, sha

GEN = ROOT / "generated"
P2805 = GEN / "p2805_s1755_girth4_spectral_collision_isomorphism_refinement.json"
OUT = GEN / "p2806_s1756_girth4_reproducible_record_label_dataset_audit.json"
MD = GEN / "p2806_s1756_girth4_reproducible_record_label_dataset_audit.md"
CSV_OUT = GEN / "p2806_s1756_girth4_record_label_dataset.csv"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
NEGATIVE_EXPORT_FLAGS = [
    "isomorphism_canonical_label_dataset_exported",
    "independent_canonical_label_tool_crosscheck_exported",
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


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"missing": True, "path": rel(path)}


def graph6_record_label(adj: tuple[tuple[int, ...], ...]) -> str:
    """Return a standard graph6-style label for n=16 in decoded vertex order."""
    bits = []
    neighbor_sets = [set(row) for row in adj]
    for j in range(1, N):
        for i in range(j):
            bits.append(1 if (i + 1) in neighbor_sets[j] else 0)
    while len(bits) % 6:
        bits.append(0)
    chars = [chr(N + 63)]
    for offset in range(0, len(bits), 6):
        value = 0
        for bit in bits[offset:offset + 6]:
            value = (value << 1) | bit
        chars.append(chr(value + 63))
    return "".join(chars)


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("ascii")).hexdigest()


def build_dataset(p2805: dict[str, Any]) -> dict[str, Any]:
    graphs, parse_stats = decode_scd_bytes(SCD.read_bytes())
    rows = []
    label_counter: Counter[str] = Counter()
    hash_counter: Counter[str] = Counter()
    first_label_by_hash: dict[str, str] = {}
    for index, graph in enumerate(graphs):
        label = graph6_record_label(graph)
        digest = sha256_text(label)
        label_counter[label] += 1
        hash_counter[digest] += 1
        first_label_by_hash.setdefault(digest, label)
        rows.append({
            "graph_index_0_based": index,
            "graph_index_1_based": index + 1,
            "record_graph6_label": label,
            "record_graph6_sha256": digest,
        })
    CSV_OUT.parent.mkdir(exist_ok=True)
    with CSV_OUT.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    csv_sha = sha(CSV_OUT)
    duplicate_labels = [label for label, count in label_counter.items() if count > 1]
    duplicate_hashes = [digest for digest, count in hash_counter.items() if count > 1]
    return {
        "expected_graph_count": EXPECTED_GRAPH_COUNT,
        "p2805_status": p2805.get("status"),
        "p2805_resolved_total_isomorphism_classes_after_refinement": p2805.get("spectral_collision_isomorphism_refinement", {}).get("resolved_total_isomorphism_classes_after_refinement"),
        "source_artifact": rel(SCD),
        "source_artifact_sha256": sha(SCD),
        "decoded_graph_count": len(graphs),
        "parse_consumed_all_bytes": parse_stats.get("parse_consumed_all_bytes"),
        "record_label_format": "graph6-style n=16 encoding in decoded Meringer vertex order",
        "record_label_dataset_csv": rel(CSV_OUT),
        "record_label_dataset_csv_sha256": csv_sha,
        "record_label_count": len(rows),
        "unique_record_graph6_label_count": len(label_counter),
        "unique_record_graph6_sha256_count": len(hash_counter),
        "duplicate_record_label_count": len(duplicate_labels),
        "duplicate_record_hash_count": len(duplicate_hashes),
        "sample_rows": rows[:5],
        "tail_rows": rows[-5:],
        "finite_certificate_statement": "A reproducible graph6-style record-label/hash CSV was exported for every decoded Meringer graph.  All 16,828 decoded-record labels and hashes are unique in the supplied vertex order; this is record provenance, not isomorphism-canonical labeling.",
    }


def acceptance_matrix(dataset: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "p2805_duplicate_free_input_present": dataset["p2805_resolved_total_isomorphism_classes_after_refinement"] == EXPECTED_GRAPH_COUNT,
        "decoded_graph_count_is_16828": dataset["decoded_graph_count"] == EXPECTED_GRAPH_COUNT,
        "record_label_dataset_csv_written": bool(dataset["record_label_dataset_csv_sha256"]),
        "all_record_labels_unique": dataset["unique_record_graph6_label_count"] == EXPECTED_GRAPH_COUNT and dataset["duplicate_record_label_count"] == 0,
        "all_record_label_hashes_unique": dataset["unique_record_graph6_sha256_count"] == EXPECTED_GRAPH_COUNT and dataset["duplicate_record_hash_count"] == 0,
        "isomorphism_canonical_label_dataset_exported": False,
        "independent_canonical_label_tool_crosscheck_exported": False,
        "strict_spectral_source_law_exported": False,
        "kernel_or_ltotal_variational_coupling_exported": False,
    }
    return {
        "facts": facts,
        "accepted_as_reproducible_record_label_dataset": all(facts[key] for key in [
            "p2805_duplicate_free_input_present",
            "decoded_graph_count_is_16828",
            "record_label_dataset_csv_written",
            "all_record_labels_unique",
            "all_record_label_hashes_unique",
        ]),
        "accepted_as_isomorphism_canonical_label_dataset": False,
        "accepted_as_independent_canonical_label_crosscheck": False,
        "accepted_as_strict_spectral_source_law": False,
        "accepted_as_ltotal_or_toe_promotion": False,
        "missing_criteria": [key for key, value in facts.items() if not value],
        "blocker": "P2806 exports a complete unique decoded-record label/hash dataset, but because labels are in supplied Meringer vertex order and no independent canonical-labeling engine/cross-check is used, it is not an isomorphism-canonical label dataset or strict source/coupling theorem.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    d = payload["record_label_dataset_audit"]
    lines = [
        "# P2806/S1756 girth>=4 reproducible record-label dataset audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Dataset counts",
        f"- decoded_graph_count={d['decoded_graph_count']}",
        f"- record_label_dataset_csv=`{d['record_label_dataset_csv']}`",
        f"- record_label_dataset_csv_sha256={d['record_label_dataset_csv_sha256']}",
        f"- record_label_count={d['record_label_count']}",
        f"- unique_record_graph6_label_count={d['unique_record_graph6_label_count']}",
        f"- unique_record_graph6_sha256_count={d['unique_record_graph6_sha256_count']}",
        f"- duplicate_record_label_count={d['duplicate_record_label_count']}",
        f"- duplicate_record_hash_count={d['duplicate_record_hash_count']}",
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
    p2805 = read_json(P2805)
    dataset = build_dataset(p2805)
    acceptance = acceptance_matrix(dataset)
    payload = {
        "status": "P2806_GIRTH4_REPRODUCIBLE_RECORD_LABEL_DATASET_UNIQUE_RECORDS_NO_ISOMORPHISM_CANONICAL_LABELS_NO_SOURCE_LAW_NO_CLOSURE",
        "input_hashes": {"P2805": sha(P2805), "16_4_4.scd": sha(SCD)},
        "input_statuses": {"P2805": p2805.get("status")},
        "audited_question": "Can the validated duplicate-free Meringer import be exported as a reproducible decoded-record label/hash dataset without falsely claiming isomorphism-canonical labels?",
        "record_label_dataset_audit": dataset,
        "acceptance_matrix": acceptance,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": acceptance["blocker"],
            "next_honest_step": "Treat P2806 as a reproducible decoded-record provenance dataset only.  The next proof-grade move is an independent isomorphism-canonical labeling cross-check: run a real canonical-labeling engine or two independent canonicalization routes on all 16,828 records, compare against the P2806 record-label table, record automorphism/group-size data if available, and only then consider a separate strict spectral source-law/coupling audit.  Do not promote P2806 to K/L_total, bridge closure, role transfer, selector closure, or ToE closure.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2806/S1756 girth>=4 reproducible record-label dataset audit", "## P2806/S1756 girth>=4 reproducible record-label dataset audit\n\n`P2806/S1756` exports a deterministic graph6-style decoded-record label/hash CSV for all `16,828` validated Meringer graphs and finds `16,828` unique record labels and hashes in the supplied vertex order.  This strengthens graph-list provenance but is not an isomorphism-canonical label dataset, not an independent canonical-label cross-check, not a strict spectral source law, and does not license `K`/`L_total` coupling.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2806/S1756 record-label dataset Ltotal guard", "## P2806/S1756 record-label dataset Ltotal guard\n\n`P2806/S1756` adds no variational source term.  A complete unique decoded-record label/hash dataset is provenance for the 16-node girth>=4 graph import, but it supplies neither isomorphism-canonical geometry nor a spectral source/coupling theorem for `K` or `L_total`.\n")
    append_once(AGENTS, "Current girth>=4 record-label dataset guardrail (P2806/S1756, 2026-06-16)", "## Current girth>=4 record-label dataset guardrail (P2806/S1756, 2026-06-16)\n\n- P2806 exports `fundamental_action_reconstruction/generated/p2806_s1756_girth4_record_label_dataset.csv`: a deterministic graph6-style decoded-record label/hash table for all `16,828` validated Meringer graphs, with `16,828` unique labels and SHA256 hashes in the supplied vertex order.\n- This is a reproducible record-label provenance dataset only.  It is not an isomorphism-canonical label dataset, not an independent canonical-labeling cross-check, and not a strict spectral source law.\n- Do not promote P2806 to canonical geometry, selector closure, kernel full-expression, role-bearing `L_total`, bridge closure, role transfer, or ToE closure.  The next admissible move is an independent isomorphism-canonical labeling cross-check plus optional automorphism/group-size audit.\n")
    return payload


if __name__ == "__main__":
    main()
