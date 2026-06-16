#!/usr/bin/env python3
"""P2807/S1757: pynauty canonical-label toolchain gate.

P2806 exported a deterministic decoded-record label/hash table, but explicitly
stopped short of isomorphism-canonical labels.  P2807 is the next honest
canonical-label step: check whether a real pynauty/nauty canonical-label stack is
available in this execution environment, and if it is available, compute compact
canonical certificate/hash and automorphism-size statistics for all 16,828
validated Meringer graphs.  If the stack is not importable, P2807 records a
bounded toolchain blocker rather than pretending that record labels are canonical
labels.

This file deliberately does not commit a 16,828-row canonical-label CSV.  It
writes only a compact JSON/MD certificate so review diffs stay small.
"""
from __future__ import annotations

import hashlib
import importlib
import json
from collections import Counter
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import EXPECTED_GRAPH_COUNT, SCD, decode_scd_bytes, sha
from p2806_s1756_girth4_reproducible_record_label_dataset_audit import graph6_record_label, sha256_text

GEN = ROOT / "generated"
P2806 = GEN / "p2806_s1756_girth4_reproducible_record_label_dataset_audit.json"
OUT = GEN / "p2807_s1757_pynauty_canonical_label_toolchain_gate.json"
MD = GEN / "p2807_s1757_pynauty_canonical_label_toolchain_gate.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
NEGATIVE_EXPORT_FLAGS = [
    "pynauty_stack_available_in_current_environment",
    "isomorphism_canonical_label_dataset_exported",
    "canonical_geometry_source_exported",
    "strict_spectral_source_law_exported",
    "kernel_geometry_closure_exported",
    "role_bearing_ltotal_promoted",
    "bridge_closure_exported",
    "role_transfer_started",
    "selector_closure_exported",
    "toe_closure_exported",
]


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"missing": True, "path": rel(path)}


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def pynauty_probe() -> dict[str, Any]:
    try:
        module = importlib.import_module("pynauty")
    except Exception as exc:  # pragma: no cover - exact exception depends on environment
        return {
            "available": False,
            "import_error_type": type(exc).__name__,
            "import_error": str(exc),
            "available_symbols": [],
        }
    wanted = ["Graph", "certificate", "canon_label", "autgrp"]
    return {
        "available": True,
        "version": getattr(module, "__version__", "unknown"),
        "available_symbols": [name for name in wanted if hasattr(module, name)],
        "missing_symbols": [name for name in wanted if not hasattr(module, name)],
    }


def to_pynauty_graph(module: Any, adj: tuple[tuple[int, ...], ...]) -> Any:
    graph = module.Graph(number_of_vertices=len(adj), directed=False)
    for vertex, row in enumerate(adj):
        graph.connect_vertex(vertex, [neighbor - 1 for neighbor in row if neighbor - 1 > vertex])
    return graph


def compute_pynauty_compact_audit(module: Any) -> dict[str, Any]:
    graphs, _ = decode_scd_bytes(SCD.read_bytes())
    certificate_counter: Counter[str] = Counter()
    record_to_certificate_agreement_counter: Counter[bool] = Counter()
    automorphism_size_counter: Counter[str] = Counter()
    sample_rows = []
    for index, adj in enumerate(graphs):
        graph = to_pynauty_graph(module, adj)
        certificate = module.certificate(graph)
        certificate_sha = sha256_bytes(certificate)
        certificate_counter[certificate_sha] += 1
        record_hash = sha256_text(graph6_record_label(adj))
        record_to_certificate_agreement_counter[record_hash == certificate_sha] += 1
        autgrp = module.autgrp(graph)
        automorphism_size = str(autgrp[1]) if isinstance(autgrp, tuple) and len(autgrp) > 1 else "unknown"
        automorphism_size_counter[automorphism_size] += 1
        if index < 10:
            sample_rows.append({
                "graph_index_1_based": index + 1,
                "pynauty_certificate_sha256": certificate_sha,
                "automorphism_group_size": automorphism_size,
            })
    duplicate_certificate_classes = sum(1 for count in certificate_counter.values() if count > 1)
    return {
        "decoded_graph_count": len(graphs),
        "pynauty_certificate_hash_class_count": len(certificate_counter),
        "pynauty_duplicate_certificate_class_count": duplicate_certificate_classes,
        "pynauty_certificate_max_class_size": max(certificate_counter.values()) if certificate_counter else 0,
        "record_hash_equals_certificate_hash_histogram": dict(record_to_certificate_agreement_counter),
        "automorphism_group_size_histogram": dict(sorted(automorphism_size_counter.items(), key=lambda item: item[0])),
        "sample_rows": sample_rows,
    }


def build_gate(p2806: dict[str, Any]) -> dict[str, Any]:
    probe = pynauty_probe()
    p2806_audit = p2806.get("record_label_dataset_audit", {})
    gate: dict[str, Any] = {
        "expected_graph_count": EXPECTED_GRAPH_COUNT,
        "p2806_status": p2806.get("status"),
        "p2806_record_label_count": p2806_audit.get("record_label_count"),
        "p2806_unique_record_graph6_label_count": p2806_audit.get("unique_record_graph6_label_count"),
        "source_artifact": rel(SCD),
        "source_artifact_sha256": sha(SCD),
        "pynauty_probe": probe,
        "diff_size_policy": "compact JSON/MD only; no 16,828-row canonical CSV committed",
    }
    if probe["available"] and not probe.get("missing_symbols"):
        module = importlib.import_module("pynauty")
        gate["pynauty_compact_audit"] = compute_pynauty_compact_audit(module)
        gate["finite_certificate_statement"] = "pynauty is available; compact canonical certificate/hash and automorphism-size statistics were computed without committing a large row-level CSV."
    else:
        gate["pynauty_compact_audit"] = None
        gate["finite_certificate_statement"] = "pynauty/nauty canonical-label stack is not importable with the required Graph/certificate/canon_label/autgrp symbols in this environment; canonical-label closure remains blocked instead of being inferred from P2806 record labels."
    return gate


def acceptance_matrix(gate: dict[str, Any]) -> dict[str, Any]:
    probe = gate["pynauty_probe"]
    audit = gate.get("pynauty_compact_audit")
    facts = {
        "p2806_unique_record_dataset_present": gate["p2806_record_label_count"] == EXPECTED_GRAPH_COUNT and gate["p2806_unique_record_graph6_label_count"] == EXPECTED_GRAPH_COUNT,
        "pynauty_stack_available_with_required_symbols": probe["available"] and not probe.get("missing_symbols"),
        "pynauty_audit_covers_16828_graphs": bool(audit) and audit.get("decoded_graph_count") == EXPECTED_GRAPH_COUNT,
        "pynauty_certificate_hashes_unique": bool(audit) and audit.get("pynauty_certificate_hash_class_count") == EXPECTED_GRAPH_COUNT,
        "row_level_canonical_csv_committed": False,
        "strict_spectral_source_law_exported": False,
        "kernel_or_ltotal_variational_coupling_exported": False,
    }
    accepted = all(facts[key] for key in [
        "p2806_unique_record_dataset_present",
        "pynauty_stack_available_with_required_symbols",
        "pynauty_audit_covers_16828_graphs",
        "pynauty_certificate_hashes_unique",
    ])
    return {
        "facts": facts,
        "accepted_as_pynauty_canonical_toolchain_gate": facts["p2806_unique_record_dataset_present"],
        "accepted_as_compact_pynauty_canonical_certificate_audit": accepted,
        "accepted_as_full_row_level_canonical_label_dataset": False,
        "accepted_as_strict_spectral_source_law": False,
        "accepted_as_ltotal_or_toe_promotion": False,
        "missing_criteria": [key for key, value in facts.items() if not value],
        "blocker": "P2807 checks the pynauty canonical-label stack and preserves a compact-diff boundary.  If pynauty is unavailable, canonical labeling remains a toolchain blocker; if available, the compact certificate audit is still not a strict spectral source/coupling theorem or K/L_total closure.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    gate = payload["pynauty_canonical_label_toolchain_gate"]
    probe = gate["pynauty_probe"]
    lines = [
        "# P2807/S1757 pynauty canonical-label toolchain gate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Toolchain probe",
        f"- pynauty_available={probe['available']}",
        f"- available_symbols={probe.get('available_symbols', [])}",
        f"- missing_symbols={probe.get('missing_symbols', [])}",
        f"- import_error_type={probe.get('import_error_type')}",
        "",
        "## Decision",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    audit = gate.get("pynauty_compact_audit")
    if audit:
        lines[8:8] = [
            "",
            "## Compact pynauty audit counts",
            f"- decoded_graph_count={audit['decoded_graph_count']}",
            f"- pynauty_certificate_hash_class_count={audit['pynauty_certificate_hash_class_count']}",
            f"- pynauty_duplicate_certificate_class_count={audit['pynauty_duplicate_certificate_class_count']}",
            f"- pynauty_certificate_max_class_size={audit['pynauty_certificate_max_class_size']}",
        ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    p2806 = read_json(P2806)
    gate = build_gate(p2806)
    acceptance = acceptance_matrix(gate)
    status = "P2807_PYNAUTY_CANONICAL_LABEL_TOOLCHAIN_GATE_BLOCKED_NO_CANONICAL_LABELS_NO_SOURCE_LAW_NO_CLOSURE"
    if acceptance["accepted_as_compact_pynauty_canonical_certificate_audit"]:
        status = "P2807_PYNAUTY_COMPACT_CANONICAL_CERTIFICATE_AUDIT_NO_ROW_CSV_NO_SOURCE_LAW_NO_CLOSURE"
    payload = {
        "status": status,
        "input_hashes": {"P2806": sha(P2806), "16_4_4.scd": sha(SCD)},
        "input_statuses": {"P2806": p2806.get("status")},
        "audited_question": "Is a real pynauty/nauty canonical-label stack available for a compact canonical-certificate audit without committing another 16,828-row diff?",
        "pynauty_canonical_label_toolchain_gate": gate,
        "acceptance_matrix": acceptance,
        "decision": {
            "negative_export_flags": {
                flag: (bool(acceptance["facts"]["pynauty_stack_available_with_required_symbols"]) if flag == "pynauty_stack_available_in_current_environment" else False)
                for flag in NEGATIVE_EXPORT_FLAGS
            },
            "reason": acceptance["blocker"],
            "next_honest_step": ("Pynauty is available and P2807 produced unique compact canonical certificate hashes, so the next honest step is P2808: a small canonical digest/automorphism manifest with any row-level CSV kept out of the review diff, followed only then by a separate strict spectral source-law/coupling audit.  Do not promote P2807 to K/L_total, bridge closure, role transfer, selector closure, or ToE closure." if acceptance["accepted_as_compact_pynauty_canonical_certificate_audit"] else "If pynauty remains blocked, obtain a review-safe canonical-label toolchain artifact outside the repository diff (wheel/cache/container or nauty executable) and rerun P2807.  If pynauty becomes available and P2807 produces unique canonical certificate hashes, the next step is a small row-level canonical digest manifest plus automorphism-size table, followed only then by a separate strict spectral source-law/coupling audit.  Do not promote P2807 to K/L_total, bridge closure, role transfer, selector closure, or ToE closure."),
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2807/S1757 pynauty canonical-label toolchain gate", "## P2807/S1757 pynauty canonical-label toolchain gate\n\n`P2807/S1757` checks whether a real `pynauty`/nauty canonical-label stack is importable for the validated `16,828`-graph Meringer dataset while preserving a compact-diff policy.  In the current environment this is a toolchain gate rather than a source-law theorem: if `pynauty` is unavailable, canonical labels remain blocked; if it is available, the resulting compact certificate audit still does not license `K`/`L_total` coupling without a separate spectral source/coupling theorem.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2807/S1757 pynauty canonical-label Ltotal guard", "## P2807/S1757 pynauty canonical-label Ltotal guard\n\n`P2807/S1757` adds no variational source term.  A `pynauty` canonical-label/certificate audit, even when the toolchain is available, is graph-isomorphism provenance only; it is not a spectral action/source law and cannot promote `K`, `L_total`, bridge closure, role transfer, selector closure, or ToE closure.\n")
    append_once(AGENTS, "Current pynauty canonical-label toolchain guardrail (P2807/S1757, 2026-06-16)", "## Current pynauty canonical-label toolchain guardrail (P2807/S1757, 2026-06-16)\n\n- P2807 is the review-size-safe canonical-label toolchain gate after P2806: it checks for an importable `pynauty`/nauty stack and, only if present, computes compact canonical certificate/hash and automorphism-size statistics without committing another `16,828`-row CSV.\n- If `pynauty` is unavailable, do not infer canonical labels from P2806 record labels.  If `pynauty` is available, treat the output as graph-isomorphism provenance only, not a strict spectral source law.\n- Do not promote P2807 to canonical geometry, selector closure, kernel full-expression, role-bearing `L_total`, bridge closure, role transfer, or ToE closure.  The next admissible move is to supply the missing canonical-label toolchain artifact or, after a successful compact audit, a small canonical digest/automorphism manifest.\n")
    return payload


if __name__ == "__main__":
    main()
