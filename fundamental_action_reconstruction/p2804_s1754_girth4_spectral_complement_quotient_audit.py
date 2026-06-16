#!/usr/bin/env python3
"""P2804/S1754: exact spectral/complement quotient audit for Meringer 16_4_4.

P2803 validated the supplied Meringer SCD import/decode layer.  P2804 performs
the next bounded computation over all 16,828 decoded graphs: exact integer
adjacency characteristic polynomials, exact integer complement characteristic
polynomials, spectral quotient collision tables, and complement-pair quotient
statistics.

This is a spectral quotient/witness matrix only.  It is not a canonical
isomorphism quotient, not a strict spectral source law, and not K/L_total or ToE
closure.
"""
from __future__ import annotations

import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import numpy as np

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import EXPECTED_GRAPH_COUNT, N, SCD, decode_scd_bytes, sha

GEN = ROOT / "generated"
P2803 = GEN / "p2803_s1753_meringer_scd_import_decode_validation_certificate.json"
OUT = GEN / "p2804_s1754_girth4_spectral_complement_quotient_audit.json"
MD = GEN / "p2804_s1754_girth4_spectral_complement_quotient_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
NEGATIVE_EXPORT_FLAGS = [
    "canonical_isomorphism_quotient_completed",
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


def adjacency_matrix(adj: tuple[tuple[int, ...], ...]) -> np.ndarray:
    matrix = np.zeros((N, N), dtype=np.int64)
    for row_index, neighbors in enumerate(adj):
        for neighbor in neighbors:
            matrix[row_index, neighbor - 1] = 1
    return matrix


def charpoly_coefficients_from_matrix(matrix: np.ndarray) -> tuple[int, ...]:
    """Return monic characteristic-polynomial coefficients exactly.

    Coefficients are ordered as `(1, c1, ..., c16)` for
    `lambda**16 + c1*lambda**15 + ... + c16`.  Newton identities are applied to
    exact integer traces of powers of a 16x16 integer matrix.
    """
    power = matrix.copy()
    traces: list[int] = []
    for exponent in range(1, N + 1):
        if exponent > 1:
            power = power @ matrix
        traces.append(int(np.trace(power)))

    elementary = [1]
    for k in range(1, N + 1):
        numerator = sum(((-1) ** (i - 1)) * elementary[k - i] * traces[i - 1] for i in range(1, k + 1))
        if numerator % k != 0:
            raise ArithmeticError(f"Newton identity numerator not divisible at k={k}: {numerator}")
        elementary.append(numerator // k)
    return tuple(((-1) ** k) * elementary[k] for k in range(0, N + 1))


def complement_matrix(matrix: np.ndarray) -> np.ndarray:
    return np.ones((N, N), dtype=np.int64) - np.eye(N, dtype=np.int64) - matrix


def digest_coefficients(coefficients: tuple[int, ...]) -> str:
    return hashlib.sha256(json.dumps(coefficients, separators=(",", ":")).encode("utf-8")).hexdigest()


def class_summary(counter: Counter[tuple[int, ...]], examples: dict[tuple[int, ...], list[int]]) -> dict[str, Any]:
    sizes = sorted(counter.values(), reverse=True)
    collision_classes = [poly for poly, count in counter.items() if count > 1]
    top = sorted(counter.items(), key=lambda item: (-item[1], digest_coefficients(item[0])))[:10]
    return {
        "class_count": len(counter),
        "singleton_class_count": sum(1 for size in sizes if size == 1),
        "collision_class_count": len(collision_classes),
        "max_class_size": sizes[0] if sizes else 0,
        "top_collision_classes": [
            {
                "size": count,
                "example_indices_0_based": examples[poly][:8],
                "charpoly_sha256": digest_coefficients(poly),
                "charpoly_coefficients": list(poly),
            }
            for poly, count in top
            if count > 1
        ],
    }


def build_audit(p2803: dict[str, Any]) -> dict[str, Any]:
    graphs, parse_stats = decode_scd_bytes(SCD.read_bytes())
    adjacency_counter: Counter[tuple[int, ...]] = Counter()
    complement_counter: Counter[tuple[int, ...]] = Counter()
    pair_counter: Counter[tuple[tuple[int, ...], tuple[int, ...]]] = Counter()
    adjacency_examples: dict[tuple[int, ...], list[int]] = defaultdict(list)
    complement_examples: dict[tuple[int, ...], list[int]] = defaultdict(list)
    pair_examples: dict[tuple[tuple[int, ...], tuple[int, ...]], list[int]] = defaultdict(list)
    complement_degree_histogram: Counter[int] = Counter()

    for index, graph in enumerate(graphs):
        matrix = adjacency_matrix(graph)
        comp = complement_matrix(matrix)
        adjacency_poly = charpoly_coefficients_from_matrix(matrix)
        complement_poly = charpoly_coefficients_from_matrix(comp)
        pair = (adjacency_poly, complement_poly)
        adjacency_counter[adjacency_poly] += 1
        complement_counter[complement_poly] += 1
        pair_counter[pair] += 1
        if len(adjacency_examples[adjacency_poly]) < 8:
            adjacency_examples[adjacency_poly].append(index)
        if len(complement_examples[complement_poly]) < 8:
            complement_examples[complement_poly].append(index)
        if len(pair_examples[pair]) < 8:
            pair_examples[pair].append(index)
        complement_degree_histogram.update(int(deg) for deg in comp.sum(axis=1))

    pair_sizes = sorted(pair_counter.values(), reverse=True)
    top_pairs = sorted(pair_counter.items(), key=lambda item: (-item[1], digest_coefficients(item[0][0]), digest_coefficients(item[0][1])))[:10]
    return {
        "expected_graph_count": EXPECTED_GRAPH_COUNT,
        "p2803_status": p2803.get("status"),
        "p2803_accepts_import_decode": p2803.get("acceptance_matrix", {}).get("accepted_as_meringer_scd_import_decode_validation"),
        "source_artifact": rel(SCD),
        "source_artifact_sha256": sha(SCD),
        "decoded_graph_count": len(graphs),
        "p2803_parse_decoded_entry_count": parse_stats.get("decoded_entry_count"),
        "adjacency_spectral_quotient": class_summary(adjacency_counter, adjacency_examples),
        "complement_spectral_quotient": class_summary(complement_counter, complement_examples),
        "adjacency_complement_pair_quotient": {
            "class_count": len(pair_counter),
            "singleton_class_count": sum(1 for size in pair_sizes if size == 1),
            "collision_class_count": sum(1 for size in pair_sizes if size > 1),
            "max_class_size": pair_sizes[0] if pair_sizes else 0,
            "top_collision_classes": [
                {
                    "size": count,
                    "example_indices_0_based": pair_examples[pair][:8],
                    "adjacency_charpoly_sha256": digest_coefficients(pair[0]),
                    "complement_charpoly_sha256": digest_coefficients(pair[1]),
                    "adjacency_charpoly_coefficients": list(pair[0]),
                    "complement_charpoly_coefficients": list(pair[1]),
                }
                for pair, count in top_pairs
                if count > 1
            ],
        },
        "complement_degree_histogram": dict(sorted(complement_degree_histogram.items())),
        "finite_certificate_statement": "All 16,828 decoded graphs were processed through exact integer adjacency and complement characteristic polynomials.  The resulting spectral quotient is finite and reproducible, but it is not an isomorphism quotient and does not by itself export a strict spectral source law.",
    }


def acceptance_matrix(audit: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "p2803_validated_import_decode_input": audit["p2803_accepts_import_decode"] is True,
        "decoded_graph_count_is_16828": audit["decoded_graph_count"] == EXPECTED_GRAPH_COUNT,
        "adjacency_charpoly_rows_processed": sum(
            row["size"] for row in audit["adjacency_spectral_quotient"]["top_collision_classes"]
        ) >= 0 and audit["adjacency_spectral_quotient"]["class_count"] > 0,
        "complement_charpoly_rows_processed": audit["complement_spectral_quotient"]["class_count"] > 0,
        "pair_quotient_rows_processed": audit["adjacency_complement_pair_quotient"]["class_count"] > 0,
        "complement_degree_is_11_for_all_vertices": audit["complement_degree_histogram"] == {11: EXPECTED_GRAPH_COUNT * N} or audit["complement_degree_histogram"] == {"11": EXPECTED_GRAPH_COUNT * N},
        "canonical_isomorphism_quotient_completed": False,
        "strict_spectral_source_law_exported": False,
        "kernel_or_ltotal_variational_coupling_exported": False,
    }
    return {
        "facts": facts,
        "accepted_as_exact_spectral_complement_quotient_audit": all(facts[key] for key in [
            "p2803_validated_import_decode_input",
            "decoded_graph_count_is_16828",
            "adjacency_charpoly_rows_processed",
            "complement_charpoly_rows_processed",
            "pair_quotient_rows_processed",
            "complement_degree_is_11_for_all_vertices",
        ]),
        "accepted_as_canonical_isomorphism_quotient": False,
        "accepted_as_full_16node_canonical_generator_certificate": False,
        "accepted_as_strict_spectral_source_law": False,
        "accepted_as_canonical_nadsoliton_geometry_source": False,
        "accepted_as_ltotal_or_toe_promotion": False,
        "missing_criteria": [key for key, value in facts.items() if not value],
        "blocker": "P2804 computes exact spectral/complement quotient data over the validated 16,828 graph import, but spectral collisions remain and no canonical isomorphism quotient, strict source law, or K/L_total coupling theorem is exported.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["girth4_spectral_complement_quotient_audit"]
    adj = audit["adjacency_spectral_quotient"]
    comp = audit["complement_spectral_quotient"]
    pair = audit["adjacency_complement_pair_quotient"]
    lines = [
        "# P2804/S1754 girth>=4 spectral/complement quotient audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Exact quotient counts",
        f"- decoded_graph_count={audit['decoded_graph_count']}",
        f"- adjacency_charpoly_class_count={adj['class_count']}",
        f"- adjacency_charpoly_collision_class_count={adj['collision_class_count']}",
        f"- adjacency_charpoly_max_class_size={adj['max_class_size']}",
        f"- complement_charpoly_class_count={comp['class_count']}",
        f"- complement_charpoly_collision_class_count={comp['collision_class_count']}",
        f"- complement_charpoly_max_class_size={comp['max_class_size']}",
        f"- adjacency_complement_pair_class_count={pair['class_count']}",
        f"- adjacency_complement_pair_collision_class_count={pair['collision_class_count']}",
        f"- adjacency_complement_pair_max_class_size={pair['max_class_size']}",
        f"- complement_degree_histogram={audit['complement_degree_histogram']}",
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
    p2803 = read_json(P2803)
    audit = build_audit(p2803)
    acceptance = acceptance_matrix(audit)
    payload = {
        "status": "P2804_GIRTH4_SPECTRAL_COMPLEMENT_QUOTIENT_AUDIT_NO_ISOMORPHISM_NO_SOURCE_LAW_NO_CLOSURE",
        "input_hashes": {"P2803": sha(P2803), "16_4_4.scd": sha(SCD)},
        "input_statuses": {"P2803": p2803.get("status")},
        "audited_question": "What exact adjacency/complement characteristic-polynomial quotient and collision matrix is induced by the validated 16,828-entry Meringer girth>=4 graph import?",
        "girth4_spectral_complement_quotient_audit": audit,
        "acceptance_matrix": acceptance,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": acceptance["blocker"],
            "next_honest_step": "Treat P2804 as an exact spectral/complement quotient witness, not canonical geometry.  The next proof-grade move is a canonical isomorphism quotient or collision-refinement audit for the residual spectral collision classes: run a graph-isomorphism/canonical-label toolchain on each non-singleton adjacency-complement pair class, record canonical labels, automorphism/group-size data if available, and produce a collision-resolved quotient table.  Do not promote the spectral quotient itself to a strict spectral source law, K/L_total coupling, role transfer, bridge closure, selector closure, or ToE closure.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2804/S1754 girth>=4 spectral/complement quotient audit", "## P2804/S1754 girth>=4 spectral/complement quotient audit\n\n`P2804/S1754` computes exact integer adjacency and complement characteristic-polynomial quotients over all `16,828` decoded Meringer `16_4_4.scd` graphs validated by P2803.  The audit records finite spectral class counts, collision counts, maximum collision sizes, adjacency/complement pair classes, and confirms complement degree `11` on every vertex.  This is a spectral/complement quotient witness only; it is not a canonical isomorphism quotient, not a strict spectral source law, and does not license `K`/`L_total` coupling.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2804/S1754 spectral/complement quotient Ltotal guard", "## P2804/S1754 spectral/complement quotient Ltotal guard\n\n`P2804/S1754` adds no variational source term.  Exact characteristic-polynomial quotient data over the validated girth>=4 graph import is spectral evidence, but without a strict source law and coupling theorem it cannot promote `K`/`L_total`, role transfer, bridge closure, or ToE closure.\n")
    append_once(AGENTS, "Current girth>=4 spectral/complement quotient guardrail (P2804/S1754, 2026-06-16)", "## Current girth>=4 spectral/complement quotient guardrail (P2804/S1754, 2026-06-16)\n\n- P2804 computes exact integer adjacency and complement characteristic-polynomial quotient data over all `16,828` decoded Meringer girth>=4 graphs validated by P2803, including collision classes and adjacency/complement pair classes.\n- This is not a canonical isomorphism quotient and not a strict spectral source law; spectral collisions and missing source/coupling theorems remain blockers.\n- Do not promote P2804 to canonical geometry, selector closure, kernel full-expression, role-bearing `L_total`, bridge closure, role transfer, or ToE closure.  The next admissible move is canonical-label/isomorphism refinement of the residual spectral collision classes, or a separate strict source-law/coupling theorem.\n")
    return payload


if __name__ == "__main__":
    main()
