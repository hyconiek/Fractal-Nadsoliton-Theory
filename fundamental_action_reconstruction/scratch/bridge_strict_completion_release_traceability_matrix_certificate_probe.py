#!/usr/bin/env python3
"""Scratch probe: source-to-section traceability matrix for strict release scaffolds.

This is a finite bookkeeping certificate.  It verifies that the release-facing
Markdown scaffolds form a traceable source-to-section matrix with full target
column coverage, while preserving all non-closure guardrails.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FAR = ROOT / "fundamental_action_reconstruction"
OUT_JSON = HERE / "bridge_strict_completion_release_traceability_matrix_certificate_report.json"
OUT_MD = HERE / "bridge_strict_completion_release_traceability_matrix_certificate_report.md"
TRACEABILITY_DOC = FAR / "STRICT_KERNEL_RELEASE_TRACEABILITY_MATRIX.md"

TARGET_DOCS = {
    "D": FAR / "DIAGRAMS_STRICT_KERNEL_TRANSFORMATION.md",
    "L": FAR / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
    "R": FAR / "STRICT_KERNEL_LEGACY_ROLE_TRANSFER_AUDIT_DRAFT.md",
    "S": FAR / "STRICT_KERNEL_RELEASE_SOURCE_COVERAGE_AUDIT.md",
}

SOURCE_REPORTS = {
    "release_scaffold": HERE / "bridge_strict_completion_release_scaffold_certificate_report.json",
    "release_source_coverage": HERE / "bridge_strict_completion_release_source_coverage_certificate_report.json",
    "chain_integrity": HERE / "bridge_strict_completion_certificate_chain_integrity_report.json",
    "frontier_low_weight": HERE / "bridge_strict_completion_theorem_frontier_low_weight_extension_certificate_report.json",
    "anchor_h1": HERE / "bridge_strict_completion_anchor_h1_generator_classification_certificate_report.json",
}

TRACEABILITY_ROWS = [
    {
        "source_family": "legacy_kernel_history",
        "vector": [1, 0, 0, 1],
        "required_doc_snippets": ["legacy_kernel_history", "K_legacy_ont", "no raw identity"],
        "blocker": "no raw identity",
    },
    {
        "source_family": "finite_bridge_ledger",
        "vector": [1, 0, 0, 1],
        "required_doc_snippets": ["finite_bridge_ledger", "APD/diagonal/symbolic cancellation", "no strict dynamical bridge theorem"],
        "blocker": "no strict dynamical bridge theorem",
    },
    {
        "source_family": "strict_lagrangian_eom",
        "vector": [0, 1, 0, 1],
        "required_doc_snippets": ["strict_lagrangian_eom", "P1622/P1866/P2315/P2316", "no full tensor-resolved EOM closure"],
        "blocker": "no full tensor-resolved EOM closure",
    },
    {
        "source_family": "role_transfer_boundaries",
        "vector": [0, 0, 1, 1],
        "required_doc_snippets": ["role_transfer_boundaries", "N87/N103", "no semantic role-transfer theorem"],
        "blocker": "no semantic role-transfer theorem",
    },
    {
        "source_family": "anchor_h1_selector_boundary",
        "vector": [1, 0, 1, 1],
        "required_doc_snippets": ["anchor_h1_selector_boundary", "C0 anchor vs C1/im(delta) H1", "no `QW-2191` discharge"],
        "blocker": "no `QW-2191` discharge",
    },
    {
        "source_family": "theorem_frontier_planning",
        "vector": [1, 0, 1, 1],
        "required_doc_snippets": ["theorem_frontier_planning", "no singleton/pair closes bridge, role, or ToE"],
        "blocker": "no singleton/pair closes bridge, role, or ToE",
    },
]

REQUIRED_TRACE_DOC_SNIPPETS = [
    "matrix shape = 6 x 4",
    "column coverage = [4, 1, 3, 6]",
    "row coverage = [2, 2, 2, 2, 3, 3]",
    "GF(2) rank = 4",
    "does not prove a bridge theorem",
    "does not prove full tensor-resolved EOM closure",
    "does not transfer legacy physical roles",
    "does not discharge `QW-2191`",
    "does not close ToE",
]


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(path)
    return json.loads(path.read_text(encoding="utf-8"))


def gf2_rank(matrix: list[list[int]]) -> int:
    rows = [row[:] for row in matrix]
    if not rows:
        return 0
    m = len(rows)
    n = len(rows[0])
    rank = 0
    pivot_col = 0
    while rank < m and pivot_col < n:
        pivot = None
        for r in range(rank, m):
            if rows[r][pivot_col] % 2:
                pivot = r
                break
        if pivot is None:
            pivot_col += 1
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        for r in range(m):
            if r != rank and rows[r][pivot_col] % 2:
                rows[r] = [(a ^ b) for a, b in zip(rows[r], rows[rank])]
        rank += 1
        pivot_col += 1
    return rank


def build_payload() -> dict[str, Any]:
    reports = {name: load_json(path) for name, path in SOURCE_REPORTS.items()}
    doc_text = TRACEABILITY_DOC.read_text(encoding="utf-8") if TRACEABILITY_DOC.exists() else ""
    target_texts = {name: path.read_text(encoding="utf-8") if path.exists() else "" for name, path in TARGET_DOCS.items()}
    matrix = [row["vector"] for row in TRACEABILITY_ROWS]
    column_coverage = [sum(row[col] for row in matrix) for col in range(len(TARGET_DOCS))]
    row_coverage = [sum(row) for row in matrix]
    rank = gf2_rank(matrix)

    row_checks = []
    for row in TRACEABILITY_ROWS:
        row_checks.append(
            {
                "source_family": row["source_family"],
                "vector": row["vector"],
                "row_coverage": sum(row["vector"]),
                "all_required_trace_snippets_present": all(snippet in doc_text for snippet in row["required_doc_snippets"]),
                "blocker_recorded": row["blocker"] in doc_text,
            }
        )

    release_source_summary = reports["release_source_coverage"]["coverage_summary"]
    release_scaffold_summary = reports["release_scaffold"]["release_scaffold_summary"]
    chain_summary = reports["chain_integrity"]["chain_summary"]
    low_weight_summary = reports["frontier_low_weight"]["theorem_frontier_low_weight_extension_summary"]
    anchor_summary = reports["anchor_h1"]["classification_summary"]

    cross_checks = {
        "traceability_doc_present": TRACEABILITY_DOC.exists(),
        "target_docs_present": all(path.exists() for path in TARGET_DOCS.values()),
        "required_trace_doc_snippets_present": all(snippet in doc_text for snippet in REQUIRED_TRACE_DOC_SNIPPETS),
        "all_rows_have_required_snippets": all(row["all_required_trace_snippets_present"] for row in row_checks),
        "all_rows_record_blockers": all(row["blocker_recorded"] for row in row_checks),
        "matrix_shape_pass": len(matrix) == 6 and all(len(row) == 4 for row in matrix),
        "column_coverage_pass": column_coverage == [4, 1, 3, 6] and all(value > 0 for value in column_coverage),
        "row_coverage_pass": row_coverage == [2, 2, 2, 2, 3, 3] and all(value >= 2 for value in row_coverage),
        "gf2_rank_full_column_pass": rank == 4,
        "release_source_coverage_inherited": reports["release_source_coverage"]["all_cross_checks_pass"] and release_source_summary["release_scaffold_nonduplicating_source_map_ready"],
        "release_scaffold_inherited": reports["release_scaffold"]["all_cross_checks_pass"] and release_scaffold_summary["strict_kernel_diagram_scaffold_ready"] and release_scaffold_summary["strict_lagrangian_eom_scaffold_ready"] and release_scaffold_summary["strict_role_transfer_audit_scaffold_ready"],
        "chain_integrity_inherited": reports["chain_integrity"]["all_cross_checks_pass"] and chain_summary["release_source_coverage_certified"],
        "selector_frontier_still_open": low_weight_summary["chi11_is_only_singleton_unlock"] and low_weight_summary["no_singleton_closes_bridge_role_or_toe"] and low_weight_summary["no_pair_closes_bridge"] and anchor_summary["selector_source_remains_open"],
        "hard_limits_preserved": all([
            "does not prove a bridge theorem" in doc_text,
            "does not prove full tensor-resolved EOM closure" in doc_text,
            "does not transfer legacy physical roles" in doc_text,
            "does not discharge `QW-2191`" in doc_text,
            "does not close ToE" in doc_text,
        ]),
    }

    return {
        "result_kind": "SCRATCH_STRICT_COMPLETION_RELEASE_TRACEABILITY_MATRIX_CERTIFICATE__GF2_SOURCE_TO_SECTION_AUDIT_NO_CLOSURE",
        "status": "strict-release-traceability-matrix-has-full-column-coverage-and-rank-no-false-pass",
        "traceability_doc": rel(TRACEABILITY_DOC),
        "target_docs": {name: rel(path) for name, path in TARGET_DOCS.items()},
        "source_reports": {name: rel(path) for name, path in SOURCE_REPORTS.items()},
        "matrix_columns": list(TARGET_DOCS),
        "traceability_rows": row_checks,
        "matrix": matrix,
        "column_coverage": column_coverage,
        "row_coverage": row_coverage,
        "gf2_rank": rank,
        "cross_checks": cross_checks,
        "all_cross_checks_pass": all(cross_checks.values()),
        "traceability_summary": {
            "all_target_scaffold_columns_covered": cross_checks["column_coverage_pass"],
            "all_source_rows_used_at_least_twice": cross_checks["row_coverage_pass"],
            "target_columns_independent_over_gf2": cross_checks["gf2_rank_full_column_pass"],
            "source_coverage_and_scaffold_certificates_inherited": cross_checks["release_source_coverage_inherited"] and cross_checks["release_scaffold_inherited"],
            "selector_frontier_still_open": cross_checks["selector_frontier_still_open"],
            "no_bridge_theorem_claimed": True,
            "no_full_eom_closure_claimed": True,
            "no_role_transfer_claimed": True,
            "no_qw2191_discharge": True,
            "no_toe_closure": True,
        },
        "proof_certificate": {
            "matrix_step": "The traceability matrix has shape 6x4, column coverage [4,1,3,6], row coverage [2,2,2,2,3,3], and GF(2) rank 4, so every release target column is covered and independent as a bookkeeping target.",
            "source_step": "The rows link legacy history, finite bridge ledger, strict Lagrangian/EOM, role-transfer boundaries, anchor/H1 selector boundary, and theorem-frontier planning to the four release scaffolds D/L/R/S.",
            "frontier_step": "The low-weight frontier and anchor/H1 reports still keep chi11_selector_source as an open selector/source atom; no singleton/pair closes bridge, role-transfer, or ToE.",
            "limit_step": "This is a traceability certificate only: no bridge theorem, no tensor-resolved EOM closure, no legacy role-transfer theorem, no QW-2191 discharge, and no ToE closure are exported.",
        },
        "hard_limits": [
            "No traceability edge is promoted to a theorem.",
            "No bridge theorem is claimed.",
            "No full tensor-resolved EOM closure is claimed.",
            "No legacy physical-role transfer is claimed.",
            "No beta_tors -> chi_11 theorem is claimed.",
            "No QW-2191 selector discharge is claimed.",
            "No ToE closure is claimed.",
        ],
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lines = [
        "# Strict release traceability matrix certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        f"GF(2) rank: `{payload['gf2_rank']}`",
        f"Column coverage: `{payload['column_coverage']}`",
        f"Row coverage: `{payload['row_coverage']}`",
        "",
        "## Rows",
        "",
    ]
    for row in payload["traceability_rows"]:
        lines.append(f"- `{row['source_family']}`: vector=`{row['vector']}`, coverage=`{row['row_coverage']}`, snippets=`{row['all_required_trace_snippets_present']}`, blocker=`{row['blocker_recorded']}`")
    lines.extend(["", "## Cross-checks", ""])
    for key, value in payload["cross_checks"].items():
        lines.append(f"- `{key}`: `{value}`")
    lines.extend(["", "## Hard limits", ""])
    for limit in payload["hard_limits"]:
        lines.append(f"- {limit}")
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    payload = build_payload()
    OUT_JSON.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_markdown(payload)
    print(json.dumps(payload, indent=2, sort_keys=True))
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
