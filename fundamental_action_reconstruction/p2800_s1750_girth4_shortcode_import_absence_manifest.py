#!/usr/bin/env python3
"""P2800/S1750: girth>=4 shortcode import absence/readiness manifest.

P2799 closed the exact local spectral/orbit table for the six currently present
girth>=4 witnesses.  The honest next proof-grade object would be the actual
linked 16,828-class shortcode/graph-list artifact.  P2800 therefore does not
invent another subclass: it performs a finite repository-local import-readiness
and absence audit for that exact artifact obligation.

This is deliberately negative.  It records that the external girth>=4 target is
known, but no local graph-list/shortcode artifact with 16,828 rows is present in
the current repository, so no full subtarget quotient or spectral audit can be
claimed yet.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
P2796 = GEN / "p2796_s1746_full_generator_artifact_absence_provenance_gate.json"
P2798 = GEN / "p2798_s1748_external_girth4_subtarget_local_girth_filter_certificate.json"
P2799 = GEN / "p2799_s1749_local_girth4_spectral_orbit_table_certificate.json"
OUT = GEN / "p2800_s1750_girth4_shortcode_import_absence_manifest.json"
MD = GEN / "p2800_s1750_girth4_shortcode_import_absence_manifest.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

EXPECTED_GIRTH4_CLASS_COUNT = 16_828
NEGATIVE_EXPORT_FLAGS = [
    "girth4_shortcode_graph_list_imported",
    "girth4_16828_class_quotient_completed",
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
CANDIDATE_SUFFIXES = {".g6", ".graph6", ".short", ".shortcode", ".sct", ".txt", ".json", ".jsonl"}
SEARCH_ROOTS = [ROOT / "generated", ROOT / "data", ROOT / "artifacts", ROOT / "external"]


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"missing": True, "path": rel(path)}


def graph6_like(line: str) -> bool:
    s = line.strip()
    if not s or s.startswith("#") or len(s) < 2:
        return False
    return all(63 <= ord(ch) <= 126 for ch in s)


def scan_candidate_files() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for root in SEARCH_ROOTS:
        if not root.exists():
            rows.append({"root": rel(root), "exists": False, "candidate_file_count": 0, "max_graph6_like_line_count": 0})
            continue
        candidates = [p for p in root.rglob("*") if p.is_file() and p.suffix.lower() in CANDIDATE_SUFFIXES]
        max_count = 0
        max_path = None
        exact_hits = []
        for path in candidates:
            try:
                text = path.read_text(encoding="utf-8", errors="ignore")
            except OSError:
                continue
            count = sum(1 for line in text.splitlines() if graph6_like(line))
            if count > max_count:
                max_count = count
                max_path = path
            if count == EXPECTED_GIRTH4_CLASS_COUNT:
                exact_hits.append(rel(path))
        rows.append({
            "root": rel(root),
            "exists": True,
            "candidate_file_count": len(candidates),
            "max_graph6_like_line_count": max_count,
            "max_graph6_like_line_file": rel(max_path) if max_path else None,
            "exact_16828_line_candidate_files": exact_hits,
        })
    return rows


def build_manifest(p2796: dict[str, Any], p2798: dict[str, Any], p2799: dict[str, Any]) -> dict[str, Any]:
    scan_rows = scan_candidate_files()
    exact_hits = [hit for row in scan_rows for hit in row.get("exact_16828_line_candidate_files", [])]
    local_rows = p2799.get("local_girth4_spectral_orbit_witness", {}).get("table_row_count")
    target = p2798.get("external_girth4_subtarget_witness", {}).get("external_connected_16node_4regular_girth4_class_count")
    return {
        "expected_girth4_shortcode_class_count": EXPECTED_GIRTH4_CLASS_COUNT,
        "p2798_external_target_count": target,
        "p2799_local_girth4_table_row_count": local_rows,
        "p2796_full_generator_required_artifact_present": p2796.get("full_generator_artifact_absence_witness", {}).get("required_artifact_present"),
        "searched_roots": [rel(root) for root in SEARCH_ROOTS],
        "candidate_suffixes": sorted(CANDIDATE_SUFFIXES),
        "candidate_scan_rows": scan_rows,
        "exact_16828_line_candidate_file_count": len(exact_hits),
        "exact_16828_line_candidate_files": exact_hits,
        "required_shortcode_artifact_present": len(exact_hits) > 0,
        "subtarget_gap_if_no_import": EXPECTED_GIRTH4_CLASS_COUNT - (local_rows or 0),
        "finite_certificate_statement": "Repository-local scan finds no candidate graph-list/shortcode file with 16,828 graph6-like rows; P2799 remains a six-row local table, not a girth>=4 subtarget import.",
    }


def acceptance_matrix(manifest: dict[str, Any], p2796: dict[str, Any], p2798: dict[str, Any], p2799: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "p2796_absence_gate_present": p2796.get("status") == "P2796_FULL_GENERATOR_ARTIFACT_ABSENCE_PROVENANCE_GATE_NO_CLOSURE",
        "p2798_girth4_target_present": p2798.get("status") == "P2798_EXTERNAL_GIRTH4_SUBTARGET_LOCAL_GIRTH_FILTER_CERTIFICATE_NO_CLOSURE",
        "p2799_local_table_present": p2799.get("status") == "P2799_LOCAL_GIRTH4_SPECTRAL_ORBIT_TABLE_CERTIFICATE_NO_CLOSURE",
        "external_girth4_target_is_16828": manifest["p2798_external_target_count"] == EXPECTED_GIRTH4_CLASS_COUNT,
        "local_table_is_six_rows": manifest["p2799_local_girth4_table_row_count"] == 6,
        "no_16828_row_shortcode_candidate_found": manifest["exact_16828_line_candidate_file_count"] == 0,
        "required_shortcode_artifact_present": manifest["required_shortcode_artifact_present"],
        "strict_nadsoliton_spectral_source_law_exported": False,
        "kernel_or_ltotal_variational_coupling_exported": False,
    }
    return {
        "facts": facts,
        "accepted_as_girth4_shortcode_import_absence_manifest": all(facts[key] for key in [
            "p2796_absence_gate_present",
            "p2798_girth4_target_present",
            "p2799_local_table_present",
            "external_girth4_target_is_16828",
            "local_table_is_six_rows",
            "no_16828_row_shortcode_candidate_found",
        ]),
        "accepted_as_girth4_shortcode_graph_list_import": False,
        "accepted_as_full_16node_canonical_generator_certificate": False,
        "accepted_as_strict_spectral_source_law": False,
        "accepted_as_canonical_nadsoliton_geometry_source": False,
        "accepted_as_ltotal_or_toe_promotion": False,
        "missing_criteria": [key for key, value in facts.items() if not value],
        "blocker": "P2800 proves only repository-local absence/readiness for the requested 16,828-class girth>=4 shortcode import; it does not import the external graph list or export a strict spectral source law.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    m = payload["girth4_shortcode_import_absence_manifest"]
    lines = [
        "# P2800/S1750 girth>=4 shortcode import absence manifest",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite repository-local scan",
        f"- expected_girth4_shortcode_class_count={m['expected_girth4_shortcode_class_count']}",
        f"- p2799_local_girth4_table_row_count={m['p2799_local_girth4_table_row_count']}",
        f"- exact_16828_line_candidate_file_count={m['exact_16828_line_candidate_file_count']}",
        f"- required_shortcode_artifact_present={m['required_shortcode_artifact_present']}",
        f"- subtarget_gap_if_no_import={m['subtarget_gap_if_no_import']}",
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
    p2796 = read_json(P2796)
    p2798 = read_json(P2798)
    p2799 = read_json(P2799)
    manifest = build_manifest(p2796, p2798, p2799)
    acceptance = acceptance_matrix(manifest, p2796, p2798, p2799)
    payload = {
        "status": "P2800_GIRTH4_SHORTCODE_IMPORT_ABSENCE_MANIFEST_NO_CLOSURE",
        "input_hashes": {"P2796": sha(P2796), "P2798": sha(P2798), "P2799": sha(P2799)},
        "input_statuses": {"P2796": p2796.get("status"), "P2798": p2798.get("status"), "P2799": p2799.get("status")},
        "audited_question": "Is the actual 16,828-class girth>=4 shortcode/graph-list artifact already present locally so that the next full subtarget audit can honestly run?",
        "girth4_shortcode_import_absence_manifest": manifest,
        "acceptance_matrix": acceptance,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": acceptance["blocker"],
            "next_honest_step": "Use P2800 only as the repository-local absence/readiness manifest.  The next proof-grade move is to add or fetch the actual 16,828-class girth>=4 shortcode/graph-list artifact with a stable source URL, retrieval date, byte size, SHA256, row-count validation, and graph decoding smoke test; only after that should the exact quotient/charpoly/complement/orbit audit run.  If that artifact cannot be supplied, pivot to the full 8,037,418-class generator/toolchain import or export a strict spectral source law fixing the admissible class and K/L_total coupling.  Otherwise preserve the P2697-P2800 no-canonical-geometry/no-closure certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2800/S1750 girth4 shortcode import absence manifest", "## P2800/S1750 girth>=4 shortcode import absence manifest\n\n`P2800/S1750` audits whether the actual `16,828`-class girth>=4 shortcode/graph-list artifact is already present locally after the P2799 six-row table.  The finite repository scan finds no candidate file with `16,828` graph6-like rows, so the required subtarget import remains absent and the gap remains `16,822` beyond the six local rows.  This is an import-readiness/absence certificate only; it does not import the shortcode list, does not complete the girth>=4 quotient, and does not export a strict spectral source law or `K`/`L_total` coupling.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2800/S1750 girth4 shortcode import Ltotal guard", "## P2800/S1750 girth>=4 shortcode import Ltotal guard\n\n`P2800/S1750` adds no variational source term.  A repository-local absence/readiness manifest for the missing 16,828-class graph-list import cannot source a spectral action or license any `K`/`L_total` promotion.\n")
    append_once(AGENTS, "Current girth4 shortcode import absence guardrail (P2800/S1750, 2026-06-16)", "## Current girth>=4 shortcode import absence guardrail (P2800/S1750, 2026-06-16)\n\n- P2800 scans the repository for a local `16,828`-row girth>=4 shortcode/graph-list artifact after P2799 and finds no such artifact; the current evidence remains a six-row local girth>=4 table.\n- This is an absence/readiness manifest only; it is not the linked 16,828-class shortcode graph-list import, not a completed girth>=4 quotient, and not a strict spectral source law.\n- Do not promote P2800 to canonical geometry, strict spectral source law, selector closure, kernel full-expression, role-bearing `L_total`, bridge closure, role transfer, or ToE closure.  The next admissible move is adding/fetching the actual graph-list artifact with hash and row-count provenance, importing the full generator/toolchain, or exporting a strict spectral source law.\n")
    return payload


if __name__ == "__main__":
    main()
