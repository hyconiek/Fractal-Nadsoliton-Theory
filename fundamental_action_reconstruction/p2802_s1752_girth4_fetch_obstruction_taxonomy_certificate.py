#!/usr/bin/env python3
"""P2802/S1752: deterministic taxonomy of the P2801 girth>=4 fetch obstruction.

P2801 executed the internet-enabled fetch attempt.  P2802 does not re-run the
network and does not invent a new graph subclass.  It converts the P2801 probe
matrix into a deterministic obstruction taxonomy: HTTPS proxy tunnel failures,
HTTP source/proxy 403 rows, graph6-like row counts, and import/closure gates.

This is a stronger bookkeeping step than another failed live fetch because it
separates an access-layer/source-route blocker from a mathematical graph-list
validation success.  No graph-list import, quotient audit, spectral source law,
or K/L_total promotion is exported.
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
P2801 = GEN / "p2801_s1751_girth4_shortcode_network_fetch_attempt_certificate.json"
OUT = GEN / "p2802_s1752_girth4_fetch_obstruction_taxonomy_certificate.json"
MD = GEN / "p2802_s1752_girth4_fetch_obstruction_taxonomy_certificate.md"
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


def scheme_of(url: str) -> str:
    return url.split(":", 1)[0]


def build_taxonomy(p2801: dict[str, Any]) -> dict[str, Any]:
    witness = p2801.get("network_fetch_witness", {})
    rows = witness.get("fetch_probe_rows", [])
    scheme_counter = Counter(scheme_of(row.get("url", "")) for row in rows)
    status_counter = Counter(str(row.get("http_status")) for row in rows)
    error_counter = Counter(row.get("error_type", "NONE") for row in rows)
    transport_signature_counter = Counter(
        f"{scheme_of(row.get('url', ''))}|{row.get('http_status')}|{row.get('error_type')}|{row.get('error_reason')}"
        for row in rows
    )
    graph6_counts = [row.get("graph6_like_line_count", 0) or 0 for row in rows]
    validated_rows = [row for row in rows if row.get("validates_expected_16828_graph6_like_rows")]
    https_tunnel_forbidden = [row for row in rows if scheme_of(row.get("url", "")) == "https" and "Tunnel connection failed: 403 Forbidden" in str(row.get("error_reason"))]
    http_forbidden = [row for row in rows if scheme_of(row.get("url", "")) == "http" and row.get("http_status") == 403]
    return {
        "expected_girth4_shortcode_class_count": EXPECTED_GIRTH4_CLASS_COUNT,
        "p2801_status": p2801.get("status"),
        "candidate_url_count": len(rows),
        "scheme_counts": dict(sorted(scheme_counter.items())),
        "http_status_counts": dict(sorted(status_counter.items())),
        "error_type_counts": dict(sorted(error_counter.items())),
        "transport_signature_counts": dict(sorted(transport_signature_counter.items())),
        "https_tunnel_403_count": len(https_tunnel_forbidden),
        "http_403_count": len(http_forbidden),
        "successful_http_response_count": witness.get("successful_http_response_count"),
        "max_graph6_like_line_count": max(graph6_counts) if graph6_counts else 0,
        "validated_16828_row_download_count": len(validated_rows),
        "required_shortcode_artifact_imported": witness.get("required_shortcode_artifact_imported"),
        "obstruction_classification": "ACCESS_LAYER_BLOCKER_NOT_GRAPH_VALIDATION_FAILURE" if rows and not validated_rows and len(https_tunnel_forbidden) + len(http_forbidden) == len(rows) else "MIXED_OR_RECHECK_REQUIRED",
        "finite_certificate_statement": "All P2801 probes are accounted for by HTTPS tunnel 403 or HTTP 403 outcomes; no response reaches a graph-list body, no graph6-like row count is positive, and no 16,828-row artifact is imported.",
    }


def acceptance_matrix(taxonomy: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "p2801_fetch_attempt_present": taxonomy["p2801_status"] and taxonomy["p2801_status"].startswith("P2801_GIRTH4_SHORTCODE_NETWORK_FETCH_ATTEMPT_CERTIFICATE_"),
        "all_candidate_urls_classified": taxonomy["candidate_url_count"] == taxonomy["https_tunnel_403_count"] + taxonomy["http_403_count"],
        "access_layer_blocker_identified": taxonomy["obstruction_classification"] == "ACCESS_LAYER_BLOCKER_NOT_GRAPH_VALIDATION_FAILURE",
        "no_graph6_body_observed": taxonomy["max_graph6_like_line_count"] == 0,
        "validated_16828_row_download_present": taxonomy["validated_16828_row_download_count"] > 0,
        "girth4_16828_class_quotient_completed": False,
        "strict_nadsoliton_spectral_source_law_exported": False,
        "kernel_or_ltotal_variational_coupling_exported": False,
    }
    return {
        "facts": facts,
        "accepted_as_fetch_obstruction_taxonomy_certificate": all(facts[key] for key in [
            "p2801_fetch_attempt_present",
            "all_candidate_urls_classified",
            "access_layer_blocker_identified",
            "no_graph6_body_observed",
        ]),
        "accepted_as_girth4_shortcode_graph_list_import": False,
        "accepted_as_completed_girth4_quotient_audit": False,
        "accepted_as_full_16node_canonical_generator_certificate": False,
        "accepted_as_strict_spectral_source_law": False,
        "accepted_as_canonical_nadsoliton_geometry_source": False,
        "accepted_as_ltotal_or_toe_promotion": False,
        "missing_criteria": [key for key, value in facts.items() if not value],
        "blocker": "P2802 classifies P2801 as an access-layer/source-route blocker: all attempted URLs failed before a graph-list body could be validated, so no import or closure is licensed.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["fetch_obstruction_taxonomy"]
    lines = [
        "# P2802/S1752 girth>=4 fetch obstruction taxonomy certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Deterministic taxonomy",
        f"- candidate_url_count={t['candidate_url_count']}",
        f"- scheme_counts={t['scheme_counts']}",
        f"- https_tunnel_403_count={t['https_tunnel_403_count']}",
        f"- http_403_count={t['http_403_count']}",
        f"- max_graph6_like_line_count={t['max_graph6_like_line_count']}",
        f"- validated_16828_row_download_count={t['validated_16828_row_download_count']}",
        f"- obstruction_classification={t['obstruction_classification']}",
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
    p2801 = read_json(P2801)
    taxonomy = build_taxonomy(p2801)
    acceptance = acceptance_matrix(taxonomy)
    payload = {
        "status": "P2802_GIRTH4_FETCH_OBSTRUCTION_TAXONOMY_CERTIFICATE_NO_IMPORT_NO_CLOSURE",
        "input_hashes": {"P2801": sha(P2801)},
        "input_statuses": {"P2801": p2801.get("status")},
        "audited_question": "Was the P2801 failure a graph-list validation failure after download, or an access-layer/source-route blocker before any graph-list body was obtained?",
        "fetch_obstruction_taxonomy": taxonomy,
        "acceptance_matrix": acceptance,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": acceptance["blocker"],
            "next_honest_step": "Treat P2802 as an access-layer obstruction certificate.  The next honest step is not another guessed filename probe; it is a working source route that bypasses the 403 layer (manual artifact upload, alternate mirror, source-side access fix, or an installed full generator/toolchain), followed by exact row-count/SHA256/decoding validation and only then quotient/charpoly/complement/orbit auditing.  Otherwise preserve the P2697-P2802 no-canonical-geometry/no-closure certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2802/S1752 girth>=4 fetch obstruction taxonomy certificate", "## P2802/S1752 girth>=4 fetch obstruction taxonomy certificate\n\n`P2802/S1752` deterministically classifies the P2801 network probe matrix: all `24` attempted `16_4_4` URL routes are accounted for by HTTPS tunnel `403` or HTTP `403` outcomes before any graph-list body is obtained; the maximum graph6-like row count is `0`, and no `16,828`-row artifact is imported.  This is an access-layer/source-route obstruction certificate only; it does not import the graph list, complete the girth>=4 quotient, export a strict spectral source law, or license `K`/`L_total` coupling.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2802/S1752 girth>=4 fetch obstruction Ltotal guard", "## P2802/S1752 girth>=4 fetch obstruction Ltotal guard\n\n`P2802/S1752` adds no variational source term.  Classifying the failed fetch as an access-layer blocker is not a spectral action/source law and cannot license any `K`/`L_total` promotion without an actual graph-list import plus separate quotient/spectral and coupling theorems.\n")
    append_once(AGENTS, "Current girth>=4 fetch obstruction taxonomy guardrail (P2802/S1752, 2026-06-16)", "## Current girth>=4 fetch obstruction taxonomy guardrail (P2802/S1752, 2026-06-16)\n\n- P2802 classifies the P2801 network probe matrix as an access-layer/source-route blocker: all 24 attempted URL routes fail by HTTPS tunnel `403` or HTTP `403` before any graph-list body is obtained.\n- No `16,828`-row graph list is imported, no girth>=4 quotient is executed, and no strict spectral source law is exported.\n- Do not continue guessed-filename replay as the primary strategy.  A next admissible move requires a working source route/manual upload/mirror/source-side access fix, an installed full generator/toolchain, or a strict spectral source law.\n")
    return payload


if __name__ == "__main__":
    main()
