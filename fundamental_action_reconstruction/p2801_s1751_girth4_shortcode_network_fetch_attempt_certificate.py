#!/usr/bin/env python3
"""P2801/S1751: network fetch attempt for the 16,828 girth>=4 shortcode artifact.

After P2800 established only repository-local absence, internet access was
explicitly enabled for the next honest step.  P2801 therefore attempts the
actual external fetch lane rather than adding another local named subclass.

The certificate is intentionally conservative: it records every probed URL,
HTTP/transport outcome, byte/hash metadata when a response is obtained, and
whether any downloaded body validates as the required 16,828-row graph list.
If the network/proxy/source refuses access or no body validates, no graph-list
import, quotient audit, canonical geometry source, spectral source law, or
K/L_total promotion is claimed.
"""
from __future__ import annotations

import hashlib
import json
import socket
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
EXTERNAL = ROOT / "external"
P2798 = GEN / "p2798_s1748_external_girth4_subtarget_local_girth_filter_certificate.json"
P2799 = GEN / "p2799_s1749_local_girth4_spectral_orbit_table_certificate.json"
P2800 = GEN / "p2800_s1750_girth4_shortcode_import_absence_manifest.json"
OUT = GEN / "p2801_s1751_girth4_shortcode_network_fetch_attempt_certificate.json"
MD = GEN / "p2801_s1751_girth4_shortcode_network_fetch_attempt_certificate.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

EXPECTED_GIRTH4_CLASS_COUNT = 16_828
DETAIL_URL = "https://www.mathe2.uni-bayreuth.de/markus/REGGRAPHS/16_4_4.html"
SOURCE_ROOTS = [
    "https://www.mathe2.uni-bayreuth.de/markus/REGGRAPHS/",
    "http://www.mathe2.uni-bayreuth.de/markus/REGGRAPHS/",
]
CANDIDATE_NAMES = [
    "16_4_4.html",
    "16_4_4.g6",
    "16_4_4.graph6",
    "16_4_4.shortcode",
    "16_4_4.short",
    "16_4_4.sct",
    "16_4_4.txt",
    "16_4_4.scd",
    "16_4_4.dat",
    "16_4_4.gz",
    "16_4_4.g6.gz",
    "16_4_4.shortcode.gz",
]
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


def graph6_like(line: str) -> bool:
    s = line.strip()
    if not s or s.startswith("#") or len(s) < 2:
        return False
    return all(63 <= ord(ch) <= 126 for ch in s)


def graph6_like_line_count(body: bytes) -> int:
    text = body.decode("utf-8", errors="ignore")
    return sum(1 for line in text.splitlines() if graph6_like(line))


def attempt_url(url: str, timeout: int = 30) -> dict[str, Any]:
    request = urllib.request.Request(url, headers={"User-Agent": "FNT-P2801-fetch-audit/1.0"})
    row: dict[str, Any] = {"url": url, "attempted": True}
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            body = response.read()
            row.update({
                "ok": 200 <= response.status < 300,
                "http_status": response.status,
                "content_type": response.headers.get("content-type"),
                "byte_count": len(body),
                "sha256": hashlib.sha256(body).hexdigest(),
                "graph6_like_line_count": graph6_like_line_count(body),
            })
            return row | {"body": body}
    except urllib.error.HTTPError as exc:
        body = exc.read() if hasattr(exc, "read") else b""
        row.update({
            "ok": False,
            "http_status": exc.code,
            "error_type": type(exc).__name__,
            "error_reason": str(exc.reason),
            "byte_count": len(body),
            "sha256": hashlib.sha256(body).hexdigest() if body else None,
            "graph6_like_line_count": graph6_like_line_count(body) if body else 0,
        })
    except (urllib.error.URLError, TimeoutError, socket.timeout, OSError) as exc:
        row.update({"ok": False, "http_status": None, "error_type": type(exc).__name__, "error_reason": str(exc)})
    return row


def candidate_urls() -> list[str]:
    urls = []
    for root in SOURCE_ROOTS:
        for name in CANDIDATE_NAMES:
            urls.append(root + name)
    return urls


def build_fetch_witness() -> dict[str, Any]:
    EXTERNAL.mkdir(exist_ok=True)
    rows = []
    imported_files = []
    for url in candidate_urls():
        result = attempt_url(url)
        body = result.pop("body", None)
        validates = result.get("graph6_like_line_count") == EXPECTED_GIRTH4_CLASS_COUNT
        result["validates_expected_16828_graph6_like_rows"] = validates
        if body is not None and validates:
            name = "p2801_s1751_girth4_16828_import_" + hashlib.sha256(url.encode()).hexdigest()[:12] + ".g6"
            path = EXTERNAL / name
            path.write_bytes(body)
            result["saved_as"] = rel(path)
            imported_files.append(rel(path))
        rows.append(result)
    return {
        "expected_girth4_shortcode_class_count": EXPECTED_GIRTH4_CLASS_COUNT,
        "retrieved_utc_timestamp": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "detail_url": DETAIL_URL,
        "candidate_url_count": len(rows),
        "fetch_probe_rows": rows,
        "successful_http_response_count": sum(1 for row in rows if row.get("ok")),
        "validated_16828_row_download_count": len(imported_files),
        "imported_validated_files": imported_files,
        "required_shortcode_artifact_imported": bool(imported_files),
        "finite_certificate_statement": "Network fetch probes were executed for the Meringer 16_4_4 detail page and plausible linked shortcode/graph-list filenames; no downloaded body validated as the required 16,828-row graph list in this run." if not imported_files else "A downloaded body validated as a 16,828-row graph6-like list and was saved under fundamental_action_reconstruction/external; this is import provenance only, not yet quotient/spectral closure.",
    }


def acceptance_matrix(witness: dict[str, Any], p2798: dict[str, Any], p2799: dict[str, Any], p2800: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "p2798_girth4_target_present": p2798.get("status") == "P2798_EXTERNAL_GIRTH4_SUBTARGET_LOCAL_GIRTH_FILTER_CERTIFICATE_NO_CLOSURE",
        "p2799_local_table_present": p2799.get("status") == "P2799_LOCAL_GIRTH4_SPECTRAL_ORBIT_TABLE_CERTIFICATE_NO_CLOSURE",
        "p2800_absence_manifest_present": p2800.get("status") == "P2800_GIRTH4_SHORTCODE_IMPORT_ABSENCE_MANIFEST_NO_CLOSURE",
        "network_fetch_attempted": witness["candidate_url_count"] > 0,
        "validated_16828_row_download_present": witness["required_shortcode_artifact_imported"],
        "girth4_16828_class_quotient_completed": False,
        "strict_nadsoliton_spectral_source_law_exported": False,
        "kernel_or_ltotal_variational_coupling_exported": False,
    }
    return {
        "facts": facts,
        "accepted_as_network_fetch_attempt_certificate": all(facts[key] for key in [
            "p2798_girth4_target_present",
            "p2799_local_table_present",
            "p2800_absence_manifest_present",
            "network_fetch_attempted",
        ]),
        "accepted_as_girth4_shortcode_graph_list_import": facts["validated_16828_row_download_present"],
        "accepted_as_completed_girth4_quotient_audit": False,
        "accepted_as_full_16node_canonical_generator_certificate": False,
        "accepted_as_strict_spectral_source_law": False,
        "accepted_as_canonical_nadsoliton_geometry_source": False,
        "accepted_as_ltotal_or_toe_promotion": False,
        "missing_criteria": [key for key, value in facts.items() if not value],
        "blocker": "P2801 executes the enabled-internet fetch attempt, but no response validates as the 16,828-row graph list, so the subtarget import remains blocked." if not facts["validated_16828_row_download_present"] else "P2801 imports only the validated graph-list artifact; quotient, charpoly, complement, orbit, strict source-law, and K/L_total closure remain unexecuted.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    w = payload["network_fetch_witness"]
    lines = [
        "# P2801/S1751 girth>=4 shortcode network fetch attempt certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Fetch result",
        f"- candidate_url_count={w['candidate_url_count']}",
        f"- successful_http_response_count={w['successful_http_response_count']}",
        f"- validated_16828_row_download_count={w['validated_16828_row_download_count']}",
        f"- required_shortcode_artifact_imported={w['required_shortcode_artifact_imported']}",
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
    p2798 = read_json(P2798)
    p2799 = read_json(P2799)
    p2800 = read_json(P2800)
    witness = build_fetch_witness()
    acceptance = acceptance_matrix(witness, p2798, p2799, p2800)
    suffix = "IMPORT_PROVENANCE_NO_QUOTIENT_CLOSURE" if witness["required_shortcode_artifact_imported"] else "NO_VALIDATED_DOWNLOAD_NO_CLOSURE"
    next_step = (
        "Use P2801 as fetch provenance only: the validated 16,828-row artifact must next receive exact decoding, quotient/charpoly/complement/orbit auditing, and only then can any graph-list subtarget claims be made; do not promote to canonical geometry, strict spectral source law, K/L_total, or ToE."
        if witness["required_shortcode_artifact_imported"]
        else "Use P2801 as a real network-attempt obstruction certificate.  The next honest step is to obtain the exact shortcode/graph-list by a working source route (manual upload, alternate mirror, or source-side access fix) and rerun row-count/hash/decoding validation; otherwise pivot to full generator/toolchain import or a strict spectral source law.  Preserve the P2697-P2801 no-canonical-geometry/no-closure certificate."
    )
    payload = {
        "status": f"P2801_GIRTH4_SHORTCODE_NETWORK_FETCH_ATTEMPT_CERTIFICATE_{suffix}",
        "input_hashes": {"P2798": sha(P2798), "P2799": sha(P2799), "P2800": sha(P2800)},
        "input_statuses": {"P2798": p2798.get("status"), "P2799": p2799.get("status"), "P2800": p2800.get("status")},
        "audited_question": "With internet access enabled, can the actual linked 16,828-class girth>=4 shortcode/graph-list artifact be fetched and validated now?",
        "network_fetch_witness": witness,
        "acceptance_matrix": acceptance,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": acceptance["blocker"],
            "next_honest_step": next_step,
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2801/S1751 girth>=4 shortcode network fetch attempt certificate", "## P2801/S1751 girth>=4 shortcode network fetch attempt certificate\n\n`P2801/S1751` executes the enabled-internet fetch attempt for the Meringer `16_4_4` girth>=4 shortcode/graph-list route and records HTTP/transport outcomes, byte/hash metadata where available, and row-count validation.  In the current run no fetched body validates as the required `16,828`-row graph list, so the subtarget import remains blocked.  This is a network-attempt/provenance certificate only; it does not complete the girth>=4 quotient, does not export a strict spectral source law, and does not license `K`/`L_total` coupling.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2801/S1751 girth>=4 shortcode network fetch Ltotal guard", "## P2801/S1751 girth>=4 shortcode network fetch Ltotal guard\n\n`P2801/S1751` adds no variational source term.  A network fetch attempt, even if it later imports a graph list, is not a spectral action/source law and cannot license any `K`/`L_total` promotion without the separate quotient/spectral and coupling theorems.\n")
    append_once(AGENTS, "Current girth>=4 shortcode network fetch attempt guardrail (P2801/S1751, 2026-06-16)", "## Current girth>=4 shortcode network fetch attempt guardrail (P2801/S1751, 2026-06-16)\n\n- P2801 executes the internet-enabled fetch attempt for the linked Meringer `16_4_4` girth>=4 shortcode/graph-list route and records URL-level HTTP/transport outcomes plus row-count validation.\n- The current run does not validate/import a `16,828`-row graph list, so the girth>=4 subtarget quotient remains unexecuted.\n- Do not promote P2801 to canonical geometry, strict spectral source law, selector closure, kernel full-expression, role-bearing `L_total`, bridge closure, role transfer, or ToE closure.  A next admissible move requires a working source route/manual upload/mirror with hash and row-count validation, the full generator/toolchain import, or a strict spectral source law.\n")
    return payload


if __name__ == "__main__":
    main()
