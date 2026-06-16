#!/usr/bin/env python3
"""P2808/S1758: pynauty canonical digest manifest.

P2807 became executable once `pynauty` was available and produced a compact
canonical-certificate audit.  P2808 is the review-size-safe follow-up: generate
an ignored row-level canonical certificate digest CSV locally, then commit only a
small manifest containing the CSV hash, uniqueness counts, sample rows, and
automorphism-size histogram.  This is still graph-isomorphism provenance, not a
strict spectral source law or K/L_total coupling theorem.
"""
from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import EXPECTED_GRAPH_COUNT, SCD, decode_scd_bytes, sha
from p2807_s1757_pynauty_canonical_label_toolchain_gate import pynauty_probe, sha256_bytes, to_pynauty_graph

GEN = ROOT / "generated"
P2807 = GEN / "p2807_s1757_pynauty_canonical_label_toolchain_gate.json"
CSV_PATH = GEN / "p2808_s1758_pynauty_canonical_digest_rows.csv"
OUT = GEN / "p2808_s1758_pynauty_canonical_digest_manifest.json"
MD = GEN / "p2808_s1758_pynauty_canonical_digest_manifest.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def compute_rows() -> dict[str, Any]:
    probe = pynauty_probe()
    if not probe["available"] or probe.get("missing_symbols"):
        return {"probe": probe, "rows_written": False, "blocker": "pynauty unavailable or missing required symbols"}
    import pynauty  # type: ignore

    graphs, _ = decode_scd_bytes(SCD.read_bytes())
    cert_counter: Counter[str] = Counter()
    aut_counter: Counter[str] = Counter()
    sample_rows: list[dict[str, Any]] = []
    ordered_cert_hashes: list[str] = []
    with CSV_PATH.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["graph_index_1_based", "pynauty_certificate_sha256", "automorphism_group_size"])
        writer.writeheader()
        for index, adj in enumerate(graphs, start=1):
            graph = to_pynauty_graph(pynauty, adj)
            cert_hash = sha256_bytes(pynauty.certificate(graph))
            autgrp = pynauty.autgrp(graph)
            aut_size = str(autgrp[1]) if isinstance(autgrp, tuple) and len(autgrp) > 1 else "unknown"
            row = {"graph_index_1_based": index, "pynauty_certificate_sha256": cert_hash, "automorphism_group_size": aut_size}
            writer.writerow(row)
            cert_counter[cert_hash] += 1
            aut_counter[aut_size] += 1
            ordered_cert_hashes.append(cert_hash)
            if len(sample_rows) < 10:
                sample_rows.append(row)
    duplicate_classes = sum(1 for count in cert_counter.values() if count > 1)
    return {
        "probe": probe,
        "rows_written": True,
        "row_level_csv_path": rel(CSV_PATH),
        "row_level_csv_sha256": sha(CSV_PATH),
        "ordered_certificate_hash_stream_sha256": sha256_bytes("\n".join(ordered_cert_hashes).encode("utf-8")),
        "decoded_graph_count": len(graphs),
        "canonical_certificate_hash_class_count": len(cert_counter),
        "duplicate_certificate_class_count": duplicate_classes,
        "canonical_certificate_max_class_size": max(cert_counter.values()) if cert_counter else 0,
        "automorphism_group_size_histogram": dict(sorted(aut_counter.items(), key=lambda item: item[0])),
        "sample_rows": sample_rows,
    }


def acceptance_matrix(manifest: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "p2807_compact_audit_accepted": manifest["p2807_status"] == "P2807_PYNAUTY_COMPACT_CANONICAL_CERTIFICATE_AUDIT_NO_ROW_CSV_NO_SOURCE_LAW_NO_CLOSURE",
        "row_level_digest_csv_generated_but_ignored": manifest["canonical_digest_manifest"].get("rows_written") is True,
        "manifest_covers_16828_graphs": manifest["canonical_digest_manifest"].get("decoded_graph_count") == EXPECTED_GRAPH_COUNT,
        "canonical_certificate_hashes_unique": manifest["canonical_digest_manifest"].get("canonical_certificate_hash_class_count") == EXPECTED_GRAPH_COUNT,
        "strict_spectral_source_law_exported": False,
        "kernel_or_ltotal_variational_coupling_exported": False,
    }
    return {
        "facts": facts,
        "accepted_as_review_safe_canonical_digest_manifest": all(facts[key] for key in [
            "p2807_compact_audit_accepted",
            "row_level_digest_csv_generated_but_ignored",
            "manifest_covers_16828_graphs",
            "canonical_certificate_hashes_unique",
        ]),
        "accepted_as_strict_spectral_source_law": False,
        "accepted_as_ltotal_or_toe_promotion": False,
        "missing_criteria": [key for key, value in facts.items() if not value],
    }


def write_markdown(payload: dict[str, Any]) -> None:
    manifest = payload["canonical_digest_manifest"]
    lines = [
        "# P2808/S1758 pynauty canonical digest manifest",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Counts",
        f"- decoded_graph_count={manifest.get('decoded_graph_count')}",
        f"- canonical_certificate_hash_class_count={manifest.get('canonical_certificate_hash_class_count')}",
        f"- duplicate_certificate_class_count={manifest.get('duplicate_certificate_class_count')}",
        f"- canonical_certificate_max_class_size={manifest.get('canonical_certificate_max_class_size')}",
        f"- row_level_csv_sha256={manifest.get('row_level_csv_sha256')}",
        f"- ordered_certificate_hash_stream_sha256={manifest.get('ordered_certificate_hash_stream_sha256')}",
        "",
        "## Boundary",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    p2807 = read_json(P2807)
    manifest = {
        "status": "P2808_PYNAUTY_CANONICAL_DIGEST_MANIFEST_NO_SOURCE_LAW_NO_CLOSURE",
        "p2807_status": p2807.get("status"),
        "input_hashes": {"P2807": sha(P2807), "16_4_4.scd": sha(SCD)},
        "canonical_digest_manifest": compute_rows(),
        "decision": {
            "negative_export_flags": {
                "strict_spectral_source_law_exported": False,
                "kernel_geometry_closure_exported": False,
                "role_bearing_ltotal_promoted": False,
                "bridge_closure_exported": False,
                "role_transfer_started": False,
                "selector_closure_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2808 exports a compact digest over local pynauty canonical certificates and an ignored row-level CSV hash.  It proves duplicate-free canonical-certificate provenance for the validated 16,828 graph list, but it is not a strict spectral source law, not a variational coupling theorem, and not K/L_total or ToE closure.",
            "next_honest_step": "Run a narrow source-law/coupling admissibility audit: test whether the now-canonical 16,828-class girth>=4 graph quotient is actually connected to any exported strict spectral action/source functional.  If no such coupling theorem exists, emit a bounded no-source-law certificate rather than promoting canonical graph provenance to K/L_total or ToE.",
        },
    }
    manifest["acceptance_matrix"] = acceptance_matrix(manifest)
    OUT.write_text(json.dumps(manifest, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(manifest)
    append_once(STRICT_EQUATION_SHEET, "P2808/S1758 pynauty canonical digest manifest", "## P2808/S1758 pynauty canonical digest manifest\n\n`P2808/S1758` uses the now-available `pynauty` stack to generate a review-safe canonical-certificate digest manifest for all `16,828` validated Meringer graphs.  The row-level CSV is reproducible and ignored; the committed manifest records its hash, the ordered certificate-hash stream hash, uniqueness counts, and automorphism-size histogram.  This is graph-isomorphism provenance only, not a strict spectral source law or `K`/`L_total` coupling theorem.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2808/S1758 canonical digest Ltotal guard", "## P2808/S1758 canonical digest Ltotal guard\n\n`P2808/S1758` adds no variational term and no Euler-Lagrange source.  Canonical `pynauty` certificate uniqueness over the Meringer girth>=4 list strengthens graph provenance but cannot by itself promote `K`, role-bearing `L_total`, bridge closure, role transfer, selector closure, or ToE closure.\n")
    append_once(AGENTS, "Current pynauty canonical digest guardrail (P2808/S1758, 2026-06-16)", "## Current pynauty canonical digest guardrail (P2808/S1758, 2026-06-16)\n\n- P2808 runs after `pynauty` is available and records a review-safe canonical-certificate digest for all `16,828` validated Meringer girth>=4 graphs: unique certificate hashes, row-level CSV hash, ordered hash-stream digest, and automorphism-size histogram.\n- Treat P2808 as canonical graph-list provenance only.  It does not export a strict spectral source law, a coupling theorem from the graph quotient to `K`/`L_total`, role transfer, bridge closure, selector closure, or ToE closure.\n- The next admissible move is a narrow source-law/coupling admissibility audit for the canonical graph quotient; if no exported strict coupling exists, preserve a bounded no-source-law certificate.\n")
    return manifest


if __name__ == "__main__":
    main()
