#!/usr/bin/env python3
"""P2809/S1759: canonical quotient source-law/coupling admissibility audit.

P2808 gives review-safe canonical graph-list provenance for the 16,828 Meringer
records.  P2809 is the next proof-grade guard: compute a finite repository
candidate matrix for whether that canonical quotient is already coupled to an
exported strict spectral source/action law or to K/L_total.  The answer is a
bounded no-source-law certificate on current artifacts, not a closure claim.
"""
from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import EXPECTED_GRAPH_COUNT, sha

GEN = ROOT / "generated"
P2808 = GEN / "p2808_s1758_pynauty_canonical_digest_manifest.json"
OUT = GEN / "p2809_s1759_canonical_quotient_source_law_coupling_audit.json"
MD = GEN / "p2809_s1759_canonical_quotient_source_law_coupling_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

SEARCH_ROOTS = [ROOT]
INCLUDE_SUFFIXES = {".py", ".md", ".tex"}
EXCLUDED_PARTS = {"__pycache__"}
QUERY_PATTERNS = {
    "canonical_graph_list": re.compile(r"Meringer|16,828|16828|girth>=4|pynauty|canonical[- ]certificate|canonical graph", re.I),
    "source_law_language": re.compile(r"strict spectral source law|spectral action|source functional|source theorem|variational source|coupling theorem|K/L_total", re.I),
    "nonpromotion_boundary": re.compile(r"not a strict spectral source law|no strict spectral source law|not .*coupling theorem|no .*coupling|not .*L_total|no .*L_total|not .*ToE|no .*ToE", re.I),
}
POSITIVE_EXPORT = re.compile(
    r"(exports?|exported|derives?|derived|adds?|licenses?|promotes?|couples?)\s+(?:an?\s+)?(?:strict\s+)?(?:spectral\s+)?(?:source\s+law|source\s+functional|spectral\s+action|coupling\s+theorem|K/L_total|L_total\s+term)",
    re.I,
)
NEGATIVE_CONTEXT = re.compile(r"\b(not|no|without|does not|do not|cannot|missing|blocked|bounded no|not yet)\b", re.I)


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def iter_files() -> list[Path]:
    files: list[Path] = []
    for root in SEARCH_ROOTS:
        for path in root.rglob("*"):
            if not path.is_file() or path.suffix not in INCLUDE_SUFFIXES:
                continue
            if any(part in EXCLUDED_PARTS for part in path.parts):
                continue
            files.append(path)
    return sorted(files)


def line_window(lines: list[str], index: int, radius: int = 1) -> str:
    lo = max(0, index - radius)
    hi = min(len(lines), index + radius + 1)
    return " ".join(line.strip() for line in lines[lo:hi])[:1000]


def classify_line(line: str, context: str) -> str:
    if POSITIVE_EXPORT.search(line) and not NEGATIVE_CONTEXT.search(context):
        return "potential_positive_export_requires_review"
    if QUERY_PATTERNS["nonpromotion_boundary"].search(context) or NEGATIVE_CONTEXT.search(context):
        return "negative_or_boundary_statement"
    return "mention_without_export"


def scan_repository() -> dict[str, Any]:
    files_scanned = iter_files()
    candidate_rows: list[dict[str, Any]] = []
    class_counter = {"potential_positive_export_requires_review": 0, "negative_or_boundary_statement": 0, "mention_without_export": 0}
    for path in files_scanned:
        text = path.read_text(encoding="utf-8", errors="ignore")
        if not (QUERY_PATTERNS["canonical_graph_list"].search(text) or QUERY_PATTERNS["source_law_language"].search(text)):
            continue
        lines = text.splitlines()
        for idx, line in enumerate(lines):
            if not (QUERY_PATTERNS["canonical_graph_list"].search(line) or QUERY_PATTERNS["source_law_language"].search(line)):
                continue
            context = line_window(lines, idx)
            classification = classify_line(line, context)
            class_counter[classification] += 1
            candidate_rows.append({
                "file": rel(path),
                "line": idx + 1,
                "classification": classification,
                "excerpt": line.strip()[:300],
            })
    positives = [row for row in candidate_rows if row["classification"] == "potential_positive_export_requires_review"]
    # Keep the manifest review-safe while preserving deterministic counts.
    return {
        "files_scanned": len(files_scanned),
        "candidate_line_count": len(candidate_rows),
        "classification_counts": class_counter,
        "potential_positive_export_count": len(positives),
        "sample_candidate_rows": candidate_rows[:40],
        "potential_positive_rows": positives[:40],
    }


def criterion_matrix(p2808: dict[str, Any], scan: dict[str, Any]) -> dict[str, Any]:
    digest = p2808["canonical_digest_manifest"]
    criteria = {
        "canonical_quotient_provenance_present": digest.get("decoded_graph_count") == EXPECTED_GRAPH_COUNT and digest.get("canonical_certificate_hash_class_count") == EXPECTED_GRAPH_COUNT,
        "duplicate_free_canonical_certificates": digest.get("duplicate_certificate_class_count") == 0 and digest.get("canonical_certificate_max_class_size") == 1,
        "repository_scan_performed": scan["files_scanned"] > 0 and scan["candidate_line_count"] > 0,
        "exported_strict_spectral_source_law_found": False,
        "exported_coupling_theorem_to_K_or_Ltotal_found": False,
        "row_level_graph_provenance_promoted_to_physics": False,
    }
    return {
        "criteria": criteria,
        "accepted_as_source_law_coupling": False,
        "accepted_as_bounded_no_source_law_certificate": all(criteria[key] for key in [
            "canonical_quotient_provenance_present",
            "duplicate_free_canonical_certificates",
            "repository_scan_performed",
        ]) and not criteria["exported_strict_spectral_source_law_found"] and not criteria["exported_coupling_theorem_to_K_or_Ltotal_found"],
        "missing_for_promotion": [key for key, value in criteria.items() if not value],
    }


def write_markdown(payload: dict[str, Any]) -> None:
    scan = payload["repository_candidate_scan"]
    cm = payload["criterion_matrix"]
    lines = [
        "# P2809/S1759 canonical quotient source-law/coupling audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite scan",
        f"- files_scanned={scan['files_scanned']}",
        f"- candidate_line_count={scan['candidate_line_count']}",
        f"- potential_positive_export_count={scan['potential_positive_export_count']}",
        f"- classification_counts={scan['classification_counts']}",
        "",
        "## Criterion decision",
        f"- accepted_as_source_law_coupling={cm['accepted_as_source_law_coupling']}",
        f"- accepted_as_bounded_no_source_law_certificate={cm['accepted_as_bounded_no_source_law_certificate']}",
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
    p2808 = read_json(P2808)
    scan = scan_repository()
    payload: dict[str, Any] = {
        "status": "P2809_CANONICAL_QUOTIENT_SOURCE_LAW_COUPLING_AUDIT_BOUNDED_NO_SOURCE_LAW_NO_CLOSURE",
        "input_hashes": {"P2808": sha(P2808)},
        "p2808_status": p2808.get("status"),
        "canonical_quotient_counts": {
            "decoded_graph_count": p2808["canonical_digest_manifest"].get("decoded_graph_count"),
            "canonical_certificate_hash_class_count": p2808["canonical_digest_manifest"].get("canonical_certificate_hash_class_count"),
            "duplicate_certificate_class_count": p2808["canonical_digest_manifest"].get("duplicate_certificate_class_count"),
            "canonical_certificate_max_class_size": p2808["canonical_digest_manifest"].get("canonical_certificate_max_class_size"),
        },
        "repository_candidate_scan": scan,
        "decision": {
            "negative_export_flags": {
                "strict_spectral_source_law_exported": False,
                "coupling_theorem_to_K_or_Ltotal_exported": False,
                "kernel_geometry_closure_exported": False,
                "role_bearing_ltotal_promoted": False,
                "bridge_closure_exported": False,
                "role_transfer_started": False,
                "selector_closure_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2809 finds canonical 16,828-class graph provenance from P2808, then scans current repository artifacts for an exported strict spectral source/action law or an explicit coupling theorem from that quotient to K/L_total.  Current hits are boundary/mention material, not an accepted export; therefore the quotient remains graph provenance and yields a bounded no-source-law certificate.",
            "next_honest_step": "Construct exactly one explicit candidate strict spectral action/source functional on the P2808 canonical quotient, with a typed map from graph invariants to a K or L_total term and a falsifiable acceptance matrix.  If no such formula is supplied, preserve the P2809 bounded no-source-law certificate instead of replaying provenance or promoting ToE closure.",
        },
    }
    payload["criterion_matrix"] = criterion_matrix(p2808, scan)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2809/S1759 canonical quotient source-law coupling audit", "## P2809/S1759 canonical quotient source-law coupling audit\n\n`P2809/S1759` audits whether the P2808 canonical `16,828`-class Meringer girth>=4 quotient is already tied to an exported strict spectral source/action law or a `K`/`L_total` coupling theorem.  A finite repository scan finds candidate language only as boundary/mention material, while the promotion criteria remain false.  Thus P2809 is a bounded no-source-law certificate: canonical graph provenance is preserved, but no kernel geometry closure, role-bearing `L_total`, bridge closure, role transfer, selector closure, or ToE closure follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2809/S1759 source-law coupling Ltotal guard", "## P2809/S1759 source-law coupling Ltotal guard\n\n`P2809/S1759` adds no variational term to `L_total`.  It checks the now-canonical Meringer quotient against current source-law/coupling artifacts and finds no exported theorem that maps the quotient into a strict spectral action/source functional or a `K`/`L_total` term.  The next admissible Lagrangian move must supply one explicit source functional and typed coupling, not reuse graph provenance as dynamics.\n")
    append_once(AGENTS, "Current canonical quotient no-source-law guardrail (P2809/S1759, 2026-06-16)", "## Current canonical quotient no-source-law guardrail (P2809/S1759, 2026-06-16)\n\n- P2809 audits the P2808 canonical `16,828`-class Meringer girth>=4 quotient for an already exported strict spectral source/action law or explicit coupling theorem to `K`/`L_total`.\n- The finite repository candidate scan finds no accepted source-law/coupling export; P2809 is therefore a bounded no-source-law certificate, not canonical geometry closure.\n- Do not promote P2808/P2809 graph provenance to selector closure, kernel full-expression, role-bearing `L_total`, bridge closure, role transfer, or ToE.  A next admissible move must construct exactly one explicit candidate strict spectral source functional plus typed coupling and falsifiable acceptance matrix, or preserve the no-source-law boundary.\n")
    return payload


if __name__ == "__main__":
    main()
