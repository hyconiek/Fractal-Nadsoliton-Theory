#!/usr/bin/env python3
"""P2900/S1850: defect-placement source-law inventory no-go.

P2899 isolated the exact missing object: a strict law sourcing defect
placement/basepoint/polarity with coupling to the 9/5 variational density.  P2900
performs the next honest intake step: scan the current generated JSON artifacts
for that exact combined export rather than replaying torsor, scalar, relation, or
defect constructions.

This is a current-artifact inventory certificate.  It does not prove no future
law can exist; it proves that the required coupled source-law artifact is not
currently exported in the generated artifact corpus.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import sha
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN

P2899 = GEN / "p2899_s1849_post_defect_potential_readiness_matrix.json"
OUT = GEN / "p2900_s1850_defect_placement_source_law_inventory_no_go.json"
MD = GEN / "p2900_s1850_defect_placement_source_law_inventory_no_go.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

SOURCE_TERMS = (
    "strict_defect_placement_source_law_exported",
    "nonimported_single_defect_source_exported",
    "defect_placement_source_law_exported",
    "nonimported_basepoint_or_polarity_law_exported",
    "strict_phase_origin_source_artifact_exported",
)
COUPLING_TERMS = (
    "coupling_to_9_over_5_variational_density_exported",
    "exports_unique_coupling_to_9_over_5_carrier",
    "nonimported_9_over_5_variational_chain_rule_exported",
)
CLOSURE_TERMS = (
    "localized_action_density_exported",
    "ltotal_exported",
    "eom_closure_exported",
    "hamiltonian_closure_exported",
    "toe_closure_exported",
)
QUARANTINE_PATH_MARKERS = ("negative_export_flags", "no_", "not_", "blocked", "accepted_as_toe_potential_evidence")


def json_files() -> list[Path]:
    return sorted(path for path in GEN.glob("*.json") if path.is_file() and path != OUT)


def walk(value: Any, prefix: str = "") -> list[tuple[str, Any]]:
    rows: list[tuple[str, Any]] = []
    if isinstance(value, dict):
        for key, child in value.items():
            rows.extend(walk(child, f"{prefix}.{key}" if prefix else str(key)))
    elif isinstance(value, list):
        for index, child in enumerate(value):
            rows.extend(walk(child, f"{prefix}[{index}]"))
    else:
        rows.append((prefix, value))
    return rows


def classify_hit(path: str, value: Any) -> str:
    lower = path.lower()
    if any(marker in lower for marker in QUARANTINE_PATH_MARKERS):
        return "quarantined_boundary_or_negative_flag"
    if value is True:
        return "positive_like_true"
    return "non_true_occurrence"


def scan_generated_json() -> dict[str, Any]:
    hits = []
    source_positive = []
    coupling_positive = []
    closure_positive = []
    coupled_positive_artifacts: dict[str, set[str]] = {}
    for file_path in json_files():
        try:
            payload = read_json(file_path)
        except Exception:
            continue
        rel = str(file_path.relative_to(ROOT))
        file_positive_kinds: set[str] = set()
        for path, value in walk(payload):
            lower = path.lower()
            term_kind = None
            if any(term in lower for term in SOURCE_TERMS):
                term_kind = "source"
            elif any(term in lower for term in COUPLING_TERMS):
                term_kind = "coupling"
            elif any(term in lower for term in CLOSURE_TERMS):
                term_kind = "closure"
            if term_kind is None:
                continue
            classification = classify_hit(path, value)
            hit = {"file": rel, "path": path, "value": value, "term_kind": term_kind, "classification": classification}
            hits.append(hit)
            if classification == "positive_like_true":
                file_positive_kinds.add(term_kind)
                if term_kind == "source":
                    source_positive.append(hit)
                elif term_kind == "coupling":
                    coupling_positive.append(hit)
                elif term_kind == "closure":
                    closure_positive.append(hit)
        if file_positive_kinds:
            coupled_positive_artifacts[rel] = file_positive_kinds
    coupled_source_and_coupling = sorted(file for file, kinds in coupled_positive_artifacts.items() if {"source", "coupling"}.issubset(kinds))
    coupled_source_coupling_and_closure = sorted(file for file, kinds in coupled_positive_artifacts.items() if {"source", "coupling", "closure"}.issubset(kinds))
    return {
        "generated_json_file_count": len(json_files()),
        "term_occurrence_count": len(hits),
        "positive_source_hit_count": len(source_positive),
        "positive_coupling_hit_count": len(coupling_positive),
        "positive_closure_hit_count": len(closure_positive),
        "coupled_source_and_coupling_artifact_count": len(coupled_source_and_coupling),
        "coupled_source_coupling_and_closure_artifact_count": len(coupled_source_coupling_and_closure),
        "coupled_source_and_coupling_artifacts": coupled_source_and_coupling[:20],
        "coupled_source_coupling_and_closure_artifacts": coupled_source_coupling_and_closure[:20],
        "sample_hits": hits[:30],
        "positive_source_hits": source_positive[:20],
        "positive_coupling_hits": coupling_positive[:20],
        "positive_closure_hits": closure_positive[:20],
    }


def build_payload(p2899: dict[str, Any]) -> dict[str, Any]:
    scan = scan_generated_json()
    exact_missing_object_found = scan["coupled_source_and_coupling_artifact_count"] > 0
    return {
        "status": "P2900_DEFECT_PLACEMENT_SOURCE_LAW_INVENTORY_NO_GO_NO_CLOSURE" if not exact_missing_object_found else "P2900_DEFECT_PLACEMENT_SOURCE_LAW_REQUIRES_REVIEW",
        "input_hashes": {"P2899": sha(P2899)},
        "defect_placement_source_law_inventory": {
            "input_status_rechecked": p2899.get("status"),
            "candidate_class": "current generated-artifact corpus scan for strict defect-placement/basepoint/polarity source law coupled to the 9/5 variational density",
            "scan": scan,
        },
        "acceptance_matrix": {
            "p2899_rechecked": p2899.get("status") == "P2899_POST_DEFECT_POTENTIAL_READINESS_MATRIX_NO_CLOSURE",
            "exact_missing_object_found": exact_missing_object_found,
            "positive_source_hits_found": scan["positive_source_hit_count"] > 0,
            "positive_coupling_hits_found": scan["positive_coupling_hit_count"] > 0,
            "accepted_as_current_artifact_inventory_no_go": not exact_missing_object_found,
        },
        "decision": {
            "negative_export_flags": {
                "strict_defect_placement_source_law_exported": False,
                "nonimported_basepoint_or_polarity_law_exported": False,
                "coupling_to_9_over_5_variational_density_exported": False,
                "nonimported_9_over_5_variational_chain_rule_exported": False,
                "localized_action_density_exported": False,
                "strict_damping_compression_bridge_exported": False,
                "full_kernel_bridge_exported": False,
                "role_transfer_exported": False,
                "ltotal_exported": False,
                "eom_closure_exported": False,
                "hamiltonian_closure_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2900 scans the current generated JSON corpus for the exact P2899 missing object: a strict defect-placement/basepoint/polarity source law coupled to the 9/5 variational density.  The scan finds no unquarantined artifact with both positive source and coupling exports, and no artifact with source+coupling+closure exports.  Therefore the current corpus does not supply the missing law, localized action density, L_total, EOM, Hamiltonian, role transfer, or ToE closure.",
            "next_honest_step": "Do not replay torsor capacity, scalar scores, circulant relations, one-edge defects, generated-inventory keyword hits, or conditional potential language as closure.  A next proof-grade move must either introduce a new explicit formula/artifact for a strict defect-placement source law with computed nonconventional sign/phase and a coupling theorem to the 9/5 variational density, or pivot to a genuinely different typed object outside torsor/basepoint/scalar-score/relation/defect/support/orbit/Fourier/inventory constructions; otherwise preserve no-new-live-frontier.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    scan = payload["defect_placement_source_law_inventory"]["scan"]
    lines = [
        "# P2900/S1850 defect-placement source-law inventory no-go",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Generated JSON scan",
        f"- generated JSON files scanned: `{scan['generated_json_file_count']}`",
        f"- source/coupling/closure term occurrences: `{scan['term_occurrence_count']}`",
        f"- positive source hits: `{scan['positive_source_hit_count']}`",
        f"- positive coupling hits: `{scan['positive_coupling_hit_count']}`",
        f"- positive closure hits: `{scan['positive_closure_hit_count']}`",
        f"- artifacts with source+coupling positives: `{scan['coupled_source_and_coupling_artifact_count']}`",
        f"- artifacts with source+coupling+closure positives: `{scan['coupled_source_coupling_and_closure_artifact_count']}`",
        "",
        "## Boundary",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload(read_json(P2899))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2900/S1850 defect-placement source-law inventory no-go",
        "## P2900/S1850 defect-placement source-law inventory no-go\n\n"
        "`P2900/S1850` scans the current generated JSON corpus for the exact `P2899` missing object: a strict defect-placement/basepoint/polarity source law coupled to the `9/5` variational density.  The scan finds no unquarantined artifact with both positive source and coupling exports, and no artifact with source+coupling+closure exports.  The current corpus therefore does not supply the missing law, localized action density, strict damping bridge, `L_total`, EOM, Hamiltonian, role transfer, or ToE closure.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2900/S1850 defect-placement source-law inventory `L_total` guard",
        "## P2900/S1850 defect-placement source-law inventory `L_total` guard\n\n"
        "`P2900/S1850` is a current-artifact inventory certificate, not a strict action construction.  It adds no strict defect-placement source law, no localized unit-bearing density, no coupling theorem to the `9/5` carrier, and no variational chain rule into nonproxy `L_total`, EOM, or Hamiltonian source.\n",
    )
    append_once(
        AGENTS,
        "Current defect-placement source-law inventory guardrail (P2900/S1850, 2026-06-19)",
        "## Current defect-placement source-law inventory guardrail (P2900/S1850, 2026-06-19)\n\n"
        "- P2900 scans the current generated JSON corpus for the exact P2899 missing object: a strict defect-placement/basepoint/polarity source law coupled to the `9/5` variational density.\n"
        "- The scan finds no unquarantined artifact with both positive source and coupling exports, and no source+coupling+closure export; inventory language does not supply the missing law.\n"
        "- Do not promote torsor capacity, scalar scores, circulant relations, one-edge defects, generated-inventory keyword hits, conditional potential language, support/orbit/Fourier data, or inventory hits to strict phase/origin sourcehood, strict damping/compression bridge, role transfer, `L_total`, EOM, Hamiltonian, or ToE closure.\n"
        "- A next admissible proof-grade move must introduce a new explicit formula/artifact for a strict defect-placement source law with computed nonconventional sign/phase and coupling to the `9/5` variational density, pivot to a genuinely different typed object outside torsor/basepoint/scalar-score/relation/defect/support/orbit/Fourier/inventory constructions, or preserve no-new-live-frontier.\n",
    )
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
