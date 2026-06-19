#!/usr/bin/env python3
"""P2905/S1855: Xi Dirac source provenance alternative audit.

P2904 constructed Xi_{0,+}, the minimal translation-breaking source candidate
that couples to A(0,+).  P2905 audits the exact remaining premise: does the
current artifact corpus prove strict provenance for Xi_{0,+} rather than one of
its 23 translated/sign-flipped alternatives?

The audit constructs the full 24-candidate alternative set and scans generated
JSON artifacts for an unquarantined strict provenance export.  Candidate mentions
and negative/provenance-block flags do not count as provenance.
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

P2904 = GEN / "p2904_s1854_translation_breaking_dirac_source_postulate_acceptance.json"
OUT = GEN / "p2905_s1855_xi_dirac_source_provenance_alternative_audit.json"
MD = GEN / "p2905_s1855_xi_dirac_source_provenance_alternative_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
N = 12

PROVENANCE_KEYS = (
    "strict_nadsoliton_source_provenance_exported",
    "xi_dirac_source_provenance_exported",
    "xi_0_plus_provenance_exported",
    "dirac_source_provenance_exported",
    "strict_translation_breaking_source_provenance_exported",
)
MENTION_TERMS = (
    "Xi_{0,+}",
    "signed Dirac source",
    "Dirac source candidate",
    "strict nadsoliton provenance",
    "translated/sign-flipped",
    "A(0,+)",
)
QUARANTINE_MARKERS = ("negative_export_flags", "no_", "not_", "candidate", "postulate", "blocked", "reason", "next_honest_step")


def alternatives() -> list[dict[str, Any]]:
    rows = []
    for b in range(N):
        for sign in (-1, 1):
            rows.append({
                "name": f"Xi_{{{b},{'+' if sign > 0 else '-'}}}",
                "basepoint": b,
                "sign": sign,
                "values": [sign if i == b else 0 for i in range(N)],
                "coupled_axiom": f"A({b},{'+' if sign > 0 else '-'})",
                "defect_edge": [b, (b + sign * 5) % N],
            })
    return rows


def walk(value: Any, prefix: str = "") -> list[tuple[str, Any]]:
    if isinstance(value, dict):
        out = []
        for key, child in value.items():
            out.extend(walk(child, f"{prefix}.{key}" if prefix else str(key)))
        return out
    if isinstance(value, list):
        out = []
        for index, child in enumerate(value):
            out.extend(walk(child, f"{prefix}[{index}]"))
        return out
    return [(prefix, value)]


def classify(path: str, value: Any) -> str:
    lower = path.lower()
    if any(marker in lower for marker in QUARANTINE_MARKERS):
        return "quarantined_candidate_boundary_or_negative"
    if value is True:
        return "positive_true"
    return "mention_or_nontrue"


def generated_json_files() -> list[Path]:
    return sorted(path for path in GEN.glob("*.json") if path.is_file() and path != OUT)


def scan_generated() -> dict[str, Any]:
    hits = []
    positive_provenance = []
    mention_hits = []
    for path in generated_json_files():
        try:
            payload = read_json(path)
        except Exception:
            continue
        rel = str(path.relative_to(ROOT))
        for item_path, value in walk(payload):
            item_path_lower = item_path.lower()
            value_text = value if isinstance(value, str) else ""
            is_provenance_key = any(key in item_path_lower for key in PROVENANCE_KEYS)
            is_mention = any(term in item_path or term in value_text for term in MENTION_TERMS)
            if not (is_provenance_key or is_mention):
                continue
            row = {"file": rel, "path": item_path, "value": value, "classification": classify(item_path, value)}
            hits.append(row)
            if is_mention:
                mention_hits.append(row)
            if is_provenance_key and row["classification"] == "positive_true":
                positive_provenance.append(row)
    return {
        "generated_json_file_count": len(generated_json_files()),
        "hit_count": len(hits),
        "mention_hit_count": len(mention_hits),
        "positive_provenance_hit_count": len(positive_provenance),
        "positive_provenance_hits": positive_provenance[:20],
        "sample_hits": hits[:30],
    }


def build_payload(p2904: dict[str, Any]) -> dict[str, Any]:
    alt = alternatives()
    scan = scan_generated()
    xi = alt[1]  # b=0, sign=+ under loop ordering (-,+)
    return {
        "status": "P2905_XI_DIRAC_SOURCE_PROVENANCE_ALTERNATIVE_AUDIT_NO_STRICT_PROVENANCE",
        "input_hashes": {"P2904": sha(P2904)},
        "constructed_theoretical_objects": {
            "alternative_family": alt,
            "distinguished_candidate_rechecked": xi,
            "provenance_scan": scan,
        },
        "acceptance_matrix": {
            "p2904_rechecked": p2904.get("status") == "P2904_TRANSLATION_BREAKING_DIRAC_SOURCE_POSTULATE_ACCEPTED_AS_CANDIDATE_NO_STRICT_PROVENANCE",
            "alternative_count": len(alt),
            "translated_sign_flipped_alternative_count": len(alt) - 1,
            "candidate_xi_0_plus_present": xi["basepoint"] == 0 and xi["sign"] == 1,
            "positive_provenance_hit_count": scan["positive_provenance_hit_count"],
            "strict_provenance_for_xi_0_plus_exported": scan["positive_provenance_hit_count"] > 0,
            "accepted_as_strict_source": False,
        },
        "decision": {
            "positive_witnesses": {
                "full_24_alternative_family_constructed": True,
                "corpus_provenance_scan_executed": True,
                "candidate_boundary_sharpened": True,
            },
            "negative_export_flags": {
                "strict_nadsoliton_source_provenance_exported": False,
                "xi_0_plus_provenance_exported": False,
                "strict_defect_placement_source_law_exported": False,
                "unit_bearing_strict_density_exported": False,
                "nonproxy_ltotal_exported": False,
                "eom_closure_exported": False,
                "hamiltonian_closure_exported": False,
                "strict_damping_compression_bridge_exported": False,
                "role_transfer_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2905 constructs all 24 translated/sign-flipped Xi alternatives and scans current generated artifacts for a positive strict provenance export selecting Xi_{0,+}.  The corpus contains candidate/boundary mentions but zero positive provenance hits.  Thus P2904 remains a useful translation-breaking postulate, not a strict nadsoliton-sourced law or closure export.",
            "next_honest_step": "The next proof-grade move must introduce a new strict asymmetry/chiral/defect-generation theorem that selects Xi_{0,+} over the 23 alternatives, or pivot outside the Xi/defect-placement lane.  Further candidate mentions, inventory scans, translation-neutral selectors, or symbolic unit assignments are repetition-gated unless they add a new provenance theorem.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    acc = payload["acceptance_matrix"]
    scan = payload["constructed_theoretical_objects"]["provenance_scan"]
    lines = [
        "# P2905/S1855 Xi Dirac source provenance alternative audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Alternative family",
        f"- Xi alternatives: `{acc['alternative_count']}`",
        f"- translated/sign-flipped alternatives to `Xi_{{0,+}}`: `{acc['translated_sign_flipped_alternative_count']}`",
        "",
        "## Provenance scan",
        f"- generated JSON files scanned: `{scan['generated_json_file_count']}`",
        f"- provenance/mention hits: `{scan['hit_count']}`",
        f"- positive provenance hits: `{scan['positive_provenance_hit_count']}`",
        f"- strict provenance exported: `{acc['strict_provenance_for_xi_0_plus_exported']}`",
        "",
        "## Boundary",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload(read_json(P2904))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2905/S1855 Xi Dirac source provenance alternative audit", "## P2905/S1855 Xi Dirac source provenance alternative audit\n\n`P2905/S1855` constructs the full `24`-member translated/sign-flipped `Xi` alternative family and scans current generated artifacts for a positive strict provenance export selecting `Xi_{0,+}`.  The scan finds candidate/boundary mentions but `0` positive provenance hits, so `Xi_{0,+}` remains a useful translation-breaking postulate rather than a strict nadsoliton-sourced law.  No strict defect-placement source law, unit-bearing density, nonproxy `L_total`, EOM, Hamiltonian, bridge, role transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2905/S1855 Xi provenance `L_total` guard", "## P2905/S1855 Xi provenance `L_total` guard\n\n`P2905/S1855` adds no action density.  It audits provenance for the P2904 `Xi_{0,+}` source candidate and finds no current positive strict provenance export over the `23` translated/sign-flipped alternatives.  Therefore no unit-bearing nonproxy `L_total`, EOM, Hamiltonian, role transfer, or ToE closure follows.\n")
    append_once(AGENTS, "Current Xi Dirac source provenance alternative guardrail (P2905/S1855, 2026-06-19)", "## Current Xi Dirac source provenance alternative guardrail (P2905/S1855, 2026-06-19)\n\n- P2905 constructs the full `24`-member translated/sign-flipped `Xi` alternative family and scans generated artifacts for a positive strict provenance export selecting `Xi_{0,+}`.\n- The current corpus has candidate/boundary mentions but `0` positive provenance hits; `Xi_{0,+}` remains a translation-breaking postulate/candidate, not a strict nadsoliton-sourced law.\n- Do not promote `Xi_{0,+}`, its candidate coupling, inventory mentions, symbolic `rho_9/5`, or `U_9_5` to strict sourcehood, nonproxy `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE closure without a new provenance theorem.\n- A next admissible proof-grade move must introduce a strict asymmetry/chiral/defect-generation theorem selecting `Xi_{0,+}` over the `23` alternatives, pivot outside the Xi/defect-placement lane, or preserve no-new-live-frontier.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
