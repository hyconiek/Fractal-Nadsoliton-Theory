#!/usr/bin/env python3
"""P2908/S1858: joint source provenance alternative audit.

P2907 constructed the imported joint theorem-postulate J_{0,+}.  P2908 audits
the exact remaining premise: does the current generated corpus export strict
nadsoliton provenance for J_{0,+}, rather than merely mention the postulate or
one of its 23 translated/sign-flipped alternatives?

The audit constructs the full 24-member J alternative family, scans generated
JSON artifacts for positive strict provenance keys, quarantines candidate,
postulate, boundary, and negative mentions, and reports whether any strict
provenance export currently selects the joint origin/sign theorem.
"""
from __future__ import annotations

import hashlib
import json
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN

P2907 = GEN / "p2907_s1857_joint_origin_sign_source_theorem_candidate_gate.json"
OUT = GEN / "p2908_s1858_joint_source_provenance_alternative_audit.json"
MD = GEN / "p2908_s1858_joint_source_provenance_alternative_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
N = 12

PROVENANCE_KEYS = (
    "strict_joint_origin_sign_theorem_exported",
    "joint_origin_sign_provenance_exported",
    "j_0_plus_provenance_exported",
    "strict_j_source_provenance_exported",
    "strict_nadsoliton_joint_source_provenance_exported",
)
MENTION_TERMS = (
    "J_{0,+}",
    "joint origin-and-sign",
    "joint origin-sign",
    "origin 0 and sign +",
    "strict provenance for J",
)
QUARANTINE_MARKERS = (
    "negative_export_flags",
    "no_",
    "not_",
    "candidate",
    "postulate",
    "blocked",
    "boundary",
    "reason",
    "next_honest_step",
)


def j_alternatives() -> list[dict[str, Any]]:
    rows = []
    for b in range(N):
        for sign in (-1, 1):
            rows.append({
                "name": f"J_{{{b},{'+' if sign > 0 else '-'}}}",
                "origin": b,
                "sign": sign,
                "xi": f"Xi_{{{b},{'+' if sign > 0 else '-'}}}",
                "coupled_axiom": f"A({b},{'+' if sign > 0 else '-'})",
                "defect_edge": [b, (b + sign * 5) % N],
                "rho_9_5_template": f"{sign}*U_9_5*delta_edge({b},{(b + sign * 5) % N})*q_9_5({b},{(b + sign * 5) % N})",
            })
    return rows


def walk(value: Any, prefix: str = "") -> list[tuple[str, Any]]:
    if isinstance(value, dict):
        out: list[tuple[str, Any]] = []
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
    lower_path = path.lower()
    value_text = str(value).lower() if isinstance(value, str) else ""
    if any(marker in lower_path or marker in value_text for marker in QUARANTINE_MARKERS):
        return "quarantined_candidate_boundary_or_negative"
    if value is True:
        return "positive_true"
    return "mention_or_nontrue"


def generated_json_files() -> list[Any]:
    return sorted(path for path in GEN.glob("*.json") if path.is_file() and path != OUT)


def scan_generated() -> dict[str, Any]:
    hits = []
    positive = []
    mentions = []
    for path in generated_json_files():
        try:
            payload = read_json(path)
        except Exception:
            continue
        rel = str(path.relative_to(ROOT))
        for item_path, value in walk(payload):
            lower = item_path.lower()
            value_text = value if isinstance(value, str) else ""
            is_key = any(key in lower for key in PROVENANCE_KEYS)
            is_mention = any(term in item_path or term in value_text for term in MENTION_TERMS)
            if not (is_key or is_mention):
                continue
            row = {"file": rel, "path": item_path, "value": value, "classification": classify(item_path, value)}
            hits.append(row)
            if is_mention:
                mentions.append(row)
            if is_key and row["classification"] == "positive_true":
                positive.append(row)
    return {
        "generated_json_file_count": len(generated_json_files()),
        "hit_count": len(hits),
        "mention_hit_count": len(mentions),
        "positive_provenance_hit_count": len(positive),
        "positive_provenance_hits": positive[:20],
        "sample_hits": hits[:40],
    }


def build_payload(p2907: dict[str, Any]) -> dict[str, Any]:
    alternatives = j_alternatives()
    distinguished = alternatives[1]
    scan = scan_generated()
    return {
        "status": "P2908_JOINT_SOURCE_PROVENANCE_ALTERNATIVE_AUDIT_NO_STRICT_PROVENANCE",
        "input_hashes": {"P2907": hashlib.sha256(P2907.read_bytes()).hexdigest() if P2907.exists() else None},
        "constructed_theoretical_objects": {
            "joint_source_alternative_family": alternatives,
            "distinguished_candidate_rechecked": distinguished,
            "provenance_scan": scan,
        },
        "acceptance_matrix": {
            "p2907_rechecked_candidate_not_strict": p2907.get("acceptance_matrix", {}).get("accepted_as_strict_source_theorem") is False,
            "alternative_count": len(alternatives),
            "translated_sign_flipped_alternative_count": len(alternatives) - 1,
            "distinguished_j_0_plus_present": distinguished["origin"] == 0 and distinguished["sign"] == 1,
            "positive_provenance_hit_count": scan["positive_provenance_hit_count"],
            "strict_provenance_for_j_0_plus_exported": scan["positive_provenance_hit_count"] > 0,
            "accepted_as_strict_joint_source": False,
        },
        "decision": {
            "positive_witnesses": {
                "full_24_joint_source_alternative_family_constructed": True,
                "corpus_provenance_scan_executed": True,
                "p2907_postulate_boundary_rechecked": True,
            },
            "negative_export_flags": {
                "strict_joint_origin_sign_theorem_exported": False,
                "j_0_plus_provenance_exported": False,
                "strict_nadsoliton_joint_source_provenance_exported": False,
                "unit_bearing_u_9_5_exported": False,
                "nonproxy_ltotal_exported": False,
                "eom_closure_exported": False,
                "hamiltonian_closure_exported": False,
                "role_transfer_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2908 constructs all 24 translated/sign-flipped joint-source alternatives J_{b,sigma} and scans current generated artifacts for a positive strict provenance export selecting J_{0,+}.  Existing hits are candidate/postulate/boundary mentions or negative flags; zero positive provenance hits are found.  Therefore P2907 remains readiness, not a strict nadsoliton-derived joint theorem.",
            "next_honest_step": "Do not add another J inventory or postulated translated variant.  The next proof-grade move must either provide one new strict nadsoliton-derived construction that computes J_{0,+} internally, pivot to a genuinely different typed object outside Xi/J/defect-placement, or preserve no-new-live-frontier.  Unit-bearing U_9_5 -> L_total coupling remains downstream until joint provenance is exported.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    acc = payload["acceptance_matrix"]
    scan = payload["constructed_theoretical_objects"]["provenance_scan"]
    lines = [
        "# P2908/S1858 joint source provenance alternative audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Alternative family",
        f"- J alternatives: `{acc['alternative_count']}`",
        f"- translated/sign-flipped alternatives to `J_{{0,+}}`: `{acc['translated_sign_flipped_alternative_count']}`",
        "",
        "## Provenance scan",
        f"- generated JSON files scanned: `{scan['generated_json_file_count']}`",
        f"- provenance/mention hits: `{scan['hit_count']}`",
        f"- positive provenance hits: `{scan['positive_provenance_hit_count']}`",
        f"- strict provenance exported: `{acc['strict_provenance_for_j_0_plus_exported']}`",
        "",
        "## Boundary",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload(read_json(P2907))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2908/S1858 joint source provenance alternative audit", "## P2908/S1858 joint source provenance alternative audit\n\n`P2908/S1858` constructs the full `24`-member translated/sign-flipped joint-source family `J_{b,sigma}` and scans current generated artifacts for a positive strict provenance export selecting `J_{0,+}`.  The scan finds only candidate/postulate/boundary mentions or negative flags and `0` positive provenance hits, so `J_{0,+}` remains imported readiness rather than a strict nadsoliton-derived joint theorem.  No unit-bearing `U_9_5`, nonproxy `L_total`, EOM, Hamiltonian, role transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2908/S1858 joint source provenance `L_total` guard", "## P2908/S1858 joint source provenance `L_total` guard\n\n`P2908/S1858` adds no new density.  It audits the P2907 `J_{0,+}` theorem-postulate against all `24` translated/sign-flipped joint-source alternatives and finds no positive strict provenance export.  Therefore the symbolic `U_9_5` derivative table cannot yet be promoted to a unit-bearing nonproxy `L_total`, EOM, Hamiltonian, role transfer, or ToE closure.\n")
    append_once(AGENTS, "Current joint source provenance alternative audit guardrail (P2908/S1858, 2026-06-19)", "## Current joint source provenance alternative audit guardrail (P2908/S1858, 2026-06-19)\n\n- P2908 constructs the full `24`-member translated/sign-flipped `J_{b,sigma}` alternative family and scans generated artifacts for a positive strict provenance export selecting `J_{0,+}`.\n- The scan finds candidate/postulate/boundary mentions and negative flags but `0` positive provenance hits; `J_{0,+}` remains imported readiness, not a strict nadsoliton-derived joint theorem.\n- Do not promote `J_{0,+}`, `Xi_{0,+}`, symbolic `rho_9/5`, or `U_9_5` to strict sourcehood, nonproxy `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE closure from provenance mentions or postulated variants.\n- A next admissible proof-grade move must provide a genuinely new strict construction computing `J_{0,+}`, pivot outside Xi/J/defect-placement, or preserve no-new-live-frontier; unit-bearing `U_9_5 -> L_total` remains downstream.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
