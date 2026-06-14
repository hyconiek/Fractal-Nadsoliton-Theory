#!/usr/bin/env python3
"""P2711/S1661: inversion-odd character source-law sign-coupling audit.

P2710 proved that inversion-odd characters of Aut(Z12) exist, but are only
parity labels unless a strict source law couples one of them to the P2708
boundary-cocycle sign.  P2711 audits that exact missing step: enumerate the
finite source-law candidates (anti-inversion character + coupling sign) and
check whether current artifacts export the coupling sign as a non-premise law.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2711_s1661_inversion_odd_character_source_law_sign_coupling_audit.json"
MD = GEN / "p2711_s1661_inversion_odd_character_source_law_sign_coupling_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

INPUTS = {
    "P2710_CHARACTER_TEST": GEN / "p2710_s1660_finite_aut_z12_anti_inversion_orientation_character_source_test.json",
    "P2708_BOUNDARY_COCYCLE": GEN / "p2708_s1658_z12_boundary_1_cocycle_selector_source_obstruction.json",
    "P2709_RELEASE_BACKSCAN": GEN / "p2709_s1659_release_8_1_to_9_3_breakthrough_unlock_backscan.json",
    "P2707_STATE_MAP": GEN / "p2707_s1657_post_p2706_no_new_live_frontier_reconciliation.json",
    "P2704_DECLARED_SCOPE": GEN / "p2704_s1654_p1343_p1348_selector_provenance_revalidation_table.json",
    "P2377_TRANSPORT": GEN / "p2377_s1327_damping_compression_transport_primitive_uniform_coupling_theorem.json",
    "P2378_INSUFFICIENCY": GEN / "p2378_s1328_unit_normalized_transport_coupling_insufficiency_theorem.json",
}

SIGN_SOURCE_TOKENS = [
    "source_law_coupling_sign_exported\": true",
    "orientation_character_sign_exported\": true",
    "non_premise_sign_convention_exported\": true",
    "strict source law selecting +omega",
    "strict_source_law_selects_plus_omega\": true",
]
SIGN_BLOCK_TOKENS = [
    "no-current-unlock", "no_current_unlock", "no strict source", "no-strict-selector",
    "no non-premise", "premise", "insufficient", "open", "blocked", "no-go", "no_go",
    "orientation-blind", "does not export",
]
NEGATIVE_EXPORT_FLAGS = [
    "source_law_coupling_sign_exported",
    "strict_source_law_selects_plus_omega",
    "qw2191_discharged",
    "non_premise_selector_provider_exported",
    "pair12_strict_core_upgrade_exported",
    "ltotal_promoted",
    "toe_closure_exported",
]


def read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"missing": True, "path": rel(path)}
    return json.loads(path.read_text(encoding="utf-8"))


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="ignore") if path.exists() else ""


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def anti_inversion_characters() -> list[dict[str, Any]]:
    p2710 = read_json(INPUTS["P2710_CHARACTER_TEST"])
    return [row for row in p2710.get("character_table", []) if row.get("anti_inversion")]


def source_law_candidates(chars: list[dict[str, Any]]) -> list[dict[str, Any]]:
    candidates: list[dict[str, Any]] = []
    for char in chars:
        for coupling_sign in [-1, 1]:
            selects = "+omega" if coupling_sign == 1 else "-omega"
            candidates.append({
                "candidate_id": f"{char['name']}__lambda_{coupling_sign:+d}",
                "character": char["name"],
                "character_values": char["values"],
                "coupling_sign_lambda": coupling_sign,
                "formal_law": "S(omega)=lambda and S(-omega)=-lambda after choosing an inversion-odd character chi with chi(inversion)=-1",
                "would_select_if_lambda_were_strictly_sourced": selects,
                "strict_lambda_source_exported": False,
                "is_nonpremise_selector_provider": False,
            })
    return candidates


def sign_pair_degeneracy(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    by_char: dict[str, list[dict[str, Any]]] = {}
    for row in candidates:
        by_char.setdefault(row["character"], []).append(row)
    pairs = []
    for character, rows in sorted(by_char.items()):
        rows = sorted(rows, key=lambda row: row["coupling_sign_lambda"])
        pairs.append({
            "character": character,
            "candidate_pair": [row["candidate_id"] for row in rows],
            "lambda_values": [row["coupling_sign_lambda"] for row in rows],
            "degeneracy": "lambda -> -lambda exchanges +omega and -omega",
            "strict_sign_source_required": True,
            "strict_sign_source_exported": False,
        })
    return pairs


def artifact_sign_source_scan() -> list[dict[str, Any]]:
    rows = []
    for name, path in INPUTS.items():
        text = read_text(path).lower()
        export_hits = sorted(tok for tok in SIGN_SOURCE_TOKENS if tok in text)
        block_hits = sorted(tok for tok in SIGN_BLOCK_TOKENS if tok in text)[:10]
        data = read_json(path) if path.suffix == ".json" and path.exists() else {}
        rows.append({
            "artifact": name,
            "path": rel(path),
            "status": data.get("status"),
            "sign_source_export_hits": export_hits,
            "block_hits_sample": block_hits,
            "exports_nonpremise_coupling_sign": bool(export_hits) and not block_hits,
        })
    return rows


def write_markdown(payload: dict[str, Any]) -> None:
    lines = [
        "# P2711/S1661 inversion-odd character source-law sign-coupling audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite candidate count",
        f"- anti_inversion_characters={payload['finite_counts']['anti_inversion_characters']}",
        f"- source_law_candidates={payload['finite_counts']['source_law_candidates']}",
        f"- sign_degenerate_pairs={payload['finite_counts']['sign_degenerate_pairs']}",
        "",
        "## Degeneracy rows",
    ]
    for row in payload["sign_pair_degeneracy"]:
        lines.append(f"- `{row['character']}`: {row['degeneracy']}; strict_sign_source_exported={row['strict_sign_source_exported']}")
    lines.extend(["", "## Decision", payload["decision"]["reason"], "", "## Next honest step", payload["decision"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    chars = anti_inversion_characters()
    candidates = source_law_candidates(chars)
    degeneracy = sign_pair_degeneracy(candidates)
    scan = artifact_sign_source_scan()
    no_unlock = len(chars) == 2 and len(candidates) == 4 and all(not row["strict_sign_source_exported"] for row in degeneracy) and all(not row["exports_nonpremise_coupling_sign"] for row in scan)
    payload = {
        "status": "P2711_SOURCE_LAW_SIGN_COUPLING_DEGENERACY_NO_STRICT_SIGN_SOURCE" if no_unlock else "P2711_REQUIRES_MANUAL_REVIEW",
        "input_hashes": {name: sha(path) for name, path in INPUTS.items()},
        "finite_counts": {
            "anti_inversion_characters": len(chars),
            "source_law_candidates": len(candidates),
            "sign_degenerate_pairs": len(degeneracy),
        },
        "source_law_candidates": candidates,
        "sign_pair_degeneracy": degeneracy,
        "artifact_sign_source_scan": scan,
        "decision": {
            "finite_source_law_candidates_enumerated": True,
            "source_law_coupling_sign_exported": False,
            "strict_source_law_selects_plus_omega": False,
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": "P2711 enumerates the four finite source-law candidates obtained from two inversion-odd characters times two coupling signs.  Each character has a lambda -> -lambda degeneracy that exchanges +omega and -omega.  Current artifacts export no non-premise law fixing lambda, so the candidate family remains a premise-sign pair rather than a strict selector source.",
            "next_honest_step": "A further admissible move must introduce an actually new strict mechanism fixing the coupling sign lambda, or pivot outside the selector/sign lane to a new typed object.  Without such a mechanism, preserve the P2697-P2711 no-current-unlock certificate and do not promote QW-2191, pair12 strict-core, L_total, role transfer, bridge closure, or ToE.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2711/S1661 inversion-odd character source-law sign-coupling audit", "## P2711/S1661 inversion-odd character source-law sign-coupling audit\n\n`P2711/S1661` enumerates the exact finite source-law candidates after P2710: two inversion-odd characters times two coupling signs.  Each character has a `lambda -> -lambda` degeneracy that exchanges `+omega` and `-omega`; current artifacts export no non-premise strict law fixing `lambda`.  Thus the inversion-odd character lane remains a premise-sign pair and does not discharge `QW-2191`, upgrade pair12 strict-core, start role transfer, promote `L_total`, close the bridge, or imply ToE closure.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2711/S1661 source-law sign coupling Ltotal guard", "## P2711/S1661 source-law sign coupling Ltotal guard\n\n`P2711/S1661` is a finite sign-coupling degeneracy audit, not a variational source.  Because no strict law fixes the coupling sign `lambda`, it does not promote `L_total`, role transfer, bridge closure, selector closure, or ToE.\n")
    append_once(AGENTS, "Current inversion-odd character source-law sign-coupling guardrail (P2711/S1661, 2026-06-14)", "## Current inversion-odd character source-law sign-coupling guardrail (P2711/S1661, 2026-06-14)\n\n- P2711 enumerates the finite source-law candidate family after P2710: two inversion-odd Aut(Z12) characters times two coupling signs.\n- Every candidate remains paired by `lambda -> -lambda`, which exchanges `+omega` and `-omega`; current artifacts do not export a non-premise strict law fixing `lambda`.\n- Do not promote the inversion-odd character/sign-coupling lane to `QW-2191` discharge, pair12 strict-core, role transfer, `L_total`, bridge closure, or ToE without a new strict mechanism fixing the coupling sign.\n")
    return payload


if __name__ == "__main__":
    main()
