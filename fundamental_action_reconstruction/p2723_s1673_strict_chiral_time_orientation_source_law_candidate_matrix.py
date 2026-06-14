#!/usr/bin/env python3
"""P2723/S1673: strict chiral/time-orientation source-law candidate matrix.

P2722 showed that a moving clock hand is not enough: the next admissible move
requires an actual strict chiral/time-orientation source law with a nonzero
signed value coupled to the P2721 polarity pair.  P2723 audits the currently
available concrete candidates against that exact standard and refuses to turn a
clock-direction intuition into selector closure.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2723_s1673_strict_chiral_time_orientation_source_law_candidate_matrix.json"
MD = GEN / "p2723_s1673_strict_chiral_time_orientation_source_law_candidate_matrix.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

INPUTS = {
    "P2722_CLOCK_PHASE_SOURCE": GEN / "p2722_s1672_clock_phase_polarity_selection_source_law_audit.json",
    "P2721_SIGN_TORSOR_COUPLING": GEN / "p2721_s1671_chiral_bispectrum_sign_to_orientation_torsor_coupling_audit.json",
    "P2711_CHARACTER_SIGN_COUPLING": GEN / "p2711_s1661_inversion_odd_character_source_law_sign_coupling_audit.json",
    "P2708_BOUNDARY_COCYCLE": GEN / "p2708_s1658_z12_boundary_1_cocycle_selector_source_obstruction.json",
}

STRICT_CRITERIA = [
    "strict_artifact_exports_source_law",
    "nonzero_signed_value_exported",
    "sign_not_convention_or_time_arrow_premise",
    "coupled_to_p2721_polarity_pair",
    "selects_exactly_one_polarity",
]

CANDIDATES = [
    {
        "name": "external_time_arrow_or_clockwise_premise",
        "strict_artifact_exports_source_law": False,
        "nonzero_signed_value_exported": True,
        "sign_not_convention_or_time_arrow_premise": False,
        "coupled_to_p2721_polarity_pair": False,
        "selects_exactly_one_polarity": True,
        "blocker": "It can choose a direction only by premise/convention, not by a strict exported internal source law.",
    },
    {
        "name": "chiral_bispectrum_temporal_phase_velocity",
        "strict_artifact_exports_source_law": False,
        "nonzero_signed_value_exported": False,
        "sign_not_convention_or_time_arrow_premise": True,
        "coupled_to_p2721_polarity_pair": False,
        "selects_exactly_one_polarity": False,
        "blocker": "Current P2718-P2722 artifacts contain static source/orientation rows, not an exported time-series/EOM law for d phase/dt with a nonzero signed value.",
    },
    {
        "name": "inversion_odd_character_lambda_sign_law",
        "strict_artifact_exports_source_law": False,
        "nonzero_signed_value_exported": False,
        "sign_not_convention_or_time_arrow_premise": True,
        "coupled_to_p2721_polarity_pair": False,
        "selects_exactly_one_polarity": False,
        "blocker": "P2711 leaves the lambda sign paired; current artifacts do not export a strict law choosing one sign and coupling it to P2721.",
    },
    {
        "name": "boundary_cocycle_orientation_flow",
        "strict_artifact_exports_source_law": False,
        "nonzero_signed_value_exported": True,
        "sign_not_convention_or_time_arrow_premise": False,
        "coupled_to_p2721_polarity_pair": False,
        "selects_exactly_one_polarity": False,
        "blocker": "P2708/P2714 provide a real orientable torsor, but no non-premise source for choosing one of its two signs.",
    },
]

NEGATIVE_EXPORT_FLAGS = [
    "strict_chiral_time_orientation_source_law_exported",
    "nonzero_signed_time_orientation_value_exported",
    "canonical_coupling_polarity_selected",
    "strict_mechanism_fixing_lambda_exported",
    "qw2191_discharged",
    "pair12_strict_core_upgrade_exported",
    "bridge_closure_exported",
    "role_transfer_started",
    "ltotal_promoted",
    "toe_closure_exported",
]


def read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"missing": True, "path": rel(path)}
    return json.loads(path.read_text(encoding="utf-8"))


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def evaluate_candidate(candidate: dict[str, Any]) -> dict[str, Any]:
    missing = [criterion for criterion in STRICT_CRITERIA if not candidate[criterion]]
    return {
        **candidate,
        "missing_criteria": missing,
        "accepted_as_strict_chiral_time_source": not missing,
    }


def source_law_matrix() -> dict[str, Any]:
    rows = [evaluate_candidate(candidate) for candidate in CANDIDATES]
    accepted = [row for row in rows if row["accepted_as_strict_chiral_time_source"]]
    return {
        "criteria": STRICT_CRITERIA,
        "candidate_rows": rows,
        "accepted_candidate_count": len(accepted),
        "strict_chiral_time_orientation_source_law_exported": len(accepted) > 0,
        "blocker": "No current candidate exports all required data: an internal source law, a nonzero signed value, a nonconventional sign, coupling to the P2721 polarity pair, and selection of exactly one polarity.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lines = [
        "# P2723/S1673 strict chiral/time-orientation source-law candidate matrix",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Matrix result",
        f"- accepted_candidate_count={payload['source_law_matrix']['accepted_candidate_count']}",
        payload["source_law_matrix"]["blocker"],
        "",
        "## Decision",
        payload["decision"]["reason"],
        "",
        "## Next honest step",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    inputs = {key: read_json(path) for key, path in INPUTS.items()}
    matrix = source_law_matrix()
    no_unlock = not matrix["strict_chiral_time_orientation_source_law_exported"]
    payload = {
        "status": "P2723_STRICT_CHIRAL_TIME_ORIENTATION_SOURCE_LAW_MATRIX_NO_ACCEPTED_SOURCE" if no_unlock else "P2723_REQUIRES_MANUAL_REVIEW",
        "input_hashes": {key: sha(path) for key, path in INPUTS.items()},
        "input_statuses": {key: value.get("status") for key, value in inputs.items()},
        "audited_question": "Does any current artifact export a strict chiral/time-orientation source law with nonzero signed value coupled to the P2721 polarity pair?",
        "source_law_matrix": matrix,
        "decision": {
            "strict_chiral_time_orientation_source_law_exported": False,
            "nonzero_signed_time_orientation_value_exported": False,
            "canonical_coupling_polarity_selected": False,
            "strict_mechanism_fixing_lambda_exported": False,
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": "P2723 audits the concrete chiral/time-orientation source-law candidates demanded by P2722.  External time-arrow direction is premise/convention; chiral-bispectrum phase velocity lacks an exported time-series/EOM source; inversion-odd character lambda remains sign-paired; and the boundary-cocycle orientation flow lacks a non-premise sign source.  No candidate selects the P2721 polarity or fixes lambda.",
            "next_honest_step": "Do not continue naming time-arrow/chiral candidates without an explicit exported source law and nonzero signed value.  A next admissible move must introduce one new strict dynamic/chiral source artifact with a computable signed value coupled to P2721, or preserve the P2697-P2723 no-new-live-frontier certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2723/S1673 strict chiral time-orientation source-law matrix", "## P2723/S1673 strict chiral time-orientation source-law matrix\n\n`P2723/S1673` audits candidate strict chiral/time-orientation source laws after P2722: external time-arrow premise, chiral-bispectrum temporal phase velocity, inversion-odd character `lambda` sign law, and boundary-cocycle orientation flow.  None exports all required data: internal source law, nonzero signed value, nonconventional sign, coupling to the P2721 polarity pair, and exact polarity selection.  No strict `lambda` fixing, `QW-2191` discharge, pair12 strict-core upgrade, bridge closure, role transfer, `L_total`, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2723/S1673 chiral time-orientation source Ltotal guard", "## P2723/S1673 chiral time-orientation source Ltotal guard\n\n`P2723/S1673` is a source-law candidate matrix, not a variational source construction.  Since no strict dynamic/chiral source artifact with a computable nonzero signed value is exported, no `L_total`, selector closure, pair12 strict-core, role transfer, bridge closure, or ToE promotion follows.\n")
    append_once(AGENTS, "Current strict chiral/time-orientation source-law matrix guardrail (P2723/S1673, 2026-06-14)", "## Current strict chiral/time-orientation source-law matrix guardrail (P2723/S1673, 2026-06-14)\n\n- P2723 audits the current concrete chiral/time-orientation source-law candidates demanded by P2722: external time-arrow premise, chiral-bispectrum temporal phase velocity, inversion-odd character `lambda` sign law, and boundary-cocycle orientation flow.\n- No candidate exports all required strict data: internal source law, nonzero signed value, nonconventional sign, coupling to the P2721 polarity pair, and exact polarity selection; therefore `lambda` remains unfixed and `QW-2191` remains open.\n- Do not continue naming time-arrow/chiral candidates as closure evidence without an explicit exported source law and nonzero signed value.  A next admissible move must introduce one new strict dynamic/chiral source artifact with a computable signed value coupled to P2721, or preserve the P2697-P2723 no-new-live-frontier certificate.\n")
    return payload


if __name__ == "__main__":
    main()
