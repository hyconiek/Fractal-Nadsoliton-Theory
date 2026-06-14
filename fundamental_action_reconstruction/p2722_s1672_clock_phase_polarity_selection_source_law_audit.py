#!/usr/bin/env python3
"""P2722/S1672: clock-phase polarity-selection source-law audit.

P2721 left two Aut(Z12)-equivariant sign-to-orientation couplings but no strict
law selecting one polarity.  P2722 tests the natural "clock hand keeps moving"
idea: can a continuously changing phase/position source select the missing
polarity without merely choosing a clock face origin, time orientation, or
external convention?
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2722_s1672_clock_phase_polarity_selection_source_law_audit.json"
MD = GEN / "p2722_s1672_clock_phase_polarity_selection_source_law_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

INPUTS = {
    "P2721_SIGN_TORSOR_COUPLING": GEN / "p2721_s1671_chiral_bispectrum_sign_to_orientation_torsor_coupling_audit.json",
    "P2720_TRANSLATION_ORBIT_LOCALIZER": GEN / "p2720_s1670_chiral_bispectrum_translation_orbit_phase_origin_localizer_no_go.json",
    "P2719_PHASE_ORIGIN_AUDIT": GEN / "p2719_s1669_chiral_bispectrum_phase_origin_source_localizer_theorem_audit.json",
    "P2718_CHIRAL_BISPECTRUM_FORMULA": GEN / "p2718_s1668_chiral_bispectrum_explicit_signed_formula_torsor_coupling_audit.json",
}

CANDIDATE_SOURCE_LAWS = [
    {
        "name": "fixed_clock_face_zero",
        "informal_clock_reading": "declare one mark to be twelve o'clock and read all positions from it",
        "strictly_internal": False,
        "nonconventional": False,
        "origin_sensitive": True,
        "selects_polarity": True,
        "blocker": "It chooses an origin, but only as an external face-label convention; P2720 forbids importing that as a non-premise source.",
    },
    {
        "name": "moving_clock_hand_position_theta_t",
        "informal_clock_reading": "let the hand keep moving and use its current position",
        "strictly_internal": False,
        "nonconventional": False,
        "origin_sensitive": False,
        "selects_polarity": False,
        "blocker": "A changing position is origin-relative: it is meaningful only relative to a zero mark or reference event; without that anchor, a global rotation gives an equivalent description.",
    },
    {
        "name": "angular_velocity_sign_dtheta_dt",
        "informal_clock_reading": "use clockwise versus counter-clockwise motion",
        "strictly_internal": False,
        "nonconventional": False,
        "origin_sensitive": False,
        "selects_polarity": True,
        "blocker": "This would be the right kind of chiral/time-orientation datum, but current artifacts do not export it as a non-premise strict source coupled to P2721.",
    },
    {
        "name": "p2721_equivariant_coupling_pair",
        "informal_clock_reading": "map the two possible hand directions to the two possible orientation labels",
        "strictly_internal": True,
        "nonconventional": True,
        "origin_sensitive": False,
        "selects_polarity": False,
        "blocker": "P2721 exports a conditional pair of opposite polarities, not a law choosing one member of the pair.",
    },
]

NEGATIVE_EXPORT_FLAGS = [
    "strict_polarity_selection_source_law_exported",
    "clock_phase_motion_source_exported",
    "continuous_phase_origin_exported",
    "time_orientation_chiral_source_exported",
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
    accepted = all([
        candidate["strictly_internal"],
        candidate["nonconventional"],
        candidate["origin_sensitive"],
        candidate["selects_polarity"],
    ])
    return {
        **candidate,
        "accepted_as_strict_polarity_source": accepted,
        "missing_requirements": [
            requirement
            for requirement in ("strictly_internal", "nonconventional", "origin_sensitive", "selects_polarity")
            if not candidate[requirement]
        ],
    }


def clock_phase_audit(p2721: dict[str, Any]) -> dict[str, Any]:
    rows = [evaluate_candidate(candidate) for candidate in CANDIDATE_SOURCE_LAWS]
    accepted = [row for row in rows if row["accepted_as_strict_polarity_source"]]
    p2721_conditional_pair = bool(p2721.get("decision", {}).get("conditional_equivariant_coupling_exported"))
    return {
        "candidate_rows": rows,
        "accepted_candidate_count": len(accepted),
        "p2721_conditional_pair_reused": p2721_conditional_pair,
        "clock_motion_explanation": "A moving hand does not by itself solve the selector problem: if the whole clock face can be rotated, the absolute position of the hand is gauge/reference data.  Only a nonconventional origin event or a strict chiral/time-orientation source coupled to the P2721 polarity pair would select one polarity.",
        "strict_polarity_selection_source_law_exported": len(accepted) > 0,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lines = [
        "# P2722/S1672 clock-phase polarity-selection source-law audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        f"- accepted_candidate_count={payload['clock_phase_audit']['accepted_candidate_count']}",
        payload["clock_phase_audit"]["clock_motion_explanation"],
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
    audit = clock_phase_audit(inputs["P2721_SIGN_TORSOR_COUPLING"])
    no_unlock = not audit["strict_polarity_selection_source_law_exported"]
    payload = {
        "status": "P2722_CLOCK_PHASE_POLARITY_SOURCE_LAW_AUDIT_NO_UNLOCK" if no_unlock else "P2722_REQUIRES_MANUAL_REVIEW",
        "input_hashes": {key: sha(path) for key, path in INPUTS.items()},
        "input_statuses": {key: value.get("status") for key, value in inputs.items()},
        "audited_question": "Can a continuously changing clock-phase/position law select the P2721 coupling polarity non-premise?",
        "clock_phase_audit": audit,
        "decision": {
            "strict_polarity_selection_source_law_exported": False,
            "clock_phase_motion_source_exported": False,
            "continuous_phase_origin_exported": False,
            "canonical_coupling_polarity_selected": False,
            "strict_mechanism_fixing_lambda_exported": False,
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": "P2722 tests the natural moving-clock analogy as a polarity-selection source.  Continuous motion supplies change, but not an absolute non-premise phase origin or a strict polarity law: fixed clock zero is a convention, theta(t) is origin-relative, angular-velocity sign would need a new exported chiral/time-orientation source, and P2721 still leaves a two-polarity pair.  Therefore no strict lambda fixing or QW-2191 discharge is exported.",
            "next_honest_step": "Do not repeat clock-position or clock-face-origin arguments as closure.  A next admissible move must supply one actual strict chiral/time-orientation source law with a nonzero signed value coupled to the P2721 polarity pair, or pivot to a genuinely new typed object outside the closed lanes; otherwise preserve the P2697-P2722 no-new-live-frontier certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2722/S1672 clock-phase polarity-selection source-law audit", "## P2722/S1672 clock-phase polarity-selection source-law audit\n\n`P2722/S1672` tests the clock-hand analogy for the P2721 polarity problem.  A continuously changing phase/position does not by itself select a non-premise polarity: a fixed clock-face zero is convention, `theta(t)` is origin-relative, angular-velocity sign would require a new strict chiral/time-orientation source, and the P2721 coupling pair remains unselected.  No strict `lambda` fixing, `QW-2191` discharge, pair12 strict-core upgrade, bridge closure, role transfer, `L_total`, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2722/S1672 clock-phase polarity Ltotal guard", "## P2722/S1672 clock-phase polarity Ltotal guard\n\n`P2722/S1672` is a source-law audit for moving phase/clock-position intuition, not a variational source construction.  Since no strict chiral/time-orientation source law with nonzero signed value is exported, clock motion cannot promote `L_total`, selector closure, pair12 strict-core, role transfer, bridge closure, or ToE.\n")
    append_once(AGENTS, "Current clock-phase polarity-selection source-law guardrail (P2722/S1672, 2026-06-14)", "## Current clock-phase polarity-selection source-law guardrail (P2722/S1672, 2026-06-14)\n\n- P2722 tests the moving-clock/continuous-phase analogy for the P2721 coupling-polarity problem.\n- A changing hand position does not by itself select a non-premise polarity: fixed clock zero is convention, `theta(t)` is origin-relative, angular-velocity sign would require a new exported strict chiral/time-orientation source, and the P2721 two-polarity pair remains unselected.\n- Do not replay clock-position or clock-face-origin arguments as selector closure, pair12 strict-core, role transfer, bridge closure, `L_total`, or ToE.  A next admissible move must supply one actual strict chiral/time-orientation source law with a nonzero signed value coupled to the P2721 polarity pair, or pivot to a genuinely new typed object outside closed lanes; otherwise preserve the P2697-P2722 no-new-live-frontier certificate.\n")
    return payload


if __name__ == "__main__":
    main()
