#!/usr/bin/env python3
"""P2724/S1674: post-P2723 commit intake and next honest step certificate.

The user supplied commit 88eb860... and asked for the next honest research
step.  The commit already contains P2719-P2723 and ends with the requirement:
either introduce a new strict dynamic/chiral source artifact with computable
signed value coupled to P2721, or preserve the no-new-live-frontier certificate.
This packet performs that intake without manufacturing closure.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2724_s1674_post_p2723_commit_intake_next_honest_step_certificate.json"
MD = GEN / "p2724_s1674_post_p2723_commit_intake_next_honest_step_certificate.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

COMMIT = "88eb860b1658ac5b648253fa65dd83bd4abbe922"
INPUTS = {
    "P2723_STRICT_CHIRAL_TIME_SOURCE_MATRIX": GEN / "p2723_s1673_strict_chiral_time_orientation_source_law_candidate_matrix.json",
    "P2722_CLOCK_PHASE_SOURCE_LAW_AUDIT": GEN / "p2722_s1672_clock_phase_polarity_selection_source_law_audit.json",
    "P2721_SIGN_TORSOR_COUPLING_AUDIT": GEN / "p2721_s1671_chiral_bispectrum_sign_to_orientation_torsor_coupling_audit.json",
    "P2718_CHIRAL_BISPECTRUM_MARKER": GEN / "p2718_s1668_chiral_bispectrum_explicit_signed_formula_torsor_coupling_audit.json",
}

CLOSED_OR_BLOCKED_LANES = [
    "selector/QW-2191 replay",
    "damping-to-selector replay",
    "older-release transfer",
    "direct-route residual replay",
    "P2680 bridge-source atom replay",
    "Lagrangian/EOM reverse-closure replay",
    "lower-boundary recursion",
    "role transfer",
    "L_total promotion",
    "ToE closure",
    "generic clock/time-arrow naming",
    "generic pseudoscalar/chiral name enumeration",
]

REQUIRED_NEW_OBJECT = {
    "name": "strict_dynamic_chiral_source_artifact_coupled_to_P2721",
    "requirements": [
        "exported internal strict source law",
        "computable nonzero signed value",
        "nonconventional sign, not external time-arrow or clock-origin premise",
        "explicit coupling to the P2721 two-polarity sign-to-orientation torsor pair",
        "polarity-selection theorem choosing exactly one coupling polarity",
    ],
}

NEGATIVE_EXPORT_FLAGS = [
    "new_strict_dynamic_chiral_source_artifact_supplied_by_commit_intake",
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


def commit_summary() -> dict[str, Any]:
    try:
        raw = subprocess.check_output(
            ["git", "show", "--first-parent", "--name-only", "--oneline", "--no-renames", COMMIT],
            cwd=REPO,
            text=True,
        )
    except subprocess.CalledProcessError as exc:
        return {"available": False, "error": str(exc)}
    lines = raw.splitlines()
    return {
        "available": True,
        "first_line": lines[0] if lines else "",
        "mentions_p2719_to_p2723": all(token in raw.lower() for token in ["p2719", "p2720", "p2721", "p2722", "p2723"]),
        "raw_stat_line_count": len(lines),
    }


def evaluate_intake(inputs: dict[str, dict[str, Any]]) -> dict[str, Any]:
    p2723 = inputs["P2723_STRICT_CHIRAL_TIME_SOURCE_MATRIX"]
    p2723_no_source = p2723.get("status") == "P2723_STRICT_CHIRAL_TIME_ORIENTATION_SOURCE_LAW_MATRIX_NO_ACCEPTED_SOURCE"
    new_artifact_supplied = False
    return {
        "p2723_no_source_boundary_confirmed": p2723_no_source,
        "new_strict_dynamic_chiral_source_artifact_supplied_by_commit_intake": new_artifact_supplied,
        "required_new_object": REQUIRED_NEW_OBJECT,
        "closed_or_blocked_lanes": CLOSED_OR_BLOCKED_LANES,
        "admissible_next_move": "Supply exactly one explicit strict dynamic/chiral source artifact satisfying all five listed requirements, then run a bounded acceptance/witness test; absent that artifact, preserve the P2697-P2724 no-new-live-frontier certificate.",
        "recommended_next_honest_step": "Do not open another replay lane.  Draft or locate one concrete internal law for a signed dynamic/chiral quantity and a theorem coupling its sign to the P2721 polarity pair; if no such law is supplied, stop at the no-new-live-frontier certificate.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lines = [
        "# P2724/S1674 post-P2723 commit intake and next honest step certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Supplied commit intake",
        f"- commit={COMMIT}",
        f"- summary={payload['commit_summary']['first_line']}",
        f"- mentions P2719-P2723={payload['commit_summary']['mentions_p2719_to_p2723']}",
        "",
        "## Decision",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["intake"]["recommended_next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    inputs = {key: read_json(path) for key, path in INPUTS.items()}
    intake = evaluate_intake(inputs)
    no_unlock = not intake["new_strict_dynamic_chiral_source_artifact_supplied_by_commit_intake"]
    payload = {
        "status": "P2724_POST_P2723_COMMIT_INTAKE_NO_NEW_LIVE_FRONTIER_CERTIFICATE" if no_unlock else "P2724_REQUIRES_MANUAL_REVIEW",
        "commit": COMMIT,
        "commit_summary": commit_summary(),
        "input_hashes": {key: sha(path) for key, path in INPUTS.items()},
        "input_statuses": {key: value.get("status") for key, value in inputs.items()},
        "intake": intake,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": "The supplied commit already reaches P2723: current clock/chiral/time-orientation candidates do not export a strict source law with a nonzero signed value coupled to the P2721 polarity pair.  No additional strict dynamic/chiral source artifact is supplied by the intake itself, so the honest continuation is a no-new-live-frontier certificate plus a precise requirement for the next admissible object.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2724/S1674 post-P2723 commit intake", "## P2724/S1674 post-P2723 commit intake\n\n`P2724/S1674` ingests commit `88eb860b1658ac5b648253fa65dd83bd4abbe922` and preserves the P2723 boundary: no new strict dynamic/chiral source artifact with computable signed value coupled to the P2721 polarity pair is supplied.  No `QW-2191` discharge, pair12 strict-core upgrade, bridge closure, role transfer, `L_total`, or ToE closure follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2724/S1674 no new Ltotal source", "## P2724/S1674 no new Ltotal source\n\n`P2724/S1674` is an intake/recommendation certificate.  Because it supplies no variational source term or signed dynamic/chiral law, it does not promote `L_total` and does not reopen the frozen Lagrangian/EOM lane.\n")
    append_once(AGENTS, "Current post-P2723 commit-intake no-new-live-frontier guardrail (P2724/S1674, 2026-06-14)", "## Current post-P2723 commit-intake no-new-live-frontier guardrail (P2724/S1674, 2026-06-14)\n\n- P2724 ingests commit `88eb860b1658ac5b648253fa65dd83bd4abbe922`, which already contains the P2719-P2723 localizer/coupling/clock/chiral-time audits.\n- The intake supplies no new strict dynamic/chiral source artifact with an exported internal law, computable nonzero signed value, nonconventional sign, coupling to the P2721 polarity pair, and exact polarity-selection theorem.\n- Preserve the P2697-P2724 no-new-live-frontier certificate: do not reopen selector/sign replay, damping-to-selector replay, older-release transfer, direct-route residual replay, bridge-source replay, Lagrangian/EOM reverse closure, role transfer, `L_total`, or ToE.  A next admissible move must supply exactly one such strict dynamic/chiral source artifact and run only a bounded acceptance/witness test.\n")
    return payload


if __name__ == "__main__":
    main()
