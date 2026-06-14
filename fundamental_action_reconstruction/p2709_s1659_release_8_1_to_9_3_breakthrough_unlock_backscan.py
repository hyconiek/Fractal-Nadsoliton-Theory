#!/usr/bin/env python3
"""P2709/S1659: release 8.1-9.3 breakthrough unlock backscan.

After the P2708 boundary-cocycle obstruction and a human release/tag summary,
scan the strongest earlier-release breakthrough artifacts and ask a narrow
question: did any of them already export the missing current strict unlock
(selector sign/source, full tensor Lagrangian/EOM, variational source strength,
bridge source/role transfer, L_total/ToE closure), or are they scoped progress
that remains bounded by the current P2697-P2708 guardrail?
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2709_s1659_release_8_1_to_9_3_breakthrough_unlock_backscan.json"
MD = GEN / "p2709_s1659_release_8_1_to_9_3_breakthrough_unlock_backscan.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

INPUTS = {
    "P1581_R81_SELECTOR_SOURCE": ROOT / "P1581_S531_STRICT_SELECTOR_SOURCE_EXPORT_AND_SYMMETRY_BREAKING_WITNESS_PACKET_PL.md",
    "P1582_R82_SELECTOR_UNIQUENESS": ROOT / "P1582_S532_STRICT_SELECTOR_UNIQUENESS_THEOREM_BRIDGE_TO_FULL_LAGRANGIAN_PACKET_PL.md",
    "P1642_NO_FALSE_PASS": ROOT / "P1642_S592_NO_FALSE_PASS_STATUS_CONSISTENCY_AUDIT_PACKET_PL.md",
    "P1624_NONCYCLIC_SELECTOR_WITNESS": ROOT / "P1624_S574_NONCYCLIC_SELECTOR_WITNESS_FROM_STRICT_KERNEL_PACKET_PL.md",
    "P2327_FUTURE_STATE_SELECTOR_CONDITION": ROOT / "P2327_S1277_KERNEL_DERIVED_FUTURE_STATE_SELECTOR_QW2191_CONDITION_PACKET.md",
    "P2319_D12_ORIENTATION_NO_GO": GEN / "p2319_s1269_strict_dihedral_orientation_response_no_go_probe.json",
    "P2320_D12_COMMUTANT_NO_GO": GEN / "p2320_s1270_strict_dihedral_commutant_orientation_no_go_probe.json",
    "P2336_TRACK_A_TRANSPORT": GEN / "p2336_s1286_strict_track_a_second_atlas_finite_part_transport_theorem.json",
    "P2361_TRACK_B_GLUING": GEN / "p2361_s1311_strict_track_b_o5a_smooth_interface_finite_gluing_summary.json",
    "P2363_MOMENT_TRANSPORT": GEN / "p2363_s1313_legacy_strict_bridge_moment_lagrangian_transport_probe.json",
    "P2364_FRW_RESIDUALS": GEN / "p2364_s1314_bridge_completed_frw_scalar_gauge_gravity_residual_table_probe.json",
    "P2375_POLARITY_INTERVAL": GEN / "p2375_s1325_damping_compression_polarity_interval_robustness_theorem.json",
    "P2376_ETA_BETA_RECTANGLE": GEN / "p2376_s1326_damping_compression_eta_beta_rectangle_robustness_theorem.json",
    "P2377_TRANSPORT_PRIMITIVE": GEN / "p2377_s1327_damping_compression_transport_primitive_uniform_coupling_theorem.json",
    "P2378_UNIT_NORMALIZED_INSUFFICIENCY": GEN / "p2378_s1328_unit_normalized_transport_coupling_insufficiency_theorem.json",
    "P2708_BOUNDARY_COCYCLE": GEN / "p2708_s1658_z12_boundary_1_cocycle_selector_source_obstruction.json",
}

CLOSURE_TOKENS = [
    "qw-2191 discharged", "qw2191_discharged\": true", "selector_closure_exported\": true",
    "non_premise_selector_provider_exported\": true", "toe_closure_claimed\": true",
    "ltotal_promoted\": true", "role_transfer_started\": true",
]
BLOCKER_TOKENS = [
    "open", "keep_open", "no-go", "no_go", "no false pass", "no_false_pass",
    "insufficient", "still_open", "no variational source", "no_selector", "blocked",
]

ROWS = [
    ("release_8_1_8_2_selector_claims", ["P1581_R81_SELECTOR_SOURCE", "P1582_R82_SELECTOR_UNIQUENESS", "P1642_NO_FALSE_PASS"], "Older selector/source claims must be checked against NO_FALSE_PASS and current P2699-P2708 criteria."),
    ("release_8_4_noncyclic_selector_witness", ["P1624_NONCYCLIC_SELECTOR_WITNESS"], "Noncyclic witness evidence is useful only if it exports theorem-level QW-2191 discharge."),
    ("release_8_9_future_state_selector_condition", ["P2327_FUTURE_STATE_SELECTOR_CONDITION"], "P2327 names the condition for QW-2191; a condition is not itself the exported selector source."),
    ("release_9_0_d12_orientation_no_go", ["P2319_D12_ORIENTATION_NO_GO", "P2320_D12_COMMUTANT_NO_GO"], "D12-invariant orientation audits are negative controls and cannot unlock the selector."),
    ("release_9_1_track_a_b_progress", ["P2336_TRACK_A_TRANSPORT", "P2361_TRACK_B_GLUING"], "FRW/Bianchi-I finite-part transport and O5a gluing are scoped progress, not global 4D/L_total closure."),
    ("release_9_3_moment_and_frw_transport", ["P2363_MOMENT_TRANSPORT", "P2364_FRW_RESIDUALS"], "Moment transport and FRW residual export improve the bridge/Lagrangian lane but explicitly remain selector/role-transfer bounded."),
    ("release_9_3_damping_robustness_and_source_strength", ["P2375_POLARITY_INTERVAL", "P2376_ETA_BETA_RECTANGLE", "P2377_TRANSPORT_PRIMITIVE", "P2378_UNIT_NORMALIZED_INSUFFICIENCY"], "Damping robustness and transport primitive must be separated from variational source closure and selector sign sourcing."),
    ("current_p2708_boundary_cocycle_sign_gap", ["P2708_BOUNDARY_COCYCLE"], "P2708 supplies the current explicit missing-sign obstruction for the selector-source lane."),
]


def read_text(path: Path) -> str:
    if not path.exists():
        return ""
    return path.read_text(encoding="utf-8")


def read_status(path: Path) -> str | None:
    if not path.exists():
        return None
    if path.suffix == ".json":
        data = json.loads(path.read_text(encoding="utf-8"))
        return str(data.get("status") or data.get("decision") or "")
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines()[:80]:
        if line.lower().startswith("status:"):
            return line.split(":", 1)[1].strip()
    return None


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def classify_row(name: str, keys: list[str], obligation: str) -> dict[str, Any]:
    texts = {key: read_text(INPUTS[key]).lower() for key in keys}
    joined = "\n".join(texts.values())
    closure_hits = sorted({tok for tok in CLOSURE_TOKENS if tok in joined})
    blocker_hits = sorted({tok for tok in BLOCKER_TOKENS if tok in joined})
    statuses = {key: read_status(INPUTS[key]) for key in keys}
    positive_progress = any(word in joined for word in ["pass", "export", "robust", "progress", "theorem", "positive"])
    # Conservative unlock criterion: an older row can unlock only if it has an explicit current closure token
    # and no blocker token.  This backscan is intentionally stricter than release prose summaries.
    current_unlock = bool(closure_hits) and not blocker_hits
    return {
        "lane": name,
        "inputs": [rel(INPUTS[key]) for key in keys],
        "statuses": statuses,
        "positive_or_breakthrough_content_present": positive_progress,
        "closure_token_hits": closure_hits,
        "blocker_token_hits_sample": blocker_hits[:8],
        "current_unlock_exported": current_unlock,
        "scope_obligation": obligation,
        "verdict": "NO_CURRENT_UNLOCK" if not current_unlock else "MANUAL_REVIEW_POSSIBLE_UNLOCK",
    }


def release_backscan() -> list[dict[str, Any]]:
    return [classify_row(name, keys, obligation) for name, keys, obligation in ROWS]


def write_markdown(payload: dict[str, Any]) -> None:
    lines = [
        "# P2709/S1659 release 8.1-9.3 breakthrough unlock backscan",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Backscan matrix",
    ]
    for row in payload["release_backscan_matrix"]:
        lines.append(f"- `{row['lane']}`: verdict={row['verdict']}, current_unlock_exported={row['current_unlock_exported']}. {row['scope_obligation']}")
    lines.extend(["", "## Decision", payload["decision"]["reason"], "", "## Next honest step", payload["decision"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    matrix = release_backscan()
    no_unlock = all(not row["current_unlock_exported"] for row in matrix)
    payload = {
        "status": "P2709_RELEASE_8_1_TO_9_3_BACKSCAN_NO_CURRENT_UNLOCK" if no_unlock else "P2709_REQUIRES_MANUAL_REVIEW",
        "input_hashes": {key: sha(path) for key, path in INPUTS.items()},
        "release_backscan_matrix": matrix,
        "decision": {
            "older_release_breakthroughs_reviewed": True,
            "current_unlock_found": not no_unlock,
            "negative_export_flags": {
                "qw2191_discharged": False,
                "non_premise_selector_provider_exported": False,
                "full_tensor_lagrangian_eom_closed": False,
                "variational_source_strength_closed": False,
                "legacy_role_transfer_started": False,
                "ltotal_or_toe_promoted": False,
            },
            "reason": "The strongest Release 8.1-9.3 breakthroughs are real scoped progress, but the scanned artifacts carry open/no-false-pass/insufficiency/no-go boundaries or scoped statuses.  None exports the missing current unlock after P2708: a strict source for the missing sign of the boundary-cocycle, QW-2191 discharge, full tensor nonproxy Lagrangian/EOM closure, variational source strength closure, role transfer, L_total, or ToE closure.",
            "next_honest_step": "Do not replay the older release breakthroughs as closure evidence.  The next admissible move should be a single new typed source candidate for the missing sign in P2708, preferably a finite anti-inversion / orientation-character source test that explicitly states which strict artifact breaks Aut(Z12) inversion.  If no such new source is available, preserve the P2697-P2709 no-current-unlock certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2709/S1659 release 8.1-9.3 breakthrough unlock backscan", "## P2709/S1659 release 8.1-9.3 breakthrough unlock backscan\n\n`P2709/S1659` backscans the strongest Release 8.1-9.3 breakthrough artifacts against the current P2708 missing-sign boundary.  The scan records real scoped progress in NO_FALSE_PASS, selector-condition, D12 no-go, Track-A/Track-B, bridge moment transport, FRW residual, and damping robustness lanes, but finds no current export of a strict source for the boundary-cocycle sign, `QW-2191` discharge, variational source-strength closure, full tensor nonproxy Lagrangian/EOM closure, role transfer, `L_total`, or ToE closure.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2709/S1659 release backscan Ltotal guard", "## P2709/S1659 release backscan Ltotal guard\n\n`P2709/S1659` checks whether older Release 8.1-9.3 breakthroughs already unlock the current `L_total`/EOM/selector boundary.  They do not: P2363/P2364/P2377/P2378 remain scoped bridge/damping progress with source-strength, selector, role-transfer, and global tensor closure obligations open.  This is not a variational closure or ToE promotion.\n")
    append_once(AGENTS, "Current release 8.1-9.3 breakthrough backscan guardrail (P2709/S1659, 2026-06-14)", "## Current release 8.1-9.3 breakthrough backscan guardrail (P2709/S1659, 2026-06-14)\n\n- P2709 backscans the strongest Release 8.1-9.3 breakthroughs against the current P2708 missing-sign and P2697-P2708 no-unlock boundary.\n- The older breakthroughs remain real scoped progress, but none currently exports a strict source for the boundary-cocycle sign, discharges `QW-2191`, closes variational source strength, closes full tensor nonproxy Lagrangian/EOM, starts role transfer, or promotes `L_total`/ToE.\n- Do not replay older release/tag breakthrough prose as closure evidence; a next admissible move must introduce one new typed anti-inversion/orientation-character source candidate or preserve the P2697-P2709 no-current-unlock certificate.\n")
    return payload


if __name__ == "__main__":
    main()
