#!/usr/bin/env python3
"""P2713/S1663: post-P2712 new-typed-object intake gate certificate.

P2712 froze the selector/sign lane unless a genuinely new strict mechanism
fixing lambda, or a different new typed object outside the closed lanes, is
supplied.  P2713 is deliberately not another replay audit.  It formalizes the
intake gate for the next research move and checks the current session's
candidate list.  With no new object supplied, the only honest result is to
preserve the P2697-P2712 no-new-live-frontier certificate.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2713_s1663_post_p2712_new_typed_object_intake_gate_certificate.json"
MD = GEN / "p2713_s1663_post_p2712_new_typed_object_intake_gate_certificate.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

INPUTS = {
    "P2697_STATE_MAP": GEN / "p2697_s1647_post_direct_route_state_map_no_new_live_frontier_reconciliation.json",
    "P2707_POST_DAMPING_STATE_MAP": GEN / "p2707_s1657_post_p2706_no_new_live_frontier_reconciliation.json",
    "P2712_SIGN_LANE_STATE_MAP": GEN / "p2712_s1662_post_p2711_sign_lane_state_map_no_new_live_frontier_certificate.json",
}

CLOSED_REPLAY_LANES = [
    "selector_sign_replay_or_qw2191_discharge_without_new_source",
    "damping_to_selector_replay_from_p2377_p2378",
    "older_release_transfer_from_p1343_p1348_or_release_prose",
    "direct_route_residual_replay_without_new_provider",
    "p2680_bridge_source_atom_replay",
    "lagrangian_eom_reverse_closure_without_new_anisotropic_source",
    "lower_boundary_recursion_without_new_seed_or_provider",
    "role_transfer_before_bridge_source_closure",
    "ltotal_or_toe_promotion",
]

# No new typed object was supplied in the user's request beyond the current
# P2708-P2712 sign-lane artifacts.  Keeping this as data makes the gate explicit
# and prevents the script from manufacturing a pseudo-candidate.
SUPPLIED_CANDIDATES: list[dict[str, Any]] = []

ACCEPTANCE_CRITERIA = [
    {
        "criterion": "strict_source_or_typed_object_is_new",
        "required": True,
        "description": "The object must not be one of the replayed P2708-P2712 boundary cocycle, character-label, or lambda-pair candidates.",
    },
    {
        "criterion": "outside_closed_lane_or_new_blocker_cut",
        "required": True,
        "description": "It must either be outside the closed lanes or introduce a genuinely new blocker-cut/provider class.",
    },
    {
        "criterion": "non_premise_non_convention_export",
        "required": True,
        "description": "It must be exported by a strict artifact, not assumed as a branch convention, premise selector, or older release prose.",
    },
    {
        "criterion": "no_closure_promotion_before_witness",
        "required": True,
        "description": "Even an admitted candidate only unlocks a bounded test; it does not promote QW-2191, pair12 strict-core, bridge, role-transfer, L_total, or ToE.",
    },
]

NEGATIVE_EXPORT_FLAGS = {
    "new_live_frontier_found": False,
    "strict_mechanism_fixing_lambda_exported": False,
    "different_new_typed_object_exported": False,
    "non_premise_selector_provider_exported": False,
    "qw2191_discharged": False,
    "pair12_strict_core_upgrade_exported": False,
    "bridge_closure_exported": False,
    "role_transfer_started": False,
    "ltotal_promoted": False,
    "toe_closure_exported": False,
}


def read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"missing": True, "path": rel(path)}
    return json.loads(path.read_text(encoding="utf-8"))


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def prior_certificates_hold(loaded: dict[str, dict[str, Any]]) -> bool:
    p2712 = loaded["P2712_SIGN_LANE_STATE_MAP"].get("decision", {})
    return bool(p2712.get("no_new_live_frontier_certificate")) and bool(p2712.get("selector_sign_lane_replay_frozen"))


def candidate_intake_rows() -> list[dict[str, Any]]:
    if not SUPPLIED_CANDIDATES:
        return [{
            "candidate_id": "NO_NEW_TYPED_OBJECT_SUPPLIED",
            "admitted_for_next_test": False,
            "reason": "The current request supplies no new strict mechanism fixing lambda and no different new typed object outside the closed lanes.",
            "criteria_checked": [row["criterion"] for row in ACCEPTANCE_CRITERIA],
        }]
    rows = []
    for candidate in SUPPLIED_CANDIDATES:
        checks = candidate.get("checks", {})
        admitted = all(checks.get(row["criterion"]) is True for row in ACCEPTANCE_CRITERIA)
        rows.append({
            "candidate_id": candidate.get("candidate_id", "unnamed_candidate"),
            "admitted_for_next_test": admitted,
            "reason": "candidate satisfies all intake criteria" if admitted else "candidate fails at least one intake criterion",
            "criteria_checked": checks,
        })
    return rows


def write_markdown(payload: dict[str, Any]) -> None:
    lines = [
        "# P2713/S1663 post-P2712 new-typed-object intake gate certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Closed replay lanes",
    ]
    for lane in payload["closed_replay_lanes"]:
        lines.append(f"- `{lane}`")
    lines.extend(["", "## Candidate intake"])
    for row in payload["candidate_intake"]:
        lines.append(f"- `{row['candidate_id']}`: admitted={row['admitted_for_next_test']}. {row['reason']}")
    lines.extend(["", "## Decision", payload["decision"]["reason"], "", "## Next honest step", payload["decision"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    loaded = {key: read_json(path) for key, path in INPUTS.items()}
    prior_hold = prior_certificates_hold(loaded)
    intake = candidate_intake_rows()
    admitted = [row for row in intake if row["admitted_for_next_test"]]
    no_new_frontier = prior_hold and not admitted
    payload = {
        "status": "P2713_NEW_TYPED_OBJECT_INTAKE_GATE_NO_CANDIDATE_CERTIFICATE" if no_new_frontier else "P2713_CANDIDATE_REQUIRES_BOUNDED_TEST",
        "input_hashes": {key: sha(path) for key, path in INPUTS.items()},
        "prior_certificates_hold": prior_hold,
        "closed_replay_lanes": CLOSED_REPLAY_LANES,
        "acceptance_criteria": ACCEPTANCE_CRITERIA,
        "candidate_intake": intake,
        "admitted_candidate_count": len(admitted),
        "decision": {
            "no_new_live_frontier_certificate_preserved": no_new_frontier,
            "replay_lanes_frozen": no_new_frontier,
            "negative_export_flags": NEGATIVE_EXPORT_FLAGS,
            "reason": "P2713 applies the post-P2712 intake gate rather than replaying closed lanes.  No new strict mechanism fixing lambda and no different new typed object were supplied, so the P2697-P2712 no-new-live-frontier certificate is preserved.",
            "next_honest_step": "Supply exactly one genuinely new strict typed object/source/mechanism outside the closed replay lanes, or a strict mechanism fixing lambda, and then run only a bounded acceptance/witness test for that object.  If no such object is supplied, stop at the P2697-P2713 no-new-live-frontier certificate and do not manufacture closure.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2713/S1663 post-P2712 new-typed-object intake gate certificate", "## P2713/S1663 post-P2712 new-typed-object intake gate certificate\n\n`P2713/S1663` applies the post-P2712 intake gate for any next research move.  With no new strict mechanism fixing `lambda` and no different genuinely new typed object supplied, it preserves the `P2697-P2712` no-new-live-frontier certificate.  Closed replay lanes remain frozen, and any future candidate may only unlock a bounded acceptance/witness test, not `QW-2191`, pair12 strict-core, bridge closure, role transfer, `L_total`, or ToE promotion.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2713/S1663 intake gate Ltotal guard", "## P2713/S1663 intake gate Ltotal guard\n\n`P2713/S1663` is an intake-gate/no-candidate certificate, not a variational source construction.  It preserves the no-new-live-frontier state unless a genuinely new strict typed object or strict mechanism fixing `lambda` is supplied, and it does not promote `L_total`, selector closure, pair12 strict-core, role transfer, bridge closure, or ToE.\n")
    append_once(AGENTS, "Current post-P2712 new-typed-object intake gate guardrail (P2713/S1663, 2026-06-14)", "## Current post-P2712 new-typed-object intake gate guardrail (P2713/S1663, 2026-06-14)\n\n- P2713 applies the post-P2712 admissibility gate rather than replaying closed lanes: no new strict mechanism fixing `lambda` and no different genuinely new typed object were supplied in the current artifacts.\n- Preserve the P2697-P2713 no-new-live-frontier certificate; selector/sign replay, damping-to-selector replay, older-release transfer, direct-route residual replay, bridge-source replay, Lagrangian/EOM reverse closure, lower-boundary recursion, role transfer, `L_total`, and ToE promotion remain blocked.\n- A next admissible move must supply exactly one genuinely new strict typed object/source/mechanism or a strict mechanism fixing `lambda`, then run only a bounded acceptance/witness test; otherwise do not manufacture closure.\n")
    return payload


if __name__ == "__main__":
    main()
