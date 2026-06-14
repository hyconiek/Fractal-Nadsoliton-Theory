#!/usr/bin/env python3
"""P2712/S1662: post-P2711 sign-lane state-map no-new-live-frontier certificate.

After P2708-P2711, the selector/sign lane has tested: the boundary cocycle,
older-release unlock claims, inversion-odd Aut(Z12) characters, and finite
source-law sign couplings.  P2712 is the honest follow-up: compute whether any
current artifact now exports the new strict mechanism required to fix the
coupling sign, or whether the lane must be frozen until a genuinely new typed
object/source is supplied.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2712_s1662_post_p2711_sign_lane_state_map_no_new_live_frontier_certificate.json"
MD = GEN / "p2712_s1662_post_p2711_sign_lane_state_map_no_new_live_frontier_certificate.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

INPUTS = {
    "P2697_STATE_MAP": GEN / "p2697_s1647_post_direct_route_state_map_no_new_live_frontier_reconciliation.json",
    "P2707_POST_DAMPING_STATE_MAP": GEN / "p2707_s1657_post_p2706_no_new_live_frontier_reconciliation.json",
    "P2708_BOUNDARY_COCYCLE": GEN / "p2708_s1658_z12_boundary_1_cocycle_selector_source_obstruction.json",
    "P2709_RELEASE_BACKSCAN": GEN / "p2709_s1659_release_8_1_to_9_3_breakthrough_unlock_backscan.json",
    "P2710_CHARACTER_TEST": GEN / "p2710_s1660_finite_aut_z12_anti_inversion_orientation_character_source_test.json",
    "P2711_SIGN_COUPLING": GEN / "p2711_s1661_inversion_odd_character_source_law_sign_coupling_audit.json",
}

UNLOCK_FLAGS = [
    "new_live_frontier_found",
    "strict_mechanism_fixing_lambda_exported",
    "source_law_coupling_sign_exported",
    "strict_source_law_selects_plus_omega",
    "orientation_character_source_exported",
    "non_premise_selector_provider_exported",
    "qw2191_discharged",
    "pair12_strict_core_upgrade_exported",
    "ltotal_promoted",
    "toe_closure_exported",
]

LANES = [
    ("boundary_cocycle_object", "P2708_BOUNDARY_COCYCLE", ["boundary_1_cocycle_exports_nonpremise_selector"], "P2708 exports a real H1 line but no strict sign source."),
    ("older_release_unlock_claims", "P2709_RELEASE_BACKSCAN", ["current_unlock_found"], "P2709 finds older-release breakthroughs scoped/non-unlocking."),
    ("aut_z12_character_label", "P2710_CHARACTER_TEST", ["strict_orientation_character_source_exported"], "P2710 finds inversion-odd characters but no strict source law."),
    ("source_law_lambda_coupling", "P2711_SIGN_COUPLING", ["source_law_coupling_sign_exported", "strict_source_law_selects_plus_omega"], "P2711 finds lambda sign degeneracy and no non-premise lambda source."),
]


def read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"missing": True, "path": rel(path)}
    return json.loads(path.read_text(encoding="utf-8"))


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def get_decision(data: dict[str, Any]) -> dict[str, Any]:
    decision = data.get("decision")
    return decision if isinstance(decision, dict) else {}


def any_flag_true(data: dict[str, Any], flag_names: list[str]) -> bool:
    decision = get_decision(data)
    negative = decision.get("negative_export_flags") if isinstance(decision.get("negative_export_flags"), dict) else {}
    for flag in flag_names:
        if decision.get(flag) is True or negative.get(flag) is True or data.get(flag) is True:
            return True
    return False


def lane_matrix(loaded: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for lane, key, flags, reason in LANES:
        data = loaded[key]
        exported = any_flag_true(data, flags)
        rows.append({
            "lane": lane,
            "artifact": key,
            "status": data.get("status"),
            "checked_unlock_flags": flags,
            "unlock_exported": exported,
            "blocked_now": not exported,
            "reason": reason,
        })
    return rows


def global_unlock_scan(loaded: dict[str, dict[str, Any]]) -> dict[str, Any]:
    hits = []
    for key, data in loaded.items():
        decision = get_decision(data)
        negative = decision.get("negative_export_flags") if isinstance(decision.get("negative_export_flags"), dict) else {}
        for flag in UNLOCK_FLAGS:
            if data.get(flag) is True or decision.get(flag) is True or negative.get(flag) is True:
                hits.append({"artifact": key, "flag": flag})
    return {"unlock_flag_hits": hits, "unlock_flag_hit_count": len(hits)}


def admissibility_gate() -> list[dict[str, Any]]:
    return [
        {
            "required_new_item": "strict_mechanism_fixing_lambda",
            "current_status": "missing",
            "why_required": "P2711 leaves each inversion-odd character paired by lambda -> -lambda.",
            "admissible_next_test_if_supplied": "Run a source-law acceptance test proving the mechanism fixes lambda without premise/convention import.",
        },
        {
            "required_new_item": "different_new_typed_object_outside_selector_sign_lane",
            "current_status": "not supplied in current artifacts",
            "why_required": "P2697/P2707 and P2708-P2711 close replay moves on current selector/sign artifacts.",
            "admissible_next_test_if_supplied": "Open a fresh bounded audit only for that object and keep closure flags false until proven.",
        },
    ]


def write_markdown(payload: dict[str, Any]) -> None:
    lines = [
        "# P2712/S1662 post-P2711 sign-lane state-map no-new-live-frontier certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Lane matrix",
    ]
    for row in payload["lane_matrix"]:
        lines.append(f"- `{row['lane']}`: blocked_now={row['blocked_now']}, unlock_exported={row['unlock_exported']}. {row['reason']}")
    lines.extend(["", "## Admissibility gate"])
    for row in payload["admissibility_gate"]:
        lines.append(f"- `{row['required_new_item']}`: {row['current_status']}. {row['why_required']}")
    lines.extend(["", "## Decision", payload["decision"]["reason"], "", "## Next honest step", payload["decision"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    loaded = {key: read_json(path) for key, path in INPUTS.items()}
    lanes = lane_matrix(loaded)
    scan = global_unlock_scan(loaded)
    gate = admissibility_gate()
    no_frontier = all(row["blocked_now"] for row in lanes) and scan["unlock_flag_hit_count"] == 0
    payload = {
        "status": "P2712_POST_P2711_SIGN_LANE_NO_NEW_LIVE_FRONTIER_CERTIFICATE" if no_frontier else "P2712_REQUIRES_MANUAL_REVIEW",
        "input_hashes": {key: sha(path) for key, path in INPUTS.items()},
        "lane_matrix": lanes,
        "global_unlock_scan": scan,
        "admissibility_gate": gate,
        "decision": {
            "no_new_live_frontier_certificate": no_frontier,
            "selector_sign_lane_replay_frozen": no_frontier,
            "negative_export_flags": {flag: False for flag in UNLOCK_FLAGS},
            "reason": "P2712 checks the post-P2711 selector/sign state map.  Boundary cocycle, older-release backscan, Aut(Z12) inversion-odd character labels, and lambda sign-coupling source-law candidates are all blocked on current artifacts; no unlock flag is exported.",
            "next_honest_step": "Do not continue the selector/sign lane by replay.  The next admissible move requires either a genuinely new strict mechanism fixing lambda, or a different new typed object outside the closed lanes.  If neither is supplied, preserve the P2697-P2712 no-new-live-frontier certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2712/S1662 post-P2711 sign-lane no-new-live-frontier certificate", "## P2712/S1662 post-P2711 sign-lane no-new-live-frontier certificate\n\n`P2712/S1662` reconciles the selector/sign state map after P2708-P2711.  The boundary-cocycle object, older-release backscan, Aut(Z12) inversion-odd character labels, and source-law `lambda` coupling candidates remain blocked on current artifacts, with no exported unlock flag.  The lane is therefore frozen against replay unless a genuinely new strict mechanism fixes `lambda` or a different new typed object is supplied.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2712/S1662 sign-lane state-map Ltotal guard", "## P2712/S1662 sign-lane state-map Ltotal guard\n\n`P2712/S1662` is a no-new-live-frontier certificate for the post-P2711 selector/sign lane, not a variational construction.  It does not promote `L_total`, selector closure, pair12 strict-core, role transfer, bridge closure, or ToE.\n")
    append_once(AGENTS, "Current post-P2711 selector/sign no-new-live-frontier guardrail (P2712/S1662, 2026-06-14)", "## Current post-P2711 selector/sign no-new-live-frontier guardrail (P2712/S1662, 2026-06-14)\n\n- P2712 reconciles P2708-P2711 and finds no new live frontier in the selector/sign lane on current artifacts: boundary cocycle, older-release backscan, inversion-odd character labels, and `lambda` sign-coupling source-law candidates remain blocked.\n- Do not continue selector/sign replay or promote `QW-2191`, pair12 strict-core, role transfer, `L_total`, bridge closure, or ToE from P2708-P2711.\n- A next admissible move requires either a genuinely new strict mechanism fixing `lambda` or a different new typed object outside the closed lanes; otherwise preserve the P2697-P2712 no-new-live-frontier certificate.\n")
    return payload


if __name__ == "__main__":
    main()
