#!/usr/bin/env python3
"""P2714/S1664: Z12 orientation-torsor global-section obstruction.

P2713 made the intake rule explicit: do not replay closed lanes unless one new
strict typed object/source/mechanism is supplied.  P2714 supplies and tests one
narrow typed candidate: the orientation torsor of the P2708 boundary-cocycle
line under Aut(Z12).  If this torsor had an Aut-compatible global section, it
could fix the missing lambda/sign without importing a premise selector.  The
finite calculation below checks that exact condition.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2714_s1664_z12_orientation_torsor_global_section_obstruction.json"
MD = GEN / "p2714_s1664_z12_orientation_torsor_global_section_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

UNITS = [1, 5, 7, 11]
ORIENTATIONS = [-1, 1]

INPUTS = {
    "P2708_BOUNDARY_COCYCLE": GEN / "p2708_s1658_z12_boundary_1_cocycle_selector_source_obstruction.json",
    "P2710_CHARACTER_TABLE": GEN / "p2710_s1660_finite_aut_z12_anti_inversion_orientation_character_source_test.json",
    "P2711_SIGN_COUPLING": GEN / "p2711_s1661_inversion_odd_character_source_law_sign_coupling_audit.json",
    "P2712_SIGN_LANE_STATE_MAP": GEN / "p2712_s1662_post_p2711_sign_lane_state_map_no_new_live_frontier_certificate.json",
    "P2713_INTAKE_GATE": GEN / "p2713_s1663_post_p2712_new_typed_object_intake_gate_certificate.json",
}

NEGATIVE_EXPORT_FLAGS = [
    "orientation_torsor_global_section_exported",
    "strict_mechanism_fixing_lambda_exported",
    "non_premise_selector_provider_exported",
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


def aut_orientation_action(unit: int, orientation: int) -> int:
    """Aut(Z12) acts on the cycle orientation by unit parity around the 12-cycle.

    Units 1 and 5 preserve the chosen generator orientation; units 7 and 11
    reverse it.  In particular, inversion 11 sends omega to -omega.
    """
    if unit in (1, 5):
        return orientation
    if unit in (7, 11):
        return -orientation
    raise ValueError(f"not a Z12 unit: {unit}")


def action_table() -> list[dict[str, Any]]:
    rows = []
    for unit in UNITS:
        image_plus = aut_orientation_action(unit, 1)
        rows.append({
            "unit": unit,
            "image_of_plus_omega": "+omega" if image_plus == 1 else "-omega",
            "orientation_preserving": image_plus == 1,
            "orientation_reversing": image_plus == -1,
        })
    return rows


def global_section_candidates() -> list[dict[str, Any]]:
    rows = []
    for chosen in ORIENTATIONS:
        failures = []
        for unit in UNITS:
            image = aut_orientation_action(unit, chosen)
            if image != chosen:
                failures.append({
                    "unit": unit,
                    "chosen_section": "+omega" if chosen == 1 else "-omega",
                    "image": "+omega" if image == 1 else "-omega",
                    "failure": "Aut action does not preserve the proposed section.",
                })
        rows.append({
            "candidate_section": "+omega" if chosen == 1 else "-omega",
            "aut_compatible": not failures,
            "failure_count": len(failures),
            "failures": failures,
        })
    return rows


def torsor_summary(rows: list[dict[str, Any]]) -> dict[str, Any]:
    admitted = [row for row in rows if row["aut_compatible"]]
    return {
        "typed_object": "orientation_torsor_of_P2708_boundary_cocycle_line",
        "torsor_fiber": ["+omega", "-omega"],
        "aut_group": "Aut(Z12)=U(12)={1,5,7,11}",
        "global_section_count": len(admitted),
        "has_aut_compatible_global_section": bool(admitted),
        "obstruction": "orientation-reversing automorphisms 7 and 11 exchange the two torsor points",
    }


def prior_boundary_check() -> dict[str, Any]:
    loaded = {key: read_json(path) for key, path in INPUTS.items()}
    p2713_decision = loaded["P2713_INTAKE_GATE"].get("decision", {})
    return {
        "p2713_no_new_frontier_preserved_before_this_candidate": bool(p2713_decision.get("no_new_live_frontier_certificate_preserved")),
        "p2713_replay_lanes_frozen": bool(p2713_decision.get("replay_lanes_frozen")),
        "input_statuses": {key: value.get("status") for key, value in loaded.items()},
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lines = [
        "# P2714/S1664 Z12 orientation-torsor global-section obstruction",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## New typed candidate",
        f"- `{payload['torsor_summary']['typed_object']}` with fiber `{payload['torsor_summary']['torsor_fiber']}`",
        "",
        "## Aut action",
    ]
    for row in payload["action_table"]:
        lines.append(f"- unit `{row['unit']}` maps `+omega` to `{row['image_of_plus_omega']}`; reversing={row['orientation_reversing']}")
    lines.extend(["", "## Global section candidates"])
    for row in payload["global_section_candidates"]:
        lines.append(f"- `{row['candidate_section']}`: aut_compatible={row['aut_compatible']}, failure_count={row['failure_count']}")
    lines.extend(["", "## Decision", payload["decision"]["reason"], "", "## Next honest step", payload["decision"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    actions = action_table()
    sections = global_section_candidates()
    summary = torsor_summary(sections)
    no_unlock = not summary["has_aut_compatible_global_section"]
    payload = {
        "status": "P2714_ORIENTATION_TORSOR_GLOBAL_SECTION_OBSTRUCTION_NO_STRICT_LAMBDA_FIX" if no_unlock else "P2714_REQUIRES_MANUAL_REVIEW",
        "input_hashes": {key: sha(path) for key, path in INPUTS.items()},
        "prior_boundary_check": prior_boundary_check(),
        "action_table": actions,
        "global_section_candidates": sections,
        "torsor_summary": summary,
        "decision": {
            "new_typed_candidate_tested": True,
            "orientation_torsor_global_section_exported": False,
            "strict_mechanism_fixing_lambda_exported": False,
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": "P2714 tests one new typed candidate admitted after P2713: the orientation torsor of the P2708 boundary-cocycle line.  The finite Aut(Z12) action has orientation-reversing units 7 and 11, so neither +omega nor -omega is an Aut-compatible global section.  The torsor is real, but it remains a two-point premise-sign torsor rather than a strict mechanism fixing lambda.",
            "next_honest_step": "Do not replay the orientation-torsor, character, lambda-coupling, damping-to-selector, older-release, or direct-route lanes.  A next admissible move must supply a new strict source that breaks the orientation torsor by an internal non-premise law, or a different genuinely new typed object outside the closed lanes.  Without that, preserve the P2697-P2714 no-new-live-frontier certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2714/S1664 Z12 orientation-torsor global-section obstruction", "## P2714/S1664 Z12 orientation-torsor global-section obstruction\n\n`P2714/S1664` supplies one new typed candidate after the P2713 intake gate: the orientation torsor of the P2708 boundary-cocycle line.  The finite `Aut(Z12)` action sends `+omega` to `-omega` for orientation-reversing units `7` and `11`, so neither torsor point is an Aut-compatible global section.  This is a real obstruction object, but it exports no strict mechanism fixing `lambda` and does not discharge `QW-2191`, upgrade pair12 strict-core, close the bridge, start role transfer, promote `L_total`, or imply ToE.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2714/S1664 orientation-torsor Ltotal guard", "## P2714/S1664 orientation-torsor Ltotal guard\n\n`P2714/S1664` is a finite orientation-torsor/global-section obstruction, not a variational source construction.  Because no Aut-compatible global section is exported, it does not fix `lambda` and does not promote `L_total`, selector closure, pair12 strict-core, role transfer, bridge closure, or ToE.\n")
    append_once(AGENTS, "Current Z12 orientation-torsor global-section guardrail (P2714/S1664, 2026-06-14)", "## Current Z12 orientation-torsor global-section guardrail (P2714/S1664, 2026-06-14)\n\n- P2714 tests one new typed candidate after the P2713 intake gate: the orientation torsor of the P2708 boundary-cocycle line under `Aut(Z12)`.\n- The finite action has orientation-reversing units `7` and `11`; neither `+omega` nor `-omega` is an Aut-compatible global section, so the torsor does not export a strict mechanism fixing `lambda` or discharge `QW-2191`.\n- Do not promote the orientation-torsor candidate to pair12 strict-core, role transfer, bridge closure, `L_total`, or ToE; a next admissible move requires a new internal non-premise law breaking the torsor or a different genuinely new typed object outside closed lanes.\n")
    return payload


if __name__ == "__main__":
    main()
