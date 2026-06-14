#!/usr/bin/env python3
"""P2743/S1693: affine-frame transition character no-go.

After P2742 closes affine weighting of source-triple chirality, this tests a
separate signed observable candidate: inversion-odd characters of the affine
frame transition unit in U(12).  The characters are real finite sign data, but
the audit asks whether their transition ensemble exports a non-premise,
orbit-safe signed value coupled to one P2721 polarity.
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
OUT = GEN / "p2743_s1693_affine_frame_transition_character_no_go.json"
MD = GEN / "p2743_s1693_affine_frame_transition_character_no_go.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
INPUTS = {"P2742_AFFINE_WEIGHT_NO_GO": GEN / "p2742_s1692_source_triple_affine_weighted_signed_aggregate_no_go.json"}

CONTENT_PATTERNS = {
    "post_p2742_pivot_obligation": r"P2742|genuinely different strict signed observable|orbit-safe signed value|source-triple chirality lane",
    "aut_character_boundary": r"inversion-odd character|Aut\(Z12\).*character|U\(12\)|character sign|lambda.*sign",
    "affine_frame_boundary": r"affine frame|transition unit|affine transition|frame transition|source frame",
    "p2721_closure_boundary": r"P2721 polarity|lambda/P2721|QW-2191|selector closure|L_total|ToE closure",
}
NEGATIVE_EXPORT_FLAGS = [
    "affine_frame_character_exported_as_strict_source",
    "nonpremise_transition_unit_selected",
    "orbit_safe_signed_value_nonzero",
    "p2721_polarity_selected",
    "lambda_fixed",
    "qw2191_discharged",
    "selector_closure_exported",
    "role_transfer_started",
    "ltotal_promoted",
    "toe_closure_exported",
]
Z12 = tuple(range(12))
UNITS = (1, 5, 7, 11)
INV = {1: 1, 5: 5, 7: 7, 11: 11}
INVERSION_ODD_CHARACTERS = {
    "chi_5_even_7_odd": {1: 1, 5: 1, 7: -1, 11: -1},
    "chi_5_odd_7_even": {1: 1, 5: -1, 7: 1, 11: -1},
}


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"missing": True, "path": rel(path)}


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def run_rg(pattern: str) -> list[str]:
    cmd = ["rg", "-n", "--glob", "!generated/*.json", pattern, "AGENTS.md", "fundamental_action_reconstruction"]
    proc = subprocess.run(cmd, cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    if proc.returncode not in (0, 1):
        raise RuntimeError(f"rg failed for {pattern!r}: {proc.stderr}")
    return [line for line in proc.stdout.splitlines() if line.strip()]


def evidence_scan() -> dict[str, Any]:
    rows = []
    for name, pattern in CONTENT_PATTERNS.items():
        hits = run_rg(pattern)
        rows.append({"content_lane": name, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:8]})
    return {
        "content_pattern_count": len(CONTENT_PATTERNS),
        "rows": rows,
        "hit_counts": {row["content_lane"]: row["hit_count"] for row in rows},
        "all_patterns_have_hits": all(row["hit_count"] > 0 for row in rows),
    }


def transition_unit(source_unit: int, target_unit: int) -> int:
    return (target_unit * INV[source_unit]) % 12


def transition_character_audit() -> dict[str, Any]:
    frames = [(base, unit) for base in Z12 for unit in UNITS]
    transitions = []
    for source in frames:
        for target in frames:
            unit = transition_unit(source[1], target[1])
            shift = (target[0] - unit * source[0]) % 12
            transitions.append({"source": source, "target": target, "unit": unit, "shift": shift})
    unit_counts = {unit: sum(1 for row in transitions if row["unit"] == unit) for unit in UNITS}
    character_rows = []
    for name, values in INVERSION_ODD_CHARACTERS.items():
        signs = [values[row["unit"]] for row in transitions]
        unit_sums = {unit: unit_counts[unit] * values[unit] for unit in UNITS}
        character_rows.append(
            {
                "character": name,
                "values_on_units": values,
                "positive_transition_count": signs.count(1),
                "negative_transition_count": signs.count(-1),
                "global_signed_sum": sum(signs),
                "unit_signed_contributions": unit_sums,
                "unit_11_flips_character_sign": all(values[(11 * unit) % 12] == -values[unit] for unit in UNITS),
            }
        )
    return {
        "typed_candidate": "inversion-odd Aut(Z12)=U(12) character evaluated on affine-frame transition units",
        "frame_count": len(frames),
        "transition_count": len(transitions),
        "unit_counts": unit_counts,
        "transition_unit_orbit_count_under_simultaneous_affine_action": len(UNITS),
        "inversion_odd_character_count": len(INVERSION_ODD_CHARACTERS),
        "characters_with_nonzero_global_signed_sum": sum(1 for row in character_rows if row["global_signed_sum"] != 0),
        "characters_with_balanced_positive_negative_transitions": sum(1 for row in character_rows if row["positive_transition_count"] == row["negative_transition_count"]),
        "all_unit_11_sign_flip_witnesses_pass": all(row["unit_11_flips_character_sign"] for row in character_rows),
        "character_rows": character_rows,
        "finite_theorem": "The two inversion-odd U(12) characters are valid pointwise signs on affine-frame transition units, but the 48x48 transition ensemble contains each unit exactly 576 times.  Hence every inversion-odd character has 1152 positive and 1152 negative transitions and global signed sum zero; multiplication by unit 11 pairs character signs.  Without a strict source selecting one transition unit or one character polarity, this frame-character observable does not export a non-premise orbit-safe signed value or fix P2721/lambda.",
    }


def acceptance_matrix(scan: dict[str, Any], audit: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "content_grep_detects_post_p2742_pivot": scan["all_patterns_have_hits"],
        "candidate_has_pointwise_inversion_odd_signs": audit["inversion_odd_character_count"] == 2,
        "orbit_safe_signed_value_nonzero": audit["characters_with_nonzero_global_signed_sum"] > 0,
        "strict_transition_unit_source_exported": False,
        "p2721_polarity_coupling_exported": False,
    }
    return {
        "facts": facts,
        "accepted_as_strict_signed_source": all(facts.values()),
        "missing_criteria": [key for key, value in facts.items() if not value],
        "blocker": "Affine-frame transition characters are real pointwise signs, but all inversion-odd character sums are globally balanced and no strict source selects a transition unit, character polarity, or P2721 coupling.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["transition_character_audit"]
    lines = [
        "# P2743/S1693 affine-frame transition character no-go",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite transition-character audit",
        f"- frame_count={audit['frame_count']}",
        f"- transition_count={audit['transition_count']}",
        f"- unit_counts={audit['unit_counts']}",
        f"- transition_unit_orbit_count_under_simultaneous_affine_action={audit['transition_unit_orbit_count_under_simultaneous_affine_action']}",
        f"- inversion_odd_character_count={audit['inversion_odd_character_count']}",
        f"- characters_with_nonzero_global_signed_sum={audit['characters_with_nonzero_global_signed_sum']}",
        f"- characters_with_balanced_positive_negative_transitions={audit['characters_with_balanced_positive_negative_transitions']}",
        f"- all_unit_11_sign_flip_witnesses_pass={audit['all_unit_11_sign_flip_witnesses_pass']}",
        "",
        "## Theorem statement",
        audit["finite_theorem"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    inputs = {key: read_json(path) for key, path in INPUTS.items()}
    scan = evidence_scan()
    audit = transition_character_audit()
    acceptance = acceptance_matrix(scan, audit)
    payload = {
        "status": "P2743_AFFINE_FRAME_TRANSITION_CHARACTER_NO_GO" if not acceptance["accepted_as_strict_signed_source"] else "P2743_REQUIRES_MANUAL_REVIEW",
        "input_hashes": {key: sha(path) for key, path in INPUTS.items()},
        "input_statuses": {key: value.get("status") for key, value in inputs.items()},
        "audited_candidate_class": "inversion-odd affine-frame transition unit characters as a different signed observable after P2742",
        "content_evidence_scan": scan,
        "transition_character_audit": audit,
        "acceptance_matrix": acceptance,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": acceptance["blocker"],
            "next_honest_step": "Do not promote affine-frame transition characters: P2743 shows the two inversion-odd U(12) character signs are pointwise real but globally balanced across all 2304 frame transitions, and current artifacts export no strict transition-unit source or P2721 polarity coupling.  The next proof-grade move must supply a strict source law selecting a nonzero transition-unit/character polarity with P2721 coupling, or pivot outside finite Z12 sign-character/frame observables; otherwise preserve the P2697-P2743 no-new-live-frontier certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2743/S1693 affine-frame transition character no-go", "## P2743/S1693 affine-frame transition character no-go\n\n`P2743/S1693` pivots outside the source-triple chirality lane to affine-frame transition characters.  There are `48` affine frames and `2304` transitions; each unit in `U(12)={1,5,7,11}` appears `576` times.  The two inversion-odd characters are real pointwise signs, but each has `1152` positive and `1152` negative transitions, global signed sum `0`, and unit `11` pairs opposite signs.  No strict transition-unit source, `P2721` coupling theorem, `lambda/P2721` fixing, `QW-2191` discharge, selector closure, role transfer, `L_total`, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2743/S1693 affine-frame transition character Ltotal guard", "## P2743/S1693 affine-frame transition character Ltotal guard\n\n`P2743/S1693` adds no variational source term: affine-frame transition characters are globally balanced and lack a strict transition-unit source plus `P2721` coupling theorem.  Therefore this does not promote role-bearing `L_total`, selector closure, bridge closure, role transfer, or ToE closure.\n")
    append_once(AGENTS, "Current affine-frame transition character no-go guardrail (P2743/S1693, 2026-06-14)", "## Current affine-frame transition character no-go guardrail (P2743/S1693, 2026-06-14)\n\n- P2743 pivots outside the source-triple chirality lane to a different signed observable candidate: inversion-odd characters of affine-frame transition units in `U(12)`.\n- The finite computation finds `48` affine frames and `2304` transitions; each unit appears `576` times, so both inversion-odd characters have `1152` positive and `1152` negative transitions and global signed sum `0`, with unit `11` pairing opposite signs.\n- Do not promote affine-frame transition characters to `lambda/P2721` fixing, `QW-2191` discharge, selector closure, bridge closure, role transfer, `L_total`, or ToE closure without a strict transition-unit source and `P2721` coupling theorem.  A next admissible move must supply such a source law, pivot outside finite `Z12` sign-character/frame observables, or preserve the P2697-P2743 no-new-live-frontier certificate.\n")
    return payload


if __name__ == "__main__":
    main()
