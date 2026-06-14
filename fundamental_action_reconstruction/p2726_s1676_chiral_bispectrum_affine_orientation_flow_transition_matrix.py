#!/usr/bin/env python3
"""P2726/S1676: affine orientation-flow transition matrix for Im(B_{1,5}).

P2725 closed the pure translation-flow velocity of the P2718 marker: all
orientation-preserving Z12 source flows have zero signed velocity.  P2726 tests
the next finite dynamic extension without pretending it is sourced: allow the
transition to preserve or flip the orientation torsor while translating the
source, then enumerate the complete 24-state affine transition matrix.

The calculation separates two facts: orientation-preserving flows remain zero,
while orientation-flipping transitions have nonzero +/-4 jumps.  The nonzero
jumps are real but they are conditional on importing an orientation-flip branch;
current artifacts do not export a strict law selecting that flip direction or a
P2721 coupling polarity.  Therefore this matrix is evidence of where a source
would have to enter, not selector closure.
"""
from __future__ import annotations

import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2726_s1676_chiral_bispectrum_affine_orientation_flow_transition_matrix.json"
MD = GEN / "p2726_s1676_chiral_bispectrum_affine_orientation_flow_transition_matrix.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

Z12 = 12
VELOCITIES = tuple(range(Z12))
ORIENTATION_MULTIPLIERS = (1, -1)
INPUTS = {
    "P2725_TRANSLATION_FLOW_NO_GO": GEN / "p2725_s1675_chiral_bispectrum_translation_flow_signed_velocity_no_go.json",
    "P2721_SIGN_TORSOR_COUPLING": GEN / "p2721_s1671_chiral_bispectrum_sign_to_orientation_torsor_coupling_audit.json",
    "P2718_CHIRAL_BISPECTRUM_MARKER": GEN / "p2718_s1668_chiral_bispectrum_explicit_signed_formula_torsor_coupling_audit.json",
}

NEGATIVE_EXPORT_FLAGS = [
    "strict_orientation_flip_source_law_exported",
    "nonpremise_orientation_flip_direction_selected",
    "p2721_coupling_polarity_selected",
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


def marker_table(p2718: dict[str, Any]) -> dict[tuple[int, int], float]:
    return {
        (int(row["orientation"]), int(row["source"])): float(row["marker_imag"])
        for row in p2718.get("marker_rows", [])
    }


def transition_matrix(table: dict[tuple[int, int], float]) -> dict[str, Any]:
    rows = []
    by_multiplier: dict[int, list[float]] = defaultdict(list)
    by_rule: Counter[str] = Counter()
    for orientation in (-1, 1):
        for source in range(Z12):
            current_value = table[(orientation, source)]
            for multiplier in ORIENTATION_MULTIPLIERS:
                for velocity in VELOCITIES:
                    next_orientation = multiplier * orientation
                    next_source = (source + velocity) % Z12
                    next_value = table[(next_orientation, next_source)]
                    delta = round(next_value - current_value, 9)
                    if delta == -0.0:
                        delta = 0.0
                    rule = "orientation_preserving_translation" if multiplier == 1 else "orientation_flipping_affine_transition"
                    rows.append({
                        "rule": rule,
                        "orientation": orientation,
                        "source": source,
                        "orientation_multiplier": multiplier,
                        "velocity": velocity,
                        "next_orientation": next_orientation,
                        "next_source": next_source,
                        "marker_imag": current_value,
                        "next_marker_imag": next_value,
                        "delta_marker_imag": delta,
                        "nonzero_signed_delta": abs(delta) > 1e-9,
                    })
                    by_multiplier[multiplier].append(delta)
                    by_rule[rule] += 1

    preserving = [row for row in rows if row["orientation_multiplier"] == 1]
    flipping = [row for row in rows if row["orientation_multiplier"] == -1]
    preserving_nonzero = [row for row in preserving if row["nonzero_signed_delta"]]
    flipping_nonzero = [row for row in flipping if row["nonzero_signed_delta"]]
    return {
        "checked_transition_rows": len(rows),
        "rule_counts": dict(by_rule),
        "velocity_count_including_zero": len(VELOCITIES),
        "orientation_multiplier_values": list(ORIENTATION_MULTIPLIERS),
        "orientation_preserving_delta_values": sorted(set(by_multiplier[1])),
        "orientation_flipping_delta_values": sorted(set(by_multiplier[-1])),
        "orientation_preserving_nonzero_count": len(preserving_nonzero),
        "orientation_flipping_nonzero_count": len(flipping_nonzero),
        "orientation_flipping_delta_histogram": dict(sorted(Counter(row["delta_marker_imag"] for row in flipping).items())),
        "all_preserving_deltas_zero": len(preserving_nonzero) == 0,
        "all_flipping_deltas_nonzero": len(flipping_nonzero) == len(flipping),
        "sample_flipping_rows": flipping[:4],
        "obstruction": "The only nonzero signed jumps in the complete affine transition matrix occur when the transition flips the orientation torsor.  Since current artifacts do not export a non-premise law selecting such a flip direction or its P2721 coupling polarity, the nonzero +/-4 jumps are conditional evidence, not a strict source.",
    }


def acceptance_matrix(matrix: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "complete_affine_transition_matrix_enumerated": matrix["checked_transition_rows"] == 576,
        "pure_translation_subcase_matches_p2725_zero": matrix["all_preserving_deltas_zero"],
        "orientation_flip_has_nonzero_signed_jump": matrix["all_flipping_deltas_nonzero"],
        "strict_orientation_flip_source_law_exported": False,
        "nonpremise_orientation_flip_direction_selected": False,
        "p2721_coupling_polarity_selected": False,
    }
    criteria = list(facts)
    missing = [criterion for criterion in criteria if not facts[criterion]]
    return {
        "criteria": criteria,
        "facts": facts,
        "missing_criteria": missing,
        "accepted_as_strict_dynamic_chiral_source": not missing,
        "blocker": "The matrix exposes a nonzero orientation-flip jump, but a flip is exactly the unsourced torsor branch choice.  Without a strict law selecting that branch and the P2721 polarity, it cannot fix lambda or discharge QW-2191.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    matrix = payload["transition_matrix"]
    lines = [
        "# P2726/S1676 chiral-bispectrum affine orientation-flow transition matrix",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite transition result",
        f"- checked_transition_rows={matrix['checked_transition_rows']}",
        f"- orientation_preserving_delta_values={matrix['orientation_preserving_delta_values']}",
        f"- orientation_flipping_delta_values={matrix['orientation_flipping_delta_values']}",
        f"- orientation_flipping_delta_histogram={matrix['orientation_flipping_delta_histogram']}",
        matrix["obstruction"],
        "",
        "## Acceptance",
        f"- accepted_as_strict_dynamic_chiral_source={payload['acceptance_matrix']['accepted_as_strict_dynamic_chiral_source']}",
        f"- missing={payload['acceptance_matrix']['missing_criteria']}",
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    inputs = {key: read_json(path) for key, path in INPUTS.items()}
    matrix = transition_matrix(marker_table(inputs["P2718_CHIRAL_BISPECTRUM_MARKER"]))
    acceptance = acceptance_matrix(matrix)
    bounded_no_go = matrix["all_preserving_deltas_zero"] and matrix["all_flipping_deltas_nonzero"] and not acceptance["accepted_as_strict_dynamic_chiral_source"]
    payload = {
        "status": "P2726_AFFINE_ORIENTATION_FLOW_NONZERO_BUT_UNSOURCED_NO_STRICT_DYNAMIC_CHIRAL_SOURCE" if bounded_no_go else "P2726_REQUIRES_MANUAL_REVIEW",
        "input_hashes": {key: sha(path) for key, path in INPUTS.items()},
        "input_statuses": {key: value.get("status") for key, value in inputs.items()},
        "audited_dynamic_candidate": "complete affine source-translation plus optional orientation-flip transition matrix for P2718 Im(B_{1,5})",
        "transition_matrix": matrix,
        "acceptance_matrix": acceptance,
        "decision": {
            "strict_orientation_flip_source_law_exported": False,
            "nonpremise_orientation_flip_direction_selected": False,
            "p2721_coupling_polarity_selected": False,
            "strict_mechanism_fixing_lambda_exported": False,
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": "P2726 enumerates all 576 affine transitions on the 24 P2718 marker states.  Orientation-preserving translation reproduces the P2725 zero-velocity no-go.  Orientation-flipping transitions give nonzero signed jumps (+/-4), but only by importing the very orientation torsor flip that lacks a non-premise source law and P2721 polarity selection.",
            "next_honest_step": "The next admissible step is no longer to search translation-flow variants.  It must either export a strict orientation-flip/chiral-flow source law that selects one branch and one P2721 polarity, then rerun this matrix as an acceptance test, or preserve the P2697-P2726 no-new-live-frontier certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2726/S1676 chiral-bispectrum affine orientation-flow transition matrix", "## P2726/S1676 chiral-bispectrum affine orientation-flow transition matrix\n\n`P2726/S1676` enumerates the complete affine transition matrix for the P2718 marker under source translation plus optional orientation-torsor flip.  The 288 orientation-preserving rows have zero signed finite difference, matching P2725; the 288 orientation-flipping rows have nonzero `+/-4` jumps, but only by importing an unsourced flip branch.  No strict orientation-flip source law, P2721 polarity selection, `lambda` fixing, `QW-2191` discharge, role transfer, `L_total`, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2726/S1676 unsourced orientation-flip Ltotal guard", "## P2726/S1676 unsourced orientation-flip Ltotal guard\n\n`P2726/S1676` finds nonzero jumps only in orientation-flipping affine transitions of `Im(B_{1,5})`.  Because the flip branch is not sourced by a strict law and no P2721 polarity is selected, those jumps cannot be promoted to a variational source term or `L_total`.\n")
    append_once(AGENTS, "Current affine orientation-flow transition matrix guardrail (P2726/S1676, 2026-06-14)", "## Current affine orientation-flow transition matrix guardrail (P2726/S1676, 2026-06-14)\n\n- P2726 exhausts the finite affine transition matrix of the P2718 marker under Z12 source translation plus optional orientation-torsor flip: 576 rows total.\n- The orientation-preserving half reproduces the P2725 zero-flow no-go; the orientation-flipping half has nonzero `+/-4` marker jumps, but only by importing the unsourced torsor branch choice.\n- Do not promote these conditional flip jumps to `QW-2191` discharge, selector closure, pair12 strict-core upgrade, role transfer, bridge closure, `L_total`, or ToE.  A next admissible move must export a strict orientation-flip/chiral-flow source law selecting one branch and one P2721 polarity, or preserve the P2697-P2726 no-new-live-frontier certificate.\n")
    return payload


if __name__ == "__main__":
    main()
