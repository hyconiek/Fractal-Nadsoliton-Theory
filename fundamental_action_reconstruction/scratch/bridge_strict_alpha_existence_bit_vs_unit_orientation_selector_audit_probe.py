#!/usr/bin/env python3
"""Scratch probe: existence/nonexistence bit is not the unit-orientation bit.

User intuition audited here: perhaps the missing first selector bit is the
information that distinguishes existence from nothingness, since the nadsoliton
is primordial information and information describes state.

Finite audit:

    E = 0/1  is an existence-vs-nothingness gate.
    U in {A1, A5} is the unit-orientation branch after existence is nonempty.

The full Aut(Z_12) unit action fixes E and swaps A1 with A5.  Therefore E can
honestly gate empty versus nonempty/orbit existence, but it cannot select the
singleton d5/A5 branch unless an additional unit-axis bit is supplied.

No false pass: this does not deny the philosophical usefulness of an existence
bit.  It only proves that, in this finite selector model, existence/nonexistence
and unit orientation are different bits with different symmetry types.
"""
from __future__ import annotations

import json
from collections import Counter
from itertools import combinations
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
OUT_JSON = HERE / "bridge_strict_alpha_existence_bit_vs_unit_orientation_selector_audit_report.json"
OUT_MD = HERE / "bridge_strict_alpha_existence_bit_vs_unit_orientation_selector_audit_report.md"

N = 12
ACTIVE_COUNT = 5
AUT_UNITS = [1, 5, 7, 11]
TARGET_Q_POWER = "256/243"
TARGET_ETA = "9/5"
A1 = "A1_k1_contiguous"
A5 = "A5_k5_d5"
AXES = [A1, A5]
EXISTENCE_STATES = [0, 1]
EMPTY = "empty_or_nothingness"
D5_HISTOGRAM = [0, 3, 2, 1, 4, 0]
CONTIGUOUS_HISTOGRAM = [4, 3, 2, 1, 0, 0]


def all_supports() -> list[tuple[int, ...]]:
    return [tuple(support) for support in combinations(range(N), ACTIVE_COUNT)]


def folded(value: int) -> int:
    residue = value % N
    return min(residue, N - residue)


def distance_histogram(support: tuple[int, ...]) -> tuple[int, int, int, int, int, int]:
    counts = [0] * (N // 2)
    for left, right in combinations(support, 2):
        counts[folded(right - left) - 1] += 1
    return tuple(counts)  # type: ignore[return-value]


def unit_action_on_axis(axis: str, unit: int) -> str:
    if axis == A1:
        return A1 if folded(unit * 1) == 1 else A5
    if axis == A5:
        return A5 if folded(unit * 5) == 5 else A1
    if axis == EMPTY:
        return EMPTY
    raise ValueError(axis)


def is_invariant_output(output: set[str]) -> bool:
    return all(unit_action_on_axis(axis, unit) in output for axis in output for unit in AUT_UNITS)


def output_name(output: set[str]) -> str:
    if not output:
        return "empty_output"
    if output == {EMPTY}:
        return "nothingness_only"
    if output == {A1}:
        return "singleton_A1_contiguous"
    if output == {A5}:
        return "singleton_A5_d5"
    if output == {A1, A5}:
        return "existence_unoriented_A1_A5_orbit"
    if output == {EMPTY, A1, A5}:
        return "nothingness_plus_unoriented_orbit"
    return "+".join(sorted(output))


def allowed_outputs() -> list[set[str]]:
    base = [EMPTY, A1, A5]
    return [{axis for index, axis in enumerate(base) if mask & (1 << index)} for mask in range(1 << len(base))]


def existence_gate_functions() -> list[dict[str, Any]]:
    # Functions from E={0,1} to relation-valued outputs over {EMPTY,A1,A5}.
    outputs = allowed_outputs()
    rows = []
    for out0 in outputs:
        for out1 in outputs:
            full_aut_invariant = is_invariant_output(out0) and is_invariant_output(out1)
            existence_gate_well_formed = out0 <= {EMPTY} and A1 in out1 and A5 in out1
            singleton_d5_when_exists = out1 == {A5}
            rows.append(
                {
                    "E0_output": output_name(out0),
                    "E1_output": output_name(out1),
                    "full_aut_invariant": full_aut_invariant,
                    "honest_existence_gate": existence_gate_well_formed and full_aut_invariant,
                    "selects_singleton_d5_when_E1": singleton_d5_when_exists,
                    "classification": classify_gate(out0, out1, full_aut_invariant),
                }
            )
    return rows


def classify_gate(out0: set[str], out1: set[str], full_aut_invariant: bool) -> str:
    if out0 == {EMPTY} and out1 == {A1, A5} and full_aut_invariant:
        return "canonical_existence_bit_gates_nothingness_to_unoriented_orbit"
    if out1 == {A5}:
        return "forbidden_singleton_d5_imports_unit_orientation_bit"
    if full_aut_invariant:
        return "aut_invariant_but_not_canonical_existence_gate"
    return "not_full_aut_invariant"


def bit_character_table() -> list[dict[str, Any]]:
    return [
        {
            "datum": "existence_bit_E",
            "values": EXISTENCE_STATES,
            "aut_action": "trivial/fixed",
            "can_gate_nonempty_state": True,
            "can_select_A5_over_A1_without_extra_bit": False,
            "symmetry_type": "Aut-invariant scalar",
        },
        {
            "datum": "unit_orientation_bit_U",
            "values": ["{1,11}-d5 side", "{5,7}-mirror side"],
            "aut_action": "changes under unit mirror 5 or 7",
            "can_gate_nonempty_state": False,
            "can_select_A5_over_A1_without_extra_bit": True,
            "symmetry_type": "Aut-breaking character",
        },
    ]


def exact_proof_certificate() -> dict[str, str]:
    return {
        "existence_bit_action": "Aut(Z_12) fixes E=0 and E=1; E is a scalar yes/no datum.",
        "unit_branch_action": "Aut units 5 and 7 swap A1_k1_contiguous with A5_k5_d5.",
        "type_mismatch": "A fixed scalar bit can distinguish empty vs nonempty, but cannot choose a member of a swapped two-point orbit equivariantly.",
        "allowed_equation": "E=0 -> empty; E=1 -> {A1,A5} is Aut-invariant and well formed.",
        "blocked_equation": "E=1 -> A5/d5 is not Aut-invariant; it requires a second, Aut-breaking unit-orientation bit.",
        "information_reading": "The nadsoliton may be primordial information, but an information-bearing distinction must still have the right symmetry type to solve a selector problem.",
    }


def main() -> None:
    supports = all_supports()
    histogram_counter = Counter(distance_histogram(support) for support in supports)
    gate_rows = existence_gate_functions()
    canonical_gates = [row for row in gate_rows if row["classification"] == "canonical_existence_bit_gates_nothingness_to_unoriented_orbit"]
    invariant_singleton_d5 = [row for row in gate_rows if row["full_aut_invariant"] and row["selects_singleton_d5_when_E1"]]

    report: dict[str, Any] = {
        "result_kind": "SCRATCH_STRICT_ALPHA_EXISTENCE_BIT_VS_UNIT_ORIENTATION_SELECTOR_AUDIT_PROBE__NOT_A_THEOREM",
        "status": "existence-bit-can-gate-unoriented-orbit-not-select-d5-unit-branch",
        "finite_model": {
            "ring": "Z_12",
            "active_count": ACTIVE_COUNT,
            "support_count": len(supports),
            "histogram_class_count": len(histogram_counter),
            "automorphism_units": AUT_UNITS,
            "existence_states": EXISTENCE_STATES,
            "survivor_axes": [
                {"name": A1, "mode": 1, "distance_histogram_d1_to_d6": CONTIGUOUS_HISTOGRAM},
                {"name": A5, "mode": 5, "distance_histogram_d1_to_d6": D5_HISTOGRAM},
            ],
            "target_q_power": TARGET_Q_POWER,
            "target_eta": TARGET_ETA,
        },
        "bit_character_table": bit_character_table(),
        "existence_gate_function_audit": {
            "function_count": len(gate_rows),
            "rows": gate_rows,
            "canonical_existence_gate_count": len(canonical_gates),
            "full_aut_invariant_singleton_d5_gate_count": len(invariant_singleton_d5),
            "canonical_existence_gate_example": canonical_gates[0] if canonical_gates else None,
        },
        "exact_proof_certificate": exact_proof_certificate(),
        "interpretation": {
            "direct_correction": "The existence/nothingness bit is a plausible first scalar distinction, but it is not the same symmetry type as the missing unit-orientation bit.",
            "what_it_can_do": "It can support an equation E=0 -> empty and E=1 -> nonempty strict survivor orbit {A1,A5}.",
            "what_it_cannot_do": "It cannot by itself support E=1 -> A5/d5 without silently adding the unit-orientation bit.",
            "information_principle": "Being information is not enough; the information must transform as the selector obstruction requires.",
        },
        "ontology_guardrail": {
            "allowed_reading": "The nadsoliton itself may be read as primordial information whose internal state includes an existence/nonexistence distinction.",
            "forbidden_reading": "No separate informational layer underneath the nadsoliton is introduced to produce the existence bit or unit-orientation bit.",
            "preferred_order_preserved": "nadsoliton -> light -> matter -> emergent observer",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is asserted or used.",
            "No legacy physical-role or matter-generation claim is transferred onto K_strict_gate.",
            "No theorem derives the existence bit from strict nadsoliton geometry.",
            "No theorem derives the required one-bit unit-axis record from strict core.",
            "The existence bit does not discharge QW-2191 because it is Aut-invariant rather than Aut-breaking.",
            "Without an Aut-breaking unit-orientation bit, singleton d5 branch selection remains blocked.",
            "No endpoint, arrow orientation, ledger selector, positive lambda action, cycle metric source, anti-Nyquist source, fifth-mode source, future-probability source, future-value source, future-path source, matter-bit source, existence-bit source, or legacy-strict bridge theorem is claimed.",
            "No QW-2191 discharge and no strict-core selector closure are claimed.",
            "No ToE closure is claimed.",
        ],
        "next_honest_step": "If pursuing the origin intuition, search for a strict theorem that turns an existence/nonexistence distinction into an Aut-breaking unit character; otherwise keep them as two different bits.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha existence bit vs unit-orientation selector audit\n\n"
        "Status: existence bit gates nonempty orbit but does not select d5 branch.\n\n"
        f"- Supports scanned: `{len(supports)}`; histogram classes: `{len(histogram_counter)}`.\n"
        f"- Existence gate functions audited: `{len(gate_rows)}`.\n"
        f"- Canonical existence gates E=0->empty, E=1->{{A1,A5}}: `{len(canonical_gates)}`.\n"
        f"- Full-Aut invariant singleton-d5 gates from E=1: `{len(invariant_singleton_d5)}`.\n"
        "- Direct correction: existence/nonexistence is an Aut-invariant scalar bit; unit orientation is an Aut-breaking bit.\n"
        "- Honest equation: `E=0 -> empty`, `E=1 -> {A1,A5}`.\n"
        "- Forbidden shortcut: `E=1 -> A5/d5` silently imports the unit-orientation bit.\n"
        f"- Target replay kept conditional: `q^5={TARGET_Q_POWER}`, eta `{TARGET_ETA}`.\n"
        "- No false pass: no existence-bit source theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
