#!/usr/bin/env python3
"""Scratch probe: what bridge equation is possible without the unit bit?

Question audited here: can we write equations describing a transition from the
legacy ontological/effective kernel to the strict working kernel before deriving
the one-bit unit-orientation selector?

Answer in this finite selector model:

* yes, an honest *unoriented / quotient / relation-valued* bridge ansatz can be
  written; it maps the legacy side to the Aut(Z_12) orbit {A1, A5} and keeps the
  parameter-translation slots explicit;
* no, a deterministic single-branch equation selecting the strict d5 branch is
  not available without the one-bit unit-axis datum.

No false pass: this probe does not assert K_legacy_ont == K_strict_gate, does
not transfer legacy physical roles onto K_strict_gate, and does not discharge
QW-2191.  It only classifies which equation types are well-formed before the
selector bit is supplied.
"""
from __future__ import annotations

import json
from collections import Counter
from itertools import combinations
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
OUT_JSON = HERE / "bridge_legacy_strict_unoriented_equation_without_selector_bit_audit_report.json"
OUT_MD = HERE / "bridge_legacy_strict_unoriented_equation_without_selector_bit_audit_report.md"

N = 12
ACTIVE_COUNT = 5
AUT_UNITS = [1, 5, 7, 11]
TARGET_Q_POWER = "256/243"
TARGET_ETA = "9/5"
A1 = "A1_k1_contiguous"
A5 = "A5_k5_d5"
AXES = [A1, A5]
D5_HISTOGRAM = [0, 3, 2, 1, 4, 0]
CONTIGUOUS_HISTOGRAM = [4, 3, 2, 1, 0, 0]
LEGACY_KERNEL = "K_legacy_ont(d)=alpha_geo*cos(omega_L*d+phi_L)/(1+beta_tors*d)"
STRICT_KERNEL = "K_strict_gate(d)=cos(omega_S*d+phi_S)/(1+beta*d^eta)"


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
    raise ValueError(axis)


def subset_name(subset: set[str]) -> str:
    if not subset:
        return "empty_relation"
    if subset == {A1}:
        return "singleton_A1_contiguous"
    if subset == {A5}:
        return "singleton_A5_d5"
    if subset == set(AXES):
        return "unoriented_orbit_A1_A5"
    return "+".join(sorted(subset))


def is_full_aut_invariant_subset(subset: set[str]) -> bool:
    return all(unit_action_on_axis(axis, unit) in subset for axis in subset for unit in AUT_UNITS)


def subset_rows() -> list[dict[str, Any]]:
    rows = []
    for mask in range(1 << len(AXES)):
        subset = {axis for index, axis in enumerate(AXES) if mask & (1 << index)}
        rows.append(
            {
                "name": subset_name(subset),
                "selected_axes": sorted(subset),
                "full_aut_invariant": is_full_aut_invariant_subset(subset),
                "can_be_output_of_unoriented_bridge_relation": is_full_aut_invariant_subset(subset),
                "is_singleton_d5": subset == {A5},
            }
        )
    return rows


def oriented_equation_rows() -> list[dict[str, Any]]:
    return [
        {
            "selector_bit_b": 0,
            "unit_axis_character": "{5,7}-side / mirror of d5 convention not selected",
            "oriented_branch": A1,
            "distance_histogram_d1_to_d6": CONTIGUOUS_HISTOGRAM,
            "status": "conditional branch after supplying an orientation convention",
        },
        {
            "selector_bit_b": 1,
            "unit_axis_character": "{1,11}-stabilized d5-side selected",
            "oriented_branch": A5,
            "distance_histogram_d1_to_d6": D5_HISTOGRAM,
            "status": "conditional d5 branch after supplying the missing one-bit selector",
        },
    ]


def bridge_equation_classification() -> list[dict[str, Any]]:
    return [
        {
            "equation_type": "quotient_or_relation_bridge_without_selector_bit",
            "schematic_equation": "Bridge_0(K_legacy_ont) -> Orbit_strict({A1_k1_contiguous,A5_k5_d5}) with amplitude/damping/phase maps left explicit",
            "well_formed_without_selector_bit": True,
            "selects_singleton_d5": False,
            "honest_status": "allowed as a bridge ansatz/classification target, not as strict kernel identity",
        },
        {
            "equation_type": "deterministic_strict_d5_branch_bridge_without_selector_bit",
            "schematic_equation": "Bridge_d5(K_legacy_ont) -> K_strict_gate on A5_k5_d5",
            "well_formed_without_selector_bit": False,
            "selects_singleton_d5": True,
            "honest_status": "blocked by full Aut(Z_12) unit-mirror symmetry / QW-2191 selector obstruction",
        },
        {
            "equation_type": "conditional_oriented_bridge_with_selector_bit",
            "schematic_equation": "Bridge_b(K_legacy_ont, b) -> K_strict_gate on A5 iff b=1, otherwise mirror branch A1",
            "well_formed_without_selector_bit": False,
            "selects_singleton_d5": "only conditionally when b=1 is supplied",
            "honest_status": "non-strict/axiom-augmented unless a strict source theorem derives b",
        },
    ]


def missing_bridge_slots() -> list[dict[str, str]]:
    return [
        {
            "slot": "amplitude_normalization",
            "legacy_side": "explicit alpha_geo in K_legacy_ont",
            "strict_side": "no explicit alpha_geo in K_strict_gate",
            "without_selector_bit_status": "can be written as an unknown normalization map, but not solved here",
        },
        {
            "slot": "damping_renormalization",
            "legacy_side": "linear beta_tors*d damping",
            "strict_side": "nonlinear beta*d^eta damping",
            "without_selector_bit_status": "can be written as an unknown beta_tors -> (beta, eta) map; target replay q^5=256/243 and eta=9/5 remain conditional data in this lane",
        },
        {
            "slot": "phase_frequency_translation",
            "legacy_side": "legacy omega_L, phi_L layer",
            "strict_side": "strict omega_S, phi_S layer",
            "without_selector_bit_status": "can be represented as an unknown phase/frequency map, not derived here",
        },
        {
            "slot": "unit_orientation_selector",
            "legacy_side": "no strict source currently exported",
            "strict_side": "A5/d5 branch requires one bit separating {1,11} from {5,7}",
            "without_selector_bit_status": "cannot be eliminated for singleton d5; only quotient/orbit equation is available",
        },
    ]


def exact_proof_certificate() -> dict[str, str]:
    return {
        "domain_assumption": "Before the one-bit selector is supplied, the legacy-side bridge input is treated as full-Aut invariant for the A1/A5 unit-mirror question.",
        "codomain_action": "Aut units 5 and 7 swap A1_k1_contiguous with A5_k5_d5; units 1 and 11 preserve each branch.",
        "fixed_point_fact": "The two-point branch set {A1,A5} has no singleton fixed point under the full unit action.",
        "deterministic_map_obstruction": "An equivariant deterministic map from a fixed input point to the two-point branch set would need a fixed output point; none exists.",
        "relation_map_permission": "An equivariant relation-valued map may output the whole orbit {A1,A5}; this is the maximal honest bridge equation without the bit.",
        "selector_bit_role": "Supplying b chooses one member of the orbit; b=1 is the conditional d5 convention, not a strict derivation by itself.",
    }


def main() -> None:
    supports = all_supports()
    histogram_counter = Counter(distance_histogram(support) for support in supports)
    subsets = subset_rows()
    full_aut_names = [row["name"] for row in subsets if row["full_aut_invariant"]]

    report: dict[str, Any] = {
        "result_kind": "SCRATCH_LEGACY_STRICT_UNORIENTED_EQUATION_WITHOUT_SELECTOR_BIT_AUDIT_PROBE__NOT_A_THEOREM",
        "status": "quotient-bridge-equation-possible-singleton-d5-bridge-blocked-without-unit-bit",
        "finite_model": {
            "ring": "Z_12",
            "active_count": ACTIVE_COUNT,
            "support_count": len(supports),
            "histogram_class_count": len(histogram_counter),
            "automorphism_units": AUT_UNITS,
            "survivor_axes": [
                {"name": A1, "mode": 1, "distance_histogram_d1_to_d6": CONTIGUOUS_HISTOGRAM},
                {"name": A5, "mode": 5, "distance_histogram_d1_to_d6": D5_HISTOGRAM},
            ],
            "target_q_power": TARGET_Q_POWER,
            "target_eta": TARGET_ETA,
        },
        "kernel_forms_kept_separate": {
            "legacy_kernel_form": LEGACY_KERNEL,
            "strict_kernel_form": STRICT_KERNEL,
            "identity_asserted": False,
            "reason": "Current repo guardrails require the legacy and strict kernel objects to remain split unless an explicit bridge theorem is supplied.",
        },
        "aut_invariant_output_subset_audit": {
            "selector_subset_count": len(subsets),
            "rows": subsets,
            "full_aut_invariant_output_names": full_aut_names,
            "singleton_d5_full_aut_invariant": any(row["is_singleton_d5"] and row["full_aut_invariant"] for row in subsets),
            "unoriented_orbit_output_available": "unoriented_orbit_A1_A5" in full_aut_names,
        },
        "bridge_equation_classification": bridge_equation_classification(),
        "conditional_oriented_equation_after_bit": {
            "requires_external_or_derived_bit": True,
            "bit_meaning": "one binary unit-axis record separating {1,11} from {5,7}",
            "branches": oriented_equation_rows(),
        },
        "missing_bridge_slots": missing_bridge_slots(),
        "exact_proof_certificate": exact_proof_certificate(),
        "interpretation": {
            "direct_answer": "Without the one-bit selector we can write only an unoriented/quotient bridge equation, not a deterministic equation selecting the strict d5 branch.",
            "what_can_be_done_now": "Formulate Bridge_0 as a relation from K_legacy_ont to the strict A1/A5 orbit while keeping amplitude, damping, phase/frequency, and selector slots explicit.",
            "what_cannot_be_done_now": "Do not write K_legacy_ont -> K_strict_gate(d5) as a strict theorem; that silently imports the missing bit and would falsely discharge QW-2191.",
            "matter_shannon_relation": "Matter/Shannon self-replication may later carry an already supplied bit, but it does not remove the need for a strict source of the bit in the bridge equation.",
        },
        "ontology_guardrail": {
            "allowed_reading": "The bridge ansatz may describe internal projections of the informational nadsoliton route.",
            "forbidden_reading": "No separate informational layer underneath the nadsoliton is introduced to provide the selector bit.",
            "preferred_order_preserved": "nadsoliton -> light -> matter -> emergent observer",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is asserted or used.",
            "No legacy physical-role or matter-generation claim is transferred onto K_strict_gate.",
            "No theorem derives amplitude normalization, damping renormalization, or phase/frequency translation maps.",
            "No theorem derives the required one-bit unit-axis record from strict core.",
            "Without the bit, singleton d5 branch selection is blocked by full Aut(Z_12) symmetry.",
            "The quotient/orbit bridge equation is a NOT_A_THEOREM ansatz/classification, not ToE closure.",
            "No endpoint, arrow orientation, ledger selector, positive lambda action, cycle metric source, anti-Nyquist source, fifth-mode source, future-probability source, future-value source, future-path source, matter-bit source, or legacy-strict bridge theorem is claimed.",
            "No QW-2191 discharge and no strict-core selector closure are claimed.",
            "No ToE closure is claimed.",
        ],
        "next_honest_step": "Develop either a strict source theorem for the unit-orientation bit or a formal non-bridge theorem proving that only quotient/orbit legacy->strict equations are available without extra selector data.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch legacy->strict unoriented bridge equation without selector bit audit\n\n"
        "Status: quotient bridge equation possible; singleton d5 bridge blocked without unit bit.\n\n"
        f"- Supports scanned: `{len(supports)}`; histogram classes: `{len(histogram_counter)}`.\n"
        f"- Full-Aut invariant bridge outputs: `{full_aut_names}`.\n"
        f"- Singleton d5 full-Aut invariant: `{report['aut_invariant_output_subset_audit']['singleton_d5_full_aut_invariant']}`.\n"
        f"- Unoriented orbit output available: `{report['aut_invariant_output_subset_audit']['unoriented_orbit_output_available']}`.\n"
        "- Direct answer: without the one-bit selector, only `Bridge_0(K_legacy_ont) -> {A1,A5}` is honest.\n"
        "- Forbidden shortcut: `K_legacy_ont -> K_strict_gate(d5)` would silently import the selector bit.\n"
        f"- Target replay kept as conditional data: `q^5={TARGET_Q_POWER}`, eta `{TARGET_ETA}`.\n"
        "- No false pass: no K_legacy_ont == K_strict_gate theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
