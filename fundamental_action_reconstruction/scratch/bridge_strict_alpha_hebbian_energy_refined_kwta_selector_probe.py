#!/usr/bin/env python3
"""Scratch probe: energy-refined Hebbian kWTA selector audit.

Previous Hebbian kWTA work found a real tie-boundary obstruction: a raw
centered zero-self Hebbian kWTA update reaches d5 from most states, but exact
boundary ties leave non-d5 cycles or tie-sensitive branches.  This probe tests
the next natural finite refinement:

    first apply kWTA; if the kWTA boundary is tied, retain only those candidate
    next states that maximize the same learned Hebbian energy E(T)=T^T W T.

The refinement is more internal than a lexicographic tie break because it reuses
W learned from the d5 self-record/teacher trace.  It is still a conditional
premise: this file does not derive the d5 teacher trace, centered Hebbian rule,
kWTA rule, or energy-refined selector from strict nadsoliton geometry.

Honest finite result on all 792 five-active-node supports of Z_12:
- raw kWTA branch counts are {1:624, 2:120, 3:24, 10:12, 21:12};
- energy refinement reduces them to {1:720, 2:60, 4:12};
- every retained energy-maximizing branch reaches the 12 d5 fixed points within
  three closure layers;
- every non-d5 state has strictly positive learned-energy ascent under every
  retained energy-maximizing branch, with minimum exact increase 4.

Thus the previous tie obstruction can be closed by an explicit energy-refined
Hebbian selector premise, but that premise is not exported here as strict-core
source data and does not discharge QW-2191 by itself.
"""
from __future__ import annotations

import json
from collections import Counter
from fractions import Fraction
from itertools import combinations
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
OUT_JSON = HERE / "bridge_strict_alpha_hebbian_energy_refined_kwta_selector_report.json"
OUT_MD = HERE / "bridge_strict_alpha_hebbian_energy_refined_kwta_selector_report.md"

N = 12
ACTIVE_COUNT = 5
TARGET_MODE = 5
RHO = Fraction(ACTIVE_COUNT, N)
TARGET_Q_POWER = "256/243"
TARGET_ETA = "9/5"

Support = tuple[int, ...]
Matrix = list[list[Fraction]]


def all_supports() -> list[Support]:
    return list(combinations(range(N), ACTIVE_COUNT))


def canonical_support(nodes: list[int] | tuple[int, ...]) -> Support:
    return tuple(sorted(int(node) % N for node in nodes))


def d5_teacher_orbit() -> list[Support]:
    return sorted(
        {
            canonical_support([(start + TARGET_MODE * index) % N for index in range(ACTIVE_COUNT)])
            for start in range(N)
        }
    )


def centered_activity(support: Support) -> list[Fraction]:
    active = set(support)
    return [Fraction(1) - RHO if node in active else -RHO for node in range(N)]


def build_centered_zero_self_hebbian(teacher: list[Support]) -> Matrix:
    weights = [[Fraction(0) for _col in range(N)] for _row in range(N)]
    for support in teacher:
        activity = centered_activity(support)
        for row in range(N):
            for col in range(N):
                weights[row][col] += activity[row] * activity[col]
    for index in range(N):
        weights[index][index] = Fraction(0)
    return weights


def fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def support_text(support: Support) -> str:
    return " ".join(map(str, support))


def score_nodes(support: Support, weights: Matrix) -> list[tuple[Fraction, int]]:
    return [(sum(weights[node][active] for active in support), node) for node in range(N)]


def raw_kwta_branches(support: Support, weights: Matrix) -> set[Support]:
    scores = score_nodes(support, weights)
    values = sorted({score for score, _node in scores}, reverse=True)
    selected: list[int] = []
    remaining = ACTIVE_COUNT
    for value in values:
        group = [node for score, node in scores if score == value]
        if len(group) < remaining:
            selected.extend(group)
            remaining -= len(group)
        elif len(group) == remaining:
            return {canonical_support(selected + group)}
        else:
            return {canonical_support(selected + list(choice)) for choice in combinations(group, remaining)}
    return {canonical_support(selected)}


def learned_energy(support: Support, weights: Matrix) -> Fraction:
    return sum(weights[row][col] for row in support for col in support)


def energy_refined_branches(support: Support, weights: Matrix) -> set[Support]:
    candidates = raw_kwta_branches(support, weights)
    best_energy = max(learned_energy(candidate, weights) for candidate in candidates)
    return {candidate for candidate in candidates if learned_energy(candidate, weights) == best_energy}


def closure_layers(
    supports: list[Support],
    target_set: set[Support],
    next_map: dict[Support, set[Support]],
) -> tuple[set[Support], dict[Support, int]]:
    guaranteed = set(target_set)
    layers = {support: 0 for support in target_set}
    changed = True
    layer = 0
    while changed:
        changed = False
        layer += 1
        for support in set(supports) - guaranteed:
            if next_map[support].issubset(guaranteed):
                guaranteed.add(support)
                layers[support] = layer
                changed = True
    return guaranteed, layers


def deterministic_steps(
    supports: list[Support],
    target_set: set[Support],
    next_map: dict[Support, set[Support]],
) -> dict[str, Any]:
    step_counter: Counter[int] = Counter()
    basin_counter: Counter[Support] = Counter()
    failures: list[dict[str, Any]] = []
    for initial in supports:
        current = initial
        seen = {current: 0}
        for step in range(20):
            if current in target_set:
                step_counter[step] += 1
                basin_counter[current] += 1
                break
            current = sorted(next_map[current])[0]
            if current in seen:
                failures.append(
                    {
                        "initial": list(initial),
                        "repeat_state": list(current),
                        "cycle_start_step": seen[current],
                        "repeat_step": step + 1,
                    }
                )
                break
            seen[current] = step + 1
    return {
        "tie_within_energy_refined_rule": "lexicographic only among equal-energy retained branches for this deterministic replay",
        "reached_d5_count": sum(step_counter.values()),
        "all_initial_states_reach_d5": sum(step_counter.values()) == len(supports),
        "failure_count": len(failures),
        "failure_examples": failures[:10],
        "steps_to_d5_histogram": {str(step): count for step, count in sorted(step_counter.items())},
        "d5_basin_counts": {support_text(support): count for support, count in sorted(basin_counter.items())},
    }


def energy_ascent_audit(
    supports: list[Support],
    target_set: set[Support],
    weights: Matrix,
    next_map: dict[Support, set[Support]],
) -> dict[str, Any]:
    delta_counter: Counter[Fraction] = Counter()
    non_d5_deltas: list[Fraction] = []
    violations: list[dict[str, Any]] = []
    for support in supports:
        source_energy = learned_energy(support, weights)
        for nxt in next_map[support]:
            delta = learned_energy(nxt, weights) - source_energy
            delta_counter[delta] += 1
            if support not in target_set:
                non_d5_deltas.append(delta)
                if delta <= 0:
                    violations.append(
                        {
                            "support": list(support),
                            "next": list(nxt),
                            "delta": fraction_text(delta),
                            "source_energy": fraction_text(source_energy),
                            "next_energy": fraction_text(learned_energy(nxt, weights)),
                        }
                    )
    minimum_non_d5_delta = min(non_d5_deltas)
    return {
        "energy_function": "E(S)=sum_{i,j in S} W_ij = x_S^T W x_S with diag(W)=0",
        "strict_positive_ascent_for_every_retained_non_d5_branch": not violations and minimum_non_d5_delta > 0,
        "minimum_non_d5_delta": fraction_text(minimum_non_d5_delta),
        "delta_histogram_over_retained_branches": {
            fraction_text(delta): count for delta, count in sorted(delta_counter.items())
        },
        "violation_count": len(violations),
        "violation_examples": violations[:10],
    }


def branch_refinement_audit(supports: list[Support], weights: Matrix) -> dict[str, Any]:
    raw_map = {support: raw_kwta_branches(support, weights) for support in supports}
    refined_map = {support: energy_refined_branches(support, weights) for support in supports}
    changed = [support for support in supports if raw_map[support] != refined_map[support]]
    return {
        "raw_kwta_branch_count_histogram": {
            str(size): count for size, count in sorted(Counter(len(branches) for branches in raw_map.values()).items())
        },
        "energy_refined_branch_count_histogram": {
            str(size): count for size, count in sorted(Counter(len(branches) for branches in refined_map.values()).items())
        },
        "supports_changed_by_energy_refinement_count": len(changed),
        "supports_still_multibranch_after_energy_refinement_count": sum(1 for branches in refined_map.values() if len(branches) > 1),
        "changed_examples": [
            {
                "support": list(support),
                "raw_branch_count": len(raw_map[support]),
                "energy_refined_branch_count": len(refined_map[support]),
                "raw_branches_sample": [list(row) for row in sorted(raw_map[support])[:8]],
                "energy_refined_branches": [list(row) for row in sorted(refined_map[support])],
            }
            for support in changed[:8]
        ],
        "raw_next_map": raw_map,
        "energy_refined_next_map": refined_map,
    }


def main() -> None:
    supports = all_supports()
    teacher = d5_teacher_orbit()
    target_set = set(teacher)
    weights = build_centered_zero_self_hebbian(teacher)
    branch_audit = branch_refinement_audit(supports, weights)
    refined_map = branch_audit["energy_refined_next_map"]
    guaranteed, layers = closure_layers(supports, target_set, refined_map)
    deterministic = deterministic_steps(supports, target_set, refined_map)
    energy_audit = energy_ascent_audit(supports, target_set, weights, refined_map)

    serializable_branch_audit = {key: value for key, value in branch_audit.items() if not key.endswith("_next_map")}
    layer_histogram = {str(layer): count for layer, count in sorted(Counter(layers.values()).items())}

    report: dict[str, Any] = {
        "result_kind": "SCRATCH_STRICT_ALPHA_HEBBIAN_ENERGY_REFINED_KWTA_SELECTOR_PROBE__NOT_A_THEOREM",
        "status": "conditional-energy-refined-hebbian-selector-closes-finite-d5-basin-not-strict-source",
        "finite_model": {
            "ring": "Z_12",
            "active_count": ACTIVE_COUNT,
            "support_count": len(supports),
            "teacher_record": "12 translated fifth-step/d5 supports",
            "learning_rule": "centered Hebbian covariance W=sum (x-rho)(x-rho)^T with rho=5/12 and diag(W)=0",
            "primary_update": "set-valued 5-winner-take-all at boundary ties",
            "tie_refinement": "retain only kWTA branches maximizing learned energy E(T)=T^T W T",
            "target_q_power": TARGET_Q_POWER,
            "target_eta": TARGET_ETA,
        },
        "exact_weight_replay": {
            "first_row_by_cyclic_distance": [fraction_text(value) for value in weights[0]],
            "circulant_check": all(weights[row][col] == weights[0][(col - row) % N] for row in range(N) for col in range(N)),
            "diagonal_zero_check": all(weights[index][index] == 0 for index in range(N)),
        },
        "branch_refinement_audit": serializable_branch_audit,
        "energy_refined_all_branch_closure_certificate": {
            "guaranteed_all_retained_branches_reach_d5_count": len(guaranteed),
            "guaranteed_all_retained_branches_reach_d5_fraction": f"{len(guaranteed)}/{len(supports)}",
            "closure_layer_histogram": layer_histogram,
            "max_closure_layer": max(layers.values()),
            "unclosed_count": len(set(supports) - guaranteed),
            "unclosed_examples": [list(row) for row in sorted(set(supports) - guaranteed)[:10]],
        },
        "deterministic_replay_with_lex_only_inside_equal_energy_ties": deterministic,
        "energy_ascent_certificate": energy_audit,
        "impact_on_previous_kwta_nonclosure": {
            "positive_gain": "The exact 6/792 deterministic non-d5 cycles and 24/792 raw tie-sensitive blockers disappear once ties are refined by learned Hebbian energy; every retained branch reaches d5.",
            "remaining_limit": "The closure is conditional on a supplied d5 teacher trace plus an explicit energy-refined kWTA selector premise; neither premise is derived from strict geometry here.",
            "selector_status": "candidate internal selector mechanism, still non-strict and not a QW-2191 discharge",
        },
        "ontology_guardrail": {
            "allowed_reading": "The nadsoliton itself may contain/update finite self-recorded resonance patterns and use their learned energy internally.",
            "forbidden_reading": "The Hebbian network, kWTA rule, and energy refinement are not a separate informational layer underneath the nadsoliton.",
            "preferred_order_preserved": "nadsoliton -> light -> matter -> emergent observer",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is asserted or used.",
            "No legacy physical-role claim is transferred onto K_strict_gate.",
            "No theorem derives the d5 teacher/self-record trace from strict nadsoliton geometry.",
            "No theorem derives the centered zero-self Hebbian rule, kWTA rule, or energy-refined tie selector from strict geometry.",
            "The result closes the finite basin only after the energy-refinement premise is added; it is not strict-core selector closure.",
            "No endpoint, arrow orientation, ledger selector, positive lambda action, or fifth-mode source theorem is claimed.",
            "No QW-2191 discharge and no strict-core selector closure are claimed.",
            "No ToE closure is claimed.",
        ],
        "next_honest_step": "Try to derive the energy-refined kWTA selector and d5 teacher trace from strict nadsoliton self-record structure; otherwise record this as a conditional non-strict closure mechanism rather than a strict theorem.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha Hebbian energy-refined kWTA selector probe\n\n"
        "Status: energy-refined Hebbian kWTA closes the finite d5 basin conditionally; not a strict source theorem.\n\n"
        f"- Supports scanned: `{len(supports)}` five-active-node states on `Z_12`.\n"
        f"- Raw kWTA branch histogram: `{serializable_branch_audit['raw_kwta_branch_count_histogram']}`.\n"
        f"- Energy-refined branch histogram: `{serializable_branch_audit['energy_refined_branch_count_histogram']}`.\n"
        f"- All retained energy-refined branches reach d5: `{len(guaranteed)}/{len(supports)}`.\n"
        f"- Closure layers: `{layer_histogram}`; max layer `{max(layers.values())}`.\n"
        f"- Deterministic replay reaches d5: `{deterministic['reached_d5_count']}/{len(supports)}` with histogram `{deterministic['steps_to_d5_histogram']}`.\n"
        f"- Minimum non-d5 learned-energy ascent over retained branches: `{energy_audit['minimum_non_d5_delta']}`.\n"
        f"- Target replay: `q^5={TARGET_Q_POWER}`, eta `{TARGET_ETA}`.\n"
        "- Honest read: the previous kWTA tie obstruction is closed only after adding an explicit learned-energy refinement premise.\n"
        "- No false pass: no d5-source theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
