#!/usr/bin/env python3
"""Scratch probe: centered Hebbian kWTA attractor-basin certificate.

This is the next finite/proof-oriented step after the Hebbian resonance-learning
spectral audit.  The previous audit showed that Hebbian learning can amplify a
supplied fifth-mode/d5 teacher trace, but does not derive that trace from an
unbiased start.  This probe asks the stronger dynamical question:

    If the nadsoliton already contains a d5 self-record/teacher trace, does a
    concrete Hebbian update actually drive arbitrary 5-node states toward the
    d5 resonance attractor family?

Finite model:
- state: a 5-active-node binary support S on Z_12;
- teacher record: the 12 translated fifth-step/d5 supports;
- learning rule: centered Hebbian covariance, W=sum (x-rho)(x-rho)^T, rho=5/12;
- no self-couplings: diag(W)=0;
- update rule: k-winner-take-all (k=5) from scores y_i=sum_j W_ij x_j.

Honest finite result:
- Exact deterministic lexicographic tie resolution reaches the d5 attractor
  family from 786/792 states within at most four updates, and exposes 6 exact
  deterministic tie-artifact non-d5 cycles.
- The stronger tie-robust claim is also not fully true: if all boundary ties are
  treated set-valuedly, 768/792 states are guaranteed to reach d5 for every tie
  branch, while 24/792 states are tie-sensitive.  All 792 states still have at
  least one tie branch reaching d5.
- Therefore this is a conditional dynamical mechanism/nonclosure certificate, not a strict
  derivation of the d5 teacher trace, not a canonical tie/source theorem, and
  not a QW-2191 discharge.
"""
from __future__ import annotations

import json
from collections import Counter, deque
from fractions import Fraction
from itertools import combinations
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
OUT_JSON = HERE / "bridge_strict_alpha_hebbian_kwta_attractor_basin_certificate_report.json"
OUT_MD = HERE / "bridge_strict_alpha_hebbian_kwta_attractor_basin_certificate_report.md"

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
    support_set = set(support)
    return [Fraction(1) - RHO if node in support_set else -RHO for node in range(N)]


def build_centered_zero_self_hebbian(teacher: list[Support]) -> Matrix:
    weights = [[Fraction(0) for _ in range(N)] for _ in range(N)]
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


def first_row_certificate(weights: Matrix) -> list[str]:
    return [fraction_text(value) for value in weights[0]]


def score_nodes(support: Support, weights: Matrix) -> list[tuple[Fraction, int]]:
    return [(sum(weights[node][active] for active in support), node) for node in range(N)]


def deterministic_kwta_update(support: Support, weights: Matrix) -> Support:
    ranked = sorted(score_nodes(support, weights), key=lambda row: (-row[0], row[1]))
    return canonical_support([node for _score, node in ranked[:ACTIVE_COUNT]])


def boundary_tie_nodes(support: Support, weights: Matrix) -> list[int]:
    ranked = sorted(score_nodes(support, weights), key=lambda row: (-row[0], row[1]))
    boundary_score = ranked[ACTIVE_COUNT - 1][0]
    next_score = ranked[ACTIVE_COUNT][0]
    if boundary_score != next_score:
        return []
    return [node for score, node in ranked if score == boundary_score]


def set_valued_kwta_updates(support: Support, weights: Matrix) -> set[Support]:
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


def deterministic_attractor_audit(supports: list[Support], teacher_set: set[Support], weights: Matrix) -> dict[str, Any]:
    step_counter: Counter[int] = Counter()
    basin_counter: Counter[Support] = Counter()
    fixed_failures: list[dict[str, Any]] = []
    non_d5_attractors: list[dict[str, Any]] = []
    max_steps = 0

    for teacher_support in sorted(teacher_set):
        updated = deterministic_kwta_update(teacher_support, weights)
        if updated != teacher_support:
            fixed_failures.append({"support": list(teacher_support), "update": list(updated)})

    for initial in supports:
        seen: dict[Support, int] = {initial: 0}
        current = initial
        for step in range(0, 20):
            if current in teacher_set:
                step_counter[step] += 1
                basin_counter[current] += 1
                max_steps = max(max_steps, step)
                break
            nxt = deterministic_kwta_update(current, weights)
            if nxt in seen:
                non_d5_attractors.append(
                    {
                        "initial": list(initial),
                        "cycle_start_step": seen[nxt],
                        "repeat_step": step + 1,
                        "repeated_state": list(nxt),
                    }
                )
                break
            seen[nxt] = step + 1
            current = nxt

    return {
        "all_teacher_patterns_fixed": not fixed_failures,
        "teacher_fixed_failures": fixed_failures,
        "all_initial_states_reach_d5": sum(step_counter.values()) == len(supports),
        "support_count": len(supports),
        "reached_d5_count": sum(step_counter.values()),
        "non_d5_attractor_count": len(non_d5_attractors),
        "max_steps_to_d5": max_steps,
        "steps_to_d5_histogram": {str(step): count for step, count in sorted(step_counter.items())},
        "d5_basin_counts": {" ".join(map(str, support)): count for support, count in sorted(basin_counter.items())},
        "non_d5_attractor_examples": non_d5_attractors[:10],
    }


def tie_branch_closure_audit(supports: list[Support], teacher_set: set[Support], weights: Matrix) -> dict[str, Any]:
    support_set = set(supports)
    next_map = {support: set_valued_kwta_updates(support, weights) for support in supports}
    branch_counter = Counter(len(nexts) for nexts in next_map.values())

    guaranteed = set(teacher_set)
    changed = True
    while changed:
        changed = False
        for support in support_set - guaranteed:
            if next_map[support].issubset(guaranteed):
                guaranteed.add(support)
                changed = True

    existential = set(teacher_set)
    changed = True
    while changed:
        changed = False
        for support in support_set - existential:
            if next_map[support] & existential:
                existential.add(support)
                changed = True

    tie_sensitive = sorted(support_set - guaranteed)
    boundary_tie_counter: Counter[int] = Counter()
    boundary_tie_examples: list[dict[str, Any]] = []
    for support in supports:
        tie_nodes = boundary_tie_nodes(support, weights)
        if tie_nodes:
            boundary_tie_counter[len(tie_nodes)] += 1
            if len(boundary_tie_examples) < 12:
                boundary_tie_examples.append(
                    {
                        "support": list(support),
                        "boundary_tie_nodes": tie_nodes,
                        "branch_count": len(next_map[support]),
                        "example_updates": [list(row) for row in sorted(next_map[support])[:6]],
                    }
                )

    return {
        "set_valued_rule": "At a kWTA boundary tie, keep every top-k completion instead of using lexicographic tie resolution.",
        "support_count": len(supports),
        "branch_count_histogram": {str(branches): count for branches, count in sorted(branch_counter.items())},
        "boundary_tie_support_count": sum(boundary_tie_counter.values()),
        "boundary_tie_size_histogram": {str(size): count for size, count in sorted(boundary_tie_counter.items())},
        "guaranteed_all_tie_branches_reach_d5_count": len(guaranteed),
        "guaranteed_all_tie_branches_reach_d5_fraction": f"{len(guaranteed)}/{len(supports)}",
        "existential_some_tie_branch_reaches_d5_count": len(existential),
        "existential_some_tie_branch_reaches_d5_fraction": f"{len(existential)}/{len(supports)}",
        "tie_sensitive_blocker_count": len(tie_sensitive),
        "tie_sensitive_blocker_fraction": f"{len(tie_sensitive)}/{len(supports)}",
        "tie_sensitive_blocker_examples": [list(row) for row in tie_sensitive[:12]],
        "boundary_tie_examples": boundary_tie_examples,
    }


def main() -> None:
    supports = all_supports()
    teacher = d5_teacher_orbit()
    teacher_set = set(teacher)
    weights = build_centered_zero_self_hebbian(teacher)
    deterministic = deterministic_attractor_audit(supports, teacher_set, weights)
    tie_audit = tie_branch_closure_audit(supports, teacher_set, weights)

    report: dict[str, Any] = {
        "result_kind": "SCRATCH_STRICT_ALPHA_HEBBIAN_KWTA_ATTRACTOR_BASIN_CERTIFICATE_PROBE__NOT_A_THEOREM",
        "status": "conditional-d5-teacher-hebbian-dynamics-certificate-with-exact-tie-nonclosure",
        "finite_model": {
            "ring": "Z_12",
            "active_count": ACTIVE_COUNT,
            "support_count": len(supports),
            "teacher_record": "12 translated fifth-step/d5 supports",
            "learning_rule": "centered Hebbian covariance W=sum (x-rho)(x-rho)^T with rho=5/12 and diag(W)=0",
            "update_rule": "5-winner-take-all from y_i=sum_j W_ij x_j",
            "target_q_power": TARGET_Q_POWER,
            "target_eta": TARGET_ETA,
        },
        "exact_weight_certificate": {
            "first_row_by_cyclic_distance": first_row_certificate(weights),
            "circulant_check": all(weights[row][col] == weights[0][(col - row) % N] for row in range(N) for col in range(N)),
            "diagonal_zero_check": all(weights[index][index] == 0 for index in range(N)),
            "positive_d5_neighbor_weight": "23/12 at cyclic distances 5 and 7",
            "negative_opposite_parity_weight": "-25/12 at cyclic distances 1, 6, and 11",
        },
        "deterministic_lexicographic_kwta_certificate": deterministic,
        "set_valued_tie_branch_audit": tie_audit,
        "impact_on_hebbian_hypothesis": {
            "positive_gain": "Given a centered zero-self Hebbian d5 self-record/teacher trace, exact kWTA dynamics reaches the d5 resonance attractor family from 786/792 states and all 12 teacher patterns are fixed.",
            "honest_limit": "The exact all-state deterministic claim is false: 6/792 states enter non-d5 tie-artifact cycles, and if boundary ties are treated as unresolved selector data, only 768/792 states are guaranteed for every tie branch.",
            "meaning_for_highest_resonance": "Hebbian learning can be upgraded from a spectral amplification story to a concrete conditional attractor mechanism for most states and an existential branch for all states, but the d5 teacher trace and tie selector are still extra premises.",
            "selector_status": "conditional non-strict mechanism certificate; no QW-2191 discharge",
        },
        "ontology_guardrail": {
            "allowed_reading": "The nadsoliton itself may contain and update finite self-recorded formation/resonance patterns as primordial information in a solitonic state.",
            "forbidden_reading": "The Hebbian network/dynamics is not a separate informational substrate under the nadsoliton.",
            "preferred_order_preserved": "nadsoliton -> light -> matter -> emergent observer",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is asserted or used.",
            "No legacy physical-role claim is transferred onto K_strict_gate.",
            "No theorem derives the d5 teacher/self-record trace from strict nadsoliton geometry.",
            "No theorem derives the centered zero-self Hebbian rule, kWTA rule, or lexicographic tie selector from strict geometry.",
            "The exact deterministic all-state version is not fully closed: 6/792 states enter non-d5 tie-artifact cycles.",
            "The tie-robust all-branches version is not fully closed: 24/792 states are tie-sensitive blockers.",
            "No endpoint, arrow orientation, ledger selector, positive lambda action, or fifth-mode source theorem is claimed.",
            "No QW-2191 discharge and no strict-core selector closure are claimed.",
            "No ToE closure is claimed.",
        ],
        "next_honest_step": "Either seek a strict internal source for the d5 teacher trace plus a canonical tie selector, or record the non-bridge result that Hebbian kWTA dynamics is only a conditional attractor mechanism once those premises are supplied.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha Hebbian kWTA attractor-basin certificate probe\n\n"
        "Status: conditional d5-teacher Hebbian dynamics certificate with exact tie nonclosure.\n\n"
        f"- Supports scanned: `{len(supports)}` five-active-node states on `Z_12`.\n"
        f"- Exact first row of centered zero-self Hebbian d5 weights: `{first_row_certificate(weights)}`.\n"
        f"- Deterministic lexicographic kWTA reaches d5 count: `{deterministic['reached_d5_count']}/{len(supports)}`.\n"
        f"- Deterministic non-d5 attractor/cycle count: `{deterministic['non_d5_attractor_count']}`.\n"
        f"- Deterministic max steps to d5 among reached states: `{deterministic['max_steps_to_d5']}` with histogram `{deterministic['steps_to_d5_histogram']}`.\n"
        f"- Set-valued all-tie-branches guaranteed basin: `{tie_audit['guaranteed_all_tie_branches_reach_d5_fraction']}`.\n"
        f"- Set-valued some-branch/existential basin: `{tie_audit['existential_some_tie_branch_reaches_d5_fraction']}`.\n"
        f"- Tie-sensitive blockers: `{tie_audit['tie_sensitive_blocker_fraction']}`.\n"
        f"- Target replay: `q^5={TARGET_Q_POWER}`, eta `{TARGET_ETA}`.\n"
        "- Honest read: Hebbian kWTA gives a strong but nonclosed conditional attractor mechanism after a d5 teacher trace is supplied; exact ties still require extra selector data.\n"
        "- No false pass: no d5-source theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
