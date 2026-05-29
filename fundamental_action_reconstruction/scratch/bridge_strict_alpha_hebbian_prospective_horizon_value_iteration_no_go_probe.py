#!/usr/bin/env python3
"""Scratch probe: finite-horizon prospective value iteration still needs an orientation bit.

The previous prospective-state probability probe showed that a single full-Aut
invariant future probability cannot break the A1/A5 unit-mirror tie.  This probe
checks the stronger temporal/lookahead version: can finite-horizon "highest
possible future resonance" value iteration create a d5 selector from symmetric
prospective dynamics?

Finite model audited here:

    A1 = (k=1, contiguous histogram)
    A5 = (k=5, d5 histogram)

with the nontrivial Aut(Z_12) action swapping A1 and A5.  Every exact rational
full-Aut-equivariant two-state Markov kernel has the form

    T_a = [[a, 1-a], [1-a, a]].

No false pass: invariant terminal/future value remains tied for every enumerated
kernel and every finite horizon.  A d5-favoring terminal value is allowed as a
conditional future-resonance premise, but it is exactly the missing unit-axis bit
in value-language, not a strict-core derivation.
"""
from __future__ import annotations

import json
from collections import Counter
from fractions import Fraction
from itertools import combinations, product
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
OUT_JSON = HERE / "bridge_strict_alpha_hebbian_prospective_horizon_value_iteration_no_go_report.json"
OUT_MD = HERE / "bridge_strict_alpha_hebbian_prospective_horizon_value_iteration_no_go_report.md"

N = 12
ACTIVE_COUNT = 5
AUT_UNITS = [1, 5, 7, 11]
TARGET_Q_POWER = "256/243"
TARGET_ETA = "9/5"
GRID_DENOMINATOR = 12
MAX_HORIZON = 8
A1 = "A1_k1_contiguous"
A5 = "A5_k5_d5"
AXES = [A1, A5]
D5_HISTOGRAM = [0, 3, 2, 1, 4, 0]
CONTIGUOUS_HISTOGRAM = [4, 3, 2, 1, 0, 0]
Vector = tuple[Fraction, Fraction]
Kernel = tuple[tuple[Fraction, Fraction], tuple[Fraction, Fraction]]


def fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def vector_json(vector: Vector) -> dict[str, Any]:
    return {
        "A1": fraction_text(vector[0]),
        "A5": fraction_text(vector[1]),
        "tied": vector[0] == vector[1],
        "d5_dominant": vector[1] > vector[0],
    }


def kernel_json(kernel: Kernel) -> dict[str, Any]:
    return {
        "from_A1": {"to_A1": fraction_text(kernel[0][0]), "to_A5": fraction_text(kernel[0][1])},
        "from_A5": {"to_A1": fraction_text(kernel[1][0]), "to_A5": fraction_text(kernel[1][1])},
    }


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


def swap(vector: Vector) -> Vector:
    return vector[1], vector[0]


def is_invariant_vector(vector: Vector) -> bool:
    return vector == swap(vector)


def equivariant_kernel(a: Fraction) -> Kernel:
    return ((a, 1 - a), (1 - a, a))


def all_equivariant_kernels(denominator: int) -> list[Kernel]:
    return [equivariant_kernel(Fraction(numerator, denominator)) for numerator in range(denominator + 1)]


def apply_kernel(kernel: Kernel, vector: Vector) -> Vector:
    return (
        kernel[0][0] * vector[0] + kernel[0][1] * vector[1],
        kernel[1][0] * vector[0] + kernel[1][1] * vector[1],
    )


def iterate_kernel(kernel: Kernel, vector: Vector, horizon: int) -> Vector:
    result = vector
    for _ in range(horizon):
        result = apply_kernel(kernel, result)
    return result


def bellman_max_step(kernels: list[Kernel], vector: Vector) -> Vector:
    next_values = [apply_kernel(kernel, vector) for kernel in kernels]
    return max(candidate[0] for candidate in next_values), max(candidate[1] for candidate in next_values)


def bellman_max_trace(kernels: list[Kernel], terminal: Vector, horizon: int) -> list[Vector]:
    trace = [terminal]
    current = terminal
    for _ in range(horizon):
        current = bellman_max_step(kernels, current)
        trace.append(current)
    return trace


def reward_grid(denominator: int) -> list[Vector]:
    return [
        (Fraction(a1_reward, denominator), Fraction(a5_reward, denominator))
        for a1_reward, a5_reward in product(range(denominator + 1), repeat=2)
    ]


def reward_grid_audit(rewards: list[Vector]) -> dict[str, Any]:
    invariant_rewards = [reward for reward in rewards if is_invariant_vector(reward)]
    d5_dominant_rewards = [reward for reward in rewards if reward[1] > reward[0]]
    return {
        "denominator": GRID_DENOMINATOR,
        "reward_count": len(rewards),
        "full_aut_invariant_reward_count": len(invariant_rewards),
        "d5_dominant_reward_count": len(d5_dominant_rewards),
        "full_aut_invariant_d5_dominant_reward_count": sum(1 for reward in invariant_rewards if reward[1] > reward[0]),
        "sample_invariant_rewards": [vector_json(reward) for reward in [invariant_rewards[0], invariant_rewards[len(invariant_rewards) // 2], invariant_rewards[-1]]],
        "sample_d5_dominant_rewards": [vector_json(reward) for reward in d5_dominant_rewards[:3]],
    }


def fixed_kernel_horizon_audit(kernels: list[Kernel], invariant_rewards: list[Vector]) -> dict[str, Any]:
    rows = []
    violation_count = 0
    for horizon in range(MAX_HORIZON + 1):
        tied_count = 0
        d5_dominant_count = 0
        for kernel in kernels:
            for reward in invariant_rewards:
                value = iterate_kernel(kernel, reward, horizon)
                tied_count += int(value[0] == value[1])
                d5_dominant_count += int(value[1] > value[0])
        expected = len(kernels) * len(invariant_rewards)
        if tied_count != expected or d5_dominant_count != 0:
            violation_count += 1
        rows.append(
            {
                "horizon": horizon,
                "cases_checked": expected,
                "tied_value_count": tied_count,
                "d5_dominant_value_count": d5_dominant_count,
            }
        )
    return {
        "max_horizon": MAX_HORIZON,
        "rows": rows,
        "invariant_reward_value_tie_violations": violation_count,
        "all_invariant_rewards_remain_tied": violation_count == 0,
    }


def biased_reward_horizon_samples(kernels: list[Kernel]) -> list[dict[str, Any]]:
    d5_terminal = (Fraction(0), Fraction(1))
    samples = []
    for a_numerator in [0, GRID_DENOMINATOR // 2, GRID_DENOMINATOR]:
        kernel = equivariant_kernel(Fraction(a_numerator, GRID_DENOMINATOR))
        samples.append(
            {
                "kernel_parameter_a": fraction_text(Fraction(a_numerator, GRID_DENOMINATOR)),
                "kernel": kernel_json(kernel),
                "values_by_horizon": [vector_json(iterate_kernel(kernel, d5_terminal, horizon)) for horizon in range(4)],
                "interpretation": "d5 terminal value is non-invariant input; dynamics transports or averages that imported bit.",
            }
        )
    controlled_trace = bellman_max_trace(kernels, d5_terminal, 3)
    samples.append(
        {
            "controlled_action_set": "all denominator-12 full-Aut-equivariant kernels T_a",
            "terminal_value": vector_json(d5_terminal),
            "bellman_max_values_by_horizon": [vector_json(value) for value in controlled_trace],
            "interpretation": "max over possible equivariant futures lets both axes reach the d5-labelled future, so it does not select a unique current axis.",
        }
    )
    return samples


def exact_proof_certificate() -> dict[str, Any]:
    return {
        "swap_operator": "J(A1,A5)=(A5,A1)",
        "equivariance_condition": "T J = J T",
        "two_state_solution": "T_a = [[a, 1-a], [1-a, a]]",
        "invariant_value_condition": "J v = v iff v_A1 = v_A5",
        "closure_step": "If T J = J T and J v = v, then J(T v)=T(J v)=T v; hence T v is tied.",
        "finite_horizon_induction": "Repeat the closure step for any finite horizon h; no h creates singleton d5 from invariant terminal data.",
        "controlled_bellman_note": "The max over an Aut-closed action family also maps invariant value vectors to invariant value vectors.",
    }


def main() -> None:
    supports = all_supports()
    histogram_counter = Counter(distance_histogram(support) for support in supports)
    kernels = all_equivariant_kernels(GRID_DENOMINATOR)
    rewards = reward_grid(GRID_DENOMINATOR)
    invariant_rewards = [reward for reward in rewards if is_invariant_vector(reward)]
    reward_audit = reward_grid_audit(rewards)
    fixed_audit = fixed_kernel_horizon_audit(kernels, invariant_rewards)
    invariant_terminal = (Fraction(1, 2), Fraction(1, 2))
    controlled_invariant_trace = bellman_max_trace(kernels, invariant_terminal, MAX_HORIZON)

    report: dict[str, Any] = {
        "result_kind": "SCRATCH_STRICT_ALPHA_HEBBIAN_PROSPECTIVE_HORIZON_VALUE_ITERATION_NO_GO_PROBE__NOT_A_THEOREM",
        "status": "finite-horizon-aut-equivariant-lookahead-cannot-create-d5-selector-from-invariant-future-value",
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
        "equivariant_transition_family": {
            "denominator": GRID_DENOMINATOR,
            "kernel_count": len(kernels),
            "closed_form": "T_a = [[a, 1-a], [1-a, a]]",
            "commutes_with_unit_mirror_swap": True,
            "sample_kernels": [kernel_json(kernels[index]) for index in [0, len(kernels) // 2, -1]],
        },
        "future_reward_grid_audit": reward_audit,
        "fixed_kernel_finite_horizon_audit": fixed_audit,
        "controlled_possible_future_audit": {
            "action_family": "all denominator-12 full-Aut-equivariant kernels",
            "invariant_terminal_value": vector_json(invariant_terminal),
            "invariant_terminal_bellman_trace": [vector_json(value) for value in controlled_invariant_trace],
            "all_controlled_invariant_values_remain_tied": all(value[0] == value[1] for value in controlled_invariant_trace),
            "d5_labelled_terminal_samples": biased_reward_horizon_samples(kernels),
        },
        "exact_proof_certificate": exact_proof_certificate(),
        "selector_interpretation": {
            "question_tested": "Can highest possible future resonance / faith-like finite lookahead select d5 without an added orientation datum?",
            "negative_result": "For invariant terminal future values and full-Aut-equivariant prospective dynamics, every finite horizon remains tied on A1/A5.",
            "conditional_positive_result": "A d5-favoring terminal value or reward can bias the result only because it already names the d5 axis, i.e. imports the missing one-bit unit-axis record.",
            "controlled_future_result": "Maximizing over possible equivariant futures does not help: invariant values stay invariant, and a d5-labelled future is reachable from both mirror axes under the Aut-closed action family.",
            "honest_limit": "This is a finite symmetry/value-iteration no-go, not a theorem deriving prospective value, Hebbian learning, or the unit-axis bit from strict nadsoliton geometry.",
        },
        "ontology_guardrail": {
            "allowed_reading": "The nadsoliton itself may contain prospective self-value about its own possible futures.",
            "forbidden_reading": "No separate informational layer underneath the nadsoliton is introduced to store future values or future probabilities.",
            "preferred_order_preserved": "nadsoliton -> light -> matter -> emergent observer",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is asserted or used.",
            "No legacy physical-role claim is transferred onto K_strict_gate.",
            "No theorem derives a faith-like prospective value functional from strict nadsoliton geometry.",
            "No theorem derives a Hebbian learning law as strict-core dynamics.",
            "No theorem derives the required one-bit unit-axis record from strict core.",
            "Full Aut(Z_12)-equivariant finite-horizon lookahead still forbids singleton d5 selection from invariant future values.",
            "Any d5-favoring future value is classified as a non-strict selector premise unless a new bridge/source theorem is supplied.",
            "No endpoint, arrow orientation, ledger selector, positive lambda action, cycle metric source, anti-Nyquist source, fifth-mode source, future-probability source, or future-value source theorem is claimed.",
            "No QW-2191 discharge and no strict-core selector closure are claimed.",
            "No ToE closure is claimed.",
        ],
        "next_honest_step": "Search for a strict internal source of a non-Aut-invariant future-value functional; otherwise keep highest-future-resonance selection explicitly axiom-augmented/non-strict.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha Hebbian prospective horizon value-iteration no-go probe\n\n"
        "Status: finite-horizon full-Aut-equivariant lookahead cannot create singleton d5 from invariant future value.\n\n"
        f"- Supports scanned: `{len(supports)}`; histogram classes: `{len(histogram_counter)}`.\n"
        f"- Equivariant kernels audited: `{len(kernels)}` with `T_a = [[a, 1-a], [1-a, a]]`.\n"
        f"- Future reward grid: `{reward_audit['reward_count']}` rewards; invariant rewards: `{reward_audit['full_aut_invariant_reward_count']}`; "
        f"invariant+d5-dominant rewards: `{reward_audit['full_aut_invariant_d5_dominant_reward_count']}`.\n"
        f"- Horizons checked: `0..{MAX_HORIZON}`; invariant reward tie violations: `{fixed_audit['invariant_reward_value_tie_violations']}`.\n"
        f"- Controlled invariant Bellman trace remains tied: `{report['controlled_possible_future_audit']['all_controlled_invariant_values_remain_tied']}`.\n"
        f"- Target replay: `q^5={TARGET_Q_POWER}`, eta `{TARGET_ETA}`.\n"
        "- Honest read: highest-future-resonance selects d5 only if the future value already names d5.\n"
        "- No false pass: no future-value source theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
