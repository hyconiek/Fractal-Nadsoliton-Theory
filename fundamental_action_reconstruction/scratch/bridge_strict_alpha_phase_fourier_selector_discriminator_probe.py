#!/usr/bin/env python3
"""Scratch probe: phase/Fourier selector discriminator for strict alpha -> eta.

The previous limma audit accepted the arithmetic identity

    q^5 = 256/243 = 2^8/3^5

and rejected the over-strong claim that this is the 12-fifth closure residual.
This packet takes the next honest step: if phase-labelled resonance is used as
the missing selector, what does a concrete Fourier criterion actually select?

Result: phase labelling is useful, but "highest resonance" is not a unique
strict selector.  Parseval gives an exact discriminator: for fixed total binary
weight 8 over five active positive branches, total non-DC Fourier power is a
monotone function of sum(e_i^2).  Therefore a maximum-ripple/highest-power
criterion selects the most unbalanced canonical ledger (4,1,1,1,1), while a
minimum-ripple/coherence-compression criterion selects the balanced ledger
(2,2,2,1,1).  Single-mode scans are placement-sensitive and do not provide an
Aut(Z_12)-invariant selector.

So this is more proof-like than a metaphor, but still not a theorem deriving
eta=9/5: one must explicitly prove the *minimum Fourier ripple / coherence
compression* principle from strict nadsoliton geometry, not merely invoke
resonance language.
"""
from __future__ import annotations

import cmath
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PHASE_AUDIT = HERE / "bridge_strict_alpha_phase_resonance_limmatic_audit_report.json"
ENTROPY_AUDIT = HERE / "bridge_strict_alpha_entropy_selector_discriminator_report.json"
OUT_JSON = HERE / "bridge_strict_alpha_phase_fourier_selector_discriminator_report.json"
OUT_MD = HERE / "bridge_strict_alpha_phase_fourier_selector_discriminator_report.md"

Z12_NODE_COUNT = 12
FIVE_BRANCH_DFT_SIZE = 5
TARGET_BINARY_EXPONENT = 8
STRICT_TARGET_ETA = 9.0 / 5.0
NAD12_SUPPORT_SIZE = 12
ALPHA_SCALE = 4.0
LIMMA = Fraction(256, 243)
FIFTH_STEP = 7
CANONICAL_LEDGERS = [
    (4, 1, 1, 1, 1),
    (3, 2, 1, 1, 1),
    (2, 2, 2, 1, 1),
]


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def eta_from_correction(correction: float) -> float:
    return math.log(NAD12_SUPPORT_SIZE * correction) / math.log(ALPHA_SCALE)


def unique_permutations(values: tuple[int, ...]) -> list[tuple[int, ...]]:
    return sorted(set(itertools.permutations(values)))


def dft_power(values: list[float]) -> list[float]:
    size = len(values)
    rows = []
    for mode in range(size):
        amplitude = sum(values[index] * cmath.exp(-2j * math.pi * mode * index / size) for index in range(size))
        rows.append(float(abs(amplitude) ** 2))
    return rows


def z12_sparse_power(nodes: tuple[int, ...], weights: tuple[int, ...]) -> list[float]:
    values = [0.0] * Z12_NODE_COUNT
    for node, weight in zip(nodes, weights):
        values[node] = float(weight)
    return dft_power(values)


def total_non_dc_power(size: int, ledger: tuple[int, ...]) -> int:
    # Parseval: sum_k |F_k|^2 = N * sum_j |x_j|^2, while the DC power is
    # |sum_j x_j|^2.  Zero-padding to a larger phase ring changes only N.
    return size * sum(value * value for value in ledger) - sum(ledger) ** 2


def ledger_row(ledger: tuple[int, ...]) -> dict[str, Any]:
    product = Fraction(2 ** sum(ledger), 3 ** len(ledger))
    correction = float(product) ** (1.0 / len(ledger))
    five_power = dft_power([float(value) for value in ledger])
    return {
        "ledger": list(ledger),
        "sum": sum(ledger),
        "sum_squares": sum(value * value for value in ledger),
        "five_branch_total_non_dc_power_parseval": total_non_dc_power(FIVE_BRANCH_DFT_SIZE, ledger),
        "z12_zero_padded_total_non_dc_power_parseval": total_non_dc_power(Z12_NODE_COUNT, ledger),
        "five_branch_non_dc_powers": five_power[1:],
        "five_branch_max_single_non_dc_power": max(five_power[1:]),
        "product_fraction": f"{product.numerator}/{product.denominator}",
        "eta_from_product": eta_from_correction(correction),
        "eta_residual_vs_9_5": eta_from_correction(correction) - STRICT_TARGET_ETA,
    }


def fixed_fifth_window_nodes() -> tuple[int, ...]:
    current = 0
    out = []
    for _ in range(5):
        out.append(current)
        current = (current + FIFTH_STEP) % Z12_NODE_COUNT
    return tuple(out)


def fixed_window_single_mode_summary(ledger: tuple[int, ...]) -> dict[str, Any]:
    nodes = fixed_fifth_window_nodes()
    maxima = []
    best_rows = []
    for weights in unique_permutations(ledger):
        powers = z12_sparse_power(nodes, weights)
        non_dc = powers[1:]
        max_power = max(non_dc)
        mode = 1 + non_dc.index(max_power)
        maxima.append(max_power)
        best_rows.append({"weights_on_fifth_window": list(weights), "mode": mode, "max_single_non_dc_power": max_power})
    min_row = min(best_rows, key=lambda row: row["max_single_non_dc_power"])
    max_row = max(best_rows, key=lambda row: row["max_single_non_dc_power"])
    return {
        "ledger": list(ledger),
        "fifth_window_nodes": list(nodes),
        "permutation_count": len(maxima),
        "min_over_label_permutations": min(maxima),
        "max_over_label_permutations": max(maxima),
        "mean_over_label_permutations": sum(maxima) / len(maxima),
        "minimizing_assignment": min_row,
        "maximizing_assignment": max_row,
    }


def bounded_z12_placement_scan(ledger: tuple[int, ...]) -> dict[str, Any]:
    max_single_values = []
    total_rows = 0
    best = None
    worst = None
    for nodes in itertools.combinations(range(Z12_NODE_COUNT), len(ledger)):
        for weights in unique_permutations(ledger):
            powers = z12_sparse_power(nodes, weights)
            non_dc = powers[1:]
            max_power = max(non_dc)
            mode = 1 + non_dc.index(max_power)
            row = {"nodes": list(nodes), "weights": list(weights), "mode": mode, "max_single_non_dc_power": max_power}
            max_single_values.append(max_power)
            total_rows += 1
            if best is None or max_power < best["max_single_non_dc_power"]:
                best = row
            if worst is None or max_power > worst["max_single_non_dc_power"]:
                worst = row
    return {
        "ledger": list(ledger),
        "rows_scanned": total_rows,
        "min_possible_max_single_non_dc_power": best["max_single_non_dc_power"],
        "max_possible_max_single_non_dc_power": worst["max_single_non_dc_power"],
        "best_min_ripple_placement": best,
        "worst_highest_resonance_placement": worst,
        "distinct_rounded_max_single_values": len({round(value, 12) for value in max_single_values}),
    }


def argmin(rows: list[dict[str, Any]], key: str) -> dict[str, Any]:
    return min(rows, key=lambda row: row[key])


def argmax(rows: list[dict[str, Any]], key: str) -> dict[str, Any]:
    return max(rows, key=lambda row: row[key])


def main() -> None:
    phase_report = load_json(PHASE_AUDIT)
    entropy_report = load_json(ENTROPY_AUDIT)
    ledger_rows = [ledger_row(ledger) for ledger in CANONICAL_LEDGERS]
    fixed_window_rows = [fixed_window_single_mode_summary(ledger) for ledger in CANONICAL_LEDGERS]
    placement_scan_rows = [bounded_z12_placement_scan(ledger) for ledger in CANONICAL_LEDGERS]

    min_ripple_winner = argmin(ledger_rows, "five_branch_total_non_dc_power_parseval")
    max_ripple_winner = argmax(ledger_rows, "five_branch_total_non_dc_power_parseval")
    z12_min_ripple_winner = argmin(ledger_rows, "z12_zero_padded_total_non_dc_power_parseval")
    z12_max_ripple_winner = argmax(ledger_rows, "z12_zero_padded_total_non_dc_power_parseval")

    report = {
        "status": "OPEN_STRICT_ALPHA_PHASE_FOURIER_SELECTOR_DISCRIMINATOR_NO_SELECTOR_THEOREM",
        "result_kind": "SCRATCH_STRICT_ALPHA_PHASE_FOURIER_SELECTOR_DISCRIMINATOR_PROBE__NOT_A_THEOREM",
        "source_reports": {
            "phase_resonance_limmatic_audit": str(PHASE_AUDIT.relative_to(ROOT)),
            "entropy_selector_discriminator": str(ENTROPY_AUDIT.relative_to(ROOT)),
        },
        "target_identity_replay": {
            "q_radical_power_5": f"{LIMMA.numerator}/{LIMMA.denominator}",
            "all_canonical_ledgers_share_product": "2^8/3^5 = 256/243",
            "eta_target": STRICT_TARGET_ETA,
            "phase_audit_result_kind": phase_report["result_kind"],
            "entropy_audit_result_kind": entropy_report["result_kind"],
        },
        "parseval_discriminator": {
            "identity": "total_non_dc_power = N*sum(e_i^2) - (sum e_i)^2",
            "fixed_total_binary_units": TARGET_BINARY_EXPONENT,
            "ledger_rows": ledger_rows,
            "minimum_fourier_ripple_winner": min_ripple_winner,
            "maximum_fourier_ripple_winner": max_ripple_winner,
            "z12_zero_padded_minimum_fourier_ripple_winner": z12_min_ripple_winner,
            "z12_zero_padded_maximum_fourier_ripple_winner": z12_max_ripple_winner,
            "verdict": "BALANCED_LEDGER_IS_MINIMUM_RIPPLE_NOT_HIGHEST_TOTAL_RESONANCE_POWER",
        },
        "fixed_fifth_window_single_mode_scan": {
            "interpretation": "Use the first five nodes of the Z12 fifth-step orbit as phase-labelled branch slots; scan all label permutations.",
            "rows": fixed_window_rows,
            "verdict": "SINGLE_MODE_CRITERIA_DEPEND_ON_ASSIGNMENT_AND_NEED_AN_EXTRA_PLACEMENT_OR_COHERENCE_PRINCIPLE",
        },
        "bounded_z12_placement_scan": {
            "interpretation": "Scan every 5-node support in Z12 and every distinct assignment of each ledger; this is a finite obstruction check, not a physical theorem.",
            "rows": placement_scan_rows,
            "all_ledgers_can_reach_dc_sized_single_mode_power_64": all(abs(row["max_possible_max_single_non_dc_power"] - 64.0) < 1e-9 for row in placement_scan_rows),
            "verdict": "HIGHEST_SINGLE_MODE_RESONANCE_IS_PLACEMENT_DOMINATED_AND_DOES_NOT_SELECT_THE_BALANCED_LEDGER",
        },
        "candidate_interpretation": {
            "supported_by_this_probe": True,
            "content": "Fourier arithmetic gives a precise selector discriminator: balanced (2,2,2,1,1) is selected by minimum ripple / coherence compression, not by naive maximum resonance power.",
            "why_this_is_more_proof_like": "The Parseval identity is exact and independent of numerical fitting; bounded scans show that single-mode resonance language is underdetermined by phase labels alone.",
            "why_this_is_not_enough": "No strict theorem derives the minimum-ripple/coherence-compression variational principle, the branch support, or the labelled placement convention from strict nadsoliton geometry.",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is exported.",
            "No legacy physical-role transfer is licensed.",
            "No theorem derives a strict-core phase/Fourier selector from nadsoliton geometry.",
            "Naive highest total Fourier resonance selects (4,1,1,1,1), not (2,2,2,1,1).",
            "Balanced (2,2,2,1,1) is supported only after adding a minimum-ripple/coherence-compression premise.",
            "Single-mode resonance scans are placement-sensitive and do not discharge Aut(Z_12)/N462/QW-2191 selector obstructions.",
            "No theorem derives eta=9/5 without adopting branch count, denominator/quotient, phase labelling, and minimum-ripple selector premises.",
            "No QW-2191 selector discharge and no ToE closure are claimed.",
        ],
        "next_honest_step": "Try to derive the minimum Fourier ripple / coherence-compression variational principle as a strict nadsoliton theorem; otherwise record it as an explicit non-strict selector premise.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha phase Fourier selector discriminator probe\n\n"
        "Status: Fourier selector discriminator for eta=9/5; no strict selector theorem.\n\n"
        "- Exact Parseval discriminator: `total_non_dc_power = N*sum(e_i^2) - (sum e_i)^2`.\n"
        f"- Minimum-ripple winner: `{min_ripple_winner['ledger']}`; this matches the balanced ledger.\n"
        f"- Maximum total-resonance-power winner: `{max_ripple_winner['ledger']}`; this is not the balanced ledger.\n"
        "- Bounded Z12 placement scan: every canonical ledger can reach single-mode power `64`, so naive highest single-mode resonance is placement-dominated.\n"
        "- Honest read: phase/Fourier language supports `(2,2,2,1,1)` only as a minimum-ripple/coherence-compression selector, not as an automatic highest-resonance selector.\n"
        "- No false pass: no phase/Fourier strict selector theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
