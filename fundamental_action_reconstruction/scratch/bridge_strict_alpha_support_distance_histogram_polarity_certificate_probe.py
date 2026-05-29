#!/usr/bin/env python3
"""Scratch probe: distance-histogram polarity certificate for strict-alpha support.

The Dirichlet support obstruction showed that the plain strict-gate Dirichlet
candidate clusters 5-node support instead of selecting the fifth-step resonance
support.  This packet makes the obstruction algebraic rather than only a scan.

For a 5-node support S on Z_12, define the internal distance histogram

    h_d(S) = #{{i,j} subset S : dist_Z12(i,j)=d},  d=1..6.

For any translation-invariant indicator Dirichlet action with distance weights
w_d, the support-only energy is

    E_ind(S) = 5*T - 2*sum_d h_d(S)*w_d,
    T = 2*(w_1+...+w_5)+w_6.

Thus support selection is a linear functional of h(S).  For the current strict
weights, the cyclic-contiguous support has h=(4,3,2,1,0,0), while the fifth-step
support has h=(0,3,2,1,4,0).  Their boundary-energy difference is exactly

    E_ind(fifth) - E_ind(contiguous) = 8*(w_1-w_5) > 0.

So the obstruction is not numerical noise: any attractive decreasing weighting
with w_1>w_5 pushes the plain Dirichlet minimizer toward contiguous clustering,
not fifth-step phase support.  To use the fifth-step support one needs a polarity
change (maximize this boundary action / use repulsive or band-pass support cost)
or an additional strict support theorem.
"""
from __future__ import annotations

import itertools
import json
import math
from collections import Counter
from fractions import Fraction
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SUPPORT_OBSTRUCTION = HERE / "bridge_strict_alpha_dirichlet_support_obstruction_certificate_report.json"
OUT_JSON = HERE / "bridge_strict_alpha_support_distance_histogram_polarity_certificate_report.json"
OUT_MD = HERE / "bridge_strict_alpha_support_distance_histogram_polarity_certificate_report.md"

Z12_NODE_COUNT = 12
SUPPORT_SIZE = 5
STRICT_OMEGA = 0.18575
STRICT_PHI = 0.16250
STRICT_BETA = 1.0
STRICT_ETA = 9.0 / 5.0
DENOMINATOR = 3
TARGET_BINARY_EXPONENT = 8
NAD12_SUPPORT_SIZE = 12
ALPHA_SCALE = 4.0
CONTIGUOUS_SUPPORT = (0, 1, 2, 3, 4)
FIFTH_STEP_SUPPORT_ORDERED = (0, 7, 2, 9, 4)
FIFTH_STEP_SUPPORT = tuple(sorted(FIFTH_STEP_SUPPORT_ORDERED))


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def strict_gate_weight(distance: int) -> float:
    return math.cos(STRICT_OMEGA * distance + STRICT_PHI) / (1.0 + STRICT_BETA * (distance**STRICT_ETA))


def cyclic_distance(left: int, right: int) -> int:
    raw = abs(left - right) % Z12_NODE_COUNT
    return min(raw, Z12_NODE_COUNT - raw)


def distance_histogram(support: tuple[int, ...]) -> tuple[int, ...]:
    counts = [0] * (Z12_NODE_COUNT // 2)
    for left, right in itertools.combinations(support, 2):
        counts[cyclic_distance(left, right) - 1] += 1
    return tuple(counts)


def gap_signature(support: tuple[int, ...]) -> tuple[int, ...]:
    nodes = sorted(support)
    return tuple(sorted((nodes[(index + 1) % len(nodes)] - node) % Z12_NODE_COUNT for index, node in enumerate(nodes)))


def eta_from_product(product: Fraction, branch_count: int) -> float:
    correction = float(product) ** (1.0 / branch_count)
    return math.log(NAD12_SUPPORT_SIZE * correction) / math.log(ALPHA_SCALE)


def weight_vector() -> tuple[float, ...]:
    return tuple(strict_gate_weight(distance) for distance in range(1, Z12_NODE_COUNT // 2 + 1))


def internal_weighted_score(histogram: tuple[int, ...], weights: tuple[float, ...]) -> float:
    return sum(count * weight for count, weight in zip(histogram, weights))


def total_node_weight(weights: tuple[float, ...]) -> float:
    return 2.0 * sum(weights[:5]) + weights[5]


def indicator_energy_from_histogram(histogram: tuple[int, ...], weights: tuple[float, ...]) -> float:
    return SUPPORT_SIZE * total_node_weight(weights) - 2.0 * internal_weighted_score(histogram, weights)


def all_support_rows() -> list[dict[str, Any]]:
    weights = weight_vector()
    rows = []
    for support in itertools.combinations(range(Z12_NODE_COUNT), SUPPORT_SIZE):
        hist = distance_histogram(support)
        rows.append(
            {
                "support": list(support),
                "gap_signature": list(gap_signature(support)),
                "distance_histogram": list(hist),
                "internal_weighted_score": internal_weighted_score(hist, weights),
                "indicator_energy_from_histogram": indicator_energy_from_histogram(hist, weights),
            }
        )
    return rows


def histogram_class_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[tuple[int, ...], list[dict[str, Any]]] = {}
    for row in rows:
        grouped.setdefault(tuple(row["distance_histogram"]), []).append(row)
    out = []
    for histogram, members in grouped.items():
        out.append(
            {
                "distance_histogram": list(histogram),
                "support_count": len(members),
                "sample_support": members[0]["support"],
                "sample_gap_signature": members[0]["gap_signature"],
                "internal_weighted_score": members[0]["internal_weighted_score"],
                "indicator_energy_from_histogram": members[0]["indicator_energy_from_histogram"],
            }
        )
    return sorted(out, key=lambda row: (row["indicator_energy_from_histogram"], row["distance_histogram"]))


def fifth_step_orbit_supports() -> list[list[int]]:
    supports = []
    for start in range(Z12_NODE_COUNT):
        support = sorted((start + offset * 7) % Z12_NODE_COUNT for offset in range(SUPPORT_SIZE))
        if support not in supports:
            supports.append(support)
    return supports


def main() -> None:
    support_report = load_json(SUPPORT_OBSTRUCTION)
    weights = weight_vector()
    rows = all_support_rows()
    classes = histogram_class_rows(rows)
    contiguous_hist = distance_histogram(CONTIGUOUS_SUPPORT)
    fifth_hist = distance_histogram(FIFTH_STEP_SUPPORT)
    contiguous_energy = indicator_energy_from_histogram(contiguous_hist, weights)
    fifth_energy = indicator_energy_from_histogram(fifth_hist, weights)
    min_energy = min(row["indicator_energy_from_histogram"] for row in rows)
    max_energy = max(row["indicator_energy_from_histogram"] for row in rows)
    product = Fraction(2 ** TARGET_BINARY_EXPONENT, DENOMINATOR ** SUPPORT_SIZE)
    histogram_counts = Counter(tuple(row["distance_histogram"]) for row in rows)

    report = {
        "status": "OPEN_STRICT_ALPHA_SUPPORT_DISTANCE_HISTOGRAM_POLARITY_CERTIFICATE_NO_SUPPORT_SELECTOR_THEOREM",
        "result_kind": "SCRATCH_STRICT_ALPHA_SUPPORT_DISTANCE_HISTOGRAM_POLARITY_CERTIFICATE_PROBE__NOT_A_THEOREM",
        "source_reports": {
            "dirichlet_support_obstruction_certificate": str(SUPPORT_OBSTRUCTION.relative_to(ROOT)),
        },
        "previous_support_obstruction_replay": {
            "result_kind": support_report["result_kind"],
            "indicator_fifth_step_is_maximum": support_report["obstruction_summary"]["indicator_fifth_step_is_maximum"],
            "balanced_fifth_step_is_maximum_of_support_minima": support_report["obstruction_summary"]["balanced_fifth_step_is_maximum_of_support_minima"],
        },
        "strict_gate_parameters": {
            "kernel": "K_strict_gate(d)=cos(omega*d+phi)/(1+beta*d^eta)",
            "omega": STRICT_OMEGA,
            "phi": STRICT_PHI,
            "beta": STRICT_BETA,
            "eta": STRICT_ETA,
            "note": "strict-side working tuple only; no K_legacy_ont identification is used",
        },
        "target_identity_replay": {
            "q_power_product": f"{product.numerator}/{product.denominator}",
            "eta_residual_vs_9_5": eta_from_product(product, SUPPORT_SIZE) - STRICT_ETA,
        },
        "linear_histogram_identity": {
            "support_size": SUPPORT_SIZE,
            "distance_weights": {f"w_{index + 1}": weight for index, weight in enumerate(weights)},
            "total_node_weight_T": total_node_weight(weights),
            "indicator_energy_formula": "E_ind(S)=5*T-2*sum_d h_d(S)*w_d",
            "histogram_class_count": len(classes),
            "support_count": len(rows),
        },
        "contiguous_vs_fifth_certificate": {
            "contiguous_support": list(CONTIGUOUS_SUPPORT),
            "contiguous_histogram": list(contiguous_hist),
            "contiguous_indicator_energy": contiguous_energy,
            "contiguous_is_global_indicator_minimum": abs(contiguous_energy - min_energy) < 1e-12,
            "fifth_step_support_ordered": list(FIFTH_STEP_SUPPORT_ORDERED),
            "fifth_step_support_canonical": list(FIFTH_STEP_SUPPORT),
            "fifth_step_histogram": list(fifth_hist),
            "fifth_step_indicator_energy": fifth_energy,
            "fifth_step_is_global_indicator_maximum": abs(fifth_energy - max_energy) < 1e-12,
            "histogram_difference_fifth_minus_contiguous": [fifth - contiguous for fifth, contiguous in zip(fifth_hist, contiguous_hist)],
            "exact_two_orbit_gap_formula": "E_ind(fifth)-E_ind(contiguous)=8*(w_1-w_5)",
            "computed_gap": fifth_energy - contiguous_energy,
            "formula_gap": 8.0 * (weights[0] - weights[4]),
            "polarity_condition_for_plain_minimization_to_prefer_fifth_over_contiguous": "requires w_5>w_1 for this pairwise comparison; current strict weights have w_1>w_5",
        },
        "extremal_histogram_classes": {
            "min_indicator_energy": min_energy,
            "max_indicator_energy": max_energy,
            "minimizer_classes": [row for row in classes if abs(row["indicator_energy_from_histogram"] - min_energy) < 1e-12],
            "maximizer_classes": [row for row in classes if abs(row["indicator_energy_from_histogram"] - max_energy) < 1e-12],
            "fifth_step_orbit_supports": fifth_step_orbit_supports(),
            "fifth_step_histogram_support_count": histogram_counts[fifth_hist],
            "contiguous_histogram_support_count": histogram_counts[contiguous_hist],
        },
        "candidate_interpretation": {
            "supported_by_this_probe": True,
            "content": "The support obstruction is exactly a distance-histogram polarity problem: the attractive strict-gate Dirichlet minimizer rewards short internal distances and therefore clusters support.",
            "why_this_is_more_proof_like": "It derives the support-only energy as a linear histogram functional and gives the exact gap E(fifth)-E(contiguous)=8*(w_1-w_5), matching the exhaustive scan.",
            "why_this_is_not_enough": "It does not derive a new support selector; it says the plain attractive Dirichlet polarity must be changed or supplemented before fifth-step resonance support can be selected.",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is exported.",
            "No legacy physical-role transfer is licensed.",
            "This packet explains a support obstruction for the strict-gate Dirichlet candidate; it does not derive a replacement support/phase-placement theorem.",
            "The polarity condition is a model-internal support-functional statement, not a strict nadsoliton action theorem.",
            "No theorem derives branch count m=5, total exponent n=8, denominator 3, or binary-rescale quotient.",
            "No Aut(Z_12)/N462/QW-2191 selector obstruction is discharged.",
            "No theorem derives eta=9/5 without the branch model, denominator/quotient convention, convex ledger selector, and support/phase-placement premise.",
            "No QW-2191 selector discharge and no ToE closure are claimed.",
        ],
        "next_honest_step": "Search for a strict source of the required support polarity change, such as a repulsive or band-pass phase action, or keep fifth-step support as an explicit non-strict support premise.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha support distance-histogram polarity certificate probe\n\n"
        "Status: distance-histogram polarity certificate; no support selector theorem.\n\n"
        f"- Supports scanned: `{len(rows)}`; histogram classes: `{len(classes)}`.\n"
        f"- Contiguous histogram `{list(contiguous_hist)}` has indicator energy `{contiguous_energy:.12f}`.\n"
        f"- Fifth-step histogram `{list(fifth_hist)}` has indicator energy `{fifth_energy:.12f}`.\n"
        f"- Exact gap: `E_fifth-E_contiguous = 8*(w_1-w_5) = {8.0 * (weights[0] - weights[4]):.12f}`.\n"
        f"- Target replay: `q^5={product.numerator}/{product.denominator}`, eta residual `{eta_from_product(product, SUPPORT_SIZE)-STRICT_ETA:.3e}`.\n"
        "- Honest read: plain attractive Dirichlet polarity clusters support; fifth-step support requires a polarity-changing source or explicit support premise.\n"
        "- No false pass: no support/phase-placement theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
