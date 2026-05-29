#!/usr/bin/env python3
"""Scratch probe: Hebbian resonance-learning audit for highest-resonance claims.

The new suggestion is that the nadsoliton may behave like a self-learning
Hebbian network and that this could affect the recent highest-resonance/mode
lock probes.  This file makes the suggestion finite and checkable on Z_12.

Finite model used here:
- a 5-node support S is a binary activity pattern on 12 nodes;
- Hebbian learning stores W_D = sum_{S in D} x_S x_S^T (or its centered
  covariance; for non-DC Fourier modes the same spectral ranking is obtained);
- because translated orbit data make W_D circulant, the learned non-DC
  resonance channels are exactly ranked by averaged coherent powers P_k(S).

Honest finite result:
- Hebbian learning can stabilize/amplify a resonance channel already present in
  the formative record.
- A d5/fifth-step training orbit learns k=5 as the unique leading non-DC
  channel; a contiguous orbit learns k=1; a parity/Nyquist orbit learns k=6.
- Unbiased training on all 792 five-node supports learns no channel at all: all
  non-DC Fourier channels tie.
- Across individual supports, Hebbian winner basins are broad; only 84/792
  supports have unique k=5 dominance and 144/792 include k=5 among tied or
  unique winners.  Therefore Hebbian learning is not, by itself, a strict proof
  of a fifth-mode lock or a discharge of QW-2191.
"""
from __future__ import annotations

import json
import math
from collections import Counter
from itertools import combinations
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
OUT_JSON = HERE / "bridge_strict_alpha_hebbian_resonance_learning_audit_report.json"
OUT_MD = HERE / "bridge_strict_alpha_hebbian_resonance_learning_audit_report.md"

N = 12
SUPPORT_SIZE = 5
TARGET_MODE = 5
INDEPENDENT_MODES = [1, 2, 3, 4, 5, 6]
TARGET_Q_POWER = "256/243"
TARGET_ETA = "9/5"
FLOAT_TOL = 1e-9


def support_list() -> list[tuple[int, ...]]:
    return list(combinations(range(N), SUPPORT_SIZE))


def canonical_support(support: tuple[int, ...] | list[int]) -> tuple[int, ...]:
    return tuple(sorted(int(x) % N for x in support))


def orbit_supports(step: int) -> list[tuple[int, ...]]:
    return sorted(
        {
            canonical_support([(start + step * index) % N for index in range(SUPPORT_SIZE)])
            for start in range(N)
        }
    )


def parity_nyquist_supports(supports: list[tuple[int, ...]]) -> list[tuple[int, ...]]:
    return [support for support in supports if all(node % 2 == support[0] % 2 for node in support)]


def coherent_power(support: tuple[int, ...], mode: int) -> float:
    real = sum(math.cos(2.0 * math.pi * mode * node / N) for node in support)
    imag = sum(math.sin(2.0 * math.pi * mode * node / N) for node in support)
    return real * real + imag * imag


def rounded(value: float) -> float:
    return round(value, 12)


def averaged_mode_scores(dataset: list[tuple[int, ...]]) -> dict[int, float]:
    return {
        mode: sum(coherent_power(support, mode) for support in dataset) / len(dataset)
        for mode in INDEPENDENT_MODES
    }


def leading_modes(scores: dict[int, float]) -> list[int]:
    maximum = max(scores.values())
    return [mode for mode in INDEPENDENT_MODES if abs(scores[mode] - maximum) <= FLOAT_TOL]


def hebbian_dataset_certificate(name: str, dataset: list[tuple[int, ...]]) -> dict[str, Any]:
    scores = averaged_mode_scores(dataset)
    winners = leading_modes(scores)
    return {
        "dataset": name,
        "pattern_count": len(dataset),
        "hebbian_rule": "W_D=sum_{S in D} x_S x_S^T; non-DC Fourier ranking equals average P_k(S)",
        "centered_covariance_note": "Subtracting the mean activity only changes the DC channel, so the non-DC winner list is unchanged in this Z_12 orbit audit.",
        "average_non_dc_mode_scores": {str(mode): rounded(scores[mode]) for mode in INDEPENDENT_MODES},
        "leading_non_dc_modes": winners,
        "unique_leading_mode": winners[0] if len(winners) == 1 else None,
        "sample_patterns": [list(row) for row in dataset[:12]],
    }


def individual_support_basin_audit(supports: list[tuple[int, ...]]) -> dict[str, Any]:
    winner_counter: Counter[tuple[int, ...]] = Counter()
    unique_counter: Counter[int] = Counter()
    tied_counter: Counter[tuple[int, ...]] = Counter()
    examples_by_winner: dict[tuple[int, ...], list[list[int]]] = {}

    for support in supports:
        powers = {mode: coherent_power(support, mode) for mode in INDEPENDENT_MODES}
        maximum = max(powers.values())
        winners = tuple(mode for mode in INDEPENDENT_MODES if abs(powers[mode] - maximum) <= FLOAT_TOL)
        winner_counter[winners] += 1
        examples_by_winner.setdefault(winners, [])
        if len(examples_by_winner[winners]) < 5:
            examples_by_winner[winners].append(list(support))
        if len(winners) == 1:
            unique_counter[winners[0]] += 1
        else:
            tied_counter[winners] += 1

    target_including_count = sum(count for winners, count in winner_counter.items() if TARGET_MODE in winners)
    return {
        "support_count": len(supports),
        "winner_class_count": len(winner_counter),
        "unique_winner_count": sum(unique_counter.values()),
        "tied_winner_count": sum(tied_counter.values()),
        "unique_winner_counts_by_mode": {str(mode): unique_counter[mode] for mode in INDEPENDENT_MODES},
        "tied_winner_classes": {"/".join(map(str, key)): value for key, value in sorted(tied_counter.items())},
        "all_winner_classes": {"/".join(map(str, key)): value for key, value in sorted(winner_counter.items())},
        "target_mode_unique_count": unique_counter[TARGET_MODE],
        "target_mode_including_ties_count": target_including_count,
        "target_mode_unique_fraction": f"{unique_counter[TARGET_MODE]}/{len(supports)}",
        "target_mode_including_ties_fraction": f"{target_including_count}/{len(supports)}",
        "examples_by_winner_class": {"/".join(map(str, key)): value for key, value in sorted(examples_by_winner.items())},
        "interpretation": "Single-pattern Hebbian learning has broad basins; k=5 is one basin among several, not a strict selector by itself.",
    }


def main() -> None:
    supports = support_list()
    d5_orbit = orbit_supports(TARGET_MODE)
    contiguous_orbit = orbit_supports(1)
    nyquist_orbit = parity_nyquist_supports(supports)

    all_certificate = hebbian_dataset_certificate("all_5_node_supports_unbiased", supports)
    d5_certificate = hebbian_dataset_certificate("translated_fifth_step_d5_orbit", d5_orbit)
    contiguous_certificate = hebbian_dataset_certificate("translated_contiguous_orbit", contiguous_orbit)
    nyquist_certificate = hebbian_dataset_certificate("parity_nyquist_orbit", nyquist_orbit)
    basin_audit = individual_support_basin_audit(supports)

    report: dict[str, Any] = {
        "result_kind": "SCRATCH_STRICT_ALPHA_HEBBIAN_RESONANCE_LEARNING_AUDIT_PROBE__NOT_A_THEOREM",
        "status": "hebbian-amplification-is-conditional-not-a-fifth-mode-source-theorem",
        "finite_model": {
            "ring": "Z_12",
            "support_size": SUPPORT_SIZE,
            "support_count": len(supports),
            "activity_pattern": "binary vector x_S in {0,1}^12 with five active nodes",
            "hebbian_weight": "W_D=sum_{S in D} x_S x_S^T; centered covariance has the same non-DC Fourier ranking",
            "resonance_score": "P_k(S)=|sum_{x in S} exp(2*pi*i*k*x/12)|^2",
            "independent_non_dc_modes": INDEPENDENT_MODES,
            "target_q_power": TARGET_Q_POWER,
            "target_eta": TARGET_ETA,
        },
        "finite_hebbian_spectral_lemma": {
            "statement": "For a translated orbit dataset on Z_12, the Hebbian matrix is circulant, so Fourier modes diagonalize the learned weights; non-DC learned-channel ranking is the dataset-average coherent power ranking.",
            "proof_status": "finite computational certificate plus standard circulant diagonalization; not a strict nadsoliton theorem",
            "why_this_matters": "Hebbian learning can explain persistence/amplification of a recorded resonance channel, but cannot choose which teacher/formative channel exists without extra input.",
        },
        "dataset_certificates": {
            "unbiased_all_supports": all_certificate,
            "fifth_step_d5_orbit": d5_certificate,
            "contiguous_orbit": contiguous_certificate,
            "parity_nyquist_orbit": nyquist_certificate,
        },
        "individual_support_basin_audit": basin_audit,
        "impact_on_recent_highest_resonance_work": {
            "positive_impact": "If a self-recorded formative trace already supplies the fifth-step/d5 orbit or a k=5 teacher channel, Hebbian resonance learning gives a concrete finite mechanism that locks/amplifies k=5.",
            "negative_result": "Hebbian learning does not create the fifth-mode premise from an unbiased start: all five-node supports tie across non-DC modes on average, and other teacher orbits learn k=1 or k=6 instead.",
            "relation_to_mode_lock_probe": "The previous highest-resonance mode-lock certificate can be re-read as a Hebbian readout after a fifth-mode teacher record, not as a derivation of that teacher record.",
            "selector_status": "conditional mechanism candidate; no strict selector theorem and no QW-2191 discharge",
        },
        "ontology_guardrail": {
            "allowed_reading": "The nadsoliton itself is treated as primordial information in a solitonic state that can contain finite self-recorded formation patterns and update/retain resonant correlations.",
            "forbidden_reading": "The Hebbian network is not introduced as a separate informational layer underneath the nadsoliton.",
            "preferred_order_preserved": "nadsoliton -> light -> matter -> emergent observer",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is asserted or used.",
            "No legacy physical-role claim is transferred onto K_strict_gate.",
            "Hebbian learning is modeled only as an internal self-record/update mechanism candidate of the nadsoliton, not as a sub-nadsoliton ontology.",
            "No theorem derives the fifth-mode teacher trace, d5 orbit, endpoint, arrow orientation, ledger selector, or positive lambda action from strict geometry.",
            "Unbiased Hebbian training over all five-node supports gives a complete non-DC tie, not a k=5 selector.",
            "Other plausible teacher records select other modes: contiguous selects k=1 and parity/Nyquist selects k=6.",
            "No QW-2191 discharge and no strict-core selector closure are claimed.",
            "No ToE closure is claimed.",
        ],
        "next_honest_step": "Look for a strict-side source of the fifth-mode teacher/self-record trace; absent that, record a non-bridge theorem that Hebbian learning only amplifies already supplied resonance data.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha Hebbian resonance-learning audit probe\n\n"
        "Status: Hebbian learning is a conditional amplification mechanism, not a fifth-mode source theorem.\n\n"
        f"- Supports scanned: `{len(supports)}` five-node supports on `Z_12`.\n"
        f"- Unbiased all-support Hebbian leading modes: `{all_certificate['leading_non_dc_modes']}` with equal score `{all_certificate['average_non_dc_mode_scores']['1']}`.\n"
        f"- d5-orbit Hebbian leading modes: `{d5_certificate['leading_non_dc_modes']}`.\n"
        f"- Contiguous-orbit Hebbian leading modes: `{contiguous_certificate['leading_non_dc_modes']}`.\n"
        f"- Parity/Nyquist-orbit Hebbian leading modes: `{nyquist_certificate['leading_non_dc_modes']}`.\n"
        f"- Individual support unique k=5 basin: `{basin_audit['target_mode_unique_fraction']}`; k=5 including ties: `{basin_audit['target_mode_including_ties_fraction']}`.\n"
        f"- Target replay: `q^5={TARGET_Q_POWER}`, eta `{TARGET_ETA}`.\n"
        "- Honest read: Hebbian learning can lock/amplify a self-recorded fifth-mode teacher trace, but does not derive that trace.\n"
        "- No false pass: no fifth-mode theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
