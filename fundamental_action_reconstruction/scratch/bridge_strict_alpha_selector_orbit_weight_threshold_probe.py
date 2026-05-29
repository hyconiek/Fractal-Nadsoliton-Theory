#!/usr/bin/env python3
"""Scratch probe: orbit-weight threshold for strict alpha -> eta selector.

The entropy discriminator and Robin-Hood certificate showed two conditional facts:

1. fixed-labelled multiplicity W(e)=8!/prod_i(e_i!) selects the balanced ledger
   (2,2,2,1,1), and
2. full unlabelled orbit aggregation O(e)*W(e) selects (3,2,1,1,1).

This probe inserts a one-parameter selector family between those conventions:

    Score_gamma(e) = log W(e) + gamma * log O(e),

where O(e) is the canonical orbit size of the multiset e.  The exact switch
threshold between the balanced ledger B=(2,2,2,1,1) and the orbit-favoured ledger
C=(3,2,1,1,1) is

    gamma_c = log(W_B/W_C) / log(O_C/O_B) = log(3/2) / log(2).

Thus fixed-labelled selection is not merely different from orbit aggregation;
it is quantitatively stable only for gamma < gamma_c ~= 0.58496.  For gamma >
gamma_c, the selector flips to C.  At gamma=gamma_c the selector is tied.

This is not a strict selector theorem.  It sharpens the missing proof obligation:
derive either gamma=0 / fixed-labelled channels, or prove an intrinsic strict
orbit-weight exponent below gamma_c.  Otherwise the eta=9/5 balanced-ledger route
remains selector-convention sensitive.
"""
from __future__ import annotations

import json
import math
from collections import Counter
from fractions import Fraction
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
ROBIN_HOOD = HERE / "bridge_strict_alpha_entropy_robin_hood_certificate_report.json"
OUT_JSON = HERE / "bridge_strict_alpha_selector_orbit_weight_threshold_report.json"
OUT_MD = HERE / "bridge_strict_alpha_selector_orbit_weight_threshold_report.md"

LEDGER_ORBIT_LOW = (4, 1, 1, 1, 1)
LEDGER_ORBIT_FAVOURED = (3, 2, 1, 1, 1)
LEDGER_BALANCED = (2, 2, 2, 1, 1)
CANONICAL_LEDGERS = [LEDGER_ORBIT_LOW, LEDGER_ORBIT_FAVOURED, LEDGER_BALANCED]
STRICT_TARGET_ETA = 9.0 / 5.0
NAD12_SUPPORT_SIZE = 12
ALPHA_SCALE = 4.0


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def labelled_multinomial_count(ledger: tuple[int, ...]) -> int:
    numerator = math.factorial(sum(ledger))
    denominator = math.prod(math.factorial(e) for e in ledger)
    return numerator // denominator


def orbit_size(ledger: tuple[int, ...]) -> int:
    multiplicities = Counter(ledger).values()
    denominator = math.prod(math.factorial(m) for m in multiplicities)
    return math.factorial(len(ledger)) // denominator


def selector_score(ledger: tuple[int, ...], gamma: float) -> float:
    return math.log(labelled_multinomial_count(ledger)) + gamma * math.log(orbit_size(ledger))


def product_fraction(ledger: tuple[int, ...]) -> Fraction:
    return Fraction(2 ** sum(ledger), 3 ** len(ledger))


def eta_from_ledger(ledger: tuple[int, ...]) -> float:
    correction = float(product_fraction(ledger)) ** (1.0 / len(ledger))
    return math.log(NAD12_SUPPORT_SIZE * correction) / math.log(ALPHA_SCALE)


def row_for_ledger(ledger: tuple[int, ...], gamma_values: list[float]) -> dict[str, Any]:
    weights = labelled_multinomial_count(ledger)
    orbit = orbit_size(ledger)
    product = product_fraction(ledger)
    return {
        "ledger": list(ledger),
        "labelled_multinomial_count_W": weights,
        "orbit_size_O": orbit,
        "full_orbit_aggregate_O_times_W": orbit * weights,
        "product_fraction": f"{product.numerator}/{product.denominator}",
        "eta_from_shared_product": eta_from_ledger(ledger),
        "eta_residual_vs_9_5": eta_from_ledger(ledger) - STRICT_TARGET_ETA,
        "scores_by_gamma": {f"{gamma:.6f}": selector_score(ledger, gamma) for gamma in gamma_values},
    }


def winner_for_gamma(gamma: float) -> tuple[int, ...]:
    return max(CANONICAL_LEDGERS, key=lambda ledger: selector_score(ledger, gamma))


def score_gap_balanced_minus_orbit_favoured(gamma: float) -> float:
    return selector_score(LEDGER_BALANCED, gamma) - selector_score(LEDGER_ORBIT_FAVOURED, gamma)


def main() -> None:
    robin_hood_report = load_json(ROBIN_HOOD)
    gamma_critical = math.log(Fraction(3, 2)) / math.log(2)
    probe_gammas = [0.0, 0.5, gamma_critical, 0.75, 1.0]
    rows = [row_for_ledger(ledger, probe_gammas) for ledger in CANONICAL_LEDGERS]
    winners = {f"{gamma:.6f}": list(winner_for_gamma(gamma)) for gamma in probe_gammas}

    just_below = gamma_critical - 1e-6
    just_above = gamma_critical + 1e-6
    report = {
        "status": "OPEN_STRICT_ALPHA_SELECTOR_ORBIT_WEIGHT_THRESHOLD_NO_STRICT_SELECTOR_THEOREM",
        "result_kind": "SCRATCH_STRICT_ALPHA_SELECTOR_ORBIT_WEIGHT_THRESHOLD_PROBE__NOT_A_THEOREM",
        "source_reports": {
            "entropy_robin_hood_certificate": str(ROBIN_HOOD.relative_to(ROOT)),
        },
        "selector_family": {
            "score_definition": "Score_gamma(e)=log W(e)+gamma*log O(e)",
            "W": "fixed-labelled bit-assignment count 8!/prod_i(e_i!)",
            "O": "canonical orbit size of the exponent multiset",
            "gamma_0_reading": "fixed-labelled entropy/multinomial selector",
            "gamma_1_reading": "full unlabelled orbit-aggregate selector",
        },
        "candidate_rows": rows,
        "critical_threshold": {
            "balanced_ledger_B": list(LEDGER_BALANCED),
            "orbit_favoured_ledger_C": list(LEDGER_ORBIT_FAVOURED),
            "W_B_over_W_C": "3/2",
            "O_C_over_O_B": "2/1",
            "gamma_critical_closed_form": "log(3/2)/log(2)",
            "gamma_critical_numeric": gamma_critical,
            "score_gap_B_minus_C_at_gamma_0": score_gap_balanced_minus_orbit_favoured(0.0),
            "score_gap_B_minus_C_at_gamma_critical": score_gap_balanced_minus_orbit_favoured(gamma_critical),
            "score_gap_B_minus_C_at_gamma_1": score_gap_balanced_minus_orbit_favoured(1.0),
            "winner_just_below_threshold": list(winner_for_gamma(just_below)),
            "winner_just_above_threshold": list(winner_for_gamma(just_above)),
        },
        "winner_scan": {
            "probe_gammas": probe_gammas,
            "winners_by_gamma": winners,
            "balanced_wins_for_gamma_less_than_gamma_c": True,
            "orbit_favoured_wins_for_gamma_greater_than_gamma_c": True,
            "tie_at_gamma_c_between_B_and_C": abs(score_gap_balanced_minus_orbit_favoured(gamma_critical)) < 1e-12,
        },
        "candidate_interpretation": {
            "supported_by_this_probe": bool(
                winner_for_gamma(0.0) == LEDGER_BALANCED
                and winner_for_gamma(1.0) == LEDGER_ORBIT_FAVOURED
                and winner_for_gamma(just_below) == LEDGER_BALANCED
                and winner_for_gamma(just_above) == LEDGER_ORBIT_FAVOURED
                and abs(score_gap_balanced_minus_orbit_favoured(gamma_critical)) < 1e-12
                and robin_hood_report["certificate_summary"]["balanced_terminal"] == list(LEDGER_BALANCED)
            ),
            "content": "The balanced-ledger selector is stable only below the exact orbit-weight threshold gamma_c=log(3/2)/log(2).",
            "why_this_is_more_proof_like": "It turns the labelled-vs-orbit selector ambiguity into an exact bifurcation threshold for a one-parameter selector family.",
            "why_this_is_not_enough": "No strict theorem derives gamma=0 or any intrinsic orbit-weight exponent gamma<gamma_c from nadsoliton geometry.",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is exported.",
            "No legacy physical-role transfer is licensed.",
            "No theorem derives fixed-labelled branch entropy as the strict-core selector.",
            "No theorem derives an intrinsic orbit-weight exponent gamma below log(3/2)/log(2).",
            "For gamma>log(3/2)/log(2), the selector chooses (3,2,1,1,1), not the balanced eta-ledger.",
            "No theorem derives eta=9/5 without adopting the branch count, ternary normalization, and a below-threshold selector convention as extra premises.",
            "No QW-2191 selector discharge and no ToE closure are claimed.",
        ],
        "next_honest_step": "Try to derive fixed-labelled channels or a below-threshold orbit-weight exponent from strict nadsoliton geometry; otherwise record the balanced eta-ledger as selector-sensitive.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha selector orbit-weight threshold probe\n\n"
        "Status: exact selector-threshold discriminator for eta=9/5; no strict selector theorem.\n\n"
        "- Selector family: `Score_gamma(e)=log W(e)+gamma*log O(e)`, interpolating fixed-labelled (`gamma=0`) and orbit-aggregate (`gamma=1`) conventions.\n"
        f"- Critical threshold: `gamma_c=log(3/2)/log(2)={gamma_critical:.12f}`.\n"
        f"- Winners: below threshold `{list(winner_for_gamma(just_below))}`; above threshold `{list(winner_for_gamma(just_above))}`; at `gamma=1`, `{list(winner_for_gamma(1.0))}`.\n"
        "- Honest read: the balanced eta-ledger is stable only for a below-threshold selector convention; this is a quantitative selector proof obligation, not closure.\n"
        "- No false pass: no strict derivation of gamma=0 or gamma<gamma_c, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
