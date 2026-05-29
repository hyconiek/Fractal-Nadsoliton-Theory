#!/usr/bin/env python3
"""Scratch probe: minimal branch-count obstruction for strict alpha -> eta.

The selector-threshold probes made the selector convention explicit, but the
branch count itself was still an input.  This packet asks what the exact radical
correction already forces in a narrow binary/ternary branch model.

Assume only a denominator-3 branch model with m branches and integer binary
numerator exponent total n:

    q_model = (2^n / 3^m)^(1/m) = 2^(n/m) / 3.

Matching the radical target

    q_radical = 2^(8/5) / 3

requires n/m = 8/5, equivalently 5n = 8m.  Since gcd(8,5)=1, every exact branch
model has

    m = 5t, n = 8t.

So five branches are not derived physically here, but they are the irreducible
minimal branch count in this denominator-3 integer-binary model.  Any larger
exact model is a t-fold refinement of the same 5:8 exponent ratio.
"""
from __future__ import annotations

import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RATIONAL_CERT = HERE / "bridge_strict_alpha_selector_rational_threshold_certificate_report.json"
OUT_JSON = HERE / "bridge_strict_alpha_minimal_branch_count_obstruction_report.json"
OUT_MD = HERE / "bridge_strict_alpha_minimal_branch_count_obstruction_report.md"

TARGET_EXPONENT_RATIO = Fraction(8, 5)
TARGET_BRANCH_DENOMINATOR = 3
STRICT_TARGET_ETA = Fraction(9, 5)
NAD12_SUPPORT_SIZE = 12
ALPHA_SCALE = 4.0
MAX_BRANCHES = 40
MAX_TOTAL_BINARY_EXPONENT = 80


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def binomial(n: int, k: int) -> int:
    if k < 0 or k > n:
        return 0
    return math.comb(n, k)


def exact_match(branch_count: int, total_binary_exponent: int) -> bool:
    return Fraction(total_binary_exponent, branch_count) == TARGET_EXPONENT_RATIO


def eta_from_branch_model(branch_count: int, total_binary_exponent: int) -> float:
    q_model = (2 ** total_binary_exponent / TARGET_BRANCH_DENOMINATOR ** branch_count) ** (1.0 / branch_count)
    return math.log(NAD12_SUPPORT_SIZE * q_model) / math.log(ALPHA_SCALE)


def row_for_pair(branch_count: int, total_binary_exponent: int) -> dict[str, Any]:
    ratio = Fraction(total_binary_exponent, branch_count)
    eta = eta_from_branch_model(branch_count, total_binary_exponent)
    return {
        "branch_count_m": branch_count,
        "total_binary_exponent_n": total_binary_exponent,
        "n_over_m": f"{ratio.numerator}/{ratio.denominator}",
        "refinement_factor_t": branch_count // TARGET_EXPONENT_RATIO.denominator if exact_match(branch_count, total_binary_exponent) else None,
        "eta_from_model": eta,
        "eta_residual_vs_9_5": eta - float(STRICT_TARGET_ETA),
        "positive_labelled_compositions_count": binomial(total_binary_exponent - 1, branch_count - 1),
        "positive_compositions_admissible": total_binary_exponent >= branch_count,
    }


def exact_match_rows(max_branches: int, max_total_binary_exponent: int) -> list[dict[str, Any]]:
    rows = []
    for branch_count in range(1, max_branches + 1):
        for total_binary_exponent in range(1, max_total_binary_exponent + 1):
            if exact_match(branch_count, total_binary_exponent):
                rows.append(row_for_pair(branch_count, total_binary_exponent))
    return rows


def near_miss_rows(max_branches: int, max_total_binary_exponent: int, limit: int = 16) -> list[dict[str, Any]]:
    rows = []
    for branch_count in range(1, max_branches + 1):
        for total_binary_exponent in range(1, max_total_binary_exponent + 1):
            if exact_match(branch_count, total_binary_exponent):
                continue
            ratio = Fraction(total_binary_exponent, branch_count)
            distance = abs(float(ratio - TARGET_EXPONENT_RATIO))
            rows.append(
                {
                    **row_for_pair(branch_count, total_binary_exponent),
                    "abs_ratio_distance_vs_8_5": distance,
                }
            )
    return sorted(rows, key=lambda row: row["abs_ratio_distance_vs_8_5"])[:limit]


def main() -> None:
    rational_cert_report = load_json(RATIONAL_CERT)
    matches = exact_match_rows(MAX_BRANCHES, MAX_TOTAL_BINARY_EXPONENT)
    near_misses = near_miss_rows(MAX_BRANCHES, MAX_TOTAL_BINARY_EXPONENT)
    minimal = matches[0]
    all_match_multiples = all(row["branch_count_m"] % TARGET_EXPONENT_RATIO.denominator == 0 for row in matches)
    all_match_total_exponent_multiples = all(row["total_binary_exponent_n"] % TARGET_EXPONENT_RATIO.numerator == 0 for row in matches)

    report = {
        "status": "OPEN_STRICT_ALPHA_MINIMAL_BRANCH_COUNT_OBSTRUCTION_NO_PHYSICAL_BRANCH_THEOREM",
        "result_kind": "SCRATCH_STRICT_ALPHA_MINIMAL_BRANCH_COUNT_OBSTRUCTION_PROBE__NOT_A_THEOREM",
        "source_reports": {
            "selector_rational_threshold_certificate": str(RATIONAL_CERT.relative_to(ROOT)),
        },
        "exact_arithmetic_statement": {
            "branch_model": "q_model=(2^n/3^m)^(1/m)=2^(n/m)/3",
            "target": "q_radical=2^(8/5)/3",
            "required_ratio": "n/m=8/5",
            "diophantine_equation": "5*n=8*m",
            "gcd_8_5": math.gcd(8, 5),
            "general_solution": "m=5*t, n=8*t for positive integer t",
            "minimal_solution": {"branch_count_m": 5, "total_binary_exponent_n": 8, "t": 1},
            "physical_five_branch_derivation_claimed": False,
        },
        "bounded_exact_match_scan": {
            "max_branches": MAX_BRANCHES,
            "max_total_binary_exponent": MAX_TOTAL_BINARY_EXPONENT,
            "exact_match_count": len(matches),
            "exact_match_rows": matches,
            "minimal_exact_match_row": minimal,
            "all_branch_counts_are_multiples_of_5": all_match_multiples,
            "all_total_binary_exponents_are_multiples_of_8": all_match_total_exponent_multiples,
            "near_miss_rows_by_ratio_distance": near_misses,
        },
        "minimal_branch_model_replay": {
            "minimal_branch_count": 5,
            "minimal_total_binary_exponent": 8,
            "minimal_labelled_positive_compositions_count": binomial(7, 4),
            "minimal_canonical_positive_ledgers": [[4, 1, 1, 1, 1], [3, 2, 1, 1, 1], [2, 2, 2, 1, 1]],
            "balanced_ledger_requires_selector_after_minimal_count": True,
        },
        "selector_context_replay": {
            "rational_certificate_result_kind": rational_cert_report["result_kind"],
            "selector_still_needed_after_minimal_branch_count": True,
            "reason": "m=5,n=8 fixes the minimal count and total exponent but still leaves three canonical positive ledgers; a selector is still required to choose (2,2,2,1,1).",
        },
        "candidate_interpretation": {
            "supported_by_this_probe": bool(
                minimal["branch_count_m"] == 5
                and minimal["total_binary_exponent_n"] == 8
                and all_match_multiples
                and all_match_total_exponent_multiples
                and abs(minimal["eta_residual_vs_9_5"]) < 1e-12
            ),
            "content": "In the denominator-3 integer-binary branch model, exact eta=9/5 forces branch counts to be multiples of five; five is the irreducible minimal count.",
            "why_this_is_more_proof_like": "It converts the branch-count question into the Diophantine equation 5n=8m and verifies the bounded exact-match scan.",
            "why_this_is_not_enough": "It does not derive that strict nadsoliton geometry must instantiate the denominator-3 integer-binary branch model, nor does it select the balanced ledger among the minimal canonical ledgers.",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is exported.",
            "No legacy physical-role transfer is licensed.",
            "No theorem derives the denominator-3 integer-binary branch model from strict nadsoliton geometry.",
            "No theorem derives physical five-branch channels; five is only the irreducible arithmetic count inside the stated model.",
            "No selector theorem chooses (2,2,2,1,1) from the minimal m=5,n=8 canonical ledgers.",
            "No theorem derives eta=9/5 without adopting the branch model and a selector convention as extra premises.",
            "No QW-2191 selector discharge and no ToE closure are claimed.",
        ],
        "next_honest_step": "Try to derive the denominator-3 integer-binary branch model from strict nadsoliton geometry; if that is supplied, the arithmetic forces minimal m=5,n=8 but still leaves the selector obligation open.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha minimal branch-count obstruction probe\n\n"
        "Status: irreducible branch-count certificate for eta=9/5; no physical branch theorem.\n\n"
        "- Model rule: `q_model=(2^n/3^m)^(1/m)=2^(n/m)/3`; exact target requires `n/m=8/5`, equivalently `5n=8m`.\n"
        f"- General solution: `m=5*t`, `n=8*t`; bounded scan found `{len(matches)}` exact matches through `m<={MAX_BRANCHES}`, `n<={MAX_TOTAL_BINARY_EXPONENT}`.\n"
        f"- Minimal exact match: `m={minimal['branch_count_m']}`, `n={minimal['total_binary_exponent_n']}`, eta residual `{minimal['eta_residual_vs_9_5']:.3e}`.\n"
        "- Honest read: five branches are arithmetically irreducible inside this branch model, but the model itself and the balanced-ledger selector are still not derived.\n"
        "- No false pass: no denominator-3 branch theorem, no selector discharge, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
