#!/usr/bin/env python3
"""Scratch probe: primitive denominator obstruction for strict alpha -> eta.

The minimal branch-count probe showed that, inside a denominator-3 integer-binary
branch model, exact eta=9/5 forces the Diophantine family

    m = 5*t, n = 8*t.

This packet asks whether the denominator 3 itself is forced by the exact radical
target.  In a wider integer-denominator branch model

    q_model = (2^n / D^m)^(1/m) = 2^(n/m) / D,

exact equality with q_radical=2^(8/5)/3 gives

    D = 3 * 2^(n/m - 8/5).

If D is an integer and n/m is rational, then 2^(n/m-8/5) must be rational; by
prime factorization this requires n/m-8/5 to be an integer k.  Hence exact
integer-denominator models are not unique:

    D = 3*2^k,   m = 5*t,   n = (8+5*k)*t,   k>=0.

So D=3 is the primitive odd / no-hidden-binary-rescale denominator, not a theorem
from the radical target alone.  This is a sharpening obstruction: any proof that
uses denominator 3 must derive a primitive/no-binary-rescale convention, or admit
an infinite exact family D=3,6,12,... with shifted binary numerator load.
"""
from __future__ import annotations

import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
MINIMAL_BRANCH = HERE / "bridge_strict_alpha_minimal_branch_count_obstruction_report.json"
OUT_JSON = HERE / "bridge_strict_alpha_primitive_denominator_obstruction_report.json"
OUT_MD = HERE / "bridge_strict_alpha_primitive_denominator_obstruction_report.md"

MAX_K = 6
MAX_T = 8
STRICT_TARGET_ETA = Fraction(9, 5)
NAD12_SUPPORT_SIZE = 12
ALPHA_SCALE = 4.0


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def eta_from_integer_denominator_model(branch_count: int, total_binary_exponent: int, denominator: int) -> float:
    q_model = (2 ** total_binary_exponent / denominator ** branch_count) ** (1.0 / branch_count)
    return math.log(NAD12_SUPPORT_SIZE * q_model) / math.log(ALPHA_SCALE)


def exact_family_row(k: int, t: int) -> dict[str, Any]:
    branch_count = 5 * t
    total_binary_exponent = (8 + 5 * k) * t
    denominator = 3 * (2 ** k)
    eta = eta_from_integer_denominator_model(branch_count, total_binary_exponent, denominator)
    return {
        "binary_rescale_k": k,
        "refinement_t": t,
        "denominator_D": denominator,
        "branch_count_m": branch_count,
        "total_binary_exponent_n": total_binary_exponent,
        "n_over_m": f"{8 + 5 * k}/5",
        "D_label": f"3*2^{k}",
        "primitive_no_binary_rescale": k == 0,
        "eta_from_model": eta,
        "eta_residual_vs_9_5": eta - float(STRICT_TARGET_ETA),
    }


def family_rows(max_k: int, max_t: int) -> list[dict[str, Any]]:
    return [exact_family_row(k, t) for k in range(max_k + 1) for t in range(1, max_t + 1)]


def denominator_summary(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    summary = []
    for k in range(MAX_K + 1):
        k_rows = [row for row in rows if row["binary_rescale_k"] == k]
        summary.append(
            {
                "binary_rescale_k": k,
                "denominator_D": 3 * (2 ** k),
                "n_over_m": f"{8 + 5 * k}/5",
                "exact_rows_in_scan": len(k_rows),
                "minimal_t_row": k_rows[0],
            }
        )
    return summary


def non_family_scan(max_denominator: int, max_m: int, max_n: int) -> list[dict[str, Any]]:
    """Return accidental exact-looking rows not in D=3*2^k family; should be empty.

    This uses exact rational arithmetic on the necessary condition:
    n/m - 8/5 must be a nonnegative integer k and D must be 3*2^k.
    """
    rows = []
    for denominator in range(1, max_denominator + 1):
        for branch_count in range(1, max_m + 1):
            for total_binary_exponent in range(1, max_n + 1):
                exponent_shift = Fraction(total_binary_exponent, branch_count) - Fraction(8, 5)
                if exponent_shift.denominator == 1 and exponent_shift >= 0:
                    k = exponent_shift.numerator
                    exact_family = denominator == 3 * (2 ** k)
                    if not exact_family and denominator <= max_denominator:
                        # This row satisfies the rational-exponent side but not the integer denominator identity.
                        continue
    return rows


def main() -> None:
    minimal_branch_report = load_json(MINIMAL_BRANCH)
    rows = family_rows(MAX_K, MAX_T)
    summary = denominator_summary(rows)
    non_family_rows = non_family_scan(max_denominator=3 * (2 ** MAX_K), max_m=5 * MAX_T, max_n=(8 + 5 * MAX_K) * MAX_T)
    primitive_rows = [row for row in rows if row["primitive_no_binary_rescale"]]
    nonprimitive_rows = [row for row in rows if not row["primitive_no_binary_rescale"]]

    report = {
        "status": "OPEN_STRICT_ALPHA_PRIMITIVE_DENOMINATOR_OBSTRUCTION_NO_DENOMINATOR_THEOREM",
        "result_kind": "SCRATCH_STRICT_ALPHA_PRIMITIVE_DENOMINATOR_OBSTRUCTION_PROBE__NOT_A_THEOREM",
        "source_reports": {
            "minimal_branch_count_obstruction": str(MINIMAL_BRANCH.relative_to(ROOT)),
        },
        "exact_arithmetic_statement": {
            "wide_model": "q_model=(2^n/D^m)^(1/m)=2^(n/m)/D",
            "target": "q_radical=2^(8/5)/3",
            "required_denominator_form": "D=3*2^(n/m-8/5)",
            "integer_D_condition": "n/m-8/5 must be a nonnegative integer k",
            "general_solution": "D=3*2^k, m=5*t, n=(8+5*k)*t for integers k>=0,t>=1",
            "primitive_solution": "k=0 gives D=3, m=5*t, n=8*t",
            "denominator_3_derived_from_target_alone": False,
        },
        "bounded_family_scan": {
            "max_binary_rescale_k": MAX_K,
            "max_refinement_t": MAX_T,
            "exact_family_row_count": len(rows),
            "primitive_row_count": len(primitive_rows),
            "nonprimitive_binary_rescale_row_count": len(nonprimitive_rows),
            "denominator_summary": summary,
            "all_rows_eta_exact": all(abs(row["eta_residual_vs_9_5"]) < 1e-12 for row in rows),
            "non_family_exact_rows_found": non_family_rows,
        },
        "minimal_branch_context_replay": {
            "minimal_branch_result_kind": minimal_branch_report["result_kind"],
            "minimal_branch_model_used_D_3": True,
            "primitive_denominator_needed_to_recover_prior_minimal_n": True,
            "without_primitive_convention": "D=6,n/m=13/5; D=12,n/m=18/5; etc. are exact eta=9/5 rescaled families.",
        },
        "candidate_interpretation": {
            "supported_by_this_probe": bool(
                len(rows) == (MAX_K + 1) * MAX_T
                and len(primitive_rows) == MAX_T
                and len(nonprimitive_rows) == MAX_K * MAX_T
                and all(abs(row["eta_residual_vs_9_5"]) < 1e-12 for row in rows)
                and not non_family_rows
            ),
            "content": "The radical target permits an infinite denominator family D=3*2^k; D=3 is primitive/no-binary-rescale, not forced by the target alone.",
            "why_this_is_more_proof_like": "It isolates the exact integer-denominator family and distinguishes primitive denominator selection from eta arithmetic.",
            "why_this_is_not_enough": "It does not derive the primitive/no-hidden-binary-rescale convention from strict nadsoliton geometry, and it does not solve the balanced-ledger selector obligation.",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is exported.",
            "No legacy physical-role transfer is licensed.",
            "No theorem derives the primitive denominator D=3 from strict nadsoliton geometry.",
            "The exact target admits nonprimitive binary-rescaled denominators D=3*2^k for k>0.",
            "No selector theorem chooses (2,2,2,1,1) from the minimal canonical ledgers.",
            "No theorem derives eta=9/5 without adopting a branch model, a primitive-denominator convention, and a selector convention as extra premises.",
            "No QW-2191 selector discharge and no ToE closure are claimed.",
        ],
        "next_honest_step": "Try to derive a primitive/no-hidden-binary-rescale denominator convention from strict nadsoliton geometry; absent that, denominator-3 remains a primitive gauge choice inside the branch model.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha primitive denominator obstruction probe\n\n"
        "Status: primitive-denominator certificate for eta=9/5; no denominator theorem.\n\n"
        "- Wide model: `q_model=(2^n/D^m)^(1/m)`; exact target requires `D=3*2^(n/m-8/5)`.\n"
        f"- Exact integer family: `D=3*2^k`, `m=5*t`, `n=(8+5*k)*t`; scanned `k<= {MAX_K}`, `t<= {MAX_T}` with `{len(rows)}` exact rows.\n"
        "- Primitive branch: `k=0` recovers `D=3,m=5*t,n=8*t`; nonprimitive branches `k>0` are hidden-binary-rescale exact families.\n"
        "- Honest read: denominator 3 is primitive inside the model, but the primitive/no-rescale convention is still not derived.\n"
        "- No false pass: no denominator theorem, no selector discharge, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
