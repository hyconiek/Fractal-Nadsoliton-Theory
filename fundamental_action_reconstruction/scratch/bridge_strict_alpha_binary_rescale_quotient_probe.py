#!/usr/bin/env python3
"""Scratch probe: binary-rescale quotient for strict alpha -> eta.

The primitive-denominator obstruction found an exact integer-denominator family

    D = 3*2^k,  m = 5*t,  n = (8+5*k)*t.

This packet asks whether the nonprimitive rows are genuinely new or just a gauge
copy under a binary-rescale equivalence.  In the model

    q_model^m = 2^n / D^m,

if D is even then

    2^n / D^m = 2^(n-m) / (D/2)^m.

So the reduction step

    (D, m, n) -> (D/2, m, n-m)

preserves the model exactly whenever D is even and n>=m.  Every family row with
D=3*2^k reduces in exactly k steps to the primitive representative

    (D, m, n) = (3, 5*t, 8*t).

This is a quotient/canonicalization certificate, not a physical theorem.  It
says D=3 is the unique odd representative under binary-rescale equivalence, but
it does not derive that strict nadsoliton geometry must quotient by this gauge or
instantiate the branch model in the first place.
"""
from __future__ import annotations

import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PRIMITIVE_DENOMINATOR = HERE / "bridge_strict_alpha_primitive_denominator_obstruction_report.json"
OUT_JSON = HERE / "bridge_strict_alpha_binary_rescale_quotient_report.json"
OUT_MD = HERE / "bridge_strict_alpha_binary_rescale_quotient_report.md"

MAX_K = 6
MAX_T = 8
STRICT_TARGET_ETA = Fraction(9, 5)
NAD12_SUPPORT_SIZE = 12
ALPHA_SCALE = 4.0


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def eta_from_model(denominator: int, branch_count: int, total_binary_exponent: int) -> float:
    q_model = (2 ** total_binary_exponent / denominator ** branch_count) ** (1.0 / branch_count)
    return math.log(NAD12_SUPPORT_SIZE * q_model) / math.log(ALPHA_SCALE)


def invariant_fraction(denominator: int, branch_count: int, total_binary_exponent: int) -> Fraction:
    return Fraction(2 ** total_binary_exponent, denominator ** branch_count)


def reduce_once(denominator: int, branch_count: int, total_binary_exponent: int) -> tuple[int, int, int]:
    if denominator % 2 != 0:
        raise ValueError("binary-rescale reduction requires even denominator")
    if total_binary_exponent < branch_count:
        raise ValueError("binary-rescale reduction requires n>=m")
    return denominator // 2, branch_count, total_binary_exponent - branch_count


def canonicalize(denominator: int, branch_count: int, total_binary_exponent: int) -> dict[str, Any]:
    states = []
    current = (denominator, branch_count, total_binary_exponent)
    initial_invariant = invariant_fraction(*current)
    while current[0] % 2 == 0:
        before = current
        after = reduce_once(*before)
        states.append(
            {
                "before": {"D": before[0], "m": before[1], "n": before[2]},
                "after": {"D": after[0], "m": after[1], "n": after[2]},
                "invariant_before": f"{invariant_fraction(*before).numerator}/{invariant_fraction(*before).denominator}",
                "invariant_after": f"{invariant_fraction(*after).numerator}/{invariant_fraction(*after).denominator}",
                "invariant_preserved": invariant_fraction(*before) == invariant_fraction(*after),
            }
        )
        current = after
    final_invariant = invariant_fraction(*current)
    return {
        "initial": {"D": denominator, "m": branch_count, "n": total_binary_exponent},
        "final": {"D": current[0], "m": current[1], "n": current[2]},
        "steps": states,
        "step_count": len(states),
        "initial_invariant": f"{initial_invariant.numerator}/{initial_invariant.denominator}",
        "final_invariant": f"{final_invariant.numerator}/{final_invariant.denominator}",
        "invariant_preserved_globally": initial_invariant == final_invariant,
        "final_denominator_odd": current[0] % 2 == 1,
        "eta_initial": eta_from_model(denominator, branch_count, total_binary_exponent),
        "eta_final": eta_from_model(*current),
        "eta_residual_initial_vs_9_5": eta_from_model(denominator, branch_count, total_binary_exponent) - float(STRICT_TARGET_ETA),
        "eta_residual_final_vs_9_5": eta_from_model(*current) - float(STRICT_TARGET_ETA),
    }


def family_row(k: int, t: int) -> dict[str, Any]:
    denominator = 3 * (2 ** k)
    branch_count = 5 * t
    total_binary_exponent = (8 + 5 * k) * t
    certificate = canonicalize(denominator, branch_count, total_binary_exponent)
    return {
        "binary_rescale_k": k,
        "refinement_t": t,
        "expected_primitive_final": {"D": 3, "m": 5 * t, "n": 8 * t},
        "canonicalization": certificate,
        "reduction_steps_equal_k": certificate["step_count"] == k,
        "final_matches_expected_primitive": certificate["final"] == {"D": 3, "m": 5 * t, "n": 8 * t},
    }


def main() -> None:
    primitive_report = load_json(PRIMITIVE_DENOMINATOR)
    rows = [family_row(k, t) for k in range(MAX_K + 1) for t in range(1, MAX_T + 1)]
    nonprimitive_rows = [row for row in rows if row["binary_rescale_k"] > 0]
    primitive_rows = [row for row in rows if row["binary_rescale_k"] == 0]
    sample_paths = [family_row(k, 1) for k in range(MAX_K + 1)]

    report = {
        "status": "OPEN_STRICT_ALPHA_BINARY_RESCALE_QUOTIENT_NO_PHYSICAL_GAUGE_THEOREM",
        "result_kind": "SCRATCH_STRICT_ALPHA_BINARY_RESCALE_QUOTIENT_PROBE__NOT_A_THEOREM",
        "source_reports": {
            "primitive_denominator_obstruction": str(PRIMITIVE_DENOMINATOR.relative_to(ROOT)),
        },
        "equivalence_rule": {
            "model_invariant": "q_model^m = 2^n/D^m",
            "reduction_step": "if D even and n>=m, (D,m,n)->(D/2,m,n-m)",
            "invariant_identity": "2^n/D^m = 2^(n-m)/(D/2)^m",
            "canonical_representative_condition": "D odd",
            "physical_gauge_quotient_derived_here": False,
        },
        "bounded_family_quotient_scan": {
            "max_binary_rescale_k": MAX_K,
            "max_refinement_t": MAX_T,
            "row_count": len(rows),
            "primitive_row_count": len(primitive_rows),
            "nonprimitive_row_count": len(nonprimitive_rows),
            "all_invariants_preserved": all(row["canonicalization"]["invariant_preserved_globally"] for row in rows),
            "all_final_denominators_odd": all(row["canonicalization"]["final_denominator_odd"] for row in rows),
            "all_reduction_steps_equal_k": all(row["reduction_steps_equal_k"] for row in rows),
            "all_finals_match_expected_primitive": all(row["final_matches_expected_primitive"] for row in rows),
            "all_eta_residuals_exact": all(
                abs(row["canonicalization"]["eta_residual_initial_vs_9_5"]) < 1e-12
                and abs(row["canonicalization"]["eta_residual_final_vs_9_5"]) < 1e-12
                for row in rows
            ),
            "sample_t1_paths_by_k": sample_paths,
        },
        "primitive_denominator_context_replay": {
            "primitive_denominator_result_kind": primitive_report["result_kind"],
            "primitive_report_nonprimitive_rows": primitive_report["bounded_family_scan"]["nonprimitive_binary_rescale_row_count"],
            "quotient_collapses_nonprimitive_rows_to_odd_representatives": True,
            "remaining_obligation": "derive that strict nadsoliton branch descriptions should quotient binary rescale rather than treat D=3*2^k rows as distinct sectors",
        },
        "candidate_interpretation": {
            "supported_by_this_probe": bool(
                len(rows) == (MAX_K + 1) * MAX_T
                and all(row["canonicalization"]["invariant_preserved_globally"] for row in rows)
                and all(row["final_matches_expected_primitive"] for row in rows)
                and all(row["reduction_steps_equal_k"] for row in rows)
            ),
            "content": "All D=3*2^k exact families are binary-rescale equivalent to a unique odd representative D=3 with n=8*t.",
            "why_this_is_more_proof_like": "It provides an explicit invariant-preserving reduction algorithm and bounded certificate paths collapsing nonprimitive denominators to the primitive representative.",
            "why_this_is_not_enough": "It does not derive the physical/gauge principle that binary-rescale equivalent branch descriptions must be identified in strict nadsoliton geometry.",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is exported.",
            "No legacy physical-role transfer is licensed.",
            "No theorem derives binary-rescale quotienting as a strict nadsoliton gauge principle.",
            "No theorem derives the denominator branch model itself from strict nadsoliton geometry.",
            "No selector theorem chooses (2,2,2,1,1) from the minimal canonical ledgers.",
            "No theorem derives eta=9/5 without adopting the branch model, binary-rescale quotient, and selector convention as extra premises.",
            "No QW-2191 selector discharge and no ToE closure are claimed.",
        ],
        "next_honest_step": "Try to derive binary-rescale quotienting as an internal strict gauge/canonicalization principle; if successful, D=3 becomes the canonical odd representative, but the balanced-ledger selector remains open.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha binary-rescale quotient probe\n\n"
        "Status: quotient/canonicalization certificate for eta=9/5; no physical gauge theorem.\n\n"
        "- Reduction rule: `(D,m,n)->(D/2,m,n-m)` preserves `2^n/D^m` when `D` is even and `n>=m`.\n"
        f"- Bounded scan: `{len(rows)}` rows through `k<={MAX_K}`, `t<={MAX_T}`; all reduce to odd representatives with steps equal to `k`.\n"
        "- Canonical result: every `D=3*2^k,m=5*t,n=(8+5*k)*t` row reduces to `D=3,m=5*t,n=8*t`.\n"
        "- Honest read: D=3 is canonical under binary-rescale quotienting, but the quotient/gauge principle itself is not derived.\n"
        "- No false pass: no binary-rescale gauge theorem, no selector discharge, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
