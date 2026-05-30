#!/usr/bin/env python3
"""Scratch probe: gap-certificate proof of the full-Aut exact-cover UNSAT result.

The previous full-Aut clause-closure audit enumerated all C(12,5)=792 Boolean
supports and found that forbidding folded shells d1,d5,d6 has zero 5-support
solutions.  This probe adds a more proof-like finite certificate: the graph on
Z_12 with forbidden distances {1,5,6} has independence number exactly 3.

Proof idea.  Any independent k-support has cyclic gaps summing to 12.  The d1
ban forces every gap to be at least 2.  Enumerating only the integer gap
necklaces for k=4 and k=5, each necklace contains a consecutive interval whose
folded length is 5 or 6, hence a forbidden d5/d6 pair.  Since explicit 3-support
witnesses exist, alpha(G)=3.  Therefore the five-active full-Aut exact-cover
closure is UNSAT for a structural reason, not merely because a brute-force scan
returned no rows.

No false pass: this is a graph/gap certificate for the finite clause closure,
not a strict derivation of the chi_11 kernel, not a QW-2191 discharge, and not
ToE closure.
"""
from __future__ import annotations

import json
from itertools import combinations, product
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
OUT_JSON = HERE / "bridge_strict_alpha_full_aut_gap_independence_certificate_report.json"
OUT_MD = HERE / "bridge_strict_alpha_full_aut_gap_independence_certificate_report.md"

N = 12
UNITS = [1, 5, 7, 11]
FORBIDDEN_SHELLS = [1, 5, 6]
TARGET_ACTIVE_COUNT = 5
TARGET_Q_POWER = "256/243"
TARGET_ETA = "9/5"


def folded(value: int) -> int:
    residue = value % N
    return min(residue, N - residue)


def canonical_necklace(gaps: tuple[int, ...]) -> tuple[int, ...]:
    rotations = [gaps[index:] + gaps[:index] for index in range(len(gaps))]
    reversals = [tuple(reversed(rotation)) for rotation in rotations]
    return min(rotations + reversals)


def cyclic_gap_tuple(support: tuple[int, ...]) -> tuple[int, ...]:
    cyclic = sorted(support)
    gaps = [cyclic[index + 1] - cyclic[index] for index in range(len(cyclic) - 1)]
    gaps.append(N + cyclic[0] - cyclic[-1])
    return tuple(gaps)


def support_satisfies(support: tuple[int, ...], forbidden_shells: set[int]) -> bool:
    return all(folded(right - left) not in forbidden_shells for left, right in combinations(support, 2))


def independent_supports(size: int) -> list[tuple[int, ...]]:
    forbidden = set(FORBIDDEN_SHELLS)
    return [support for support in combinations(range(N), size) if support_satisfies(support, forbidden)]


def gap_necklaces_with_min_gap(size: int, min_gap: int = 2) -> list[tuple[int, ...]]:
    necklaces = {
        canonical_necklace(tuple(gaps))
        for gaps in product(range(min_gap, N + 1), repeat=size)
        if sum(gaps) == N
    }
    return sorted(necklaces)


def first_forbidden_interval(gaps: tuple[int, ...]) -> dict[str, Any] | None:
    size = len(gaps)
    for start in range(size):
        interval_sum = 0
        for length in range(1, size):
            interval_sum += gaps[(start + length - 1) % size]
            interval_folded = folded(interval_sum)
            if interval_folded in {5, 6}:
                return {
                    "start_gap_index": start,
                    "consecutive_gap_count": length,
                    "oriented_interval_sum": interval_sum,
                    "folded_distance": interval_folded,
                    "forbidden_shell_hit": f"d{interval_folded}",
                }
    return None


def gap_elimination_rows(size: int) -> list[dict[str, Any]]:
    rows = []
    for necklace in gap_necklaces_with_min_gap(size):
        obstruction = first_forbidden_interval(necklace)
        rows.append(
            {
                "active_count": size,
                "gap_necklace": list(necklace),
                "gap_sum": sum(necklace),
                "all_gaps_at_least_2_from_d1_ban": all(gap >= 2 for gap in necklace),
                "forbidden_interval_certificate": obstruction,
                "survives_d5_d6_bans": obstruction is None,
            }
        )
    return rows


def independence_profile() -> list[dict[str, Any]]:
    rows = []
    for size in range(0, TARGET_ACTIVE_COUNT + 1):
        supports = independent_supports(size)
        necklaces = sorted({canonical_necklace(cyclic_gap_tuple(support)) for support in supports if support})
        rows.append(
            {
                "active_count": size,
                "independent_support_count": len(supports),
                "representative_supports": [list(support) for support in supports[:8]],
                "gap_necklaces": [list(necklace) for necklace in necklaces],
                "exists": bool(supports),
            }
        )
    return rows


def exact_proof_certificate(profile: list[dict[str, Any]], gap_rows_4: list[dict[str, Any]], gap_rows_5: list[dict[str, Any]]) -> dict[str, str]:
    size3 = next(row for row in profile if row["active_count"] == 3)
    size4 = next(row for row in profile if row["active_count"] == 4)
    size5 = next(row for row in profile if row["active_count"] == 5)
    return {
        "graph_definition": "G has vertex set Z_12 and forbidden edges at folded distances d1,d5,d6, the full-Aut closure of the earlier d1+d6 exact-cover clauses.",
        "gap_reduction": "For any independent support, the d1 ban forces every cyclic gap to be at least 2; all possible k-support gap certificates are therefore integer necklaces summing to 12.",
        "four_support_elimination": f"All {len(gap_rows_4)} min-gap necklaces for k=4 contain a consecutive interval of folded length 5 or 6, so no 4-support survives.",
        "five_support_elimination": f"All {len(gap_rows_5)} min-gap necklaces for k=5 contain a consecutive interval of folded length 5 or 6, so no 5-support survives.",
        "tightness": f"There are {size3['independent_support_count']} independent 3-supports, while k=4 count is {size4['independent_support_count']} and k=5 count is {size5['independent_support_count']}; hence alpha(G)=3.",
        "selector_consequence": "The target exact-cover cardinality 5 is impossible after full-Aut clause closure; the successful d1+d6 selector remains chi_11-conditional.",
    }


def build_payload() -> dict[str, Any]:
    profile = independence_profile()
    gap_rows_4 = gap_elimination_rows(4)
    gap_rows_5 = gap_elimination_rows(5)
    maximum_independent_size = max(row["active_count"] for row in profile if row["exists"])
    return {
        "result_kind": "SCRATCH_STRICT_ALPHA_FULL_AUT_GAP_INDEPENDENCE_CERTIFICATE_PROBE__NO_GO_NOT_A_THEOREM",
        "status": "full-aut-d1-d5-d6-closure-has-independence-number-3-so-five-active-exact-cover-is-unsat",
        "finite_model": {
            "ring": "Z_12",
            "automorphism_units": UNITS,
            "forbidden_shells_full_Aut_closed": FORBIDDEN_SHELLS,
            "target_active_count": TARGET_ACTIVE_COUNT,
            "target_q_power": TARGET_Q_POWER,
            "target_eta": TARGET_ETA,
        },
        "independence_profile": profile,
        "maximum_independent_size": maximum_independent_size,
        "target_five_active_is_unsat": maximum_independent_size < TARGET_ACTIVE_COUNT,
        "gap_elimination_certificate": {
            "k4_rows": gap_rows_4,
            "k5_rows": gap_rows_5,
            "all_k4_necklaces_eliminated_by_d5_or_d6": all(not row["survives_d5_d6_bans"] for row in gap_rows_4),
            "all_k5_necklaces_eliminated_by_d5_or_d6": all(not row["survives_d5_d6_bans"] for row in gap_rows_5),
        },
        "exact_proof_certificate": exact_proof_certificate(profile, gap_rows_4, gap_rows_5),
        "interpretation": {
            "what_was_added": "A compact gap-necklace certificate proves the previous full-Aut UNSAT result and strengthens it to alpha(G)=3.",
            "why_more_proof_like": "Only 8 k=4 and 3 k=5 min-gap necklaces need elimination after the d1 gap reduction, instead of relying only on all 792 Boolean assignments.",
            "honest_limit": "This proves a finite full-Aut closure obstruction, not the origin of chi_11, not QW-2191 discharge, and not ToE closure.",
        },
        "ontology_guardrail": {
            "allowed_reading": "The nadsoliton itself is the primordial information in a solitonic state.",
            "forbidden_reading": "No separate informational layer underneath the nadsoliton is introduced by this graph certificate.",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is used or claimed.",
            "No legacy role transfer to K_strict_gate, alpha_geo, beta_tors, or D_f is made.",
            "No theorem derives the chi_11-kernel, shell-label d1 vs d5, or unit-axis bit from strict nadsoliton geometry.",
            "The result is a finite full-Aut graph/gap no-go certificate, not a selector-origin theorem.",
            "No QW-2191 discharge.",
            "No ToE closure.",
        ],
    }


def write_markdown(payload: dict[str, Any]) -> str:
    proof = payload["exact_proof_certificate"]
    profile = payload["independence_profile"]
    lines = [
        "# Strict alpha full-Aut gap independence certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite model",
        "",
        f"- Ring: `{payload['finite_model']['ring']}`",
        f"- Full-Aut closed forbidden shells: `{payload['finite_model']['forbidden_shells_full_Aut_closed']}`",
        f"- Target active count: `{payload['finite_model']['target_active_count']}`",
        f"- Maximum independent size: `{payload['maximum_independent_size']}`",
        f"- Target five-active UNSAT: `{payload['target_five_active_is_unsat']}`",
        "",
        "## Independence profile",
        "",
    ]
    for row in profile:
        lines.append(
            f"- k=`{row['active_count']}`: count=`{row['independent_support_count']}`, "
            f"gap necklaces=`{row['gap_necklaces']}`"
        )
    lines.extend([
        "",
        "## Gap elimination certificate",
        "",
        f"- k=4 necklaces eliminated: `{payload['gap_elimination_certificate']['all_k4_necklaces_eliminated_by_d5_or_d6']}`",
        f"- k=5 necklaces eliminated: `{payload['gap_elimination_certificate']['all_k5_necklaces_eliminated_by_d5_or_d6']}`",
        f"- k=4 necklace count after d1 gap reduction: `{len(payload['gap_elimination_certificate']['k4_rows'])}`",
        f"- k=5 necklace count after d1 gap reduction: `{len(payload['gap_elimination_certificate']['k5_rows'])}`",
        "",
        "## Proof certificate",
        "",
    ])
    for key, value in proof.items():
        lines.append(f"- `{key}`: {value}")
    lines.extend(["", "## Hard limits", ""])
    lines.extend(f"- {item}" for item in payload["hard_limits"])
    lines.append("")
    return "\n".join(lines)


def main() -> None:
    payload = build_payload()
    OUT_JSON.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    OUT_MD.write_text(write_markdown(payload), encoding="utf-8")
    print(json.dumps(payload, indent=2, sort_keys=True))
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
