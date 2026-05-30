#!/usr/bin/env python3
"""Scratch probe: bounded linear histogram functionals and the chi_11 source gap.

The previous candidate-source transform audit checked a small hand-picked list of
branch records and found no allowed strict chi_11 source without importing the
d1-vs-d5 shell label.  This probe makes that source audit more proof-like by
exhaustively enumerating a finite basis class: every linear functional

    L_w(h) = sum_{i=1..6} w_i h_i,  w_i in {-1,0,1},

on the folded distance histogram h=(d1,...,d6) of the four unit branches
A1,A5,A7,A11.

Result: among all 3^6=729 bounded shell-linear functionals, 54 transform with
the required chi_11 sign on the branch pairs, but every such functional has
w1 != w5 and therefore imports the d1-vs-d5 shell label.  Among the 3^5=243
full-Aut-invariant functionals (w1=w5), zero distinguish A1 from A5 and zero
export chi_11.  This is a finite no-go for this whole bounded linear class.

No false pass: this is not an exhaustive theorem over all possible strict
nadsoliton sources; it is only an exhaustive certificate for the declared
bounded shell-linear histogram class.  It does not discharge QW-2191 and does
not close ToE.
"""
from __future__ import annotations

import json
from itertools import combinations, product
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
OUT_JSON = HERE / "bridge_strict_alpha_histogram_linear_functional_chi11_source_no_go_report.json"
OUT_MD = HERE / "bridge_strict_alpha_histogram_linear_functional_chi11_source_no_go_report.md"

N = 12
ACTIVE_COUNT = 5
UNITS = [1, 5, 7, 11]
BRANCH_MODES = [1, 5, 7, 11]
COEFFICIENT_DOMAIN = [-1, 0, 1]
TARGET_Q_POWER = "256/243"
TARGET_ETA = "9/5"
A1_HISTOGRAM = (4, 3, 2, 1, 0, 0)
A5_HISTOGRAM = (0, 3, 2, 1, 4, 0)
CONTIGUOUS_PAIR = "contiguous_pair_A1_A11"
D5_PAIR = "d5_pair_A5_A7"


def folded(value: int) -> int:
    residue = value % N
    return min(residue, N - residue)


def unit_support(mode: int) -> tuple[int, ...]:
    return tuple(sorted((mode * index) % N for index in range(ACTIVE_COUNT)))


def distance_histogram(support: tuple[int, ...]) -> tuple[int, int, int, int, int, int]:
    counts = [0] * (N // 2)
    for left, right in combinations(support, 2):
        counts[folded(right - left) - 1] += 1
    return tuple(counts)  # type: ignore[return-value]


def branch_name(mode: int) -> str:
    return f"A{mode}_k{mode}"


def branch_pair(name: str) -> str:
    if name in {"A1_k1", "A11_k11"}:
        return CONTIGUOUS_PAIR
    if name in {"A5_k5", "A7_k7"}:
        return D5_PAIR
    raise ValueError(name)


def required_chi11_value(name: str) -> int:
    return 1 if branch_pair(name) == D5_PAIR else -1


def dot(weights: tuple[int, ...], histogram: tuple[int, ...]) -> int:
    return sum(weight * count for weight, count in zip(weights, histogram))


def is_full_aut_invariant_weight(weights: tuple[int, ...]) -> bool:
    # Full Aut shell orbits identify d1 and d5; linear invariance forces w1=w5.
    return weights[0] == weights[4]


def is_chi11_covariant_weight(weights: tuple[int, ...]) -> bool:
    score_a1 = dot(weights, A1_HISTOGRAM)
    score_a5 = dot(weights, A5_HISTOGRAM)
    return score_a5 == -score_a1 and score_a5 != 0


def support_size(weights: tuple[int, ...]) -> int:
    return sum(1 for weight in weights if weight != 0)


def weight_row(weights: tuple[int, ...]) -> dict[str, Any]:
    score_a1 = dot(weights, A1_HISTOGRAM)
    score_a5 = dot(weights, A5_HISTOGRAM)
    full_aut_invariant = is_full_aut_invariant_weight(weights)
    chi11_covariant = is_chi11_covariant_weight(weights)
    return {
        "weights_d1_to_d6": list(weights),
        "score_A1_A11": score_a1,
        "score_A5_A7": score_a5,
        "score_sum_A1_plus_A5": score_a1 + score_a5,
        "score_difference_A5_minus_A1": score_a5 - score_a1,
        "support_size": support_size(weights),
        "full_Aut_invariant_weight": full_aut_invariant,
        "pair_distinguishing": score_a1 != score_a5,
        "numeric_chi11_covariant": chi11_covariant,
        "imports_d1_vs_d5_shell_label": weights[0] != weights[4],
        "allowed_strict_chi11_source": chi11_covariant and full_aut_invariant,
    }


def all_weight_rows() -> list[dict[str, Any]]:
    return [weight_row(tuple(weights)) for weights in product(COEFFICIENT_DOMAIN, repeat=6)]


def branch_rows() -> list[dict[str, Any]]:
    return [
        {
            "name": branch_name(mode),
            "mode": mode,
            "support": list(unit_support(mode)),
            "branch_pair": branch_pair(branch_name(mode)),
            "required_chi11_value": required_chi11_value(branch_name(mode)),
            "distance_histogram_d1_to_d6": list(distance_histogram(unit_support(mode))),
        }
        for mode in BRANCH_MODES
    ]


def symbolic_certificate() -> dict[str, str]:
    return {
        "score_A1": "For A1/A11, h=[4,3,2,1,0,0], so L(A1)=4*w1+3*w2+2*w3+w4.",
        "score_A5": "For A5/A7, h=[0,3,2,1,4,0], so L(A5)=3*w2+2*w3+w4+4*w5.",
        "pair_difference": "L(A5)-L(A1)=4*(w5-w1); therefore any full-Aut invariant linear functional with w1=w5 cannot distinguish A1 from A5.",
        "chi11_condition": "To transform as chi_11 on the two branch pairs, L(A5)=-L(A1) and L(A5) must be nonzero.",
        "full_aut_no_go": "Combining w1=w5 with the pair-difference formula gives L(A5)=L(A1), incompatible with nonzero chi_11 covariance.",
    }


def build_payload() -> dict[str, Any]:
    rows = all_weight_rows()
    chi_rows = [row for row in rows if row["numeric_chi11_covariant"]]
    full_rows = [row for row in rows if row["full_Aut_invariant_weight"]]
    full_pair_rows = [row for row in full_rows if row["pair_distinguishing"]]
    allowed_rows = [row for row in rows if row["allowed_strict_chi11_source"]]
    minimal_chi_rows = sorted(chi_rows, key=lambda row: (row["support_size"], row["weights_d1_to_d6"]))[:12]
    return {
        "result_kind": "SCRATCH_STRICT_ALPHA_HISTOGRAM_LINEAR_FUNCTIONAL_CHI11_SOURCE_NO_GO_PROBE__NO_GO_NOT_A_THEOREM",
        "status": "bounded-shell-linear-histogram-functionals-export-chi11-only-by-importing-d1-vs-d5-label",
        "finite_model": {
            "ring": "Z_12",
            "active_count": ACTIVE_COUNT,
            "automorphism_units": UNITS,
            "branch_modes": BRANCH_MODES,
            "coefficient_domain": COEFFICIENT_DOMAIN,
            "weight_count": len(rows),
            "target_q_power": TARGET_Q_POWER,
            "target_eta": TARGET_ETA,
        },
        "branch_rows": branch_rows(),
        "symbolic_certificate": symbolic_certificate(),
        "linear_functional_summary": {
            "total_weight_count": len(rows),
            "full_Aut_invariant_weight_count": len(full_rows),
            "pair_distinguishing_weight_count": sum(row["pair_distinguishing"] for row in rows),
            "full_Aut_pair_distinguishing_weight_count": len(full_pair_rows),
            "numeric_chi11_covariant_weight_count": len(chi_rows),
            "numeric_chi11_covariant_full_Aut_invariant_weight_count": sum(row["full_Aut_invariant_weight"] for row in chi_rows),
            "allowed_strict_chi11_source_weight_count": len(allowed_rows),
            "chi11_covariant_imports_shell_label_count": sum(row["imports_d1_vs_d5_shell_label"] for row in chi_rows),
            "chi11_covariant_support_size_histogram": {
                str(size): sum(1 for row in chi_rows if row["support_size"] == size)
                for size in range(7)
                if any(row["support_size"] == size for row in chi_rows)
            },
        },
        "minimal_chi11_covariant_examples": minimal_chi_rows,
        "allowed_strict_chi11_source_rows": allowed_rows,
        "exact_proof_certificate": {
            "finite_domain": "All 3^6=729 weights w_i in {-1,0,1} for shell-linear histogram functionals are enumerated exactly.",
            "symbolic_no_go": symbolic_certificate()["pair_difference"] + " " + symbolic_certificate()["full_aut_no_go"],
            "enumerated_no_go": f"Full-Aut invariant pair-distinguishing weights: {len(full_pair_rows)}; allowed strict chi_11 source weights: {len(allowed_rows)}.",
            "import_boundary": f"There are {len(chi_rows)} chi_11-covariant weights, and {sum(row['imports_d1_vs_d5_shell_label'] for row in chi_rows)} of them have w1!=w5, i.e. they import the d1-vs-d5 shell label.",
            "scope_limit": "This is exhaustive only for bounded shell-linear histogram functionals, not for all possible strict nadsoliton source mechanisms.",
        },
        "interpretation": {
            "honest_positive": "The finite linear class contains many chi_11-covariant witnesses, with minimal examples such as d5-d1.",
            "honest_negative": "Every chi_11-covariant witness in this class breaks full Aut by w1!=w5, so none supplies chi_11 as an unlabelled strict source.",
            "relation_to_previous_probe": "This upgrades the hand-picked candidate-source audit to an exhaustive bounded linear histogram audit.",
        },
        "ontology_guardrail": {
            "allowed_reading": "The nadsoliton itself remains the only admissible place where a future non-linear or non-histogram strict chi_11 source could be derived.",
            "forbidden_reading": "No separate informational layer underneath the nadsoliton is introduced to supply chi_11.",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is used or claimed.",
            "No legacy role transfer to K_strict_gate, alpha_geo, beta_tors, or D_f is made.",
            "No theorem derives the chi_11-kernel, shell-label d1 vs d5, unit-axis bit, exact-cover clauses, or cardinality 5 from strict nadsoliton geometry.",
            "The result is a finite bounded shell-linear histogram no-go, not an exhaustive strict-source theorem.",
            "No QW-2191 discharge.",
            "No ToE closure.",
        ],
    }


def write_markdown(payload: dict[str, Any]) -> str:
    summary = payload["linear_functional_summary"]
    proof = payload["exact_proof_certificate"]
    lines = [
        "# Strict alpha bounded shell-linear chi_11 source no-go",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Linear functional summary",
        "",
        f"- Total weights: `{summary['total_weight_count']}`",
        f"- Full-Aut invariant weights: `{summary['full_Aut_invariant_weight_count']}`",
        f"- Pair-distinguishing weights: `{summary['pair_distinguishing_weight_count']}`",
        f"- Full-Aut pair-distinguishing weights: `{summary['full_Aut_pair_distinguishing_weight_count']}`",
        f"- chi_11-covariant weights: `{summary['numeric_chi11_covariant_weight_count']}`",
        f"- Full-Aut invariant chi_11-covariant weights: `{summary['numeric_chi11_covariant_full_Aut_invariant_weight_count']}`",
        f"- Allowed strict chi_11 source weights: `{summary['allowed_strict_chi11_source_weight_count']}`",
        f"- chi_11-covariant shell-label imports: `{summary['chi11_covariant_imports_shell_label_count']}`",
        f"- chi_11 support-size histogram: `{summary['chi11_covariant_support_size_histogram']}`",
        "",
        "## Symbolic certificate",
        "",
    ]
    for key, value in payload["symbolic_certificate"].items():
        lines.append(f"- `{key}`: {value}")
    lines.extend(["", "## Proof certificate", ""])
    for key, value in proof.items():
        lines.append(f"- `{key}`: {value}")
    lines.extend(["", "## Minimal chi_11-covariant examples", ""])
    for row in payload["minimal_chi11_covariant_examples"][:6]:
        lines.append(
            f"- weights=`{row['weights_d1_to_d6']}`, score_A1=`{row['score_A1_A11']}`, "
            f"score_A5=`{row['score_A5_A7']}`, imports_label=`{row['imports_d1_vs_d5_shell_label']}`"
        )
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
