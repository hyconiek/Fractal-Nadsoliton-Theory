#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
from fractions import Fraction
from itertools import permutations
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"

OUT = GEN / "p2369_s1319_self_recorded_ledger_closed_form_uniqueness_theorem.json"
MD = GEN / "p2369_s1319_self_recorded_ledger_closed_form_uniqueness_theorem.md"

SOURCE_FILES = {
    "P2368_SELF_RECORDED_ENDPOINT_ANCHOR": GEN
    / "p2368_s1318_self_recorded_endpoint_anchor_selector_candidate_probe.json",
    "P2367_ADMISSIBILITY_NO_GO": GEN
    / "p2367_s1317_selector_phase_origin_admissibility_no_go_probe.json",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

PARTS = 5
TOTAL = 8
BASE_VALUE = TOTAL // PARTS
REMAINDER = TOTAL % PARTS
BALANCED_LEDGER = (2, 2, 2, 1, 1)
STRICT_TARGET_ETA = Fraction(9, 5)
NAD12_SUPPORT_SIZE = 12
ALPHA_SCALE = 4
DENOMINATOR = 3


def rel(path: Path) -> str:
    return path.relative_to(REPO).as_posix()


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": rel(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode(
            "utf-8"
        )
    ).hexdigest()


def positive_compositions(total: int, parts: int) -> list[tuple[int, ...]]:
    if parts == 1:
        return [(total,)] if total >= 1 else []
    rows: list[tuple[int, ...]] = []
    for head in range(1, total - parts + 2):
        for tail in positive_compositions(total - head, parts - 1):
            rows.append((head,) + tail)
    return rows


def ripple(ledger: tuple[int, ...]) -> int:
    return sum(value * value for value in ledger)


def arrow(ledger: tuple[int, ...]) -> int:
    return sum(max(0, right - left) for left, right in zip(ledger, ledger[1:]))


def pairwise_smoothing_delta(a: int, b: int) -> int:
    """Positive drop in sum of squares after moving one unit from a to b when a>=b+2."""
    return a * a + b * b - ((a - 1) * (a - 1) + (b + 1) * (b + 1))


def closed_form_ripple_certificate() -> dict[str, Any]:
    lower_bound = (PARTS - REMAINDER) * BASE_VALUE * BASE_VALUE + REMAINDER * (BASE_VALUE + 1) ** 2
    balanced_multiset = tuple(sorted([BASE_VALUE + 1] * REMAINDER + [BASE_VALUE] * (PARTS - REMAINDER), reverse=True))
    minimizer_permutations = sorted({perm for perm in permutations(balanced_multiset)})
    all_rows = positive_compositions(TOTAL, PARTS)
    brute_min = min(ripple(row) for row in all_rows)
    brute_minimizers = sorted(row for row in all_rows if ripple(row) == brute_min)
    smoothing_deltas = sorted(
        {
            pairwise_smoothing_delta(a, b)
            for a in range(1, TOTAL + 1)
            for b in range(1, TOTAL + 1)
            if a >= b + 2
        }
    )
    return {
        "integer_convexity_statement": "For fixed positive integer sum, sum(e_i^2) strictly decreases whenever a pair differs by at least 2 and one unit is moved from the larger to the smaller.",
        "smoothing_delta_formula": "a^2+b^2-((a-1)^2+(b+1)^2)=2*(a-b-1)>0 for a>=b+2",
        "observed_positive_smoothing_deltas": smoothing_deltas,
        "base_value": BASE_VALUE,
        "remainder": REMAINDER,
        "closed_form_lower_bound": lower_bound,
        "balanced_multiset_descending": list(balanced_multiset),
        "permutation_minimizer_count": len(minimizer_permutations),
        "permutation_minimizers": [list(row) for row in minimizer_permutations],
        "brute_force_composition_count_cross_check": len(all_rows),
        "brute_force_minimum_cross_check": brute_min,
        "brute_force_minimizers_cross_check": [list(row) for row in brute_minimizers],
        "closed_form_matches_bruteforce": lower_bound == brute_min
        and [list(row) for row in minimizer_permutations] == [list(row) for row in brute_minimizers],
    }


def arrow_tiebreak_certificate(minimizers: list[list[int]]) -> dict[str, Any]:
    rows = [
        {
            "ledger": row,
            "arrow_penalty": arrow(tuple(row)),
            "is_nonincreasing": all(left >= right for left, right in zip(row, row[1:])),
        }
        for row in minimizers
    ]
    minimum_arrow = min(row["arrow_penalty"] for row in rows)
    winners = [row["ledger"] for row in rows if row["arrow_penalty"] == minimum_arrow]
    return {
        "arrow_zero_iff_nonincreasing_on_minimizer_permutations": all(
            (row["arrow_penalty"] == 0) == row["is_nonincreasing"] for row in rows
        ),
        "minimum_arrow": minimum_arrow,
        "winner_count": len(winners),
        "winners": winners,
        "unique_winner_is_balanced_ledger": winners == [list(BALANCED_LEDGER)],
        "minimizer_tiebreak_rows": rows,
    }


def endpoint_anchor_closed_form_certificate() -> dict[str, Any]:
    endpoints = [BALANCED_LEDGER[0], BALANCED_LEDGER[-1]]
    return {
        "ordered_ledger": list(BALANCED_LEDGER),
        "endpoint_values": endpoints,
        "endpoint_values_are_distinct": endpoints[0] != endpoints[1],
        "source_endpoint_rule": "Choose the endpoint with value 2 and orient the d5 path toward the endpoint with value 1.",
        "source_endpoint_index": 0,
        "terminus_endpoint_index": PARTS - 1,
        "reversal_effect": "Reversing the ordered ledger swaps the endpoint source and flips orientation.",
        "why_internal_order_is_required": "The unordered value multiset {2,2,2,1,1} has no endpoint labels and cannot by itself choose a source endpoint.",
    }


def eta_identity_certificate() -> dict[str, Any]:
    product = Fraction(2**TOTAL, DENOMINATOR**PARTS)
    lhs_power = Fraction(NAD12_SUPPORT_SIZE**PARTS) * product
    rhs_power = Fraction(ALPHA_SCALE) ** (STRICT_TARGET_ETA * PARTS)
    return {
        "q_power_product": f"{product.numerator}/{product.denominator}",
        "identity_checked": "(12^5)*(2^8/3^5)=4^9, hence log_4(12*(2^8/3^5)^(1/5))=9/5",
        "lhs_power_integer": lhs_power.numerator,
        "rhs_power_integer": rhs_power.numerator,
        "exact_eta_identity_holds": lhs_power == rhs_power,
        "floating_eta_residual_vs_9_5": math.log(float(NAD12_SUPPORT_SIZE * product ** Fraction(1, PARTS)), ALPHA_SCALE)
        - float(STRICT_TARGET_ETA),
    }


def main() -> None:
    GEN.mkdir(exist_ok=True)
    artifacts = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    docs = {name: path.read_text(encoding="utf-8") if path.exists() else "" for name, path in DOC_FILES.items()}

    ripple_cert = closed_form_ripple_certificate()
    arrow_cert = arrow_tiebreak_certificate(ripple_cert["permutation_minimizers"])
    endpoint_cert = endpoint_anchor_closed_form_certificate()
    eta_cert = eta_identity_certificate()

    theorem_export = {
        "theorem_name": "P2369 closed-form self-recorded ledger uniqueness theorem",
        "claim": (
            "For five positive integer slots with total 8, convex smoothing proves the "
            "closed-form ripple lower bound 14 and shows that the only ripple minimizers "
            "are the ten permutations of the multiset {2,2,2,1,1}. The self-recorded "
            "arrow penalty is zero exactly on the nonincreasing permutation, so the "
            "lexicographic (ripple, arrow) action uniquely selects (2,2,2,1,1). Its "
            "distinct endpoint values then define a conditional endpoint-source anchor, "
            "but no strict derivation of the ordered d5 support or arrow action is claimed."
        ),
        "positive_content": [
            "Closed-form convexity replaces raw enumeration as the primary proof of the ripple minimizer set.",
            "A brute-force 35-composition scan is retained only as a cross-check.",
            "Arrow tiebreak uniquely selects the nonincreasing ordered ledger (2,2,2,1,1).",
            "The endpoint values 2 and 1 are distinct, so the ordered d5 path conditionally self-records a source endpoint.",
            "The eta identity associated with q^5=2^8/3^5 is checked exactly as (12^5)*(2^8/3^5)=4^9.",
        ],
        "not_licensed": [
            "strict derivation of ordered d5 support",
            "strict derivation of the arrow action from nadsoliton dynamics",
            "strict derivation of endpoint-source convention",
            "strict selector closure",
            "QW-2191 discharge",
            "beta_tors -> chi11 theorem",
            "legacy physical-role transfer",
            "ToE closure",
        ],
    }

    probe = {
        "probe_id": "P2369_S1319_SELF_RECORDED_LEDGER_CLOSED_FORM_UNIQUENESS_THEOREM",
        "source_artifacts": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "docs": {name: rel(path) for name, path in DOC_FILES.items()},
        "p2368_replay": {
            "packet_id": artifacts["P2368_SELF_RECORDED_ENDPOINT_ANCHOR"].get("packet_id"),
            "ledger_unique_winner": artifacts["P2368_SELF_RECORDED_ENDPOINT_ANCHOR"]
            .get("gatekeeper_checks", {})
            .get("ledger_unique_winner"),
            "all_24_anchors_recovered": artifacts["P2368_SELF_RECORDED_ENDPOINT_ANCHOR"]
            .get("gatekeeper_checks", {})
            .get("all_24_anchors_recovered"),
        },
        "closed_form_ripple_certificate": ripple_cert,
        "arrow_tiebreak_certificate": arrow_cert,
        "endpoint_anchor_closed_form_certificate": endpoint_cert,
        "eta_identity_certificate": eta_cert,
        "theorem_export": theorem_export,
    }
    probe["theorem_fingerprint_sha256"] = sha256_json(theorem_export)

    gatekeeper_checks = {
        "source_artifacts_present": all("_missing" not in artifact for artifact in artifacts.values()),
        "p2368_loaded": probe["p2368_replay"]["packet_id"] == "P2368",
        "closed_form_matches_bruteforce": ripple_cert["closed_form_matches_bruteforce"],
        "ripple_minimizer_count_is_ten": ripple_cert["permutation_minimizer_count"] == 10,
        "arrow_tiebreak_unique": arrow_cert["unique_winner_is_balanced_ledger"],
        "endpoint_anchor_distinct": endpoint_cert["endpoint_values_are_distinct"],
        "eta_identity_exact": eta_cert["exact_eta_identity_holds"],
        "docs_updated_with_p2369_theorem": all("P2369/S1319" in text for text in docs.values()),
        "no_qw2191_discharge_claimed": True,
    }

    payload = {
        "schema_version": "p2369_s1319_v1",
        "packet_id": "P2369",
        "stage_id": "S1319",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-06-01T00:00:00+00:00",
        "status": "OPEN_PROGRESS_CLOSED_FORM_LEDGER_UNIQUENESS_NO_QW2191_DISCHARGE",
        "result_kind": "SELF_RECORDED_LEDGER_CLOSED_FORM_UNIQUENESS_THEOREM",
        "self_recorded_ledger_closed_form_uniqueness_theorem": probe,
        "recommended_next_honest_step": (
            "Move from the now closed finite ledger uniqueness proof to the real missing premise: "
            "derive or refute the ordered d5 support and arrow action from bridge-completed "
            "nadsoliton dynamics, without importing beta_tors or any legacy physical role as a selector."
        ),
        "gatekeeper_checks": gatekeeper_checks,
        "global_status": "OPEN_PROGRESS_WITH_TRACE_NO_SELECTOR_OR_TOE_CLOSURE",
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2369 S1319: self-recorded ledger closed-form uniqueness theorem",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Result",
                "",
                "P2369/S1319 replaces the P2368 ledger enumeration with a closed-form integer proof: convex smoothing gives the ripple lower bound, and the arrow penalty uniquely selects the nonincreasing ordered minimizer.",
                "",
                "## Certificate",
                "",
                f"- Ripple lower bound: `{ripple_cert['closed_form_lower_bound']}`.",
                f"- Ripple minimizer permutations: `{ripple_cert['permutation_minimizer_count']}`.",
                f"- Brute-force cross-check over 35 compositions: `{ripple_cert['closed_form_matches_bruteforce']}`.",
                f"- Arrow tiebreak unique winner: `{arrow_cert['unique_winner_is_balanced_ledger']}` with `{arrow_cert['winners']}`.",
                f"- Endpoint values distinct: `{endpoint_cert['endpoint_values_are_distinct']}`.",
                f"- Exact eta identity holds: `{eta_cert['exact_eta_identity_holds']}`.",
                "",
                "## Hard limits",
                "",
                "- This proves finite ledger uniqueness, not strict derivation of the ordered d5 support or arrow action.",
                "- No `beta_tors -> chi11` theorem, legacy role transfer, QW-2191 discharge, selector closure, or ToE closure is claimed.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
