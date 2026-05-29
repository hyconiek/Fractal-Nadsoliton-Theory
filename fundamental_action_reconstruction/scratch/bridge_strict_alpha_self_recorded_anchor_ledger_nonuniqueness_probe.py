#!/usr/bin/env python3
"""Scratch probe: self-recorded endpoint anchor ledger non-uniqueness audit.

The self-recorded d5 endpoint certificate showed that the balanced ordered ledger
(2,2,2,1,1) can carry an internal endpoint anchor once d5 support and that
ledger are already present.  This probe asks the next honest question: does the
existence of such a self-record select the balanced ledger, or is it a more
common property of positive d5 ledgers?

Finite answer: it is not unique to the balanced ledger.  For a five-node d5 path
with total binary exponent 8, there are 35 positive ordered ledgers.  A simple
endpoint convention can orient/source exactly those with unequal endpoint values;
22 of 35 pass.  Every canonical partition class has endpoint-distinct examples,
including (4,1,1,1,1), (3,2,1,1,1), and (2,2,2,1,1).  Thus self-recording is a
useful decoder/consistency condition, not a ledger selector by itself.
"""
from __future__ import annotations

import json
import math
from collections import Counter
from fractions import Fraction
from itertools import permutations
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SELF_RECORDED = HERE / "bridge_strict_alpha_self_recorded_d5_endpoint_anchor_certificate_report.json"
OUT_JSON = HERE / "bridge_strict_alpha_self_recorded_anchor_ledger_nonuniqueness_report.json"
OUT_MD = HERE / "bridge_strict_alpha_self_recorded_anchor_ledger_nonuniqueness_report.md"

Z12_NODE_COUNT = 12
SUPPORT_SIZE = 5
DISTANCE_SELECTED = 5
TARGET_BINARY_EXPONENT = 8
DENOMINATOR = 3
STRICT_TARGET_ETA = 9.0 / 5.0
NAD12_SUPPORT_SIZE = 12
ALPHA_SCALE = 4.0
BALANCED_LEDGER = (2, 2, 2, 1, 1)
CANONICAL_PARTITIONS = ((4, 1, 1, 1, 1), (3, 2, 1, 1, 1), (2, 2, 2, 1, 1))


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def eta_from_product(product: Fraction, branch_count: int) -> float:
    correction = float(product) ** (1.0 / branch_count)
    return math.log(NAD12_SUPPORT_SIZE * correction) / math.log(ALPHA_SCALE)


def positive_compositions(total: int, parts: int) -> list[tuple[int, ...]]:
    if parts == 1:
        return [(total,)] if total >= 1 else []
    rows: list[tuple[int, ...]] = []
    for head in range(1, total - parts + 2):
        for tail in positive_compositions(total - head, parts - 1):
            rows.append((head,) + tail)
    return rows


def canonical_partition(ledger: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sorted(ledger, reverse=True))


def endpoint_anchor_status(ledger: tuple[int, ...]) -> dict[str, Any]:
    left, right = ledger[0], ledger[-1]
    endpoint_distinct = left != right
    if endpoint_distinct:
        inferred_source_endpoint = "left" if left > right else "right"
        inferred_orientation_sign = 1 if inferred_source_endpoint == "left" else -1
    else:
        inferred_source_endpoint = None
        inferred_orientation_sign = None
    return {
        "endpoint_values": [left, right],
        "endpoint_values_distinct": endpoint_distinct,
        "endpoint_value_gap_abs": abs(left - right),
        "higher_endpoint_source_convention_recovers_orientation": endpoint_distinct,
        "inferred_source_endpoint_under_higher_endpoint_convention": inferred_source_endpoint,
        "inferred_orientation_sign_under_higher_endpoint_convention": inferred_orientation_sign,
        "ordered_values_if_walked_from_higher_endpoint": list(ledger if left > right else tuple(reversed(ledger))) if endpoint_distinct else None,
    }


def ledger_rows() -> list[dict[str, Any]]:
    rows = []
    for ledger in positive_compositions(TARGET_BINARY_EXPONENT, SUPPORT_SIZE):
        status = endpoint_anchor_status(ledger)
        rows.append(
            {
                "ordered_ledger": list(ledger),
                "canonical_partition": list(canonical_partition(ledger)),
                **status,
            }
        )
    return rows


def partition_summary(rows: list[dict[str, Any]]) -> dict[str, Any]:
    by_partition: dict[tuple[int, ...], list[dict[str, Any]]] = {}
    for row in rows:
        by_partition.setdefault(tuple(row["canonical_partition"]), []).append(row)

    summary: dict[str, Any] = {}
    for partition in sorted(by_partition, reverse=True):
        partition_rows = by_partition[partition]
        passing = [row for row in partition_rows if row["endpoint_values_distinct"]]
        failing = [row for row in partition_rows if not row["endpoint_values_distinct"]]
        summary[",".join(map(str, partition))] = {
            "canonical_partition": list(partition),
            "ordered_permutation_count": len(partition_rows),
            "endpoint_distinct_count": len(passing),
            "endpoint_equal_count": len(failing),
            "has_endpoint_anchor_examples": bool(passing),
            "example_endpoint_anchor_ledgers": [row["ordered_ledger"] for row in passing[:5]],
            "example_endpoint_equal_ledgers": [row["ordered_ledger"] for row in failing[:5]],
        }
    return summary


def unique_permutation_count(partition: tuple[int, ...]) -> int:
    return len(set(permutations(partition)))


def canonical_partition_replay() -> dict[str, Any]:
    replay = {}
    for partition in CANONICAL_PARTITIONS:
        perms = sorted(set(permutations(partition)))
        passing = [perm for perm in perms if perm[0] != perm[-1]]
        replay[",".join(map(str, partition))] = {
            "canonical_partition": list(partition),
            "unique_permutation_count": unique_permutation_count(partition),
            "endpoint_distinct_permutation_count": len(passing),
            "sorted_descending_order_endpoint_anchor_status": endpoint_anchor_status(partition),
        }
    return replay


def main() -> None:
    self_recorded = load_json(SELF_RECORDED)
    product = Fraction(2 ** TARGET_BINARY_EXPONENT, DENOMINATOR ** SUPPORT_SIZE)
    rows = ledger_rows()
    passing = [row for row in rows if row["endpoint_values_distinct"]]
    failing = [row for row in rows if not row["endpoint_values_distinct"]]
    partition_counts = Counter(tuple(row["canonical_partition"]) for row in rows)

    report = {
        "status": "OPEN_STRICT_ALPHA_SELF_RECORDED_ANCHOR_LEDGER_NONUNIQUENESS_NO_STRICT_SELECTOR_DISCHARGE",
        "result_kind": "SCRATCH_STRICT_ALPHA_SELF_RECORDED_ANCHOR_LEDGER_NONUNIQUENESS_PROBE__NOT_A_THEOREM",
        "source_reports": {
            "self_recorded_d5_endpoint_anchor_certificate": str(SELF_RECORDED.relative_to(ROOT)),
        },
        "previous_self_recorded_anchor_replay": {
            "result_kind": self_recorded["result_kind"],
            "all_sources_recovered": self_recorded["anchor_reconstruction_scan"]["all_sources_recovered"],
            "all_orientations_recovered": self_recorded["anchor_reconstruction_scan"]["all_orientations_recovered"],
            "d12_equivariance_checked_cases": self_recorded["d12_equivariance_audit"]["checked_cases"],
            "candidate_status": self_recorded["candidate_interpretation"]["status"],
        },
        "target_identity_replay": {
            "q_power_product": f"{product.numerator}/{product.denominator}",
            "eta_residual_vs_9_5": eta_from_product(product, SUPPORT_SIZE) - STRICT_TARGET_ETA,
            "balanced_ledger": list(BALANCED_LEDGER),
        },
        "endpoint_anchor_criterion": {
            "premise": "A five-node d5 path with a positive ordered ledger of total exponent 8 is already present.",
            "criterion": "The endpoint values are unequal; with a higher-endpoint-as-source convention, the path has a self-recorded orientation/source.",
            "failure_mode": "Equal endpoint values leave the two path directions endpoint-symmetric under this convention.",
            "why_not_selector": "The criterion identifies ledgers that can carry an endpoint anchor, but it does not rank or derive which ledger should exist.",
        },
        "ordered_ledger_scan": {
            "total_positive_ordered_ledgers": len(rows),
            "endpoint_distinct_self_recording_count": len(passing),
            "endpoint_equal_ambiguous_count": len(failing),
            "rows": rows,
        },
        "canonical_partition_summary": partition_summary(rows),
        "canonical_partition_replay": canonical_partition_replay(),
        "selector_consequence": {
            "what_is_gained": "Self-recorded endpoint anchoring is a real internal decoding property for endpoint-distinct d5 path ledgers.",
            "what_is_ruled_out": "Self-recording alone does not select the balanced ledger, because all canonical partition classes have endpoint-anchor examples.",
            "remaining_selector_obligation": "A strict derivation still needs a ledger selector such as min-ripple/majorization or another strict self-record formation theorem that singles out (2,2,2,1,1).",
        },
        "candidate_interpretation": {
            "supported_by_this_probe": True,
            "content": "Endpoint self-recording is nonunique across positive d5 ledgers with total exponent 8.",
            "why_this_is_more_proof_like": "The probe exhaustively enumerates all 35 positive ordered ledgers and partitions them by the exact endpoint-distinct criterion.",
            "why_this_is_not_enough": "It shows the self-record idea is a decoder, not a complete ledger selector or strict QW-2191 discharge.",
            "status": "candidate-supported-but-self-record-anchor-not-ledger-selector",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is exported.",
            "No legacy physical-role transfer is licensed.",
            "The nadsoliton is treated as information carrying finite self-record patterns, but no separate sub-nadsoliton informational layer is introduced.",
            "Endpoint self-recording is not unique to the balanced ledger; it is not a ledger selector theorem.",
            "The higher-endpoint-as-source convention is audited as a finite criterion, not derived from strict nadsoliton geometry.",
            "No theorem derives the d5 support, balanced ledger, endpoint-source convention, or strict source-localizing term from strict nadsoliton geometry.",
            "No theorem derives branch count m=5, total exponent n=8, denominator 3, or binary-rescale quotient.",
            "No Aut(Z_12)/N462/QW-2191 selector obstruction is discharged by this endpoint-anchor nonuniqueness audit.",
            "No theorem derives eta=9/5 without the branch model, denominator/quotient convention, support premise, ledger selector, orientation premise, and endpoint-anchor premise.",
            "No QW-2191 discharge and no ToE closure are claimed.",
            "No ToE closure is claimed.",
        ],
        "next_honest_step": "Try to derive why the self-recorded formation ledger should be the balanced endpoint-valued ledger rather than another endpoint-distinct positive d5 ledger, or prove a no-go for such a derivation under current strict assumptions.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha self-recorded anchor ledger nonuniqueness probe\n\n"
        "Status: endpoint self-record nonuniqueness audit; no strict selector discharge.\n\n"
        f"- Positive ordered ledgers: `{len(rows)}`; endpoint-distinct self-recording: `{len(passing)}`; endpoint-equal ambiguous: `{len(failing)}`.\n"
        f"- Canonical partition counts: `{ {','.join(map(str, key)): value for key, value in partition_counts.items()} }`.\n"
        "- Honest read: endpoint self-recording is a decoder/consistency property, not a unique ledger selector.\n"
        f"- Target replay: `q^5={product.numerator}/{product.denominator}`, eta residual `{eta_from_product(product, SUPPORT_SIZE)-STRICT_TARGET_ETA:.3e}`.\n"
        "- No false pass: no strict ledger theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
