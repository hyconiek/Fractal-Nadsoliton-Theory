#!/usr/bin/env python3
"""Scratch probe: matter-like self-replication can carry, not originate, the unit bit.

This probe checks the tempting interpretation raised after the prospective
no-go probes: perhaps the missing one-bit unit orientation is connected with
matter because matter carries self-replicating information.  We audit the finite
binary version of that idea with exact combinatorics and a Shannon-style record
model.

Model:

    B=0  means the A1/contiguous side of the unit-mirror orbit,
    B=1  means the A5/d5 side of the unit-mirror orbit.

A matter-like record is an L-bit copy word R in {0,1}^L.  Unit mirror flips B
and every bit of R.  A replication channel is admissible only when it is mirror
equivariant: P(R=r|B=b)=P(R=not(r)|B=1-b).

No false pass: such a self-replicating record can preserve or amplify an already
chosen orientation bit, but with a symmetric prior it cannot create a singleton
d5 selector.  Any decoder that calls the 1-side "d5" imports exactly the missing
unit-axis label.
"""
from __future__ import annotations

import json
from collections import Counter
from fractions import Fraction
from itertools import combinations, product
from math import comb
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
OUT_JSON = HERE / "bridge_strict_alpha_hebbian_matter_shannon_self_replication_selector_audit_report.json"
OUT_MD = HERE / "bridge_strict_alpha_hebbian_matter_shannon_self_replication_selector_audit_report.md"

N = 12
ACTIVE_COUNT = 5
AUT_UNITS = [1, 5, 7, 11]
TARGET_Q_POWER = "256/243"
TARGET_ETA = "9/5"
MAX_RECORD_LENGTH = 8
ENUMERATED_DECODER_MAX_LENGTH = 4
ERROR_GRID_DENOMINATOR = 8
A1 = "A1_k1_contiguous"
A5 = "A5_k5_d5"
D5_HISTOGRAM = [0, 3, 2, 1, 4, 0]
CONTIGUOUS_HISTOGRAM = [4, 3, 2, 1, 0, 0]
Record = tuple[int, ...]


def fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def all_supports() -> list[tuple[int, ...]]:
    return [tuple(support) for support in combinations(range(N), ACTIVE_COUNT)]


def folded(value: int) -> int:
    residue = value % N
    return min(residue, N - residue)


def distance_histogram(support: tuple[int, ...]) -> tuple[int, int, int, int, int, int]:
    counts = [0] * (N // 2)
    for left, right in combinations(support, 2):
        counts[folded(right - left) - 1] += 1
    return tuple(counts)  # type: ignore[return-value]


def flip_record(record: Record) -> Record:
    return tuple(1 - bit for bit in record)


def all_records(length: int) -> list[Record]:
    return list(product([0, 1], repeat=length))


def hamming_weight(record: Record) -> int:
    return sum(record)


def hamming_distance_to_seed(record: Record, seed: int) -> int:
    return sum(bit != seed for bit in record)


def channel_probability(record: Record, seed: int, error: Fraction) -> Fraction:
    errors = hamming_distance_to_seed(record, seed)
    return (error**errors) * ((1 - error) ** (len(record) - errors))


def orbit_representative(record: Record) -> Record:
    return min(record, flip_record(record))


def record_orbits(length: int) -> list[tuple[Record, Record]]:
    seen: set[Record] = set()
    orbits = []
    for record in all_records(length):
        if record in seen:
            continue
        flipped = flip_record(record)
        orbit = tuple(sorted([record, flipped]))  # type: ignore[assignment]
        seen.update(orbit)
        orbits.append(orbit)  # type: ignore[arg-type]
    return orbits


def majority_success_probability(length: int, error: Fraction) -> Fraction | None:
    if length % 2 == 0:
        return None
    threshold = length // 2 + 1
    success = Fraction(0)
    for correct_copies in range(threshold, length + 1):
        success += Fraction(comb(length, correct_copies)) * ((1 - error) ** correct_copies) * (error ** (length - correct_copies))
    return success


def exact_channel_rows() -> list[dict[str, Any]]:
    rows = []
    for length in range(1, MAX_RECORD_LENGTH + 1):
        records = all_records(length)
        orbits = record_orbits(length)
        for error_numerator in range((ERROR_GRID_DENOMINATOR // 2) + 1):
            error = Fraction(error_numerator, ERROR_GRID_DENOMINATOR)
            equivariance_violations = 0
            symmetric_marginal_violations = 0
            for record in records:
                if channel_probability(record, 0, error) != channel_probability(flip_record(record), 1, error):
                    equivariance_violations += 1
                marginal = Fraction(1, 2) * channel_probability(record, 0, error) + Fraction(1, 2) * channel_probability(record, 1, error)
                flipped_marginal = Fraction(1, 2) * channel_probability(flip_record(record), 0, error) + Fraction(
                    1, 2
                ) * channel_probability(flip_record(record), 1, error)
                if marginal != flipped_marginal:
                    symmetric_marginal_violations += 1
            success = majority_success_probability(length, error)
            rows.append(
                {
                    "record_length": length,
                    "error_probability": fraction_text(error),
                    "record_count": len(records),
                    "record_mirror_orbit_count": len(orbits),
                    "channel_equivariance_violations": equivariance_violations,
                    "symmetric_prior_marginal_violations": symmetric_marginal_violations,
                    "majority_decoder_success_given_seed": None if success is None else fraction_text(success),
                    "majority_decoder_interpretation": (
                        "tie-prone-even-length-record" if success is None else "amplifies/preserves supplied seed bit"
                    ),
                }
            )
    return rows


def invariant_decoder_event_counts(length: int) -> dict[str, Any]:
    orbits = record_orbits(length)
    if length <= ENUMERATED_DECODER_MAX_LENGTH:
        invariant_events = []
        for mask in range(1 << len(orbits)):
            event: set[Record] = set()
            for index, orbit in enumerate(orbits):
                if mask & (1 << index):
                    event.update(orbit)
            invariant_events.append(event)
        d5_dominant = 0
        nonzero_bias = 0
        for event in invariant_events:
            ones_bias = sum(1 for record in event if hamming_weight(record) > length / 2)
            zeros_bias = sum(1 for record in event if hamming_weight(record) < length / 2)
            if ones_bias > zeros_bias:
                d5_dominant += 1
            if ones_bias != zeros_bias:
                nonzero_bias += 1
        return {
            "record_length": length,
            "record_mirror_orbit_count": len(orbits),
            "invariant_decoder_event_count": len(invariant_events),
            "d5_majority_dominant_invariant_event_count": d5_dominant,
            "nonzero_majority_bias_invariant_event_count": nonzero_bias,
        }
    return {
        "record_length": length,
        "record_mirror_orbit_count": len(orbits),
        "invariant_decoder_event_count_formula": f"2^{len(orbits)}",
        "enumeration_skipped_reason": "decoder event count too large; mirror-pair proof used instead",
        "d5_majority_dominant_invariant_event_count": None,
        "nonzero_majority_bias_invariant_event_count": None,
    }


def decoder_rows() -> list[dict[str, Any]]:
    return [invariant_decoder_event_counts(length) for length in range(1, MAX_RECORD_LENGTH + 1)]


def opinion_audit() -> dict[str, Any]:
    return {
        "overall_verdict": "reject-bridge-closure-claim; accept-one-bit-measurement-and-conditional-carrier-claim",
        "claims": [
            {
                "claim": "The Legacy-Strict bridge is fully formed / 99 percent closed.",
                "verdict": "rejected_by_current_repo_guardrails",
                "reason": "Current FAR guardrails and probe hard-limits still say no identity K_legacy_ont == K_strict_gate, no strict-core selector closure, no QW-2191 discharge, and no ToE closure.",
            },
            {
                "claim": "The missing selector has been isolated to one bit separating {1,11} from {5,7}.",
                "verdict": "accepted_as_measured_finite_obstruction",
                "reason": "The unit-orientation probe measures the missing datum as one binary unit-axis record, but does not derive it from strict core.",
            },
            {
                "claim": "Hebbian/self-recorded weights can carry the required subgroup.",
                "verdict": "accepted_conditionally",
                "reason": "The weight-stabilizer self-record audit says learned d5 weights can carry the {1,11} subgroup conditionally while the d5 teacher trace and Hebbian law remain supplied inputs.",
            },
            {
                "claim": "Matter-like self-replicating information may explain the bit source.",
                "verdict": "new_probe_result_conditional_carrier_not_origin",
                "reason": "Mirror-equivariant self-replication preserves/amplifies a supplied bit but cannot create a singleton d5 label from symmetric prior data.",
            },
        ],
    }


def exact_proof_certificate() -> dict[str, str]:
    return {
        "orientation_variable": "B in {0,1}, with B=1 labelled as A5/d5 only after a unit-axis convention is supplied.",
        "record_space": "R_L={0,1}^L, with mirror J(r)=bitwise-not(r).",
        "equivariant_replication_channel": "P(R=r|B=b)=P(R=Jr|B=1-b).",
        "symmetric_prior_consequence": "If P(B=0)=P(B=1)=1/2, then P(R=r)=P(R=Jr); matter records are mirror-balanced in the marginal.",
        "invariant_decoder_consequence": "A full-Aut-invariant decoder/event is a union of record mirror pairs and cannot prefer the 1/d5 side.",
        "shannon_reading": "Replication can increase reliable mutual record of an already supplied B, but Shannon information is correlation with a seed, not creation of the seed label.",
        "selector_consequence": "Calling majority-1 records d5 is a non-invariant decoder unless the missing unit-axis bit has already labelled 1 as d5.",
    }


def main() -> None:
    supports = all_supports()
    histogram_counter = Counter(distance_histogram(support) for support in supports)
    channel_rows = exact_channel_rows()
    decoder_audit_rows = decoder_rows()
    enumerated_decoder_rows = [row for row in decoder_audit_rows if row.get("invariant_decoder_event_count") is not None]
    odd_error_free = [
        row
        for row in channel_rows
        if row["error_probability"] == "0" and row["record_length"] % 2 == 1 and row["majority_decoder_success_given_seed"] == "1"
    ]

    report: dict[str, Any] = {
        "result_kind": "SCRATCH_STRICT_ALPHA_HEBBIAN_MATTER_SHANNON_SELF_REPLICATION_SELECTOR_AUDIT_PROBE__NOT_A_THEOREM",
        "status": "matter-like-self-replication-can-carry-not-originate-unit-orientation-bit",
        "finite_model": {
            "ring": "Z_12",
            "active_count": ACTIVE_COUNT,
            "support_count": len(supports),
            "histogram_class_count": len(histogram_counter),
            "automorphism_units": AUT_UNITS,
            "survivor_axes": [
                {"name": A1, "mode": 1, "distance_histogram_d1_to_d6": CONTIGUOUS_HISTOGRAM},
                {"name": A5, "mode": 5, "distance_histogram_d1_to_d6": D5_HISTOGRAM},
            ],
            "target_q_power": TARGET_Q_POWER,
            "target_eta": TARGET_ETA,
        },
        "opinion_audit": opinion_audit(),
        "matter_replication_channel_audit": {
            "max_record_length": MAX_RECORD_LENGTH,
            "error_grid_denominator": ERROR_GRID_DENOMINATOR,
            "channel_case_count": len(channel_rows),
            "rows": channel_rows,
            "total_channel_equivariance_violations": sum(row["channel_equivariance_violations"] for row in channel_rows),
            "total_symmetric_prior_marginal_violations": sum(row["symmetric_prior_marginal_violations"] for row in channel_rows),
            "error_free_odd_lengths_with_perfect_seed_recovery": [row["record_length"] for row in odd_error_free],
        },
        "invariant_decoder_event_audit": {
            "enumerated_decoder_max_length": ENUMERATED_DECODER_MAX_LENGTH,
            "rows": decoder_audit_rows,
            "enumerated_d5_majority_dominant_total": sum(
                row["d5_majority_dominant_invariant_event_count"] for row in enumerated_decoder_rows
            ),
            "enumerated_nonzero_majority_bias_total": sum(
                row["nonzero_majority_bias_invariant_event_count"] for row in enumerated_decoder_rows
            ),
        },
        "exact_proof_certificate": exact_proof_certificate(),
        "selector_interpretation": {
            "question_tested": "Can matter, read as self-replicating Shannon information, generate the missing d5 unit-orientation bit?",
            "negative_result": "No mirror-equivariant replication channel creates a d5-labelled singleton from a symmetric prior and invariant decoder.",
            "conditional_positive_result": "Matter-like replication can carry and stabilize a supplied orientation bit; in that conditional sense it is a plausible downstream self-record medium.",
            "legacy_strict_bridge_warning": "This does not close the legacy->strict bridge or transfer legacy matter-generation claims onto K_strict_gate.",
            "honest_limit": "This probe is a finite Shannon/replication audit, not a theorem deriving matter, Hebbian learning, or the one-bit source from strict nadsoliton geometry.",
        },
        "ontology_guardrail": {
            "allowed_reading": "Matter may be studied as a downstream self-replicating record internal to the informational nadsoliton route.",
            "forbidden_reading": "No separate informational layer underneath the nadsoliton is introduced, and matter is not placed before the nadsoliton or light in the preferred order.",
            "preferred_order_preserved": "nadsoliton -> light -> matter -> emergent observer",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is asserted or used.",
            "No legacy physical-role or matter-generation claim is transferred onto K_strict_gate.",
            "No theorem derives a matter/self-replication source for the unit-axis bit from strict nadsoliton geometry.",
            "No theorem derives a Hebbian learning law as strict-core dynamics.",
            "No theorem derives the required one-bit unit-axis record from strict core.",
            "Full Aut(Z_12)-equivariant matter-like replication still forbids singleton d5 selection from symmetric prior data.",
            "Any majority-1 -> d5 decoder is classified as a non-strict selector premise unless a new bridge/source theorem is supplied.",
            "No endpoint, arrow orientation, ledger selector, positive lambda action, cycle metric source, anti-Nyquist source, fifth-mode source, future-probability source, future-value source, future-path source, or matter-bit source theorem is claimed.",
            "No QW-2191 discharge and no strict-core selector closure are claimed.",
            "No ToE closure is claimed.",
        ],
        "next_honest_step": "If pursuing matter, search for a strict non-Aut-invariant matter-bit source theorem; otherwise treat matter as conditional carrier/amplifier of a supplied selector bit.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha Hebbian matter/Shannon self-replication selector audit probe\n\n"
        "Status: matter-like self-replication can carry, not originate, the unit-orientation bit.\n\n"
        f"- Supports scanned: `{len(supports)}`; histogram classes: `{len(histogram_counter)}`.\n"
        f"- Replication channel cases audited: `{len(channel_rows)}`; equivariance violations: "
        f"`{report['matter_replication_channel_audit']['total_channel_equivariance_violations']}`; symmetric marginal violations: "
        f"`{report['matter_replication_channel_audit']['total_symmetric_prior_marginal_violations']}`.\n"
        f"- Invariant decoder events enumerated through record length `{ENUMERATED_DECODER_MAX_LENGTH}`; d5-majority dominant invariant events: "
        f"`{report['invariant_decoder_event_audit']['enumerated_d5_majority_dominant_total']}`.\n"
        f"- Error-free odd record lengths with perfect supplied-seed recovery: "
        f"`{report['matter_replication_channel_audit']['error_free_odd_lengths_with_perfect_seed_recovery']}`.\n"
        f"- Opinion audit verdict: `{report['opinion_audit']['overall_verdict']}`.\n"
        f"- Target replay: `q^5={TARGET_Q_POWER}`, eta `{TARGET_ETA}`.\n"
        "- Honest read: matter/self-replication may stabilize a supplied bit but does not derive the bit.\n"
        "- No false pass: no matter-bit source theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
