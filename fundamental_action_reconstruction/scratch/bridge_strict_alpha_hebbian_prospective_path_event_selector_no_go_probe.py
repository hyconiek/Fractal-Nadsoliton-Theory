#!/usr/bin/env python3
"""Scratch probe: Aut-invariant prospective path events cannot select current d5.

The previous prospective probes checked one-step probabilities and finite-horizon
value iteration on the two-axis survivor orbit

    A1 = (k=1, contiguous histogram),
    A5 = (k=5, d5 histogram).

This probe audits a more expressive "future possible state" object: an entire
finite future path.  A path is a binary word over {A1, A5}; the unit-mirror
Aut(Z_12) action flips every symbol.  A future-event predicate is full-Aut
invariant exactly when it is a union of mirror-paired path orbits.

No false pass: every such invariant event contains equally many paths starting
at A1 and A5, and equally many paths ending at A1 and A5.  Therefore even a rich
faith-like/prospective path record cannot select singleton current d5 unless it
already imports the missing unit-axis bit that chooses one side of every mirror
pair.
"""
from __future__ import annotations

import json
from collections import Counter
from fractions import Fraction
from itertools import combinations, product
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
OUT_JSON = HERE / "bridge_strict_alpha_hebbian_prospective_path_event_selector_no_go_report.json"
OUT_MD = HERE / "bridge_strict_alpha_hebbian_prospective_path_event_selector_no_go_report.md"

N = 12
ACTIVE_COUNT = 5
AUT_UNITS = [1, 5, 7, 11]
TARGET_Q_POWER = "256/243"
TARGET_ETA = "9/5"
MAX_HORIZON = 8
ENUMERATED_EVENT_MAX_HORIZON = 4
A1 = "A1_k1_contiguous"
A5 = "A5_k5_d5"
AXES = [A1, A5]
D5_HISTOGRAM = [0, 3, 2, 1, 4, 0]
CONTIGUOUS_HISTOGRAM = [4, 3, 2, 1, 0, 0]
PathWord = tuple[int, ...]


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


def axis_name(bit: int) -> str:
    return A5 if bit else A1


def path_name(path: PathWord) -> str:
    return "->".join(axis_name(bit) for bit in path)


def mirror_path(path: PathWord) -> PathWord:
    return tuple(1 - bit for bit in path)


def all_paths(horizon: int) -> list[PathWord]:
    return list(product([0, 1], repeat=horizon + 1))


def path_orbits(horizon: int) -> list[tuple[PathWord, PathWord]]:
    seen: set[PathWord] = set()
    orbits = []
    for path in all_paths(horizon):
        if path in seen:
            continue
        mirror = mirror_path(path)
        orbit = tuple(sorted([path, mirror]))  # type: ignore[assignment]
        seen.update(orbit)
        orbits.append(orbit)  # type: ignore[arg-type]
    return orbits


def event_from_orbit_mask(orbits: list[tuple[PathWord, PathWord]], mask: int) -> set[PathWord]:
    event: set[PathWord] = set()
    for index, orbit in enumerate(orbits):
        if mask & (1 << index):
            event.update(orbit)
    return event


def event_counts(event: set[PathWord]) -> dict[str, int]:
    current_a1 = sum(1 for path in event if path[0] == 0)
    current_a5 = sum(1 for path in event if path[0] == 1)
    terminal_a1 = sum(1 for path in event if path[-1] == 0)
    terminal_a5 = sum(1 for path in event if path[-1] == 1)
    return {
        "event_size": len(event),
        "current_A1_path_count": current_a1,
        "current_A5_path_count": current_a5,
        "current_A5_minus_A1": current_a5 - current_a1,
        "terminal_A1_path_count": terminal_a1,
        "terminal_A5_path_count": terminal_a5,
        "terminal_A5_minus_A1": terminal_a5 - terminal_a1,
    }


def is_invariant_event(event: set[PathWord]) -> bool:
    return all(mirror_path(path) in event for path in event)


def cylinder_current_a5(horizon: int) -> set[PathWord]:
    return {path for path in all_paths(horizon) if path[0] == 1}


def cylinder_terminal_a5(horizon: int) -> set[PathWord]:
    return {path for path in all_paths(horizon) if path[-1] == 1}


def exact_horizon_rows() -> list[dict[str, Any]]:
    rows = []
    for horizon in range(MAX_HORIZON + 1):
        paths = all_paths(horizon)
        orbits = path_orbits(horizon)
        current_cylinder = cylinder_current_a5(horizon)
        terminal_cylinder = cylinder_terminal_a5(horizon)
        row: dict[str, Any] = {
            "horizon": horizon,
            "path_length": horizon + 1,
            "path_count": len(paths),
            "mirror_orbit_count": len(orbits),
            "fixed_path_count_under_mirror": sum(1 for path in paths if mirror_path(path) == path),
            "full_aut_invariant_event_count_formula": f"2^{len(orbits)}",
            "current_A5_cylinder_invariant": is_invariant_event(current_cylinder),
            "terminal_A5_cylinder_invariant": is_invariant_event(terminal_cylinder),
            "current_A5_cylinder_counts": event_counts(current_cylinder),
            "terminal_A5_cylinder_counts": event_counts(terminal_cylinder),
            "global_unit_axis_bits_required": 1,
        }
        if horizon <= ENUMERATED_EVENT_MAX_HORIZON:
            invariant_events = [event_from_orbit_mask(orbits, mask) for mask in range(1 << len(orbits))]
            current_biases = [event_counts(event)["current_A5_minus_A1"] for event in invariant_events]
            terminal_biases = [event_counts(event)["terminal_A5_minus_A1"] for event in invariant_events]
            row.update(
                {
                    "enumerated_invariant_event_count": len(invariant_events),
                    "enumerated_current_d5_dominant_event_count": sum(1 for bias in current_biases if bias > 0),
                    "enumerated_terminal_d5_dominant_event_count": sum(1 for bias in terminal_biases if bias > 0),
                    "enumerated_nonzero_current_bias_count": sum(1 for bias in current_biases if bias != 0),
                    "enumerated_nonzero_terminal_bias_count": sum(1 for bias in terminal_biases if bias != 0),
                }
            )
        else:
            row.update(
                {
                    "enumerated_invariant_event_count": None,
                    "enumeration_skipped_reason": "event count too large; exact orbit-pair proof used instead",
                    "enumerated_current_d5_dominant_event_count": None,
                    "enumerated_terminal_d5_dominant_event_count": None,
                    "enumerated_nonzero_current_bias_count": None,
                    "enumerated_nonzero_terminal_bias_count": None,
                }
            )
        rows.append(row)
    return rows


def sample_orbit_certificate(horizon: int) -> dict[str, Any]:
    orbits = path_orbits(horizon)
    samples = []
    for orbit in [orbits[0], orbits[len(orbits) // 2], orbits[-1]]:
        samples.append(
            {
                "paths": [path_name(path) for path in orbit],
                "current_symbols": [axis_name(path[0]) for path in orbit],
                "terminal_symbols": [axis_name(path[-1]) for path in orbit],
                "orbit_has_one_current_A1_and_one_current_A5": sorted(path[0] for path in orbit) == [0, 1],
                "orbit_has_one_terminal_A1_and_one_terminal_A5": sorted(path[-1] for path in orbit) == [0, 1],
            }
        )
    return {"horizon": horizon, "sample_orbits": samples}


def exact_proof_certificate() -> dict[str, str]:
    return {
        "path_space": "Omega_H = {A1,A5}^{H+1}",
        "mirror_action": "J flips every path symbol A1<->A5 and has no fixed path for finite H.",
        "orbit_pair_fact": "Each orbit is {w,Jw}; w and Jw have opposite current symbols and opposite terminal symbols.",
        "invariant_event_condition": "A full-Aut-invariant future event is a union of whole mirror pairs.",
        "current_balance_consequence": "Every invariant event contains the same number of paths with current A1 as with current A5.",
        "terminal_balance_consequence": "Every invariant event contains the same number of paths with terminal A1 as with terminal A5.",
        "selector_consequence": "The current-A5 cylinder chooses exactly one member from every mirror pair, so it is the missing unit-axis bit, not an invariant event.",
    }


def main() -> None:
    supports = all_supports()
    histogram_counter = Counter(distance_histogram(support) for support in supports)
    horizon_rows = exact_horizon_rows()
    enumerated_rows = [row for row in horizon_rows if row["enumerated_invariant_event_count"] is not None]

    report: dict[str, Any] = {
        "result_kind": "SCRATCH_STRICT_ALPHA_HEBBIAN_PROSPECTIVE_PATH_EVENT_SELECTOR_NO_GO_PROBE__NOT_A_THEOREM",
        "status": "aut-invariant-prospective-path-events-cannot-select-current-d5",
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
        "path_event_horizon_audit": {
            "max_horizon": MAX_HORIZON,
            "enumerated_event_max_horizon": ENUMERATED_EVENT_MAX_HORIZON,
            "rows": horizon_rows,
            "enumerated_current_d5_dominant_total": sum(
                row["enumerated_current_d5_dominant_event_count"] for row in enumerated_rows
            ),
            "enumerated_terminal_d5_dominant_total": sum(
                row["enumerated_terminal_d5_dominant_event_count"] for row in enumerated_rows
            ),
            "all_exact_rows_have_noninvariant_current_A5_cylinder": all(
                not row["current_A5_cylinder_invariant"] for row in horizon_rows
            ),
            "all_exact_rows_have_noninvariant_terminal_A5_cylinder": all(
                not row["terminal_A5_cylinder_invariant"] for row in horizon_rows
            ),
        },
        "sample_orbit_certificate": sample_orbit_certificate(3),
        "exact_proof_certificate": exact_proof_certificate(),
        "selector_interpretation": {
            "question_tested": "Can a rich prospective record of possible future paths select d5 without an added unit orientation?",
            "negative_result": "No full-Aut-invariant future-path event has a current-d5 or terminal-d5 count advantage; invariant path predicates are balanced mirror-pair unions.",
            "conditional_positive_result": "The current-d5 path cylinder selects one side of every mirror pair and therefore works only as the same missing one-bit unit-axis premise.",
            "relation_to_highest_future_resonance": "A highest-future-resonance event built from invariant path data remains balanced; if it names the d5 side, the naming is extra selector input.",
            "honest_limit": "This is a path-space symmetry no-go, not a theorem deriving future path records, Hebbian learning, or the unit-axis bit from strict nadsoliton geometry.",
        },
        "ontology_guardrail": {
            "allowed_reading": "The nadsoliton itself may contain information about possible future paths of its own state.",
            "forbidden_reading": "No separate informational layer underneath the nadsoliton is introduced to store path futures.",
            "preferred_order_preserved": "nadsoliton -> light -> matter -> emergent observer",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is asserted or used.",
            "No legacy physical-role claim is transferred onto K_strict_gate.",
            "No theorem derives a faith-like prospective path-record functional from strict nadsoliton geometry.",
            "No theorem derives a Hebbian learning law as strict-core dynamics.",
            "No theorem derives the required one-bit unit-axis record from strict core.",
            "Full Aut(Z_12)-invariant finite future-path events still forbid singleton current-d5 selection.",
            "Any current-d5 future-path cylinder is classified as a non-strict selector premise unless a new bridge/source theorem is supplied.",
            "No endpoint, arrow orientation, ledger selector, positive lambda action, cycle metric source, anti-Nyquist source, fifth-mode source, future-probability source, future-value source, or future-path source theorem is claimed.",
            "No QW-2191 discharge and no strict-core selector closure are claimed.",
            "No ToE closure is claimed.",
        ],
        "next_honest_step": "Search for a strict source of a non-Aut-invariant path-orientation character; otherwise keep future-path/highest-resonance selection explicitly axiom-augmented/non-strict.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha Hebbian prospective path-event selector no-go probe\n\n"
        "Status: full-Aut-invariant future-path events cannot select current d5.\n\n"
        f"- Supports scanned: `{len(supports)}`; histogram classes: `{len(histogram_counter)}`.\n"
        f"- Horizons audited exactly: `0..{MAX_HORIZON}`; invariant events enumerated through horizon `{ENUMERATED_EVENT_MAX_HORIZON}`.\n"
        f"- Enumerated current-d5-dominant invariant events: `{report['path_event_horizon_audit']['enumerated_current_d5_dominant_total']}`.\n"
        f"- Enumerated terminal-d5-dominant invariant events: `{report['path_event_horizon_audit']['enumerated_terminal_d5_dominant_total']}`.\n"
        f"- Current-A5 cylinder invariant for any audited horizon: `{not report['path_event_horizon_audit']['all_exact_rows_have_noninvariant_current_A5_cylinder']}`.\n"
        f"- Target replay: `q^5={TARGET_Q_POWER}`, eta `{TARGET_ETA}`.\n"
        "- Honest read: future-path selection of d5 requires choosing one side of every mirror pair.\n"
        "- No false pass: no future-path source theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
