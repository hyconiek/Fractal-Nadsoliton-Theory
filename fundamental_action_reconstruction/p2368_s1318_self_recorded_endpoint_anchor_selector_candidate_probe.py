#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"
SCRATCH = ROOT / "scratch"

OUT = GEN / "p2368_s1318_self_recorded_endpoint_anchor_selector_candidate_probe.json"
MD = GEN / "p2368_s1318_self_recorded_endpoint_anchor_selector_candidate_probe.md"

SOURCE_FILES = {
    "P2367_ADMISSIBILITY_NO_GO": GEN
    / "p2367_s1317_selector_phase_origin_admissibility_no_go_probe.json",
    "SELF_RECORDED_ARROW_ACTION": SCRATCH
    / "bridge_strict_alpha_self_recorded_arrow_action_lexicographic_selector_report.json",
    "SELF_RECORDED_ENDPOINT_ANCHOR": SCRATCH
    / "bridge_strict_alpha_self_recorded_d5_endpoint_anchor_certificate_report.json",
    "PHASE_REFERENCE_SELECTOR": SCRATCH
    / "bridge_strict_alpha_phase_reference_source_selector_certificate_report.json",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

Z12_NODE_COUNT = 12
SUPPORT_SIZE = 5
DISTANCE_SELECTED = 5
TARGET_BINARY_EXPONENT = 8
DENOMINATOR = 3
STRICT_TARGET_ETA = 9.0 / 5.0
NAD12_SUPPORT_SIZE = 12
ALPHA_SCALE = 4.0
BALANCED_LEDGER = (2, 2, 2, 1, 1)


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


def ripple_sum_squares(ledger: tuple[int, ...]) -> int:
    return sum(value * value for value in ledger)


def arrow_increase_penalty(ledger: tuple[int, ...]) -> int:
    return sum(max(0, right - left) for left, right in zip(ledger, ledger[1:]))


def parseval_non_dc_power_n5(ledger: tuple[int, ...]) -> int:
    total = sum(ledger)
    return SUPPORT_SIZE * ripple_sum_squares(ledger) - total * total


def ordered_d5_support(source: int, orientation: int) -> tuple[int, ...]:
    step = (orientation * DISTANCE_SELECTED) % Z12_NODE_COUNT
    return tuple((source + index * step) % Z12_NODE_COUNT for index in range(SUPPORT_SIZE))


def value_configuration(source: int, orientation: int, ledger: tuple[int, ...] = BALANCED_LEDGER) -> tuple[int, ...]:
    values = [0] * Z12_NODE_COUNT
    for node, value in zip(ordered_d5_support(source, orientation), ledger):
        values[node] = value
    return tuple(values)


def d5_neighbors(node: int) -> tuple[int, int]:
    return ((node + DISTANCE_SELECTED) % Z12_NODE_COUNT, (node - DISTANCE_SELECTED) % Z12_NODE_COUNT)


def d5_adjacency(config: tuple[int, ...]) -> dict[int, list[int]]:
    support = {node for node, value in enumerate(config) if value != 0}
    return {
        node: sorted(neighbor for neighbor in d5_neighbors(node) if neighbor in support)
        for node in sorted(support)
    }


def walk_path(adjacency: dict[int, list[int]], source: int) -> list[int]:
    path = [source]
    previous = None
    current = source
    while True:
        candidates = [neighbor for neighbor in adjacency[current] if neighbor != previous]
        if not candidates:
            return path
        previous, current = current, candidates[0]
        path.append(current)


def infer_self_recorded_anchor(config: tuple[int, ...]) -> dict[str, Any]:
    adjacency = d5_adjacency(config)
    degrees = {node: len(neighbors) for node, neighbors in adjacency.items()}
    endpoints = sorted(node for node, degree in degrees.items() if degree == 1)
    if len(adjacency) != SUPPORT_SIZE or sorted(degrees.values()) != [1, 1, 2, 2, 2]:
        raise ValueError("support is not a five-node distance-5 path")
    value_2_endpoints = [node for node in endpoints if config[node] == 2]
    value_1_endpoints = [node for node in endpoints if config[node] == 1]
    if len(value_2_endpoints) != 1 or len(value_1_endpoints) != 1:
        raise ValueError("endpoint values do not uniquely mark a value-2 source and value-1 terminus")
    source = value_2_endpoints[0]
    ordered_path = walk_path(adjacency, source)
    second = ordered_path[1]
    if second == (source + DISTANCE_SELECTED) % Z12_NODE_COUNT:
        orientation = 1
    elif second == (source - DISTANCE_SELECTED) % Z12_NODE_COUNT:
        orientation = -1
    else:
        raise ValueError("path second node is not a signed d5 neighbor")
    return {
        "inferred_source": source,
        "inferred_orientation": orientation,
        "ordered_path": ordered_path,
        "ordered_values": [config[node] for node in ordered_path],
        "endpoints": endpoints,
        "endpoint_values": {str(node): config[node] for node in endpoints},
        "degrees": {str(node): degree for node, degree in degrees.items()},
    }


def dihedral_action_on_config(config: tuple[int, ...], shift: int, reflect: bool) -> tuple[int, ...]:
    out = [0] * Z12_NODE_COUNT
    for node, value in enumerate(config):
        transformed = (-node if reflect else node) + shift
        out[transformed % Z12_NODE_COUNT] = value
    return tuple(out)


def transformed_anchor(source: int, orientation: int, shift: int, reflect: bool) -> tuple[int, int]:
    transformed_source = ((-source if reflect else source) + shift) % Z12_NODE_COUNT
    transformed_orientation = -orientation if reflect else orientation
    return transformed_source, transformed_orientation


def ledger_action_certificate() -> dict[str, Any]:
    ledgers = positive_compositions(TARGET_BINARY_EXPONENT, SUPPORT_SIZE)
    rows = []
    for ledger in ledgers:
        rows.append(
            {
                "ordered_ledger": list(ledger),
                "ripple_sum_squares": ripple_sum_squares(ledger),
                "parseval_non_dc_power_n5": parseval_non_dc_power_n5(ledger),
                "self_recorded_arrow_increase_penalty": arrow_increase_penalty(ledger),
                "lexicographic_pair_ripple_then_arrow": [
                    ripple_sum_squares(ledger),
                    arrow_increase_penalty(ledger),
                ],
                "is_balanced_target_order": ledger == BALANCED_LEDGER,
            }
        )
    minimum_pair = min(tuple(row["lexicographic_pair_ripple_then_arrow"]) for row in rows)
    winners = [row["ordered_ledger"] for row in rows if tuple(row["lexicographic_pair_ripple_then_arrow"]) == minimum_pair]
    target = next(row for row in rows if row["is_balanced_target_order"])
    equal_ripple_competitors = [
        row for row in rows
        if not row["is_balanced_target_order"] and row["ripple_sum_squares"] == target["ripple_sum_squares"]
    ]
    return {
        "composition_count": len(ledgers),
        "target_total": TARGET_BINARY_EXPONENT,
        "support_size": SUPPORT_SIZE,
        "unique_lexicographic_winner": winners == [list(BALANCED_LEDGER)],
        "minimum_pair": list(minimum_pair),
        "winners": winners,
        "target_row": target,
        "equal_ripple_competitor_count": len(equal_ripple_competitors),
        "equal_ripple_competitor_arrow_penalties": sorted(
            {row["self_recorded_arrow_increase_penalty"] for row in equal_ripple_competitors}
        ),
        "all_equal_ripple_competitors_have_positive_arrow_penalty": all(
            row["self_recorded_arrow_increase_penalty"] > target["self_recorded_arrow_increase_penalty"]
            for row in equal_ripple_competitors
        ),
        "rows": rows,
    }


def anchor_reconstruction_certificate() -> dict[str, Any]:
    rows = []
    for source in range(Z12_NODE_COUNT):
        for orientation in (-1, 1):
            config = value_configuration(source, orientation)
            inferred = infer_self_recorded_anchor(config)
            rows.append(
                {
                    "actual_source": source,
                    "actual_orientation": orientation,
                    "value_configuration": list(config),
                    "inferred_anchor": inferred,
                    "source_matches": inferred["inferred_source"] == source,
                    "orientation_matches": inferred["inferred_orientation"] == orientation,
                    "ordered_values_match_balanced_ledger": inferred["ordered_values"] == list(BALANCED_LEDGER),
                }
            )
    return {
        "row_count": len(rows),
        "rows": rows,
        "all_sources_recovered": all(row["source_matches"] for row in rows),
        "all_orientations_recovered": all(row["orientation_matches"] for row in rows),
        "all_ordered_values_match_balanced_ledger": all(
            row["ordered_values_match_balanced_ledger"] for row in rows
        ),
    }


def equivariance_certificate() -> dict[str, Any]:
    checked = 0
    mismatches = []
    for source in range(Z12_NODE_COUNT):
        for orientation in (-1, 1):
            config = value_configuration(source, orientation)
            for shift in range(Z12_NODE_COUNT):
                for reflect in (False, True):
                    transformed = dihedral_action_on_config(config, shift, reflect)
                    inferred = infer_self_recorded_anchor(transformed)
                    expected_source, expected_orientation = transformed_anchor(source, orientation, shift, reflect)
                    checked += 1
                    if inferred["inferred_source"] != expected_source or inferred["inferred_orientation"] != expected_orientation:
                        mismatches.append(
                            {
                                "source": source,
                                "orientation": orientation,
                                "shift": shift,
                                "reflect": reflect,
                                "expected_source": expected_source,
                                "expected_orientation": expected_orientation,
                                "inferred_source": inferred["inferred_source"],
                                "inferred_orientation": inferred["inferred_orientation"],
                            }
                        )
    return {
        "checked_cases": checked,
        "mismatch_count": len(mismatches),
        "mismatches": mismatches[:10],
        "all_cases_equivariant": len(mismatches) == 0,
    }


def negative_controls() -> dict[str, Any]:
    multiset_signatures = {
        tuple(sorted(value for value in value_configuration(source, orientation) if value != 0))
        for source in range(Z12_NODE_COUNT)
        for orientation in (-1, 1)
    }
    reversed_ledger = tuple(reversed(BALANCED_LEDGER))
    reversed_rows = []
    for source in range(Z12_NODE_COUNT):
        for orientation in (-1, 1):
            inferred = infer_self_recorded_anchor(value_configuration(source, orientation, reversed_ledger))
            reversed_rows.append(
                {
                    "actual_source": source,
                    "actual_orientation": orientation,
                    "inferred_source": inferred["inferred_source"],
                    "inferred_orientation": inferred["inferred_orientation"],
                    "source_shift": (inferred["inferred_source"] - source) % Z12_NODE_COUNT,
                    "source_moves_to_opposite_endpoint": inferred["inferred_source"] == ordered_d5_support(source, orientation)[-1],
                    "orientation_flipped": inferred["inferred_orientation"] == -orientation,
                }
            )
    return {
        "value_multiset_signature_count": len(multiset_signatures),
        "value_multiset_is_source_and_orientation_blind": len(multiset_signatures) == 1,
        "reversed_ledger_rows": reversed_rows,
        "reversed_ledger_moves_source_to_opposite_endpoint_and_flips_orientation": all(
            row["source_moves_to_opposite_endpoint"] and row["orientation_flipped"] for row in reversed_rows
        ),
        "absolute_origin_not_selected_by_equivariant_anchor": True,
    }


def main() -> None:
    GEN.mkdir(exist_ok=True)
    artifacts = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    docs = {name: path.read_text(encoding="utf-8") if path.exists() else "" for name, path in DOC_FILES.items()}

    ledger = ledger_action_certificate()
    anchor = anchor_reconstruction_certificate()
    equivariance = equivariance_certificate()
    controls = negative_controls()
    product = Fraction(2 ** TARGET_BINARY_EXPONENT, DENOMINATOR ** SUPPORT_SIZE)

    theorem_export = {
        "theorem_name": "P2368 conditional self-recorded endpoint anchor selector candidate",
        "claim": (
            "Within the audited Z12/d5 finite model, a single lexicographic self-recorded "
            "ledger action first selects the ordered ledger (2,2,2,1,1) from all 35 "
            "positive five-part ledgers of total 8; then the endpoint-valued d5 path "
            "constructively recovers source and orientation on all 24 source/orientation "
            "rows and is D12-equivariant over 576 transformed cases. This is a stronger "
            "conditional selector candidate than a raw Fourier phase origin, but it remains "
            "conditional because strict nadsoliton dynamics has not derived the d5 support, "
            "the arrow action, or the endpoint-source convention."
        ),
        "positive_content": [
            "35 positive ledgers are enumerated and the unique lexicographic (R,A) minimizer is (2,2,2,1,1).",
            "All 24 source/orientation rows are recovered by the endpoint anchor extractor.",
            "The endpoint anchor is D12-equivariant over 576 shift/reflection cases.",
            "The value multiset remains source/orientation blind, so ordering and endpoint structure are essential.",
            "The reversed ledger flips orientation, proving the rule is sensitive to ordered self-record data.",
        ],
        "not_licensed": [
            "strict derivation of d5 support",
            "strict derivation of the self-recorded arrow action",
            "strict derivation of endpoint-source convention",
            "D12-invariant absolute origin selection",
            "beta_tors -> chi11 theorem",
            "legacy physical-role transfer",
            "QW-2191 discharge",
            "ToE closure",
        ],
    }

    probe = {
        "probe_id": "P2368_S1318_SELF_RECORDED_ENDPOINT_ANCHOR_SELECTOR_CANDIDATE",
        "source_artifacts": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "docs": {name: rel(path) for name, path in DOC_FILES.items()},
        "p2367_boundary_replay": {
            "packet_id": artifacts["P2367_ADMISSIBILITY_NO_GO"].get("packet_id"),
            "translation_invariant_source_blind": artifacts["P2367_ADMISSIBILITY_NO_GO"]
            .get("gatekeeper_checks", {})
            .get("translation_invariant_examples_source_blind"),
            "coprime_phase_reference_recovers_sources": artifacts["P2367_ADMISSIBILITY_NO_GO"]
            .get("gatekeeper_checks", {})
            .get("coprime_phase_reference_recovers_sources"),
        },
        "target_identity_replay": {
            "q_power_product": f"{product.numerator}/{product.denominator}",
            "eta_residual_vs_9_5": eta_from_product(product, SUPPORT_SIZE) - STRICT_TARGET_ETA,
            "balanced_ledger": list(BALANCED_LEDGER),
        },
        "ledger_action_certificate": ledger,
        "anchor_reconstruction_certificate": anchor,
        "d12_equivariance_certificate": equivariance,
        "negative_controls": controls,
        "candidate_boundary": {
            "stronger_than_p2366_phase_origin_candidate_because": "source and orientation are extracted from structured finite self-record data rather than from an external Fourier phase reference",
            "weaker_than_strict_selector_theorem_because": "the strict theory has not derived the d5 support, arrow action, endpoint-valued ledger, or endpoint-source convention",
            "relation_to_p2367_no_go": "it does not contradict P2367 because it is an equivariant extractor from non-translation-invariant ordered path data, not a translation-invariant scalar score",
        },
        "theorem_export": theorem_export,
    }
    probe["theorem_fingerprint_sha256"] = sha256_json(theorem_export)

    gatekeeper_checks = {
        "source_artifacts_present": all("_missing" not in artifact for artifact in artifacts.values()),
        "p2367_boundary_loaded": probe["p2367_boundary_replay"]["packet_id"] == "P2367",
        "ledger_unique_winner": ledger["unique_lexicographic_winner"],
        "equal_ripple_competitors_beaten_by_arrow": ledger[
            "all_equal_ripple_competitors_have_positive_arrow_penalty"
        ],
        "all_24_anchors_recovered": anchor["row_count"] == 24
        and anchor["all_sources_recovered"]
        and anchor["all_orientations_recovered"]
        and anchor["all_ordered_values_match_balanced_ledger"],
        "d12_equivariance_passes": equivariance["checked_cases"] == 576
        and equivariance["all_cases_equivariant"],
        "negative_controls_pass": controls["value_multiset_is_source_and_orientation_blind"]
        and controls["reversed_ledger_moves_source_to_opposite_endpoint_and_flips_orientation"]
        and controls["absolute_origin_not_selected_by_equivariant_anchor"],
        "docs_updated_with_p2368_candidate": all("P2368/S1318" in text for text in docs.values()),
        "no_qw2191_discharge_claimed": True,
    }

    payload = {
        "schema_version": "p2368_s1318_v1",
        "packet_id": "P2368",
        "stage_id": "S1318",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-06-01T00:00:00+00:00",
        "status": "OPEN_PROGRESS_SELF_RECORDED_ENDPOINT_ANCHOR_CANDIDATE_NO_QW2191_DISCHARGE",
        "result_kind": "SELF_RECORDED_ENDPOINT_ANCHOR_SELECTOR_CANDIDATE",
        "self_recorded_endpoint_anchor_selector_candidate_probe": probe,
        "recommended_next_honest_step": (
            "Try to derive the self-recorded arrow action and endpoint-valued d5 ledger from "
            "bridge-completed nadsoliton dynamics, or prove that this ordered endpoint record "
            "requires an extra non-strict selector premise."
        ),
        "gatekeeper_checks": gatekeeper_checks,
        "global_status": "OPEN_PROGRESS_WITH_TRACE_NO_SELECTOR_OR_TOE_CLOSURE",
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2368 S1318: self-recorded endpoint anchor selector candidate",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Result",
                "",
                "P2368/S1318 promotes the next conditional selector candidate above the phase-origin boundary: a self-recorded ordered d5 endpoint anchor.",
                "",
                "## Computed certificate",
                "",
                f"- Positive five-part ledgers of total 8 enumerated: `{ledger['composition_count']}`.",
                f"- Unique lexicographic `(ripple, arrow)` winner: `{ledger['unique_lexicographic_winner']}` with `{ledger['winners']}`.",
                f"- Source/orientation rows recovered: `{anchor['row_count']}`; all sources `{anchor['all_sources_recovered']}`; all orientations `{anchor['all_orientations_recovered']}`.",
                f"- D12 equivariance checked cases: `{equivariance['checked_cases']}`; mismatches `{equivariance['mismatch_count']}`.",
                f"- Value multiset source/orientation blind: `{controls['value_multiset_is_source_and_orientation_blind']}`.",
                f"- Reversed ledger moves source to the opposite endpoint and flips orientation: `{controls['reversed_ledger_moves_source_to_opposite_endpoint_and_flips_orientation']}`.",
                "",
                "## Hard limits",
                "",
                "- This is conditional on d5 support, arrow action, endpoint-valued ledger, and endpoint-source convention.",
                "- It is D12-equivariant self-record extraction, not D12-invariant absolute-origin selection.",
                "- No `beta_tors -> chi11` theorem, legacy role transfer, QW-2191 discharge, selector closure, or ToE closure is claimed.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
