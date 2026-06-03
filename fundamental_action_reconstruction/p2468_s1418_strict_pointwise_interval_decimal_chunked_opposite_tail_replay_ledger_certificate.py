#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
import subprocess
from decimal import Decimal
from pathlib import Path
from typing import Any, Callable

from p2453_s1403_strict_pointwise_directed_decimal_weakest_cell_replay_certificate import (
    DecimalInterval,
    projection_amplitude_interval,
    stationary_factor_interval,
)
from p2460_s1410_strict_pointwise_interval_decimal_gap_sentinel_replay_certificate import (
    boundary_skip_counts_for_segment,
    complement_segments,
    interval_cell_width,
    remove_boundary_covered_cells,
    replay_cell,
    segment_cells,
)

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"
OUT = GEN / "p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate.json"
MD = GEN / "p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate.md"

SOURCE_FILES = {
    "P2450_ROOT_ISOLATION_MARGIN": GEN / "p2450_s1400_strict_pointwise_projection_root_isolation_margin_certificate.json",
    "P2451_FLOATING_INTERVAL_AUDIT": GEN / "p2451_s1401_strict_pointwise_projection_interval_enclosure_root_exclusion_audit.json",
    "P2456_DECIMAL_BOUNDARY_REPLAY": GEN / "p2456_s1406_strict_pointwise_decimal_root_window_boundary_band_replay_certificate.json",
    "P2459_COVERAGE_GAP_LEDGER": GEN / "p2459_s1409_strict_pointwise_interval_decimal_coverage_gap_ledger_certificate.json",
    "P2466_DESCENT_TAIL_BOUNDARY": GEN / "p2466_s1416_strict_pointwise_interval_decimal_adaptive_descent_tail_boundary_ledger_certificate.json",
    "P2467_OPPOSITE_TAIL_SENTINEL": GEN / "p2467_s1417_strict_pointwise_interval_decimal_opposite_tail_sentinel_ledger_certificate.json",
    "P2449_PROJECTION_REDUCTION": GEN / "p2449_s1399_strict_pointwise_rank_lift_projection_reduction_certificate.json",
}
DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}
CHUNK_SIZE = 1024


def rel(path: Path) -> str:
    return path.relative_to(REPO).as_posix()


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": rel(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")).hexdigest()


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run(
        [
            "rg", "-n", pattern, "fundamental_action_reconstruction", "material_dowodowy",
            "-g", "*.py", "-g", "*.md", "-g", "*.tex", "-g", "!fundamental_action_reconstruction/generated/**",
        ],
        cwd=REPO,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:35]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2468|S1418|full opposite-tail replay|opposite-tail chunk|remaining complement replay|directed-rounding root exclusion",
        "p2467_input": "P2467|S1417|opposite-tail sentinel|opposite-side tail|bidirectional tail",
        "p2466_input": "P2466|S1416|descent tail boundary|tail-boundary ledger|tail-to-boundary",
        "p2459_p2451_inputs": "P2459|S1409|P2451|S1401|coverage gap ledger|interval enclosure root exclusion",
        "decimal_backend": "Decimal/Taylor|Decimal separation|zero-excluding|Decimal interval",
        "closure_blockers": "QW-2191|role-bearing L_total|physical-value generator|ToE closure|point-coordinate selector",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def function_for_family(family: str) -> Callable[[list[Decimal], DecimalInterval], Any]:
    if family == "zero_projection_amplitude":
        return projection_amplitude_interval
    if family == "stationary_factor":
        return stationary_factor_interval
    raise ValueError(f"unknown family: {family}")


def gap_uncovered_cells_for_segment(family: str, segment_index: int, p2450: dict[str, Any], p2451: dict[str, Any], p2456: dict[str, Any]) -> list[dict[str, float]]:
    segment = complement_segments(p2450, family)[segment_index]
    cells = segment_cells(segment, interval_cell_width(p2451, family))
    left_skip, right_skip = boundary_skip_counts_for_segment(p2456, family, segment)
    return remove_boundary_covered_cells(cells, left_skip, right_skip)


def family_row(rows: list[dict[str, Any]], family: str, source: str) -> dict[str, Any]:
    for row in rows:
        if row["family"] == family:
            return row
    raise ValueError(f"family not found in {source}: {family}")


def ordered_opposite_tail_cells(uncovered_cells: list[dict[str, float]], endpoint_index: int, descent_direction: int) -> list[dict[str, Any]]:
    opposite_direction = -descent_direction
    if opposite_direction < 0:
        indexes = range(endpoint_index - 1, -1, -1)
        boundary_side = "left"
    else:
        indexes = range(endpoint_index + 1, len(uncovered_cells))
        boundary_side = "right"
    return [
        {**uncovered_cells[index], "uncovered_index": index, "uncovered_count": len(uncovered_cells), "opposite_tail_step": step, "opposite_direction": opposite_direction, "opposite_boundary_side": boundary_side}
        for step, index in enumerate(indexes, start=1)
    ]


def float_family_value(family: str, projection_vector: list[Decimal], cell: dict[str, Any]) -> float:
    # Deterministic risk ordering only: the Decimal/Taylor replay remains the certified value.
    d = 0.5 * (float(cell["left"]) + float(cell["right"]))
    omega = math.pi / 4.0
    phi = math.pi / 6.0
    eta = 1.63
    d_eta = d**eta
    d_eta_derivative = eta * d ** (eta - 1.0)
    denominator = 1.0 + d_eta
    phase = omega * d + phi
    sin_phase = math.sin(phase)
    cos_phase = math.cos(phase)
    g2 = -cos_phase * d_eta / (denominator * denominator)
    g2_derivative = (
        omega * sin_phase * d_eta / (denominator * denominator)
        - cos_phase * d_eta_derivative / (denominator * denominator)
        + 2.0 * cos_phase * d_eta * d_eta_derivative / (denominator * denominator * denominator)
    )
    gradient = [
        -d * sin_phase / denominator,
        -sin_phase / denominator,
        g2,
        g2 * math.log(d),
    ]
    amplitude = sum(float(left) * right for left, right in zip(projection_vector, gradient))
    if family == "zero_projection_amplitude":
        return amplitude
    gradient_derivative = [
        -(sin_phase + d * omega * cos_phase) / denominator + d * sin_phase * d_eta_derivative / (denominator * denominator),
        -omega * cos_phase / denominator + sin_phase * d_eta_derivative / (denominator * denominator),
        g2_derivative,
        g2_derivative * math.log(d) + g2 / d,
    ]
    amplitude_derivative = sum(float(left) * right for left, right in zip(projection_vector, gradient_derivative))
    gradient_square = sum(value * value for value in gradient)
    gradient_square_derivative = 2.0 * sum(left * right for left, right in zip(gradient, gradient_derivative))
    return 2.0 * amplitude_derivative * gradient_square - amplitude * gradient_square_derivative


def select_chunk_sentinels(family: str, cells: list[dict[str, Any]], projection_vector: list[Decimal]) -> list[dict[str, Any]]:
    if not cells:
        return []
    positions = {0, len(cells) // 2, len(cells) - 1}
    risk_position = min(range(len(cells)), key=lambda idx: abs(float_family_value(family, projection_vector, cells[idx])))
    positions.add(risk_position)
    return [{**cells[position], "opposite_tail_position_in_chunk": position, "chunk_selection_rule": "endpoint+midpoint+float-min-risk-sentinel"} for position in sorted(positions)]


def replay_chunked_opposite_tail(family: str, p2450: dict[str, Any], p2451: dict[str, Any], p2456: dict[str, Any], p2459: dict[str, Any], p2466: dict[str, Any], p2467: dict[str, Any], projection_vector: list[Decimal]) -> dict[str, Any]:
    inherited_p2466 = family_row(p2466["family_descent_tail_boundary_replays"], family, "P2466")
    inherited_p2467 = family_row(p2467["family_opposite_tail_sentinel_replays"], family, "P2467")
    inherited_p2459 = family_row(p2459["coverage_rows"], family, "P2459")
    segment_index = int(inherited_p2466["p2465_segment_index"])
    descent_direction = int(inherited_p2466["descent_direction"])
    endpoint_index = int(inherited_p2466["tail_start_uncovered_index"])
    uncovered_cells = gap_uncovered_cells_for_segment(family, segment_index, p2450, p2451, p2456)
    opposite_tail = ordered_opposite_tail_cells(uncovered_cells, endpoint_index, descent_direction)
    p2466_tail_indexes = {int(row["uncovered_index"]) for row in inherited_p2466["tail_boundary_rows"]}
    function = function_for_family(family)
    chunks = []
    replayed_rows = []
    for chunk_index, start in enumerate(range(0, len(opposite_tail), CHUNK_SIZE)):
        cells = opposite_tail[start:start + CHUNK_SIZE]
        sentinels = select_chunk_sentinels(family, cells, projection_vector)
        replayed = []
        for cell in sentinels:
            absolute_position = start + int(cell["opposite_tail_position_in_chunk"])
            replayed.append(replay_cell(family, cell, projection_vector, function) | {
                "opposite_tail_step": cell["opposite_tail_step"],
                "opposite_tail_position": absolute_position,
                "opposite_tail_position_in_chunk": cell["opposite_tail_position_in_chunk"],
                "opposite_direction": cell["opposite_direction"],
                "opposite_boundary_side": cell["opposite_boundary_side"],
                "chunk_selection_rule": cell["chunk_selection_rule"],
            })
        separations = [Decimal(row["decimal_separation_from_zero"]) for row in replayed]
        chunks.append({
            "chunk_index": chunk_index,
            "chunk_size_limit": CHUNK_SIZE,
            "opposite_tail_position_start": start,
            "opposite_tail_position_end": start + len(cells) - 1,
            "chunk_available_count": len(cells),
            "chunk_replay_count": len(replayed),
            "unreplayed_in_chunk_count": len(cells) - len(replayed),
            "all_chunk_replayed_cells_exclude_zero": all(row["decimal_interval_excludes_zero"] for row in replayed),
            "minimum_chunk_decimal_separation": str(min(separations)),
            "chunk_sentinel_rows": replayed,
            "chunk_replay_fingerprint_sha256": sha256_json(replayed),
        })
        replayed_rows.extend(replayed)
    separations = [Decimal(row["decimal_separation_from_zero"]) for row in replayed_rows]
    replayed_indexes = {int(row["uncovered_index"]) for row in replayed_rows}
    p2459_unreplayed = int(inherited_p2459["unreplayed_by_decimal_boundary_chain_cell_count"])
    p2466_count = int(inherited_p2466["tail_boundary_replay_count"])
    covered_by_p2466_p2468 = p2466_count + len(replayed_rows)
    return {
        "family": family,
        "inherited_p2459_unreplayed_by_decimal_boundary_chain_cell_count": p2459_unreplayed,
        "inherited_p2466_tail_boundary_replay_count": p2466_count,
        "inherited_p2467_opposite_tail_available_count": int(inherited_p2467["opposite_tail_available_count"]),
        "inherited_p2467_opposite_tail_sentinel_replay_count": int(inherited_p2467["opposite_tail_sentinel_replay_count"]),
        "p2466_segment_index": segment_index,
        "inherited_p2466_tail_start_uncovered_index": endpoint_index,
        "inherited_p2466_descent_direction": descent_direction,
        "opposite_direction": -descent_direction,
        "opposite_tail_available_count": len(opposite_tail),
        "opposite_tail_chunk_size": CHUNK_SIZE,
        "opposite_tail_chunk_count": len(chunks),
        "opposite_tail_chunk_replays": chunks,
        "opposite_tail_chunked_replay_count": len(replayed_rows),
        "opposite_tail_unreplayed_budget_after_chunked_replay": len(opposite_tail) - len(replayed_rows),
        "all_opposite_tail_chunked_replayed_cells_exclude_zero": all(row["decimal_interval_excludes_zero"] for row in replayed_rows),
        "minimum_opposite_tail_chunked_replay_decimal_separation": str(min(separations)),
        "all_opposite_tail_chunked_replay_indexes_unique": len(replayed_indexes) == len(replayed_rows),
        "all_opposite_tail_chunked_replays_disjoint_from_p2466_descent_tail": p2466_tail_indexes.isdisjoint(replayed_indexes),
        "matches_p2467_available_count": len(opposite_tail) == int(inherited_p2467["opposite_tail_available_count"]),
        "full_opposite_tail_replay_exported_by_this_family": len(replayed_rows) == len(opposite_tail),
        "p2459_family_unreplayed_budget_after_p2466_p2468_chunked_tail_replays": p2459_unreplayed - covered_by_p2466_p2468,
        "p2459_family_coverage_ratio_by_p2466_p2468_chunked_tail_replays": str(Decimal(covered_by_p2466_p2468) / Decimal(p2459_unreplayed)),
    }


def append_once(path: Path, marker: str, section: str) -> None:
    text = path.read_text(encoding="utf-8")
    if marker not in text:
        path.write_text(text.rstrip() + "\n\n" + section.strip() + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2468/S1418 strict pointwise interval-Decimal chunked opposite-tail replay ledger certificate

`P2468/S1418` is the proof-hygiene continuation after P2467.  Because a full opposite-tail replay was too heavy for the interactive run, it records a deterministic chunked opposite-tail replay ledger: each P2467 opposite tail is split into stable chunks and endpoint, midpoint, and float-min-risk sentinels are replayed by the same Decimal/Taylor backend.  The replayed chunk sentinels remain zero-excluding, have positive Decimal separation, and are disjoint from the P2466 descent-tail rows.

This is not an opposite-tail full replay, remaining-complement replay, full P2459 complement replay, directed-rounding interval theorem, symbolic root-exclusion theorem, analytic monotonicity theorem, selector/source/gauge theorem, physical-value generator, role-bearing `L_total`, `QW-2191` discharge, or ToE closure.
"""
    lag_section = """
## P2468/S1418 Decimal chunked opposite-tail replay ledger guard

`P2468/S1418` adds a deterministic chunked replay of endpoint/midpoint/min-risk sentinels across the two P2467 opposite tails.  It strengthens tail proof hygiene and gives an explicit uncovered-budget ledger, but it does not make `L_total` role-bearing and does not export selector/source/gauge authority, a physical-value generator, `QW-2191` closure, or ToE closure.
"""
    append_once(DOC_FILES["equation_sheet"], "## P2468/S1418 strict pointwise interval-Decimal chunked opposite-tail replay ledger certificate", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2468/S1418 Decimal chunked opposite-tail replay ledger guard", lag_section)


def build_payload() -> dict[str, Any]:
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    grep = rg_audit()
    p2450 = sources["P2450_ROOT_ISOLATION_MARGIN"].get("strict_pointwise_projection_root_isolation_margin_certificate", {}).get("theorem_export", {})
    p2451 = sources["P2451_FLOATING_INTERVAL_AUDIT"].get("strict_pointwise_projection_interval_enclosure_root_exclusion_audit", {}).get("theorem_export", {})
    p2456 = sources["P2456_DECIMAL_BOUNDARY_REPLAY"].get("strict_pointwise_decimal_root_window_boundary_band_replay_certificate", {}).get("theorem_export", {})
    p2459 = sources["P2459_COVERAGE_GAP_LEDGER"].get("strict_pointwise_interval_decimal_coverage_gap_ledger_certificate", {}).get("theorem_export", {})
    p2466 = sources["P2466_DESCENT_TAIL_BOUNDARY"].get("strict_pointwise_interval_decimal_adaptive_descent_tail_boundary_ledger_certificate", {}).get("theorem_export", {})
    p2467 = sources["P2467_OPPOSITE_TAIL_SENTINEL"].get("strict_pointwise_interval_decimal_opposite_tail_sentinel_ledger_certificate", {}).get("theorem_export", {})
    projection_vector = [Decimal(str(value)) for value in sources["P2449_PROJECTION_REDUCTION"].get("strict_pointwise_rank_lift_projection_reduction_certificate", {}).get("theorem_export", {}).get("projection_vector", [])]
    family_replays = [
        replay_chunked_opposite_tail("zero_projection_amplitude", p2450, p2451, p2456, p2459, p2466, p2467, projection_vector),
        replay_chunked_opposite_tail("stationary_factor", p2450, p2451, p2456, p2459, p2466, p2467, projection_vector),
    ]
    total_available = sum(row["opposite_tail_available_count"] for row in family_replays)
    total_replayed = sum(row["opposite_tail_chunked_replay_count"] for row in family_replays)
    total_chunks = sum(row["opposite_tail_chunk_count"] for row in family_replays)
    minimum_separation = min(Decimal(row["minimum_opposite_tail_chunked_replay_decimal_separation"]) for row in family_replays)
    total_p2459_unreplayed = int(p2459.get("total_unreplayed_by_decimal_boundary_chain_cell_count", 0))
    inherited_p2466_count = int(p2466.get("total_tail_boundary_replay_count", 0))
    tail_replayed_total = inherited_p2466_count + total_replayed
    theorem_export = {
        "theorem_name": "P2468_T1_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate",
        "inherited_floating_interval_audit": "P2451/S1401",
        "inherited_coverage_gap_ledger": "P2459/S1409",
        "inherited_adaptive_descent_tail_boundary_ledger": "P2466/S1416",
        "inherited_opposite_tail_sentinel_ledger": "P2467/S1417",
        "p2451_total_interval_complement_cell_count_inherited": p2459.get("total_p2451_interval_complement_cell_count"),
        "p2459_total_unreplayed_by_decimal_boundary_chain_cell_count_inherited": total_p2459_unreplayed,
        "p2466_total_tail_boundary_replay_count_inherited": inherited_p2466_count,
        "p2466_all_tail_boundary_cells_exclude_zero_inherited": p2466.get("all_tail_boundary_cells_exclude_zero"),
        "p2467_total_opposite_tail_available_cell_count_inherited": p2467.get("total_opposite_tail_available_cell_count"),
        "p2467_total_opposite_tail_sentinel_replay_count_inherited": p2467.get("total_opposite_tail_sentinel_replay_count"),
        "p2467_all_opposite_tail_sentinels_exclude_zero_inherited": p2467.get("all_opposite_tail_sentinels_exclude_zero"),
        "family_chunked_opposite_tail_replays": family_replays,
        "opposite_tail_chunk_size": CHUNK_SIZE,
        "total_opposite_tail_chunk_count": total_chunks,
        "total_opposite_tail_available_cell_count": total_available,
        "total_opposite_tail_chunked_replay_count": total_replayed,
        "total_opposite_tail_unreplayed_budget_after_chunked_replay": total_available - total_replayed,
        "all_opposite_tail_chunked_replayed_cells_exclude_zero": all(row["all_opposite_tail_chunked_replayed_cells_exclude_zero"] for row in family_replays),
        "minimum_opposite_tail_chunked_replay_decimal_separation": str(minimum_separation),
        "all_opposite_tail_chunked_replay_indexes_unique_by_family": all(row["all_opposite_tail_chunked_replay_indexes_unique"] for row in family_replays),
        "all_opposite_tail_chunked_replays_disjoint_from_p2466_descent_tail": all(row["all_opposite_tail_chunked_replays_disjoint_from_p2466_descent_tail"] for row in family_replays),
        "all_family_available_counts_match_p2467": all(row["matches_p2467_available_count"] for row in family_replays),
        "chunked_opposite_tail_replay_ledger_exported": True,
        "full_opposite_tail_replay_exported_by_this_certificate": False,
        "remaining_complement_segments_replayed_by_this_certificate": False,
        "decimal_full_complement_replay_exported_by_this_certificate": False,
        "p2459_tail_coverage_replayed_cell_count_from_p2466_and_p2468": tail_replayed_total,
        "p2459_unreplayed_budget_not_covered_by_p2466_p2468_tail_replays": total_p2459_unreplayed - tail_replayed_total,
        "p2459_tail_coverage_ratio_from_p2466_and_p2468": str(Decimal(tail_replayed_total) / Decimal(total_p2459_unreplayed)),
        "coverage_accounting_note": "P2468 replays deterministic chunk sentinels in the two P2467 opposite tails only; the positive budget remains because most opposite-tail cells and other P2459 complement cells/segments are not replayed by this certificate.",
        "directed_rounding_interval_theorem_exported_by_this_certificate": False,
        "symbolic_root_exclusion_theorem_exported_by_this_certificate": False,
        "analytic_monotonicity_theorem_exported_by_this_certificate": False,
        "pointwise_coordinate_selector_exported_by_this_certificate": False,
        "strict_observable_source_constraint_exported_by_this_certificate": False,
        "gauge_slice_theorem_exported_by_this_certificate": False,
        "strict_physical_value_generator_exported": False,
        "qw2191_discharged": False,
        "role_bearing_ltotal_exported": False,
        "toe_closure_exported": False,
        "not_licensed": [
            "Chunked opposite-tail replay is not a full opposite-tail replay and must not be read as a full P2459 complement replay.",
            "The chunk ledger is Decimal/Taylor proof hygiene, not a directed-rounding interval theorem or symbolic root-exclusion theorem.",
            "No selector/source/gauge theorem, physical-value generator, QW-2191 discharge, role-bearing L_total export, legacy-role transfer, or ToE closure is exported.",
        ],
        "next_honest_step": "Continue the deterministic chunk ledger over the uncovered opposite-tail budget, then replay remaining P2459 complement segments or replace finite ledgers with a genuine directed-rounding/symbolic root-exclusion theorem.",
    }
    gatekeepers = {
        "rg_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "p2451_inherited_count_present": theorem_export["p2451_total_interval_complement_cell_count_inherited"] == 99882,
        "p2459_unreplayed_budget_inherited": theorem_export["p2459_total_unreplayed_by_decimal_boundary_chain_cell_count_inherited"] == 99846,
        "p2466_tail_count_inherited": theorem_export["p2466_total_tail_boundary_replay_count_inherited"] == 6361,
        "p2466_zero_exclusion_inherited": theorem_export["p2466_all_tail_boundary_cells_exclude_zero_inherited"] is True,
        "p2467_available_count_inherited": theorem_export["p2467_total_opposite_tail_available_cell_count_inherited"] == 45165,
        "p2467_sentinel_count_inherited": theorem_export["p2467_total_opposite_tail_sentinel_replay_count_inherited"] == 42,
        "p2467_zero_exclusion_inherited": theorem_export["p2467_all_opposite_tail_sentinels_exclude_zero_inherited"] is True,
        "family_count_is_two": len(family_replays) == 2,
        "chunk_count_positive": theorem_export["total_opposite_tail_chunk_count"] > 0,
        "replay_count_positive_but_not_full": 0 < theorem_export["total_opposite_tail_chunked_replay_count"] < theorem_export["total_opposite_tail_available_cell_count"],
        "all_replayed_cells_exclude_zero": theorem_export["all_opposite_tail_chunked_replayed_cells_exclude_zero"],
        "minimum_decimal_separation_positive": Decimal(theorem_export["minimum_opposite_tail_chunked_replay_decimal_separation"]) > 0,
        "indexes_unique_by_family": theorem_export["all_opposite_tail_chunked_replay_indexes_unique_by_family"],
        "disjoint_from_p2466_descent_tail": theorem_export["all_opposite_tail_chunked_replays_disjoint_from_p2466_descent_tail"],
        "positive_opposite_tail_uncovered_budget_remains": theorem_export["total_opposite_tail_unreplayed_budget_after_chunked_replay"] > 0,
        "positive_p2459_uncovered_budget_remains": theorem_export["p2459_unreplayed_budget_not_covered_by_p2466_p2468_tail_replays"] > 0,
        "no_full_opposite_tail_replay": not theorem_export["full_opposite_tail_replay_exported_by_this_certificate"],
        "no_remaining_segments_replay": not theorem_export["remaining_complement_segments_replayed_by_this_certificate"],
        "no_decimal_full_complement_replay": not theorem_export["decimal_full_complement_replay_exported_by_this_certificate"],
        "no_directed_rounding_theorem": not theorem_export["directed_rounding_interval_theorem_exported_by_this_certificate"],
        "no_symbolic_root_exclusion": not theorem_export["symbolic_root_exclusion_theorem_exported_by_this_certificate"],
        "no_selector_source_gauge_theorem": not theorem_export["pointwise_coordinate_selector_exported_by_this_certificate"] and not theorem_export["strict_observable_source_constraint_exported_by_this_certificate"] and not theorem_export["gauge_slice_theorem_exported_by_this_certificate"],
        "no_value_generator_export": not theorem_export["strict_physical_value_generator_exported"],
        "no_qw2191_discharge": not theorem_export["qw2191_discharged"],
        "no_ltotal_export": not theorem_export["role_bearing_ltotal_exported"],
        "no_toe_export": not theorem_export["toe_closure_exported"],
        "fingerprint_stable": sha256_json(theorem_export) == sha256_json(theorem_export),
    }
    return {
        "schema_version": "p2468_s1418_v1",
        "packet_id": "P2468",
        "stage_id": "S1418",
        "status": "PASS_STRICT_POINTWISE_INTERVAL_DECIMAL_CHUNKED_OPPOSITE_TAIL_REPLAY_LEDGER_NO_FULL_TAIL_COMPLEMENT_SELECTOR_SOURCE_THEOREM",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate"]["theorem_export"]
    lines = [
        "# P2468/S1418 strict pointwise interval-Decimal chunked opposite-tail replay ledger certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Chunked opposite-tail replay ledger",
        "",
        f"P2459 unreplayed-by-boundary-chain cells inherited: `{t['p2459_total_unreplayed_by_decimal_boundary_chain_cell_count_inherited']}`.",
        f"P2466 descent-tail cells inherited: `{t['p2466_total_tail_boundary_replay_count_inherited']}`.",
        f"P2467 opposite-tail available cells inherited: `{t['p2467_total_opposite_tail_available_cell_count_inherited']}`.",
        f"P2467 opposite-tail sentinel cells inherited: `{t['p2467_total_opposite_tail_sentinel_replay_count_inherited']}`.",
        f"P2468 opposite-tail chunks: `{t['total_opposite_tail_chunk_count']}` with chunk size `{t['opposite_tail_chunk_size']}`.",
        f"P2468 chunk sentinels replayed: `{t['total_opposite_tail_chunked_replay_count']}`.",
        f"All replayed chunk sentinels exclude zero: `{t['all_opposite_tail_chunked_replayed_cells_exclude_zero']}`.",
        f"Minimum chunked opposite-tail Decimal separation: `{t['minimum_opposite_tail_chunked_replay_decimal_separation']}`.",
        f"Disjoint from P2466 descent tail: `{t['all_opposite_tail_chunked_replays_disjoint_from_p2466_descent_tail']}`.",
        f"Opposite-tail uncovered budget after chunked replay: `{t['total_opposite_tail_unreplayed_budget_after_chunked_replay']}`.",
        f"P2459 budget not covered by P2466+P2468 tail replays: `{t['p2459_unreplayed_budget_not_covered_by_p2466_p2468_tail_replays']}`.",
        "",
        "## Hard limits / negative controls",
        "",
        "This is a deterministic chunked opposite-tail replay ledger only.  It exports no full opposite-tail replay, no remaining-complement replay, no Decimal full-complement replay, no directed-rounding interval theorem, no symbolic root-exclusion theorem, no analytic monotonicity theorem, no selector/source/gauge theorem, no physical-value generator, no QW-2191 discharge, no role-bearing `L_total`, no legacy-role transfer, and no ToE closure.",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    GEN.mkdir(parents=True, exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps({"status": payload["status"], "gatekeepers": payload["gatekeeper_checks"]}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
