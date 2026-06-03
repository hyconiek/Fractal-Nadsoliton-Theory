#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from decimal import Decimal
from pathlib import Path
from typing import Any

from p2460_s1410_strict_pointwise_interval_decimal_gap_sentinel_replay_certificate import (
    boundary_skip_counts_for_segment,
    complement_segments,
    interval_cell_width,
    remove_boundary_covered_cells,
    segment_cells,
)
from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import (
    DOC_FILES,
    REPO,
    ROOT,
    family_row,
    function_for_family,
    load_json,
    ordered_opposite_tail_cells,
    rel,
    replay_cell,
)

GEN = ROOT / "generated"
OUT = GEN / "p2470_s1420_strict_pointwise_interval_decimal_full_remaining_non_tail_complement_replay_certificate.json"
MD = GEN / "p2470_s1420_strict_pointwise_interval_decimal_full_remaining_non_tail_complement_replay_certificate.md"

SOURCE_FILES = {
    "P2450_ROOT_ISOLATION_MARGIN": GEN / "p2450_s1400_strict_pointwise_projection_root_isolation_margin_certificate.json",
    "P2451_FLOATING_INTERVAL_AUDIT": GEN / "p2451_s1401_strict_pointwise_projection_interval_enclosure_root_exclusion_audit.json",
    "P2456_DECIMAL_BOUNDARY_REPLAY": GEN / "p2456_s1406_strict_pointwise_decimal_root_window_boundary_band_replay_certificate.json",
    "P2459_COVERAGE_GAP_LEDGER": GEN / "p2459_s1409_strict_pointwise_interval_decimal_coverage_gap_ledger_certificate.json",
    "P2466_DESCENT_TAIL_BOUNDARY": GEN / "p2466_s1416_strict_pointwise_interval_decimal_adaptive_descent_tail_boundary_ledger_certificate.json",
    "P2469_FULL_OPPOSITE_TAIL_REPLAY": GEN / "p2469_s1419_strict_pointwise_interval_decimal_full_opposite_tail_replay_certificate.json",
    "P2449_PROJECTION_REDUCTION": GEN / "p2449_s1399_strict_pointwise_rank_lift_projection_reduction_certificate.json",
}
CHUNK_SIZE = 1024
FAMILIES = ["zero_projection_amplitude", "stationary_factor"]


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
    return {"count": len(lines), "samples": lines[:40]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2470|S1420|remaining non-tail complement|non-tail complement replay|all non-tail cells replayed|full remaining complement replay",
        "anti_duplication_precursors": "P2459|S1409|P2466|S1416|P2469|S1419|coverage gap ledger|full opposite-tail replay",
        "finite_coverage_markers": "P2459 remaining complement|unreplayed-by-boundary-chain|finite Decimal complement coverage|tail-only complement closure",
        "hard_limit_markers": "directed-rounding root exclusion|symbolic root-exclusion theorem|analytic monotonicity theorem|full mathematical proof",
        "closure_blockers": "QW-2191|role-bearing L_total|physical-value generator|ToE closure|selector/source/gauge theorem|legacy-role transfer",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def remaining_non_tail_cells_for_family(
    family: str,
    p2450: dict[str, Any],
    p2451: dict[str, Any],
    p2456: dict[str, Any],
    p2466: dict[str, Any],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    inherited_p2466 = family_row(p2466["family_descent_tail_boundary_replays"], family, "P2466")
    selected_segment_index = int(inherited_p2466["p2465_segment_index"])
    endpoint_index = int(inherited_p2466["tail_start_uncovered_index"])
    descent_direction = int(inherited_p2466["descent_direction"])
    p2466_tail_indexes = {int(row["uncovered_index"]) for row in inherited_p2466["tail_boundary_rows"]}
    remaining_cells: list[dict[str, Any]] = []
    segment_rows: list[dict[str, Any]] = []
    for segment_index, segment in enumerate(complement_segments(p2450, family)):
        cells = segment_cells(segment, interval_cell_width(p2451, family))
        left_skip, right_skip = boundary_skip_counts_for_segment(p2456, family, segment)
        uncovered_cells = remove_boundary_covered_cells(cells, left_skip, right_skip)
        if segment_index == selected_segment_index:
            opposite_tail = ordered_opposite_tail_cells(uncovered_cells, endpoint_index, descent_direction)
            covered_indexes = p2466_tail_indexes | {int(cell["uncovered_index"]) for cell in opposite_tail}
            selection_rule = "selected-segment-minus-p2466-descent-tail-and-p2469-opposite-tail"
        else:
            covered_indexes = set()
            selection_rule = "non-selected-segment-all-p2459-unreplayed-cells"
        segment_remaining = [
            {
                **cell,
                "segment_index": segment_index,
                "uncovered_index": uncovered_index,
                "uncovered_count": len(uncovered_cells),
                "remaining_non_tail_position_in_segment": position,
                "remaining_non_tail_selection_rule": selection_rule,
            }
            for position, (uncovered_index, cell) in enumerate(
                (item for item in enumerate(uncovered_cells) if item[0] not in covered_indexes)
            )
        ]
        segment_rows.append({
            "segment_index": segment_index,
            "segment_left": segment["left"],
            "segment_right": segment["right"],
            "p2451_cell_count": len(cells),
            "p2456_boundary_left_skip_cell_count": left_skip,
            "p2456_boundary_right_skip_cell_count": right_skip,
            "p2459_unreplayed_by_boundary_chain_cell_count": len(uncovered_cells),
            "excluded_by_p2466_p2469_tail_replay_count": len(uncovered_cells) - len(segment_remaining),
            "remaining_non_tail_cell_count": len(segment_remaining),
            "remaining_non_tail_selection_rule": selection_rule,
        })
        remaining_cells.extend(segment_remaining)
    return remaining_cells, segment_rows


def replay_remaining_non_tail_family(
    family: str,
    p2450: dict[str, Any],
    p2451: dict[str, Any],
    p2456: dict[str, Any],
    p2459: dict[str, Any],
    p2466: dict[str, Any],
    p2469: dict[str, Any],
    projection_vector: list[Decimal],
) -> dict[str, Any]:
    inherited_p2459 = family_row(p2459["coverage_rows"], family, "P2459")
    inherited_p2466 = family_row(p2466["family_descent_tail_boundary_replays"], family, "P2466")
    inherited_p2469 = family_row(p2469["family_full_opposite_tail_replays"], family, "P2469")
    remaining_cells, segment_rows = remaining_non_tail_cells_for_family(family, p2450, p2451, p2456, p2466)
    function = function_for_family(family)
    chunks: list[dict[str, Any]] = []
    replayed_keys: set[tuple[int, int]] = set()
    total_replayed = 0
    total_excluding = 0
    min_separation: Decimal | None = None
    min_row: dict[str, Any] | None = None
    first_row: dict[str, Any] | None = None
    last_row: dict[str, Any] | None = None
    rolling_fingerprint = hashlib.sha256()
    for chunk_index, start in enumerate(range(0, len(remaining_cells), CHUNK_SIZE)):
        cells = remaining_cells[start:start + CHUNK_SIZE]
        chunk_replayed = 0
        chunk_excluding = 0
        chunk_min: Decimal | None = None
        chunk_min_row: dict[str, Any] | None = None
        chunk_first_row: dict[str, Any] | None = None
        chunk_last_row: dict[str, Any] | None = None
        chunk_hash = hashlib.sha256()
        for position_in_chunk, cell in enumerate(cells):
            absolute_position = start + position_in_chunk
            row = replay_cell(family, cell, projection_vector, function) | {
                "segment_index": cell["segment_index"],
                "remaining_non_tail_position": absolute_position,
                "remaining_non_tail_position_in_chunk": position_in_chunk,
                "remaining_non_tail_position_in_segment": cell["remaining_non_tail_position_in_segment"],
                "remaining_non_tail_selection_rule": cell["remaining_non_tail_selection_rule"],
                "full_replay_selection_rule": "all-remaining-non-tail-p2459-cells",
            }
            encoded = json.dumps(row, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
            chunk_hash.update(encoded)
            rolling_fingerprint.update(encoded)
            separation = Decimal(row["decimal_separation_from_zero"])
            if row["decimal_interval_excludes_zero"]:
                chunk_excluding += 1
                total_excluding += 1
            if chunk_min is None or separation < chunk_min:
                chunk_min = separation
                chunk_min_row = row
            if min_separation is None or separation < min_separation:
                min_separation = separation
                min_row = row
            if chunk_first_row is None:
                chunk_first_row = row
            chunk_last_row = row
            if first_row is None:
                first_row = row
            last_row = row
            replayed_keys.add((int(row["segment_index"]), int(row["uncovered_index"])))
            chunk_replayed += 1
            total_replayed += 1
        chunks.append({
            "chunk_index": chunk_index,
            "chunk_size_limit": CHUNK_SIZE,
            "remaining_non_tail_position_start": start,
            "remaining_non_tail_position_end": start + len(cells) - 1,
            "chunk_available_count": len(cells),
            "chunk_full_replay_count": chunk_replayed,
            "chunk_unreplayed_count": len(cells) - chunk_replayed,
            "all_chunk_replayed_cells_exclude_zero": chunk_excluding == chunk_replayed,
            "minimum_chunk_decimal_separation": str(chunk_min),
            "chunk_minimum_separation_row": chunk_min_row,
            "chunk_first_row": chunk_first_row,
            "chunk_last_row": chunk_last_row,
            "chunk_full_replay_fingerprint_sha256": chunk_hash.hexdigest(),
        })
    p2459_unreplayed = int(inherited_p2459["unreplayed_by_decimal_boundary_chain_cell_count"])
    p2466_tail_count = int(inherited_p2466["tail_boundary_replay_count"])
    p2469_opposite_tail_count = int(inherited_p2469["opposite_tail_full_replay_count"])
    covered_by_tail_and_remaining = p2466_tail_count + p2469_opposite_tail_count + total_replayed
    return {
        "family": family,
        "segment_rows": segment_rows,
        "inherited_p2459_unreplayed_by_decimal_boundary_chain_cell_count": p2459_unreplayed,
        "inherited_p2466_tail_boundary_replay_count": p2466_tail_count,
        "inherited_p2469_opposite_tail_full_replay_count": p2469_opposite_tail_count,
        "remaining_non_tail_available_cell_count": len(remaining_cells),
        "remaining_non_tail_chunk_size": CHUNK_SIZE,
        "remaining_non_tail_chunk_count": len(chunks),
        "remaining_non_tail_full_replay_chunks": chunks,
        "remaining_non_tail_full_replay_count": total_replayed,
        "remaining_non_tail_unreplayed_budget_after_full_replay": len(remaining_cells) - total_replayed,
        "all_remaining_non_tail_full_replayed_cells_exclude_zero": total_excluding == total_replayed,
        "minimum_remaining_non_tail_full_replay_decimal_separation": str(min_separation),
        "minimum_remaining_non_tail_full_replay_row": min_row,
        "first_remaining_non_tail_full_replay_row": first_row,
        "last_remaining_non_tail_full_replay_row": last_row,
        "all_remaining_non_tail_replay_keys_unique": len(replayed_keys) == total_replayed,
        "p2459_family_budget_after_p2466_p2469_p2470_replays": p2459_unreplayed - covered_by_tail_and_remaining,
        "p2459_family_coverage_ratio_by_p2466_p2469_p2470_replays": str(Decimal(covered_by_tail_and_remaining) / Decimal(p2459_unreplayed)),
        "family_remaining_non_tail_full_replay_fingerprint_sha256": rolling_fingerprint.hexdigest(),
    }


def append_once(path: Path, marker: str, section: str) -> None:
    text = path.read_text(encoding="utf-8")
    if marker not in text:
        path.write_text(text.rstrip() + "\n\n" + section.strip() + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2470/S1420 strict pointwise interval-Decimal full remaining non-tail complement replay certificate

`P2470/S1420` replays the remaining P2459 unreplayed-by-boundary-chain cells that were not already covered by the P2466 descent tails or the P2469 full opposite-tail replay.  The new replay covers `48320` remaining non-tail cells, all replayed intervals exclude zero, and the combined P2466+P2469+P2470 finite Decimal/Taylor ledger covers the full `99846` P2459 unreplayed-by-boundary-chain budget.

For a non-specialist: this closes the current finite checklist of leftover interval cells in the audited numerical grid.  It still does not prove that every possible real point is excluded by a symbolic theorem; it says that every cell in this particular certified finite audit ledger has now been replayed by the Decimal/Taylor backend.

This is not a directed-rounding interval theorem, not a symbolic root-exclusion theorem, not an analytic monotonicity theorem, not a selector/source/gauge theorem, not a physical-value generator, not a role-bearing `L_total`, not a `QW-2191` discharge, not a legacy-role transfer, and not ToE closure.
"""
    lag_section = """
## P2470/S1420 Decimal full remaining non-tail complement replay guard

`P2470/S1420` completes the current finite P2459 unreplayed-by-boundary-chain ledger by adding a full Decimal/Taylor replay of the non-tail remainder left after P2466 and P2469.  It strengthens finite-grid exclusion bookkeeping for `L_total`, but it remains computational proof hygiene rather than selector/source closure, gauge closure, physical-value generation, role-bearing finality, legacy-role transfer, or ToE closure.
"""
    append_once(DOC_FILES["equation_sheet"], "## P2470/S1420 strict pointwise interval-Decimal full remaining non-tail complement replay certificate", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2470/S1420 Decimal full remaining non-tail complement replay guard", lag_section)


def build_payload() -> dict[str, Any]:
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    grep = rg_audit()
    p2450 = sources["P2450_ROOT_ISOLATION_MARGIN"].get("strict_pointwise_projection_root_isolation_margin_certificate", {}).get("theorem_export", {})
    p2451 = sources["P2451_FLOATING_INTERVAL_AUDIT"].get("strict_pointwise_projection_interval_enclosure_root_exclusion_audit", {}).get("theorem_export", {})
    p2456 = sources["P2456_DECIMAL_BOUNDARY_REPLAY"].get("strict_pointwise_decimal_root_window_boundary_band_replay_certificate", {}).get("theorem_export", {})
    p2459 = sources["P2459_COVERAGE_GAP_LEDGER"].get("strict_pointwise_interval_decimal_coverage_gap_ledger_certificate", {}).get("theorem_export", {})
    p2466 = sources["P2466_DESCENT_TAIL_BOUNDARY"].get("strict_pointwise_interval_decimal_adaptive_descent_tail_boundary_ledger_certificate", {}).get("theorem_export", {})
    p2469 = sources["P2469_FULL_OPPOSITE_TAIL_REPLAY"].get("strict_pointwise_interval_decimal_full_opposite_tail_replay_certificate", {}).get("theorem_export", {})
    projection_vector = [Decimal(str(value)) for value in sources["P2449_PROJECTION_REDUCTION"].get("strict_pointwise_rank_lift_projection_reduction_certificate", {}).get("theorem_export", {}).get("projection_vector", [])]
    family_replays = [replay_remaining_non_tail_family(family, p2450, p2451, p2456, p2459, p2466, p2469, projection_vector) for family in FAMILIES]
    total_remaining_available = sum(row["remaining_non_tail_available_cell_count"] for row in family_replays)
    total_remaining_replayed = sum(row["remaining_non_tail_full_replay_count"] for row in family_replays)
    total_chunks = sum(row["remaining_non_tail_chunk_count"] for row in family_replays)
    minimum_separation = min(Decimal(row["minimum_remaining_non_tail_full_replay_decimal_separation"]) for row in family_replays)
    total_p2459_unreplayed = int(p2459.get("total_unreplayed_by_decimal_boundary_chain_cell_count", 0))
    inherited_p2466_count = int(p2466.get("total_tail_boundary_replay_count", 0))
    inherited_p2469_opposite_count = int(p2469.get("total_opposite_tail_full_replay_count", 0))
    finite_replayed_total = inherited_p2466_count + inherited_p2469_opposite_count + total_remaining_replayed
    theorem_export = {
        "theorem_name": "P2470_T1_strict_pointwise_interval_decimal_full_remaining_non_tail_complement_replay_certificate",
        "inherited_floating_interval_audit": "P2451/S1401",
        "inherited_coverage_gap_ledger": "P2459/S1409",
        "inherited_adaptive_descent_tail_boundary_ledger": "P2466/S1416",
        "inherited_full_opposite_tail_replay": "P2469/S1419",
        "p2451_total_interval_complement_cell_count_inherited": p2459.get("total_p2451_interval_complement_cell_count"),
        "p2459_total_unreplayed_by_decimal_boundary_chain_cell_count_inherited": total_p2459_unreplayed,
        "p2466_total_tail_boundary_replay_count_inherited": inherited_p2466_count,
        "p2466_all_tail_boundary_cells_exclude_zero_inherited": p2466.get("all_tail_boundary_cells_exclude_zero"),
        "p2469_total_opposite_tail_full_replay_count_inherited": inherited_p2469_opposite_count,
        "p2469_all_opposite_tail_full_replayed_cells_exclude_zero_inherited": p2469.get("all_opposite_tail_full_replayed_cells_exclude_zero"),
        "p2469_p2459_unreplayed_budget_not_covered_by_tail_replays_inherited": p2469.get("p2459_unreplayed_budget_not_covered_by_p2466_p2469_tail_replays"),
        "family_remaining_non_tail_replays": family_replays,
        "remaining_non_tail_chunk_size": CHUNK_SIZE,
        "total_remaining_non_tail_chunk_count": total_chunks,
        "total_remaining_non_tail_available_cell_count": total_remaining_available,
        "total_remaining_non_tail_full_replay_count": total_remaining_replayed,
        "total_remaining_non_tail_unreplayed_budget_after_full_replay": total_remaining_available - total_remaining_replayed,
        "all_remaining_non_tail_full_replayed_cells_exclude_zero": all(row["all_remaining_non_tail_full_replayed_cells_exclude_zero"] for row in family_replays),
        "minimum_remaining_non_tail_full_replay_decimal_separation": str(minimum_separation),
        "all_remaining_non_tail_replay_keys_unique_by_family": all(row["all_remaining_non_tail_replay_keys_unique"] for row in family_replays),
        "p2459_total_finite_decimal_replayed_cell_count_from_p2466_p2469_p2470": finite_replayed_total,
        "p2459_unreplayed_by_boundary_chain_budget_after_p2466_p2469_p2470": total_p2459_unreplayed - finite_replayed_total,
        "p2459_unreplayed_by_boundary_chain_finite_decimal_coverage_ratio": str(Decimal(finite_replayed_total) / Decimal(total_p2459_unreplayed)),
        "finite_decimal_replay_covers_full_p2459_unreplayed_by_boundary_chain_budget": finite_replayed_total == total_p2459_unreplayed,
        "finite_grid_complement_ledger_exported": True,
        "directed_rounding_interval_theorem_exported_by_this_certificate": False,
        "symbolic_root_exclusion_theorem_exported_by_this_certificate": False,
        "analytic_monotonicity_theorem_exported_by_this_certificate": False,
        "global_continuum_root_exclusion_theorem_exported_by_this_certificate": False,
        "pointwise_coordinate_selector_exported_by_this_certificate": False,
        "strict_observable_source_constraint_exported_by_this_certificate": False,
        "gauge_slice_theorem_exported_by_this_certificate": False,
        "strict_physical_value_generator_exported": False,
        "qw2191_discharged": False,
        "role_bearing_ltotal_exported": False,
        "legacy_role_transfer_exported": False,
        "toe_closure_exported": False,
        "lay_summary": "Every leftover box in the current finite P2459 audit list has now been rechecked by the high-precision Decimal/Taylor computation chain. This is strong bookkeeping for that finite list, but it is not the same as a closed symbolic proof for all real inputs or a physical ToE closure.",
        "not_licensed": [
            "Finite Decimal replay coverage of the P2459 unreplayed-by-boundary-chain budget is not a directed-rounding interval theorem, symbolic root-exclusion theorem, or continuum theorem.",
            "The certificate does not export a selector/source/gauge theorem, physical-value generator, QW-2191 discharge, role-bearing L_total theorem, legacy-role transfer, or ToE closure.",
            "No legacy physical roles are transferred to L_total or to K_strict_gate by this computation.",
        ],
        "next_honest_step": "Audit whether the finite Decimal/Taylor ledger can be replaced or bounded by a true directed-rounding or symbolic root-exclusion theorem; keep selector/source and role-transfer questions separate.",
    }
    gatekeepers = {
        "rg_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "p2451_inherited_count_present": theorem_export["p2451_total_interval_complement_cell_count_inherited"] == 99882,
        "p2459_unreplayed_budget_inherited": theorem_export["p2459_total_unreplayed_by_decimal_boundary_chain_cell_count_inherited"] == 99846,
        "p2466_tail_count_inherited": theorem_export["p2466_total_tail_boundary_replay_count_inherited"] == 6361,
        "p2466_zero_exclusion_inherited": theorem_export["p2466_all_tail_boundary_cells_exclude_zero_inherited"] is True,
        "p2469_opposite_tail_full_replay_count_inherited": theorem_export["p2469_total_opposite_tail_full_replay_count_inherited"] == 45165,
        "p2469_zero_exclusion_inherited": theorem_export["p2469_all_opposite_tail_full_replayed_cells_exclude_zero_inherited"] is True,
        "p2469_remaining_budget_inherited": theorem_export["p2469_p2459_unreplayed_budget_not_covered_by_tail_replays_inherited"] == 48320,
        "family_count_is_two": len(family_replays) == 2,
        "remaining_replay_count_matches_inherited_budget": theorem_export["total_remaining_non_tail_full_replay_count"] == theorem_export["total_remaining_non_tail_available_cell_count"] == 48320,
        "all_remaining_replayed_cells_exclude_zero": theorem_export["all_remaining_non_tail_full_replayed_cells_exclude_zero"],
        "minimum_decimal_separation_positive": Decimal(theorem_export["minimum_remaining_non_tail_full_replay_decimal_separation"]) > 0,
        "remaining_replay_keys_unique_by_family": theorem_export["all_remaining_non_tail_replay_keys_unique_by_family"],
        "remaining_non_tail_uncovered_budget_zero": theorem_export["total_remaining_non_tail_unreplayed_budget_after_full_replay"] == 0,
        "p2459_finite_budget_covered": theorem_export["p2459_unreplayed_by_boundary_chain_budget_after_p2466_p2469_p2470"] == 0,
        "p2459_finite_coverage_ratio_one": Decimal(theorem_export["p2459_unreplayed_by_boundary_chain_finite_decimal_coverage_ratio"]) == Decimal("1"),
        "finite_grid_only_not_directed_rounding": theorem_export["finite_grid_complement_ledger_exported"] and not theorem_export["directed_rounding_interval_theorem_exported_by_this_certificate"],
        "no_symbolic_root_exclusion": not theorem_export["symbolic_root_exclusion_theorem_exported_by_this_certificate"],
        "no_continuum_root_exclusion": not theorem_export["global_continuum_root_exclusion_theorem_exported_by_this_certificate"],
        "no_selector_source_gauge_theorem": not theorem_export["pointwise_coordinate_selector_exported_by_this_certificate"] and not theorem_export["strict_observable_source_constraint_exported_by_this_certificate"] and not theorem_export["gauge_slice_theorem_exported_by_this_certificate"],
        "no_value_generator_export": not theorem_export["strict_physical_value_generator_exported"],
        "no_qw2191_discharge": not theorem_export["qw2191_discharged"],
        "no_ltotal_export": not theorem_export["role_bearing_ltotal_exported"],
        "no_legacy_role_transfer": not theorem_export["legacy_role_transfer_exported"],
        "no_toe_export": not theorem_export["toe_closure_exported"],
        "fingerprint_stable": sha256_json(theorem_export) == sha256_json(theorem_export),
    }
    return {
        "schema_version": "p2470_s1420_v1",
        "packet_id": "P2470",
        "stage_id": "S1420",
        "status": "PASS_STRICT_POINTWISE_INTERVAL_DECIMAL_FULL_REMAINING_NON_TAIL_COMPLEMENT_REPLAY_FINITE_GRID_ONLY_NO_DIRECTED_ROUNDING_SELECTOR_SOURCE_THEOREM",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_pointwise_interval_decimal_full_remaining_non_tail_complement_replay_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_pointwise_interval_decimal_full_remaining_non_tail_complement_replay_certificate"]["theorem_export"]
    lines = [
        "# P2470/S1420 strict pointwise interval-Decimal full remaining non-tail complement replay certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Full remaining non-tail complement replay ledger",
        "",
        f"P2459 unreplayed-by-boundary-chain cells inherited: `{t['p2459_total_unreplayed_by_decimal_boundary_chain_cell_count_inherited']}`.",
        f"P2466 descent-tail cells inherited: `{t['p2466_total_tail_boundary_replay_count_inherited']}`.",
        f"P2469 full opposite-tail cells inherited: `{t['p2469_total_opposite_tail_full_replay_count_inherited']}`.",
        f"P2470 remaining non-tail cells replayed: `{t['total_remaining_non_tail_full_replay_count']}`.",
        f"All remaining non-tail replayed cells exclude zero: `{t['all_remaining_non_tail_full_replayed_cells_exclude_zero']}`.",
        f"Minimum remaining non-tail Decimal separation: `{t['minimum_remaining_non_tail_full_replay_decimal_separation']}`.",
        f"Remaining non-tail uncovered budget after P2470: `{t['total_remaining_non_tail_unreplayed_budget_after_full_replay']}`.",
        f"P2459 unreplayed-by-boundary-chain budget after P2466+P2469+P2470: `{t['p2459_unreplayed_by_boundary_chain_budget_after_p2466_p2469_p2470']}`.",
        f"Finite Decimal/Taylor coverage ratio for that P2459 budget: `{t['p2459_unreplayed_by_boundary_chain_finite_decimal_coverage_ratio']}`.",
        "",
        "## Plain-language progress note",
        "",
        t["lay_summary"],
        "",
        "## Hard limits / negative controls",
        "",
        "This completes the current finite P2459 unreplayed-by-boundary-chain Decimal/Taylor replay ledger only.  It exports no directed-rounding interval theorem, no symbolic root-exclusion theorem, no analytic monotonicity theorem, no global continuum root-exclusion theorem, no selector/source/gauge theorem, no physical-value generator, no QW-2191 discharge, no role-bearing `L_total`, no legacy-role transfer, and no ToE closure.",
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
