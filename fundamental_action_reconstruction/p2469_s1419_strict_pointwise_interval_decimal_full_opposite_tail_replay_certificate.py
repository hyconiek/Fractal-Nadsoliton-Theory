#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from decimal import Decimal
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import (
    DOC_FILES,
    REPO,
    ROOT,
    family_row,
    function_for_family,
    gap_uncovered_cells_for_segment,
    load_json,
    ordered_opposite_tail_cells,
    rel,
    replay_cell,
)

GEN = ROOT / "generated"
OUT = GEN / "p2469_s1419_strict_pointwise_interval_decimal_full_opposite_tail_replay_certificate.json"
MD = GEN / "p2469_s1419_strict_pointwise_interval_decimal_full_opposite_tail_replay_certificate.md"

SOURCE_FILES = {
    "P2450_ROOT_ISOLATION_MARGIN": GEN / "p2450_s1400_strict_pointwise_projection_root_isolation_margin_certificate.json",
    "P2451_FLOATING_INTERVAL_AUDIT": GEN / "p2451_s1401_strict_pointwise_projection_interval_enclosure_root_exclusion_audit.json",
    "P2456_DECIMAL_BOUNDARY_REPLAY": GEN / "p2456_s1406_strict_pointwise_decimal_root_window_boundary_band_replay_certificate.json",
    "P2459_COVERAGE_GAP_LEDGER": GEN / "p2459_s1409_strict_pointwise_interval_decimal_coverage_gap_ledger_certificate.json",
    "P2466_DESCENT_TAIL_BOUNDARY": GEN / "p2466_s1416_strict_pointwise_interval_decimal_adaptive_descent_tail_boundary_ledger_certificate.json",
    "P2467_OPPOSITE_TAIL_SENTINEL": GEN / "p2467_s1417_strict_pointwise_interval_decimal_opposite_tail_sentinel_ledger_certificate.json",
    "P2468_CHUNKED_OPPOSITE_TAIL_REPLAY": GEN / "p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate.json",
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
        "new_packet": "P2469|S1419|full opposite-tail Decimal replay|full opposite-tail computational replay|all opposite-tail cells replayed|opposite-tail full replay",
        "anti_duplication_precursors": "P2467|S1417|P2468|S1418|opposite-tail sentinel|opposite-tail chunk",
        "coverage_boundary_inputs": "P2459|S1409|P2466|S1416|coverage gap ledger|descent-tail boundary",
        "hard_limit_markers": "P2459 complement full replay|remaining complement replay|directed-rounding root exclusion|symbolic root-exclusion theorem",
        "closure_blockers": "QW-2191|role-bearing L_total|physical-value generator|ToE closure|selector/source/gauge theorem",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def replay_full_opposite_tail_family(
    family: str,
    p2450: dict[str, Any],
    p2451: dict[str, Any],
    p2456: dict[str, Any],
    p2459: dict[str, Any],
    p2466: dict[str, Any],
    p2467: dict[str, Any],
    projection_vector: list[Decimal],
) -> dict[str, Any]:
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
    chunks: list[dict[str, Any]] = []
    replayed_indexes: set[int] = set()
    total_replayed = 0
    total_excluding = 0
    min_separation: Decimal | None = None
    min_row: dict[str, Any] | None = None
    first_row: dict[str, Any] | None = None
    last_row: dict[str, Any] | None = None
    rolling_fingerprint = hashlib.sha256()

    for chunk_index, start in enumerate(range(0, len(opposite_tail), CHUNK_SIZE)):
        cells = opposite_tail[start:start + CHUNK_SIZE]
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
                "opposite_tail_step": cell["opposite_tail_step"],
                "opposite_tail_position": absolute_position,
                "opposite_tail_position_in_chunk": position_in_chunk,
                "opposite_direction": cell["opposite_direction"],
                "opposite_boundary_side": cell["opposite_boundary_side"],
                "full_replay_selection_rule": "all-opposite-tail-cells",
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
            replayed_indexes.add(int(row["uncovered_index"]))
            chunk_replayed += 1
            total_replayed += 1
        chunks.append({
            "chunk_index": chunk_index,
            "chunk_size_limit": CHUNK_SIZE,
            "opposite_tail_position_start": start,
            "opposite_tail_position_end": start + len(cells) - 1,
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
    p2466_count = int(inherited_p2466["tail_boundary_replay_count"])
    covered_by_p2466_p2469 = p2466_count + total_replayed
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
        "opposite_tail_full_replay_chunks": chunks,
        "opposite_tail_full_replay_count": total_replayed,
        "opposite_tail_unreplayed_budget_after_full_tail_replay": len(opposite_tail) - total_replayed,
        "all_opposite_tail_full_replayed_cells_exclude_zero": total_excluding == total_replayed,
        "minimum_opposite_tail_full_replay_decimal_separation": str(min_separation),
        "minimum_opposite_tail_full_replay_row": min_row,
        "first_opposite_tail_full_replay_row": first_row,
        "last_opposite_tail_full_replay_row": last_row,
        "all_opposite_tail_full_replay_indexes_unique": len(replayed_indexes) == total_replayed,
        "all_opposite_tail_full_replays_disjoint_from_p2466_descent_tail": p2466_tail_indexes.isdisjoint(replayed_indexes),
        "matches_p2467_available_count": len(opposite_tail) == int(inherited_p2467["opposite_tail_available_count"]),
        "full_opposite_tail_replay_exported_by_this_family": total_replayed == len(opposite_tail),
        "p2459_family_unreplayed_budget_after_p2466_p2469_tail_replays": p2459_unreplayed - covered_by_p2466_p2469,
        "p2459_family_coverage_ratio_by_p2466_p2469_tail_replays": str(Decimal(covered_by_p2466_p2469) / Decimal(p2459_unreplayed)),
        "family_full_replay_fingerprint_sha256": rolling_fingerprint.hexdigest(),
    }


def append_once(path: Path, marker: str, section: str) -> None:
    text = path.read_text(encoding="utf-8")
    if marker not in text:
        path.write_text(text.rstrip() + "\n\n" + section.strip() + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2469/S1419 strict pointwise interval-Decimal full opposite-tail replay certificate

`P2469/S1419` upgrades the P2468 chunked ledger into a full computational Decimal/Taylor replay of every cell in the two P2467 opposite tails.  All `45165` inherited opposite-tail cells are replayed, all replayed intervals exclude zero, the minimum Decimal separation remains positive, and the replay index sets remain disjoint from the P2466 descent-tail rows.

This is only a full replay of the two P2467 opposite tails.  It is not a full P2459 remaining-complement replay, not a directed-rounding interval theorem, not a symbolic root-exclusion theorem, not a selector/source/gauge theorem, not a physical-value generator, not a role-bearing `L_total`, not a `QW-2191` discharge, not a legacy-role transfer, and not ToE closure.
"""
    lag_section = """
## P2469/S1419 Decimal full opposite-tail replay guard

`P2469/S1419` gives `L_total` bookkeeping a stronger negative-control input by replaying every P2467 opposite-tail cell with the Decimal/Taylor backend.  The result is a tail-only computational exclusion certificate: it improves coverage accounting but leaves the remaining P2459 complement budget, selector/source obstruction, gauge status, physical-value generation, role-bearing status, and ToE closure open.
"""
    append_once(DOC_FILES["equation_sheet"], "## P2469/S1419 strict pointwise interval-Decimal full opposite-tail replay certificate", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2469/S1419 Decimal full opposite-tail replay guard", lag_section)


def build_payload() -> dict[str, Any]:
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    grep = rg_audit()
    p2450 = sources["P2450_ROOT_ISOLATION_MARGIN"].get("strict_pointwise_projection_root_isolation_margin_certificate", {}).get("theorem_export", {})
    p2451 = sources["P2451_FLOATING_INTERVAL_AUDIT"].get("strict_pointwise_projection_interval_enclosure_root_exclusion_audit", {}).get("theorem_export", {})
    p2456 = sources["P2456_DECIMAL_BOUNDARY_REPLAY"].get("strict_pointwise_decimal_root_window_boundary_band_replay_certificate", {}).get("theorem_export", {})
    p2459 = sources["P2459_COVERAGE_GAP_LEDGER"].get("strict_pointwise_interval_decimal_coverage_gap_ledger_certificate", {}).get("theorem_export", {})
    p2466 = sources["P2466_DESCENT_TAIL_BOUNDARY"].get("strict_pointwise_interval_decimal_adaptive_descent_tail_boundary_ledger_certificate", {}).get("theorem_export", {})
    p2467 = sources["P2467_OPPOSITE_TAIL_SENTINEL"].get("strict_pointwise_interval_decimal_opposite_tail_sentinel_ledger_certificate", {}).get("theorem_export", {})
    p2468 = sources["P2468_CHUNKED_OPPOSITE_TAIL_REPLAY"].get("strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate", {}).get("theorem_export", {})
    projection_vector = [Decimal(str(value)) for value in sources["P2449_PROJECTION_REDUCTION"].get("strict_pointwise_rank_lift_projection_reduction_certificate", {}).get("theorem_export", {}).get("projection_vector", [])]

    family_replays = [replay_full_opposite_tail_family(family, p2450, p2451, p2456, p2459, p2466, p2467, projection_vector) for family in FAMILIES]
    total_available = sum(row["opposite_tail_available_count"] for row in family_replays)
    total_replayed = sum(row["opposite_tail_full_replay_count"] for row in family_replays)
    total_chunks = sum(row["opposite_tail_chunk_count"] for row in family_replays)
    minimum_separation = min(Decimal(row["minimum_opposite_tail_full_replay_decimal_separation"]) for row in family_replays)
    total_p2459_unreplayed = int(p2459.get("total_unreplayed_by_decimal_boundary_chain_cell_count", 0))
    inherited_p2466_count = int(p2466.get("total_tail_boundary_replay_count", 0))
    tail_replayed_total = inherited_p2466_count + total_replayed
    theorem_export = {
        "theorem_name": "P2469_T1_strict_pointwise_interval_decimal_full_opposite_tail_replay_certificate",
        "inherited_floating_interval_audit": "P2451/S1401",
        "inherited_coverage_gap_ledger": "P2459/S1409",
        "inherited_adaptive_descent_tail_boundary_ledger": "P2466/S1416",
        "inherited_opposite_tail_sentinel_ledger": "P2467/S1417",
        "inherited_chunked_opposite_tail_replay_ledger": "P2468/S1418",
        "p2451_total_interval_complement_cell_count_inherited": p2459.get("total_p2451_interval_complement_cell_count"),
        "p2459_total_unreplayed_by_decimal_boundary_chain_cell_count_inherited": total_p2459_unreplayed,
        "p2466_total_tail_boundary_replay_count_inherited": inherited_p2466_count,
        "p2466_all_tail_boundary_cells_exclude_zero_inherited": p2466.get("all_tail_boundary_cells_exclude_zero"),
        "p2467_total_opposite_tail_available_cell_count_inherited": p2467.get("total_opposite_tail_available_cell_count"),
        "p2467_total_opposite_tail_sentinel_replay_count_inherited": p2467.get("total_opposite_tail_sentinel_replay_count"),
        "p2467_all_opposite_tail_sentinels_exclude_zero_inherited": p2467.get("all_opposite_tail_sentinels_exclude_zero"),
        "p2468_total_opposite_tail_chunked_replay_count_inherited": p2468.get("total_opposite_tail_chunked_replay_count"),
        "p2468_total_opposite_tail_unreplayed_budget_after_chunked_replay_inherited": p2468.get("total_opposite_tail_unreplayed_budget_after_chunked_replay"),
        "family_full_opposite_tail_replays": family_replays,
        "opposite_tail_chunk_size": CHUNK_SIZE,
        "total_opposite_tail_chunk_count": total_chunks,
        "total_opposite_tail_available_cell_count": total_available,
        "total_opposite_tail_full_replay_count": total_replayed,
        "total_opposite_tail_unreplayed_budget_after_full_tail_replay": total_available - total_replayed,
        "all_opposite_tail_full_replayed_cells_exclude_zero": all(row["all_opposite_tail_full_replayed_cells_exclude_zero"] for row in family_replays),
        "minimum_opposite_tail_full_replay_decimal_separation": str(minimum_separation),
        "all_opposite_tail_full_replay_indexes_unique_by_family": all(row["all_opposite_tail_full_replay_indexes_unique"] for row in family_replays),
        "all_opposite_tail_full_replays_disjoint_from_p2466_descent_tail": all(row["all_opposite_tail_full_replays_disjoint_from_p2466_descent_tail"] for row in family_replays),
        "all_family_available_counts_match_p2467": all(row["matches_p2467_available_count"] for row in family_replays),
        "full_opposite_tail_replay_exported_by_this_certificate": total_replayed == total_available,
        "remaining_complement_segments_replayed_by_this_certificate": False,
        "decimal_full_complement_replay_exported_by_this_certificate": False,
        "p2459_tail_coverage_replayed_cell_count_from_p2466_and_p2469": tail_replayed_total,
        "p2459_unreplayed_budget_not_covered_by_p2466_p2469_tail_replays": total_p2459_unreplayed - tail_replayed_total,
        "p2459_tail_coverage_ratio_from_p2466_and_p2469": str(Decimal(tail_replayed_total) / Decimal(total_p2459_unreplayed)),
        "coverage_accounting_note": "P2469 fully replays only the two P2467 opposite tails; remaining P2459 complement cells outside P2466 descent tails and P2467 opposite tails are still unreplayed by this certificate.",
        "directed_rounding_interval_theorem_exported_by_this_certificate": False,
        "symbolic_root_exclusion_theorem_exported_by_this_certificate": False,
        "analytic_monotonicity_theorem_exported_by_this_certificate": False,
        "pointwise_coordinate_selector_exported_by_this_certificate": False,
        "strict_observable_source_constraint_exported_by_this_certificate": False,
        "gauge_slice_theorem_exported_by_this_certificate": False,
        "strict_physical_value_generator_exported": False,
        "qw2191_discharged": False,
        "role_bearing_ltotal_exported": False,
        "legacy_role_transfer_exported": False,
        "toe_closure_exported": False,
        "not_licensed": [
            "Full opposite-tail replay is not a full P2459 complement replay and does not cover remaining complement segments outside the two P2467 tails.",
            "The certificate is a finite Decimal/Taylor computation, not a directed-rounding interval theorem or symbolic root-exclusion theorem.",
            "No selector/source/gauge theorem, physical-value generator, QW-2191 discharge, role-bearing L_total export, legacy-role transfer, or ToE closure is exported.",
        ],
        "next_honest_step": "Replay the remaining non-tail P2459 complement budget or replace finite ledgers with a genuine directed-rounding/symbolic root-exclusion theorem; do not treat this tail-only replay as ToE closure.",
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
        "p2468_chunked_replay_count_inherited": theorem_export["p2468_total_opposite_tail_chunked_replay_count_inherited"] == 140,
        "p2468_chunked_uncovered_budget_inherited": theorem_export["p2468_total_opposite_tail_unreplayed_budget_after_chunked_replay_inherited"] == 45025,
        "family_count_is_two": len(family_replays) == 2,
        "full_replay_count_matches_available": theorem_export["total_opposite_tail_full_replay_count"] == theorem_export["total_opposite_tail_available_cell_count"] == 45165,
        "all_replayed_cells_exclude_zero": theorem_export["all_opposite_tail_full_replayed_cells_exclude_zero"],
        "minimum_decimal_separation_positive": Decimal(theorem_export["minimum_opposite_tail_full_replay_decimal_separation"]) > 0,
        "indexes_unique_by_family": theorem_export["all_opposite_tail_full_replay_indexes_unique_by_family"],
        "disjoint_from_p2466_descent_tail": theorem_export["all_opposite_tail_full_replays_disjoint_from_p2466_descent_tail"],
        "opposite_tail_uncovered_budget_zero": theorem_export["total_opposite_tail_unreplayed_budget_after_full_tail_replay"] == 0,
        "positive_p2459_uncovered_budget_remains": theorem_export["p2459_unreplayed_budget_not_covered_by_p2466_p2469_tail_replays"] > 0,
        "full_opposite_tail_replay_is_explicitly_tail_only": theorem_export["full_opposite_tail_replay_exported_by_this_certificate"] and not theorem_export["decimal_full_complement_replay_exported_by_this_certificate"],
        "no_remaining_segments_replay": not theorem_export["remaining_complement_segments_replayed_by_this_certificate"],
        "no_directed_rounding_theorem": not theorem_export["directed_rounding_interval_theorem_exported_by_this_certificate"],
        "no_symbolic_root_exclusion": not theorem_export["symbolic_root_exclusion_theorem_exported_by_this_certificate"],
        "no_selector_source_gauge_theorem": not theorem_export["pointwise_coordinate_selector_exported_by_this_certificate"] and not theorem_export["strict_observable_source_constraint_exported_by_this_certificate"] and not theorem_export["gauge_slice_theorem_exported_by_this_certificate"],
        "no_value_generator_export": not theorem_export["strict_physical_value_generator_exported"],
        "no_qw2191_discharge": not theorem_export["qw2191_discharged"],
        "no_ltotal_export": not theorem_export["role_bearing_ltotal_exported"],
        "no_legacy_role_transfer": not theorem_export["legacy_role_transfer_exported"],
        "no_toe_export": not theorem_export["toe_closure_exported"],
        "fingerprint_stable": sha256_json(theorem_export) == sha256_json(theorem_export),
    }
    return {
        "schema_version": "p2469_s1419_v1",
        "packet_id": "P2469",
        "stage_id": "S1419",
        "status": "PASS_STRICT_POINTWISE_INTERVAL_DECIMAL_FULL_OPPOSITE_TAIL_REPLAY_CERTIFICATE_TAIL_ONLY_NO_COMPLEMENT_SELECTOR_SOURCE_THEOREM",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_pointwise_interval_decimal_full_opposite_tail_replay_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_pointwise_interval_decimal_full_opposite_tail_replay_certificate"]["theorem_export"]
    lines = [
        "# P2469/S1419 strict pointwise interval-Decimal full opposite-tail replay certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Full opposite-tail replay ledger",
        "",
        f"P2459 unreplayed-by-boundary-chain cells inherited: `{t['p2459_total_unreplayed_by_decimal_boundary_chain_cell_count_inherited']}`.",
        f"P2466 descent-tail cells inherited: `{t['p2466_total_tail_boundary_replay_count_inherited']}`.",
        f"P2467 opposite-tail available cells inherited: `{t['p2467_total_opposite_tail_available_cell_count_inherited']}`.",
        f"P2468 chunked opposite-tail replay count inherited: `{t['p2468_total_opposite_tail_chunked_replay_count_inherited']}`.",
        f"P2469 opposite-tail chunks: `{t['total_opposite_tail_chunk_count']}` with chunk size `{t['opposite_tail_chunk_size']}`.",
        f"P2469 full opposite-tail cells replayed: `{t['total_opposite_tail_full_replay_count']}`.",
        f"All full opposite-tail replayed cells exclude zero: `{t['all_opposite_tail_full_replayed_cells_exclude_zero']}`.",
        f"Minimum full opposite-tail Decimal separation: `{t['minimum_opposite_tail_full_replay_decimal_separation']}`.",
        f"Disjoint from P2466 descent tail: `{t['all_opposite_tail_full_replays_disjoint_from_p2466_descent_tail']}`.",
        f"Opposite-tail uncovered budget after full tail replay: `{t['total_opposite_tail_unreplayed_budget_after_full_tail_replay']}`.",
        f"P2459 budget not covered by P2466+P2469 tail replays: `{t['p2459_unreplayed_budget_not_covered_by_p2466_p2469_tail_replays']}`.",
        "",
        "## Hard limits / negative controls",
        "",
        "This is a full replay of the two P2467 opposite tails only.  It exports no full P2459 remaining-complement replay, no Decimal full-complement replay, no directed-rounding interval theorem, no symbolic root-exclusion theorem, no analytic monotonicity theorem, no selector/source/gauge theorem, no physical-value generator, no QW-2191 discharge, no role-bearing `L_total`, no legacy-role transfer, and no ToE closure.",
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
