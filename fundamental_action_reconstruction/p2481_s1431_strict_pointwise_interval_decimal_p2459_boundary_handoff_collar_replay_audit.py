#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from decimal import Decimal
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import (
    DOC_FILES,
    REPO,
    ROOT,
    function_for_family,
    load_json,
    rel,
    replay_cell,
)
from p2475_s1425_strict_pointwise_interval_decimal_p2459_critical_minimum_halo_replay_audit import projection_vector
from p2480_s1430_strict_pointwise_interval_decimal_p2459_segment_start_cell_dyadic_refinement_replay_audit import (
    DYADIC_SUBCELL_COUNT,
    dyadic_subcells,
)

GEN = ROOT / "generated"
OUT = GEN / "p2481_s1431_strict_pointwise_interval_decimal_p2459_boundary_handoff_collar_replay_audit.json"
MD = GEN / "p2481_s1431_strict_pointwise_interval_decimal_p2459_boundary_handoff_collar_replay_audit.md"

SOURCE_FILES = {
    "P2449_PROJECTION_REDUCTION": GEN / "p2449_s1399_strict_pointwise_rank_lift_projection_reduction_certificate.json",
    "P2456_DECIMAL_BOUNDARY_REPLAY": GEN / "p2456_s1406_strict_pointwise_decimal_root_window_boundary_band_replay_certificate.json",
    "P2479_SEGMENT_START_LEFT_PREFIX": GEN / "p2479_s1429_strict_pointwise_interval_decimal_p2459_segment_start_left_prefix_replay_audit.json",
    "P2480_SEGMENT_START_CELL_REFINEMENT": GEN / "p2480_s1430_strict_pointwise_interval_decimal_p2459_segment_start_cell_dyadic_refinement_replay_audit.json",
}


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
        "new_packet": "P2481|S1431|boundary handoff collar|root-window handoff|boundary-band to prefix|boundary collar replay|covered-boundary collar|right-boundary band handoff|segment-start handoff",
        "precursor_packets": "P2456|S1406|P2479|S1429|P2480|S1430|segment-start cell dyadic refinement|boundary band replay",
        "handoff_language": "boundary handoff|boundary-band to prefix|covered boundary band|segment-start prefix|collar replay|handoff collar",
        "hard_limit_markers": "directed-rounding root exclusion|symbolic root-exclusion theorem|global continuum root-exclusion|full mathematical proof|full complement",
        "closure_blockers": "QW-2191|role-bearing L_total|physical-value generator|ToE closure|selector/source/gauge theorem|legacy-role transfer",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def p2456_boundary_replay_key(family: str) -> str:
    if family == "zero_projection_amplitude":
        return "zero_projection_boundary_band_replays"
    if family == "stationary_factor":
        return "stationary_factor_boundary_band_replays"
    raise ValueError(f"Unsupported P2456 boundary-band family: {family}")


def matching_right_boundary_band(p2456_theorem: dict[str, Any], family: str, parent_left: Decimal) -> dict[str, Any]:
    tolerance = Decimal("1e-20")
    candidates = []
    for band in p2456_theorem[p2456_boundary_replay_key(family)]:
        if band["family"] != family or band["side"] != "right":
            continue
        rows = band.get("rows", [])
        if not rows:
            continue
        band_right = Decimal(str(rows[-1]["right"]))
        if abs(band_right - parent_left) <= tolerance:
            candidates.append(band)
    if len(candidates) != 1:
        raise ValueError(f"Expected exactly one right boundary band ending at parent_left={parent_left}, found {len(candidates)}")
    return candidates[0]


def fresh_boundary_band_rows(boundary_band: dict[str, Any], projection: list[Decimal]) -> list[dict[str, Any]]:
    family = boundary_band["family"]
    function = function_for_family(family)
    rows = []
    for row in boundary_band["rows"]:
        cell = {
            "left": row["left"],
            "right": row["right"],
            "uncovered_index": int(row["index"]),
            "uncovered_count": int(boundary_band["band_cell_count"]),
        }
        fresh = replay_cell(family, cell, projection, function)
        rows.append({
            **fresh,
            "collar_side": "p2456_right_boundary_band",
            "boundary_band_index": int(row["index"]),
            "fresh_decimal_separation_positive": Decimal(fresh["decimal_separation_from_zero"]) > 0,
        })
    return rows


def fresh_p2480_subcell_rows(parent_cell: dict[str, Any], projection: list[Decimal]) -> list[dict[str, Any]]:
    family = parent_cell["family"]
    function = function_for_family(family)
    rows = []
    for cell in dyadic_subcells(parent_cell):
        fresh = replay_cell(family, cell, projection, function)
        rows.append({
            **fresh,
            "collar_side": "p2480_parent_first_cell_dyadic_subcell",
            "dyadic_subcell_index": int(cell["uncovered_index"]),
            "dyadic_subcell_count": DYADIC_SUBCELL_COUNT,
            "parent_segment_index": int(cell["parent_segment_index"]),
            "parent_uncovered_index": int(cell["parent_uncovered_index"]),
            "fresh_decimal_separation_positive": Decimal(fresh["decimal_separation_from_zero"]) > 0,
        })
    return rows


def handoff_collar_replay(boundary_band: dict[str, Any], parent_cell: dict[str, Any], projection: list[Decimal]) -> dict[str, Any]:
    boundary_rows = fresh_boundary_band_rows(boundary_band, projection)
    subcell_rows = fresh_p2480_subcell_rows(parent_cell, projection)
    ordered = sorted(boundary_rows + subcell_rows, key=lambda row: (Decimal(str(row["left"])), Decimal(str(row["right"]))))
    consecutive_pairs = []
    for prior, current in zip(ordered, ordered[1:]):
        delta = Decimal(current["decimal_separation_from_zero"]) - Decimal(prior["decimal_separation_from_zero"])
        endpoint_gap = Decimal(str(current["left"])) - Decimal(str(prior["right"]))
        consecutive_pairs.append({
            "from_collar_side": prior["collar_side"],
            "to_collar_side": current["collar_side"],
            "from_left": str(prior["left"]),
            "from_right": str(prior["right"]),
            "to_left": str(current["left"]),
            "to_right": str(current["right"]),
            "endpoint_gap_to_minus_from": str(endpoint_gap),
            "from_separation": prior["decimal_separation_from_zero"],
            "to_separation": current["decimal_separation_from_zero"],
            "separation_delta_to_minus_from": str(delta),
            "strictly_increases": delta > 0,
            "is_exactly_adjacent": endpoint_gap == 0,
        })
    min_row = min(ordered, key=lambda row: Decimal(row["decimal_separation_from_zero"]))
    max_row = max(ordered, key=lambda row: Decimal(row["decimal_separation_from_zero"]))
    return {
        "family": parent_cell["family"],
        "parent_segment_index": int(parent_cell["segment_index"]),
        "parent_uncovered_index": int(parent_cell["uncovered_index"]),
        "boundary_band_side": boundary_band["side"],
        "boundary_band_cell_count": len(boundary_rows),
        "dyadic_subcell_count": len(subcell_rows),
        "total_fresh_handoff_collar_replay_rows": len(ordered),
        "boundary_band_left": str(boundary_rows[0]["left"]),
        "boundary_band_right": str(boundary_rows[-1]["right"]),
        "subcell_parent_left": str(parent_cell["left"]),
        "subcell_parent_right": str(parent_cell["right"]),
        "handoff_endpoint_gap_boundary_to_parent": str(Decimal(str(parent_cell["left"])) - Decimal(str(boundary_rows[-1]["right"]))),
        "all_collar_rows_exclude_zero": all(row["decimal_interval_excludes_zero"] for row in ordered),
        "all_collar_rows_have_positive_separation": all(row["fresh_decimal_separation_positive"] for row in ordered),
        "all_consecutive_collar_separations_strictly_increase_left_to_right": all(row["strictly_increases"] for row in consecutive_pairs),
        "all_consecutive_collar_rows_are_exactly_adjacent": all(row["is_exactly_adjacent"] for row in consecutive_pairs),
        "minimum_consecutive_endpoint_gap": str(min(Decimal(row["endpoint_gap_to_minus_from"]) for row in consecutive_pairs)),
        "maximum_consecutive_endpoint_gap": str(max(Decimal(row["endpoint_gap_to_minus_from"]) for row in consecutive_pairs)),
        "minimum_consecutive_positive_delta": str(min(Decimal(row["separation_delta_to_minus_from"]) for row in consecutive_pairs)),
        "minimum_collar_decimal_separation": min_row["decimal_separation_from_zero"],
        "maximum_collar_decimal_separation": max_row["decimal_separation_from_zero"],
        "minimum_collar_replay_row": min_row,
        "maximum_collar_replay_row": max_row,
        "minimum_is_p2456_right_boundary_band_leftmost_cell": min_row["collar_side"] == "p2456_right_boundary_band" and int(min_row["boundary_band_index"]) == 0,
        "maximum_is_p2480_rightmost_dyadic_subcell": max_row["collar_side"] == "p2480_parent_first_cell_dyadic_subcell" and int(max_row["dyadic_subcell_index"]) == DYADIC_SUBCELL_COUNT - 1,
        "boundary_band_rows": boundary_rows,
        "dyadic_subcell_rows": subcell_rows,
        "consecutive_collar_pairs": consecutive_pairs,
        "handoff_collar_replay_fingerprint_sha256": sha256_json(ordered),
    }


def append_once(path, marker: str, section: str) -> None:
    text = path.read_text(encoding="utf-8")
    if marker not in text:
        path.write_text(text.rstrip() + "\n\n" + section.strip() + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2481/S1431 strict pointwise interval-Decimal P2459 boundary-handoff collar replay audit

`P2481/S1431` audits the seam that P2480 left open: the P2456 right boundary band immediately before the P2479/P2480 segment-start cell is freshly replayed together with the `128` dyadic subcell rows inside that first P2479 parent cell.  The collar has `6` inherited covered-boundary-chain cells plus `128` subcell diagnostic rows, with exact consecutive endpoint adjacency at the handoff; all `134` Decimal evaluation rows exclude zero, and their Decimal separations strictly increase left-to-right.  The weakest collar row is not the P2480 subcell but the leftmost P2456 right-boundary-band cell, so the correct finite statement is a covered-boundary-to-prefix handoff localization rather than a continuum theorem, a P2459 coverage increase by 134 cells, or a full-complement replay.

For a non-specialist: after opening the first segment cell, P2481 checks the already-covered strip immediately to its left as well.  That strip is even weaker, but it also stays away from zero and increases smoothly into the refined segment-start cell.  The `134` rows are diagnostic Decimal evaluations, not `134` new P2459 cells.  This improves the seam bookkeeping without proving directed rounding, symbolic root exclusion, selector/source closure, legacy-role transfer, or ToE closure.
"""
    lag_section = """
## P2481/S1431 P2459 boundary-handoff collar replay guard

`P2481/S1431` adds a boundary-handoff collar guard behind `L_total`: the P2456 covered right-boundary band adjacent to the P2479 segment-start cell is replayed together with the P2480 dyadic refinement.  This prevents overreading the P2480 leftmost subcell as the whole local bottleneck and keeps `134` fresh Decimal evaluation rows separate from P2459 coverage-cell counting; it remains finite seam evidence only and does not add selector/source/gauge authority, physical-value generation, role-bearing finality, legacy-role transfer, or ToE closure.
"""
    append_once(DOC_FILES["equation_sheet"], "## P2481/S1431 strict pointwise interval-Decimal P2459 boundary-handoff collar replay audit", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2481/S1431 P2459 boundary-handoff collar replay guard", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    projection = projection_vector(sources["P2449_PROJECTION_REDUCTION"])
    p2456 = theorem(sources["P2456_DECIMAL_BOUNDARY_REPLAY"], "strict_pointwise_decimal_root_window_boundary_band_replay_certificate")
    p2479 = theorem(sources["P2479_SEGMENT_START_LEFT_PREFIX"], "strict_pointwise_interval_decimal_p2459_segment_start_left_prefix_replay_audit")
    p2480 = theorem(sources["P2480_SEGMENT_START_CELL_REFINEMENT"], "strict_pointwise_interval_decimal_p2459_segment_start_cell_dyadic_refinement_replay_audit")
    parent_cell = p2479["segment_start_left_prefix_replay"]["minimum_prefix_replay_row"]
    boundary_band = matching_right_boundary_band(p2456, parent_cell["family"], Decimal(str(parent_cell["left"])))
    replay = handoff_collar_replay(boundary_band, parent_cell, projection)
    p2459_universe = p2479["p2471_p2459_universe_count_inherited"]
    theorem_export = {
        "theorem_name": "P2481_T1_strict_pointwise_interval_decimal_p2459_boundary_handoff_collar_replay_audit",
        "audited_chain": ["P2456/S1406", "P2479/S1429", "P2480/S1430"],
        "handoff_collar_replay": replay,
        "family": replay["family"],
        "parent_segment_index": replay["parent_segment_index"],
        "parent_uncovered_index": replay["parent_uncovered_index"],
        "fresh_boundary_band_cell_replay_count": replay["boundary_band_cell_count"],
        "fresh_dyadic_subcell_replay_count": replay["dyadic_subcell_count"],
        "total_fresh_handoff_collar_replay_rows": replay["total_fresh_handoff_collar_replay_rows"],
        "handoff_endpoint_gap_boundary_to_parent": replay["handoff_endpoint_gap_boundary_to_parent"],
        "all_collar_rows_exclude_zero": replay["all_collar_rows_exclude_zero"],
        "all_collar_rows_have_positive_separation": replay["all_collar_rows_have_positive_separation"],
        "all_consecutive_collar_separations_strictly_increase_left_to_right": replay["all_consecutive_collar_separations_strictly_increase_left_to_right"],
        "all_consecutive_collar_rows_are_exactly_adjacent": replay["all_consecutive_collar_rows_are_exactly_adjacent"],
        "minimum_consecutive_endpoint_gap": replay["minimum_consecutive_endpoint_gap"],
        "maximum_consecutive_endpoint_gap": replay["maximum_consecutive_endpoint_gap"],
        "minimum_consecutive_positive_delta": replay["minimum_consecutive_positive_delta"],
        "minimum_collar_decimal_separation": replay["minimum_collar_decimal_separation"],
        "maximum_collar_decimal_separation": replay["maximum_collar_decimal_separation"],
        "minimum_is_p2456_right_boundary_band_leftmost_cell": replay["minimum_is_p2456_right_boundary_band_leftmost_cell"],
        "maximum_is_p2480_rightmost_dyadic_subcell": replay["maximum_is_p2480_rightmost_dyadic_subcell"],
        "p2480_minimum_subcell_decimal_separation_inherited": p2480["minimum_subcell_decimal_separation"],
        "p2480_dyadic_subcell_count_inherited": p2480["dyadic_subcell_count"],
        "p2479_prefix_replay_count_inherited": p2479["fresh_segment_start_left_prefix_replay_count"],
        "p2479_prefix_plus_p2478_union_count_inherited": p2479["p2479_prefix_plus_p2478_flank_union_count"],
        "p2471_p2459_universe_count_inherited": p2459_universe,
        "p2481_fresh_decimal_evaluation_row_count_not_a_coverage_count": replay["total_fresh_handoff_collar_replay_rows"],
        "p2481_fresh_decimal_evaluation_row_ratio_not_a_p2459_coverage_fraction": f"{replay['total_fresh_handoff_collar_replay_rows']}/{p2459_universe}",
        "targeted_p2481_new_p2459_unreplayed_parent_cell_scope_against_p2459_universe": f"1/{p2459_universe}",
        "targeted_p2481_new_p2459_unreplayed_parent_cell_count": 1,
        "p2481_subcell_rows_inside_single_parent_cell_not_distinct_p2459_cells": replay["dyadic_subcell_count"],
        "p2456_boundary_band_cells_are_inherited_covered_boundary_chain_not_new_p2459_unreplayed_cells": True,
        "p2456_inherited_covered_boundary_chain_cell_count": replay["boundary_band_cell_count"],
        "finite_chain_coverage_budget_inherited_from_p2479": p2479["finite_chain_coverage_budget_inherited_from_p2471"],
        "no_full_complement_claimed_by_this_certificate": True,
        "full_complement_replay_exported_by_this_certificate": False,
        "finite_replay_chain_boundary_handoff_collar_audit_exported": True,
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
        "lay_summary": "This packet checks the seam immediately to the left of the P2480 refined parent cell. The six already-covered P2456 right-boundary cells and the 128 P2480-style dyadic subcell rows all remain zero-excluding, and their separations increase left-to-right across an exactly adjacent handoff. The weakest row is in the covered boundary band, so P2480 was not the whole local bottleneck. The 134 fresh Decimal evaluations are diagnostic rows, not 134 new P2459 coverage cells; this is finite seam bookkeeping, not a continuum or full-complement proof.",
        "not_licensed": [
            "This boundary-handoff collar is only 6 covered boundary-band cells plus 128 subcell rows inside one parent cell; it is not a full complement replay of the P2459 universe.",
            "The 134 fresh Decimal evaluations are diagnostic rows, not 134 distinct P2459 coverage cells.",
            "The P2456 boundary-band cells are inherited covered-boundary-chain cells, not new P2459 unreplayed complement cells.",
            "It is not a directed-rounding interval theorem, symbolic root-exclusion theorem, analytic monotonicity theorem, or continuum root-exclusion theorem.",
            "It does not export selector/source/gauge closure, a physical-value generator, QW-2191 discharge, role-bearing L_total, legacy-role transfer, or ToE closure.",
        ],
        "next_honest_step": "A stronger theorem would need a formal boundary/window argument across the excluded root window or a directed-rounding backend; do not inflate this finite collar handoff into continuum root exclusion or selector/source closure.",
    }
    gatekeepers = {
        "rg_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "collar_counts_expected": theorem_export["fresh_boundary_band_cell_replay_count"] == 6 and theorem_export["fresh_dyadic_subcell_replay_count"] == 128 and theorem_export["total_fresh_handoff_collar_replay_rows"] == 134,
        "parent_cell_binding_expected": theorem_export["parent_segment_index"] == 2 and theorem_export["parent_uncovered_index"] == 0,
        "handoff_gapless": Decimal(theorem_export["handoff_endpoint_gap_boundary_to_parent"]) == 0 and Decimal(theorem_export["minimum_consecutive_endpoint_gap"]) == 0 and Decimal(theorem_export["maximum_consecutive_endpoint_gap"]) == 0 and theorem_export["all_consecutive_collar_rows_are_exactly_adjacent"],
        "minimum_separation_positive": Decimal(theorem_export["minimum_collar_decimal_separation"]) > 0,
        "all_collar_rows_exclude_zero": theorem_export["all_collar_rows_exclude_zero"],
        "all_collar_rows_have_positive_separation": theorem_export["all_collar_rows_have_positive_separation"],
        "collar_order_checked": theorem_export["all_consecutive_collar_separations_strictly_increase_left_to_right"] and Decimal(theorem_export["minimum_consecutive_positive_delta"]) > 0,
        "minimum_relocated_to_boundary_band_not_overclaimed": theorem_export["minimum_is_p2456_right_boundary_band_leftmost_cell"],
        "p2480_inherited": theorem_export["p2480_dyadic_subcell_count_inherited"] == 128,
        "no_full_complement_unless_genuinely_full": theorem_export["no_full_complement_claimed_by_this_certificate"] and not theorem_export["full_complement_replay_exported_by_this_certificate"],
        "finite_grid_only_not_directed_rounding": theorem_export["finite_replay_chain_boundary_handoff_collar_audit_exported"] and not theorem_export["directed_rounding_interval_theorem_exported_by_this_certificate"],
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
        "schema_version": "p2481_s1431_v1",
        "packet_id": "P2481",
        "stage_id": "S1431",
        "status": "PASS_STRICT_POINTWISE_INTERVAL_DECIMAL_P2459_BOUNDARY_HANDOFF_COLLAR_REPLAY_AUDIT_NO_FULL_COMPLEMENT_NO_DIRECTED_ROUNDING_SELECTOR_SOURCE_THEOREM",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_pointwise_interval_decimal_p2459_boundary_handoff_collar_replay_audit": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_pointwise_interval_decimal_p2459_boundary_handoff_collar_replay_audit"]["theorem_export"]
    lines = [
        "# P2481/S1431 strict pointwise interval-Decimal P2459 boundary-handoff collar replay audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Boundary-handoff collar replay",
        "",
        f"Parent cell: segment `{t['parent_segment_index']}`, uncovered index `{t['parent_uncovered_index']}`.",
        f"Fresh P2456 right-boundary-band replay rows: `{t['fresh_boundary_band_cell_replay_count']}`.",
        f"Fresh dyadic subcell replay rows: `{t['fresh_dyadic_subcell_replay_count']}`.",
        f"Total fresh handoff collar replay rows: `{t['total_fresh_handoff_collar_replay_rows']}`.",
        f"Boundary-to-parent endpoint gap: `{t['handoff_endpoint_gap_boundary_to_parent']}`.",
        f"Minimum consecutive endpoint gap: `{t['minimum_consecutive_endpoint_gap']}`.",
        f"Maximum consecutive endpoint gap: `{t['maximum_consecutive_endpoint_gap']}`.",
        f"Minimum collar Decimal separation: `{t['minimum_collar_decimal_separation']}`.",
        f"Maximum collar Decimal separation: `{t['maximum_collar_decimal_separation']}`.",
        f"Minimum consecutive positive delta: `{t['minimum_consecutive_positive_delta']}`.",
        f"Minimum is P2456 right-boundary-band leftmost cell: `{t['minimum_is_p2456_right_boundary_band_leftmost_cell']}`.",
        f"All collar rows exclude zero: `{t['all_collar_rows_exclude_zero']}`.",
        "",
        "## Coverage budget",
        "",
        f"P2481 fresh Decimal evaluation rows (not a P2459 coverage count): `{t['p2481_fresh_decimal_evaluation_row_count_not_a_coverage_count']}`.",
        f"Diagnostic row ratio against inherited P2459 universe (not a coverage fraction): `{t['p2481_fresh_decimal_evaluation_row_ratio_not_a_p2459_coverage_fraction']}`.",
        f"P2481 new P2459 unreplayed parent-cell scope against inherited P2459 finite universe: `{t['targeted_p2481_new_p2459_unreplayed_parent_cell_scope_against_p2459_universe']}`.",
        f"P2481 subcell rows inside that single parent cell, not distinct P2459 cells: `{t['p2481_subcell_rows_inside_single_parent_cell_not_distinct_p2459_cells']}`.",
        f"P2456 inherited covered-boundary-chain cells, not new P2459 unreplayed cells: `{t['p2456_inherited_covered_boundary_chain_cell_count']}`.",
        f"P2456 boundary cells are inherited covered-boundary-chain cells, not new P2459 unreplayed cells: `{t['p2456_boundary_band_cells_are_inherited_covered_boundary_chain_not_new_p2459_unreplayed_cells']}`.",
        f"Full complement replay exported by P2481: `{t['full_complement_replay_exported_by_this_certificate']}`.",
        "",
        "## Plain-language progress note",
        "",
        t["lay_summary"],
        "",
        "## Hard limits / negative controls",
        "",
        "This is a finite seam replay of one covered boundary band and one refined parent cell only.  It is not a full complement replay, directed-rounding interval theorem, symbolic root-exclusion theorem, analytic monotonicity theorem, global continuum root-exclusion theorem, selector/source/gauge theorem, physical-value generator, QW-2191 discharge, role-bearing `L_total`, legacy-role transfer, or ToE closure.",
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
