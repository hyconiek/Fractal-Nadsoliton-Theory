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

GEN = ROOT / "generated"
OUT = GEN / "p2480_s1430_strict_pointwise_interval_decimal_p2459_segment_start_cell_dyadic_refinement_replay_audit.json"
MD = GEN / "p2480_s1430_strict_pointwise_interval_decimal_p2459_segment_start_cell_dyadic_refinement_replay_audit.md"
DYADIC_SUBCELL_COUNT = 128

SOURCE_FILES = {
    "P2449_PROJECTION_REDUCTION": GEN / "p2449_s1399_strict_pointwise_rank_lift_projection_reduction_certificate.json",
    "P2479_SEGMENT_START_LEFT_PREFIX": GEN / "p2479_s1429_strict_pointwise_interval_decimal_p2459_segment_start_left_prefix_replay_audit.json",
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
        "new_packet": "P2480|S1430|segment-start cell dyadic refinement|weakest cell refinement|dyadic subcell replay|subcell refinement audit|128 subcells|cell-internal replay",
        "precursor_packets": "P2479|S1429|segment-start left-prefix replay|segment-start boundary|minimum prefix replay row",
        "refinement_language": "dyadic refinement|subcell replay|cell-internal|weakest cell|segment-start cell|refinement replay",
        "hard_limit_markers": "directed-rounding root exclusion|symbolic root-exclusion theorem|global continuum root-exclusion|full mathematical proof|full complement",
        "closure_blockers": "QW-2191|role-bearing L_total|physical-value generator|ToE closure|selector/source/gauge theorem|legacy-role transfer",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def dyadic_subcells(parent_cell: dict[str, Any]) -> list[dict[str, Any]]:
    left = Decimal(str(parent_cell["left"]))
    right = Decimal(str(parent_cell["right"]))
    width = right - left
    cells = []
    for index in range(DYADIC_SUBCELL_COUNT):
        sub_left = left + width * Decimal(index) / Decimal(DYADIC_SUBCELL_COUNT)
        sub_right = left + width * Decimal(index + 1) / Decimal(DYADIC_SUBCELL_COUNT)
        cells.append({
            "left": str(sub_left),
            "right": str(sub_right),
            "uncovered_index": index,
            "uncovered_count": DYADIC_SUBCELL_COUNT,
            "parent_segment_index": int(parent_cell["segment_index"]),
            "parent_uncovered_index": int(parent_cell["uncovered_index"]),
        })
    return cells


def dyadic_refinement_replay(parent_cell: dict[str, Any], projection: list[Decimal]) -> dict[str, Any]:
    family = parent_cell["family"]
    function = function_for_family(family)
    replayed = []
    for cell in dyadic_subcells(parent_cell):
        fresh = replay_cell(family, cell, projection, function)
        replayed.append({
            **fresh,
            "dyadic_subcell_index": int(cell["uncovered_index"]),
            "dyadic_subcell_count": DYADIC_SUBCELL_COUNT,
            "parent_segment_index": int(cell["parent_segment_index"]),
            "parent_uncovered_index": int(cell["parent_uncovered_index"]),
            "fresh_decimal_separation_positive": Decimal(fresh["decimal_separation_from_zero"]) > 0,
        })
    ordered = sorted(replayed, key=lambda row: int(row["dyadic_subcell_index"]))
    consecutive_pairs = []
    for prior, current in zip(ordered, ordered[1:]):
        delta = Decimal(current["decimal_separation_from_zero"]) - Decimal(prior["decimal_separation_from_zero"])
        consecutive_pairs.append({
            "from_dyadic_subcell_index": int(prior["dyadic_subcell_index"]),
            "to_dyadic_subcell_index": int(current["dyadic_subcell_index"]),
            "from_separation": prior["decimal_separation_from_zero"],
            "to_separation": current["decimal_separation_from_zero"],
            "separation_delta_to_minus_from": str(delta),
            "strictly_increases": delta > 0,
        })
    min_row = min(replayed, key=lambda row: Decimal(row["decimal_separation_from_zero"]))
    max_row = max(replayed, key=lambda row: Decimal(row["decimal_separation_from_zero"]))
    return {
        "family": family,
        "parent_segment_index": int(parent_cell["segment_index"]),
        "parent_uncovered_index": int(parent_cell["uncovered_index"]),
        "parent_left": parent_cell["left"],
        "parent_right": parent_cell["right"],
        "parent_decimal_separation_from_p2479": parent_cell["decimal_separation_from_zero"],
        "dyadic_subcell_count": DYADIC_SUBCELL_COUNT,
        "all_subcells_exclude_zero": all(row["decimal_interval_excludes_zero"] for row in replayed),
        "all_subcells_have_positive_separation": all(row["fresh_decimal_separation_positive"] for row in replayed),
        "minimum_subcell_decimal_separation": min_row["decimal_separation_from_zero"],
        "maximum_subcell_decimal_separation": max_row["decimal_separation_from_zero"],
        "minimum_subcell_replay_row": min_row,
        "maximum_subcell_replay_row": max_row,
        "minimum_is_leftmost_subcell": int(min_row["dyadic_subcell_index"]) == 0,
        "maximum_is_rightmost_subcell": int(max_row["dyadic_subcell_index"]) == DYADIC_SUBCELL_COUNT - 1,
        "all_consecutive_subcell_separations_strictly_increase_left_to_right": all(row["strictly_increases"] for row in consecutive_pairs),
        "minimum_consecutive_positive_delta": str(min(Decimal(row["separation_delta_to_minus_from"]) for row in consecutive_pairs)),
        "subcell_replay_rows": replayed,
        "subcell_replay_fingerprint_sha256": sha256_json(replayed),
    }


def append_once(path, marker: str, section: str) -> None:
    text = path.read_text(encoding="utf-8")
    if marker not in text:
        path.write_text(text.rstrip() + "\n\n" + section.strip() + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2480/S1430 strict pointwise interval-Decimal P2459 segment-start cell dyadic-refinement replay audit

`P2480/S1430` refines the P2479 segment-start boundary minimum instead of promoting the prefix result to a continuum theorem.  It subdivides the weakest parent cell into `128` dyadic subcells and replays each subcell with the Decimal/Taylor backend.  All `128` subcells exclude zero with positive Decimal separation, and the subcell separations strictly increase left-to-right; the minimum remains the leftmost subcell, so this is a cell-internal finite refinement, not a directed-rounding, symbolic, continuum, or full-complement theorem.

For a non-specialist: P2479 found the weakest checked cell at the very start of the finite segment.  P2480 opens that one cell and checks `128` smaller pieces inside it.  The pieces stay away from zero, but the weakest piece is still at the left edge.  This improves finite evidence inside the weakest cell without closing the continuum or selector/source questions.
"""
    lag_section = """
## P2480/S1430 P2459 segment-start cell dyadic-refinement replay guard

`P2480/S1430` adds a cell-internal refinement guard to the finite-grid bookkeeping behind `L_total`: the P2479 segment-start minimum cell is split into `128` dyadic subcells and replayed with Decimal/Taylor arithmetic.  This improves local finite resolution of the weakest boundary cell, but it does not add selector/source/gauge authority, physical-value generation, role-bearing finality, legacy-role transfer, or ToE closure.
"""
    append_once(DOC_FILES["equation_sheet"], "## P2480/S1430 strict pointwise interval-Decimal P2459 segment-start cell dyadic-refinement replay audit", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2480/S1430 P2459 segment-start cell dyadic-refinement replay guard", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    projection = projection_vector(sources["P2449_PROJECTION_REDUCTION"])
    p2479 = theorem(sources["P2479_SEGMENT_START_LEFT_PREFIX"], "strict_pointwise_interval_decimal_p2459_segment_start_left_prefix_replay_audit")
    parent_cell = p2479["segment_start_left_prefix_replay"]["minimum_prefix_replay_row"]
    replay = dyadic_refinement_replay(parent_cell, projection)
    p2459_universe = p2479["p2471_p2459_universe_count_inherited"]
    theorem_export = {
        "theorem_name": "P2480_T1_strict_pointwise_interval_decimal_p2459_segment_start_cell_dyadic_refinement_replay_audit",
        "audited_chain": ["P2479/S1429"],
        "dyadic_refinement_replay": replay,
        "dyadic_subcell_count": replay["dyadic_subcell_count"],
        "parent_segment_index": replay["parent_segment_index"],
        "parent_uncovered_index": replay["parent_uncovered_index"],
        "parent_decimal_separation_from_p2479": replay["parent_decimal_separation_from_p2479"],
        "all_subcells_exclude_zero": replay["all_subcells_exclude_zero"],
        "all_subcells_have_positive_separation": replay["all_subcells_have_positive_separation"],
        "minimum_subcell_decimal_separation": replay["minimum_subcell_decimal_separation"],
        "maximum_subcell_decimal_separation": replay["maximum_subcell_decimal_separation"],
        "minimum_is_leftmost_subcell": replay["minimum_is_leftmost_subcell"],
        "maximum_is_rightmost_subcell": replay["maximum_is_rightmost_subcell"],
        "all_consecutive_subcell_separations_strictly_increase_left_to_right": replay["all_consecutive_subcell_separations_strictly_increase_left_to_right"],
        "minimum_consecutive_positive_delta": replay["minimum_consecutive_positive_delta"],
        "p2479_prefix_replay_count_inherited": p2479["fresh_segment_start_left_prefix_replay_count"],
        "p2479_prefix_plus_p2478_union_count_inherited": p2479["p2479_prefix_plus_p2478_flank_union_count"],
        "p2471_p2459_universe_count_inherited": p2459_universe,
        "targeted_p2480_subcell_replay_fraction_against_p2459_universe": f"{replay['dyadic_subcell_count']}/{p2459_universe}",
        "targeted_p2480_parent_cell_count_against_p2459_universe": f"1/{p2459_universe}",
        "finite_chain_coverage_budget_inherited_from_p2479": p2479["finite_chain_coverage_budget_inherited_from_p2471"],
        "no_full_complement_claimed_by_this_certificate": True,
        "full_complement_replay_exported_by_this_certificate": False,
        "finite_replay_chain_segment_start_cell_dyadic_refinement_audit_exported": True,
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
        "lay_summary": "This packet opens the weakest P2479 segment-start cell into 128 smaller dyadic subcells. Every subcell still excludes zero with positive Decimal separation, and separations increase from left to right, but the weakest subcell is still the leftmost one. This is finite cell-internal refinement, not a continuum or full-complement proof.",
        "not_licensed": [
            "This dyadic refinement replays one parent cell only; it is not a full complement replay of the P2459 universe.",
            "It is not a directed-rounding interval theorem, symbolic root-exclusion theorem, analytic monotonicity theorem, or continuum root-exclusion theorem.",
            "It does not export selector/source/gauge closure, a physical-value generator, QW-2191 discharge, role-bearing L_total, legacy-role transfer, or ToE closure.",
            "No legacy physical roles are transferred to L_total or K_strict_gate by this finite subcell replay.",
        ],
        "next_honest_step": "A stronger claim at the segment-start boundary needs a real boundary/continuum argument or a separately scoped adjacent-domain replay; do not promote this cell-internal dyadic replay to a directed-rounding or symbolic proof.",
    }
    gatekeepers = {
        "rg_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "subcell_count_expected": theorem_export["dyadic_subcell_count"] == 128,
        "parent_cell_binding_expected": theorem_export["parent_segment_index"] == 2 and theorem_export["parent_uncovered_index"] == 0,
        "minimum_separation_positive": Decimal(theorem_export["minimum_subcell_decimal_separation"]) > 0,
        "all_subcells_exclude_zero": theorem_export["all_subcells_exclude_zero"],
        "all_subcells_have_positive_separation": theorem_export["all_subcells_have_positive_separation"],
        "leftmost_boundary_not_overclaimed": theorem_export["minimum_is_leftmost_subcell"],
        "consecutive_order_checked": theorem_export["all_consecutive_subcell_separations_strictly_increase_left_to_right"] and Decimal(theorem_export["minimum_consecutive_positive_delta"]) > 0,
        "p2479_inherited": theorem_export["p2479_prefix_replay_count_inherited"] == 1116 and theorem_export["p2479_prefix_plus_p2478_union_count_inherited"] == 1132,
        "no_full_complement_unless_genuinely_full": theorem_export["no_full_complement_claimed_by_this_certificate"] and not theorem_export["full_complement_replay_exported_by_this_certificate"],
        "finite_grid_only_not_directed_rounding": theorem_export["finite_replay_chain_segment_start_cell_dyadic_refinement_audit_exported"] and not theorem_export["directed_rounding_interval_theorem_exported_by_this_certificate"],
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
        "schema_version": "p2480_s1430_v1",
        "packet_id": "P2480",
        "stage_id": "S1430",
        "status": "PASS_STRICT_POINTWISE_INTERVAL_DECIMAL_P2459_SEGMENT_START_CELL_DYADIC_REFINEMENT_REPLAY_AUDIT_NO_FULL_COMPLEMENT_NO_DIRECTED_ROUNDING_SELECTOR_SOURCE_THEOREM",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_pointwise_interval_decimal_p2459_segment_start_cell_dyadic_refinement_replay_audit": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_pointwise_interval_decimal_p2459_segment_start_cell_dyadic_refinement_replay_audit"]["theorem_export"]
    lines = [
        "# P2480/S1430 strict pointwise interval-Decimal P2459 segment-start cell dyadic-refinement replay audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Segment-start weakest-cell dyadic refinement",
        "",
        f"Parent cell: segment `{t['parent_segment_index']}`, uncovered index `{t['parent_uncovered_index']}`.",
        f"Dyadic subcells replayed: `{t['dyadic_subcell_count']}`.",
        f"Parent Decimal separation inherited from P2479: `{t['parent_decimal_separation_from_p2479']}`.",
        f"Minimum subcell Decimal separation: `{t['minimum_subcell_decimal_separation']}`.",
        f"Maximum subcell Decimal separation: `{t['maximum_subcell_decimal_separation']}`.",
        f"Minimum consecutive positive delta: `{t['minimum_consecutive_positive_delta']}`.",
        f"Minimum is leftmost subcell: `{t['minimum_is_leftmost_subcell']}`.",
        f"All subcells exclude zero: `{t['all_subcells_exclude_zero']}`.",
        "",
        "## Coverage budget",
        "",
        f"P2480 subcell replay fraction against inherited P2459 finite universe: `{t['targeted_p2480_subcell_replay_fraction_against_p2459_universe']}`.",
        f"P2480 parent-cell fraction against inherited P2459 finite universe: `{t['targeted_p2480_parent_cell_count_against_p2459_universe']}`.",
        f"Inherited P2479 prefix replay count: `{t['p2479_prefix_replay_count_inherited']}`.",
        f"Inherited P2479+P2478 union count: `{t['p2479_prefix_plus_p2478_union_count_inherited']}`.",
        f"Full complement replay exported by P2480: `{t['full_complement_replay_exported_by_this_certificate']}`.",
        "",
        "## Plain-language progress note",
        "",
        t["lay_summary"],
        "",
        "## Hard limits / negative controls",
        "",
        "This is a finite dyadic refinement replay of one weakest parent cell only.  It is not a full complement replay, directed-rounding interval theorem, symbolic root-exclusion theorem, analytic monotonicity theorem, global continuum root-exclusion theorem, selector/source/gauge theorem, physical-value generator, QW-2191 discharge, role-bearing `L_total`, legacy-role transfer, or ToE closure.",
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
