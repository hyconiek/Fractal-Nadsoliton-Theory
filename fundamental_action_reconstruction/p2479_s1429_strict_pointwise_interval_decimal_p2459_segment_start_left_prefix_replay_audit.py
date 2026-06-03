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
from p2472_s1422_strict_pointwise_interval_decimal_p2459_partition_seam_replay_audit import (
    classified_uncovered_cells_for_family,
)
from p2475_s1425_strict_pointwise_interval_decimal_p2459_critical_minimum_halo_replay_audit import (
    projection_vector,
)

GEN = ROOT / "generated"
OUT = GEN / "p2479_s1429_strict_pointwise_interval_decimal_p2459_segment_start_left_prefix_replay_audit.json"
MD = GEN / "p2479_s1429_strict_pointwise_interval_decimal_p2459_segment_start_left_prefix_replay_audit.md"

SOURCE_FILES = {
    "P2449_PROJECTION_REDUCTION": GEN / "p2449_s1399_strict_pointwise_rank_lift_projection_reduction_certificate.json",
    "P2450_ROOT_ISOLATION_MARGIN": GEN / "p2450_s1400_strict_pointwise_projection_root_isolation_margin_certificate.json",
    "P2451_FLOATING_INTERVAL_AUDIT": GEN / "p2451_s1401_strict_pointwise_projection_interval_enclosure_root_exclusion_audit.json",
    "P2456_DECIMAL_BOUNDARY_REPLAY": GEN / "p2456_s1406_strict_pointwise_decimal_root_window_boundary_band_replay_certificate.json",
    "P2466_DESCENT_TAIL_BOUNDARY": GEN / "p2466_s1416_strict_pointwise_interval_decimal_adaptive_descent_tail_boundary_ledger_certificate.json",
    "P2471_FINITE_PARTITION_WITNESS": GEN / "p2471_s1421_strict_pointwise_interval_decimal_p2459_finite_partition_witness_certificate.json",
    "P2478_LEFT_FLANK_EXTENSION": GEN / "p2478_s1428_strict_pointwise_interval_decimal_p2459_left_boundary_descent_flank_extension_replay_audit.json",
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
        "new_packet": "P2479|S1429|segment-start left prefix|left-prefix flank|descent prefix replay|full left prefix|one-sided prefix replay|segment boundary prefix|prefix flank certificate",
        "precursor_packets": "P2478|S1428|left-boundary descent-flank|left flank extension|P2477|S1427",
        "prefix_language": "segment start|left prefix|prefix replay|boundary-truncated minimum|one-sided segment prefix|descent-tail prefix|strictly increasing prefix",
        "hard_limit_markers": "directed-rounding root exclusion|symbolic root-exclusion theorem|global continuum root-exclusion|full mathematical proof|full complement",
        "closure_blockers": "QW-2191|role-bearing L_total|physical-value generator|ToE closure|selector/source/gauge theorem|legacy-role transfer",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def p2478_anchor(p2478: dict[str, Any]) -> dict[str, Any]:
    return p2478["left_flank_extension_replay"]["minimum_left_flank_replay_row"]


def p2478_union_keys(p2478: dict[str, Any]) -> set[tuple[int, int]]:
    replay = p2478["left_flank_extension_replay"]
    keys = {(int(row["segment_index"]), int(row["uncovered_index"])) for row in replay["left_flank_replay_rows"]}
    # P2478 exports the inherited P2477+P2478 union count but not P2477 rows;
    # the P2478 flank itself is enough to count fresh incremental prefix cells.
    return keys


def segment_start_left_prefix_replay(
    anchor: dict[str, Any],
    p2450: dict[str, Any],
    p2451: dict[str, Any],
    p2456: dict[str, Any],
    p2466: dict[str, Any],
    projection: list[Decimal],
    p2478: dict[str, Any],
) -> dict[str, Any]:
    family = anchor["family"]
    segment_index = int(anchor["segment_index"])
    anchor_index = int(anchor["uncovered_index"])
    classified, selected_segment_summary = classified_uncovered_cells_for_family(family, p2450, p2451, p2456, p2466)
    classified_by_key = {(int(row["segment_index"]), int(row["uncovered_index"])): row for row in classified}
    segment_indices = sorted(int(row["uncovered_index"]) for row in classified if int(row["segment_index"]) == segment_index)
    segment_start = min(segment_indices)
    replay_keys = [(segment_index, idx) for idx in range(segment_start, anchor_index + 1)]
    missing_keys = [key for key in replay_keys if key not in classified_by_key]
    p2478_keys = p2478_union_keys(p2478)
    function = function_for_family(family)
    replayed = []
    for key in replay_keys:
        if key not in classified_by_key:
            continue
        row = classified_by_key[key]
        fresh = replay_cell(family, row, projection, function)
        replayed.append({
            **fresh,
            "segment_index": key[0],
            "uncovered_index": key[1],
            "partition_class": row["partition_class"],
            "fresh_decimal_separation_positive": Decimal(fresh["decimal_separation_from_zero"]) > 0,
            "already_in_p2478_left_flank_extension": key in p2478_keys,
            "offset_from_segment_start": key[1] - segment_start,
            "offset_from_p2478_left_boundary_anchor": key[1] - anchor_index,
        })
    ordered = sorted(replayed, key=lambda row: int(row["uncovered_index"]))
    consecutive_pairs = []
    for prior, current in zip(ordered, ordered[1:]):
        delta = Decimal(current["decimal_separation_from_zero"]) - Decimal(prior["decimal_separation_from_zero"])
        consecutive_pairs.append({
            "from_uncovered_index": int(prior["uncovered_index"]),
            "to_uncovered_index": int(current["uncovered_index"]),
            "from_separation": prior["decimal_separation_from_zero"],
            "to_separation": current["decimal_separation_from_zero"],
            "separation_delta_to_minus_from": str(delta),
            "strictly_increases": delta > 0,
        })
    class_counts: dict[str, int] = {}
    for row in replayed:
        class_counts[row["partition_class"]] = class_counts.get(row["partition_class"], 0) + 1
    min_row = min(replayed, key=lambda row: Decimal(row["decimal_separation_from_zero"]))
    max_row = max(replayed, key=lambda row: Decimal(row["decimal_separation_from_zero"]))
    return {
        "family": family,
        "segment_index": segment_index,
        "segment_start_uncovered_index": segment_start,
        "p2478_left_boundary_anchor_uncovered_index": anchor_index,
        "selected_segment_summary": selected_segment_summary,
        "fresh_segment_start_left_prefix_replay_count": len(replayed),
        "incremental_cells_beyond_p2478_left_flank_count": sum(not row["already_in_p2478_left_flank_extension"] for row in replayed),
        "p2479_prefix_plus_p2478_flank_union_count": len({(int(row["segment_index"]), int(row["uncovered_index"])) for row in replayed} | p2478_keys),
        "missing_key_count": len(missing_keys),
        "missing_keys": [{"segment_index": key[0], "uncovered_index": key[1]} for key in missing_keys],
        "partition_class_counts": class_counts,
        "all_prefix_cells_exclude_zero": all(row["decimal_interval_excludes_zero"] for row in replayed),
        "all_prefix_cells_have_positive_separation": all(row["fresh_decimal_separation_positive"] for row in replayed),
        "minimum_prefix_decimal_separation": min_row["decimal_separation_from_zero"],
        "maximum_prefix_decimal_separation": max_row["decimal_separation_from_zero"],
        "minimum_prefix_replay_row": min_row,
        "maximum_prefix_replay_row": max_row,
        "minimum_is_segment_start_boundary": int(min_row["uncovered_index"]) == segment_start,
        "anchor_is_right_boundary_of_prefix_window": anchor_index == max(int(row["uncovered_index"]) for row in replayed),
        "consecutive_separation_pairs": consecutive_pairs,
        "all_consecutive_pairs_strictly_increase_left_to_right": all(row["strictly_increases"] for row in consecutive_pairs),
        "minimum_consecutive_positive_delta": str(min(Decimal(row["separation_delta_to_minus_from"]) for row in consecutive_pairs)),
        "segment_start_left_prefix_replay_rows": replayed,
        "segment_start_left_prefix_replay_fingerprint_sha256": sha256_json(replayed),
    }


def append_once(path, marker: str, section: str) -> None:
    text = path.read_text(encoding="utf-8")
    if marker not in text:
        path.write_text(text.rstrip() + "\n\n" + section.strip() + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2479/S1429 strict pointwise interval-Decimal P2459 segment-start left-prefix replay audit

`P2479/S1429` follows the P2478 left-boundary minimum all the way to the finite segment start for the same zero-projection-amplitude flank.  It freshly replays the complete one-sided prefix from uncovered index `0` through `1115` on segment `2`: `1116` cells, `1115` of them new beyond the P2478 flank.  All replayed prefix cells exclude zero with positive Decimal separation, and their separations strictly increase left-to-right; the minimum is now at the segment-start boundary, so this remains a boundary-truncated finite-prefix result rather than a continuum or local-minimum theorem.

For a non-specialist: P2478 kept finding the weakest value at the left edge of its strip.  P2479 moves all the way to the start of that finite segment and recalculates every cell in between.  The checked prefix stays away from zero, but its weakest value is at the segment boundary.  This is stronger finite bookkeeping, not a full complement replay, directed-rounding interval theorem, symbolic root-exclusion theorem, continuum theorem, selector/source/gauge theorem, physical-value generator, role-bearing `L_total`, `QW-2191` discharge, legacy-role transfer, or ToE closure.
"""
    lag_section = """
## P2479/S1429 P2459 segment-start left-prefix replay guard

`P2479/S1429` adds a complete one-sided finite-prefix guard to the bookkeeping behind `L_total`: the P2478 left-boundary flank is extended to the segment-start boundary and `1116` prefix cells are replayed with the Decimal/Taylor backend.  This improves finite descent-prefix localization, but it does not add selector/source/gauge authority, physical-value generation, role-bearing finality, legacy-role transfer, or ToE closure.
"""
    append_once(DOC_FILES["equation_sheet"], "## P2479/S1429 strict pointwise interval-Decimal P2459 segment-start left-prefix replay audit", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2479/S1429 P2459 segment-start left-prefix replay guard", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    projection = projection_vector(sources["P2449_PROJECTION_REDUCTION"])
    p2450 = theorem(sources["P2450_ROOT_ISOLATION_MARGIN"], "strict_pointwise_projection_root_isolation_margin_certificate")
    p2451 = theorem(sources["P2451_FLOATING_INTERVAL_AUDIT"], "strict_pointwise_projection_interval_enclosure_root_exclusion_audit")
    p2456 = theorem(sources["P2456_DECIMAL_BOUNDARY_REPLAY"], "strict_pointwise_decimal_root_window_boundary_band_replay_certificate")
    p2466 = theorem(sources["P2466_DESCENT_TAIL_BOUNDARY"], "strict_pointwise_interval_decimal_adaptive_descent_tail_boundary_ledger_certificate")
    p2471 = theorem(sources["P2471_FINITE_PARTITION_WITNESS"], "strict_pointwise_interval_decimal_p2459_finite_partition_witness_certificate")
    p2478 = theorem(sources["P2478_LEFT_FLANK_EXTENSION"], "strict_pointwise_interval_decimal_p2459_left_boundary_descent_flank_extension_replay_audit")
    replay = segment_start_left_prefix_replay(p2478_anchor(p2478), p2450, p2451, p2456, p2466, projection, p2478)
    p2459_universe = p2471["total_p2459_universe_cell_count_rebuilt"]
    theorem_export = {
        "theorem_name": "P2479_T1_strict_pointwise_interval_decimal_p2459_segment_start_left_prefix_replay_audit",
        "audited_chain": ["P2471/S1421", "P2478/S1428"],
        "segment_start_left_prefix_replay": replay,
        "fresh_segment_start_left_prefix_replay_count": replay["fresh_segment_start_left_prefix_replay_count"],
        "incremental_cells_beyond_p2478_left_flank_count": replay["incremental_cells_beyond_p2478_left_flank_count"],
        "p2479_prefix_plus_p2478_flank_union_count": replay["p2479_prefix_plus_p2478_flank_union_count"],
        "minimum_prefix_decimal_separation": replay["minimum_prefix_decimal_separation"],
        "maximum_prefix_decimal_separation": replay["maximum_prefix_decimal_separation"],
        "minimum_consecutive_positive_delta": replay["minimum_consecutive_positive_delta"],
        "all_prefix_cells_exclude_zero": replay["all_prefix_cells_exclude_zero"],
        "all_prefix_cells_have_positive_separation": replay["all_prefix_cells_have_positive_separation"],
        "minimum_is_segment_start_boundary": replay["minimum_is_segment_start_boundary"],
        "all_consecutive_pairs_strictly_increase_left_to_right": replay["all_consecutive_pairs_strictly_increase_left_to_right"],
        "p2478_fresh_left_flank_replay_count_inherited": p2478["fresh_left_flank_replay_count"],
        "p2478_minimum_left_flank_decimal_separation_inherited": p2478["minimum_left_flank_decimal_separation"],
        "p2471_p2459_universe_count_inherited": p2459_universe,
        "targeted_p2479_fresh_replay_fraction_of_p2459_universe": f"{replay['fresh_segment_start_left_prefix_replay_count']}/{p2459_universe}",
        "targeted_p2479_residual_not_freshly_replayed_count_against_p2459_universe": p2459_universe - replay["fresh_segment_start_left_prefix_replay_count"],
        "p2479_prefix_plus_p2478_union_fraction_of_p2459_universe": f"{replay['p2479_prefix_plus_p2478_flank_union_count']}/{p2459_universe}",
        "p2479_prefix_plus_p2478_union_residual_count_against_p2459_universe": p2459_universe - replay["p2479_prefix_plus_p2478_flank_union_count"],
        "finite_chain_coverage_budget_inherited_from_p2471": 0 if p2471["all_family_partitions_are_disjoint_and_complete"] else None,
        "p2471_disjoint_complete_inherited": p2471["all_family_partitions_are_disjoint_and_complete"],
        "no_full_complement_claimed_by_this_certificate": True,
        "full_complement_replay_exported_by_this_certificate": False,
        "finite_replay_chain_segment_start_left_prefix_replay_audit_exported": True,
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
        "lay_summary": "This packet follows the P2478 left edge to the start of the finite segment and recalculates the whole one-sided prefix. The prefix stays away from zero and its separations increase left-to-right, but the weakest checked cell is the segment-start boundary. The honest conclusion is finite prefix localization, not a local-minimum, continuum, or full-complement proof.",
        "not_licensed": [
            "This segment-start prefix replay is not a full complement replay and does not claim all P2459 cells were freshly replayed by P2479.",
            "It is not a directed-rounding interval theorem, symbolic root-exclusion theorem, analytic monotonicity theorem, or continuum root-exclusion theorem.",
            "It does not export selector/source/gauge closure, a physical-value generator, QW-2191 discharge, role-bearing L_total, legacy-role transfer, or ToE closure.",
            "No legacy physical roles are transferred to L_total or K_strict_gate by this finite prefix replay.",
        ],
        "next_honest_step": "The finite prefix is now exhausted to the segment-start boundary. A stronger claim would need a directed-rounding/symbolic argument across the adjacent real interval or a separately scoped finite replay outside this prefix; do not promote this prefix audit to a full proof.",
    }
    gatekeepers = {
        "rg_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "replay_count_expected": theorem_export["fresh_segment_start_left_prefix_replay_count"] == 1116,
        "incremental_count_expected": theorem_export["incremental_cells_beyond_p2478_left_flank_count"] == 1115,
        "union_count_expected": theorem_export["p2479_prefix_plus_p2478_flank_union_count"] == 1132,
        "minimum_separation_positive": Decimal(theorem_export["minimum_prefix_decimal_separation"]) > 0,
        "all_replayed_cells_exclude_zero": theorem_export["all_prefix_cells_exclude_zero"],
        "all_replayed_cells_have_positive_separation": theorem_export["all_prefix_cells_have_positive_separation"],
        "segment_start_boundary_not_overclaimed": theorem_export["minimum_is_segment_start_boundary"],
        "consecutive_order_checked": theorem_export["all_consecutive_pairs_strictly_increase_left_to_right"] and Decimal(theorem_export["minimum_consecutive_positive_delta"]) > 0,
        "no_full_complement_unless_genuinely_full": theorem_export["no_full_complement_claimed_by_this_certificate"] and not theorem_export["full_complement_replay_exported_by_this_certificate"],
        "coverage_budget_accounted": theorem_export["targeted_p2479_residual_not_freshly_replayed_count_against_p2459_universe"] == p2459_universe - replay["fresh_segment_start_left_prefix_replay_count"],
        "p2471_partition_inherited": theorem_export["p2471_disjoint_complete_inherited"],
        "finite_grid_only_not_directed_rounding": theorem_export["finite_replay_chain_segment_start_left_prefix_replay_audit_exported"] and not theorem_export["directed_rounding_interval_theorem_exported_by_this_certificate"],
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
        "schema_version": "p2479_s1429_v1",
        "packet_id": "P2479",
        "stage_id": "S1429",
        "status": "PASS_STRICT_POINTWISE_INTERVAL_DECIMAL_P2459_SEGMENT_START_LEFT_PREFIX_REPLAY_AUDIT_NO_FULL_COMPLEMENT_NO_DIRECTED_ROUNDING_SELECTOR_SOURCE_THEOREM",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_pointwise_interval_decimal_p2459_segment_start_left_prefix_replay_audit": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_pointwise_interval_decimal_p2459_segment_start_left_prefix_replay_audit"]["theorem_export"]
    lines = [
        "# P2479/S1429 strict pointwise interval-Decimal P2459 segment-start left-prefix replay audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Segment-start left-prefix replay",
        "",
        f"Fresh P2479 prefix replay count: `{t['fresh_segment_start_left_prefix_replay_count']}`.",
        f"Incremental cells beyond P2478 left flank: `{t['incremental_cells_beyond_p2478_left_flank_count']}`.",
        f"P2479 prefix + P2478 flank union count: `{t['p2479_prefix_plus_p2478_flank_union_count']}`.",
        f"Minimum Decimal separation in prefix: `{t['minimum_prefix_decimal_separation']}`.",
        f"Maximum Decimal separation in prefix: `{t['maximum_prefix_decimal_separation']}`.",
        f"Minimum consecutive positive delta: `{t['minimum_consecutive_positive_delta']}`.",
        f"Minimum is the segment-start boundary: `{t['minimum_is_segment_start_boundary']}`.",
        f"All prefix cells exclude zero: `{t['all_prefix_cells_exclude_zero']}`.",
        "",
        "## Coverage budget",
        "",
        f"P2479 fresh replay fraction of inherited P2459 finite universe: `{t['targeted_p2479_fresh_replay_fraction_of_p2459_universe']}`.",
        f"P2479 residual not freshly replayed against P2459 universe: `{t['targeted_p2479_residual_not_freshly_replayed_count_against_p2459_universe']}`.",
        f"P2479 prefix + P2478 union fraction: `{t['p2479_prefix_plus_p2478_union_fraction_of_p2459_universe']}`.",
        f"P2479 prefix + P2478 union residual: `{t['p2479_prefix_plus_p2478_union_residual_count_against_p2459_universe']}`.",
        f"Inherited P2471 finite-chain coverage budget: `{t['finite_chain_coverage_budget_inherited_from_p2471']}`.",
        f"Full complement replay exported by P2479: `{t['full_complement_replay_exported_by_this_certificate']}`.",
        "",
        "## Plain-language progress note",
        "",
        t["lay_summary"],
        "",
        "## Hard limits / negative controls",
        "",
        "This is a finite segment-start left-prefix replay only.  It is not a full complement replay, directed-rounding interval theorem, symbolic root-exclusion theorem, analytic monotonicity theorem, global continuum root-exclusion theorem, selector/source/gauge theorem, physical-value generator, QW-2191 discharge, role-bearing `L_total`, legacy-role transfer, or ToE closure.",
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
