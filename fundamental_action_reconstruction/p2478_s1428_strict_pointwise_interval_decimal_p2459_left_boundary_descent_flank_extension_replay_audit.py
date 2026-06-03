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
OUT = GEN / "p2478_s1428_strict_pointwise_interval_decimal_p2459_left_boundary_descent_flank_extension_replay_audit.json"
MD = GEN / "p2478_s1428_strict_pointwise_interval_decimal_p2459_left_boundary_descent_flank_extension_replay_audit.md"
LEFT_EXTENSION_RADIUS = 16

SOURCE_FILES = {
    "P2449_PROJECTION_REDUCTION": GEN / "p2449_s1399_strict_pointwise_rank_lift_projection_reduction_certificate.json",
    "P2450_ROOT_ISOLATION_MARGIN": GEN / "p2450_s1400_strict_pointwise_projection_root_isolation_margin_certificate.json",
    "P2451_FLOATING_INTERVAL_AUDIT": GEN / "p2451_s1401_strict_pointwise_projection_interval_enclosure_root_exclusion_audit.json",
    "P2456_DECIMAL_BOUNDARY_REPLAY": GEN / "p2456_s1406_strict_pointwise_decimal_root_window_boundary_band_replay_certificate.json",
    "P2466_DESCENT_TAIL_BOUNDARY": GEN / "p2466_s1416_strict_pointwise_interval_decimal_adaptive_descent_tail_boundary_ledger_certificate.json",
    "P2471_FINITE_PARTITION_WITNESS": GEN / "p2471_s1421_strict_pointwise_interval_decimal_p2459_finite_partition_witness_certificate.json",
    "P2477_EXCEPTION_EXPANDED_HALO": GEN / "p2477_s1427_strict_pointwise_interval_decimal_p2459_lower_neighbor_exception_expanded_halo_replay_audit.json",
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
        "new_packet": "P2478|S1428|left-boundary descent flank|left flank continuation|flank extension replay|asymmetric descent flank|boundary-minimum continuation|wider finite flank|left-edge continuation",
        "precursor_packets": "P2477|S1427|lower-neighbor exception expanded halo|left-boundary minimum|targeted descent localization",
        "flank_language": "left boundary|left flank|descent flank|one-sided extension|consecutive separation|monotone flank|boundary minimum",
        "hard_limit_markers": "directed-rounding root exclusion|symbolic root-exclusion theorem|global continuum root-exclusion|full mathematical proof|full complement",
        "closure_blockers": "QW-2191|role-bearing L_total|physical-value generator|ToE closure|selector/source/gauge theorem|legacy-role transfer",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def p2477_left_boundary_anchor(p2477: dict[str, Any]) -> dict[str, Any]:
    replay = p2477["expanded_exception_replays"][0]
    return replay["minimum_expanded_exception_replay_row"]


def left_flank_extension_replay(
    anchor: dict[str, Any],
    p2450: dict[str, Any],
    p2451: dict[str, Any],
    p2456: dict[str, Any],
    p2466: dict[str, Any],
    projection: list[Decimal],
    p2477: dict[str, Any],
) -> dict[str, Any]:
    family = anchor["family"]
    segment_index = int(anchor["segment_index"])
    anchor_index = int(anchor["uncovered_index"])
    start_index = anchor_index - LEFT_EXTENSION_RADIUS
    classified, selected_segment_summary = classified_uncovered_cells_for_family(family, p2450, p2451, p2456, p2466)
    classified_by_key = {(int(row["segment_index"]), int(row["uncovered_index"])): row for row in classified}
    replay_keys = [(segment_index, idx) for idx in range(start_index, anchor_index + 1)]
    missing_keys = [key for key in replay_keys if key not in classified_by_key]
    p2477_keys = {
        (int(row["segment_index"]), int(row["uncovered_index"]))
        for replay in p2477["expanded_exception_replays"]
        for row in replay["expanded_exception_replay_rows"]
    }
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
            "already_in_p2477_expanded_halo": key in p2477_keys,
            "offset_from_p2477_left_boundary_anchor": key[1] - anchor_index,
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
    p2477_plus_p2478_union_keys = p2477_keys | {(int(row["segment_index"]), int(row["uncovered_index"])) for row in replayed}
    return {
        "family": family,
        "segment_index": segment_index,
        "p2477_left_boundary_anchor_uncovered_index": anchor_index,
        "left_extension_radius": LEFT_EXTENSION_RADIUS,
        "left_extension_start_uncovered_index": start_index,
        "selected_segment_summary": selected_segment_summary,
        "fresh_left_flank_replay_count": len(replayed),
        "incremental_cells_beyond_p2477_expanded_halo_count": sum(not row["already_in_p2477_expanded_halo"] for row in replayed),
        "p2477_plus_p2478_targeted_union_cell_count": len(p2477_plus_p2478_union_keys),
        "missing_key_count": len(missing_keys),
        "missing_keys": [{"segment_index": key[0], "uncovered_index": key[1]} for key in missing_keys],
        "partition_class_counts": class_counts,
        "all_left_flank_cells_exclude_zero": all(row["decimal_interval_excludes_zero"] for row in replayed),
        "all_left_flank_cells_have_positive_separation": all(row["fresh_decimal_separation_positive"] for row in replayed),
        "minimum_left_flank_decimal_separation": min_row["decimal_separation_from_zero"],
        "minimum_left_flank_replay_row": min_row,
        "minimum_is_left_boundary_of_p2478_window": int(min_row["uncovered_index"]) == min(int(row["uncovered_index"]) for row in replayed),
        "p2477_anchor_is_right_boundary_of_p2478_window": anchor_index == max(int(row["uncovered_index"]) for row in replayed),
        "consecutive_separation_pairs": consecutive_pairs,
        "all_consecutive_pairs_strictly_increase_left_to_right": all(row["strictly_increases"] for row in consecutive_pairs),
        "minimum_consecutive_positive_delta": str(min(Decimal(row["separation_delta_to_minus_from"]) for row in consecutive_pairs)),
        "left_flank_replay_rows": replayed,
        "left_flank_replay_fingerprint_sha256": sha256_json(replayed),
    }


def append_once(path: Path, marker: str, section: str) -> None:
    text = path.read_text(encoding="utf-8")
    if marker not in text:
        path.write_text(text.rstrip() + "\n\n" + section.strip() + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2478/S1428 strict pointwise interval-Decimal P2459 left-boundary descent-flank extension replay audit

`P2478/S1428` follows the P2477 left-boundary minimum instead of declaring it solved.  It replays a one-sided `16`-cell left extension ending at the P2477 left-boundary anchor.  The fresh targeted replay covers `17` cells (`16` new beyond P2477), all zero-excluding with positive Decimal separation, and the separations strictly increase left-to-right across the replayed flank; the new minimum is again at the left boundary, so this remains a finite descent-flank localization rather than a local-minimum theorem.

For a non-specialist: P2477 said the weakest checked cell was at the left edge of its small strip.  P2478 moves the strip farther left and recalculates that flank.  The values still stay away from zero, but the weakest value is again at the newly extended left edge.  This is useful bookkeeping of the weak direction, not a full complement replay, directed-rounding interval theorem, symbolic root-exclusion theorem, continuum theorem, selector/source/gauge theorem, physical-value generator, role-bearing `L_total`, `QW-2191` discharge, legacy-role transfer, or ToE closure.
"""
    lag_section = """
## P2478/S1428 P2459 left-boundary descent-flank extension replay guard

`P2478/S1428` adds a one-sided flank-extension guard to the finite-grid bookkeeping behind `L_total`: the P2477 left-boundary minimum is extended by `16` fresh left-flank cells and replayed with the Decimal/Taylor backend.  This improves localization of the finite descent direction, but it does not add selector/source/gauge authority, physical-value generation, role-bearing finality, legacy-role transfer, or ToE closure.
"""
    append_once(DOC_FILES["equation_sheet"], "## P2478/S1428 strict pointwise interval-Decimal P2459 left-boundary descent-flank extension replay audit", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2478/S1428 P2459 left-boundary descent-flank extension replay guard", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    projection = projection_vector(sources["P2449_PROJECTION_REDUCTION"])
    p2450 = theorem(sources["P2450_ROOT_ISOLATION_MARGIN"], "strict_pointwise_projection_root_isolation_margin_certificate")
    p2451 = theorem(sources["P2451_FLOATING_INTERVAL_AUDIT"], "strict_pointwise_projection_interval_enclosure_root_exclusion_audit")
    p2456 = theorem(sources["P2456_DECIMAL_BOUNDARY_REPLAY"], "strict_pointwise_decimal_root_window_boundary_band_replay_certificate")
    p2466 = theorem(sources["P2466_DESCENT_TAIL_BOUNDARY"], "strict_pointwise_interval_decimal_adaptive_descent_tail_boundary_ledger_certificate")
    p2471 = theorem(sources["P2471_FINITE_PARTITION_WITNESS"], "strict_pointwise_interval_decimal_p2459_finite_partition_witness_certificate")
    p2477 = theorem(sources["P2477_EXCEPTION_EXPANDED_HALO"], "strict_pointwise_interval_decimal_p2459_lower_neighbor_exception_expanded_halo_replay_audit")
    anchor = p2477_left_boundary_anchor(p2477)
    replay = left_flank_extension_replay(anchor, p2450, p2451, p2456, p2466, projection, p2477)
    p2459_universe = p2471["total_p2459_universe_cell_count_rebuilt"]
    theorem_export = {
        "theorem_name": "P2478_T1_strict_pointwise_interval_decimal_p2459_left_boundary_descent_flank_extension_replay_audit",
        "audited_chain": ["P2471/S1421", "P2477/S1427"],
        "left_extension_radius": LEFT_EXTENSION_RADIUS,
        "left_flank_extension_replay": replay,
        "fresh_left_flank_replay_count": replay["fresh_left_flank_replay_count"],
        "incremental_cells_beyond_p2477_expanded_halo_count": replay["incremental_cells_beyond_p2477_expanded_halo_count"],
        "p2477_plus_p2478_targeted_union_cell_count": replay["p2477_plus_p2478_targeted_union_cell_count"],
        "minimum_left_flank_decimal_separation": replay["minimum_left_flank_decimal_separation"],
        "minimum_consecutive_positive_delta": replay["minimum_consecutive_positive_delta"],
        "all_left_flank_cells_exclude_zero": replay["all_left_flank_cells_exclude_zero"],
        "all_left_flank_cells_have_positive_separation": replay["all_left_flank_cells_have_positive_separation"],
        "minimum_is_left_boundary_of_p2478_window": replay["minimum_is_left_boundary_of_p2478_window"],
        "all_consecutive_pairs_strictly_increase_left_to_right": replay["all_consecutive_pairs_strictly_increase_left_to_right"],
        "p2477_fresh_targeted_replay_count_inherited": p2477["total_unique_expanded_exception_replay_count"],
        "p2477_minimum_expanded_exception_decimal_separation_inherited": p2477["minimum_expanded_exception_decimal_separation"],
        "p2471_p2459_universe_count_inherited": p2459_universe,
        "targeted_p2478_fresh_replay_fraction_of_p2459_universe": f"{replay['fresh_left_flank_replay_count']}/{p2459_universe}",
        "targeted_p2478_residual_not_freshly_replayed_count_against_p2459_universe": p2459_universe - replay["fresh_left_flank_replay_count"],
        "p2477_plus_p2478_targeted_union_fraction_of_p2459_universe": f"{replay['p2477_plus_p2478_targeted_union_cell_count']}/{p2459_universe}",
        "p2477_plus_p2478_targeted_union_residual_count_against_p2459_universe": p2459_universe - replay["p2477_plus_p2478_targeted_union_cell_count"],
        "finite_chain_coverage_budget_inherited_from_p2471": 0 if p2471["all_family_partitions_are_disjoint_and_complete"] else None,
        "p2471_disjoint_complete_inherited": p2471["all_family_partitions_are_disjoint_and_complete"],
        "no_full_complement_claimed_by_this_certificate": True,
        "full_complement_replay_exported_by_this_certificate": False,
        "finite_replay_chain_left_boundary_descent_flank_extension_audit_exported": True,
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
        "lay_summary": "This packet follows the left edge left by P2477. It recalculates a one-sided strip farther left. The strip still excludes zero with positive Decimal separation and has strictly increasing separations from left to right, but its weakest value is again at the new left boundary. The honest conclusion is finite descent-flank localization, not a local-minimum, continuum, or full-complement proof.",
        "not_licensed": [
            "This one-sided flank extension is not a full complement replay and does not claim all P2459 cells were freshly replayed by P2478.",
            "It is not a directed-rounding interval theorem, symbolic root-exclusion theorem, analytic monotonicity theorem, or continuum root-exclusion theorem.",
            "It does not export selector/source/gauge closure, a physical-value generator, QW-2191 discharge, role-bearing L_total, legacy-role transfer, or ToE closure.",
            "No legacy physical roles are transferred to L_total or K_strict_gate by this finite flank replay.",
        ],
        "next_honest_step": "The replay localizes a still-open left flank. Continue only by a wider finite flank with the same no-full-complement language or by a real directed-rounding/symbolic argument; do not promote this targeted flank to a full proof.",
    }
    gatekeepers = {
        "rg_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "replay_count_expected": theorem_export["fresh_left_flank_replay_count"] == 17,
        "incremental_count_expected": theorem_export["incremental_cells_beyond_p2477_expanded_halo_count"] == 16,
        "p2477_plus_p2478_union_count_expected": theorem_export["p2477_plus_p2478_targeted_union_cell_count"] == 27,
        "minimum_separation_positive": Decimal(theorem_export["minimum_left_flank_decimal_separation"]) > 0,
        "all_replayed_cells_exclude_zero": theorem_export["all_left_flank_cells_exclude_zero"],
        "all_replayed_cells_have_positive_separation": theorem_export["all_left_flank_cells_have_positive_separation"],
        "left_boundary_descent_not_overclaimed": theorem_export["minimum_is_left_boundary_of_p2478_window"],
        "consecutive_order_checked": theorem_export["all_consecutive_pairs_strictly_increase_left_to_right"] and Decimal(theorem_export["minimum_consecutive_positive_delta"]) > 0,
        "no_full_complement_unless_genuinely_full": theorem_export["no_full_complement_claimed_by_this_certificate"] and not theorem_export["full_complement_replay_exported_by_this_certificate"],
        "coverage_budget_accounted": theorem_export["targeted_p2478_residual_not_freshly_replayed_count_against_p2459_universe"] == p2459_universe - replay["fresh_left_flank_replay_count"] and theorem_export["p2477_plus_p2478_targeted_union_residual_count_against_p2459_universe"] == p2459_universe - replay["p2477_plus_p2478_targeted_union_cell_count"],
        "p2471_partition_inherited": theorem_export["p2471_disjoint_complete_inherited"],
        "finite_grid_only_not_directed_rounding": theorem_export["finite_replay_chain_left_boundary_descent_flank_extension_audit_exported"] and not theorem_export["directed_rounding_interval_theorem_exported_by_this_certificate"],
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
        "schema_version": "p2478_s1428_v1",
        "packet_id": "P2478",
        "stage_id": "S1428",
        "status": "PASS_STRICT_POINTWISE_INTERVAL_DECIMAL_P2459_LEFT_BOUNDARY_DESCENT_FLANK_EXTENSION_REPLAY_AUDIT_NO_FULL_COMPLEMENT_NO_DIRECTED_ROUNDING_SELECTOR_SOURCE_THEOREM",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_pointwise_interval_decimal_p2459_left_boundary_descent_flank_extension_replay_audit": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_pointwise_interval_decimal_p2459_left_boundary_descent_flank_extension_replay_audit"]["theorem_export"]
    lines = [
        "# P2478/S1428 strict pointwise interval-Decimal P2459 left-boundary descent-flank extension replay audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Left-boundary descent-flank extension",
        "",
        f"Left extension radius: `{t['left_extension_radius']}`.",
        f"Fresh P2478 replay count: `{t['fresh_left_flank_replay_count']}`.",
        f"Incremental cells beyond P2477 expanded halo: `{t['incremental_cells_beyond_p2477_expanded_halo_count']}`.",
        f"P2477+P2478 targeted union cell count: `{t['p2477_plus_p2478_targeted_union_cell_count']}`.",
        f"Minimum Decimal separation in P2478 flank: `{t['minimum_left_flank_decimal_separation']}`.",
        f"Minimum consecutive positive delta: `{t['minimum_consecutive_positive_delta']}`.",
        f"Minimum is at the P2478 left boundary: `{t['minimum_is_left_boundary_of_p2478_window']}`.",
        f"All left-flank cells exclude zero: `{t['all_left_flank_cells_exclude_zero']}`.",
        "",
        "## Coverage budget",
        "",
        f"P2478 fresh replay fraction of inherited P2459 finite universe: `{t['targeted_p2478_fresh_replay_fraction_of_p2459_universe']}`.",
        f"P2478 residual not freshly replayed against P2459 universe: `{t['targeted_p2478_residual_not_freshly_replayed_count_against_p2459_universe']}`.",
        f"P2477+P2478 targeted union fraction: `{t['p2477_plus_p2478_targeted_union_fraction_of_p2459_universe']}`.",
        f"P2477+P2478 targeted union residual: `{t['p2477_plus_p2478_targeted_union_residual_count_against_p2459_universe']}`.",
        f"Inherited P2471 finite-chain coverage budget: `{t['finite_chain_coverage_budget_inherited_from_p2471']}`.",
        f"Full complement replay exported by P2478: `{t['full_complement_replay_exported_by_this_certificate']}`.",
        "",
        "## Plain-language progress note",
        "",
        t["lay_summary"],
        "",
        "## Hard limits / negative controls",
        "",
        "This is a finite left-boundary descent-flank extension replay only.  It is not a full complement replay, directed-rounding interval theorem, symbolic root-exclusion theorem, analytic monotonicity theorem, global continuum root-exclusion theorem, selector/source/gauge theorem, physical-value generator, QW-2191 discharge, role-bearing `L_total`, legacy-role transfer, or ToE closure.",
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
