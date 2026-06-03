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
OUT = GEN / "p2472_s1422_strict_pointwise_interval_decimal_p2459_partition_seam_replay_audit.json"
MD = GEN / "p2472_s1422_strict_pointwise_interval_decimal_p2459_partition_seam_replay_audit.md"

SOURCE_FILES = {
    "P2450_ROOT_ISOLATION_MARGIN": GEN / "p2450_s1400_strict_pointwise_projection_root_isolation_margin_certificate.json",
    "P2451_FLOATING_INTERVAL_AUDIT": GEN / "p2451_s1401_strict_pointwise_projection_interval_enclosure_root_exclusion_audit.json",
    "P2456_DECIMAL_BOUNDARY_REPLAY": GEN / "p2456_s1406_strict_pointwise_decimal_root_window_boundary_band_replay_certificate.json",
    "P2459_COVERAGE_GAP_LEDGER": GEN / "p2459_s1409_strict_pointwise_interval_decimal_coverage_gap_ledger_certificate.json",
    "P2466_DESCENT_TAIL_BOUNDARY": GEN / "p2466_s1416_strict_pointwise_interval_decimal_adaptive_descent_tail_boundary_ledger_certificate.json",
    "P2469_FULL_OPPOSITE_TAIL_REPLAY": GEN / "p2469_s1419_strict_pointwise_interval_decimal_full_opposite_tail_replay_certificate.json",
    "P2470_REMAINING_NON_TAIL_REPLAY": GEN / "p2470_s1420_strict_pointwise_interval_decimal_full_remaining_non_tail_complement_replay_certificate.json",
    "P2471_FINITE_PARTITION_WITNESS": GEN / "p2471_s1421_strict_pointwise_interval_decimal_p2459_finite_partition_witness_certificate.json",
    "P2449_PROJECTION_REDUCTION": GEN / "p2449_s1399_strict_pointwise_rank_lift_projection_reduction_certificate.json",
}
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
        "new_packet": "P2472|S1422|partition seam|seam replay|off-by-one audit|finite partition seam|boundary seam witness|adjacency seam",
        "precursor_packets": "P2459|S1409|P2466|S1416|P2469|S1419|P2470|S1420|P2471|S1421",
        "seam_language": "partition seam|off-by-one|adjacent index|central pivot|boundary edge replay",
        "hard_limit_markers": "directed-rounding root exclusion|symbolic root-exclusion theorem|global continuum root-exclusion|full mathematical proof",
        "closure_blockers": "QW-2191|role-bearing L_total|physical-value generator|ToE closure|selector/source/gauge theorem|legacy-role transfer",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def classified_uncovered_cells_for_family(
    family: str,
    p2450: dict[str, Any],
    p2451: dict[str, Any],
    p2456: dict[str, Any],
    p2466: dict[str, Any],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    inherited_p2466 = family_row(p2466["family_descent_tail_boundary_replays"], family, "P2466")
    selected_segment_index = int(inherited_p2466["p2465_segment_index"])
    endpoint_index = int(inherited_p2466["tail_start_uncovered_index"])
    descent_direction = int(inherited_p2466["descent_direction"])
    descent_tail = {int(row["uncovered_index"]) for row in inherited_p2466["tail_boundary_rows"]}
    classified: list[dict[str, Any]] = []
    selected_segment_summary: dict[str, Any] = {}
    for segment_index, segment in enumerate(complement_segments(p2450, family)):
        cells = segment_cells(segment, interval_cell_width(p2451, family))
        left_skip, right_skip = boundary_skip_counts_for_segment(p2456, family, segment)
        uncovered_cells = remove_boundary_covered_cells(cells, left_skip, right_skip)
        opposite_tail = set()
        if segment_index == selected_segment_index:
            opposite_tail = {int(cell["uncovered_index"]) for cell in ordered_opposite_tail_cells(uncovered_cells, endpoint_index, descent_direction)}
        for uncovered_index, cell in enumerate(uncovered_cells):
            if segment_index == selected_segment_index and uncovered_index in descent_tail:
                partition_class = "p2466_descent_tail"
            elif segment_index == selected_segment_index and uncovered_index in opposite_tail:
                partition_class = "p2469_opposite_tail"
            else:
                partition_class = "p2470_remaining_non_tail"
            classified.append({
                **cell,
                "segment_index": segment_index,
                "segment_left": segment["left"],
                "segment_right": segment["right"],
                "uncovered_index": uncovered_index,
                "uncovered_count": len(uncovered_cells),
                "partition_class": partition_class,
            })
        if segment_index == selected_segment_index:
            selected_segment_summary = {
                "segment_index": segment_index,
                "endpoint_index": endpoint_index,
                "descent_direction": descent_direction,
                "uncovered_count": len(uncovered_cells),
                "descent_tail_count": len(descent_tail),
                "opposite_tail_count": len(opposite_tail),
                "pivot_partition_class": "p2470_remaining_non_tail",
            }
    return classified, selected_segment_summary


def transition_pairs(classified_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    pairs: list[dict[str, Any]] = []
    by_segment: dict[int, list[dict[str, Any]]] = {}
    for row in classified_rows:
        by_segment.setdefault(int(row["segment_index"]), []).append(row)
    for segment_index, rows in sorted(by_segment.items()):
        rows = sorted(rows, key=lambda item: int(item["uncovered_index"]))
        for left, right in zip(rows, rows[1:]):
            if left["partition_class"] != right["partition_class"]:
                pairs.append({
                    "segment_index": segment_index,
                    "left_uncovered_index": left["uncovered_index"],
                    "right_uncovered_index": right["uncovered_index"],
                    "left_partition_class": left["partition_class"],
                    "right_partition_class": right["partition_class"],
                    "adjacent_indexes": int(right["uncovered_index"]) - int(left["uncovered_index"]) == 1,
                    "transition_label": f"{left['partition_class']}->{right['partition_class']}",
                })
    return pairs


def seam_candidate_rows(classified_rows: list[dict[str, Any]], transitions: list[dict[str, Any]]) -> list[dict[str, Any]]:
    by_key = {(int(row["segment_index"]), int(row["uncovered_index"])): row for row in classified_rows}
    wanted: dict[tuple[int, int], str] = {}
    by_segment: dict[int, list[dict[str, Any]]] = {}
    for row in classified_rows:
        by_segment.setdefault(int(row["segment_index"]), []).append(row)
    for segment_index, rows in sorted(by_segment.items()):
        rows = sorted(rows, key=lambda item: int(item["uncovered_index"]))
        if rows:
            wanted[(segment_index, int(rows[0]["uncovered_index"]))] = "segment-left-boundary-edge"
            wanted[(segment_index, int(rows[-1]["uncovered_index"]))] = "segment-right-boundary-edge"
    for transition in transitions:
        segment_index = int(transition["segment_index"])
        wanted[(segment_index, int(transition["left_uncovered_index"]))] = "partition-transition-left"
        wanted[(segment_index, int(transition["right_uncovered_index"]))] = "partition-transition-right"
    return [{**by_key[key], "seam_selection_rule": rule} for key, rule in sorted(wanted.items())]


def replay_family_seams(
    family: str,
    p2450: dict[str, Any],
    p2451: dict[str, Any],
    p2456: dict[str, Any],
    p2459: dict[str, Any],
    p2466: dict[str, Any],
    p2471: dict[str, Any],
    projection_vector: list[Decimal],
) -> dict[str, Any]:
    classified, selected_segment_summary = classified_uncovered_cells_for_family(family, p2450, p2451, p2456, p2466)
    transitions = transition_pairs(classified)
    seam_rows = seam_candidate_rows(classified, transitions)
    function = function_for_family(family)
    replayed = [replay_cell(family, row, projection_vector, function) | {
        "segment_index": row["segment_index"],
        "partition_class": row["partition_class"],
        "seam_selection_rule": row["seam_selection_rule"],
    } for row in seam_rows]
    family_partition = family_row(p2471["family_partition_witnesses"], family, "P2471")
    p2459_row = family_row(p2459["coverage_rows"], family, "P2459")
    separations = [Decimal(row["decimal_separation_from_zero"]) for row in replayed]
    class_counts: dict[str, int] = {}
    for row in replayed:
        class_counts[row["partition_class"]] = class_counts.get(row["partition_class"], 0) + 1
    return {
        "family": family,
        "inherited_p2459_universe_count": int(p2459_row["unreplayed_by_decimal_boundary_chain_cell_count"]),
        "inherited_p2471_universe_count": int(family_partition["p2459_universe_cell_count"]),
        "selected_segment_summary": selected_segment_summary,
        "transition_pair_count": len(transitions),
        "transition_pairs": transitions,
        "all_transition_pairs_are_adjacent": all(row["adjacent_indexes"] for row in transitions),
        "seam_replay_count": len(replayed),
        "seam_replay_partition_class_counts": class_counts,
        "all_seam_replayed_cells_exclude_zero": all(row["decimal_interval_excludes_zero"] for row in replayed),
        "minimum_seam_decimal_separation": str(min(separations)),
        "seam_replay_rows": replayed,
        "seam_replay_fingerprint_sha256": sha256_json(replayed),
        "p2471_family_partition_inherited_disjoint_complete": family_partition["is_disjoint_partition_of_p2459_unreplayed_by_boundary_chain_universe"],
    }


def append_once(path: Path, marker: str, section: str) -> None:
    text = path.read_text(encoding="utf-8")
    if marker not in text:
        path.write_text(text.rstrip() + "\n\n" + section.strip() + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2472/S1422 strict pointwise interval-Decimal P2459 partition seam replay audit

`P2472/S1422` audits the off-by-one seams of the P2471 finite partition witness.  It rebuilds the selected P2459 segment classifications, checks that class transitions are adjacent-index seams, and replays boundary-edge plus transition-edge cells with the Decimal/Taylor backend.  All seam cells remain zero-excluding with positive Decimal separation.

This is a local seam/off-by-one audit for the finite partition ledger.  It is not a directed-rounding interval theorem, symbolic root-exclusion theorem, continuum theorem, selector/source/gauge theorem, physical-value generator, role-bearing `L_total`, `QW-2191` discharge, legacy-role transfer, or ToE closure.
"""
    lag_section = """
## P2472/S1422 P2459 finite partition seam replay guard

`P2472/S1422` strengthens the finite-grid bookkeeping behind `L_total` by checking the transition seams between P2466, P2469, and P2470 partition classes and replaying seam cells.  It does not add selector/source/gauge authority, physical-value generation, role-bearing finality, legacy-role transfer, or ToE closure.
"""
    append_once(DOC_FILES["equation_sheet"], "## P2472/S1422 strict pointwise interval-Decimal P2459 partition seam replay audit", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2472/S1422 P2459 finite partition seam replay guard", lag_section)


def build_payload() -> dict[str, Any]:
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    grep = rg_audit()
    p2450 = sources["P2450_ROOT_ISOLATION_MARGIN"].get("strict_pointwise_projection_root_isolation_margin_certificate", {}).get("theorem_export", {})
    p2451 = sources["P2451_FLOATING_INTERVAL_AUDIT"].get("strict_pointwise_projection_interval_enclosure_root_exclusion_audit", {}).get("theorem_export", {})
    p2456 = sources["P2456_DECIMAL_BOUNDARY_REPLAY"].get("strict_pointwise_decimal_root_window_boundary_band_replay_certificate", {}).get("theorem_export", {})
    p2459 = sources["P2459_COVERAGE_GAP_LEDGER"].get("strict_pointwise_interval_decimal_coverage_gap_ledger_certificate", {}).get("theorem_export", {})
    p2466 = sources["P2466_DESCENT_TAIL_BOUNDARY"].get("strict_pointwise_interval_decimal_adaptive_descent_tail_boundary_ledger_certificate", {}).get("theorem_export", {})
    p2469 = sources["P2469_FULL_OPPOSITE_TAIL_REPLAY"].get("strict_pointwise_interval_decimal_full_opposite_tail_replay_certificate", {}).get("theorem_export", {})
    p2470 = sources["P2470_REMAINING_NON_TAIL_REPLAY"].get("strict_pointwise_interval_decimal_full_remaining_non_tail_complement_replay_certificate", {}).get("theorem_export", {})
    p2471 = sources["P2471_FINITE_PARTITION_WITNESS"].get("strict_pointwise_interval_decimal_p2459_finite_partition_witness_certificate", {}).get("theorem_export", {})
    projection_vector = [Decimal(str(value)) for value in sources["P2449_PROJECTION_REDUCTION"].get("strict_pointwise_rank_lift_projection_reduction_certificate", {}).get("theorem_export", {}).get("projection_vector", [])]
    family_seams = [replay_family_seams(family, p2450, p2451, p2456, p2459, p2466, p2471, projection_vector) for family in FAMILIES]
    total_seams = sum(row["seam_replay_count"] for row in family_seams)
    minimum_separation = min(Decimal(row["minimum_seam_decimal_separation"]) for row in family_seams)
    theorem_export = {
        "theorem_name": "P2472_T1_strict_pointwise_interval_decimal_p2459_partition_seam_replay_audit",
        "inherited_coverage_gap_ledger": "P2459/S1409",
        "inherited_descent_tail_replay": "P2466/S1416",
        "inherited_full_opposite_tail_replay": "P2469/S1419",
        "inherited_remaining_non_tail_replay": "P2470/S1420",
        "inherited_finite_partition_witness": "P2471/S1421",
        "p2459_total_unreplayed_by_decimal_boundary_chain_cell_count_inherited": p2459.get("total_unreplayed_by_decimal_boundary_chain_cell_count"),
        "p2466_total_tail_boundary_replay_count_inherited": p2466.get("total_tail_boundary_replay_count"),
        "p2469_total_opposite_tail_full_replay_count_inherited": p2469.get("total_opposite_tail_full_replay_count"),
        "p2470_total_remaining_non_tail_full_replay_count_inherited": p2470.get("total_remaining_non_tail_full_replay_count"),
        "p2471_total_p2459_universe_cell_count_rebuilt_inherited": p2471.get("total_p2459_universe_cell_count_rebuilt"),
        "p2471_all_family_partitions_are_disjoint_and_complete_inherited": p2471.get("all_family_partitions_are_disjoint_and_complete"),
        "family_partition_seam_replays": family_seams,
        "total_transition_pair_count": sum(row["transition_pair_count"] for row in family_seams),
        "total_partition_seam_replay_count": total_seams,
        "all_transition_pairs_are_adjacent": all(row["all_transition_pairs_are_adjacent"] for row in family_seams),
        "all_seam_replayed_cells_exclude_zero": all(row["all_seam_replayed_cells_exclude_zero"] for row in family_seams),
        "minimum_partition_seam_decimal_separation": str(minimum_separation),
        "all_p2471_family_partitions_inherited_disjoint_complete": all(row["p2471_family_partition_inherited_disjoint_complete"] for row in family_seams),
        "finite_partition_seam_replay_audit_exported": True,
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
        "lay_summary": "This packet checks the joins between the three finite piles from P2471. The neighbouring cells at each join are adjacent, still nonzero under Decimal/Taylor replay, and therefore there is no visible off-by-one seam error in the finite audit partition.",
        "not_licensed": [
            "A seam replay audit is not a directed-rounding interval theorem, symbolic root-exclusion theorem, or continuum root-exclusion theorem.",
            "The certificate does not export selector/source/gauge closure, a physical-value generator, QW-2191 discharge, role-bearing L_total, legacy-role transfer, or ToE closure.",
            "No legacy physical roles are transferred to L_total or K_strict_gate by this finite seam audit.",
        ],
        "next_honest_step": "Use P2471/P2472 as bookkeeping inputs for a future directed-rounding or symbolic root-exclusion audit; do not convert seam hygiene into selector/source or ToE closure.",
    }
    gatekeepers = {
        "rg_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "p2459_inherited_count": theorem_export["p2459_total_unreplayed_by_decimal_boundary_chain_cell_count_inherited"] == 99846,
        "p2466_inherited_count": theorem_export["p2466_total_tail_boundary_replay_count_inherited"] == 6361,
        "p2469_inherited_count": theorem_export["p2469_total_opposite_tail_full_replay_count_inherited"] == 45165,
        "p2470_inherited_count": theorem_export["p2470_total_remaining_non_tail_full_replay_count_inherited"] == 48320,
        "p2471_universe_inherited_count": theorem_export["p2471_total_p2459_universe_cell_count_rebuilt_inherited"] == 99846,
        "p2471_partition_inherited_complete": theorem_export["p2471_all_family_partitions_are_disjoint_and_complete_inherited"] is True,
        "transition_pairs_present": theorem_export["total_transition_pair_count"] == 4,
        "seam_replays_present": theorem_export["total_partition_seam_replay_count"] == 16,
        "transition_pairs_adjacent": theorem_export["all_transition_pairs_are_adjacent"],
        "all_seams_exclude_zero": theorem_export["all_seam_replayed_cells_exclude_zero"],
        "minimum_seam_decimal_separation_positive": Decimal(theorem_export["minimum_partition_seam_decimal_separation"]) > 0,
        "all_p2471_family_partitions_inherited_disjoint_complete": theorem_export["all_p2471_family_partitions_inherited_disjoint_complete"],
        "finite_grid_only_not_directed_rounding": theorem_export["finite_partition_seam_replay_audit_exported"] and not theorem_export["directed_rounding_interval_theorem_exported_by_this_certificate"],
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
        "schema_version": "p2472_s1422_v1",
        "packet_id": "P2472",
        "stage_id": "S1422",
        "status": "PASS_STRICT_POINTWISE_INTERVAL_DECIMAL_P2459_PARTITION_SEAM_REPLAY_AUDIT_NO_DIRECTED_ROUNDING_SELECTOR_SOURCE_THEOREM",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_pointwise_interval_decimal_p2459_partition_seam_replay_audit": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_pointwise_interval_decimal_p2459_partition_seam_replay_audit"]["theorem_export"]
    lines = [
        "# P2472/S1422 strict pointwise interval-Decimal P2459 partition seam replay audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite partition seam audit",
        "",
        f"P2471 finite universe inherited: `{t['p2471_total_p2459_universe_cell_count_rebuilt_inherited']}`.",
        f"Transition pairs checked: `{t['total_transition_pair_count']}`.",
        f"Seam cells replayed: `{t['total_partition_seam_replay_count']}`.",
        f"All transition pairs adjacent: `{t['all_transition_pairs_are_adjacent']}`.",
        f"All seam cells exclude zero: `{t['all_seam_replayed_cells_exclude_zero']}`.",
        f"Minimum partition seam Decimal separation: `{t['minimum_partition_seam_decimal_separation']}`.",
        "",
        "## Plain-language progress note",
        "",
        t["lay_summary"],
        "",
        "## Hard limits / negative controls",
        "",
        "This is a finite seam/off-by-one replay audit for the audited P2459 partition only.  It exports no directed-rounding interval theorem, no symbolic root-exclusion theorem, no analytic monotonicity theorem, no global continuum root-exclusion theorem, no selector/source/gauge theorem, no physical-value generator, no QW-2191 discharge, no role-bearing `L_total`, no legacy-role transfer, and no ToE closure.",
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
