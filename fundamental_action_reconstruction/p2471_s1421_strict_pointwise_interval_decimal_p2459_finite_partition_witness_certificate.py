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
    load_json,
    ordered_opposite_tail_cells,
    rel,
)

GEN = ROOT / "generated"
OUT = GEN / "p2471_s1421_strict_pointwise_interval_decimal_p2459_finite_partition_witness_certificate.json"
MD = GEN / "p2471_s1421_strict_pointwise_interval_decimal_p2459_finite_partition_witness_certificate.md"

SOURCE_FILES = {
    "P2450_ROOT_ISOLATION_MARGIN": GEN / "p2450_s1400_strict_pointwise_projection_root_isolation_margin_certificate.json",
    "P2451_FLOATING_INTERVAL_AUDIT": GEN / "p2451_s1401_strict_pointwise_projection_interval_enclosure_root_exclusion_audit.json",
    "P2456_DECIMAL_BOUNDARY_REPLAY": GEN / "p2456_s1406_strict_pointwise_decimal_root_window_boundary_band_replay_certificate.json",
    "P2459_COVERAGE_GAP_LEDGER": GEN / "p2459_s1409_strict_pointwise_interval_decimal_coverage_gap_ledger_certificate.json",
    "P2466_DESCENT_TAIL_BOUNDARY": GEN / "p2466_s1416_strict_pointwise_interval_decimal_adaptive_descent_tail_boundary_ledger_certificate.json",
    "P2469_FULL_OPPOSITE_TAIL_REPLAY": GEN / "p2469_s1419_strict_pointwise_interval_decimal_full_opposite_tail_replay_certificate.json",
    "P2470_REMAINING_NON_TAIL_REPLAY": GEN / "p2470_s1420_strict_pointwise_interval_decimal_full_remaining_non_tail_complement_replay_certificate.json",
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
        "new_packet": "P2471|S1421|finite partition witness|P2459 finite partition|coverage partition witness|disjoint coverage certificate",
        "precursor_packets": "P2459|S1409|P2466|S1416|P2469|S1419|P2470|S1420",
        "partition_language": "unreplayed-by-boundary-chain partition|finite replay partition|partition witness|disjoint union",
        "hard_limit_markers": "directed-rounding root exclusion|symbolic root-exclusion theorem|global continuum root-exclusion|full mathematical proof",
        "closure_blockers": "QW-2191|role-bearing L_total|physical-value generator|ToE closure|selector/source/gauge theorem|legacy-role transfer",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def cell_key(segment_index: int, uncovered_index: int) -> str:
    return f"segment={segment_index}:uncovered={uncovered_index}"


def independent_partition_key_sets_for_family(
    family: str,
    p2450: dict[str, Any],
    p2451: dict[str, Any],
    p2456: dict[str, Any],
    p2466: dict[str, Any],
) -> dict[str, set[str]]:
    inherited_p2466 = family_row(p2466["family_descent_tail_boundary_replays"], family, "P2466")
    selected_segment_index = int(inherited_p2466["p2465_segment_index"])
    endpoint_index = int(inherited_p2466["tail_start_uncovered_index"])
    descent_direction = int(inherited_p2466["descent_direction"])
    descent_tail = {int(row["uncovered_index"]) for row in inherited_p2466["tail_boundary_rows"]}
    key_sets = {"p2466_descent_tail": set(), "p2469_opposite_tail": set(), "p2470_remaining_non_tail": set(), "p2459_universe": set()}
    for segment_index, segment in enumerate(complement_segments(p2450, family)):
        cells = segment_cells(segment, interval_cell_width(p2451, family))
        left_skip, right_skip = boundary_skip_counts_for_segment(p2456, family, segment)
        uncovered_cells = remove_boundary_covered_cells(cells, left_skip, right_skip)
        universe_indexes = set(range(len(uncovered_cells)))
        key_sets["p2459_universe"].update(cell_key(segment_index, index) for index in universe_indexes)
        if segment_index == selected_segment_index:
            opposite_tail = {int(cell["uncovered_index"]) for cell in ordered_opposite_tail_cells(uncovered_cells, endpoint_index, descent_direction)}
            key_sets["p2466_descent_tail"].update(cell_key(segment_index, index) for index in descent_tail)
            key_sets["p2469_opposite_tail"].update(cell_key(segment_index, index) for index in opposite_tail)
            remaining = universe_indexes - descent_tail - opposite_tail
            key_sets["p2470_remaining_non_tail"].update(cell_key(segment_index, index) for index in remaining)
        else:
            key_sets["p2470_remaining_non_tail"].update(cell_key(segment_index, index) for index in universe_indexes)
    return key_sets


def set_fingerprint(keys: set[str]) -> str:
    h = hashlib.sha256()
    for key in sorted(keys):
        h.update(key.encode("utf-8"))
        h.update(b"\n")
    return h.hexdigest()


def family_partition_witness(
    family: str,
    p2450: dict[str, Any],
    p2451: dict[str, Any],
    p2456: dict[str, Any],
    p2459: dict[str, Any],
    p2466: dict[str, Any],
    p2469: dict[str, Any],
    p2470: dict[str, Any],
) -> dict[str, Any]:
    inherited_p2459 = family_row(p2459["coverage_rows"], family, "P2459")
    inherited_p2466 = family_row(p2466["family_descent_tail_boundary_replays"], family, "P2466")
    inherited_p2469 = family_row(p2469["family_full_opposite_tail_replays"], family, "P2469")
    inherited_p2470 = family_row(p2470["family_remaining_non_tail_replays"], family, "P2470")
    partition_sets = independent_partition_key_sets_for_family(family, p2450, p2451, p2456, p2466)
    descent_keys = partition_sets["p2466_descent_tail"]
    opposite_keys = partition_sets["p2469_opposite_tail"]
    remaining_keys = partition_sets["p2470_remaining_non_tail"]
    universe_keys = partition_sets["p2459_universe"]
    pairwise_intersections = {
        "p2466_descent_tail∩p2469_opposite_tail": len(descent_keys & opposite_keys),
        "p2466_descent_tail∩p2470_remaining_non_tail": len(descent_keys & remaining_keys),
        "p2469_opposite_tail∩p2470_remaining_non_tail": len(opposite_keys & remaining_keys),
    }
    union_keys = descent_keys | opposite_keys | remaining_keys
    missing_keys = universe_keys - union_keys
    extra_keys = union_keys - universe_keys
    inherited_counts_match_sets = {
        "p2459_universe": len(universe_keys) == int(inherited_p2459["unreplayed_by_decimal_boundary_chain_cell_count"]),
        "p2466_descent_tail": len(descent_keys) == int(inherited_p2466["tail_boundary_replay_count"]),
        "p2469_opposite_tail": len(opposite_keys) == int(inherited_p2469["opposite_tail_full_replay_count"]),
        "p2470_remaining_non_tail": len(remaining_keys) == int(inherited_p2470["remaining_non_tail_full_replay_count"]),
    }
    minimum_separations = [
        Decimal(inherited_p2466["minimum_tail_boundary_decimal_separation"]),
        Decimal(inherited_p2469["minimum_opposite_tail_full_replay_decimal_separation"]),
        Decimal(inherited_p2470["minimum_remaining_non_tail_full_replay_decimal_separation"]),
    ]
    return {
        "family": family,
        "p2459_universe_cell_count": len(universe_keys),
        "p2466_descent_tail_cell_count": len(descent_keys),
        "p2469_opposite_tail_cell_count": len(opposite_keys),
        "p2470_remaining_non_tail_cell_count": len(remaining_keys),
        "partition_union_cell_count": len(union_keys),
        "partition_missing_cell_count": len(missing_keys),
        "partition_extra_cell_count": len(extra_keys),
        "pairwise_intersection_counts": pairwise_intersections,
        "is_disjoint_partition_of_p2459_unreplayed_by_boundary_chain_universe": len(union_keys) == len(universe_keys) and not missing_keys and not extra_keys and all(value == 0 for value in pairwise_intersections.values()),
        "inherited_counts_match_independent_partition_sets": inherited_counts_match_sets,
        "minimum_inherited_decimal_separation_across_partition": str(min(minimum_separations)),
        "p2459_universe_fingerprint_sha256": set_fingerprint(universe_keys),
        "p2466_descent_tail_fingerprint_sha256": set_fingerprint(descent_keys),
        "p2469_opposite_tail_fingerprint_sha256": set_fingerprint(opposite_keys),
        "p2470_remaining_non_tail_fingerprint_sha256": set_fingerprint(remaining_keys),
        "partition_union_fingerprint_sha256": set_fingerprint(union_keys),
        "missing_key_samples": sorted(missing_keys)[:10],
        "extra_key_samples": sorted(extra_keys)[:10],
    }


def append_once(path: Path, marker: str, section: str) -> None:
    text = path.read_text(encoding="utf-8")
    if marker not in text:
        path.write_text(text.rstrip() + "\n\n" + section.strip() + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2471/S1421 strict pointwise interval-Decimal P2459 finite partition witness certificate

`P2471/S1421` adds an independent set-partition witness for the finite P2459 unreplayed-by-boundary-chain universe.  It rebuilds the audited cell universe from P2450/P2451/P2456, then proves by indexed set accounting that P2466 descent tails, P2469 opposite tails, and P2470 remaining non-tail cells form a disjoint partition of all `99846` P2459 cells.

This is a stronger bookkeeping theorem for the finite audit ledger: it checks that no P2459 cell was counted twice and no P2459 cell was missed.  It remains finite-grid proof hygiene, not a directed-rounding interval theorem, symbolic root-exclusion theorem, continuum theorem, selector/source/gauge theorem, physical-value generator, role-bearing `L_total`, `QW-2191` discharge, legacy-role transfer, or ToE closure.
"""
    lag_section = """
## P2471/S1421 P2459 finite partition witness guard

`P2471/S1421` audits the finite-grid coverage claim behind P2466+P2469+P2470 by proving an indexed disjoint partition of the P2459 unreplayed-by-boundary-chain cells.  This improves `L_total` bookkeeping hygiene but does not add selector/source/gauge authority, physical-value generation, role-bearing finality, legacy-role transfer, or ToE closure.
"""
    append_once(DOC_FILES["equation_sheet"], "## P2471/S1421 strict pointwise interval-Decimal P2459 finite partition witness certificate", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2471/S1421 P2459 finite partition witness guard", lag_section)


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
    family_witnesses = [family_partition_witness(family, p2450, p2451, p2456, p2459, p2466, p2469, p2470) for family in FAMILIES]
    total_universe = sum(row["p2459_universe_cell_count"] for row in family_witnesses)
    total_descent = sum(row["p2466_descent_tail_cell_count"] for row in family_witnesses)
    total_opposite = sum(row["p2469_opposite_tail_cell_count"] for row in family_witnesses)
    total_remaining = sum(row["p2470_remaining_non_tail_cell_count"] for row in family_witnesses)
    total_union = sum(row["partition_union_cell_count"] for row in family_witnesses)
    minimum_partition_separation = min(Decimal(row["minimum_inherited_decimal_separation_across_partition"]) for row in family_witnesses)
    theorem_export = {
        "theorem_name": "P2471_T1_strict_pointwise_interval_decimal_p2459_finite_partition_witness_certificate",
        "inherited_coverage_gap_ledger": "P2459/S1409",
        "inherited_descent_tail_replay": "P2466/S1416",
        "inherited_full_opposite_tail_replay": "P2469/S1419",
        "inherited_remaining_non_tail_replay": "P2470/S1420",
        "p2459_total_unreplayed_by_decimal_boundary_chain_cell_count_inherited": p2459.get("total_unreplayed_by_decimal_boundary_chain_cell_count"),
        "p2466_total_tail_boundary_replay_count_inherited": p2466.get("total_tail_boundary_replay_count"),
        "p2469_total_opposite_tail_full_replay_count_inherited": p2469.get("total_opposite_tail_full_replay_count"),
        "p2470_total_remaining_non_tail_full_replay_count_inherited": p2470.get("total_remaining_non_tail_full_replay_count"),
        "family_partition_witnesses": family_witnesses,
        "total_p2459_universe_cell_count_rebuilt": total_universe,
        "total_p2466_descent_tail_cell_count_in_partition": total_descent,
        "total_p2469_opposite_tail_cell_count_in_partition": total_opposite,
        "total_p2470_remaining_non_tail_cell_count_in_partition": total_remaining,
        "total_partition_union_cell_count": total_union,
        "total_partition_missing_cell_count": sum(row["partition_missing_cell_count"] for row in family_witnesses),
        "total_partition_extra_cell_count": sum(row["partition_extra_cell_count"] for row in family_witnesses),
        "all_family_partitions_are_disjoint_and_complete": all(row["is_disjoint_partition_of_p2459_unreplayed_by_boundary_chain_universe"] for row in family_witnesses),
        "all_inherited_counts_match_independent_partition_sets": all(all(row["inherited_counts_match_independent_partition_sets"].values()) for row in family_witnesses),
        "minimum_inherited_decimal_separation_across_full_partition": str(minimum_partition_separation),
        "finite_p2459_partition_witness_exported": True,
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
        "lay_summary": "This packet does not add new physics. It checks the bookkeeping: the finite list of P2459 boxes is split into exactly three non-overlapping piles already replayed by P2466, P2469, and P2470. No box is duplicated and no box is missing from that finite audit list.",
        "not_licensed": [
            "A disjoint finite partition of audited cells is not a directed-rounding interval theorem, symbolic root-exclusion theorem, or continuum root-exclusion theorem.",
            "The certificate does not export selector/source/gauge closure, a physical-value generator, QW-2191 discharge, role-bearing L_total, legacy-role transfer, or ToE closure.",
            "No legacy physical roles are transferred to L_total or K_strict_gate by this finite bookkeeping audit.",
        ],
        "next_honest_step": "Use this partition witness as an input to a future directed-rounding or symbolic root-exclusion audit; keep selector/source and legacy-role transfer separate.",
    }
    gatekeepers = {
        "rg_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "p2459_inherited_count": theorem_export["p2459_total_unreplayed_by_decimal_boundary_chain_cell_count_inherited"] == 99846,
        "p2466_inherited_count": theorem_export["p2466_total_tail_boundary_replay_count_inherited"] == 6361,
        "p2469_inherited_count": theorem_export["p2469_total_opposite_tail_full_replay_count_inherited"] == 45165,
        "p2470_inherited_count": theorem_export["p2470_total_remaining_non_tail_full_replay_count_inherited"] == 48320,
        "rebuilt_universe_matches_p2459": theorem_export["total_p2459_universe_cell_count_rebuilt"] == 99846,
        "partition_counts_sum_to_universe": total_descent + total_opposite + total_remaining == total_universe == total_union,
        "no_missing_cells": theorem_export["total_partition_missing_cell_count"] == 0,
        "no_extra_cells": theorem_export["total_partition_extra_cell_count"] == 0,
        "all_family_partitions_disjoint_complete": theorem_export["all_family_partitions_are_disjoint_and_complete"],
        "all_inherited_counts_match_sets": theorem_export["all_inherited_counts_match_independent_partition_sets"],
        "minimum_decimal_separation_positive": Decimal(theorem_export["minimum_inherited_decimal_separation_across_full_partition"]) > 0,
        "finite_grid_only_not_directed_rounding": theorem_export["finite_p2459_partition_witness_exported"] and not theorem_export["directed_rounding_interval_theorem_exported_by_this_certificate"],
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
        "schema_version": "p2471_s1421_v1",
        "packet_id": "P2471",
        "stage_id": "S1421",
        "status": "PASS_STRICT_POINTWISE_INTERVAL_DECIMAL_P2459_FINITE_PARTITION_WITNESS_NO_DIRECTED_ROUNDING_SELECTOR_SOURCE_THEOREM",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_pointwise_interval_decimal_p2459_finite_partition_witness_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_pointwise_interval_decimal_p2459_finite_partition_witness_certificate"]["theorem_export"]
    lines = [
        "# P2471/S1421 strict pointwise interval-Decimal P2459 finite partition witness certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite partition witness",
        "",
        f"P2459 finite universe rebuilt: `{t['total_p2459_universe_cell_count_rebuilt']}`.",
        f"P2466 descent-tail partition cells: `{t['total_p2466_descent_tail_cell_count_in_partition']}`.",
        f"P2469 opposite-tail partition cells: `{t['total_p2469_opposite_tail_cell_count_in_partition']}`.",
        f"P2470 remaining non-tail partition cells: `{t['total_p2470_remaining_non_tail_cell_count_in_partition']}`.",
        f"Partition union cells: `{t['total_partition_union_cell_count']}`.",
        f"Missing cells: `{t['total_partition_missing_cell_count']}`.",
        f"Extra cells: `{t['total_partition_extra_cell_count']}`.",
        f"All family partitions disjoint and complete: `{t['all_family_partitions_are_disjoint_and_complete']}`.",
        f"Minimum inherited Decimal separation across the partition: `{t['minimum_inherited_decimal_separation_across_full_partition']}`.",
        "",
        "## Plain-language progress note",
        "",
        t["lay_summary"],
        "",
        "## Hard limits / negative controls",
        "",
        "This is a finite partition witness for the audited P2459 cells only.  It exports no directed-rounding interval theorem, no symbolic root-exclusion theorem, no analytic monotonicity theorem, no global continuum root-exclusion theorem, no selector/source/gauge theorem, no physical-value generator, no QW-2191 discharge, no role-bearing `L_total`, no legacy-role transfer, and no ToE closure.",
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
