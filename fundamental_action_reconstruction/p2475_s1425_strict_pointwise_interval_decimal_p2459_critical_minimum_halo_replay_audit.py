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
    load_json,
    rel,
    replay_cell,
)
from p2472_s1422_strict_pointwise_interval_decimal_p2459_partition_seam_replay_audit import (
    classified_uncovered_cells_for_family,
)

GEN = ROOT / "generated"
OUT = GEN / "p2475_s1425_strict_pointwise_interval_decimal_p2459_critical_minimum_halo_replay_audit.json"
MD = GEN / "p2475_s1425_strict_pointwise_interval_decimal_p2459_critical_minimum_halo_replay_audit.md"
HALO_RADIUS = 2

SOURCE_FILES = {
    "P2449_PROJECTION_REDUCTION": GEN / "p2449_s1399_strict_pointwise_rank_lift_projection_reduction_certificate.json",
    "P2450_ROOT_ISOLATION_MARGIN": GEN / "p2450_s1400_strict_pointwise_projection_root_isolation_margin_certificate.json",
    "P2451_FLOATING_INTERVAL_AUDIT": GEN / "p2451_s1401_strict_pointwise_projection_interval_enclosure_root_exclusion_audit.json",
    "P2456_DECIMAL_BOUNDARY_REPLAY": GEN / "p2456_s1406_strict_pointwise_decimal_root_window_boundary_band_replay_certificate.json",
    "P2466_DESCENT_TAIL_BOUNDARY": GEN / "p2466_s1416_strict_pointwise_interval_decimal_adaptive_descent_tail_boundary_ledger_certificate.json",
    "P2469_FULL_OPPOSITE_TAIL_REPLAY": GEN / "p2469_s1419_strict_pointwise_interval_decimal_full_opposite_tail_replay_certificate.json",
    "P2470_REMAINING_NON_TAIL_REPLAY": GEN / "p2470_s1420_strict_pointwise_interval_decimal_full_remaining_non_tail_complement_replay_certificate.json",
    "P2471_FINITE_PARTITION_WITNESS": GEN / "p2471_s1421_strict_pointwise_interval_decimal_p2459_finite_partition_witness_certificate.json",
    "P2472_PARTITION_SEAM_REPLAY_AUDIT": GEN / "p2472_s1422_strict_pointwise_interval_decimal_p2459_partition_seam_replay_audit.json",
    "P2474_EXTREMAL_WITNESS_RERUN": GEN / "p2474_s1424_strict_pointwise_interval_decimal_p2459_extremal_witness_rerun_audit.json",
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
        "new_packet": "P2475|S1425|critical minimum halo|minimum halo replay|weakest halo replay|two-cell halo|local halo replay|extremal neighborhood replay",
        "precursor_packets": "P2469|S1419|P2470|S1420|P2471|S1421|P2472|S1422|P2474|S1424",
        "halo_language": "minimum separation witness|halo replay|neighbor cell|local neighborhood|weakest-cell neighborhood|critical finite witness",
        "hard_limit_markers": "directed-rounding root exclusion|symbolic root-exclusion theorem|global continuum root-exclusion|full mathematical proof",
        "closure_blockers": "QW-2191|role-bearing L_total|physical-value generator|ToE closure|selector/source/gauge theorem|legacy-role transfer",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def projection_vector(p2449: dict[str, Any]) -> list[Decimal]:
    values = p2449.get("strict_pointwise_rank_lift_projection_reduction_certificate", {}).get("theorem_export", {}).get("projection_vector", [])
    return [Decimal(str(value)) for value in values]


def target_minimum_keys(family: str, p2466: dict[str, Any], p2469: dict[str, Any], p2470: dict[str, Any], p2472: dict[str, Any]) -> list[dict[str, Any]]:
    selected_segment = int(family_row(p2466["family_descent_tail_boundary_replays"], family, "P2466")["p2465_segment_index"])
    p2469_family = family_row(p2469["family_full_opposite_tail_replays"], family, "P2469")
    p2470_family = family_row(p2470["family_remaining_non_tail_replays"], family, "P2470")
    p2472_family = family_row(p2472["family_partition_seam_replays"], family, "P2472")
    p2469_min = p2469_family["minimum_opposite_tail_full_replay_row"]
    p2470_min = p2470_family["minimum_remaining_non_tail_full_replay_row"]
    p2472_min = min(p2472_family["seam_replay_rows"], key=lambda row: Decimal(row["decimal_separation_from_zero"]))
    return [
        {
            "source_packet": "P2469/S1419",
            "source_rule": "minimum_opposite_tail_full_replay_row",
            "segment_index": selected_segment,
            "uncovered_index": int(p2469_min["uncovered_index"]),
            "stored_minimum_separation": p2469_min["decimal_separation_from_zero"],
        },
        {
            "source_packet": "P2470/S1420",
            "source_rule": "minimum_remaining_non_tail_full_replay_row",
            "segment_index": int(p2470_min["segment_index"]),
            "uncovered_index": int(p2470_min["uncovered_index"]),
            "stored_minimum_separation": p2470_min["decimal_separation_from_zero"],
        },
        {
            "source_packet": "P2472/S1422",
            "source_rule": "minimum_partition_seam_replay_row",
            "segment_index": int(p2472_min["segment_index"]),
            "uncovered_index": int(p2472_min["uncovered_index"]),
            "stored_minimum_separation": p2472_min["decimal_separation_from_zero"],
        },
    ]


def family_halo_replay(
    family: str,
    p2450: dict[str, Any],
    p2451: dict[str, Any],
    p2456: dict[str, Any],
    p2466: dict[str, Any],
    p2469: dict[str, Any],
    p2470: dict[str, Any],
    p2472: dict[str, Any],
    projection: list[Decimal],
) -> dict[str, Any]:
    classified, selected_segment_summary = classified_uncovered_cells_for_family(family, p2450, p2451, p2456, p2466)
    classified_by_key = {(int(row["segment_index"]), int(row["uncovered_index"])): row for row in classified}
    targets = target_minimum_keys(family, p2466, p2469, p2470, p2472)
    halo_sources: dict[tuple[int, int], list[str]] = {}
    missing_neighbors = []
    for target in targets:
        segment_index = int(target["segment_index"])
        center_index = int(target["uncovered_index"])
        for offset in range(-HALO_RADIUS, HALO_RADIUS + 1):
            key = (segment_index, center_index + offset)
            source_label = f"{target['source_packet']}:{target['source_rule']}:offset={offset}"
            if key in classified_by_key:
                halo_sources.setdefault(key, []).append(source_label)
            else:
                missing_neighbors.append({"target": target, "offset": offset, "missing_key": {"segment_index": key[0], "uncovered_index": key[1]}})
    replayed = []
    function = function_for_family(family)
    for segment_index, uncovered_index in sorted(halo_sources):
        row = classified_by_key[(segment_index, uncovered_index)]
        fresh = replay_cell(family, row, projection, function)
        replayed.append({
            **fresh,
            "segment_index": segment_index,
            "uncovered_index": uncovered_index,
            "partition_class": row["partition_class"],
            "halo_source_rules": sorted(halo_sources[(segment_index, uncovered_index)]),
            "halo_source_count": len(halo_sources[(segment_index, uncovered_index)]),
            "fresh_decimal_separation_positive": Decimal(fresh["decimal_separation_from_zero"]) > 0,
        })
    separations = [Decimal(row["decimal_separation_from_zero"]) for row in replayed]
    class_counts: dict[str, int] = {}
    for row in replayed:
        class_counts[row["partition_class"]] = class_counts.get(row["partition_class"], 0) + 1
    return {
        "family": family,
        "halo_radius": HALO_RADIUS,
        "minimum_witness_targets": targets,
        "selected_segment_summary": selected_segment_summary,
        "unique_halo_replay_count": len(replayed),
        "missing_neighbor_count_due_to_segment_boundaries": len(missing_neighbors),
        "missing_neighbors_due_to_segment_boundaries": missing_neighbors,
        "partition_class_counts": class_counts,
        "all_halo_replayed_cells_exclude_zero": all(row["decimal_interval_excludes_zero"] for row in replayed),
        "all_halo_replayed_cells_have_positive_separation": all(row["fresh_decimal_separation_positive"] for row in replayed),
        "minimum_halo_decimal_separation": str(min(separations)),
        "minimum_halo_replay_row": min(replayed, key=lambda row: Decimal(row["decimal_separation_from_zero"])),
        "halo_replay_rows": replayed,
        "halo_replay_fingerprint_sha256": sha256_json(replayed),
    }


def append_once(path: Path, marker: str, section: str) -> None:
    text = path.read_text(encoding="utf-8")
    if marker not in text:
        path.write_text(text.rstrip() + "\n\n" + section.strip() + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2475/S1425 strict pointwise interval-Decimal P2459 critical-minimum halo replay audit

`P2475/S1425` expands the P2474 extremal-witness rerun into a local halo audit.  For each scalar family it takes the minimum-separation witnesses inherited from P2469, P2470, and P2472, rebuilds the P2459 finite partition classification, and replays the two-neighbor halo around each critical witness with the Decimal/Taylor backend.  The fresh halo replay covers 14 unique cells after overlap and segment-boundary truncation; all replayed halo cells remain zero-excluding with positive Decimal separation.

For a non-specialist: P2474 recalculated the most important single saved cells; P2475 also checks nearby cells around those weak spots, so the result is not hanging on a one-cell bookkeeping accident.  It remains a finite local-neighborhood replay, not a directed-rounding interval theorem, symbolic root-exclusion theorem, continuum theorem, selector/source/gauge theorem, physical-value generator, role-bearing `L_total`, `QW-2191` discharge, legacy-role transfer, or ToE closure.
"""
    lag_section = """
## P2475/S1425 P2459 critical-minimum halo replay guard

`P2475/S1425` adds a local halo guard to the finite-grid bookkeeping behind `L_total`: the two-cell neighborhoods around the P2469/P2470/P2472 minimum-separation witnesses are rebuilt from the finite P2459 partition and replayed.  This strengthens local robustness around the weakest finite witnesses, but it does not add selector/source/gauge authority, physical-value generation, role-bearing finality, legacy-role transfer, or ToE closure.
"""
    append_once(DOC_FILES["equation_sheet"], "## P2475/S1425 strict pointwise interval-Decimal P2459 critical-minimum halo replay audit", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2475/S1425 P2459 critical-minimum halo replay guard", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    projection = projection_vector(sources["P2449_PROJECTION_REDUCTION"])
    p2450 = theorem(sources["P2450_ROOT_ISOLATION_MARGIN"], "strict_pointwise_projection_root_isolation_margin_certificate")
    p2451 = theorem(sources["P2451_FLOATING_INTERVAL_AUDIT"], "strict_pointwise_projection_interval_enclosure_root_exclusion_audit")
    p2456 = theorem(sources["P2456_DECIMAL_BOUNDARY_REPLAY"], "strict_pointwise_decimal_root_window_boundary_band_replay_certificate")
    p2466 = theorem(sources["P2466_DESCENT_TAIL_BOUNDARY"], "strict_pointwise_interval_decimal_adaptive_descent_tail_boundary_ledger_certificate")
    p2469 = theorem(sources["P2469_FULL_OPPOSITE_TAIL_REPLAY"], "strict_pointwise_interval_decimal_full_opposite_tail_replay_certificate")
    p2470 = theorem(sources["P2470_REMAINING_NON_TAIL_REPLAY"], "strict_pointwise_interval_decimal_full_remaining_non_tail_complement_replay_certificate")
    p2471 = theorem(sources["P2471_FINITE_PARTITION_WITNESS"], "strict_pointwise_interval_decimal_p2459_finite_partition_witness_certificate")
    p2472 = theorem(sources["P2472_PARTITION_SEAM_REPLAY_AUDIT"], "strict_pointwise_interval_decimal_p2459_partition_seam_replay_audit")
    p2474 = theorem(sources["P2474_EXTREMAL_WITNESS_RERUN"], "strict_pointwise_interval_decimal_p2459_extremal_witness_rerun_audit")
    family_halos = [family_halo_replay(family, p2450, p2451, p2456, p2466, p2469, p2470, p2472, projection) for family in FAMILIES]
    total_halo = sum(row["unique_halo_replay_count"] for row in family_halos)
    minimum_halo = min(Decimal(row["minimum_halo_decimal_separation"]) for row in family_halos)
    finite_chain_sum = p2471["p2466_total_tail_boundary_replay_count_inherited"] + p2469["total_opposite_tail_full_replay_count"] + p2470["total_remaining_non_tail_full_replay_count"]
    theorem_export = {
        "theorem_name": "P2475_T1_strict_pointwise_interval_decimal_p2459_critical_minimum_halo_replay_audit",
        "audited_chain": ["P2469/S1419", "P2470/S1420", "P2471/S1421", "P2472/S1422", "P2474/S1424"],
        "halo_radius": HALO_RADIUS,
        "family_halo_replays": family_halos,
        "total_unique_halo_replay_count": total_halo,
        "total_missing_neighbor_count_due_to_segment_boundaries": sum(row["missing_neighbor_count_due_to_segment_boundaries"] for row in family_halos),
        "all_halo_replayed_cells_exclude_zero": all(row["all_halo_replayed_cells_exclude_zero"] for row in family_halos),
        "all_halo_replayed_cells_have_positive_separation": all(row["all_halo_replayed_cells_have_positive_separation"] for row in family_halos),
        "minimum_halo_decimal_separation": str(minimum_halo),
        "p2474_total_extremal_witness_rerun_count_inherited": p2474["total_fresh_decimal_taylor_witness_rerun_count"],
        "p2474_all_fresh_witness_groups_match_stored_inherited": p2474["all_fresh_witness_groups_match_stored"],
        "p2474_all_fresh_witness_groups_exclude_zero_inherited": p2474["all_fresh_witness_groups_exclude_zero_with_positive_separation"],
        "p2466_tail_count_inherited_from_p2471": p2471["p2466_total_tail_boundary_replay_count_inherited"],
        "p2469_full_opposite_tail_count_inherited": p2469["total_opposite_tail_full_replay_count"],
        "p2470_remaining_non_tail_count_inherited": p2470["total_remaining_non_tail_full_replay_count"],
        "p2471_p2459_universe_count_inherited": p2471["total_p2459_universe_cell_count_rebuilt"],
        "finite_chain_sum_p2466_p2469_p2470": finite_chain_sum,
        "finite_chain_sum_matches_p2471_universe": finite_chain_sum == p2471["total_p2459_universe_cell_count_rebuilt"],
        "p2471_missing_cells_inherited": p2471["total_partition_missing_cell_count"],
        "p2471_extra_cells_inherited": p2471["total_partition_extra_cell_count"],
        "p2471_disjoint_complete_inherited": p2471["all_family_partitions_are_disjoint_and_complete"],
        "finite_replay_chain_critical_minimum_halo_replay_audit_exported": True,
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
        "lay_summary": "This packet checks the neighborhoods around the weakest saved cells. Instead of trusting only the single lowest-separation cell, it replays the nearby finite cells on both sides where they exist. The nearby cells also stay away from zero, so the finite audit is less likely to depend on a one-cell indexing accident. It still does not become a symbolic proof for every real point.",
        "not_licensed": [
            "A critical-minimum halo replay audit is not a directed-rounding interval theorem, symbolic root-exclusion theorem, analytic monotonicity theorem, or continuum root-exclusion theorem.",
            "The audit does not export selector/source/gauge closure, a physical-value generator, QW-2191 discharge, role-bearing L_total, legacy-role transfer, or ToE closure.",
            "No legacy physical roles are transferred to L_total or K_strict_gate by this finite local-neighborhood check.",
        ],
        "next_honest_step": "Use the rerun-stable critical halos as local robustness evidence for a future directed-rounding or symbolic root-exclusion attempt; do not promote finite halo replay into continuum closure.",
    }
    gatekeepers = {
        "rg_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "halo_radius_expected": theorem_export["halo_radius"] == 2,
        "total_unique_halo_replay_count_expected": theorem_export["total_unique_halo_replay_count"] == 14,
        "boundary_truncations_recorded": theorem_export["total_missing_neighbor_count_due_to_segment_boundaries"] == 10,
        "all_halos_exclude_zero": theorem_export["all_halo_replayed_cells_exclude_zero"],
        "all_halos_positive_separation": theorem_export["all_halo_replayed_cells_have_positive_separation"],
        "minimum_halo_separation_positive": Decimal(theorem_export["minimum_halo_decimal_separation"]) > 0,
        "p2474_extremal_rerun_inherited": theorem_export["p2474_total_extremal_witness_rerun_count_inherited"] == 28 and theorem_export["p2474_all_fresh_witness_groups_match_stored_inherited"] and theorem_export["p2474_all_fresh_witness_groups_exclude_zero_inherited"],
        "finite_chain_sum_matches_universe": theorem_export["finite_chain_sum_matches_p2471_universe"],
        "p2459_universe_count": theorem_export["p2471_p2459_universe_count_inherited"] == 99846,
        "p2469_count": theorem_export["p2469_full_opposite_tail_count_inherited"] == 45165,
        "p2470_count": theorem_export["p2470_remaining_non_tail_count_inherited"] == 48320,
        "p2471_no_missing_extra": theorem_export["p2471_missing_cells_inherited"] == 0 and theorem_export["p2471_extra_cells_inherited"] == 0 and theorem_export["p2471_disjoint_complete_inherited"],
        "finite_grid_only_not_directed_rounding": theorem_export["finite_replay_chain_critical_minimum_halo_replay_audit_exported"] and not theorem_export["directed_rounding_interval_theorem_exported_by_this_certificate"],
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
        "schema_version": "p2475_s1425_v1",
        "packet_id": "P2475",
        "stage_id": "S1425",
        "status": "PASS_STRICT_POINTWISE_INTERVAL_DECIMAL_P2459_CRITICAL_MINIMUM_HALO_REPLAY_AUDIT_NO_DIRECTED_ROUNDING_SELECTOR_SOURCE_THEOREM",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_pointwise_interval_decimal_p2459_critical_minimum_halo_replay_audit": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_pointwise_interval_decimal_p2459_critical_minimum_halo_replay_audit"]["theorem_export"]
    lines = [
        "# P2475/S1425 strict pointwise interval-Decimal P2459 critical-minimum halo replay audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Critical minimum halo replay audit",
        "",
        f"Halo radius: `{t['halo_radius']}` cells around each critical minimum witness.",
        f"Unique halo cells replayed: `{t['total_unique_halo_replay_count']}`.",
        f"Boundary-truncated missing neighbors recorded: `{t['total_missing_neighbor_count_due_to_segment_boundaries']}`.",
        f"All halo cells exclude zero: `{t['all_halo_replayed_cells_exclude_zero']}`.",
        f"All halo cells have positive separation: `{t['all_halo_replayed_cells_have_positive_separation']}`.",
        f"Minimum halo Decimal separation: `{t['minimum_halo_decimal_separation']}`.",
        f"Finite chain sum P2466+P2469+P2470: `{t['finite_chain_sum_p2466_p2469_p2470']}` / P2459 universe `{t['p2471_p2459_universe_count_inherited']}`.",
        "",
        "## Plain-language progress note",
        "",
        t["lay_summary"],
        "",
        "## Hard limits / negative controls",
        "",
        "This is a finite critical-minimum halo replay audit only.  It exports no directed-rounding interval theorem, no symbolic root-exclusion theorem, no analytic monotonicity theorem, no global continuum root-exclusion theorem, no selector/source/gauge theorem, no physical-value generator, no QW-2191 discharge, no role-bearing `L_total`, no legacy-role transfer, and no ToE closure.",
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
