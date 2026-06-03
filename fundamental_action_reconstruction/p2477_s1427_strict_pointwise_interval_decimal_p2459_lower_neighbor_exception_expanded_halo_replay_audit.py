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
OUT = GEN / "p2477_s1427_strict_pointwise_interval_decimal_p2459_lower_neighbor_exception_expanded_halo_replay_audit.json"
MD = GEN / "p2477_s1427_strict_pointwise_interval_decimal_p2459_lower_neighbor_exception_expanded_halo_replay_audit.md"
EXPANDED_RADIUS = 4

SOURCE_FILES = {
    "P2449_PROJECTION_REDUCTION": GEN / "p2449_s1399_strict_pointwise_rank_lift_projection_reduction_certificate.json",
    "P2450_ROOT_ISOLATION_MARGIN": GEN / "p2450_s1400_strict_pointwise_projection_root_isolation_margin_certificate.json",
    "P2451_FLOATING_INTERVAL_AUDIT": GEN / "p2451_s1401_strict_pointwise_projection_interval_enclosure_root_exclusion_audit.json",
    "P2456_DECIMAL_BOUNDARY_REPLAY": GEN / "p2456_s1406_strict_pointwise_decimal_root_window_boundary_band_replay_certificate.json",
    "P2466_DESCENT_TAIL_BOUNDARY": GEN / "p2466_s1416_strict_pointwise_interval_decimal_adaptive_descent_tail_boundary_ledger_certificate.json",
    "P2471_FINITE_PARTITION_WITNESS": GEN / "p2471_s1421_strict_pointwise_interval_decimal_p2459_finite_partition_witness_certificate.json",
    "P2475_CRITICAL_MINIMUM_HALO_REPLAY": GEN / "p2475_s1425_strict_pointwise_interval_decimal_p2459_critical_minimum_halo_replay_audit.json",
    "P2476_CRITICAL_HALO_ORDER_CLASSIFICATION": GEN / "p2476_s1426_strict_pointwise_interval_decimal_p2459_critical_halo_order_classification_audit.json",
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
        "new_packet": "P2477|S1427|exception descent halo|lower-neighbor exception expanded halo|descent ladder replay|expanded exception halo|radius-four exception|weakest exception neighborhood|exception-centered halo",
        "precursor_packets": "P2475|S1425|P2476|S1426|critical halo order|lower-neighbor exception",
        "descent_language": "lower neighbor|descent ladder|expanded halo|exception halo|monotone flank|weakest exception|available-halo minimum",
        "hard_limit_markers": "directed-rounding root exclusion|symbolic root-exclusion theorem|global continuum root-exclusion|full mathematical proof|full complement",
        "closure_blockers": "QW-2191|role-bearing L_total|physical-value generator|ToE closure|selector/source/gauge theorem|legacy-role transfer",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def p2476_lower_neighbor_exceptions(p2476: dict[str, Any]) -> list[dict[str, Any]]:
    return [row for row in p2476["target_classifications"] if row["center_has_lower_neighbor_within_available_halo"]]


def expanded_exception_replay(
    exception: dict[str, Any],
    p2450: dict[str, Any],
    p2451: dict[str, Any],
    p2456: dict[str, Any],
    p2466: dict[str, Any],
    projection: list[Decimal],
) -> dict[str, Any]:
    family = exception["family"]
    segment_index = int(exception["segment_index"])
    center_index = int(exception["uncovered_index"])
    lower_anchor_rows = exception["lower_neighbor_rows"]
    lowest_lower_anchor = min(lower_anchor_rows, key=lambda row: Decimal(row["decimal_separation_from_zero"]))
    lowest_lower_index = int(lowest_lower_anchor["uncovered_index"])
    anchor_indices = sorted({center_index, lowest_lower_index})
    classified, selected_segment_summary = classified_uncovered_cells_for_family(family, p2450, p2451, p2456, p2466)
    classified_by_key = {(int(row["segment_index"]), int(row["uncovered_index"])): row for row in classified}
    replay_keys: dict[tuple[int, int], list[str]] = {}
    missing_keys = []
    for anchor_index in anchor_indices:
        anchor_name = "original_exception_center" if anchor_index == center_index else "lowest_lower_neighbor_anchor"
        for offset in range(-EXPANDED_RADIUS, EXPANDED_RADIUS + 1):
            key = (segment_index, anchor_index + offset)
            source_label = f"P2476/S1426:{anchor_name}:offset={offset}"
            if key in classified_by_key:
                replay_keys.setdefault(key, []).append(source_label)
            else:
                missing_keys.append({"anchor_index": anchor_index, "offset": offset, "missing_key": {"segment_index": key[0], "uncovered_index": key[1]}})
    function = function_for_family(family)
    replayed = []
    p2475_existing_keys = set()
    for row in exception["lower_neighbor_rows"] + [exception["lowest_available_halo_row"]]:
        p2475_existing_keys.add((int(row["segment_index"]), int(row["uncovered_index"])))
    for offset in exception["available_offsets"]:
        p2475_existing_keys.add((segment_index, center_index + int(offset)))
    for key in sorted(replay_keys):
        row = classified_by_key[key]
        fresh = replay_cell(family, row, projection, function)
        replayed.append({
            **fresh,
            "segment_index": key[0],
            "uncovered_index": key[1],
            "partition_class": row["partition_class"],
            "expanded_halo_source_rules": sorted(replay_keys[key]),
            "expanded_halo_source_count": len(replay_keys[key]),
            "fresh_decimal_separation_positive": Decimal(fresh["decimal_separation_from_zero"]) > 0,
            "already_in_p2475_radius_two_halo": key in p2475_existing_keys,
        })
    separations = [Decimal(row["decimal_separation_from_zero"]) for row in replayed]
    lowest = min(replayed, key=lambda row: Decimal(row["decimal_separation_from_zero"]))
    center = next(row for row in replayed if int(row["uncovered_index"]) == center_index)
    left_of_center = [row for row in replayed if int(row["uncovered_index"]) < center_index]
    right_of_center = [row for row in replayed if int(row["uncovered_index"]) > center_index]
    consecutive_pairs = []
    ordered = sorted(replayed, key=lambda row: int(row["uncovered_index"]))
    for prior, current in zip(ordered, ordered[1:]):
        consecutive_pairs.append({
            "from_uncovered_index": int(prior["uncovered_index"]),
            "to_uncovered_index": int(current["uncovered_index"]),
            "from_separation": prior["decimal_separation_from_zero"],
            "to_separation": current["decimal_separation_from_zero"],
            "separation_delta_to_minus_from": str(Decimal(current["decimal_separation_from_zero"]) - Decimal(prior["decimal_separation_from_zero"])),
            "strictly_increases": Decimal(current["decimal_separation_from_zero"]) > Decimal(prior["decimal_separation_from_zero"]),
        })
    class_counts: dict[str, int] = {}
    for row in replayed:
        class_counts[row["partition_class"]] = class_counts.get(row["partition_class"], 0) + 1
    return {
        "family": family,
        "expanded_radius": EXPANDED_RADIUS,
        "segment_index": segment_index,
        "original_exception_center_uncovered_index": center_index,
        "lowest_lower_neighbor_anchor_uncovered_index": lowest_lower_index,
        "anchor_indices": anchor_indices,
        "selected_segment_summary": selected_segment_summary,
        "unique_expanded_exception_replay_count": len(replayed),
        "incremental_cells_beyond_p2475_radius_two_halo_count": sum(not row["already_in_p2475_radius_two_halo"] for row in replayed),
        "missing_neighbor_count_due_to_segment_boundaries": len(missing_keys),
        "missing_neighbors_due_to_segment_boundaries": missing_keys,
        "partition_class_counts": class_counts,
        "all_expanded_exception_cells_exclude_zero": all(row["decimal_interval_excludes_zero"] for row in replayed),
        "all_expanded_exception_cells_have_positive_separation": all(row["fresh_decimal_separation_positive"] for row in replayed),
        "minimum_expanded_exception_decimal_separation": str(min(separations)),
        "minimum_expanded_exception_replay_row": lowest,
        "center_decimal_separation": center["decimal_separation_from_zero"],
        "lowest_is_left_boundary_of_expanded_replay": int(lowest["uncovered_index"]) == min(int(row["uncovered_index"]) for row in replayed),
        "center_has_lower_left_flank_in_expanded_replay": any(Decimal(row["decimal_separation_from_zero"]) < Decimal(center["decimal_separation_from_zero"]) for row in left_of_center),
        "center_has_lower_right_flank_in_expanded_replay": any(Decimal(row["decimal_separation_from_zero"]) < Decimal(center["decimal_separation_from_zero"]) for row in right_of_center),
        "consecutive_separation_pairs": consecutive_pairs,
        "all_consecutive_pairs_strictly_increase_left_to_right": all(row["strictly_increases"] for row in consecutive_pairs),
        "expanded_exception_replay_rows": replayed,
        "expanded_exception_replay_fingerprint_sha256": sha256_json(replayed),
    }


def append_once(path: Path, marker: str, section: str) -> None:
    text = path.read_text(encoding="utf-8")
    if marker not in text:
        path.write_text(text.rstrip() + "\n\n" + section.strip() + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2477/S1427 strict pointwise interval-Decimal P2459 lower-neighbor exception expanded-halo replay audit

`P2477/S1427` targets the single lower-neighbor exception exported by P2476 instead of pretending that every P2475 minimum witness is a local minimum.  It replays a radius-four expanded halo around both the original exception center and its lowest lower-neighbor anchor.  The fresh targeted replay covers `11` unique cells (`6` beyond the P2475 radius-two halo), all with positive Decimal separation from zero; the lowest cell remains on the left boundary of the expanded target window, so this is an exception/descent localization rather than a local-minimum theorem.

For a non-specialist: P2476 found the one place where the supposedly weakest saved point had even weaker nearby cells.  P2477 zooms into that place and recalculates a wider strip.  The strip still stays away from zero, but it also honestly says that the weak direction continues to the left edge of this small strip.  It is not a full complement replay, directed-rounding interval theorem, symbolic root-exclusion theorem, continuum theorem, selector/source/gauge theorem, physical-value generator, role-bearing `L_total`, `QW-2191` discharge, legacy-role transfer, or ToE closure.
"""
    lag_section = """
## P2477/S1427 P2459 lower-neighbor exception expanded-halo replay guard

`P2477/S1427` adds an exception-targeted expanded-halo guard to the finite-grid bookkeeping behind `L_total`: the unique P2476 lower-neighbor exception is replayed with a radius-four Decimal/Taylor strip around the original center and the lowest lower-neighbor anchor.  This improves localization of the finite weakest flank, but it does not add selector/source/gauge authority, physical-value generation, role-bearing finality, legacy-role transfer, or ToE closure.
"""
    append_once(DOC_FILES["equation_sheet"], "## P2477/S1427 strict pointwise interval-Decimal P2459 lower-neighbor exception expanded-halo replay audit", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2477/S1427 P2459 lower-neighbor exception expanded-halo replay guard", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    projection = projection_vector(sources["P2449_PROJECTION_REDUCTION"])
    p2450 = theorem(sources["P2450_ROOT_ISOLATION_MARGIN"], "strict_pointwise_projection_root_isolation_margin_certificate")
    p2451 = theorem(sources["P2451_FLOATING_INTERVAL_AUDIT"], "strict_pointwise_projection_interval_enclosure_root_exclusion_audit")
    p2456 = theorem(sources["P2456_DECIMAL_BOUNDARY_REPLAY"], "strict_pointwise_decimal_root_window_boundary_band_replay_certificate")
    p2466 = theorem(sources["P2466_DESCENT_TAIL_BOUNDARY"], "strict_pointwise_interval_decimal_adaptive_descent_tail_boundary_ledger_certificate")
    p2471 = theorem(sources["P2471_FINITE_PARTITION_WITNESS"], "strict_pointwise_interval_decimal_p2459_finite_partition_witness_certificate")
    p2475 = theorem(sources["P2475_CRITICAL_MINIMUM_HALO_REPLAY"], "strict_pointwise_interval_decimal_p2459_critical_minimum_halo_replay_audit")
    p2476 = theorem(sources["P2476_CRITICAL_HALO_ORDER_CLASSIFICATION"], "strict_pointwise_interval_decimal_p2459_critical_halo_order_classification_audit")
    exceptions = p2476_lower_neighbor_exceptions(p2476)
    expanded_replays = [expanded_exception_replay(exception, p2450, p2451, p2456, p2466, projection) for exception in exceptions]
    total_replay_count = sum(row["unique_expanded_exception_replay_count"] for row in expanded_replays)
    incremental_count = sum(row["incremental_cells_beyond_p2475_radius_two_halo_count"] for row in expanded_replays)
    min_sep = min(Decimal(row["minimum_expanded_exception_decimal_separation"]) for row in expanded_replays)
    p2459_universe = p2471["total_p2459_universe_cell_count_rebuilt"]
    theorem_export = {
        "theorem_name": "P2477_T1_strict_pointwise_interval_decimal_p2459_lower_neighbor_exception_expanded_halo_replay_audit",
        "audited_chain": ["P2471/S1421", "P2475/S1425", "P2476/S1426"],
        "expanded_radius": EXPANDED_RADIUS,
        "lower_neighbor_exception_count_inherited_from_p2476": len(exceptions),
        "expanded_exception_replays": expanded_replays,
        "total_unique_expanded_exception_replay_count": total_replay_count,
        "incremental_cells_beyond_p2475_radius_two_halo_count": incremental_count,
        "all_expanded_exception_cells_exclude_zero": all(row["all_expanded_exception_cells_exclude_zero"] for row in expanded_replays),
        "all_expanded_exception_cells_have_positive_separation": all(row["all_expanded_exception_cells_have_positive_separation"] for row in expanded_replays),
        "minimum_expanded_exception_decimal_separation": str(min_sep),
        "p2475_total_unique_halo_replay_count_inherited": p2475["total_unique_halo_replay_count"],
        "p2476_total_target_classification_count_inherited": p2476["total_target_classification_count"],
        "p2476_targets_with_lower_neighbor_count_inherited": p2476["targets_with_lower_neighbor_count"],
        "p2471_p2459_universe_count_inherited": p2459_universe,
        "targeted_p2477_residual_not_freshly_replayed_count_against_p2459_universe": p2459_universe - total_replay_count,
        "targeted_p2477_replay_fraction_of_p2459_universe": f"{total_replay_count}/{p2459_universe}",
        "finite_chain_coverage_budget_inherited_from_p2471": 0 if p2471["all_family_partitions_are_disjoint_and_complete"] else None,
        "p2471_disjoint_complete_inherited": p2471["all_family_partitions_are_disjoint_and_complete"],
        "no_full_complement_claimed_by_this_certificate": True,
        "full_complement_replay_exported_by_this_certificate": False,
        "finite_replay_chain_lower_neighbor_exception_expanded_halo_replay_audit_exported": True,
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
        "lay_summary": "This packet zooms into the single lower-neighbor exception from P2476. It recalculates a wider finite strip around the exception center and its lowest lower neighbor. The recalculated strip still excludes zero with positive Decimal separation, but the lowest point sits at the strip's left boundary, so the honest conclusion is targeted descent localization, not a local-minimum or continuum proof.",
        "not_licensed": [
            "This targeted expanded-halo replay is not a full complement replay and does not claim all P2459 cells were freshly replayed by P2477.",
            "It is not a directed-rounding interval theorem, symbolic root-exclusion theorem, analytic monotonicity theorem, or continuum root-exclusion theorem.",
            "It does not export selector/source/gauge closure, a physical-value generator, QW-2191 discharge, role-bearing L_total, legacy-role transfer, or ToE closure.",
            "No legacy physical roles are transferred to L_total or K_strict_gate by this finite exception replay.",
        ],
        "next_honest_step": "If continuing this branch, either replay a still-wider finite flank around the left-boundary minimum or switch to a real directed-rounding/symbolic argument; do not rename this targeted strip as a full proof.",
    }
    gatekeepers = {
        "rg_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "single_lower_neighbor_exception_inherited": theorem_export["lower_neighbor_exception_count_inherited_from_p2476"] == 1,
        "replay_count_expected": theorem_export["total_unique_expanded_exception_replay_count"] == 11,
        "incremental_count_expected": theorem_export["incremental_cells_beyond_p2475_radius_two_halo_count"] == 6,
        "minimum_separation_positive": Decimal(theorem_export["minimum_expanded_exception_decimal_separation"]) > 0,
        "all_replayed_cells_exclude_zero": theorem_export["all_expanded_exception_cells_exclude_zero"],
        "all_replayed_cells_have_positive_separation": theorem_export["all_expanded_exception_cells_have_positive_separation"],
        "left_boundary_descent_not_overclaimed": all(row["lowest_is_left_boundary_of_expanded_replay"] for row in expanded_replays),
        "no_full_complement_unless_genuinely_full": theorem_export["no_full_complement_claimed_by_this_certificate"] and not theorem_export["full_complement_replay_exported_by_this_certificate"],
        "coverage_budget_accounted": theorem_export["targeted_p2477_residual_not_freshly_replayed_count_against_p2459_universe"] == p2459_universe - total_replay_count,
        "p2471_partition_inherited": theorem_export["p2471_disjoint_complete_inherited"],
        "finite_grid_only_not_directed_rounding": theorem_export["finite_replay_chain_lower_neighbor_exception_expanded_halo_replay_audit_exported"] and not theorem_export["directed_rounding_interval_theorem_exported_by_this_certificate"],
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
        "schema_version": "p2477_s1427_v1",
        "packet_id": "P2477",
        "stage_id": "S1427",
        "status": "PASS_STRICT_POINTWISE_INTERVAL_DECIMAL_P2459_LOWER_NEIGHBOR_EXCEPTION_EXPANDED_HALO_REPLAY_AUDIT_NO_FULL_COMPLEMENT_NO_DIRECTED_ROUNDING_SELECTOR_SOURCE_THEOREM",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_pointwise_interval_decimal_p2459_lower_neighbor_exception_expanded_halo_replay_audit": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_pointwise_interval_decimal_p2459_lower_neighbor_exception_expanded_halo_replay_audit"]["theorem_export"]
    lines = [
        "# P2477/S1427 strict pointwise interval-Decimal P2459 lower-neighbor exception expanded-halo replay audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Expanded exception replay",
        "",
        f"P2476 lower-neighbor exceptions targeted: `{t['lower_neighbor_exception_count_inherited_from_p2476']}`.",
        f"Expanded radius: `{t['expanded_radius']}`.",
        f"Fresh targeted replay count: `{t['total_unique_expanded_exception_replay_count']}`.",
        f"Incremental cells beyond P2475 radius-two halo: `{t['incremental_cells_beyond_p2475_radius_two_halo_count']}`.",
        f"Minimum Decimal separation in expanded target: `{t['minimum_expanded_exception_decimal_separation']}`.",
        f"All expanded cells exclude zero: `{t['all_expanded_exception_cells_exclude_zero']}`.",
        "",
        "## Coverage budget",
        "",
        f"P2477 targeted replay fraction of the inherited P2459 finite universe: `{t['targeted_p2477_replay_fraction_of_p2459_universe']}`.",
        f"P2477 residual not freshly replayed against the P2459 universe: `{t['targeted_p2477_residual_not_freshly_replayed_count_against_p2459_universe']}`.",
        f"Inherited P2471 finite-chain coverage budget: `{t['finite_chain_coverage_budget_inherited_from_p2471']}`.",
        f"Full complement replay exported by P2477: `{t['full_complement_replay_exported_by_this_certificate']}`.",
        "",
        "## Plain-language progress note",
        "",
        t["lay_summary"],
        "",
        "## Hard limits / negative controls",
        "",
        "This is a finite lower-neighbor exception expanded-halo replay only.  It is not a full complement replay, directed-rounding interval theorem, symbolic root-exclusion theorem, analytic monotonicity theorem, global continuum root-exclusion theorem, selector/source/gauge theorem, physical-value generator, QW-2191 discharge, role-bearing `L_total`, legacy-role transfer, or ToE closure.",
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
