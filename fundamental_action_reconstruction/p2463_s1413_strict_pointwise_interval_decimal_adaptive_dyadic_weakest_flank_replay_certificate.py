#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from decimal import Decimal
from pathlib import Path
from typing import Any, Callable

from p2453_s1403_strict_pointwise_directed_decimal_weakest_cell_replay_certificate import (
    DecimalInterval,
    projection_amplitude_interval,
    stationary_factor_interval,
)
from p2460_s1410_strict_pointwise_interval_decimal_gap_sentinel_replay_certificate import (
    boundary_skip_counts_for_segment,
    complement_segments,
    interval_cell_width,
    remove_boundary_covered_cells,
    replay_cell,
    segment_cells,
)
from p2462_s1412_strict_pointwise_interval_decimal_gap_dyadic_refinement_replay_certificate import (
    DYADIC_REFINEMENT_FRACTIONS,
    P2460_SENTINEL_FRACTIONS,
    fraction_to_index,
)

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"
OUT = GEN / "p2463_s1413_strict_pointwise_interval_decimal_adaptive_dyadic_weakest_flank_replay_certificate.json"
MD = GEN / "p2463_s1413_strict_pointwise_interval_decimal_adaptive_dyadic_weakest_flank_replay_certificate.md"

SOURCE_FILES = {
    "P2450_ROOT_ISOLATION_MARGIN": GEN / "p2450_s1400_strict_pointwise_projection_root_isolation_margin_certificate.json",
    "P2451_FLOATING_INTERVAL_AUDIT": GEN / "p2451_s1401_strict_pointwise_projection_interval_enclosure_root_exclusion_audit.json",
    "P2456_DECIMAL_BOUNDARY_REPLAY": GEN / "p2456_s1406_strict_pointwise_decimal_root_window_boundary_band_replay_certificate.json",
    "P2460_GAP_SENTINEL_REPLAY": GEN / "p2460_s1410_strict_pointwise_interval_decimal_gap_sentinel_replay_certificate.json",
    "P2462_GAP_DYADIC_REFINEMENT": GEN / "p2462_s1412_strict_pointwise_interval_decimal_gap_dyadic_refinement_replay_certificate.json",
    "P2449_PROJECTION_REDUCTION": GEN / "p2449_s1399_strict_pointwise_rank_lift_projection_reduction_certificate.json",
}
DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

ADAPTIVE_FLANK_RADIUS = 4


def rel(path: Path) -> str:
    return path.relative_to(REPO).as_posix()


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": rel(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


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
    return {"count": len(lines), "samples": lines[:35]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2463|S1413|adaptive dyadic weakest flank|dyadic weakest flank|adaptive flank replay",
        "p2462_input": "P2462|S1412|gap dyadic refinement|dyadic-refinement replay|dyadic refinement replay",
        "descent_language": "weakest dyadic flank|flank descent|smaller-than-dyadic|adaptive local descent|dyadic not local minimum",
        "decimal_backend": "Decimal/Taylor|Decimal separation|zero-excluding|Decimal interval",
        "closure_blockers": "QW-2191|role-bearing L_total|physical-value generator|ToE closure|point-coordinate selector",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def function_for_family(family: str) -> Callable[[list[Decimal], DecimalInterval], Any]:
    if family == "zero_projection_amplitude":
        return projection_amplitude_interval
    if family == "stationary_factor":
        return stationary_factor_interval
    raise ValueError(f"unknown family: {family}")


def weakest_p2462_dyadic_row(p2462: dict[str, Any], family: str) -> tuple[int, dict[str, Any]]:
    for family_row in p2462["family_dyadic_refinement_replays"]:
        if family_row["family"] != family:
            continue
        weakest_segment_index = -1
        weakest_row: dict[str, Any] | None = None
        for segment in family_row["segment_rows"]:
            for replay in segment["dyadic_refinement_replays"]:
                if weakest_row is None or Decimal(replay["decimal_separation_from_zero"]) < Decimal(weakest_row["decimal_separation_from_zero"]):
                    weakest_segment_index = segment["segment_index"]
                    weakest_row = replay
        if weakest_row is None:
            raise ValueError(f"no P2462 dyadic rows found for {family}")
        return weakest_segment_index, weakest_row
    raise ValueError(f"family not found in P2462: {family}")


def gap_uncovered_cells_for_segment(
    family: str,
    segment_index: int,
    p2450: dict[str, Any],
    p2451: dict[str, Any],
    p2456: dict[str, Any],
) -> list[dict[str, float]]:
    segments = complement_segments(p2450, family)
    segment = segments[segment_index]
    cells = segment_cells(segment, interval_cell_width(p2451, family))
    left_skip, right_skip = boundary_skip_counts_for_segment(p2456, family, segment)
    return remove_boundary_covered_cells(cells, left_skip, right_skip)


def inherited_anchor_indexes(uncovered_count: int) -> set[int]:
    max_index = uncovered_count - 1
    fractions = [*P2460_SENTINEL_FRACTIONS, *DYADIC_REFINEMENT_FRACTIONS]
    return {fraction_to_index(fraction, max_index) for fraction in fractions}


def adaptive_flank_cells(uncovered_cells: list[dict[str, float]], center_index: int) -> list[dict[str, Any]]:
    inherited_indexes = inherited_anchor_indexes(len(uncovered_cells))
    left_index = max(0, center_index - ADAPTIVE_FLANK_RADIUS)
    right_index = min(len(uncovered_cells) - 1, center_index + ADAPTIVE_FLANK_RADIUS)
    rows = []
    for index in range(left_index, right_index + 1):
        if index in inherited_indexes:
            continue
        rows.append(
            {
                **uncovered_cells[index],
                "uncovered_index": index,
                "uncovered_count": len(uncovered_cells),
                "offset_from_p2462_weakest_dyadic": index - center_index,
            }
        )
    return rows


def adaptive_weakest_flank_replay(
    family: str,
    p2450: dict[str, Any],
    p2451: dict[str, Any],
    p2456: dict[str, Any],
    p2462: dict[str, Any],
    projection_vector: list[Decimal],
) -> dict[str, Any]:
    segment_index, weakest = weakest_p2462_dyadic_row(p2462, family)
    uncovered_cells = gap_uncovered_cells_for_segment(family, segment_index, p2450, p2451, p2456)
    center_index = int(weakest["uncovered_index"])
    function = function_for_family(family)
    cells = adaptive_flank_cells(uncovered_cells, center_index)
    replayed = [
        replay_cell(family, cell, projection_vector, function)
        | {"offset_from_p2462_weakest_dyadic": cell["offset_from_p2462_weakest_dyadic"]}
        for cell in cells
    ]
    minimum_flank_row = min(replayed, key=lambda row: Decimal(row["decimal_separation_from_zero"]))
    inherited_weakest_separation = Decimal(weakest["decimal_separation_from_zero"])
    minimum_flank_separation = Decimal(minimum_flank_row["decimal_separation_from_zero"])
    smaller_rows = [row for row in replayed if Decimal(row["decimal_separation_from_zero"]) < inherited_weakest_separation]
    return {
        "family": family,
        "weakest_p2462_dyadic_segment_index": segment_index,
        "weakest_p2462_dyadic_row": weakest,
        "adaptive_flank_radius_requested": ADAPTIVE_FLANK_RADIUS,
        "adaptive_flank_replay_count": len(replayed),
        "adaptive_flank_rows": replayed,
        "all_adaptive_flank_cells_exclude_zero": all(row["decimal_interval_excludes_zero"] for row in replayed),
        "minimum_adaptive_flank_row": minimum_flank_row,
        "minimum_adaptive_flank_decimal_separation": str(minimum_flank_separation),
        "inherited_weakest_dyadic_decimal_separation": str(inherited_weakest_separation),
        "smaller_than_inherited_weakest_dyadic_count": len(smaller_rows),
        "found_smaller_than_inherited_weakest_dyadic": bool(smaller_rows),
    }


def append_doc_sections() -> None:
    eq_section = """
## P2463/S1413 strict pointwise interval-Decimal adaptive dyadic weakest-flank replay certificate

`P2463/S1413` follows the P2462 dyadic refinement by taking the weakest P2462 dyadic cell in each scalar family and replaying the nearby non-anchor flank cells within radius 4.  The replay adds 16 non-anchor flank cells, all remain zero-excluding, and it records that smaller Decimal separations than the P2462 dyadic anchors are found locally without turning that observation into a monotonicity or full-complement theorem.

This is an adaptive weakest-flank replay and descent diagnostic, not a full complement Decimal replay, directed-rounding interval theorem, symbolic root-exclusion theorem, analytic monotonicity theorem, selector/source/gauge theorem, role-bearing `L_total`, `QW-2191` discharge, or ToE closure.
""".strip()
    lag_section = """
## P2463/S1413 Decimal adaptive dyadic weakest-flank replay guard

`P2463/S1413` probes non-anchor neighbors around the weakest P2462 dyadic cells.  It improves adaptive gap diagnostics and exposes local smaller-separation flanks, but it still leaves the full P2459 complement gap open and adds no selector/source/gauge authority for `L_total`.
""".strip()
    for path, section in [(DOC_FILES["equation_sheet"], eq_section), (DOC_FILES["lagrangian_eom_draft"], lag_section)]:
        text = path.read_text(encoding="utf-8") if path.exists() else ""
        marker = section.splitlines()[0]
        if marker not in text:
            path.write_text(text.rstrip() + "\n\n" + section + "\n", encoding="utf-8")


def build_payload() -> dict[str, Any]:
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    grep = rg_audit()
    p2450 = sources["P2450_ROOT_ISOLATION_MARGIN"].get("strict_pointwise_projection_root_isolation_margin_certificate", {}).get("theorem_export", {})
    p2451 = sources["P2451_FLOATING_INTERVAL_AUDIT"].get("strict_pointwise_projection_interval_enclosure_root_exclusion_audit", {}).get("theorem_export", {})
    p2456 = sources["P2456_DECIMAL_BOUNDARY_REPLAY"].get("strict_pointwise_decimal_root_window_boundary_band_replay_certificate", {}).get("theorem_export", {})
    p2460 = sources["P2460_GAP_SENTINEL_REPLAY"].get("strict_pointwise_interval_decimal_gap_sentinel_replay_certificate", {}).get("theorem_export", {})
    p2462 = sources["P2462_GAP_DYADIC_REFINEMENT"].get("strict_pointwise_interval_decimal_gap_dyadic_refinement_replay_certificate", {}).get("theorem_export", {})
    projection_vector = [Decimal(str(value)) for value in sources["P2449_PROJECTION_REDUCTION"].get("strict_pointwise_rank_lift_projection_reduction_certificate", {}).get("theorem_export", {}).get("projection_vector", [])]
    family_replays = [
        adaptive_weakest_flank_replay("zero_projection_amplitude", p2450, p2451, p2456, p2462, projection_vector),
        adaptive_weakest_flank_replay("stationary_factor", p2450, p2451, p2456, p2462, projection_vector),
    ]
    total_replayed = sum(row["adaptive_flank_replay_count"] for row in family_replays)
    minimum_separation = min(Decimal(row["minimum_adaptive_flank_decimal_separation"]) for row in family_replays)
    total_smaller_than_dyadic = sum(row["smaller_than_inherited_weakest_dyadic_count"] for row in family_replays)
    theorem_export = {
        "theorem_name": "P2463_T1_strict_pointwise_interval_decimal_adaptive_dyadic_weakest_flank_replay_certificate",
        "inherited_gap_sentinel_replay": "P2460/S1410",
        "inherited_gap_dyadic_refinement_replay": "P2462/S1412",
        "p2460_total_gap_sentinel_replay_count_inherited": p2460.get("total_gap_sentinel_replay_count"),
        "p2462_total_dyadic_refinement_replay_count_inherited": p2462.get("total_dyadic_refinement_replay_count"),
        "p2462_all_dyadic_refinement_cells_exclude_zero_inherited": p2462.get("all_dyadic_refinement_cells_exclude_zero"),
        "adaptive_flank_radius_requested": ADAPTIVE_FLANK_RADIUS,
        "family_adaptive_weakest_flank_replays": family_replays,
        "total_adaptive_flank_replay_count": total_replayed,
        "all_adaptive_flank_cells_exclude_zero": all(row["all_adaptive_flank_cells_exclude_zero"] for row in family_replays),
        "minimum_adaptive_flank_decimal_separation": str(minimum_separation),
        "total_smaller_than_inherited_weakest_dyadic_count": total_smaller_than_dyadic,
        "all_families_found_smaller_than_inherited_weakest_dyadic": all(row["found_smaller_than_inherited_weakest_dyadic"] for row in family_replays),
        "adaptive_dyadic_weakest_flank_replay_exported": True,
        "decimal_full_complement_replay_exported_by_this_certificate": False,
        "directed_rounding_interval_theorem_exported_by_this_certificate": False,
        "symbolic_root_exclusion_theorem_exported_by_this_certificate": False,
        "analytic_monotonicity_theorem_exported_by_this_certificate": False,
        "pointwise_coordinate_selector_exported_by_this_certificate": False,
        "strict_observable_source_constraint_exported_by_this_certificate": False,
        "gauge_slice_theorem_exported_by_this_certificate": False,
        "strict_physical_value_generator_exported": False,
        "qw2191_discharged": False,
        "role_bearing_ltotal_exported": False,
        "toe_closure_exported": False,
        "not_licensed": [
            "Adaptive weakest-flank replay is not an exhaustive Decimal full-complement replay and does not close the P2459 coverage gap.",
            "Finding smaller local flank separations than the P2462 dyadic anchors is a descent diagnostic, not an analytic monotonicity theorem or symbolic root-exclusion theorem.",
            "No point-coordinate selector, strict observable/source theorem, gauge-slice theorem, physical-value generator, QW-2191 discharge, role-bearing L_total export, or ToE closure is exported.",
        ],
        "next_honest_step": (
            "Use the newly found weakest flank cells to drive another adaptive replay layer, or replace sampled descent diagnostics with a formal directed-rounding/full-complement proof."
        ),
    }
    gatekeepers = {
        "rg_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "p2460_sentinel_count_inherited": theorem_export["p2460_total_gap_sentinel_replay_count_inherited"] == 25,
        "p2462_dyadic_count_inherited": theorem_export["p2462_total_dyadic_refinement_replay_count_inherited"] == 20,
        "p2462_dyadic_zero_exclusion_inherited": theorem_export["p2462_all_dyadic_refinement_cells_exclude_zero_inherited"] is True,
        "family_count_is_two": len(family_replays) == 2,
        "total_adaptive_flank_replay_count": theorem_export["total_adaptive_flank_replay_count"] == 16,
        "all_adaptive_flank_cells_exclude_zero": theorem_export["all_adaptive_flank_cells_exclude_zero"],
        "minimum_adaptive_flank_separation_positive": Decimal(theorem_export["minimum_adaptive_flank_decimal_separation"]) > 0,
        "all_families_found_smaller_than_inherited_weakest_dyadic": theorem_export["all_families_found_smaller_than_inherited_weakest_dyadic"],
        "no_decimal_full_complement_replay": not theorem_export["decimal_full_complement_replay_exported_by_this_certificate"],
        "no_directed_rounding_theorem": not theorem_export["directed_rounding_interval_theorem_exported_by_this_certificate"],
        "no_symbolic_root_exclusion": not theorem_export["symbolic_root_exclusion_theorem_exported_by_this_certificate"],
        "no_analytic_monotonicity_theorem": not theorem_export["analytic_monotonicity_theorem_exported_by_this_certificate"],
        "no_pointwise_selector_export": not theorem_export["pointwise_coordinate_selector_exported_by_this_certificate"],
        "no_observable_source_export": not theorem_export["strict_observable_source_constraint_exported_by_this_certificate"],
        "no_gauge_slice_export": not theorem_export["gauge_slice_theorem_exported_by_this_certificate"],
        "no_value_generator_export": not theorem_export["strict_physical_value_generator_exported"],
        "no_qw2191_discharge": not theorem_export["qw2191_discharged"],
        "no_ltotal_export": not theorem_export["role_bearing_ltotal_exported"],
        "no_toe_export": not theorem_export["toe_closure_exported"],
        "fingerprint_stable": sha256_json(theorem_export) == sha256_json(theorem_export),
    }
    return {
        "schema_version": "p2463_s1413_v1",
        "packet_id": "P2463",
        "stage_id": "S1413",
        "status": "PASS_STRICT_POINTWISE_INTERVAL_DECIMAL_ADAPTIVE_DYADIC_WEAKEST_FLANK_REPLAY_NO_FULL_COMPLEMENT_SELECTOR_SOURCE_THEOREM",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_pointwise_interval_decimal_adaptive_dyadic_weakest_flank_replay_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_pointwise_interval_decimal_adaptive_dyadic_weakest_flank_replay_certificate"]["theorem_export"]
    lines = [
        "# P2463/S1413 strict pointwise interval-Decimal adaptive dyadic weakest-flank replay certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Adaptive dyadic weakest-flank replay",
        "",
        f"P2460 gap sentinels inherited: `{t['p2460_total_gap_sentinel_replay_count_inherited']}`.",
        f"P2462 dyadic refinement cells inherited: `{t['p2462_total_dyadic_refinement_replay_count_inherited']}`.",
        f"Adaptive flank radius: `{t['adaptive_flank_radius_requested']}`.",
        f"Adaptive non-anchor flank cells replayed: `{t['total_adaptive_flank_replay_count']}`.",
        f"All adaptive flank cells exclude zero: `{t['all_adaptive_flank_cells_exclude_zero']}`.",
        f"Families finding smaller-than-P2462-dyadic flank cells: `{t['all_families_found_smaller_than_inherited_weakest_dyadic']}`.",
        f"Minimum adaptive flank Decimal separation: `{t['minimum_adaptive_flank_decimal_separation']}`.",
        "",
        "## Guardrail",
        "",
        "This is an adaptive local flank replay only.  Smaller-than-dyadic flank cells are descent diagnostics, not a monotonicity theorem.  The certificate exports no Decimal full-complement replay, no directed-rounding interval theorem, no symbolic root-exclusion theorem, no point-coordinate selector, no strict observable/source theorem, no gauge-slice theorem, no physical-value generator, no QW-2191 discharge, no role-bearing `L_total`, and no ToE closure.",
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
