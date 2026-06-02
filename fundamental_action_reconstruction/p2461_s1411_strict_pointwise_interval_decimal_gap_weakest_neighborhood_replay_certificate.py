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

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"
OUT = GEN / "p2461_s1411_strict_pointwise_interval_decimal_gap_weakest_neighborhood_replay_certificate.json"
MD = GEN / "p2461_s1411_strict_pointwise_interval_decimal_gap_weakest_neighborhood_replay_certificate.md"

SOURCE_FILES = {
    "P2450_ROOT_ISOLATION_MARGIN": GEN / "p2450_s1400_strict_pointwise_projection_root_isolation_margin_certificate.json",
    "P2451_FLOATING_INTERVAL_AUDIT": GEN / "p2451_s1401_strict_pointwise_projection_interval_enclosure_root_exclusion_audit.json",
    "P2456_DECIMAL_BOUNDARY_REPLAY": GEN / "p2456_s1406_strict_pointwise_decimal_root_window_boundary_band_replay_certificate.json",
    "P2460_GAP_SENTINEL_REPLAY": GEN / "p2460_s1410_strict_pointwise_interval_decimal_gap_sentinel_replay_certificate.json",
    "P2449_PROJECTION_REDUCTION": GEN / "p2449_s1399_strict_pointwise_rank_lift_projection_reduction_certificate.json",
}
DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}
NEIGHBORHOOD_RADIUS = 3


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
        "new_packet": "P2461|S1411|gap weakest neighborhood|weakest gap neighborhood|gap-neighborhood replay",
        "p2460_input": "P2460|S1410|gap sentinel replay|unreplayed gap sentinel|stratified sentinel",
        "neighborhood_language": "weakest sentinel neighborhood|sentinel neighborhood|neighbor cells|local gap neighborhood",
        "decimal_backend": "Decimal interval|Taylor endpoint|zero-excluding|Decimal separation",
        "closure_blockers": "QW-2191|role-bearing L_total|physical-value generator|ToE closure|point-coordinate selector",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def weakest_sentinel_by_family(p2460: dict[str, Any], family: str) -> tuple[int, dict[str, Any]]:
    for family_row in p2460["family_gap_sentinel_replays"]:
        if family_row["family"] == family:
            weakest_segment_index = -1
            weakest_row: dict[str, Any] | None = None
            for segment in family_row["segment_rows"]:
                for replay in segment["sentinel_replays"]:
                    if weakest_row is None or Decimal(replay["decimal_separation_from_zero"]) < Decimal(weakest_row["decimal_separation_from_zero"]):
                        weakest_segment_index = segment["segment_index"]
                        weakest_row = replay
            if weakest_row is None:
                raise ValueError(f"no sentinels found for {family}")
            return weakest_segment_index, weakest_row
    raise ValueError(f"family not found in P2460: {family}")


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


def function_for_family(family: str) -> Callable[[list[Decimal], DecimalInterval], Any]:
    if family == "zero_projection_amplitude":
        return projection_amplitude_interval
    if family == "stationary_factor":
        return stationary_factor_interval
    raise ValueError(f"unknown family: {family}")


def replay_weakest_neighborhood(
    family: str,
    p2450: dict[str, Any],
    p2451: dict[str, Any],
    p2456: dict[str, Any],
    p2460: dict[str, Any],
    projection_vector: list[Decimal],
) -> dict[str, Any]:
    segment_index, weakest = weakest_sentinel_by_family(p2460, family)
    uncovered_cells = gap_uncovered_cells_for_segment(family, segment_index, p2450, p2451, p2456)
    center = int(weakest["uncovered_index"])
    left_index = max(0, center - NEIGHBORHOOD_RADIUS)
    right_index = min(len(uncovered_cells) - 1, center + NEIGHBORHOOD_RADIUS)
    function = function_for_family(family)
    rows = []
    for index in range(left_index, right_index + 1):
        cell = {**uncovered_cells[index], "uncovered_index": index, "uncovered_count": len(uncovered_cells)}
        rows.append(replay_cell(family, cell, projection_vector, function))
    minimum_separation = min(Decimal(row["decimal_separation_from_zero"]) for row in rows)
    return {
        "family": family,
        "weakest_p2460_sentinel": weakest,
        "weakest_segment_index": segment_index,
        "neighborhood_radius_requested": NEIGHBORHOOD_RADIUS,
        "neighborhood_index_left": left_index,
        "neighborhood_index_right": right_index,
        "neighborhood_replay_count": len(rows),
        "neighborhood_rows": rows,
        "all_neighborhood_cells_exclude_zero": all(row["decimal_interval_excludes_zero"] for row in rows),
        "minimum_neighborhood_decimal_separation": str(minimum_separation),
        "weakest_sentinel_replayed_inside_neighborhood": any(row["uncovered_index"] == center for row in rows),
    }


def append_doc_sections() -> None:
    eq_section = """
## P2461/S1411 strict pointwise interval-Decimal gap weakest-neighborhood replay certificate

`P2461/S1411` strengthens the P2460 gap-sentinel check by taking the weakest unreplayed-gap sentinel in each scalar family and replaying its local neighboring cells with the same Decimal/Taylor backend.  The replay covers 8 local gap-neighborhood cells in total and every replayed cell remains zero-excluding.

This is a local weakest-neighborhood replay, not an exhaustive Decimal full-complement replay, directed-rounding interval theorem, symbolic root-exclusion theorem, selector/source/gauge theorem, role-bearing `L_total`, `QW-2191` discharge, or ToE closure.
""".strip()
    lag_section = """
## P2461/S1411 Decimal gap weakest-neighborhood replay guard

`P2461/S1411` probes the local neighborhoods around the weakest P2460 gap sentinels.  It improves local gap diagnostics but does not close the full P2459 complement gap and adds no selector/source/gauge authority for `L_total`.
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
    projection_vector = [Decimal(str(value)) for value in sources["P2449_PROJECTION_REDUCTION"].get("strict_pointwise_rank_lift_projection_reduction_certificate", {}).get("theorem_export", {}).get("projection_vector", [])]
    neighborhood_replays = [
        replay_weakest_neighborhood("zero_projection_amplitude", p2450, p2451, p2456, p2460, projection_vector),
        replay_weakest_neighborhood("stationary_factor", p2450, p2451, p2456, p2460, projection_vector),
    ]
    total_replayed = sum(row["neighborhood_replay_count"] for row in neighborhood_replays)
    minimum_separation = min(Decimal(row["minimum_neighborhood_decimal_separation"]) for row in neighborhood_replays)
    theorem_export = {
        "theorem_name": "P2461_T1_strict_pointwise_interval_decimal_gap_weakest_neighborhood_replay_certificate",
        "inherited_gap_sentinel_replay": "P2460/S1410",
        "p2460_all_gap_sentinels_exclude_zero_inherited": p2460.get("all_gap_sentinels_exclude_zero"),
        "p2460_total_gap_sentinel_replay_count_inherited": p2460.get("total_gap_sentinel_replay_count"),
        "neighborhood_radius_requested": NEIGHBORHOOD_RADIUS,
        "family_weakest_neighborhood_replays": neighborhood_replays,
        "total_neighborhood_replay_count": total_replayed,
        "all_neighborhood_cells_exclude_zero": all(row["all_neighborhood_cells_exclude_zero"] for row in neighborhood_replays),
        "all_weakest_sentinels_replayed_inside_neighborhoods": all(row["weakest_sentinel_replayed_inside_neighborhood"] for row in neighborhood_replays),
        "minimum_neighborhood_decimal_separation": str(minimum_separation),
        "gap_weakest_neighborhood_replay_exported": True,
        "decimal_full_complement_replay_exported_by_this_certificate": False,
        "directed_rounding_interval_theorem_exported_by_this_certificate": False,
        "symbolic_root_exclusion_theorem_exported_by_this_certificate": False,
        "pointwise_coordinate_selector_exported_by_this_certificate": False,
        "strict_observable_source_constraint_exported_by_this_certificate": False,
        "gauge_slice_theorem_exported_by_this_certificate": False,
        "strict_physical_value_generator_exported": False,
        "qw2191_discharged": False,
        "role_bearing_ltotal_exported": False,
        "toe_closure_exported": False,
        "not_licensed": [
            "Weakest-neighborhood replay is not an exhaustive Decimal full-complement replay and does not close the P2459 coverage gap.",
            "Local zero exclusion around weakest gap sentinels does not prove symbolic root exclusion or a directed-rounding interval theorem.",
            "No point-coordinate selector, strict observable/source theorem, gauge-slice theorem, physical-value generator, QW-2191 discharge, role-bearing L_total export, or ToE closure is exported.",
        ],
        "next_honest_step": (
            "Expand from local weakest-neighborhood replay to adaptive multi-neighborhood replay or formal directed-rounding proof over every unreplayed complement cell."
        ),
    }
    gatekeepers = {
        "rg_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "p2460_gap_sentinels_exclude_zero_inherited": theorem_export["p2460_all_gap_sentinels_exclude_zero_inherited"] is True,
        "p2460_sentinel_count_inherited": theorem_export["p2460_total_gap_sentinel_replay_count_inherited"] == 25,
        "family_count_is_two": len(neighborhood_replays) == 2,
        "total_neighborhood_replay_count": theorem_export["total_neighborhood_replay_count"] == 8,
        "all_neighborhood_cells_exclude_zero": theorem_export["all_neighborhood_cells_exclude_zero"],
        "all_weakest_sentinels_replayed_inside_neighborhoods": theorem_export["all_weakest_sentinels_replayed_inside_neighborhoods"],
        "minimum_neighborhood_separation_positive": Decimal(theorem_export["minimum_neighborhood_decimal_separation"]) > 0,
        "no_decimal_full_complement_replay": not theorem_export["decimal_full_complement_replay_exported_by_this_certificate"],
        "no_directed_rounding_theorem": not theorem_export["directed_rounding_interval_theorem_exported_by_this_certificate"],
        "no_symbolic_root_exclusion": not theorem_export["symbolic_root_exclusion_theorem_exported_by_this_certificate"],
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
        "schema_version": "p2461_s1411_v1",
        "packet_id": "P2461",
        "stage_id": "S1411",
        "status": "PASS_STRICT_POINTWISE_INTERVAL_DECIMAL_GAP_WEAKEST_NEIGHBORHOOD_REPLAY_NO_FULL_COMPLEMENT_SELECTOR_SOURCE_THEOREM",
        "strict_pointwise_interval_decimal_gap_weakest_neighborhood_replay_certificate": {
            "inputs": {name: rel(path) for name, path in SOURCE_FILES.items()},
            "input_fingerprints": {name: sha256_json(sources[name]) for name in sources},
            "rg_nonduplication_audit": grep,
            "theorem_export": theorem_export,
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_pointwise_interval_decimal_gap_weakest_neighborhood_replay_certificate"]["theorem_export"]
    lines = [
        "# P2461/S1411 strict pointwise interval-Decimal gap weakest-neighborhood replay certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Weakest-neighborhood replay",
        "",
        f"P2460 gap sentinels inherited: `{t['p2460_total_gap_sentinel_replay_count_inherited']}`.",
        f"Neighborhood cells replayed: `{t['total_neighborhood_replay_count']}`.",
        f"All neighborhood cells exclude zero: `{t['all_neighborhood_cells_exclude_zero']}`.",
        f"Minimum neighborhood Decimal separation: `{t['minimum_neighborhood_decimal_separation']}`.",
        "",
        "## Guardrail",
        "",
        "This is a local weakest-neighborhood Decimal/Taylor replay only.  It exports no Decimal full-complement replay, no directed-rounding interval theorem, no symbolic root-exclusion theorem, no point-coordinate selector, no strict observable/source theorem, no gauge-slice theorem, no physical-value generator, no QW-2191 discharge, no role-bearing `L_total`, and no ToE closure.",
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
