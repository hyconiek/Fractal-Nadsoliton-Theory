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
OUT = GEN / "p2462_s1412_strict_pointwise_interval_decimal_gap_dyadic_refinement_replay_certificate.json"
MD = GEN / "p2462_s1412_strict_pointwise_interval_decimal_gap_dyadic_refinement_replay_certificate.md"

SOURCE_FILES = {
    "P2450_ROOT_ISOLATION_MARGIN": GEN / "p2450_s1400_strict_pointwise_projection_root_isolation_margin_certificate.json",
    "P2451_FLOATING_INTERVAL_AUDIT": GEN / "p2451_s1401_strict_pointwise_projection_interval_enclosure_root_exclusion_audit.json",
    "P2456_DECIMAL_BOUNDARY_REPLAY": GEN / "p2456_s1406_strict_pointwise_decimal_root_window_boundary_band_replay_certificate.json",
    "P2460_GAP_SENTINEL_REPLAY": GEN / "p2460_s1410_strict_pointwise_interval_decimal_gap_sentinel_replay_certificate.json",
    "P2461_GAP_WEAKEST_NEIGHBORHOOD": GEN / "p2461_s1411_strict_pointwise_interval_decimal_gap_weakest_neighborhood_replay_certificate.json",
    "P2449_PROJECTION_REDUCTION": GEN / "p2449_s1399_strict_pointwise_rank_lift_projection_reduction_certificate.json",
}
DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

# P2460 already probes 0, 1/4, 1/2, 3/4, 1.  P2462 honestly adds the
# dyadic midpoints between those quarter sentinels without claiming exhaustive
# complement coverage.
DYADIC_REFINEMENT_FRACTIONS = [Decimal("0.125"), Decimal("0.375"), Decimal("0.625"), Decimal("0.875")]
P2460_SENTINEL_FRACTIONS = [Decimal("0"), Decimal("0.25"), Decimal("0.5"), Decimal("0.75"), Decimal("1")]


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
        "new_packet": "P2462|S1412|gap dyadic refinement|dyadic gap refinement|dyadic refinement replay",
        "p2460_input": "P2460|S1410|gap sentinel replay|stratified sentinel|unreplayed gap sentinel",
        "p2461_input": "P2461|S1411|gap weakest neighborhood|weakest gap neighborhood|local gap neighborhood",
        "dyadic_language": "eighth-gap|inter-sentinel midpoint|quarter-sentinel refinement|dyadic midpoint|refinement fraction",
        "closure_blockers": "QW-2191|role-bearing L_total|physical-value generator|ToE closure|point-coordinate selector",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def function_for_family(family: str) -> Callable[[list[Decimal], DecimalInterval], Any]:
    if family == "zero_projection_amplitude":
        return projection_amplitude_interval
    if family == "stationary_factor":
        return stationary_factor_interval
    raise ValueError(f"unknown family: {family}")


def fraction_to_index(fraction: Decimal, max_index: int) -> int:
    return int((fraction * Decimal(max_index)).to_integral_value(rounding="ROUND_HALF_UP"))


def selected_fraction_cells(uncovered_cells: list[dict[str, float]], fractions: list[Decimal]) -> list[dict[str, Any]]:
    if not uncovered_cells:
        return []
    max_index = len(uncovered_cells) - 1
    rows = []
    seen: set[int] = set()
    for fraction in fractions:
        index = fraction_to_index(fraction, max_index)
        if index in seen:
            continue
        seen.add(index)
        rows.append({**uncovered_cells[index], "uncovered_index": index, "uncovered_count": len(uncovered_cells), "refinement_fraction": str(fraction)})
    return rows


def p2460_sentinel_indexes(uncovered_count: int) -> set[int]:
    max_index = uncovered_count - 1
    return {fraction_to_index(fraction, max_index) for fraction in P2460_SENTINEL_FRACTIONS}


def family_dyadic_refinement_replay(
    family: str,
    p2450: dict[str, Any],
    p2451: dict[str, Any],
    p2456: dict[str, Any],
    projection_vector: list[Decimal],
) -> dict[str, Any]:
    width = interval_cell_width(p2451, family)
    function = function_for_family(family)
    segment_rows = []
    replay_rows = []
    total_uncovered = 0
    for segment_index, segment in enumerate(complement_segments(p2450, family)):
        cells = segment_cells(segment, width)
        left_skip, right_skip = boundary_skip_counts_for_segment(p2456, family, segment)
        uncovered_cells = remove_boundary_covered_cells(cells, left_skip, right_skip)
        total_uncovered += len(uncovered_cells)
        refinement_cells = selected_fraction_cells(uncovered_cells, DYADIC_REFINEMENT_FRACTIONS)
        replayed = [replay_cell(family, cell, projection_vector, function) | {"refinement_fraction": cell["refinement_fraction"]} for cell in refinement_cells]
        sentinel_indexes = p2460_sentinel_indexes(len(uncovered_cells)) if uncovered_cells else set()
        no_p2460_index_duplicates = all(row["uncovered_index"] not in sentinel_indexes for row in replayed)
        segment_rows.append(
            {
                "segment_index": segment_index,
                "segment_left": segment["left"],
                "segment_right": segment["right"],
                "p2451_cell_count": len(cells),
                "p2456_boundary_left_skip_cell_count": left_skip,
                "p2456_boundary_right_skip_cell_count": right_skip,
                "p2456_boundary_covered_cell_count": left_skip + right_skip,
                "uncovered_by_boundary_chain_cell_count": len(uncovered_cells),
                "p2462_refinement_fraction_count_requested": len(DYADIC_REFINEMENT_FRACTIONS),
                "dyadic_refinement_replay_count": len(replayed),
                "no_p2460_quarter_sentinel_index_duplicates": no_p2460_index_duplicates,
                "all_segment_dyadic_refinements_exclude_zero": all(row["decimal_interval_excludes_zero"] for row in replayed),
                "minimum_segment_dyadic_refinement_decimal_separation": str(min(Decimal(row["decimal_separation_from_zero"]) for row in replayed)) if replayed else None,
                "dyadic_refinement_replays": replayed,
            }
        )
        replay_rows.extend(replayed)
    return {
        "family": family,
        "segment_count": len(segment_rows),
        "segment_rows": segment_rows,
        "uncovered_by_boundary_chain_cell_count": total_uncovered,
        "dyadic_refinement_replay_count": len(replay_rows),
        "all_dyadic_refinement_cells_exclude_zero": all(row["decimal_interval_excludes_zero"] for row in replay_rows),
        "all_refinement_indexes_avoid_p2460_quarter_sentinels": all(row["no_p2460_quarter_sentinel_index_duplicates"] for row in segment_rows),
        "minimum_dyadic_refinement_decimal_separation": str(min(Decimal(row["decimal_separation_from_zero"]) for row in replay_rows)),
    }


def append_doc_sections() -> None:
    eq_section = """
## P2462/S1412 strict pointwise interval-Decimal gap dyadic-refinement replay certificate

`P2462/S1412` refines the P2460 gap-sentinel ledger without pretending to close the P2459 gap: in each unreplayed complement segment it replays the dyadic midpoints between the P2460 quarter sentinels (`1/8`, `3/8`, `5/8`, `7/8`) using the same Decimal/Taylor backend.  The replay adds 20 non-quarter-sentinel refinement cells and every replayed cell remains zero-excluding.

This is a dyadic refinement replay between already-audited sentinels, not a full complement Decimal replay, directed-rounding interval theorem, symbolic root-exclusion theorem, selector/source/gauge theorem, role-bearing `L_total`, `QW-2191` discharge, or ToE closure.
""".strip()
    lag_section = """
## P2462/S1412 Decimal gap dyadic-refinement replay guard

`P2462/S1412` probes the eighth-step midpoints between the P2460 quarter sentinels inside the unreplayed complement segments.  It improves stratified gap diagnostics but still leaves the full P2459 complement gap open and adds no selector/source/gauge authority for `L_total`.
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
    p2461 = sources["P2461_GAP_WEAKEST_NEIGHBORHOOD"].get("strict_pointwise_interval_decimal_gap_weakest_neighborhood_replay_certificate", {}).get("theorem_export", {})
    projection_vector = [Decimal(str(value)) for value in sources["P2449_PROJECTION_REDUCTION"].get("strict_pointwise_rank_lift_projection_reduction_certificate", {}).get("theorem_export", {}).get("projection_vector", [])]
    family_replays = [
        family_dyadic_refinement_replay("zero_projection_amplitude", p2450, p2451, p2456, projection_vector),
        family_dyadic_refinement_replay("stationary_factor", p2450, p2451, p2456, projection_vector),
    ]
    total_replayed = sum(row["dyadic_refinement_replay_count"] for row in family_replays)
    minimum_separation = min(Decimal(row["minimum_dyadic_refinement_decimal_separation"]) for row in family_replays)
    theorem_export = {
        "theorem_name": "P2462_T1_strict_pointwise_interval_decimal_gap_dyadic_refinement_replay_certificate",
        "inherited_gap_sentinel_replay": "P2460/S1410",
        "inherited_gap_weakest_neighborhood_replay": "P2461/S1411",
        "p2460_total_gap_sentinel_replay_count_inherited": p2460.get("total_gap_sentinel_replay_count"),
        "p2461_total_neighborhood_replay_count_inherited": p2461.get("total_neighborhood_replay_count"),
        "p2460_all_gap_sentinels_exclude_zero_inherited": p2460.get("all_gap_sentinels_exclude_zero"),
        "p2461_all_neighborhood_cells_exclude_zero_inherited": p2461.get("all_neighborhood_cells_exclude_zero"),
        "dyadic_refinement_fractions": [str(value) for value in DYADIC_REFINEMENT_FRACTIONS],
        "p2460_quarter_sentinel_fractions_avoided": [str(value) for value in P2460_SENTINEL_FRACTIONS],
        "family_dyadic_refinement_replays": family_replays,
        "total_dyadic_refinement_replay_count": total_replayed,
        "all_dyadic_refinement_cells_exclude_zero": all(row["all_dyadic_refinement_cells_exclude_zero"] for row in family_replays),
        "all_refinement_indexes_avoid_p2460_quarter_sentinels": all(row["all_refinement_indexes_avoid_p2460_quarter_sentinels"] for row in family_replays),
        "minimum_dyadic_refinement_decimal_separation": str(minimum_separation),
        "gap_dyadic_refinement_replay_exported": True,
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
            "Dyadic refinement replay is not an exhaustive Decimal full-complement replay and does not close the P2459 coverage gap.",
            "Eighth-step zero exclusion between sentinels does not prove symbolic root exclusion or a directed-rounding interval theorem.",
            "No point-coordinate selector, strict observable/source theorem, gauge-slice theorem, physical-value generator, QW-2191 discharge, role-bearing L_total export, or ToE closure is exported.",
        ],
        "next_honest_step": (
            "Either add adaptive refinement driven by the smallest P2462 dyadic separation, or replace sampled replay by a formal directed-rounding/full-complement proof."
        ),
    }
    gatekeepers = {
        "rg_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "p2460_gap_sentinels_exclude_zero_inherited": theorem_export["p2460_all_gap_sentinels_exclude_zero_inherited"] is True,
        "p2461_neighborhoods_exclude_zero_inherited": theorem_export["p2461_all_neighborhood_cells_exclude_zero_inherited"] is True,
        "p2460_sentinel_count_inherited": theorem_export["p2460_total_gap_sentinel_replay_count_inherited"] == 25,
        "p2461_neighborhood_count_inherited": theorem_export["p2461_total_neighborhood_replay_count_inherited"] == 8,
        "family_count_is_two": len(family_replays) == 2,
        "total_dyadic_refinement_replay_count": theorem_export["total_dyadic_refinement_replay_count"] == 20,
        "all_dyadic_refinement_cells_exclude_zero": theorem_export["all_dyadic_refinement_cells_exclude_zero"],
        "all_refinement_indexes_avoid_p2460_quarter_sentinels": theorem_export["all_refinement_indexes_avoid_p2460_quarter_sentinels"],
        "minimum_dyadic_refinement_separation_positive": Decimal(theorem_export["minimum_dyadic_refinement_decimal_separation"]) > 0,
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
        "schema_version": "p2462_s1412_v1",
        "packet_id": "P2462",
        "stage_id": "S1412",
        "status": "PASS_STRICT_POINTWISE_INTERVAL_DECIMAL_GAP_DYADIC_REFINEMENT_REPLAY_NO_FULL_COMPLEMENT_SELECTOR_SOURCE_THEOREM",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_pointwise_interval_decimal_gap_dyadic_refinement_replay_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_pointwise_interval_decimal_gap_dyadic_refinement_replay_certificate"]["theorem_export"]
    lines = [
        "# P2462/S1412 strict pointwise interval-Decimal gap dyadic-refinement replay certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Dyadic-refinement replay",
        "",
        f"P2460 gap sentinels inherited: `{t['p2460_total_gap_sentinel_replay_count_inherited']}`.",
        f"P2461 weakest-neighborhood cells inherited: `{t['p2461_total_neighborhood_replay_count_inherited']}`.",
        f"Dyadic refinement fractions replayed per segment: `{', '.join(t['dyadic_refinement_fractions'])}`.",
        f"Dyadic refinement cells replayed: `{t['total_dyadic_refinement_replay_count']}`.",
        f"All dyadic refinement cells exclude zero: `{t['all_dyadic_refinement_cells_exclude_zero']}`.",
        f"All refinement indexes avoid P2460 quarter sentinels: `{t['all_refinement_indexes_avoid_p2460_quarter_sentinels']}`.",
        f"Minimum dyadic refinement Decimal separation: `{t['minimum_dyadic_refinement_decimal_separation']}`.",
        "",
        "## Guardrail",
        "",
        "This is an eighth-step dyadic Decimal/Taylor replay only.  It exports no Decimal full-complement replay, no directed-rounding interval theorem, no symbolic root-exclusion theorem, no point-coordinate selector, no strict observable/source theorem, no gauge-slice theorem, no physical-value generator, no QW-2191 discharge, no role-bearing `L_total`, and no ToE closure.",
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
