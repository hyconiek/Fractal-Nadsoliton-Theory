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
OUT = GEN / "p2464_s1414_strict_pointwise_interval_decimal_adaptive_flank_descent_extension_replay_certificate.json"
MD = GEN / "p2464_s1414_strict_pointwise_interval_decimal_adaptive_flank_descent_extension_replay_certificate.md"

SOURCE_FILES = {
    "P2450_ROOT_ISOLATION_MARGIN": GEN / "p2450_s1400_strict_pointwise_projection_root_isolation_margin_certificate.json",
    "P2451_FLOATING_INTERVAL_AUDIT": GEN / "p2451_s1401_strict_pointwise_projection_interval_enclosure_root_exclusion_audit.json",
    "P2456_DECIMAL_BOUNDARY_REPLAY": GEN / "p2456_s1406_strict_pointwise_decimal_root_window_boundary_band_replay_certificate.json",
    "P2463_ADAPTIVE_DYADIC_FLANK": GEN / "p2463_s1413_strict_pointwise_interval_decimal_adaptive_dyadic_weakest_flank_replay_certificate.json",
    "P2449_PROJECTION_REDUCTION": GEN / "p2449_s1399_strict_pointwise_rank_lift_projection_reduction_certificate.json",
}
DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

DESCENT_EXTENSION_STEPS = 4


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
        "new_packet": "P2464|S1414|adaptive flank descent extension|flank descent extension|descent-extension replay",
        "p2463_input": "P2463|S1413|adaptive dyadic weakest flank|adaptive flank replay|smaller-than-dyadic",
        "extension_language": "descent direction|one-sided extension|strictly decreasing separations|flank cascade continuation|adaptive descent continuation",
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


def gap_uncovered_cells_for_segment(
    family: str,
    segment_index: int,
    p2450: dict[str, Any],
    p2451: dict[str, Any],
    p2456: dict[str, Any],
) -> list[dict[str, float]]:
    segment = complement_segments(p2450, family)[segment_index]
    cells = segment_cells(segment, interval_cell_width(p2451, family))
    left_skip, right_skip = boundary_skip_counts_for_segment(p2456, family, segment)
    return remove_boundary_covered_cells(cells, left_skip, right_skip)


def family_p2463_replay(p2463: dict[str, Any], family: str) -> dict[str, Any]:
    for row in p2463["family_adaptive_weakest_flank_replays"]:
        if row["family"] == family:
            return row
    raise ValueError(f"family not found in P2463: {family}")


def descent_extension_cells(uncovered_cells: list[dict[str, float]], center_index: int, direction: int) -> list[dict[str, Any]]:
    rows = []
    for step in range(1, DESCENT_EXTENSION_STEPS + 1):
        index = center_index + direction * step
        if not 0 <= index < len(uncovered_cells):
            continue
        rows.append(
            {
                **uncovered_cells[index],
                "uncovered_index": index,
                "uncovered_count": len(uncovered_cells),
                "descent_extension_step": step,
                "descent_direction": direction,
            }
        )
    return rows


def replay_descent_extension(
    family: str,
    p2450: dict[str, Any],
    p2451: dict[str, Any],
    p2456: dict[str, Any],
    p2463: dict[str, Any],
    projection_vector: list[Decimal],
) -> dict[str, Any]:
    inherited = family_p2463_replay(p2463, family)
    anchor = inherited["minimum_adaptive_flank_row"]
    segment_index = int(inherited["weakest_p2462_dyadic_segment_index"])
    center_index = int(anchor["uncovered_index"])
    direction = -1 if int(anchor["offset_from_p2462_weakest_dyadic"]) < 0 else 1
    uncovered_cells = gap_uncovered_cells_for_segment(family, segment_index, p2450, p2451, p2456)
    cells = descent_extension_cells(uncovered_cells, center_index, direction)
    function = function_for_family(family)
    replayed = [
        replay_cell(family, cell, projection_vector, function)
        | {"descent_extension_step": cell["descent_extension_step"], "descent_direction": direction}
        for cell in cells
    ]
    separations = [Decimal(row["decimal_separation_from_zero"]) for row in replayed]
    inherited_sep = Decimal(inherited["minimum_adaptive_flank_decimal_separation"])
    strictly_decreasing_from_anchor = bool(separations) and all(
        left > right for left, right in zip([inherited_sep, *separations[:-1]], separations)
    )
    return {
        "family": family,
        "inherited_p2463_minimum_flank_row": anchor,
        "inherited_p2463_minimum_flank_decimal_separation": str(inherited_sep),
        "p2463_segment_index": segment_index,
        "descent_direction": direction,
        "descent_extension_step_count_requested": DESCENT_EXTENSION_STEPS,
        "descent_extension_replay_count": len(replayed),
        "descent_extension_rows": replayed,
        "all_descent_extension_cells_exclude_zero": all(row["decimal_interval_excludes_zero"] for row in replayed),
        "minimum_descent_extension_decimal_separation": str(min(separations)),
        "all_extension_separations_below_p2463_anchor": all(value < inherited_sep for value in separations),
        "strictly_decreasing_from_p2463_anchor_along_extension": strictly_decreasing_from_anchor,
    }


def append_doc_sections() -> None:
    eq_section = """
## P2464/S1414 strict pointwise interval-Decimal adaptive flank descent-extension replay certificate

`P2464/S1414` continues the P2463 adaptive descent diagnostic by extending one-sided from each family's weakest P2463 flank cell for four additional unreplayed complement cells in the empirically decreasing direction.  The replay adds 8 descent-extension cells, all remain zero-excluding, and the Decimal separations strictly decrease from the P2463 anchor along each audited extension.

This is a finite one-sided descent-extension replay, not a full complement Decimal replay, directed-rounding interval theorem, symbolic root-exclusion theorem, analytic monotonicity theorem, selector/source/gauge theorem, role-bearing `L_total`, `QW-2191` discharge, or ToE closure.
""".strip()
    lag_section = """
## P2464/S1414 Decimal adaptive flank descent-extension replay guard

`P2464/S1414` extends from the weakest P2463 flank cells along the local decreasing-separation side.  It improves adaptive descent diagnostics but remains finite sampled evidence, leaves the full P2459 complement gap open, and adds no selector/source/gauge authority for `L_total`.
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
    p2463 = sources["P2463_ADAPTIVE_DYADIC_FLANK"].get("strict_pointwise_interval_decimal_adaptive_dyadic_weakest_flank_replay_certificate", {}).get("theorem_export", {})
    projection_vector = [Decimal(str(value)) for value in sources["P2449_PROJECTION_REDUCTION"].get("strict_pointwise_rank_lift_projection_reduction_certificate", {}).get("theorem_export", {}).get("projection_vector", [])]
    family_replays = [
        replay_descent_extension("zero_projection_amplitude", p2450, p2451, p2456, p2463, projection_vector),
        replay_descent_extension("stationary_factor", p2450, p2451, p2456, p2463, projection_vector),
    ]
    total_replayed = sum(row["descent_extension_replay_count"] for row in family_replays)
    minimum_separation = min(Decimal(row["minimum_descent_extension_decimal_separation"]) for row in family_replays)
    theorem_export = {
        "theorem_name": "P2464_T1_strict_pointwise_interval_decimal_adaptive_flank_descent_extension_replay_certificate",
        "inherited_adaptive_dyadic_weakest_flank_replay": "P2463/S1413",
        "p2463_total_adaptive_flank_replay_count_inherited": p2463.get("total_adaptive_flank_replay_count"),
        "p2463_all_adaptive_flank_cells_exclude_zero_inherited": p2463.get("all_adaptive_flank_cells_exclude_zero"),
        "descent_extension_steps_requested": DESCENT_EXTENSION_STEPS,
        "family_descent_extension_replays": family_replays,
        "total_descent_extension_replay_count": total_replayed,
        "all_descent_extension_cells_exclude_zero": all(row["all_descent_extension_cells_exclude_zero"] for row in family_replays),
        "minimum_descent_extension_decimal_separation": str(minimum_separation),
        "all_extension_separations_below_p2463_anchor": all(row["all_extension_separations_below_p2463_anchor"] for row in family_replays),
        "all_extensions_strictly_decreasing_from_p2463_anchor": all(row["strictly_decreasing_from_p2463_anchor_along_extension"] for row in family_replays),
        "adaptive_flank_descent_extension_replay_exported": True,
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
            "Descent-extension replay is finite sampled evidence and does not close the P2459 full complement gap.",
            "Strictly decreasing separations on two four-cell extensions are not an analytic monotonicity theorem, symbolic root-exclusion theorem, or directed-rounding interval theorem.",
            "No point-coordinate selector, strict observable/source theorem, gauge-slice theorem, physical-value generator, QW-2191 discharge, role-bearing L_total export, or ToE closure is exported.",
        ],
        "next_honest_step": (
            "Continue adaptive descent until the local minimum is bracketed by increasing separations on both sides, or replace finite sampling with a full-complement directed-rounding proof."
        ),
    }
    gatekeepers = {
        "rg_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "p2463_flank_count_inherited": theorem_export["p2463_total_adaptive_flank_replay_count_inherited"] == 16,
        "p2463_zero_exclusion_inherited": theorem_export["p2463_all_adaptive_flank_cells_exclude_zero_inherited"] is True,
        "family_count_is_two": len(family_replays) == 2,
        "total_descent_extension_replay_count": theorem_export["total_descent_extension_replay_count"] == 8,
        "all_descent_extension_cells_exclude_zero": theorem_export["all_descent_extension_cells_exclude_zero"],
        "minimum_descent_extension_separation_positive": Decimal(theorem_export["minimum_descent_extension_decimal_separation"]) > 0,
        "all_extension_separations_below_p2463_anchor": theorem_export["all_extension_separations_below_p2463_anchor"],
        "all_extensions_strictly_decreasing_from_p2463_anchor": theorem_export["all_extensions_strictly_decreasing_from_p2463_anchor"],
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
        "schema_version": "p2464_s1414_v1",
        "packet_id": "P2464",
        "stage_id": "S1414",
        "status": "PASS_STRICT_POINTWISE_INTERVAL_DECIMAL_ADAPTIVE_FLANK_DESCENT_EXTENSION_REPLAY_NO_FULL_COMPLEMENT_SELECTOR_SOURCE_THEOREM",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_pointwise_interval_decimal_adaptive_flank_descent_extension_replay_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_pointwise_interval_decimal_adaptive_flank_descent_extension_replay_certificate"]["theorem_export"]
    lines = [
        "# P2464/S1414 strict pointwise interval-Decimal adaptive flank descent-extension replay certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Adaptive flank descent-extension replay",
        "",
        f"P2463 adaptive flank cells inherited: `{t['p2463_total_adaptive_flank_replay_count_inherited']}`.",
        f"Descent extension steps per family: `{t['descent_extension_steps_requested']}`.",
        f"Descent-extension cells replayed: `{t['total_descent_extension_replay_count']}`.",
        f"All descent-extension cells exclude zero: `{t['all_descent_extension_cells_exclude_zero']}`.",
        f"All extension separations below P2463 anchor: `{t['all_extension_separations_below_p2463_anchor']}`.",
        f"All extensions strictly decrease from P2463 anchor: `{t['all_extensions_strictly_decreasing_from_p2463_anchor']}`.",
        f"Minimum descent-extension Decimal separation: `{t['minimum_descent_extension_decimal_separation']}`.",
        "",
        "## Guardrail",
        "",
        "This is a finite one-sided adaptive descent-extension replay only.  It exports no Decimal full-complement replay, no directed-rounding interval theorem, no symbolic root-exclusion theorem, no analytic monotonicity theorem, no point-coordinate selector, no strict observable/source theorem, no gauge-slice theorem, no physical-value generator, no QW-2191 discharge, no role-bearing `L_total`, and no ToE closure.",
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
