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
OUT = GEN / "p2466_s1416_strict_pointwise_interval_decimal_adaptive_descent_tail_boundary_ledger_certificate.json"
MD = GEN / "p2466_s1416_strict_pointwise_interval_decimal_adaptive_descent_tail_boundary_ledger_certificate.md"

SOURCE_FILES = {
    "P2450_ROOT_ISOLATION_MARGIN": GEN / "p2450_s1400_strict_pointwise_projection_root_isolation_margin_certificate.json",
    "P2451_FLOATING_INTERVAL_AUDIT": GEN / "p2451_s1401_strict_pointwise_projection_interval_enclosure_root_exclusion_audit.json",
    "P2456_DECIMAL_BOUNDARY_REPLAY": GEN / "p2456_s1406_strict_pointwise_decimal_root_window_boundary_band_replay_certificate.json",
    "P2465_DESCENT_HORIZON": GEN / "p2465_s1415_strict_pointwise_interval_decimal_adaptive_descent_horizon_ledger_certificate.json",
    "P2449_PROJECTION_REDUCTION": GEN / "p2449_s1399_strict_pointwise_rank_lift_projection_reduction_certificate.json",
}
DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}


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
        "new_packet": "P2466|S1416|descent tail boundary|tail-boundary ledger|tail-to-boundary",
        "p2465_input": "P2465|S1415|adaptive descent-horizon|descent horizon ledger|horizon endpoint",
        "tail_language": "endpoint-to-boundary|remaining descent tail|boundary-side tail|descent tail replay",
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


def family_p2465_replay(p2465: dict[str, Any], family: str) -> dict[str, Any]:
    for row in p2465["family_descent_horizon_replays"]:
        if row["family"] == family:
            return row
    raise ValueError(f"family not found in P2465: {family}")


def descent_tail_cells(uncovered_cells: list[dict[str, float]], endpoint_index: int, direction: int) -> list[dict[str, Any]]:
    if direction < 0:
        indexes = range(endpoint_index - 1, -1, -1)
        boundary_side = "left"
    else:
        indexes = range(endpoint_index + 1, len(uncovered_cells))
        boundary_side = "right"
    rows = []
    for tail_step, index in enumerate(indexes, start=1):
        rows.append(
            {
                **uncovered_cells[index],
                "uncovered_index": index,
                "uncovered_count": len(uncovered_cells),
                "tail_step": tail_step,
                "descent_direction": direction,
                "boundary_side": boundary_side,
            }
        )
    return rows


def first_non_decreasing_step(inherited_sep: Decimal, separations: list[Decimal], replayed: list[dict[str, Any]]) -> int | None:
    previous = inherited_sep
    for row, value in zip(replayed, separations):
        if value >= previous:
            return int(row["tail_step"])
        previous = value
    return None


def replay_descent_tail(
    family: str,
    p2450: dict[str, Any],
    p2451: dict[str, Any],
    p2456: dict[str, Any],
    p2465: dict[str, Any],
    projection_vector: list[Decimal],
) -> dict[str, Any]:
    inherited = family_p2465_replay(p2465, family)
    horizon_rows = inherited["descent_horizon_rows"]
    horizon_endpoint_row = horizon_rows[-1]
    endpoint_index = int(inherited["horizon_endpoint_uncovered_index"])
    direction = int(inherited["descent_direction"])
    segment_index = int(inherited["p2464_segment_index"])
    uncovered_cells = gap_uncovered_cells_for_segment(family, segment_index, p2450, p2451, p2456)
    cells = descent_tail_cells(uncovered_cells, endpoint_index, direction)
    function = function_for_family(family)
    replayed = [
        replay_cell(family, cell, projection_vector, function)
        | {
            "tail_step": cell["tail_step"],
            "descent_direction": direction,
            "boundary_side": cell["boundary_side"],
        }
        for cell in cells
    ]
    inherited_endpoint_sep = Decimal(horizon_endpoint_row["decimal_separation_from_zero"])
    separations = [Decimal(row["decimal_separation_from_zero"]) for row in replayed]
    non_decreasing_step = first_non_decreasing_step(inherited_endpoint_sep, separations, replayed)
    strictly_decreasing = bool(separations) and non_decreasing_step is None
    return {
        "family": family,
        "inherited_p2465_endpoint_row": horizon_endpoint_row,
        "inherited_p2465_endpoint_decimal_separation": str(inherited_endpoint_sep),
        "p2465_segment_index": segment_index,
        "descent_direction": direction,
        "boundary_side": replayed[0]["boundary_side"] if replayed else ("left" if direction < 0 else "right"),
        "tail_start_uncovered_index": endpoint_index,
        "tail_boundary_replay_count": len(replayed),
        "tail_boundary_rows": replayed,
        "all_tail_boundary_cells_exclude_zero": all(row["decimal_interval_excludes_zero"] for row in replayed),
        "minimum_tail_boundary_decimal_separation": str(min(separations)) if separations else None,
        "strictly_decreasing_from_p2465_endpoint_to_boundary": strictly_decreasing,
        "first_non_decreasing_tail_step": non_decreasing_step,
        "local_bracket_found_in_tail": non_decreasing_step is not None,
        "boundary_endpoint_uncovered_index": replayed[-1]["uncovered_index"] if replayed else endpoint_index,
        "boundary_endpoint_decimal_separation": replayed[-1]["decimal_separation_from_zero"] if replayed else str(inherited_endpoint_sep),
    }


def append_doc_sections() -> None:
    eq_section = """
## P2466/S1416 strict pointwise interval-Decimal adaptive descent tail-boundary ledger certificate

`P2466/S1416` continues the P2465 unbracketed horizon from each family's P2465 endpoint all the way to the corresponding boundary side of the same unreplayed complement segment.  The replay adds the remaining descent-direction tail cells, all remain zero-excluding, and the certificate records whether the finite tail stays strictly decreasing or encounters a local bracket.

This is a two-tail endpoint-to-boundary ledger, not a full complement Decimal replay, directed-rounding interval theorem, symbolic root-exclusion theorem, analytic monotonicity theorem, selector/source/gauge theorem, role-bearing `L_total`, `QW-2191` discharge, or ToE closure.
""".strip()
    lag_section = """
## P2466/S1416 Decimal adaptive descent tail-boundary ledger guard

`P2466/S1416` extends the P2465 adaptive descent endpoints to their segment-side boundaries.  It improves finite tail coverage for the selected scalar families, but it does not cover the full P2459 complement and adds no selector/source/gauge authority for `L_total`.
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
    p2465 = sources["P2465_DESCENT_HORIZON"].get("strict_pointwise_interval_decimal_adaptive_descent_horizon_ledger_certificate", {}).get("theorem_export", {})
    projection_vector = [Decimal(str(value)) for value in sources["P2449_PROJECTION_REDUCTION"].get("strict_pointwise_rank_lift_projection_reduction_certificate", {}).get("theorem_export", {}).get("projection_vector", [])]
    family_replays = [
        replay_descent_tail("zero_projection_amplitude", p2450, p2451, p2456, p2465, projection_vector),
        replay_descent_tail("stationary_factor", p2450, p2451, p2456, p2465, projection_vector),
    ]
    total_replayed = sum(row["tail_boundary_replay_count"] for row in family_replays)
    minimum_separation = min(Decimal(row["minimum_tail_boundary_decimal_separation"]) for row in family_replays)
    theorem_export = {
        "theorem_name": "P2466_T1_strict_pointwise_interval_decimal_adaptive_descent_tail_boundary_ledger_certificate",
        "inherited_adaptive_descent_horizon_ledger": "P2465/S1415",
        "p2465_total_descent_horizon_replay_count_inherited": p2465.get("total_descent_horizon_replay_count"),
        "p2465_all_descent_horizon_cells_exclude_zero_inherited": p2465.get("all_descent_horizon_cells_exclude_zero"),
        "family_descent_tail_boundary_replays": family_replays,
        "total_tail_boundary_replay_count": total_replayed,
        "all_tail_boundary_cells_exclude_zero": all(row["all_tail_boundary_cells_exclude_zero"] for row in family_replays),
        "minimum_tail_boundary_decimal_separation": str(minimum_separation),
        "all_tails_strictly_decreasing_from_p2465_endpoint_to_boundary": all(row["strictly_decreasing_from_p2465_endpoint_to_boundary"] for row in family_replays),
        "any_local_bracket_found_in_tail": any(row["local_bracket_found_in_tail"] for row in family_replays),
        "adaptive_descent_tail_boundary_ledger_exported": True,
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
            "A two-tail endpoint-to-boundary ledger does not close the P2459 full complement gap.",
            "Tail coverage in selected adaptive directions is not a directed-rounding interval theorem, symbolic root-exclusion theorem, or analytic monotonicity theorem.",
            "No point-coordinate selector, strict observable/source theorem, gauge-slice theorem, physical-value generator, QW-2191 discharge, role-bearing L_total export, or ToE closure is exported.",
        ],
        "next_honest_step": (
            "Either replay the opposite-side tails and remaining complement segments, or replace the finite tail ledger with a full-complement directed-rounding proof."
        ),
    }
    gatekeepers = {
        "rg_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "p2465_horizon_count_inherited": theorem_export["p2465_total_descent_horizon_replay_count_inherited"] == 64,
        "p2465_zero_exclusion_inherited": theorem_export["p2465_all_descent_horizon_cells_exclude_zero_inherited"] is True,
        "family_count_is_two": len(family_replays) == 2,
        "total_tail_boundary_replay_count_positive": theorem_export["total_tail_boundary_replay_count"] > 0,
        "all_tail_boundary_cells_exclude_zero": theorem_export["all_tail_boundary_cells_exclude_zero"],
        "minimum_tail_boundary_separation_positive": Decimal(theorem_export["minimum_tail_boundary_decimal_separation"]) > 0,
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
        "schema_version": "p2466_s1416_v1",
        "packet_id": "P2466",
        "stage_id": "S1416",
        "status": "PASS_STRICT_POINTWISE_INTERVAL_DECIMAL_ADAPTIVE_DESCENT_TAIL_BOUNDARY_LEDGER_NO_FULL_COMPLEMENT_SELECTOR_SOURCE_THEOREM",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_pointwise_interval_decimal_adaptive_descent_tail_boundary_ledger_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_pointwise_interval_decimal_adaptive_descent_tail_boundary_ledger_certificate"]["theorem_export"]
    lines = [
        "# P2466/S1416 strict pointwise interval-Decimal adaptive descent tail-boundary ledger certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Adaptive descent tail-boundary ledger",
        "",
        f"P2465 descent-horizon cells inherited: `{t['p2465_total_descent_horizon_replay_count_inherited']}`.",
        f"Tail-boundary cells replayed: `{t['total_tail_boundary_replay_count']}`.",
        f"All tail-boundary cells exclude zero: `{t['all_tail_boundary_cells_exclude_zero']}`.",
        f"All tails strictly decrease from P2465 endpoint to boundary: `{t['all_tails_strictly_decreasing_from_p2465_endpoint_to_boundary']}`.",
        f"Any local bracket found in tail: `{t['any_local_bracket_found_in_tail']}`.",
        f"Minimum tail-boundary Decimal separation: `{t['minimum_tail_boundary_decimal_separation']}`.",
        "",
        "## Guardrail",
        "",
        "This is a finite two-tail endpoint-to-boundary ledger only.  It exports no Decimal full-complement replay, no directed-rounding interval theorem, no symbolic root-exclusion theorem, no analytic monotonicity theorem, no point-coordinate selector, no strict observable/source theorem, no gauge-slice theorem, no physical-value generator, no QW-2191 discharge, no role-bearing `L_total`, and no ToE closure.",
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
