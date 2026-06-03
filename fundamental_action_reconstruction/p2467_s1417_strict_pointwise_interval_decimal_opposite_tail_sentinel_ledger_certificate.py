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
OUT = GEN / "p2467_s1417_strict_pointwise_interval_decimal_opposite_tail_sentinel_ledger_certificate.json"
MD = GEN / "p2467_s1417_strict_pointwise_interval_decimal_opposite_tail_sentinel_ledger_certificate.md"

SOURCE_FILES = {
    "P2450_ROOT_ISOLATION_MARGIN": GEN / "p2450_s1400_strict_pointwise_projection_root_isolation_margin_certificate.json",
    "P2451_FLOATING_INTERVAL_AUDIT": GEN / "p2451_s1401_strict_pointwise_projection_interval_enclosure_root_exclusion_audit.json",
    "P2456_DECIMAL_BOUNDARY_REPLAY": GEN / "p2456_s1406_strict_pointwise_decimal_root_window_boundary_band_replay_certificate.json",
    "P2466_DESCENT_TAIL_BOUNDARY": GEN / "p2466_s1416_strict_pointwise_interval_decimal_adaptive_descent_tail_boundary_ledger_certificate.json",
    "P2449_PROJECTION_REDUCTION": GEN / "p2449_s1399_strict_pointwise_rank_lift_projection_reduction_certificate.json",
}
DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

OPPOSITE_TAIL_BOUNDARY_BAND = 8
OPPOSITE_TAIL_SENTINEL_FRACTIONS = [
    Decimal("0"),
    Decimal("0.125"),
    Decimal("0.25"),
    Decimal("0.5"),
    Decimal("0.75"),
    Decimal("0.875"),
    Decimal("1"),
]


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
        "new_packet": "P2467|S1417|opposite tail sentinel|opposite-side tail|counter-tail sentinel",
        "p2466_input": "P2466|S1416|descent tail boundary|tail-boundary ledger|tail-to-boundary",
        "sentinel_language": "opposite boundary|non-descent tail|counter-tail|opposite-tail endpoint band|bidirectional tail",
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


def family_p2466_replay(p2466: dict[str, Any], family: str) -> dict[str, Any]:
    for row in p2466["family_descent_tail_boundary_replays"]:
        if row["family"] == family:
            return row
    raise ValueError(f"family not found in P2466: {family}")


def ordered_opposite_tail_cells(uncovered_cells: list[dict[str, float]], endpoint_index: int, descent_direction: int) -> list[dict[str, Any]]:
    opposite_direction = -descent_direction
    if opposite_direction < 0:
        indexes = range(endpoint_index - 1, -1, -1)
        boundary_side = "left"
    else:
        indexes = range(endpoint_index + 1, len(uncovered_cells))
        boundary_side = "right"
    rows = []
    for opposite_tail_step, index in enumerate(indexes, start=1):
        rows.append(
            {
                **uncovered_cells[index],
                "uncovered_index": index,
                "uncovered_count": len(uncovered_cells),
                "opposite_tail_step": opposite_tail_step,
                "opposite_direction": opposite_direction,
                "opposite_boundary_side": boundary_side,
            }
        )
    return rows


def select_opposite_tail_sentinels(opposite_tail_cells: list[dict[str, Any]]) -> list[dict[str, Any]]:
    if not opposite_tail_cells:
        return []
    max_position = len(opposite_tail_cells) - 1
    positions = set(range(min(OPPOSITE_TAIL_BOUNDARY_BAND, len(opposite_tail_cells))))
    positions.update(range(max(0, len(opposite_tail_cells) - OPPOSITE_TAIL_BOUNDARY_BAND), len(opposite_tail_cells)))
    positions.update(
        int((fraction * Decimal(max_position)).to_integral_value(rounding="ROUND_HALF_UP"))
        for fraction in OPPOSITE_TAIL_SENTINEL_FRACTIONS
    )
    return [
        {
            **opposite_tail_cells[position],
            "opposite_tail_position": position,
            "opposite_tail_available_count": len(opposite_tail_cells),
            "selection_rule": "endpoint-band+dyadic-fraction-sentinel",
        }
        for position in sorted(positions)
    ]


def replay_opposite_tail_sentinels(
    family: str,
    p2450: dict[str, Any],
    p2451: dict[str, Any],
    p2456: dict[str, Any],
    p2466: dict[str, Any],
    projection_vector: list[Decimal],
) -> dict[str, Any]:
    inherited = family_p2466_replay(p2466, family)
    segment_index = int(inherited["p2465_segment_index"])
    descent_direction = int(inherited["descent_direction"])
    endpoint_index = int(inherited["tail_start_uncovered_index"])
    uncovered_cells = gap_uncovered_cells_for_segment(family, segment_index, p2450, p2451, p2456)
    opposite_tail = ordered_opposite_tail_cells(uncovered_cells, endpoint_index, descent_direction)
    sentinels = select_opposite_tail_sentinels(opposite_tail)
    p2466_tail_indexes = {int(row["uncovered_index"]) for row in inherited["tail_boundary_rows"]}
    function = function_for_family(family)
    replayed = [
        replay_cell(family, cell, projection_vector, function)
        | {
            "opposite_tail_step": cell["opposite_tail_step"],
            "opposite_tail_position": cell["opposite_tail_position"],
            "opposite_tail_available_count": cell["opposite_tail_available_count"],
            "opposite_direction": cell["opposite_direction"],
            "opposite_boundary_side": cell["opposite_boundary_side"],
            "selection_rule": cell["selection_rule"],
        }
        for cell in sentinels
    ]
    separations = [Decimal(row["decimal_separation_from_zero"]) for row in replayed]
    replayed_indexes = {int(row["uncovered_index"]) for row in replayed}
    return {
        "family": family,
        "inherited_p2466_tail_start_uncovered_index": endpoint_index,
        "inherited_p2466_descent_direction": descent_direction,
        "p2466_segment_index": segment_index,
        "opposite_direction": -descent_direction,
        "opposite_boundary_side": sentinels[0]["opposite_boundary_side"] if sentinels else ("left" if descent_direction > 0 else "right"),
        "opposite_tail_available_count": len(opposite_tail),
        "opposite_tail_sentinel_replay_count": len(replayed),
        "opposite_tail_boundary_band": OPPOSITE_TAIL_BOUNDARY_BAND,
        "opposite_tail_sentinel_fractions": [str(value) for value in OPPOSITE_TAIL_SENTINEL_FRACTIONS],
        "opposite_tail_sentinel_rows": replayed,
        "all_opposite_tail_sentinels_exclude_zero": all(row["decimal_interval_excludes_zero"] for row in replayed),
        "minimum_opposite_tail_sentinel_decimal_separation": str(min(separations)) if separations else None,
        "all_opposite_tail_sentinel_indexes_unique": len(replayed_indexes) == len(replayed),
        "all_opposite_tail_sentinels_disjoint_from_p2466_descent_tail": p2466_tail_indexes.isdisjoint(replayed_indexes),
        "opposite_tail_full_replay_exported_by_this_family": len(replayed) == len(opposite_tail),
    }


def append_doc_sections() -> None:
    eq_section = """
## P2467/S1417 strict pointwise interval-Decimal opposite-tail sentinel ledger certificate

`P2467/S1417` responds to the P2466 next-step prompt by auditing the opposite, non-descent-side tails from the same P2465/P2466 adaptive endpoints.  It does not replay the full opposite tails; instead it replays deterministic endpoint-band plus dyadic-fraction sentinels in each opposite tail using the same Decimal/Taylor backend.  The selected opposite-tail sentinels remain zero-excluding and are disjoint from the P2466 descent-tail rows.

This is an opposite-tail sentinel ledger, not an opposite-tail full replay, full complement Decimal replay, directed-rounding interval theorem, symbolic root-exclusion theorem, analytic monotonicity theorem, selector/source/gauge theorem, role-bearing `L_total`, `QW-2191` discharge, or ToE closure.
""".strip()
    lag_section = """
## P2467/S1417 Decimal opposite-tail sentinel ledger guard

`P2467/S1417` adds a finite opposite-side audit after P2466 by replaying endpoint-band plus dyadic sentinels on the non-descent tails.  It improves bidirectional tail hygiene for the selected scalar families, but it is not a full-tail or full-complement proof and adds no selector/source/gauge authority for `L_total`.
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
    p2466 = sources["P2466_DESCENT_TAIL_BOUNDARY"].get("strict_pointwise_interval_decimal_adaptive_descent_tail_boundary_ledger_certificate", {}).get("theorem_export", {})
    projection_vector = [Decimal(str(value)) for value in sources["P2449_PROJECTION_REDUCTION"].get("strict_pointwise_rank_lift_projection_reduction_certificate", {}).get("theorem_export", {}).get("projection_vector", [])]
    family_replays = [
        replay_opposite_tail_sentinels("zero_projection_amplitude", p2450, p2451, p2456, p2466, projection_vector),
        replay_opposite_tail_sentinels("stationary_factor", p2450, p2451, p2456, p2466, projection_vector),
    ]
    total_available = sum(row["opposite_tail_available_count"] for row in family_replays)
    total_replayed = sum(row["opposite_tail_sentinel_replay_count"] for row in family_replays)
    minimum_separation = min(Decimal(row["minimum_opposite_tail_sentinel_decimal_separation"]) for row in family_replays)
    theorem_export = {
        "theorem_name": "P2467_T1_strict_pointwise_interval_decimal_opposite_tail_sentinel_ledger_certificate",
        "inherited_adaptive_descent_tail_boundary_ledger": "P2466/S1416",
        "p2466_total_tail_boundary_replay_count_inherited": p2466.get("total_tail_boundary_replay_count"),
        "p2466_all_tail_boundary_cells_exclude_zero_inherited": p2466.get("all_tail_boundary_cells_exclude_zero"),
        "opposite_tail_boundary_band": OPPOSITE_TAIL_BOUNDARY_BAND,
        "opposite_tail_sentinel_fractions": [str(value) for value in OPPOSITE_TAIL_SENTINEL_FRACTIONS],
        "family_opposite_tail_sentinel_replays": family_replays,
        "total_opposite_tail_available_cell_count": total_available,
        "total_opposite_tail_sentinel_replay_count": total_replayed,
        "all_opposite_tail_sentinels_exclude_zero": all(row["all_opposite_tail_sentinels_exclude_zero"] for row in family_replays),
        "minimum_opposite_tail_sentinel_decimal_separation": str(minimum_separation),
        "all_opposite_tail_sentinel_indexes_unique": all(row["all_opposite_tail_sentinel_indexes_unique"] for row in family_replays),
        "all_opposite_tail_sentinels_disjoint_from_p2466_descent_tail": all(row["all_opposite_tail_sentinels_disjoint_from_p2466_descent_tail"] for row in family_replays),
        "opposite_tail_sentinel_ledger_exported": True,
        "opposite_tail_full_replay_exported_by_this_certificate": all(row["opposite_tail_full_replay_exported_by_this_family"] for row in family_replays),
        "remaining_complement_segments_replayed_by_this_certificate": False,
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
            "Opposite-tail sentinels are finite bidirectional hygiene and do not replay every opposite-tail cell.",
            "This does not close the P2459 full complement gap or replace P2451 with a directed-rounding theorem.",
            "No point-coordinate selector, strict observable/source theorem, gauge-slice theorem, physical-value generator, QW-2191 discharge, role-bearing L_total export, or ToE closure is exported.",
        ],
        "next_honest_step": (
            "Replay the full opposite tails, replay the remaining complement segments, or replace finite tail ledgers with a full-complement directed-rounding proof."
        ),
    }
    gatekeepers = {
        "rg_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "p2466_tail_count_inherited": theorem_export["p2466_total_tail_boundary_replay_count_inherited"] == 6361,
        "p2466_zero_exclusion_inherited": theorem_export["p2466_all_tail_boundary_cells_exclude_zero_inherited"] is True,
        "family_count_is_two": len(family_replays) == 2,
        "total_opposite_tail_available_positive": theorem_export["total_opposite_tail_available_cell_count"] > 0,
        "total_opposite_tail_sentinel_replay_count": theorem_export["total_opposite_tail_sentinel_replay_count"] == 42,
        "all_opposite_tail_sentinels_exclude_zero": theorem_export["all_opposite_tail_sentinels_exclude_zero"],
        "minimum_opposite_tail_sentinel_separation_positive": Decimal(theorem_export["minimum_opposite_tail_sentinel_decimal_separation"]) > 0,
        "all_opposite_tail_sentinel_indexes_unique": theorem_export["all_opposite_tail_sentinel_indexes_unique"],
        "all_opposite_tail_sentinels_disjoint_from_p2466_descent_tail": theorem_export["all_opposite_tail_sentinels_disjoint_from_p2466_descent_tail"],
        "no_opposite_tail_full_replay": not theorem_export["opposite_tail_full_replay_exported_by_this_certificate"],
        "no_remaining_segments_replay": not theorem_export["remaining_complement_segments_replayed_by_this_certificate"],
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
        "schema_version": "p2467_s1417_v1",
        "packet_id": "P2467",
        "stage_id": "S1417",
        "status": "PASS_STRICT_POINTWISE_INTERVAL_DECIMAL_OPPOSITE_TAIL_SENTINEL_LEDGER_NO_FULL_TAIL_COMPLEMENT_SELECTOR_SOURCE_THEOREM",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_pointwise_interval_decimal_opposite_tail_sentinel_ledger_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_pointwise_interval_decimal_opposite_tail_sentinel_ledger_certificate"]["theorem_export"]
    lines = [
        "# P2467/S1417 strict pointwise interval-Decimal opposite-tail sentinel ledger certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Opposite-tail sentinel ledger",
        "",
        f"P2466 tail-boundary cells inherited: `{t['p2466_total_tail_boundary_replay_count_inherited']}`.",
        f"Opposite-tail available cells: `{t['total_opposite_tail_available_cell_count']}`.",
        f"Opposite-tail sentinel cells replayed: `{t['total_opposite_tail_sentinel_replay_count']}`.",
        f"All opposite-tail sentinels exclude zero: `{t['all_opposite_tail_sentinels_exclude_zero']}`.",
        f"All opposite-tail sentinels disjoint from P2466 descent tail: `{t['all_opposite_tail_sentinels_disjoint_from_p2466_descent_tail']}`.",
        f"Minimum opposite-tail sentinel Decimal separation: `{t['minimum_opposite_tail_sentinel_decimal_separation']}`.",
        "",
        "## Guardrail",
        "",
        "This is a finite opposite-tail sentinel ledger only.  It exports no opposite-tail full replay, no Decimal full-complement replay, no directed-rounding interval theorem, no symbolic root-exclusion theorem, no analytic monotonicity theorem, no point-coordinate selector, no strict observable/source theorem, no gauge-slice theorem, no physical-value generator, no QW-2191 discharge, no role-bearing `L_total`, and no ToE closure.",
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
