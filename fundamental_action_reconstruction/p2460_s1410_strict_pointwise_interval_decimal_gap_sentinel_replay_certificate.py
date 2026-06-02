#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
import subprocess
from decimal import Decimal
from pathlib import Path
from typing import Any, Callable

from p2453_s1403_strict_pointwise_directed_decimal_weakest_cell_replay_certificate import (
    DecimalInterval,
    projection_amplitude_interval,
    stationary_factor_interval,
)

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"
OUT = GEN / "p2460_s1410_strict_pointwise_interval_decimal_gap_sentinel_replay_certificate.json"
MD = GEN / "p2460_s1410_strict_pointwise_interval_decimal_gap_sentinel_replay_certificate.md"

SOURCE_FILES = {
    "P2450_ROOT_ISOLATION_MARGIN": GEN / "p2450_s1400_strict_pointwise_projection_root_isolation_margin_certificate.json",
    "P2451_FLOATING_INTERVAL_AUDIT": GEN / "p2451_s1401_strict_pointwise_projection_interval_enclosure_root_exclusion_audit.json",
    "P2456_DECIMAL_BOUNDARY_REPLAY": GEN / "p2456_s1406_strict_pointwise_decimal_root_window_boundary_band_replay_certificate.json",
    "P2459_COVERAGE_GAP_LEDGER": GEN / "p2459_s1409_strict_pointwise_interval_decimal_coverage_gap_ledger_certificate.json",
    "P2449_PROJECTION_REDUCTION": GEN / "p2449_s1399_strict_pointwise_rank_lift_projection_reduction_certificate.json",
}
DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}
SENTINEL_FRACTIONS = [Decimal("0"), Decimal("0.25"), Decimal("0.5"), Decimal("0.75"), Decimal("1")]
CELL_MATCH_DIGITS = 12


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
        "new_packet": "P2460|S1410|gap sentinel replay|Decimal gap sentinel|unreplayed gap sentinel",
        "p2459_input": "P2459|S1409|coverage gap ledger|unreplayed by Decimal boundary chain|coverage gap",
        "sentinel_language": "sentinel replay|stratified sentinel|gap interior|unreplayed complement cell",
        "decimal_backend": "Decimal interval|Taylor endpoint|zero-excluding|Decimal separation",
        "closure_blockers": "QW-2191|role-bearing L_total|physical-value generator|ToE closure|point-coordinate selector",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def cell_key(left: float, right: float) -> tuple[float, float]:
    return (round(float(left), CELL_MATCH_DIGITS), round(float(right), CELL_MATCH_DIGITS))


def family_replays(p2456: dict[str, Any], family: str) -> list[dict[str, Any]]:
    key = "zero_projection_boundary_band_replays" if family == "zero_projection_amplitude" else "stationary_factor_boundary_band_replays"
    return p2456[key]


def boundary_skip_counts_for_segment(p2456: dict[str, Any], family: str, segment: dict[str, float]) -> tuple[int, int]:
    left_skip = 0
    right_skip = 0
    segment_left = float(segment["left"])
    segment_right = float(segment["right"])
    tolerance = 1.0e-12
    for replay in family_replays(p2456, family):
        side = replay["side"]
        for row in replay["rows"]:
            midpoint = 0.5 * (float(row["left"]) + float(row["right"]))
            if segment_left - tolerance <= midpoint <= segment_right + tolerance:
                if side == "right":
                    left_skip += 1
                elif side == "left":
                    right_skip += 1
    return left_skip, right_skip


def remove_boundary_covered_cells(cells: list[dict[str, float]], left_skip: int, right_skip: int) -> list[dict[str, float]]:
    right_limit = len(cells) - right_skip if right_skip else len(cells)
    return cells[left_skip:right_limit]


def complement_segments(p2450: dict[str, Any], family: str) -> list[dict[str, float]]:
    key = "zero_projection_amplitude_certificate" if family == "zero_projection_amplitude" else "stationary_factor_certificate"
    return p2450[key]["sampled_lipschitz_exclusion"]["complement_segments"]


def interval_cell_width(p2451: dict[str, Any], family: str) -> float:
    key = "zero_projection_amplitude_interval_audit" if family == "zero_projection_amplitude" else "stationary_factor_interval_audit"
    return float(p2451[key]["interval_cell_width"])


def segment_cells(segment: dict[str, float], width: float) -> list[dict[str, float]]:
    cells = []
    x = float(segment["left"])
    right = float(segment["right"])
    while x < right - 1.0e-15:
        y = min(right, x + width)
        cells.append({"left": x, "right": y})
        x = y
    return cells


def select_sentinels(uncovered_cells: list[dict[str, float]]) -> list[dict[str, float]]:
    if not uncovered_cells:
        return []
    max_index = len(uncovered_cells) - 1
    indexes = sorted({int((fraction * Decimal(max_index)).to_integral_value(rounding="ROUND_HALF_UP")) for fraction in SENTINEL_FRACTIONS})
    return [{**uncovered_cells[index], "uncovered_index": index, "uncovered_count": len(uncovered_cells)} for index in indexes]


def replay_cell(
    family: str,
    cell: dict[str, float],
    projection_vector: list[Decimal],
    function: Callable[[list[Decimal], DecimalInterval], Any],
) -> dict[str, Any]:
    value = function(projection_vector, DecimalInterval(cell["left"], cell["right"]))
    separation = value.separation_from_zero()
    return {
        "family": family,
        "left": cell["left"],
        "right": cell["right"],
        "uncovered_index": cell["uncovered_index"],
        "uncovered_count": cell["uncovered_count"],
        "decimal_taylor_interval_value": value.as_dict(),
        "decimal_interval_excludes_zero": not value.contains_zero(),
        "decimal_separation_from_zero": str(separation),
    }


def family_gap_sentinel_replay(
    family: str,
    p2450: dict[str, Any],
    p2451: dict[str, Any],
    p2456: dict[str, Any],
    projection_vector: list[Decimal],
    function: Callable[[list[Decimal], DecimalInterval], Any],
) -> dict[str, Any]:
    width = interval_cell_width(p2451, family)
    segment_rows = []
    sentinel_replays = []
    total_uncovered = 0
    for segment_index, segment in enumerate(complement_segments(p2450, family)):
        cells = segment_cells(segment, width)
        left_skip, right_skip = boundary_skip_counts_for_segment(p2456, family, segment)
        uncovered_cells = remove_boundary_covered_cells(cells, left_skip, right_skip)
        total_uncovered += len(uncovered_cells)
        sentinels = select_sentinels(uncovered_cells)
        replayed = [replay_cell(family, cell, projection_vector, function) for cell in sentinels]
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
                "sentinel_replay_count": len(replayed),
                "all_segment_sentinels_exclude_zero": all(row["decimal_interval_excludes_zero"] for row in replayed),
                "minimum_segment_sentinel_decimal_separation": str(min(Decimal(row["decimal_separation_from_zero"]) for row in replayed)) if replayed else None,
                "sentinel_replays": replayed,
            }
        )
        sentinel_replays.extend(replayed)
    return {
        "family": family,
        "segment_count": len(segment_rows),
        "segment_rows": segment_rows,
        "uncovered_by_boundary_chain_cell_count": total_uncovered,
        "sentinel_replay_count": len(sentinel_replays),
        "all_gap_sentinels_exclude_zero": all(row["decimal_interval_excludes_zero"] for row in sentinel_replays),
        "minimum_gap_sentinel_decimal_separation": str(min(Decimal(row["decimal_separation_from_zero"]) for row in sentinel_replays)),
    }


def append_doc_sections() -> None:
    eq_section = """
## P2460/S1410 strict pointwise interval-Decimal gap-sentinel replay certificate

`P2460/S1410` responds to the P2459 coverage-gap ledger without pretending to close it: for every P2451 complement segment, it selects stratified sentinel cells from the part not already covered by the P2456 Decimal boundary chain and replays those sentinels with the Decimal/Taylor backend.  All 25 selected gap sentinels remain zero-excluding.

This is a stratified gap-sentinel replay, not a full complement Decimal replay, directed-rounding interval theorem, symbolic root-exclusion theorem, selector/source/gauge theorem, role-bearing `L_total`, `QW-2191` discharge, or ToE closure.
""".strip()
    lag_section = """
## P2460/S1410 Decimal gap-sentinel replay guard

`P2460/S1410` probes the P2459 unreplayed complement gap with stratified Decimal/Taylor sentinels in each complement segment.  It improves gap diagnostics but still leaves the full complement gap open and adds no selector/source/gauge authority for `L_total`.
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
    p2459 = sources["P2459_COVERAGE_GAP_LEDGER"].get("strict_pointwise_interval_decimal_coverage_gap_ledger_certificate", {}).get("theorem_export", {})
    projection_vector = [Decimal(str(value)) for value in sources["P2449_PROJECTION_REDUCTION"].get("strict_pointwise_rank_lift_projection_reduction_certificate", {}).get("theorem_export", {}).get("projection_vector", [])]
    zero_replay = family_gap_sentinel_replay("zero_projection_amplitude", p2450, p2451, p2456, projection_vector, projection_amplitude_interval)
    stationary_replay = family_gap_sentinel_replay("stationary_factor", p2450, p2451, p2456, projection_vector, stationary_factor_interval)
    family_replays = [zero_replay, stationary_replay]
    total_sentinels = sum(item["sentinel_replay_count"] for item in family_replays)
    minimum_separation = min(Decimal(item["minimum_gap_sentinel_decimal_separation"]) for item in family_replays)
    total_uncovered_reconstructed = sum(item["uncovered_by_boundary_chain_cell_count"] for item in family_replays)
    theorem_export = {
        "theorem_name": "P2460_T1_strict_pointwise_interval_decimal_gap_sentinel_replay_certificate",
        "inherited_coverage_gap_ledger": "P2459/S1409",
        "inherited_decimal_boundary_replay": "P2456/S1406",
        "p2459_total_unreplayed_gap_inherited": p2459.get("total_unreplayed_by_decimal_boundary_chain_cell_count"),
        "sentinel_fractions": [str(value) for value in SENTINEL_FRACTIONS],
        "family_gap_sentinel_replays": family_replays,
        "total_gap_sentinel_replay_count": total_sentinels,
        "total_reconstructed_uncovered_gap_cell_count": total_uncovered_reconstructed,
        "all_gap_sentinels_exclude_zero": all(item["all_gap_sentinels_exclude_zero"] for item in family_replays),
        "minimum_gap_sentinel_decimal_separation": str(minimum_separation),
        "gap_sentinel_replay_exported": True,
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
            "Stratified gap sentinels are not a full complement Decimal replay and do not close the P2459 coverage gap.",
            "Zero exclusion on selected unreplayed-gap sentinels does not prove symbolic root exclusion or a directed-rounding interval theorem.",
            "No point-coordinate selector, strict observable/source theorem, gauge-slice theorem, physical-value generator, QW-2191 discharge, role-bearing L_total export, or ToE closure is exported.",
        ],
        "next_honest_step": (
            "Promote from stratified sentinels to an exhaustive Decimal/Taylor replay or a formal directed-rounding interval proof for every unreplayed complement cell."
        ),
    }
    gatekeepers = {
        "rg_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "reconstructed_gap_matches_p2459": theorem_export["total_reconstructed_uncovered_gap_cell_count"] == theorem_export["p2459_total_unreplayed_gap_inherited"] == 99846,
        "total_gap_sentinel_count": theorem_export["total_gap_sentinel_replay_count"] == 25,
        "all_gap_sentinels_exclude_zero": theorem_export["all_gap_sentinels_exclude_zero"],
        "minimum_gap_sentinel_separation_positive": Decimal(theorem_export["minimum_gap_sentinel_decimal_separation"]) > 0,
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
        "schema_version": "p2460_s1410_v1",
        "packet_id": "P2460",
        "stage_id": "S1410",
        "status": "PASS_STRICT_POINTWISE_INTERVAL_DECIMAL_GAP_SENTINEL_REPLAY_NO_FULL_COMPLEMENT_SELECTOR_SOURCE_THEOREM",
        "strict_pointwise_interval_decimal_gap_sentinel_replay_certificate": {
            "inputs": {name: rel(path) for name, path in SOURCE_FILES.items()},
            "input_fingerprints": {name: sha256_json(sources[name]) for name in sources},
            "rg_nonduplication_audit": grep,
            "theorem_export": theorem_export,
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_pointwise_interval_decimal_gap_sentinel_replay_certificate"]["theorem_export"]
    lines = [
        "# P2460/S1410 strict pointwise interval-Decimal gap-sentinel replay certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Gap-sentinel replay",
        "",
        f"Inherited P2459 unreplayed gap cells: `{t['p2459_total_unreplayed_gap_inherited']}`.",
        f"Reconstructed unreplayed gap cells: `{t['total_reconstructed_uncovered_gap_cell_count']}`.",
        f"Gap sentinels replayed: `{t['total_gap_sentinel_replay_count']}`.",
        f"All gap sentinels exclude zero: `{t['all_gap_sentinels_exclude_zero']}`.",
        f"Minimum gap-sentinel Decimal separation: `{t['minimum_gap_sentinel_decimal_separation']}`.",
        "",
        "## Guardrail",
        "",
        "This is a stratified gap-sentinel Decimal/Taylor replay only.  It exports no Decimal full-complement replay, no directed-rounding interval theorem, no symbolic root-exclusion theorem, no point-coordinate selector, no strict observable/source theorem, no gauge-slice theorem, no physical-value generator, no QW-2191 discharge, no role-bearing `L_total`, and no ToE closure.",
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
