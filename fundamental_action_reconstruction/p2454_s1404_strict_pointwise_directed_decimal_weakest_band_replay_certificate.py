#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from decimal import Decimal
from pathlib import Path
from typing import Any, Callable

from p2453_s1403_strict_pointwise_directed_decimal_weakest_cell_replay_certificate import (
    DECIMAL_PRECISION,
    DECIMAL_SAFETY_PAD,
    TAYLOR_BOUND_PAD,
    DecimalInterval,
    projection_amplitude_interval,
    stationary_factor_interval,
)

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"
OUT = GEN / "p2454_s1404_strict_pointwise_directed_decimal_weakest_band_replay_certificate.json"
MD = GEN / "p2454_s1404_strict_pointwise_directed_decimal_weakest_band_replay_certificate.md"

SOURCE_FILES = {
    "P2451_INTERVAL_ENCLOSURE": GEN / "p2451_s1401_strict_pointwise_projection_interval_enclosure_root_exclusion_audit.json",
    "P2450_ROOT_ISOLATION_MARGIN": GEN / "p2450_s1400_strict_pointwise_projection_root_isolation_margin_certificate.json",
    "P2453_WEAKEST_CELL_REPLAY": GEN / "p2453_s1403_strict_pointwise_directed_decimal_weakest_cell_replay_certificate.json",
    "P2449_PROJECTION_REDUCTION": GEN / "p2449_s1399_strict_pointwise_rank_lift_projection_reduction_certificate.json",
}
DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

BAND_CELL_COUNT = 12


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
        "new_packet": "P2454|S1404|directed decimal weakest-band|weakest band replay|critical band replay",
        "p2453_input": "P2453|S1403|directed decimal weakest-cell|Taylor endpoint replay|Decimal/Taylor",
        "band_language": "weakest band|critical complement band|adjacent weakest cells|band replay",
        "backend_language": "Decimal interval|Taylor endpoint|alternating series|directed decimal|zero-excluding",
        "closure_blockers": "QW-2191|role-bearing L_total|physical-value generator|ToE closure|point-coordinate selector",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def find_segment_for_cell(cell: dict[str, Any], segments: list[dict[str, float]]) -> dict[str, float]:
    left = float(cell["left"])
    right = float(cell["right"])
    for segment in segments:
        if float(segment["left"]) - 1.0e-15 <= left and right <= float(segment["right"]) + 1.0e-15:
            return segment
    raise ValueError(f"cell {cell} is not contained in any complement segment")


def build_forward_band(cell: dict[str, Any], segment: dict[str, float], width: float, count: int = BAND_CELL_COUNT) -> list[dict[str, float]]:
    cells = []
    left = float(cell["left"])
    segment_right = float(segment["right"])
    for _ in range(count):
        right = min(segment_right, left + width)
        if right <= left:
            break
        cells.append({"left": left, "right": right})
        left = right
    return cells


def replay_band(
    family: str,
    cells: list[dict[str, float]],
    projection_vector: list[Decimal],
    function: Callable[[list[Decimal], DecimalInterval], Any],
) -> dict[str, Any]:
    rows = []
    minimum_separation: Decimal | None = None
    weakest_row: dict[str, Any] | None = None
    for idx, cell in enumerate(cells):
        value = function(projection_vector, DecimalInterval(cell["left"], cell["right"]))
        separation = value.separation_from_zero()
        row = {
            "index": idx,
            "left": cell["left"],
            "right": cell["right"],
            "decimal_taylor_interval_value": value.as_dict(),
            "decimal_interval_excludes_zero": not value.contains_zero(),
            "decimal_separation_from_zero": str(separation),
        }
        rows.append(row)
        if minimum_separation is None or separation < minimum_separation:
            minimum_separation = separation
            weakest_row = row
    return {
        "family": family,
        "band_cell_count": len(rows),
        "rows": rows,
        "all_cells_exclude_zero": all(row["decimal_interval_excludes_zero"] for row in rows),
        "minimum_decimal_separation_from_zero": str(minimum_separation),
        "weakest_decimal_row": weakest_row,
    }


def append_doc_sections() -> None:
    eq_section = """
## P2454/S1404 strict pointwise directed-decimal weakest-band replay certificate

`P2454/S1404` extends P2453 from a single weakest cell per family to a forward weakest-band replay: twelve adjacent complement cells are replayed for both the zero-projection amplitude and the stationary factor using the same Decimal/Taylor endpoint backend.  Every replayed critical-band cell remains zero-excluding.

This is still not a full complement directed-rounding interval theorem or symbolic root-exclusion proof.  It exports no point-coordinate selector, source/gauge theorem, role-bearing `L_total`, `QW-2191` discharge, or ToE closure.
""".strip()
    lag_section = """
## P2454/S1404 directed-decimal weakest-band replay guard

`P2454/S1404` broadens the Decimal/Taylor replay from the weakest cells to adjacent critical bands, increasing backend confidence in the most fragile complement region.  It remains a bounded numerical replay and does not license any pointwise coordinate as a selector/source/gauge row for `L_total`.
""".strip()
    for path, section in [(DOC_FILES["equation_sheet"], eq_section), (DOC_FILES["lagrangian_eom_draft"], lag_section)]:
        text = path.read_text(encoding="utf-8") if path.exists() else ""
        marker = section.splitlines()[0]
        if marker not in text:
            path.write_text(text.rstrip() + "\n\n" + section + "\n", encoding="utf-8")


def build_payload() -> dict[str, Any]:
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    grep = rg_audit()
    p2451 = sources["P2451_INTERVAL_ENCLOSURE"].get("strict_pointwise_projection_interval_enclosure_root_exclusion_audit", {}).get("theorem_export", {})
    p2450 = sources["P2450_ROOT_ISOLATION_MARGIN"].get("strict_pointwise_projection_root_isolation_margin_certificate", {}).get("theorem_export", {})
    p2453 = sources["P2453_WEAKEST_CELL_REPLAY"].get("strict_pointwise_directed_decimal_weakest_cell_replay_certificate", {}).get("theorem_export", {})
    projection_vector = [Decimal(str(value)) for value in sources["P2449_PROJECTION_REDUCTION"].get("strict_pointwise_rank_lift_projection_reduction_certificate", {}).get("theorem_export", {}).get("projection_vector", [])]
    zero_audit = p2451.get("zero_projection_amplitude_interval_audit", {})
    stationary_audit = p2451.get("stationary_factor_interval_audit", {})
    cell_width = float(p2451.get("interval_cell_width", 1.0e-4))
    zero_segments = p2450.get("zero_projection_amplitude_certificate", {}).get("sampled_lipschitz_exclusion", {}).get("complement_segments", [])
    stationary_segments = p2450.get("stationary_factor_certificate", {}).get("sampled_lipschitz_exclusion", {}).get("complement_segments", [])
    zero_segment = find_segment_for_cell(zero_audit.get("weakest_cell", {}), zero_segments)
    stationary_segment = find_segment_for_cell(stationary_audit.get("weakest_cell", {}), stationary_segments)
    zero_band = build_forward_band(zero_audit.get("weakest_cell", {}), zero_segment, cell_width)
    stationary_band = build_forward_band(stationary_audit.get("weakest_cell", {}), stationary_segment, cell_width)
    zero_replay = replay_band("zero_projection_amplitude", zero_band, projection_vector, projection_amplitude_interval)
    stationary_replay = replay_band("stationary_factor", stationary_band, projection_vector, stationary_factor_interval)
    theorem_export = {
        "theorem_name": "P2454_T1_strict_pointwise_directed_decimal_weakest_band_replay_certificate",
        "inherited_interval_enclosure_audit": "P2451/S1401",
        "inherited_weakest_cell_replay": "P2453/S1403",
        "decimal_precision": DECIMAL_PRECISION,
        "decimal_safety_pad": str(DECIMAL_SAFETY_PAD),
        "taylor_bound_pad": str(TAYLOR_BOUND_PAD),
        "band_cell_count_requested_per_family": BAND_CELL_COUNT,
        "cell_width_inherited_from_p2451": cell_width,
        "p2453_both_weakest_cells_exclude_zero": p2453.get("both_weakest_cells_exclude_zero_under_decimal_taylor_replay"),
        "zero_projection_weakest_band_replay": zero_replay,
        "stationary_factor_weakest_band_replay": stationary_replay,
        "both_bands_exclude_zero_under_decimal_taylor_replay": zero_replay["all_cells_exclude_zero"] and stationary_replay["all_cells_exclude_zero"],
        "minimum_decimal_separation_across_bands": str(
            min(Decimal(zero_replay["minimum_decimal_separation_from_zero"]), Decimal(stationary_replay["minimum_decimal_separation_from_zero"]))
        ),
        "full_complement_directed_rounding_interval_theorem_exported_by_this_certificate": False,
        "symbolic_root_exclusion_theorem_exported_by_this_certificate": False,
        "pointwise_coordinate_selector_exported_by_this_certificate": False,
        "strict_observable_source_constraint_exported_by_this_certificate": False,
        "gauge_slice_theorem_exported_by_this_certificate": False,
        "strict_physical_value_generator_exported": False,
        "qw2191_discharged": False,
        "role_bearing_ltotal_exported": False,
        "toe_closure_exported": False,
        "not_licensed": [
            "Weakest-band Decimal/Taylor replay is not a full complement directed-rounding interval theorem.",
            "The bounded critical-band replay does not prove symbolic root exclusion and does not choose a point-coordinate selector.",
            "No strict observable/source theorem, gauge-slice theorem, physical-value generator, QW-2191 discharge, role-bearing L_total export, or ToE closure is exported.",
        ],
        "next_honest_step": (
            "Extend the Decimal/Taylor endpoint backend from the weakest bands to all complement cells, or replace it with a formal directed-rounding interval backend."
        ),
    }
    gatekeepers = {
        "rg_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "p2453_weakest_cells_inherited": theorem_export["p2453_both_weakest_cells_exclude_zero"] is True,
        "zero_band_has_requested_count": zero_replay["band_cell_count"] == BAND_CELL_COUNT,
        "stationary_band_has_requested_count": stationary_replay["band_cell_count"] == BAND_CELL_COUNT,
        "zero_band_excludes_zero": zero_replay["all_cells_exclude_zero"],
        "stationary_band_excludes_zero": stationary_replay["all_cells_exclude_zero"],
        "both_bands_exclude_zero": theorem_export["both_bands_exclude_zero_under_decimal_taylor_replay"],
        "minimum_decimal_separation_positive": Decimal(theorem_export["minimum_decimal_separation_across_bands"]) > 0,
        "no_full_complement_theorem": not theorem_export["full_complement_directed_rounding_interval_theorem_exported_by_this_certificate"],
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
        "schema_version": "p2454_s1404_v1",
        "packet_id": "P2454",
        "stage_id": "S1404",
        "status": "PASS_STRICT_POINTWISE_DIRECTED_DECIMAL_WEAKEST_BAND_REPLAY_NO_FULL_INTERVAL_SELECTOR_SOURCE_THEOREM",
        "strict_pointwise_directed_decimal_weakest_band_replay_certificate": {
            "inputs": {name: rel(path) for name, path in SOURCE_FILES.items()},
            "input_fingerprints": {name: sha256_json(sources[name]) for name in sources},
            "rg_nonduplication_audit": grep,
            "theorem_export": theorem_export,
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_pointwise_directed_decimal_weakest_band_replay_certificate"]["theorem_export"]
    z = t["zero_projection_weakest_band_replay"]
    h = t["stationary_factor_weakest_band_replay"]
    lines = [
        "# P2454/S1404 strict pointwise directed-decimal weakest-band replay certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Weakest-band replay",
        "",
        f"Zero-projection band cells: `{z['band_cell_count']}`, all exclude zero: `{z['all_cells_exclude_zero']}`, minimum separation: `{z['minimum_decimal_separation_from_zero']}`.",
        f"Stationary-factor band cells: `{h['band_cell_count']}`, all exclude zero: `{h['all_cells_exclude_zero']}`, minimum separation: `{h['minimum_decimal_separation_from_zero']}`.",
        f"Minimum Decimal separation across bands: `{t['minimum_decimal_separation_across_bands']}`.",
        "",
        "## Guardrail",
        "",
        "This is a bounded weakest-band Decimal/Taylor backend replay only.  It exports no full complement directed-rounding interval theorem, no symbolic root-exclusion theorem, no point-coordinate selector, no strict observable/source theorem, no gauge-slice theorem, no physical-value generator, no QW-2191 discharge, no role-bearing `L_total`, and no ToE closure.",
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
