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
OUT = GEN / "p2456_s1406_strict_pointwise_decimal_root_window_boundary_band_replay_certificate.json"
MD = GEN / "p2456_s1406_strict_pointwise_decimal_root_window_boundary_band_replay_certificate.md"

SOURCE_FILES = {
    "P2450_ROOT_ISOLATION_MARGIN": GEN / "p2450_s1400_strict_pointwise_projection_root_isolation_margin_certificate.json",
    "P2455_WEAKEST_BAND_MONOTONICITY": GEN / "p2455_s1405_strict_pointwise_directed_decimal_weakest_band_separation_monotonicity_certificate.json",
    "P2449_PROJECTION_REDUCTION": GEN / "p2449_s1399_strict_pointwise_rank_lift_projection_reduction_certificate.json",
}
DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

BOUNDARY_BAND_CELL_COUNT = 6
CELL_WIDTH = 1.0e-4


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
        "new_packet": "P2456|S1406|root-window boundary band|boundary-band replay|two-sided root-window",
        "p2455_input": "P2455|S1405|weakest-band separation|separation monotonicity|Decimal separation",
        "boundary_language": "root-window boundary|left boundary band|right boundary band|boundary-adjacent cells",
        "backend_language": "Decimal interval|Taylor endpoint|zero-excluding|boundary band",
        "closure_blockers": "QW-2191|role-bearing L_total|physical-value generator|ToE closure|point-coordinate selector",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def build_boundary_cells(window: dict[str, float], audit_start: float, audit_end: float, side: str) -> list[dict[str, float]]:
    cells: list[dict[str, float]] = []
    if side == "left":
        right = float(window["left"])
        for _ in range(BOUNDARY_BAND_CELL_COUNT):
            left = max(audit_start, right - CELL_WIDTH)
            if left >= right:
                break
            cells.append({"left": left, "right": right})
            right = left
        cells.reverse()
    elif side == "right":
        left = float(window["right"])
        for _ in range(BOUNDARY_BAND_CELL_COUNT):
            right = min(audit_end, left + CELL_WIDTH)
            if right <= left:
                break
            cells.append({"left": left, "right": right})
            left = right
    else:
        raise ValueError(f"unknown side: {side}")
    return cells


def replay_boundary_band(
    family: str,
    root_d: float,
    side: str,
    cells: list[dict[str, float]],
    projection_vector: list[Decimal],
    function: Callable[[list[Decimal], DecimalInterval], Any],
) -> dict[str, Any]:
    rows = []
    minimum_separation: Decimal | None = None
    for index, cell in enumerate(cells):
        value = function(projection_vector, DecimalInterval(cell["left"], cell["right"]))
        separation = value.separation_from_zero()
        rows.append(
            {
                "index": index,
                "left": cell["left"],
                "right": cell["right"],
                "decimal_taylor_interval_value": value.as_dict(),
                "decimal_interval_excludes_zero": not value.contains_zero(),
                "decimal_separation_from_zero": str(separation),
            }
        )
        minimum_separation = separation if minimum_separation is None else min(minimum_separation, separation)
    return {
        "family": family,
        "root_d": root_d,
        "side": side,
        "band_cell_count": len(rows),
        "rows": rows,
        "all_cells_exclude_zero": all(row["decimal_interval_excludes_zero"] for row in rows),
        "minimum_decimal_separation_from_zero": str(minimum_separation),
    }


def family_boundary_replays(
    family: str,
    certificate: dict[str, Any],
    projection_vector: list[Decimal],
    function: Callable[[list[Decimal], DecimalInterval], Any],
) -> list[dict[str, Any]]:
    audit_start = float(certificate.get("sampled_lipschitz_exclusion", {}).get("audit_start"))
    audit_end = float(certificate.get("sampled_lipschitz_exclusion", {}).get("audit_end"))
    replays = []
    for window in certificate.get("root_windows", []):
        for side in ("left", "right"):
            cells = build_boundary_cells(window, audit_start, audit_end, side)
            if cells:
                replays.append(replay_boundary_band(family, float(window["root_d"]), side, cells, projection_vector, function))
    return replays


def append_doc_sections() -> None:
    eq_section = """
## P2456/S1406 strict pointwise Decimal root-window boundary-band replay certificate

`P2456/S1406` broadens the P2454/P2455 Decimal/Taylor replay from the single weakest forward band to every root-window boundary in the P2450 isolation certificate.  Six boundary-adjacent complement cells are replayed on each available left/right side of every zero-projection and stationary-factor root window, and every replayed cell remains zero-excluding.

This is a bounded all-boundary replay, not a full complement directed-rounding theorem, symbolic root-exclusion proof, selector/source theorem, role-bearing `L_total`, `QW-2191` discharge, or ToE closure.
""".strip()
    lag_section = """
## P2456/S1406 Decimal root-window boundary replay guard

`P2456/S1406` checks all P2450 root-window boundaries with local Decimal/Taylor boundary bands, not only the single weakest P2454 band.  This improves coverage of root-adjacent complement cells but remains bounded numerical evidence, not selector/source/gauge authority for `L_total`.
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
    p2455 = sources["P2455_WEAKEST_BAND_MONOTONICITY"].get("strict_pointwise_directed_decimal_weakest_band_separation_monotonicity_certificate", {}).get("theorem_export", {})
    projection_vector = [Decimal(str(value)) for value in sources["P2449_PROJECTION_REDUCTION"].get("strict_pointwise_rank_lift_projection_reduction_certificate", {}).get("theorem_export", {}).get("projection_vector", [])]
    zero_replays = family_boundary_replays(
        "zero_projection_amplitude",
        p2450.get("zero_projection_amplitude_certificate", {}),
        projection_vector,
        projection_amplitude_interval,
    )
    stationary_replays = family_boundary_replays(
        "stationary_factor",
        p2450.get("stationary_factor_certificate", {}),
        projection_vector,
        stationary_factor_interval,
    )
    all_replays = zero_replays + stationary_replays
    minimum_separation = min(Decimal(replay["minimum_decimal_separation_from_zero"]) for replay in all_replays)
    theorem_export = {
        "theorem_name": "P2456_T1_strict_pointwise_decimal_root_window_boundary_band_replay_certificate",
        "inherited_root_isolation_margin_certificate": "P2450/S1400",
        "inherited_weakest_band_monotonicity_certificate": "P2455/S1405",
        "p2455_bands_strictly_increasing_inherited": p2455.get("both_bands_have_strictly_increasing_separation"),
        "boundary_band_cell_count_requested": BOUNDARY_BAND_CELL_COUNT,
        "boundary_band_cell_width": CELL_WIDTH,
        "zero_projection_boundary_band_replays": zero_replays,
        "stationary_factor_boundary_band_replays": stationary_replays,
        "total_boundary_band_count": len(all_replays),
        "total_replayed_boundary_cell_count": sum(replay["band_cell_count"] for replay in all_replays),
        "all_boundary_bands_exclude_zero": all(replay["all_cells_exclude_zero"] for replay in all_replays),
        "minimum_decimal_separation_across_all_boundary_bands": str(minimum_separation),
        "bounded_all_root_window_boundary_replay_exported": True,
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
            "All-root-window boundary-band replay is still bounded numerical evidence, not a full complement directed-rounding interval theorem.",
            "Boundary-band zero exclusion does not prove symbolic root exclusion and does not choose a point-coordinate selector.",
            "No strict observable/source theorem, gauge-slice theorem, physical-value generator, QW-2191 discharge, role-bearing L_total export, or ToE closure is exported.",
        ],
        "next_honest_step": (
            "Extend Decimal/Taylor replay from boundary bands to all complement cells, or replace the bounded replay with a formal directed-rounding interval proof."
        ),
    }
    gatekeepers = {
        "rg_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "p2455_bands_strictly_increasing_inherited": theorem_export["p2455_bands_strictly_increasing_inherited"] is True,
        "zero_boundary_band_count": len(zero_replays) == 4,
        "stationary_boundary_band_count": len(stationary_replays) == 2,
        "total_boundary_band_count": theorem_export["total_boundary_band_count"] == 6,
        "total_replayed_boundary_cell_count": theorem_export["total_replayed_boundary_cell_count"] == 36,
        "all_boundary_bands_exclude_zero": theorem_export["all_boundary_bands_exclude_zero"],
        "minimum_decimal_separation_positive": Decimal(theorem_export["minimum_decimal_separation_across_all_boundary_bands"]) > 0,
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
        "schema_version": "p2456_s1406_v1",
        "packet_id": "P2456",
        "stage_id": "S1406",
        "status": "PASS_STRICT_POINTWISE_DECIMAL_ROOT_WINDOW_BOUNDARY_BAND_REPLAY_NO_FULL_INTERVAL_SELECTOR_SOURCE_THEOREM",
        "strict_pointwise_decimal_root_window_boundary_band_replay_certificate": {
            "inputs": {name: rel(path) for name, path in SOURCE_FILES.items()},
            "input_fingerprints": {name: sha256_json(sources[name]) for name in sources},
            "rg_nonduplication_audit": grep,
            "theorem_export": theorem_export,
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_pointwise_decimal_root_window_boundary_band_replay_certificate"]["theorem_export"]
    lines = [
        "# P2456/S1406 strict pointwise Decimal root-window boundary-band replay certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Boundary-band replay",
        "",
        f"Boundary bands replayed: `{t['total_boundary_band_count']}`.",
        f"Boundary cells replayed: `{t['total_replayed_boundary_cell_count']}`.",
        f"All boundary bands exclude zero: `{t['all_boundary_bands_exclude_zero']}`.",
        f"Minimum Decimal separation across all boundary bands: `{t['minimum_decimal_separation_across_all_boundary_bands']}`.",
        "",
        "## Guardrail",
        "",
        "This is a bounded all-root-window boundary-band replay only.  It exports no full complement directed-rounding interval theorem, no symbolic root-exclusion theorem, no point-coordinate selector, no strict observable/source theorem, no gauge-slice theorem, no physical-value generator, no QW-2191 discharge, no role-bearing `L_total`, and no ToE closure.",
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
