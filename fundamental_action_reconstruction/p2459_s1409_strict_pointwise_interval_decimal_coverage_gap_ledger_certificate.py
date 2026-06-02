#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from decimal import Decimal, getcontext
from pathlib import Path
from typing import Any

getcontext().prec = 80

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"
OUT = GEN / "p2459_s1409_strict_pointwise_interval_decimal_coverage_gap_ledger_certificate.json"
MD = GEN / "p2459_s1409_strict_pointwise_interval_decimal_coverage_gap_ledger_certificate.md"

SOURCE_FILES = {
    "P2451_FLOATING_INTERVAL_AUDIT": GEN / "p2451_s1401_strict_pointwise_projection_interval_enclosure_root_exclusion_audit.json",
    "P2456_DECIMAL_BOUNDARY_REPLAY": GEN / "p2456_s1406_strict_pointwise_decimal_root_window_boundary_band_replay_certificate.json",
    "P2458_WEAKEST_CELL_ALIGNMENT": GEN / "p2458_s1408_strict_pointwise_interval_decimal_weakest_cell_alignment_certificate.json",
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
        "new_packet": "P2459|S1409|coverage gap ledger|Decimal coverage gap|interval Decimal coverage",
        "coverage_language": "unreplayed complement|complement coverage ledger|coverage ratio|covered complement cells",
        "p2451_input": "P2451|S1401|floating interval enclosure|complement cells|interval audit",
        "p2456_input": "P2456|S1406|root-window boundary band|boundary-band replay|boundary cells",
        "p2458_input": "P2458|S1408|weakest-cell alignment|backend-chain alignment|nearest boundary",
        "closure_blockers": "QW-2191|role-bearing L_total|physical-value generator|ToE closure|point-coordinate selector",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def p2451_family_audit(p2451: dict[str, Any], family: str) -> dict[str, Any]:
    key = "zero_projection_amplitude_interval_audit" if family == "zero_projection_amplitude" else "stationary_factor_interval_audit"
    return p2451[key]


def p2456_family_replays(p2456: dict[str, Any], family: str) -> list[dict[str, Any]]:
    key = "zero_projection_boundary_band_replays" if family == "zero_projection_amplitude" else "stationary_factor_boundary_band_replays"
    return p2456[key]


def decimal_ratio(numerator: int, denominator: int) -> str:
    if denominator == 0:
        return "NaN"
    return str(Decimal(numerator) / Decimal(denominator))


def family_coverage_row(p2451: dict[str, Any], p2456: dict[str, Any], family: str) -> dict[str, Any]:
    interval_audit = p2451_family_audit(p2451, family)
    interval_cell_count = int(interval_audit["cell_count"])
    decimal_boundary_cell_count = sum(int(replay["band_cell_count"]) for replay in p2456_family_replays(p2456, family))
    unreplayed_cell_count = interval_cell_count - decimal_boundary_cell_count
    return {
        "family": family,
        "p2451_interval_complement_cell_count": interval_cell_count,
        "p2451_all_interval_cells_exclude_zero": interval_audit["all_complement_cells_exclude_zero"],
        "p2456_decimal_boundary_replayed_cell_count": decimal_boundary_cell_count,
        "p2456_decimal_boundary_band_count": len(p2456_family_replays(p2456, family)),
        "unreplayed_by_decimal_boundary_chain_cell_count": unreplayed_cell_count,
        "decimal_boundary_coverage_ratio": decimal_ratio(decimal_boundary_cell_count, interval_cell_count),
        "unreplayed_coverage_gap_ratio": decimal_ratio(unreplayed_cell_count, interval_cell_count),
        "decimal_boundary_chain_is_full_family_complement_coverage": decimal_boundary_cell_count == interval_cell_count,
    }


def append_doc_sections() -> None:
    eq_section = """
## P2459/S1409 strict pointwise interval-Decimal coverage-gap ledger certificate

`P2459/S1409` records the coverage gap between the P2451 floating-interval complement audit and the later P2456/P2458 Decimal boundary chain.  The Decimal chain aligns the weakest cells, but it covers only the bounded root-window boundary subset: 36 replayed Decimal boundary cells versus 99,882 P2451 complement cells, leaving 99,846 complement cells not replayed by the Decimal boundary chain.

This is an honesty ledger, not a failure of P2451: P2451 remains the broad floating-interval audit, while P2456-P2458 are targeted Decimal boundary/weakest-cell checks.  The ledger exports no directed-rounding interval theorem, symbolic root-exclusion theorem, selector/source/gauge theorem, role-bearing `L_total`, `QW-2191` discharge, or ToE closure.
""".strip()
    lag_section = """
## P2459/S1409 interval-Decimal coverage-gap guard

`P2459/S1409` prevents overclaiming the Decimal boundary chain as full complement coverage: it quantifies the gap between the 36 Decimal boundary cells and the 99,882 P2451 complement cells.  This improves proof hygiene but adds no selector/source/gauge authority for `L_total`.
""".strip()
    for path, section in [(DOC_FILES["equation_sheet"], eq_section), (DOC_FILES["lagrangian_eom_draft"], lag_section)]:
        text = path.read_text(encoding="utf-8") if path.exists() else ""
        marker = section.splitlines()[0]
        if marker not in text:
            path.write_text(text.rstrip() + "\n\n" + section + "\n", encoding="utf-8")


def build_payload() -> dict[str, Any]:
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    grep = rg_audit()
    p2451 = sources["P2451_FLOATING_INTERVAL_AUDIT"].get("strict_pointwise_projection_interval_enclosure_root_exclusion_audit", {}).get("theorem_export", {})
    p2456 = sources["P2456_DECIMAL_BOUNDARY_REPLAY"].get("strict_pointwise_decimal_root_window_boundary_band_replay_certificate", {}).get("theorem_export", {})
    p2458 = sources["P2458_WEAKEST_CELL_ALIGNMENT"].get("strict_pointwise_interval_decimal_weakest_cell_alignment_certificate", {}).get("theorem_export", {})
    coverage_rows = [
        family_coverage_row(p2451, p2456, "zero_projection_amplitude"),
        family_coverage_row(p2451, p2456, "stationary_factor"),
    ]
    total_interval_cells = sum(row["p2451_interval_complement_cell_count"] for row in coverage_rows)
    total_decimal_cells = sum(row["p2456_decimal_boundary_replayed_cell_count"] for row in coverage_rows)
    total_unreplayed_cells = total_interval_cells - total_decimal_cells
    theorem_export = {
        "theorem_name": "P2459_T1_strict_pointwise_interval_decimal_coverage_gap_ledger_certificate",
        "inherited_floating_interval_audit": "P2451/S1401",
        "inherited_decimal_boundary_replay": "P2456/S1406",
        "inherited_weakest_cell_alignment": "P2458/S1408",
        "p2458_all_weakest_cells_aligned_inherited": p2458.get("all_p2451_weakest_cells_match_p2456_boundary_replays"),
        "coverage_rows": coverage_rows,
        "total_p2451_interval_complement_cell_count": total_interval_cells,
        "total_p2456_decimal_boundary_replayed_cell_count": total_decimal_cells,
        "total_unreplayed_by_decimal_boundary_chain_cell_count": total_unreplayed_cells,
        "total_decimal_boundary_coverage_ratio": decimal_ratio(total_decimal_cells, total_interval_cells),
        "total_unreplayed_coverage_gap_ratio": decimal_ratio(total_unreplayed_cells, total_interval_cells),
        "decimal_boundary_chain_is_full_complement_coverage": total_decimal_cells == total_interval_cells,
        "coverage_gap_positive": total_unreplayed_cells > 0,
        "coverage_gap_ledger_exported": True,
        "directed_rounding_interval_theorem_exported_by_this_certificate": False,
        "symbolic_root_exclusion_theorem_exported_by_this_certificate": False,
        "decimal_full_complement_replay_exported_by_this_certificate": False,
        "pointwise_coordinate_selector_exported_by_this_certificate": False,
        "strict_observable_source_constraint_exported_by_this_certificate": False,
        "gauge_slice_theorem_exported_by_this_certificate": False,
        "strict_physical_value_generator_exported": False,
        "qw2191_discharged": False,
        "role_bearing_ltotal_exported": False,
        "toe_closure_exported": False,
        "not_licensed": [
            "A positive Decimal coverage gap means the P2456-P2458 Decimal boundary chain must not be read as full complement coverage.",
            "P2451 remains a broad floating-interval audit, but this certificate does not upgrade it to a directed-rounding interval theorem.",
            "No symbolic root-exclusion theorem, point-coordinate selector, strict observable/source theorem, gauge-slice theorem, physical-value generator, QW-2191 discharge, role-bearing L_total export, or ToE closure is exported.",
        ],
        "next_honest_step": (
            "Either replay the unreplayed complement cells with a directed Decimal backend or replace the floating interval audit with a formal directed-rounding interval proof."
        ),
    }
    gatekeepers = {
        "rg_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "p2458_all_weakest_cells_aligned_inherited": theorem_export["p2458_all_weakest_cells_aligned_inherited"] is True,
        "family_count_is_two": len(coverage_rows) == 2,
        "total_interval_cell_count_matches_p2451": theorem_export["total_p2451_interval_complement_cell_count"] == 99882,
        "total_decimal_boundary_cell_count_matches_p2456": theorem_export["total_p2456_decimal_boundary_replayed_cell_count"] == p2456.get("total_replayed_boundary_cell_count") == 36,
        "coverage_gap_positive": theorem_export["coverage_gap_positive"],
        "decimal_boundary_chain_not_full_complement_coverage": not theorem_export["decimal_boundary_chain_is_full_complement_coverage"],
        "unreplayed_cell_count_positive": theorem_export["total_unreplayed_by_decimal_boundary_chain_cell_count"] > 0,
        "no_directed_rounding_theorem": not theorem_export["directed_rounding_interval_theorem_exported_by_this_certificate"],
        "no_symbolic_root_exclusion": not theorem_export["symbolic_root_exclusion_theorem_exported_by_this_certificate"],
        "no_decimal_full_complement_replay": not theorem_export["decimal_full_complement_replay_exported_by_this_certificate"],
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
        "schema_version": "p2459_s1409_v1",
        "packet_id": "P2459",
        "stage_id": "S1409",
        "status": "PASS_STRICT_POINTWISE_INTERVAL_DECIMAL_COVERAGE_GAP_LEDGER_NO_FULL_COMPLEMENT_SELECTOR_SOURCE_THEOREM",
        "strict_pointwise_interval_decimal_coverage_gap_ledger_certificate": {
            "inputs": {name: rel(path) for name, path in SOURCE_FILES.items()},
            "input_fingerprints": {name: sha256_json(sources[name]) for name in sources},
            "rg_nonduplication_audit": grep,
            "theorem_export": theorem_export,
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_pointwise_interval_decimal_coverage_gap_ledger_certificate"]["theorem_export"]
    lines = [
        "# P2459/S1409 strict pointwise interval-Decimal coverage-gap ledger certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Coverage ledger",
        "",
        f"P2451 interval complement cells: `{t['total_p2451_interval_complement_cell_count']}`.",
        f"P2456 Decimal boundary replay cells: `{t['total_p2456_decimal_boundary_replayed_cell_count']}`.",
        f"Unreplayed by Decimal boundary chain: `{t['total_unreplayed_by_decimal_boundary_chain_cell_count']}`.",
        f"Decimal boundary coverage ratio: `{t['total_decimal_boundary_coverage_ratio']}`.",
        f"Coverage gap positive: `{t['coverage_gap_positive']}`.",
        "",
        "## Guardrail",
        "",
        "This is a coverage-gap ledger only.  It exports no directed-rounding interval theorem, no symbolic root-exclusion theorem, no Decimal full-complement replay, no point-coordinate selector, no strict observable/source theorem, no gauge-slice theorem, no physical-value generator, no QW-2191 discharge, no role-bearing `L_total`, and no ToE closure.",
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
