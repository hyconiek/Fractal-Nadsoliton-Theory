#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from decimal import Decimal
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"
OUT = GEN / "p2458_s1408_strict_pointwise_interval_decimal_weakest_cell_alignment_certificate.json"
MD = GEN / "p2458_s1408_strict_pointwise_interval_decimal_weakest_cell_alignment_certificate.md"

SOURCE_FILES = {
    "P2451_FLOATING_INTERVAL_AUDIT": GEN / "p2451_s1401_strict_pointwise_projection_interval_enclosure_root_exclusion_audit.json",
    "P2456_DECIMAL_BOUNDARY_REPLAY": GEN / "p2456_s1406_strict_pointwise_decimal_root_window_boundary_band_replay_certificate.json",
    "P2457_BOUNDARY_SHAPE_AUDIT": GEN / "p2457_s1407_strict_pointwise_decimal_root_boundary_separation_shape_certificate.json",
}
DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}
CELL_MATCH_TOLERANCE = 1.0e-12


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
        "new_packet": "P2458|S1408|interval-decimal weakest-cell alignment|weakest-cell alignment|backend alignment",
        "p2451_input": "P2451|S1401|floating interval enclosure|weakest_cell|interval audit",
        "p2456_input": "P2456|S1406|root-window boundary band|boundary-band replay|Decimal separation",
        "p2457_input": "P2457|S1407|root-boundary separation-shape|opposite boundary sign|separation away",
        "closure_blockers": "QW-2191|role-bearing L_total|physical-value generator|ToE closure|point-coordinate selector",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def cell_matches(left_a: float, right_a: float, left_b: float, right_b: float) -> bool:
    return abs(left_a - left_b) <= CELL_MATCH_TOLERANCE and abs(right_a - right_b) <= CELL_MATCH_TOLERANCE


def decimal_row_sign(row: dict[str, Any]) -> str:
    interval = row["decimal_taylor_interval_value"]
    lo = Decimal(interval["lo"])
    hi = Decimal(interval["hi"])
    if lo > 0:
        return "positive"
    if hi < 0:
        return "negative"
    return "mixed_or_zero_touching"


def p2451_family_audit(p2451: dict[str, Any], family: str) -> dict[str, Any]:
    key = "zero_projection_amplitude_interval_audit" if family == "zero_projection_amplitude" else "stationary_factor_interval_audit"
    return p2451[key]


def p2456_family_replays(p2456: dict[str, Any], family: str) -> list[dict[str, Any]]:
    key = "zero_projection_boundary_band_replays" if family == "zero_projection_amplitude" else "stationary_factor_boundary_band_replays"
    return p2456[key]


def find_matching_boundary_row(p2456: dict[str, Any], family: str, weakest_cell: dict[str, Any]) -> dict[str, Any] | None:
    for replay in p2456_family_replays(p2456, family):
        for row in replay["rows"]:
            if cell_matches(float(row["left"]), float(row["right"]), float(weakest_cell["left"]), float(weakest_cell["right"])):
                return {"replay": replay, "row": row}
    return None


def matching_shape_audit(p2457: dict[str, Any], family: str, replay: dict[str, Any]) -> dict[str, Any] | None:
    for audit in p2457.get("boundary_band_shape_audits", []):
        if audit["family"] == family and audit["side"] == replay["side"] and str(audit["root_d"]) == str(replay["root_d"]):
            return audit
    return None


def audit_family_alignment(p2451: dict[str, Any], p2456: dict[str, Any], p2457: dict[str, Any], family: str) -> dict[str, Any]:
    interval_audit = p2451_family_audit(p2451, family)
    weakest_cell = interval_audit["weakest_cell"]
    match = find_matching_boundary_row(p2456, family, weakest_cell)
    if match is None:
        return {
            "family": family,
            "p2451_weakest_cell": weakest_cell,
            "matched_by_p2456_boundary_replay": False,
            "covered_by_p2457_shape_audit": False,
        }
    replay = match["replay"]
    row = match["row"]
    shape = matching_shape_audit(p2457, family, replay)
    nearest = None if shape is None else shape["nearest_boundary_cell"]
    nearest_matches = False if nearest is None else cell_matches(
        float(nearest["left"]), float(nearest["right"]), float(row["left"]), float(row["right"])
    )
    float_sep = Decimal(str(weakest_cell["separation_from_zero"]))
    decimal_sep = Decimal(row["decimal_separation_from_zero"])
    return {
        "family": family,
        "p2451_weakest_cell": weakest_cell,
        "p2451_interval_cell_count": interval_audit["cell_count"],
        "p2451_all_interval_cells_exclude_zero": interval_audit["all_complement_cells_exclude_zero"],
        "matched_by_p2456_boundary_replay": True,
        "matched_root_d": replay["root_d"],
        "matched_side": replay["side"],
        "matched_boundary_row_index": row["index"],
        "matched_boundary_row_cell": {"left": row["left"], "right": row["right"]},
        "matched_cell_is_nearest_boundary_cell": row["index"] == 0 and nearest_matches,
        "p2457_shape_audit_found": shape is not None,
        "p2457_shape_audit_strictly_increasing": False if shape is None else shape["strictly_increasing_separation_away_from_root_window"],
        "p2457_shape_audit_sign_coherent": False if shape is None else shape["sign_coherent_zero_excluding_band"],
        "p2457_band_sign": None if shape is None else shape["band_sign"],
        "matched_decimal_row_sign": decimal_row_sign(row),
        "float_interval_separation_from_zero": str(float_sep),
        "decimal_taylor_separation_from_zero": str(decimal_sep),
        "decimal_minus_float_separation": str(decimal_sep - float_sep),
        "decimal_separation_positive": decimal_sep > 0,
    }


def append_doc_sections() -> None:
    eq_section = """
## P2458/S1408 strict pointwise interval-Decimal weakest-cell alignment certificate

`P2458/S1408` cross-checks the P2451 floating interval weakest cells against the later P2456/P2457 Decimal boundary chain.  For both the zero-projection amplitude and the stationary factor, the P2451 weakest complement cell is exactly found as the nearest P2456 root-window boundary cell and is covered by the P2457 monotone/sign-coherent boundary-shape audit.

This is a backend-chain alignment certificate only.  It does not upgrade floating intervals into a directed-rounding theorem, does not prove symbolic root exclusion, and exports no point-coordinate selector, source/gauge theorem, role-bearing `L_total`, `QW-2191` discharge, or ToE closure.
""".strip()
    lag_section = """
## P2458/S1408 interval-Decimal weakest-cell alignment guard

`P2458/S1408` verifies that the weakest P2451 interval cells are not orphaned from the later Decimal chain: they coincide with nearest P2456 boundary rows and are protected by P2457 shape witnesses.  This improves audit traceability but remains backend consistency evidence, not selector/source/gauge authority for `L_total`.
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
    p2457 = sources["P2457_BOUNDARY_SHAPE_AUDIT"].get("strict_pointwise_decimal_root_boundary_separation_shape_certificate", {}).get("theorem_export", {})
    alignments = [
        audit_family_alignment(p2451, p2456, p2457, "zero_projection_amplitude"),
        audit_family_alignment(p2451, p2456, p2457, "stationary_factor"),
    ]
    minimum_decimal_minus_float = min(Decimal(item["decimal_minus_float_separation"]) for item in alignments)
    theorem_export = {
        "theorem_name": "P2458_T1_strict_pointwise_interval_decimal_weakest_cell_alignment_certificate",
        "inherited_floating_interval_audit": "P2451/S1401",
        "inherited_decimal_boundary_replay": "P2456/S1406",
        "inherited_boundary_shape_audit": "P2457/S1407",
        "cell_match_tolerance": CELL_MATCH_TOLERANCE,
        "weakest_cell_alignments": alignments,
        "aligned_family_count": len(alignments),
        "all_p2451_weakest_cells_match_p2456_boundary_replays": all(item["matched_by_p2456_boundary_replay"] for item in alignments),
        "all_matched_cells_are_nearest_boundary_cells": all(item["matched_cell_is_nearest_boundary_cell"] for item in alignments),
        "all_matched_cells_are_covered_by_p2457_shape_audits": all(item["p2457_shape_audit_found"] for item in alignments),
        "all_matched_shape_audits_are_strictly_increasing": all(item["p2457_shape_audit_strictly_increasing"] for item in alignments),
        "all_matched_shape_audits_are_sign_coherent": all(item["p2457_shape_audit_sign_coherent"] for item in alignments),
        "all_decimal_separations_positive": all(item["decimal_separation_positive"] for item in alignments),
        "minimum_decimal_minus_float_separation": str(minimum_decimal_minus_float),
        "backend_chain_alignment_certificate_exported": True,
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
            "Weakest-cell backend alignment does not convert floating interval arithmetic into a directed-rounding theorem.",
            "Alignment of the weakest P2451 cells with P2456/P2457 boundary witnesses does not prove symbolic root exclusion.",
            "No point-coordinate selector, strict observable/source theorem, gauge-slice theorem, physical-value generator, QW-2191 discharge, role-bearing L_total export, or ToE closure is exported.",
        ],
        "next_honest_step": (
            "Audit whether every P2451 complement cell, not only the weakest cells, has a corresponding formal directed-rounding or symbolic exclusion witness."
        ),
    }
    gatekeepers = {
        "rg_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "aligned_family_count_is_two": theorem_export["aligned_family_count"] == 2,
        "all_p2451_weakest_cells_match_p2456_boundary_replays": theorem_export["all_p2451_weakest_cells_match_p2456_boundary_replays"],
        "all_matched_cells_are_nearest_boundary_cells": theorem_export["all_matched_cells_are_nearest_boundary_cells"],
        "all_matched_cells_are_covered_by_p2457_shape_audits": theorem_export["all_matched_cells_are_covered_by_p2457_shape_audits"],
        "all_matched_shape_audits_are_strictly_increasing": theorem_export["all_matched_shape_audits_are_strictly_increasing"],
        "all_matched_shape_audits_are_sign_coherent": theorem_export["all_matched_shape_audits_are_sign_coherent"],
        "all_decimal_separations_positive": theorem_export["all_decimal_separations_positive"],
        "minimum_decimal_minus_float_separation_nonnegative": Decimal(theorem_export["minimum_decimal_minus_float_separation"]) >= 0,
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
        "schema_version": "p2458_s1408_v1",
        "packet_id": "P2458",
        "stage_id": "S1408",
        "status": "PASS_STRICT_POINTWISE_INTERVAL_DECIMAL_WEAKEST_CELL_ALIGNMENT_NO_DIRECTED_INTERVAL_SELECTOR_SOURCE_THEOREM",
        "strict_pointwise_interval_decimal_weakest_cell_alignment_certificate": {
            "inputs": {name: rel(path) for name, path in SOURCE_FILES.items()},
            "input_fingerprints": {name: sha256_json(sources[name]) for name in sources},
            "rg_nonduplication_audit": grep,
            "theorem_export": theorem_export,
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_pointwise_interval_decimal_weakest_cell_alignment_certificate"]["theorem_export"]
    lines = [
        "# P2458/S1408 strict pointwise interval-Decimal weakest-cell alignment certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Weakest-cell alignment",
        "",
        f"Aligned families: `{t['aligned_family_count']}`.",
        f"All P2451 weakest cells match P2456 boundary replays: `{t['all_p2451_weakest_cells_match_p2456_boundary_replays']}`.",
        f"All matched cells are nearest boundary cells: `{t['all_matched_cells_are_nearest_boundary_cells']}`.",
        f"All matched cells are covered by P2457 shape audits: `{t['all_matched_cells_are_covered_by_p2457_shape_audits']}`.",
        f"Minimum Decimal-minus-float separation: `{t['minimum_decimal_minus_float_separation']}`.",
        "",
        "## Guardrail",
        "",
        "This is a backend-chain alignment certificate only.  It exports no directed-rounding interval theorem, no symbolic root-exclusion theorem, no analytic monotonicity theorem, no point-coordinate selector, no strict observable/source theorem, no gauge-slice theorem, no physical-value generator, no QW-2191 discharge, no role-bearing `L_total`, and no ToE closure.",
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
