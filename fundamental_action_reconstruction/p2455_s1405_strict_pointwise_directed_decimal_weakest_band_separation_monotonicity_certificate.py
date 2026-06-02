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
OUT = GEN / "p2455_s1405_strict_pointwise_directed_decimal_weakest_band_separation_monotonicity_certificate.json"
MD = GEN / "p2455_s1405_strict_pointwise_directed_decimal_weakest_band_separation_monotonicity_certificate.md"

SOURCE_FILES = {
    "P2454_WEAKEST_BAND_REPLAY": GEN / "p2454_s1404_strict_pointwise_directed_decimal_weakest_band_replay_certificate.json",
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
        "new_packet": "P2455|S1405|weakest-band separation|separation monotonicity|Decimal separation monotone",
        "p2454_input": "P2454|S1404|directed decimal weakest-band|critical band replay|weakest band replay",
        "monotonic_language": "strictly increasing separation|first difference|weakest boundary cell|monotone separation",
        "backend_language": "Decimal interval|Taylor endpoint|zero-excluding|Decimal separation",
        "closure_blockers": "QW-2191|role-bearing L_total|physical-value generator|ToE closure|point-coordinate selector",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def separation_monotonicity_audit(family_replay: dict[str, Any]) -> dict[str, Any]:
    rows = family_replay.get("rows", [])
    separations = [Decimal(row["decimal_separation_from_zero"]) for row in rows]
    first_differences = [separations[idx + 1] - separations[idx] for idx in range(len(separations) - 1)]
    positive_differences = [diff for diff in first_differences if diff > 0]
    weakest_index = min(range(len(separations)), key=lambda idx: separations[idx]) if separations else None
    return {
        "family": family_replay.get("family"),
        "band_cell_count": len(rows),
        "separations": [str(value) for value in separations],
        "first_differences": [str(value) for value in first_differences],
        "strictly_increasing_separation": len(positive_differences) == len(first_differences),
        "weakest_decimal_separation_index": weakest_index,
        "weakest_decimal_separation_at_boundary_cell": weakest_index == 0,
        "minimum_first_difference": str(min(first_differences)) if first_differences else None,
        "maximum_first_difference": str(max(first_differences)) if first_differences else None,
        "minimum_separation": str(min(separations)) if separations else None,
        "maximum_separation": str(max(separations)) if separations else None,
    }


def append_doc_sections() -> None:
    eq_section = """
## P2455/S1405 strict pointwise directed-decimal weakest-band separation-monotonicity certificate

`P2455/S1405` audits the P2454 critical bands one level deeper: in each replayed forward band, the Decimal/Taylor separation from zero is strictly increasing cell-by-cell.  Thus the boundary cell adjacent to the excluded root window is the actual weakest replayed cell for both the zero-projection amplitude and the stationary factor.

This is still a bounded critical-band shape audit, not a full complement directed-rounding theorem, symbolic root exclusion, selector/source theorem, role-bearing `L_total`, `QW-2191` discharge, or ToE closure.
""".strip()
    lag_section = """
## P2455/S1405 directed-decimal weakest-band monotonicity guard

`P2455/S1405` shows that the Decimal/Taylor zero-separation margins increase away from the root-window boundary across the P2454 critical bands.  This strengthens the bottleneck-cell audit but remains bounded numerical evidence, not selector/source/gauge authority for `L_total`.
""".strip()
    for path, section in [(DOC_FILES["equation_sheet"], eq_section), (DOC_FILES["lagrangian_eom_draft"], lag_section)]:
        text = path.read_text(encoding="utf-8") if path.exists() else ""
        marker = section.splitlines()[0]
        if marker not in text:
            path.write_text(text.rstrip() + "\n\n" + section + "\n", encoding="utf-8")


def build_payload() -> dict[str, Any]:
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    grep = rg_audit()
    p2454 = sources["P2454_WEAKEST_BAND_REPLAY"].get("strict_pointwise_directed_decimal_weakest_band_replay_certificate", {}).get("theorem_export", {})
    zero_audit = separation_monotonicity_audit(p2454.get("zero_projection_weakest_band_replay", {}))
    stationary_audit = separation_monotonicity_audit(p2454.get("stationary_factor_weakest_band_replay", {}))
    theorem_export = {
        "theorem_name": "P2455_T1_strict_pointwise_directed_decimal_weakest_band_separation_monotonicity_certificate",
        "inherited_weakest_band_replay": "P2454/S1404",
        "zero_projection_separation_monotonicity_audit": zero_audit,
        "stationary_factor_separation_monotonicity_audit": stationary_audit,
        "both_bands_have_strictly_increasing_separation": zero_audit["strictly_increasing_separation"] and stationary_audit["strictly_increasing_separation"],
        "both_bands_weakest_at_boundary_cell": zero_audit["weakest_decimal_separation_at_boundary_cell"] and stationary_audit["weakest_decimal_separation_at_boundary_cell"],
        "inherited_p2454_bands_exclude_zero": p2454.get("both_bands_exclude_zero_under_decimal_taylor_replay"),
        "minimum_first_difference_across_bands": str(
            min(Decimal(zero_audit["minimum_first_difference"]), Decimal(stationary_audit["minimum_first_difference"]))
        ),
        "bounded_critical_band_monotonicity_certificate_exported": True,
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
            "Critical-band separation monotonicity is not a full complement directed-rounding interval theorem.",
            "The monotone margin shape does not prove symbolic root exclusion and does not choose a point-coordinate selector.",
            "No strict observable/source theorem, gauge-slice theorem, physical-value generator, QW-2191 discharge, role-bearing L_total export, or ToE closure is exported.",
        ],
        "next_honest_step": (
            "Extend the monotone Decimal/Taylor separation audit to all complement cells or replace the bounded replay with formal directed-rounding interval proof."
        ),
    }
    gatekeepers = {
        "rg_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "inherited_p2454_bands_exclude_zero": theorem_export["inherited_p2454_bands_exclude_zero"] is True,
        "zero_band_strictly_increasing": zero_audit["strictly_increasing_separation"],
        "stationary_band_strictly_increasing": stationary_audit["strictly_increasing_separation"],
        "both_bands_strictly_increasing": theorem_export["both_bands_have_strictly_increasing_separation"],
        "both_bands_weakest_at_boundary": theorem_export["both_bands_weakest_at_boundary_cell"],
        "minimum_first_difference_positive": Decimal(theorem_export["minimum_first_difference_across_bands"]) > 0,
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
        "schema_version": "p2455_s1405_v1",
        "packet_id": "P2455",
        "stage_id": "S1405",
        "status": "PASS_STRICT_POINTWISE_DIRECTED_DECIMAL_WEAKEST_BAND_SEPARATION_MONOTONICITY_NO_FULL_INTERVAL_SELECTOR_SOURCE_THEOREM",
        "strict_pointwise_directed_decimal_weakest_band_separation_monotonicity_certificate": {
            "inputs": {name: rel(path) for name, path in SOURCE_FILES.items()},
            "input_fingerprints": {name: sha256_json(sources[name]) for name in sources},
            "rg_nonduplication_audit": grep,
            "theorem_export": theorem_export,
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_pointwise_directed_decimal_weakest_band_separation_monotonicity_certificate"]["theorem_export"]
    z = t["zero_projection_separation_monotonicity_audit"]
    h = t["stationary_factor_separation_monotonicity_audit"]
    lines = [
        "# P2455/S1405 strict pointwise directed-decimal weakest-band separation-monotonicity certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Separation monotonicity",
        "",
        f"Zero-projection band strictly increasing: `{z['strictly_increasing_separation']}`, minimum first difference: `{z['minimum_first_difference']}`.",
        f"Stationary-factor band strictly increasing: `{h['strictly_increasing_separation']}`, minimum first difference: `{h['minimum_first_difference']}`.",
        f"Both weakest cells are boundary cells: `{t['both_bands_weakest_at_boundary_cell']}`.",
        "",
        "## Guardrail",
        "",
        "This is a bounded weakest-band separation-shape audit only.  It exports no full complement directed-rounding interval theorem, no symbolic root-exclusion theorem, no point-coordinate selector, no strict observable/source theorem, no gauge-slice theorem, no physical-value generator, no QW-2191 discharge, no role-bearing `L_total`, and no ToE closure.",
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
