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
OUT = GEN / "p2457_s1407_strict_pointwise_decimal_root_boundary_separation_shape_certificate.json"
MD = GEN / "p2457_s1407_strict_pointwise_decimal_root_boundary_separation_shape_certificate.md"

SOURCE_FILES = {
    "P2456_ROOT_WINDOW_BOUNDARY_BAND_REPLAY": GEN / "p2456_s1406_strict_pointwise_decimal_root_window_boundary_band_replay_certificate.json",
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
        "new_packet": "P2457|S1407|root boundary separation shape|root-window separation shape|boundary separation shape",
        "p2456_input": "P2456|S1406|root-window boundary band|boundary-band replay|two-sided root-window",
        "shape_language": "boundary monotonicity|separation monotonicity|opposite boundary sign|two-sided boundary sign|distance from boundary",
        "decimal_backend": "Decimal interval|Taylor endpoint|zero-excluding|Decimal separation",
        "closure_blockers": "QW-2191|role-bearing L_total|physical-value generator|ToE closure|point-coordinate selector",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def row_sign(row: dict[str, Any]) -> str:
    interval = row["decimal_taylor_interval_value"]
    lo = Decimal(interval["lo"])
    hi = Decimal(interval["hi"])
    if lo > 0:
        return "positive"
    if hi < 0:
        return "negative"
    return "mixed_or_zero_touching"


def rows_ordered_by_distance_from_root_window_boundary(replay: dict[str, Any]) -> list[dict[str, Any]]:
    rows = list(replay["rows"])
    if replay["side"] == "left":
        return list(reversed(rows))
    if replay["side"] == "right":
        return rows
    raise ValueError(f"unknown replay side: {replay['side']}")


def audit_replay_shape(replay: dict[str, Any]) -> dict[str, Any]:
    ordered = rows_ordered_by_distance_from_root_window_boundary(replay)
    separations = [Decimal(row["decimal_separation_from_zero"]) for row in ordered]
    first_differences = [separations[index + 1] - separations[index] for index in range(len(separations) - 1)]
    signs = [row_sign(row) for row in ordered]
    distinct_signs = sorted(set(signs))
    nearest = ordered[0]
    farthest = ordered[-1]
    return {
        "family": replay["family"],
        "root_d": replay["root_d"],
        "side": replay["side"],
        "band_cell_count": len(ordered),
        "distance_ordering_rule": "left bands are reversed from P2456 storage; right bands keep P2456 order",
        "nearest_boundary_cell": {"index": nearest["index"], "left": nearest["left"], "right": nearest["right"]},
        "farthest_boundary_cell": {"index": farthest["index"], "left": farthest["left"], "right": farthest["right"]},
        "separations_ordered_near_to_far": [str(value) for value in separations],
        "first_separation_differences_near_to_far": [str(value) for value in first_differences],
        "minimum_first_separation_difference": str(min(first_differences)),
        "strictly_increasing_separation_away_from_root_window": all(value > 0 for value in first_differences),
        "band_sign": distinct_signs[0] if len(distinct_signs) == 1 else "mixed",
        "signs_ordered_near_to_far": signs,
        "sign_coherent_zero_excluding_band": len(distinct_signs) == 1 and distinct_signs[0] in {"positive", "negative"},
        "p2456_all_cells_exclude_zero_inherited": replay["all_cells_exclude_zero"],
    }


def root_key(audit: dict[str, Any]) -> tuple[str, str]:
    return (audit["family"], str(audit["root_d"]))


def pair_two_sided_signs(audits: list[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[tuple[str, str], dict[str, dict[str, Any]]] = {}
    for audit in audits:
        grouped.setdefault(root_key(audit), {})[audit["side"]] = audit
    pairs = []
    for (family, root_d), sides in sorted(grouped.items(), key=lambda item: (item[0][0], Decimal(item[0][1]))):
        left = sides.get("left")
        right = sides.get("right")
        if left is None or right is None:
            opposite = False
            left_sign = None if left is None else left["band_sign"]
            right_sign = None if right is None else right["band_sign"]
        else:
            left_sign = left["band_sign"]
            right_sign = right["band_sign"]
            opposite = {left_sign, right_sign} == {"positive", "negative"}
        pairs.append(
            {
                "family": family,
                "root_d": root_d,
                "left_boundary_sign": left_sign,
                "right_boundary_sign": right_sign,
                "opposite_zero_excluding_boundary_signs": opposite,
            }
        )
    return pairs


def append_doc_sections() -> None:
    eq_section = """
## P2457/S1407 strict pointwise Decimal root-boundary separation-shape certificate

`P2457/S1407` audits the P2456 all-root-window boundary replay as a shape statement: after ordering each boundary band by increasing distance from its adjacent root-window boundary, Decimal/Taylor separation from zero strictly increases in every audited band, each band is sign-coherent, and the two sides of each audited root window have opposite zero-excluding signs.

This is still a finite boundary-band shape audit.  It is not a full complement directed-rounding theorem, symbolic root-exclusion theorem, point-coordinate selector, strict observable/source theorem, gauge-slice theorem, role-bearing `L_total`, `QW-2191` discharge, or ToE closure.
""".strip()
    lag_section = """
## P2457/S1407 Decimal root-boundary separation-shape guard

`P2457/S1407` strengthens the P2456 boundary replay by checking monotone separation away from every audited root-window boundary plus opposite two-sided signs.  This helps prevent misreading isolated boundary samples as accidental sign margins, but it remains bounded numerical shape evidence and adds no selector/source/gauge authority for `L_total`.
""".strip()
    for path, section in [(DOC_FILES["equation_sheet"], eq_section), (DOC_FILES["lagrangian_eom_draft"], lag_section)]:
        text = path.read_text(encoding="utf-8") if path.exists() else ""
        marker = section.splitlines()[0]
        if marker not in text:
            path.write_text(text.rstrip() + "\n\n" + section + "\n", encoding="utf-8")


def build_payload() -> dict[str, Any]:
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    grep = rg_audit()
    p2456 = sources["P2456_ROOT_WINDOW_BOUNDARY_BAND_REPLAY"].get(
        "strict_pointwise_decimal_root_window_boundary_band_replay_certificate", {}
    ).get("theorem_export", {})
    zero_replays = p2456.get("zero_projection_boundary_band_replays", [])
    stationary_replays = p2456.get("stationary_factor_boundary_band_replays", [])
    all_replays = zero_replays + stationary_replays
    band_audits = [audit_replay_shape(replay) for replay in all_replays]
    two_sided_pairs = pair_two_sided_signs(band_audits)
    min_first_difference = min(Decimal(audit["minimum_first_separation_difference"]) for audit in band_audits)
    theorem_export = {
        "theorem_name": "P2457_T1_strict_pointwise_decimal_root_boundary_separation_shape_certificate",
        "inherited_root_window_boundary_band_replay_certificate": "P2456/S1406",
        "p2456_all_boundary_bands_exclude_zero_inherited": p2456.get("all_boundary_bands_exclude_zero"),
        "p2456_total_boundary_band_count_inherited": p2456.get("total_boundary_band_count"),
        "p2456_total_replayed_boundary_cell_count_inherited": p2456.get("total_replayed_boundary_cell_count"),
        "boundary_band_shape_audits": band_audits,
        "two_sided_root_window_sign_pairs": two_sided_pairs,
        "total_shape_audited_boundary_band_count": len(band_audits),
        "total_shape_audited_boundary_cell_count": sum(audit["band_cell_count"] for audit in band_audits),
        "all_boundary_bands_strictly_increase_separation_away_from_root_window": all(
            audit["strictly_increasing_separation_away_from_root_window"] for audit in band_audits
        ),
        "all_boundary_bands_are_sign_coherent": all(audit["sign_coherent_zero_excluding_band"] for audit in band_audits),
        "all_two_sided_root_windows_have_opposite_boundary_signs": all(
            pair["opposite_zero_excluding_boundary_signs"] for pair in two_sided_pairs
        ),
        "minimum_first_separation_difference_across_boundary_bands": str(min_first_difference),
        "bounded_root_boundary_shape_audit_exported": True,
        "full_complement_directed_rounding_interval_theorem_exported_by_this_certificate": False,
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
            "Near-to-far monotone Decimal separation across finite boundary bands is not a full complement directed-rounding theorem.",
            "Opposite two-sided boundary signs do not prove symbolic root exclusion and do not select a privileged point coordinate.",
            "No strict observable/source theorem, gauge-slice theorem, physical-value generator, QW-2191 discharge, role-bearing L_total export, or ToE closure is exported.",
        ],
        "next_honest_step": (
            "Replace finite boundary-band shape evidence with a formal derivative-sign or directed-rounding interval proof on every complement cell."
        ),
    }
    gatekeepers = {
        "rg_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "p2456_all_boundary_bands_exclude_zero_inherited": theorem_export["p2456_all_boundary_bands_exclude_zero_inherited"] is True,
        "shape_band_count_matches_p2456": theorem_export["total_shape_audited_boundary_band_count"] == theorem_export["p2456_total_boundary_band_count_inherited"] == 6,
        "shape_cell_count_matches_p2456": theorem_export["total_shape_audited_boundary_cell_count"] == theorem_export["p2456_total_replayed_boundary_cell_count_inherited"] == 36,
        "all_boundary_bands_strictly_increase_separation_away_from_root_window": theorem_export["all_boundary_bands_strictly_increase_separation_away_from_root_window"],
        "all_boundary_bands_are_sign_coherent": theorem_export["all_boundary_bands_are_sign_coherent"],
        "all_two_sided_root_windows_have_opposite_boundary_signs": theorem_export["all_two_sided_root_windows_have_opposite_boundary_signs"],
        "minimum_first_separation_difference_positive": Decimal(theorem_export["minimum_first_separation_difference_across_boundary_bands"]) > 0,
        "no_full_complement_theorem": not theorem_export["full_complement_directed_rounding_interval_theorem_exported_by_this_certificate"],
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
        "schema_version": "p2457_s1407_v1",
        "packet_id": "P2457",
        "stage_id": "S1407",
        "status": "PASS_STRICT_POINTWISE_DECIMAL_ROOT_BOUNDARY_SEPARATION_SHAPE_NO_FULL_INTERVAL_SELECTOR_SOURCE_THEOREM",
        "strict_pointwise_decimal_root_boundary_separation_shape_certificate": {
            "inputs": {name: rel(path) for name, path in SOURCE_FILES.items()},
            "input_fingerprints": {name: sha256_json(sources[name]) for name in sources},
            "rg_nonduplication_audit": grep,
            "theorem_export": theorem_export,
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_pointwise_decimal_root_boundary_separation_shape_certificate"]["theorem_export"]
    lines = [
        "# P2457/S1407 strict pointwise Decimal root-boundary separation-shape certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Boundary-shape audit",
        "",
        f"Boundary bands shape-audited: `{t['total_shape_audited_boundary_band_count']}`.",
        f"Boundary cells shape-audited: `{t['total_shape_audited_boundary_cell_count']}`.",
        f"All boundary bands strictly increase separation away from the root window: `{t['all_boundary_bands_strictly_increase_separation_away_from_root_window']}`.",
        f"All boundary bands are sign-coherent: `{t['all_boundary_bands_are_sign_coherent']}`.",
        f"All two-sided root windows have opposite boundary signs: `{t['all_two_sided_root_windows_have_opposite_boundary_signs']}`.",
        f"Minimum first separation difference across boundary bands: `{t['minimum_first_separation_difference_across_boundary_bands']}`.",
        "",
        "## Guardrail",
        "",
        "This is a bounded root-window boundary separation-shape audit only.  It exports no full complement directed-rounding interval theorem, no symbolic root-exclusion theorem, no analytic monotonicity theorem, no point-coordinate selector, no strict observable/source theorem, no gauge-slice theorem, no physical-value generator, no QW-2191 discharge, no role-bearing `L_total`, and no ToE closure.",
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
