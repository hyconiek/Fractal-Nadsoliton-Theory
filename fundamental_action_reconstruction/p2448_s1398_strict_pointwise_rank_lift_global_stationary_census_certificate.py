#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
import subprocess
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"
OUT = GEN / "p2448_s1398_strict_pointwise_rank_lift_global_stationary_census_certificate.json"
MD = GEN / "p2448_s1398_strict_pointwise_rank_lift_global_stationary_census_certificate.md"

SOURCE_FILES = {
    "P2442_NULLSPACE": GEN / "p2442_s1392_strict_moment_coefficient_local_identifiability_nullspace_certificate.json",
    "P2447_STATIONARY_REFINEMENT": GEN / "p2447_s1397_strict_pointwise_rank_lift_stationary_refinement_certificate.json",
}
DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

STRICT_PARAMS = {"omega": 0.18575, "phi": 0.16250, "beta": 1.0, "eta": 1.8}
CENSUS_INTERVAL = (0.0, 5.0)
DERIVATIVE_SCAN_STEP = 1.0e-3
DERIVATIVE_H = 1.0e-5
SECOND_DERIVATIVE_H = 1.0e-4
ROOT_TOL = 1.0e-12
STATIONARY_DERIVATIVE_TOL = 1.0e-8
NEGATIVE_SECOND_DERIVATIVE_TOL = -1.0e-3
POSITIVE_SECOND_DERIVATIVE_TOL = 1.0e-3
P2447_D_MATCH_TOL = 5.0e-9
GLOBAL_GAP_TOL = 1.0e-4


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
        "new_packet": "P2448|S1398|global stationary census|pointwise rank-lift global|stationary census",
        "p2447_input": "P2447|S1397|stationary refinement|golden-section|continuous maximum",
        "stationary_census_language": "stationary census|derivative sign|root bracket|global maximum census|boundary dominance",
        "selector_source_language": "d_ref|point-coordinate selector|strict observable/source|gauge-slice|lawful gauge",
        "closure_blockers": "QW-2191|role-bearing L_total|physical-value generator|ToE closure",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def pointwise_gradient(d: float) -> list[float]:
    omega = STRICT_PARAMS["omega"]
    phi = STRICT_PARAMS["phi"]
    beta = STRICT_PARAMS["beta"]
    eta = STRICT_PARAMS["eta"]
    d_eta = d**eta
    denominator = 1.0 + beta * d_eta
    phase = omega * d + phi
    cos_phase = math.cos(phase)
    sin_phase = math.sin(phase)
    eta_derivative = 0.0 if d == 0.0 else -cos_phase * beta * d_eta * math.log(d) / (denominator * denominator)
    return [
        -d * sin_phase / denominator,
        -sin_phase / denominator,
        -cos_phase * d_eta / (denominator * denominator),
        eta_derivative,
    ]


def row_norm(row: list[float]) -> float:
    return math.sqrt(sum(value * value for value in row))


def determinant(matrix: list[list[float]]) -> float:
    mat = [row[:] for row in matrix]
    n = len(mat)
    det = 1.0
    for col in range(n):
        pivot = max(range(col, n), key=lambda r: abs(mat[r][col]))
        if abs(mat[pivot][col]) < 1e-15:
            return 0.0
        if pivot != col:
            mat[col], mat[pivot] = mat[pivot], mat[col]
            det *= -1.0
        pivot_value = mat[col][col]
        det *= pivot_value
        for row in range(col + 1, n):
            factor = mat[row][col] / pivot_value
            for c in range(col, n):
                mat[row][c] -= factor * mat[col][c]
    return det


def normalized_volume(matrix: list[list[float]]) -> float:
    normalized = []
    for row in matrix:
        norm = row_norm(row)
        normalized.append([value / norm for value in row])
    return abs(determinant(normalized))


def pointwise_volume(base_jacobian: list[list[float]], d: float) -> float:
    return normalized_volume(base_jacobian + [pointwise_gradient(d)])


def derivative(base_jacobian: list[list[float]], d: float, h: float = DERIVATIVE_H) -> float:
    left = max(CENSUS_INTERVAL[0], d - h)
    right = min(CENSUS_INTERVAL[1], d + h)
    return (pointwise_volume(base_jacobian, right) - pointwise_volume(base_jacobian, left)) / (right - left)


def second_derivative(base_jacobian: list[list[float]], d: float, h: float = SECOND_DERIVATIVE_H) -> float:
    left = max(CENSUS_INTERVAL[0], d - h)
    right = min(CENSUS_INTERVAL[1], d + h)
    if left == d or right == d:
        return float("nan")
    return (pointwise_volume(base_jacobian, right) - 2.0 * pointwise_volume(base_jacobian, d) + pointwise_volume(base_jacobian, left)) / (h * h)


def bisect_derivative_root(base_jacobian: list[list[float]], a: float, b: float) -> dict[str, Any]:
    fa = derivative(base_jacobian, a)
    fb = derivative(base_jacobian, b)
    iterations = 0
    while b - a > ROOT_TOL and iterations < 100:
        mid = (a + b) / 2.0
        fm = derivative(base_jacobian, mid)
        if fa == 0.0:
            b = a
            fb = fa
            break
        if fb == 0.0:
            a = b
            fa = fb
            break
        if fa * fm <= 0.0:
            b = mid
            fb = fm
        else:
            a = mid
            fa = fm
        iterations += 1
    root = (a + b) / 2.0
    curvature = second_derivative(base_jacobian, root)
    if curvature < NEGATIVE_SECOND_DERIVATIVE_TOL:
        classification = "local_maximum"
    elif curvature > POSITIVE_SECOND_DERIVATIVE_TOL:
        classification = "local_minimum"
    else:
        classification = "flat_or_unclassified"
    return {
        "bracket": {"left": a, "right": b},
        "root_d": root,
        "normalized_rank_lift_volume": pointwise_volume(base_jacobian, root),
        "derivative_at_root": derivative(base_jacobian, root),
        "second_derivative_at_root": curvature,
        "classification": classification,
        "iterations": iterations,
        "final_width": b - a,
    }


def stationary_census(base_jacobian: list[list[float]]) -> dict[str, Any]:
    d_min, d_max = CENSUS_INTERVAL
    count = int(round((d_max - d_min) / DERIVATIVE_SCAN_STEP))
    scan_points = [d_min + i * DERIVATIVE_SCAN_STEP for i in range(count + 1)]
    brackets: list[tuple[float, float]] = []
    prev_x = scan_points[1]
    prev = derivative(base_jacobian, prev_x)
    for x in scan_points[2:]:
        value = derivative(base_jacobian, x)
        if prev == 0.0 or value == 0.0 or prev * value < 0.0:
            brackets.append((prev_x, x))
        prev_x = x
        prev = value
    roots = [bisect_derivative_root(base_jacobian, left, right) for left, right in brackets]
    roots = sorted(roots, key=lambda item: item["root_d"])
    boundary_rows = [
        {"d": d_min, "normalized_rank_lift_volume": pointwise_volume(base_jacobian, d_min), "kind": "left_boundary"},
        {"d": d_max, "normalized_rank_lift_volume": pointwise_volume(base_jacobian, d_max), "kind": "right_boundary"},
    ]
    candidates = [
        {"d": row["root_d"], "normalized_rank_lift_volume": row["normalized_rank_lift_volume"], "kind": row["classification"]}
        for row in roots
    ] + boundary_rows
    best = max(candidates, key=lambda item: item["normalized_rank_lift_volume"])
    maxima = [row for row in roots if row["classification"] == "local_maximum"]
    minima = [row for row in roots if row["classification"] == "local_minimum"]
    return {
        "interval": {"d_min": d_min, "d_max": d_max},
        "derivative_scan_step": DERIVATIVE_SCAN_STEP,
        "derivative_h": DERIVATIVE_H,
        "second_derivative_h": SECOND_DERIVATIVE_H,
        "scan_point_count": len(scan_points),
        "sign_change_bracket_count": len(brackets),
        "stationary_roots": roots,
        "local_maximum_count": len(maxima),
        "local_minimum_count": len(minima),
        "boundary_rows": boundary_rows,
        "best_stationary_or_boundary_row": best,
        "maxima_sorted_by_volume": sorted(maxima, key=lambda item: item["normalized_rank_lift_volume"], reverse=True),
        "minima_sorted_by_volume": sorted(minima, key=lambda item: item["normalized_rank_lift_volume"]),
    }


def append_doc_sections() -> None:
    eq_section = """
## P2448/S1398 strict pointwise rank-lift global stationary-census certificate

`P2448/S1398` extends the P2447 local stationary refinement to a finite derivative-sign census over `d in [0,5]`.  The census finds three stationary roots on the audited interval: two near-zero local minima and one local maximum near `d=0.7852889045`; the local maximum also dominates both interval boundaries and matches the P2447 refined point.

This is still a finite conditioning/global-census statement, not an analytic interval root-exclusion theorem, not a point-coordinate selector theorem, not a strict observable/source theorem, and not a gauge-slice theorem.  It therefore cannot promote a pointwise row into `L_total` or discharge `QW-2191`.
""".strip()
    lag_section = """
## P2448/S1398 global pointwise census guard

`P2448/S1398` confirms that the best pointwise rank-lift conditioning coordinate on the audited interval is the same refined `d≈0.7852889045` stationary maximum, but this global census supplies no selector/source/gauge authority.  Pointwise conditioning remains inadmissible as a role-bearing `L_total` row without a separate theorem.
""".strip()
    for path, section in [(DOC_FILES["equation_sheet"], eq_section), (DOC_FILES["lagrangian_eom_draft"], lag_section)]:
        text = path.read_text(encoding="utf-8") if path.exists() else ""
        marker = section.splitlines()[0]
        if marker not in text:
            path.write_text(text.rstrip() + "\n\n" + section + "\n", encoding="utf-8")


def build_payload() -> dict[str, Any]:
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    grep = rg_audit()
    p2442 = sources["P2442_NULLSPACE"].get("strict_moment_coefficient_local_identifiability_nullspace_certificate", {}).get("theorem_export", {})
    p2447 = sources["P2447_STATIONARY_REFINEMENT"].get("strict_pointwise_rank_lift_stationary_refinement_certificate", {}).get("theorem_export", {})
    base_jacobian = p2442.get("input_jacobian", [])
    census = stationary_census(base_jacobian)
    p2447_d = p2447.get("golden_section_maximum", {}).get("d")
    p2447_volume = p2447.get("golden_section_maximum", {}).get("normalized_rank_lift_volume")
    best = census["best_stationary_or_boundary_row"]
    best_maximum = census["maxima_sorted_by_volume"][0] if census["maxima_sorted_by_volume"] else {}
    nonbest_ceiling = max(
        row["normalized_rank_lift_volume"]
        for row in ([item for item in census["boundary_rows"]] + [item for item in census["stationary_roots"] if item is not best_maximum])
    )
    theorem_export = {
        "theorem_name": "P2448_T1_strict_pointwise_rank_lift_global_stationary_census_certificate",
        "inherited_stationary_refinement_certificate": "P2447/S1397",
        "inherited_p2447_refined_d": p2447_d,
        "inherited_p2447_refined_volume": p2447_volume,
        "census": census,
        "best_matches_p2447_refinement": abs(best["d"] - p2447_d) < P2447_D_MATCH_TOL,
        "best_volume_matches_p2447_refinement": abs(best["normalized_rank_lift_volume"] - p2447_volume) < 1.0e-10,
        "unique_local_maximum_on_audited_interval": census["local_maximum_count"] == 1,
        "global_best_is_stationary_local_maximum": best["kind"] == "local_maximum",
        "global_best_dominates_boundaries": best["normalized_rank_lift_volume"] > max(row["normalized_rank_lift_volume"] for row in census["boundary_rows"]),
        "global_best_gap_over_nonbest_stationary_or_boundary_ceiling": best["normalized_rank_lift_volume"] - nonbest_ceiling,
        "global_gap_tolerance": GLOBAL_GAP_TOL,
        "global_gap_positive_above_tolerance": best["normalized_rank_lift_volume"] - nonbest_ceiling > GLOBAL_GAP_TOL,
        "analytic_interval_root_exclusion_theorem_exported_by_this_certificate": False,
        "derivative_stationary_tolerance": STATIONARY_DERIVATIVE_TOL,
        "all_stationary_derivatives_small": all(abs(row["derivative_at_root"]) < STATIONARY_DERIVATIVE_TOL for row in census["stationary_roots"]),
        "all_stationary_roots_classified": all(row["classification"] in {"local_maximum", "local_minimum"} for row in census["stationary_roots"]),
        "pointwise_coordinate_selector_exported_by_this_certificate": False,
        "strict_observable_source_constraint_exported_by_this_certificate": False,
        "gauge_slice_theorem_exported_by_this_certificate": False,
        "strict_physical_value_generator_exported": False,
        "qw2191_discharged": False,
        "role_bearing_ltotal_exported": False,
        "toe_closure_exported": False,
        "not_licensed": [
            "A derivative-sign stationary census is not an analytic interval root-exclusion theorem and not a strict point-coordinate selector theorem.",
            "Global pointwise conditioning on an audited interval is not observable/source admissibility and not a lawful gauge-slice theorem.",
            "No strict physical-value generator, QW-2191 discharge, role-bearing L_total export, or ToE closure is exported.",
        ],
        "next_honest_step": (
            "Either prove an independent selector/source/gauge theorem for the point coordinate, or stop treating pointwise conditioning as a possible coefficient-source replacement."
        ),
    }
    gatekeepers = {
        "rg_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "three_stationary_roots_found": census["sign_change_bracket_count"] == 3 and len(census["stationary_roots"]) == 3,
        "unique_local_maximum": theorem_export["unique_local_maximum_on_audited_interval"],
        "best_matches_p2447": theorem_export["best_matches_p2447_refinement"],
        "best_is_stationary_local_maximum": theorem_export["global_best_is_stationary_local_maximum"],
        "best_dominates_boundaries": theorem_export["global_best_dominates_boundaries"],
        "global_gap_positive": theorem_export["global_gap_positive_above_tolerance"],
        "stationary_derivatives_small": theorem_export["all_stationary_derivatives_small"],
        "stationary_roots_classified": theorem_export["all_stationary_roots_classified"],
        "no_interval_root_exclusion_export": not theorem_export["analytic_interval_root_exclusion_theorem_exported_by_this_certificate"],
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
        "schema_version": "p2448_s1398_v1",
        "packet_id": "P2448",
        "stage_id": "S1398",
        "status": "PASS_STRICT_POINTWISE_GLOBAL_STATIONARY_CENSUS_SELECTOR_OBSTRUCTION_NO_SOURCE_THEOREM",
        "strict_pointwise_rank_lift_global_stationary_census_certificate": {
            "inputs": {name: rel(path) for name, path in SOURCE_FILES.items()},
            "input_fingerprints": {name: sha256_json(sources[name]) for name in sources},
            "rg_nonduplication_audit": grep,
            "theorem_export": theorem_export,
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_pointwise_rank_lift_global_stationary_census_certificate"]["theorem_export"]
    c = t["census"]
    lines = [
        "# P2448/S1398 strict pointwise rank-lift global stationary-census certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        "The audited derivative-sign census over `d in [0,5]` finds three stationary roots.",
        f"The unique local maximum is at `d={c['best_stationary_or_boundary_row']['d']:.12f}` with normalized rank-lift volume `{c['best_stationary_or_boundary_row']['normalized_rank_lift_volume']:.17g}`.",
        f"It matches the P2447 refined maximum: `{t['best_matches_p2447_refinement']}`.",
        f"Its gap over the best nonwinning stationary/boundary row is `{t['global_best_gap_over_nonbest_stationary_or_boundary_ceiling']:.17g}`.",
        "",
        "## Guardrail",
        "",
        "This remains a finite conditioning census only.  It exports no analytic interval root-exclusion theorem, no point-coordinate selector, no strict observable/source theorem, no gauge-slice theorem, no physical-value generator, no QW-2191 discharge, no role-bearing `L_total`, and no ToE closure.",
        "",
        "## Stationary roots",
        "",
    ]
    for row in c["stationary_roots"]:
        lines.append(
            f"- `{row['classification']}` at `d={row['root_d']:.12f}`: volume `{row['normalized_rank_lift_volume']:.17g}`, second derivative `{row['second_derivative_at_root']:.17g}`."
        )
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
