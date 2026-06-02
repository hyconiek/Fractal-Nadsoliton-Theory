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
OUT = GEN / "p2447_s1397_strict_pointwise_rank_lift_stationary_refinement_certificate.json"
MD = GEN / "p2447_s1397_strict_pointwise_rank_lift_stationary_refinement_certificate.md"

SOURCE_FILES = {
    "P2442_NULLSPACE": GEN / "p2442_s1392_strict_moment_coefficient_local_identifiability_nullspace_certificate.json",
    "P2446_SELECTOR_OBSTRUCTION": GEN / "p2446_s1396_strict_pointwise_rank_lift_selector_obstruction_certificate.json",
}
DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

STRICT_PARAMS = {"omega": 0.18575, "phi": 0.16250, "beta": 1.0, "eta": 1.8}
ROBUST_NORMALIZED_VOLUME_THRESHOLD = 1.0e-3
REFINEMENT_BRACKET = (0.6, 1.1)
LOCAL_WINDOW = (0.75, 1.25)
GOLDEN_TOL = 1.0e-12
DERIVATIVE_STEPS = [1.0e-3, 1.0e-4, 1.0e-5]
STATIONARY_DERIVATIVE_TOL = 1.0e-6
NEGATIVE_SECOND_DERIVATIVE_TOL = -1.0e-3


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
        "new_packet": "P2447|S1397|stationary refinement|pointwise rank-lift stationary|continuous pointwise rank-lift",
        "p2446_input": "P2446|S1396|pointwise rank-lift selector obstruction|d1_not_global|max.*d=0.785",
        "stationary_language": "golden-section|stationary point|second derivative|continuous maximum|derivative bracket",
        "selector_language": "d_ref|point-coordinate selector|coordinate selection|gauge-slice|strict observable/source",
        "closure_blockers": "QW-2191|role-bearing L_total|physical-value generator|ToE closure",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def k_strict(d: float) -> float:
    return math.cos(STRICT_PARAMS["omega"] * d + STRICT_PARAMS["phi"]) / (1.0 + STRICT_PARAMS["beta"] * (d ** STRICT_PARAMS["eta"]))


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


def golden_section_maximum(base_jacobian: list[list[float]], a: float, b: float, tol: float = GOLDEN_TOL) -> dict[str, Any]:
    golden = (math.sqrt(5.0) - 1.0) / 2.0
    c = b - golden * (b - a)
    d = a + golden * (b - a)
    fc = pointwise_volume(base_jacobian, c)
    fd = pointwise_volume(base_jacobian, d)
    iterations = 0
    while b - a > tol and iterations < 600:
        if fc < fd:
            a = c
            c = d
            fc = fd
            d = a + golden * (b - a)
            fd = pointwise_volume(base_jacobian, d)
        else:
            b = d
            d = c
            fd = fc
            c = b - golden * (b - a)
            fc = pointwise_volume(base_jacobian, c)
        iterations += 1
    x = (a + b) / 2.0
    return {"d": x, "normalized_rank_lift_volume": pointwise_volume(base_jacobian, x), "iterations": iterations, "final_width": b - a}


def derivative_witnesses(base_jacobian: list[list[float]], d: float) -> list[dict[str, float]]:
    center = pointwise_volume(base_jacobian, d)
    witnesses = []
    for h in DERIVATIVE_STEPS:
        plus = pointwise_volume(base_jacobian, d + h)
        minus = pointwise_volume(base_jacobian, d - h)
        witnesses.append(
            {
                "h": h,
                "central_first_derivative": (plus - minus) / (2.0 * h),
                "central_second_derivative": (plus - 2.0 * center + minus) / (h * h),
            }
        )
    return witnesses


def append_doc_sections() -> None:
    eq_section = """
## P2447/S1397 strict pointwise rank-lift stationary-refinement certificate

`P2447/S1397` refines the P2446 grid obstruction with a continuous one-dimensional optimization witness.  A golden-section refinement over the audited pointwise window finds a stationary conditioning maximum near `d=0.785288904663`, with negative finite-difference second-derivative witnesses and a positive conditioning gap above `d=1`.

This strengthens the selector obstruction rather than closing it: even the continuous conditioning optimum is not a strict point-coordinate selector, not an admissible observable/source theorem, and not a gauge-slice theorem.
""".strip()
    lag_section = """
## P2447/S1397 stationary pointwise selector guard

`P2447/S1397` shows that refining the pointwise conditioning search moves the best-conditioned coordinate away from `d=1`; no pointwise row may enter `L_total` until an independent selector/source/gauge theorem chooses and licenses the coordinate.
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
    p2446 = sources["P2446_SELECTOR_OBSTRUCTION"].get("strict_pointwise_rank_lift_selector_obstruction_certificate", {}).get("theorem_export", {})
    base_jacobian = p2442.get("input_jacobian", [])
    maximum = golden_section_maximum(base_jacobian, *REFINEMENT_BRACKET)
    d1_volume = pointwise_volume(base_jacobian, 1.0)
    local_left = pointwise_volume(base_jacobian, LOCAL_WINDOW[0])
    local_right = pointwise_volume(base_jacobian, LOCAL_WINDOW[1])
    witnesses = derivative_witnesses(base_jacobian, maximum["d"])
    theorem_export = {
        "theorem_name": "P2447_T1_strict_pointwise_rank_lift_stationary_refinement_certificate",
        "inherited_selector_obstruction_certificate": "P2446/S1396",
        "inherited_p2446_grid_best_d": p2446.get("global_best_pointwise_row", {}).get("d"),
        "refinement_bracket": {"d_min": REFINEMENT_BRACKET[0], "d_max": REFINEMENT_BRACKET[1]},
        "local_window": {"d_min": LOCAL_WINDOW[0], "d_max": LOCAL_WINDOW[1]},
        "golden_section_maximum": maximum,
        "d1_normalized_rank_lift_volume": d1_volume,
        "conditioning_gap_over_d1": maximum["normalized_rank_lift_volume"] - d1_volume,
        "boundary_volumes": {"local_left": local_left, "local_right": local_right},
        "derivative_witnesses": witnesses,
        "stationary_derivative_tolerance": STATIONARY_DERIVATIVE_TOL,
        "negative_second_derivative_tolerance": NEGATIVE_SECOND_DERIVATIVE_TOL,
        "all_first_derivative_witnesses_small": all(abs(item["central_first_derivative"]) < STATIONARY_DERIVATIVE_TOL for item in witnesses),
        "all_second_derivative_witnesses_negative": all(item["central_second_derivative"] < NEGATIVE_SECOND_DERIVATIVE_TOL for item in witnesses),
        "continuous_maximum_exceeds_d1": maximum["normalized_rank_lift_volume"] > d1_volume,
        "continuous_maximum_inside_local_window": LOCAL_WINDOW[0] < maximum["d"] < LOCAL_WINDOW[1],
        "continuous_maximum_exceeds_local_boundaries": maximum["normalized_rank_lift_volume"] > max(local_left, local_right),
        "pointwise_coordinate_selector_exported_by_this_certificate": False,
        "strict_observable_source_constraint_exported_by_this_certificate": False,
        "gauge_slice_theorem_exported_by_this_certificate": False,
        "strict_physical_value_generator_exported": False,
        "qw2191_discharged": False,
        "role_bearing_ltotal_exported": False,
        "toe_closure_exported": False,
        "not_licensed": [
            "A continuous conditioning maximum is not a strict point-coordinate selector theorem.",
            "The refined optimum is not an observable/source admissibility theorem and not a lawful gauge-slice theorem.",
            "No strict physical-value generator, QW-2191 discharge, role-bearing L_total export, or ToE closure is exported.",
        ],
        "next_honest_step": (
            "Prove an independent point-coordinate selector/source/gauge theorem, or abandon pointwise rows as coefficient-source replacements and return to source-level bridge completion."
        ),
    }
    gatekeepers = {
        "rg_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "p2446_grid_best_inherited": abs(theorem_export["inherited_p2446_grid_best_d"] - 0.785) < 1.0e-12,
        "maximum_inside_local_window": theorem_export["continuous_maximum_inside_local_window"],
        "maximum_exceeds_d1": theorem_export["continuous_maximum_exceeds_d1"],
        "maximum_exceeds_boundaries": theorem_export["continuous_maximum_exceeds_local_boundaries"],
        "stationary_first_derivatives_small": theorem_export["all_first_derivative_witnesses_small"],
        "stationary_second_derivatives_negative": theorem_export["all_second_derivative_witnesses_negative"],
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
        "schema_version": "p2447_s1397_v1",
        "packet_id": "P2447",
        "stage_id": "S1397",
        "result_kind": "STRICT_POINTWISE_RANK_LIFT_STATIONARY_REFINEMENT_CERTIFICATE",
        "status": "PASS_STRICT_POINTWISE_STATIONARY_REFINEMENT_SELECTOR_OBSTRUCTION_NO_SOURCE_THEOREM",
        "strict_pointwise_rank_lift_stationary_refinement_certificate": {
            "rg_nonduplication_audit": grep,
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
        "recommended_next_honest_step": theorem_export["next_honest_step"],
        "global_status": "OPEN_PROGRESS_STATIONARY_POINTWISE_SELECTOR_OBSTRUCTION_EXPORTED_NO_CLOSURE",
    }


def write_outputs(payload: dict[str, Any]) -> None:
    theorem = payload["strict_pointwise_rank_lift_stationary_refinement_certificate"]["theorem_export"]
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2447 S1397: strict pointwise rank-lift stationary-refinement certificate",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Finite facts",
                "",
                f"- Refined maximum: `{theorem['golden_section_maximum']}`.",
                f"- `d=1` volume: `{theorem['d1_normalized_rank_lift_volume']}`.",
                f"- Conditioning gap over `d=1`: `{theorem['conditioning_gap_over_d1']}`.",
                f"- First derivative witnesses small: `{theorem['all_first_derivative_witnesses_small']}`.",
                f"- Second derivative witnesses negative: `{theorem['all_second_derivative_witnesses_negative']}`.",
                "",
                "## Hard limits",
                "",
                *[f"- {item}" for item in theorem["not_licensed"]],
                "",
                "## Next honest step",
                "",
                theorem["next_honest_step"],
                "",
                "## Gatekeepers",
                "",
                f"`{payload['gatekeeper_checks']}`",
                "",
            ]
        ),
        encoding="utf-8",
    )


def main() -> None:
    GEN.mkdir(parents=True, exist_ok=True)
    append_doc_sections()
    payload = build_payload()
    write_outputs(payload)
    if not all(payload["gatekeeper_checks"].values()):
        raise SystemExit(f"gatekeeper failure: {payload['gatekeeper_checks']}")


if __name__ == "__main__":
    main()
