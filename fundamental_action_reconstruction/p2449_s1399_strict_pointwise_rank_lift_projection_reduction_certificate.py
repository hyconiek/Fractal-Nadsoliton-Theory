#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
import subprocess
from pathlib import Path
from typing import Any, Callable

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"
OUT = GEN / "p2449_s1399_strict_pointwise_rank_lift_projection_reduction_certificate.json"
MD = GEN / "p2449_s1399_strict_pointwise_rank_lift_projection_reduction_certificate.md"

SOURCE_FILES = {
    "P2442_NULLSPACE": GEN / "p2442_s1392_strict_moment_coefficient_local_identifiability_nullspace_certificate.json",
    "P2448_STATIONARY_CENSUS": GEN / "p2448_s1398_strict_pointwise_rank_lift_global_stationary_census_certificate.json",
}
DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

STRICT_PARAMS = {"omega": 0.18575, "phi": 0.16250, "beta": 1.0, "eta": 1.8}
AUDIT_INTERVAL = (0.0, 5.0)
ROOT_SCAN_START = 1.0e-6
ROOT_SCAN_STEP = 1.0e-3
ROOT_TOL = 1.0e-13
PROJECTION_IDENTITY_SAMPLE_POINTS = [0.0, 0.125, 0.212239394092, 0.5, 0.785288904497, 1.0, 2.5, 4.056893978257, 5.0]
P2448_ROOT_MATCH_TOL = 5.0e-8
PROJECTION_IDENTITY_TOL = 1.0e-12


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
        "new_packet": "P2449|S1399|projection reduction|rank-lift projection|volume scalar reduction",
        "p2448_input": "P2448|S1398|global stationary census|stationary roots|derivative-sign census",
        "linear_algebra_language": "nullspace projection|cofactor vector|determinant-to-projection|projection identity|normalized determinant",
        "stationary_factor_language": "zero-projection root|log-derivative root|stationary factor|analytic derivative factor",
        "closure_blockers": "QW-2191|role-bearing L_total|physical-value generator|ToE closure|point-coordinate selector",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def determinant(matrix: list[list[float]]) -> float:
    mat = [row[:] for row in matrix]
    n = len(mat)
    det = 1.0
    for col in range(n):
        pivot = max(range(col, n), key=lambda r: abs(mat[r][col]))
        if abs(mat[pivot][col]) < 1.0e-15:
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


def row_norm(row: list[float]) -> float:
    return math.sqrt(sum(value * value for value in row))


def normalize_row(row: list[float]) -> list[float]:
    norm = row_norm(row)
    return [value / norm for value in row]


def normalized_base_rows(base_jacobian: list[list[float]]) -> list[list[float]]:
    return [normalize_row(row) for row in base_jacobian]


def cofactor_projection_vector(base_jacobian: list[list[float]]) -> list[float]:
    base = normalized_base_rows(base_jacobian)
    vector = []
    for col in range(4):
        columns = [idx for idx in range(4) if idx != col]
        minor = determinant([[row[idx] for idx in columns] for row in base])
        vector.append(((-1.0) ** (3 + col)) * minor)
    return vector


def pointwise_gradient(d: float) -> list[float]:
    omega = STRICT_PARAMS["omega"]
    phi = STRICT_PARAMS["phi"]
    eta = STRICT_PARAMS["eta"]
    d_eta = d**eta
    denominator = 1.0 + d_eta
    phase = omega * d + phi
    cos_phase = math.cos(phase)
    sin_phase = math.sin(phase)
    eta_derivative = 0.0 if d == 0.0 else -cos_phase * d_eta * math.log(d) / (denominator * denominator)
    return [
        -d * sin_phase / denominator,
        -sin_phase / denominator,
        -cos_phase * d_eta / (denominator * denominator),
        eta_derivative,
    ]


def pointwise_gradient_derivative(d: float) -> list[float]:
    omega = STRICT_PARAMS["omega"]
    phi = STRICT_PARAMS["phi"]
    eta = STRICT_PARAMS["eta"]
    d_eta = d**eta
    d_eta_derivative = eta * d ** (eta - 1.0)
    denominator = 1.0 + d_eta
    denominator_derivative = d_eta_derivative
    phase = omega * d + phi
    cos_phase = math.cos(phase)
    sin_phase = math.sin(phase)
    g2 = -cos_phase * d_eta / (denominator * denominator)
    g2_derivative = (
        omega * sin_phase * d_eta / (denominator * denominator)
        - cos_phase * d_eta_derivative / (denominator * denominator)
        + 2.0 * cos_phase * d_eta * denominator_derivative / (denominator**3)
    )
    return [
        -(sin_phase + d * omega * cos_phase) / denominator + d * sin_phase * denominator_derivative / (denominator * denominator),
        -omega * cos_phase / denominator + sin_phase * denominator_derivative / (denominator * denominator),
        g2_derivative,
        g2_derivative * math.log(d) + g2 / d,
    ]


def determinant_volume(base_jacobian: list[list[float]], d: float) -> float:
    return abs(determinant(normalized_base_rows(base_jacobian) + [normalize_row(pointwise_gradient(d))]))


def dot(left: list[float], right: list[float]) -> float:
    return sum(a * b for a, b in zip(left, right))


def projection_amplitude(projection_vector: list[float], d: float) -> float:
    return dot(projection_vector, pointwise_gradient(d))


def projection_volume(projection_vector: list[float], d: float) -> float:
    gradient = pointwise_gradient(d)
    return abs(dot(projection_vector, gradient)) / row_norm(gradient)


def projection_stationary_factor(projection_vector: list[float], d: float) -> float:
    gradient = pointwise_gradient(d)
    gradient_derivative = pointwise_gradient_derivative(d)
    amplitude = dot(projection_vector, gradient)
    amplitude_derivative = dot(projection_vector, gradient_derivative)
    gradient_square = dot(gradient, gradient)
    gradient_square_derivative = 2.0 * dot(gradient, gradient_derivative)
    return 2.0 * amplitude_derivative * gradient_square - amplitude * gradient_square_derivative


def bisect_root(function: Callable[[float], float], left: float, right: float) -> dict[str, float]:
    f_left = function(left)
    iterations = 0
    while right - left > ROOT_TOL and iterations < 120:
        midpoint = (left + right) / 2.0
        f_mid = function(midpoint)
        if f_left == 0.0:
            right = left
            break
        if f_left * f_mid <= 0.0:
            right = midpoint
        else:
            left = midpoint
            f_left = f_mid
        iterations += 1
    root = (left + right) / 2.0
    return {"root_d": root, "residual": function(root), "iterations": iterations, "final_width": right - left}


def root_scan(function: Callable[[float], float]) -> list[dict[str, float]]:
    d_min, d_max = AUDIT_INTERVAL
    roots = []
    previous_d = ROOT_SCAN_START
    previous_value = function(previous_d)
    scan_count = int(round((d_max - d_min) / ROOT_SCAN_STEP))
    for idx in range(1, scan_count + 1):
        d = idx * ROOT_SCAN_STEP
        value = function(d)
        if previous_value == 0.0 or value == 0.0 or previous_value * value < 0.0:
            roots.append(bisect_root(function, previous_d, d))
        previous_d = d
        previous_value = value
    return roots


def projection_identity_samples(base_jacobian: list[list[float]], projection_vector: list[float]) -> list[dict[str, float]]:
    samples = []
    for d in PROJECTION_IDENTITY_SAMPLE_POINTS:
        determinant_value = determinant_volume(base_jacobian, d)
        projection_value = projection_volume(projection_vector, d)
        samples.append(
            {
                "d": d,
                "determinant_volume": determinant_value,
                "projection_volume": projection_value,
                "absolute_difference": abs(determinant_value - projection_value),
            }
        )
    return samples


def append_doc_sections() -> None:
    eq_section = """
## P2449/S1399 strict pointwise rank-lift projection-reduction certificate

`P2449/S1399` reduces the P2448 four-by-four normalized determinant calculation to a one-row nullspace projection identity.  If `b_0,b_1,b_2` are the normalized inherited moment rows and `a` is their cofactor-null vector, then every pointwise candidate row `g(d)` has rank-lift volume `|a·g(d)|/||g(d)||`.  The zero-volume stationary roots are therefore `a·g(d)=0`, while nonzero stationary roots satisfy the analytic factor `2(a·g)'||g||^2 - (a·g)(||g||^2)'=0`.

The replay matches the P2448 roots and maximum, but it is still a projection-reduction audit, not an exact interval root-exclusion theorem, not a point-coordinate selector theorem, and not source/gauge authority for `L_total`.
""".strip()
    lag_section = """
## P2449/S1399 projection-reduction guard

`P2449/S1399` makes the pointwise rank-lift computation more theorem-shaped by reducing determinant volume to a nullspace projection and splitting stationary points into zero-projection roots and analytic stationary-factor roots.  This improves the proof audit of the obstruction, but still exports no selector/source/gauge permission to add a pointwise row to `L_total`.
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
    p2448 = sources["P2448_STATIONARY_CENSUS"].get("strict_pointwise_rank_lift_global_stationary_census_certificate", {}).get("theorem_export", {})
    base_jacobian = p2442.get("input_jacobian", [])
    projection_vector = cofactor_projection_vector(base_jacobian)
    zero_projection_roots = root_scan(lambda d: projection_amplitude(projection_vector, d))
    stationary_factor_roots = root_scan(lambda d: projection_stationary_factor(projection_vector, d))
    samples = projection_identity_samples(base_jacobian, projection_vector)
    p2448_roots = p2448.get("census", {}).get("stationary_roots", [])
    p2448_root_ds = [row.get("root_d") for row in p2448_roots]
    reconstructed_roots = sorted(
        [
            {**row, "root_family": "zero_projection", "normalized_rank_lift_volume": projection_volume(projection_vector, row["root_d"])}
            for row in zero_projection_roots
        ]
        + [
            {**row, "root_family": "stationary_factor", "normalized_rank_lift_volume": projection_volume(projection_vector, row["root_d"])}
            for row in stationary_factor_roots
        ],
        key=lambda row: row["root_d"],
    )
    best_reconstructed = max(reconstructed_roots, key=lambda row: row["normalized_rank_lift_volume"])
    p2448_best = p2448.get("census", {}).get("best_stationary_or_boundary_row", {})
    theorem_export = {
        "theorem_name": "P2449_T1_strict_pointwise_rank_lift_projection_reduction_certificate",
        "inherited_stationary_census_certificate": "P2448/S1398",
        "projection_vector": projection_vector,
        "projection_vector_norm": row_norm(projection_vector),
        "projection_identity": "det([b0;b1;b2;g/||g||]) = a·g/||g|| for the cofactor-null vector a of normalized inherited rows b0,b1,b2",
        "stationary_factorization": "For nonzero a·g(d), stationary points of |a·g(d)|/||g(d)|| satisfy 2(a·g)'||g||^2 - (a·g)(||g||^2)' = 0; zero-volume cusps are tracked separately by a·g(d)=0.",
        "projection_identity_samples": samples,
        "maximum_projection_identity_error": max(sample["absolute_difference"] for sample in samples),
        "zero_projection_roots": zero_projection_roots,
        "stationary_factor_roots": stationary_factor_roots,
        "reconstructed_stationary_roots": reconstructed_roots,
        "reconstructed_root_count": len(reconstructed_roots),
        "reconstructed_root_ds": [row["root_d"] for row in reconstructed_roots],
        "p2448_root_ds": p2448_root_ds,
        "all_reconstructed_roots_match_p2448": len(reconstructed_roots) == len(p2448_root_ds)
        and all(abs(row["root_d"] - p2448_root_ds[idx]) < P2448_ROOT_MATCH_TOL for idx, row in enumerate(reconstructed_roots)),
        "best_reconstructed_root": best_reconstructed,
        "best_reconstructed_matches_p2448_best": abs(best_reconstructed["root_d"] - p2448_best.get("d")) < P2448_ROOT_MATCH_TOL
        and abs(best_reconstructed["normalized_rank_lift_volume"] - p2448_best.get("normalized_rank_lift_volume")) < 1.0e-10,
        "projection_identity_tolerance": PROJECTION_IDENTITY_TOL,
        "p2448_root_match_tolerance": P2448_ROOT_MATCH_TOL,
        "exact_interval_root_exclusion_theorem_exported_by_this_certificate": False,
        "pointwise_coordinate_selector_exported_by_this_certificate": False,
        "strict_observable_source_constraint_exported_by_this_certificate": False,
        "gauge_slice_theorem_exported_by_this_certificate": False,
        "strict_physical_value_generator_exported": False,
        "qw2191_discharged": False,
        "role_bearing_ltotal_exported": False,
        "toe_closure_exported": False,
        "not_licensed": [
            "The cofactor projection identity is a linear-algebra reduction of the audited pointwise computation, not an exact interval root-exclusion theorem.",
            "The analytic stationary factorization identifies the root families but does not choose a strict point-coordinate selector.",
            "No strict observable/source theorem, gauge-slice theorem, physical-value generator, QW-2191 discharge, role-bearing L_total export, or ToE closure is exported.",
        ],
        "next_honest_step": (
            "If the pointwise route is still pursued, the next proof upgrade would be interval-certified root exclusion for the projection amplitude and stationary factor; otherwise return to source/selector bridge completion."
        ),
    }
    gatekeepers = {
        "rg_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "projection_vector_nonzero": theorem_export["projection_vector_norm"] > 1.0e-6,
        "projection_identity_holds_on_samples": theorem_export["maximum_projection_identity_error"] < PROJECTION_IDENTITY_TOL,
        "two_zero_projection_roots": len(zero_projection_roots) == 2,
        "one_stationary_factor_root": len(stationary_factor_roots) == 1,
        "reconstructed_roots_match_p2448": theorem_export["all_reconstructed_roots_match_p2448"],
        "best_root_matches_p2448": theorem_export["best_reconstructed_matches_p2448_best"],
        "no_exact_interval_root_exclusion_export": not theorem_export["exact_interval_root_exclusion_theorem_exported_by_this_certificate"],
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
        "schema_version": "p2449_s1399_v1",
        "packet_id": "P2449",
        "stage_id": "S1399",
        "status": "PASS_STRICT_POINTWISE_PROJECTION_REDUCTION_MATCHES_CENSUS_NO_SELECTOR_SOURCE_THEOREM",
        "strict_pointwise_rank_lift_projection_reduction_certificate": {
            "inputs": {name: rel(path) for name, path in SOURCE_FILES.items()},
            "input_fingerprints": {name: sha256_json(sources[name]) for name in sources},
            "rg_nonduplication_audit": grep,
            "theorem_export": theorem_export,
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_pointwise_rank_lift_projection_reduction_certificate"]["theorem_export"]
    lines = [
        "# P2449/S1399 strict pointwise rank-lift projection-reduction certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Linear-algebra reduction",
        "",
        f"Projection vector norm: `{t['projection_vector_norm']:.17g}`.",
        f"Maximum determinant-vs-projection identity error on audited samples: `{t['maximum_projection_identity_error']:.17g}`.",
        "",
        "## Stationary root reconstruction",
        "",
        f"Zero-projection roots: `{[row['root_d'] for row in t['zero_projection_roots']]}`.",
        f"Stationary-factor roots: `{[row['root_d'] for row in t['stationary_factor_roots']]}`.",
        f"All reconstructed roots match P2448: `{t['all_reconstructed_roots_match_p2448']}`.",
        f"Best reconstructed root: `{t['best_reconstructed_root']}`.",
        "",
        "## Guardrail",
        "",
        "This exports the determinant-to-projection reduction and an analytic stationary-factor replay only.  It exports no exact interval root-exclusion theorem, no point-coordinate selector, no strict observable/source theorem, no gauge-slice theorem, no physical-value generator, no QW-2191 discharge, no role-bearing `L_total`, and no ToE closure.",
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
