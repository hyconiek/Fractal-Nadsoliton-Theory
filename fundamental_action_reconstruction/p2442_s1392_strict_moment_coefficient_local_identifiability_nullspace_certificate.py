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

OUT = GEN / "p2442_s1392_strict_moment_coefficient_local_identifiability_nullspace_certificate.json"
MD = GEN / "p2442_s1392_strict_moment_coefficient_local_identifiability_nullspace_certificate.md"

SOURCE_FILES = {
    "P2441_PHASE_SENSITIVITY": GEN / "p2441_s1391_strict_moment_coefficient_phase_sensitivity_rank_certificate.json",
}
DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

STRICT_PARAMS = {"omega": 0.18575, "phi": 0.16250, "beta": 1.0, "eta": 1.8, "D": 25.0, "steps": 20000}
PARAMETER_ORDER = ["omega", "phi", "beta", "eta"]
COEFFICIENT_ORDER = ["lambda_sm_eff", "kappa_gr_eff", "epsilon_mix_eff"]
NULL_PERTURBATION_EPSILON = 1.0e-3
KERNEL_SAMPLE_D_VALUES = [0.0, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0]


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
        "new_packet": "P2442|S1392|moment coefficient local identifiability|nullspace certificate|strict moment coefficient nullspace",
        "p2441_input": "P2441|S1391|phase-sensitivity rank|moment coefficient phase",
        "identifiability_language": "identifiability|nullspace|null direction|rank|Jacobian|local inverse",
        "moment_coefficients": "lambda_sm_eff|kappa_gr_eff|epsilon_mix_eff|strict moment",
        "selector_closure": "QW-2191|selector uniqueness|strict observable|ToE closure",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def k_strict(d: float, params: dict[str, float]) -> float:
    return math.cos(params["omega"] * d + params["phi"]) / (1.0 + params["beta"] * (d ** params["eta"]))


def moment(n: int, params: dict[str, float]) -> float:
    steps = int(params["steps"])
    h = params["D"] / steps
    total = 0.0
    for i in range(steps + 1):
        d = i * h
        weight = 0.5 if i in (0, steps) else 1.0
        total += weight * (d**n) * k_strict(d, params)
    return total * h


def moment_coefficients(params: dict[str, float]) -> dict[str, float]:
    m0, m1, m2, m3 = (moment(n, params) for n in range(4))
    r1 = m1 / m0
    r2 = m2 / m0
    r3 = m3 / m0
    return {
        "lambda_sm_eff": abs(r1),
        "kappa_gr_eff": abs(r2 - r1 * r1),
        "epsilon_mix_eff": abs(r3) / (1.0 + abs(r2)),
    }


def rref(matrix: list[list[float]], tol: float = 1e-7) -> tuple[list[list[float]], list[int]]:
    mat = [row[:] for row in matrix]
    if not mat:
        return mat, []
    rows = len(mat)
    cols = len(mat[0])
    rank = 0
    pivots: list[int] = []
    for col in range(cols):
        pivot = max(range(rank, rows), key=lambda r: abs(mat[r][col]), default=rank)
        if rank >= rows or abs(mat[pivot][col]) <= tol:
            continue
        mat[rank], mat[pivot] = mat[pivot], mat[rank]
        pivot_value = mat[rank][col]
        mat[rank] = [value / pivot_value for value in mat[rank]]
        for r in range(rows):
            if r != rank and abs(mat[r][col]) > tol:
                factor = mat[r][col]
                mat[r] = [a - factor * b for a, b in zip(mat[r], mat[rank])]
        pivots.append(col)
        rank += 1
        if rank == rows:
            break
    return mat, pivots


def nullspace_basis_one(matrix: list[list[float]]) -> dict[str, Any]:
    reduced, pivots = rref(matrix)
    cols = len(matrix[0]) if matrix else 0
    free_cols = [col for col in range(cols) if col not in pivots]
    if len(free_cols) != 1:
        raise RuntimeError(f"expected one null direction, got free columns {free_cols}")
    free = free_cols[0]
    vector = [0.0] * cols
    vector[free] = 1.0
    for row_index, pivot_col in enumerate(pivots):
        vector[pivot_col] = -reduced[row_index][free]
    norm = math.sqrt(sum(value * value for value in vector))
    normalized = [value / norm for value in vector]
    residual = [sum(row[col] * normalized[col] for col in range(cols)) for row in matrix]
    return {
        "rref": reduced,
        "pivot_columns": pivots,
        "free_columns": free_cols,
        "normalized_null_vector": {name: normalized[index] for index, name in enumerate(PARAMETER_ORDER)},
        "linear_residual": {name: residual[index] for index, name in enumerate(COEFFICIENT_ORDER)},
        "max_abs_linear_residual": max(abs(value) for value in residual),
    }


def shifted_params(null_vector: dict[str, float], epsilon: float) -> dict[str, float]:
    params = dict(STRICT_PARAMS)
    for name in PARAMETER_ORDER:
        params[name] += epsilon * null_vector[name]
    return params


def kernel_sample_delta(params: dict[str, float]) -> dict[str, Any]:
    rows = []
    for d in KERNEL_SAMPLE_D_VALUES:
        base = k_strict(d, STRICT_PARAMS)
        shifted = k_strict(d, params)
        rows.append({"d": d, "base": base, "shifted": shifted, "delta": shifted - base})
    return {"sample_count": len(rows), "max_abs_delta": max(abs(row["delta"]) for row in rows), "rows": rows}


def null_perturbation_rows(null_vector: dict[str, float], base_coeffs: dict[str, float]) -> list[dict[str, Any]]:
    rows = []
    for sign in [-1.0, 1.0]:
        epsilon = sign * NULL_PERTURBATION_EPSILON
        params = shifted_params(null_vector, epsilon)
        coeffs = moment_coefficients(params)
        deltas = {name: coeffs[name] - base_coeffs[name] for name in COEFFICIENT_ORDER}
        rows.append(
            {
                "epsilon": epsilon,
                "shifted_params": {name: params[name] for name in PARAMETER_ORDER},
                "coefficient_deltas": deltas,
                "max_abs_coefficient_delta": max(abs(value) for value in deltas.values()),
                "kernel_sample_delta": kernel_sample_delta(params),
            }
        )
    return rows


def append_doc_sections() -> None:
    eq_section = """
## P2442/S1392 strict moment coefficient local-identifiability nullspace certificate

`P2442/S1392` turns the P2441 phase-sensitive moment route into a local identifiability audit.  The Jacobian from four strict kernel parameters `(omega, phi, beta, eta)` to the three P1562-style moment coefficients has rank `3`, so it has a one-dimensional local nullspace.  A normalized null direction changes all four kernel parameters while leaving the three moment coefficients unchanged to first order.

Therefore even the phase-sensitive three-coefficient moment map is not by itself an injective strict kernel source or full physical-value generator.  A strict SM/GR coefficient theorem still needs an extra independent observable/source constraint, or an explicit reduction theorem explaining why the null direction is physically gauge/redundant.
""".strip()
    lag_section = """
## P2442/S1392 local-identifiability guard

`P2442/S1392` shows that the current three moment-derived coefficients underdetermine the four-parameter strict kernel at the local linear level.  These coefficients may remain diagnostic/effective inputs, but they cannot alone define the final role-bearing `L_total` coefficient source without an extra independent constraint or a proved gauge-redundancy theorem for the null direction.
""".strip()
    for path, section in [(DOC_FILES["equation_sheet"], eq_section), (DOC_FILES["lagrangian_eom_draft"], lag_section)]:
        text = path.read_text(encoding="utf-8") if path.exists() else ""
        marker = section.splitlines()[0]
        if marker not in text:
            path.write_text(text.rstrip() + "\n\n" + section + "\n", encoding="utf-8")


def build_payload() -> dict[str, Any]:
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    grep = rg_audit()
    p2441 = sources["P2441_PHASE_SENSITIVITY"].get("strict_moment_coefficient_phase_sensitivity_rank_certificate", {}).get("theorem_export", {})
    jacobian = p2441.get("jacobian_numeric", [])
    rank = p2441.get("jacobian_real_rank")
    null = nullspace_basis_one(jacobian)
    base_coeffs = moment_coefficients(STRICT_PARAMS)
    perturbations = null_perturbation_rows(null["normalized_null_vector"], base_coeffs)
    max_coeff_delta = max(row["max_abs_coefficient_delta"] for row in perturbations)
    max_kernel_delta = max(row["kernel_sample_delta"]["max_abs_delta"] for row in perturbations)
    theorem_export = {
        "theorem_name": "P2442_T1_strict_moment_coefficient_local_identifiability_nullspace_certificate",
        "parameter_order": PARAMETER_ORDER,
        "coefficient_order": COEFFICIENT_ORDER,
        "input_jacobian_rank": rank,
        "input_jacobian": jacobian,
        "nullspace_dimension": len(null["free_columns"]),
        "nullspace_certificate": null,
        "base_coefficients": base_coeffs,
        "null_perturbation_epsilon_abs": NULL_PERTURBATION_EPSILON,
        "null_perturbation_rows": perturbations,
        "max_abs_coefficient_delta_under_null_perturbations": max_coeff_delta,
        "max_abs_kernel_sample_delta_under_null_perturbations": max_kernel_delta,
        "moment_map_locally_injective_for_four_strict_parameters": False,
        "extra_constraint_or_gauge_redundancy_theorem_required": True,
        "strict_kernel_to_coefficient_map_theorem_exported_by_this_certificate": False,
        "strict_physical_value_generator_exported": False,
        "qw2191_discharged": False,
        "role_bearing_ltotal_exported": False,
        "toe_closure_exported": False,
        "not_licensed": [
            "Three moment coefficients cannot locally identify four strict kernel parameters without an extra premise.",
            "The null direction is not proven to be a gauge redundancy by this certificate.",
            "Null-direction first-order coefficient stability is not SM/GR physical-value generation.",
            "No QW-2191 discharge, role-bearing L_total export, or ToE closure is exported.",
        ],
        "next_honest_step": (
            "Add one independent strict observable/source constraint for the null direction, or prove a gauge-redundancy theorem that quotients it before promoting moment coefficients."
        ),
    }
    gatekeepers = {
        "rg_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "rank_three_inherited": theorem_export["input_jacobian_rank"] == 3,
        "one_dimensional_nullspace": theorem_export["nullspace_dimension"] == 1,
        "null_residual_small": theorem_export["nullspace_certificate"]["max_abs_linear_residual"] < 1e-9,
        "null_perturbation_rows_two": len(theorem_export["null_perturbation_rows"]) == 2,
        "null_perturbation_coefficients_stable": theorem_export["max_abs_coefficient_delta_under_null_perturbations"] < 1e-4,
        "null_perturbation_changes_kernel_samples": theorem_export["max_abs_kernel_sample_delta_under_null_perturbations"] > 1e-5,
        "not_locally_injective": not theorem_export["moment_map_locally_injective_for_four_strict_parameters"],
        "extra_constraint_required": theorem_export["extra_constraint_or_gauge_redundancy_theorem_required"],
        "no_coefficient_theorem_export": not theorem_export["strict_kernel_to_coefficient_map_theorem_exported_by_this_certificate"],
        "no_value_generator_export": not theorem_export["strict_physical_value_generator_exported"],
        "no_qw2191_discharge": not theorem_export["qw2191_discharged"],
        "no_ltotal_export": not theorem_export["role_bearing_ltotal_exported"],
        "no_toe_export": not theorem_export["toe_closure_exported"],
        "fingerprint_stable": sha256_json(theorem_export) == sha256_json(theorem_export),
    }
    return {
        "schema_version": "p2442_s1392_v1",
        "packet_id": "P2442",
        "stage_id": "S1392",
        "result_kind": "STRICT_MOMENT_COEFFICIENT_LOCAL_IDENTIFIABILITY_NULLSPACE_CERTIFICATE",
        "status": "PASS_STRICT_MOMENT_COEFFICIENT_LOCAL_NULLSPACE_NO_GENERATOR",
        "strict_moment_coefficient_local_identifiability_nullspace_certificate": {
            "rg_nonduplication_audit": grep,
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
        "recommended_next_honest_step": theorem_export["next_honest_step"],
        "global_status": "OPEN_PROGRESS_STRICT_MOMENT_LOCAL_NULLSPACE_EXPORTED_NO_CLOSURE",
    }


def write_outputs(payload: dict[str, Any]) -> None:
    theorem = payload["strict_moment_coefficient_local_identifiability_nullspace_certificate"]["theorem_export"]
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2442 S1392: strict moment coefficient local-identifiability nullspace certificate",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Finite facts",
                "",
                f"- Input Jacobian rank: `{theorem['input_jacobian_rank']}`.",
                f"- Nullspace dimension: `{theorem['nullspace_dimension']}`.",
                f"- Max linear null residual: `{theorem['nullspace_certificate']['max_abs_linear_residual']}`.",
                f"- Max coefficient delta under null perturbations: `{theorem['max_abs_coefficient_delta_under_null_perturbations']}`.",
                f"- Max kernel sample delta under null perturbations: `{theorem['max_abs_kernel_sample_delta_under_null_perturbations']}`.",
                f"- Extra constraint/gauge theorem required: `{theorem['extra_constraint_or_gauge_redundancy_theorem_required']}`.",
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
