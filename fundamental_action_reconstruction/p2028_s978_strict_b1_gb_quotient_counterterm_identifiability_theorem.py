#!/usr/bin/env python3
"""P2028 S978: strict B1 GB-quotient counterterm identifiability theorem.

This consumes P2027 and turns the rank-3/null observation into an explicit
quotient theorem.  The theorem is deliberately local: it proves the B1 scalar
counterterm class modulo the Gauss-Bonnet null direction.  It does not identify
an independent a_GB and does not replace a tensor-resolved projection.
"""

from __future__ import annotations

import hashlib
import json
import platform
from pathlib import Path
from typing import Any

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2028_s978_strict_b1_gb_quotient_counterterm_identifiability_theorem.json"
MD = GEN / "p2028_s978_strict_b1_gb_quotient_counterterm_identifiability_theorem.md"

SCHEMA_VERSION = "p2028_s978_v1"
TIMESTAMP_UTC = "2026-05-25T00:00:00+00:00"
CHANNELS = ("R2", "Ric2", "Riem2", "GB")
QUOTIENT_CHANNELS = ("R2_bar", "Ric2_bar", "Riem2_bar")
NULL_VECTOR = np.array([1.0, -4.0, 1.0, -1.0], dtype=float)
TAU_SAMPLES = (-3.0, -1.0, -1.0 / 19.0, 0.0, 1.0, 3.0)


def load(name: str) -> dict[str, Any]:
    path = GEN / name
    if not path.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def file_sha256(path: Path) -> str:
    if not path.exists():
        return "MISSING"
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows_to_matrix(rows: list[dict[str, Any]]) -> np.ndarray:
    index = {name: i for i, name in enumerate(CHANNELS)}
    mat = np.zeros((len(CHANNELS), len(CHANNELS)), dtype=float)
    for row in rows:
        mat[index[row["row_channel"]], index[row["col_channel"]]] = float(row["value"])
    return mat


def rows_to_rhs(rows: list[dict[str, Any]]) -> np.ndarray:
    index = {name: i for i, name in enumerate(CHANNELS)}
    rhs = np.zeros(len(CHANNELS), dtype=float)
    for row in rows:
        rhs[index[row["channel"]]] = float(row["value"])
    return rhs


def sympy_rows(mat: sp.Matrix) -> list[list[str]]:
    return [[sp.sstr(sp.simplify(mat[i, j])) for j in range(mat.cols)] for i in range(mat.rows)]


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2027 = load("p2027_s977_strict_b1_rank3_gb_null_adaptive_quadrature_witness.json")

    basis = p2027.get("basis_decomposition") or {}
    quotient = p2027.get("quotient_identifiability") or {}
    quadrature = p2027.get("adaptive_quadrature") or {}
    full_family = p2027.get("full_four_channel_min_norm_family") or {}

    matrix_rows = quadrature.get("assembled_full_matrix_entries") or []
    rhs_rows = quadrature.get("assembled_full_rhs_entries") or []
    gram_full = rows_to_matrix(matrix_rows)
    rhs_full = rows_to_rhs(rhs_rows)

    # T maps a four-channel coefficient vector to the unique rank-3 scalar B1
    # quotient coefficients after using GB = Riem2 - 4*Ric2 + R2.
    t_np = np.array(
        [
            [1.0, 0.0, 0.0, 1.0],
            [0.0, 1.0, 0.0, -4.0],
            [0.0, 0.0, 1.0, 1.0],
        ],
        dtype=float,
    )
    section_np = np.array(
        [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, 0.0],
        ],
        dtype=float,
    )
    projector_np = section_np @ t_np
    target_q = np.array(
        [
            float((quotient.get("rank3_quotient_coefficients") or {}).get("a_R2", np.nan)),
            float((quotient.get("rank3_quotient_coefficients") or {}).get("a_Ric2", np.nan)),
            float((quotient.get("rank3_quotient_coefficients") or {}).get("a_Riem2", np.nan)),
        ],
        dtype=float,
    )
    canonical_rep = section_np @ target_q

    tau_rows = []
    max_residual_linf = 0.0
    max_q_gap_linf = 0.0
    for tau in TAU_SAMPLES:
        coeff = canonical_rep + tau * NULL_VECTOR
        residual = gram_full @ coeff - rhs_full
        q = t_np @ coeff
        residual_linf = float(np.linalg.norm(residual, ord=np.inf))
        q_gap_linf = float(np.linalg.norm(q - target_q, ord=np.inf))
        max_residual_linf = max(max_residual_linf, residual_linf)
        max_q_gap_linf = max(max_q_gap_linf, q_gap_linf)
        tau_rows.append(
            {
                "tau": float(tau),
                "coefficients_R2_Ric2_Riem2_GB": [float(v) for v in coeff],
                "quotient_coefficients": [float(v) for v in q],
                "quotient_gap_linf": q_gap_linf,
                "full_system_residual_linf": residual_linf,
            }
        )

    min_norm_coeff = np.array(
        [
            float((full_family.get("minimum_norm_coefficients_from_full_rank_deficient_lstsq") or {}).get("a_R2", np.nan)),
            float((full_family.get("minimum_norm_coefficients_from_full_rank_deficient_lstsq") or {}).get("a_Ric2", np.nan)),
            float((full_family.get("minimum_norm_coefficients_from_full_rank_deficient_lstsq") or {}).get("a_Riem2", np.nan)),
            float((full_family.get("minimum_norm_coefficients_from_full_rank_deficient_lstsq") or {}).get("a_GB", np.nan)),
        ],
        dtype=float,
    )
    min_norm_q = t_np @ min_norm_coeff
    min_norm_q_gap = float(np.linalg.norm(min_norm_q - target_q, ord=np.inf))

    t = sp.Matrix([[1, 0, 0, 1], [0, 1, 0, -4], [0, 0, 1, 1]])
    section = sp.Matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1], [0, 0, 0]])
    projector = section * t
    null_vec = sp.Matrix([1, -4, 1, -1])
    exact_checks = {
        "T_matrix_rows": sympy_rows(t),
        "section_matrix_rows": sympy_rows(section),
        "projector_matrix_rows": sympy_rows(projector),
        "T_null_vector": [sp.sstr(x) for x in t * null_vec],
        "T_section_equals_identity": bool(t * section == sp.eye(3)),
        "projector_idempotent": bool(projector * projector == projector),
        "projector_null_vector_zero": bool(projector * null_vec == sp.zeros(4, 1)),
        "projector_kernel_direction": [sp.sstr(x) for x in null_vec],
        "quotient_rule": "T(a_R2,a_Ric2,a_Riem2,a_GB)=(a_R2+a_GB, a_Ric2-4*a_GB, a_Riem2+a_GB)",
        "section_rule": "s(q_R2,q_Ric2,q_Riem2)=(q_R2,q_Ric2,q_Riem2,0)",
    }

    numerical_tolerance = max(
        1.0e-10,
        float(quadrature.get("adaptive_abs_gap_primary_vs_check", 0.0)) * 100.0,
        float(quotient.get("full_row_residual_linf_for_gb_zero_representative", 0.0)) * 100.0,
    )
    quotient_theorem_pass = (
        p2027.get("local_verdict") == "PASS_RANK3_QUOTIENT_IDENTIFIABLE_ON_SCALAR_B1_WITH_GB_NULL_TRACE"
        and basis.get("full_gram_rank") == 3
        and basis.get("full_gram_nullity") == 1
        and exact_checks["T_section_equals_identity"]
        and exact_checks["projector_idempotent"]
        and exact_checks["projector_null_vector_zero"]
        and max_residual_linf <= numerical_tolerance
        and max_q_gap_linf <= numerical_tolerance
        and min_norm_q_gap <= numerical_tolerance
    )

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2028",
        "stage_id": "S978",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "result_kind": (
            "PASS_B1_SCALAR_GB_QUOTIENT_COUNTERTERM_CLASS_THEOREM__TENSOR_EXTENSION_OPEN"
            if quotient_theorem_pass
            else "OPEN_B1_SCALAR_GB_QUOTIENT_THEOREM_OBSTRUCTION_WITH_TRACE"
        ),
        "local_verdict": "PASS_QUOTIENT_THEOREM_WITH_TRACE" if quotient_theorem_pass else "OPEN_QUOTIENT_THEOREM_WITH_TRACE",
        "route": "strict_only",
        "strict_lane_assumptions": [
            "strict_kernel_only",
            "no_legacy_transfer",
            "B1_scalar_projection_only",
            "GB_direction_quotiented_not_independently_identified",
        ],
        "depends_on": {
            "p2027_present": p2027.get("_missing") is None,
            "p2027_local_verdict": p2027.get("local_verdict"),
            "p2027_full_gram_rank": basis.get("full_gram_rank"),
            "p2027_full_gram_nullity": basis.get("full_gram_nullity"),
        },
        "input_hashes": {
            "p2027_json_sha256": file_sha256(GEN / "p2027_s977_strict_b1_rank3_gb_null_adaptive_quadrature_witness.json"),
        },
        "quotient_space": {
            "full_channel_order": list(CHANNELS),
            "quotient_channel_order": list(QUOTIENT_CHANNELS),
            "null_direction_R2_Ric2_Riem2_GB": [float(v) for v in NULL_VECTOR],
            "target_quotient_coefficients": {
                "R2_bar": float(target_q[0]),
                "Ric2_bar": float(target_q[1]),
                "Riem2_bar": float(target_q[2]),
            },
            "canonical_section_representative_R2_Ric2_Riem2_GB": [float(v) for v in canonical_rep],
            "minimum_norm_representative_quotient_coefficients": [float(v) for v in min_norm_q],
            "minimum_norm_quotient_gap_linf": min_norm_q_gap,
        },
        "exact_linear_algebra": exact_checks,
        "tau_family_invariance": {
            "coefficient_family_rule": "a(tau)=section(target_q)+tau*(1,-4,1,-1)",
            "tau_samples": tau_rows,
            "max_quotient_gap_linf": max_q_gap_linf,
            "max_full_system_residual_linf": max_residual_linf,
            "numerical_tolerance": numerical_tolerance,
            "tau_family_invariance_pass": max_q_gap_linf <= numerical_tolerance and max_residual_linf <= numerical_tolerance,
        },
        "gatekeeper_checks": {
            "p2027_rank3_null1": basis.get("full_gram_rank") == 3 and basis.get("full_gram_nullity") == 1,
            "exact_quotient_map_valid": exact_checks["T_section_equals_identity"],
            "exact_projector_valid": exact_checks["projector_idempotent"] and exact_checks["projector_null_vector_zero"],
            "tau_family_invariance_pass": max_q_gap_linf <= numerical_tolerance and max_residual_linf <= numerical_tolerance,
            "minimum_norm_rep_maps_to_same_quotient": min_norm_q_gap <= numerical_tolerance,
            "independent_a_GB_identified": False,
            "full_tensor_projection_claimed": False,
            "toe_closure_claimed": False,
        },
        "theorem_scope": {
            "licensed_statement": "On the strict B1 scalar projection, counterterm coefficients are identifiable as the quotient class T(a)=(a_R2+a_GB, a_Ric2-4*a_GB, a_Riem2+a_GB) modulo span(1,-4,1,-1).",
            "not_licensed": [
                "independent scalar B1 identification of a_GB",
                "four-channel coefficient uniqueness",
                "tensor-resolved curvature-operator projection",
                "background-family globalization beyond B1",
                "BRST/Cutkosky unitarity closure",
                "QW-2191 selector closure",
                "ToE closure",
            ],
        },
        "false_pass_guard": "This is a quotient-class theorem only.  The apparent four-coefficient minimum-norm vector is a gauge/section choice inside the GB-null family, not an independent a_GB measurement.",
        "next_honest_step": "Attempt a tensor-resolved B1 projection that separates the GB/topological direction at the operator level, or explicitly carry this quotient object into the Task-1 ledger without upgrading it to four-channel uniqueness.",
        "lay_explanation": "W skalarnej projekcji B1 da sie jednoznacznie okreslic tylko klase trzech kombinacji wspolczynnikow. Przesuniecie wzdluz kierunku Gaussa-Bonneta zmienia cztery liczby, ale nie zmienia tego, co widzi B1.",
        "environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "sympy": sp.__version__,
        },
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    md = f"""# P2028 S978 Strict B1 GB-Quotient Counterterm Identifiability Theorem

Status: `{payload['status']}`

Local verdict: `{payload['local_verdict']}`

P2027 showed that the strict B1 scalar operator profiles have rank `3` with
one Gauss-Bonnet null direction:

`n_GB = (1,-4,1,-1)`.

P2028 exports the quotient map

`T(a_R2,a_Ric2,a_Riem2,a_GB) = (a_R2+a_GB, a_Ric2-4*a_GB, a_Riem2+a_GB)`.

The canonical section is

`s(q_R2,q_Ric2,q_Riem2) = (q_R2,q_Ric2,q_Riem2,0)`.

Exact checks:
- `T(n_GB)=0`
- `T*s=I_3`
- `s*T` is idempotent
- `s*T(n_GB)=0`

Target quotient coefficients:
- R2_bar={float(target_q[0]):.16e}
- Ric2_bar={float(target_q[1]):.16e}
- Riem2_bar={float(target_q[2]):.16e}

The whole family

`a(tau)=s(T(a)) + tau*(1,-4,1,-1)`

has the same B1 scalar quotient class.  Across sampled tau values the maximum
full-system residual L-infinity norm is `{max_residual_linf:.3e}`, with
tolerance `{numerical_tolerance:.3e}`.

## Honest Interpretation

This is a local quotient theorem.  It licenses the rank-3 B1 scalar
counterterm class modulo the GB null direction.  It does not license an
independent `a_GB`, four-channel uniqueness, a tensor-resolved projection,
background-global renormalization, unitarity, selector closure, or ToE closure.
"""
    MD.write_text(md, encoding="utf-8")
    print(OUT)
    print(MD)


if __name__ == "__main__":
    main()
