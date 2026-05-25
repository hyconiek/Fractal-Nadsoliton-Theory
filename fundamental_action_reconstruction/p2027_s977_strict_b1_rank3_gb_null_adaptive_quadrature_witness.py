#!/usr/bin/env python3
"""P2027 S977: p1950-next rank-3/GB-null adaptive quadrature witness.

This is the next honest step after P1950 exposed the B1 scalar-profile
rank deficiency.  It keeps the strict-only lane, separates the independent
rank-3 profile quotient from the explicit Gauss-Bonnet null direction, and
replaces the Simpson replay gate with adaptive endpoint-transformed
quadrature for all Kdd**2-bearing profile products.

Scope discipline: this can identify the B1 scalar quotient class of the local
counterterm profile, but it cannot identify an independent four-channel GB
coefficient and does not prove full tensor-level renormalization closure.
"""

from __future__ import annotations

import hashlib
import json
import platform
from pathlib import Path
from typing import Any, Callable

import numpy as np
import scipy.integrate as si
import sympy as sp

from p1950_s900_strict_renormalization_exact_integration import (
    _backend_channel_defs,
    strict_kernel,
)

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2027_s977_strict_b1_rank3_gb_null_adaptive_quadrature_witness.json"
MD = GEN / "p2027_s977_strict_b1_rank3_gb_null_adaptive_quadrature_witness_theorem.md"

SCHEMA_VERSION = "p2027_s977_v1"
TIMESTAMP_UTC = "2026-05-25T00:00:00+00:00"
EPS_UV = 1.0e-6
QUAD_EPSABS = 1.0e-10
QUAD_EPSREL = 1.0e-10
QUAD_LIMIT = 800
PRIMARY_LEFT_POWER = 5.0
CHECK_LEFT_POWER = 7.0
ADAPTIVE_STABILITY_TOL = 1.0e-8
CHANNEL_ORDER = ("R2", "Ric2", "Riem2", "GB")
INDEPENDENT_CHANNELS = ("R2", "Ric2", "Riem2")
FULL_COEFFICIENT_NAMES = ("a_R2", "a_Ric2", "a_Riem2", "a_GB")
NULL_VECTOR_EXACT = np.array([1.0, -4.0, 1.0, -1.0], dtype=float)


def load(name: str) -> dict[str, Any]:
    path = GEN / name
    if not path.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def file_sha256(path: Path) -> str:
    if not path.exists():
        return "MISSING"
    return hashlib.sha256(path.read_bytes()).hexdigest()


def as_float(value: Any) -> float:
    arr = np.asarray(value, dtype=float)
    if arr.shape == ():
        return float(arr)
    return float(arr.reshape(-1)[0])


def scalar_fn(fn: Callable[[Any], Any]) -> Callable[[float], float]:
    return lambda x: as_float(fn(float(x)))


def left_endpoint_quad(
    fn: Callable[[float], float],
    a: float,
    b: float,
    power: float,
) -> tuple[float, float]:
    """Integrate fn on [a,b] with x=a+(b-a)t**power."""

    width = b - a

    def transformed(t: float) -> float:
        if t <= 0.0:
            return 0.0
        x = a + width * (t**power)
        jac = width * power * (t ** (power - 1.0))
        return float(fn(x) * jac)

    val, err = si.quad(
        transformed,
        0.0,
        1.0,
        epsabs=QUAD_EPSABS,
        epsrel=QUAD_EPSREL,
        limit=QUAD_LIMIT,
    )
    return float(val), float(err)


def native_quad(fn: Callable[[float], float], a: float, b: float) -> tuple[float, float]:
    val, err = si.quad(
        fn,
        a,
        b,
        epsabs=QUAD_EPSABS,
        epsrel=QUAD_EPSREL,
        limit=QUAD_LIMIT,
    )
    return float(val), float(err)


def uses_kdd_profile(channel: str) -> bool:
    return channel in {"Riem2", "GB"}


def integrate_profile_product(
    fn: Callable[[float], float],
    left_channel: str,
    right_channel: str,
    power: float,
) -> tuple[float, float, str]:
    if uses_kdd_profile(left_channel) or uses_kdd_profile(right_channel):
        val, err = left_endpoint_quad(fn, EPS_UV, 1.0, power)
        return val, err, f"left_endpoint_power_{power:g}"
    val, err = native_quad(fn, EPS_UV, 1.0)
    return val, err, "native_adaptive_quad"


def build_quadrature_objects(
    basis_fns: list[Callable[[float], float]],
    target_fn: Callable[[float], float],
    channel_names: list[str],
    power: float,
) -> tuple[np.ndarray, np.ndarray, list[dict[str, Any]], list[dict[str, Any]]]:
    n = len(channel_names)
    gram = np.zeros((n, n), dtype=float)
    rhs = np.zeros(n, dtype=float)
    gram_rows: list[dict[str, Any]] = []
    rhs_rows: list[dict[str, Any]] = []

    for i, ci in enumerate(channel_names):
        for j, cj in enumerate(channel_names):
            def prod(x: float, ii: int = i, jj: int = j) -> float:
                return basis_fns[ii](x) * basis_fns[jj](x)

            val, err, transform = integrate_profile_product(prod, ci, cj, power)
            gram[i, j] = val
            gram_rows.append(
                {
                    "row_channel": ci,
                    "col_channel": cj,
                    "value": val,
                    "quad_error_estimate": err,
                    "transform": transform,
                }
            )

        def rhs_prod(x: float, ii: int = i) -> float:
            return basis_fns[ii](x) * target_fn(x)

        val, err, transform = integrate_profile_product(rhs_prod, ci, "R2_target_K2", power)
        rhs[i] = val
        rhs_rows.append(
            {
                "channel": ci,
                "value": val,
                "quad_error_estimate": err,
                "transform": transform,
            }
        )

    return gram, rhs, gram_rows, rhs_rows


def full_from_rank3(gram3: np.ndarray, rhs3: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    transform = np.array(
        [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, -4.0, 1.0],
        ],
        dtype=float,
    )
    return transform @ gram3 @ transform.T, transform @ rhs3


def assembled_full_rows(gram: np.ndarray, rhs: np.ndarray) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    gram_rows = []
    rhs_rows = []
    for i, ci in enumerate(CHANNEL_ORDER):
        for j, cj in enumerate(CHANNEL_ORDER):
            gram_rows.append(
                {
                    "row_channel": ci,
                    "col_channel": cj,
                    "value": float(gram[i, j]),
                    "assembly": "assembled_from_rank3_integrals_and_GB_identity",
                }
            )
        rhs_rows.append(
            {
                "channel": ci,
                "value": float(rhs[i]),
                "assembly": "assembled_from_rank3_integrals_and_GB_identity",
            }
        )
    return gram_rows, rhs_rows


def matrix_rank_with_tol(mat: np.ndarray) -> tuple[int, float, list[float], float | None]:
    singular_values = np.linalg.svd(mat, compute_uv=False)
    max_singular = float(np.max(singular_values)) if singular_values.size else 0.0
    min_singular = float(np.min(singular_values)) if singular_values.size else 0.0
    tol = max(mat.shape) * max_singular * np.finfo(float).eps * 100.0 if max_singular > 0.0 else 0.0
    rank = int(np.sum(singular_values > tol)) if max_singular > 0.0 else 0
    cond = None if min_singular <= tol else float(max_singular / min_singular)
    return rank, tol, [float(v) for v in singular_values], cond


def solve_quotient(gram: np.ndarray, rhs: np.ndarray) -> dict[str, Any]:
    independent_indices = [CHANNEL_ORDER.index(ch) for ch in INDEPENDENT_CHANNELS]
    gram3 = gram[np.ix_(independent_indices, independent_indices)]
    rhs3 = rhs[independent_indices]
    rank3, tol3, sv3, cond3 = matrix_rank_with_tol(gram3)
    coeff3 = np.linalg.solve(gram3, rhs3)
    coeff4 = np.array([coeff3[0], coeff3[1], coeff3[2], 0.0], dtype=float)
    residual4 = gram @ coeff4 - rhs
    residual3 = gram3 @ coeff3 - rhs3
    return {
        "independent_channels": list(INDEPENDENT_CHANNELS),
        "dropped_dependent_channel": "GB",
        "dependency_rule": "GB = Riem2 - 4*Ric2 + R2",
        "rank3_gram_rank": rank3,
        "rank3_gram_rank_tolerance": tol3,
        "rank3_gram_singular_values": sv3,
        "rank3_gram_condition_number": cond3,
        "rank3_quotient_coefficients": {
            "a_R2": float(coeff3[0]),
            "a_Ric2": float(coeff3[1]),
            "a_Riem2": float(coeff3[2]),
            "a_GB_representative": 0.0,
        },
        "rank3_residual_l2": float(np.linalg.norm(residual3, ord=2)),
        "rank3_residual_linf": float(np.linalg.norm(residual3, ord=np.inf)),
        "full_row_residual_l2_for_gb_zero_representative": float(np.linalg.norm(residual4, ord=2)),
        "full_row_residual_linf_for_gb_zero_representative": float(np.linalg.norm(residual4, ord=np.inf)),
    }


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p1848 = load("p1848_s798_strict_gravity_componentwise_variation_and_counterterm_witness_checkpoint.json")
    p1950 = load("p1950_s900_strict_renormalization_exact_integration_probe.json")

    omega = sp.Float("0.18575")
    phi = sp.Float("0.16250")
    beta = sp.Float("1.0")
    eta = sp.Rational(9, 5)
    d = sp.symbols("d", positive=True)
    k_expr = strict_kernel(d, omega, phi, beta, eta)
    kd_expr = sp.diff(k_expr, d)
    kdd_expr = sp.diff(kd_expr, d)
    channel_defs, backend_meta = _backend_channel_defs(d, k_expr, kd_expr, kdd_expr)

    delta_names = [row[0] for row in channel_defs]
    channel_names = [row[1] for row in channel_defs]
    if tuple(channel_names) != CHANNEL_ORDER:
        raise ValueError(f"Unexpected channel order from backend: {channel_names}")

    profile_exprs = [row[2] for row in channel_defs]
    profile_by_channel = dict(zip(channel_names, profile_exprs))
    gb_relation_expr = (
        profile_by_channel["R2"]
        - 4 * profile_by_channel["Ric2"]
        + profile_by_channel["Riem2"]
        - profile_by_channel["GB"]
    )
    symbolic_gb_relation_zero = bool(gb_relation_expr == 0)

    independent_indices = [CHANNEL_ORDER.index(ch) for ch in INDEPENDENT_CHANNELS]
    independent_exprs = [profile_exprs[idx] for idx in independent_indices]
    independent_fns = [scalar_fn(sp.lambdify(d, expr, "numpy")) for expr in independent_exprs]
    target_fn = scalar_fn(sp.lambdify(d, sp.simplify(k_expr**2), "numpy"))

    gram3, rhs3, independent_gram_rows, independent_rhs_rows = build_quadrature_objects(
        independent_fns,
        target_fn,
        list(INDEPENDENT_CHANNELS),
        PRIMARY_LEFT_POWER,
    )
    gram3_check, rhs3_check, _, _ = build_quadrature_objects(
        independent_fns,
        target_fn,
        list(INDEPENDENT_CHANNELS),
        CHECK_LEFT_POWER,
    )
    gram, rhs = full_from_rank3(gram3, rhs3)
    gram_check, rhs_check = full_from_rank3(gram3_check, rhs3_check)
    gram_rows, rhs_rows = assembled_full_rows(gram, rhs)

    adaptive_abs_gap = max(
        float(np.max(np.abs(gram - gram_check))),
        float(np.max(np.abs(rhs - rhs_check))),
    )
    adaptive_scale = max(float(np.max(np.abs(gram))), float(np.max(np.abs(rhs))), 1.0)
    adaptive_relative_gap = adaptive_abs_gap / adaptive_scale
    adaptive_transform_gate_pass = adaptive_relative_gap <= ADAPTIVE_STABILITY_TOL

    rank4, tol4, sv4, cond4 = matrix_rank_with_tol(gram)
    nullity4 = int(gram.shape[1] - rank4)
    u_full, s_full, vh_full = np.linalg.svd(gram, full_matrices=True)
    del u_full, s_full
    null_vec = np.array(vh_full[-1, :], dtype=float)
    null_vec /= max(np.linalg.norm(null_vec, ord=2), np.finfo(float).tiny)
    exact_null_norm = NULL_VECTOR_EXACT / np.linalg.norm(NULL_VECTOR_EXACT, ord=2)
    if float(np.dot(null_vec, exact_null_norm)) < 0.0:
        null_vec = -null_vec
    null_alignment = float(np.dot(null_vec, exact_null_norm))
    exact_null_residual = gram @ exact_null_norm
    exact_null_residual_l2 = float(np.linalg.norm(exact_null_residual, ord=2))
    exact_null_residual_linf = float(np.linalg.norm(exact_null_residual, ord=np.inf))

    rcond = (tol4 / max(sv4[0], np.finfo(float).tiny)) if sv4 else None
    full_min_norm_coeff, *_ = np.linalg.lstsq(gram, rhs, rcond=rcond)
    full_min_norm_residual = gram @ full_min_norm_coeff - rhs

    quotient = solve_quotient(gram, rhs)
    quotient_gate_pass = (
        quotient["rank3_gram_rank"] == 3
        and quotient["full_row_residual_linf_for_gb_zero_representative"] <= max(1.0e-8, adaptive_abs_gap * 20.0)
        and symbolic_gb_relation_zero
    )

    exact_null_square_norm = float(np.dot(NULL_VECTOR_EXACT, NULL_VECTOR_EXACT))
    canonical_rep = np.array([1.0, 0.0, 0.0, 0.0], dtype=float)
    min_norm_tau_from_canonical = -float(np.dot(canonical_rep, NULL_VECTOR_EXACT)) / exact_null_square_norm
    analytic_min_norm_from_canonical = canonical_rep + min_norm_tau_from_canonical * NULL_VECTOR_EXACT

    previous_metrics = p1950.get("aggregate_metrics") or {}
    previous_quadrature_relative_bound = previous_metrics.get("quadrature_relative_discretization_bound")
    previous_rank = previous_metrics.get("operator_profile_gram_rank")

    local_verdict = (
        "PASS_RANK3_QUOTIENT_IDENTIFIABLE_ON_SCALAR_B1_WITH_GB_NULL_TRACE"
        if quotient_gate_pass and adaptive_transform_gate_pass and rank4 == 3 and nullity4 == 1
        else "OPEN_QUOTIENT_OR_QUADRATURE_OBSTRUCTION_WITH_TRACE"
    )
    status = "OPEN_OBSTRUCTION_WITH_TRACE"
    result_kind = (
        "PARTIAL_P1950_NEXT_RANK3_QUOTIENT_PASS__FOUR_CHANNEL_GB_AND_TENSOR_PROJECTION_OPEN"
        if local_verdict.startswith("PASS")
        else "P1950_NEXT_OPEN_OBSTRUCTION_WITH_TRACE"
    )

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2027",
        "stage_id": "S977",
        "p1950_next_alias": "p1950-next",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": status,
        "result_kind": result_kind,
        "local_verdict": local_verdict,
        "route": "strict_only",
        "strict_lane_assumptions": [
            "strict_kernel_only",
            "no_legacy_transfer",
            "background_family_B1_backend_projection_channels",
            "GB_exported_as_tensor_identity_not_independent_surrogate",
        ],
        "depends_on": {
            "p1848_present": "gravity_operator_profiles_B1" in p1848,
            "p1950_probe_present": p1950.get("_missing") is None,
            "p1950_previous_rank": previous_rank,
            "p1950_previous_quadrature_relative_discretization_bound": previous_quadrature_relative_bound,
        },
        "input_hashes": {
            "p1848_json_sha256": file_sha256(GEN / "p1848_s798_strict_gravity_componentwise_variation_and_counterterm_witness_checkpoint.json"),
            "p1950_probe_json_sha256": file_sha256(GEN / "p1950_s900_strict_renormalization_exact_integration_probe.json"),
        },
        "domain": {
            "background_family": "B1",
            "epsilon_uv": EPS_UV,
            "kernel": "K_strict(d)=cos(omega*d+phi)/(1+beta*d^eta)",
            "kernel_params": {
                "omega": float(omega),
                "phi": float(phi),
                "beta": float(beta),
                "eta": float(eta),
            },
            "backend_counterterm_basis": [
                {"delta_label": delta_label, "channel": channel}
                for delta_label, channel in zip(delta_names, channel_names)
            ],
            "backend_direct_operator_profile_source_kind": backend_meta["direct_operator_profile_source_kind"],
            "backend_direct_operator_profile_generation_rule": backend_meta["direct_operator_profile_generation_rule"],
        },
        "basis_decomposition": {
            "full_channels": list(CHANNEL_ORDER),
            "full_gram_rank": rank4,
            "full_gram_nullity": nullity4,
            "full_gram_rank_tolerance": tol4,
            "full_gram_singular_values": sv4,
            "full_gram_condition_number": cond4,
            "gb_dependency_rule": "GB = Riem2 - 4*Ric2 + R2",
            "symbolic_null_relation": "R2 - 4*Ric2 + Riem2 - GB = 0",
            "symbolic_gb_relation_zero": symbolic_gb_relation_zero,
            "symbolic_gb_relation_residual": sp.sstr(gb_relation_expr),
            "exact_null_vector_R2_Ric2_Riem2_GB": [float(v) for v in NULL_VECTOR_EXACT],
            "svd_null_vector_R2_Ric2_Riem2_GB": [float(v) for v in null_vec],
            "null_vector_cosine_alignment_with_exact": null_alignment,
            "exact_null_gram_residual_l2": exact_null_residual_l2,
            "exact_null_gram_residual_linf": exact_null_residual_linf,
        },
        "quotient_identifiability": {
            **quotient,
            "quotient_gate_pass": quotient_gate_pass,
            "meaning": "The scalar B1 profile determines the rank-3 quotient representative after dropping the dependent GB direction; it does not determine an independent GB coefficient.",
        },
        "full_four_channel_min_norm_family": {
            "minimum_norm_coefficients_from_full_rank_deficient_lstsq": {
                name: float(value) for name, value in zip(FULL_COEFFICIENT_NAMES, full_min_norm_coeff)
            },
            "minimum_norm_residual_l2": float(np.linalg.norm(full_min_norm_residual, ord=2)),
            "minimum_norm_residual_linf": float(np.linalg.norm(full_min_norm_residual, ord=np.inf)),
            "canonical_rank3_representative": {
                "a_R2": 1.0,
                "a_Ric2": 0.0,
                "a_Riem2": 0.0,
                "a_GB": 0.0,
            },
            "coefficient_family_rule": "a_full = canonical_rank3_representative + tau*(1,-4,1,-1)",
            "minimum_norm_tau_from_canonical": min_norm_tau_from_canonical,
            "analytic_minimum_norm_from_canonical": {
                name: float(value) for name, value in zip(FULL_COEFFICIENT_NAMES, analytic_min_norm_from_canonical)
            },
            "gb_coefficient_identified": False,
        },
        "adaptive_quadrature": {
            "policy": "scipy.integrate.quad; entries involving Riem2 or GB use x=epsilon+(1-epsilon)*t**power to resolve Kdd**2 near the UV endpoint",
            "epsabs": QUAD_EPSABS,
            "epsrel": QUAD_EPSREL,
            "limit": QUAD_LIMIT,
            "primary_left_endpoint_power": PRIMARY_LEFT_POWER,
            "check_left_endpoint_power": CHECK_LEFT_POWER,
            "adaptive_abs_gap_primary_vs_check": adaptive_abs_gap,
            "adaptive_scale": adaptive_scale,
            "adaptive_relative_gap_primary_vs_check": adaptive_relative_gap,
            "adaptive_relative_stability_tolerance": ADAPTIVE_STABILITY_TOL,
            "adaptive_transform_gate_pass": adaptive_transform_gate_pass,
            "independent_matrix_entries_directly_integrated": independent_gram_rows,
            "independent_rhs_entries_directly_integrated": independent_rhs_rows,
            "assembled_full_matrix_entries": gram_rows,
            "assembled_full_rhs_entries": rhs_rows,
        },
        "gatekeeper_checks": {
            "p1848_direct_profiles_present": backend_meta["direct_operator_profiles_present"],
            "symbolic_gb_relation_zero": symbolic_gb_relation_zero,
            "full_gram_rank_eq_3": rank4 == 3,
            "full_gram_nullity_eq_1": nullity4 == 1,
            "rank3_quotient_gate_pass": quotient_gate_pass,
            "adaptive_transform_gate_pass": adaptive_transform_gate_pass,
            "no_four_channel_closure_claimed": True,
            "no_toe_closure_claimed": True,
        },
        "theorem_scope": {
            "licensed_statement": "On the B1 scalar projection, the direct p1848 profiles split into an identifiable rank-3 quotient plus the explicit GB null direction R2-4*Ric2+Riem2-GB.",
            "not_licensed": [
                "independent identification of a_GB on the scalar B1 projection",
                "four-channel counterterm uniqueness",
                "full tensor projection closure",
                "background-family globalization beyond B1",
                "BRST/Cutkosky unitarity closure",
                "QW-2191 selector closure",
                "ToE closure",
            ],
        },
        "false_pass_guard": "The rank-3 quotient pass is not full renormalization closure: GB remains a null/topological direction in scalar B1, and a full tensor projection is still required if the goal is four-channel coefficient identification.",
        "next_honest_step": "Lift the p1848 profile projection from scalar B1 to a tensor-resolved background family, or explicitly formulate the renormalization theorem in the quotient by the GB/topological null direction.",
        "lay_explanation": "Po usunieciu zaleznego kierunku GB trzy niezalezne profile B1 sa identyfikowalne, a nowa kwadratura stabilizuje osobliwy profil Kdd**2. Nadal nie ma dowodu na niezalezny czwarty wspolczynnik GB ani na pelna tensorowa renormalizacje.",
        "environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": __import__("scipy").__version__,
            "sympy": sp.__version__,
        },
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    md = f"""# P2027 S977 p1950-next Rank-3 / GB-Null Adaptive Quadrature Witness

Status: `{payload['status']}`

Local verdict: `{payload['local_verdict']}`

This packet consumes the direct P1848 B1 operator profiles and keeps the strict-only
kernel lane.  The four scalar profiles obey the explicit relation

`R2 - 4*Ric2 + Riem2 - GB = 0`.

Therefore the full B1 scalar Gram matrix has rank `{rank4}` and nullity `{nullity4}`.
The exact null vector in `(R2,Ric2,Riem2,GB)` coordinates is
`(1,-4,1,-1)`.

## Rank-3 Quotient

Independent representative channels: `R2, Ric2, Riem2`.

Rank-3 quotient coefficients:
- a_R2={quotient['rank3_quotient_coefficients']['a_R2']:.16e}
- a_Ric2={quotient['rank3_quotient_coefficients']['a_Ric2']:.16e}
- a_Riem2={quotient['rank3_quotient_coefficients']['a_Riem2']:.16e}
- a_GB_representative=0.0000000000000000e+00

Full-row residual for this GB-zero representative:
`{quotient['full_row_residual_linf_for_gb_zero_representative']:.3e}` in L-infinity norm.

The equivalent four-channel family is

`a_full = (1,0,0,0) + tau*(1,-4,1,-1)`.

The minimum-norm family point has `tau={min_norm_tau_from_canonical:.16e}`.

## Adaptive Quadrature

P1950 used a Simpson refinement replay and exposed a large discretization
bound.  P2027 replaces that gate with adaptive `scipy.integrate.quad`.
Any matrix/RHS entry involving `Riem2` or `GB` uses the left-endpoint transform

`x = epsilon + (1-epsilon)*t**power`

to resolve the `Kdd**2` endpoint profile.

Primary power: `{PRIMARY_LEFT_POWER:g}`.
Check power: `{CHECK_LEFT_POWER:g}`.
Relative primary/check gap: `{adaptive_relative_gap:.3e}`.
Gate tolerance: `{ADAPTIVE_STABILITY_TOL:.3e}`.
Adaptive gate pass: `{adaptive_transform_gate_pass}`.

## Honest Interpretation

The scalar B1 quotient is locally identifiable after removing the GB null
direction.  This does not identify an independent `a_GB`, does not prove
four-channel counterterm uniqueness, and does not replace the need for a
tensor-resolved projection if the theorem target requires the full curvature
operator basis.
"""
    MD.write_text(md, encoding="utf-8")
    print(OUT)
    print(MD)


if __name__ == "__main__":
    main()
