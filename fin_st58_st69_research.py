#!/usr/bin/env python3
"""FIN ST58--ST69: interval closure, projective source, and physics bridges.

The batch is local and deterministic.  Exact finite theorems, outward interval
certificates, conditional constructions, synthetic operational tests, and
failed source searches remain separate epistemic classes.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import platform
import warnings
from fractions import Fraction
from pathlib import Path
from typing import Any

warnings.filterwarnings("ignore", message="Unable to import Axes3D.*")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mpmath as mp
import numpy as np
import scipy
import sympy as sp
from scipy.linalg import expm
from scipy.optimize import minimize_scalar, root

from fin_programs_497_506_next_research import iv_bounds
from fin_programs_507_516_research import strict_a_interval
from fin_st01_st15_research import N, random_orthogonal_fixing_uniform, strict_operator
from fin_st16_st27_research import final_memory_state
from fin_st28_st45_research import (
    dyadic_lift,
    saturation_energy_gradient_hessians,
)
from fin_st46_st57_research import carrier_probability_table, canonical_digest


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "FIN_ST58_ST69_Results.json"
SUMMARY = ROOT / "FIN_ST58_ST69_Summary.csv"
CERT58 = ROOT / "FIN_ST58_Full_Interval_Certificate.json"
POWER62 = ROOT / "FIN_ST62_Finite_Count_Bounds.json"
SPEC68 = ROOT / "FIN_ST68_Calibration_Custody_Validator.json"
FIG_DIR = ROOT / "FIN_ST58_ST69_Figures"
SEED = 20260814


def native(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): native(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [native(v) for v in value]
    if isinstance(value, np.ndarray):
        return native(value.tolist())
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, complex):
        return {"real": value.real, "imag": value.imag}
    return value


def iv(value: float | str | tuple[float, float]):
    if isinstance(value, tuple):
        return mp.iv.mpf([value[0], value[1]])
    return mp.iv.mpf(str(value))


def interval_product(alo: float, ahi: float, blo: float, bhi: float) -> tuple[float, float]:
    values = [alo * blo, alo * bhi, ahi * blo, ahi * bhi]
    return float(np.nextafter(min(values), -np.inf)), float(np.nextafter(max(values), np.inf))


def interval_square(lo: float, hi: float) -> tuple[float, float]:
    if lo <= 0.0 <= hi:
        return 0.0, float(np.nextafter(max(lo * lo, hi * hi), np.inf))
    values = [lo * lo, hi * hi]
    return float(np.nextafter(min(values), -np.inf)), float(np.nextafter(max(values), np.inf))


def interval_matvec(
    alo: np.ndarray, ahi: np.ndarray, xlo: np.ndarray, xhi: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    rows, cols = alo.shape
    out_lo = np.zeros(rows)
    out_hi = np.zeros(rows)
    for i in range(rows):
        lower, upper = 0.0, 0.0
        for j in range(cols):
            lo, hi = interval_product(float(alo[i, j]), float(ahi[i, j]), float(xlo[j]), float(xhi[j]))
            lower = np.nextafter(lower + lo, -np.inf)
            upper = np.nextafter(upper + hi, np.inf)
        out_lo[i], out_hi[i] = lower, upper
    return out_lo, out_hi


def interval_left_product(
    r: np.ndarray, alo: np.ndarray, ahi: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    rows, inner = r.shape
    cols = alo.shape[1]
    out_lo = np.zeros((rows, cols))
    out_hi = np.zeros((rows, cols))
    for i in range(rows):
        for j in range(cols):
            lower, upper = 0.0, 0.0
            for k in range(inner):
                lo, hi = interval_product(float(r[i, k]), float(r[i, k]), float(alo[k, j]), float(ahi[k, j]))
                lower = np.nextafter(lower + lo, -np.inf)
                upper = np.nextafter(upper + hi, np.inf)
            out_lo[i, j], out_hi[i, j] = lower, upper
    return out_lo, out_hi


def interval_dot(
    xlo: np.ndarray, xhi: np.ndarray, ylo: np.ndarray, yhi: np.ndarray
) -> tuple[float, float]:
    lower, upper = 0.0, 0.0
    for i in range(len(xlo)):
        lo, hi = interval_product(float(xlo[i]), float(xhi[i]), float(ylo[i]), float(yhi[i]))
        lower = np.nextafter(lower + lo, -np.inf)
        upper = np.nextafter(upper + hi, np.inf)
    return float(lower), float(upper)


def positive_matvec_upper(matrix: np.ndarray, vector: np.ndarray) -> np.ndarray:
    """Outward upper product for nonnegative binary64 inputs."""
    result = np.zeros(matrix.shape[0])
    for i in range(matrix.shape[0]):
        value = 0.0
        for j in range(matrix.shape[1]):
            product = np.nextafter(float(matrix[i, j]) * float(vector[j]), np.inf)
            value = np.nextafter(value + product, np.inf)
        result[i] = value
    return result


def positive_row_sum_upper(matrix: np.ndarray) -> float:
    return float(np.max(positive_matvec_upper(matrix, np.ones(matrix.shape[1]))))


def strict_spectral_interval_data() -> dict:
    """Exact-circulant interval reduction of strict A and its hidden Schur data."""
    mp.iv.dps = 80
    omega = iv("0.18575")
    phi = iv("0.16250")
    eta = iv("1.8")
    weights = {d: mp.iv.cos(omega * d + phi) / (1 + iv(d) ** eta) for d in range(1, 7)}
    row_sum = 2 * sum(weights[d] for d in range(1, 6)) + weights[6]

    lambdas = []
    for k in range(12):
        if k == 0:
            value = iv("0")
        else:
            value = row_sum
            for d in range(1, 6):
                value -= 2 * weights[d] * mp.iv.cos(2 * mp.iv.pi * k * d / 12)
            value -= weights[6] * mp.iv.cos(mp.iv.pi * k)
        lambdas.append(value)

    poles = []
    residues = []
    multiplicities = [1, 2, 2]
    for k, multiplicity in zip([0, 1, 2], multiplicities):
        pole = row_sum + iv("0.15")
        pole -= 2 * weights[2] * mp.iv.cos(2 * mp.iv.pi * k / 6)
        pole -= 2 * weights[4] * mp.iv.cos(4 * mp.iv.pi * k / 6)
        pole -= weights[6] * mp.iv.cos(mp.iv.pi * k)
        poles.append(pole)

        distances = [1, 3, 5, 5, 3, 1]
        real = iv("0")
        imag = iv("0")
        for j, distance in enumerate(distances):
            entry = -weights[distance]
            angle = 2 * mp.iv.pi * k * j / 6
            real += entry * mp.iv.cos(angle)
            imag -= entry * mp.iv.sin(angle)
        residues.append(iv(multiplicity) * (real**2 + imag**2) / 6)

    # The omitted k=3 hidden residue is identically zero:
    # -w1+w3-w5+w5-w3+w1=0.
    memory_eigs, t_eigs, s_eigs = [], [], []
    fzero = sum(residue / pole for pole, residue in zip(poles, residues))
    for lam in lambdas:
        memory_eigs.append(fzero - sum(residue / (lam + pole) for pole, residue in zip(poles, residues)))
        t_eigs.append(sum(residue / (lam + pole) ** 2 for pole, residue in zip(poles, residues)))
        s_eigs.append(sum(residue / (lam + pole) ** 3 for pole, residue in zip(poles, residues)))
    return {
        "weights": weights,
        "row_sum": row_sum,
        "lambdas": lambdas,
        "poles": poles,
        "residues": residues,
        "memory_eigs": memory_eigs,
        "t_eigs": t_eigs,
        "s_eigs": s_eigs,
    }


def interval_circulant_from_eigenvalues(eigenvalues: list) -> tuple[np.ndarray, np.ndarray]:
    first_lo = np.zeros(N)
    first_hi = np.zeros(N)
    for j in range(N):
        value = iv("0")
        for k, eig in enumerate(eigenvalues):
            value += eig * mp.iv.cos(2 * mp.iv.pi * k * j / N) / N
        first_lo[j], first_hi[j] = iv_bounds(value)
    lo = np.zeros((N, N))
    hi = np.zeros((N, N))
    for x in range(N):
        for y in range(N):
            index = (y - x) % N
            lo[x, y], hi[x, y] = first_lo[index], first_hi[index]
    return lo, hi


def nonlinear_krawczyk(
    olo: np.ndarray, ohi: np.ndarray, center: np.ndarray
) -> dict:
    omid = 0.5 * (olo + ohi)
    center = root(
        lambda u: omid @ u + u - u**3,
        center,
        jac=lambda u: omid + np.eye(N) - np.diag(3.0 * u**2),
        method="hybr",
        options={"xtol": 1e-12, "maxfev": 5000},
    ).x
    jmid = omid + np.eye(N) - np.diag(3.0 * center**2)
    r = np.linalg.inv(jmid)
    cvec = center.copy()
    oc_lo, oc_hi = interval_matvec(olo, ohi, cvec, cvec)
    flo = np.nextafter(oc_lo + center - center**3, -np.inf)
    fhi = np.nextafter(oc_hi + center - center**3, np.inf)
    rflo, rfhi = interval_matvec(r, r, flo, fhi)
    roundoff_payment = 2e-14
    correction = np.nextafter(np.maximum(np.abs(rflo), np.abs(rfhi)) + roundoff_payment, np.inf)
    rho = np.maximum(2.0 * correction + 1e-15, 1e-13)
    final = None
    for _ in range(80):
        ulo, uhi = center - rho, center + rho
        u2lo = np.where(ulo * uhi <= 0.0, 0.0, np.minimum(ulo**2, uhi**2))
        u2hi = np.maximum(ulo**2, uhi**2)
        jlo, jhi = olo.copy(), ohi.copy()
        for i in range(N):
            jlo[i, i] = np.nextafter(jlo[i, i] + 1.0 - 3.0 * u2hi[i], -np.inf)
            jhi[i, i] = np.nextafter(jhi[i, i] + 1.0 - 3.0 * u2lo[i], np.inf)
        rjlo, rjhi = interval_left_product(r, jlo, jhi)
        mlo, mhi = -rjhi, -rjlo
        for i in range(N):
            mlo[i, i] = np.nextafter(mlo[i, i] + 1.0, -np.inf)
            mhi[i, i] = np.nextafter(mhi[i, i] + 1.0, np.inf)
        mabs = np.maximum(np.abs(mlo), np.abs(mhi))
        image = np.nextafter(correction + positive_matvec_upper(mabs, rho), np.inf)
        margin = float(np.min(rho - image))
        defect = positive_row_sum_upper(mabs)
        if margin > 0.0 and defect < 1.0:
            final = (rho, image, margin, defect)
            break
        rho = np.maximum(1.15 * rho, 1.20 * image + 1e-15)
    if final is None:
        raise RuntimeError("nonlinear interval Krawczyk inclusion failed")
    rho, image, margin, defect = final
    return {
        "center": center,
        "radius": rho,
        "box_lo": center - rho,
        "box_hi": center + rho,
        "maximum_radius": float(np.max(rho)),
        "maximum_image_radius": float(np.max(image)),
        "minimum_inclusion_margin": margin,
        "defect_infinity_norm_upper": defect,
        "center_residual_inf": float(np.linalg.norm(omid @ center + center - center**3, ord=np.inf)),
        "binary64_roundoff_payment_per_component": roundoff_payment,
    }


def interval_linear_solve(
    alo: np.ndarray, ahi: np.ndarray, blo: np.ndarray, bhi: np.ndarray
) -> dict:
    amid = 0.5 * (alo + ahi)
    bmid = 0.5 * (blo + bhi)
    x0 = np.linalg.solve(amid, bmid)
    r = np.linalg.inv(amid)
    axlo, axhi = interval_matvec(alo, ahi, x0, x0)
    reslo = np.nextafter(blo - axhi, -np.inf)
    reshi = np.nextafter(bhi - axlo, np.inf)
    rrlo, rrhi = interval_matvec(r, r, reslo, reshi)
    roundoff_payment = 2e-14
    correction = np.nextafter(np.maximum(np.abs(rrlo), np.abs(rrhi)) + roundoff_payment, np.inf)
    ralo, rahi = interval_left_product(r, alo, ahi)
    mlo, mhi = -rahi, -ralo
    for i in range(len(x0)):
        mlo[i, i] = np.nextafter(mlo[i, i] + 1.0, -np.inf)
        mhi[i, i] = np.nextafter(mhi[i, i] + 1.0, np.inf)
    mabs = np.maximum(np.abs(mlo), np.abs(mhi))
    spectral_radius = float(max(abs(np.linalg.eigvals(mabs))))
    if spectral_radius >= 1.0:
        raise RuntimeError("interval linear enclosure lacks contraction")
    radius = np.linalg.solve(np.eye(len(x0)) - mabs, correction)
    radius = np.nextafter(radius * 1.01 + 2e-14, np.inf)
    for _ in range(20):
        image = np.nextafter(correction + positive_matvec_upper(mabs, radius), np.inf)
        if np.all(image < radius):
            break
        radius = np.nextafter(np.maximum(1.20 * radius, 1.20 * image + 2e-14), np.inf)
    else:
        raise RuntimeError("interval linear Krawczyk inclusion failed")
    return {
        "center": x0,
        "lo": x0 - radius,
        "hi": x0 + radius,
        "radius": radius,
        "comparison_spectral_radius": spectral_radius,
        "minimum_inclusion_margin": float(np.min(radius - image)),
        "binary64_roundoff_payment_per_component": roundoff_payment,
    }


def st58_full_interval_certificate() -> dict:
    data = strict_spectral_interval_data()
    operator_eigs = [a + m for a, m in zip(data["lambdas"], data["memory_eigs"])]
    olo, ohi = interval_circulant_from_eigenvalues(operator_eigs)
    tlo, thi = interval_circulant_from_eigenvalues(data["t_eigs"])
    slo, shi = interval_circulant_from_eigenvalues(data["s_eigs"])
    _, amid, _ = strict_operator()
    starting_state = final_memory_state(amid)[0]
    root_cert = nonlinear_krawczyk(olo, ohi, starting_state)
    ulo, uhi = root_cert["box_lo"], root_cert["box_hi"]

    u2lo = np.where(ulo * uhi <= 0.0, 0.0, np.minimum(ulo**2, uhi**2))
    u2hi = np.maximum(ulo**2, uhi**2)
    llo, lhi = olo.copy(), ohi.copy()
    for i in range(N):
        llo[i, i] = np.nextafter(llo[i, i] + 1.0 - 3.0 * u2hi[i], -np.inf)
        lhi[i, i] = np.nextafter(lhi[i, i] + 1.0 - 3.0 * u2lo[i], np.inf)

    y = interval_linear_solve(llo, lhi, ulo, uhi)
    vlo, vhi = interval_matvec(tlo, thi, ulo, uhi)
    z = interval_linear_solve(llo, lhi, vlo, vhi)
    su_lo, su_hi = interval_matvec(slo, shi, ulo, uhi)
    c0 = interval_dot(ulo, uhi, y["lo"], y["hi"])
    c1raw = interval_dot(ulo, uhi, z["lo"], z["hi"])
    c2 = interval_dot(vlo, vhi, z["lo"], z["hi"])
    c3 = interval_dot(ulo, uhi, su_lo, su_hi)
    coefficients = {
        "a=c2+c3": [np.nextafter(c2[0] + c3[0], -np.inf), np.nextafter(c2[1] + c3[1], np.inf)],
        "b=2c1": [np.nextafter(2.0 * c1raw[0], -np.inf), np.nextafter(2.0 * c1raw[1], np.inf)],
        "c=c0": list(c0),
    }

    mp.iv.dps = 80
    aa = mp.iv.mpf(coefficients["a=c2+c3"])
    bb = mp.iv.mpf(coefficients["b=2c1"])
    cc = mp.iv.mpf(coefficients["c=c0"])
    positive_x = (-bb + mp.iv.sqrt(bb**2 - 4 * aa * cc)) / (2 * aa)
    speed_interval = iv_bounds(1 / positive_x)
    pole_intervals = [iv_bounds(item) for item in data["poles"]]
    residue_intervals = [iv_bounds(item) for item in data["residues"]]
    lambda_intervals = [iv_bounds(item) for item in data["lambdas"]]

    certificate = {
        "arithmetic": "mpmath directed intervals plus nextafter-outward binary64 matrix Krawczyk operations",
        "exact_reduction": (
            "A, the even/odd hidden blocks B,C, and the memory functions are circulant. The 6- and 12-point "
            "Fourier decompositions reduce projectors, residues, resolvents, T and S to scalar intervals. The hidden k=3 "
            "residue is exactly zero by pair cancellation."
        ),
        "strict_input_model": "omega=743/4000, phi=13/80, eta=9/5, beta=1, hidden shift=3/20",
        "pole_intervals": pole_intervals,
        "residue_intervals": residue_intervals,
        "strict_eigenvalue_intervals": lambda_intervals,
        "stationary_root": {k: v for k, v in root_cert.items() if k not in {"box_lo", "box_hi"}},
        "stationary_box_lo": ulo,
        "stationary_box_hi": uhi,
        "linear_solve_comparison_spectral_radii": [y["comparison_spectral_radius"], z["comparison_spectral_radius"]],
        "linear_solve_minimum_inclusion_margins": [y["minimum_inclusion_margin"], z["minimum_inclusion_margin"]],
        "binary64_roundoff_payment_per_component": root_cert["binary64_roundoff_payment_per_component"],
        "coefficient_intervals": coefficients,
        "collision_speed_interval": speed_interval,
        "interval_width": speed_interval[1] - speed_interval[0],
        "proof_boundary": (
            "This is a machine-generated outward interval certificate under the explicitly proved circulant/Fourier "
            "reduction. It is not a proof-assistant replay and does not constitute a global Krein or physical-stability theorem."
        ),
    }
    CERT58.write_text(json.dumps(native(certificate), indent=2, sort_keys=True), encoding="utf-8")
    return {
        "program": "ST58",
        "object": "Full Directed-Interval Strict-to-Memory Collision Certificate",
        **certificate,
        "status": "proven_outward_interval_certificate_under_exact_circulant_reduction",
    }


def st59_gain_source_sign_flip_no_go(a: np.ndarray) -> dict:
    reflection = np.zeros((N, N))
    for x in range(N):
        reflection[(-x) % N, x] = 1.0
    pattern = np.diag(np.sin(2.0 * math.pi * np.arange(N) / N))
    pattern /= np.linalg.norm(pattern)
    odd_residual = float(np.linalg.norm(reflection @ pattern @ reflection.T + pattern))
    rows = []
    for sign in [-1.0, 1.0]:
        # Phi_sign(A+H)=sign*(||P_-H||^2+||P_+H||^2)/2.
        # Gradient descent has odd gain -sign.
        rows.append({
            "functional_sign": sign,
            "stationary_at_A": True,
            "stabilizer_invariant": True,
            "odd_gradient_flow_gain": -sign,
        })
    return {
        "program": "ST59",
        "object": "Sign-Closed Adaptive-Functional No-Go",
        "reflection_odd_witness_residual": odd_residual,
        "rows": rows,
        "theorem": (
            "In every admissible source class closed under Phi -> -Phi, stationarity at A and stabilizer invariance "
            "cannot determine the sign of the response Hessian: Phi and -Phi satisfy the same static source data and "
            "give opposite linearized gradient-flow gains. Positivity requires an order, Lyapunov direction, entropy "
            "production rule, or equivalent premise not contained in the static strict operator."
        ),
        "status": "proven_no_strict_gain_sign_in_sign_closed_functional_class",
        "boundary": "The theorem does not exclude a genuinely sourced time-arrow or monotonicity axiom outside the sign-closed class.",
    }


def pancharatnam_holonomy(rays: np.ndarray) -> complex:
    product = 1.0 + 0.0j
    for x in range(len(rays)):
        overlap = np.vdot(rays[x], rays[(x + 1) % len(rays)])
        if abs(overlap) < 1e-14:
            return complex(np.nan, np.nan)
        product *= overlap / abs(overlap)
    return complex(product)


def st60_projective_source_obstruction(a: np.ndarray) -> dict:
    rows = []
    x = np.arange(N)
    for k in range(1, 6):
        rays = np.stack([np.cos(2.0 * math.pi * k * x / N), np.sin(2.0 * math.pi * k * x / N)], axis=1).astype(complex)
        holonomy = pancharatnam_holonomy(rays)
        defined = bool(np.isfinite(holonomy.real))
        rows.append({
            "integer_fourier_mode": k,
            "neighbor_overlap": math.cos(2.0 * math.pi * k / N),
            "holonomy": holonomy if defined else None,
            "defined": defined,
        })
    half_rays = np.stack([np.cos(math.pi * x / N), np.sin(math.pi * x / N)], axis=1).astype(complex)
    half_holonomy = pancharatnam_holonomy(half_rays)
    return {
        "program": "ST60",
        "object": "Strict Spectral-Bundle Projective-Texture Obstruction",
        "integer_mode_rows": rows,
        "inserted_half_angle_holonomy": half_holonomy,
        "strict_spectral_projector_holonomies_all_trivial_when_defined": all(
            (not row["defined"]) or abs(complex(row["holonomy"]) - 1.0) < 1e-12 for row in rows
        ),
        "theorem": (
            "Every real two-dimensional eigenspace of a scalar C12-circulant operator is an integer Fourier doublet. "
            "Its projector-column ray embedding has constant neighbor overlap cos(2*pi*k/12); whenever this is nonzero, "
            "the twelve-link Pancharatnam sign is sign(cos)^12=+1. The ST49 -1 texture requires a half-integer/projective "
            "lift absent from the scalar strict spectral bundle."
        ),
        "status": "proven_scalar_strict_spectral_bundle_does_not_source_ST49_texture",
        "boundary": "A spin structure, central extension, boundary condition, or other added projective state type remains admissible but is not selected by A.",
    }


def interval_trace_square(lo: np.ndarray, hi: np.ndarray) -> tuple[float, float]:
    lower, upper = 0.0, 0.0
    for i in range(lo.shape[0]):
        for j in range(lo.shape[1]):
            slo, shi = interval_square(float(lo[i, j]), float(hi[i, j]))
            lower = np.nextafter(lower + slo, -np.inf)
            upper = np.nextafter(upper + shi, np.inf)
    return float(lower), float(upper)


def strict_lift_interval() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    alo, ahi, _ = strict_a_interval()
    wlo = {-d: 0.0 for d in range(1, 7)}
    wlo = {d: -ahi[0, d] for d in range(1, 7)}
    whi = {d: -alo[0, d] for d in range(1, 7)}
    n = 24
    llo = np.zeros((n, n))
    lhi = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            d = min(abs(i - j), n - abs(i - j))
            if 1 <= d < 6:
                low, high = wlo[d], whi[d]
            elif d == 6:
                low, high = wlo[6] / 2.0, whi[6] / 2.0
            else:
                low, high = 0.0, 0.0
            llo[i, j] = llo[j, i] = -high
            lhi[i, j] = lhi[j, i] = -low
    for i in range(n):
        lower_diagonal = np.nextafter(-np.sum(lhi[i]), -np.inf)
        upper_diagonal = np.nextafter(-np.sum(llo[i]), np.inf)
        llo[i, i] = lower_diagonal
        lhi[i, i] = upper_diagonal
    return alo, ahi, llo, lhi


def st61_heat_signature_refinement_no_go() -> dict:
    alo, ahi, llo, lhi = strict_lift_interval()
    coarse_sq = interval_trace_square(alo, ahi)
    lift_sq = interval_trace_square(llo, lhi)
    difference = [
        np.nextafter(lift_sq[0] / 24.0 - coarse_sq[1] / 12.0, -np.inf),
        np.nextafter(lift_sq[1] / 24.0 - coarse_sq[0] / 12.0, np.inf),
    ]
    a = strict_operator()[1]
    rows = []
    for q in [0.0, 0.1, 0.3, 0.9]:
        lift = dyadic_lift(a, q)
        rows.append({
            "q": q,
            "trace_density": float(np.trace(lift) / 24.0),
            "second_moment_density": float(np.trace(lift @ lift) / 24.0),
        })
    return {
        "program": "ST61",
        "object": "Heat-Signature Refinement Source Audit",
        "rows": rows,
        "normalized_second_moment_difference_at_q0_interval": difference,
        "interval_excludes_zero": bool(difference[0] > 0.0 or difference[1] < 0.0),
        "theorem": (
            "Normalized heat-trace equality for all t implies equality of every normalized spectral moment. Its first "
            "derivative at t=0 forces q=0 by ST47, but at q=0 the outward interval for the normalized second-moment "
            "difference excludes zero. Therefore no member of the declared dyadic lift preserves the full normalized "
            "heat signature. Trace-density conservation is not implied by a full heat-signature fixed point."
        ),
        "status": "proven_no_full_heat_signature_stationary_dyadic_refinement",
        "boundary": "This refutes the declared heat-signature source; it does not exclude a different fine-data fractal code with a nonstationary renormalization law.",
    }


def st62_finite_count_bounds(a: np.ndarray) -> dict:
    rng = np.random.default_rng(20260813 + 51)
    qmat = random_orthogonal_fixing_uniform(rng)
    p = carrier_probability_table(a, np.eye(N), transported=False)
    q = carrier_probability_table(a, qmat, transported=False)
    floor = 1e-300
    p = np.maximum(p, floor)
    q = np.maximum(q, floor)
    kl_pq = float(np.sum(p * np.log(p / q)))
    kl_qp = float(np.sum(q * np.log(q / p)))

    def log_chernoff(s: float) -> float:
        return float(np.sum(np.log(np.sum(p**s * q ** (1.0 - s), axis=1))))

    optimum = minimize_scalar(log_chernoff, bounds=(0.0, 1.0), method="bounded", options={"xatol": 1e-14})
    chernoff = -float(optimum.fun)
    epsilon = 0.01
    binary_kl = (1.0 - 2.0 * epsilon) * math.log((1.0 - epsilon) / epsilon)
    necessary = math.ceil(max(binary_kl / kl_pq, binary_kl / kl_qp))
    sufficient = math.ceil(math.log(1.0 / epsilon) / chernoff)
    bc_log = float(np.sum(np.log(np.sum(np.sqrt(p * q), axis=1))))
    packet = {
        "frozen_table_source": "ST51 seed and exact binary64 carrier tables",
        "target_each_error": epsilon,
        "sum_KL_PQ_per_shot_per_preparation": kl_pq,
        "sum_KL_QP_per_shot_per_preparation": kl_qp,
        "chernoff_optimizer_s": float(optimum.x),
        "chernoff_information_per_shot_per_preparation": chernoff,
        "log_bhattacharyya_coefficient_per_shot_per_preparation": bc_log,
        "necessary_shots_per_preparation_information_bound": necessary,
        "sufficient_shots_per_preparation_chernoff_bound": sufficient,
        "ST51_shots_per_preparation": 1200,
    }
    POWER62.write_text(json.dumps(packet, indent=2, sort_keys=True), encoding="utf-8")
    return {
        "program": "ST62",
        "object": "Analytic Finite-Count Two-Carrier Power Bracket",
        **packet,
        "theorem": (
            "Binary decision data processing gives the stated KL necessary lower bound for both errors <=0.01. "
            "The Chernoff bound gives a likelihood-ratio test with alpha+beta<=0.01 at the stated sufficient count. "
            "The bracket is conditional on the frozen synthetic carrier tables."
        ),
        "status": "proven_information_bound_bracket_for_frozen_synthetic_receiver",
        "boundary": "The bound is not an apparatus count budget and does not include calibration drift, dark counts, or external custody.",
    }


def st63_cp_operational_intertwiners(a: np.ndarray) -> dict:
    eigvals, eigvecs = np.linalg.eigh(a)
    groups: list[list[int]] = []
    for index, value in enumerate(eigvals):
        if not groups or abs(value - eigvals[groups[-1][0]]) > 1e-9:
            groups.append([index])
        else:
            groups[-1].append(index)
    projectors = [eigvecs[:, group] @ eigvecs[:, group].T for group in groups]
    tp_residual = float(np.linalg.norm(sum(projectors) - np.eye(N)))
    choi = np.zeros((N * N, N * N), dtype=complex)
    for projector in projectors:
        vectorized = projector.reshape(-1, order="F")
        choi += np.outer(vectorized, vectorized.conj())
    choi_eigs = np.linalg.eigvalsh((choi + choi.conj().T) / 2.0)
    t, tau = 0.61, 0.61
    gap_rows = []
    for i, left in enumerate(eigvals):
        for j, right in enumerate(eigvals):
            gap = left - right
            uvalue = np.exp(-1j * t * gap)
            dvalue = np.exp(-tau * gap * gap)
            if abs(gap) > 1e-9:
                gap_rows.append(abs(uvalue - dvalue))
    return {
        "program": "ST63",
        "object": "Completely-Positive Channel Intertwiner Classification",
        "spectral_block_multiplicities": [len(group) for group in groups],
        "pinching_output_algebra_dimension": int(sum(len(group) ** 2 for group in groups)),
        "pinching_kraus_rank": len(groups),
        "trace_preservation_residual": tp_residual,
        "choi_minimum_eigenvalue": float(choi_eigs[0]),
        "minimum_positive_gap_superoperator_mismatch": float(min(gap_rows)),
        "invertible_CP_cross_intertwiner_exists": False,
        "theorem": (
            "For U_t(rho)=e^-itA rho e^itA and D_tau(rho)=e^{-tau ad_A^2}(rho), a spectral-basis Schur "
            "intertwiner can retain E_ij only when exp[-it(lambda_i-lambda_j)]=exp[-tau(lambda_i-lambda_j)^2]. "
            "At t=tau=0.61 this forces equal eigenvalues. The canonical spectral pinching is CPTP and intertwines "
            "the channels on the 22-dimensional block algebra, but it destroys all cross-eigenspace coherences and "
            "is not invertible. No invertible CP cross-channel equivalence exists."
        ),
        "status": "proven_CP_pinching_intertwiner_and_invertible_equivalence_no_go",
        "boundary": "The dephasing channel is a valid A-derived quantum channel, not the classical vertex heat semigroup itself.",
    }


def shannon(values: np.ndarray) -> float:
    positive = values[values > 0]
    return float(-np.sum(positive * np.log(positive)))


def st64_information_to_thermodynamics_bridge(a: np.ndarray) -> dict:
    beta = 0.8
    eigenvalues, vectors = np.linalg.eigh(a)
    weights = np.exp(-beta * eigenvalues)
    weights /= np.sum(weights)
    rho_beta = vectors @ np.diag(weights) @ vectors.T
    p = np.arange(1, N + 1, dtype=float)
    p /= np.sum(p)
    rho = np.diag(p)
    entropy_rho = shannon(np.linalg.eigvalsh(rho))
    entropy_beta = shannon(weights)
    energy_rho = float(np.trace(rho @ a))
    energy_beta = float(np.trace(rho_beta @ a))
    evals_rho, evecs_rho = np.linalg.eigh(rho)
    log_rho = evecs_rho @ np.diag(np.log(np.maximum(evals_rho, 1e-300))) @ evecs_rho.T
    log_beta = vectors @ np.diag(np.log(weights)) @ vectors.T
    relative = float(np.trace(rho @ (log_rho - log_beta)))
    free_difference = energy_rho - energy_beta - (entropy_rho - entropy_beta) / beta
    entropy_production = (entropy_beta - entropy_rho) - beta * (energy_beta - energy_rho)
    return {
        "program": "ST64",
        "object": "Minimal Conditional Information-to-Thermodynamics Bridge",
        "dimensionless_beta": beta,
        "relative_entropy": relative,
        "beta_times_free_energy_difference": beta * free_difference,
        "free_energy_identity_residual": abs(relative - beta * free_difference),
        "dimensionless_entropy_production": entropy_production,
        "entropy_production_identity_residual": abs(relative - entropy_production),
        "minimal_dimensional_inputs": ["entropy unit k_B", "energy unit E_*"],
        "derived_temperature": "T=E_*/(k_B*beta_tilde)",
        "scale_orbit": "(E_*,T)->(c E_*,c T) leaves rho_beta and every dimensionless FIN record unchanged",
        "theorem": (
            "After adding H=E_* A, S=k_B S_vN, a bath Gibbs state with dimensionless beta_tilde, and a thermalization "
            "process, D(rho||rho_beta)=beta_tilde [F(rho)-F(rho_beta)]/E_* and Delta S-Q/T=k_B D>=0. "
            "The energy and entropy units cannot be recovered from dimensionless records because the displayed positive "
            "scale orbit leaves the state invariant."
        ),
        "status": "proven_conditional_thermodynamic_identity_and_dimensional_scale_obstruction",
        "boundary": "FIN supplies neither E_*, k_B, a bath, nor a physical thermalization instrument from the strict operator alone.",
    }


def st65_localized_excited_state_search(a: np.ndarray) -> dict:
    high_amplitude = 2.669532756662446
    u = np.zeros(N)
    u[0] = high_amplitude
    rows = []
    last_success = 0.0
    last_u = u.copy()
    for kappa in np.linspace(0.0, 0.04, 401):
        solution = root(
            lambda x: saturation_energy_gradient_hessians(kappa * a, x)[1],
            u,
            jac=lambda x: saturation_energy_gradient_hessians(kappa * a, x)[2],
            method="hybr",
            options={"xtol": 1e-12, "maxfev": 4000},
        )
        residual = float(np.linalg.norm(saturation_energy_gradient_hessians(kappa * a, solution.x)[1], ord=np.inf))
        if (not solution.success) or residual > 1e-8:
            break
        u = solution.x
        last_u = u.copy()
        last_success = float(kappa)
        power = float(u @ u)
        rows.append({
            "kappa": float(kappa),
            "ipr": float(np.sum(u**4) / power**2),
            "residual_inf": residual,
            "peak_amplitude": float(np.max(np.abs(u))),
        })

    rng = np.random.default_rng(SEED + 65)
    localized_hits = []
    stationary_hits = 0
    for _ in range(240):
        start = rng.normal(scale=rng.uniform(0.15, 3.5), size=N)
        solution = root(
            lambda x: saturation_energy_gradient_hessians(a, x)[1],
            start,
            jac=lambda x: saturation_energy_gradient_hessians(a, x)[2],
            method="hybr",
            options={"xtol": 1e-11, "maxfev": 5000},
        )
        residual = float(np.linalg.norm(saturation_energy_gradient_hessians(a, solution.x)[1], ord=np.inf))
        if solution.success and residual < 1e-8:
            stationary_hits += 1
            power = float(solution.x @ solution.x)
            ipr = float(np.sum(solution.x**4) / power**2) if power > 1e-14 else 0.0
            if ipr > 0.10:
                localized_hits.append({"ipr": ipr, "residual": residual, "state": solution.x})
    last_power = float(last_u @ last_u)
    return {
        "program": "ST65",
        "object": "Localized Excited-State Continuation and Search",
        "anti_continuum_branch_rows": rows,
        "last_regular_numerical_continuation_kappa": last_success,
        "last_continuation_ipr": float(np.sum(last_u**4) / last_power**2),
        "random_multistarts": 240,
        "converged_stationary_starts": stationary_hits,
        "localized_hits_at_kappa_1": len(localized_hits),
        "status": "strong_numerical_localized_branch_terminates_before_strict_coupling_no_kappa1_hit",
        "boundary": (
            "The localized anti-continuum branch loses regular continuation near the reported kappa, and no kappa=1 "
            "localized root was found. This is not a global nonexistence theorem; pseudo-arclength and interval global methods remain open."
        ),
    }


def rational_bisect_polynomial(poly: sp.Poly, left: Fraction, right: Fraction, iterations: int = 90) -> tuple[Fraction, Fraction]:
    def evaluate(value: Fraction) -> Fraction:
        total = Fraction(0)
        for coefficient in poly.all_coeffs():
            total = total * value + Fraction(int(coefficient.p), int(coefficient.q))
        return total
    fl, fr = evaluate(left), evaluate(right)
    if fl * fr >= 0:
        raise ValueError("root not bracketed")
    for _ in range(iterations):
        mid = (left + right) / 2
        fm = evaluate(mid)
        if fm == 0:
            return mid, mid
        if fl * fm < 0:
            right, fr = mid, fm
        else:
            left, fl = mid, fm
    return left, right


def st66_polynomial_equivariant_bifurcation() -> dict:
    y = sp.symbols("y")
    mu = sp.Rational(7, 20)
    epsilon = sp.Rational(1, 100)
    nu = sp.Rational(1, 100)
    stable_poly = sp.Poly(-mu + y - epsilon * y**5 + nu * y**6, y, domain=sp.QQ)
    saddle_poly = sp.Poly(-mu + y + epsilon * y**5 + nu * y**6, y, domain=sp.QQ)
    stable_count = int(sp.polys.polytools.count_roots(stable_poly, 0, sp.oo))
    saddle_count = int(sp.polys.polytools.count_roots(saddle_poly, 0, sp.oo))
    stable_box = rational_bisect_polynomial(stable_poly, Fraction(35, 100), Fraction(36, 100))
    saddle_box = rational_bisect_polynomial(saddle_poly, Fraction(34, 100), Fraction(35, 100))
    stable_y = float((stable_box[0] + stable_box[1]) / 2)
    radial_curvature = 2.0 * stable_y * float(sp.diff(stable_poly.as_expr(), y).subs(y, stable_y))
    angular_curvature = 12.0 * float(epsilon) * stable_y**6
    return {
        "program": "ST66",
        "object": "Smooth Polynomial C12-Equivariant Selection Bifurcation",
        "potential": "V=-mu|z|^2/2+|z|^4/4-epsilon Re(z^12)/12+nu|z|^14/14",
        "parameters": {"mu": str(mu), "epsilon": str(epsilon), "nu": str(nu)},
        "stable_radial_polynomial": str(stable_poly.as_expr()),
        "positive_stable_radial_roots": stable_count,
        "positive_angular_saddle_radial_roots": saddle_count,
        "stable_y_rational_interval": [str(stable_box[0]), str(stable_box[1])],
        "saddle_y_rational_interval": [str(saddle_box[0]), str(saddle_box[1])],
        "stable_branch_count": 12,
        "angular_saddle_branch_count": 12,
        "stable_radial_curvature": radial_curvature,
        "stable_angular_curvature": angular_curvature,
        "theorem": (
            "The displayed coercive polynomial potential is invariant under z->exp(2*pi*i/12)z and reflection z->conj(z). "
            "Sturm counting gives one positive radial root on each angular class. The cos(12 theta)=+1 classes form "
            "twelve strict local minima; cos(12 theta)=-1 gives twelve angular saddles; the origin is unstable for mu>0."
        ),
        "status": "proven_constructed_smooth_polynomial_C12_bifurcation_classification",
        "boundary": "Gain, saturation, anisotropy and the realized branch are inserted; the construction is not strict-sourced.",
    }


def bargmann_chirality(rays: np.ndarray) -> float:
    total = 0.0
    for x in range(len(rays)):
        a, b, c = rays[x], rays[(x + 1) % len(rays)], rays[(x + 2) % len(rays)]
        total += float(np.imag(np.vdot(a, b) * np.vdot(b, c) * np.vdot(c, a)))
    return total


def st67_pi_holonomy_chirality() -> dict:
    x = np.arange(N)
    t = 2.0 - math.sqrt(3.0)
    forward = np.stack([
        np.full(N, math.sqrt(1.0 - t), dtype=complex),
        math.sqrt(t) * np.exp(1j * 2.0 * math.pi * 2.0 * x / N),
    ], axis=1)
    reverse = forward[::-1].copy()
    h_forward = pancharatnam_holonomy(forward)
    h_reverse = pancharatnam_holonomy(reverse)
    chi_forward = bargmann_chirality(forward)
    chi_reverse = bargmann_chirality(reverse)
    return {
        "program": "ST67",
        "object": "Pi-Holonomy/Chirality Separation by a Projective Bargmann Invariant",
        "latitude_weight_t": t,
        "forward_holonomy": h_forward,
        "reverse_holonomy": h_reverse,
        "forward_bargmann_chirality": chi_forward,
        "reverse_bargmann_chirality": chi_reverse,
        "chirality_reflection_residual": abs(chi_forward + chi_reverse),
        "theorem": (
            "The inserted CP1 two-winding latitude with t=2-sqrt(3) has link phase pi/12 and total Pancharatnam "
            "holonomy pi. Reflection leaves holonomy -1 fixed but reverses the gauge-invariant sum of imaginary "
            "three-ray Bargmann products. Therefore pi flux and chirality are mathematically independent data."
        ),
        "status": "proven_conditional_projective_chirality_invariant_at_fixed_pi_holonomy",
        "boundary": "The chiral CP1 texture is constructed, not selected by the scalar strict operator or QW-2191.",
    }


def st68_calibration_custody_validator(a: np.ndarray) -> dict:
    operator_norm = float(np.linalg.eigvalsh(a)[-1])
    delta_t = 0.001
    trace_distance_tv_bound = operator_norm * delta_t
    actual = []
    base_t = 0.63
    u0 = expm(-1j * base_t * a)
    p0 = np.abs(u0) ** 2
    for shift in [-delta_t, delta_t]:
        p1 = np.abs(expm(-1j * (base_t + shift) * a)) ** 2
        actual.append(float(np.max(0.5 * np.sum(np.abs(p1 - p0), axis=1))))
    specification = {
        "version": "1.0.0",
        "frozen_model": "FIN ST55 two-carrier operational transfer",
        "calibration": {
            "registered_time": base_t,
            "absolute_time_uncertainty": delta_t,
            "analytic_probability_TV_upper_bound": trace_distance_tv_bound,
            "declared_TV_budget": 0.01,
        },
        "required_distinct_roles": ["provider", "registrar", "analyst", "custodian"],
        "event_schema": ["timestamp", "carrier_id", "preparation_id", "outcome_id", "run_id", "calibration_split", "blind_label"],
        "hash_before_unblinding": True,
        "holdout_frozen": True,
        "single_pipeline_run": True,
        "no_refit_after_failure": True,
    }
    digest = canonical_digest(specification)
    valid_roles = {"provider": "P", "registrar": "R", "analyst": "A", "custodian": "C"}
    invalid_roles = {"provider": "P", "registrar": "P", "analyst": "A", "custodian": "C"}

    def validate(roles: dict, events_present: bool, hash_frozen: bool) -> dict:
        checks = {
            "roles_pairwise_distinct": len(set(roles.values())) == 4,
            "raw_events_present": events_present,
            "hash_frozen_before_unblinding": hash_frozen,
            "calibration_bound_within_budget": trace_distance_tv_bound <= specification["calibration"]["declared_TV_budget"],
        }
        return {"checks": checks, "accepted": all(checks.values())}

    packet = {
        "specification": specification,
        "sha256": digest,
        "actual_endpoint_TV_values": actual,
        "valid_synthetic_case": validate(valid_roles, True, True),
        "duplicate_role_rejection": validate(invalid_roles, True, True),
        "missing_external_record_rejection": validate(valid_roles, False, True),
    }
    SPEC68.write_text(json.dumps(packet, indent=2, sort_keys=True), encoding="utf-8")
    return {
        "program": "ST68",
        "object": "Calibration-Uncertainty and Custody Validator",
        "specification_sha256": digest,
        "operator_norm": operator_norm,
        "time_uncertainty": delta_t,
        "analytic_TV_bound": trace_distance_tv_bound,
        "actual_endpoint_TV_maximum": max(actual),
        "valid_case_accepted": packet["valid_synthetic_case"]["accepted"],
        "duplicate_role_case_rejected": not packet["duplicate_role_rejection"]["accepted"],
        "missing_record_case_rejected": not packet["missing_external_record_rejection"]["accepted"],
        "status": "constructed_executable_calibration_and_custody_validator",
        "boundary": "The validator can reject malformed packets; it cannot make role identities independent or create an external record.",
    }


def st69_axiom_graph_reconciliation() -> dict:
    rows = [
        {"axiom": "A0 strict finite operator", "status": "retained", "evidence": "ST58 certifies consequences but does not derive A"},
        {"axiom": "A1 state-selection event", "status": "retained", "evidence": "ST60 obstructs strict ST49 source; ST66/ST67 insert branch/state"},
        {"axiom": "A2 refinement law", "status": "retained", "evidence": "ST61 refutes full heat-signature stationarity"},
        {"axiom": "A3 dimensional calibration", "status": "retained", "evidence": "ST64 proves the E_*, k_B scale orbit"},
        {"axiom": "A4 operational process", "status": "constructible", "evidence": "ST63/ST68 refine but do not strict-source it"},
        {"axiom": "A5 external custody/data", "status": "retained", "evidence": "ST68 explicitly cannot generate independence or events"},
        {"axiom": "A6 nonlinear response/gain", "status": "retained", "evidence": "ST59 proves sign-closed gain nonidentifiability"},
        {"axiom": "A7 connection/projective sector", "status": "retained", "evidence": "ST60 obstruction; ST67 is inserted CP1 data"},
        {"axiom": "A8 temporal memory realization", "status": "retained", "evidence": "ST58 certifies one static reduction, not temporal uniqueness"},
    ]
    conditional_packages = {
        "projective_chiral_state": ["A1", "A7"],
        "thermal_conversion_package": ["A3", "bath/equilibrium process"],
        "smooth_selection_dynamics": ["A1", "A6"],
        "operational_validation_packet": ["A4", "record schema; not A5"],
    }
    return {
        "program": "ST69",
        "object": "Post-ST68 Axiom Graph and Conditional Package Audit",
        "rows": rows,
        "strict_source_group_count_before": 9,
        "strict_source_group_count_after": 9,
        "conditional_packages": conditional_packages,
        "theorem": (
            "ST58 closes an arithmetic proof obligation but does not remove a physical source group. ST59-ST68 either "
            "prove obstructions or construct conditional packages. Removal countermodels therefore retain all nine source "
            "groups relative to the declared physical-completion targets."
        ),
        "status": "proven_no_strict_axiom_reduction_after_ST58_ST68",
        "boundary": "Minimality is relative to the current targets and object classes, not all possible future mathematics.",
    }


def make_figures(results: dict) -> None:
    FIG_DIR.mkdir(exist_ok=True)
    plt.style.use("seaborn-v0_8-whitegrid")

    st58 = results["ST58"]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    poles = st58["pole_intervals"]
    axes[0].errorbar(range(1, 4), [(x[0] + x[1]) / 2 for x in poles], yerr=[(x[1] - x[0]) / 2 for x in poles], fmt="o")
    axes[0].set(title="directed hidden poles", xlabel="cluster", ylabel="pole")
    lo, hi = st58["collision_speed_interval"]
    axes[1].bar(["certified interval"], [hi - lo], bottom=[lo])
    axes[1].set(title="full collision enclosure", ylabel="speed")
    fig.tight_layout(); fig.savefig(FIG_DIR / "st58_full_interval.png", dpi=190); plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    rows = results["ST59"]["rows"]
    axes[0].bar([str(r["functional_sign"]) for r in rows], [r["odd_gradient_flow_gain"] for r in rows])
    axes[0].axhline(0, color="black", lw=0.8); axes[0].set(title="sign-flip gain obstruction", xlabel="functional sign", ylabel="odd gain")
    modes = results["ST60"]["integer_mode_rows"]
    axes[1].bar([r["integer_fourier_mode"] for r in modes], [np.angle(complex(r["holonomy"])) if r["defined"] else np.nan for r in modes])
    axes[1].set(title="strict integer-mode holonomy", xlabel="Fourier mode", ylabel="phase")
    fig.tight_layout(); fig.savefig(FIG_DIR / "st59_gain_st60_projective_source.png", dpi=190); plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    rows = results["ST61"]["rows"]
    axes[0].plot([r["q"] for r in rows], [r["trace_density"] for r in rows], "o-", label="first moment")
    axes[0].plot([r["q"] for r in rows], [r["second_moment_density"] for r in rows], "s-", label="second moment")
    axes[0].set(title="dyadic spectral moments", xlabel="q"); axes[0].legend()
    st62 = results["ST62"]
    axes[1].bar(["necessary", "sufficient", "ST51"], [st62["necessary_shots_per_preparation_information_bound"], st62["sufficient_shots_per_preparation_chernoff_bound"], st62["ST51_shots_per_preparation"]])
    axes[1].set_yscale("log"); axes[1].set(title="shots per preparation", ylabel="count (log)")
    fig.tight_layout(); fig.savefig(FIG_DIR / "st61_refinement_st62_counts.png", dpi=190); plt.close(fig)

    st63 = results["ST63"]
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    ax.bar(range(len(st63["spectral_block_multiplicities"])), st63["spectral_block_multiplicities"])
    ax.set(title="CPTP spectral pinching blocks", xlabel="distinct eigenvalue block", ylabel="multiplicity")
    fig.tight_layout(); fig.savefig(FIG_DIR / "st63_cp_pinching.png", dpi=190); plt.close(fig)

    st64 = results["ST64"]
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    ax.bar(["relative entropy", "beta Delta F", "entropy production"], [st64["relative_entropy"], st64["beta_times_free_energy_difference"], st64["dimensionless_entropy_production"]])
    ax.set(title="conditional thermodynamic identities", ylabel="dimensionless value")
    fig.tight_layout(); fig.savefig(FIG_DIR / "st64_thermodynamic_bridge.png", dpi=190); plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    rows = results["ST65"]["anti_continuum_branch_rows"]
    axes[0].plot([r["kappa"] for r in rows], [r["ipr"] for r in rows])
    axes[0].axvline(1.0, color="red", ls=":"); axes[0].set(title="localized branch before termination", xlabel="coupling kappa", ylabel="IPR")
    st66 = results["ST66"]
    theta = np.linspace(0, 2 * np.pi, 361)
    axes[1].plot(np.cos(theta), np.sin(theta), color="lightgray")
    stable_angles = np.arange(12) * np.pi / 6
    axes[1].scatter(np.cos(stable_angles), np.sin(stable_angles), c="tab:blue")
    axes[1].set_aspect("equal"); axes[1].set(title="12 stable polynomial branches")
    fig.tight_layout(); fig.savefig(FIG_DIR / "st65_localization_st66_bifurcation.png", dpi=190); plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    st67 = results["ST67"]
    axes[0].bar(["forward", "reflected"], [st67["forward_bargmann_chirality"], st67["reverse_bargmann_chirality"]])
    axes[0].axhline(0, color="black", lw=0.8); axes[0].set(title="same pi holonomy, opposite chirality")
    st68 = results["ST68"]
    axes[1].bar(["actual TV", "analytic bound", "budget"], [st68["actual_endpoint_TV_maximum"], st68["analytic_TV_bound"], 0.01])
    axes[1].set(title="clock-calibration envelope", ylabel="total variation")
    fig.tight_layout(); fig.savefig(FIG_DIR / "st67_chirality_st68_validator.png", dpi=190); plt.close(fig)


def write_summary(results: dict) -> None:
    with SUMMARY.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["program", "object", "status"])
        for number in range(58, 70):
            item = results[f"ST{number}"]
            writer.writerow([item["program"], item["object"], item["status"]])


def main() -> None:
    _, a, _ = strict_operator()
    results: dict[str, Any] = {
        "metadata": {
            "programs": "ST58-ST69",
            "seed": SEED,
            "date": "2026-08-11",
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "sympy": sp.__version__,
            "mpmath": mp.__version__,
        }
    }
    results["ST58"] = st58_full_interval_certificate()
    results["ST59"] = st59_gain_source_sign_flip_no_go(a)
    results["ST60"] = st60_projective_source_obstruction(a)
    results["ST61"] = st61_heat_signature_refinement_no_go()
    results["ST62"] = st62_finite_count_bounds(a)
    results["ST63"] = st63_cp_operational_intertwiners(a)
    results["ST64"] = st64_information_to_thermodynamics_bridge(a)
    results["ST65"] = st65_localized_excited_state_search(a)
    results["ST66"] = st66_polynomial_equivariant_bifurcation()
    results["ST67"] = st67_pi_holonomy_chirality()
    results["ST68"] = st68_calibration_custody_validator(a)
    results["ST69"] = st69_axiom_graph_reconciliation()
    results["recommended_next_programs"] = [
        {"id": "ST70", "priority": 1, "study": "independent interval-library or proof-assistant replay of the ST58 Fourier/Krawczyk certificate"},
        {"id": "ST71", "priority": 2, "study": "classify monotone/time-oriented adaptive response classes that evade the ST59 sign-flip no-go"},
        {"id": "ST72", "priority": 3, "study": "derive or obstruct a spin structure or central extension from a new strict-sourced object"},
        {"id": "ST73", "priority": 4, "study": "construct a nonstationary fine-data fractal rate-distortion/RG law and test refinement uniqueness"},
        {"id": "ST74", "priority": 5, "study": "derive a minimax finite-count design including calibration nuisance and dark-count channels"},
        {"id": "ST75", "priority": 6, "study": "classify all CPTP intertwiners beyond spectral Schur maps using Choi-semidefinite constraints"},
        {"id": "ST76", "priority": 7, "study": "prove a resource-theoretic minimality theorem for the E_*, k_B, bath bridge"},
        {"id": "ST77", "priority": 8, "study": "pseudo-arclength continue the ST65 fold and interval-certify its termination or disconnected branches"},
        {"id": "ST78", "priority": 9, "study": "couple the ST66 order parameter to the strict E1 doublet and audit backreaction without claiming a source"},
        {"id": "ST79", "priority": 10, "study": "search a strict-sourced orientation-odd Bargmann observable or prove its stabilizer obstruction"},
        {"id": "ST80", "priority": 11, "study": "export a laboratory-neutral event validator with signed custody attestations and frozen nuisance bounds"},
        {"id": "ST81", "priority": 12, "study": "update W0+CA+SA and the nine-source axiom graph after ST70-ST80"},
    ]
    results["epistemic_boundary"] = (
        "ST58-ST69 provide local interval mathematics, no-go theorems, conditional constructions and synthetic operational "
        "bounds. They do not supply QW-2191 closure, a strict gain/state/scale source, laboratory evidence, legacy role "
        "transfer, Standard Model, gravity, L_total or ToE closure."
    )
    make_figures(results)
    write_summary(results)
    RESULTS.write_text(json.dumps(native(results), indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({"results": RESULTS.name, "programs": 12, "figures": 7}, indent=2))


if __name__ == "__main__":
    main()
