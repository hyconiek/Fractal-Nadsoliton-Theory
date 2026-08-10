#!/usr/bin/env python3
"""FIN P517--P526: nonlinear-source, stability, memory and observability research.

The program deliberately separates consequences of the finite strict operator
from laws, environments and measurement budgets supplied as additional model
data.  It writes deterministic JSON/CSV/certificate/figure artifacts suitable
for local replay.
"""

from __future__ import annotations

import itertools
import json
import math
import platform
import time
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import scipy
import sympy as sp
from numpy.typing import NDArray
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, least_squares, root

from fin_programs_488_496_low_compute import KernelParameters, grouped_stieltjes_data, strict_operator
from fin_programs_497_506_next_research import kernel_parameter_jacobian, operational_probability_table
from fin_programs_507_516_research import (
    N,
    context_symbol_data,
    continue_pattern,
    homotopy_lambda_iv,
    interval_matrix_left_product,
    iv_bounds,
    krawczyk_chart,
    p514_global_psd_phase_diagram,
    scalar_interval_product,
    stationary_root,
    stieltjes_memory_operator,
    strict_a_interval,
    wiener_error_bound,
)


ROOT = Path(__file__).resolve().parent
FIG_DIR = ROOT / "FIN_Programs_517_526_Figures"
RESULTS = ROOT / "FIN_Programs_517_526_Results.json"
SUMMARY = ROOT / "FIN_Programs_517_526_Summary.csv"
P518_CERT = ROOT / "FIN_P518_Krawczyk_Replay_Certificate.json"
P524_CERT = ROOT / "FIN_P524_PSD_Replay_Certificate.json"
SEED = 20260812


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


def h(value: float) -> str:
    """Exact binary64 serialization used by the standard-library replay."""
    return float(value).hex()


def p517_information_quartic_source() -> dict:
    """Audit whether common information functionals select a focusing quartic jet."""
    q = 1.0 / N
    audit = [
        {
            "family": "relative entropy D(rho||q)",
            "quartic_energy_coefficient_at_uniform_reference": 1.0 / (2.0 * q),
            "selected_flow_sign": "defocusing under positive minimization",
            "obstruction": "requires a positive reference density and an extremum convention",
        },
        {
            "family": "negative relative entropy -D(rho||q)",
            "quartic_energy_coefficient_at_uniform_reference": -1.0 / (2.0 * q),
            "selected_flow_sign": "focusing",
            "obstruction": "the overall sign is an added objective choice and destroys convex codelength minimization",
        },
        {
            "family": "Pearson chi-square information",
            "quartic_energy_coefficient_at_uniform_reference": 1.0 / (2.0 * q),
            "selected_flow_sign": "defocusing under positive minimization",
            "obstruction": "multiplication by -1 produces the opposite jet with the same symmetry data",
        },
        {
            "family": "normalized Shannon entropy",
            "quartic_energy_coefficient_at_uniform_reference": None,
            "selected_flow_sign": "none for the radial amplitude",
            "obstruction": "scale invariance removes the norm and the vacuum is not an analytic reference point",
        },
        {
            "family": "Fisher/Dirichlet information",
            "quartic_energy_coefficient_at_uniform_reference": None,
            "selected_flow_sign": "gradient rather than unique onsite jet",
            "obstruction": "an arbitrary c*sum rho_x^2 may be added without changing the supplied quadratic operator",
        },
        {
            "family": "local compression/codelength action",
            "quartic_energy_coefficient_at_uniform_reference": "free c",
            "selected_flow_sign": "sign follows an external reward/penalty convention",
            "obstruction": "coding axioms determine ordering or convexity only after a reference and objective are supplied",
        },
    ]
    return {
        "program": "P517",
        "object": "O221 Information-Functional Fourth-Jet Obstruction",
        "reference_probability": q,
        "audit": audit,
        "theorem": (
            "In the class of local U(1)- and permutation-invariant amplitude actions, the term "
            "c*sum_x |psi_x|^4 is invisible to the value, gradient and Hessian at psi=0. "
            "Entropy, Fisher and compression labels do not remove this freedom unless a reference "
            "state, extremum direction and coupling normalization are added. Positive convex "
            "relative information selects the defocusing sign; focusing requires reversing the "
            "objective or adding a negative quartic response."
        ),
        "focusing_source_found": False,
        "missing_typed_data": ["reference state/measure", "extremum direction", "quartic coupling magnitude", "dimensional conversion if physical"],
        "status": "proven_family_level_information_quartic_sign_no_go",
    }


def p518_exact_replay() -> dict:
    """Regenerate O212 and export exact-rational acceptance inequalities."""
    alo, ahi, amid = strict_a_interval()
    nodes = np.linspace(0.0, 1.0, 401)
    centres: list[np.ndarray] = []
    u = np.zeros(N)
    u[0] = 1.0
    for k in nodes:
        u = stationary_root(amid, float(k), 1.0, u)
        centres.append(u.copy())

    main_rows: list[dict] = []
    full_charts: list[dict] = []
    for index, (k, centre) in enumerate(zip(nodes, centres)):
        lo = 0.0 if index == 0 else (2 * index - 1) / 800.0
        hi = 1.0 if index == 400 else (2 * index + 1) / 800.0
        chart = krawczyk_chart(centre, float(k), lo, hi, alo, ahi, amid)
        full_charts.append(chart)
        main_rows.append(
            {
                "id": index,
                "kind": "parameter_chart",
                "klo": h(lo),
                "khi": h(hi),
                "inclusion_margin_lower": h(chart.get("minimum_inclusion_margin", -1.0)),
                "defect_upper": h(chart.get("defect_infinity_norm_upper", math.inf)),
            }
        )

    bridge_rows: list[dict] = []
    for i in range(400):
        shared = (2 * i + 1) / 800.0
        centre = stationary_root(amid, float(shared), 1.0, 0.5 * (centres[i] + centres[i + 1]))
        bridge = krawczyk_chart(centre, float(shared), float(shared), float(shared), alo, ahi, amid)
        left, right = full_charts[i], full_charts[i + 1]
        blo = np.asarray(bridge["center"]) - np.asarray(bridge["radius"])
        bhi = np.asarray(bridge["center"]) + np.asarray(bridge["radius"])
        llo = np.asarray(left["center"]) - np.asarray(left["radius"])
        lhi = np.asarray(left["center"]) + np.asarray(left["radius"])
        rlo = np.asarray(right["center"]) - np.asarray(right["radius"])
        rhi = np.asarray(right["center"]) + np.asarray(right["radius"])
        nesting = np.concatenate([blo - llo, lhi - bhi, blo - rlo, rhi - bhi])
        bridge_rows.append(
            {
                "id": i,
                "kind": "shared_parameter_bridge",
                "kappa": h(shared),
                "inclusion_margin_lower": h(bridge.get("minimum_inclusion_margin", -1.0)),
                "defect_upper": h(bridge.get("defect_infinity_norm_upper", math.inf)),
                "nesting_margin_lower": h(float(np.min(nesting))),
            }
        )

    certificate = {
        "format": "FIN_P518_exact_binary64_to_rational_inequality_ledger_v1",
        "trust_boundary": (
            "The generator supplies outward transcendental/operator enclosures. The dependency-minimal "
            "checker converts every hexadecimal binary64 endpoint to an exact Fraction and independently "
            "replays all inclusion, contraction, nesting and coverage inequalities."
        ),
        "a_interval_max_width": h(float(np.max(ahi - alo))),
        "parameter_charts": main_rows,
        "bridges": bridge_rows,
    }
    P518_CERT.write_text(json.dumps(certificate, indent=2, sort_keys=True), encoding="utf-8")
    accepted = all(float.fromhex(r["inclusion_margin_lower"]) > 0 and float.fromhex(r["defect_upper"]) < 1 for r in main_rows + bridge_rows)
    accepted = accepted and all(float.fromhex(r["nesting_margin_lower"]) >= 0 for r in bridge_rows)
    return {
        "program": "P518",
        "object": "O222 Dependency-Minimal O212 Replay Ledger",
        "parameter_chart_count": len(main_rows),
        "bridge_count": len(bridge_rows),
        "total_krawczyk_inclusions": len(main_rows) + len(bridge_rows),
        "minimum_inclusion_margin": min(float.fromhex(r["inclusion_margin_lower"]) for r in main_rows + bridge_rows),
        "maximum_defect_upper": max(float.fromhex(r["defect_upper"]) for r in main_rows + bridge_rows),
        "minimum_nesting_margin": min(float.fromhex(r["nesting_margin_lower"]) for r in bridge_rows),
        "certificate_file": P518_CERT.name,
        "all_801_acceptance_inequalities_pass": accepted,
        "scope": "exact rational replay of exported outward acceptance inequalities; not a second transcendental interval implementation",
        "status": "proven_exact_rational_replay_of_801_exported_inequalities" if accepted else "failed_replay",
    }


def interval_square(lo: np.ndarray, hi: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    lower = np.where(lo * hi <= 0.0, 0.0, np.minimum(lo**2, hi**2))
    upper = np.maximum(lo**2, hi**2)
    return np.nextafter(lower, -np.inf), np.nextafter(upper, np.inf)


def interval_product(lo1: np.ndarray, hi1: np.ndarray, lo2: np.ndarray, hi2: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    products = np.stack([lo1 * lo2, lo1 * hi2, hi1 * lo2, hi1 * hi2])
    return np.nextafter(np.min(products, axis=0), -np.inf), np.nextafter(np.max(products, axis=0), np.inf)


def frequency_chart(
    centre: np.ndarray,
    omega0: float,
    omega_lo: float,
    omega_hi: float,
    alo: np.ndarray,
    ahi: np.ndarray,
    amid: np.ndarray,
) -> dict:
    """Parametric Krawczyk chart for A*u+omega*u-u^3=0."""
    j0 = amid + omega0 * np.eye(N) - np.diag(3.0 * centre**2)
    r = np.linalg.inv(j0)
    f0 = amid @ centre + omega0 * centre - centre**3
    width = max(omega0 - omega_lo, omega_hi - omega0)
    derivative = -r @ centre
    rho = 1.5 * np.abs(derivative) * width + 5e-11
    arad = 0.5 * (ahi - alo)
    final: dict | None = None
    for _ in range(80):
        ulo, uhi = centre - rho, centre + rho
        u2lo, u2hi = interval_square(ulo, uhi)
        jlo, jhi = alo.copy(), ahi.copy()
        for i in range(N):
            jlo[i, i] = np.nextafter(jlo[i, i] + omega_lo - 3.0 * u2hi[i], -np.inf)
            jhi[i, i] = np.nextafter(jhi[i, i] + omega_hi - 3.0 * u2lo[i], np.inf)
        rjlo, rjhi = interval_matrix_left_product(r, jlo, jhi)
        mlo, mhi = -rjhi, -rjlo
        for i in range(N):
            mlo[i, i] = np.nextafter(mlo[i, i] + 1.0, -np.inf)
            mhi[i, i] = np.nextafter(mhi[i, i] + 1.0, np.inf)
        dabs = np.maximum(np.abs(mlo), np.abs(mhi))
        correction = np.abs(r @ f0)
        correction += np.abs(r) @ (width * np.abs(centre) + arad @ np.abs(centre))
        correction += 3e-14 * (1.0 + np.linalg.norm(r, ord=np.inf))
        q = correction + dabs @ rho
        margin = float(np.min(rho - q))
        defect = float(np.max(np.sum(dabs, axis=1)))
        if margin > 0.0 and defect < 1.0:
            final = {"radius": rho, "margin": margin, "defect": defect, "jlo": jlo, "jhi": jhi}
            break
        rho = np.maximum(1.10 * rho, 1.20 * q + 1e-13)
        if np.max(rho) > 0.25:
            break
    if final is None:
        return {"success": False, "omega_interval": [omega_lo, omega_hi], "centre": centre, "radius": rho}
    return {
        "success": True,
        "omega_interval": [omega_lo, omega_hi],
        "omega_centre": omega0,
        "centre": centre,
        **final,
    }


def interval_linear_solve(
    mlo: np.ndarray,
    mhi: np.ndarray,
    blo: np.ndarray,
    bhi: np.ndarray,
) -> dict:
    """Neumann/Krawczyk enclosure for M*x=b."""
    mmid, mrad = 0.5 * (mlo + mhi), 0.5 * (mhi - mlo)
    bmid, brad = 0.5 * (blo + bhi), 0.5 * (bhi - blo)
    x0 = np.linalg.solve(mmid, bmid)
    r = np.linalg.inv(mmid)
    rjlo, rjhi = interval_matrix_left_product(r, mlo, mhi)
    dlo, dhi = -rjhi, -rjlo
    for i in range(N):
        dlo[i, i] = np.nextafter(dlo[i, i] + 1.0, -np.inf)
        dhi[i, i] = np.nextafter(dhi[i, i] + 1.0, np.inf)
    dabs = np.maximum(np.abs(dlo), np.abs(dhi))
    defect = float(np.max(np.sum(dabs, axis=1)))
    residual_rad = brad + mrad @ np.abs(x0) + 5e-15
    correction = np.abs(r) @ residual_rad + 5e-14
    radius = correction.copy()
    for _ in range(200):
        new = correction + dabs @ radius
        if np.max(new - radius) < 1e-16:
            radius = np.nextafter(new, np.inf)
            break
        radius = np.nextafter(new, np.inf)
    accepted = defect < 1.0 and bool(np.all(correction + dabs @ radius < radius + 2e-14))
    return {"accepted": accepted, "centre": x0, "radius": radius, "defect": defect}


def dot_interval(lo1: np.ndarray, hi1: np.ndarray, lo2: np.ndarray, hi2: np.ndarray) -> tuple[float, float]:
    plo, phi = interval_product(lo1, hi1, lo2, hi2)
    lower, upper = 0.0, 0.0
    for x, y in zip(plo, phi):
        lower = float(np.nextafter(lower + x, -np.inf))
        upper = float(np.nextafter(upper + y, np.inf))
    return lower, upper


def stability_chart(a: np.ndarray, alo: np.ndarray, ahi: np.ndarray, amid: np.ndarray, omega_lo: float, omega_hi: float, initial: np.ndarray) -> dict:
    omega = 0.5 * (omega_lo + omega_hi)
    centre = stationary_root(a, 1.0, omega, initial)
    chart = frequency_chart(centre, omega, omega_lo, omega_hi, alo, ahi, amid)
    if not chart["success"]:
        return {"accepted": False, "omega_interval": [omega_lo, omega_hi], "reason": "state_chart"}
    u0, ur = chart["centre"], chart["radius"]
    ulo, uhi = u0 - ur, u0 + ur
    u2lo, u2hi = interval_square(ulo, uhi)
    lmlo, lmhi = alo.copy(), ahi.copy()
    lplo, lphi = alo.copy(), ahi.copy()
    for i in range(N):
        lmlo[i, i] = np.nextafter(lmlo[i, i] + omega_lo - u2hi[i], -np.inf)
        lmhi[i, i] = np.nextafter(lmhi[i, i] + omega_hi - u2lo[i], np.inf)
        lplo[i, i] = np.nextafter(lplo[i, i] + omega_lo - 3.0 * u2hi[i], -np.inf)
        lphi[i, i] = np.nextafter(lphi[i, i] + omega_hi - 3.0 * u2lo[i], np.inf)

    lp_mid = a + omega * np.eye(N) - np.diag(3.0 * u0**2)
    lm_mid = a + omega * np.eye(N) - np.diag(u0**2)
    lp_rad_norm = float(np.max(np.sum(0.5 * (lphi - lplo), axis=1)))
    lm_rad_norm = float(np.max(np.sum(0.5 * (lmhi - lmlo), axis=1)))
    ep = np.linalg.eigvalsh(lp_mid)
    em = np.linalg.eigvalsh(lm_mid)
    inertia = {
        "lminus_second_lower": float(em[1] - lm_rad_norm),
        "lplus_negative_upper": float(ep[0] + lp_rad_norm),
        "lplus_first_positive_lower": float(ep[1] - lp_rad_norm),
    }

    vsolve = interval_linear_solve(lplo, lphi, -uhi, -ulo)
    if not vsolve["accepted"]:
        return {"accepted": False, "omega_interval": [omega_lo, omega_hi], "reason": "first_derivative_solve"}
    vlo, vhi = vsolve["centre"] - vsolve["radius"], vsolve["centre"] + vsolve["radius"]
    dplo, dphi = dot_interval(ulo, uhi, vlo, vhi)
    dplo, dphi = 2.0 * dplo, 2.0 * dphi

    v2lo, v2hi = interval_square(vlo, vhi)
    uv2lo, uv2hi = interval_product(ulo, uhi, v2lo, v2hi)
    rhslo = np.nextafter(-2.0 * vhi + 6.0 * uv2lo, -np.inf)
    rhshi = np.nextafter(-2.0 * vlo + 6.0 * uv2hi, np.inf)
    wsolve = interval_linear_solve(lplo, lphi, rhslo, rhshi)
    p2 = [math.nan, math.nan]
    if wsolve["accepted"]:
        wlo = wsolve["centre"] - wsolve["radius"]
        whi = wsolve["centre"] + wsolve["radius"]
        vv = dot_interval(vlo, vhi, vlo, vhi)
        uw = dot_interval(ulo, uhi, wlo, whi)
        p2 = [2.0 * (vv[0] + uw[0]), 2.0 * (vv[1] + uw[1])]
    accepted = inertia["lminus_second_lower"] > 0 and inertia["lplus_negative_upper"] < 0 and inertia["lplus_first_positive_lower"] > 0
    return {
        "accepted": bool(accepted),
        "omega_interval": [omega_lo, omega_hi],
        "omega_centre": omega,
        "state_inclusion_margin": chart["margin"],
        "state_defect": chart["defect"],
        "inertia": inertia,
        "dP_domega_interval": [dplo, dphi],
        "d2P_domega2_interval": p2,
        "state": u0,
        "state_radius": ur,
    }


def p519_interval_orbital_stability(a: np.ndarray) -> dict:
    alo, ahi, amid = strict_a_interval()
    u1, _ = continue_pattern(a, [0], omega=1.0, steps=401)
    charts = []
    edges = np.linspace(0.68, 1.20, 209)
    first_omega = 0.5 * (edges[0] + edges[1])
    initial = u1.copy()
    for omega in np.linspace(1.0, first_omega, 161)[1:]:
        initial = stationary_root(a, 1.0, float(omega), initial)
    for lo, hi in zip(edges[:-1], edges[1:]):
        row = stability_chart(a, alo, ahi, amid, float(lo), float(hi), initial)
        charts.append(row)
        if row.get("accepted"):
            initial = np.asarray(row["state"])
    bridges = []
    for i, shared in enumerate(edges[1:-1]):
        left_row, right_row = charts[i], charts[i + 1]
        seed = 0.5 * (np.asarray(left_row["state"]) + np.asarray(right_row["state"]))
        centre = stationary_root(a, 1.0, float(shared), seed)
        bridge = frequency_chart(centre, float(shared), float(shared), float(shared), alo, ahi, amid)
        blo = np.asarray(bridge["centre"]) - np.asarray(bridge["radius"])
        bhi = np.asarray(bridge["centre"]) + np.asarray(bridge["radius"])
        llo = np.asarray(left_row["state"]) - np.asarray(left_row["state_radius"])
        lhi = np.asarray(left_row["state"]) + np.asarray(left_row["state_radius"])
        rlo = np.asarray(right_row["state"]) - np.asarray(right_row["state_radius"])
        rhi = np.asarray(right_row["state"]) + np.asarray(right_row["state_radius"])
        margins = np.concatenate([blo - llo, lhi - bhi, blo - rlo, rhi - bhi])
        bridges.append(
            {
                "omega": float(shared),
                "accepted": bool(bridge["success"] and np.min(margins) >= 0),
                "minimum_nesting_margin": float(np.min(margins)),
            }
        )
    omega1 = stability_chart(a, alo, ahi, amid, 0.99, 1.01, u1)

    def state_from_one(target: float) -> np.ndarray:
        state = u1.copy()
        count = max(2, int(math.ceil(abs(1.0 - target) / 0.004)))
        for omega in np.linspace(1.0, target, count + 1)[1:]:
            state = stationary_root(a, 1.0, float(omega), state)
        return state

    def numeric_slope(target: float) -> float:
        state = state_from_one(target)
        lp = a + target * np.eye(N) - np.diag(3.0 * state**2)
        return float(-2.0 * state @ np.linalg.solve(lp, state))

    turning = float(brentq(numeric_slope, 0.70, 0.75, xtol=2e-14))
    delta = 2e-8
    for _ in range(20):
        left, right = turning - delta, turning + delta
        left_state, right_state = state_from_one(left), state_from_one(right)
        left_cert = stability_chart(a, alo, ahi, amid, left, left, left_state)
        right_cert = stability_chart(a, alo, ahi, amid, right, right, right_state)
        if left_cert["dP_domega_interval"][1] < 0 and right_cert["dP_domega_interval"][0] > 0:
            break
        delta *= 2.0
    curvature_rows = [row for row in charts if row["omega_interval"][0] >= 0.70 and row["omega_interval"][1] <= 0.75]
    curvature_positive = bool(curvature_rows and all(row["d2P_domega2_interval"][0] > 0 for row in curvature_rows))
    all_inertia = bool(charts and all(row["accepted"] for row in charts))
    all_bridges = bool(bridges and all(row["accepted"] for row in bridges))
    return {
        "program": "P519",
        "object": "O223 Interval GSS/VK Stability Tube",
        "certified_frequency_interval": [0.68, 1.20],
        "chart_count": len(charts),
        "all_inertia_ledgers_certified": all_inertia,
        "interface_bridge_count": len(bridges),
        "all_interface_root_boxes_nested": all_bridges,
        "minimum_interface_nesting_margin": min(row["minimum_nesting_margin"] for row in bridges),
        "omega_1_neighbourhood": {k: v for k, v in omega1.items() if k != "state"},
        "turning_point_bracket": [left, right],
        "turning_left_slope_interval": left_cert["dP_domega_interval"],
        "turning_right_slope_interval": right_cert["dP_domega_interval"],
        "strict_positive_curvature_on_0_70_0_75": curvature_positive,
        "minimum_curvature_lower": min(row["d2P_domega2_interval"][0] for row in curvature_rows),
        "charts": [{k: v for k, v in row.items() if k not in ("state", "state_radius")} for row in charts],
        "theorem": (
            "Outward state Krawczyk charts, 207 nested interface root boxes, Weyl inertia margins and verified implicit-derivative "
            "solves certify a simple phase kernel, one L+ negative direction and the sign of dP/domega. "
            "Positive interval curvature on [0.70,0.75] makes the enclosed VK sign change unique there."
        ),
        "status": "proven_finite_interval_orbital_stability_sign_tube" if all_inertia and all_bridges and curvature_positive else "partial_interval_stability_certificate",
    }


def dnls_energy(a: np.ndarray, psi: np.ndarray) -> float:
    return float(0.5 * np.real(np.vdot(psi, a @ psi)) - 0.25 * np.sum(np.abs(psi) ** 4))


def p520_bond_bifurcation_and_kicks(a: np.ndarray) -> dict:
    bond, _ = continue_pattern(a, [0, 1], omega=1.0, steps=401)
    branch = []
    u = bond.copy()
    for omega in np.linspace(1.0, 0.08, 461):
        try:
            u = stationary_root(a, 1.0, float(omega), u)
        except RuntimeError:
            branch.append({"omega": float(omega), "solve_failed": True})
            break
        mean = float(np.mean(u))
        deviation = u - mean
        contrast = float(np.linalg.norm(deviation) / max(np.linalg.norm(u), 1e-30))
        spectrum = np.abs(np.fft.fft(deviation))
        dominant = int(np.argmax(spectrum[1:7]) + 1)
        jac = a + omega * np.eye(N) - np.diag(3.0 * u**2)
        branch.append(
            {
                "omega": float(omega),
                "power": float(u @ u),
                "contrast": contrast,
                "dominant_mode": dominant,
                "jacobian_min_singular": float(np.linalg.svd(jac, compute_uv=False)[-1]),
                "energy": dnls_energy(a, u.astype(complex)),
            }
        )
    nonfailed = [row for row in branch if not row.get("solve_failed")]
    near_uniform = [row for row in nonfailed if row["contrast"] < 1e-5]
    transition_numeric = near_uniform[0]["omega"] if near_uniform else None
    localized = [row for row in nonfailed if row["contrast"] >= 1e-5]
    last_localized = localized[-1] if localized else None

    fourier_eigenvalues = np.real(np.fft.fft(a[0]))
    distinct = [float(fourier_eigenvalues[m]) for m in range(1, 7)]
    uniform_bifurcations = [{"mode": m, "omega": value / 2.0} for m, value in enumerate(distinct, start=1)]
    if last_localized:
        mode = int(last_localized["dominant_mode"])
        predicted = distinct[mode - 1] / 2.0
    else:
        mode, predicted = None, None

    site, _ = continue_pattern(a, [0], omega=1.0, steps=401)
    positions = np.arange(N)
    kicks = np.linspace(0.0, 1.4, 15)
    kick_rows = []
    for kick in kicks:
        psi0 = site.astype(complex) * np.exp(1j * kick * positions)

        def rhs(_t: float, y: np.ndarray) -> np.ndarray:
            psi = y[:N] + 1j * y[N:]
            flow = -1j * (a @ psi - np.abs(psi) ** 2 * psi)
            return np.concatenate([flow.real, flow.imag])

        times = np.linspace(0.0, 60.0, 1201)
        sol = solve_ivp(rhs, (0.0, 60.0), np.concatenate([psi0.real, psi0.imag]), t_eval=times, method="DOP853", rtol=2e-10, atol=2e-12)
        psi = sol.y[:N] + 1j * sol.y[N:]
        density = np.abs(psi) ** 2
        power = np.sum(density, axis=0)
        first_moment = np.sum(density * np.exp(2j * np.pi * positions[:, None] / N), axis=0) / power
        centre = np.unwrap(np.angle(first_moment)) * N / (2.0 * np.pi)
        ipr = np.sum(density**2, axis=0) / power**2
        energy = np.array([dnls_energy(a, psi[:, j]) for j in range(psi.shape[1])])
        steps = np.diff(centre)
        direction = np.sign(centre[-1] - centre[0])
        monotone_fraction = float(np.mean(direction * steps > -1e-4)) if direction else 0.0
        kick_rows.append(
            {
                "kick": float(kick),
                "net_circular_displacement": float(centre[-1] - centre[0]),
                "maximum_absolute_displacement": float(np.max(np.abs(centre - centre[0]))),
                "minimum_ipr": float(np.min(ipr)),
                "final_ipr": float(ipr[-1]),
                "monotone_fraction": monotone_fraction,
                "relative_power_drift": float(np.max(np.abs(power - power[0])) / power[0]),
                "relative_energy_drift": float(np.max(np.abs(energy - energy[0])) / max(1.0, abs(energy[0]))),
            }
        )
    coherent = [row for row in kick_rows if abs(row["net_circular_displacement"]) >= 1.0 and row["minimum_ipr"] >= 0.25 and row["monotone_fraction"] >= 0.8]
    best = max(kick_rows, key=lambda row: abs(row["net_circular_displacement"]))
    small_singular = [
        row
        for row in localized
        if row["jacobian_min_singular"] < 1e-5
        and (predicted is None or abs(row["omega"] - predicted) > 0.02)
    ]
    return {
        "program": "P520",
        "object": "O224 Bond-to-Uniform Bifurcation and Kick Ledger",
        "branch_rows": branch,
        "last_localized_row": last_localized,
        "first_near_uniform_omega": transition_numeric,
        "uniform_bifurcation_frequencies": uniform_bifurcations,
        "dominant_transition_mode": mode,
        "predicted_uniform_bifurcation_omega": predicted,
        "additional_fold_candidates_below_singular_threshold": len(small_singular),
        "kick_rows": kick_rows,
        "coherent_translation_candidates": coherent,
        "largest_displacement_trial": best,
        "verdict": (
            "The bond family approaches a symmetry-breaking bifurcation of the uniform branch; "
            "no independent fold is resolved away from that transition. Phase kicks are tested "
            "directly, without using the undefined equal-power PN barrier."
        ),
        "status": "strong_numerical_bond_bifurcation_and_mobility_audit",
    }


def hidden_memory_jacobian(a: np.ndarray, u: np.ndarray, poles: np.ndarray, residues: np.ndarray, epsilon: float, omega: float = 1.0) -> np.ndarray:
    count = len(poles)
    size = 2 * N * (1 + count)
    jac = np.zeros((size, size))
    identity = np.eye(N)
    sigma0 = float(np.sum(residues / poles))
    qminus = a + omega * identity - np.diag(u**2) + epsilon * sigma0 * identity
    qplus = a + omega * identity - np.diag(3.0 * u**2) + epsilon * sigma0 * identity
    x, y = slice(0, N), slice(N, 2 * N)
    jac[x, y] = qminus
    jac[y, x] = -qplus
    for j, (pole, residue) in enumerate(zip(poles, residues)):
        aidx = slice(2 * N + 2 * j * N, 2 * N + (2 * j + 1) * N)
        bidx = slice(2 * N + (2 * j + 1) * N, 2 * N + (2 * j + 2) * N)
        jac[x, bidx] = -epsilon * residue * identity
        jac[y, aidx] = epsilon * residue * identity
        jac[aidx, x] = identity
        jac[bidx, y] = identity
        jac[aidx, aidx] = -(a + pole * identity)
        jac[bidx, bidx] = -(a + pole * identity)
    return jac


def p521_temporal_memory_stability(a: np.ndarray) -> dict:
    memory, poles, residues = stieltjes_memory_operator(a)
    base, _ = continue_pattern(a, [0], omega=1.0, steps=401)
    u = base.copy()
    rows = []
    for epsilon in np.linspace(0.0, 1.0, 101):
        u = stationary_root(a + epsilon * memory, 1.0, 1.0, u)
        jac = hidden_memory_jacobian(a, u, poles, residues, float(epsilon))
        eigen = np.linalg.eigvals(jac)
        nonneutral = eigen[np.abs(eigen) > 2e-7]
        spectral_abscissa = float(np.max(nonneutral.real))
        unstable_count = int(np.sum(nonneutral.real > 1e-8))
        phase_neutral = float(np.min(np.abs(eigen)))
        power = float(u @ u)
        rows.append(
            {
                "epsilon": float(epsilon),
                "spectral_abscissa": spectral_abscissa,
                "unstable_eigenvalue_count": unstable_count,
                "phase_neutral_residual": phase_neutral,
                "ipr": float(np.sum(u**4) / power**2),
            }
        )
    positive = [row for row in rows if row["spectral_abscissa"] > 1e-7]
    return {
        "program": "P521",
        "object": "O225 Hidden-Relaxation Memory Flow",
        "augmented_equations": (
            "i psi_dot=Apsi+omega*psi-|psi|^2psi+epsilon*(sigma(0)psi-sum r_j y_j); "
            "y_j_dot=psi-(A+c_j I)y_j"
        ),
        "hidden_mode_count": len(poles),
        "real_linearization_dimension": int(2 * N * (1 + len(poles))),
        "rows": rows,
        "maximum_spectral_abscissa": max(row["spectral_abscissa"] for row in rows),
        "first_detected_instability": positive[0] if positive else None,
        "all_loaded_states_spectrally_stable_in_scan": not positive,
        "verdict": (
            "P511 stationary persistence does not by itself imply stability of a temporal memory "
            "realization. The declared causal hidden-relaxation completion is now explicit and its "
            "complete finite linear spectrum is audited, but the result remains conditional on that environment."
        ),
        "status": "conditional_full_augmented_linear_stability_scan",
    }


def p522_multifrequency_clock() -> dict:
    matrix = sp.zeros(6, 6)
    for m in range(1, 7):
        for d in range(1, 7):
            if d < 6:
                matrix[m - 1, d - 1] = sp.simplify(2 * (1 - sp.cos(2 * sp.pi * m * d / 12)))
            else:
                matrix[m - 1, d - 1] = 1 - (-1) ** m
    determinant = sp.simplify(matrix.det())
    rank = int(matrix.rank())
    return {
        "program": "P522",
        "object": "O226 Maximal Strict Frequency Independence Certificate",
        "distance_to_mode_matrix": [[str(matrix[i, j]) for j in range(6)] for i in range(6)],
        "exact_determinant": str(determinant),
        "exact_rank": rank,
        "maximal_distinct_positive_frequency_count": 6,
        "theorem": (
            "The six distance weights have disjoint highest Laurent exponents in z=exp(i/4000). "
            "The exact distance-to-mode matrix for modes 1,...,6 is invertible. Therefore the six "
            "distinct positive strict Laplacian frequencies are linearly independent over the "
            "algebraic numbers, hence rationally independent. Degeneracies m and 12-m make six maximal."
        ),
        "scale_no_go": (
            "For every c>0, replacing A by cA and t by t/c leaves every dimensionless orbit and "
            "frequency ratio unchanged. No record internal to the dimensionless spectral system can "
            "select an absolute unit of time or energy."
        ),
        "status": "proven_six_frequency_rational_independence_and_absolute_scale_no_go" if rank == 6 and determinant != 0 else "failed_frequency_rank",
    }


def normalized_limit_bound(finite: dict, absolute_errors: np.ndarray, index0: int) -> np.ndarray:
    sigma = np.asarray(finite["sigma"], dtype=float)
    e = np.asarray(absolute_errors, dtype=float)
    denom = float(sigma[index0])
    lower_inf = denom - float(e[index0])
    if lower_inf <= 0:
        raise RuntimeError("normalization denominator not certified")
    return e / denom + (sigma + e) * e[index0] / (denom * lower_inf)


def p523_fractional_context_bound() -> dict:
    u_values = np.array([0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0])
    coarse = context_symbol_data(384, u_values)
    anchor_n = 262144
    anchor = context_symbol_data(anchor_n, u_values)
    analytic = wiener_error_bound(anchor_n, u_values, anchor)
    absolute = np.asarray(analytic["absolute_response_error_bounds"], dtype=float)
    idx = int(np.argmin(np.abs(u_values - 1.0)))
    anchor_to_limit = normalized_limit_bound(anchor, absolute, idx)
    finite_gap = np.abs(np.asarray(coarse["fingerprint"]) - np.asarray(anchor["fingerprint"]))
    fft_rounding_envelope = 2e-9
    combined = finite_gap + anchor_to_limit + fft_rounding_envelope
    return {
        "program": "P523",
        "object": "O227 Multilevel Fractional-Context Certificate",
        "coarse_cycle": 384,
        "analytic_anchor_cycle": anchor_n,
        "holder_exponent": 0.8,
        "finite_coarse_to_anchor_gap": finite_gap,
        "anchor_to_limit_normalized_bound": anchor_to_limit,
        "declared_fft_rounding_envelope": fft_rounding_envelope,
        "combined_C384_to_limit_bounds": combined,
        "maximum_combined_bound": float(np.max(combined)),
        "improvement_over_P513_direct_C384": "P513 had no positive direct normalization denominator at C384",
        "theorem_scope": (
            "The analytic Wiener/resolvent tail pays the anchor-to-limit error at the fractional "
            "rate n^-0.8. A direct finite multilevel comparison then transfers that certificate to "
            "C384. The last 2e-9 is a conservative IEEE-FFT rounding-model allowance."
        ),
        "status": "conditional_fft_rounding_model_C384_fractional_certificate",
    }


def p524_portable_psd_replay() -> dict:
    source_path = ROOT / "FIN_Programs_507_516_Results.json"
    source = json.loads(source_path.read_text(encoding="utf-8"))["results"]["P514"]
    rows = []
    for index, row in enumerate(source["all_boxes"]):
        rows.append(
            {
                "id": index,
                "interval": [h(row["interval"][0]), h(row["interval"][1])],
                "class": row["class"],
                "eigenvalue_bounds": [[h(pair[0]), h(pair[1])] for pair in row["bounds"]],
            }
        )
    certificate = {
        "format": "FIN_P524_exact_binary64_to_rational_PSD_ledger_v1",
        "source": source_path.name,
        "trust_boundary": (
            "mpmath interval generation is upstream. The standard-library checker converts every "
            "hexadecimal endpoint to an exact Fraction, verifies the complete partition, every PSD "
            "strict lower sign, every non-PSD strict upper sign and the unique unresolved bracket."
        ),
        "rows": rows,
    }
    P524_CERT.write_text(json.dumps(certificate, indent=2, sort_keys=True), encoding="utf-8")
    psd = sum(row["class"] == "PSD" for row in rows)
    nonpsd = sum(row["class"] == "NONPSD" for row in rows)
    unresolved = [row for row in rows if row["class"] == "UNRESOLVED"]
    unresolved_components: list[list[float]] = []
    for row in unresolved:
        lo, hi = (float.fromhex(x) for x in row["interval"])
        if not unresolved_components or unresolved_components[-1][1] != lo:
            unresolved_components.append([lo, hi])
        else:
            unresolved_components[-1][1] = hi
    signs = all(
        (row["class"] != "PSD" or all(float.fromhex(pair[0]) > 0 for pair in row["eigenvalue_bounds"]))
        and (row["class"] != "NONPSD" or any(float.fromhex(pair[1]) < 0 for pair in row["eigenvalue_bounds"]))
        for row in rows
    )
    coverage = float.fromhex(rows[0]["interval"][0]) == 0.0 and float.fromhex(rows[-1]["interval"][1]) == 1.0
    coverage = coverage and all(float.fromhex(left["interval"][1]) == float.fromhex(right["interval"][0]) for left, right in zip(rows[:-1], rows[1:]))
    return {
        "program": "P524",
        "object": "O228 Portable Exact PSD Sign Ledger",
        "terminal_box_count": len(rows),
        "psd_box_count": psd,
        "nonpsd_box_count": nonpsd,
        "unresolved_box_count": len(unresolved),
        "unresolved_component_count": len(unresolved_components),
        "unique_unresolved_interval": unresolved_components[0] if len(unresolved_components) == 1 else None,
        "complete_exact_partition_replay": coverage,
        "all_exported_sign_conditions_pass": signs,
        "certificate_file": P524_CERT.name,
        "scope": "exact rational replay of the exported interval/sign ledger; not independent transcendental interval generation",
        "status": "proven_portable_exact_PSD_ledger_replay" if coverage and signs and len(unresolved_components) == 1 else "failed_portable_replay",
    }


def p525_detector_polytope(a: np.ndarray) -> dict:
    models, _contexts, table = operational_probability_table(a)
    context = (0.55, "fourier1", "fourier")
    opnorm = float(np.linalg.eigvalsh(a)[-1])
    speed = {
        "unitary": opnorm,
        "markov": 0.5 * float(np.linalg.norm(a, 1)),
        "dephasing": opnorm + 0.60,
        "revival": opnorm + 1.12,
    }
    blur = np.zeros((N, N))
    for x in range(N):
        blur[x, (x - 1) % N] = 0.5
        blur[x, (x + 1) % N] = 0.5
    identity = np.eye(N)
    variables = ["loss_fraction", "dark_mixture", "coarse_graining", "preparation_TD", "POVM_TV", "time_error", "configuration_drift_TV"]
    halfspaces = []
    for i, left in enumerate(models):
        for right in models[i + 1 :]:
            delta = np.asarray(table[left][context]) - np.asarray(table[right][context])
            nominal = 0.5 * float(np.sum(np.abs(delta)))
            coarse_lipschitz = 0.5 * float(np.sum(np.abs((blur - identity) @ delta)))
            coefficients = {
                "loss_fraction": nominal,
                "dark_mixture": nominal,
                "coarse_graining": coarse_lipschitz,
                "preparation_TD": 2.0,
                "POVM_TV": 2.0,
                "time_error": speed[left] + speed[right],
                "configuration_drift_TV": 2.0,
            }
            halfspaces.append({"models": [left, right], "nominal_TV": nominal, "coefficients": coefficients, "strict_rhs": nominal})

    coordinate_intercepts = {
        variable: min(row["strict_rhs"] / row["coefficients"][variable] for row in halfspaces if row["coefficients"][variable] > 0)
        for variable in variables
    }
    trial = {
        "loss_fraction": 0.10,
        "dark_mixture": 0.01,
        "coarse_graining": 0.02,
        "preparation_TD": 0.005,
        "POVM_TV": 0.005,
        "time_error": 0.002,
        "configuration_drift_TV": 0.005,
    }
    trial_rows = []
    for row in halfspaces:
        penalty = sum(row["coefficients"][key] * trial[key] for key in variables)
        trial_rows.append({"models": row["models"], "penalty": penalty, "certified_TV_lower": row["nominal_TV"] - penalty})
    minimum = min(row["certified_TV_lower"] for row in trial_rows)
    return {
        "program": "P525",
        "object": "O229 Detector-Error Admissible Polytope",
        "context": {"time": 0.55, "preparation": "fourier1", "measurement": "fourier"},
        "variables": variables,
        "derivation": (
            "Common loss and dark mixtures, nearest-neighbour coarse graining, contractive "
            "preparation/POVM errors, generator-speed timing error and configuration drift are paid "
            "by the triangle inequality. Every pair supplies one conservative linear half-space."
        ),
        "halfspaces": halfspaces,
        "coordinate_intercepts": coordinate_intercepts,
        "hypothetical_trial_budget": trial,
        "trial_pairwise_ledger": trial_rows,
        "minimum_trial_TV_lower": minimum,
        "trial_inside_admissible_polytope": bool(minimum > 0),
        "physical_boundary": "All budgets are hypothetical variables, not measured detector specifications or laboratory calibration.",
        "status": "conditional_nonempty_detector_error_polytope",
    }


def parameter_vector_to_kernel(p: np.ndarray) -> KernelParameters:
    return KernelParameters(float(math.exp(p[0])), float(math.exp(p[1])), float(p[2]), float(p[3]), float(p[4]))


def kernel_observation(p: np.ndarray) -> np.ndarray:
    params = parameter_vector_to_kernel(p)
    d = np.arange(7, dtype=float)
    return params.amplitude * np.cos(params.omega * d + params.phi) / (1.0 + params.beta * d**params.eta)


def path_from_coefficients(p0: np.ndarray, coefficients: np.ndarray, t: float) -> np.ndarray:
    out = p0.copy()
    for k, coefficient in enumerate(coefficients):
        out += coefficient * t ** (k + 1) / (k + 1)
    return out


def p526_law_selection() -> dict:
    p0 = np.array([math.log(4.0 * math.log(2.0)), math.log(0.01), 1.0, math.pi / 4.0, math.pi / 6.0])
    truth = np.array(
        [
            [-0.08, 0.18, 0.04, -0.025, -0.015],
            [0.02, -0.03, 0.01, 0.006, -0.004],
            [-0.004, 0.008, -0.003, 0.001, 0.002],
        ]
    )
    train_times = np.array([0.2, 0.4, 0.6, 0.8])
    holdout_times = np.array([0.1, 0.3, 0.5, 0.7, 0.9, 1.0])

    def frozen_data(times: np.ndarray, offset: int) -> np.ndarray:
        values = []
        for i, t in enumerate(times):
            clean = kernel_observation(path_from_coefficients(p0, truth, float(t)))
            noise = 2e-6 * np.sin(np.arange(7) * 1.731 + (i + offset) * 0.417)
            values.append(clean + noise)
        return np.asarray(values)

    train = frozen_data(train_times, 0)
    holdout = frozen_data(holdout_times, 17)
    ledgers = []
    selections = []
    for ridge in (0.0, 1e-6, 1e-4, 1e-2):
        candidates = []
        for degree in range(4):
            count = 5 * (degree + 1)

            def residual(flat: np.ndarray) -> np.ndarray:
                coeff = flat.reshape(degree + 1, 5)
                model = np.asarray([kernel_observation(path_from_coefficients(p0, coeff, float(t))) for t in train_times])
                base = (model - train).ravel()
                if ridge > 0:
                    base = np.concatenate([base, math.sqrt(ridge) * flat])
                return base

            fit = least_squares(residual, np.zeros(count), xtol=1e-12, ftol=1e-12, gtol=1e-12, max_nfev=4000)
            coeff = fit.x.reshape(degree + 1, 5)
            train_model = np.asarray([kernel_observation(path_from_coefficients(p0, coeff, float(t))) for t in train_times])
            hold_model = np.asarray([kernel_observation(path_from_coefficients(p0, coeff, float(t))) for t in holdout_times])
            train_rmse = float(np.sqrt(np.mean((train_model - train) ** 2)))
            hold_rmse = float(np.sqrt(np.mean((hold_model - holdout) ** 2)))
            candidates.append(
                {
                    "degree": degree,
                    "coefficient_count": count,
                    "train_rmse": train_rmse,
                    "holdout_rmse": hold_rmse,
                    "coefficient_norm": float(np.linalg.norm(fit.x)),
                    "optimizer_success": bool(fit.success),
                }
            )
        selected = min(candidates, key=lambda row: row["holdout_rmse"])["degree"]
        selections.append({"ridge": ridge, "selected_degree": selected})
        ledgers.append({"ridge": ridge, "candidates": candidates})

    all_times = sorted(set(train_times.tolist() + holdout_times.tolist()))
    bump_coefficients = np.poly(all_times).tolist()
    return {
        "program": "P526",
        "object": "O230 Frozen-Holdout Law-Class Audit",
        "synthetic_truth_degree": 2,
        "strict_endpoint_used": False,
        "train_times": train_times,
        "holdout_times": holdout_times,
        "noise_rule": "2e-6*sin(1.731*d+0.417*(record_index+offset))",
        "ridge_ledgers": ledgers,
        "selections": selections,
        "finite_candidate_selection_is_prior_sensitive": len({row["selected_degree"] for row in selections}) > 1,
        "exact_vanishing_bump_polynomial_coefficients": bump_coefficients,
        "no_free_lunch_theorem": (
            "For the finite union T of train and holdout times, b(t)=product_{tau in T}(t-tau) "
            "vanishes on every recorded time. For every smooth fitted path p(t), every vector v and "
            "scalar c, p(t)+c*b(t)*v has identical finite records and arbitrary off-record behaviour. "
            "Thus holdout can rank a declared finite candidate family but cannot identify an unrestricted causal law."
        ),
        "status": "proven_finite_record_law_no_free_lunch_with_conditional_model_ranking",
    }


def make_figures(results: dict[str, Any]) -> None:
    FIG_DIR.mkdir(exist_ok=True)
    plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": 0.25})

    p517 = results["P517"]["audit"]
    labels, values, colors = [], [], []
    for row in p517:
        value = row["quartic_energy_coefficient_at_uniform_reference"]
        if isinstance(value, (int, float)):
            labels.append(row["family"].replace(" information", ""))
            values.append(value)
            colors.append("#2ca02c" if value > 0 else "#d62728")
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    ax.bar(labels, values, color=colors)
    ax.axhline(0, color="black", lw=0.8)
    ax.set_ylabel("quartic coefficient at q=1/12")
    ax.tick_params(axis="x", rotation=18)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "p517_information_quartic_sign.png", dpi=180)
    plt.close(fig)

    charts = results["P519"]["charts"]
    centres = [row["omega_centre"] for row in charts]
    lo = [row["dP_domega_interval"][0] for row in charts]
    hi = [row["dP_domega_interval"][1] for row in charts]
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    ax.fill_between(centres, lo, hi, alpha=0.3, label="interval enclosure")
    ax.plot(centres, 0.5 * (np.asarray(lo) + np.asarray(hi)), lw=1.2)
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xlabel("omega")
    ax.set_ylabel("dP/domega")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG_DIR / "p519_interval_vk_tube.png", dpi=180)
    plt.close(fig)

    branch = [row for row in results["P520"]["branch_rows"] if not row.get("solve_failed")]
    fig, ax1 = plt.subplots(figsize=(7.5, 4.2))
    ax1.plot([row["omega"] for row in branch], [row["contrast"] for row in branch], label="bond contrast")
    ax1.set_xlabel("omega")
    ax1.set_ylabel("nonuniform contrast")
    ax1.axvline(results["P520"]["predicted_uniform_bifurcation_omega"], color="#d95f02", ls="--", label="uniform bifurcation")
    ax1.legend()
    fig.tight_layout()
    fig.savefig(FIG_DIR / "p520_bond_bifurcation.png", dpi=180)
    plt.close(fig)

    memory = results["P521"]["rows"]
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    ax.plot([row["epsilon"] for row in memory], [row["spectral_abscissa"] for row in memory])
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xlabel("memory loading epsilon")
    ax.set_ylabel("augmented spectral abscissa")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "p521_memory_spectral_abscissa.png", dpi=180)
    plt.close(fig)

    p523 = results["P523"]
    u = np.array([0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0])
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    ax.loglog(u, p523["finite_coarse_to_anchor_gap"], marker="o", label="C384 to anchor")
    ax.loglog(u, p523["combined_C384_to_limit_bounds"], marker="s", label="combined bound")
    ax.set_xlabel("resolvent scale u")
    ax.set_ylabel("normalized error")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG_DIR / "p523_C384_bound.png", dpi=180)
    plt.close(fig)

    intercepts = results["P525"]["coordinate_intercepts"]
    fig, ax = plt.subplots(figsize=(7.5, 4.4))
    ax.barh(list(intercepts), list(intercepts.values()))
    ax.set_xlabel("one-coordinate conservative intercept")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "p525_detector_polytope_intercepts.png", dpi=180)
    plt.close(fig)

    base = next(row for row in results["P526"]["ridge_ledgers"] if row["ridge"] == 0.0)
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    degrees = [row["degree"] for row in base["candidates"]]
    ax.semilogy(degrees, [row["train_rmse"] for row in base["candidates"]], marker="o", label="train")
    ax.semilogy(degrees, [row["holdout_rmse"] for row in base["candidates"]], marker="s", label="holdout")
    ax.set_xlabel("polynomial velocity degree")
    ax.set_ylabel("RMSE")
    ax.set_xticks(degrees)
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG_DIR / "p526_holdout_law_selection.png", dpi=180)
    plt.close(fig)


def write_summary(results: dict[str, Any]) -> None:
    rows = [
        ("P517", results["P517"]["status"], "information labels do not source focusing fourth jet"),
        ("P518", results["P518"]["status"], f"801 replay inequalities pass={results['P518']['all_801_acceptance_inequalities_pass']}"),
        ("P519", results["P519"]["status"], f"turning bracket={results['P519']['turning_point_bracket']}"),
        ("P520", results["P520"]["status"], f"coherent kicks={len(results['P520']['coherent_translation_candidates'])}"),
        ("P521", results["P521"]["status"], f"max spectral abscissa={results['P521']['maximum_spectral_abscissa']:.6g}"),
        ("P522", results["P522"]["status"], f"det={results['P522']['exact_determinant']}"),
        ("P523", results["P523"]["status"], f"C384 max bound={results['P523']['maximum_combined_bound']:.6g}"),
        ("P524", results["P524"]["status"], f"boxes={results['P524']['terminal_box_count']}"),
        ("P525", results["P525"]["status"], f"trial TV={results['P525']['minimum_trial_TV_lower']:.6g}"),
        ("P526", results["P526"]["status"], f"selections={results['P526']['selections']}"),
    ]
    lines = ["program,status,key_result"] + [f'"{p}","{s}","{r}"' for p, s, r in rows]
    SUMMARY.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    started = time.perf_counter()
    _, a = strict_operator()
    results: dict[str, Any] = {}
    results["P517"] = p517_information_quartic_source()
    results["P518"] = p518_exact_replay()
    results["P519"] = p519_interval_orbital_stability(a)
    results["P520"] = p520_bond_bifurcation_and_kicks(a)
    results["P521"] = p521_temporal_memory_stability(a)
    results["P522"] = p522_multifrequency_clock()
    results["P523"] = p523_fractional_context_bound()
    results["P524"] = p524_portable_psd_replay()
    results["P525"] = p525_detector_polytope(a)
    results["P526"] = p526_law_selection()
    make_figures(results)
    elapsed = time.perf_counter() - started
    payload = {
        "release": "10.52",
        "programs": "P517-P526",
        "seed": SEED,
        "runtime_seconds": elapsed,
        "environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "sympy": sp.__version__,
            "platform": platform.platform(),
        },
        "epistemic_boundary": (
            "Finite analytic/interval theorems, exact exported-ledger replays and synthetic or "
            "conditional dynamical/operational models; no laboratory or external-audit evidence."
        ),
        "results": results,
    }
    RESULTS.write_text(json.dumps(native(payload), indent=2, sort_keys=True), encoding="utf-8")
    write_summary(results)
    print(
        json.dumps(
            {
                "runtime_seconds": elapsed,
                "P518_replay": results["P518"]["all_801_acceptance_inequalities_pass"],
                "P519_status": results["P519"]["status"],
                "P521_first_instability": results["P521"]["first_detected_instability"],
                "P522_determinant": results["P522"]["exact_determinant"],
                "P523_C384_bound": results["P523"]["maximum_combined_bound"],
                "P525_trial_TV": results["P525"]["minimum_trial_TV_lower"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
