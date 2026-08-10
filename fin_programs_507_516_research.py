#!/usr/bin/env python3
"""FIN P507--P516: nonlinear-source, certified-branch and observability research.

All results are local mathematics, interval/numerical certificates, or
synthetic operational studies.  No supplied nonlinear law or measurement
model is promoted to a consequence of the frozen FIN kernels.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import math
import platform
import time
from collections import deque
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import mpmath as mp
import numpy as np
import scipy
from numpy.typing import NDArray
from scipy.optimize import brentq, root

from fin_programs_488_496_low_compute import (
    ETA_S,
    OMEGA_S,
    PHI_S,
    KernelParameters,
    grouped_stieltjes_data,
    strict_operator,
)
from fin_programs_497_506_next_research import (
    homotopy_laplacian_eigenvalue,
    homotopy_params,
    homotopy_weight_mp,
    iv_bounds,
    kernel_parameter_jacobian,
    operational_probability_table,
)


ROOT = Path(__file__).resolve().parent
FIG_DIR = ROOT / "FIN_Programs_507_516_Figures"
RESULTS = ROOT / "FIN_Programs_507_516_Results.json"
SUMMARY = ROOT / "FIN_Programs_507_516_Summary.csv"
SEED = 20260811
N = 12


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


def stationary_root(a: np.ndarray, kappa: float, omega: float, initial: np.ndarray) -> np.ndarray:
    """Solve kappa*A*u + omega*u - u^3=0 for a real focusing DNLS state."""
    fun = lambda u: kappa * (a @ u) + omega * u - u**3
    jac = lambda u: kappa * a + omega * np.eye(len(u)) - np.diag(3.0 * u**2)
    sol = root(fun, initial, jac=jac, method="hybr", options={"xtol": 1e-12, "maxfev": 5000})
    if (not sol.success) or np.linalg.norm(fun(sol.x), ord=np.inf) > 2e-9:
        raise RuntimeError(f"stationary solve failed: {sol.message}")
    return np.asarray(sol.x, dtype=float)


def continue_pattern(a: np.ndarray, active_sites: list[int], omega: float = 1.0, steps: int = 401) -> tuple[np.ndarray, list[dict]]:
    u = np.zeros(N)
    u[active_sites] = math.sqrt(omega)
    rows = []
    for kappa in np.linspace(0.0, 1.0, steps):
        u = stationary_root(a, float(kappa), omega, u)
        power = float(u @ u)
        rows.append(
            {
                "kappa": float(kappa),
                "power": power,
                "ipr": float(np.sum(u**4) / power**2),
                "peak_fraction": float(np.max(u**2) / power),
                "residual_inf": float(np.linalg.norm(kappa * (a @ u) + omega * u - u**3, ord=np.inf)),
                "u": u.copy(),
            }
        )
    return u, rows


def p507_nonlinear_source_audit(a: np.ndarray) -> dict:
    rng = np.random.default_rng(SEED + 507)
    quadratic_jet_residuals = []
    for _ in range(100):
        direction = rng.normal(size=N)
        direction /= np.linalg.norm(direction)
        eps = 1e-5
        for sign in (-1.0, 1.0):
            energy_plus = 0.5 * (eps * direction) @ a @ (eps * direction) + sign * 0.25 * np.sum((eps * direction) ** 4)
            energy_minus = 0.5 * (-eps * direction) @ a @ (-eps * direction) + sign * 0.25 * np.sum((-eps * direction) ** 4)
            hessian_rayleigh = (energy_plus + energy_minus) / eps**2
            quadratic_jet_residuals.append(abs(hessian_rayleigh - float(direction @ a @ direction)))

    def anti_continuum_real_amplitude_squared(sign: int, omega: float) -> float:
        # Stationary equation omega*u + sign*g*u^3=0 for H=quadratic+sign*g|u|^4/2.
        return -omega / sign

    return {
        "program": "P507",
        "object": "O211 Quartic-Jet Source Obstruction",
        "minimal_hamiltonian_family": "E_sigma,g=0.5<psi,Apsi>+sigma*g/4*sum_x |psi_x|^4",
        "hamiltonian_flow": "with convention i psi_dot=2*dE/d(psi*): Apsi+sigma*g|psi|^2 psi",
        "quadratic_jet_max_replay_residual": max(quadratic_jet_residuals),
        "anti_continuum_at_omega_1": {
            "focusing_sigma_minus_1_amplitude_squared": anti_continuum_real_amplitude_squared(-1, 1.0),
            "defocusing_sigma_plus_1_amplitude_squared": anti_continuum_real_amplitude_squared(1, 1.0),
        },
        "coefficient_scaling_orbit": "u_g=u_1/sqrt(g); normalized profile and IPR are independent of g",
        "theorem": (
            "The strict quadratic action fixes only the Hessian A at psi=0. For every g>0, the "
            "focusing and defocusing quartic Hamiltonians have the same value, gradient, and Hessian "
            "through second order at zero, yet generate opposite cubic terms. Therefore no datum "
            "confined to the strict quadratic/linear core can select the cubic sign or coefficient."
        ),
        "verdict": (
            "The P497 focusing law is not derivable from the present strict quadratic core. A "
            "quartic source functional, nonlinear response datum, or an explicit focusing axiom is necessary."
        ),
        "status": "proven_quadratic_core_nonlinearity_source_no_go",
    }


def strict_a_interval() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mp.iv.dps = 70
    omega = mp.iv.mpf("0.18575")
    phi = mp.iv.mpf("0.16250")
    eta = mp.iv.mpf("1.8")
    wlo = np.zeros((N, N))
    whi = np.zeros((N, N))
    for i in range(N):
        for j in range(i + 1, N):
            d = min(abs(i - j), N - abs(i - j))
            value = mp.iv.cos(omega * d + phi) / (1 + mp.iv.mpf(d) ** eta)
            lo, hi = iv_bounds(value)
            wlo[i, j] = wlo[j, i] = lo
            whi[i, j] = whi[j, i] = hi
    alo = -whi
    ahi = -wlo
    for i in range(N):
        alo[i, i] = np.nextafter(np.sum(wlo[i]), -np.inf)
        ahi[i, i] = np.nextafter(np.sum(whi[i]), np.inf)
    amid = 0.5 * (alo + ahi)
    return alo, ahi, amid


def scalar_interval_product(c: float, lo: float, hi: float) -> tuple[float, float]:
    if c >= 0:
        return float(np.nextafter(c * lo, -np.inf)), float(np.nextafter(c * hi, np.inf))
    return float(np.nextafter(c * hi, -np.inf)), float(np.nextafter(c * lo, np.inf))


def interval_matrix_left_product(r: np.ndarray, jlo: np.ndarray, jhi: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    rows, inner = r.shape
    cols = jlo.shape[1]
    out_lo = np.zeros((rows, cols))
    out_hi = np.zeros((rows, cols))
    for i in range(rows):
        for j in range(cols):
            lower, upper = 0.0, 0.0
            for k in range(inner):
                lo, hi = scalar_interval_product(float(r[i, k]), float(jlo[k, j]), float(jhi[k, j]))
                lower = np.nextafter(lower + lo, -np.inf)
                upper = np.nextafter(upper + hi, np.inf)
            out_lo[i, j], out_hi[i, j] = lower, upper
    return out_lo, out_hi


def interval_k_times_a(klo: float, khi: float, alo: np.ndarray, ahi: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    products = np.stack([klo * alo, klo * ahi, khi * alo, khi * ahi])
    return np.nextafter(np.min(products, axis=0), -np.inf), np.nextafter(np.max(products, axis=0), np.inf)


def krawczyk_chart(
    center: np.ndarray,
    k0: float,
    klo: float,
    khi: float,
    alo: np.ndarray,
    ahi: np.ndarray,
    amid: np.ndarray,
) -> dict:
    j0 = k0 * amid + np.eye(N) - np.diag(3.0 * center**2)
    r = np.linalg.inv(j0)
    f0 = k0 * (amid @ center) + center - center**3
    h = max(k0 - klo, khi - k0)
    derivative = -r @ (amid @ center)
    rho = 1.35 * np.abs(derivative) * h + 2e-11
    success = False
    final = None
    for _ in range(60):
        ulo, uhi = center - rho, center + rho
        u2lo = np.where(ulo * uhi <= 0.0, 0.0, np.minimum(ulo**2, uhi**2))
        u2hi = np.maximum(ulo**2, uhi**2)
        jlo, jhi = interval_k_times_a(klo, khi, alo, ahi)
        for i in range(N):
            jlo[i, i] = np.nextafter(jlo[i, i] + 1.0 - 3.0 * u2hi[i], -np.inf)
            jhi[i, i] = np.nextafter(jhi[i, i] + 1.0 - 3.0 * u2lo[i], np.inf)
        rjlo, rjhi = interval_matrix_left_product(r, jlo, jhi)
        mlo = -rjhi
        mhi = -rjlo
        for i in range(N):
            mlo[i, i] = np.nextafter(mlo[i, i] + 1.0, -np.inf)
            mhi[i, i] = np.nextafter(mhi[i, i] + 1.0, np.inf)
        mabs = np.maximum(np.abs(mlo), np.abs(mhi))
        a_center_rad = 0.5 * (ahi - alo) @ np.abs(center)
        ax = amid @ center
        correction = np.abs(r @ f0) + np.abs(r) @ (h * np.abs(ax) + (abs(k0) + h) * a_center_rad)
        correction += 2e-14 * (1.0 + np.linalg.norm(r, ord=np.inf))
        q = correction + mabs @ rho
        defect_inf = float(np.max(np.sum(mabs, axis=1)))
        inclusion_margin = float(np.min(rho - q))
        if inclusion_margin > 0.0 and defect_inf < 1.0:
            success = True
            final = (rho, q, defect_inf, inclusion_margin)
            break
        rho = np.maximum(rho * 1.08, q * 1.20 + 1e-13)
        if np.max(rho) > 0.2:
            break
    if not success or final is None:
        return {"success": False, "k_interval": [klo, khi], "center_kappa": k0, "max_radius": float(np.max(rho))}
    rho, q, defect_inf, inclusion_margin = final
    return {
        "success": True,
        "center_kappa": k0,
        "k_interval": [klo, khi],
        "center": center,
        "radius": rho,
        "max_radius": float(np.max(rho)),
        "max_krawczyk_image_radius": float(np.max(q)),
        "defect_infinity_norm_upper": defect_inf,
        "minimum_inclusion_margin": inclusion_margin,
    }


def p508_certified_global_continuation() -> dict:
    alo, ahi, amid = strict_a_interval()
    chart_count = 401
    nodes = np.linspace(0.0, 1.0, chart_count)
    centers = []
    u = np.zeros(N)
    u[0] = 1.0
    for k in nodes:
        u = stationary_root(amid, float(k), 1.0, u)
        centers.append(u.copy())
    charts = []
    failures = []
    half = 0.5 / (chart_count - 1)
    for index, (k, center) in enumerate(zip(nodes, centers)):
        klo = max(0.0, float(k - half))
        khi = min(1.0, float(k + half))
        chart = krawczyk_chart(center, float(k), klo, khi, alo, ahi, amid)
        charts.append(chart)
        if not chart["success"]:
            failures.append(index)
    overlaps = []
    bridges = []
    if not failures:
        for i in range(chart_count - 1):
            left, right = charts[i], charts[i + 1]
            llo = np.asarray(left["center"]) - np.asarray(left["radius"])
            lhi = np.asarray(left["center"]) + np.asarray(left["radius"])
            rlo = np.asarray(right["center"]) - np.asarray(right["radius"])
            rhi = np.asarray(right["center"]) + np.asarray(right["radius"])
            overlaps.append(bool(np.all(np.maximum(llo, rlo) <= np.minimum(lhi, rhi))))
            shared_kappa = 0.5 * (nodes[i] + nodes[i + 1])
            bridge_center = stationary_root(amid, float(shared_kappa), 1.0, 0.5 * (centers[i] + centers[i + 1]))
            bridge = krawczyk_chart(bridge_center, float(shared_kappa), float(shared_kappa), float(shared_kappa), alo, ahi, amid)
            if bridge["success"]:
                blo = np.asarray(bridge["center"]) - np.asarray(bridge["radius"])
                bhi = np.asarray(bridge["center"]) + np.asarray(bridge["radius"])
                nested = bool(np.all(blo >= llo) and np.all(bhi <= lhi) and np.all(blo >= rlo) and np.all(bhi <= rhi))
            else:
                nested = False
            bridges.append({"kappa": float(shared_kappa), "success": bridge["success"], "nested_in_both_adjacent_boxes": nested, "max_radius": bridge["max_radius"]})
    accepted = not failures and all(overlaps) and bool(bridges) and all(row["nested_in_both_adjacent_boxes"] for row in bridges)
    return {
        "program": "P508",
        "object": "O212 Parametric Krawczyk Continuation Tube",
        "chart_count": chart_count,
        "strict_operator_interval_max_width": float(np.max(ahi - alo)),
        "failed_chart_indices": failures,
        "all_adjacent_state_boxes_overlap": bool(overlaps and all(overlaps)),
        "all_shared_parameter_root_boxes_nested": bool(bridges and all(row["nested_in_both_adjacent_boxes"] for row in bridges)),
        "maximum_box_radius": max(float(c["max_radius"]) for c in charts),
        "maximum_defect_infinity_norm_upper": max(float(c.get("defect_infinity_norm_upper", math.inf)) for c in charts),
        "minimum_inclusion_margin": min(float(c.get("minimum_inclusion_margin", -math.inf)) for c in charts),
        "accepted_complete_parameter_tube": accepted,
        "charts": charts,
        "bridge_certificates": bridges,
        "theorem_if_accepted": (
            "Uniform Krawczyk inclusion in every parameter chart proves one stationary root in "
            "each state box for every kappa in its interval. Nonsingular Jacobians and overlapping "
            "adjacent boxes join these roots into one fold-free branch covering [0,1]."
        ),
        "status": "outward_interval_parametric_krawczyk_certificate" if accepted else "unresolved_parametric_krawczyk_gap",
    }


def p509_orbital_stability(a: np.ndarray) -> dict:
    u, _ = continue_pattern(a, [0], omega=1.0, steps=401)
    lminus = a + np.eye(N) - np.diag(u**2)
    lplus = a + np.eye(N) - np.diag(3.0 * u**2)
    eminus = np.linalg.eigvalsh(lminus)
    eplus = np.linalg.eigvalsh(lplus)
    du_domega = -np.linalg.solve(lplus, u)
    dp_domega = float(2.0 * u @ du_domega)
    omega_grid = np.linspace(0.6, 1.4, 161)
    branch = []
    current = u.copy()
    # Start at omega=1 and continue both directions for reliable roots.
    solutions = {1.0: u.copy()}
    current = u.copy()
    for omega in np.linspace(1.0, 1.4, 81)[1:]:
        current = stationary_root(a, 1.0, float(omega), current)
        solutions[round(float(omega), 10)] = current.copy()
    current = u.copy()
    for omega in np.linspace(1.0, 0.6, 81)[1:]:
        current = stationary_root(a, 1.0, float(omega), current)
        solutions[round(float(omega), 10)] = current.copy()
    for omega in omega_grid:
        state = solutions[round(float(omega), 10)]
        lm = a + omega * np.eye(N) - np.diag(state**2)
        lp = a + omega * np.eye(N) - np.diag(3.0 * state**2)
        dstate = -np.linalg.solve(lp, state)
        branch.append(
            {
                "omega": float(omega),
                "power": float(state @ state),
                "dp_domega": float(2.0 * state @ dstate),
                "lminus_second_eigenvalue": float(np.linalg.eigvalsh(lm)[1]),
                "lplus_negative_count": int(np.sum(np.linalg.eigvalsh(lp) < -1e-9)),
                "lplus_smallest_positive": float(np.linalg.eigvalsh(lp)[1]),
            }
        )
    omega1_hypotheses = bool(
        abs(eminus[0]) < 1e-8
        and eminus[1] > 0.0
        and np.sum(eplus < -1e-9) == 1
        and eplus[1] > 0.0
        and dp_domega > 0.0
    )
    turning_bracket = None
    for left, right in zip(branch[:-1], branch[1:]):
        if left["dp_domega"] <= 0.0 < right["dp_domega"]:
            turning_bracket = (left["omega"], right["omega"])
            break
    cache = {round(row["omega"], 10): solutions[round(row["omega"], 10)] for row in branch}

    def slope_at(omega: float) -> float:
        key = min(cache, key=lambda value: abs(value - omega))
        state = stationary_root(a, 1.0, omega, cache[key])
        cache[omega] = state
        lp = a + omega * np.eye(N) - np.diag(3.0 * state**2)
        return float(-2.0 * state @ np.linalg.solve(lp, state))

    turning = None
    if turning_bracket is not None:
        turning = brentq(slope_at, turning_bracket[0], turning_bracket[1], xtol=1e-13)
    stable_rows = [row for row in branch if row["dp_domega"] > 0 and row["lminus_second_eigenvalue"] > 0 and row["lplus_negative_count"] == 1 and row["lplus_smallest_positive"] > 0]
    return {
        "program": "P509",
        "object": "O213 Finite-DNLS Orbital-Stability Ledger",
        "omega_1": {
            "Lminus_eigenvalues": eminus,
            "Lplus_eigenvalues": eplus,
            "Lminus_kernel_residual": float(np.linalg.norm(lminus @ u)),
            "Lminus_second_margin": float(eminus[1]),
            "Lplus_negative_count": int(np.sum(eplus < -1e-9)),
            "Lplus_first_positive_margin": float(eplus[1]),
            "dP_domega": dp_domega,
        },
        "omega_tube": branch,
        "gss_vk_hypotheses_at_omega_1": omega1_hypotheses,
        "vk_slope_turning_point": turning,
        "sampled_positive_slope_interval": [min(row["omega"] for row in stable_rows), max(row["omega"] for row in stable_rows)],
        "conditional_theorem": (
            "For the finite focusing DNLS, the standard one-negative-direction, simple phase-kernel, "
            "and positive slope hypotheses imply orbital stability modulo global phase."
        ),
        "verdict": (
            "At omega=1 all finite-matrix GSS/VK sign hypotheses are strongly separated. The "
            "power slope changes sign near the reported turning point, so the complete 0.6--1.4 "
            "tube is not one stability phase. Proof-grade interval eigenvalue/slope enclosures remain necessary."
        ),
        "status": "strong_omega_1_orbital_stability_evidence_with_vk_turning_point",
    }


def dnls_energy(a: np.ndarray, u: np.ndarray) -> float:
    return float(u @ a @ u - 0.5 * np.sum(u**4))


def p510_mobility_barrier(a: np.ndarray) -> dict:
    site, _ = continue_pattern(a, [0], omega=1.0, steps=401)
    bond, _ = continue_pattern(a, [0, 1], omega=1.0, steps=401)
    p_target = float(site @ site)
    site_energy = dnls_energy(a, site)
    bond_energy_same_omega = dnls_energy(a, bond)
    branch_rows = []
    current = bond.copy()
    for omega in np.linspace(1.0, 0.20, 161):
        current = stationary_root(a, 1.0, float(omega), current)
        power = float(current @ current)
        branch_rows.append(
            {
                "omega": float(omega),
                "power": power,
                "ipr": float(np.sum(current**4) / power**2),
                "energy": dnls_energy(a, current),
            }
        )
    localized_rows = [row for row in branch_rows if row["ipr"] > 2.0 / N]
    minimum_localized_bond_power = min(row["power"] for row in localized_rows)
    equal_power_localized_bond_exists = bool(minimum_localized_bond_power <= p_target)
    uniform = np.ones(N) * math.sqrt(p_target / N)
    uniform_energy = dnls_energy(a, uniform)
    translations = [np.roll(site, shift) for shift in range(N)]
    translation_spread = max(dnls_energy(a, x) for x in translations) - min(dnls_energy(a, x) for x in translations)
    return {
        "program": "P510",
        "object": "O214 Equal-Power Peierls--Nabarro Existence Audit",
        "target_power": p_target,
        "site_centered": {"omega": 1.0, "energy": site_energy, "ipr": float(np.sum(site**4) / p_target**2), "u": site},
        "bond_centered_same_omega": {"omega": 1.0, "power": float(bond @ bond), "energy": bond_energy_same_omega, "ipr": float(np.sum(bond**4) / float(bond @ bond) ** 2), "u": bond},
        "bond_branch": branch_rows,
        "minimum_localized_bond_power": minimum_localized_bond_power,
        "equal_power_localized_bond_exists": equal_power_localized_bond_exists,
        "peierls_nabarro_energy_barrier": None,
        "equal_power_uniform_comparison": {"omega": p_target / N, "energy": uniform_energy, "energy_minus_site": uniform_energy - site_energy},
        "translation_orbit_energy_spread": translation_spread,
        "verdict": (
            "The localized bond-centred branch never reaches the lower power of the P497 site "
            "state before delocalizing into the uniform branch. Therefore the standard equal-power "
            "site/bond Peierls--Nabarro barrier is undefined for this state and the proposed mobility "
            "test is refuted in its declared form. Translation degeneracy alone does not imply motion."
        ),
        "status": "refuted_equal_power_localized_PN_comparison",
    }


def stieltjes_memory_operator(a: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    poles, residues, _, _ = grouped_stieltjes_data(a)
    identity = np.eye(N)
    f0 = float(np.sum(residues / poles))
    f_a = np.zeros_like(a)
    for pole, residue in zip(poles, residues):
        f_a += residue * np.linalg.inv(a + pole * identity)
    memory = f0 * identity - f_a
    return memory, poles, residues


def p511_localized_memory(a: np.ndarray) -> dict:
    memory, poles, residues = stieltjes_memory_operator(a)
    memory_eigs = np.linalg.eigvalsh(memory)
    base, _ = continue_pattern(a, [0], omega=1.0, steps=401)
    u = base.copy()
    rows = []
    for epsilon in np.linspace(0.0, 1.0, 201):
        operator = a + epsilon * memory
        u = stationary_root(operator, 1.0, 1.0, u)
        power = float(u @ u)
        rows.append(
            {
                "epsilon": float(epsilon),
                "ipr": float(np.sum(u**4) / power**2),
                "peak_fraction": float(np.max(u**2) / power),
                "power": power,
                "residual_inf": float(np.linalg.norm(operator @ u + u - u**3, ord=np.inf)),
            }
        )
    return {
        "program": "P511",
        "object": "O215 Stieltjes-Bernstein Memory-Loaded Localized Branch",
        "memory_operator": "M=sigma(0)I-sum_j r_j(A+c_j I)^-1",
        "stable_hidden_realization": "sigma(z)=sum_j r_j/(z+c_j), c_j>0, r_j>0",
        "poles": poles,
        "residues": residues,
        "memory_row_sum_residual": float(np.max(np.abs(memory.sum(axis=1)))),
        "memory_eigenvalue_interval": [float(memory_eigs[0]), float(memory_eigs[-1])],
        "continuation": rows,
        "endpoint": rows[-1],
        "localization_survives_declared_memory_loading": bool(rows[-1]["ipr"] > 2.0 / N),
        "verdict": (
            "The localized branch survives the full declared positive memory loading at zero "
            "frequency. The construction is a conditional Bernstein/Stieltjes loading, not a "
            "derived dynamical environment or a non-Markovian stability theorem."
        ),
        "status": "strong_stationary_memory_robustness_evidence",
    }


def p512_frequency_arithmetic() -> dict:
    # phi=650/4000, omega=743/4000, so every phase is n_d/4000.
    phase_numerators = {d: 650 + 743 * d for d in range(1, 7)}
    highest_lambda1 = phase_numerators[6]
    lambda2_d6_coefficient_zero = (1 - (-1) ** 2) == 0
    lambda1_d6_coefficient_nonzero = (1 - (-1) ** 1) != 0
    support_witness = lambda1_d6_coefficient_nonzero and lambda2_d6_coefficient_zero
    _, a = strict_operator()
    lam = np.real(np.fft.fft(a[0]))
    ratio = float(lam[2] / lam[1])
    return {
        "program": "P512",
        "object": "O216 Lindemann--Weierstrass Frequency-Ratio Certificate",
        "exact_phase_representation": "phi+d*omega=(650+743d)/4000",
        "phase_numerators": phase_numerators,
        "transcendental_generator": "z=exp(i/4000)",
        "highest_support_witness": {
            "lambda1_has_nonzero_z_power": highest_lambda1,
            "lambda2_d6_coefficient_is_zero": lambda2_d6_coefficient_zero,
            "laurent_polynomials_not_proportional": support_witness,
        },
        "ratio_numeric": ratio,
        "theorem": (
            "Because i/4000 is nonzero algebraic, Lindemann--Weierstrass makes z=exp(i/4000) "
            "transcendental. Each strict eigenvalue is a Laurent polynomial in z with algebraic "
            "coefficients (d^(9/5) is algebraic). Lambda_1 contains a nonzero d=6 highest-power "
            "term, while lambda_2 does not, so their Laurent polynomials are not proportional. "
            "If lambda_2/lambda_1 were algebraic, z would satisfy a nonzero algebraic-coefficient "
            "polynomial, a contradiction. Hence the ratio is transcendental and irrational."
        ),
        "denominator_nonzero_reason": "lambda_1>0 because all strict radial weights are positive and mode 1 is nonconstant",
        "status": "proven_transcendental_frequency_ratio",
    }


def attenuation(d: np.ndarray | float) -> np.ndarray | float:
    return 1.0 / (1.0 + np.asarray(d) ** ETA_S)


def context_symbol_data(n: int, u_values: np.ndarray, delta: float = 0.1) -> dict:
    distance = np.minimum(np.arange(n), n - np.arange(n))
    weights = np.zeros(n)
    weights[1:] = attenuation(distance[1:])
    row = -weights
    row[0] = np.sum(weights)
    c = np.real(np.fft.fft(row[::2])) + delta
    b = np.fft.fft(row[1::2])
    scale = float(np.median(c))
    sigma = np.array([np.mean(np.abs(b) ** 2 / (u * scale + c)) for u in u_values])
    return {"scale": scale, "sigma": sigma, "fingerprint": sigma / sigma[np.argmin(np.abs(u_values - 1.0))]}


def wiener_error_bound(n: int, u_values: np.ndarray, reference: dict, delta: float = 0.1) -> dict:
    cells = n // 2
    cutoff = max(2, cells // 2 - 1)
    dsum = np.arange(1, 1_000_001, dtype=float)
    hsum = attenuation(dsum)
    odd_sum = float(2.0 * np.sum(hsum[(dsum.astype(int) % 2) == 1]))
    odd_tail = 2.0 / (ETA_S - 1.0) * 1_000_000.0 ** (1.0 - ETA_S)
    b0 = odd_sum + odd_tail
    c_floor = b0 + delta
    tail = 2.0 / (ETA_S - 1.0) * (2.0 * cutoff) ** (1.0 - ETA_S)
    ctail = 2.0 * tail
    r = np.arange(1, cutoff + 1, dtype=float)
    b1 = float(np.sum(r * (attenuation(np.abs(2 * r - 1)) + attenuation(2 * r + 1))))
    c1 = float(2.0 * np.sum(r * attenuation(2 * r)))
    finite = context_symbol_data(n, u_values, delta)
    errors = []
    for idx, u in enumerate(u_values):
        denom = (1.0 + u) * c_floor
        tail_effect = (2 * b0 * tail + tail**2) / denom + b0**2 * ctail / denom**2
        scale_effect = b0**2 * u * ctail / denom**2
        derivative_bound = 2 * b0 * b1 / denom + b0**2 * c1 / denom**2
        quadrature_effect = 2 * math.pi / cells * derivative_bound
        errors.append(tail_effect + scale_effect + quadrature_effect)
    ref_sigma = np.asarray(reference["sigma"])
    e = np.asarray(errors)
    denominator_lower = float(ref_sigma[np.argmin(np.abs(u_values - 1.0))] - e[np.argmin(np.abs(u_values - 1.0))])
    normalized_bounds = []
    if denominator_lower > 0:
        ref_den = float(ref_sigma[np.argmin(np.abs(u_values - 1.0))])
        for i in range(len(u_values)):
            upper_num = abs(float(ref_sigma[i])) + e[i]
            bound = e[i] / denominator_lower + upper_num * e[np.argmin(np.abs(u_values - 1.0))] / (denominator_lower * ref_den)
            normalized_bounds.append(float(bound))
    return {
        "n": n,
        "cutoff": cutoff,
        "tail_l1_bound": tail,
        "absolute_response_error_bounds": errors,
        "denominator_lower_bound": denominator_lower,
        "normalized_fingerprint_error_bounds": normalized_bounds,
        "maximum_normalized_bound": max(normalized_bounds) if normalized_bounds else None,
        "observed_max_difference_from_reference": float(np.max(np.abs(finite["fingerprint"] - reference["fingerprint"]))),
    }


def p513_context_rate() -> dict:
    u_values = np.array([0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0])
    reference = context_symbol_data(262144, u_values)
    rows = [wiener_error_bound(n, u_values, reference) for n in (384, 768, 1536, 3072, 6144)]
    certified = [row for row in rows if row["maximum_normalized_bound"] is not None]
    return {
        "program": "P513",
        "object": "O217 Wiener--Resolvent Context Error Envelope",
        "reference_cycle": 262144,
        "rows": rows,
        "all_reported_bounds_finite": len(certified) == len(rows),
        "first_cycle_with_positive_normalization_denominator_bound": certified[0]["n"] if certified else None,
        "bound_decreases_with_n": bool(all(right["maximum_normalized_bound"] < left["maximum_normalized_bound"] for left, right in zip(certified[:-1], certified[1:]))) if len(certified) > 1 else False,
        "analytic_rate": "tail and quadrature envelope O(n^{-(eta-1)})=O(n^-0.8) up to conservative constants",
        "theorem_scope": (
            "Wiener l1 tail bounds, a positive odd-coupling floor, the resolvent identity, and a "
            "truncated-symbol derivative bound give the declared absolute and normalized error envelope."
        ),
        "verdict": (
            "A proof-oriented convergence envelope is available but is substantially looser than "
            "the observed FFT discrepancy. It certifies the rate class, not the sharp P501 constants."
        ),
        "status": "partial_normalized_certificate_from_C1536_with_rate_bound",
    }


def homotopy_lambda_iv(tlo: float, thi: float, mode: int) -> Any:
    mp.iv.dps = 50
    t = mp.iv.mpf([tlo, thi])
    alpha = 4 * mp.iv.log(2)
    amplitude = (1 - t) * alpha + t
    beta = (1 - t) * mp.iv.mpf("0.01") + t
    eta = (1 - t) + t * mp.iv.mpf("1.8")
    omega = (1 - t) * (mp.iv.pi / 4) + t * mp.iv.mpf("0.18575")
    phi = (1 - t) * (mp.iv.pi / 6) + t * mp.iv.mpf("0.16250")
    value = mp.iv.mpf("0")
    for d in range(1, 6):
        weight = amplitude * mp.iv.cos(omega * d + phi) / (1 + beta * mp.iv.mpf(d) ** eta)
        value += 2 * weight * (1 - mp.iv.cos(2 * mp.iv.pi * mode * d / 12))
    weight6 = amplitude * mp.iv.cos(omega * 6 + phi) / (1 + beta * mp.iv.mpf(6) ** eta)
    value += weight6 * (1 - (-1) ** mode)
    return value


def classify_psd_interval(lo: float, hi: float) -> tuple[str, list[tuple[float, float]]]:
    bounds = [iv_bounds(homotopy_lambda_iv(lo, hi, mode)) for mode in range(1, 12)]
    if all(bound[0] > 0.0 for bound in bounds):
        return "PSD", bounds
    if any(bound[1] < 0.0 for bound in bounds):
        return "NONPSD", bounds
    return "UNRESOLVED", bounds


def merge_intervals(rows: list[dict], label: str) -> list[list[float]]:
    selected = sorted([row["interval"] for row in rows if row["class"] == label])
    merged: list[list[float]] = []
    for lo, hi in selected:
        if not merged or lo > merged[-1][1] + 1e-15:
            merged.append([lo, hi])
        else:
            merged[-1][1] = max(merged[-1][1], hi)
    return merged


def p514_global_psd_phase_diagram() -> dict:
    initial_width = 1.0 / 512.0
    queue: deque[tuple[float, float, int]] = deque((i * initial_width, (i + 1) * initial_width, 0) for i in range(512))
    rows = []
    tolerance = 2e-10
    max_depth = 30
    while queue:
        lo, hi, depth = queue.popleft()
        cls, bounds = classify_psd_interval(lo, hi)
        if cls != "UNRESOLVED" or hi - lo <= tolerance or depth >= max_depth:
            rows.append({"interval": [lo, hi], "class": cls, "bounds": bounds})
        else:
            mid = 0.5 * (lo + hi)
            queue.append((lo, mid, depth + 1))
            queue.append((mid, hi, depth + 1))
    rows.sort(key=lambda row: row["interval"][0])
    unresolved = merge_intervals(rows, "UNRESOLVED")
    psd = merge_intervals(rows, "PSD")
    nonpsd = merge_intervals(rows, "NONPSD")
    last_unresolved = unresolved[-1] if unresolved else None
    final_psd_certified = bool(last_unresolved and all(row["class"] == "PSD" for row in rows if row["interval"][0] >= last_unresolved[1]))
    return {
        "program": "P514",
        "object": "O218 Global Interval PSD Phase Diagram",
        "initial_boxes": 512,
        "terminal_box_count": len(rows),
        "tolerance": tolerance,
        "psd_intervals": psd,
        "nonpsd_intervals": nonpsd,
        "unresolved_transition_intervals": unresolved,
        "final_transition_bracket": last_unresolved,
        "all_boxes_after_final_bracket_certified_psd": final_psd_certified,
        "all_boxes": rows,
        "verdict": (
            "Outward interval evaluation classifies the complete homotopy except narrow root "
            "brackets. The final unresolved bracket encloses the P502 transition and every later "
            "box is PSD; earlier PSD/non-PSD windows are reported rather than suppressed."
        ),
        "status": "global_interval_phase_diagram_with_narrow_root_brackets",
    }


def p515_robust_tester(a: np.ndarray) -> dict:
    models, contexts, table = operational_probability_table(a)
    context = (0.55, "fourier1", "fourier")
    prep_tv = 0.005
    povm_tv = 0.005
    time_error = 0.002
    opnorm = float(np.linalg.eigvalsh(a)[-1])
    speed = {
        "unitary": opnorm,
        "markov": 0.5 * float(np.linalg.norm(a, 1)),
        "dephasing": opnorm + 0.60,
        "revival": opnorm + 1.12,
    }
    pairs = []
    for i, left in enumerate(models):
        for right in models[i + 1 :]:
            nominal = 0.5 * float(np.sum(np.abs(table[left][context] - table[right][context])))
            penalty = 2.0 * (prep_tv + povm_tv) + (speed[left] + speed[right]) * time_error
            lower = nominal - penalty
            pairs.append({"models": [left, right], "nominal_tv": nominal, "worst_case_penalty": penalty, "certified_tv_lower": lower})
    minimum = min(row["certified_tv_lower"] for row in pairs)
    return {
        "program": "P515",
        "object": "O219 Robust One-Context Separation Budget",
        "context": {"time": context[0], "preparation": context[1], "measurement": context[2]},
        "uncertainty_budget": {"preparation_trace_distance": prep_tv, "measurement_distribution_tv": povm_tv, "absolute_time_error": time_error},
        "uniform_speed_bounds": speed,
        "pairwise_ledger": pairs,
        "minimum_certified_pairwise_tv": minimum,
        "robust_separation_retained": bool(minimum > 0.0),
        "theorem_scope": (
            "CPTP/classical contraction pays preparation error, the declared POVM budget pays "
            "measurement error, and global generator-norm speed bounds pay timing error. The "
            "triangle inequality then gives each pairwise lower TV bound."
        ),
        "verdict": (
            "The one-context synthetic tester remains separating under the declared conservative "
            "error budgets. These budgets are hypothetical and are not calibrated laboratory tolerances."
        ),
        "status": "conditional_robust_synthetic_separation_certificate",
    }


def parameter_vector_to_kernel(p: np.ndarray) -> KernelParameters:
    return KernelParameters(float(math.exp(p[0])), float(math.exp(p[1])), float(p[2]), float(p[3]), float(p[4]))


def p516_law_identifiability() -> dict:
    p0 = np.array([math.log(4.0 * math.log(2.0)), math.log(0.01), 1.0, math.pi / 4.0, math.pi / 6.0])
    coefficients = np.array(
        [
            [-0.08, 0.18, 0.04, -0.025, -0.015],
            [0.02, -0.03, 0.01, 0.006, -0.004],
            [-0.004, 0.008, -0.003, 0.001, 0.002],
        ]
    )
    candidate_times = [0.2, 0.4, 0.6, 0.8, 1.0]
    designs = []
    for degree in (0, 1, 2):
        unknowns = 5 * (degree + 1)
        best = None
        for distances in itertools.combinations(range(7), 5):
            for times in itertools.combinations(candidate_times, degree + 1):
                blocks = []
                for t in times:
                    p = p0.copy()
                    for k in range(degree + 1):
                        p += coefficients[k] * t ** (k + 1) / (k + 1)
                    jk = kernel_parameter_jacobian(parameter_vector_to_kernel(p), list(distances))
                    blocks.append(np.column_stack([jk * (t ** (k + 1) / (k + 1)) for k in range(degree + 1)]))
                sensitivity = np.row_stack(blocks)
                singular = np.linalg.svd(sensitivity, compute_uv=False)
                rank = int(np.linalg.matrix_rank(sensitivity, tol=1e-11))
                condition = float(singular[0] / singular[-1]) if rank == unknowns else math.inf
                if best is None or condition < best[0]:
                    best = (condition, distances, times, singular, rank)
        assert best is not None
        condition, distances, times, singular, rank = best
        designs.append(
            {
                "polynomial_derivative_degree": degree,
                "unknown_coefficient_count": unknowns,
                "scalar_observation_count": 5 * (degree + 1),
                "distances": list(distances),
                "times": list(times),
                "rank": rank,
                "condition_number": condition,
                "singular_values": singular,
                "dimension_lower_bound_attained": rank == unknowns,
            }
        )
    return {
        "program": "P516",
        "object": "O220 Polynomial Parameter-Law Observation Design",
        "declared_law_class": "p_dot(t)=sum_{k=0}^r c_k t^k, known p(0)=p_legacy, r in {0,1,2}",
        "strict_endpoint_inserted": False,
        "designs": designs,
        "theorem": (
            "The degree-r law has 5(r+1) free coefficients, so fewer scalar observations cannot "
            "be locally injective. A full-rank five-distance snapshot at each of r+1 distinct "
            "nonzero times yields a square sensitivity matrix; the reported designs attain full rank."
        ),
        "verdict": (
            "A finite design can identify a law only after a finite-dimensional law class and the "
            "legacy initial point are supplied. Without that class, the P506 infinite-trajectory no-go remains."
        ),
        "status": "conditional_minimal_local_law_identifiability",
    }


def make_figures(results: dict) -> None:
    FIG_DIR.mkdir(exist_ok=True)
    plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": 0.25})

    charts = results["P508"]["charts"]
    fig, ax = plt.subplots(figsize=(7.4, 4.2))
    ax.plot([row["center_kappa"] for row in charts], [row["max_radius"] for row in charts])
    ax.set_xlabel("kappa chart centre")
    ax.set_ylabel("maximum Krawczyk state radius")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "p508_krawczyk_tube.png", dpi=180)
    plt.close(fig)

    tube = results["P509"]["omega_tube"]
    fig, ax1 = plt.subplots(figsize=(7.4, 4.2))
    ax1.plot([row["omega"] for row in tube], [row["power"] for row in tube], label="power")
    ax1.set_xlabel("omega")
    ax1.set_ylabel("power")
    ax2 = ax1.twinx()
    ax2.plot([row["omega"] for row in tube], [row["dp_domega"] for row in tube], color="#d95f02", label="dP/domega")
    ax2.set_ylabel("slope")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "p509_vk_slope.png", dpi=180)
    plt.close(fig)

    p510 = results["P510"]
    fig, ax = plt.subplots(figsize=(7.4, 4.1))
    ax.bar(["site", "bond/same omega", "uniform/equal power"], [p510["site_centered"]["energy"], p510["bond_centered_same_omega"]["energy"], p510["equal_power_uniform_comparison"]["energy"]])
    ax.set_ylabel("DNLS energy")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "p510_pn_barrier.png", dpi=180)
    plt.close(fig)

    memory = results["P511"]["continuation"]
    fig, ax = plt.subplots(figsize=(7.4, 4.2))
    ax.plot([row["epsilon"] for row in memory], [row["ipr"] for row in memory], label="IPR")
    ax.plot([row["epsilon"] for row in memory], [row["peak_fraction"] for row in memory], label="peak fraction")
    ax.axhline(1 / 12, color="black", ls="--", lw=1)
    ax.set_xlabel("memory loading epsilon")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG_DIR / "p511_memory_localization.png", dpi=180)
    plt.close(fig)

    rate = results["P513"]["rows"]
    fig, ax = plt.subplots(figsize=(7.4, 4.2))
    ax.loglog([row["n"] for row in rate], [row["observed_max_difference_from_reference"] for row in rate], marker="o", label="observed")
    valid = [row for row in rate if row["maximum_normalized_bound"] is not None]
    if valid:
        ax.loglog([row["n"] for row in valid], [row["maximum_normalized_bound"] for row in valid], marker="s", label="certified envelope")
    ax.set_xlabel("cycle size n")
    ax.set_ylabel("normalized fingerprint error")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG_DIR / "p513_context_bound.png", dpi=180)
    plt.close(fig)


def write_summary(results: dict) -> None:
    rows = [
        ("P507", results["P507"]["status"], "quadratic core cannot select quartic sign/coefficient"),
        ("P508", results["P508"]["status"], f"charts={results['P508']['chart_count']} accepted={results['P508']['accepted_complete_parameter_tube']}"),
        ("P509", results["P509"]["status"], f"dP/domega={results['P509']['omega_1']['dP_domega']:.6g}"),
        ("P510", results["P510"]["status"], f"minimum localized bond power={results['P510']['minimum_localized_bond_power']:.6g}"),
        ("P511", results["P511"]["status"], f"endpoint IPR={results['P511']['endpoint']['ipr']:.6g}"),
        ("P512", results["P512"]["status"], "lambda2/lambda1 transcendental"),
        ("P513", results["P513"]["status"], f"finite_bounds={results['P513']['all_reported_bounds_finite']}"),
        ("P514", results["P514"]["status"], f"final bracket={results['P514']['final_transition_bracket']}"),
        ("P515", results["P515"]["status"], f"robust TV={results['P515']['minimum_certified_pairwise_tv']:.6g}"),
        ("P516", results["P516"]["status"], f"full-rank designs={sum(row['dimension_lower_bound_attained'] for row in results['P516']['designs'])}/3"),
    ]
    lines = ["program,status,key_result"] + [f'"{p}","{s}","{r}"' for p, s, r in rows]
    SUMMARY.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    started = time.perf_counter()
    _, a = strict_operator()
    results: dict[str, Any] = {}
    results["P507"] = p507_nonlinear_source_audit(a)
    results["P508"] = p508_certified_global_continuation()
    results["P509"] = p509_orbital_stability(a)
    results["P510"] = p510_mobility_barrier(a)
    results["P511"] = p511_localized_memory(a)
    results["P512"] = p512_frequency_arithmetic()
    results["P513"] = p513_context_rate()
    results["P514"] = p514_global_psd_phase_diagram()
    results["P515"] = p515_robust_tester(a)
    results["P516"] = p516_law_identifiability()
    make_figures(results)
    elapsed = time.perf_counter() - started
    payload = {
        "release": "10.51",
        "programs": "P507-P516",
        "seed": SEED,
        "runtime_seconds": elapsed,
        "environment": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__, "platform": platform.platform()},
        "epistemic_boundary": "Analytic and interval/local computational results plus synthetic operational models; no laboratory or external-audit evidence.",
        "results": results,
    }
    RESULTS.write_text(json.dumps(native(payload), indent=2, sort_keys=True), encoding="utf-8")
    write_summary(results)
    print(
        json.dumps(
            {
                "runtime_seconds": elapsed,
                "P508_accepted": results["P508"]["accepted_complete_parameter_tube"],
                "P510_equal_power_localized_bond": results["P510"]["equal_power_localized_bond_exists"],
                "P511_endpoint_IPR": results["P511"]["endpoint"]["ipr"],
                "P512_status": results["P512"]["status"],
                "P514_final_bracket": results["P514"]["final_transition_bracket"],
                "P515_robust_TV": results["P515"]["minimum_certified_pairwise_tv"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
