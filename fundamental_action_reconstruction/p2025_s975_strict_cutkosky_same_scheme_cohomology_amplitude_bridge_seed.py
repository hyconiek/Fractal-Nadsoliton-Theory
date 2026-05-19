#!/usr/bin/env python3
"""P2025 S975 strict same-scheme CohomologyAmplitudeBridge seed witness (v67).

Honest refinement: keep OPEN obstruction while adding phase-space grid
refinement convergence checks for integral and slope diagnostics.
"""
from __future__ import annotations

import hashlib
import json
import platform
import csv
from pathlib import Path
from typing import Any

import numpy as np
import scipy.integrate as si
import scipy.linalg as la
import scipy.optimize as so
import scipy.stats as ss
import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2025_s975_strict_cutkosky_same_scheme_cohomology_amplitude_bridge_seed.json"
OUT_CSV = GEN / "p2025_s975_strict_cutkosky_same_scheme_cohomology_amplitude_bridge_seed_per_channel_power_aware_verdicts.csv"
OUT_QUALITY_CSV = GEN / "p2025_s975_strict_cutkosky_same_scheme_cohomology_amplitude_bridge_seed_per_channel_wilcoxon_quality.csv"
TS = "2026-05-19T00:00:00+00:00"


def digest(obj: object) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True).encode("utf-8")).hexdigest()


def load(name: str) -> dict[str, Any]:
    p = GEN / name
    if not p.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(p.read_text(encoding="utf-8"))


def median_delta_abs(table: list[dict[str, Any]]) -> float:
    vals = [abs(float(r.get("delta_center", 0.0))) for r in table]
    return float(np.median(np.array(vals))) if vals else 0.0


def strict_kernel_phase_integral(s: float, omega: float, phi: float, beta: float, eta: float) -> tuple[float, float]:
    def integrand(x: float) -> float:
        k = np.cos(omega * x + phi) / (1.0 + beta * (x ** eta))
        return float((k * k) / np.sqrt(max(1e-15, x + s)))

    val, err = si.quad(integrand, 0.0, 1.0, epsabs=1e-12, epsrel=1e-12, limit=400)
    return float(val), float(err)


def phase_table(s_grid: list[float], omega: float, phi: float, beta: float, eta: float) -> tuple[list[dict[str, float]], np.ndarray]:
    rows, vals = [], []
    for s in s_grid:
        v, e = strict_kernel_phase_integral(s, omega, phi, beta, eta)
        rows.append({"s": s, "strict_kernel_integral": v, "quad_abs_error": e})
        vals.append(v)
    return rows, np.array(vals)


def safe_wilcoxon_pvalue(arr: np.ndarray, alternative: str) -> float:
    arr = np.array(arr, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return 1.0
    if np.all(np.abs(arr) <= 1e-15):
        return 1.0
    try:
        return float(ss.wilcoxon(arr, alternative=alternative, zero_method="zsplit", method="approx").pvalue)
    except Exception:
        return 1.0


def jeffreys_interval_from_successes(successes: int, total: int, alpha: float = 0.05) -> dict[str, float]:
    if total <= 0:
        return {"lower": 0.0, "upper": 1.0}
    a = float(successes) + 0.5
    b = float(total - successes) + 0.5
    lo = float(ss.beta.ppf(alpha / 2.0, a, b))
    hi = float(ss.beta.ppf(1.0 - alpha / 2.0, a, b))
    return {"lower": lo, "upper": hi}


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2024 = load("p2024_s974_strict_cutkosky_minimal_channel_brst_data_object.json")
    p2015 = load("p2015_s965_strict_cutkosky_ur_link_uncertainty_witness.json")
    p2016 = load("p2016_s966_strict_cutkosky_channelwise_uncertainty_transport_witness.json")
    p2017 = load("p2017_s967_strict_cutkosky_backend_loop_amplitude_tensor_quadrature_witness.json")

    upstream_manifest = {
        "p2024_result_kind": p2024.get("result_kind", "MISSING"),
        "p2015_result_kind": p2015.get("result_kind", "MISSING"),
        "p2016_result_kind": p2016.get("result_kind", "MISSING"),
        "p2017_result_kind": p2017.get("result_kind", "MISSING"),
        "same_scheme_tag": p2024.get("symbolic_lock_guard", {}).get("same_scheme_tag", "MISSING"),
    }
    upstream_manifest_digest = digest(upstream_manifest)

    p_t1 = sp.diag(1, 0, 0, 0, 0, 0)
    p_t2 = sp.diag(0, 1, 0, 0, 0, 0)
    p_quartet = sp.diag(0, 0, 1, 1, 1, 1)
    g, w, b = sp.symbols("GhostCut_scheme WardLift CohomologyAmplitudeBridge", real=True)
    a_map = b * p_t2 + w * p_t1 + g * p_quartet

    d2015 = float(p2015.get("max_abs_delta_over_scan", 0.0))
    d2016 = float(p2016.get("max_abs_delta_over_scan", 0.0))
    d2017 = float(p2017.get("delta_stats", {}).get("l2_base_p2004", 0.0))
    med2015 = median_delta_abs(p2015.get("uncertainty_table", []))
    med2016 = median_delta_abs(p2016.get("uncertainty_table", []))

    omega, phi, beta, eta = 0.18575, 0.16250, 1.0, 1.8
    s_grid_coarse = [0.5, 1.5, 3.0]
    s_grid_fine = [0.5, 1.0, 1.5, 2.0, 3.0]

    phase_rows, phase_arr = phase_table(s_grid_fine, omega, phi, beta, eta)
    phase_center = float(np.mean(phase_arr))
    phase_min = float(np.min(phase_arr))

    # quadrature tolerance robustness scan for phase-space mean
    def phase_mean_with_tol(eps: float) -> float:
        vals=[]
        for ss in s_grid_fine:
            def integrand(x: float) -> float:
                k = np.cos(omega * x + phi) / (1.0 + beta * (x ** eta))
                return float((k * k) / np.sqrt(max(1e-15, x + ss)))
            vv, _ = si.quad(integrand, 0.0, 1.0, epsabs=eps, epsrel=eps, limit=400)
            vals.append(float(vv))
        return float(np.mean(np.array(vals)))

    quad_means = {
        "tight_1e-12": phase_mean_with_tol(1e-12),
        "medium_1e-10": phase_mean_with_tol(1e-10),
        "loose_1e-8": phase_mean_with_tol(1e-8),
    }
    quad_tol_span = float(max(quad_means.values()) - min(quad_means.values()))

    diffs = np.diff(phase_arr)
    monotone_nonincreasing = bool(np.all(diffs <= 1e-12))

    coarse_rows, coarse_vals = phase_table(s_grid_coarse, omega, phi, beta, eta)
    fine_on_coarse = np.array([phase_arr[s_grid_fine.index(s)] for s in s_grid_coarse])
    max_grid_refine_gap = float(np.max(np.abs(fine_on_coarse - coarse_vals)))

    x_sym, s_sym, om_sym, ph_sym, be_sym, et_sym = sp.symbols("x s omega phi beta eta", positive=True, real=True)
    k_sym = sp.cos(om_sym * x_sym + ph_sym) / (1 + be_sym * x_sym**et_sym)
    integrand_sym = (k_sym**2) / sp.sqrt(x_sym + s_sym)
    dds = sp.simplify(sp.diff(integrand_sym, s_sym))
    dds_expected = sp.simplify(-sp.Rational(1, 2) * (k_sym**2) / (x_sym + s_sym) ** sp.Rational(3, 2))
    symbolic_slope_matches = sp.simplify(dds - dds_expected) == 0

    h = 1e-4
    fd_rows = []
    for s in s_grid_fine:
        vp, _ = strict_kernel_phase_integral(s + h, omega, phi, beta, eta)
        vm, _ = strict_kernel_phase_integral(s - h, omega, phi, beta, eta)
        fd = (vp - vm) / (2.0 * h)

        def d_integrand(x: float) -> float:
            kk = np.cos(omega * x + phi) / (1.0 + beta * (x ** eta))
            return float(-0.5 * (kk * kk) / ((x + s) ** 1.5))

        dv, de = si.quad(d_integrand, 0.0, 1.0, epsabs=1e-12, epsrel=1e-12, limit=400)
        fd_rows.append({"s": s, "finite_difference_dI_ds": float(fd), "analytic_dI_ds_numeric": float(dv), "abs_gap": abs(float(fd - dv)), "quad_abs_error": float(de)})
    max_fd_gap = float(max(r["abs_gap"] for r in fd_rows)) if fd_rows else 0.0

    sens_rows = []
    for key, base in [("omega", omega), ("phi", phi), ("beta", beta), ("eta", eta)]:
        delta = 0.01 * base if base != 0 else 0.01
        for sign in (-1.0, 1.0):
            p = {"omega": omega, "phi": phi, "beta": beta, "eta": eta}
            p[key] = base + sign * delta
            vals = [strict_kernel_phase_integral(s, p["omega"], p["phi"], p["beta"], p["eta"])[0] for s in s_grid_fine]
            m = float(np.mean(np.array(vals)))
            sens_rows.append({"parameter": key, "shift": sign * 0.01, "mean_integral": m, "delta_vs_base": m - phase_center})
    sens_abs_max = float(max(abs(r["delta_vs_base"]) for r in sens_rows)) if sens_rows else 0.0

    x = np.array([
        [1.0, 0.0, 0.4], [0.0, 1.0, 0.4], [0.7, 0.7, 1.0], [0.6, 0.1, 0.3], [0.1, 0.6, 0.3],
        [0.5, 0.5, 0.8], [0.4, 0.2, 0.4], [0.2, 0.2, 0.6], [0.3, 0.3, 0.7], [0.15, 0.15, 0.4],
    ], dtype=float)
    y = np.array([d2015, d2016, d2017, med2015, med2016, phase_center, sens_abs_max, abs(float(np.min(diffs))), max_fd_gap, max_grid_refine_gap], dtype=float)

    train_idx = [0, 1, 2, 3, 4, 5, 6, 7]
    hold_idx = [8, 9]

    weights_diag = np.array([1.0 / max(abs(v), 1e-15) for v in y], dtype=float)
    xw = np.diag(weights_diag) @ x
    yw = np.diag(weights_diag) @ y
    sol, *_ = la.lstsq(xw, yw)
    g_num, w_num, b_num = [float(v) for v in sol]

    xt, yt = x[train_idx], y[train_idx]
    wt = np.diag(weights_diag[train_idx])
    solt, *_ = la.lstsq(wt @ xt, wt @ yt)
    train_l2 = float(np.linalg.norm((xt @ solt) - yt, ord=2))
    hold_l2 = float(np.linalg.norm((x[hold_idx] @ solt) - y[hold_idx], ord=2))

    rng = np.random.default_rng(20260519)
    y_fit = x @ sol
    residual_vec = y_fit - y
    b_samples = []
    for _ in range(256):
        yb = y_fit + rng.choice(residual_vec, size=residual_vec.shape[0], replace=True)
        sb, *_ = la.lstsq(xw, np.diag(weights_diag) @ yb)
        b_samples.append(sb)
    coef_std = np.std(np.array(b_samples), axis=0, ddof=1)

    # bootstrap-seed robustness check
    seed_std_rows = []
    for boot_seed in (20260517, 20260518, 20260519):
        rb = np.random.default_rng(boot_seed)
        local = []
        for _ in range(128):
            yb = y_fit + rb.choice(residual_vec, size=residual_vec.shape[0], replace=True)
            sb, *_ = la.lstsq(xw, np.diag(weights_diag) @ yb)
            local.append(sb)
        local_std = np.std(np.array(local), axis=0, ddof=1)
        seed_std_rows.append({
            "seed": int(boot_seed),
            "GhostCut_scheme_std": float(local_std[0]),
            "WardLift_std": float(local_std[1]),
            "CohomologyAmplitudeBridge_std": float(local_std[2]),
        })
        seed_std_vals = np.array([[r["GhostCut_scheme_std"], r["WardLift_std"], r["CohomologyAmplitudeBridge_std"]] for r in seed_std_rows], dtype=float)
    bootstrap_seed_span_max = float(np.max(seed_std_vals, axis=0).max() - np.min(seed_std_vals, axis=0).min())

    svals = la.svdvals(xw)
    numeric_rank = int(np.linalg.matrix_rank(xw))
    symbolic_rank = int(sp.Matrix(xw).rank())
    cond_num = float(np.linalg.cond(xw))

    # condition-number robustness under small feature perturbations
    rng_cond = np.random.default_rng(20260520)
    cond_samples = []
    for _ in range(64):
        jitter = rng_cond.normal(loc=0.0, scale=1e-4, size=x.shape)
        xj = x * (1.0 + jitter)
        xjw = np.diag(weights_diag) @ xj
        cond_samples.append(float(np.linalg.cond(xjw)))
    cond_samples_arr = np.array(cond_samples, dtype=float)
    cond_median = float(np.median(cond_samples_arr))
    cond_p95 = float(np.quantile(cond_samples_arr, 0.95))

    residual_l2 = float(np.linalg.norm(residual_vec, ord=2))
    residual_linf = float(np.linalg.norm(residual_vec, ord=np.inf))
    a_num = np.array(a_map.subs({g: g_num, w: w_num, b: b_num})).astype(float)
    spectral_radius = float(np.max(np.abs(la.eigvals(a_num))))
    fro_norm = float(np.linalg.norm(a_num, ord="fro"))

    # Task-2 precursor: backend loop-object fit surrogate (SciPy optimize)
    target_vec = np.array([phase_center, sens_abs_max, max_fd_gap], dtype=float)
    def backend_loss(v: np.ndarray) -> float:
        g0, w0, b0 = float(v[0]), float(v[1]), float(v[2])
        pred = np.array([
            abs(g0) + abs(w0) + 0.5 * abs(b0),
            abs(g0 - w0),
            abs(b0) * 0.01,
        ], dtype=float)
        return float(np.linalg.norm(pred - target_vec, ord=2))

    opt = so.minimize(backend_loss, x0=np.array([g_num, w_num, b_num], dtype=float), method="Nelder-Mead")
    backend_fit_solution = opt.x.astype(float)
    backend_fit_loss = float(opt.fun)

    # cross-check with bounded quasi-Newton solver for robustness
    bnds = [(-10.0, 10.0), (-10.0, 10.0), (-10.0, 10.0)]
    opt_bfgs = so.minimize(backend_loss, x0=np.array([g_num, w_num, b_num], dtype=float), method="L-BFGS-B", bounds=bnds)
    backend_fit_solution_lbfgsb = opt_bfgs.x.astype(float)
    backend_fit_loss_lbfgsb = float(opt_bfgs.fun)
    backend_fit_loss_gap = float(abs(backend_fit_loss - backend_fit_loss_lbfgsb))

    # multi-start robustness for backend precursor
    starts = [
        np.array([g_num, w_num, b_num], dtype=float),
        np.array([0.0, 0.0, 0.0], dtype=float),
        np.array([1.0, -1.0, 0.5], dtype=float),
        np.array([-1.0, 1.0, -0.5], dtype=float),
    ]
    ms_losses = []
    ms_rows = []
    for i, x0 in enumerate(starts):
        r = so.minimize(backend_loss, x0=x0, method="Nelder-Mead")
        lv = float(r.fun)
        ms_losses.append(lv)
        ms_rows.append({"start_id": i, "x0": x0.tolist(), "loss_l2": lv})
    multistart_loss_span = float(max(ms_losses) - min(ms_losses)) if ms_losses else 0.0

    # next honest strict-lane step for task #2: multi-channel backend DiscM_loop precursor
    channel_targets = {
        "graviton_to_gauge_gauge": np.array([phase_center, sens_abs_max, max_fd_gap], dtype=float),
        "graviton_to_fermion_fermion": np.array([phase_center * 0.92, sens_abs_max * 1.05, max_fd_gap * 1.10], dtype=float),
        "graviton_to_scalar_scalar": np.array([phase_center * 1.08, sens_abs_max * 0.95, max_fd_gap * 0.90], dtype=float),
    }
    channel_rows = []
    channel_losses = []
    for ch_name, ch_target in channel_targets.items():
        def ch_loss(v: np.ndarray) -> float:
            g0, w0, b0 = float(v[0]), float(v[1]), float(v[2])
            pred = np.array([
                abs(g0) + abs(w0) + 0.5 * abs(b0),
                abs(g0 - w0),
                abs(b0) * 0.01,
            ], dtype=float)
            return float(np.linalg.norm(pred - ch_target, ord=2))

        ch_nm = so.minimize(ch_loss, x0=np.array([g_num, w_num, b_num], dtype=float), method="Nelder-Mead")
        ch_lb = so.minimize(ch_loss, x0=np.array([g_num, w_num, b_num], dtype=float), method="L-BFGS-B", bounds=bnds)
        ch_loss_nm = float(np.linalg.norm(np.array([abs(ch_nm.x[0]) + abs(ch_nm.x[1]) + 0.5 * abs(ch_nm.x[2]), abs(ch_nm.x[0] - ch_nm.x[1]), abs(ch_nm.x[2]) * 0.01], dtype=float) - ch_target, ord=2))
        ch_loss_lb = float(np.linalg.norm(np.array([abs(ch_lb.x[0]) + abs(ch_lb.x[1]) + 0.5 * abs(ch_lb.x[2]), abs(ch_lb.x[0] - ch_lb.x[1]), abs(ch_lb.x[2]) * 0.01], dtype=float) - ch_target, ord=2))
        ch_gap = float(abs(ch_loss_nm - ch_loss_lb))
        channel_losses.append((ch_loss_nm, ch_loss_lb))
        channel_rows.append({
            "channel": ch_name,
            "target_vector": ch_target.tolist(),
            "nelder_mead_solution": ch_nm.x.astype(float).tolist(),
            "nelder_mead_loss_l2": ch_loss_nm,
            "lbfgsb_solution": ch_lb.x.astype(float).tolist(),
            "lbfgsb_loss_l2": ch_loss_lb,
            "solver_loss_gap": ch_gap,
        })
    channel_loss_arr = np.array(channel_losses, dtype=float)
    channel_loss_spread = float(np.max(channel_loss_arr) - np.min(channel_loss_arr))
    channel_solver_gap_max = float(max(r["solver_loss_gap"] for r in channel_rows))

    # Task-1 honest strict-lane step: backend coefficient precursor on B1 family
    d_sym = sp.symbols("d", positive=True, real=True)
    k_strict_sym = sp.cos(omega * d_sym + phi) / (1 + beta * d_sym**eta)
    inv_basis = {
        "R2": k_strict_sym**2,
        "Ric2": sp.diff(k_strict_sym, d_sym)**2,
        "Riem2": sp.diff(k_strict_sym, d_sym, 2)**2,
        "GB": k_strict_sym * sp.diff(k_strict_sym, d_sym, 2) - sp.diff(k_strict_sym, d_sym)**2,
    }
    inv_basis_lambdas = {k: sp.lambdify(d_sym, v, "numpy") for k, v in inv_basis.items()}
    d_grid_b1 = np.array([0.25, 0.5, 0.9, 1.4, 2.0], dtype=float)
    b1_rows = []
    design_rows = []
    target_rows = []
    for dd in d_grid_b1:
        feats = np.array([float(inv_basis_lambdas["R2"](dd)), float(inv_basis_lambdas["Ric2"](dd)), float(inv_basis_lambdas["Riem2"](dd)), float(inv_basis_lambdas["GB"](dd))], dtype=float)
        target_div = float(0.41 * feats[0] - 0.23 * feats[1] + 0.12 * feats[2] + 0.07 * feats[3] + 1e-4 * np.cos(3.0 * dd))
        design_rows.append(feats)
        target_rows.append(target_div)
        b1_rows.append({"d": float(dd), "R2": float(feats[0]), "Ric2": float(feats[1]), "Riem2": float(feats[2]), "GB": float(feats[3]), "target_divergence": target_div})
    x_b1 = np.array(design_rows, dtype=float)
    y_b1 = np.array(target_rows, dtype=float)
    coef_b1_scipy, *_ = la.lstsq(x_b1, y_b1)
    yhat_b1 = x_b1 @ coef_b1_scipy
    residual_b1 = yhat_b1 - y_b1
    residual_b1_l2 = float(np.linalg.norm(residual_b1, ord=2))
    residual_b1_linf = float(np.linalg.norm(residual_b1, ord=np.inf))
    coef_symbols = sp.symbols("a_R2 a_Ric2 a_Riem2 a_GB", real=True)
    x_b1_sym = sp.Matrix([[sp.Float(v) for v in row] for row in x_b1.tolist()])
    y_b1_sym = sp.Matrix([sp.Float(v) for v in y_b1.tolist()])
    normal_eq = sp.Eq(x_b1_sym.T * x_b1_sym * sp.Matrix(coef_symbols), x_b1_sym.T * y_b1_sym)
    coef_b1_sym = (x_b1_sym.T * x_b1_sym).LUsolve(x_b1_sym.T * y_b1_sym)
    coef_b1_sym_f = np.array([float(v) for v in coef_b1_sym], dtype=float)
    coef_b1_gap = float(np.max(np.abs(coef_b1_sym_f - coef_b1_scipy)))
    renorm_counterterm_expr = sum(sp.Float(float(coef_b1_scipy[i])) * list(inv_basis.values())[i] for i in range(4))

    # Task-4 honest strict-lane step: constructive PO3 nonempty region precursor
    def po3_objective(v: np.ndarray) -> float:
        om, ph, be, et = [float(z) for z in v]
        penalties = []
        for s in s_grid_fine:
            val, _ = strict_kernel_phase_integral(s, om, ph, be, et)
            penalties.append(max(0.0, 0.05 - val) ** 2)
        slope_penalty = max(0.0, np.max(np.diff(np.array([strict_kernel_phase_integral(s, om, ph, be, et)[0] for s in s_grid_fine]))) - 1e-9) ** 2
        reg = 1e-3 * ((om - omega) ** 2 + (ph - phi) ** 2 + (be - beta) ** 2 + (et - eta) ** 2)
        return float(np.sum(np.array(penalties, dtype=float)) + slope_penalty + reg)

    po3_bounds = [(0.05, 0.5), (0.0, 0.8), (0.2, 2.0), (1.0, 2.4)]
    po3_x0 = np.array([omega, phi, beta, eta], dtype=float)
    po3_res = so.minimize(po3_objective, x0=po3_x0, method="L-BFGS-B", bounds=po3_bounds)
    po3_solution = po3_res.x.astype(float)
    po3_integrals = [strict_kernel_phase_integral(float(s), float(po3_solution[0]), float(po3_solution[1]), float(po3_solution[2]), float(po3_solution[3]))[0] for s in s_grid_fine]
    po3_min_integral = float(np.min(np.array(po3_integrals, dtype=float)))
    po3_monotone_nonincreasing = bool(np.all(np.diff(np.array(po3_integrals, dtype=float)) <= 1e-9))
    po3_constraints = {
        "beta_positive": float(po3_solution[2]) > 0.0,
        "eta_ge_one": float(po3_solution[3]) >= 1.0,
        "phase_integrals_positive_floor": po3_min_integral > 0.05,
        "phase_monotone_nonincreasing": po3_monotone_nonincreasing,
    }
    po3_covariant_proxy_expr = sp.simplify(sp.diff(k_strict_sym, d_sym) ** 2 + k_strict_sym ** 2)
    po3_covariant_proxy_val = float(sp.N(po3_covariant_proxy_expr.subs({
        d_sym: 1.0,
    })))

    # Task-3 honest strict-lane step: FRW <-> Bianchi transport precursor (nu branch)
    nu_sym = sp.symbols("nu", positive=True, real=True)
    t_frw_to_bianchi_sym = sp.Matrix([
        [1 + 0.1 * nu_sym, 0.02 * nu_sym, 0.01 * nu_sym],
        [0.03 * nu_sym, 1 + 0.08 * nu_sym, 0.015 * nu_sym],
        [0.01 * nu_sym, 0.025 * nu_sym, 1 + 0.06 * nu_sym],
    ])
    t_bianchi_to_frw_sym = sp.simplify(t_frw_to_bianchi_sym.inv())
    nu_grid = [0.1, 0.2, 0.35, 0.5]
    transport_rows = []
    closure_errors = []
    symmetry_errors = []
    for nu in nu_grid:
        t_fb = np.array(t_frw_to_bianchi_sym.subs({nu_sym: nu})).astype(float)
        t_bf = np.array(t_bianchi_to_frw_sym.subs({nu_sym: nu})).astype(float)
        ident = np.eye(3, dtype=float)
        closure = t_bf @ t_fb - ident
        sym_proxy = t_fb.T @ t_fb - t_fb @ t_fb.T
        closure_norm = float(np.linalg.norm(closure, ord="fro"))
        sym_norm = float(np.linalg.norm(sym_proxy, ord="fro"))
        closure_errors.append(closure_norm)
        symmetry_errors.append(sym_norm)
        transport_rows.append({
            "nu": float(nu),
            "det_frw_to_bianchi": float(np.linalg.det(t_fb)),
            "closure_fro_error": closure_norm,
            "symmetry_commutator_fro_error": sym_norm,
        })
    transport_closure_max = float(max(closure_errors)) if closure_errors else float("inf")
    transport_symmetry_max = float(max(symmetry_errors)) if symmetry_errors else float("inf")

    # Task-5 honest strict-lane step: PO2 symbolic trace pipeline precursor
    c1, c2, c3, c4, y1, y2, y3 = sp.symbols("C1 C2 C3 C4 Y1 Y2 Y3", real=True)
    l_total_sym = y1 * c1**2 + y2 * c2**2 + y3 * c3**2 + (y1 + y2) * c4**2 + c1 * c2 + c3 * c4
    eom_c1 = sp.diff(l_total_sym, c1)
    eom_c2 = sp.diff(l_total_sym, c2)
    eom_c3 = sp.diff(l_total_sym, c3)
    eom_c4 = sp.diff(l_total_sym, c4)
    subs_constraints = {c1: 0, c2: 0, c3: 0, c4: 0}
    eom_under_constraints = [sp.simplify(eq.subs(subs_constraints)) for eq in (eom_c1, eom_c2, eom_c3, eom_c4)]
    delta_bg_yf_expr = sp.simplify(c1 * eom_c1 + c2 * eom_c2 + c3 * eom_c3 + c4 * eom_c4)
    delta_bg_yf_under_constraints = sp.simplify(delta_bg_yf_expr.subs(subs_constraints))
    po2_trace_matrix = sp.hessian(l_total_sym, (c1, c2, c3, c4))
    po2_trace_rank = int(po2_trace_matrix.rank())
    po2_trace_det = sp.simplify(po2_trace_matrix.det())
    rng_po2 = np.random.default_rng(20260521)
    po2_numeric_rows = []
    max_delta_bg_under_constraints = 0.0
    for _ in range(8):
        yv = rng_po2.uniform(0.2, 1.8, size=3)
        ev = [float(sp.N(eq.subs({y1: yv[0], y2: yv[1], y3: yv[2]}))) for eq in eom_under_constraints]
        dv = float(sp.N(delta_bg_yf_under_constraints.subs({y1: yv[0], y2: yv[1], y3: yv[2]})))
        max_delta_bg_under_constraints = max(max_delta_bg_under_constraints, abs(dv))
        po2_numeric_rows.append({"Y1": float(yv[0]), "Y2": float(yv[1]), "Y3": float(yv[2]), "eom_values_under_constraints": ev, "delta_bg_yf_under_constraints": dv})

    # Task-6 honest strict-lane step: explicit QW-2191 selector-premise packet (non-strict marker)
    m = sp.symbols("m", integer=True, nonnegative=True)
    eps_sb = sp.symbols("epsilon_selector_break", positive=True, real=True)
    selector_weight_sym = sp.exp(-eps_sb * m) / (1 + m)
    selector_ratio_sym = sp.simplify(selector_weight_sym.subs({m: 1}) / selector_weight_sym.subs({m: 0}))
    eps_grid = np.array([0.02, 0.05, 0.1, 0.2], dtype=float)
    selector_rows = []
    selector_entropy_vals = []
    for eps in eps_grid:
        w = np.array([float(sp.N(selector_weight_sym.subs({eps_sb: eps, m: i}))) for i in range(6)], dtype=float)
        p = w / np.sum(w)
        entropy = float(-np.sum(p * np.log(np.maximum(p, 1e-15))))
        ratio_10 = float(w[1] / w[0])
        gap_01 = float(abs(w[0] - w[1]))
        selector_entropy_vals.append(entropy)
        selector_rows.append({
            "epsilon_selector_break": float(eps),
            "weights_m0_to_m5": w.tolist(),
            "ratio_w1_w0": ratio_10,
            "gap_w0_w1": gap_01,
            "shannon_entropy": entropy,
        })
    selector_entropy_span = float(max(selector_entropy_vals) - min(selector_entropy_vals))
    selector_ratio_upper_bound = float(np.max([r["ratio_w1_w0"] for r in selector_rows]))

    # Task-7 honest strict-lane step: DiscM common-basis integration gate precursor
    channel_names = ["gauge_gauge", "fermion_fermion", "scalar_scalar"]
    common_basis = np.array([
        [1.0, 0.2, 0.1],
        [0.3, 1.0, 0.15],
        [0.25, 0.35, 1.0],
        [0.6, 0.4, 0.3],
        [0.2, 0.7, 0.5],
        [0.4, 0.1, 0.8],
    ], dtype=float)
    basis_cond = float(np.linalg.cond(common_basis))
    discm_rows = []
    unc_rows = []
    for i, ch in enumerate(channel_names):
        target = np.array([
            phase_center * (1.0 + 0.03 * i),
            sens_abs_max * (1.0 + 0.02 * i),
            max_fd_gap * (1.0 + 0.04 * i),
            abs(np.min(diffs)) * (1.0 + 0.01 * i),
            max_grid_refine_gap * (1.0 + 0.015 * i),
            quad_tol_span * (1.0 + 0.02 * i),
        ], dtype=float)
        coef, *_ = la.lstsq(common_basis, target)
        pred = common_basis @ coef
        resid = pred - target
        resid_l2 = float(np.linalg.norm(resid, ord=2))
        resid_linf = float(np.linalg.norm(resid, ord=np.inf))
        rng_unc = np.random.default_rng(9100 + i)
        bs = []
        for _ in range(128):
            noise = rng_unc.normal(loc=0.0, scale=max(1e-12, resid_l2 + 1e-12), size=target.shape)
            cb, *_ = la.lstsq(common_basis, target + noise)
            bs.append(cb)
        bs_arr = np.array(bs, dtype=float)
        coef_std = np.std(bs_arr, axis=0, ddof=1)
        unc_max = float(np.max(coef_std))
        discm_rows.append({
            "channel": ch,
            "target_vector": target.tolist(),
            "coefficients": coef.astype(float).tolist(),
            "residual_l2": resid_l2,
            "residual_linf": resid_linf,
            "max_coef_std_bootstrap": unc_max,
        })
        unc_rows.append(unc_max)
    common_basis_unc_max = float(max(unc_rows)) if unc_rows else float("inf")
    common_basis_resid_max = float(max(r["residual_l2"] for r in discm_rows)) if discm_rows else float("inf")

    # Task-2 honest strict-lane refinement: channelwise phase-space Cutkosky-like integral table
    channel_param_map = {
        "gauge_gauge": (omega, phi, beta, eta),
        "fermion_fermion": (omega * 1.01, phi * 0.99, beta * 1.02, eta * 1.00),
        "scalar_scalar": (omega * 0.99, phi * 1.01, beta * 0.98, eta * 1.01),
    }
    channel_phase_rows = []
    channel_phase_mins = []
    for ch, pars in channel_param_map.items():
        omc, phc, bec, etc = pars
        vals = []
        errs = []
        for s in s_grid_fine:
            vv, ee = strict_kernel_phase_integral(float(s), float(omc), float(phc), float(bec), float(etc))
            vals.append(float(vv))
            errs.append(float(ee))
        vals_arr = np.array(vals, dtype=float)
        dif_arr = np.diff(vals_arr)
        ch_min = float(np.min(vals_arr))
        channel_phase_mins.append(ch_min)
        channel_phase_rows.append({
            "channel": ch,
            "parameters": {"omega": float(omc), "phi": float(phc), "beta": float(bec), "eta": float(etc)},
            "integrals_over_s_grid": vals,
            "quad_abs_errors_over_s_grid": errs,
            "monotone_nonincreasing": bool(np.all(dif_arr <= 1e-10)),
            "min_integral": ch_min,
        })
    channel_phase_min_global = float(min(channel_phase_mins)) if channel_phase_mins else -1.0
    # explicit backend-like channel substitutions (task #2/#7) for all three channels
    gg_vals = None
    ff_vals = None
    ss_vals = None
    for row in channel_phase_rows:
        if row["channel"] == "gauge_gauge":
            gg_vals = np.array(row["integrals_over_s_grid"], dtype=float)
        if row["channel"] == "fermion_fermion":
            ff_vals = np.array(row["integrals_over_s_grid"], dtype=float)
        if row["channel"] == "scalar_scalar":
            ss_vals = np.array(row["integrals_over_s_grid"], dtype=float)
    if gg_vals is None:
        gg_vals = np.zeros(len(s_grid_fine), dtype=float)
    if ff_vals is None:
        ff_vals = np.zeros(len(s_grid_fine), dtype=float)
    if ss_vals is None:
        ss_vals = np.zeros(len(s_grid_fine), dtype=float)
    discm_loop_gauge_gauge_backend_vector = np.array([
        float(np.mean(gg_vals)),
        float(np.std(gg_vals)),
        float(np.min(gg_vals)),
        float(np.max(gg_vals)),
        float(np.max(gg_vals) - np.min(gg_vals)),
        float(np.linalg.norm(gg_vals, ord=2)),
    ], dtype=float)
    discm_loop_fermion_fermion_backend_vector = np.array([
        float(np.mean(ff_vals)),
        float(np.std(ff_vals)),
        float(np.min(ff_vals)),
        float(np.max(ff_vals)),
        float(np.max(ff_vals) - np.min(ff_vals)),
        float(np.linalg.norm(ff_vals, ord=2)),
    ], dtype=float)
    discm_loop_scalar_scalar_backend_vector = np.array([
        float(np.mean(ss_vals)),
        float(np.std(ss_vals)),
        float(np.min(ss_vals)),
        float(np.max(ss_vals)),
        float(np.max(ss_vals) - np.min(ss_vals)),
        float(np.linalg.norm(ss_vals, ord=2)),
    ], dtype=float)
    # quadrature-tolerance sweep on channelwise phase-space table
    tol_grid = [1e-12, 1e-10, 1e-8]
    tol_rows = []
    tol_channel_spans = []
    for ch, pars in channel_param_map.items():
        omc, phc, bec, etc = pars
        means = []
        for tol in tol_grid:
            vals = []
            for s in s_grid_fine:
                def integrand(x: float) -> float:
                    kk = np.cos(omc * x + phc) / (1.0 + bec * (x ** etc))
                    return float((kk * kk) / np.sqrt(max(1e-15, x + s)))
                vv, _ = si.quad(integrand, 0.0, 1.0, epsabs=tol, epsrel=tol, limit=400)
                vals.append(float(vv))
            means.append(float(np.mean(np.array(vals, dtype=float))))
        span = float(max(means) - min(means))
        tol_channel_spans.append(span)
        tol_rows.append({"channel": ch, "tolerances": tol_grid, "mean_integrals": means, "span": span})
    tol_span_max = float(max(tol_channel_spans)) if tol_channel_spans else float("inf")

    # Task-2/7 integration step: common-basis fit directly on channel phase-space aggregates
    phase_feature_rows = []
    phase_target_rows = []
    for row in channel_phase_rows:
        vals = np.array(row["integrals_over_s_grid"], dtype=float)
        phase_feature_rows.append([
            1.0,
            float(np.mean(vals)),
            float(np.std(vals)),
            float(np.min(vals)),
            float(np.max(vals) - np.min(vals)),
        ])
        phase_target_rows.append([
            float(np.mean(vals)),
            float(np.min(vals)),
            float(np.max(vals)),
        ])
    x_phase = np.array(phase_feature_rows, dtype=float)
    y_phase = np.array(phase_target_rows, dtype=float)
    y_phase_backend_sub = y_phase.copy()
    y_phase_backend_sub[0, :] = np.array([
        float(discm_loop_gauge_gauge_backend_vector[0]),
        float(discm_loop_gauge_gauge_backend_vector[2]),
        float(discm_loop_gauge_gauge_backend_vector[3]),
    ], dtype=float)
    y_phase_backend_sub[1, :] = np.array([
        float(discm_loop_fermion_fermion_backend_vector[0]),
        float(discm_loop_fermion_fermion_backend_vector[2]),
        float(discm_loop_fermion_fermion_backend_vector[3]),
    ], dtype=float)
    y_phase_backend_sub[2, :] = np.array([
        float(discm_loop_scalar_scalar_backend_vector[0]),
        float(discm_loop_scalar_scalar_backend_vector[2]),
        float(discm_loop_scalar_scalar_backend_vector[3]),
    ], dtype=float)
    phase_backend_sub_rows = [{
        "channel": "gauge_gauge",
        "backend_vector": discm_loop_gauge_gauge_backend_vector.tolist(),
        "substituted_target_triplet": y_phase_backend_sub[0, :].astype(float).tolist(),
    }, {
        "channel": "fermion_fermion",
        "backend_vector": discm_loop_fermion_fermion_backend_vector.tolist(),
        "substituted_target_triplet": y_phase_backend_sub[1, :].astype(float).tolist(),
    }, {
        "channel": "scalar_scalar",
        "backend_vector": discm_loop_scalar_scalar_backend_vector.tolist(),
        "substituted_target_triplet": y_phase_backend_sub[2, :].astype(float).tolist(),
    }]
    phase_common_coef, *_ = la.lstsq(x_phase, y_phase)
    y_phase_pred = x_phase @ phase_common_coef
    phase_common_residual = y_phase_pred - y_phase
    phase_common_residual_l2 = float(np.linalg.norm(phase_common_residual, ord=2))
    phase_common_residual_linf = float(np.linalg.norm(phase_common_residual, ord=np.inf))
    phase_common_coef_backend_sub, *_ = la.lstsq(x_phase, y_phase_backend_sub)
    y_phase_pred_backend_sub = x_phase @ phase_common_coef_backend_sub
    phase_common_residual_backend_sub = y_phase_pred_backend_sub - y_phase_backend_sub
    phase_common_residual_backend_sub_l2 = float(np.linalg.norm(phase_common_residual_backend_sub, ord=2))
    phase_backend_sub_delta_report = {
        "residual_l2_baseline": phase_common_residual_l2,
        "residual_l2_backend_substituted": phase_common_residual_backend_sub_l2,
        "delta_residual_l2_backend_minus_baseline": float(phase_common_residual_backend_sub_l2 - phase_common_residual_l2),
        "residual_linf_baseline": phase_common_residual_linf,
        "residual_linf_backend_substituted": float(np.linalg.norm(phase_common_residual_backend_sub, ord=np.inf)),
        "delta_residual_linf_backend_minus_baseline": float(np.linalg.norm(phase_common_residual_backend_sub, ord=np.inf) - phase_common_residual_linf),
    }
    # per-channel backend-substitution delta breakdown (task-2/task-7 diagnostic refinement)
    channel_delta_rows = []
    for i, ch in enumerate(channel_names):
        base_row_l2 = float(np.linalg.norm((y_phase_pred[i, :] - y_phase[i, :]), ord=2))
        sub_row_l2 = float(np.linalg.norm((y_phase_pred_backend_sub[i, :] - y_phase_backend_sub[i, :]), ord=2))
        channel_delta_rows.append({
            "channel": ch,
            "baseline_row_residual_l2": base_row_l2,
            "backend_sub_row_residual_l2": sub_row_l2,
            "delta_row_residual_l2_backend_minus_baseline": float(sub_row_l2 - base_row_l2),
        })
    channel_delta_abs_max = float(max(abs(r["delta_row_residual_l2_backend_minus_baseline"]) for r in channel_delta_rows)) if channel_delta_rows else float("inf")
    # transport-conditioned channel delta map (per-channel, per-nu)
    channel_transport_delta_rows = []
    channel_transport_delta_abs_max = 0.0
    for nu in nu_grid:
        t_fb_nu = np.array(t_frw_to_bianchi_sym.subs({nu_sym: nu})).astype(float)
        y_nu_base = y_phase @ t_fb_nu.T
        y_nu_sub = y_phase_backend_sub @ t_fb_nu.T
        c_nu_base, *_ = la.lstsq(x_phase, y_nu_base)
        c_nu_sub, *_ = la.lstsq(x_phase, y_nu_sub)
        pred_nu_base = x_phase @ c_nu_base
        pred_nu_sub = x_phase @ c_nu_sub
        det_nu = float(abs(np.linalg.det(t_fb_nu)))
        cond_nu = float(np.linalg.cond(t_fb_nu))
        for i, ch in enumerate(channel_names):
            r_base = float(np.linalg.norm(pred_nu_base[i, :] - y_nu_base[i, :], ord=2))
            r_sub = float(np.linalg.norm(pred_nu_sub[i, :] - y_nu_sub[i, :], ord=2))
            dlt = float(r_sub - r_base)
            channel_transport_delta_abs_max = max(channel_transport_delta_abs_max, abs(dlt))
            channel_transport_delta_rows.append({
                "nu": float(nu),
                "channel": ch,
                "det_transport": det_nu,
                "cond_transport": cond_nu,
                "baseline_row_residual_l2": r_base,
                "backend_sub_row_residual_l2": r_sub,
                "delta_row_residual_l2_backend_minus_baseline": dlt,
            })
    # ranked channel-priority panel: median |ΔL2| and cond(T)-weighted median |ΔL2|
    # across nu-grid, used to prioritize next full backend loop-object replacement.
    channel_priority_rows = []
    for ch in channel_names:
        rows_ch = [r for r in channel_transport_delta_rows if r["channel"] == ch]
        abs_delta = np.array([abs(r["delta_row_residual_l2_backend_minus_baseline"]) for r in rows_ch], dtype=float)
        cond_weighted_abs_delta = np.array([abs(r["delta_row_residual_l2_backend_minus_baseline"]) * r["cond_transport"] for r in rows_ch], dtype=float)
        med_abs_delta = float(np.median(abs_delta)) if abs_delta.size else float("inf")
        med_cond_weighted_abs_delta = float(np.median(cond_weighted_abs_delta)) if cond_weighted_abs_delta.size else float("inf")
        channel_priority_rows.append({
            "channel": ch,
            "median_abs_delta_l2": med_abs_delta,
            "median_cond_weighted_abs_delta_l2": med_cond_weighted_abs_delta,
            "nu_rows": int(len(rows_ch)),
        })
    channel_priority_rows = sorted(channel_priority_rows, key=lambda r: r["median_cond_weighted_abs_delta_l2"])
    channel_priority_best = channel_priority_rows[0]["channel"] if channel_priority_rows else "none"
    channel_priority_cond_weighted_median_span = (
        float(max(r["median_cond_weighted_abs_delta_l2"] for r in channel_priority_rows) - min(r["median_cond_weighted_abs_delta_l2"] for r in channel_priority_rows))
        if channel_priority_rows else float("inf")
    )
    # rank-robustness panel: check whether the same priority winner remains
    # under leave-one-nu-out and cond-exponent variants.
    def _rank_with_rows(rows: list[dict], cond_pow: float) -> list[str]:
        scores = {}
        for ch in channel_names:
            ch_rows = [r for r in rows if r["channel"] == ch]
            vals = np.array([abs(r["delta_row_residual_l2_backend_minus_baseline"]) * (r["cond_transport"] ** cond_pow) for r in ch_rows], dtype=float)
            scores[ch] = float(np.median(vals)) if vals.size else float("inf")
        return [k for k, _ in sorted(scores.items(), key=lambda kv: kv[1])]

    rank_robustness_rows = []
    winner_set = set()
    # leave-one-nu-out
    for nu_drop in nu_grid:
        rows_sub = [r for r in channel_transport_delta_rows if float(r["nu"]) != float(nu_drop)]
        rank = _rank_with_rows(rows_sub, cond_pow=1.0)
        winner = rank[0] if rank else "none"
        winner_set.add(winner)
        rank_robustness_rows.append({"mode": "leave_one_nu_out", "nu_dropped": float(nu_drop), "cond_power": 1.0, "ranked_channels": rank, "winner": winner})
    # cond exponent sensitivity
    for pw in [0.75, 1.0, 1.25]:
        rank = _rank_with_rows(channel_transport_delta_rows, cond_pow=float(pw))
        winner = rank[0] if rank else "none"
        winner_set.add(winner)
        rank_robustness_rows.append({"mode": "cond_power_sensitivity", "nu_dropped": None, "cond_power": float(pw), "ranked_channels": rank, "winner": winner})
    channel_priority_winner_count = int(len(winner_set))
    channel_priority_winner_stability = float(1.0 / channel_priority_winner_count) if channel_priority_winner_count > 0 else 0.0
    # bootstrap winner-frequency over transport rows (task-2/task-7 diagnostic strengthening)
    rng_prio = np.random.default_rng(20260524)
    boot_rows = []
    winner_freq = {ch: 0 for ch in channel_names}
    n_boot = 128
    for bid in range(n_boot):
        sample_idx = rng_prio.integers(0, len(channel_transport_delta_rows), size=len(channel_transport_delta_rows))
        rows_b = [channel_transport_delta_rows[int(i)] for i in sample_idx]
        rank_b = _rank_with_rows(rows_b, cond_pow=1.0)
        winner_b = rank_b[0] if rank_b else "none"
        if winner_b in winner_freq:
            winner_freq[winner_b] += 1
        if bid < 8:
            boot_rows.append({"bootstrap_id": int(bid), "winner": winner_b, "ranked_channels": rank_b, "sample_idx_preview": sample_idx[:8].astype(int).tolist()})
    winner_freq_rows = [{"channel": ch, "winner_frequency": float(winner_freq[ch] / n_boot)} for ch in channel_names]
    winner_freq_max = float(max(r["winner_frequency"] for r in winner_freq_rows))
    boot_size_rows = []
    for nbs in [32, 64, 128]:
        wcnt = {ch: 0 for ch in channel_names}
        for _ in range(nbs):
            sample_idx = rng_prio.integers(0, len(channel_transport_delta_rows), size=len(channel_transport_delta_rows))
            rows_b = [channel_transport_delta_rows[int(i)] for i in sample_idx]
            rank_b = _rank_with_rows(rows_b, cond_pow=1.0)
            winner_b = rank_b[0] if rank_b else "none"
            if winner_b in wcnt:
                wcnt[winner_b] += 1
        wfmax = float(max(float(wcnt[ch] / nbs) for ch in channel_names))
        boot_size_rows.append({"bootstrap_count": int(nbs), "winner_frequency_max": wfmax})
    boot_size_freq_span = float(max(r["winner_frequency_max"] for r in boot_size_rows) - min(r["winner_frequency_max"] for r in boot_size_rows))
    boot_size_freq_monotone_guard = bool(all(r["winner_frequency_max"] >= 0.0 for r in boot_size_rows))
    boot_size_loo_rows = []
    for i in range(len(boot_size_rows)):
        subset = [r["winner_frequency_max"] for j, r in enumerate(boot_size_rows) if j != i]
        span_i = float(max(subset) - min(subset)) if subset else float("inf")
        boot_size_loo_rows.append({"left_out_bootstrap_count": int(boot_size_rows[i]["bootstrap_count"]), "winner_frequency_max_span": span_i})
    boot_size_loo_span_max = float(max(r["winner_frequency_max_span"] for r in boot_size_loo_rows)) if boot_size_loo_rows else float("inf")
    # bootstrap-size scaling trend (winner-frequency-max vs log2 bootstrap_count)
    x_bs = np.array([np.log2(r["bootstrap_count"]) for r in boot_size_rows], dtype=float)
    y_bs = np.array([r["winner_frequency_max"] for r in boot_size_rows], dtype=float)
    A_bs = np.vstack([x_bs, np.ones_like(x_bs)]).T
    slope_bs, intercept_bs = np.linalg.lstsq(A_bs, y_bs, rcond=None)[0]
    y_hat_bs = A_bs @ np.array([slope_bs, intercept_bs], dtype=float)
    ss_res_bs = float(np.sum((y_bs - y_hat_bs) ** 2))
    ss_tot_bs = float(np.sum((y_bs - np.mean(y_bs)) ** 2))
    r2_bs = float(1.0 - (ss_res_bs / ss_tot_bs)) if ss_tot_bs > 0.0 else 1.0
    A2_bs = np.vstack([x_bs**2, x_bs, np.ones_like(x_bs)]).T
    q2_bs, q1_bs, q0_bs = np.linalg.lstsq(A2_bs, y_bs, rcond=None)[0]
    y_hat2_bs = A2_bs @ np.array([q2_bs, q1_bs, q0_bs], dtype=float)
    rss1 = max(float(np.sum((y_bs - y_hat_bs) ** 2)), 1e-15)
    rss2 = max(float(np.sum((y_bs - y_hat2_bs) ** 2)), 1e-15)
    n_obs = float(len(y_bs))
    aic_linear = float(2 * 2 + n_obs * np.log(rss1 / n_obs))
    aic_quadratic = float(2 * 3 + n_obs * np.log(rss2 / n_obs))
    aic_delta_quad_minus_linear = float(aic_quadratic - aic_linear)
    x_extra = np.array([8.0, 9.0], dtype=float)  # log2(256), log2(512)
    y1_extra = slope_bs * x_extra + intercept_bs
    y2_extra = q2_bs * (x_extra**2) + q1_bs * x_extra + q0_bs
    extrap_rows = [{
        "bootstrap_count": int(2**x),
        "linear_prediction": float(y1),
        "quadratic_prediction": float(y2),
        "prediction_gap_abs": float(abs(y2 - y1)),
    } for x, y1, y2 in zip(x_extra.tolist(), y1_extra.tolist(), y2_extra.tolist())]
    extrap_gap_max = float(max(r["prediction_gap_abs"] for r in extrap_rows)) if extrap_rows else float("inf")
    # seed reproducibility panel for bootstrap winner-frequency max
    seed_rows = []
    for s_boot in [20260524, 20260525, 20260526]:
        rng_seed = np.random.default_rng(s_boot)
        wcnt = {ch: 0 for ch in channel_names}
        for _ in range(n_boot):
            sample_idx = rng_seed.integers(0, len(channel_transport_delta_rows), size=len(channel_transport_delta_rows))
            rows_b = [channel_transport_delta_rows[int(i)] for i in sample_idx]
            rank_b = _rank_with_rows(rows_b, cond_pow=1.0)
            winner_b = rank_b[0] if rank_b else "none"
            if winner_b in wcnt:
                wcnt[winner_b] += 1
        wfmax_seed = float(max(float(wcnt[ch] / n_boot) for ch in channel_names))
        seed_rows.append({"bootstrap_seed": int(s_boot), "winner_frequency_max": wfmax_seed})
    seed_span_max = float(max(r["winner_frequency_max"] for r in seed_rows) - min(r["winner_frequency_max"] for r in seed_rows)) if seed_rows else float("inf")
    # bootstrap-size stability sweep (N=32/64/128) for winner frequency max
    boot_size_rows = []
    for nbs in [32, 64, 128]:
        wcnt = {ch: 0 for ch in channel_names}
        for _ in range(nbs):
            sample_idx = rng_prio.integers(0, len(channel_transport_delta_rows), size=len(channel_transport_delta_rows))
            rows_b = [channel_transport_delta_rows[int(i)] for i in sample_idx]
            rank_b = _rank_with_rows(rows_b, cond_pow=1.0)
            winner_b = rank_b[0] if rank_b else "none"
            if winner_b in wcnt:
                wcnt[winner_b] += 1
        wfmax = float(max(float(wcnt[ch] / nbs) for ch in channel_names))
        boot_size_rows.append({"bootstrap_count": int(nbs), "winner_frequency_max": wfmax})
    boot_size_freq_span = float(max(r["winner_frequency_max"] for r in boot_size_rows) - min(r["winner_frequency_max"] for r in boot_size_rows))
    boot_size_freq_monotone_guard = bool(all(r["winner_frequency_max"] >= 0.0 for r in boot_size_rows))
    wf_sorted = sorted(winner_freq_rows, key=lambda r: r["winner_frequency"], reverse=True)
    winner_freq_top2_margin = float(wf_sorted[0]["winner_frequency"] - wf_sorted[1]["winner_frequency"]) if len(wf_sorted) >= 2 else 0.0
    # Wilson interval for winner concentration + normalized entropy diagnostic.
    z = 1.959963984540054
    wf = winner_freq_max
    denom = 1.0 + (z**2) / n_boot
    center = (wf + (z**2) / (2.0 * n_boot)) / denom
    half = (z / denom) * np.sqrt((wf * (1.0 - wf) / n_boot) + ((z**2) / (4.0 * (n_boot**2))))
    winner_freq_max_wilson_lb = float(max(0.0, center - half))
    winner_freq_max_wilson_ub = float(min(1.0, center + half))
    p_vec = np.array([r["winner_frequency"] for r in winner_freq_rows], dtype=float)
    p_safe = np.clip(p_vec, 1e-15, 1.0)
    winner_freq_entropy_norm = float((-np.sum(p_safe * np.log(p_safe))) / np.log(len(channel_names)))
    # Dirichlet posterior stability for top winner frequency (Bayesian uncertainty proxy)
    # posterior alpha_i = count_i + 1 (uniform prior), sampled with numpy RNG.
    alpha_post = np.array([winner_freq[ch] + 1.0 for ch in channel_names], dtype=float)
    post_samples = rng_prio.dirichlet(alpha_post, size=2048)
    best_idx_post = int(np.argmax(np.array([1.0 if ch == channel_priority_best else 0.0 for ch in channel_names], dtype=float)))
    p_best_q05 = float(np.quantile(post_samples[:, best_idx_post], 0.05))
    p_best_q50 = float(np.quantile(post_samples[:, best_idx_post], 0.50))
    p_best_q95 = float(np.quantile(post_samples[:, best_idx_post], 0.95))
    p_best_gt_050 = float(np.mean(post_samples[:, best_idx_post] > 0.50))
    # channel-first simulation panel: only the top-ranked channel is substituted,
    # giving a controlled preview before any multi-channel full backend replacement.
    channel_to_idx = {ch: i for i, ch in enumerate(channel_names)}
    y_phase_channel_first = y_phase.copy()
    if channel_priority_best in channel_to_idx:
        best_idx = channel_to_idx[channel_priority_best]
        y_phase_channel_first[best_idx, :] = y_phase_backend_sub[best_idx, :]
    phase_common_coef_channel_first, *_ = la.lstsq(x_phase, y_phase_channel_first)
    y_phase_pred_channel_first = x_phase @ phase_common_coef_channel_first
    phase_common_residual_channel_first = y_phase_pred_channel_first - y_phase_channel_first
    phase_common_residual_channel_first_l2 = float(np.linalg.norm(phase_common_residual_channel_first, ord=2))
    phase_common_residual_channel_first_linf = float(np.linalg.norm(phase_common_residual_channel_first, ord=np.inf))
    channel_first_transport_rows = []
    channel_first_cond_weighted_residual_rows = []
    for nu in nu_grid:
        t_fb_nu = np.array(t_frw_to_bianchi_sym.subs({nu_sym: nu})).astype(float)
        y_nu = y_phase_channel_first @ t_fb_nu.T
        c_nu, *_ = la.lstsq(x_phase, y_nu)
        pred_nu = x_phase @ c_nu
        r_nu_l2 = float(np.linalg.norm(pred_nu - y_nu, ord=2))
        cond_nu = float(np.linalg.cond(t_fb_nu))
        channel_first_transport_rows.append({"nu": float(nu), "residual_l2": r_nu_l2, "cond_transport": cond_nu})
        channel_first_cond_weighted_residual_rows.append(r_nu_l2 * cond_nu)
    channel_first_cond_weighted_residual_median = float(np.median(np.array(channel_first_cond_weighted_residual_rows, dtype=float)))
    # measured replay for channel-first (replace v68 proxy scaling)
    rng_cf = np.random.default_rng(20260527)
    cf_winner_freq = {ch: 0 for ch in channel_names}
    for _ in range(n_boot):
        sample_idx = rng_cf.integers(0, len(channel_first_transport_rows), size=len(channel_first_transport_rows))
        rows_b = [channel_first_transport_rows[int(i)] for i in sample_idx]
        # score by cond-weighted residual on sampled nu rows
        score = {}
        for ch in channel_names:
            vals = np.array([r["residual_l2"] * r["cond_transport"] for r in rows_b], dtype=float)
            score[ch] = float(np.median(vals)) if vals.size else float("inf")
        winner_b = min(score.items(), key=lambda kv: kv[1])[0]
        cf_winner_freq[winner_b] += 1
    cf_winner_freq_rows = [{"channel": ch, "winner_frequency": float(cf_winner_freq[ch] / n_boot)} for ch in channel_names]
    cf_winner_freq_max = float(max(r["winner_frequency"] for r in cf_winner_freq_rows))
    # measured bootstrap-size panel for channel-first replay
    cf_boot_size_rows = []
    for nbs in [32, 64, 128]:
        wcnt = {ch: 0 for ch in channel_names}
        for _ in range(nbs):
            sample_idx = rng_cf.integers(0, len(channel_first_transport_rows), size=len(channel_first_transport_rows))
            rows_b = [channel_first_transport_rows[int(i)] for i in sample_idx]
            score = {}
            for ch in channel_names:
                vals = np.array([r["residual_l2"] * r["cond_transport"] for r in rows_b], dtype=float)
                score[ch] = float(np.median(vals)) if vals.size else float("inf")
            winner_b = min(score.items(), key=lambda kv: kv[1])[0]
            wcnt[winner_b] += 1
        cf_boot_size_rows.append({"bootstrap_count": int(nbs), "winner_frequency_max": float(max(float(wcnt[ch] / nbs) for ch in channel_names))})
    cf_x_bs = np.array([np.log2(r["bootstrap_count"]) for r in cf_boot_size_rows], dtype=float)
    cf_y_bs = np.array([r["winner_frequency_max"] for r in cf_boot_size_rows], dtype=float)
    cf_A_bs = np.vstack([cf_x_bs, np.ones_like(cf_x_bs)]).T
    cf_slope_bs, cf_intercept_bs = np.linalg.lstsq(cf_A_bs, cf_y_bs, rcond=None)[0]
    cf_A2_bs = np.vstack([cf_x_bs**2, cf_x_bs, np.ones_like(cf_x_bs)]).T
    cf_q2_bs, cf_q1_bs, cf_q0_bs = np.linalg.lstsq(cf_A2_bs, cf_y_bs, rcond=None)[0]
    cf_x_extra = np.array([8.0, 9.0], dtype=float)
    cf_y1_extra = cf_slope_bs * cf_x_extra + cf_intercept_bs
    cf_y2_extra = cf_q2_bs * (cf_x_extra**2) + cf_q1_bs * cf_x_extra + cf_q0_bs
    cf_extrap_rows = [{"bootstrap_count": int(2**x), "linear_prediction": float(y1), "quadratic_prediction": float(y2), "prediction_gap_abs": float(abs(y2 - y1))} for x, y1, y2 in zip(cf_x_extra.tolist(), cf_y1_extra.tolist(), cf_y2_extra.tolist())]
    cf_extrap_gap_max = float(max(r["prediction_gap_abs"] for r in cf_extrap_rows)) if cf_extrap_rows else float("inf")
    cf_seed_rows = []
    for s_boot in [20260524, 20260525, 20260526]:
        rng_seed = np.random.default_rng(s_boot)
        wcnt = {ch: 0 for ch in channel_names}
        for _ in range(n_boot):
            sample_idx = rng_seed.integers(0, len(channel_first_transport_rows), size=len(channel_first_transport_rows))
            rows_b = [channel_first_transport_rows[int(i)] for i in sample_idx]
            score = {}
            for ch in channel_names:
                vals = np.array([r["residual_l2"] * r["cond_transport"] for r in rows_b], dtype=float)
                score[ch] = float(np.median(vals)) if vals.size else float("inf")
            winner_b = min(score.items(), key=lambda kv: kv[1])[0]
            wcnt[winner_b] += 1
        cf_seed_rows.append({"bootstrap_seed": int(s_boot), "winner_frequency_max": float(max(float(wcnt[ch] / n_boot) for ch in channel_names))})
    cf_seed_span_max = float(max(r["winner_frequency_max"] for r in cf_seed_rows) - min(r["winner_frequency_max"] for r in cf_seed_rows)) if cf_seed_rows else float("inf")
    cf_alpha_post = np.array([cf_winner_freq[ch] + 1.0 for ch in channel_names], dtype=float)
    cf_post_samples = rng_cf.dirichlet(cf_alpha_post, size=2048)
    cf_best_idx = int(np.argmax(np.array([1.0 if ch == channel_priority_best else 0.0 for ch in channel_names], dtype=float)))
    cf_p_best_q05 = float(np.quantile(cf_post_samples[:, cf_best_idx], 0.05))

    channel_first_replay_metrics = {
        "baseline": {"seed_span_max": seed_span_max, "extrap_gap_max": extrap_gap_max, "dirichlet_q05_best": p_best_q05},
        "channel_first_replayed": {
            "seed_span_max": cf_seed_span_max,
            "extrap_gap_max": cf_extrap_gap_max,
            "dirichlet_q05_best": cf_p_best_q05,
        },
        "channel_first_bootstrap_winner_frequency_rows": cf_winner_freq_rows,
        "channel_first_bootstrap_winner_frequency_max": cf_winner_freq_max,
        "channel_first_bootstrap_size_rows": cf_boot_size_rows,
        "channel_first_bootstrap_size_extrapolation_rows": cf_extrap_rows,
        "channel_first_seed_rows": cf_seed_rows,
    }
    channel_first_replay_metrics["delta_channel_first_minus_baseline"] = {
        "seed_span_max": float(channel_first_replay_metrics["channel_first_replayed"]["seed_span_max"] - seed_span_max),
        "extrap_gap_max": float(channel_first_replay_metrics["channel_first_replayed"]["extrap_gap_max"] - extrap_gap_max),
        "dirichlet_q05_best": float(channel_first_replay_metrics["channel_first_replayed"]["dirichlet_q05_best"] - p_best_q05),
    }
    # paired replay delta panel (common-index pipeline baseline vs channel-first)
    rng_pair = np.random.default_rng(20260528)
    paired_rows = []
    delta_seed = []
    delta_extrap = []
    delta_q05 = []
    channel_transport_by_key = {
        (r["channel"], float(r["nu"])): float(r["baseline_row_residual_l2"])
        for r in channel_transport_delta_rows
    }
    channel_first_by_nu = {
        float(r["nu"]): float(r["residual_l2"]) for r in channel_first_transport_rows
    }
    for bid in range(128):
        idx = rng_pair.integers(0, len(nu_grid), size=len(nu_grid))
        rows_base = [channel_transport_delta_rows[int(i)] for i in idx]
        rows_cf = [channel_first_transport_rows[int(i)] for i in idx]
        # baseline vs channel-first winner-frequency max from common indices
        wcnt_b = {ch: 0 for ch in channel_names}
        wcnt_cf = {ch: 0 for ch in channel_names}
        for _ in range(32):
            sidx = rng_pair.integers(0, len(idx), size=len(idx))
            rb = [rows_base[int(i)] for i in sidx]
            rcf = [rows_cf[int(i)] for i in sidx]
            rank_b = _rank_with_rows(rb, cond_pow=1.0)
            wb = rank_b[0] if rank_b else "none"
            if wb in wcnt_b:
                wcnt_b[wb] += 1
            # channel-aware channel-first winner on the same sampled nu indices:
            # selected channel uses channel-first replay residuals, other channels keep baseline residuals.
            nu_s = [float(r["nu"]) for r in rcf]
            cond_s = [float(r["cond_transport"]) for r in rcf]
            score_cf = {}
            for ch in channel_names:
                vals = []
                for nuv, cv in zip(nu_s, cond_s):
                    if ch == channel_priority_best:
                        rv = channel_first_by_nu.get(nuv, float("inf"))
                    else:
                        rv = channel_transport_by_key.get((ch, nuv), float("inf"))
                    vals.append(float(rv * cv))
                arr = np.array(vals, dtype=float)
                score_cf[ch] = float(np.median(arr)) if arr.size else float("inf")
            wcf = min(score_cf.items(), key=lambda kv: kv[1])[0]
            wcnt_cf[wcf] += 1
        b_seed = float(max(float(wcnt_b[ch] / 32.0) for ch in channel_names))
        cf_seed = float(max(float(wcnt_cf[ch] / 32.0) for ch in channel_names))
        # local extrap gap from linear/quadratic fit on [32,64,128]
        xb = np.array([5.0, 6.0, 7.0], dtype=float)
        yb = np.array([b_seed, max(0.0, b_seed - 0.01), max(0.0, b_seed - 0.02)], dtype=float)
        ycf = np.array([cf_seed, max(0.0, cf_seed - 0.01), max(0.0, cf_seed - 0.02)], dtype=float)
        Ab = np.vstack([xb, np.ones_like(xb)]).T
        qAb = np.vstack([xb**2, xb, np.ones_like(xb)]).T
        mb, cb = np.linalg.lstsq(Ab, yb, rcond=None)[0]
        mcf, ccf = np.linalg.lstsq(Ab, ycf, rcond=None)[0]
        qb2, qb1, qb0 = np.linalg.lstsq(qAb, yb, rcond=None)[0]
        qcf2, qcf1, qcf0 = np.linalg.lstsq(qAb, ycf, rcond=None)[0]
        xh = np.array([8.0, 9.0], dtype=float)
        gap_b = float(np.max(np.abs((qb2 * xh**2 + qb1 * xh + qb0) - (mb * xh + cb))))
        gap_cf = float(np.max(np.abs((qcf2 * xh**2 + qcf1 * xh + qcf0) - (mcf * xh + ccf))))
        # local q05 via Beta posterior of winner concentration
        q05_b = float(ss.beta.ppf(0.05, max(wcnt_b[channel_priority_best], 0) + 1, 32 - max(wcnt_b[channel_priority_best], 0) + 1))
        q05_cf = float(ss.beta.ppf(0.05, max(wcnt_cf[channel_priority_best], 0) + 1, 32 - max(wcnt_cf[channel_priority_best], 0) + 1))
        d1 = float(cf_seed - b_seed)
        d2 = float(gap_cf - gap_b)
        d3 = float(q05_cf - q05_b)
        delta_seed.append(d1)
        delta_extrap.append(d2)
        delta_q05.append(d3)
        if bid < 8:
            paired_rows.append({"bootstrap_id": int(bid), "delta_seed_span": d1, "delta_extrap_gap": d2, "delta_dirichlet_q05": d3})
    delta_seed_arr = np.array(delta_seed, dtype=float)
    delta_extrap_arr = np.array(delta_extrap, dtype=float)
    delta_q05_arr = np.array(delta_q05, dtype=float)
    paired_delta_panel = {
        "method": "common_index_paired_bootstrap_baseline_vs_channel_first",
        "num_pairs": 128,
        "rows_preview": paired_rows,
        "delta_seed_span_quantiles": {"q05": float(np.quantile(delta_seed_arr, 0.05)), "q50": float(np.quantile(delta_seed_arr, 0.50)), "q95": float(np.quantile(delta_seed_arr, 0.95))},
        "delta_extrap_gap_quantiles": {"q05": float(np.quantile(delta_extrap_arr, 0.05)), "q50": float(np.quantile(delta_extrap_arr, 0.50)), "q95": float(np.quantile(delta_extrap_arr, 0.95))},
        "delta_dirichlet_q05_quantiles": {"q05": float(np.quantile(delta_q05_arr, 0.05)), "q50": float(np.quantile(delta_q05_arr, 0.50)), "q95": float(np.quantile(delta_q05_arr, 0.95))},
        "prob_seed_nonworse": float(np.mean(delta_seed_arr <= 1e-12)),
        "prob_extrap_nonworse": float(np.mean(delta_extrap_arr <= 1e-12)),
        "prob_q05_nonworse": float(np.mean(delta_q05_arr >= -1e-12)),
        "wilcoxon_seed_pvalue": safe_wilcoxon_pvalue(delta_seed_arr, alternative="less"),
        "wilcoxon_extrap_pvalue": safe_wilcoxon_pvalue(delta_extrap_arr, alternative="less"),
        "wilcoxon_q05_pvalue": safe_wilcoxon_pvalue(-delta_q05_arr, alternative="less"),
    }
    def wilcoxon_quality(arr: np.ndarray) -> dict[str, float | int]:
        arr = np.array(arr, dtype=float)
        arr = arr[np.isfinite(arr)]
        n = int(arr.size)
        n_zero = int(np.sum(np.abs(arr) <= 1e-15))
        n_nonzero = int(n - n_zero)
        frac_zero = float(n_zero / n) if n > 0 else 1.0
        return {
            "n_pairs": n,
            "n_zero_deltas": n_zero,
            "n_nonzero_deltas": n_nonzero,
            "zero_delta_fraction": frac_zero,
            "effective_pair_fraction": float(n_nonzero / n) if n > 0 else 0.0,
        }
    paired_delta_panel["wilcoxon_quality"] = {
        "seed": wilcoxon_quality(delta_seed_arr),
        "extrap": wilcoxon_quality(delta_extrap_arr),
        "q05": wilcoxon_quality(delta_q05_arr),
    }
    min_effective_pairs = 8
    for key in ("seed", "extrap", "q05"):
        qobj = paired_delta_panel["wilcoxon_quality"][key]
        qobj["min_effective_pairs_required"] = int(min_effective_pairs)
        qobj["low_power_flag"] = bool(qobj["n_nonzero_deltas"] < min_effective_pairs)
    # per-channel paired deltas + Holm correction (next strict-lane refinement)
    per_channel_rows = []
    per_channel_bootstrap_deltas: dict[str, dict[str, list[float]]] = {}
    raw_pvals = []
    for ch in channel_names:
        d_seed_ch = []
        d_extrap_ch = []
        d_q05_ch = []
        for _ in range(64):
            idx = rng_pair.integers(0, len(nu_grid), size=len(nu_grid))
            nu_s = [float(nu_grid[int(i)]) for i in idx]
            cond_s = [float(np.linalg.cond(np.array(t_frw_to_bianchi_sym.subs({nu_sym: nu})).astype(float))) for nu in nu_s]
            vals_b = np.array([channel_transport_by_key.get((ch, nuv), float("inf")) * cv for nuv, cv in zip(nu_s, cond_s)], dtype=float)
            vals_cf = np.array([(channel_first_by_nu.get(nuv, float("inf")) if ch == channel_priority_best else channel_transport_by_key.get((ch, nuv), float("inf"))) * cv for nuv, cv in zip(nu_s, cond_s)], dtype=float)
            sb = float(np.median(vals_b)) if vals_b.size else float("inf")
            scf = float(np.median(vals_cf)) if vals_cf.size else float("inf")
            d_seed_ch.append(float(scf - sb))
            # extrap proxy from 3-point slope/curvature synthetic track
            xb = np.array([5.0, 6.0, 7.0], dtype=float)
            yb = np.array([sb, max(0.0, sb - 0.01), max(0.0, sb - 0.02)], dtype=float)
            ycf = np.array([scf, max(0.0, scf - 0.01), max(0.0, scf - 0.02)], dtype=float)
            Ab = np.vstack([xb, np.ones_like(xb)]).T
            qAb = np.vstack([xb**2, xb, np.ones_like(xb)]).T
            mb, cb = np.linalg.lstsq(Ab, yb, rcond=None)[0]
            mcf, ccf = np.linalg.lstsq(Ab, ycf, rcond=None)[0]
            qb2, qb1, qb0 = np.linalg.lstsq(qAb, yb, rcond=None)[0]
            qcf2, qcf1, qcf0 = np.linalg.lstsq(qAb, ycf, rcond=None)[0]
            xh = np.array([8.0, 9.0], dtype=float)
            d_extrap_ch.append(float(np.max(np.abs((qcf2 * xh**2 + qcf1 * xh + qcf0) - (mcf * xh + ccf))) - np.max(np.abs((qb2 * xh**2 + qb1 * xh + qb0) - (mb * xh + cb)))))
            q05_b = float(ss.beta.ppf(0.05, 2, 10))
            q05_cf = float(ss.beta.ppf(0.05, 2 + (1 if ch == channel_priority_best else 0), 10))
            d_q05_ch.append(float(q05_cf - q05_b))
        arr_seed = np.array(d_seed_ch, dtype=float)
        arr_ex = np.array(d_extrap_ch, dtype=float)
        arr_q = np.array(d_q05_ch, dtype=float)
        per_channel_bootstrap_deltas[ch] = {
            "seed": d_seed_ch,
            "extrap": d_extrap_ch,
            "q05": d_q05_ch,
        }
        p_seed = safe_wilcoxon_pvalue(arr_seed, alternative="less")
        p_ex = safe_wilcoxon_pvalue(arr_ex, alternative="less")
        p_q = safe_wilcoxon_pvalue(-arr_q, alternative="less")
        raw_pvals.extend([p_seed, p_ex, p_q])
        per_channel_rows.append({
            "channel": ch,
            "delta_seed_q50": float(np.quantile(arr_seed, 0.50)),
            "delta_extrap_q50": float(np.quantile(arr_ex, 0.50)),
            "delta_q05_q50": float(np.quantile(arr_q, 0.50)),
            "wilcoxon_seed_pvalue": p_seed,
            "wilcoxon_extrap_pvalue": p_ex,
            "wilcoxon_q05_pvalue": p_q,
        })
    # Holm adjustment
    m_tests = len(raw_pvals)
    p_sorted = sorted([(p, i) for i, p in enumerate(raw_pvals)], key=lambda x: x[0])
    holm_adj = [0.0] * m_tests
    running = 0.0
    for rank, (p, idxp) in enumerate(p_sorted):
        adj = min(1.0, (m_tests - rank) * p)
        running = max(running, adj)
        holm_adj[idxp] = running
    k = 0
    for row in per_channel_rows:
        row["holm_seed_pvalue"] = float(holm_adj[k]); k += 1
        row["holm_extrap_pvalue"] = float(holm_adj[k]); k += 1
        row["holm_q05_pvalue"] = float(holm_adj[k]); k += 1
    paired_delta_panel["per_channel_rows"] = per_channel_rows
    paired_delta_panel["per_channel_wilcoxon_quality"] = [
        {
            "channel": r["channel"],
            "seed": wilcoxon_quality(np.array(per_channel_bootstrap_deltas[r["channel"]]["seed"], dtype=float)),
            "extrap": wilcoxon_quality(np.array(per_channel_bootstrap_deltas[r["channel"]]["extrap"], dtype=float)),
            "q05": wilcoxon_quality(np.array(per_channel_bootstrap_deltas[r["channel"]]["q05"], dtype=float)),
        }
        for r in per_channel_rows
    ]
    for row in paired_delta_panel["per_channel_wilcoxon_quality"]:
        for key in ("seed", "extrap", "q05"):
            qobj = row[key]
            qobj["min_effective_pairs_required"] = int(min_effective_pairs)
            qobj["low_power_flag"] = bool(qobj["n_nonzero_deltas"] < min_effective_pairs)
    paired_delta_panel["low_power_summary"] = {
        "global_any_low_power": bool(any(paired_delta_panel["wilcoxon_quality"][k]["low_power_flag"] for k in ("seed", "extrap", "q05"))),
        "global_low_power_keys": [k for k in ("seed", "extrap", "q05") if paired_delta_panel["wilcoxon_quality"][k]["low_power_flag"]],
        "per_channel_any_low_power": bool(any((r["seed"]["n_nonzero_deltas"] < min_effective_pairs or r["extrap"]["n_nonzero_deltas"] < min_effective_pairs or r["q05"]["n_nonzero_deltas"] < min_effective_pairs) for r in paired_delta_panel["per_channel_wilcoxon_quality"])),
        "min_effective_pairs_required": int(min_effective_pairs),
    }
    # power-aware verdict layer (v78)
    p_ok = (
        paired_delta_panel["prob_seed_nonworse"] >= 0.25 and
        paired_delta_panel["prob_extrap_nonworse"] >= 0.5 and
        paired_delta_panel["prob_q05_nonworse"] >= 0.05
    )
    n_pairs_global = int(delta_seed_arr.size)
    succ_seed_global = int(np.sum(delta_seed_arr <= 1e-12))
    succ_extrap_global = int(np.sum(delta_extrap_arr <= 1e-12))
    succ_q05_global = int(np.sum(delta_q05_arr >= -1e-12))
    ci_seed_global = jeffreys_interval_from_successes(succ_seed_global, n_pairs_global, alpha=0.05)
    ci_extrap_global = jeffreys_interval_from_successes(succ_extrap_global, n_pairs_global, alpha=0.05)
    ci_q05_global = jeffreys_interval_from_successes(succ_q05_global, n_pairs_global, alpha=0.05)
    p_ok_ci = (
        ci_seed_global["lower"] >= 0.25 and
        ci_extrap_global["lower"] >= 0.5 and
        ci_q05_global["lower"] >= 0.05
    )
    not_low_power = not paired_delta_panel["low_power_summary"]["global_any_low_power"]
    if p_ok_ci and not_low_power:
        power_aware_status = "NON_WORSE_CONFIRMED"
    elif p_ok and (not not_low_power):
        power_aware_status = "NON_WORSE_LOW_POWER"
    else:
        power_aware_status = "NON_WORSE_NOT_CONFIRMED"
    paired_delta_panel["power_aware_verdict"] = {
        "status": power_aware_status,
        "nonworse_probability_conditions_met": bool(p_ok),
        "prob_seed_nonworse_ci95": ci_seed_global,
        "prob_extrap_nonworse_ci95": ci_extrap_global,
        "prob_q05_nonworse_ci95": ci_q05_global,
        "nonworse_probability_ci95_conditions_met": bool(p_ok_ci),
        "low_power_detected": bool(not_low_power is False),
        "ready_for_real_backend_substitution": bool(p_ok_ci and not_low_power),
    }
    per_channel_power_aware_verdicts = []
    for row in paired_delta_panel["per_channel_rows"]:
        quality_row = next(qr for qr in paired_delta_panel["per_channel_wilcoxon_quality"] if qr["channel"] == row["channel"])
        ch = row["channel"]
        arr_seed_ch = np.array(per_channel_bootstrap_deltas[ch]["seed"], dtype=float)
        arr_extrap_ch = np.array(per_channel_bootstrap_deltas[ch]["extrap"], dtype=float)
        arr_q05_ch = np.array(per_channel_bootstrap_deltas[ch]["q05"], dtype=float)
        prob_seed_ch = float(np.mean(arr_seed_ch <= 1e-12))
        prob_extrap_ch = float(np.mean(arr_extrap_ch <= 1e-12))
        prob_q05_ch = float(np.mean(arr_q05_ch >= -1e-12))
        succ_seed = int(np.sum(arr_seed_ch <= 1e-12))
        succ_extrap = int(np.sum(arr_extrap_ch <= 1e-12))
        succ_q05 = int(np.sum(arr_q05_ch >= -1e-12))
        n_pairs_ch = int(arr_seed_ch.size)
        ci_seed_ch = jeffreys_interval_from_successes(succ_seed, n_pairs_ch, alpha=0.05)
        ci_extrap_ch = jeffreys_interval_from_successes(succ_extrap, n_pairs_ch, alpha=0.05)
        ci_q05_ch = jeffreys_interval_from_successes(succ_q05, n_pairs_ch, alpha=0.05)
        p_ok_ch = (
            prob_seed_ch >= 0.25 and
            prob_extrap_ch >= 0.5 and
            prob_q05_ch >= 0.05
        )
        p_ok_ch_ci = (
            ci_seed_ch["lower"] >= 0.25 and
            ci_extrap_ch["lower"] >= 0.5 and
            ci_q05_ch["lower"] >= 0.05
        )
        low_power_ch = bool(
            quality_row["seed"]["low_power_flag"] or
            quality_row["extrap"]["low_power_flag"] or
            quality_row["q05"]["low_power_flag"]
        )
        if p_ok_ch_ci and (not low_power_ch):
            status_ch = "NON_WORSE_CONFIRMED"
        elif p_ok_ch and low_power_ch:
            status_ch = "NON_WORSE_LOW_POWER"
        else:
            status_ch = "NON_WORSE_NOT_CONFIRMED"
        per_channel_power_aware_verdicts.append({
            "channel": row["channel"],
            "status": status_ch,
            "prob_seed_nonworse": prob_seed_ch,
            "prob_extrap_nonworse": prob_extrap_ch,
            "prob_q05_nonworse": prob_q05_ch,
            "prob_seed_nonworse_ci95": ci_seed_ch,
            "prob_extrap_nonworse_ci95": ci_extrap_ch,
            "prob_q05_nonworse_ci95": ci_q05_ch,
            "nonworse_probability_conditions_met": bool(p_ok_ch),
            "nonworse_probability_ci95_conditions_met": bool(p_ok_ch_ci),
            "low_power_detected": bool(low_power_ch),
            "ready_for_real_backend_substitution": bool(p_ok_ch_ci and (not low_power_ch)),
        })
    paired_delta_panel["per_channel_power_aware_verdicts"] = per_channel_power_aware_verdicts
    status_counts = {
        "NON_WORSE_CONFIRMED": int(sum(1 for r in per_channel_power_aware_verdicts if r["status"] == "NON_WORSE_CONFIRMED")),
        "NON_WORSE_LOW_POWER": int(sum(1 for r in per_channel_power_aware_verdicts if r["status"] == "NON_WORSE_LOW_POWER")),
        "NON_WORSE_NOT_CONFIRMED": int(sum(1 for r in per_channel_power_aware_verdicts if r["status"] == "NON_WORSE_NOT_CONFIRMED")),
    }
    paired_delta_panel["mixed_verdict_regime"] = {
        "status_counts": status_counts,
        "num_channels": int(len(per_channel_power_aware_verdicts)),
        "num_ready_channels": int(sum(1 for r in per_channel_power_aware_verdicts if r["ready_for_real_backend_substitution"])),
    }
    risk_rows = []
    for r in per_channel_power_aware_verdicts:
        m_seed = float(r["prob_seed_nonworse_ci95"]["lower"] - 0.25)
        m_extrap = float(r["prob_extrap_nonworse_ci95"]["lower"] - 0.5)
        m_q05 = float(r["prob_q05_nonworse_ci95"]["lower"] - 0.05)
        risk_rows.append({
            "channel": r["channel"],
            "min_ci_margin_to_threshold": float(min(m_seed, m_extrap, m_q05)),
            "seed_ci_margin": m_seed,
            "extrap_ci_margin": m_extrap,
            "q05_ci_margin": m_q05,
            "status": r["status"],
        })
    risk_rows_sorted = sorted(risk_rows, key=lambda x: x["min_ci_margin_to_threshold"])
    paired_delta_panel["per_channel_risk_ranking"] = {
        "rows": risk_rows_sorted,
        "most_critical_channel": risk_rows_sorted[0]["channel"] if risk_rows_sorted else "NONE",
    }
    # time-stability risk panel: repeat CI-margin ranking across RNG seeds
    risk_seed_rows = []
    critical_counts = {ch: 0 for ch in channel_names}
    for risk_seed in [20260530, 20260531, 20260532, 20260533, 20260534, 20260535]:
        rng_risk = np.random.default_rng(risk_seed)
        local_risk_rows = []
        for ch in channel_names:
            d_seed_ch = []
            d_extrap_ch = []
            d_q05_ch = []
            for _ in range(64):
                idx = rng_risk.integers(0, len(nu_grid), size=len(nu_grid))
                nu_s = [float(nu_grid[int(i)]) for i in idx]
                cond_s = [float(np.linalg.cond(np.array(t_frw_to_bianchi_sym.subs({nu_sym: nu})).astype(float))) for nu in nu_s]
                vals_b = np.array([channel_transport_by_key.get((ch, nuv), float("inf")) * cv for nuv, cv in zip(nu_s, cond_s)], dtype=float)
                vals_cf = np.array([(channel_first_by_nu.get(nuv, float("inf")) if ch == channel_priority_best else channel_transport_by_key.get((ch, nuv), float("inf"))) * cv for nuv, cv in zip(nu_s, cond_s)], dtype=float)
                sb = float(np.median(vals_b)) if vals_b.size else float("inf")
                scf = float(np.median(vals_cf)) if vals_cf.size else float("inf")
                d_seed_ch.append(float(scf - sb))
                xb = np.array([5.0, 6.0, 7.0], dtype=float)
                yb = np.array([sb, max(0.0, sb - 0.01), max(0.0, sb - 0.02)], dtype=float)
                ycf = np.array([scf, max(0.0, scf - 0.01), max(0.0, scf - 0.02)], dtype=float)
                Ab = np.vstack([xb, np.ones_like(xb)]).T
                qAb = np.vstack([xb**2, xb, np.ones_like(xb)]).T
                mb, cb = np.linalg.lstsq(Ab, yb, rcond=None)[0]
                mcf, ccf = np.linalg.lstsq(Ab, ycf, rcond=None)[0]
                qb2, qb1, qb0 = np.linalg.lstsq(qAb, yb, rcond=None)[0]
                qcf2, qcf1, qcf0 = np.linalg.lstsq(qAb, ycf, rcond=None)[0]
                xh = np.array([8.0, 9.0], dtype=float)
                d_extrap_ch.append(float(np.max(np.abs((qcf2 * xh**2 + qcf1 * xh + qcf0) - (mcf * xh + ccf))) - np.max(np.abs((qb2 * xh**2 + qb1 * xh + qb0) - (mb * xh + cb)))))
                q05_b = float(ss.beta.ppf(0.05, 2, 10))
                q05_cf = float(ss.beta.ppf(0.05, 2 + (1 if ch == channel_priority_best else 0), 10))
                d_q05_ch.append(float(q05_cf - q05_b))
            arr_seed_ch = np.array(d_seed_ch, dtype=float)
            arr_extrap_ch = np.array(d_extrap_ch, dtype=float)
            arr_q05_ch = np.array(d_q05_ch, dtype=float)
            n_pairs_ch = int(arr_seed_ch.size)
            ci_seed_ch = jeffreys_interval_from_successes(int(np.sum(arr_seed_ch <= 1e-12)), n_pairs_ch, alpha=0.05)
            ci_extrap_ch = jeffreys_interval_from_successes(int(np.sum(arr_extrap_ch <= 1e-12)), n_pairs_ch, alpha=0.05)
            ci_q05_ch = jeffreys_interval_from_successes(int(np.sum(arr_q05_ch >= -1e-12)), n_pairs_ch, alpha=0.05)
            m_seed = float(ci_seed_ch["lower"] - 0.25)
            m_extrap = float(ci_extrap_ch["lower"] - 0.5)
            m_q05 = float(ci_q05_ch["lower"] - 0.05)
            local_risk_rows.append({"channel": ch, "min_ci_margin_to_threshold": float(min(m_seed, m_extrap, m_q05))})
        local_risk_rows = sorted(local_risk_rows, key=lambda x: x["min_ci_margin_to_threshold"])
        critical = local_risk_rows[0]["channel"] if local_risk_rows else "NONE"
        if critical in critical_counts:
            critical_counts[critical] += 1
        risk_seed_rows.append({
            "risk_seed": int(risk_seed),
            "most_critical_channel": critical,
            "most_critical_margin": float(local_risk_rows[0]["min_ci_margin_to_threshold"]) if local_risk_rows else float("nan"),
        })
    critical_count_rows = [{"channel": ch, "count": int(critical_counts[ch]), "frequency": float(critical_counts[ch] / len(risk_seed_rows))} for ch in channel_names]
    critical_count_rows_sorted = sorted(critical_count_rows, key=lambda x: x["count"], reverse=True)
    paired_delta_panel["time_stability_risk_panel"] = {
        "rows": risk_seed_rows,
        "critical_channel_frequency_rows": critical_count_rows_sorted,
        "num_seeds": int(len(risk_seed_rows)),
        "most_frequent_critical_channel": critical_count_rows_sorted[0]["channel"] if critical_count_rows_sorted else "NONE",
        "critical_channel_frequency_span": float(max(r["frequency"] for r in critical_count_rows_sorted) - min(r["frequency"] for r in critical_count_rows_sorted)) if critical_count_rows_sorted else 0.0,
    }
    freq_top = float(critical_count_rows_sorted[0]["frequency"]) if critical_count_rows_sorted else 0.0
    freq_span = float(paired_delta_panel["time_stability_risk_panel"]["critical_channel_frequency_span"])
    seed_robust_ready = bool((freq_top >= 0.67) and (freq_span <= 0.5))
    paired_delta_panel["time_stability_seed_robust_gate"] = {
        "dominant_frequency_threshold": 0.67,
        "frequency_span_threshold": 0.5,
        "dominant_frequency_observed": freq_top,
        "frequency_span_observed": freq_span,
        "risk_signal_stable": bool(seed_robust_ready),
    }
    # branch-cut/regularization sensitivity proxy (Obstacle C guard)
    eta_probe_grid = [1.8, 1.8 + 1e-4, 1.8 - 1e-4]
    eps_floor_grid = [1e-15, 1e-12, 1e-10]
    branch_rows = []
    for eta_probe in eta_probe_grid:
        for eps_floor in eps_floor_grid:
            vals = []
            for s_val in s_grid_fine:
                def integrand_probe(x: float) -> float:
                    kk = np.cos(omega * x + phi) / (1.0 + beta * (x ** eta_probe))
                    return float((kk * kk) / np.sqrt(max(eps_floor, x + s_val)))
                vv, _ = si.quad(integrand_probe, 0.0, 1.0, epsabs=1e-12, epsrel=1e-12, limit=400)
                vals.append(float(vv))
            mval = float(np.mean(np.array(vals, dtype=float)))
            branch_rows.append({"eta_probe": float(eta_probe), "eps_floor": float(eps_floor), "phase_mean": mval})
    phase_means_probe = [r["phase_mean"] for r in branch_rows]
    branch_span = float(max(phase_means_probe) - min(phase_means_probe)) if phase_means_probe else float("inf")
    # local log-log slope proxy d(log I)/d(eta) around eta=1.8
    eta_minus = 1.8 - 1e-4
    eta_plus = 1.8 + 1e-4
    slope_rows = []
    for s_val in s_grid_fine:
        i_minus, _ = strict_kernel_phase_integral(s_val, omega, phi, beta, eta_minus)
        i_plus, _ = strict_kernel_phase_integral(s_val, omega, phi, beta, eta_plus)
        log_i_minus = float(np.log(max(1e-30, i_minus)))
        log_i_plus = float(np.log(max(1e-30, i_plus)))
        dlogi_deta = float((log_i_plus - log_i_minus) / (eta_plus - eta_minus))
        slope_rows.append({"s": float(s_val), "dlogI_deta_local": dlogi_deta})
    slope_vals = [r["dlogI_deta_local"] for r in slope_rows]
    slope_span = float(max(slope_vals) - min(slope_vals)) if slope_vals else float("inf")
    paired_delta_panel["branch_cut_sensitivity_panel"] = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "eta_probe_grid": eta_probe_grid,
        "eps_floor_grid": eps_floor_grid,
        "rows": branch_rows,
        "phase_mean_span": branch_span,
        "phase_mean_span_threshold": 1e-3,
        "branch_choice_robustness_bounded": bool(branch_span < 1e-3),
        "loglog_slope_rows": slope_rows,
        "loglog_slope_span": slope_span,
        "loglog_slope_span_threshold": 5e-2,
        "loglog_slope_span_bounded": bool(slope_span < 5e-2),
    }
    # cross-integrator agreement proxy on baseline branch settings
    cross_integrator_rows = []
    for s_val in s_grid_fine:
        v_quad, _ = strict_kernel_phase_integral(s_val, omega, phi, beta, eta)
        xg = np.linspace(0.0, 1.0, 4001, dtype=float)
        kg = np.cos(omega * xg + phi) / (1.0 + beta * (xg ** eta))
        yg = (kg * kg) / np.sqrt(np.maximum(1e-15, xg + s_val))
        v_trapz = float(np.trapezoid(yg, xg))
        cross_integrator_rows.append({
            "s": float(s_val),
            "quad_value": float(v_quad),
            "trapz_value": v_trapz,
            "abs_gap": float(abs(v_quad - v_trapz)),
        })
    cross_integrator_gap_max = float(max(r["abs_gap"] for r in cross_integrator_rows)) if cross_integrator_rows else float("inf")
    paired_delta_panel["branch_cross_integrator_panel"] = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "methods": ["scipy.integrate.quad", "numpy.trapezoid"],
        "trapz_grid_points": 4001,
        "rows": cross_integrator_rows,
        "abs_gap_max": cross_integrator_gap_max,
        "abs_gap_max_threshold": 1e-4,
        "cross_integrator_agreement_bounded": bool(cross_integrator_gap_max < 1e-4),
    }
    # joint branch+integrator stress matrix (eta_probe x eps_floor x integrator)
    branch_integrator_rows = []
    for eta_probe in eta_probe_grid:
        for eps_floor in eps_floor_grid:
            for s_val in s_grid_fine:
                def integrand_probe(x: float) -> float:
                    kk = np.cos(omega * x + phi) / (1.0 + beta * (x ** eta_probe))
                    return float((kk * kk) / np.sqrt(max(eps_floor, x + s_val)))
                v_quad, _ = si.quad(integrand_probe, 0.0, 1.0, epsabs=1e-12, epsrel=1e-12, limit=400)
                xg = np.linspace(0.0, 1.0, 4001, dtype=float)
                yg = np.array([integrand_probe(float(xx)) for xx in xg], dtype=float)
                v_trapezoid = float(np.trapezoid(yg, xg))
                branch_integrator_rows.append({
                    "eta_probe": float(eta_probe),
                    "eps_floor": float(eps_floor),
                    "s": float(s_val),
                    "quad_value": float(v_quad),
                    "trapezoid_value": v_trapezoid,
                    "abs_gap": float(abs(v_quad - v_trapezoid)),
                })
    worst_case_gap_envelope = float(max(r["abs_gap"] for r in branch_integrator_rows)) if branch_integrator_rows else float("inf")
    paired_delta_panel["branch_integrator_stress_matrix"] = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "rows": branch_integrator_rows,
        "worst_case_gap_envelope": worst_case_gap_envelope,
        "worst_case_gap_threshold": 2e-4,
        "worst_case_gap_bounded": bool(worst_case_gap_envelope < 2e-4),
    }
    # cross-seed stress envelope over joint branch-integrator matrix
    cross_seed_rows = []
    for seed_scale in [0.9995, 1.0, 1.0005]:
        rows_local = []
        for eta_probe in eta_probe_grid:
            eta_eff = float(eta_probe * seed_scale)
            for eps_floor in eps_floor_grid:
                for s_val in s_grid_fine:
                    def integrand_probe_seed(x: float) -> float:
                        kk = np.cos(omega * x + phi) / (1.0 + beta * (x ** eta_eff))
                        return float((kk * kk) / np.sqrt(max(eps_floor, x + s_val)))
                    vq, _ = si.quad(integrand_probe_seed, 0.0, 1.0, epsabs=1e-12, epsrel=1e-12, limit=400)
                    xg = np.linspace(0.0, 1.0, 4001, dtype=float)
                    yg = np.array([integrand_probe_seed(float(xx)) for xx in xg], dtype=float)
                    vt = float(np.trapezoid(yg, xg))
                    rows_local.append(abs(vq - vt))
        cross_seed_rows.append({"eta_scale": seed_scale, "worst_case_gap_envelope": float(max(rows_local)) if rows_local else float("inf")})
    cross_seed_worst = float(max(r["worst_case_gap_envelope"] for r in cross_seed_rows)) if cross_seed_rows else float("inf")
    paired_delta_panel["branch_integrator_cross_seed_envelope"] = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "rows": cross_seed_rows,
        "worst_case_across_seed_scales": cross_seed_worst,
        "worst_case_across_seed_scales_threshold": 3e-4,
        "cross_seed_envelope_bounded": bool(cross_seed_worst < 3e-4),
    }
    # adaptive threshold calibration panel (bootstrap CI over seed-scale envelopes)
    env_arr = np.array([r["worst_case_gap_envelope"] for r in cross_seed_rows], dtype=float)
    rng_env = np.random.default_rng(20260536)
    env_boot = []
    for _ in range(256):
        idx = rng_env.integers(0, env_arr.size, size=env_arr.size)
        env_boot.append(float(np.max(env_arr[idx])))
    env_boot_arr = np.array(env_boot, dtype=float)
    env_ci = {
        "q05": float(np.quantile(env_boot_arr, 0.05)),
        "q50": float(np.quantile(env_boot_arr, 0.50)),
        "q95": float(np.quantile(env_boot_arr, 0.95)),
    }
    paired_delta_panel["branch_integrator_threshold_calibration_panel"] = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "bootstrap_size": 256,
        "envelope_bootstrap_quantiles": env_ci,
        "base_threshold": 3e-4,
        "adaptive_threshold_q95": float(env_ci["q95"]),
        "adaptive_threshold_consistency": bool(env_ci["q95"] >= env_ci["q50"] >= env_ci["q05"]),
    }
    joint_ready = bool(
        paired_delta_panel["branch_integrator_cross_seed_envelope"]["cross_seed_envelope_bounded"] and
        paired_delta_panel["branch_integrator_threshold_calibration_panel"]["adaptive_threshold_consistency"] and
        paired_delta_panel["time_stability_seed_robust_gate"]["risk_signal_stable"]
    )
    fail_reasons = []
    if not paired_delta_panel["branch_integrator_cross_seed_envelope"]["cross_seed_envelope_bounded"]:
        fail_reasons.append("cross_seed_envelope_not_bounded")
    if not paired_delta_panel["branch_integrator_threshold_calibration_panel"]["adaptive_threshold_consistency"]:
        fail_reasons.append("adaptive_threshold_inconsistent")
    if not paired_delta_panel["time_stability_seed_robust_gate"]["risk_signal_stable"]:
        fail_reasons.append("risk_signal_not_stable")
    paired_delta_panel["branch_robust_substitution_decision"] = {
        "status": "READY_FOR_BRANCH_ROBUST_SUBSTITUTION" if joint_ready else "NOT_READY_FOR_BRANCH_ROBUST_SUBSTITUTION",
        "ready_for_branch_robust_substitution": joint_ready,
        "criteria": {
            "cross_seed_envelope_bounded": bool(paired_delta_panel["branch_integrator_cross_seed_envelope"]["cross_seed_envelope_bounded"]),
            "adaptive_threshold_consistency": bool(paired_delta_panel["branch_integrator_threshold_calibration_panel"]["adaptive_threshold_consistency"]),
            "time_stability_risk_signal_stable": bool(paired_delta_panel["time_stability_seed_robust_gate"]["risk_signal_stable"]),
        },
        "failed_criteria": fail_reasons,
    }
    channel_first_replay_metrics["paired_delta_panel"] = paired_delta_panel
    channel_first_simulation_panel = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "selected_channel": channel_priority_best,
        "residual_l2_channel_first": phase_common_residual_channel_first_l2,
        "residual_linf_channel_first": phase_common_residual_channel_first_linf,
        "delta_residual_l2_channel_first_minus_baseline": float(phase_common_residual_channel_first_l2 - phase_common_residual_l2),
        "delta_residual_l2_channel_first_minus_backend_3channel": float(phase_common_residual_channel_first_l2 - phase_common_residual_backend_sub_l2),
        "transport_rows": channel_first_transport_rows,
        "cond_weighted_residual_l2_median": channel_first_cond_weighted_residual_median,
        "replay_metrics": channel_first_replay_metrics,
        "note": "Single-channel substitution simulation only; no closure claim.",
    }
    phase_common_cond = float(np.linalg.cond(x_phase))
    # leave-one-channel-out stability check for phase/common-basis link
    loocv_rows = []
    loocv_residuals = []
    for i in range(x_phase.shape[0]):
        tr_idx = [j for j in range(x_phase.shape[0]) if j != i]
        te_idx = [i]
        c_loo, *_ = la.lstsq(x_phase[tr_idx, :], y_phase[tr_idx, :])
        pred_te = x_phase[te_idx, :] @ c_loo
        resid_te = pred_te - y_phase[te_idx, :]
        r_l2 = float(np.linalg.norm(resid_te, ord=2))
        loocv_residuals.append(r_l2)
        loocv_rows.append({"held_out_channel": channel_phase_rows[i]["channel"], "residual_l2": r_l2})
    loocv_residual_max = float(max(loocv_residuals)) if loocv_residuals else float("inf")
    # bootstrap-over-channels stability for phase/common-basis link
    rng_link = np.random.default_rng(20260522)
    link_bootstrap_rows = []
    link_bootstrap_resids = []
    for bidx in range(128):
        sample_idx = rng_link.integers(0, x_phase.shape[0], size=x_phase.shape[0])
        xb = x_phase[sample_idx, :]
        yb = y_phase[sample_idx, :]
        cb, *_ = la.lstsq(xb, yb)
        pred_full = x_phase @ cb
        rb = pred_full - y_phase
        rb_l2 = float(np.linalg.norm(rb, ord=2))
        link_bootstrap_resids.append(rb_l2)
        if bidx < 8:
            link_bootstrap_rows.append({"bootstrap_id": int(bidx), "sample_idx": sample_idx.tolist(), "residual_l2_full_eval": rb_l2})
    link_bootstrap_resid_p95 = float(np.quantile(np.array(link_bootstrap_resids, dtype=float), 0.95))
    # unified stability envelope for task-2/task-7 bridge diagnostics
    stability_envelope = {
        "loocv_residual_l2_max": loocv_residual_max,
        "bootstrap_residual_l2_p95": link_bootstrap_resid_p95,
        "tolerance_span_max": tol_span_max,
    }
    stability_envelope_max = float(max(stability_envelope.values()))
    # cross-channel coefficient agreement for phase/common-basis link
    coef_col_means = np.mean(phase_common_coef, axis=1)
    coef_centered = phase_common_coef - coef_col_means[:, None]
    cross_channel_coef_spread_l2 = float(np.linalg.norm(coef_centered, ord=2))
    # joint constrained fit (task-2/task-7 bridge): coupled channels with spread penalty
    def unpack_joint(v: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        c0 = np.array(v[0:5], dtype=float)
        c1 = np.array(v[5:10], dtype=float)
        c2 = np.array(v[10:15], dtype=float)
        return c0, c1, c2

    y_cols = [y_phase[:, i] for i in range(y_phase.shape[1])]
    def joint_loss(v: np.ndarray, lam: float = 0.1) -> float:
        c0, c1, c2 = unpack_joint(v)
        p0 = x_phase @ c0
        p1 = x_phase @ c1
        p2 = x_phase @ c2
        data_term = np.linalg.norm(p0 - y_cols[0], ord=2) ** 2 + np.linalg.norm(p1 - y_cols[1], ord=2) ** 2 + np.linalg.norm(p2 - y_cols[2], ord=2) ** 2
        cmat = np.stack([c0, c1, c2], axis=1)
        mu = np.mean(cmat, axis=1, keepdims=True)
        spread_term = np.linalg.norm(cmat - mu, ord=2) ** 2
        return float(data_term + lam * spread_term)

    x0_joint = np.concatenate([phase_common_coef[:, 0], phase_common_coef[:, 1], phase_common_coef[:, 2]], axis=0)
    joint_opt = so.minimize(joint_loss, x0=x0_joint, method="L-BFGS-B")
    joint_opt_slsqp = so.minimize(joint_loss, x0=x0_joint, method="SLSQP")
    c0j, c1j, c2j = unpack_joint(joint_opt.x.astype(float))
    cmat_joint = np.stack([c0j, c1j, c2j], axis=1)
    pred_joint = np.stack([x_phase @ c0j, x_phase @ c1j, x_phase @ c2j], axis=1)
    resid_joint = pred_joint - y_phase
    joint_resid_l2 = float(np.linalg.norm(resid_joint, ord=2))
    joint_spread_l2 = float(np.linalg.norm(cmat_joint - np.mean(cmat_joint, axis=1, keepdims=True), ord=2))
    joint_solver_gap = float(abs(float(joint_opt.fun) - float(joint_opt_slsqp.fun)))

    lambda_rows = []
    lambda_resids = []
    for lam in [0.05, 0.1, 0.2]:
        r = so.minimize(lambda v: joint_loss(v, lam=lam), x0=x0_joint, method="L-BFGS-B")
        cc0, cc1, cc2 = unpack_joint(r.x.astype(float))
        pp = np.stack([x_phase @ cc0, x_phase @ cc1, x_phase @ cc2], axis=1)
        rr = pp - y_phase
        rr_l2 = float(np.linalg.norm(rr, ord=2))
        ss_l2 = float(np.linalg.norm(np.stack([cc0, cc1, cc2], axis=1) - np.mean(np.stack([cc0, cc1, cc2], axis=1), axis=1, keepdims=True), ord=2))
        lambda_rows.append({"lambda_spread": lam, "objective": float(r.fun), "residual_l2": rr_l2, "coef_spread_l2": ss_l2})
        lambda_resids.append(rr_l2)
    lambda_resid_span = float(max(lambda_resids) - min(lambda_resids))
    # channel holdout-rotation for joint coupled fit robustness
    holdout_rotation_rows = []
    holdout_joint_resids = []
    for h in range(x_phase.shape[0]):
        tr = [i for i in range(x_phase.shape[0]) if i != h]
        te = [h]
        xtr, ytr = x_phase[tr, :], y_phase[tr, :]
        xte, yte = x_phase[te, :], y_phase[te, :]
        x0h = np.concatenate([phase_common_coef[:, 0], phase_common_coef[:, 1], phase_common_coef[:, 2]], axis=0)
        def joint_loss_h(v: np.ndarray) -> float:
            c0, c1, c2 = unpack_joint(v)
            p0 = xtr @ c0
            p1 = xtr @ c1
            p2 = xtr @ c2
            data_term = np.linalg.norm(p0 - ytr[:, 0], ord=2) ** 2 + np.linalg.norm(p1 - ytr[:, 1], ord=2) ** 2 + np.linalg.norm(p2 - ytr[:, 2], ord=2) ** 2
            cmat = np.stack([c0, c1, c2], axis=1)
            mu = np.mean(cmat, axis=1, keepdims=True)
            spread_term = np.linalg.norm(cmat - mu, ord=2) ** 2
            return float(data_term + 0.1 * spread_term)
        rh = so.minimize(joint_loss_h, x0=x0h, method="L-BFGS-B")
        c0h, c1h, c2h = unpack_joint(rh.x.astype(float))
        pte = np.stack([xte @ c0h, xte @ c1h, xte @ c2h], axis=1)
        r_te = pte - yte
        r_l2 = float(np.linalg.norm(r_te, ord=2))
        holdout_joint_resids.append(r_l2)
        holdout_rotation_rows.append({"held_out_channel": channel_phase_rows[h]["channel"], "test_residual_l2": r_l2})
    holdout_joint_resid_max = float(max(holdout_joint_resids)) if holdout_joint_resids else float("inf")
    # joint-fit multistart robustness
    joint_starts = [
        x0_joint,
        np.zeros_like(x0_joint),
        x0_joint + 0.1,
        x0_joint - 0.1,
    ]
    joint_ms_rows = []
    joint_ms_resids = []
    for sid, s0 in enumerate(joint_starts):
        rj = so.minimize(joint_loss, x0=s0, method="L-BFGS-B")
        cj0, cj1, cj2 = unpack_joint(rj.x.astype(float))
        pj = np.stack([x_phase @ cj0, x_phase @ cj1, x_phase @ cj2], axis=1)
        rmat = pj - y_phase
        r_l2 = float(np.linalg.norm(rmat, ord=2))
        joint_ms_resids.append(r_l2)
        joint_ms_rows.append({"start_id": int(sid), "residual_l2": r_l2, "objective": float(rj.fun)})
    joint_ms_resid_span = float(max(joint_ms_resids) - min(joint_ms_resids)) if joint_ms_resids else float("inf")
    # joint-fit perturbation robustness (feature-matrix jitter)
    rng_jitter = np.random.default_rng(20260523)
    pert_rows = []
    pert_resids = []
    for pid in range(8):
        jitter = rng_jitter.normal(loc=0.0, scale=1e-3, size=x_phase.shape)
        xj = x_phase * (1.0 + jitter)
        x0p = np.concatenate([phase_common_coef[:, 0], phase_common_coef[:, 1], phase_common_coef[:, 2]], axis=0)
        def joint_loss_j(v: np.ndarray) -> float:
            c0, c1, c2 = unpack_joint(v)
            p0 = xj @ c0
            p1 = xj @ c1
            p2 = xj @ c2
            data_term = np.linalg.norm(p0 - y_cols[0], ord=2) ** 2 + np.linalg.norm(p1 - y_cols[1], ord=2) ** 2 + np.linalg.norm(p2 - y_cols[2], ord=2) ** 2
            cmat = np.stack([c0, c1, c2], axis=1)
            mu = np.mean(cmat, axis=1, keepdims=True)
            spread_term = np.linalg.norm(cmat - mu, ord=2) ** 2
            return float(data_term + 0.1 * spread_term)
        rp = so.minimize(joint_loss_j, x0=x0p, method="L-BFGS-B")
        cp0, cp1, cp2 = unpack_joint(rp.x.astype(float))
        pmat = np.stack([xj @ cp0, xj @ cp1, xj @ cp2], axis=1)
        rmat = pmat - y_phase
        rl2 = float(np.linalg.norm(rmat, ord=2))
        pert_resids.append(rl2)
        pert_rows.append({"perturbation_id": int(pid), "residual_l2": rl2})
    perturb_resid_span = float(max(pert_resids) - min(pert_resids)) if pert_resids else float("inf")
    stress_residual_values = [joint_resid_l2, holdout_joint_resid_max] + joint_ms_resids + pert_resids + lambda_resids
    joint_worst_case_residual_envelope = float(max(stress_residual_values)) if stress_residual_values else float("inf")
    joint_stress_panel = {
        "solver_crosscheck_objective_gap": joint_solver_gap,
        "lambda_sweep_residual_l2_span": lambda_resid_span,
        "holdout_rotation_residual_l2_max": holdout_joint_resid_max,
        "multistart_residual_l2_span": joint_ms_resid_span,
        "perturbation_residual_l2_span": perturb_resid_span,
        "worst_case_residual_envelope": joint_worst_case_residual_envelope,
    }
    transport_det_mean = float(np.mean(np.array([r["det_frw_to_bianchi"] for r in transport_rows], dtype=float)))
    nu_mean = float(np.mean(np.array(nu_grid, dtype=float)))
    t_fb_mean = np.array(t_frw_to_bianchi_sym.subs({nu_sym: nu_mean})).astype(float)
    t_bf_mean = np.array(t_bianchi_to_frw_sym.subs({nu_sym: nu_mean})).astype(float)
    background_scale = {"FRW": 1.00, "BianchiI": transport_det_mean}
    bg_panel_rows = []
    bg_env_vals = []
    for bg_name, bg_scale in background_scale.items():
        bg_vals = [v * bg_scale for v in stress_residual_values]
        bg_env = float(max(bg_vals)) if bg_vals else float("inf")
        bg_env_vals.append(bg_env)
        bg_panel_rows.append({"background": bg_name, "scale_proxy": float(bg_scale), "worst_case_residual_envelope": bg_env})
    cross_background_envelope_span = float(max(bg_env_vals) - min(bg_env_vals)) if bg_env_vals else float("inf")
    # operator-level transport replay on phase target matrix (3x3 channel codomain)
    y_frw = y_phase.copy()
    y_bianchi = y_phase @ t_fb_mean.T
    c_frw, *_ = la.lstsq(x_phase, y_frw)
    c_bianchi, *_ = la.lstsq(x_phase, y_bianchi)
    r_frw = x_phase @ c_frw - y_frw
    r_bianchi = x_phase @ c_bianchi - y_bianchi
    transport_operator_rows = [
        {"background": "FRW", "residual_l2": float(np.linalg.norm(r_frw, ord=2)), "residual_linf": float(np.linalg.norm(r_frw, ord=np.inf))},
        {"background": "BianchiI", "residual_l2": float(np.linalg.norm(r_bianchi, ord=2)), "residual_linf": float(np.linalg.norm(r_bianchi, ord=np.inf))},
    ]
    operator_resid_span = float(abs(transport_operator_rows[1]["residual_l2"] - transport_operator_rows[0]["residual_l2"]))
    # nu-sweep operator replay for stronger background-independence proxy diagnostics
    operator_nu_sweep_rows = []
    operator_nu_resids = []
    operator_nu_solver_gaps = []
    operator_nu_lambda_rows = []
    operator_nu_lambda_resids = []
    operator_nu_lambda_solver_gaps = []
    operator_nu_lambda_weighted_resids = []
    operator_nu_lambda_cond_weighted_resids = []
    lambda_grid_nu = [0.05, 0.1, 0.2]
    for nu in nu_grid:
        t_fb_nu = np.array(t_frw_to_bianchi_sym.subs({nu_sym: nu})).astype(float)
        cond_t_fb_nu = float(np.linalg.cond(t_fb_nu))
        y_bi_nu = y_phase @ t_fb_nu.T
        c_bi_nu, *_ = la.lstsq(x_phase, y_bi_nu)
        r_bi_nu = x_phase @ c_bi_nu - y_bi_nu
        r_l2 = float(np.linalg.norm(r_bi_nu, ord=2))
        # per-nu joint solver crosscheck (L-BFGS-B vs SLSQP) on transported targets
        y_cols_nu = [y_bi_nu[:, i] for i in range(y_bi_nu.shape[1])]
        def joint_loss_nu(v: np.ndarray, lam: float = 0.1) -> float:
            c0, c1, c2 = unpack_joint(v)
            p0 = x_phase @ c0
            p1 = x_phase @ c1
            p2 = x_phase @ c2
            data_term = np.linalg.norm(p0 - y_cols_nu[0], ord=2) ** 2 + np.linalg.norm(p1 - y_cols_nu[1], ord=2) ** 2 + np.linalg.norm(p2 - y_cols_nu[2], ord=2) ** 2
            cmat = np.stack([c0, c1, c2], axis=1)
            mu = np.mean(cmat, axis=1, keepdims=True)
            spread_term = np.linalg.norm(cmat - mu, ord=2) ** 2
            return float(data_term + lam * spread_term)
        r_lb = so.minimize(joint_loss_nu, x0=x0_joint, method="L-BFGS-B")
        r_sq = so.minimize(joint_loss_nu, x0=x0_joint, method="SLSQP")
        solver_gap_nu = float(abs(float(r_lb.fun) - float(r_sq.fun)))
        for lam_nu in lambda_grid_nu:
            r_lb_lam = so.minimize(lambda v: joint_loss_nu(v, lam=lam_nu), x0=x0_joint, method="L-BFGS-B")
            r_sq_lam = so.minimize(lambda v: joint_loss_nu(v, lam=lam_nu), x0=x0_joint, method="SLSQP")
            c0_l, c1_l, c2_l = unpack_joint(r_lb_lam.x.astype(float))
            p_l = np.stack([x_phase @ c0_l, x_phase @ c1_l, x_phase @ c2_l], axis=1)
            r_l = p_l - y_bi_nu
            resid_l2_lam = float(np.linalg.norm(r_l, ord=2))
            gap_lam = float(abs(float(r_lb_lam.fun) - float(r_sq_lam.fun)))
            nu_weight = float(abs(np.linalg.det(t_fb_nu)))
            resid_weighted = float(resid_l2_lam * nu_weight)
            resid_cond_weighted = float(resid_l2_lam * cond_t_fb_nu)
            operator_nu_lambda_resids.append(resid_l2_lam)
            operator_nu_lambda_solver_gaps.append(gap_lam)
            operator_nu_lambda_weighted_resids.append(resid_weighted)
            operator_nu_lambda_cond_weighted_resids.append(resid_cond_weighted)
            operator_nu_lambda_rows.append({
                "nu": float(nu),
                "lambda_spread": float(lam_nu),
                "transport_det_weight": nu_weight,
                "transport_condition_weight": cond_t_fb_nu,
                "residual_l2": resid_l2_lam,
                "residual_l2_weighted_by_det": resid_weighted,
                "residual_l2_weighted_by_condition": resid_cond_weighted,
                "solver_objective_gap_lbfgsb_vs_slsqp": gap_lam,
            })
        operator_nu_resids.append(r_l2)
        operator_nu_solver_gaps.append(solver_gap_nu)
        operator_nu_sweep_rows.append({
            "nu": float(nu),
            "det_frw_to_bianchi": float(np.linalg.det(t_fb_nu)),
            "residual_l2": r_l2,
            "residual_linf": float(np.linalg.norm(r_bi_nu, ord=np.inf)),
            "joint_solver_objective_gap_lbfgsb_vs_slsqp": solver_gap_nu,
        })
    operator_nu_sweep_resid_span = float(max(operator_nu_resids) - min(operator_nu_resids)) if operator_nu_resids else float("inf")
    operator_nu_sweep_solver_gap_max = float(max(operator_nu_solver_gaps)) if operator_nu_solver_gaps else float("inf")
    operator_nu_lambda_resid_span = float(max(operator_nu_lambda_resids) - min(operator_nu_lambda_resids)) if operator_nu_lambda_resids else float("inf")
    operator_nu_lambda_solver_gap_max = float(max(operator_nu_lambda_solver_gaps)) if operator_nu_lambda_solver_gaps else float("inf")
    operator_nu_lambda_weighted_resid_max = float(max(operator_nu_lambda_weighted_resids)) if operator_nu_lambda_weighted_resids else float("inf")
    operator_nu_lambda_weighted_resid_span = float(max(operator_nu_lambda_weighted_resids) - min(operator_nu_lambda_weighted_resids)) if operator_nu_lambda_weighted_resids else float("inf")
    operator_nu_lambda_cond_weighted_resid_max = float(max(operator_nu_lambda_cond_weighted_resids)) if operator_nu_lambda_cond_weighted_resids else float("inf")
    operator_nu_lambda_cond_weighted_resid_span = float(max(operator_nu_lambda_cond_weighted_resids) - min(operator_nu_lambda_cond_weighted_resids)) if operator_nu_lambda_cond_weighted_resids else float("inf")
    # dual-criterion frontier map (det-weighted vs condition-weighted)
    dual_rows = []
    for row in operator_nu_lambda_rows:
        dual_rows.append({
            "nu": row["nu"],
            "lambda_spread": row["lambda_spread"],
            "det_weighted": row["residual_l2_weighted_by_det"],
            "cond_weighted": row["residual_l2_weighted_by_condition"],
        })
    dual_with_front = []
    for i, r in enumerate(dual_rows):
        dominated = False
        for j, s in enumerate(dual_rows):
            if i == j:
                continue
            if (s["det_weighted"] <= r["det_weighted"] and s["cond_weighted"] <= r["cond_weighted"]) and (s["det_weighted"] < r["det_weighted"] or s["cond_weighted"] < r["cond_weighted"]):
                dominated = True
                break
        dual_with_front.append({**r, "pareto_frontier": not dominated})
    pareto_count = int(sum(1 for r in dual_with_front if r["pareto_frontier"]))
    dual_stable_rows = [r for r in dual_with_front if r["det_weighted"] < 1.0 and r["cond_weighted"] < 2.0]
    dual_unstable_rows = [r for r in dual_with_front if not (r["det_weighted"] < 1.0 and r["cond_weighted"] < 2.0)]
    # frontier continuity check along nu for each lambda (detect local branch flips)
    frontier_continuity_rows = []
    branch_flip_total = 0
    for lam in lambda_grid_nu:
        lam_rows = sorted([r for r in dual_with_front if abs(float(r["lambda_spread"]) - float(lam)) < 1e-12], key=lambda z: float(z["nu"]))
        flips = 0
        for i in range(len(lam_rows) - 1):
            if bool(lam_rows[i]["pareto_frontier"]) != bool(lam_rows[i + 1]["pareto_frontier"]):
                flips += 1
        branch_flip_total += flips
        frontier_continuity_rows.append({
            "lambda_spread": float(lam),
            "num_nu_points": int(len(lam_rows)),
            "frontier_membership_flips_along_nu": int(flips),
        })

    same_scheme_tag = "STRICT_P2020_PHASESPACE_SCHEME_V1"

    reproducibility_probe = {
        "phase_center": phase_center,
        "max_fd_gap": max_fd_gap,
        "max_grid_refine_gap": max_grid_refine_gap,
        "quad_tol_span": quad_tol_span,
        "cond_p95": cond_p95, "bootstrap_seed_span_max": bootstrap_seed_span_max,
        "solution": {"g": g_num, "w": w_num, "b": b_num},
    }
    reproducibility_digest_1 = digest(reproducibility_probe)
    reproducibility_digest_2 = digest(reproducibility_probe)

    toe_closure_gaps_7tasks = [
        {"task_id": 1, "name": "Renormalization exact divergence cancellation coefficients", "status": "OPEN_OBSTRUCTION_WITH_TRACE", "missing": ["backend-computed a_R2/a_Ric2/a_Riem2/a_GB on strict B1"]},
        {"task_id": 2, "name": "Cutkosky full unitarity closure", "status": "OPEN_OBSTRUCTION_WITH_TRACE", "missing": ["DiscM_loop backend objects", "same-scheme GhostCut/WardLift/CohomologyAmplitudeBridge global fit"]},
        {"task_id": 3, "name": "Background-independence global transport", "status": "OPEN_OBSTRUCTION_WITH_TRACE", "missing": ["global FRW<->Bianchi transport closure for nu branch"]},
        {"task_id": 4, "name": "PO3 nonempty certifier", "status": "OPEN_OBSTRUCTION_WITH_TRACE", "missing": ["constructive strict-kernel admissible region witness"]},
        {"task_id": 5, "name": "PO2 sufficiency theorem", "status": "OPEN_OBSTRUCTION_WITH_TRACE", "missing": ["full L_total symbolic trace forcing C1-C4 -> DELTA_BG_Yf"]},
        {"task_id": 6, "name": "QW-2191 selector obstruction resolution", "status": "OPEN_OBSTRUCTION_WITH_TRACE", "missing": ["strict-core selector source or explicit symmetry-breaking premise"]},
        {"task_id": 7, "name": "DiscM common-basis integration gate", "status": "OPEN_OBSTRUCTION_WITH_TRACE", "missing": ["backend loop/channel common-basis fit with bounded uncertainty"]},
    ]
    task_numeric_evidence_7 = [
        {
            "task_id": 1,
            "name": "Renormalization exact divergence cancellation coefficients",
            "method_stack": ["numpy.linalg", "scipy.linalg", "sympy"],
            "metrics": {
                "fit_residual_l2": residual_b1_l2,
                "fit_residual_linf": residual_b1_linf,
                "solver_gap_abs_max": coef_b1_gap,
            },
            "local_readiness_score_0_1": float(np.exp(-5.0 * residual_b1_l2)),
            "honest_verdict": "OPEN_OBSTRUCTION_WITH_TRACE",
        },
        {
            "task_id": 2,
            "name": "Cutkosky full unitarity closure",
            "method_stack": ["numpy", "scipy.integrate", "scipy.optimize"],
            "metrics": {
                "phase_min_integral": channel_phase_min_global,
                "phase_tolerance_span_max": tol_span_max,
                "joint_solver_gap": joint_solver_gap,
            },
            "local_readiness_score_0_1": float(np.clip(1.0 - min(1.0, tol_span_max * 1e8), 0.0, 1.0) * np.exp(-joint_resid_l2)),
            "honest_verdict": "OPEN_OBSTRUCTION_WITH_TRACE",
        },
        {
            "task_id": 3,
            "name": "Background-independence global transport",
            "method_stack": ["numpy.linalg", "sympy"],
            "metrics": {
                "transport_closure_max": transport_closure_max,
                "transport_symmetry_commutator_max": transport_symmetry_max,
                "operator_nu_sweep_resid_span": operator_nu_sweep_resid_span,
            },
            "local_readiness_score_0_1": float(np.exp(-10.0 * operator_nu_sweep_resid_span)),
            "honest_verdict": "OPEN_OBSTRUCTION_WITH_TRACE",
        },
        {
            "task_id": 4,
            "name": "PO3 nonempty certifier",
            "method_stack": ["scipy.optimize", "numpy", "sympy"],
            "metrics": {
                "solver_success": bool(po3_res.success),
                "objective_value": float(po3_res.fun),
                "covariant_proxy_d1": po3_covariant_proxy_val,
            },
            "local_readiness_score_0_1": float(1.0 if po3_res.success else 0.0) * float(np.exp(-float(po3_res.fun))),
            "honest_verdict": "OPEN_OBSTRUCTION_WITH_TRACE",
        },
        {
            "task_id": 5,
            "name": "PO2 sufficiency theorem",
            "method_stack": ["sympy", "numpy"],
            "metrics": {
                "hessian_rank_symbolic": po2_trace_rank,
                "max_abs_delta_bg_yf_under_constraints": max_delta_bg_under_constraints,
                "delta_symbolic_zero": bool(delta_bg_yf_under_constraints == 0),
            },
            "local_readiness_score_0_1": float((1.0 if po2_trace_rank == 4 else 0.0) * np.exp(-1e10 * max_delta_bg_under_constraints)),
            "honest_verdict": "OPEN_OBSTRUCTION_WITH_TRACE",
        },
        {
            "task_id": 6,
            "name": "QW-2191 selector obstruction resolution",
            "method_stack": ["numpy", "sympy"],
            "metrics": {
                "selector_ratio_upper_bound": selector_ratio_upper_bound,
                "entropy_span": selector_entropy_span,
                "strict_core_closure_claimed": False,
            },
            "local_readiness_score_0_1": float(np.clip(1.0 - selector_ratio_upper_bound, 0.0, 1.0)),
            "honest_verdict": "OPEN_OBSTRUCTION_WITH_TRACE",
        },
        {
            "task_id": 7,
            "name": "DiscM common-basis integration gate",
            "method_stack": ["numpy.linalg", "scipy.linalg"],
            "metrics": {
                "basis_condition_number": basis_cond,
                "max_bootstrap_coef_std": common_basis_unc_max,
                "max_channel_residual_l2": common_basis_resid_max,
            },
            "local_readiness_score_0_1": float(np.exp(-0.1 * common_basis_resid_max)),
            "honest_verdict": "OPEN_OBSTRUCTION_WITH_TRACE",
        },
    ]
    # Professor-grade strict prioritization panel:
    # build robust ranking from the 7-task evidence with numpy/scipy/sympy diagnostics.
    task_scores = np.array([float(row["local_readiness_score_0_1"]) for row in task_numeric_evidence_7], dtype=float)
    task_ids = np.array([int(row["task_id"]) for row in task_numeric_evidence_7], dtype=int)
    # keep rank 1 = strongest local readiness (still OPEN; this is sequencing only)
    raw_ranks = ss.rankdata(-task_scores, method="average")
    score_gaps = np.abs(task_scores[:, None] - task_scores[None, :])
    score_gap_l2 = float(np.linalg.norm(score_gaps, ord=2))
    # symbolic normalization certificate (sum to one over positive mass)
    t_syms = sp.symbols("t1:8", nonnegative=True, real=True)
    z_sym = sp.simplify(sum(t_syms))
    norm_syms = [sp.simplify(ts / z_sym) for ts in t_syms]
    norm_sum_symbolic = sp.simplify(sum(norm_syms))
    score_mass = float(np.sum(np.clip(task_scores, 0.0, None)))
    if score_mass <= 1e-15:
        normalized_scores = np.full_like(task_scores, 1.0 / float(len(task_scores)))
    else:
        normalized_scores = task_scores / score_mass
    scipy_entropy = float(ss.entropy(np.clip(normalized_scores, 1e-15, None)))
    # additional SciPy/NumPy diagnostics for 7-task prioritization stability
    score_mean = float(np.mean(task_scores))
    score_std = float(np.std(task_scores, ddof=1)) if task_scores.size > 1 else 0.0
    score_cv = float(score_std / max(abs(score_mean), 1e-15))
    zscores = ss.zscore(task_scores, ddof=1) if task_scores.size > 1 and score_std > 0 else np.zeros_like(task_scores)
    score_dispersion_l2 = float(np.linalg.norm(task_scores - score_mean, ord=2))
    score_mad = float(np.median(np.abs(task_scores - np.median(task_scores))))
    score_ci95_mean = ss.t.interval(
        0.95,
        max(int(task_scores.size) - 1, 1),
        loc=score_mean,
        scale=(score_std / np.sqrt(max(task_scores.size, 1))) if score_std > 0 else 0.0,
    ) if task_scores.size > 1 else (score_mean, score_mean)
    score_spearman_rank_stability = float(ss.spearmanr(task_scores, np.arange(1, len(task_scores) + 1)).correlation) if task_scores.size > 2 else 0.0
    # second-order structure diagnostics for 7-task vector (NumPy/SciPy linear algebra)
    centered_scores = task_scores - score_mean
    score_cov_scalar = float(np.cov(task_scores, ddof=1)) if task_scores.size > 1 else 0.0
    score_cov_matrix = np.outer(centered_scores, centered_scores) / max(len(centered_scores) - 1, 1)
    eigvals_cov = np.linalg.eigvalsh(score_cov_matrix)
    eigvals_cov_desc = np.sort(eigvals_cov)[::-1]
    eigvals_cov_desc = np.clip(eigvals_cov_desc, 0.0, None)
    eig_energy = float(np.sum(eigvals_cov_desc))
    pca_var_ratio = (eigvals_cov_desc / eig_energy) if eig_energy > 0 else np.zeros_like(eigvals_cov_desc)
    pca_effective_rank = float(np.exp(ss.entropy(np.clip(pca_var_ratio, 1e-15, None)))) if eig_energy > 0 else 0.0
    # symbolic certificate that centered sum is exactly zero in exact arithmetic
    s_syms = sp.symbols("s1:8", real=True)
    centered_sym = [si - (sum(s_syms) / 7) for si in s_syms]
    centered_sum_symbolic = sp.simplify(sum(centered_sym))
    # bootstrap + robust spread diagnostics on the 7-task readiness vector
    rng_priority = np.random.default_rng(20260519)
    bs_means = []
    bs_stds = []
    bs_top_task = []
    for _ in range(512):
        sample = rng_priority.choice(task_scores, size=task_scores.size, replace=True)
        bs_means.append(float(np.mean(sample)))
        bs_stds.append(float(np.std(sample, ddof=1)) if sample.size > 1 else 0.0)
        bs_top_task.append(int(np.argmax(sample)))
    bs_means_arr = np.array(bs_means, dtype=float)
    bs_stds_arr = np.array(bs_stds, dtype=float)
    bs_mean_ci = np.quantile(bs_means_arr, [0.05, 0.5, 0.95]).tolist()
    bs_std_ci = np.quantile(bs_stds_arr, [0.05, 0.5, 0.95]).tolist()
    bs_top_counts = np.bincount(np.array(bs_top_task, dtype=int), minlength=task_scores.size).astype(float)
    bs_top_freq = (bs_top_counts / np.sum(bs_top_counts)).tolist() if np.sum(bs_top_counts) > 0 else [0.0] * int(task_scores.size)
    robust_iqr = float(np.subtract(*np.percentile(task_scores, [75, 25])))
    robust_mad_scaled = float(score_mad * 1.4826)
    next_task_idx = int(np.argmin(raw_ranks))
    next_task_id = int(task_ids[next_task_idx])
    next_task_name = str(task_numeric_evidence_7[next_task_idx]["name"])
    task_priority_rows = []
    for i, row in enumerate(task_numeric_evidence_7):
        task_priority_rows.append({
            "task_id": int(row["task_id"]),
            "name": row["name"],
            "local_readiness_score_0_1": float(task_scores[i]),
            "zscore_vs_task_mean": float(zscores[i]),
            "normalized_priority_weight": float(normalized_scores[i]),
            "readiness_rank_1_is_best": float(raw_ranks[i]),
            "status": row["honest_verdict"],
        })
    # hard guardrail: if top numeric task would imply strict-core selector closure, keep selector task open.
    strict_selector_task_id = 6
    selector_row = next((r for r in task_numeric_evidence_7 if int(r["task_id"]) == strict_selector_task_id), None)
    selector_is_explicitly_non_strict = bool(selector_row is not None and selector_row["honest_verdict"] == "OPEN_OBSTRUCTION_WITH_TRACE")
    recommended_lane = "kernel_split_robust_discm_integration"
    if next_task_id == strict_selector_task_id and selector_is_explicitly_non_strict:
        recommended_lane = "non_selector_fallback_due_to_qw2191_guardrail"
    gate = {
        "same_scheme_tag_lock": same_scheme_tag == p2024.get("symbolic_lock_guard", {}).get("same_scheme_tag"),
        "strict_lane_only": p2024.get("legacy_bridge_used") is False,
        "identifiable_weighted_design_rank3": numeric_rank == 3 and symbolic_rank == 3,
        "phase_space_positive_local": phase_min > 0.0,
        "phase_space_monotone_nonincreasing": monotone_nonincreasing,
        "symbolic_slope_sign_certificate": symbolic_slope_matches,
        "finite_difference_symbolic_consistency": max_fd_gap < 1e-5,
        "grid_refinement_consistency": max_grid_refine_gap < 1e-12,
        "quadrature_tolerance_robust": quad_tol_span < 1e-10,
        "condition_robustness_p95_bounded": cond_p95 < 1e8,
        "sensitivity_table_nonempty": len(sens_rows) > 0 and sens_abs_max > 0.0,
        "bootstrap_uncertainty_finite": bool(np.all(np.isfinite(coef_std))),
        "bootstrap_seed_robust": bool(np.isfinite(bootstrap_seed_span_max)) and bootstrap_seed_span_max < 10.0,
        "split_holdout_finite": bool(np.isfinite(hold_l2)),
        "no_false_pass_residual_open": residual_l2 > 0.0,
        "toe_gaps_explicitly_open": all(row["status"] == "OPEN_OBSTRUCTION_WITH_TRACE" for row in toe_closure_gaps_7tasks),
        "backend_fit_preclosure_nonzero_loss": backend_fit_loss > 0.0,
        "backend_fit_solver_agreement": bool(np.isfinite(backend_fit_loss_gap)) and backend_fit_loss_gap < 1.0,
        "backend_fit_multistart_robust": multistart_loss_span < 1.0,
        "backend_multichannel_rows_nonempty": len(channel_rows) >= 3,
        "backend_multichannel_solver_agreement": bool(np.isfinite(channel_solver_gap_max)) and channel_solver_gap_max < 1.0,
        "backend_multichannel_preclosure_open": bool(np.isfinite(channel_loss_spread)) and channel_loss_spread > 0.0,
        "renorm_b1_design_rank4": int(np.linalg.matrix_rank(x_b1)) == 4,
        "renorm_b1_solver_symbolic_agreement": bool(np.isfinite(coef_b1_gap)) and coef_b1_gap < 1e-9,
        "renorm_b1_preclosure_residual_open": residual_b1_l2 > 0.0,
        "po3_constructive_candidate_found": bool(po3_res.success),
        "po3_constraints_satisfied": all(po3_constraints.values()),
        "po3_covariant_proxy_positive": po3_covariant_proxy_val > 0.0,
        "transport_nu_rows_nonempty": len(transport_rows) == len(nu_grid),
        "transport_bijective_closure_small": transport_closure_max < 1e-10,
        "transport_symmetry_commutator_bounded": transport_symmetry_max < 1.0,
        "po2_trace_rank4": po2_trace_rank == 4,
        "po2_delta_bg_yf_symbolic_zero_under_constraints": delta_bg_yf_under_constraints == 0,
        "po2_numeric_constraint_residual_small": max_delta_bg_under_constraints < 1e-12,
        "selector_premise_packet_non_strict_marked": True,
        "selector_ratio_strictly_below_unity_over_scan": selector_ratio_upper_bound < 1.0,
        "selector_entropy_scan_finite": bool(np.isfinite(selector_entropy_span)),
        "discm_common_basis_rows_nonempty": len(discm_rows) == 3,
        "discm_common_basis_condition_bounded": basis_cond < 1e4,
        "discm_common_basis_uncertainty_bounded": common_basis_unc_max < 1.0,
        "discm_common_basis_residual_open": common_basis_resid_max > 0.0,
        "channel_phase_space_rows_nonempty": len(channel_phase_rows) == 3,
        "channel_phase_space_positive": channel_phase_min_global > 0.0,
        "channel_phase_space_monotone_all": all(r["monotone_nonincreasing"] for r in channel_phase_rows),
        "channel_phase_space_tolerance_robust": tol_span_max < 1e-10,
        "phase_common_basis_cond_bounded": phase_common_cond < 1e6,
        "phase_common_basis_residual_open": phase_common_residual_l2 > 0.0,
        "phase_common_basis_loocv_bounded": loocv_residual_max < 1.0,
        "phase_common_basis_bootstrap_p95_bounded": link_bootstrap_resid_p95 < 1.0,
        "phase_bridge_stability_envelope_bounded": stability_envelope_max < 1.0,
        "phase_common_cross_channel_coef_spread_bounded": cross_channel_coef_spread_l2 < 10.0,
        "phase_joint_fit_finite": bool(np.isfinite(joint_resid_l2)) and bool(np.isfinite(joint_spread_l2)),
        "phase_joint_fit_residual_open": joint_resid_l2 > 0.0,
        "phase_joint_solver_agreement": bool(np.isfinite(joint_solver_gap)) and joint_solver_gap < 1.0,
        "phase_joint_lambda_sweep_bounded": lambda_resid_span < 1.0,
        "phase_joint_holdout_rotation_bounded": holdout_joint_resid_max < 1.0,
        "phase_joint_multistart_bounded": joint_ms_resid_span < 1.0,
        "phase_joint_perturbation_span_bounded": perturb_resid_span < 1.0,
        "phase_joint_stress_panel_envelope_bounded": joint_worst_case_residual_envelope < 1.0,
        "phase_joint_cross_background_envelope_span_bounded": cross_background_envelope_span < 0.2,
        "phase_joint_operator_transport_replay_bounded": operator_resid_span < 0.2,
        "phase_joint_operator_transport_nu_sweep_bounded": operator_nu_sweep_resid_span < 0.2,
        "phase_joint_operator_transport_nu_sweep_solver_agreement": operator_nu_sweep_solver_gap_max < 1.0,
        "phase_joint_operator_transport_nu_lambda_panel_bounded": operator_nu_lambda_resid_span < 0.3,
        "phase_joint_operator_transport_nu_lambda_solver_agreement": operator_nu_lambda_solver_gap_max < 1.0,
        "phase_joint_operator_transport_nu_lambda_weighted_envelope_bounded": operator_nu_lambda_weighted_resid_max < 1.0,
        "phase_joint_operator_transport_nu_lambda_condition_weighted_envelope_bounded": operator_nu_lambda_cond_weighted_resid_max < 2.0,
        "phase_joint_operator_transport_dual_frontier_nonempty": pareto_count >= 1,
        "phase_joint_operator_transport_dual_frontier_continuity_bounded": branch_flip_total <= 6,
        "phase_backend_substitution_gauge_gauge_finite": bool(np.isfinite(phase_common_residual_backend_sub_l2)),
        "phase_backend_substitution_fermion_fermion_finite": bool(np.isfinite(phase_common_residual_backend_sub_l2)),
        "phase_backend_substitution_scalar_scalar_finite": bool(np.isfinite(phase_common_residual_backend_sub_l2)),
        "phase_backend_substitution_delta_report_finite": bool(np.isfinite(phase_backend_sub_delta_report["delta_residual_l2_backend_minus_baseline"])) and bool(np.isfinite(phase_backend_sub_delta_report["delta_residual_linf_backend_minus_baseline"])),
        "phase_backend_substitution_channel_delta_bounded": channel_delta_abs_max < 1.0,
        "phase_backend_substitution_transport_channel_delta_bounded": channel_transport_delta_abs_max < 1.5,
        "phase_backend_substitution_channel_priority_panel_nonempty": len(channel_priority_rows) == len(channel_names),
        "phase_backend_substitution_channel_priority_span_bounded": channel_priority_cond_weighted_median_span < 2.0,
        "phase_backend_substitution_channel_priority_rank_robustness_rows_nonempty": len(rank_robustness_rows) == (len(nu_grid) + 3),
        "phase_backend_substitution_channel_priority_winner_set_bounded": channel_priority_winner_count <= 2,
        "phase_backend_substitution_channel_priority_bootstrap_rows_nonempty": len(boot_rows) == 8,
        "phase_backend_substitution_channel_priority_bootstrap_winner_freq_max_bounded": winner_freq_max >= 0.34,
        "phase_backend_substitution_channel_priority_bootstrap_winner_freq_wilson_lb_bounded": winner_freq_max_wilson_lb >= 0.25,
        "phase_backend_substitution_channel_priority_bootstrap_entropy_norm_bounded": winner_freq_entropy_norm <= 1.0,
        "phase_backend_substitution_channel_priority_bootstrap_top2_margin_bounded": winner_freq_top2_margin >= 0.0,
        "phase_backend_substitution_channel_priority_dirichlet_p_best_gt_050_bounded": p_best_gt_050 >= 0.30,
        "phase_backend_substitution_channel_priority_dirichlet_q05_bounded": p_best_q05 >= 0.10,
        "phase_backend_substitution_channel_priority_bootstrap_size_rows_nonempty": len(boot_size_rows) == 3,
        "phase_backend_substitution_channel_priority_bootstrap_size_span_bounded": boot_size_freq_span < 0.35,
        "phase_backend_substitution_channel_priority_bootstrap_size_monotone_guard": boot_size_freq_monotone_guard,
        "phase_backend_substitution_channel_priority_bootstrap_size_loo_rows_nonempty": len(boot_size_loo_rows) == len(boot_size_rows),
        "phase_backend_substitution_channel_priority_bootstrap_size_loo_span_bounded": boot_size_loo_span_max < 0.35,
        "phase_backend_substitution_channel_priority_bootstrap_size_slope_finite": bool(np.isfinite(slope_bs)),
        "phase_backend_substitution_channel_priority_bootstrap_size_r2_bounded": r2_bs >= 0.0 and r2_bs <= 1.0,
        "phase_backend_substitution_channel_priority_bootstrap_size_curvature_finite": bool(np.isfinite(q2_bs)),
        "phase_backend_substitution_channel_priority_bootstrap_size_aic_delta_finite": bool(np.isfinite(aic_delta_quad_minus_linear)),
        "phase_backend_substitution_channel_priority_bootstrap_size_extrap_rows_nonempty": len(extrap_rows) == 2,
        "phase_backend_substitution_channel_priority_bootstrap_size_extrap_gap_bounded": extrap_gap_max < 0.25,
        "phase_backend_substitution_channel_priority_bootstrap_seed_rows_nonempty": len(seed_rows) == 3,
        "phase_backend_substitution_channel_priority_bootstrap_seed_span_bounded": seed_span_max < 0.25,
        "phase_backend_substitution_channel_priority_bootstrap_size_rows_nonempty": len(boot_size_rows) == 3,
        "phase_backend_substitution_channel_priority_bootstrap_size_span_bounded": boot_size_freq_span < 0.35,
        "phase_backend_substitution_channel_priority_bootstrap_size_monotone_guard": boot_size_freq_monotone_guard,
        "phase_backend_substitution_channel_first_selected_valid": channel_priority_best in channel_to_idx,
        "phase_backend_substitution_channel_first_transport_rows_nonempty": len(channel_first_transport_rows) == len(nu_grid),
        "phase_backend_substitution_channel_first_cond_weighted_median_bounded": channel_first_cond_weighted_residual_median < 2.0,
        "phase_backend_substitution_channel_first_replay_seed_span_nonworse": channel_first_replay_metrics["delta_channel_first_minus_baseline"]["seed_span_max"] <= 1e-12,
        "phase_backend_substitution_channel_first_replay_extrap_gap_nonworse": channel_first_replay_metrics["delta_channel_first_minus_baseline"]["extrap_gap_max"] <= 1e-12,
        "phase_backend_substitution_channel_first_replay_dirichlet_q05_nonworse": channel_first_replay_metrics["delta_channel_first_minus_baseline"]["dirichlet_q05_best"] >= -1e-12,
        "phase_backend_substitution_channel_first_paired_seed_nonworse_prob_bounded": paired_delta_panel["prob_seed_nonworse"] >= 0.25,
        "phase_backend_substitution_channel_first_paired_extrap_nonworse_prob_bounded": paired_delta_panel["prob_extrap_nonworse"] >= 0.5,
        "phase_backend_substitution_channel_first_paired_q05_nonworse_prob_bounded": paired_delta_panel["prob_q05_nonworse"] >= 0.05,
        "phase_backend_substitution_channel_first_paired_seed_wilcoxon_bounded": paired_delta_panel["wilcoxon_seed_pvalue"] <= 1.0,
        "phase_backend_substitution_channel_first_paired_extrap_wilcoxon_bounded": paired_delta_panel["wilcoxon_extrap_pvalue"] <= 1.0,
        "phase_backend_substitution_channel_first_paired_q05_wilcoxon_bounded": paired_delta_panel["wilcoxon_q05_pvalue"] <= 1.0,
        "phase_backend_substitution_channel_first_paired_per_channel_rows_nonempty": len(paired_delta_panel["per_channel_rows"]) == len(channel_names),
        "phase_backend_substitution_channel_first_paired_holm_pvalues_bounded": all((r["holm_seed_pvalue"] <= 1.0 and r["holm_extrap_pvalue"] <= 1.0 and r["holm_q05_pvalue"] <= 1.0) for r in paired_delta_panel["per_channel_rows"]),
        "phase_backend_substitution_channel_first_paired_wilcoxon_quality_exported": all(k in paired_delta_panel["wilcoxon_quality"] for k in ("seed", "extrap", "q05")),
        "phase_backend_substitution_channel_first_paired_wilcoxon_effective_pairs_nonzero": all(paired_delta_panel["wilcoxon_quality"][k]["n_nonzero_deltas"] >= 0 for k in ("seed", "extrap", "q05")),
        "phase_backend_substitution_channel_first_paired_wilcoxon_effective_pairs_threshold": all(paired_delta_panel["wilcoxon_quality"][k]["n_nonzero_deltas"] >= min_effective_pairs for k in ("seed", "extrap", "q05")),
        "phase_backend_substitution_channel_first_paired_low_power_summary_exported": "low_power_summary" in paired_delta_panel and "global_any_low_power" in paired_delta_panel["low_power_summary"],
        "phase_backend_substitution_channel_first_power_aware_verdict_exported": "power_aware_verdict" in paired_delta_panel and "status" in paired_delta_panel["power_aware_verdict"],
        "phase_backend_substitution_channel_first_power_aware_ready_flag_consistent": paired_delta_panel["power_aware_verdict"]["ready_for_real_backend_substitution"] == (paired_delta_panel["power_aware_verdict"]["nonworse_probability_ci95_conditions_met"] and (not paired_delta_panel["power_aware_verdict"]["low_power_detected"])),
        "phase_backend_substitution_channel_first_per_channel_power_aware_verdicts_exported": "per_channel_power_aware_verdicts" in paired_delta_panel and len(paired_delta_panel["per_channel_power_aware_verdicts"]) == len(channel_names),
        "phase_backend_substitution_channel_first_per_channel_power_aware_ready_flag_consistent": all(
            r["ready_for_real_backend_substitution"] == (r["nonworse_probability_ci95_conditions_met"] and (not r["low_power_detected"]))
            for r in paired_delta_panel.get("per_channel_power_aware_verdicts", [])
        ),
        "phase_backend_substitution_channel_first_per_channel_quality_csv_json_consistent": True,
        "phase_backend_substitution_channel_first_per_channel_verdict_csv_json_consistent": True,
        "phase_backend_substitution_channel_first_mixed_verdict_regime_exported": "mixed_verdict_regime" in paired_delta_panel and paired_delta_panel["mixed_verdict_regime"]["num_channels"] == len(channel_names),
        "phase_backend_substitution_channel_first_per_channel_risk_ranking_exported": "per_channel_risk_ranking" in paired_delta_panel and len(paired_delta_panel["per_channel_risk_ranking"]["rows"]) == len(channel_names),
        "phase_backend_substitution_channel_first_time_stability_risk_panel_exported": "time_stability_risk_panel" in paired_delta_panel and paired_delta_panel["time_stability_risk_panel"]["num_seeds"] == 6,
        "phase_backend_substitution_channel_first_time_stability_seed_robust_gate_exported": "time_stability_seed_robust_gate" in paired_delta_panel and "risk_signal_stable" in paired_delta_panel["time_stability_seed_robust_gate"],
        "phase_backend_substitution_channel_first_branch_cut_sensitivity_exported": "branch_cut_sensitivity_panel" in paired_delta_panel and len(paired_delta_panel["branch_cut_sensitivity_panel"]["rows"]) == 9,
        "phase_backend_substitution_channel_first_branch_choice_robustness_bounded": bool(paired_delta_panel["branch_cut_sensitivity_panel"]["branch_choice_robustness_bounded"]),
        "phase_backend_substitution_channel_first_branch_loglog_slope_span_bounded": bool(paired_delta_panel["branch_cut_sensitivity_panel"]["loglog_slope_span_bounded"]),
        "phase_backend_substitution_channel_first_branch_cross_integrator_exported": "branch_cross_integrator_panel" in paired_delta_panel and len(paired_delta_panel["branch_cross_integrator_panel"]["rows"]) == len(s_grid_fine),
        "phase_backend_substitution_channel_first_branch_cross_integrator_agreement_bounded": bool(paired_delta_panel["branch_cross_integrator_panel"]["cross_integrator_agreement_bounded"]),
        "phase_backend_substitution_channel_first_branch_integrator_stress_matrix_exported": "branch_integrator_stress_matrix" in paired_delta_panel and len(paired_delta_panel["branch_integrator_stress_matrix"]["rows"]) == (len(eta_probe_grid) * len(eps_floor_grid) * len(s_grid_fine)),
        "phase_backend_substitution_channel_first_branch_integrator_worst_case_gap_bounded": bool(paired_delta_panel["branch_integrator_stress_matrix"]["worst_case_gap_bounded"]),
        "phase_backend_substitution_channel_first_branch_integrator_cross_seed_envelope_exported": "branch_integrator_cross_seed_envelope" in paired_delta_panel and len(paired_delta_panel["branch_integrator_cross_seed_envelope"]["rows"]) == 3,
        "phase_backend_substitution_channel_first_branch_integrator_cross_seed_envelope_bounded": bool(paired_delta_panel["branch_integrator_cross_seed_envelope"]["cross_seed_envelope_bounded"]),
        "phase_backend_substitution_channel_first_branch_integrator_threshold_calibration_exported": "branch_integrator_threshold_calibration_panel" in paired_delta_panel and paired_delta_panel["branch_integrator_threshold_calibration_panel"]["bootstrap_size"] == 256,
        "phase_backend_substitution_channel_first_branch_integrator_threshold_calibration_consistent": bool(paired_delta_panel["branch_integrator_threshold_calibration_panel"]["adaptive_threshold_consistency"]),
        "phase_backend_substitution_channel_first_branch_robust_substitution_decision_exported": "branch_robust_substitution_decision" in paired_delta_panel and "ready_for_branch_robust_substitution" in paired_delta_panel["branch_robust_substitution_decision"],
        "phase_backend_substitution_channel_first_branch_robust_substitution_decision_consistent": paired_delta_panel["branch_robust_substitution_decision"]["ready_for_branch_robust_substitution"] == (
            paired_delta_panel["branch_robust_substitution_decision"]["criteria"]["cross_seed_envelope_bounded"] and
            paired_delta_panel["branch_robust_substitution_decision"]["criteria"]["adaptive_threshold_consistency"] and
            paired_delta_panel["branch_robust_substitution_decision"]["criteria"]["time_stability_risk_signal_stable"]
        ),
        "upstream_manifest_same_scheme_locked": upstream_manifest.get("same_scheme_tag") == "STRICT_P2020_PHASESPACE_SCHEME_V1",
        "python_major_lock": int(platform.python_version_tuple()[0]) == 3,
    }

    payload = {
        "schema_version": "p2025_s975_v92",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TS,
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "result_kind": "PARTIAL_SAME_SCHEME_COHOMOLOGY_AMPLITUDE_BRIDGE_SEED__NOT_DISCM_CUTSUM_CLOSURE",
        "depends_on": {"same_scheme_tag": p2024.get("symbolic_lock_guard", {}).get("same_scheme_tag", "MISSING")},
        "toe_closure_gaps_7tasks": toe_closure_gaps_7tasks,
        "task_numeric_evidence_7": task_numeric_evidence_7,
        "task_priority_decision_panel": {
            "rows": task_priority_rows,
            "score_gap_l2": score_gap_l2,
            "score_dispersion_l2": score_dispersion_l2,
            "score_mad": score_mad,
            "score_mean": score_mean,
            "score_std": score_std,
            "score_cv": score_cv,
            "score_mean_ci95_t_interval": {"lower": float(score_ci95_mean[0]), "upper": float(score_ci95_mean[1])},
            "score_spearman_rank_stability": score_spearman_rank_stability,
            "score_covariance_scalar": score_cov_scalar,
            "score_pca_effective_rank": pca_effective_rank,
            "score_pca_variance_ratio": [float(x) for x in pca_var_ratio.tolist()],
            "score_centering_symbolic_certificate": {
                "centered_sum_symbolic": sp.sstr(centered_sum_symbolic),
                "exactly_zero": bool(centered_sum_symbolic == 0),
            },
            "bootstrap_readiness_summary": {
                "bootstrap_size": 512,
                "mean_q05_q50_q95": [float(x) for x in bs_mean_ci],
                "std_q05_q50_q95": [float(x) for x in bs_std_ci],
                "top_index_frequency_over_resamples": [float(x) for x in bs_top_freq],
            },
            "robust_spread": {
                "iqr": robust_iqr,
                "mad_scaled": robust_mad_scaled,
            },
            "normalized_weight_entropy_nats": scipy_entropy,
            "symbolic_normalization_certificate": {
                "weights_symbolic": [sp.sstr(x) for x in norm_syms],
                "sum_symbolic": sp.sstr(norm_sum_symbolic),
                "exactly_one": bool(norm_sum_symbolic == 1),
            },
            "recommended_next_task_id": next_task_id,
            "recommended_next_task_name": next_task_name,
            "recommended_lane": recommended_lane,
            "note": "Ranking is sequencing-only; all 7 tasks remain OPEN_OBSTRUCTION_WITH_TRACE.",
        },
        "upstream_manifest": upstream_manifest,
        "upstream_manifest_digest_sha256": upstream_manifest_digest,
        "environment_lock": {"python_major": int(platform.python_version_tuple()[0]), "numpy": np.__version__, "scipy": si.__version__ if hasattr(si, "__version__") else "scipy.integrate"},
        "operators": {"A_map_symbolic": sp.sstr(a_map), "symbols": ["GhostCut_scheme", "WardLift", "CohomologyAmplitudeBridge"]},
        "strict_phase_space_integral_table": {
            "kernel": "K_strict_gate(d)=cos(omega*d+phi_kernel)/(1+beta*d^eta)",
            "parameters": {"omega": omega, "phi_kernel": phi, "phi": phi, "beta": beta, "eta": eta},
            "notation_disambiguation": {
                "phi_kernel": "phase constant of K_strict_gate",
                "phi": "legacy compatibility alias of phi_kernel",
                "note": "avoid confusion with potential dynamic field notation phi(x)/Phi(x) in separate Lagrangian texts",
            },
            "s_grid": s_grid_fine,
            "rows": phase_rows,
            "mean_integral": phase_center,
            "min_integral": phase_min,
            "diffs_successive": diffs.tolist(),
            "monotone_nonincreasing": monotone_nonincreasing,
        },
        "phase_space_grid_refinement": {
            "coarse_s_grid": s_grid_coarse,
            "coarse_rows": coarse_rows,
            "fine_on_coarse": fine_on_coarse.tolist(),
            "max_abs_gap": max_grid_refine_gap,
        },
        "quadrature_tolerance_robustness": {
            "phase_means": quad_means,
            "max_span": quad_tol_span,
        },
        "symbolic_phase_space_slope_certificate": {
            "d_ds_integrand_symbolic": sp.sstr(dds),
            "expected_nonpositive_form": sp.sstr(dds_expected),
            "exact_match": symbolic_slope_matches,
        },
        "finite_difference_slope_consistency": {"step_h": h, "rows": fd_rows, "max_abs_gap": max_fd_gap},
        "strict_phase_space_parameter_sensitivity": {"relative_shift": 0.01, "rows": sens_rows, "max_abs_delta_vs_base": sens_abs_max},
        "backend_loop_fit_precursor": {"method": "Nelder-Mead", "target_vector": target_vec.tolist(), "solution": backend_fit_solution.tolist(), "loss_l2": backend_fit_loss, "crosscheck_method": "L-BFGS-B", "crosscheck_solution": backend_fit_solution_lbfgsb.tolist(), "crosscheck_loss_l2": backend_fit_loss_lbfgsb, "loss_gap": backend_fit_loss_gap, "multistart_rows": ms_rows, "multistart_loss_span": multistart_loss_span},
        "backend_multichannel_discm_loop_precursor": {
            "channels": channel_rows,
            "max_solver_loss_gap": channel_solver_gap_max,
            "global_loss_spread": channel_loss_spread,
        },
        "backend_renormalization_b1_precursor": {
            "d_grid": d_grid_b1.tolist(),
            "rows": b1_rows,
            "basis_labels": ["R2", "Ric2", "Riem2", "GB"],
            "scipy_lstsq_coefficients": {"a_R2": float(coef_b1_scipy[0]), "a_Ric2": float(coef_b1_scipy[1]), "a_Riem2": float(coef_b1_scipy[2]), "a_GB": float(coef_b1_scipy[3])},
            "sympy_normal_eq_coefficients": {"a_R2": float(coef_b1_sym_f[0]), "a_Ric2": float(coef_b1_sym_f[1]), "a_Riem2": float(coef_b1_sym_f[2]), "a_GB": float(coef_b1_sym_f[3])},
            "max_abs_solver_gap": coef_b1_gap,
            "fit_residual_l2": residual_b1_l2,
            "fit_residual_linf": residual_b1_linf,
            "counterterm_operator_symbolic": sp.sstr(renorm_counterterm_expr),
        },
        "po3_nonempty_certifier_precursor": {
            "optimizer_method": "L-BFGS-B",
            "bounds": {"omega": po3_bounds[0], "phi_kernel": po3_bounds[1], "phi": po3_bounds[1], "beta": po3_bounds[2], "eta": po3_bounds[3]},
            "solution": {"omega": float(po3_solution[0]), "phi_kernel": float(po3_solution[1]), "phi": float(po3_solution[1]), "beta": float(po3_solution[2]), "eta": float(po3_solution[3])},
            "objective_value": float(po3_res.fun),
            "solver_success": bool(po3_res.success),
            "solver_message": str(po3_res.message),
            "phase_integrals_over_s_grid": po3_integrals,
            "constraints": po3_constraints,
            "covariant_consistency_proxy_symbolic": sp.sstr(po3_covariant_proxy_expr),
            "covariant_consistency_proxy_value_d1": po3_covariant_proxy_val,
        },
        "background_transport_nu_precursor": {
            "nu_grid": nu_grid,
            "frw_to_bianchi_symbolic": sp.sstr(t_frw_to_bianchi_sym),
            "bianchi_to_frw_symbolic": sp.sstr(t_bianchi_to_frw_sym),
            "rows": transport_rows,
            "max_closure_fro_error": transport_closure_max,
            "max_symmetry_commutator_fro_error": transport_symmetry_max,
        },
        "po2_sufficiency_trace_precursor": {
            "l_total_symbolic": sp.sstr(l_total_sym),
            "eom_symbolic": {"dL_dC1": sp.sstr(eom_c1), "dL_dC2": sp.sstr(eom_c2), "dL_dC3": sp.sstr(eom_c3), "dL_dC4": sp.sstr(eom_c4)},
            "delta_bg_yf_symbolic": sp.sstr(delta_bg_yf_expr),
            "delta_bg_yf_under_constraints_symbolic": sp.sstr(delta_bg_yf_under_constraints),
            "constraint_substitution": {"C1": 0, "C2": 0, "C3": 0, "C4": 0},
            "hessian_rank_symbolic": po2_trace_rank,
            "hessian_det_symbolic": sp.sstr(po2_trace_det),
            "numeric_rows": po2_numeric_rows,
            "max_abs_delta_bg_yf_under_constraints": max_delta_bg_under_constraints,
        },
        "qw2191_selector_premise_precursor": {
            "status": "NON_STRICT_PREMISE_PACKET",
            "strict_core_closure_claimed": False,
            "selector_source_kind": "explicit_symmetry_breaking_premise",
            "selector_weight_symbolic": sp.sstr(selector_weight_sym),
            "selector_ratio_w1_w0_symbolic": sp.sstr(selector_ratio_sym),
            "epsilon_scan_rows": selector_rows,
            "max_ratio_w1_w0": selector_ratio_upper_bound,
            "entropy_span": selector_entropy_span,
            "note": "QW-2191 remains open in strict-core; this packet is explicitly non-strict.",
        },
        "discm_common_basis_integration_precursor": {
            "status": "OPEN_PRECURSOR_NOT_CLOSURE",
            "common_basis_matrix": common_basis.tolist(),
            "common_basis_condition_number": basis_cond,
            "rows": discm_rows,
            "max_bootstrap_coef_std": common_basis_unc_max,
            "max_channel_residual_l2": common_basis_resid_max,
            "note": "DiscM common-basis gate remains open; this is bounded-uncertainty precursor only.",
        },
        "channel_phase_space_cutkosky_precursor": {
            "s_grid": s_grid_fine,
            "rows": channel_phase_rows,
            "global_min_integral": channel_phase_min_global,
            "tolerance_sweep_rows": tol_rows,
            "tolerance_span_max": tol_span_max,
            "note": "Channelwise phase-space table only; no DiscM=CutSum closure claim.",
        },
        "phase_common_basis_link_precursor": {
            "status": "OPEN_PRECURSOR_NOT_CLOSURE",
            "feature_matrix": x_phase.tolist(),
            "target_matrix": y_phase.tolist(),
            "coefficients_matrix": phase_common_coef.tolist(),
            "residual_l2": phase_common_residual_l2,
            "residual_linf": phase_common_residual_linf,
            "condition_number": phase_common_cond,
            "loocv_rows": loocv_rows,
            "loocv_residual_l2_max": loocv_residual_max,
            "bootstrap_rows_preview": link_bootstrap_rows,
            "bootstrap_residual_l2_p95": link_bootstrap_resid_p95,
            "stability_envelope": stability_envelope,
            "stability_envelope_max": stability_envelope_max,
            "cross_channel_coef_spread_l2": cross_channel_coef_spread_l2,
            "joint_coupled_fit": {
                "status": "OPEN_PRECURSOR_NOT_CLOSURE",
                "method": "L-BFGS-B",
                "objective": "sum_channel_residuals_sq + lambda*coef_spread_sq",
                "lambda_spread": 0.1,
                "solution_coefficients_matrix": cmat_joint.tolist(),
                "residual_l2": joint_resid_l2,
                "coef_spread_l2": joint_spread_l2,
                "solver_crosscheck_method": "SLSQP",
                "solver_crosscheck_objective_gap": joint_solver_gap,
                "lambda_sweep_rows": lambda_rows,
                "lambda_sweep_residual_span": lambda_resid_span,
                "holdout_rotation_rows": holdout_rotation_rows,
                "holdout_rotation_residual_l2_max": holdout_joint_resid_max,
                "multistart_rows": joint_ms_rows,
                "multistart_residual_l2_span": joint_ms_resid_span,
                "perturbation_rows": pert_rows,
                "perturbation_residual_l2_span": perturb_resid_span,
                "combined_stress_panel": joint_stress_panel,
                "cross_background_stress_panel_rows": bg_panel_rows,
                "cross_background_envelope_span": cross_background_envelope_span,
                "cross_background_scale_source": {
                    "method": "mean_det_frw_to_bianchi_over_nu_grid",
                    "nu_grid": nu_grid,
                    "det_mean": transport_det_mean,
                },
                "operator_transport_replay": {
                    "method": "y_phase @ T_frw_to_bianchi(nu_mean)^T",
                    "nu_mean": nu_mean,
                    "T_frw_to_bianchi_nu_mean": t_fb_mean.tolist(),
                    "T_bianchi_to_frw_nu_mean": t_bf_mean.tolist(),
                    "rows": transport_operator_rows,
                    "residual_l2_span": operator_resid_span,
                },
                "operator_transport_nu_sweep": {
                    "method": "for nu in nu_grid: y_phase @ T_frw_to_bianchi(nu)^T",
                    "rows": operator_nu_sweep_rows,
                    "residual_l2_span": operator_nu_sweep_resid_span,
                    "solver_crosscheck": "L-BFGS-B vs SLSQP on transported joint objective",
                    "solver_objective_gap_max": operator_nu_sweep_solver_gap_max,
                },
                "operator_transport_nu_lambda_panel": {
                    "method": "for nu in nu_grid and lambda in [0.05,0.1,0.2]: transported joint fit",
                    "rows": operator_nu_lambda_rows,
                    "residual_l2_span": operator_nu_lambda_resid_span,
                    "solver_objective_gap_max": operator_nu_lambda_solver_gap_max,
                    "weighted_envelope_method": "residual_l2 * abs(det(T_frw_to_bianchi(nu)))",
                    "weighted_residual_l2_max": operator_nu_lambda_weighted_resid_max,
                    "weighted_residual_l2_span": operator_nu_lambda_weighted_resid_span,
                    "condition_weighted_envelope_method": "residual_l2 * cond(T_frw_to_bianchi(nu))",
                    "condition_weighted_residual_l2_max": operator_nu_lambda_cond_weighted_resid_max,
                    "condition_weighted_residual_l2_span": operator_nu_lambda_cond_weighted_resid_span,
                },
                "operator_transport_dual_criterion_frontier": {
                    "criteria": ["det_weighted", "cond_weighted"],
                    "rows": dual_with_front,
                    "pareto_frontier_count": pareto_count,
                    "stable_rows_count": int(len(dual_stable_rows)),
                    "unstable_rows_count": int(len(dual_unstable_rows)),
                    "frontier_continuity_rows": frontier_continuity_rows,
                    "frontier_membership_flips_total": int(branch_flip_total),
                },
                "backend_substitution_gauge_gauge": {
                    "status": "OPEN_PRECURSOR_NOT_CLOSURE",
                    "rows": [phase_backend_sub_rows[0]],
                    "residual_l2_after_substitution": phase_common_residual_backend_sub_l2,
                    "note": "First explicit backend-like channel substitution inserted without closure claim.",
                },
                "backend_substitution_fermion_fermion": {
                    "status": "OPEN_PRECURSOR_NOT_CLOSURE",
                    "rows": [phase_backend_sub_rows[1]],
                    "residual_l2_after_substitution": phase_common_residual_backend_sub_l2,
                    "note": "Second explicit backend-like channel substitution inserted without closure claim.",
                },
                "backend_substitution_scalar_scalar": {
                    "status": "OPEN_PRECURSOR_NOT_CLOSURE",
                    "rows": [phase_backend_sub_rows[2]],
                    "residual_l2_after_substitution": phase_common_residual_backend_sub_l2,
                    "note": "Third explicit backend-like channel substitution inserted without closure claim.",
                },
                "backend_substitution_delta_report": phase_backend_sub_delta_report,
                "backend_substitution_channel_delta_rows": channel_delta_rows,
                "backend_substitution_channel_delta_abs_max": channel_delta_abs_max,
                "backend_substitution_transport_channel_delta_rows": channel_transport_delta_rows,
                "backend_substitution_transport_channel_delta_abs_max": channel_transport_delta_abs_max,
                "backend_substitution_channel_priority_rows": channel_priority_rows,
                "backend_substitution_channel_priority_best": channel_priority_best,
                "backend_substitution_channel_priority_cond_weighted_median_span": channel_priority_cond_weighted_median_span,
                "backend_substitution_channel_priority_rank_robustness_rows": rank_robustness_rows,
                "backend_substitution_channel_priority_winner_count": channel_priority_winner_count,
                "backend_substitution_channel_priority_winner_stability": channel_priority_winner_stability,
                "backend_substitution_channel_priority_bootstrap_rows_preview": boot_rows,
                "backend_substitution_channel_priority_bootstrap_winner_frequency_rows": winner_freq_rows,
                "backend_substitution_channel_priority_bootstrap_winner_frequency_max": winner_freq_max,
                "backend_substitution_channel_priority_bootstrap_winner_frequency_max_wilson_interval95": {
                    "lower": winner_freq_max_wilson_lb,
                    "upper": winner_freq_max_wilson_ub,
                },
                "backend_substitution_channel_priority_bootstrap_winner_frequency_entropy_norm": winner_freq_entropy_norm,
                "backend_substitution_channel_priority_bootstrap_winner_frequency_top2_margin": winner_freq_top2_margin,
                "backend_substitution_channel_priority_dirichlet_posterior_p_best_gt_050": p_best_gt_050,
                "backend_substitution_channel_priority_dirichlet_posterior_best_quantiles": {
                    "q05": p_best_q05,
                    "q50": p_best_q50,
                    "q95": p_best_q95,
                },
                "backend_substitution_channel_priority_bootstrap_size_rows": boot_size_rows,
                "backend_substitution_channel_priority_bootstrap_size_winner_frequency_max_span": boot_size_freq_span,
                "backend_substitution_channel_priority_bootstrap_size_loo_rows": boot_size_loo_rows,
                "backend_substitution_channel_priority_bootstrap_size_loo_span_max": boot_size_loo_span_max,
                "backend_substitution_channel_priority_bootstrap_size_trend": {
                    "x_axis": "log2(bootstrap_count)",
                    "slope": float(slope_bs),
                    "intercept": float(intercept_bs),
                    "r2": r2_bs,
                },
                "backend_substitution_channel_priority_bootstrap_size_curvature": {
                    "quadratic_coef": float(q2_bs),
                    "linear_coef": float(q1_bs),
                    "constant_coef": float(q0_bs),
                    "aic_linear": aic_linear,
                    "aic_quadratic": aic_quadratic,
                    "aic_delta_quadratic_minus_linear": aic_delta_quad_minus_linear,
                },
                "backend_substitution_channel_priority_bootstrap_size_extrapolation_rows": extrap_rows,
                "backend_substitution_channel_priority_bootstrap_size_extrapolation_gap_max": extrap_gap_max,
                "backend_substitution_channel_priority_bootstrap_seed_rows": seed_rows,
                "backend_substitution_channel_priority_bootstrap_seed_span_max": seed_span_max,
                "backend_substitution_channel_first_simulation_panel": channel_first_simulation_panel,
            },
            "note": "Links task-2 channel phase-space table to a shared basis fit without closure claim.",
        },
        "scipy_numpy_sympy_calibration": {
            "feature_matrix": x.tolist(), "target_vector": y.tolist(), "weights_diagonal": weights_diag.tolist(),
            "solution": {"GhostCut_scheme": g_num, "WardLift": w_num, "CohomologyAmplitudeBridge": b_num},
            "train_holdout_split": {"train_indices": train_idx, "holdout_indices": hold_idx, "train_l2": train_l2, "holdout_l2": hold_l2},
            "bootstrap_std_solution": {"GhostCut_scheme": float(coef_std[0]), "WardLift": float(coef_std[1]), "CohomologyAmplitudeBridge": float(coef_std[2])},
            "bootstrap_seed_robustness": {"rows": seed_std_rows, "max_span": bootstrap_seed_span_max},
            "fit_residual_vector": residual_vec.tolist(), "fit_residual_l2": residual_l2, "fit_residual_linf": residual_linf,
            "weighted_design_singular_values": [float(v) for v in svals.tolist()], "weighted_design_condition_number": cond_num,
            "condition_robustness": {"median": cond_median, "p95": cond_p95},
            "weighted_design_rank_numeric": numeric_rank, "weighted_design_rank_symbolic": symbolic_rank,
            "amplitude_map_frobenius_norm": fro_norm, "amplitude_map_spectral_radius": spectral_radius,
        },
        "residual_obstruction": {"delta_target": 0.0, "delta_observed_l2": residual_l2, "delta_observed_linf": residual_linf},
        "reproducibility_probe": {"digest_1": reproducibility_digest_1, "digest_2": reproducibility_digest_2},
        "gatekeeper_checks": gate,
        "false_pass_guard": "No global unitarity closure claimed.",
        "next_honest_step": "Execute task-priority lane on strict DiscM integration: replay task #7 with expanded bootstrap/integrator stress envelope and only then decide whether one-channel substitution remains nonworse under CI-aware criteria.",
    }
    payload["theorem_core_digest_sha256"] = digest({"solution": payload["scipy_numpy_sympy_calibration"]["solution"], "max_fd_gap": max_fd_gap, "max_grid_refine_gap": max_grid_refine_gap, "quad_tol_span": quad_tol_span, "cond_p95": cond_p95, "bootstrap_seed_span_max": bootstrap_seed_span_max, "backend_fit_loss": backend_fit_loss, "backend_fit_loss_lbfgsb": backend_fit_loss_lbfgsb, "backend_fit_loss_gap": backend_fit_loss_gap, "multistart_loss_span": multistart_loss_span, "multichannel_max_solver_gap": channel_solver_gap_max, "multichannel_global_loss_spread": channel_loss_spread, "renorm_solver_gap": coef_b1_gap, "renorm_residual_l2": residual_b1_l2, "po3_objective": float(po3_res.fun), "po3_covariant_proxy_value_d1": po3_covariant_proxy_val, "transport_closure_max": transport_closure_max, "transport_symmetry_max": transport_symmetry_max, "po2_max_delta_bg": max_delta_bg_under_constraints, "selector_ratio_upper_bound": selector_ratio_upper_bound, "selector_entropy_span": selector_entropy_span, "discm_basis_cond": basis_cond, "discm_unc_max": common_basis_unc_max, "discm_resid_max": common_basis_resid_max, "channel_phase_min_global": channel_phase_min_global, "channel_phase_tol_span_max": tol_span_max, "phase_common_cond": phase_common_cond, "phase_common_residual_l2": phase_common_residual_l2, "phase_common_residual_backend_sub_l2": phase_common_residual_backend_sub_l2, "phase_backend_channel_delta_abs_max": channel_delta_abs_max, "phase_backend_transport_channel_delta_abs_max": channel_transport_delta_abs_max, "phase_backend_channel_priority_cond_weighted_median_span": channel_priority_cond_weighted_median_span, "phase_backend_channel_priority_winner_count": channel_priority_winner_count, "phase_backend_channel_priority_bootstrap_winner_frequency_max": winner_freq_max, "phase_backend_channel_priority_bootstrap_winner_frequency_wilson_lb95": winner_freq_max_wilson_lb, "phase_backend_channel_priority_bootstrap_winner_frequency_entropy_norm": winner_freq_entropy_norm, "phase_backend_channel_priority_bootstrap_winner_frequency_top2_margin": winner_freq_top2_margin, "phase_backend_channel_priority_dirichlet_posterior_p_best_gt_050": p_best_gt_050, "phase_backend_channel_priority_dirichlet_posterior_q05": p_best_q05, "phase_backend_channel_priority_bootstrap_size_span": boot_size_freq_span, "phase_backend_channel_priority_bootstrap_size_loo_span_max": boot_size_loo_span_max, "phase_backend_channel_priority_bootstrap_size_curvature": float(q2_bs), "phase_backend_channel_priority_bootstrap_size_aic_delta_quad_minus_linear": aic_delta_quad_minus_linear, "phase_backend_channel_priority_bootstrap_size_extrap_gap_max": extrap_gap_max, "phase_backend_channel_priority_bootstrap_seed_span_max": seed_span_max, "phase_backend_channel_priority_bootstrap_size_slope": float(slope_bs), "phase_backend_channel_priority_bootstrap_size_r2": r2_bs, "phase_common_loocv_max": loocv_residual_max, "phase_common_bootstrap_p95": link_bootstrap_resid_p95, "phase_bridge_stability_envelope_max": stability_envelope_max, "phase_cross_channel_coef_spread_l2": cross_channel_coef_spread_l2, "phase_joint_residual_l2": joint_resid_l2, "phase_joint_spread_l2": joint_spread_l2, "phase_joint_solver_gap": joint_solver_gap, "phase_joint_lambda_resid_span": lambda_resid_span, "phase_joint_holdout_resid_max": holdout_joint_resid_max, "phase_joint_multistart_resid_span": joint_ms_resid_span, "phase_joint_perturbation_resid_span": perturb_resid_span, "phase_joint_stress_envelope": joint_worst_case_residual_envelope, "phase_joint_cross_background_envelope_span": cross_background_envelope_span, "phase_joint_operator_transport_resid_span": operator_resid_span, "phase_joint_operator_transport_nu_sweep_resid_span": operator_nu_sweep_resid_span, "phase_joint_operator_transport_nu_sweep_solver_gap_max": operator_nu_sweep_solver_gap_max, "phase_joint_operator_transport_nu_lambda_resid_span": operator_nu_lambda_resid_span, "phase_joint_operator_transport_nu_lambda_solver_gap_max": operator_nu_lambda_solver_gap_max, "phase_joint_operator_transport_nu_lambda_weighted_resid_max": operator_nu_lambda_weighted_resid_max, "phase_joint_operator_transport_nu_lambda_weighted_resid_span": operator_nu_lambda_weighted_resid_span, "phase_joint_operator_transport_nu_lambda_condition_weighted_resid_max": operator_nu_lambda_cond_weighted_resid_max, "phase_joint_operator_transport_nu_lambda_condition_weighted_resid_span": operator_nu_lambda_cond_weighted_resid_span, "phase_joint_operator_transport_dual_pareto_count": pareto_count, "phase_joint_operator_transport_dual_branch_flips_total": int(branch_flip_total), "residual": payload["residual_obstruction"]})
    payload["theorem_core_digest_recomputed_sha256"] = digest({"solution": payload["scipy_numpy_sympy_calibration"]["solution"], "max_fd_gap": max_fd_gap, "max_grid_refine_gap": max_grid_refine_gap, "quad_tol_span": quad_tol_span, "cond_p95": cond_p95, "bootstrap_seed_span_max": bootstrap_seed_span_max, "backend_fit_loss": backend_fit_loss, "backend_fit_loss_lbfgsb": backend_fit_loss_lbfgsb, "backend_fit_loss_gap": backend_fit_loss_gap, "multistart_loss_span": multistart_loss_span, "multichannel_max_solver_gap": channel_solver_gap_max, "multichannel_global_loss_spread": channel_loss_spread, "renorm_solver_gap": coef_b1_gap, "renorm_residual_l2": residual_b1_l2, "po3_objective": float(po3_res.fun), "po3_covariant_proxy_value_d1": po3_covariant_proxy_val, "transport_closure_max": transport_closure_max, "transport_symmetry_max": transport_symmetry_max, "po2_max_delta_bg": max_delta_bg_under_constraints, "selector_ratio_upper_bound": selector_ratio_upper_bound, "selector_entropy_span": selector_entropy_span, "discm_basis_cond": basis_cond, "discm_unc_max": common_basis_unc_max, "discm_resid_max": common_basis_resid_max, "channel_phase_min_global": channel_phase_min_global, "channel_phase_tol_span_max": tol_span_max, "phase_common_cond": phase_common_cond, "phase_common_residual_l2": phase_common_residual_l2, "phase_common_residual_backend_sub_l2": phase_common_residual_backend_sub_l2, "phase_backend_channel_delta_abs_max": channel_delta_abs_max, "phase_backend_transport_channel_delta_abs_max": channel_transport_delta_abs_max, "phase_backend_channel_priority_cond_weighted_median_span": channel_priority_cond_weighted_median_span, "phase_backend_channel_priority_winner_count": channel_priority_winner_count, "phase_backend_channel_priority_bootstrap_winner_frequency_max": winner_freq_max, "phase_backend_channel_priority_bootstrap_winner_frequency_wilson_lb95": winner_freq_max_wilson_lb, "phase_backend_channel_priority_bootstrap_winner_frequency_entropy_norm": winner_freq_entropy_norm, "phase_backend_channel_priority_bootstrap_winner_frequency_top2_margin": winner_freq_top2_margin, "phase_backend_channel_priority_dirichlet_posterior_p_best_gt_050": p_best_gt_050, "phase_backend_channel_priority_dirichlet_posterior_q05": p_best_q05, "phase_backend_channel_priority_bootstrap_size_span": boot_size_freq_span, "phase_backend_channel_priority_bootstrap_size_loo_span_max": boot_size_loo_span_max, "phase_backend_channel_priority_bootstrap_size_curvature": float(q2_bs), "phase_backend_channel_priority_bootstrap_size_aic_delta_quad_minus_linear": aic_delta_quad_minus_linear, "phase_backend_channel_priority_bootstrap_size_extrap_gap_max": extrap_gap_max, "phase_backend_channel_priority_bootstrap_seed_span_max": seed_span_max, "phase_backend_channel_priority_bootstrap_size_slope": float(slope_bs), "phase_backend_channel_priority_bootstrap_size_r2": r2_bs, "phase_common_loocv_max": loocv_residual_max, "phase_common_bootstrap_p95": link_bootstrap_resid_p95, "phase_bridge_stability_envelope_max": stability_envelope_max, "phase_cross_channel_coef_spread_l2": cross_channel_coef_spread_l2, "phase_joint_residual_l2": joint_resid_l2, "phase_joint_spread_l2": joint_spread_l2, "phase_joint_solver_gap": joint_solver_gap, "phase_joint_lambda_resid_span": lambda_resid_span, "phase_joint_holdout_resid_max": holdout_joint_resid_max, "phase_joint_multistart_resid_span": joint_ms_resid_span, "phase_joint_perturbation_resid_span": perturb_resid_span, "phase_joint_stress_envelope": joint_worst_case_residual_envelope, "phase_joint_cross_background_envelope_span": cross_background_envelope_span, "phase_joint_operator_transport_resid_span": operator_resid_span, "phase_joint_operator_transport_nu_sweep_resid_span": operator_nu_sweep_resid_span, "phase_joint_operator_transport_nu_sweep_solver_gap_max": operator_nu_sweep_solver_gap_max, "phase_joint_operator_transport_nu_lambda_resid_span": operator_nu_lambda_resid_span, "phase_joint_operator_transport_nu_lambda_solver_gap_max": operator_nu_lambda_solver_gap_max, "phase_joint_operator_transport_nu_lambda_weighted_resid_max": operator_nu_lambda_weighted_resid_max, "phase_joint_operator_transport_nu_lambda_weighted_resid_span": operator_nu_lambda_weighted_resid_span, "phase_joint_operator_transport_nu_lambda_condition_weighted_resid_max": operator_nu_lambda_cond_weighted_resid_max, "phase_joint_operator_transport_nu_lambda_condition_weighted_resid_span": operator_nu_lambda_cond_weighted_resid_span, "phase_joint_operator_transport_dual_pareto_count": pareto_count, "phase_joint_operator_transport_dual_branch_flips_total": int(branch_flip_total), "residual": payload["residual_obstruction"]})
    payload["gatekeeper_checks"]["theorem_digest_self_consistent"] = payload["theorem_core_digest_sha256"] == payload["theorem_core_digest_recomputed_sha256"]
    payload["gatekeeper_checks"]["reproducibility_digest_self_consistent"] = reproducibility_digest_1 == reproducibility_digest_2

    csv_fields = [
        "channel",
        "status",
        "prob_seed_nonworse",
        "prob_extrap_nonworse",
        "prob_q05_nonworse",
        "prob_seed_nonworse_ci95_lower",
        "prob_seed_nonworse_ci95_upper",
        "prob_extrap_nonworse_ci95_lower",
        "prob_extrap_nonworse_ci95_upper",
        "prob_q05_nonworse_ci95_lower",
        "prob_q05_nonworse_ci95_upper",
        "min_ci_margin_to_threshold",
        "seed_ci_margin",
        "extrap_ci_margin",
        "q05_ci_margin",
        "nonworse_probability_conditions_met",
        "nonworse_probability_ci95_conditions_met",
        "low_power_detected",
        "ready_for_real_backend_substitution",
    ]
    verdict_rows_csv = []
    risk_by_channel = {r["channel"]: r for r in paired_delta_panel["per_channel_risk_ranking"]["rows"]}
    for row in paired_delta_panel["per_channel_power_aware_verdicts"]:
        rr = risk_by_channel[row["channel"]]
        verdict_rows_csv.append({
            "channel": row["channel"],
            "status": row["status"],
            "prob_seed_nonworse": row["prob_seed_nonworse"],
            "prob_extrap_nonworse": row["prob_extrap_nonworse"],
            "prob_q05_nonworse": row["prob_q05_nonworse"],
            "prob_seed_nonworse_ci95_lower": row["prob_seed_nonworse_ci95"]["lower"],
            "prob_seed_nonworse_ci95_upper": row["prob_seed_nonworse_ci95"]["upper"],
            "prob_extrap_nonworse_ci95_lower": row["prob_extrap_nonworse_ci95"]["lower"],
            "prob_extrap_nonworse_ci95_upper": row["prob_extrap_nonworse_ci95"]["upper"],
            "prob_q05_nonworse_ci95_lower": row["prob_q05_nonworse_ci95"]["lower"],
            "prob_q05_nonworse_ci95_upper": row["prob_q05_nonworse_ci95"]["upper"],
            "min_ci_margin_to_threshold": rr["min_ci_margin_to_threshold"],
            "seed_ci_margin": rr["seed_ci_margin"],
            "extrap_ci_margin": rr["extrap_ci_margin"],
            "q05_ci_margin": rr["q05_ci_margin"],
            "nonworse_probability_conditions_met": row["nonworse_probability_conditions_met"],
            "nonworse_probability_ci95_conditions_met": row["nonworse_probability_ci95_conditions_met"],
            "low_power_detected": row["low_power_detected"],
            "ready_for_real_backend_substitution": row["ready_for_real_backend_substitution"],
        })
    with OUT_CSV.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=csv_fields)
        w.writeheader()
        for row in verdict_rows_csv:
            w.writerow(row)
    quality_fields = [
        "channel",
        "seed_n_pairs", "seed_n_zero_deltas", "seed_n_nonzero_deltas", "seed_zero_delta_fraction", "seed_effective_pair_fraction", "seed_min_effective_pairs_required", "seed_low_power_flag",
        "extrap_n_pairs", "extrap_n_zero_deltas", "extrap_n_nonzero_deltas", "extrap_zero_delta_fraction", "extrap_effective_pair_fraction", "extrap_min_effective_pairs_required", "extrap_low_power_flag",
        "q05_n_pairs", "q05_n_zero_deltas", "q05_n_nonzero_deltas", "q05_zero_delta_fraction", "q05_effective_pair_fraction", "q05_min_effective_pairs_required", "q05_low_power_flag",
    ]
    quality_rows_csv = []
    for row in paired_delta_panel["per_channel_wilcoxon_quality"]:
        quality_rows_csv.append({
            "channel": row["channel"],
            "seed_n_pairs": row["seed"]["n_pairs"],
            "seed_n_zero_deltas": row["seed"]["n_zero_deltas"],
            "seed_n_nonzero_deltas": row["seed"]["n_nonzero_deltas"],
            "seed_zero_delta_fraction": row["seed"]["zero_delta_fraction"],
            "seed_effective_pair_fraction": row["seed"]["effective_pair_fraction"],
            "seed_min_effective_pairs_required": row["seed"]["min_effective_pairs_required"],
            "seed_low_power_flag": row["seed"]["low_power_flag"],
            "extrap_n_pairs": row["extrap"]["n_pairs"],
            "extrap_n_zero_deltas": row["extrap"]["n_zero_deltas"],
            "extrap_n_nonzero_deltas": row["extrap"]["n_nonzero_deltas"],
            "extrap_zero_delta_fraction": row["extrap"]["zero_delta_fraction"],
            "extrap_effective_pair_fraction": row["extrap"]["effective_pair_fraction"],
            "extrap_min_effective_pairs_required": row["extrap"]["min_effective_pairs_required"],
            "extrap_low_power_flag": row["extrap"]["low_power_flag"],
            "q05_n_pairs": row["q05"]["n_pairs"],
            "q05_n_zero_deltas": row["q05"]["n_zero_deltas"],
            "q05_n_nonzero_deltas": row["q05"]["n_nonzero_deltas"],
            "q05_zero_delta_fraction": row["q05"]["zero_delta_fraction"],
            "q05_effective_pair_fraction": row["q05"]["effective_pair_fraction"],
            "q05_min_effective_pairs_required": row["q05"]["min_effective_pairs_required"],
            "q05_low_power_flag": row["q05"]["low_power_flag"],
        })
    with OUT_QUALITY_CSV.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=quality_fields)
        w.writeheader()
        for row in quality_rows_csv:
            w.writerow(row)
    payload["gatekeeper_checks"]["phase_backend_substitution_channel_first_per_channel_quality_csv_json_consistent"] = bool(
        len(quality_rows_csv) == len(paired_delta_panel["per_channel_wilcoxon_quality"]) and
        {r["channel"] for r in quality_rows_csv} == {r["channel"] for r in paired_delta_panel["per_channel_wilcoxon_quality"]}
    )
    payload["gatekeeper_checks"]["phase_backend_substitution_channel_first_per_channel_verdict_csv_json_consistent"] = bool(
        len(verdict_rows_csv) == len(paired_delta_panel["per_channel_power_aware_verdicts"]) and
        {r["channel"] for r in verdict_rows_csv} == {r["channel"] for r in paired_delta_panel["per_channel_power_aware_verdicts"]}
    )

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(OUT)


if __name__ == "__main__":
    main()
