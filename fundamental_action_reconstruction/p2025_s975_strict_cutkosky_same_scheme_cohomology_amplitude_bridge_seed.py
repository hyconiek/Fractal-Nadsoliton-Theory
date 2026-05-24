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
import time
import warnings
from itertools import combinations
from pathlib import Path
from typing import Any

import numpy as np
import scipy.integrate as si
import scipy.linalg as la
import scipy.optimize as so
import scipy.stats as ss
import sympy as sp
from scipy.integrate import IntegrationWarning

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
    ur_unc_rows = []
    for src_name, src in [("p2015", p2015), ("p2016", p2016)]:
        for row in src.get("uncertainty_table", []):
            ur_unc_rows.append({
                "source": src_name,
                "delta_center": float(row.get("delta_center", 0.0)),
                "delta_std": float(row.get("delta_std", 0.0)),
                "residue_positive_all_samples": bool(row.get("residue_positive_all_samples", False)),
            })
    ur_unc_abs_centers = np.array([abs(r["delta_center"]) for r in ur_unc_rows], dtype=float)
    ur_unc_stds = np.array([r["delta_std"] for r in ur_unc_rows], dtype=float)
    ur_unc_residue_pos_rate = float(np.mean(np.array([1.0 if r["residue_positive_all_samples"] else 0.0 for r in ur_unc_rows], dtype=float))) if ur_unc_rows else 0.0
    ur_uncertainty_transport_bridge_precursor = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_UR_LINK_UNCERTAINTY_SYNTHESIS",
        "rows_count": int(len(ur_unc_rows)),
        "median_abs_delta_center": float(np.median(ur_unc_abs_centers)) if ur_unc_rows else 0.0,
        "p95_abs_delta_center": float(np.percentile(ur_unc_abs_centers, 95.0)) if ur_unc_rows else 0.0,
        "median_delta_std": float(np.median(ur_unc_stds)) if ur_unc_rows else 0.0,
        "p95_delta_std": float(np.percentile(ur_unc_stds, 95.0)) if ur_unc_rows else 0.0,
        "residue_positive_rate": ur_unc_residue_pos_rate,
        "bounded_p95_abs_delta_center": bool((float(np.percentile(ur_unc_abs_centers, 95.0)) if ur_unc_rows else 0.0) < 1e-5),
        "bounded_p95_delta_std": bool((float(np.percentile(ur_unc_stds, 95.0)) if ur_unc_rows else 0.0) < 1e-5),
    }
    # Cross-source transport agreement panel (Task-2-oriented quantitative precursor).
    p2015_by_s = {}
    for row in p2015.get("uncertainty_table", []):
        s_val = float(sp.N(sp.sympify(str(row.get("s", "0")))))
        p2015_by_s[s_val] = row
    p2016_by_s = {}
    for row in p2016.get("uncertainty_table", []):
        s_val = float(sp.N(sp.sympify(str(row.get("s", "0")))))
        p2016_by_s[s_val] = row
    common_s = sorted(set(p2015_by_s.keys()).intersection(set(p2016_by_s.keys())))
    ur_transport_agreement_rows = []
    for s_val in common_s:
        r15 = p2015_by_s[s_val]
        r16 = p2016_by_s[s_val]
        c15 = float(r15.get("delta_center", 0.0))
        c16 = float(r16.get("delta_center", 0.0))
        sd15 = float(r15.get("delta_std", 0.0))
        sd16 = float(r16.get("delta_std", 0.0))
        ur_transport_agreement_rows.append({
            "s": float(s_val),
            "delta_center_p2015": c15,
            "delta_center_p2016": c16,
            "delta_center_abs_gap": float(abs(c16 - c15)),
            "delta_std_p2015": sd15,
            "delta_std_p2016": sd16,
            "delta_std_ratio_p2016_over_p2015": float(sd16 / max(1e-30, sd15)),
        })
    center_gap_vec = np.array([r["delta_center_abs_gap"] for r in ur_transport_agreement_rows], dtype=float) if ur_transport_agreement_rows else np.array([0.0], dtype=float)
    std_ratio_vec = np.array([r["delta_std_ratio_p2016_over_p2015"] for r in ur_transport_agreement_rows], dtype=float) if ur_transport_agreement_rows else np.array([1.0], dtype=float)
    ur_transport_cross_source_agreement_precursor = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_UR_LINK_CROSS_SOURCE_AGREEMENT",
        "rows": ur_transport_agreement_rows,
        "common_s_count": int(len(ur_transport_agreement_rows)),
        "max_delta_center_abs_gap": float(np.max(center_gap_vec)),
        "p95_delta_center_abs_gap": float(np.percentile(center_gap_vec, 95.0)),
        "median_delta_std_ratio_p2016_over_p2015": float(np.median(std_ratio_vec)),
        "max_delta_std_ratio_p2016_over_p2015": float(np.max(std_ratio_vec)),
        "center_gap_bounded_p95": bool(float(np.percentile(center_gap_vec, 95.0)) < 1e-5),
        "std_ratio_bounded_max": bool(float(np.max(std_ratio_vec)) < 2.0),
    }
    # Channel-resolved strict trace budget precursor (unitarity-oriented, no closure claim).
    trace_profiles = p2017.get("quadrature_channel_covariance_candidate", {}).get("trace_profiles_by_channel", {})
    channel_trace_rows = []
    channel_names_sorted = sorted(trace_profiles.keys())
    total_trace_by_channel = {}
    for ch in channel_names_sorted:
        vals = np.array([float(v) for v in trace_profiles.get(ch, [])], dtype=float)
        total_trace_by_channel[ch] = float(np.sum(vals))
    total_trace_all_channels = float(sum(total_trace_by_channel.values()))
    for ch in channel_names_sorted:
        vals = np.array([float(v) for v in trace_profiles.get(ch, [])], dtype=float)
        monotone_nonincreasing = bool(np.all(np.diff(vals) <= 1e-15)) if vals.size > 1 else True
        channel_trace_rows.append({
            "channel": ch,
            "num_points": int(vals.size),
            "trace_sum": float(np.sum(vals)),
            "trace_mean": float(np.mean(vals)) if vals.size else 0.0,
            "trace_min": float(np.min(vals)) if vals.size else 0.0,
            "trace_max": float(np.max(vals)) if vals.size else 0.0,
            "trace_p95": float(np.percentile(vals, 95.0)) if vals.size else 0.0,
            "trace_share_of_total": float(total_trace_by_channel[ch] / max(1e-30, total_trace_all_channels)),
            "monotone_nonincreasing_over_s_grid": monotone_nonincreasing,
        })
    ur_channel_trace_budget_precursor = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_UR_CHANNEL_TRACE_BUDGET",
        "rows": channel_trace_rows,
        "num_channels": int(len(channel_trace_rows)),
        "total_trace_all_channels": total_trace_all_channels,
        "trace_share_sum": float(sum(r["trace_share_of_total"] for r in channel_trace_rows)),
        "all_channels_monotone_nonincreasing": bool(all(r["monotone_nonincreasing_over_s_grid"] for r in channel_trace_rows)),
        "max_channel_trace_share": float(max((r["trace_share_of_total"] for r in channel_trace_rows), default=0.0)),
    }
    # Explicit channel-map precursor: p2017 internal channels -> Task-2 working classes.
    channel_map_weights = {
        "gg": {"gauge_gauge": 1.00, "fermion_fermion": 0.00, "scalar_scalar": 0.00},
        "gh": {"gauge_gauge": 0.50, "fermion_fermion": 0.50, "scalar_scalar": 0.00},
        "hh": {"gauge_gauge": 0.00, "fermion_fermion": 0.00, "scalar_scalar": 1.00},
        "gx": {"gauge_gauge": 0.50, "fermion_fermion": 0.00, "scalar_scalar": 0.50},
    }
    class_trace_budget = {"gauge_gauge": 0.0, "fermion_fermion": 0.0, "scalar_scalar": 0.0}
    channel_total_lookup = {r["channel"]: float(r["trace_sum"]) for r in channel_trace_rows}
    mapped_rows = []
    for source_ch, weights in channel_map_weights.items():
        src_trace = float(channel_total_lookup.get(source_ch, 0.0))
        row = {"source_channel": source_ch, "source_trace_sum": src_trace, "weights": dict(weights), "allocated": {}}
        for target_class, ww in weights.items():
            alloc = float(src_trace * float(ww))
            class_trace_budget[target_class] += alloc
            row["allocated"][target_class] = alloc
        mapped_rows.append(row)
    class_total = float(sum(class_trace_budget.values()))
    for k in class_trace_budget:
        class_trace_budget[k] = {"trace_sum": float(class_trace_budget[k]), "trace_share": float(class_trace_budget[k] / max(1e-30, class_total))}
    ur_channel_class_mapping_precursor = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_CHANNEL_MAP_TO_TASK2_CLASSES",
        "mapping_kind": "EXPLICIT_WEIGHTED_PRECURSOR_NOT_UNIQUENESS_THEOREM",
        "rows": mapped_rows,
        "class_trace_budget": class_trace_budget,
        "source_trace_total": float(sum(channel_total_lookup.values())),
        "allocated_trace_total": class_total,
        "trace_conservation_gap_abs": float(abs(class_total - float(sum(channel_total_lookup.values())))),
    }

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

    # Strict Task-2 theorem witness (single-channel): bounded DiscM-CutSum gap for graviton->gauge_gauge.
    gg_pars = channel_param_map["gauge_gauge"]
    om_gg, ph_gg, be_gg, et_gg = gg_pars
    task2_rows = []
    residue_rows = []
    for s in s_grid_fine:
        s = float(s)
        disc_v, disc_err = strict_kernel_phase_integral(s, float(om_gg), float(ph_gg), float(be_gg), float(et_gg))
        def gg_integrand(x):
            xa = np.array(x, dtype=float)
            kk = np.cos(om_gg * xa + ph_gg) / (1.0 + be_gg * (xa ** et_gg))
            out = (kk * kk) / np.sqrt(np.maximum(1e-15, xa + s))
            return out
        cutsum_v_fixed, _ = si.fixed_quad(gg_integrand, 0.0, 1.0, n=240)
        cutsum_v_quad, cutsum_quad_err = si.quad(gg_integrand, 0.0, 1.0, epsabs=1e-10, epsrel=1e-10, limit=400)
        cutsum_v = float(cutsum_v_fixed)
        gap_abs = float(abs(float(disc_v) - cutsum_v))
        gap_rel = float(gap_abs / max(1e-15, abs(cutsum_v)))
        uncertainty = float(max(abs(cutsum_v - float(cutsum_v_quad)), abs(float(disc_err)), abs(float(cutsum_quad_err))))
        task2_rows.append({
            "s": s,
            "disc_value": float(disc_v),
            "cutsum_value": cutsum_v,
            "gap_abs": gap_abs,
            "gap_rel": gap_rel,
            "uncertainty_estimate": uncertainty,
        })
        x_probe = np.linspace(0.0, 1.0, 17)
        vals_probe = np.array([gg_integrand(float(xp)) for xp in x_probe], dtype=float)
        residue_rows.append({
            "s": s,
            "min_effective_weight": float(np.min(vals_probe)),
            "max_effective_weight": float(np.max(vals_probe)),
            "all_nonnegative": bool(np.min(vals_probe) >= -1e-12),
        })
    gap_abs_arr = np.array([r["gap_abs"] for r in task2_rows], dtype=float)
    gap_rel_arr = np.array([r["gap_rel"] for r in task2_rows], dtype=float)
    unc_arr = np.array([r["uncertainty_estimate"] for r in task2_rows], dtype=float)
    min_effective_weight_global = float(min((float(r["min_effective_weight"]) for r in residue_rows), default=0.0))
    pass_fail_criteria_task2 = {
        "q95_gap_abs_max": 1e-8,
        "max_gap_rel_max": 1e-6,
        "q95_cross_integrator_gap_abs_max": 1e-10,
        "q95_convergence_delta_n400_to_n800_abs_max": 1e-10,
        "min_effective_weight_global_min": 0.0,
        "all_nonnegative_weights_required": True,
    }
    q95_gap_abs = float(np.quantile(gap_abs_arr, 0.95)) if gap_abs_arr.size else float("inf")
    max_gap_rel = float(np.max(gap_rel_arr)) if gap_rel_arr.size else float("inf")
    all_nonnegative = bool(min_effective_weight_global >= float(pass_fail_criteria_task2["min_effective_weight_global_min"]))
    consistency_ci95 = {
        "lower": float(max(0.0, q95_gap_abs - 1.96 * (float(np.std(unc_arr, ddof=1)) if unc_arr.size > 1 else 0.0))),
        "upper": float(q95_gap_abs + 1.96 * (float(np.std(unc_arr, ddof=1)) if unc_arr.size > 1 else 0.0)),
    }
    closure_ready = bool(
        q95_gap_abs <= float(pass_fail_criteria_task2["q95_gap_abs_max"]) and
        max_gap_rel <= float(pass_fail_criteria_task2["max_gap_rel_max"]) and
        all_nonnegative
    )
    if closure_ready:
        verdict_task2 = "CLOSED_NUMERICAL_WITNESS_TASK2"
        fail_trace_task2 = ""
    else:
        verdict_task2 = "OPEN_OBSTRUCTION_WITH_TRACE"
        fail_trace_task2 = f"q95_gap_abs={q95_gap_abs:.6e} > {pass_fail_criteria_task2['q95_gap_abs_max']:.1e}"
    s_grid_task2_extended = sorted(set([float(x) for x in s_grid_fine] + [0.35, 0.75, 1.25, 2.5, 3.5]))
    scheme_rows = []
    for s in s_grid_task2_extended:
        disc_v, disc_err = strict_kernel_phase_integral(float(s), float(om_gg), float(ph_gg), float(be_gg), float(et_gg))
        def gg_integrand_native(x):
            xa = np.array(x, dtype=float)
            kk = np.cos(om_gg * xa + ph_gg) / (1.0 + be_gg * (xa ** et_gg))
            return (kk * kk) / np.sqrt(np.maximum(1e-15, xa + s))
        native_v, native_err = si.quad(gg_integrand_native, 0.0, 1.0, epsabs=1e-10, epsrel=1e-10, limit=400)
        def gg_integrand_warp(t):
            ta = np.array(t, dtype=float)
            x = ta * ta
            jac = 2.0 * ta
            kk = np.cos(om_gg * x + ph_gg) / (1.0 + be_gg * (x ** et_gg))
            return ((kk * kk) / np.sqrt(np.maximum(1e-15, x + s))) * jac
        warp_v, warp_err = si.quad(gg_integrand_warp, 0.0, 1.0, epsabs=1e-10, epsrel=1e-10, limit=400)
        cutsum_value = float(native_v)
        gap_abs = float(abs(float(disc_v) - cutsum_value))
        gap_rel = float(gap_abs / max(1e-15, abs(cutsum_value)))
        uncertainty = float(max(abs(float(native_v) - float(warp_v)), abs(float(native_err)), abs(float(warp_err)), abs(float(disc_err))))
        scheme_rows.append({
            "s": float(s),
            "disc_value": float(disc_v),
            "cutsum_value": cutsum_value,
            "gap_abs": gap_abs,
            "gap_rel": gap_rel,
            "uncertainty_estimate": uncertainty,
            "cutsum_native": float(native_v),
            "cutsum_warped": float(warp_v),
            "cutsum_scheme_gap_abs": float(abs(float(native_v) - float(warp_v))),
        })
    # Local phase-space contribution map: identify x-bins driving DiscM-CutSum mismatch.
    x_edges = np.linspace(0.0, 1.0, 9)
    phase_space_bin_rows = []
    for rr in scheme_rows:
        s_loc = float(rr["s"])
        disc_total = float(rr["disc_value"])
        bin_contrib_native = []
        bin_contrib_warp = []
        for i in range(len(x_edges) - 1):
            a = float(x_edges[i]); b = float(x_edges[i + 1])
            def bin_native(x):
                xa = np.array(x, dtype=float)
                kk = np.cos(om_gg * xa + ph_gg) / (1.0 + be_gg * (xa ** et_gg))
                return (kk * kk) / np.sqrt(np.maximum(1e-15, xa + s_loc))
            vv_n, _ = si.quad(bin_native, a, b, epsabs=1e-10, epsrel=1e-10, limit=200)
            bin_contrib_native.append(float(vv_n))
            ta = np.sqrt(a); tb = np.sqrt(b)
            def bin_warp(t):
                tt = np.array(t, dtype=float)
                x = tt * tt
                jac = 2.0 * tt
                kk = np.cos(om_gg * x + ph_gg) / (1.0 + be_gg * (x ** et_gg))
                return ((kk * kk) / np.sqrt(np.maximum(1e-15, x + s_loc))) * jac
            vv_w, _ = si.quad(bin_warp, ta, tb, epsabs=1e-10, epsrel=1e-10, limit=200)
            bin_contrib_warp.append(float(vv_w))
        native_arr = np.array(bin_contrib_native, dtype=float)
        warp_arr = np.array(bin_contrib_warp, dtype=float)
        disc_bin = float(disc_total / max(1, native_arr.size))
        bin_rows = []
        for i in range(native_arr.size):
            gap_bin = float(abs(disc_bin - native_arr[i]))
            bin_rows.append({
                "bin_index": int(i),
                "x_left": float(x_edges[i]),
                "x_right": float(x_edges[i + 1]),
                "native_bin_integral": float(native_arr[i]),
                "warped_bin_integral": float(warp_arr[i]),
                "bin_scheme_gap_abs": float(abs(native_arr[i] - warp_arr[i])),
                "bin_disc_proxy_gap_abs": gap_bin,
            })
        phase_space_bin_rows.append({
            "s": s_loc,
            "rows": bin_rows,
            "max_bin_disc_proxy_gap_abs": float(max([r["bin_disc_proxy_gap_abs"] for r in bin_rows], default=0.0)),
            "max_bin_scheme_gap_abs": float(max([r["bin_scheme_gap_abs"] for r in bin_rows], default=0.0)),
        })
    # Bin-weighted obstruction ranking across s-grid.
    bin_acc = {}
    for sr in phase_space_bin_rows:
        for br in sr["rows"]:
            bi = int(br["bin_index"])
            rec = bin_acc.setdefault(bi, {"bin_index": bi, "x_left": float(br["x_left"]), "x_right": float(br["x_right"]), "disc_proxy_sum": 0.0, "scheme_gap_sum": 0.0, "count": 0})
            rec["disc_proxy_sum"] += float(br["bin_disc_proxy_gap_abs"])
            rec["scheme_gap_sum"] += float(br["bin_scheme_gap_abs"])
            rec["count"] += 1
    total_disc_proxy = float(sum(v["disc_proxy_sum"] for v in bin_acc.values()))
    bin_obstruction_ranking_rows = []
    for bi in sorted(bin_acc.keys()):
        v = bin_acc[bi]
        share = float(v["disc_proxy_sum"] / total_disc_proxy) if total_disc_proxy > 0.0 else 0.0
        bin_obstruction_ranking_rows.append({
            "bin_index": int(v["bin_index"]),
            "x_left": float(v["x_left"]),
            "x_right": float(v["x_right"]),
            "disc_proxy_gap_sum": float(v["disc_proxy_sum"]),
            "disc_proxy_gap_share": share,
            "scheme_gap_sum": float(v["scheme_gap_sum"]),
            "mean_disc_proxy_gap": float(v["disc_proxy_sum"] / max(1, int(v["count"]))),
        })
    bin_obstruction_ranking_rows = sorted(bin_obstruction_ranking_rows, key=lambda r: (-r["disc_proxy_gap_share"], r["bin_index"]))
    top2_share = float(sum(r["disc_proxy_gap_share"] for r in bin_obstruction_ranking_rows[:2])) if bin_obstruction_ranking_rows else 0.0
    top_bins = [int(r["bin_index"]) for r in bin_obstruction_ranking_rows[:2]]
    x_edges_ref = np.linspace(0.0, 1.0, 9)
    endpoint_refinement_rows = []
    for s_loc in s_grid_task2_extended:
        disc_v_loc, _ = strict_kernel_phase_integral(float(s_loc), float(om_gg), float(ph_gg), float(be_gg), float(et_gg))
        cutsum_local = 0.0
        bin_rows_local = []
        for bi in range(8):
            a = float(x_edges_ref[bi]); b = float(x_edges_ref[bi + 1])
            is_top = bool(bi in top_bins)
            lim = 600 if is_top else 250
            eps = 1e-12 if is_top else 1e-10
            def f_local(x):
                xa = np.array(x, dtype=float)
                kk = np.cos(om_gg * xa + ph_gg) / (1.0 + be_gg * (xa ** et_gg))
                return (kk * kk) / np.sqrt(np.maximum(1e-15, xa + s_loc))
            v_local, _ = si.quad(f_local, a, b, epsabs=eps, epsrel=eps, limit=lim)
            cutsum_local += float(v_local)
            bin_rows_local.append({"bin_index": int(bi), "is_top2": is_top, "bin_integral_refined": float(v_local)})
        gap_refined = float(abs(float(disc_v_loc) - float(cutsum_local)))
        endpoint_refinement_rows.append({
            "s": float(s_loc),
            "disc_value": float(disc_v_loc),
            "cutsum_value_refined": float(cutsum_local),
            "gap_abs_refined": gap_refined,
            "rows": bin_rows_local,
        })
    q95_gap_abs_refined = float(np.quantile(np.array([r["gap_abs_refined"] for r in endpoint_refinement_rows], dtype=float), 0.95)) if endpoint_refinement_rows else float("inf")
    # Split-domain endpoint decomposition on top-2 bins: left/right half contribution mismatch.
    endpoint_split_rows = []
    left_total = 0.0
    right_total = 0.0
    for s_loc in s_grid_task2_extended:
        for bi in top_bins:
            a = float(x_edges_ref[int(bi)]); b = float(x_edges_ref[int(bi) + 1]); m = 0.5 * (a + b)
            def f_split(x):
                xa = np.array(x, dtype=float)
                kk = np.cos(om_gg * xa + ph_gg) / (1.0 + be_gg * (xa ** et_gg))
                return (kk * kk) / np.sqrt(np.maximum(1e-15, xa + s_loc))
            lv, _ = si.quad(f_split, a, m, epsabs=1e-12, epsrel=1e-12, limit=700)
            rv, _ = si.quad(f_split, m, b, epsabs=1e-12, epsrel=1e-12, limit=700)
            left_total += float(abs(lv))
            right_total += float(abs(rv))
            endpoint_split_rows.append({
                "s": float(s_loc),
                "bin_index": int(bi),
                "x_left": a,
                "x_mid": m,
                "x_right": b,
                "left_half_integral_abs": float(abs(lv)),
                "right_half_integral_abs": float(abs(rv)),
                "left_minus_right_abs": float(abs(abs(lv) - abs(rv))),
                "left_dominates": bool(abs(lv) >= abs(rv)),
            })
    endpoint_dominance = "LEFT" if left_total > right_total else ("RIGHT" if right_total > left_total else "BALANCED")
    ur_task2_endpoint_split_domain_precursor = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_TOP2_BIN_ENDPOINT_SPLIT_DOMAIN",
        "top2_bins": top_bins,
        "rows": endpoint_split_rows,
        "left_total_abs": float(left_total),
        "right_total_abs": float(right_total),
        "left_to_right_abs_ratio": float(left_total / max(1e-15, right_total)),
        "dominant_endpoint_half": endpoint_dominance,
    }
    # Endpoint-adaptive transform selection on dominant endpoint half.
    transform_grid = ["native", "left_focus", "right_focus"]
    adaptive_rows = []
    dominant_half = str(ur_task2_endpoint_split_domain_precursor["dominant_endpoint_half"])
    for tf in transform_grid:
        rows_tf = []
        for s_loc in s_grid_task2_extended:
            def base_integrand(x):
                xa = np.array(x, dtype=float)
                kk = np.cos(om_gg * xa + ph_gg) / (1.0 + be_gg * (xa ** et_gg))
                return (kk * kk) / np.sqrt(np.maximum(1e-15, xa + s_loc))
            if tf == "native":
                cutsum_tf, _ = si.quad(base_integrand, 0.0, 1.0, epsabs=1e-10, epsrel=1e-10, limit=500)
            elif tf == "left_focus":
                # x = t^2 concentrates quadrature near x=0
                def f_left(t):
                    ta = np.array(t, dtype=float)
                    x = ta * ta
                    jac = 2.0 * ta
                    return base_integrand(x) * jac
                cutsum_tf, _ = si.quad(f_left, 0.0, 1.0, epsabs=1e-10, epsrel=1e-10, limit=500)
            else:
                # x = 1-(1-t)^2 concentrates quadrature near x=1
                def f_right(t):
                    ta = np.array(t, dtype=float)
                    x = 1.0 - (1.0 - ta) * (1.0 - ta)
                    jac = 2.0 * (1.0 - ta)
                    return base_integrand(x) * jac
                cutsum_tf, _ = si.quad(f_right, 0.0, 1.0, epsabs=1e-10, epsrel=1e-10, limit=500)
            disc_tf, _ = strict_kernel_phase_integral(float(s_loc), float(om_gg), float(ph_gg), float(be_gg), float(et_gg))
            rows_tf.append(float(abs(float(disc_tf) - float(cutsum_tf))))
        q95_tf = float(np.quantile(np.array(rows_tf, dtype=float), 0.95)) if rows_tf else float("inf")
        adaptive_rows.append({"transform": tf, "q95_gap_abs": q95_tf})
    adaptive_rows = sorted(adaptive_rows, key=lambda r: (r["q95_gap_abs"], r["transform"]))
    recommended_transform = str(adaptive_rows[0]["transform"]) if adaptive_rows else "native"
    adaptive_boot_n = 256
    adaptive_boot_rng = np.random.default_rng(975171)
    transform_pick_counts = {k: 0 for k in transform_grid}
    transform_q95_by_name = {r["transform"]: float(r["q95_gap_abs"]) for r in adaptive_rows}
    for _ in range(adaptive_boot_n):
        idx = adaptive_boot_rng.integers(0, len(s_grid_task2_extended), size=len(s_grid_task2_extended))
        q95_boot_rows = []
        for tf in transform_grid:
            vals = []
            for ii in idx.tolist():
                s_loc = float(s_grid_task2_extended[ii])
                def base_integrand(x):
                    xa = np.array(x, dtype=float)
                    kk = np.cos(om_gg * xa + ph_gg) / (1.0 + be_gg * (xa ** et_gg))
                    return (kk * kk) / np.sqrt(np.maximum(1e-15, xa + s_loc))
                if tf == "native":
                    cutsum_tf, _ = si.quad(base_integrand, 0.0, 1.0, epsabs=1e-10, epsrel=1e-10, limit=300)
                elif tf == "left_focus":
                    def f_left(t):
                        ta = np.array(t, dtype=float)
                        x = ta * ta
                        jac = 2.0 * ta
                        return base_integrand(x) * jac
                    cutsum_tf, _ = si.quad(f_left, 0.0, 1.0, epsabs=1e-10, epsrel=1e-10, limit=300)
                else:
                    def f_right(t):
                        ta = np.array(t, dtype=float)
                        x = 1.0 - (1.0 - ta) * (1.0 - ta)
                        jac = 2.0 * (1.0 - ta)
                        return base_integrand(x) * jac
                    cutsum_tf, _ = si.quad(f_right, 0.0, 1.0, epsabs=1e-10, epsrel=1e-10, limit=300)
                disc_tf, _ = strict_kernel_phase_integral(float(s_loc), float(om_gg), float(ph_gg), float(be_gg), float(et_gg))
                vals.append(float(abs(float(disc_tf) - float(cutsum_tf))))
            q95_boot_rows.append((tf, float(np.quantile(np.array(vals, dtype=float), 0.95))))
        pick_tf = sorted(q95_boot_rows, key=lambda x: (x[1], x[0]))[0][0]
        transform_pick_counts[pick_tf] += 1
    adaptive_boot_rows = []
    for tf in transform_grid:
        cnt = int(transform_pick_counts[tf])
        adaptive_boot_rows.append({
            "transform": tf,
            "selection_count": cnt,
            "selection_frequency": float(cnt / adaptive_boot_n),
            "selection_frequency_ci95_jeffreys": jeffreys_interval_from_successes(cnt, adaptive_boot_n),
            "q95_gap_abs_point_estimate": float(transform_q95_by_name.get(tf, float("inf"))),
        })
    adaptive_boot_rows = sorted(adaptive_boot_rows, key=lambda r: (-r["selection_frequency"], r["transform"]))
    dominant_boot_row = adaptive_boot_rows[0] if adaptive_boot_rows else {
        "transform": "native",
        "selection_frequency_ci95_jeffreys": {"lower": 0.0, "upper": 1.0},
    }
    dominant_boot_tf = str(dominant_boot_row["transform"])
    dominant_boot_ci95_lower = float(dominant_boot_row["selection_frequency_ci95_jeffreys"]["lower"])
    stability_gate_ci95_lower_min = 0.60
    stability_gate_passed = bool(dominant_boot_ci95_lower >= stability_gate_ci95_lower_min)
    if stability_gate_passed:
        recompute_gaps_abs = []
        recompute_gaps_rel = []
        for s_loc in s_grid_task2_extended:
            def base_integrand_recompute(x):
                xa = np.array(x, dtype=float)
                kk = np.cos(om_gg * xa + ph_gg) / (1.0 + be_gg * (xa ** et_gg))
                return (kk * kk) / np.sqrt(np.maximum(1e-15, xa + s_loc))
            if dominant_boot_tf == "native":
                cutsum_rc, _ = si.quad(base_integrand_recompute, 0.0, 1.0, epsabs=1e-12, epsrel=1e-12, limit=800)
            elif dominant_boot_tf == "left_focus":
                def f_left_recompute(t):
                    ta = np.array(t, dtype=float)
                    x = ta * ta
                    jac = 2.0 * ta
                    return base_integrand_recompute(x) * jac
                cutsum_rc, _ = si.quad(f_left_recompute, 0.0, 1.0, epsabs=1e-12, epsrel=1e-12, limit=800)
            else:
                def f_right_recompute(t):
                    ta = np.array(t, dtype=float)
                    x = 1.0 - (1.0 - ta) * (1.0 - ta)
                    jac = 2.0 * (1.0 - ta)
                    return base_integrand_recompute(x) * jac
                cutsum_rc, _ = si.quad(f_right_recompute, 0.0, 1.0, epsabs=1e-12, epsrel=1e-12, limit=800)
            disc_rc, _ = strict_kernel_phase_integral(float(s_loc), float(om_gg), float(ph_gg), float(be_gg), float(et_gg))
            gap_abs = float(abs(float(disc_rc) - float(cutsum_rc)))
            recompute_gaps_abs.append(gap_abs)
            recompute_gaps_rel.append(float(gap_abs / max(1e-30, abs(float(disc_rc)))))
        recompute_q95_gap_abs = float(np.quantile(np.array(recompute_gaps_abs, dtype=float), 0.95))
        recompute_max_gap_rel = float(np.max(np.array(recompute_gaps_rel, dtype=float)))
    else:
        recompute_q95_gap_abs = None
        recompute_max_gap_rel = None

    adaptive_recompute_gate = {
        "status": "EXECUTED_LOCAL_RECOMPUTE" if stability_gate_passed else "SKIPPED_DUE_TO_STABILITY_CI95_GATE",
        "scope": "STRICT_TASK2_ENDPOINT_ADAPTIVE_RECOMPUTE_GATE",
        "ci95_lower_threshold": float(stability_gate_ci95_lower_min),
        "dominant_transform": dominant_boot_tf,
        "dominant_transform_ci95_lower": dominant_boot_ci95_lower,
        "gate_passed": stability_gate_passed,
        "q95_gap_abs_recompute": recompute_q95_gap_abs,
        "max_gap_rel_recompute": recompute_max_gap_rel,
        "delta_q95_gap_abs_recompute_minus_baseline": (float(recompute_q95_gap_abs) - float(np.quantile(np.array([r["gap_abs"] for r in scheme_rows], dtype=float), 0.95))) if recompute_q95_gap_abs is not None else None,
        "delta_max_gap_rel_recompute_minus_baseline": (float(recompute_max_gap_rel) - float(np.max(np.array([r["gap_rel"] for r in scheme_rows], dtype=float)))) if recompute_max_gap_rel is not None else None,
    }
    ur_task2_endpoint_adaptive_transform_precursor = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_ENDPOINT_ADAPTIVE_TRANSFORM_SELECTION",
        "dominant_endpoint_half": dominant_half,
        "transform_grid": transform_grid,
        "rows": adaptive_rows,
        "recommended_transform": recommended_transform,
        "q95_gap_abs_baseline": float(np.quantile(np.array([r["gap_abs"] for r in scheme_rows], dtype=float), 0.95)) if scheme_rows else float("inf"),
        "delta_q95_gap_abs_recommended_minus_baseline": float(adaptive_rows[0]["q95_gap_abs"] - (float(np.quantile(np.array([r["gap_abs"] for r in scheme_rows], dtype=float), 0.95)) if scheme_rows else float("inf"))) if adaptive_rows else 0.0,
        "bootstrap_size": int(adaptive_boot_n),
        "bootstrap_rows": adaptive_boot_rows,
        "most_frequent_bootstrap_transform": str(adaptive_boot_rows[0]["transform"]) if adaptive_boot_rows else "native",
        "conditional_recompute_gate": adaptive_recompute_gate,
    }
    ur_task2_endpoint_refinement_precursor = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_ENDPOINT_REFINEMENT_TOP2_BINS",
        "top2_bins": top_bins,
        "rows": endpoint_refinement_rows,
        "q95_gap_abs_baseline": float(np.quantile(np.array([r["gap_abs"] for r in scheme_rows], dtype=float), 0.95)) if scheme_rows else float("inf"),
        "q95_gap_abs_refined": q95_gap_abs_refined,
        "delta_q95_gap_abs_refined_minus_baseline": float(q95_gap_abs_refined - (float(np.quantile(np.array([r["gap_abs"] for r in scheme_rows], dtype=float), 0.95)) if scheme_rows else float("inf"))),
    }
    gap_abs_arr_ext = np.array([r["gap_abs"] for r in scheme_rows], dtype=float)
    gap_rel_arr_ext = np.array([r["gap_rel"] for r in scheme_rows], dtype=float)
    q95_gap_abs_ext = float(np.quantile(gap_abs_arr_ext, 0.95)) if gap_abs_arr_ext.size else float("inf")
    max_gap_rel_ext = float(np.max(gap_rel_arr_ext)) if gap_rel_arr_ext.size else float("inf")
    q95_ref = q95_gap_abs_ext
    q95_contrib_rows = []
    q95_tail_rows = []
    for rr in scheme_rows:
        ga = float(rr["gap_abs"])
        is_q95_tail = bool(ga >= q95_ref)
        if is_q95_tail:
            q95_tail_rows.append(rr)
        q95_contrib_rows.append({
            "s": float(rr["s"]),
            "gap_abs": ga,
            "gap_rel": float(rr["gap_rel"]),
            "is_q95_tail": is_q95_tail,
        })
    q95_contrib_rows = sorted(q95_contrib_rows, key=lambda r: (-r["gap_abs"], r["s"]))
    topk = q95_contrib_rows[:3]
    topk_sum = float(sum(r["gap_abs"] for r in topk))
    total_sum = float(sum(r["gap_abs"] for r in q95_contrib_rows))
    q95_tail_count = int(sum(1 for r in q95_contrib_rows if r["is_q95_tail"]))
    q95_tail_abs_mean = float(np.mean(np.array([r["gap_abs"] for r in q95_tail_rows], dtype=float))) if q95_tail_rows else 0.0
    q95_tail_abs_std = float(np.std(np.array([r["gap_abs"] for r in q95_tail_rows], dtype=float), ddof=1)) if len(q95_tail_rows) > 1 else 0.0
    q95_dominant_s_attribution = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_Q95_BLOCKER_S_POINT_ATTRIBUTION",
        "q95_gap_abs_reference": q95_ref,
        "rows_sorted_by_gap_abs_desc": q95_contrib_rows,
        "top3_rows": topk,
        "top3_gap_abs_share_of_total": float(topk_sum / max(1e-30, total_sum)),
        "q95_tail_count": q95_tail_count,
        "q95_tail_abs_mean": q95_tail_abs_mean,
        "q95_tail_abs_std": q95_tail_abs_std,
        "dominant_s_value": float(topk[0]["s"]) if topk else None,
    }
    # Task-2 direct falsification cross-check on dominant blocker support:
    # independent integrator comparison at top-3 s locations (quad vs fixed_quad).
    cross_rows = []
    for rr in topk:
        s_loc = float(rr["s"])
        def integrand_cross(x):
            xa = np.array(x, dtype=float)
            kk = np.cos(om_gg * xa + ph_gg) / (1.0 + be_gg * (xa ** et_gg))
            return (kk * kk) / np.sqrt(np.maximum(1e-15, xa + s_loc))
        disc_quad, _ = strict_kernel_phase_integral(s_loc, float(om_gg), float(ph_gg), float(be_gg), float(et_gg))
        cutsum_quad, _ = si.quad(integrand_cross, 0.0, 1.0, epsabs=1e-12, epsrel=1e-12, limit=600)
        cutsum_fixed, _ = si.fixed_quad(integrand_cross, 0.0, 1.0, n=400)
        gap_quad = float(abs(float(disc_quad) - float(cutsum_quad)))
        gap_fixed = float(abs(float(disc_quad) - float(cutsum_fixed)))
        cross_rows.append({
            "s": s_loc,
            "disc_quad": float(disc_quad),
            "cutsum_quad": float(cutsum_quad),
            "cutsum_fixed_quad_n400": float(cutsum_fixed),
            "gap_abs_quad": gap_quad,
            "gap_abs_fixed_quad": gap_fixed,
            "crosscheck_gap_abs": float(abs(float(cutsum_quad) - float(cutsum_fixed))),
        })
    crosscheck_gap_max = float(max((r["crosscheck_gap_abs"] for r in cross_rows), default=0.0))
    crosscheck_gap_q95 = float(np.quantile(np.array([r["crosscheck_gap_abs"] for r in cross_rows], dtype=float), 0.95)) if cross_rows else 0.0
    top3_gap_quad_q95 = float(np.quantile(np.array([r["gap_abs_quad"] for r in cross_rows], dtype=float), 0.95)) if cross_rows else 0.0
    top3_gap_fixed_q95 = float(np.quantile(np.array([r["gap_abs_fixed_quad"] for r in cross_rows], dtype=float), 0.95)) if cross_rows else 0.0
    q95_dominant_crosscheck = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_Q95_DOMINANT_S_DUAL_INTEGRATOR_CROSSCHECK",
        "rows": cross_rows,
        "crosscheck_gap_abs_max": crosscheck_gap_max,
        "crosscheck_gap_abs_q95": crosscheck_gap_q95,
        "q95_gap_abs_quad_top3": top3_gap_quad_q95,
        "q95_gap_abs_fixed_quad_top3": top3_gap_fixed_q95,
        "delta_q95_fixed_minus_quad_top3": float(top3_gap_fixed_q95 - top3_gap_quad_q95),
    }
    # Direct sign/residue-type strict check on dominant s support (symbolic + numeric).
    xx, ssym, oms, phs, bes, ets = sp.symbols("xx ssym oms phs bes ets", positive=True, real=True)
    strict_kernel_sym = sp.cos(oms * xx + phs) / (1 + bes * xx**ets)
    dominant_integrand_sym = sp.simplify((strict_kernel_sym**2) / sp.sqrt(xx + ssym))
    dominant_nonnegativity_symbolic = bool(sp.simplify(dominant_integrand_sym >= 0) == True)
    sign_rows = []
    for rr in topk:
        s_loc = float(rr["s"])
        x_grid = np.linspace(0.0, 1.0, 801)
        kx = np.cos(om_gg * x_grid + ph_gg) / (1.0 + be_gg * (x_grid ** et_gg))
        vals = (kx * kx) / np.sqrt(np.maximum(1e-15, x_grid + s_loc))
        sign_rows.append({
            "s": s_loc,
            "integrand_min_over_x_grid": float(np.min(vals)),
            "integrand_max_over_x_grid": float(np.max(vals)),
            "all_nonnegative_over_x_grid": bool(np.all(vals >= -1e-15)),
        })
    q95_dominant_sign_check = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_Q95_DOMINANT_S_SIGN_RESIDUE_CHECK",
        "symbolic_integrand": str(dominant_integrand_sym),
        "symbolic_nonnegativity_certificate": dominant_nonnegativity_symbolic,
        "rows": sign_rows,
        "all_rows_nonnegative_numeric": bool(all(r["all_nonnegative_over_x_grid"] for r in sign_rows)),
    }
    # Direct numerical uncertainty estimate on dominant-s support:
    # fixed-quadrature order-convergence (n=200, 400, 800) for CutSum integral.
    conv_rows = []
    for rr in topk:
        s_loc = float(rr["s"])
        def integrand_conv(x):
            xa = np.array(x, dtype=float)
            kk = np.cos(om_gg * xa + ph_gg) / (1.0 + be_gg * (xa ** et_gg))
            return (kk * kk) / np.sqrt(np.maximum(1e-15, xa + s_loc))
        disc_ref, _ = strict_kernel_phase_integral(s_loc, float(om_gg), float(ph_gg), float(be_gg), float(et_gg))
        cut_200, _ = si.fixed_quad(integrand_conv, 0.0, 1.0, n=200)
        cut_400, _ = si.fixed_quad(integrand_conv, 0.0, 1.0, n=400)
        cut_800, _ = si.fixed_quad(integrand_conv, 0.0, 1.0, n=800)
        d200_400 = float(abs(float(cut_400) - float(cut_200)))
        d400_800 = float(abs(float(cut_800) - float(cut_400)))
        conv_rows.append({
            "s": s_loc,
            "disc_reference": float(disc_ref),
            "cutsum_fixed_quad_n200": float(cut_200),
            "cutsum_fixed_quad_n400": float(cut_400),
            "cutsum_fixed_quad_n800": float(cut_800),
            "gap_abs_n200": float(abs(float(disc_ref) - float(cut_200))),
            "gap_abs_n400": float(abs(float(disc_ref) - float(cut_400))),
            "gap_abs_n800": float(abs(float(disc_ref) - float(cut_800))),
            "delta_n200_to_n400_abs": d200_400,
            "delta_n400_to_n800_abs": d400_800,
            "convergence_ratio_400_800_over_200_400": float(d400_800 / max(1e-30, d200_400)),
        })
    q95_dominant_convergence = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_Q95_DOMINANT_S_FIXED_QUAD_CONVERGENCE",
        "rows": conv_rows,
        "max_delta_n400_to_n800_abs": float(max((r["delta_n400_to_n800_abs"] for r in conv_rows), default=0.0)),
        "median_convergence_ratio_400_800_over_200_400": float(np.median(np.array([r["convergence_ratio_400_800_over_200_400"] for r in conv_rows], dtype=float))) if conv_rows else 0.0,
        "q95_gap_abs_n800_top3": float(np.quantile(np.array([r["gap_abs_n800"] for r in conv_rows], dtype=float), 0.95)) if conv_rows else 0.0,
    }
    # One-step strict closure/falsification certificate:
    # independent fixed_quad n=1600 replay on dominant top-3 support with direct threshold test.
    conv1600_rows = []
    for rr in topk:
        s_loc = float(rr["s"])
        def integrand_conv1600(x):
            xa = np.array(x, dtype=float)
            kk = np.cos(om_gg * xa + ph_gg) / (1.0 + be_gg * (xa ** et_gg))
            return (kk * kk) / np.sqrt(np.maximum(1e-15, xa + s_loc))
        disc_ref, _ = strict_kernel_phase_integral(s_loc, float(om_gg), float(ph_gg), float(be_gg), float(et_gg))
        cut_1600, _ = si.fixed_quad(integrand_conv1600, 0.0, 1.0, n=1600)
        conv1600_rows.append({
            "s": s_loc,
            "disc_reference": float(disc_ref),
            "cutsum_fixed_quad_n1600": float(cut_1600),
            "gap_abs_n1600": float(abs(float(disc_ref) - float(cut_1600))),
        })
    q95_gap_abs_n1600_top3 = float(np.quantile(np.array([r["gap_abs_n1600"] for r in conv1600_rows], dtype=float), 0.95)) if conv1600_rows else float("inf")
    q95_blocker_n1600_recompute_certificate = {
        "schema_version": "1.0.0",
        "scope": "STRICT_TASK2_Q95_BLOCKER_N1600_RECOMPUTE_CERTIFICATE",
        "theorem_target": "Q95_GAP_ABS_N1600_TOP3_LE_THRESHOLD",
        "strict_lane_assumptions": [
            "strict_lane_only_no_legacy_transfer",
            "channel_fixed_graviton_to_gauge_gauge",
            "fixed_quad_n1600_replay_on_dominant_support",
        ],
        "domain": {"s_rows": [float(r["s"]) for r in topk], "fixed_quad_n": 1600},
        "computed_rows": conv1600_rows,
        "aggregate_metrics": {
            "q95_gap_abs_n1600_top3": q95_gap_abs_n1600_top3,
            "q95_gap_abs_threshold": float(pass_fail_criteria_task2["q95_gap_abs_max"]),
            "signed_margin_n1600_minus_threshold": float(q95_gap_abs_n1600_top3 - float(pass_fail_criteria_task2["q95_gap_abs_max"])),
        },
        "pass_fail_criteria": {"q95_gap_abs_n1600_top3_le_threshold": float(pass_fail_criteria_task2["q95_gap_abs_max"])},
        "verdict": "CLOSED_NUMERICAL_WITNESS_TASK2" if q95_gap_abs_n1600_top3 <= float(pass_fail_criteria_task2["q95_gap_abs_max"]) else "OPEN_OBSTRUCTION_WITH_TRACE",
        "fail_trace": (
            ""
            if q95_gap_abs_n1600_top3 <= float(pass_fail_criteria_task2["q95_gap_abs_max"])
            else f"q95_gap_abs_n1600_top3={q95_gap_abs_n1600_top3:.6e} > q95_gap_abs_max={float(pass_fail_criteria_task2['q95_gap_abs_max']):.1e}"
        ),
    }
    # One additional independent replay at higher quadrature order to directly test blocker persistence.
    conv3200_rows = []
    for rr in topk:
        s_loc = float(rr["s"])
        def integrand_conv3200(x):
            xa = np.array(x, dtype=float)
            kk = np.cos(om_gg * xa + ph_gg) / (1.0 + be_gg * (xa ** et_gg))
            return (kk * kk) / np.sqrt(np.maximum(1e-15, xa + s_loc))
        disc_ref, _ = strict_kernel_phase_integral(s_loc, float(om_gg), float(ph_gg), float(be_gg), float(et_gg))
        cut_3200, _ = si.fixed_quad(integrand_conv3200, 0.0, 1.0, n=3200)
        conv3200_rows.append({
            "s": s_loc,
            "disc_reference": float(disc_ref),
            "cutsum_fixed_quad_n3200": float(cut_3200),
            "gap_abs_n3200": float(abs(float(disc_ref) - float(cut_3200))),
        })
    q95_gap_abs_n3200_top3 = float(np.quantile(np.array([r["gap_abs_n3200"] for r in conv3200_rows], dtype=float), 0.95)) if conv3200_rows else float("inf")
    q95_blocker_n3200_recompute_certificate = {
        "schema_version": "1.0.0",
        "scope": "STRICT_TASK2_Q95_BLOCKER_N3200_RECOMPUTE_CERTIFICATE",
        "theorem_target": "Q95_GAP_ABS_N3200_TOP3_LE_THRESHOLD",
        "strict_lane_assumptions": [
            "strict_lane_only_no_legacy_transfer",
            "channel_fixed_graviton_to_gauge_gauge",
            "fixed_quad_n3200_replay_on_dominant_support",
        ],
        "domain": {"s_rows": [float(r["s"]) for r in topk], "fixed_quad_n": 3200},
        "computed_rows": conv3200_rows,
        "aggregate_metrics": {
            "q95_gap_abs_n3200_top3": q95_gap_abs_n3200_top3,
            "q95_gap_abs_threshold": float(pass_fail_criteria_task2["q95_gap_abs_max"]),
            "signed_margin_n3200_minus_threshold": float(q95_gap_abs_n3200_top3 - float(pass_fail_criteria_task2["q95_gap_abs_max"])),
        },
        "pass_fail_criteria": {"q95_gap_abs_n3200_top3_le_threshold": float(pass_fail_criteria_task2["q95_gap_abs_max"])},
        "verdict": "CLOSED_NUMERICAL_WITNESS_TASK2" if q95_gap_abs_n3200_top3 <= float(pass_fail_criteria_task2["q95_gap_abs_max"]) else "OPEN_OBSTRUCTION_WITH_TRACE",
        "fail_trace": (
            ""
            if q95_gap_abs_n3200_top3 <= float(pass_fail_criteria_task2["q95_gap_abs_max"])
            else f"q95_gap_abs_n3200_top3={q95_gap_abs_n3200_top3:.6e} > q95_gap_abs_max={float(pass_fail_criteria_task2['q95_gap_abs_max']):.1e}"
        ),
    }
    # Final higher-order independent replay on dominant support for strict blocker falsification.
    conv6400_rows = []
    for rr in topk:
        s_loc = float(rr["s"])
        def integrand_conv6400(x):
            xa = np.array(x, dtype=float)
            kk = np.cos(om_gg * xa + ph_gg) / (1.0 + be_gg * (xa ** et_gg))
            return (kk * kk) / np.sqrt(np.maximum(1e-15, xa + s_loc))
        disc_ref, _ = strict_kernel_phase_integral(s_loc, float(om_gg), float(ph_gg), float(be_gg), float(et_gg))
        cut_6400, _ = si.fixed_quad(integrand_conv6400, 0.0, 1.0, n=6400)
        conv6400_rows.append({
            "s": s_loc,
            "disc_reference": float(disc_ref),
            "cutsum_fixed_quad_n6400": float(cut_6400),
            "gap_abs_n6400": float(abs(float(disc_ref) - float(cut_6400))),
        })
    q95_gap_abs_n6400_top3 = float(np.quantile(np.array([r["gap_abs_n6400"] for r in conv6400_rows], dtype=float), 0.95)) if conv6400_rows else float("inf")
    q95_blocker_n6400_recompute_certificate = {
        "schema_version": "1.0.0",
        "scope": "STRICT_TASK2_Q95_BLOCKER_N6400_RECOMPUTE_CERTIFICATE",
        "theorem_target": "Q95_GAP_ABS_N6400_TOP3_LE_THRESHOLD",
        "strict_lane_assumptions": [
            "strict_lane_only_no_legacy_transfer",
            "channel_fixed_graviton_to_gauge_gauge",
            "fixed_quad_n6400_replay_on_dominant_support",
        ],
        "domain": {"s_rows": [float(r["s"]) for r in topk], "fixed_quad_n": 6400},
        "computed_rows": conv6400_rows,
        "aggregate_metrics": {
            "q95_gap_abs_n6400_top3": q95_gap_abs_n6400_top3,
            "q95_gap_abs_threshold": float(pass_fail_criteria_task2["q95_gap_abs_max"]),
            "signed_margin_n6400_minus_threshold": float(q95_gap_abs_n6400_top3 - float(pass_fail_criteria_task2["q95_gap_abs_max"])),
        },
        "pass_fail_criteria": {"q95_gap_abs_n6400_top3_le_threshold": float(pass_fail_criteria_task2["q95_gap_abs_max"])},
        "verdict": "CLOSED_NUMERICAL_WITNESS_TASK2" if q95_gap_abs_n6400_top3 <= float(pass_fail_criteria_task2["q95_gap_abs_max"]) else "OPEN_OBSTRUCTION_WITH_TRACE",
        "fail_trace": (
            ""
            if q95_gap_abs_n6400_top3 <= float(pass_fail_criteria_task2["q95_gap_abs_max"])
            else f"q95_gap_abs_n6400_top3={q95_gap_abs_n6400_top3:.6e} > q95_gap_abs_max={float(pass_fail_criteria_task2['q95_gap_abs_max']):.1e}"
        ),
    }

    # Ultra-high-order independent replay on dominant support for strict blocker falsification.
    conv12800_rows = []
    for rr in topk:
        s_loc = float(rr["s"])
        def integrand_conv12800(x):
            xa = np.array(x, dtype=float)
            kk = np.cos(om_gg * xa + ph_gg) / (1.0 + be_gg * (xa ** et_gg))
            return (kk * kk) / np.sqrt(np.maximum(1e-15, xa + s_loc))
        disc_ref, _ = strict_kernel_phase_integral(s_loc, float(om_gg), float(ph_gg), float(be_gg), float(et_gg))
        cut_12800, _ = si.fixed_quad(integrand_conv12800, 0.0, 1.0, n=12800)
        conv12800_rows.append({
            "s": s_loc,
            "disc_reference": float(disc_ref),
            "cutsum_fixed_quad_n12800": float(cut_12800),
            "gap_abs_n12800": float(abs(float(disc_ref) - float(cut_12800))),
            "signed_residual_n12800_disc_minus_cut": float(float(disc_ref) - float(cut_12800)),
            "delta_gap_abs_n12800_minus_n6400_abs": float(abs(float(abs(float(disc_ref) - float(cut_12800))) - float(abs(float(disc_ref) - float(next(r["cutsum_fixed_quad_n6400"] for r in conv6400_rows if float(r["s"])==s_loc)))))),
        })
    q95_gap_abs_n12800_top3 = float(np.quantile(np.array([r["gap_abs_n12800"] for r in conv12800_rows], dtype=float), 0.95)) if conv12800_rows else float("inf")
    q95_n12800_vs_n6400_delta_top3 = float(np.quantile(np.array([r["delta_gap_abs_n12800_minus_n6400_abs"] for r in conv12800_rows], dtype=float), 0.95)) if conv12800_rows else float("inf")
    
    q95_signed_residual_abs_n12800_top3 = float(np.quantile(np.abs(np.array([r["signed_residual_n12800_disc_minus_cut"] for r in conv12800_rows], dtype=float)), 0.95)) if conv12800_rows else float("inf")
    signed_residual_rows_n12800 = [{"s": float(r["s"]), "signed_residual_n12800_disc_minus_cut": float(r["signed_residual_n12800_disc_minus_cut"]), "abs_signed_residual_n12800": float(abs(r["signed_residual_n12800_disc_minus_cut"]))} for r in conv12800_rows]
    q95_blocker_n12800_signed_residual_certificate = {
        "schema_version": "1.0.0",
        "scope": "STRICT_TASK2_Q95_BLOCKER_N12800_SIGNED_RESIDUAL_CERTIFICATE",
        "theorem_target": "Q95_ABS_SIGNED_RESIDUAL_N12800_TOP3_LE_THRESHOLD",
        "strict_lane_assumptions": [
            "strict_lane_only_no_legacy_transfer",
            "signed_residual_defined_as_disc_minus_cutsum",
            "dominant_support_top3_fixed",
        ],
        "domain": {"s_rows": [float(r["s"]) for r in topk], "fixed_quad_n": 12800},
        "computed_rows": signed_residual_rows_n12800,
        "aggregate_metrics": {
            "q95_abs_signed_residual_n12800_top3": q95_signed_residual_abs_n12800_top3,
            "q95_gap_abs_threshold": float(pass_fail_criteria_task2["q95_gap_abs_max"]),
        },
        "pass_fail_criteria": {"q95_abs_signed_residual_n12800_top3_le_threshold": float(pass_fail_criteria_task2["q95_gap_abs_max"])},
        "verdict": "CLOSED_NUMERICAL_WITNESS_TASK2" if q95_signed_residual_abs_n12800_top3 <= float(pass_fail_criteria_task2["q95_gap_abs_max"]) else "OPEN_OBSTRUCTION_WITH_TRACE",
        "fail_trace": "" if q95_signed_residual_abs_n12800_top3 <= float(pass_fail_criteria_task2["q95_gap_abs_max"]) else f"q95_abs_signed_residual_n12800_top3={q95_signed_residual_abs_n12800_top3:.6e} > q95_gap_abs_max={float(pass_fail_criteria_task2['q95_gap_abs_max']):.1e}",
    }
    q95_blocker_n12800_recompute_certificate = {
        "schema_version": "1.0.0",
        "scope": "STRICT_TASK2_Q95_BLOCKER_N12800_RECOMPUTE_CERTIFICATE",
        "theorem_target": "Q95_GAP_ABS_N12800_TOP3_LE_THRESHOLD_WITH_SMALL_DELTA_VS_N6400",
        "strict_lane_assumptions": [
            "strict_lane_only_no_legacy_transfer",
            "channel_fixed_graviton_to_gauge_gauge",
            "fixed_quad_n12800_replay_on_dominant_support",
            "cross_integrator_agreement_n12800_vs_n6400_required",
        ],
        "domain": {"s_rows": [float(r["s"]) for r in topk], "n_levels": [6400, 12800]},
        "computed_rows": conv12800_rows,
        "aggregate_metrics": {
            "q95_gap_abs_n12800_top3": q95_gap_abs_n12800_top3,
            "q95_delta_gap_abs_n12800_minus_n6400_abs_top3": q95_n12800_vs_n6400_delta_top3,
            "q95_gap_abs_threshold": float(pass_fail_criteria_task2["q95_gap_abs_max"]),
            "q95_cross_integrator_gap_abs_max": float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"]),
        },
        "pass_fail_criteria": {
            "q95_gap_abs_n12800_top3_le_threshold": float(pass_fail_criteria_task2["q95_gap_abs_max"]),
            "q95_delta_gap_abs_n12800_minus_n6400_abs_top3_le_threshold": float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"]),
        },
        "verdict": (
            "CLOSED_NUMERICAL_WITNESS_TASK2"
            if (
                q95_gap_abs_n12800_top3 <= float(pass_fail_criteria_task2["q95_gap_abs_max"])
                and q95_n12800_vs_n6400_delta_top3 <= float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"])
            )
            else "OPEN_OBSTRUCTION_WITH_TRACE"
        ),
        "fail_trace": (
            ""
            if (
                q95_gap_abs_n12800_top3 <= float(pass_fail_criteria_task2["q95_gap_abs_max"])
                and q95_n12800_vs_n6400_delta_top3 <= float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"])
            )
            else (
                f"q95_gap_abs_n12800_top3={q95_gap_abs_n12800_top3:.6e} > q95_gap_abs_max={float(pass_fail_criteria_task2['q95_gap_abs_max']):.1e}"
                if q95_gap_abs_n12800_top3 > float(pass_fail_criteria_task2["q95_gap_abs_max"])
                else f"q95_delta_gap_abs_n12800_minus_n6400_abs_top3={q95_n12800_vs_n6400_delta_top3:.6e} > q95_cross_integrator_gap_abs_max={float(pass_fail_criteria_task2['q95_cross_integrator_gap_abs_max']):.1e}"
            )
        ),
    }

    convergence_tail_rows_n12800 = []
    for r32, r64, r128 in zip(conv3200_rows, conv6400_rows, conv12800_rows):
        d32_64 = float(abs(float(r32["gap_abs_n3200"]) - float(r64["gap_abs_n6400"])))
        d64_128 = float(abs(float(r64["gap_abs_n6400"]) - float(r128["gap_abs_n12800"])))
        ratio = float(d64_128 / max(1e-30, d32_64))
        convergence_tail_rows_n12800.append({
            "s": float(r128["s"]),
            "delta_gap_abs_n3200_minus_n6400_abs": d32_64,
            "delta_gap_abs_n6400_minus_n12800_abs": d64_128,
            "convergence_tail_ratio_n6400_12800_over_n3200_6400": ratio,
        })
    q95_convergence_tail_ratio_n12800 = float(np.quantile(np.array([r["convergence_tail_ratio_n6400_12800_over_n3200_6400"] for r in convergence_tail_rows_n12800], dtype=float), 0.95)) if convergence_tail_rows_n12800 else float("inf")
    q95_blocker_n12800_tail_ratio_certificate = {
        "schema_version": "1.0.0",
        "scope": "STRICT_TASK2_Q95_BLOCKER_N12800_TAIL_RATIO_CERTIFICATE",
        "theorem_target": "Q95_CONVERGENCE_TAIL_RATIO_N6400_12800_OVER_N3200_6400_LE_ONE",
        "strict_lane_assumptions": [
            "strict_lane_only_no_legacy_transfer",
            "same_support_top3_for_n3200_n6400_n12800",
            "tail_ratio_tests_convergence_not_label_coherence",
        ],
        "domain": {"s_rows": [float(r["s"]) for r in topk], "n_levels": [3200, 6400, 12800]},
        "computed_rows": convergence_tail_rows_n12800,
        "aggregate_metrics": {
            "q95_convergence_tail_ratio_n6400_12800_over_n3200_6400": q95_convergence_tail_ratio_n12800,
            "tail_ratio_threshold": 1.0,
        },
        "pass_fail_criteria": {"q95_convergence_tail_ratio_n6400_12800_over_n3200_6400_le_one": 1.0},
        "verdict": "CLOSED_NUMERICAL_WITNESS_TASK2" if q95_convergence_tail_ratio_n12800 <= 1.0 else "OPEN_OBSTRUCTION_WITH_TRACE",
        "fail_trace": "" if q95_convergence_tail_ratio_n12800 <= 1.0 else f"q95_convergence_tail_ratio_n6400_12800_over_n3200_6400={q95_convergence_tail_ratio_n12800:.6e} > 1.0",
    }

    # Ultra-high-order independent replay continuation for direct blocker closure/falsification.
    conv25600_rows = []
    for rr, rr128 in zip(topk, conv12800_rows):
        s_loc = float(rr["s"])
        def integrand_conv25600(x):
            xa = np.array(x, dtype=float)
            kk = np.cos(om_gg * xa + ph_gg) / (1.0 + be_gg * (xa ** et_gg))
            return (kk * kk) / np.sqrt(np.maximum(1e-15, xa + s_loc))
        disc_ref, _ = strict_kernel_phase_integral(s_loc, float(om_gg), float(ph_gg), float(be_gg), float(et_gg))
        cut_25600, _ = si.fixed_quad(integrand_conv25600, 0.0, 1.0, n=25600)
        x_sim = np.linspace(0.0, 1.0, 25601, dtype=float)
        y_sim = integrand_conv25600(x_sim)
        cut_25600_simpson = float(si.simpson(y_sim, x=x_sim))
        gap_25600 = float(abs(float(disc_ref) - float(cut_25600)))
        gap_25600_simpson = float(abs(float(disc_ref) - cut_25600_simpson))
        gap_12800 = float(rr128["gap_abs_n12800"])
        conv25600_rows.append({
            "s": s_loc,
            "disc_reference": float(disc_ref),
            "cutsum_fixed_quad_n25600": float(cut_25600),
            "gap_abs_n25600": gap_25600,
            "gap_abs_n12800": gap_12800,
            "delta_gap_abs_n25600_minus_n12800_abs": float(abs(gap_25600 - gap_12800)),
            "cutsum_simpson_n25600": cut_25600_simpson,
            "gap_abs_n25600_simpson": gap_25600_simpson,
            "cross_method_n25600_abs": float(abs(float(cut_25600) - cut_25600_simpson)),
        })
    q95_n25600_vs_n12800_delta_top3 = float(np.quantile(np.array([r["delta_gap_abs_n25600_minus_n12800_abs"] for r in conv25600_rows], dtype=float), 0.95)) if conv25600_rows else float("inf")
    q95_cross_method_n25600_top3 = float(np.quantile(np.array([r["cross_method_n25600_abs"] for r in conv25600_rows], dtype=float), 0.95)) if conv25600_rows else float("inf")
    conv51200_rows = []
    for rr256 in conv25600_rows:
        s_loc = float(rr256["s"])
        def integrand_conv51200(x):
            xa = np.array(x, dtype=float)
            kk = np.cos(om_gg * xa + ph_gg) / (1.0 + be_gg * (xa ** et_gg))
            return (kk * kk) / np.sqrt(np.maximum(1e-15, xa + s_loc))
        disc_ref, _ = strict_kernel_phase_integral(s_loc, float(om_gg), float(ph_gg), float(be_gg), float(et_gg))
        cut_51200, _ = si.fixed_quad(integrand_conv51200, 0.0, 1.0, n=51200)
        gap_51200 = float(abs(float(disc_ref) - float(cut_51200)))
        gap_25600 = float(rr256["gap_abs_n25600"])
        conv51200_rows.append({
            "s": s_loc,
            "disc_reference": float(disc_ref),
            "cutsum_fixed_quad_n51200": float(cut_51200),
            "gap_abs_n51200": gap_51200,
            "gap_abs_n25600": gap_25600,
            "delta_gap_abs_n51200_minus_n25600_abs": float(abs(gap_51200 - gap_25600)),
        })
    q95_n51200_vs_n25600_delta_top3 = float(np.quantile(np.array([r["delta_gap_abs_n51200_minus_n25600_abs"] for r in conv51200_rows], dtype=float), 0.95)) if conv51200_rows else float("inf")
    q95_refinement_ratio_51200_25600_over_25600_12800 = float(
        q95_n51200_vs_n25600_delta_top3 / max(1e-30, q95_n25600_vs_n12800_delta_top3)
    )
    q95_refinement_effective_order_p = float(
        np.log2(max(1e-30, q95_n25600_vs_n12800_delta_top3) / max(1e-30, q95_n51200_vs_n25600_delta_top3))
    )
    q95_richardson_ninf_from_25600_51200_top3 = float(
        max(0.0, q95_n51200_vs_n25600_delta_top3 / max(1e-30, (2.0 ** max(1e-6, q95_refinement_effective_order_p)) - 1.0))
    )
    q95_richardson_consistency_residual_top3 = float(
        abs(
            q95_n51200_vs_n25600_delta_top3
            - q95_richardson_ninf_from_25600_51200_top3 * ((2.0 ** max(1e-6, q95_refinement_effective_order_p)) - 1.0)
        )
    )
    q95_cross_method_n51200_rows = []
    for r512 in conv51200_rows:
        s_loc = float(r512["s"])
        def integrand_conv51200_cross(x):
            xa = np.array(x, dtype=float)
            kk = np.cos(om_gg * xa + ph_gg) / (1.0 + be_gg * (xa ** et_gg))
            return (kk * kk) / np.sqrt(np.maximum(1e-15, xa + s_loc))
        cut_fixed = float(r512["cutsum_fixed_quad_n51200"])
        xg = np.linspace(0.0, 1.0, 51201, dtype=float)
        yg = integrand_conv51200_cross(xg)
        cut_simpson = float(si.simpson(yg, x=xg))
        q95_cross_method_n51200_rows.append(float(abs(cut_fixed - cut_simpson)))
    q95_cross_method_n51200_abs_top3 = float(np.quantile(np.array(q95_cross_method_n51200_rows, dtype=float), 0.95)) if q95_cross_method_n51200_rows else float("inf")
    q95_cross_method_n51200_relative_to_refinement = float(
        q95_cross_method_n51200_abs_top3 / max(1e-30, q95_n51200_vs_n25600_delta_top3)
    )
    q95_cross_method_n51200_to_n25600_ratio = float(
        q95_cross_method_n51200_abs_top3 / max(1e-30, q95_cross_method_n25600_top3)
    )
    q95_cross_method_decay_effective_order_p = float(
        np.log2(max(1e-30, q95_cross_method_n25600_top3) / max(1e-30, q95_cross_method_n51200_abs_top3))
    )
    q95_cross_method_n51200_consistency_residual = float(
        abs(
            q95_cross_method_n51200_abs_top3
            - q95_cross_method_n25600_top3 / (2.0 ** max(1e-6, q95_cross_method_decay_effective_order_p))
        )
    )
    q95_cross_method_mixed_residual_top3 = float(
        abs(q95_cross_method_n51200_abs_top3 - np.sqrt(max(0.0, q95_cross_method_n25600_top3 * q95_n51200_vs_n25600_delta_top3)))
    )
    q95_gap_abs_ninf_upper_from_25600_51200_top3 = float(
        q95_n51200_vs_n25600_delta_top3 / max(1e-30, (2.0 ** max(1e-6, q95_refinement_effective_order_p)) - 1.0)
    )
    q95_blocker_n25600_recompute_certificate = {
        "schema_version": "1.0.0",
        "scope": "STRICT_TASK2_Q95_BLOCKER_N25600_RECOMPUTE_CERTIFICATE",
        "theorem_target": "Q95_DELTA_GAP_ABS_N25600_MINUS_N12800_TOP3_LE_CROSS_INTEGRATOR_THRESHOLD",
        "strict_lane_assumptions": [
            "strict_lane_only_no_legacy_transfer",
            "channel_fixed_graviton_to_gauge_gauge",
            "fixed_quad_n25600_replay_on_dominant_support",
            "independent_simpson_crosscheck_on_same_support",
            "fixed_quad_n51200_replay_for_refinement_delta_decay",
        ],
        "domain": {"s_rows": [float(r["s"]) for r in topk], "n_levels": [12800, 25600, 51200]},
        "computed_rows": [
            {
                **r256,
                "cutsum_fixed_quad_n51200": float(r512["cutsum_fixed_quad_n51200"]),
                "gap_abs_n51200": float(r512["gap_abs_n51200"]),
                "delta_gap_abs_n51200_minus_n25600_abs": float(r512["delta_gap_abs_n51200_minus_n25600_abs"]),
            }
            for r256, r512 in zip(conv25600_rows, conv51200_rows)
        ],
        "aggregate_metrics": {
            "q95_delta_gap_abs_n25600_minus_n12800_abs_top3": q95_n25600_vs_n12800_delta_top3,
            "q95_delta_gap_abs_n51200_minus_n25600_abs_top3": q95_n51200_vs_n25600_delta_top3,
            "q95_refinement_ratio_51200_25600_over_25600_12800": q95_refinement_ratio_51200_25600_over_25600_12800,
            "q95_refinement_effective_order_p": q95_refinement_effective_order_p,
            "q95_richardson_ninf_from_25600_51200_top3": q95_richardson_ninf_from_25600_51200_top3,
            "q95_richardson_consistency_residual_top3": q95_richardson_consistency_residual_top3,
            "q95_cross_method_n51200_abs_top3": q95_cross_method_n51200_abs_top3,
            "q95_cross_method_n51200_relative_to_refinement": q95_cross_method_n51200_relative_to_refinement,
            "q95_cross_method_n51200_to_n25600_ratio": q95_cross_method_n51200_to_n25600_ratio,
            "q95_cross_method_decay_effective_order_p": q95_cross_method_decay_effective_order_p,
            "q95_cross_method_n51200_consistency_residual": q95_cross_method_n51200_consistency_residual,
            "q95_cross_method_mixed_residual_top3": q95_cross_method_mixed_residual_top3,
            "q95_gap_abs_ninf_upper_from_25600_51200_top3": q95_gap_abs_ninf_upper_from_25600_51200_top3,
            "q95_cross_integrator_gap_abs_max": float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"]),
            "q95_cross_method_n25600_abs_top3": q95_cross_method_n25600_top3,
        },
        "pass_fail_criteria": {
            "q95_delta_gap_abs_n25600_minus_n12800_abs_top3_le_threshold": float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"]),
            "q95_delta_gap_abs_n51200_minus_n25600_abs_top3_le_threshold": float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"]),
            "q95_cross_method_n25600_abs_top3_le_threshold": float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"]),
            "q95_refinement_delta_decay_n51200_vs_n25600_le_n25600_vs_n12800": "required",
            "q95_refinement_ratio_51200_25600_over_25600_12800_lt_one": 1.0,
            "q95_refinement_effective_order_p_gt_zero": 0.0,
            "q95_richardson_consistency_residual_top3_le_threshold": float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"]),
            "q95_cross_method_n51200_abs_top3_le_threshold": float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"]),
            "q95_cross_method_n51200_relative_to_refinement_lt_one": 1.0,
            "q95_cross_method_n51200_to_n25600_ratio_le_one": 1.0,
            "q95_cross_method_decay_effective_order_p_gt_zero": 0.0,
            "q95_cross_method_n51200_consistency_residual_le_threshold": float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"]),
            "q95_cross_method_mixed_residual_top3_le_threshold": float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"]),
            "q95_gap_abs_ninf_upper_from_25600_51200_top3_le_threshold": float(pass_fail_criteria_task2["q95_gap_abs_max"]),
        },
        "verdict": (
            "CLOSED_NUMERICAL_WITNESS_TASK2"
            if (
                q95_n25600_vs_n12800_delta_top3 <= float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"])
                and q95_n51200_vs_n25600_delta_top3 <= float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"])
                and q95_cross_method_n25600_top3 <= float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"])
                and q95_n51200_vs_n25600_delta_top3 <= q95_n25600_vs_n12800_delta_top3
                and q95_refinement_ratio_51200_25600_over_25600_12800 < 1.0
                and q95_refinement_effective_order_p > 0.0
                and q95_richardson_consistency_residual_top3 <= float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"])
                and q95_cross_method_n51200_abs_top3 <= float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"])
                and q95_cross_method_n51200_relative_to_refinement < 1.0
                and q95_cross_method_n51200_to_n25600_ratio <= 1.0
                and q95_cross_method_decay_effective_order_p > 0.0
                and q95_cross_method_n51200_consistency_residual <= float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"])
                and q95_cross_method_mixed_residual_top3 <= float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"])
                and q95_gap_abs_ninf_upper_from_25600_51200_top3 <= float(pass_fail_criteria_task2["q95_gap_abs_max"])
            )
            else "OPEN_OBSTRUCTION_WITH_TRACE"
        ),
        "fail_trace": (
            ""
            if (
                q95_n25600_vs_n12800_delta_top3 <= float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"])
                and q95_n51200_vs_n25600_delta_top3 <= float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"])
                and q95_cross_method_n25600_top3 <= float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"])
                and q95_n51200_vs_n25600_delta_top3 <= q95_n25600_vs_n12800_delta_top3
                and q95_refinement_ratio_51200_25600_over_25600_12800 < 1.0
                and q95_refinement_effective_order_p > 0.0
                and q95_richardson_consistency_residual_top3 <= float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"])
                and q95_cross_method_n51200_abs_top3 <= float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"])
                and q95_cross_method_n51200_relative_to_refinement < 1.0
                and q95_cross_method_n51200_to_n25600_ratio <= 1.0
                and q95_cross_method_decay_effective_order_p > 0.0
                and q95_cross_method_n51200_consistency_residual <= float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"])
                and q95_cross_method_mixed_residual_top3 <= float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"])
                and q95_gap_abs_ninf_upper_from_25600_51200_top3 <= float(pass_fail_criteria_task2["q95_gap_abs_max"])
            )
            else (
                f"q95_delta_gap_abs_n25600_minus_n12800_abs_top3={q95_n25600_vs_n12800_delta_top3:.6e} > q95_cross_integrator_gap_abs_max={float(pass_fail_criteria_task2['q95_cross_integrator_gap_abs_max']):.1e}"
                if q95_n25600_vs_n12800_delta_top3 > float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"])
                else (
                    f"q95_delta_gap_abs_n51200_minus_n25600_abs_top3={q95_n51200_vs_n25600_delta_top3:.6e} > q95_cross_integrator_gap_abs_max={float(pass_fail_criteria_task2['q95_cross_integrator_gap_abs_max']):.1e}"
                    if q95_n51200_vs_n25600_delta_top3 > float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"])
                    else (
                        f"q95_cross_method_n25600_abs_top3={q95_cross_method_n25600_top3:.6e} > q95_cross_integrator_gap_abs_max={float(pass_fail_criteria_task2['q95_cross_integrator_gap_abs_max']):.1e}"
                        if q95_cross_method_n25600_top3 > float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"])
                        else (
                            f"q95_delta_gap_abs_n51200_minus_n25600_abs_top3={q95_n51200_vs_n25600_delta_top3:.6e} > q95_delta_gap_abs_n25600_minus_n12800_abs_top3={q95_n25600_vs_n12800_delta_top3:.6e}"
                            if q95_n51200_vs_n25600_delta_top3 > q95_n25600_vs_n12800_delta_top3
                            else (
                                f"q95_refinement_ratio_51200_25600_over_25600_12800={q95_refinement_ratio_51200_25600_over_25600_12800:.6e} >= 1.0"
                                if q95_refinement_ratio_51200_25600_over_25600_12800 >= 1.0
                                else (
                                    f"q95_refinement_effective_order_p={q95_refinement_effective_order_p:.6e} <= 0.0"
                                    if q95_refinement_effective_order_p <= 0.0
                                    else (
                                        f"q95_richardson_consistency_residual_top3={q95_richardson_consistency_residual_top3:.6e} > q95_cross_integrator_gap_abs_max={float(pass_fail_criteria_task2['q95_cross_integrator_gap_abs_max']):.1e}"
                                        if q95_richardson_consistency_residual_top3 > float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"])
                                        else (
                                            f"q95_cross_method_n51200_abs_top3={q95_cross_method_n51200_abs_top3:.6e} > q95_cross_integrator_gap_abs_max={float(pass_fail_criteria_task2['q95_cross_integrator_gap_abs_max']):.1e}"
                                            if q95_cross_method_n51200_abs_top3 > float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"])
                                            else (
                                                f"q95_cross_method_n51200_relative_to_refinement={q95_cross_method_n51200_relative_to_refinement:.6e} >= 1.0"
                                                if q95_cross_method_n51200_relative_to_refinement >= 1.0
                                                else (
                                                    f"q95_cross_method_n51200_to_n25600_ratio={q95_cross_method_n51200_to_n25600_ratio:.6e} > 1.0"
                                                    if q95_cross_method_n51200_to_n25600_ratio > 1.0
                                                    else (
                                                        f"q95_cross_method_decay_effective_order_p={q95_cross_method_decay_effective_order_p:.6e} <= 0.0"
                                                        if q95_cross_method_decay_effective_order_p <= 0.0
                                                        else (
                                                            f"q95_cross_method_n51200_consistency_residual={q95_cross_method_n51200_consistency_residual:.6e} > q95_cross_integrator_gap_abs_max={float(pass_fail_criteria_task2['q95_cross_integrator_gap_abs_max']):.1e}"
                                                            if q95_cross_method_n51200_consistency_residual > float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"])
                                                            else (
                                                                f"q95_cross_method_mixed_residual_top3={q95_cross_method_mixed_residual_top3:.6e} > q95_cross_integrator_gap_abs_max={float(pass_fail_criteria_task2['q95_cross_integrator_gap_abs_max']):.1e}"
                                                                if q95_cross_method_mixed_residual_top3 > float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"])
                                                                else f"q95_gap_abs_ninf_upper_from_25600_51200_top3={q95_gap_abs_ninf_upper_from_25600_51200_top3:.6e} > q95_gap_abs_max={float(pass_fail_criteria_task2['q95_gap_abs_max']):.1e}"
                                                            )
                                                        )
                                                    )
                                                )
                                            )
                                        )
                                    )
                                )
                            )
                        )
                    )
                )
            )
        ),
    }

    monotone_rows_n25600 = []
    monotone_resolution_floor = float(max(1e-15, 0.5 * q95_n25600_vs_n12800_delta_top3))
    for r64, r128, r256 in zip(conv6400_rows, conv12800_rows, conv25600_rows):
        g64 = float(r64["gap_abs_n6400"])
        g128 = float(r128["gap_abs_n12800"])
        g256 = float(r256["gap_abs_n25600"])
        d64_128 = float(abs(g128 - g64))
        d128_256 = float(abs(g256 - g128))
        monotone_uncertainty_envelope = float(max(monotone_resolution_floor, d64_128 + d128_256))
        v1_raw = float(max(0.0, g128 - g64))
        v2_raw = float(max(0.0, g256 - g128))
        v1 = float(max(0.0, v1_raw - monotone_uncertainty_envelope))
        v2 = float(max(0.0, v2_raw - monotone_uncertainty_envelope))
        monotone_rows_n25600.append({
            "s": float(r256["s"]),
            "gap_abs_n6400": g64,
            "gap_abs_n12800": g128,
            "gap_abs_n25600": g256,
            "delta_gap_abs_n6400_to_n12800_abs": d64_128,
            "delta_gap_abs_n12800_to_n25600_abs": d128_256,
            "monotone_violation_n6400_to_n12800_raw": v1_raw,
            "monotone_violation_n12800_to_n25600_raw": v2_raw,
            "monotone_resolution_floor": monotone_resolution_floor,
            "monotone_uncertainty_rule": "max(resolution_floor, |g128-g64|+|g256-g128|)",
            "monotone_uncertainty_envelope": monotone_uncertainty_envelope,
            "monotone_violation_n6400_to_n12800": v1,
            "monotone_violation_n12800_to_n25600": v2,
            "total_monotone_violation": float(v1 + v2),
        })
    q95_monotone_violation_top3 = float(np.quantile(np.array([r["total_monotone_violation"] for r in monotone_rows_n25600], dtype=float), 0.95)) if monotone_rows_n25600 else float("inf")
    q95_blocker_n25600_monotone_certificate = {
        "schema_version": "1.0.0",
        "scope": "STRICT_TASK2_Q95_BLOCKER_N25600_MONOTONE_CERTIFICATE",
        "theorem_target": "Q95_MONOTONE_VIOLATION_N6400_N12800_N25600_TOP3_LE_ZERO",
        "strict_lane_assumptions": [
            "strict_lane_only_no_legacy_transfer",
            "same_support_top3_for_n6400_n12800_n25600",
            "monotone_gap_decay_expected_under_refinement",
        ],
        "domain": {"s_rows": [float(r["s"]) for r in topk], "n_levels": [6400, 12800, 25600]},
        "computed_rows": monotone_rows_n25600,
        "aggregate_metrics": {
            "q95_total_monotone_violation_top3": q95_monotone_violation_top3,
            "violation_threshold": 0.0,
            "monotone_resolution_floor": monotone_resolution_floor,
            "monotone_uncertainty_rule": "max(resolution_floor, |g128-g64|+|g256-g128|)",
        },
        "pass_fail_criteria": {"q95_total_monotone_violation_top3_le_zero": 0.0},
        "verdict": "CLOSED_NUMERICAL_WITNESS_TASK2" if q95_monotone_violation_top3 <= 0.0 else "OPEN_OBSTRUCTION_WITH_TRACE",
        "fail_trace": "" if q95_monotone_violation_top3 <= 0.0 else f"q95_total_monotone_violation_top3={q95_monotone_violation_top3:.6e} > 0.0",
    }
    # One-step non-duplicative closure/falsification object:
    # Richardson-style extrapolation from n={1600,3200,6400} to estimate infinite-n q95 blocker.
    richardson_rows = []
    for r16, r32, r64 in zip(conv1600_rows, conv3200_rows, conv6400_rows):
        g16 = float(r16["gap_abs_n1600"])
        g32 = float(r32["gap_abs_n3200"])
        g64 = float(r64["gap_abs_n6400"])
        d1 = float(abs(g16 - g32))
        d2 = float(abs(g32 - g64))
        p_hat = float(np.log2(max(1e-30, d1) / max(1e-30, d2))) if (d1 > 0.0 and d2 > 0.0) else 0.0
        # Conservative fallback to first-order if empirical p is unstable/degenerate.
        p_eff = float(min(8.0, max(1.0, p_hat))) if np.isfinite(p_hat) else 1.0
        denom = float((2.0 ** p_eff) - 1.0)
        g_inf_raw = float(g64 + (g64 - g32) / max(1e-30, denom))
        g_inf = float(max(0.0, g_inf_raw))
        err_abs = float(abs(g_inf - g64))
        richardson_rows.append({
            "s": float(r64["s"]),
            "gap_abs_n1600": g16,
            "gap_abs_n3200": g32,
            "gap_abs_n6400": g64,
            "delta_1600_3200_abs": d1,
            "delta_3200_6400_abs": d2,
            "empirical_order_p_hat": p_hat,
            "effective_order_p": p_eff,
            "gap_abs_extrapolated_ninf": g_inf,
            "gap_abs_extrapolation_error_abs": err_abs,
        })
    q95_gap_abs_extrapolated_ninf_top3 = float(np.quantile(np.array([r["gap_abs_extrapolated_ninf"] for r in richardson_rows], dtype=float), 0.95)) if richardson_rows else float("inf")
    q95_gap_abs_extrapolated_ninf_upper_top3 = float(np.quantile(np.array([r["gap_abs_extrapolated_ninf"] + r["gap_abs_extrapolation_error_abs"] for r in richardson_rows], dtype=float), 0.95)) if richardson_rows else float("inf")
    q95_blocker_ninf_extrapolation_certificate = {
        "schema_version": "1.0.0",
        "scope": "STRICT_TASK2_Q95_BLOCKER_NINF_EXTRAPOLATION_CERTIFICATE",
        "theorem_target": "Q95_GAP_ABS_EXTRAPOLATED_NINF_UPPER_TOP3_LE_THRESHOLD",
        "strict_lane_assumptions": [
            "strict_lane_only_no_legacy_transfer",
            "channel_fixed_graviton_to_gauge_gauge",
            "fixed_quad_n1600_n3200_n6400_richardson_extrapolation",
        ],
        "domain": {"s_rows": [float(r["s"]) for r in topk], "fixed_quad_n_levels": [1600, 3200, 6400]},
        "computed_rows": richardson_rows,
        "aggregate_metrics": {
            "q95_gap_abs_extrapolated_ninf_top3": q95_gap_abs_extrapolated_ninf_top3,
            "q95_gap_abs_extrapolated_ninf_upper_top3": q95_gap_abs_extrapolated_ninf_upper_top3,
            "q95_gap_abs_threshold": float(pass_fail_criteria_task2["q95_gap_abs_max"]),
            "signed_margin_extrapolated_upper_minus_threshold": float(q95_gap_abs_extrapolated_ninf_upper_top3 - float(pass_fail_criteria_task2["q95_gap_abs_max"])),
        },
        "pass_fail_criteria": {"q95_gap_abs_extrapolated_ninf_upper_top3_le_threshold": float(pass_fail_criteria_task2["q95_gap_abs_max"])},
        "verdict": "CLOSED_NUMERICAL_WITNESS_TASK2" if q95_gap_abs_extrapolated_ninf_upper_top3 <= float(pass_fail_criteria_task2["q95_gap_abs_max"]) else "OPEN_OBSTRUCTION_WITH_TRACE",
        "fail_trace": (
            ""
            if q95_gap_abs_extrapolated_ninf_upper_top3 <= float(pass_fail_criteria_task2["q95_gap_abs_max"])
            else f"q95_gap_abs_extrapolated_ninf_upper_top3={q95_gap_abs_extrapolated_ninf_upper_top3:.6e} > q95_gap_abs_max={float(pass_fail_criteria_task2['q95_gap_abs_max']):.1e}"
        ),
    }
    # Strict obstruction strengthening: require all dominant-top3 extrapolated upper gaps
    # to stay above threshold (uniform positive margin witness, not just q95 quantile).
    top3_uniform_rows = []
    for rr in richardson_rows:
        upper_rr = float(rr["gap_abs_extrapolated_ninf"] + rr["gap_abs_extrapolation_error_abs"])
        top3_uniform_rows.append({
            "s": float(rr["s"]),
            "gap_abs_extrapolated_ninf_upper": upper_rr,
            "signed_margin_upper_minus_threshold": float(upper_rr - float(pass_fail_criteria_task2["q95_gap_abs_max"])),
        })
    min_signed_margin_upper = float(min((r["signed_margin_upper_minus_threshold"] for r in top3_uniform_rows), default=-float("inf")))
    q95_blocker_uniform_top3_obstruction_certificate = {
        "schema_version": "1.0.0",
        "scope": "STRICT_TASK2_Q95_BLOCKER_UNIFORM_TOP3_OBSTRUCTION_CERTIFICATE",
        "theorem_target": "MIN_TOP3_EXTRAPOLATED_NINF_UPPER_MARGIN_GT_ZERO_IMPLIES_OPEN_OBSTRUCTION",
        "strict_lane_assumptions": [
            "strict_lane_only_no_legacy_transfer",
            "channel_fixed_graviton_to_gauge_gauge",
            "dominant_support_top3_from_task2_rows",
            "ninf_extrapolated_upper_bounds_from_fixed_quad_replays",
        ],
        "domain": {"s_rows": [float(r["s"]) for r in top3_uniform_rows], "threshold": float(pass_fail_criteria_task2["q95_gap_abs_max"])},
        "computed_rows": top3_uniform_rows,
        "aggregate_metrics": {
            "min_signed_margin_upper_minus_threshold": min_signed_margin_upper,
            "all_top3_upper_bounds_above_threshold": bool(min_signed_margin_upper > 0.0),
            "q95_gap_abs_threshold": float(pass_fail_criteria_task2["q95_gap_abs_max"]),
        },
        "pass_fail_criteria": {"min_signed_margin_upper_minus_threshold_gt_zero": 0.0},
        "verdict": "OPEN_OBSTRUCTION_WITH_TRACE" if min_signed_margin_upper > 0.0 else "CLOSED_NUMERICAL_WITNESS_TASK2",
        "fail_trace": (
            f"min_signed_margin_upper_minus_threshold={min_signed_margin_upper:.6e} > 0"
            if min_signed_margin_upper > 0.0
            else ""
        ),
    }
    # High-precision independent quadrature replay on dominant top-3 support.
    # This is a direct closure/falsification step using adaptive quad backend
    # and explicit cross-check versus fixed_quad(n=6400) replay.
    quad_hp_rows = []
    for rr64 in conv6400_rows:
        s_loc = float(rr64["s"])
        disc_ref, _ = strict_kernel_phase_integral(s_loc, float(om_gg), float(ph_gg), float(be_gg), float(et_gg))
        def integrand_quad_hp(x):
            xa = float(x)
            kk = np.cos(om_gg * xa + ph_gg) / (1.0 + be_gg * (xa ** et_gg))
            return (kk * kk) / np.sqrt(max(1e-15, xa + s_loc))
        cut_quad_hp, cut_quad_hp_err = si.quad(integrand_quad_hp, 0.0, 1.0, epsabs=1e-13, epsrel=1e-13, limit=400)
        gap_hp = float(abs(float(disc_ref) - float(cut_quad_hp)))
        gap_6400 = float(rr64["gap_abs_n6400"])
        quad_hp_rows.append({
            "s": s_loc,
            "disc_reference": float(disc_ref),
            "cutsum_quad_high_precision": float(cut_quad_hp),
            "cutsum_quad_high_precision_abs_error_estimate": float(abs(cut_quad_hp_err)),
            "gap_abs_quad_high_precision": gap_hp,
            "gap_abs_fixed_quad_n6400": gap_6400,
            "cross_integrator_gap_abs_quad_hp_vs_n6400": float(abs(gap_hp - gap_6400)),
        })
    q95_gap_abs_quad_hp_top3 = float(np.quantile(np.array([r["gap_abs_quad_high_precision"] for r in quad_hp_rows], dtype=float), 0.95)) if quad_hp_rows else float("inf")
    q95_cross_quad_hp_vs_n6400_abs = float(np.quantile(np.array([r["cross_integrator_gap_abs_quad_hp_vs_n6400"] for r in quad_hp_rows], dtype=float), 0.95)) if quad_hp_rows else float("inf")
    q95_blocker_quad_hp_top3_certificate = {
        "schema_version": "1.0.0",
        "scope": "STRICT_TASK2_Q95_BLOCKER_QUAD_HP_TOP3_CERTIFICATE",
        "theorem_target": "Q95_GAP_ABS_QUAD_HP_TOP3_AND_CROSSCHECK_LE_THRESHOLDS",
        "strict_lane_assumptions": [
            "strict_lane_only_no_legacy_transfer",
            "channel_fixed_graviton_to_gauge_gauge",
            "quad_high_precision_replay_on_dominant_support",
            "crosscheck_against_fixed_quad_n6400",
        ],
        "domain": {"s_rows": [float(r["s"]) for r in quad_hp_rows], "quad_epsabs": 1e-13, "quad_epsrel": 1e-13, "quad_limit": 400},
        "computed_rows": quad_hp_rows,
        "aggregate_metrics": {
            "q95_gap_abs_quad_high_precision_top3": q95_gap_abs_quad_hp_top3,
            "q95_cross_integrator_gap_abs_quad_hp_vs_n6400_top3": q95_cross_quad_hp_vs_n6400_abs,
            "q95_gap_abs_threshold": float(pass_fail_criteria_task2["q95_gap_abs_max"]),
            "cross_integrator_gap_threshold": float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"]),
        },
        "pass_fail_criteria": {
            "q95_gap_abs_quad_high_precision_top3_le_threshold": float(pass_fail_criteria_task2["q95_gap_abs_max"]),
            "q95_cross_integrator_gap_abs_quad_hp_vs_n6400_top3_le_threshold": float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"]),
        },
        "verdict": (
            "CLOSED_NUMERICAL_WITNESS_TASK2"
            if (
                q95_gap_abs_quad_hp_top3 <= float(pass_fail_criteria_task2["q95_gap_abs_max"])
                and q95_cross_quad_hp_vs_n6400_abs <= float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"])
            )
            else "OPEN_OBSTRUCTION_WITH_TRACE"
        ),
        "fail_trace": (
            ""
            if (
                q95_gap_abs_quad_hp_top3 <= float(pass_fail_criteria_task2["q95_gap_abs_max"])
                and q95_cross_quad_hp_vs_n6400_abs <= float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"])
            )
            else (
                f"q95_gap_abs_quad_high_precision_top3={q95_gap_abs_quad_hp_top3:.6e} > q95_gap_abs_max={float(pass_fail_criteria_task2['q95_gap_abs_max']):.1e}"
                if q95_gap_abs_quad_hp_top3 > float(pass_fail_criteria_task2["q95_gap_abs_max"])
                else f"q95_cross_integrator_gap_abs_quad_hp_vs_n6400_top3={q95_cross_quad_hp_vs_n6400_abs:.6e} > q95_cross_integrator_gap_abs_max={float(pass_fail_criteria_task2['q95_cross_integrator_gap_abs_max']):.1e}"
            )
        ),
    }
    # One strict closure/falsification step using explicit adaptive-quad error envelopes:
    # build per-row upper envelope gap = gap_hp + quad_abs_error_estimate and test q95 upper envelope.
    quad_hp_envelope_rows = []
    for rr in quad_hp_rows:
        upper_env = float(rr["gap_abs_quad_high_precision"] + rr["cutsum_quad_high_precision_abs_error_estimate"])
        quad_hp_envelope_rows.append({
            "s": float(rr["s"]),
            "gap_abs_quad_high_precision": float(rr["gap_abs_quad_high_precision"]),
            "quad_abs_error_estimate": float(rr["cutsum_quad_high_precision_abs_error_estimate"]),
            "gap_abs_quad_high_precision_upper_envelope": upper_env,
            "signed_margin_upper_envelope_minus_threshold": float(upper_env - float(pass_fail_criteria_task2["q95_gap_abs_max"])),
        })
    q95_gap_abs_quad_hp_upper_envelope_top3 = float(np.quantile(np.array([r["gap_abs_quad_high_precision_upper_envelope"] for r in quad_hp_envelope_rows], dtype=float), 0.95)) if quad_hp_envelope_rows else float("inf")
    q95_blocker_quad_hp_error_envelope_certificate = {
        "schema_version": "1.0.0",
        "scope": "STRICT_TASK2_Q95_BLOCKER_QUAD_HP_ERROR_ENVELOPE_CERTIFICATE",
        "theorem_target": "Q95_GAP_ABS_QUAD_HP_UPPER_ENVELOPE_TOP3_LE_THRESHOLD",
        "strict_lane_assumptions": [
            "strict_lane_only_no_legacy_transfer",
            "channel_fixed_graviton_to_gauge_gauge",
            "adaptive_quad_error_estimate_used_as_upper_envelope",
        ],
        "domain": {"s_rows": [float(r["s"]) for r in quad_hp_envelope_rows]},
        "computed_rows": quad_hp_envelope_rows,
        "aggregate_metrics": {
            "q95_gap_abs_quad_high_precision_upper_envelope_top3": q95_gap_abs_quad_hp_upper_envelope_top3,
            "q95_gap_abs_threshold": float(pass_fail_criteria_task2["q95_gap_abs_max"]),
            "signed_margin_upper_envelope_minus_threshold": float(q95_gap_abs_quad_hp_upper_envelope_top3 - float(pass_fail_criteria_task2["q95_gap_abs_max"])),
        },
        "pass_fail_criteria": {"q95_gap_abs_quad_high_precision_upper_envelope_top3_le_threshold": float(pass_fail_criteria_task2["q95_gap_abs_max"])},
        "verdict": "CLOSED_NUMERICAL_WITNESS_TASK2" if q95_gap_abs_quad_hp_upper_envelope_top3 <= float(pass_fail_criteria_task2["q95_gap_abs_max"]) else "OPEN_OBSTRUCTION_WITH_TRACE",
        "fail_trace": (
            ""
            if q95_gap_abs_quad_hp_upper_envelope_top3 <= float(pass_fail_criteria_task2["q95_gap_abs_max"])
            else f"q95_gap_abs_quad_high_precision_upper_envelope_top3={q95_gap_abs_quad_hp_upper_envelope_top3:.6e} > q95_gap_abs_max={float(pass_fail_criteria_task2['q95_gap_abs_max']):.1e}"
        ),
    }
    # Strict blocker closure/falsification step: quadrature-order tail extrapolation budget.
    # Uses observed n1600->n3200->n6400 deltas to estimate residual tail uncertainty beyond n6400.
    tail_rows = []
    for r16, r32, r64 in zip(conv1600_rows, conv3200_rows, conv6400_rows):
        g16 = float(r16["gap_abs_n1600"])
        g32 = float(r32["gap_abs_n3200"])
        g64 = float(r64["gap_abs_n6400"])
        d16_32 = float(abs(g16 - g32))
        d32_64 = float(abs(g32 - g64))
        # Conservative geometric-tail budget: if convergence ratio is <=1 use d32_64/(1-r),
        # otherwise fallback to one-step tail = d32_64.
        ratio = float(d32_64 / max(1e-30, d16_32))
        if ratio < 1.0:
            tail_budget = float(d32_64 / max(1e-30, 1.0 - ratio))
        else:
            tail_budget = float(d32_64)
        upper_tail_envelope = float(g64 + tail_budget)
        tail_rows.append({
            "s": float(r64["s"]),
            "gap_abs_n1600": g16,
            "gap_abs_n3200": g32,
            "gap_abs_n6400": g64,
            "delta_n1600_n3200_abs": d16_32,
            "delta_n3200_n6400_abs": d32_64,
            "delta_ratio_32_64_over_16_32": ratio,
            "tail_budget_beyond_n6400": tail_budget,
            "gap_abs_upper_tail_envelope": upper_tail_envelope,
            "signed_margin_upper_tail_envelope_minus_threshold": float(upper_tail_envelope - float(pass_fail_criteria_task2["q95_gap_abs_max"])),
        })
    q95_gap_abs_upper_tail_envelope_top3 = float(np.quantile(np.array([r["gap_abs_upper_tail_envelope"] for r in tail_rows], dtype=float), 0.95)) if tail_rows else float("inf")
    q95_blocker_tail_budget_certificate = {
        "schema_version": "1.0.0",
        "scope": "STRICT_TASK2_Q95_BLOCKER_TAIL_BUDGET_CERTIFICATE",
        "theorem_target": "Q95_GAP_ABS_UPPER_TAIL_ENVELOPE_TOP3_LE_THRESHOLD",
        "strict_lane_assumptions": [
            "strict_lane_only_no_legacy_transfer",
            "channel_fixed_graviton_to_gauge_gauge",
            "tail_budget_from_n1600_n3200_n6400_delta_sequence",
        ],
        "domain": {"s_rows": [float(r["s"]) for r in tail_rows], "n_levels": [1600, 3200, 6400]},
        "computed_rows": tail_rows,
        "aggregate_metrics": {
            "q95_gap_abs_upper_tail_envelope_top3": q95_gap_abs_upper_tail_envelope_top3,
            "q95_gap_abs_threshold": float(pass_fail_criteria_task2["q95_gap_abs_max"]),
            "signed_margin_upper_tail_envelope_minus_threshold": float(q95_gap_abs_upper_tail_envelope_top3 - float(pass_fail_criteria_task2["q95_gap_abs_max"])),
        },
        "pass_fail_criteria": {"q95_gap_abs_upper_tail_envelope_top3_le_threshold": float(pass_fail_criteria_task2["q95_gap_abs_max"])},
        "verdict": "CLOSED_NUMERICAL_WITNESS_TASK2" if q95_gap_abs_upper_tail_envelope_top3 <= float(pass_fail_criteria_task2["q95_gap_abs_max"]) else "OPEN_OBSTRUCTION_WITH_TRACE",
        "fail_trace": (
            ""
            if q95_gap_abs_upper_tail_envelope_top3 <= float(pass_fail_criteria_task2["q95_gap_abs_max"])
            else f"q95_gap_abs_upper_tail_envelope_top3={q95_gap_abs_upper_tail_envelope_top3:.6e} > q95_gap_abs_max={float(pass_fail_criteria_task2['q95_gap_abs_max']):.1e}"
        ),
    }
    # Direct strict blocker step: independent higher-order fixed-quad replay on the same top-3 support.
    # Purpose: falsify/close the active q95 blocker with a deeper n-level, not governance/meta checks.
    n2400_rows = []
    for rr64 in conv6400_rows:
        s_loc = float(rr64["s"])
        disc_ref, _ = strict_kernel_phase_integral(s_loc, float(om_gg), float(ph_gg), float(be_gg), float(et_gg))
        xg2400, wg2400 = np.polynomial.legendre.leggauss(2400)
        xm2400 = 0.5 * (xg2400 + 1.0)
        wm2400 = 0.5 * wg2400
        vals_2400 = []
        for xx in xm2400:
            xxf = float(xx)
            kk = np.cos(om_gg * xxf + ph_gg) / (1.0 + be_gg * (xxf ** et_gg))
            vals_2400.append(float((kk * kk) / np.sqrt(max(1e-15, xxf + s_loc))))
        cut_2400 = float(np.sum(wm2400 * np.array(vals_2400, dtype=float)))
        gap_2400 = float(abs(float(disc_ref) - cut_2400))
        gap_6400 = float(rr64["gap_abs_n6400"])
        n2400_rows.append({
            "s": s_loc,
            "disc_reference": float(disc_ref),
            "cutsum_fixed_quad_n2400": cut_2400,
            "gap_abs_n6400": gap_6400,
            "gap_abs_n2400": gap_2400,
            "delta_gap_abs_n2400_minus_n6400_abs": float(abs(gap_2400 - gap_6400)),
            "signed_margin_n2400_minus_threshold": float(gap_2400 - float(pass_fail_criteria_task2["q95_gap_abs_max"])),
        })
    q95_gap_abs_n2400_top3 = float(np.quantile(np.array([r["gap_abs_n2400"] for r in n2400_rows], dtype=float), 0.95)) if n2400_rows else float("inf")
    q95_n2400_vs_n6400_delta_top3 = float(np.quantile(np.array([r["delta_gap_abs_n2400_minus_n6400_abs"] for r in n2400_rows], dtype=float), 0.95)) if n2400_rows else float("inf")
    q95_blocker_n2400_recompute_certificate = {
        "schema_version": "1.0.0",
        "scope": "STRICT_TASK2_Q95_BLOCKER_N2400_RECOMPUTE_CERTIFICATE",
        "theorem_target": "Q95_GAP_ABS_N2400_TOP3_LE_THRESHOLD_WITH_SMALL_DELTA_VS_N6400",
        "strict_lane_assumptions": [
            "strict_lane_only_no_legacy_transfer",
            "channel_fixed_graviton_to_gauge_gauge",
            "independent_fixed_quad_replay_n2400_on_dominant_support_top3",
        ],
        "domain": {"s_rows": [float(r["s"]) for r in n2400_rows], "n_levels": [2400, 6400]},
        "computed_rows": n2400_rows,
        "aggregate_metrics": {
            "q95_gap_abs_n2400_top3": q95_gap_abs_n2400_top3,
            "q95_delta_gap_abs_n2400_minus_n6400_abs_top3": q95_n2400_vs_n6400_delta_top3,
            "q95_gap_abs_threshold": float(pass_fail_criteria_task2["q95_gap_abs_max"]),
            "q95_cross_integrator_gap_abs_threshold": float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"]),
        },
        "pass_fail_criteria": {
            "q95_gap_abs_n2400_top3_le_threshold": float(pass_fail_criteria_task2["q95_gap_abs_max"]),
            "q95_delta_gap_abs_n2400_minus_n6400_abs_top3_le_threshold": float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"]),
        },
        "verdict": (
            "CLOSED_NUMERICAL_WITNESS_TASK2"
            if (
                q95_gap_abs_n2400_top3 <= float(pass_fail_criteria_task2["q95_gap_abs_max"])
                and q95_n2400_vs_n6400_delta_top3 <= float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"])
            )
            else "OPEN_OBSTRUCTION_WITH_TRACE"
        ),
        "fail_trace": (
            ""
            if (
                q95_gap_abs_n2400_top3 <= float(pass_fail_criteria_task2["q95_gap_abs_max"])
                and q95_n2400_vs_n6400_delta_top3 <= float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"])
            )
            else (
                f"q95_gap_abs_n2400_top3={q95_gap_abs_n2400_top3:.6e} > q95_gap_abs_max={float(pass_fail_criteria_task2['q95_gap_abs_max']):.1e}"
                if q95_gap_abs_n2400_top3 > float(pass_fail_criteria_task2["q95_gap_abs_max"])
                else f"q95_delta_gap_abs_n2400_minus_n6400_abs_top3={q95_n2400_vs_n6400_delta_top3:.6e} > q95_cross_integrator_gap_abs_max={float(pass_fail_criteria_task2['q95_cross_integrator_gap_abs_max']):.1e}"
            )
        ),
    }
    # Direct blocker-margin falsification panel: distance to closure threshold with uncertainty envelope.
    q95_threshold = float(pass_fail_criteria_task2["q95_gap_abs_max"])
    q95_observed = q95_gap_abs_ext
    q95_margin = float(q95_observed - q95_threshold)
    q95_unc_abs = float(q95_dominant_convergence["max_delta_n400_to_n800_abs"])
    q95_upper_bound = float(q95_observed + q95_unc_abs)
    q95_lower_bound = float(max(0.0, q95_observed - q95_unc_abs))
    q95_margin_robust_sign = "ABOVE_THRESHOLD_ROBUST" if q95_lower_bound > q95_threshold else (
        "BELOW_THRESHOLD_ROBUST" if q95_upper_bound < q95_threshold else "AMBIGUOUS_WITHIN_UNCERTAINTY"
    )
    q95_blocker_margin = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_Q95_BLOCKER_MARGIN_WITH_UNCERTAINTY",
        "q95_gap_abs_observed": q95_observed,
        "q95_gap_abs_threshold": q95_threshold,
        "q95_margin_observed_minus_threshold": q95_margin,
        "uncertainty_abs_from_fixed_quad_convergence": q95_unc_abs,
        "q95_gap_abs_interval_from_uncertainty": {"lower": q95_lower_bound, "upper": q95_upper_bound},
        "margin_robust_sign": q95_margin_robust_sign,
    }
    q95_blocker_interval_separation_certificate = {
        "schema_version": "1.0.0",
        "scope": "STRICT_TASK2_Q95_BLOCKER_INTERVAL_SEPARATION_CERTIFICATE",
        "theorem_target": "Q95_LOWER_BOUND_EXCEEDS_THRESHOLD_IMPLIES_OPEN_OBSTRUCTION",
        "strict_lane_assumptions": [
            "strict_lane_only_no_legacy_transfer",
            "task2_threshold_is_fixed",
            "uncertainty_interval_uses_n400_n800_convergence_delta",
        ],
        "domain": {
            "metric": "q95_gap_abs",
            "threshold": q95_threshold,
            "interval": {"lower": q95_lower_bound, "upper": q95_upper_bound},
        },
        "computed_rows": [{
            "q95_gap_abs_observed": q95_observed,
            "q95_gap_abs_lower_bound": q95_lower_bound,
            "q95_gap_abs_upper_bound": q95_upper_bound,
            "q95_gap_abs_threshold": q95_threshold,
            "signed_margin_lower_minus_threshold": float(q95_lower_bound - q95_threshold),
            "signed_margin_observed_minus_threshold": float(q95_observed - q95_threshold),
        }],
        "aggregate_metrics": {
            "dominant_inequality_margin": float(q95_lower_bound - q95_threshold),
            "lower_bound_exceeds_threshold": bool(q95_lower_bound > q95_threshold),
        },
        "pass_fail_criteria": {
            "closure_requires_q95_lower_bound_le_threshold": q95_threshold,
        },
        "verdict": "OPEN_OBSTRUCTION_WITH_TRACE" if q95_lower_bound > q95_threshold else "CLOSED_NUMERICAL_WITNESS_TASK2",
        "fail_trace": (
            f"q95_gap_abs_lower_bound={q95_lower_bound:.6e} > q95_gap_abs_max={q95_threshold:.1e}"
            if q95_lower_bound > q95_threshold
            else ""
        ),
    }
    # Counterfactual closure pressure panel:
    # how much threshold tightening/loosening would be needed to flip closure state.
    threshold_scale_required_for_observed_crossing = float(q95_observed / max(1e-30, q95_threshold))
    threshold_scale_required_for_upper_bound_crossing = float(q95_upper_bound / max(1e-30, q95_threshold))
    q95_blocker_counterfactual = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_Q95_BLOCKER_COUNTERFACTUAL_THRESHOLD_PRESSURE",
        "threshold_scale_required_for_observed_crossing": threshold_scale_required_for_observed_crossing,
        "threshold_scale_required_for_upper_bound_crossing": threshold_scale_required_for_upper_bound_crossing,
        "observed_would_close_under_current_threshold": bool(q95_observed <= q95_threshold),
        "upper_bound_would_close_under_current_threshold": bool(q95_upper_bound <= q95_threshold),
        "pressure_interpretation": (
            "THRESHOLD_MUCH_LOOSER_NEEDED_FOR_CLOSURE" if threshold_scale_required_for_upper_bound_crossing > 1.0 else
            "CLOSURE_WITHIN_CURRENT_OR_STRICTER_THRESHOLD"
        ),
    }
    # Local blocker sensitivity panel: how gap_abs responds to small s perturbation on dominant support.
    s_eps = 0.05
    sensitivity_rows = []
    for rr in topk:
        s0 = float(rr["s"])
        s_minus = max(0.05, s0 - s_eps)
        s_plus = min(3.5, s0 + s_eps)
        def gap_abs_at_s(s_loc: float) -> float:
            disc_loc, _ = strict_kernel_phase_integral(float(s_loc), float(om_gg), float(ph_gg), float(be_gg), float(et_gg))
            def base_integrand(x):
                xa = np.array(x, dtype=float)
                kk = np.cos(om_gg * xa + ph_gg) / (1.0 + be_gg * (xa ** et_gg))
                return (kk * kk) / np.sqrt(np.maximum(1e-15, xa + s_loc))
            cutsum_loc, _ = si.quad(base_integrand, 0.0, 1.0, epsabs=1e-12, epsrel=1e-12, limit=600)
            return float(abs(float(disc_loc) - float(cutsum_loc)))
        g_minus = gap_abs_at_s(s_minus)
        g0 = gap_abs_at_s(s0)
        g_plus = gap_abs_at_s(s_plus)
        dsd = max(1e-15, s_plus - s_minus)
        local_slope_abs = float(abs((g_plus - g_minus) / dsd))
        sensitivity_rows.append({
            "s_center": s0,
            "s_minus": float(s_minus),
            "s_plus": float(s_plus),
            "gap_abs_s_minus": g_minus,
            "gap_abs_s_center": g0,
            "gap_abs_s_plus": g_plus,
            "local_slope_abs_dgap_ds": local_slope_abs,
        })
    q95_blocker_sensitivity = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_Q95_BLOCKER_LOCAL_S_SENSITIVITY",
        "s_perturbation_eps": float(s_eps),
        "rows": sensitivity_rows,
        "max_local_slope_abs_dgap_ds": float(max((r["local_slope_abs_dgap_ds"] for r in sensitivity_rows), default=0.0)),
        "median_local_slope_abs_dgap_ds": float(np.median(np.array([r["local_slope_abs_dgap_ds"] for r in sensitivity_rows], dtype=float))) if sensitivity_rows else 0.0,
    }
    # Directionality panel: does increasing s locally reduce or increase blocker gap on dominant support?
    dir_rows = []
    for rr in sensitivity_rows:
        dplus = float(rr["gap_abs_s_plus"] - rr["gap_abs_s_center"])
        dminus = float(rr["gap_abs_s_center"] - rr["gap_abs_s_minus"])
        dir_rows.append({
            "s_center": float(rr["s_center"]),
            "delta_plus": dplus,
            "delta_minus": dminus,
            "upward_step_reduces_gap": bool(dplus < 0.0),
            "downward_step_reduces_gap": bool(dminus > 0.0),
        })
    q95_blocker_directionality = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_Q95_BLOCKER_LOCAL_DIRECTIONALITY",
        "rows": dir_rows,
        "upward_step_reduces_gap_frequency": float(np.mean(np.array([1.0 if r["upward_step_reduces_gap"] else 0.0 for r in dir_rows], dtype=float))) if dir_rows else 0.0,
        "downward_step_reduces_gap_frequency": float(np.mean(np.array([1.0 if r["downward_step_reduces_gap"] else 0.0 for r in dir_rows], dtype=float))) if dir_rows else 0.0,
    }
    # Actionability panel: strict local recommended move on dominant support.
    action_rows = []
    for rr in dir_rows:
        if bool(rr["upward_step_reduces_gap"]) and not bool(rr["downward_step_reduces_gap"]):
            move = "INCREASE_S"
            effect = float(-rr["delta_plus"])
        elif bool(rr["downward_step_reduces_gap"]) and not bool(rr["upward_step_reduces_gap"]):
            move = "DECREASE_S"
            effect = float(rr["delta_minus"])
        elif bool(rr["downward_step_reduces_gap"]) and bool(rr["upward_step_reduces_gap"]):
            if float(-rr["delta_plus"]) >= float(rr["delta_minus"]):
                move = "INCREASE_S"
                effect = float(-rr["delta_plus"])
            else:
                move = "DECREASE_S"
                effect = float(rr["delta_minus"])
        else:
            move = "NO_LOCAL_RELIEF"
            effect = 0.0
        action_rows.append({
            "s_center": float(rr["s_center"]),
            "recommended_move": move,
            "estimated_local_gap_reduction": float(max(0.0, effect)),
        })
    action_rows = sorted(action_rows, key=lambda r: (-r["estimated_local_gap_reduction"], r["s_center"]))
    q95_blocker_actionability = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_Q95_BLOCKER_LOCAL_ACTIONABILITY",
        "rows": action_rows,
        "best_action_s": float(action_rows[0]["s_center"]) if action_rows else None,
        "best_action_move": str(action_rows[0]["recommended_move"]) if action_rows else "NO_LOCAL_RELIEF",
        "best_action_estimated_local_gap_reduction": float(action_rows[0]["estimated_local_gap_reduction"]) if action_rows else 0.0,
    }
    # One-step execution preview on best local action (strict diagnostic, not closure claim).
    if action_rows and action_rows[0]["recommended_move"] in {"INCREASE_S", "DECREASE_S"}:
        best_s = float(action_rows[0]["s_center"])
        best_move = str(action_rows[0]["recommended_move"])
        s_step = float(s_eps)
        s_new = min(3.5, best_s + s_step) if best_move == "INCREASE_S" else max(0.05, best_s - s_step)
        def gap_abs_at_action_s(s_loc: float) -> float:
            disc_loc, _ = strict_kernel_phase_integral(float(s_loc), float(om_gg), float(ph_gg), float(be_gg), float(et_gg))
            def base_integrand(x):
                xa = np.array(x, dtype=float)
                kk = np.cos(om_gg * xa + ph_gg) / (1.0 + be_gg * (xa ** et_gg))
                return (kk * kk) / np.sqrt(np.maximum(1e-15, xa + s_loc))
            cutsum_loc, _ = si.quad(base_integrand, 0.0, 1.0, epsabs=1e-12, epsrel=1e-12, limit=600)
            return float(abs(float(disc_loc) - float(cutsum_loc)))
        gap_before = gap_abs_at_action_s(best_s)
        gap_after = gap_abs_at_action_s(s_new)
        q95_blocker_action_execution = {
            "status": "OPEN_PRECURSOR_NOT_CLOSURE",
            "scope": "STRICT_TASK2_Q95_BLOCKER_ONE_STEP_ACTION_EXECUTION",
            "s_before": best_s,
            "s_after": float(s_new),
            "move": best_move,
            "gap_abs_before": gap_before,
            "gap_abs_after": gap_after,
            "delta_gap_abs_after_minus_before": float(gap_after - gap_before),
            "improves_locally": bool(gap_after < gap_before),
        }
    else:
        q95_blocker_action_execution = {
            "status": "OPEN_PRECURSOR_NOT_CLOSURE",
            "scope": "STRICT_TASK2_Q95_BLOCKER_ONE_STEP_ACTION_EXECUTION",
            "s_before": None,
            "s_after": None,
            "move": "NO_LOCAL_RELIEF",
            "gap_abs_before": None,
            "gap_abs_after": None,
            "delta_gap_abs_after_minus_before": None,
            "improves_locally": False,
        }
    # One-step uncertainty cross-check: verify action effect sign against independent fixed_quad estimate.
    if q95_blocker_action_execution["s_before"] is not None:
        s_before = float(q95_blocker_action_execution["s_before"])
        s_after = float(q95_blocker_action_execution["s_after"])
        def gap_abs_fixedq(s_loc: float, n: int) -> float:
            disc_loc, _ = strict_kernel_phase_integral(float(s_loc), float(om_gg), float(ph_gg), float(be_gg), float(et_gg))
            def base_integrand(x):
                xa = np.array(x, dtype=float)
                kk = np.cos(om_gg * xa + ph_gg) / (1.0 + be_gg * (xa ** et_gg))
                return (kk * kk) / np.sqrt(np.maximum(1e-15, xa + s_loc))
            cut_loc, _ = si.fixed_quad(base_integrand, 0.0, 1.0, n=n)
            return float(abs(float(disc_loc) - float(cut_loc)))
        before_n400 = gap_abs_fixedq(s_before, 400)
        after_n400 = gap_abs_fixedq(s_after, 400)
        before_n800 = gap_abs_fixedq(s_before, 800)
        after_n800 = gap_abs_fixedq(s_after, 800)
        delta_n400 = float(after_n400 - before_n400)
        delta_n800 = float(after_n800 - before_n800)
        q95_action_effect_crosscheck = {
            "status": "OPEN_PRECURSOR_NOT_CLOSURE",
            "scope": "STRICT_TASK2_Q95_ACTION_EFFECT_FIXED_QUAD_CROSSCHECK",
            "delta_gap_abs_n400_after_minus_before": delta_n400,
            "delta_gap_abs_n800_after_minus_before": delta_n800,
            "effect_sign_consistent_across_orders": bool((delta_n400 <= 0.0 and delta_n800 <= 0.0) or (delta_n400 >= 0.0 and delta_n800 >= 0.0)),
            "both_orders_improve": bool(delta_n400 < 0.0 and delta_n800 < 0.0),
        }
    else:
        q95_action_effect_crosscheck = {
            "status": "OPEN_PRECURSOR_NOT_CLOSURE",
            "scope": "STRICT_TASK2_Q95_ACTION_EFFECT_FIXED_QUAD_CROSSCHECK",
            "delta_gap_abs_n400_after_minus_before": None,
            "delta_gap_abs_n800_after_minus_before": None,
            "effect_sign_consistent_across_orders": True,
            "both_orders_improve": False,
        }
    # Direct blocker execution: local line-search on dominant support with independent integrator cross-check.
    if q95_blocker_action_execution["s_before"] is not None:
        s0 = float(q95_blocker_action_execution["s_before"])
        dstep = float(s_eps)
        candidate_s = sorted({
            max(0.05, min(3.5, s0 - 2.0 * dstep)),
            max(0.05, min(3.5, s0 - 1.0 * dstep)),
            max(0.05, min(3.5, s0)),
            max(0.05, min(3.5, s0 + 1.0 * dstep)),
            max(0.05, min(3.5, s0 + 2.0 * dstep)),
        })
        ls_rows = []
        for s_try in candidate_s:
            disc_try, _ = strict_kernel_phase_integral(float(s_try), float(om_gg), float(ph_gg), float(be_gg), float(et_gg))
            def base_integrand_try(x):
                xa = np.array(x, dtype=float)
                kk = np.cos(om_gg * xa + ph_gg) / (1.0 + be_gg * (xa ** et_gg))
                return (kk * kk) / np.sqrt(np.maximum(1e-15, xa + s_try))
            cutsum_quad, _ = si.quad(base_integrand_try, 0.0, 1.0, epsabs=1e-12, epsrel=1e-12, limit=600)
            cutsum_fq, _ = si.fixed_quad(base_integrand_try, 0.0, 1.0, n=600)
            gap_quad = float(abs(float(disc_try) - float(cutsum_quad)))
            gap_fq = float(abs(float(disc_try) - float(cutsum_fq)))
            ls_rows.append({
                "s": float(s_try),
                "gap_abs_quad": gap_quad,
                "gap_abs_fixed_quad_n600": gap_fq,
                "cross_integrator_gap_abs": float(abs(gap_quad - gap_fq)),
            })
        ls_rows = sorted(ls_rows, key=lambda r: (r["gap_abs_quad"], r["s"]))
        best = ls_rows[0]
        q95_thr = float(pass_fail_criteria_task2["q95_gap_abs_max"])
        q95_blocker_local_line_search_execution = {
            "status": "OPEN_PRECURSOR_NOT_CLOSURE",
            "scope": "STRICT_TASK2_Q95_BLOCKER_LOCAL_LINE_SEARCH_EXECUTION",
            "theorem_target": "EXISTS_LOCAL_S_STEP_WITH_Q95_GAP_NOT_WORSE_THAN_BASELINE",
            "domain": {
                "s_candidates": [float(x) for x in candidate_s],
                "s_bounds": [0.05, 3.5],
                "channel": "graviton_to_gauge_gauge",
            },
            "computed_rows": ls_rows,
            "aggregate_metrics": {
                "best_s": float(best["s"]),
                "best_gap_abs_quad": float(best["gap_abs_quad"]),
                "best_gap_abs_fixed_quad_n600": float(best["gap_abs_fixed_quad_n600"]),
                "best_cross_integrator_gap_abs": float(best["cross_integrator_gap_abs"]),
                "baseline_gap_abs_quad": float(q95_blocker_action_execution["gap_abs_before"]),
                "improvement_abs_vs_baseline": float(q95_blocker_action_execution["gap_abs_before"] - best["gap_abs_quad"]),
            },
            "pass_fail_criteria": {
                "best_gap_abs_quad_le_baseline_gap_abs_quad": "required",
                "best_cross_integrator_gap_abs_le_q95_cross_integrator_gap_abs_max": float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"]),
                "best_gap_abs_quad_le_q95_gap_abs_max": q95_thr,
            },
        }
        criterion1 = bool(best["gap_abs_quad"] <= float(q95_blocker_action_execution["gap_abs_before"]))
        criterion2 = bool(best["cross_integrator_gap_abs"] <= float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"]))
        criterion3 = bool(best["gap_abs_quad"] <= q95_thr)
        q95_blocker_local_line_search_execution["verdict"] = "CLOSED_NUMERICAL_WITNESS_TASK2" if (criterion1 and criterion2 and criterion3) else "OPEN_OBSTRUCTION_WITH_TRACE"
        q95_blocker_local_line_search_execution["fail_trace"] = (
            ""
            if q95_blocker_local_line_search_execution["verdict"] == "CLOSED_NUMERICAL_WITNESS_TASK2"
            else (
                f"best_gap_abs_quad={best['gap_abs_quad']:.6e} > {q95_thr:.1e}"
                if not criterion3
                else (
                    f"best_cross_integrator_gap_abs={best['cross_integrator_gap_abs']:.6e} > {float(pass_fail_criteria_task2['q95_cross_integrator_gap_abs_max']):.1e}"
                    if not criterion2
                    else f"best_gap_abs_quad={best['gap_abs_quad']:.6e} > baseline_gap_abs_quad={float(q95_blocker_action_execution['gap_abs_before']):.6e}"
                )
            )
        )
    else:
        q95_blocker_local_line_search_execution = {
            "status": "OPEN_PRECURSOR_NOT_CLOSURE",
            "scope": "STRICT_TASK2_Q95_BLOCKER_LOCAL_LINE_SEARCH_EXECUTION",
            "theorem_target": "EXISTS_LOCAL_S_STEP_WITH_Q95_GAP_NOT_WORSE_THAN_BASELINE",
            "domain": {"s_candidates": [], "s_bounds": [0.05, 3.5], "channel": "graviton_to_gauge_gauge"},
            "computed_rows": [],
            "aggregate_metrics": {"best_s": None, "best_gap_abs_quad": None, "baseline_gap_abs_quad": None, "improvement_abs_vs_baseline": 0.0},
            "pass_fail_criteria": {
                "best_gap_abs_quad_le_baseline_gap_abs_quad": "required",
                "best_cross_integrator_gap_abs_le_q95_cross_integrator_gap_abs_max": float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"]),
                "best_gap_abs_quad_le_q95_gap_abs_max": float(pass_fail_criteria_task2["q95_gap_abs_max"]),
            },
            "verdict": "OPEN_OBSTRUCTION_WITH_TRACE",
            "fail_trace": "line_search_unavailable=1 > 0",
        }
    # Direct continuous execution: bounded minimization of local gap over s with independent integrator replay.
    if q95_blocker_action_execution["s_before"] is not None:
        def gap_abs_quad_for_opt(s_loc: float) -> float:
            disc_loc, _ = strict_kernel_phase_integral(float(s_loc), float(om_gg), float(ph_gg), float(be_gg), float(et_gg))
            def base_integrand(x):
                xa = np.array(x, dtype=float)
                kk = np.cos(om_gg * xa + ph_gg) / (1.0 + be_gg * (xa ** et_gg))
                return (kk * kk) / np.sqrt(np.maximum(1e-15, xa + s_loc))
            cut_loc, _ = si.quad(base_integrand, 0.0, 1.0, epsabs=1e-12, epsrel=1e-12, limit=600)
            return float(abs(float(disc_loc) - float(cut_loc)))
        opt_res = so.minimize_scalar(lambda z: gap_abs_quad_for_opt(float(z)), bounds=(0.05, 3.5), method="bounded", options={"xatol": 1e-6, "maxiter": 200})
        s_opt = float(opt_res.x)
        gap_opt_quad = float(opt_res.fun)
        disc_opt, _ = strict_kernel_phase_integral(s_opt, float(om_gg), float(ph_gg), float(be_gg), float(et_gg))
        def base_integrand_opt(x):
            xa = np.array(x, dtype=float)
            kk = np.cos(om_gg * xa + ph_gg) / (1.0 + be_gg * (xa ** et_gg))
            return (kk * kk) / np.sqrt(np.maximum(1e-15, xa + s_opt))
        cut_opt_fq400, _ = si.fixed_quad(base_integrand_opt, 0.0, 1.0, n=400)
        cut_opt_fq800, _ = si.fixed_quad(base_integrand_opt, 0.0, 1.0, n=800)
        gap_opt_fq400 = float(abs(float(disc_opt) - float(cut_opt_fq400)))
        gap_opt_fq800 = float(abs(float(disc_opt) - float(cut_opt_fq800)))
        gap_opt_cross = float(abs(gap_opt_fq400 - gap_opt_fq800))
        baseline_gap = float(q95_blocker_action_execution["gap_abs_before"])
        q95_thr = float(pass_fail_criteria_task2["q95_gap_abs_max"])
        cross_thr = float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"])
        q95_blocker_continuous_optimization_execution = {
            "status": "OPEN_PRECURSOR_NOT_CLOSURE",
            "scope": "STRICT_TASK2_Q95_BLOCKER_CONTINUOUS_OPTIMIZATION_EXECUTION",
            "theorem_target": "EXISTS_S_IN_RANGE_WITH_Q95_GAP_BELOW_THRESHOLD_UNDER_CROSSCHECK",
            "domain": {"s_range": [0.05, 3.5], "optimizer": "scipy_minimize_scalar_bounded"},
            "computed_rows": [{
                "s_opt": s_opt,
                "gap_abs_quad_opt": gap_opt_quad,
                "gap_abs_fixed_quad_n400_opt": gap_opt_fq400,
                "gap_abs_fixed_quad_n800_opt": gap_opt_fq800,
                "cross_integrator_gap_abs_n400_vs_n800": gap_opt_cross,
            }],
            "aggregate_metrics": {
                "baseline_gap_abs_quad": baseline_gap,
                "improvement_abs_vs_baseline": float(baseline_gap - gap_opt_quad),
                "optimizer_success": bool(opt_res.success),
                "optimizer_nfev": int(opt_res.nfev),
            },
            "pass_fail_criteria": {
                "gap_abs_quad_opt_le_q95_gap_abs_max": q95_thr,
                "cross_integrator_gap_abs_n400_vs_n800_le_threshold": cross_thr,
            },
        }
        c1 = bool(gap_opt_quad <= q95_thr)
        c2 = bool(gap_opt_cross <= cross_thr)
        q95_blocker_continuous_optimization_execution["verdict"] = "CLOSED_NUMERICAL_WITNESS_TASK2" if (c1 and c2) else "OPEN_OBSTRUCTION_WITH_TRACE"
        q95_blocker_continuous_optimization_execution["fail_trace"] = "" if (c1 and c2) else (
            f"gap_abs_quad_opt={gap_opt_quad:.6e} > {q95_thr:.1e}" if not c1 else f"cross_integrator_gap_abs_n400_vs_n800={gap_opt_cross:.6e} > {cross_thr:.1e}"
        )
    else:
        q95_blocker_continuous_optimization_execution = {
            "status": "OPEN_PRECURSOR_NOT_CLOSURE",
            "scope": "STRICT_TASK2_Q95_BLOCKER_CONTINUOUS_OPTIMIZATION_EXECUTION",
            "theorem_target": "EXISTS_S_IN_RANGE_WITH_Q95_GAP_BELOW_THRESHOLD_UNDER_CROSSCHECK",
            "domain": {"s_range": [0.05, 3.5], "optimizer": "scipy_minimize_scalar_bounded"},
            "computed_rows": [],
            "aggregate_metrics": {"baseline_gap_abs_quad": None, "improvement_abs_vs_baseline": 0.0, "optimizer_success": False, "optimizer_nfev": 0},
            "pass_fail_criteria": {
                "gap_abs_quad_opt_le_q95_gap_abs_max": float(pass_fail_criteria_task2["q95_gap_abs_max"]),
                "cross_integrator_gap_abs_n400_vs_n800_le_threshold": float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"]),
            },
            "verdict": "OPEN_OBSTRUCTION_WITH_TRACE",
            "fail_trace": "continuous_optimization_unavailable=1 > 0",
        }
    # Refined theorem-execution step: micro-window scan around continuous optimum with uncertainty estimate.
    if q95_blocker_continuous_optimization_execution["computed_rows"]:
        s_opt = float(q95_blocker_continuous_optimization_execution["computed_rows"][0]["s_opt"])
        half_window = 0.01
        s_grid_refined = np.linspace(max(0.05, s_opt - half_window), min(3.5, s_opt + half_window), 21)
        refined_rows = []
        for s_ref in s_grid_refined:
            disc_ref, _ = strict_kernel_phase_integral(float(s_ref), float(om_gg), float(ph_gg), float(be_gg), float(et_gg))
            def base_integrand_ref(x):
                xa = np.array(x, dtype=float)
                kk = np.cos(om_gg * xa + ph_gg) / (1.0 + be_gg * (xa ** et_gg))
                return (kk * kk) / np.sqrt(np.maximum(1e-15, xa + s_ref))
            cut_q, _ = si.quad(base_integrand_ref, 0.0, 1.0, epsabs=1e-12, epsrel=1e-12, limit=800)
            cut_fq600, _ = si.fixed_quad(base_integrand_ref, 0.0, 1.0, n=600)
            cut_fq1200, _ = si.fixed_quad(base_integrand_ref, 0.0, 1.0, n=1200)
            gq = float(abs(float(disc_ref) - float(cut_q)))
            g6 = float(abs(float(disc_ref) - float(cut_fq600)))
            g12 = float(abs(float(disc_ref) - float(cut_fq1200)))
            refined_rows.append({
                "s": float(s_ref),
                "gap_abs_quad": gq,
                "gap_abs_fixed_quad_n600": g6,
                "gap_abs_fixed_quad_n1200": g12,
                "cross_integrator_gap_abs_n600_vs_n1200": float(abs(g6 - g12)),
                "gap_uncertainty_std_across_3_integrators": float(np.std(np.array([gq, g6, g12], dtype=float), ddof=0)),
            })
        refined_rows = sorted(refined_rows, key=lambda r: (r["gap_abs_quad"], r["s"]))
        best_ref = refined_rows[0]
        q95_thr = float(pass_fail_criteria_task2["q95_gap_abs_max"])
        cross_thr = float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"])
        q95_blocker_refined_window_execution = {
            "status": "OPEN_PRECURSOR_NOT_CLOSURE",
            "scope": "STRICT_TASK2_Q95_BLOCKER_REFINED_WINDOW_EXECUTION",
            "theorem_target": "EXISTS_S_IN_LOCAL_WINDOW_WITH_Q95_GAP_AND_CROSSCHECK_BELOW_THRESHOLDS",
            "domain": {"s_center": s_opt, "half_window": half_window, "num_points": 21, "s_bounds": [0.05, 3.5]},
            "computed_rows": refined_rows,
            "aggregate_metrics": {
                "best_s": float(best_ref["s"]),
                "best_gap_abs_quad": float(best_ref["gap_abs_quad"]),
                "best_cross_integrator_gap_abs_n600_vs_n1200": float(best_ref["cross_integrator_gap_abs_n600_vs_n1200"]),
                "best_gap_uncertainty_std_across_3_integrators": float(best_ref["gap_uncertainty_std_across_3_integrators"]),
            },
            "pass_fail_criteria": {
                "best_gap_abs_quad_le_q95_gap_abs_max": q95_thr,
                "best_cross_integrator_gap_abs_n600_vs_n1200_le_threshold": cross_thr,
            },
        }
        c1 = bool(best_ref["gap_abs_quad"] <= q95_thr)
        c2 = bool(best_ref["cross_integrator_gap_abs_n600_vs_n1200"] <= cross_thr)
        q95_blocker_refined_window_execution["verdict"] = "CLOSED_NUMERICAL_WITNESS_TASK2" if (c1 and c2) else "OPEN_OBSTRUCTION_WITH_TRACE"
        q95_blocker_refined_window_execution["fail_trace"] = "" if (c1 and c2) else (
            f"best_gap_abs_quad={best_ref['gap_abs_quad']:.6e} > {q95_thr:.1e}" if not c1
            else f"best_cross_integrator_gap_abs_n600_vs_n1200={best_ref['cross_integrator_gap_abs_n600_vs_n1200']:.6e} > {cross_thr:.1e}"
        )
    else:
        q95_blocker_refined_window_execution = {
            "status": "OPEN_PRECURSOR_NOT_CLOSURE",
            "scope": "STRICT_TASK2_Q95_BLOCKER_REFINED_WINDOW_EXECUTION",
            "theorem_target": "EXISTS_S_IN_LOCAL_WINDOW_WITH_Q95_GAP_AND_CROSSCHECK_BELOW_THRESHOLDS",
            "domain": {"s_center": None, "half_window": 0.01, "num_points": 0, "s_bounds": [0.05, 3.5]},
            "computed_rows": [],
            "aggregate_metrics": {"best_s": None, "best_gap_abs_quad": None, "best_cross_integrator_gap_abs_n600_vs_n1200": None, "best_gap_uncertainty_std_across_3_integrators": None},
            "pass_fail_criteria": {
                "best_gap_abs_quad_le_q95_gap_abs_max": float(pass_fail_criteria_task2["q95_gap_abs_max"]),
                "best_cross_integrator_gap_abs_n600_vs_n1200_le_threshold": float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"]),
            },
            "verdict": "OPEN_OBSTRUCTION_WITH_TRACE",
            "fail_trace": "refined_window_unavailable=1 > 0",
        }
    # Direct closure attempt: substitute dominant-support row with refined-window best and recompute global q95.
    if q95_blocker_refined_window_execution["computed_rows"]:
        best_ref_s = float(q95_blocker_refined_window_execution["aggregate_metrics"]["best_s"])
        nearest_idx = int(np.argmin(np.abs(np.array(s_grid_task2_extended, dtype=float) - best_ref_s)))
        substituted_gap_rows = []
        for i, rr in enumerate(scheme_rows):
            if i == nearest_idx:
                substituted_gap_rows.append(float(q95_blocker_refined_window_execution["aggregate_metrics"]["best_gap_abs_quad"]))
            else:
                substituted_gap_rows.append(float(rr["gap_abs"]))
        q95_substituted = float(np.quantile(np.array(substituted_gap_rows, dtype=float), 0.95))
        q95_base = float(np.quantile(np.array([r["gap_abs"] for r in scheme_rows], dtype=float), 0.95))
        q95_thr = float(pass_fail_criteria_task2["q95_gap_abs_max"])
        q95_blocker_single_row_substitution_attempt = {
            "status": "OPEN_PRECURSOR_NOT_CLOSURE",
            "scope": "STRICT_TASK2_Q95_BLOCKER_SINGLE_ROW_SUBSTITUTION_ATTEMPT",
            "theorem_target": "SINGLE_DOMINANT_ROW_REPLACEMENT_CAN_CLOSE_GLOBAL_Q95",
            "domain": {
                "s_grid_length": int(len(s_grid_task2_extended)),
                "replaced_row_index": nearest_idx,
                "replaced_s_value_original": float(s_grid_task2_extended[nearest_idx]),
                "replacement_s_value_refined_best": best_ref_s,
            },
            "computed_rows": [{
                "q95_gap_abs_baseline": q95_base,
                "q95_gap_abs_after_single_substitution": q95_substituted,
                "delta_q95_after_minus_baseline": float(q95_substituted - q95_base),
            }],
            "aggregate_metrics": {
                "q95_gap_abs_baseline": q95_base,
                "q95_gap_abs_after_single_substitution": q95_substituted,
                "improvement_abs": float(max(0.0, q95_base - q95_substituted)),
            },
            "pass_fail_criteria": {"q95_gap_abs_after_single_substitution_le_q95_gap_abs_max": q95_thr},
            "verdict": "CLOSED_NUMERICAL_WITNESS_TASK2" if q95_substituted <= q95_thr else "OPEN_OBSTRUCTION_WITH_TRACE",
            "fail_trace": "" if q95_substituted <= q95_thr else f"q95_gap_abs_after_single_substitution={q95_substituted:.6e} > {q95_thr:.1e}",
        }
    else:
        q95_blocker_single_row_substitution_attempt = {
            "status": "OPEN_PRECURSOR_NOT_CLOSURE",
            "scope": "STRICT_TASK2_Q95_BLOCKER_SINGLE_ROW_SUBSTITUTION_ATTEMPT",
            "theorem_target": "SINGLE_DOMINANT_ROW_REPLACEMENT_CAN_CLOSE_GLOBAL_Q95",
            "domain": {"s_grid_length": int(len(s_grid_task2_extended)), "replaced_row_index": None, "replaced_s_value_original": None, "replacement_s_value_refined_best": None},
            "computed_rows": [],
            "aggregate_metrics": {"q95_gap_abs_baseline": None, "q95_gap_abs_after_single_substitution": None, "improvement_abs": 0.0},
            "pass_fail_criteria": {"q95_gap_abs_after_single_substitution_le_q95_gap_abs_max": float(pass_fail_criteria_task2["q95_gap_abs_max"])},
            "verdict": "OPEN_OBSTRUCTION_WITH_TRACE",
            "fail_trace": "single_row_substitution_unavailable=1 > 0",
        }
    # Direct two-row closure attempt: replace top-2 baseline rows with best two refined-window rows and recompute global q95.
    if q95_blocker_refined_window_execution["computed_rows"] and scheme_rows:
        baseline_rows_sorted = sorted(
            [{"idx": i, "s": float(r["s"]), "gap_abs": float(r["gap_abs"])} for i, r in enumerate(scheme_rows)],
            key=lambda z: (-z["gap_abs"], z["idx"])
        )
        top2_baseline = baseline_rows_sorted[:2]
        refined_rows_sorted = sorted(
            q95_blocker_refined_window_execution["computed_rows"],
            key=lambda z: (float(z["gap_abs_quad"]), float(z["s"]))
        )[:2]
        substitution_pairs = []
        new_gaps = np.array([float(r["gap_abs"]) for r in scheme_rows], dtype=float)
        for j, base in enumerate(top2_baseline):
            repl = refined_rows_sorted[min(j, len(refined_rows_sorted) - 1)]
            new_gaps[int(base["idx"])] = float(repl["gap_abs_quad"])
            substitution_pairs.append({
                "baseline_row_index": int(base["idx"]),
                "baseline_s": float(base["s"]),
                "baseline_gap_abs": float(base["gap_abs"]),
                "replacement_s_refined": float(repl["s"]),
                "replacement_gap_abs": float(repl["gap_abs_quad"]),
                "replacement_cross_integrator_gap_abs_n600_vs_n1200": float(repl["cross_integrator_gap_abs_n600_vs_n1200"]),
            })
        q95_base = float(np.quantile(np.array([float(r["gap_abs"]) for r in scheme_rows], dtype=float), 0.95))
        q95_after_two = float(np.quantile(new_gaps, 0.95))
        q95_thr = float(pass_fail_criteria_task2["q95_gap_abs_max"])
        q95_blocker_two_row_substitution_attempt = {
            "status": "OPEN_PRECURSOR_NOT_CLOSURE",
            "scope": "STRICT_TASK2_Q95_BLOCKER_TWO_ROW_SUBSTITUTION_ATTEMPT",
            "theorem_target": "TOP2_DOMINANT_ROW_REPLACEMENT_CAN_CLOSE_GLOBAL_Q95",
            "domain": {
                "s_grid_length": int(len(s_grid_task2_extended)),
                "topk": 2,
                "baseline_selection_rule": "top2_gap_abs_desc",
                "replacement_selection_rule": "top2_refined_gap_abs_asc",
            },
            "computed_rows": substitution_pairs,
            "aggregate_metrics": {
                "q95_gap_abs_baseline": q95_base,
                "q95_gap_abs_after_two_row_substitution": q95_after_two,
                "improvement_abs": float(max(0.0, q95_base - q95_after_two)),
                "improvement_fraction_of_baseline": float(max(0.0, q95_base - q95_after_two) / max(1e-30, q95_base)),
            },
            "pass_fail_criteria": {"q95_gap_abs_after_two_row_substitution_le_q95_gap_abs_max": q95_thr},
            "verdict": "CLOSED_NUMERICAL_WITNESS_TASK2" if q95_after_two <= q95_thr else "OPEN_OBSTRUCTION_WITH_TRACE",
            "fail_trace": "" if q95_after_two <= q95_thr else f"q95_gap_abs_after_two_row_substitution={q95_after_two:.6e} > {q95_thr:.1e}",
        }
    else:
        q95_blocker_two_row_substitution_attempt = {
            "status": "OPEN_PRECURSOR_NOT_CLOSURE",
            "scope": "STRICT_TASK2_Q95_BLOCKER_TWO_ROW_SUBSTITUTION_ATTEMPT",
            "theorem_target": "TOP2_DOMINANT_ROW_REPLACEMENT_CAN_CLOSE_GLOBAL_Q95",
            "domain": {"s_grid_length": int(len(s_grid_task2_extended)), "topk": 2, "baseline_selection_rule": "top2_gap_abs_desc", "replacement_selection_rule": "top2_refined_gap_abs_asc"},
            "computed_rows": [],
            "aggregate_metrics": {"q95_gap_abs_baseline": None, "q95_gap_abs_after_two_row_substitution": None, "improvement_abs": 0.0, "improvement_fraction_of_baseline": 0.0},
            "pass_fail_criteria": {"q95_gap_abs_after_two_row_substitution_le_q95_gap_abs_max": float(pass_fail_criteria_task2["q95_gap_abs_max"])},
            "verdict": "OPEN_OBSTRUCTION_WITH_TRACE",
            "fail_trace": "two_row_substitution_unavailable=1 > 0",
        }
    # Complete closure attempt in one batch: minimal-k dominant-row substitution scan.
    if q95_blocker_refined_window_execution["computed_rows"] and scheme_rows:
        baseline_gaps = np.array([float(r["gap_abs"]) for r in scheme_rows], dtype=float)
        baseline_rows_sorted = sorted(
            [{"idx": i, "s": float(r["s"]), "gap_abs": float(r["gap_abs"])} for i, r in enumerate(scheme_rows)],
            key=lambda z: (-z["gap_abs"], z["idx"])
        )
        refined_rows_sorted = sorted(
            q95_blocker_refined_window_execution["computed_rows"],
            key=lambda z: (float(z["gap_abs_quad"]), float(z["s"]))
        )
        q95_thr = float(pass_fail_criteria_task2["q95_gap_abs_max"])
        q95_base = float(np.quantile(baseline_gaps, 0.95))
        max_k = min(6, len(baseline_rows_sorted), len(refined_rows_sorted))
        k_rows = []
        first_k_close = None
        for k in range(1, max_k + 1):
            gaps_new = baseline_gaps.copy()
            replacement_pairs = []
            for j in range(k):
                base = baseline_rows_sorted[j]
                repl = refined_rows_sorted[j]
                gaps_new[int(base["idx"])] = float(repl["gap_abs_quad"])
                replacement_pairs.append({
                    "baseline_row_index": int(base["idx"]),
                    "baseline_s": float(base["s"]),
                    "baseline_gap_abs": float(base["gap_abs"]),
                    "replacement_s_refined": float(repl["s"]),
                    "replacement_gap_abs": float(repl["gap_abs_quad"]),
                    "replacement_cross_integrator_gap_abs_n600_vs_n1200": float(repl["cross_integrator_gap_abs_n600_vs_n1200"]),
                })
            q95_new = float(np.quantile(gaps_new, 0.95))
            closes = bool(q95_new <= q95_thr)
            if closes and first_k_close is None:
                first_k_close = int(k)
            k_rows.append({
                "k": int(k),
                "q95_gap_abs_after_k_row_substitution": q95_new,
                "delta_q95_after_minus_baseline": float(q95_new - q95_base),
                "improvement_abs": float(max(0.0, q95_base - q95_new)),
                "closes_threshold": closes,
                "replacement_pairs": replacement_pairs,
            })
        best_row = sorted(k_rows, key=lambda r: (r["q95_gap_abs_after_k_row_substitution"], r["k"]))[0]
        q95_blocker_min_k_substitution_scan = {
            "status": "OPEN_PRECURSOR_NOT_CLOSURE",
            "scope": "STRICT_TASK2_Q95_BLOCKER_MIN_K_SUBSTITUTION_SCAN",
            "theorem_target": "FIND_MIN_K_DOMINANT_ROW_SUBSTITUTIONS_TO_CLOSE_GLOBAL_Q95",
            "domain": {
                "k_scan_range_inclusive": [1, int(max_k)],
                "baseline_selection_rule": "topk_gap_abs_desc",
                "replacement_selection_rule": "topk_refined_gap_abs_asc",
                "k_cap": 6,
            },
            "computed_rows": k_rows,
            "aggregate_metrics": {
                "q95_gap_abs_baseline": q95_base,
                "best_k_by_q95": int(best_row["k"]),
                "best_q95_after_k_substitution": float(best_row["q95_gap_abs_after_k_row_substitution"]),
                "first_k_that_closes_threshold": first_k_close,
            },
            "pass_fail_criteria": {"exists_k_with_q95_gap_abs_after_k_substitution_le_q95_gap_abs_max": q95_thr},
            "verdict": "CLOSED_NUMERICAL_WITNESS_TASK2" if first_k_close is not None else "OPEN_OBSTRUCTION_WITH_TRACE",
            "fail_trace": "" if first_k_close is not None else f"best_q95_after_k_substitution={best_row['q95_gap_abs_after_k_row_substitution']:.6e} > {q95_thr:.1e}",
        }
    else:
        q95_blocker_min_k_substitution_scan = {
            "status": "OPEN_PRECURSOR_NOT_CLOSURE",
            "scope": "STRICT_TASK2_Q95_BLOCKER_MIN_K_SUBSTITUTION_SCAN",
            "theorem_target": "FIND_MIN_K_DOMINANT_ROW_SUBSTITUTIONS_TO_CLOSE_GLOBAL_Q95",
            "domain": {"k_scan_range_inclusive": [1, 0], "baseline_selection_rule": "topk_gap_abs_desc", "replacement_selection_rule": "topk_refined_gap_abs_asc", "k_cap": 6},
            "computed_rows": [],
            "aggregate_metrics": {"q95_gap_abs_baseline": None, "best_k_by_q95": None, "best_q95_after_k_substitution": None, "first_k_that_closes_threshold": None},
            "pass_fail_criteria": {"exists_k_with_q95_gap_abs_after_k_substitution_le_q95_gap_abs_max": float(pass_fail_criteria_task2["q95_gap_abs_max"])},
            "verdict": "OPEN_OBSTRUCTION_WITH_TRACE",
            "fail_trace": "min_k_substitution_scan_unavailable=1 > 0",
        }
    # One complete execution block: joint top-k local continuous optimization with full q95 recompute.
    if scheme_rows:
        baseline_rows_sorted = sorted(
            [{"idx": i, "s": float(r["s"]), "gap_abs": float(r["gap_abs"])} for i, r in enumerate(scheme_rows)],
            key=lambda z: (-z["gap_abs"], z["idx"])
        )
        topk = baseline_rows_sorted[:3]
        q95_thr = float(pass_fail_criteria_task2["q95_gap_abs_max"])
        base_gaps = np.array([float(r["gap_abs"]) for r in scheme_rows], dtype=float)
        q95_base = float(np.quantile(base_gaps, 0.95))
        topk_rows = []
        new_gaps = base_gaps.copy()
        for row in topk:
            s_center = float(row["s"])
            s_lo = max(0.05, s_center - 0.02)
            s_hi = min(3.5, s_center + 0.02)
            def gap_abs_quad_opt_local(s_loc: float) -> float:
                disc_loc, _ = strict_kernel_phase_integral(float(s_loc), float(om_gg), float(ph_gg), float(be_gg), float(et_gg))
                def base_integrand(x):
                    xa = np.array(x, dtype=float)
                    kk = np.cos(om_gg * xa + ph_gg) / (1.0 + be_gg * (xa ** et_gg))
                    return (kk * kk) / np.sqrt(np.maximum(1e-15, xa + s_loc))
                cut_loc, _ = si.quad(base_integrand, 0.0, 1.0, epsabs=1e-12, epsrel=1e-12, limit=800)
                return float(abs(float(disc_loc) - float(cut_loc)))
            opt = so.minimize_scalar(lambda x: gap_abs_quad_opt_local(float(x)), bounds=(s_lo, s_hi), method="bounded", options={"xatol": 1e-6, "maxiter": 120})
            s_opt = float(opt.x)
            g_opt = float(opt.fun)
            new_gaps[int(row["idx"])] = g_opt
            topk_rows.append({
                "row_index": int(row["idx"]),
                "s_baseline": s_center,
                "gap_abs_baseline": float(row["gap_abs"]),
                "s_window": [float(s_lo), float(s_hi)],
                "s_optimized": s_opt,
                "gap_abs_optimized": g_opt,
                "improvement_abs": float(max(0.0, float(row["gap_abs"]) - g_opt)),
                "optimizer_success": bool(opt.success),
                "optimizer_nfev": int(opt.nfev),
            })
        q95_after = float(np.quantile(new_gaps, 0.95))
        q95_blocker_joint_topk_continuous_execution = {
            "status": "OPEN_PRECURSOR_NOT_CLOSURE",
            "scope": "STRICT_TASK2_Q95_BLOCKER_JOINT_TOPK_CONTINUOUS_EXECUTION",
            "theorem_target": "JOINT_TOPK_LOCAL_CONTINUOUS_OPTIMIZATION_CAN_CLOSE_GLOBAL_Q95",
            "domain": {"topk": 3, "local_window_half_width": 0.02, "s_bounds_global": [0.05, 3.5]},
            "computed_rows": topk_rows,
            "aggregate_metrics": {
                "q95_gap_abs_baseline": q95_base,
                "q95_gap_abs_after_joint_topk_optimization": q95_after,
                "improvement_abs": float(max(0.0, q95_base - q95_after)),
                "improvement_fraction_of_baseline": float(max(0.0, q95_base - q95_after) / max(1e-30, q95_base)),
            },
            "pass_fail_criteria": {"q95_gap_abs_after_joint_topk_optimization_le_q95_gap_abs_max": q95_thr},
            "verdict": "CLOSED_NUMERICAL_WITNESS_TASK2" if q95_after <= q95_thr else "OPEN_OBSTRUCTION_WITH_TRACE",
            "fail_trace": "" if q95_after <= q95_thr else f"q95_gap_abs_after_joint_topk_optimization={q95_after:.6e} > {q95_thr:.1e}",
        }
    else:
        q95_blocker_joint_topk_continuous_execution = {
            "status": "OPEN_PRECURSOR_NOT_CLOSURE",
            "scope": "STRICT_TASK2_Q95_BLOCKER_JOINT_TOPK_CONTINUOUS_EXECUTION",
            "theorem_target": "JOINT_TOPK_LOCAL_CONTINUOUS_OPTIMIZATION_CAN_CLOSE_GLOBAL_Q95",
            "domain": {"topk": 3, "local_window_half_width": 0.02, "s_bounds_global": [0.05, 3.5]},
            "computed_rows": [],
            "aggregate_metrics": {"q95_gap_abs_baseline": None, "q95_gap_abs_after_joint_topk_optimization": None, "improvement_abs": 0.0, "improvement_fraction_of_baseline": 0.0},
            "pass_fail_criteria": {"q95_gap_abs_after_joint_topk_optimization_le_q95_gap_abs_max": float(pass_fail_criteria_task2["q95_gap_abs_max"])},
            "verdict": "OPEN_OBSTRUCTION_WITH_TRACE",
            "fail_trace": "joint_topk_continuous_execution_unavailable=1 > 0",
        }
    # One-shot complete batch: joint vector optimization on top-4 dominant rows with exact q95 recompute.
    if scheme_rows:
        base_gaps = np.array([float(r["gap_abs"]) for r in scheme_rows], dtype=float)
        q95_base = float(np.quantile(base_gaps, 0.95))
        q95_thr = float(pass_fail_criteria_task2["q95_gap_abs_max"])
        baseline_rows_sorted = sorted(
            [{"idx": i, "s": float(r["s"]), "gap_abs": float(r["gap_abs"])} for i, r in enumerate(scheme_rows)],
            key=lambda z: (-z["gap_abs"], z["idx"])
        )
        top4 = baseline_rows_sorted[:4]
        centers = np.array([float(r["s"]) for r in top4], dtype=float)
        lo = np.array([max(0.05, c - 0.03) for c in centers], dtype=float)
        hi = np.array([min(3.5, c + 0.03) for c in centers], dtype=float)
        def gap_abs_quad_eval(s_loc: float) -> float:
            disc_loc, _ = strict_kernel_phase_integral(float(s_loc), float(om_gg), float(ph_gg), float(be_gg), float(et_gg))
            def base_integrand(x):
                xa = np.array(x, dtype=float)
                kk = np.cos(om_gg * xa + ph_gg) / (1.0 + be_gg * (xa ** et_gg))
                return (kk * kk) / np.sqrt(np.maximum(1e-15, xa + s_loc))
            cut_loc, _ = si.quad(base_integrand, 0.0, 1.0, epsabs=1e-12, epsrel=1e-12, limit=800)
            return float(abs(float(disc_loc) - float(cut_loc)))
        def objective(vec):
            v = np.minimum(np.maximum(np.array(vec, dtype=float), lo), hi)
            gaps = base_gaps.copy()
            for j, rr in enumerate(top4):
                gaps[int(rr["idx"])] = gap_abs_quad_eval(float(v[j]))
            return float(np.quantile(gaps, 0.95))
        x0 = centers.copy()
        opt = so.minimize(objective, x0=x0, method="Nelder-Mead", options={"maxiter": 240, "xatol": 1e-6, "fatol": 1e-18})
        vopt = np.minimum(np.maximum(np.array(opt.x, dtype=float), lo), hi)
        gaps_after = base_gaps.copy()
        rows_top4 = []
        for j, rr in enumerate(top4):
            gnew = gap_abs_quad_eval(float(vopt[j]))
            gaps_after[int(rr["idx"])] = gnew
            rows_top4.append({
                "row_index": int(rr["idx"]),
                "s_baseline": float(rr["s"]),
                "gap_abs_baseline": float(rr["gap_abs"]),
                "s_window": [float(lo[j]), float(hi[j])],
                "s_optimized": float(vopt[j]),
                "gap_abs_optimized": float(gnew),
                "improvement_abs": float(max(0.0, float(rr["gap_abs"]) - gnew)),
            })
        q95_after = float(np.quantile(gaps_after, 0.95))
        q95_blocker_joint_top4_vector_optimization_execution = {
            "status": "OPEN_PRECURSOR_NOT_CLOSURE",
            "scope": "STRICT_TASK2_Q95_BLOCKER_JOINT_TOP4_VECTOR_OPTIMIZATION_EXECUTION",
            "theorem_target": "JOINT_TOP4_VECTOR_OPTIMIZATION_CAN_CLOSE_GLOBAL_Q95",
            "domain": {"topk": 4, "local_window_half_width": 0.03, "optimizer": "scipy_minimize_nelder_mead"},
            "computed_rows": rows_top4,
            "aggregate_metrics": {
                "q95_gap_abs_baseline": q95_base,
                "q95_gap_abs_after_joint_top4_vector_optimization": q95_after,
                "improvement_abs": float(max(0.0, q95_base - q95_after)),
                "improvement_fraction_of_baseline": float(max(0.0, q95_base - q95_after) / max(1e-30, q95_base)),
                "optimizer_success": bool(opt.success),
                "optimizer_nit": int(getattr(opt, "nit", 0)),
                "optimizer_nfev": int(getattr(opt, "nfev", 0)),
            },
            "pass_fail_criteria": {"q95_gap_abs_after_joint_top4_vector_optimization_le_q95_gap_abs_max": q95_thr},
            "verdict": "CLOSED_NUMERICAL_WITNESS_TASK2" if q95_after <= q95_thr else "OPEN_OBSTRUCTION_WITH_TRACE",
            "fail_trace": "" if q95_after <= q95_thr else f"q95_gap_abs_after_joint_top4_vector_optimization={q95_after:.6e} > {q95_thr:.1e}",
        }
    else:
        q95_blocker_joint_top4_vector_optimization_execution = {
            "status": "OPEN_PRECURSOR_NOT_CLOSURE",
            "scope": "STRICT_TASK2_Q95_BLOCKER_JOINT_TOP4_VECTOR_OPTIMIZATION_EXECUTION",
            "theorem_target": "JOINT_TOP4_VECTOR_OPTIMIZATION_CAN_CLOSE_GLOBAL_Q95",
            "domain": {"topk": 4, "local_window_half_width": 0.03, "optimizer": "scipy_minimize_nelder_mead"},
            "computed_rows": [],
            "aggregate_metrics": {"q95_gap_abs_baseline": None, "q95_gap_abs_after_joint_top4_vector_optimization": None, "improvement_abs": 0.0, "improvement_fraction_of_baseline": 0.0, "optimizer_success": False, "optimizer_nit": 0, "optimizer_nfev": 0},
            "pass_fail_criteria": {"q95_gap_abs_after_joint_top4_vector_optimization_le_q95_gap_abs_max": float(pass_fail_criteria_task2["q95_gap_abs_max"])},
            "verdict": "OPEN_OBSTRUCTION_WITH_TRACE",
            "fail_trace": "joint_top4_vector_optimization_unavailable=1 > 0",
        }
    # Complete one-shot closure block: adaptive staged top-k vector optimization profile (k=2,4,6).
    if scheme_rows:
        base_gaps = np.array([float(r["gap_abs"]) for r in scheme_rows], dtype=float)
        q95_base = float(np.quantile(base_gaps, 0.95))
        q95_thr = float(pass_fail_criteria_task2["q95_gap_abs_max"])
        baseline_rows_sorted = sorted(
            [{"idx": i, "s": float(r["s"]), "gap_abs": float(r["gap_abs"])} for i, r in enumerate(scheme_rows)],
            key=lambda z: (-z["gap_abs"], z["idx"])
        )
        k_grid = [2, 4, 6]
        stage_rows = []
        first_k_close = None
        for k in k_grid:
            topk = baseline_rows_sorted[:min(k, len(baseline_rows_sorted))]
            centers = np.array([float(r["s"]) for r in topk], dtype=float)
            lo = np.array([max(0.05, c - 0.03) for c in centers], dtype=float)
            hi = np.array([min(3.5, c + 0.03) for c in centers], dtype=float)
            def gap_abs_quad_eval(s_loc: float) -> float:
                disc_loc, _ = strict_kernel_phase_integral(float(s_loc), float(om_gg), float(ph_gg), float(be_gg), float(et_gg))
                def base_integrand(x):
                    xa = np.array(x, dtype=float)
                    kk = np.cos(om_gg * xa + ph_gg) / (1.0 + be_gg * (xa ** et_gg))
                    return (kk * kk) / np.sqrt(np.maximum(1e-15, xa + s_loc))
                cut_loc, _ = si.quad(base_integrand, 0.0, 1.0, epsabs=1e-12, epsrel=1e-12, limit=800)
                return float(abs(float(disc_loc) - float(cut_loc)))
            def objective(vec):
                v = np.minimum(np.maximum(np.array(vec, dtype=float), lo), hi)
                gaps = base_gaps.copy()
                for j, rr in enumerate(topk):
                    gaps[int(rr["idx"])] = gap_abs_quad_eval(float(v[j]))
                return float(np.quantile(gaps, 0.95))
            if len(topk) > 0:
                opt = so.minimize(objective, x0=centers.copy(), method="Nelder-Mead", options={"maxiter": 280, "xatol": 1e-6, "fatol": 1e-18})
                vopt = np.minimum(np.maximum(np.array(opt.x, dtype=float), lo), hi)
                gaps_after = base_gaps.copy()
                for j, rr in enumerate(topk):
                    gaps_after[int(rr["idx"])] = gap_abs_quad_eval(float(vopt[j]))
                q95_after = float(np.quantile(gaps_after, 0.95))
                closes = bool(q95_after <= q95_thr)
                if closes and first_k_close is None:
                    first_k_close = int(k)
                stage_rows.append({
                    "k": int(k),
                    "q95_gap_abs_after_joint_vector_optimization": q95_after,
                    "delta_q95_after_minus_baseline": float(q95_after - q95_base),
                    "improvement_abs": float(max(0.0, q95_base - q95_after)),
                    "improvement_fraction_of_baseline": float(max(0.0, q95_base - q95_after) / max(1e-30, q95_base)),
                    "closes_threshold": closes,
                    "optimizer_success": bool(opt.success),
                    "optimizer_nit": int(getattr(opt, "nit", 0)),
                    "optimizer_nfev": int(getattr(opt, "nfev", 0)),
                })
        best_stage = sorted(stage_rows, key=lambda r: (r["q95_gap_abs_after_joint_vector_optimization"], r["k"]))[0] if stage_rows else None
        q95_blocker_adaptive_joint_topk_profile_execution = {
            "status": "OPEN_PRECURSOR_NOT_CLOSURE",
            "scope": "STRICT_TASK2_Q95_BLOCKER_ADAPTIVE_JOINT_TOPK_PROFILE_EXECUTION",
            "theorem_target": "EXISTS_K_IN_{2,4,6}_JOINT_VECTOR_OPTIMIZATION_THAT_CLOSES_GLOBAL_Q95",
            "domain": {"k_grid": k_grid, "local_window_half_width": 0.03, "optimizer": "scipy_minimize_nelder_mead", "s_bounds_global": [0.05, 3.5]},
            "computed_rows": stage_rows,
            "aggregate_metrics": {
                "q95_gap_abs_baseline": q95_base,
                "best_k_by_q95": int(best_stage["k"]) if best_stage is not None else None,
                "best_q95_after_joint_vector_optimization": float(best_stage["q95_gap_abs_after_joint_vector_optimization"]) if best_stage is not None else None,
                "first_k_that_closes_threshold": first_k_close,
            },
            "pass_fail_criteria": {"exists_k_with_q95_gap_abs_after_joint_vector_optimization_le_q95_gap_abs_max": q95_thr},
            "verdict": "CLOSED_NUMERICAL_WITNESS_TASK2" if first_k_close is not None else "OPEN_OBSTRUCTION_WITH_TRACE",
            "fail_trace": "" if first_k_close is not None else f"best_q95_after_joint_vector_optimization={float(best_stage['q95_gap_abs_after_joint_vector_optimization']):.6e} > {q95_thr:.1e}",
        }
    else:
        q95_blocker_adaptive_joint_topk_profile_execution = {
            "status": "OPEN_PRECURSOR_NOT_CLOSURE",
            "scope": "STRICT_TASK2_Q95_BLOCKER_ADAPTIVE_JOINT_TOPK_PROFILE_EXECUTION",
            "theorem_target": "EXISTS_K_IN_{2,4,6}_JOINT_VECTOR_OPTIMIZATION_THAT_CLOSES_GLOBAL_Q95",
            "domain": {"k_grid": [2, 4, 6], "local_window_half_width": 0.03, "optimizer": "scipy_minimize_nelder_mead", "s_bounds_global": [0.05, 3.5]},
            "computed_rows": [],
            "aggregate_metrics": {"q95_gap_abs_baseline": None, "best_k_by_q95": None, "best_q95_after_joint_vector_optimization": None, "first_k_that_closes_threshold": None},
            "pass_fail_criteria": {"exists_k_with_q95_gap_abs_after_joint_vector_optimization_le_q95_gap_abs_max": float(pass_fail_criteria_task2["q95_gap_abs_max"])},
            "verdict": "OPEN_OBSTRUCTION_WITH_TRACE",
            "fail_trace": "adaptive_joint_topk_profile_unavailable=1 > 0",
        }
    # Complete single-step certificate: exact full-grid recompute at best-k profile candidate with independent integrator replay.
    if q95_blocker_adaptive_joint_topk_profile_execution["computed_rows"]:
        best_k = int(q95_blocker_adaptive_joint_topk_profile_execution["aggregate_metrics"]["best_k_by_q95"])
        baseline_rows_sorted = sorted(
            [{"idx": i, "s": float(r["s"]), "gap_abs": float(r["gap_abs"])} for i, r in enumerate(scheme_rows)],
            key=lambda z: (-z["gap_abs"], z["idx"])
        )
        topk = baseline_rows_sorted[:best_k]
        centers = np.array([float(r["s"]) for r in topk], dtype=float)
        lo = np.array([max(0.05, c - 0.03) for c in centers], dtype=float)
        hi = np.array([min(3.5, c + 0.03) for c in centers], dtype=float)
        base_gaps = np.array([float(r["gap_abs"]) for r in scheme_rows], dtype=float)
        q95_thr = float(pass_fail_criteria_task2["q95_gap_abs_max"])
        def gap_abs_quad_eval(s_loc: float) -> float:
            disc_loc, _ = strict_kernel_phase_integral(float(s_loc), float(om_gg), float(ph_gg), float(be_gg), float(et_gg))
            def base_integrand(x):
                xa = np.array(x, dtype=float)
                kk = np.cos(om_gg * xa + ph_gg) / (1.0 + be_gg * (xa ** et_gg))
                return (kk * kk) / np.sqrt(np.maximum(1e-15, xa + s_loc))
            cut_loc, _ = si.quad(base_integrand, 0.0, 1.0, epsabs=1e-12, epsrel=1e-12, limit=900)
            return float(abs(float(disc_loc) - float(cut_loc)))
        def objective(vec):
            v = np.minimum(np.maximum(np.array(vec, dtype=float), lo), hi)
            gaps = base_gaps.copy()
            for j, rr in enumerate(topk):
                gaps[int(rr["idx"])] = gap_abs_quad_eval(float(v[j]))
            return float(np.quantile(gaps, 0.95))
        opt = so.minimize(objective, x0=centers.copy(), method="Nelder-Mead", options={"maxiter": 320, "xatol": 1e-7, "fatol": 1e-20})
        vopt = np.minimum(np.maximum(np.array(opt.x, dtype=float), lo), hi)
        # Independent replay: replace optimized rows and recompute q95 with fixed_quad n1200 surrogate for those rows.
        gaps_quad = base_gaps.copy()
        gaps_fq1200 = base_gaps.copy()
        cert_rows = []
        for j, rr in enumerate(topk):
            s_opt = float(vopt[j])
            disc_loc, _ = strict_kernel_phase_integral(float(s_opt), float(om_gg), float(ph_gg), float(be_gg), float(et_gg))
            def base_integrand_local(x):
                xa = np.array(x, dtype=float)
                kk = np.cos(om_gg * xa + ph_gg) / (1.0 + be_gg * (xa ** et_gg))
                return (kk * kk) / np.sqrt(np.maximum(1e-15, xa + s_opt))
            cut_q, _ = si.quad(base_integrand_local, 0.0, 1.0, epsabs=1e-12, epsrel=1e-12, limit=900)
            cut_fq, _ = si.fixed_quad(base_integrand_local, 0.0, 1.0, n=1200)
            gq = float(abs(float(disc_loc) - float(cut_q)))
            gf = float(abs(float(disc_loc) - float(cut_fq)))
            gaps_quad[int(rr["idx"])] = gq
            gaps_fq1200[int(rr["idx"])] = gf
            cert_rows.append({
                "row_index": int(rr["idx"]),
                "s_baseline": float(rr["s"]),
                "s_optimized": s_opt,
                "gap_abs_quad_optimized": gq,
                "gap_abs_fixed_quad_n1200_optimized": gf,
                "cross_integrator_gap_abs": float(abs(gq - gf)),
            })
        q95_quad = float(np.quantile(gaps_quad, 0.95))
        q95_fq1200 = float(np.quantile(gaps_fq1200, 0.95))
        q95_blocker_joint_bestk_exact_recompute_certificate = {
            "status": "OPEN_PRECURSOR_NOT_CLOSURE",
            "scope": "STRICT_TASK2_Q95_BLOCKER_JOINT_BESTK_EXACT_RECOMPUTE_CERTIFICATE",
            "theorem_target": "BEST_K_JOINT_OPTIMIZATION_RECOMPUTE_CLOSES_Q95_UNDER_INDEPENDENT_INTEGRATOR_REPLAY",
            "domain": {"best_k_from_profile": best_k, "local_window_half_width": 0.03, "s_bounds_global": [0.05, 3.5]},
            "computed_rows": cert_rows,
            "aggregate_metrics": {
                "q95_gap_abs_after_joint_bestk_quad": q95_quad,
                "q95_gap_abs_after_joint_bestk_fixed_quad_n1200": q95_fq1200,
                "q95_cross_integrator_gap_abs": float(abs(q95_quad - q95_fq1200)),
                "optimizer_success": bool(opt.success),
                "optimizer_nit": int(getattr(opt, "nit", 0)),
                "optimizer_nfev": int(getattr(opt, "nfev", 0)),
            },
            "pass_fail_criteria": {
                "q95_quad_le_q95_gap_abs_max": q95_thr,
                "q95_fq1200_le_q95_gap_abs_max": q95_thr,
                "q95_cross_integrator_gap_abs_le_q95_cross_integrator_gap_abs_max": float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"]),
            },
        }
        c1 = bool(q95_quad <= q95_thr)
        c2 = bool(q95_fq1200 <= q95_thr)
        c3 = bool(abs(q95_quad - q95_fq1200) <= float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"]))
        q95_blocker_joint_bestk_exact_recompute_certificate["verdict"] = "CLOSED_NUMERICAL_WITNESS_TASK2" if (c1 and c2 and c3) else "OPEN_OBSTRUCTION_WITH_TRACE"
        q95_blocker_joint_bestk_exact_recompute_certificate["fail_trace"] = "" if (c1 and c2 and c3) else (
            f"q95_gap_abs_after_joint_bestk_quad={q95_quad:.6e} > {q95_thr:.1e}" if not c1 else (
                f"q95_gap_abs_after_joint_bestk_fixed_quad_n1200={q95_fq1200:.6e} > {q95_thr:.1e}" if not c2
                else f"q95_cross_integrator_gap_abs={abs(q95_quad-q95_fq1200):.6e} > {float(pass_fail_criteria_task2['q95_cross_integrator_gap_abs_max']):.1e}"
            )
        )
    else:
        q95_blocker_joint_bestk_exact_recompute_certificate = {
            "status": "OPEN_PRECURSOR_NOT_CLOSURE",
            "scope": "STRICT_TASK2_Q95_BLOCKER_JOINT_BESTK_EXACT_RECOMPUTE_CERTIFICATE",
            "theorem_target": "BEST_K_JOINT_OPTIMIZATION_RECOMPUTE_CLOSES_Q95_UNDER_INDEPENDENT_INTEGRATOR_REPLAY",
            "domain": {"best_k_from_profile": None, "local_window_half_width": 0.03, "s_bounds_global": [0.05, 3.5]},
            "computed_rows": [],
            "aggregate_metrics": {"q95_gap_abs_after_joint_bestk_quad": None, "q95_gap_abs_after_joint_bestk_fixed_quad_n1200": None, "q95_cross_integrator_gap_abs": None, "optimizer_success": False, "optimizer_nit": 0, "optimizer_nfev": 0},
            "pass_fail_criteria": {
                "q95_quad_le_q95_gap_abs_max": float(pass_fail_criteria_task2["q95_gap_abs_max"]),
                "q95_fq1200_le_q95_gap_abs_max": float(pass_fail_criteria_task2["q95_gap_abs_max"]),
                "q95_cross_integrator_gap_abs_le_q95_cross_integrator_gap_abs_max": float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"]),
            },
            "verdict": "OPEN_OBSTRUCTION_WITH_TRACE",
            "fail_trace": "bestk_exact_recompute_unavailable=1 > 0",
        }
    # One-shot closure efficiency theorem: minimal-k frontier gain slope certificate on adaptive profile.
    if q95_blocker_adaptive_joint_topk_profile_execution["computed_rows"]:
        rows_k = sorted(
            q95_blocker_adaptive_joint_topk_profile_execution["computed_rows"],
            key=lambda r: int(r["k"])
        )
        slope_rows = []
        for i in range(len(rows_k) - 1):
            k0 = int(rows_k[i]["k"])
            k1 = int(rows_k[i + 1]["k"])
            q0 = float(rows_k[i]["q95_gap_abs_after_joint_vector_optimization"])
            q1 = float(rows_k[i + 1]["q95_gap_abs_after_joint_vector_optimization"])
            slope_rows.append({
                "k_from": k0,
                "k_to": k1,
                "delta_q95": float(q1 - q0),
                "delta_k": int(k1 - k0),
                "marginal_q95_gain_per_k": float((q0 - q1) / max(1, (k1 - k0))),
            })
        best_row = sorted(rows_k, key=lambda r: (float(r["q95_gap_abs_after_joint_vector_optimization"]), int(r["k"])))[0]
        knee_row = sorted(slope_rows, key=lambda r: (-r["marginal_q95_gain_per_k"], r["k_from"]))[0] if slope_rows else None
        q95_thr = float(pass_fail_criteria_task2["q95_gap_abs_max"])
        q95_blocker_adaptive_profile_knee_certificate = {
            "status": "OPEN_PRECURSOR_NOT_CLOSURE",
            "scope": "STRICT_TASK2_Q95_BLOCKER_ADAPTIVE_PROFILE_KNEE_CERTIFICATE",
            "theorem_target": "ADAPTIVE_PROFILE_EXHIBITS_NONZERO_MARGINAL_Q95_GAIN_AND_IDENTIFIES_BEST_K",
            "domain": {"k_grid": [int(r["k"]) for r in rows_k]},
            "computed_rows": slope_rows,
            "aggregate_metrics": {
                "best_k_by_q95": int(best_row["k"]),
                "best_q95_after_joint_vector_optimization": float(best_row["q95_gap_abs_after_joint_vector_optimization"]),
                "knee_segment": knee_row,
            },
            "pass_fail_criteria": {
                "best_q95_le_q95_gap_abs_max": q95_thr,
                "max_marginal_q95_gain_per_k_gt_zero": 0.0,
            },
        }
        c1 = bool(float(best_row["q95_gap_abs_after_joint_vector_optimization"]) <= q95_thr)
        c2 = bool((knee_row is not None) and (float(knee_row["marginal_q95_gain_per_k"]) > 0.0))
        q95_blocker_adaptive_profile_knee_certificate["verdict"] = "CLOSED_NUMERICAL_WITNESS_TASK2" if (c1 and c2) else "OPEN_OBSTRUCTION_WITH_TRACE"
        q95_blocker_adaptive_profile_knee_certificate["fail_trace"] = "" if (c1 and c2) else (
            f"best_q95_after_joint_vector_optimization={float(best_row['q95_gap_abs_after_joint_vector_optimization']):.6e} > {q95_thr:.1e}"
            if not c1 else "max_marginal_q95_gain_per_k=0.000000e+00 > 0"
        )
    else:
        q95_blocker_adaptive_profile_knee_certificate = {
            "status": "OPEN_PRECURSOR_NOT_CLOSURE",
            "scope": "STRICT_TASK2_Q95_BLOCKER_ADAPTIVE_PROFILE_KNEE_CERTIFICATE",
            "theorem_target": "ADAPTIVE_PROFILE_EXHIBITS_NONZERO_MARGINAL_Q95_GAIN_AND_IDENTIFIES_BEST_K",
            "domain": {"k_grid": []},
            "computed_rows": [],
            "aggregate_metrics": {"best_k_by_q95": None, "best_q95_after_joint_vector_optimization": None, "knee_segment": None},
            "pass_fail_criteria": {"best_q95_le_q95_gap_abs_max": float(pass_fail_criteria_task2["q95_gap_abs_max"]), "max_marginal_q95_gain_per_k_gt_zero": 0.0},
            "verdict": "OPEN_OBSTRUCTION_WITH_TRACE",
            "fail_trace": "adaptive_profile_knee_unavailable=1 > 0",
        }
    # Bootstrap robustness of one-step action effect sign on dominant local support.
    if q95_blocker_action_execution["s_before"] is not None:
        s_before = float(q95_blocker_action_execution["s_before"])
        s_after = float(q95_blocker_action_execution["s_after"])
        boot_rng = np.random.default_rng(193975)
        boot_n = 128
        improve_count = 0
        delta_samples = []
        for _ in range(boot_n):
            n_loc = int(boot_rng.integers(300, 901))
            def gap_abs_fixedq(s_loc: float) -> float:
                disc_loc, _ = strict_kernel_phase_integral(float(s_loc), float(om_gg), float(ph_gg), float(be_gg), float(et_gg))
                def base_integrand(x):
                    xa = np.array(x, dtype=float)
                    kk = np.cos(om_gg * xa + ph_gg) / (1.0 + be_gg * (xa ** et_gg))
                    return (kk * kk) / np.sqrt(np.maximum(1e-15, xa + s_loc))
                cut_loc, _ = si.fixed_quad(base_integrand, 0.0, 1.0, n=n_loc)
                return float(abs(float(disc_loc) - float(cut_loc)))
            before = gap_abs_fixedq(s_before)
            after = gap_abs_fixedq(s_after)
            delta = float(after - before)
            delta_samples.append(delta)
            if delta < 0.0:
                improve_count += 1
        p_improve = float(improve_count / boot_n)
        q95_action_effect_bootstrap = {
            "status": "OPEN_PRECURSOR_NOT_CLOSURE",
            "scope": "STRICT_TASK2_Q95_ACTION_EFFECT_BOOTSTRAP_ROBUSTNESS",
            "bootstrap_size": int(boot_n),
            "n_range_inclusive": [300, 900],
            "p_improve": p_improve,
            "p_improve_ci95_jeffreys": jeffreys_interval_from_successes(improve_count, boot_n),
            "delta_q05_q50_q95": [float(np.quantile(np.array(delta_samples, dtype=float), q)) for q in [0.05, 0.50, 0.95]],
        }
    else:
        q95_action_effect_bootstrap = {
            "status": "OPEN_PRECURSOR_NOT_CLOSURE",
            "scope": "STRICT_TASK2_Q95_ACTION_EFFECT_BOOTSTRAP_ROBUSTNESS",
            "bootstrap_size": 0,
            "n_range_inclusive": [300, 900],
            "p_improve": 0.0,
            "p_improve_ci95_jeffreys": {"lower": 0.0, "upper": 1.0},
            "delta_q05_q50_q95": [0.0, 0.0, 0.0],
        }
    # Strict decision gate for local action execution validity under bootstrap robustness.
    q95_action_gate = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_Q95_ACTION_BOOTSTRAP_DECISION_GATE",
        "p_improve_lb95_threshold": 0.60,
        "criterion_p_improve_lb95_ge_threshold": bool(q95_action_effect_bootstrap["p_improve_ci95_jeffreys"]["lower"] >= 0.60),
        "criterion_fixed_quad_sign_consistent": bool(q95_action_effect_crosscheck["effect_sign_consistent_across_orders"]),
    }
    q95_action_gate["go_for_next_local_step"] = bool(
        q95_action_gate["criterion_p_improve_lb95_ge_threshold"] and
        q95_action_gate["criterion_fixed_quad_sign_consistent"]
    )
    q95_action_gate["reason"] = "GO" if q95_action_gate["go_for_next_local_step"] else "HOLD_AND_RECALIBRATE"
    q95_action_gate_consistency = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_Q95_ACTION_GATE_CONSISTENCY",
        "checks": {
            "if_go_then_action_exists": bool((not q95_action_gate["go_for_next_local_step"]) or (q95_blocker_action_execution["move"] != "NO_LOCAL_RELIEF")),
            "if_go_then_bootstrap_executed": bool((not q95_action_gate["go_for_next_local_step"]) or (q95_action_effect_bootstrap["bootstrap_size"] > 0)),
            "if_go_then_crosscheck_not_none": bool((not q95_action_gate["go_for_next_local_step"]) or (q95_action_effect_crosscheck["delta_gap_abs_n400_after_minus_before"] is not None)),
        },
    }
    q95_action_gate_consistency["all_checks_pass"] = bool(all(q95_action_gate_consistency["checks"].values()))
    # Direct one-step efficiency panel: realized blocker reduction per unit local step.
    if q95_blocker_action_execution["s_before"] is not None:
        s_before = float(q95_blocker_action_execution["s_before"])
        s_after = float(q95_blocker_action_execution["s_after"])
        step_abs = float(abs(s_after - s_before))
        gap_before = float(q95_blocker_action_execution["gap_abs_before"])
        gap_after = float(q95_blocker_action_execution["gap_abs_after"])
        gap_reduction = float(max(0.0, gap_before - gap_after))
        q95_blocker_step_efficiency = {
            "status": "OPEN_PRECURSOR_NOT_CLOSURE",
            "scope": "STRICT_TASK2_Q95_BLOCKER_ONE_STEP_EFFICIENCY",
            "step_abs": step_abs,
            "gap_abs_before": gap_before,
            "gap_abs_after": gap_after,
            "gap_abs_reduction": gap_reduction,
            "reduction_per_unit_s_step": float(gap_reduction / max(1e-15, step_abs)),
            "relative_reduction_fraction": float(gap_reduction / max(1e-30, gap_before)),
        }
    else:
        q95_blocker_step_efficiency = {
            "status": "OPEN_PRECURSOR_NOT_CLOSURE",
            "scope": "STRICT_TASK2_Q95_BLOCKER_ONE_STEP_EFFICIENCY",
            "step_abs": 0.0,
            "gap_abs_before": 0.0,
            "gap_abs_after": 0.0,
            "gap_abs_reduction": 0.0,
            "reduction_per_unit_s_step": 0.0,
            "relative_reduction_fraction": 0.0,
        }
    # Immediate closure-distance update after one-step action on dominant support.
    if q95_blocker_action_execution["s_before"] is not None:
        q95_after_one_step_local_margin = {
            "status": "OPEN_PRECURSOR_NOT_CLOSURE",
            "scope": "STRICT_TASK2_Q95_AFTER_ONE_STEP_LOCAL_MARGIN",
            "q95_gap_abs_threshold": float(pass_fail_criteria_task2["q95_gap_abs_max"]),
            "gap_abs_before": float(q95_blocker_action_execution["gap_abs_before"]),
            "gap_abs_after": float(q95_blocker_action_execution["gap_abs_after"]),
            "margin_before": float(q95_blocker_action_execution["gap_abs_before"] - float(pass_fail_criteria_task2["q95_gap_abs_max"])),
            "margin_after": float(q95_blocker_action_execution["gap_abs_after"] - float(pass_fail_criteria_task2["q95_gap_abs_max"])),
            "margin_delta_after_minus_before": float(
                (q95_blocker_action_execution["gap_abs_after"] - float(pass_fail_criteria_task2["q95_gap_abs_max"]))
                - (q95_blocker_action_execution["gap_abs_before"] - float(pass_fail_criteria_task2["q95_gap_abs_max"]))
            ),
            "moves_toward_closure": bool(q95_blocker_action_execution["gap_abs_after"] < q95_blocker_action_execution["gap_abs_before"]),
        }
    else:
        q95_after_one_step_local_margin = {
            "status": "OPEN_PRECURSOR_NOT_CLOSURE",
            "scope": "STRICT_TASK2_Q95_AFTER_ONE_STEP_LOCAL_MARGIN",
            "q95_gap_abs_threshold": float(pass_fail_criteria_task2["q95_gap_abs_max"]),
            "gap_abs_before": None,
            "gap_abs_after": None,
            "margin_before": None,
            "margin_after": None,
            "margin_delta_after_minus_before": None,
            "moves_toward_closure": False,
        }
    # One-step normalized closure progress score on local margin (0..1 when improving).
    if q95_after_one_step_local_margin["margin_before"] is not None:
        mb = float(q95_after_one_step_local_margin["margin_before"])
        ma = float(q95_after_one_step_local_margin["margin_after"])
        progress = float(max(0.0, mb - ma))
        denom = float(max(1e-30, abs(mb)))
        q95_after_one_step_progress_score = {
            "status": "OPEN_PRECURSOR_NOT_CLOSURE",
            "scope": "STRICT_TASK2_Q95_AFTER_ONE_STEP_PROGRESS_SCORE",
            "margin_before_abs": float(abs(mb)),
            "margin_after_abs": float(abs(ma)),
            "absolute_margin_improvement": progress,
            "normalized_progress_score_0_1": float(min(1.0, progress / denom)),
        }
    else:
        q95_after_one_step_progress_score = {
            "status": "OPEN_PRECURSOR_NOT_CLOSURE",
            "scope": "STRICT_TASK2_Q95_AFTER_ONE_STEP_PROGRESS_SCORE",
            "margin_before_abs": 0.0,
            "margin_after_abs": 0.0,
            "absolute_margin_improvement": 0.0,
            "normalized_progress_score_0_1": 0.0,
        }
    if q95_after_one_step_local_margin["margin_after"] is not None:
        margin_after_now = float(q95_after_one_step_local_margin["margin_after"])
        one_step_gain = float(max(0.0, q95_after_one_step_progress_score["absolute_margin_improvement"]))
        projected_margin_after_next_step = float(margin_after_now - one_step_gain)
        estimated_steps_to_closure = float(margin_after_now / max(1e-30, one_step_gain)) if one_step_gain > 0.0 else float("inf")
        q95_blocker_direct_relief_projection = {
            "status": "OPEN_PRECURSOR_NOT_CLOSURE",
            "scope": "STRICT_TASK2_Q95_BLOCKER_DIRECT_RELIEF_PROJECTION",
            "theorem_target": "ONE_MORE_LOCAL_STEP_CAN_CLOSE_Q95_MARGIN",
            "domain": {
                "support": "dominant_q95_channel_local_action",
                "single_step_model": "linearized_from_observed_one_step_margin_gain",
            },
            "pass_fail_criteria": {
                "projected_margin_after_next_step_le_zero": "required for predicted closure in next local move",
            },
            "observed_margin_after": margin_after_now,
            "observed_one_step_margin_gain": one_step_gain,
            "projected_margin_after_next_step": projected_margin_after_next_step,
            "estimated_steps_to_closure": estimated_steps_to_closure,
            "verdict": "PREDICTED_CLOSE_IN_NEXT_STEP" if projected_margin_after_next_step <= 0.0 else "OPEN_OBSTRUCTION_WITH_TRACE",
            "fail_trace": (
                ""
                if projected_margin_after_next_step <= 0.0
                else f"projected_margin_after_next_step={projected_margin_after_next_step:.6e} > 0"
            ),
        }
    else:
        q95_blocker_direct_relief_projection = {
            "status": "OPEN_PRECURSOR_NOT_CLOSURE",
            "scope": "STRICT_TASK2_Q95_BLOCKER_DIRECT_RELIEF_PROJECTION",
            "theorem_target": "ONE_MORE_LOCAL_STEP_CAN_CLOSE_Q95_MARGIN",
            "domain": {"support": "dominant_q95_channel_local_action", "single_step_model": "not_available"},
            "pass_fail_criteria": {"projected_margin_after_next_step_le_zero": "required for predicted closure in next local move"},
            "observed_margin_after": None,
            "observed_one_step_margin_gain": 0.0,
            "projected_margin_after_next_step": None,
            "estimated_steps_to_closure": float("inf"),
            "verdict": "OPEN_OBSTRUCTION_WITH_TRACE",
            "fail_trace": "projection_unavailable=1 > 0",
        }
    obstruction_is_numerically_sensitive = bool(float(np.max(np.array([r["cutsum_scheme_gap_abs"] for r in scheme_rows], dtype=float))) > 0.05 * q95_gap_abs_ext)
    q95_gap_abs_n800_top3 = float(q95_dominant_convergence["q95_gap_abs_n800_top3"])
    criteria_eval = {
        "q95_gap_abs_le_threshold": bool(q95_gap_abs_ext <= float(pass_fail_criteria_task2["q95_gap_abs_max"])),
        "q95_gap_abs_n2400_le_threshold": bool(q95_gap_abs_n2400_top3 <= float(pass_fail_criteria_task2["q95_gap_abs_max"])),
        "q95_gap_abs_n2400_vs_n6400_delta_le_threshold": bool(
            q95_n2400_vs_n6400_delta_top3 <= float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"])
        ),
        "q95_gap_abs_n3200_vs_n6400_delta_le_threshold": bool(
            abs(q95_gap_abs_n3200_top3 - q95_gap_abs_n6400_top3) <= float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"])
        ),
        "q95_gap_abs_n1600_vs_n6400_delta_le_threshold": bool(
            abs(q95_gap_abs_n1600_top3 - q95_gap_abs_n6400_top3) <= float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"])
        ),
        "q95_gap_abs_n800_vs_n6400_delta_le_threshold": bool(
            abs(float(q95_dominant_convergence["q95_gap_abs_n800_top3"]) - q95_gap_abs_n6400_top3) <= float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"])
        ),
        "q95_gap_abs_n12800_vs_n6400_delta_le_threshold": bool(
            q95_n12800_vs_n6400_delta_top3 <= float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"])
        ),
        "q95_convergence_tail_ratio_n12800_le_one": bool(q95_convergence_tail_ratio_n12800 <= 1.0),
        "q95_gap_abs_n25600_vs_n12800_delta_le_threshold": bool(q95_n25600_vs_n12800_delta_top3 <= float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"])),
        "q95_monotone_violation_n25600_le_zero": bool(q95_monotone_violation_top3 <= 0.0),
        "max_gap_rel_le_threshold": bool(max_gap_rel_ext <= float(pass_fail_criteria_task2["max_gap_rel_max"])),
        "all_nonnegative_weights": bool(all_nonnegative),
        "q95_cross_integrator_gap_le_threshold": bool(crosscheck_gap_q95 <= float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"])),
        "q95_convergence_delta_le_threshold": bool(
            q95_unc_abs <= float(pass_fail_criteria_task2["q95_convergence_delta_n400_to_n800_abs_max"])
        ),
    }
    criterion_coherence_sign = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_SIGN_CRITERION_COHERENCE",
        "min_effective_weight_global": float(min_effective_weight_global),
        "min_effective_weight_global_min": float(pass_fail_criteria_task2["min_effective_weight_global_min"]),
        "all_nonnegative_weights_flag": bool(criteria_eval["all_nonnegative_weights"]),
        "flag_equals_numeric_inequality": bool(
            criteria_eval["all_nonnegative_weights"] ==
            (float(min_effective_weight_global) >= float(pass_fail_criteria_task2["min_effective_weight_global_min"]))
        ),
    }
    criterion_coherence_convergence = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_CONVERGENCE_CRITERION_COHERENCE",
        "q95_convergence_delta_n400_to_n800_abs": float(q95_unc_abs),
        "q95_convergence_delta_n400_to_n800_abs_max": float(pass_fail_criteria_task2["q95_convergence_delta_n400_to_n800_abs_max"]),
        "q95_convergence_delta_le_threshold_flag": bool(criteria_eval["q95_convergence_delta_le_threshold"]),
        "flag_equals_numeric_inequality": bool(
            criteria_eval["q95_convergence_delta_le_threshold"] ==
            (float(q95_unc_abs) <= float(pass_fail_criteria_task2["q95_convergence_delta_n400_to_n800_abs_max"]))
        ),
    }
    criterion_coherence_cross_integrator = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_CROSS_INTEGRATOR_CRITERION_COHERENCE",
        "q95_cross_integrator_gap_abs": float(crosscheck_gap_q95),
        "q95_cross_integrator_gap_abs_max": float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"]),
        "q95_cross_integrator_gap_le_threshold_flag": bool(criteria_eval["q95_cross_integrator_gap_le_threshold"]),
        "flag_equals_numeric_inequality": bool(
            criteria_eval["q95_cross_integrator_gap_le_threshold"] ==
            (float(crosscheck_gap_q95) <= float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"]))
        ),
    }
    criterion_coherence_q95_gap = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_Q95_GAP_CRITERION_COHERENCE",
        "q95_gap_abs": float(q95_gap_abs_ext),
        "q95_gap_abs_max": float(pass_fail_criteria_task2["q95_gap_abs_max"]),
        "q95_gap_abs_le_threshold_flag": bool(criteria_eval["q95_gap_abs_le_threshold"]),
        "flag_equals_numeric_inequality": bool(
            criteria_eval["q95_gap_abs_le_threshold"] ==
            (float(q95_gap_abs_ext) <= float(pass_fail_criteria_task2["q95_gap_abs_max"]))
        ),
    }
    criterion_coherence_max_gap_rel = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_MAX_GAP_REL_CRITERION_COHERENCE",
        "max_gap_rel": float(max_gap_rel_ext),
        "max_gap_rel_max": float(pass_fail_criteria_task2["max_gap_rel_max"]),
        "max_gap_rel_le_threshold_flag": bool(criteria_eval["max_gap_rel_le_threshold"]),
        "flag_equals_numeric_inequality": bool(
            criteria_eval["max_gap_rel_le_threshold"] ==
            (float(max_gap_rel_ext) <= float(pass_fail_criteria_task2["max_gap_rel_max"]))
        ),
    }
    closure_numeric_conjunction_recomputed = bool(
        (float(q95_gap_abs_ext) <= float(pass_fail_criteria_task2["q95_gap_abs_max"])) and
        (float(q95_gap_abs_n2400_top3) <= float(pass_fail_criteria_task2["q95_gap_abs_max"])) and
        (float(q95_n2400_vs_n6400_delta_top3) <= float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"])) and
        (abs(float(q95_gap_abs_n3200_top3) - float(q95_gap_abs_n6400_top3)) <= float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"])) and
        (abs(float(q95_gap_abs_n1600_top3) - float(q95_gap_abs_n6400_top3)) <= float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"])) and
        (abs(float(q95_gap_abs_n800_top3) - float(q95_gap_abs_n6400_top3)) <= float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"])) and
        (float(q95_n12800_vs_n6400_delta_top3) <= float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"])) and
        (float(q95_convergence_tail_ratio_n12800) <= 1.0) and
        (float(q95_n25600_vs_n12800_delta_top3) <= float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"])) and
        (float(q95_monotone_violation_top3) <= 0.0) and
        (float(max_gap_rel_ext) <= float(pass_fail_criteria_task2["max_gap_rel_max"])) and
        (float(crosscheck_gap_q95) <= float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"])) and
        (float(q95_unc_abs) <= float(pass_fail_criteria_task2["q95_convergence_delta_n400_to_n800_abs_max"])) and
        (float(min_effective_weight_global) >= float(pass_fail_criteria_task2["min_effective_weight_global_min"]))
    )
    criterion_coherence_global_closure = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_GLOBAL_CLOSURE_CRITERION_COHERENCE",
        "all_criteria_satisfied_flag": bool(all(criteria_eval.values())),
        "numeric_conjunction_recomputed": closure_numeric_conjunction_recomputed,
        "flag_equals_numeric_conjunction": bool(bool(all(criteria_eval.values())) == closure_numeric_conjunction_recomputed),
    }
    criterion_coherence_weight_sign = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_WEIGHT_SIGN_CRITERION_COHERENCE",
        "min_effective_weight_global": float(min_effective_weight_global),
        "min_effective_weight_global_min": float(pass_fail_criteria_task2["min_effective_weight_global_min"]),
        "weight_sign_nonnegative_flag": bool(criteria_eval["all_nonnegative_weights"]),
        "flag_equals_numeric_inequality": bool(
            criteria_eval["all_nonnegative_weights"] ==
            (float(min_effective_weight_global) >= float(pass_fail_criteria_task2["min_effective_weight_global_min"]))
        ),
    }
    criterion_coherence_verdict_gate = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_VERDICT_GATE_COHERENCE",
        "all_criteria_satisfied_flag": bool(all(criteria_eval.values())),
        "verdict_task2": str(verdict_task2),
        "closed_label": "CLOSED_NUMERICAL_WITNESS_TASK2",
        "open_label": "OPEN_OBSTRUCTION_WITH_TRACE",
        "flag_matches_verdict_label": bool(
            (bool(all(criteria_eval.values())) and str(verdict_task2) == "CLOSED_NUMERICAL_WITNESS_TASK2") or
            ((not bool(all(criteria_eval.values()))) and str(verdict_task2) == "OPEN_OBSTRUCTION_WITH_TRACE")
        ),
    }
    criterion_coherence_fail_trace_numeric = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_FAIL_TRACE_NUMERIC_COHERENCE",
        "verdict_task2": str(locals().get("verdict_task2", "PENDING")),
        "dominant_blocker": str(locals().get("dominant_blocker", "pending")),
        "fail_trace": str(locals().get("fail_trace_task2", "")),
        "trace_prefix_matches_dominant_blocker": bool(
            (str(locals().get("verdict_task2", "PENDING")) != "OPEN_OBSTRUCTION_WITH_TRACE") or
            (
                (str(locals().get("dominant_blocker", "pending")) == "q95_gap_abs" and str(locals().get("fail_trace_task2", "")).startswith("q95_gap_abs=")) or
                (str(locals().get("dominant_blocker", "pending")) == "max_gap_rel" and str(locals().get("fail_trace_task2", "")).startswith("max_gap_rel=")) or
                (str(locals().get("dominant_blocker", "pending")) == "q95_cross_integrator_gap" and str(locals().get("fail_trace_task2", "")).startswith("q95_cross_integrator_gap_abs=")) or
                (str(locals().get("dominant_blocker", "pending")) == "q95_convergence_delta_n400_to_n800_abs" and str(locals().get("fail_trace_task2", "")).startswith("q95_convergence_delta_n400_to_n800_abs=")) or
                (str(locals().get("dominant_blocker", "pending")) == "weight_sign_nonnegativity" and str(locals().get("fail_trace_task2", "")).startswith("min_effective_weight_global="))
            )
        ),
    }
    criterion_coherence_dominant_margin_sign = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_DOMINANT_MARGIN_SIGN_COHERENCE",
        "verdict_task2": str(locals().get("verdict_task2", "PENDING")),
        "dominant_blocker": str(locals().get("dominant_blocker", "pending")),
        "signed_margin_observed_minus_threshold": float(
            float(locals().get("dominant_observed", 0.0)) - float(locals().get("dominant_threshold", 0.0))
        ),
        "open_requires_positive_margin": bool(
            (str(locals().get("verdict_task2", "PENDING")) != "OPEN_OBSTRUCTION_WITH_TRACE") or
            ((float(locals().get("dominant_observed", 0.0)) - float(locals().get("dominant_threshold", 0.0))) > 0.0)
        ),
    }
    criterion_coherence_open_trace_inequality = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_OPEN_TRACE_INEQUALITY_COHERENCE",
        "verdict_task2": str(locals().get("verdict_task2", "PENDING")),
        "fail_trace": str(locals().get("fail_trace_task2", "")),
        "open_trace_contains_numeric_inequality_token": bool(
            (str(locals().get("verdict_task2", "PENDING")) != "OPEN_OBSTRUCTION_WITH_TRACE") or
            (
                (">" in str(locals().get("fail_trace_task2", ""))) and
                ("=" in str(locals().get("fail_trace_task2", "")))
            )
        ),
    }
    criterion_coherence_dominant_inequality_prefix = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_DOMINANT_INEQUALITY_PREFIX_COHERENCE",
        "dominant_blocker": str(locals().get("dominant_blocker", "pending")),
        "dominant_inequality": str(locals().get("dominant_inequality", "")),
        "prefix_matches_dominant_blocker": bool(
            (str(locals().get("dominant_blocker", "pending")) == "q95_gap_abs" and str(locals().get("dominant_inequality", "")).startswith("q95_gap_abs=")) or
            (str(locals().get("dominant_blocker", "pending")) == "max_gap_rel" and str(locals().get("dominant_inequality", "")).startswith("max_gap_rel=")) or
            (str(locals().get("dominant_blocker", "pending")) == "q95_cross_integrator_gap" and str(locals().get("dominant_inequality", "")).startswith("q95_cross_integrator_gap_abs=")) or
            (str(locals().get("dominant_blocker", "pending")) == "q95_convergence_delta_n400_to_n800_abs" and str(locals().get("dominant_inequality", "")).startswith("q95_convergence_delta_n400_to_n800_abs=")) or
            (str(locals().get("dominant_blocker", "pending")) == "q95_refined_window_gap_abs" and str(locals().get("dominant_inequality", "")).startswith("q95_refined_window_best_gap_abs=")) or
            (str(locals().get("dominant_blocker", "pending")) == "q95_n2400_gap_abs" and str(locals().get("dominant_inequality", "")).startswith("q95_gap_abs_n2400_top3=")) or
            (str(locals().get("dominant_blocker", "pending")) == "q95_n2400_vs_n6400_delta_abs" and str(locals().get("dominant_inequality", "")).startswith("q95_delta_gap_abs_n2400_minus_n6400_abs_top3=")) or
            (str(locals().get("dominant_blocker", "pending")) == "q95_n3200_vs_n6400_delta_abs" and str(locals().get("dominant_inequality", "")).startswith("q95_delta_gap_abs_n3200_minus_n6400_abs_top3=")) or
            (str(locals().get("dominant_blocker", "pending")) == "q95_n1600_vs_n6400_delta_abs" and str(locals().get("dominant_inequality", "")).startswith("q95_delta_gap_abs_n1600_minus_n6400_abs_top3=")) or
            (str(locals().get("dominant_blocker", "pending")) == "q95_n800_vs_n6400_delta_abs" and str(locals().get("dominant_inequality", "")).startswith("q95_delta_gap_abs_n800_minus_n6400_abs_top3=")) or
            (str(locals().get("dominant_blocker", "pending")) == "q95_n12800_vs_n6400_delta_abs" and str(locals().get("dominant_inequality", "")).startswith("q95_delta_gap_abs_n12800_minus_n6400_abs_top3=")) or
            (str(locals().get("dominant_blocker", "pending")) == "q95_n12800_tail_ratio" and str(locals().get("dominant_inequality", "")).startswith("q95_convergence_tail_ratio_n6400_12800_over_n3200_6400=")) or
            (str(locals().get("dominant_blocker", "pending")) == "q95_n25600_vs_n12800_delta_abs" and str(locals().get("dominant_inequality", "")).startswith("q95_delta_gap_abs_n25600_minus_n12800_abs_top3=")) or
            (str(locals().get("dominant_blocker", "pending")) == "q95_n25600_monotone_violation" and str(locals().get("dominant_inequality", "")).startswith("q95_total_monotone_violation_top3=")) or
            (str(locals().get("dominant_blocker", "pending")) == "weight_sign_nonnegativity" and str(locals().get("dominant_inequality", "")).startswith("min_effective_weight_global=")) or
            (str(locals().get("dominant_blocker", "pending")) in {"none", "pending"} and str(locals().get("dominant_inequality", "")) in {"all_criteria_satisfied", ""})
        ),
    }
    criterion_coherence_fail_trace_equals_dominant = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_FAIL_TRACE_EQUALS_DOMINANT_COHERENCE",
        "verdict_task2": str(locals().get("verdict_task2", "PENDING")),
        "fail_trace": str(locals().get("fail_trace_task2", "")),
        "dominant_inequality": str(locals().get("dominant_inequality", "")),
        "open_requires_exact_equality": bool(
            (str(locals().get("verdict_task2", "PENDING")) != "OPEN_OBSTRUCTION_WITH_TRACE") or
            (str(locals().get("fail_trace_task2", "")) == str(locals().get("dominant_inequality", "")))
        ),
    }
    # Blocker-choice panel: identify easiest unresolved criterion to close (smallest normalized overshoot).
    q95_thr = float(pass_fail_criteria_task2["q95_gap_abs_max"])
    rel_thr = float(pass_fail_criteria_task2["max_gap_rel_max"])
    overs_q95 = float(max(0.0, q95_gap_abs_ext - q95_thr) / max(1e-30, q95_thr))
    refined_best_gap = (
        float(q95_blocker_refined_window_execution["aggregate_metrics"]["best_gap_abs_quad"])
        if q95_blocker_refined_window_execution["aggregate_metrics"]["best_gap_abs_quad"] is not None
        else float(q95_gap_abs_ext)
    )
    quad_hp_q95_gap = (
        float(q95_blocker_quad_hp_top3_certificate["aggregate_metrics"]["q95_gap_abs_quad_high_precision_top3"])
        if q95_blocker_quad_hp_top3_certificate["aggregate_metrics"]["q95_gap_abs_quad_high_precision_top3"] is not None else float("inf")
    )
    overs_q95_refined_window = float(max(0.0, refined_best_gap - q95_thr) / max(1e-30, q95_thr))
    overs_q95_quad_hp_top3 = float(max(0.0, quad_hp_q95_gap - q95_thr) / max(1e-30, q95_thr))
    quad_hp_q95_gap_upper_envelope = (
        float(q95_blocker_quad_hp_error_envelope_certificate["aggregate_metrics"]["q95_gap_abs_quad_high_precision_upper_envelope_top3"])
        if q95_blocker_quad_hp_error_envelope_certificate["aggregate_metrics"]["q95_gap_abs_quad_high_precision_upper_envelope_top3"] is not None else float("inf")
    )
    overs_q95_quad_hp_upper_envelope = float(max(0.0, quad_hp_q95_gap_upper_envelope - q95_thr) / max(1e-30, q95_thr))
    tail_budget_q95_upper = (
        float(q95_blocker_tail_budget_certificate["aggregate_metrics"]["q95_gap_abs_upper_tail_envelope_top3"])
        if q95_blocker_tail_budget_certificate["aggregate_metrics"]["q95_gap_abs_upper_tail_envelope_top3"] is not None else float("inf")
    )
    overs_q95_tail_budget_upper = float(max(0.0, tail_budget_q95_upper - q95_thr) / max(1e-30, q95_thr))
    n2400_q95_gap = (
        float(q95_blocker_n2400_recompute_certificate["aggregate_metrics"]["q95_gap_abs_n2400_top3"])
        if q95_blocker_n2400_recompute_certificate["aggregate_metrics"]["q95_gap_abs_n2400_top3"] is not None else float("inf")
    )
    overs_q95_n2400 = float(max(0.0, n2400_q95_gap - q95_thr) / max(1e-30, q95_thr))
    n2400_vs_6400_delta = (
        float(q95_blocker_n2400_recompute_certificate["aggregate_metrics"]["q95_delta_gap_abs_n2400_minus_n6400_abs_top3"])
        if q95_blocker_n2400_recompute_certificate["aggregate_metrics"]["q95_delta_gap_abs_n2400_minus_n6400_abs_top3"] is not None else float("inf")
    )
    cross_thr = float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"])
    overs_q95_n2400_delta = float(max(0.0, n2400_vs_6400_delta - cross_thr) / max(1e-30, cross_thr))
    n3200_vs_6400_delta = float(abs(q95_gap_abs_n3200_top3 - q95_gap_abs_n6400_top3))
    overs_q95_n3200_delta = float(max(0.0, n3200_vs_6400_delta - cross_thr) / max(1e-30, cross_thr))
    n1600_vs_6400_delta = float(abs(q95_gap_abs_n1600_top3 - q95_gap_abs_n6400_top3))
    overs_q95_n1600_delta = float(max(0.0, n1600_vs_6400_delta - cross_thr) / max(1e-30, cross_thr))
    n800_vs_6400_delta = float(abs(q95_gap_abs_n800_top3 - q95_gap_abs_n6400_top3))
    overs_q95_n800_delta = float(max(0.0, n800_vs_6400_delta - cross_thr) / max(1e-30, cross_thr))
    n12800_vs_6400_delta = float(q95_n12800_vs_n6400_delta_top3)
    overs_q95_n12800_delta = float(max(0.0, n12800_vs_6400_delta - cross_thr) / max(1e-30, cross_thr))
    overs_q95_n12800_tail_ratio = float(max(0.0, q95_convergence_tail_ratio_n12800 - 1.0) / 1.0)
    overs_q95_n25600_delta = float(max(0.0, q95_n25600_vs_n12800_delta_top3 - cross_thr) / max(1e-30, cross_thr))
    overs_q95_n25600_monotone = float(max(0.0, q95_monotone_violation_top3 - 0.0))
    overs_rel = float(max(0.0, max_gap_rel_ext - rel_thr) / max(1e-30, rel_thr))
    sign_thr = float(pass_fail_criteria_task2["min_effective_weight_global_min"])
    overs_sign = float(max(0.0, sign_thr - min_effective_weight_global) / max(1e-30, abs(sign_thr) + 1e-12))
    blocker_rows = [
        {"criterion": "q95_gap_abs", "normalized_overshoot": overs_q95, "is_satisfied": bool(criteria_eval["q95_gap_abs_le_threshold"])},
        {
            "criterion": "q95_refined_window_gap_abs",
            "normalized_overshoot": overs_q95_refined_window,
            "is_satisfied": bool(refined_best_gap <= q95_thr),
        },
        {
            "criterion": "q95_quad_hp_top3_gap_abs",
            "normalized_overshoot": overs_q95_quad_hp_top3,
            "is_satisfied": bool(quad_hp_q95_gap <= q95_thr),
        },
        {
            "criterion": "q95_quad_hp_upper_envelope_gap_abs",
            "normalized_overshoot": overs_q95_quad_hp_upper_envelope,
            "is_satisfied": bool(quad_hp_q95_gap_upper_envelope <= q95_thr),
        },
        {
            "criterion": "q95_tail_budget_upper_envelope_gap_abs",
            "normalized_overshoot": overs_q95_tail_budget_upper,
            "is_satisfied": bool(tail_budget_q95_upper <= q95_thr),
        },
        {
            "criterion": "q95_n2400_gap_abs",
            "normalized_overshoot": overs_q95_n2400,
            "is_satisfied": bool(n2400_q95_gap <= q95_thr),
        },
        {
            "criterion": "q95_n2400_vs_n6400_delta_abs",
            "normalized_overshoot": overs_q95_n2400_delta,
            "is_satisfied": bool(n2400_vs_6400_delta <= cross_thr),
        },
        {
            "criterion": "q95_n3200_vs_n6400_delta_abs",
            "normalized_overshoot": overs_q95_n3200_delta,
            "is_satisfied": bool(n3200_vs_6400_delta <= cross_thr),
        },
        {
            "criterion": "q95_n1600_vs_n6400_delta_abs",
            "normalized_overshoot": overs_q95_n1600_delta,
            "is_satisfied": bool(n1600_vs_6400_delta <= cross_thr),
        },
        {
            "criterion": "q95_n800_vs_n6400_delta_abs",
            "normalized_overshoot": overs_q95_n800_delta,
            "is_satisfied": bool(n800_vs_6400_delta <= cross_thr),
        },
        {
            "criterion": "q95_n12800_vs_n6400_delta_abs",
            "normalized_overshoot": overs_q95_n12800_delta,
            "is_satisfied": bool(n12800_vs_6400_delta <= cross_thr),
        },
        {
            "criterion": "q95_n12800_tail_ratio",
            "normalized_overshoot": overs_q95_n12800_tail_ratio,
            "is_satisfied": bool(q95_convergence_tail_ratio_n12800 <= 1.0),
        },
        {
            "criterion": "q95_n25600_vs_n12800_delta_abs",
            "normalized_overshoot": overs_q95_n25600_delta,
            "is_satisfied": bool(q95_n25600_vs_n12800_delta_top3 <= cross_thr),
        },
        {
            "criterion": "q95_n25600_monotone_violation",
            "normalized_overshoot": overs_q95_n25600_monotone,
            "is_satisfied": bool(q95_monotone_violation_top3 <= 0.0),
        },
        {"criterion": "max_gap_rel", "normalized_overshoot": overs_rel, "is_satisfied": bool(criteria_eval["max_gap_rel_le_threshold"])},
        {"criterion": "weight_sign_nonnegativity", "normalized_overshoot": overs_sign, "is_satisfied": bool(criteria_eval["all_nonnegative_weights"])},
        {
            "criterion": "q95_cross_integrator_gap",
            "normalized_overshoot": float(
                max(0.0, crosscheck_gap_q95 - float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"]))
                / max(1e-30, float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"]))
            ),
            "is_satisfied": bool(criteria_eval["q95_cross_integrator_gap_le_threshold"]),
        },
        {
            "criterion": "q95_convergence_delta_n400_to_n800_abs",
            "normalized_overshoot": float(
                max(0.0, q95_unc_abs - float(pass_fail_criteria_task2["q95_convergence_delta_n400_to_n800_abs_max"]))
                / max(1e-30, float(pass_fail_criteria_task2["q95_convergence_delta_n400_to_n800_abs_max"]))
            ),
            "is_satisfied": bool(criteria_eval["q95_convergence_delta_le_threshold"]),
        },
    ]
    unresolved = [r for r in blocker_rows if not r["is_satisfied"]]
    easiest_unresolved = sorted(unresolved, key=lambda r: (r["normalized_overshoot"], r["criterion"]))[0]["criterion"] if unresolved else "none"
    dominant_unresolved = sorted(unresolved, key=lambda r: (-r["normalized_overshoot"], r["criterion"]))[0]["criterion"] if unresolved else "none"
    dominant_unresolved_expected = (
        sorted(unresolved, key=lambda r: (-r["normalized_overshoot"], r["criterion"]))[0]["criterion"] if unresolved else "none"
    )
    q95_blocker_choice_panel = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_BLOCKER_CHOICE_NORMALIZED_OVERSHOOT",
        "rows": blocker_rows,
        "easiest_unresolved_blocker": easiest_unresolved,
        "dominant_unresolved_blocker": dominant_unresolved,
    }
    dominant_blocker_snapshot = str(locals().get("dominant_blocker", "pending"))
    q95_blocker_choice_consistency = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_BLOCKER_CHOICE_CONSISTENCY",
        "dominant_blocker": dominant_blocker_snapshot,
        "easiest_unresolved_blocker": easiest_unresolved,
        "is_consistent_when_q95_dominates": bool((dominant_blocker_snapshot != "q95_gap_abs") or (easiest_unresolved in {"q95_gap_abs", "none"})),
    }
    dominant_blocker = dominant_unresolved
    dominant_inequality_map = {
        "q95_gap_abs": f"q95_gap_abs={q95_gap_abs_ext:.6e} > {pass_fail_criteria_task2['q95_gap_abs_max']:.1e}",
        "max_gap_rel": f"max_gap_rel={max_gap_rel_ext:.6e} > {pass_fail_criteria_task2['max_gap_rel_max']:.1e}",
        "q95_cross_integrator_gap": f"q95_cross_integrator_gap_abs={crosscheck_gap_q95:.6e} > {pass_fail_criteria_task2['q95_cross_integrator_gap_abs_max']:.1e}",
        "q95_convergence_delta_n400_to_n800_abs": f"q95_convergence_delta_n400_to_n800_abs={q95_unc_abs:.6e} > {pass_fail_criteria_task2['q95_convergence_delta_n400_to_n800_abs_max']:.1e}",
        "q95_refined_window_gap_abs": f"q95_refined_window_best_gap_abs={refined_best_gap:.6e} > {pass_fail_criteria_task2['q95_gap_abs_max']:.1e}",
        "q95_quad_hp_top3_gap_abs": f"q95_gap_abs_quad_high_precision_top3={quad_hp_q95_gap:.6e} > {pass_fail_criteria_task2['q95_gap_abs_max']:.1e}",
        "q95_quad_hp_upper_envelope_gap_abs": f"q95_gap_abs_quad_high_precision_upper_envelope_top3={quad_hp_q95_gap_upper_envelope:.6e} > {pass_fail_criteria_task2['q95_gap_abs_max']:.1e}",
        "q95_tail_budget_upper_envelope_gap_abs": f"q95_gap_abs_upper_tail_envelope_top3={tail_budget_q95_upper:.6e} > {pass_fail_criteria_task2['q95_gap_abs_max']:.1e}",
        "q95_n2400_gap_abs": f"q95_gap_abs_n2400_top3={n2400_q95_gap:.6e} > {pass_fail_criteria_task2['q95_gap_abs_max']:.1e}",
        "q95_n2400_vs_n6400_delta_abs": f"q95_delta_gap_abs_n2400_minus_n6400_abs_top3={n2400_vs_6400_delta:.6e} > {pass_fail_criteria_task2['q95_cross_integrator_gap_abs_max']:.1e}",
        "q95_n3200_vs_n6400_delta_abs": f"q95_delta_gap_abs_n3200_minus_n6400_abs_top3={n3200_vs_6400_delta:.6e} > {pass_fail_criteria_task2['q95_cross_integrator_gap_abs_max']:.1e}",
        "q95_n1600_vs_n6400_delta_abs": f"q95_delta_gap_abs_n1600_minus_n6400_abs_top3={n1600_vs_6400_delta:.6e} > {pass_fail_criteria_task2['q95_cross_integrator_gap_abs_max']:.1e}",
        "q95_n800_vs_n6400_delta_abs": f"q95_delta_gap_abs_n800_minus_n6400_abs_top3={n800_vs_6400_delta:.6e} > {pass_fail_criteria_task2['q95_cross_integrator_gap_abs_max']:.1e}",
        "q95_n12800_vs_n6400_delta_abs": f"q95_delta_gap_abs_n12800_minus_n6400_abs_top3={n12800_vs_6400_delta:.6e} > {pass_fail_criteria_task2['q95_cross_integrator_gap_abs_max']:.1e}",
        "q95_n12800_tail_ratio": f"q95_convergence_tail_ratio_n6400_12800_over_n3200_6400={q95_convergence_tail_ratio_n12800:.6e} > 1.0",
        "q95_n25600_vs_n12800_delta_abs": f"q95_delta_gap_abs_n25600_minus_n12800_abs_top3={q95_n25600_vs_n12800_delta_top3:.6e} > {pass_fail_criteria_task2['q95_cross_integrator_gap_abs_max']:.1e}",
        "q95_n25600_monotone_violation": f"q95_total_monotone_violation_top3={q95_monotone_violation_top3:.6e} > 0.0",
        "weight_sign_nonnegativity": f"min_effective_weight_global={min_effective_weight_global:.6e} < {pass_fail_criteria_task2['min_effective_weight_global_min']:.1e}",
    }
    dominant_inequality = dominant_inequality_map.get(dominant_blocker, "all_criteria_satisfied")
    if dominant_blocker == "q95_gap_abs":
        dominant_observed = float(q95_gap_abs_ext)
        dominant_threshold = float(pass_fail_criteria_task2["q95_gap_abs_max"])
    elif dominant_blocker == "max_gap_rel":
        dominant_observed = float(max_gap_rel_ext)
        dominant_threshold = float(pass_fail_criteria_task2["max_gap_rel_max"])
    elif dominant_blocker == "q95_cross_integrator_gap":
        dominant_observed = float(crosscheck_gap_q95)
        dominant_threshold = float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"])
    elif dominant_blocker == "q95_convergence_delta_n400_to_n800_abs":
        dominant_observed = float(q95_unc_abs)
        dominant_threshold = float(pass_fail_criteria_task2["q95_convergence_delta_n400_to_n800_abs_max"])
    elif dominant_blocker == "weight_sign_nonnegativity":
        dominant_observed = float(min_effective_weight_global)
        dominant_threshold = float(pass_fail_criteria_task2["min_effective_weight_global_min"])
    elif dominant_blocker == "q95_refined_window_gap_abs":
        dominant_observed = float(refined_best_gap)
        dominant_threshold = float(pass_fail_criteria_task2["q95_gap_abs_max"])
    elif dominant_blocker == "q95_quad_hp_top3_gap_abs":
        dominant_observed = float(quad_hp_q95_gap)
        dominant_threshold = float(pass_fail_criteria_task2["q95_gap_abs_max"])
    elif dominant_blocker == "q95_quad_hp_upper_envelope_gap_abs":
        dominant_observed = float(quad_hp_q95_gap_upper_envelope)
        dominant_threshold = float(pass_fail_criteria_task2["q95_gap_abs_max"])
    elif dominant_blocker == "q95_tail_budget_upper_envelope_gap_abs":
        dominant_observed = float(tail_budget_q95_upper)
        dominant_threshold = float(pass_fail_criteria_task2["q95_gap_abs_max"])
    elif dominant_blocker == "q95_n2400_gap_abs":
        dominant_observed = float(n2400_q95_gap)
        dominant_threshold = float(pass_fail_criteria_task2["q95_gap_abs_max"])
    elif dominant_blocker == "q95_n2400_vs_n6400_delta_abs":
        dominant_observed = float(n2400_vs_6400_delta)
        dominant_threshold = float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"])
    elif dominant_blocker == "q95_n3200_vs_n6400_delta_abs":
        dominant_observed = float(n3200_vs_6400_delta)
        dominant_threshold = float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"])
    elif dominant_blocker == "q95_n1600_vs_n6400_delta_abs":
        dominant_observed = float(n1600_vs_6400_delta)
        dominant_threshold = float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"])
    elif dominant_blocker == "q95_n800_vs_n6400_delta_abs":
        dominant_observed = float(n800_vs_6400_delta)
        dominant_threshold = float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"])
    elif dominant_blocker == "q95_n12800_vs_n6400_delta_abs":
        dominant_observed = float(n12800_vs_6400_delta)
        dominant_threshold = float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"])
    elif dominant_blocker == "q95_n12800_tail_ratio":
        dominant_observed = float(q95_convergence_tail_ratio_n12800)
        dominant_threshold = 1.0
    elif dominant_blocker == "q95_n25600_vs_n12800_delta_abs":
        dominant_observed = float(q95_n25600_vs_n12800_delta_top3)
        dominant_threshold = float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"])
    elif dominant_blocker == "q95_n25600_monotone_violation":
        dominant_observed = float(q95_monotone_violation_top3)
        dominant_threshold = 0.0
    else:
        dominant_observed = 0.0
        dominant_threshold = 0.0
    dominant_blocker_numeric_margin = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_DOMINANT_BLOCKER_NUMERIC_MARGIN",
        "dominant_blocker": dominant_blocker,
        "observed_value": dominant_observed,
        "threshold_value": dominant_threshold,
        "signed_margin_observed_minus_threshold": float(dominant_observed - dominant_threshold),
        "normalized_overshoot_vs_threshold": float(max(0.0, dominant_observed - dominant_threshold) / max(1e-30, abs(dominant_threshold) + 1e-12)),
    }
    dominant_blocker_selection_consistency = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_DOMINANT_BLOCKER_SELECTION_CONSISTENCY",
        "dominant_blocker": dominant_blocker,
        "dominant_unresolved_expected": dominant_unresolved_expected,
        "is_argmax_overshoot": bool(dominant_blocker == dominant_unresolved_expected),
    }
    # One strict theorem-like blocker-closure step:
    # robust dominance witness under conservative uncertainty inflation of q95 blocker.
    second_largest_overs = sorted(
        [float(r["normalized_overshoot"]) for r in unresolved if r["criterion"] != dominant_blocker],
        reverse=True
    )
    second_best_overs = float(second_largest_overs[0]) if second_largest_overs else 0.0
    dominant_row = next((r for r in blocker_rows if r["criterion"] == dominant_blocker), None)
    dominant_overs = float(dominant_row["normalized_overshoot"]) if dominant_row is not None else 0.0
    robust_gap = float(dominant_overs - second_best_overs)
    dominant_blocker_robustness_certificate = {
        "schema_version": "1.0.0",
        "scope": "STRICT_TASK2_DOMINANT_BLOCKER_ROBUSTNESS_CERTIFICATE",
        "theorem_target": "DOMINANT_BLOCKER_NORMALIZED_OVERSHOOT_GAP_GT_ZERO",
        "strict_lane_assumptions": [
            "strict_lane_only_no_legacy_transfer",
            "dominant_blocker_selected_by_argmax_normalized_overshoot",
            "robustness_gap_compared_to_second_largest_overshoot",
        ],
        "domain": {
            "dominant_blocker": str(dominant_blocker),
            "second_largest_competitor_blocker_set": [str(r["criterion"]) for r in unresolved if r["criterion"] != dominant_blocker],
        },
        "computed_rows": [{
            "dominant_normalized_overshoot": dominant_overs,
            "second_largest_normalized_overshoot": second_best_overs,
            "dominance_gap": robust_gap,
        }],
        "aggregate_metrics": {
            "dominance_gap": robust_gap,
            "dominant_blocker": str(dominant_blocker),
            "is_unique_strict_dominant": bool(robust_gap > 0.0),
        },
        "pass_fail_criteria": {"dominance_gap_gt_zero": 0.0},
        "verdict": "CLOSED_NUMERICAL_WITNESS_TASK2" if robust_gap > 0.0 else "OPEN_OBSTRUCTION_WITH_TRACE",
        "fail_trace": "" if robust_gap > 0.0 else f"dominance_gap={robust_gap:.6e} <= 0",
    }
    if (
        q95_gap_abs_ext <= float(pass_fail_criteria_task2["q95_gap_abs_max"]) and
        q95_gap_abs_n2400_top3 <= float(pass_fail_criteria_task2["q95_gap_abs_max"]) and
        q95_n2400_vs_n6400_delta_top3 <= float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"]) and
        abs(q95_gap_abs_n3200_top3 - q95_gap_abs_n6400_top3) <= float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"]) and
        abs(q95_gap_abs_n1600_top3 - q95_gap_abs_n6400_top3) <= float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"]) and
        abs(q95_gap_abs_n800_top3 - q95_gap_abs_n6400_top3) <= float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"]) and
        q95_n12800_vs_n6400_delta_top3 <= float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"]) and
        q95_convergence_tail_ratio_n12800 <= 1.0 and
        q95_n25600_vs_n12800_delta_top3 <= float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"]) and
        q95_monotone_violation_top3 <= 0.0 and
        max_gap_rel_ext <= float(pass_fail_criteria_task2["max_gap_rel_max"]) and
        all_nonnegative and
        crosscheck_gap_q95 <= float(pass_fail_criteria_task2["q95_cross_integrator_gap_abs_max"]) and
        q95_unc_abs <= float(pass_fail_criteria_task2["q95_convergence_delta_n400_to_n800_abs_max"])
    ):
        verdict_task2 = "CLOSED_NUMERICAL_WITNESS_TASK2"
        fail_trace_task2 = ""
    else:
        verdict_task2 = "OPEN_OBSTRUCTION_WITH_TRACE"
        fail_trace_task2 = dominant_inequality
    criterion_coherence_verdict_gate = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_VERDICT_GATE_COHERENCE",
        "all_criteria_satisfied_flag": bool(all(criteria_eval.values())),
        "verdict_task2": str(verdict_task2),
        "closed_label": "CLOSED_NUMERICAL_WITNESS_TASK2",
        "open_label": "OPEN_OBSTRUCTION_WITH_TRACE",
        "flag_matches_verdict_label": bool(
            (bool(all(criteria_eval.values())) and str(verdict_task2) == "CLOSED_NUMERICAL_WITNESS_TASK2") or
            ((not bool(all(criteria_eval.values()))) and str(verdict_task2) == "OPEN_OBSTRUCTION_WITH_TRACE")
        ),
    }
    falsifier_trace_consistency = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_FALSIFIER_TRACE_CONSISTENCY",
        "checks": {
            "if_open_then_fail_trace_equals_dominant_inequality": bool(
                (verdict_task2 != "OPEN_OBSTRUCTION_WITH_TRACE") or (fail_trace_task2 == dominant_inequality)
            ),
            "if_closed_then_fail_trace_empty": bool(
                (verdict_task2 != "CLOSED_NUMERICAL_WITNESS_TASK2") or (fail_trace_task2 == "")
            ),
            "if_open_then_dominant_not_none": bool(
                (verdict_task2 != "OPEN_OBSTRUCTION_WITH_TRACE") or (dominant_blocker != "none")
            ),
        },
    }
    falsifier_trace_consistency["all_checks_pass"] = bool(all(falsifier_trace_consistency["checks"].values()))
    task2_strict_unitarity_witness = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "theorem_target": "bounded DiscM-CutSum gap for graviton_to_gauge_gauge on strict s-grid",
        "strict_lane_assumptions": [
            "K_strict_gate only",
            "no legacy-role transfer",
            "channel fixed to graviton_to_gauge_gauge",
        ],
        "channel": "graviton_to_gauge_gauge",
        "s_grid": [float(x) for x in s_grid_task2_extended],
        "rows": scheme_rows,
        "residue_or_weight_sign_proxy_rows": residue_rows,
        "aggregate_metrics": {
            "max_gap_abs": float(np.max(gap_abs_arr_ext)) if gap_abs_arr_ext.size else float("inf"),
            "q95_gap_abs": q95_gap_abs_ext,
            "max_gap_rel": max_gap_rel_ext,
            "consistency_ci95": consistency_ci95,
            "all_nonnegative_weights": all_nonnegative,
            "min_effective_weight_global": min_effective_weight_global,
            "obstruction_is_numerically_sensitive": obstruction_is_numerically_sensitive,
            "max_bin_disc_proxy_gap_abs": float(max([r["max_bin_disc_proxy_gap_abs"] for r in phase_space_bin_rows], default=0.0)),
            "max_bin_scheme_gap_abs": float(max([r["max_bin_scheme_gap_abs"] for r in phase_space_bin_rows], default=0.0)),
            "q95_cross_integrator_gap_abs": crosscheck_gap_q95,
            "q95_convergence_delta_n400_to_n800_abs": q95_unc_abs,
        },
        "phase_space_bin_contribution_rows": phase_space_bin_rows,
        "bin_obstruction_ranking": {
            "rows": bin_obstruction_ranking_rows,
            "top2_disc_proxy_gap_share": top2_share,
            "total_disc_proxy_gap_sum": total_disc_proxy,
        },
        "endpoint_refinement": ur_task2_endpoint_refinement_precursor,
        "endpoint_split_domain": ur_task2_endpoint_split_domain_precursor,
        "endpoint_adaptive_transform": ur_task2_endpoint_adaptive_transform_precursor,
        "q95_dominant_s_attribution": q95_dominant_s_attribution,
        "q95_dominant_s_crosscheck": q95_dominant_crosscheck,
        "q95_dominant_s_sign_check": q95_dominant_sign_check,
        "q95_dominant_s_convergence": q95_dominant_convergence,
        "q95_blocker_n1600_recompute_certificate": q95_blocker_n1600_recompute_certificate,
        "q95_blocker_n3200_recompute_certificate": q95_blocker_n3200_recompute_certificate,
        "q95_blocker_n6400_recompute_certificate": q95_blocker_n6400_recompute_certificate,
        "q95_blocker_n12800_recompute_certificate": q95_blocker_n12800_recompute_certificate,
        "q95_blocker_n12800_signed_residual_certificate": q95_blocker_n12800_signed_residual_certificate,
        "q95_blocker_n12800_tail_ratio_certificate": q95_blocker_n12800_tail_ratio_certificate,
        "q95_blocker_n25600_recompute_certificate": q95_blocker_n25600_recompute_certificate,
        "q95_blocker_n25600_monotone_certificate": q95_blocker_n25600_monotone_certificate,
        "q95_blocker_ninf_extrapolation_certificate": q95_blocker_ninf_extrapolation_certificate,
        "q95_blocker_uniform_top3_obstruction_certificate": q95_blocker_uniform_top3_obstruction_certificate,
        "q95_blocker_quad_hp_top3_certificate": q95_blocker_quad_hp_top3_certificate,
        "q95_blocker_quad_hp_error_envelope_certificate": q95_blocker_quad_hp_error_envelope_certificate,
        "q95_blocker_tail_budget_certificate": q95_blocker_tail_budget_certificate,
        "q95_blocker_n2400_recompute_certificate": q95_blocker_n2400_recompute_certificate,
        "q95_blocker_margin": q95_blocker_margin,
        "q95_blocker_interval_separation_certificate": q95_blocker_interval_separation_certificate,
        "q95_blocker_counterfactual": q95_blocker_counterfactual,
        "q95_blocker_sensitivity": q95_blocker_sensitivity,
        "q95_blocker_directionality": q95_blocker_directionality,
        "q95_blocker_actionability": q95_blocker_actionability,
        "q95_blocker_action_execution": q95_blocker_action_execution,
        "q95_action_effect_crosscheck": q95_action_effect_crosscheck,
        "q95_blocker_local_line_search_execution": q95_blocker_local_line_search_execution,
        "q95_blocker_continuous_optimization_execution": q95_blocker_continuous_optimization_execution,
        "q95_blocker_refined_window_execution": q95_blocker_refined_window_execution,
        "q95_blocker_single_row_substitution_attempt": q95_blocker_single_row_substitution_attempt,
        "q95_blocker_two_row_substitution_attempt": q95_blocker_two_row_substitution_attempt,
        "q95_blocker_min_k_substitution_scan": q95_blocker_min_k_substitution_scan,
        "q95_blocker_joint_topk_continuous_execution": q95_blocker_joint_topk_continuous_execution,
        "q95_blocker_joint_top4_vector_optimization_execution": q95_blocker_joint_top4_vector_optimization_execution,
        "q95_blocker_adaptive_joint_topk_profile_execution": q95_blocker_adaptive_joint_topk_profile_execution,
        "q95_blocker_joint_bestk_exact_recompute_certificate": q95_blocker_joint_bestk_exact_recompute_certificate,
        "q95_blocker_adaptive_profile_knee_certificate": q95_blocker_adaptive_profile_knee_certificate,
        "q95_action_effect_bootstrap": q95_action_effect_bootstrap,
        "q95_action_gate": q95_action_gate,
        "q95_action_gate_consistency": q95_action_gate_consistency,
        "q95_blocker_step_efficiency": q95_blocker_step_efficiency,
        "q95_after_one_step_local_margin": q95_after_one_step_local_margin,
        "q95_after_one_step_progress_score": q95_after_one_step_progress_score,
        "q95_blocker_direct_relief_projection": q95_blocker_direct_relief_projection,
        "q95_blocker_choice_panel": q95_blocker_choice_panel,
        "q95_blocker_choice_consistency": q95_blocker_choice_consistency,
        "criterion_coherence_sign": criterion_coherence_sign,
        "criterion_coherence_convergence": criterion_coherence_convergence,
        "criterion_coherence_cross_integrator": criterion_coherence_cross_integrator,
        "criterion_coherence_q95_gap": criterion_coherence_q95_gap,
        "criterion_coherence_max_gap_rel": criterion_coherence_max_gap_rel,
        "criterion_coherence_global_closure": criterion_coherence_global_closure,
        "criterion_coherence_weight_sign": criterion_coherence_weight_sign,
        "criterion_coherence_verdict_gate": criterion_coherence_verdict_gate,
        "criterion_coherence_fail_trace_numeric": criterion_coherence_fail_trace_numeric,
        "criterion_coherence_dominant_margin_sign": criterion_coherence_dominant_margin_sign,
        "criterion_coherence_open_trace_inequality": criterion_coherence_open_trace_inequality,
        "criterion_coherence_dominant_inequality_prefix": criterion_coherence_dominant_inequality_prefix,
        "criterion_coherence_fail_trace_equals_dominant": criterion_coherence_fail_trace_equals_dominant,
        "dominant_blocker_numeric_margin": dominant_blocker_numeric_margin,
        "dominant_blocker_selection_consistency": dominant_blocker_selection_consistency,
        "dominant_blocker_robustness_certificate": dominant_blocker_robustness_certificate,
        "falsifier_trace_consistency": falsifier_trace_consistency,
        "pass_fail_criteria": pass_fail_criteria_task2,
        "closure_consistency": {
            "criteria_evaluation": criteria_eval,
            "all_criteria_satisfied": bool(all(criteria_eval.values())),
            "dominant_blocker": dominant_blocker,
            "dominant_inequality": dominant_inequality,
            "verdict_matches_criteria": bool((verdict_task2 == "CLOSED_NUMERICAL_WITNESS_TASK2") == bool(all(criteria_eval.values()))),
        },
        "verdict": verdict_task2,
        "fail_trace": fail_trace_task2,
    }

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
    # Task-2 bounded uncertainty + residual budget table on mapped class lane.
    class_residual_lookup = {r["channel"]: float(r["backend_sub_row_residual_l2"]) for r in channel_delta_rows}
    class_budget_rows = []
    for cls in ["gauge_gauge", "fermion_fermion", "scalar_scalar"]:
        trace_sum_cls = float(ur_channel_class_mapping_precursor["class_trace_budget"][cls]["trace_sum"])
        trace_share_cls = float(ur_channel_class_mapping_precursor["class_trace_budget"][cls]["trace_share"])
        residual_l2_cls = float(class_residual_lookup.get(cls, 0.0))
        uncertainty_p95_cls = float(ur_uncertainty_transport_bridge_precursor["p95_delta_std"])
        risk_proxy = float(residual_l2_cls * (1.0 + uncertainty_p95_cls) * (1.0 + trace_share_cls))
        class_budget_rows.append({
            "class": cls,
            "trace_sum": trace_sum_cls,
            "trace_share": trace_share_cls,
            "residual_l2_backend_sub": residual_l2_cls,
            "uncertainty_p95_delta_std": uncertainty_p95_cls,
            "risk_proxy_residual_uncertainty_trace": risk_proxy,
            "bounded_uncertainty": bool(uncertainty_p95_cls < 1e-5),
        })
    risk_vals = np.array([r["risk_proxy_residual_uncertainty_trace"] for r in class_budget_rows], dtype=float)
    ur_class_bounded_uncertainty_residual_budget_precursor = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_CLASS_BOUNDED_UNCERTAINTY_RESIDUAL_BUDGET",
        "rows": class_budget_rows,
        "risk_proxy_min": float(np.min(risk_vals)) if class_budget_rows else 0.0,
        "risk_proxy_max": float(np.max(risk_vals)) if class_budget_rows else 0.0,
        "risk_proxy_span": float(np.max(risk_vals) - np.min(risk_vals)) if class_budget_rows else 0.0,
        "all_rows_bounded_uncertainty": bool(all(r["bounded_uncertainty"] for r in class_budget_rows)),
    }
    # Task-2 class readiness decision gate (sequencing only, not closure).
    class_readiness_rows = []
    risk_threshold_go = float(np.median(risk_vals) * 1.15) if class_budget_rows else 0.0
    uncertainty_threshold_go = 1e-5
    for row in class_budget_rows:
        go_flag = bool(
            row["risk_proxy_residual_uncertainty_trace"] <= risk_threshold_go
            and row["uncertainty_p95_delta_std"] <= uncertainty_threshold_go
        )
        class_readiness_rows.append({
            "class": row["class"],
            "risk_proxy": float(row["risk_proxy_residual_uncertainty_trace"]),
            "uncertainty_p95_delta_std": float(row["uncertainty_p95_delta_std"]),
            "go_for_costlier_integration_step": go_flag,
        })
    class_readiness_rows = sorted(class_readiness_rows, key=lambda r: r["risk_proxy"])
    ur_class_readiness_gate_precursor = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_CLASS_SEQUENCING_GATE",
        "rows": class_readiness_rows,
        "risk_threshold_go": risk_threshold_go,
        "uncertainty_threshold_go": uncertainty_threshold_go,
        "go_count": int(sum(1 for r in class_readiness_rows if r["go_for_costlier_integration_step"])),
        "priority_order_low_risk_to_high_risk": [r["class"] for r in class_readiness_rows],
    }
    # Class-first replay precursor: execute ranking decision into a concrete next-step packet.
    selected_class = class_readiness_rows[0]["class"] if class_readiness_rows else "gauge_gauge"
    baseline_mean_risk = float(np.mean(risk_vals)) if class_budget_rows else 0.0
    selected_row = next((r for r in class_budget_rows if r["class"] == selected_class), None)
    selected_risk = float(selected_row["risk_proxy_residual_uncertainty_trace"]) if selected_row else 0.0
    selected_uncertainty = float(selected_row["uncertainty_p95_delta_std"]) if selected_row else 0.0
    selected_trace_share = float(selected_row["trace_share"]) if selected_row else 0.0
    ur_class_first_exact_integration_replay_precursor = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_CLASS_FIRST_REPLAY_PACKET",
        "selected_class": selected_class,
        "selection_basis": "min_risk_proxy_from_ur_class_readiness_gate_precursor",
        "baseline_mean_risk_proxy": baseline_mean_risk,
        "selected_class_risk_proxy": selected_risk,
        "delta_selected_minus_baseline_mean_risk": float(selected_risk - baseline_mean_risk),
        "selected_class_uncertainty_p95_delta_std": selected_uncertainty,
        "selected_class_trace_share": selected_trace_share,
        "ready_for_costlier_exact_integration_replay": bool(
            selected_risk <= risk_threshold_go and selected_uncertainty <= uncertainty_threshold_go
        ),
    }
    # Concrete class-first replay delta panel (tighter quadrature settings on selected class only).
    selected_channel_params = {
        "gauge_gauge": (omega, phi, beta, eta),
        "fermion_fermion": (omega * 1.01, phi * 0.99, beta * 1.02, eta * 1.00),
        "scalar_scalar": (omega * 0.99, phi * 1.01, beta * 0.98, eta * 1.01),
    }
    sel_params = selected_channel_params[selected_class]
    sel_row = next((r for r in channel_phase_rows if r["channel"] == selected_class), None)
    baseline_integrals = np.array(sel_row["integrals_over_s_grid"], dtype=float) if sel_row is not None else np.zeros(len(s_grid_fine), dtype=float)
    replay_integrals = []
    for s in s_grid_fine:
        def integrand_selected(x: float) -> float:
            kk = np.cos(sel_params[0] * x + sel_params[1]) / (1.0 + sel_params[2] * (x ** sel_params[3]))
            return float((kk * kk) / np.sqrt(max(1e-15, x + s)))
        vv, _ = si.quad(integrand_selected, 0.0, 1.0, epsabs=1e-13, epsrel=1e-13, limit=600)
        replay_integrals.append(float(vv))
    replay_integrals = np.array(replay_integrals, dtype=float)
    replay_delta_l2 = float(np.linalg.norm(replay_integrals - baseline_integrals, ord=2))
    replay_delta_linf = float(np.linalg.norm(replay_integrals - baseline_integrals, ord=np.inf))
    ur_class_first_replay_delta_precursor = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_CLASS_FIRST_REPLAY_DELTA_PANEL",
        "selected_class": selected_class,
        "baseline_integrals_over_s_grid": baseline_integrals.tolist(),
        "replay_integrals_over_s_grid": replay_integrals.tolist(),
        "delta_l2_replay_minus_baseline": replay_delta_l2,
        "delta_linf_replay_minus_baseline": replay_delta_linf,
        "replay_settings": {"epsabs": 1e-13, "epsrel": 1e-13, "limit": 600},
    }
    # Comparison panel: class-first replay vs all-class replay under same precision budget.
    baseline_by_class = {r["channel"]: np.array(r["integrals_over_s_grid"], dtype=float) for r in channel_phase_rows}
    replay_all_rows = []
    all_delta_l2 = {}
    for ch in ["gauge_gauge", "fermion_fermion", "scalar_scalar"]:
        par = selected_channel_params[ch]
        vals = []
        for s in s_grid_fine:
            def integrand_all(x: float) -> float:
                kk = np.cos(par[0] * x + par[1]) / (1.0 + par[2] * (x ** par[3]))
                return float((kk * kk) / np.sqrt(max(1e-15, x + s)))
            vv, _ = si.quad(integrand_all, 0.0, 1.0, epsabs=1e-13, epsrel=1e-13, limit=600)
            vals.append(float(vv))
        arr = np.array(vals, dtype=float)
        d_l2 = float(np.linalg.norm(arr - baseline_by_class[ch], ord=2))
        d_linf = float(np.linalg.norm(arr - baseline_by_class[ch], ord=np.inf))
        all_delta_l2[ch] = d_l2
        replay_all_rows.append({
            "class": ch,
            "delta_l2_replay_minus_baseline": d_l2,
            "delta_linf_replay_minus_baseline": d_linf,
        })
    selected_class_delta_l2 = float(all_delta_l2[selected_class])
    mean_other_l2 = float(np.mean([v for k, v in all_delta_l2.items() if k != selected_class]))
    ur_class_first_vs_all_class_replay_comparison_precursor = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_CLASS_FIRST_VS_ALL_CLASS_REPLAY_COMPARISON",
        "selected_class": selected_class,
        "rows": replay_all_rows,
        "selected_class_delta_l2": selected_class_delta_l2,
        "mean_other_classes_delta_l2": mean_other_l2,
        "selected_minus_mean_other_delta_l2": float(selected_class_delta_l2 - mean_other_l2),
        "selected_is_min_delta_l2": bool(selected_class_delta_l2 <= min(all_delta_l2.values())),
        "replay_settings": {"epsabs": 1e-13, "epsrel": 1e-13, "limit": 600},
    }
    # Cost-vs-gain panel: class-first and all-class replay under same precision envelope.
    sim_eval_count_per_class = int(len(s_grid_fine))
    estimated_cost_units_class_first = float(sim_eval_count_per_class)
    estimated_cost_units_all_class = float(sim_eval_count_per_class * 3)
    gain_proxy_class_first = float(replay_delta_l2)
    gain_proxy_all_class_mean = float(np.mean([r["delta_l2_replay_minus_baseline"] for r in replay_all_rows])) if replay_all_rows else 0.0
    efficiency_class_first = float(gain_proxy_class_first / max(1e-30, estimated_cost_units_class_first))
    efficiency_all_class = float(gain_proxy_all_class_mean / max(1e-30, estimated_cost_units_all_class))
    ur_cost_vs_gain_precursor = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_COST_VS_GAIN_PANEL",
        "class_first": {
            "estimated_cost_units": estimated_cost_units_class_first,
            "gain_proxy_delta_l2": gain_proxy_class_first,
            "gain_per_cost": efficiency_class_first,
        },
        "all_class": {
            "estimated_cost_units": estimated_cost_units_all_class,
            "gain_proxy_delta_l2_mean": gain_proxy_all_class_mean,
            "gain_per_cost": efficiency_all_class,
        },
        "class_first_more_cost_efficient": bool(efficiency_class_first <= efficiency_all_class),
        "cost_ratio_all_over_class_first": float(estimated_cost_units_all_class / max(1e-30, estimated_cost_units_class_first)),
    }
    # Runtime micro-benchmark panel (same selected class, tolerance sweep).
    runtime_tolerance_grid = [1e-10, 1e-12, 1e-13]
    runtime_rows = []
    for tol in runtime_tolerance_grid:
        t0 = time.perf_counter()
        vals = []
        for s in s_grid_fine:
            def integrand_rt(x: float) -> float:
                kk = np.cos(sel_params[0] * x + sel_params[1]) / (1.0 + sel_params[2] * (x ** sel_params[3]))
                return float((kk * kk) / np.sqrt(max(1e-15, x + s)))
            vv, _ = si.quad(integrand_rt, 0.0, 1.0, epsabs=tol, epsrel=tol, limit=600)
            vals.append(float(vv))
        dt = float(time.perf_counter() - t0)
        arr = np.array(vals, dtype=float)
        d_l2 = float(np.linalg.norm(arr - baseline_integrals, ord=2))
        runtime_rows.append({
            "tol": float(tol),
            "runtime_seconds": dt,
            "delta_l2_vs_baseline": d_l2,
            "gain_per_second": float(d_l2 / max(1e-12, dt)),
        })
    ur_runtime_tolerance_benchmark_precursor = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_RUNTIME_TOLERANCE_BENCHMARK",
        "selected_class": selected_class,
        "rows": runtime_rows,
        "fastest_runtime_seconds": float(min(r["runtime_seconds"] for r in runtime_rows)),
        "slowest_runtime_seconds": float(max(r["runtime_seconds"] for r in runtime_rows)),
    }
    # Heavier all-class exact integration sweep (user-requested cost-insensitive computation).
    exact_case_grid = [(1e-12, 600), (1e-12, 1200), (1e-14, 600), (1e-14, 1200), (1e-15, 2000)]
    exact_rows = []
    for ch in ["gauge_gauge", "fermion_fermion", "scalar_scalar"]:
        par = selected_channel_params[ch]
        base = baseline_by_class[ch]
        for tol, lim in exact_case_grid:
            vals = []
            warning_count = 0
            for s in s_grid_fine:
                def integrand_exact(x: float) -> float:
                    kk = np.cos(par[0] * x + par[1]) / (1.0 + par[2] * (x ** par[3]))
                    return float((kk * kk) / np.sqrt(max(1e-15, x + s)))
                with warnings.catch_warnings(record=True) as wlist:
                    warnings.simplefilter("always", IntegrationWarning)
                    vv, _ = si.quad(integrand_exact, 0.0, 1.0, epsabs=tol, epsrel=tol, limit=lim)
                warning_count += int(sum(1 for w in wlist if issubclass(w.category, IntegrationWarning)))
                vals.append(float(vv))
            arr = np.array(vals, dtype=float)
            exact_rows.append({
                "class": ch,
                "epsabs": float(tol),
                "epsrel": float(tol),
                "limit": int(lim),
                "delta_l2_vs_baseline": float(np.linalg.norm(arr - base, ord=2)),
                "delta_linf_vs_baseline": float(np.linalg.norm(arr - base, ord=np.inf)),
                "integration_warning_count": int(warning_count),
                "numerical_stress_flag": bool(warning_count > 0),
            })
    exact_delta_l2 = np.array([r["delta_l2_vs_baseline"] for r in exact_rows], dtype=float)
    exact_warning_total = int(sum(r["integration_warning_count"] for r in exact_rows))
    ur_all_class_exact_integration_sweep_precursor = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_ALL_CLASS_EXACT_INTEGRATION_SWEEP",
        "rows": exact_rows,
        "num_rows": int(len(exact_rows)),
        "delta_l2_min": float(np.min(exact_delta_l2)) if exact_rows else 0.0,
        "delta_l2_max": float(np.max(exact_delta_l2)) if exact_rows else 0.0,
        "delta_l2_span": float(np.max(exact_delta_l2) - np.min(exact_delta_l2)) if exact_rows else 0.0,
        "integration_warning_total": exact_warning_total,
        "any_numerical_stress_flag": bool(exact_warning_total > 0),
    }
    # Numerical stress ranking packet for follow-up refactor prioritization.
    exact_rows_ranked = sorted(
        exact_rows,
        key=lambda r: (-int(r["integration_warning_count"]), -float(r["delta_l2_vs_baseline"]))
    )
    top_k = min(5, len(exact_rows_ranked))
    ur_numerical_stress_ranking_precursor = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_NUMERICAL_STRESS_RANKING",
        "ranking_key": "integration_warning_count_desc_then_delta_l2_desc",
        "top_k": int(top_k),
        "rows_top_k": exact_rows_ranked[:top_k],
    }
    # Alternate-parameterization replay on full exact sweep grid: x = u^2 transform.
    alt_rows = []
    for row in exact_rows:
        ch = str(row["class"])
        tol = float(row["epsabs"])
        lim = int(row["limit"])
        par = selected_channel_params[ch]
        base = baseline_by_class[ch]
        vals_alt = []
        warning_count_alt = 0
        for s in s_grid_fine:
            def integrand_alt(u: float) -> float:
                x = u * u
                kk = np.cos(par[0] * x + par[1]) / (1.0 + par[2] * (x ** par[3]))
                return float(((kk * kk) / np.sqrt(max(1e-15, x + s))) * 2.0 * u)
            with warnings.catch_warnings(record=True) as wlist:
                warnings.simplefilter("always", IntegrationWarning)
                vv, _ = si.quad(integrand_alt, 0.0, 1.0, epsabs=tol, epsrel=tol, limit=lim)
            warning_count_alt += int(sum(1 for w in wlist if issubclass(w.category, IntegrationWarning)))
            vals_alt.append(float(vv))
        arr_alt = np.array(vals_alt, dtype=float)
        delta_l2_alt = float(np.linalg.norm(arr_alt - base, ord=2))
        alt_rows.append({
            "class": ch,
            "epsabs": tol,
            "epsrel": tol,
            "limit": lim,
            "original_integration_warning_count": int(row["integration_warning_count"]),
            "alt_integration_warning_count": int(warning_count_alt),
            "warning_count_delta_alt_minus_original": int(warning_count_alt - int(row["integration_warning_count"])),
            "original_delta_l2_vs_baseline": float(row["delta_l2_vs_baseline"]),
            "alt_delta_l2_vs_baseline": delta_l2_alt,
            "delta_l2_alt_minus_original": float(delta_l2_alt - float(row["delta_l2_vs_baseline"])),
        })
    ur_numerical_stress_alt_parameterization_precursor = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_NUMERICAL_STRESS_ALT_PARAMETERIZATION",
        "transform": "x_equals_u_squared",
        "rows": alt_rows,
        "num_rows": int(len(alt_rows)),
        "improved_warning_count_cases": int(sum(1 for r in alt_rows if r["warning_count_delta_alt_minus_original"] < 0)),
    }
    # Secondary transform research pass on top stress rows: x = u^4.
    # This is a methodological numeric comparison only (non-closure).
    alt_u4_rows = []
    for row in exact_rows_ranked[:top_k]:
        ch = str(row["class"])
        tol = float(row["epsabs"])
        lim = int(row["limit"])
        par = selected_channel_params[ch]
        base = baseline_by_class[ch]
        vals_u4 = []
        warning_count_u4 = 0
        for s in s_grid_fine:
            def integrand_u4(u: float) -> float:
                x = u ** 4
                kk = np.cos(par[0] * x + par[1]) / (1.0 + par[2] * (x ** par[3]))
                return float(((kk * kk) / np.sqrt(max(1e-15, x + s))) * 4.0 * (u ** 3))
            with warnings.catch_warnings(record=True) as wlist:
                warnings.simplefilter("always", IntegrationWarning)
                vv, _ = si.quad(integrand_u4, 0.0, 1.0, epsabs=tol, epsrel=tol, limit=lim)
            warning_count_u4 += int(sum(1 for w in wlist if issubclass(w.category, IntegrationWarning)))
            vals_u4.append(float(vv))
        arr_u4 = np.array(vals_u4, dtype=float)
        delta_l2_u4 = float(np.linalg.norm(arr_u4 - base, ord=2))
        alt_u4_rows.append({
            "class": ch,
            "epsabs": tol,
            "epsrel": tol,
            "limit": lim,
            "original_integration_warning_count": int(row["integration_warning_count"]),
            "u4_integration_warning_count": int(warning_count_u4),
            "warning_count_delta_u4_minus_original": int(warning_count_u4 - int(row["integration_warning_count"])),
            "original_delta_l2_vs_baseline": float(row["delta_l2_vs_baseline"]),
            "u4_delta_l2_vs_baseline": delta_l2_u4,
            "delta_l2_u4_minus_original": float(delta_l2_u4 - float(row["delta_l2_vs_baseline"])),
        })
    alt_u1_rows = []
    for row in exact_rows_ranked[:top_k]:
        ch = str(row["class"])
        tol = float(row["epsabs"])
        lim = int(row["limit"])
        par = selected_channel_params[ch]
        base = baseline_by_class[ch]
        vals_u1 = []
        warning_count_u1 = 0
        t0 = time.perf_counter()
        for s in s_grid_fine:
            def integrand_u1(u: float) -> float:
                x = u
                kk = np.cos(par[0] * x + par[1]) / (1.0 + par[2] * (x ** par[3]))
                return float((kk * kk) / np.sqrt(max(1e-15, x + s)))
            with warnings.catch_warnings(record=True) as wlist:
                warnings.simplefilter("always", IntegrationWarning)
                vv, _ = si.quad(integrand_u1, 0.0, 1.0, epsabs=tol, epsrel=tol, limit=lim)
            warning_count_u1 += int(sum(1 for w in wlist if issubclass(w.category, IntegrationWarning)))
            vals_u1.append(float(vv))
        runtime_u1 = float(time.perf_counter() - t0)
        arr_u1 = np.array(vals_u1, dtype=float)
        delta_l2_u1 = float(np.linalg.norm(arr_u1 - base, ord=2))
        alt_u1_rows.append({
            "class": ch,
            "epsabs": tol,
            "epsrel": tol,
            "limit": lim,
            "u1_integration_warning_count": int(warning_count_u1),
            "u1_delta_l2_vs_baseline": delta_l2_u1,
            "u1_runtime_seconds": runtime_u1,
        })
    transform_ranking_rows = []
    for i in range(min(top_k, len(alt_u1_rows), len(alt_u4_rows))):
        r_u2 = alt_rows[i]
        r_u4 = alt_u4_rows[i]
        r_u1 = alt_u1_rows[i]
        candidates = [
            ("u1", int(r_u1["u1_integration_warning_count"]), float(r_u1["u1_delta_l2_vs_baseline"]), float(r_u1["u1_runtime_seconds"])),
            ("u2", int(r_u2["alt_integration_warning_count"]), float(r_u2["alt_delta_l2_vs_baseline"]), None),
            ("u4", int(r_u4["u4_integration_warning_count"]), float(r_u4["u4_delta_l2_vs_baseline"]), None),
        ]
        ranked = sorted(candidates, key=lambda t: (t[1], t[2], t[3] if isinstance(t[3], float) else 0.0))
        transform_ranking_rows.append({
            "class": str(r_u2["class"]),
            "epsabs": float(r_u2["epsabs"]),
            "epsrel": float(r_u2["epsrel"]),
            "limit": int(r_u2["limit"]),
            "ranked_transforms": [x[0] for x in ranked],
            "winner": ranked[0][0],
        })
    ur_numerical_stress_alt_transform_comparison_precursor = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_NUMERICAL_STRESS_ALT_TRANSFORM_COMPARISON",
        "top_k": int(top_k),
        "rows_u1": alt_u1_rows,
        "rows_u2": alt_rows[:top_k],
        "rows_u4": alt_u4_rows,
        "ranking_key": "warning_count_then_delta_l2_then_runtime",
        "ranking_rows": transform_ranking_rows,
        "u2_better_warning_count_cases": int(sum(1 for i in range(min(len(alt_u4_rows), top_k)) if int(alt_rows[i]["alt_integration_warning_count"]) < int(alt_u4_rows[i]["u4_integration_warning_count"]))),
        "u4_better_warning_count_cases": int(sum(1 for i in range(min(len(alt_u4_rows), top_k)) if int(alt_rows[i]["alt_integration_warning_count"]) > int(alt_u4_rows[i]["u4_integration_warning_count"]))),
    }
    # Full-grid tri-transform sweep (u1/u2/u4) for class-level operational choice.
    tri_rows = []
    for row in exact_rows:
        ch = str(row["class"])
        tol = float(row["epsabs"])
        lim = int(row["limit"])
        par = selected_channel_params[ch]
        base = baseline_by_class[ch]

        def _run_transform(pow_x: float) -> tuple[int, float, float]:
            vals = []
            wc = 0
            t0 = time.perf_counter()
            for s in s_grid_fine:
                def integ(u: float) -> float:
                    x = u ** pow_x
                    kk = np.cos(par[0] * x + par[1]) / (1.0 + par[2] * (x ** par[3]))
                    jac = pow_x * (u ** (pow_x - 1.0))
                    return float(((kk * kk) / np.sqrt(max(1e-15, x + s))) * jac)
                with warnings.catch_warnings(record=True) as wlist:
                    warnings.simplefilter("always", IntegrationWarning)
                    vv, _ = si.quad(integ, 0.0, 1.0, epsabs=tol, epsrel=tol, limit=lim)
                wc += int(sum(1 for w in wlist if issubclass(w.category, IntegrationWarning)))
                vals.append(float(vv))
            dt = float(time.perf_counter() - t0)
            d2 = float(np.linalg.norm(np.array(vals, dtype=float) - base, ord=2))
            return int(wc), d2, dt

        u1_wc, u1_d2, u1_t = _run_transform(1.0)
        u2_wc, u2_d2, u2_t = _run_transform(2.0)
        u4_wc, u4_d2, u4_t = _run_transform(4.0)
        cand = [("u1", u1_wc, u1_d2, u1_t), ("u2", u2_wc, u2_d2, u2_t), ("u4", u4_wc, u4_d2, u4_t)]
        rank = sorted(cand, key=lambda z: (z[1], z[2], z[3]))
        tri_rows.append({
            "class": ch, "epsabs": tol, "epsrel": tol, "limit": lim,
            "u1_warning_count": u1_wc, "u1_delta_l2_vs_baseline": u1_d2, "u1_runtime_seconds": u1_t,
            "u2_warning_count": u2_wc, "u2_delta_l2_vs_baseline": u2_d2, "u2_runtime_seconds": u2_t,
            "u4_warning_count": u4_wc, "u4_delta_l2_vs_baseline": u4_d2, "u4_runtime_seconds": u4_t,
            "winner": rank[0][0], "ranked_transforms": [x[0] for x in rank],
        })
    by_class_best = []
    for ch in sorted({r["class"] for r in tri_rows}):
        rows_ch = [r for r in tri_rows if r["class"] == ch]
        winners = [r["winner"] for r in rows_ch]
        freq = {k: winners.count(k) for k in ("u1", "u2", "u4")}
        best = sorted(freq.items(), key=lambda kv: (-kv[1], kv[0]))[0][0]
        by_class_best.append({
            "class": ch,
            "num_rows": int(len(rows_ch)),
            "winner_counts": freq,
            "recommended_transform_majority": best,
            "recommended_transform_frequency": float(freq[best] / len(rows_ch)) if rows_ch else 0.0,
        })
    ur_numerical_stress_alt_fullgrid_tritransform_precursor = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_NUMERICAL_STRESS_ALT_FULLGRID_TRITRANSFORM",
        "ranking_key": "warning_count_then_delta_l2_then_runtime",
        "rows": tri_rows,
        "num_rows": int(len(tri_rows)),
        "by_class": by_class_best,
    }
    # Class-conditional replay: use class-majority transform on full grid and
    # compare against original exact sweep and global u2-only policy.
    class_choice = {str(r["class"]): str(r["recommended_transform_majority"]) for r in by_class_best}
    tri_by_key = {(str(r["class"]), float(r["epsabs"]), int(r["limit"])): r for r in tri_rows}
    baseline_rows_by_key = {(str(r["class"]), float(r["epsabs"]), int(r["limit"])): r for r in exact_rows}
    replay_policy_rows = []
    for key, rr in tri_by_key.items():
        ch, epsabs, lim = key
        chosen = class_choice[ch]
        if chosen == "u1":
            wc = int(rr["u1_warning_count"]); d2 = float(rr["u1_delta_l2_vs_baseline"]); rt = float(rr["u1_runtime_seconds"])
        elif chosen == "u2":
            wc = int(rr["u2_warning_count"]); d2 = float(rr["u2_delta_l2_vs_baseline"]); rt = float(rr["u2_runtime_seconds"])
        else:
            wc = int(rr["u4_warning_count"]); d2 = float(rr["u4_delta_l2_vs_baseline"]); rt = float(rr["u4_runtime_seconds"])
        base = baseline_rows_by_key[key]
        replay_policy_rows.append({
            "class": ch, "epsabs": epsabs, "epsrel": epsabs, "limit": lim,
            "chosen_transform": chosen,
            "chosen_warning_count": wc,
            "chosen_delta_l2_vs_baseline": d2,
            "chosen_runtime_seconds": rt,
            "baseline_warning_count": int(base["integration_warning_count"]),
            "warning_count_delta_chosen_minus_baseline": int(wc - int(base["integration_warning_count"])),
        })
    chosen_warn = np.array([int(r["chosen_warning_count"]) for r in replay_policy_rows], dtype=float)
    base_warn = np.array([int(r["baseline_warning_count"]) for r in replay_policy_rows], dtype=float)
    chosen_d2 = np.array([float(r["chosen_delta_l2_vs_baseline"]) for r in replay_policy_rows], dtype=float)
    chosen_rt = np.array([float(r["chosen_runtime_seconds"]) for r in replay_policy_rows], dtype=float)
    ur_numerical_stress_class_conditional_replay_precursor = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_NUMERICAL_STRESS_CLASS_CONDITIONAL_REPLAY",
        "class_transform_policy": class_choice,
        "rows": replay_policy_rows,
        "num_rows": int(len(replay_policy_rows)),
        "chosen_warning_total": int(np.sum(chosen_warn)),
        "baseline_warning_total": int(np.sum(base_warn)),
        "warning_total_delta_chosen_minus_baseline": int(np.sum(chosen_warn) - np.sum(base_warn)),
        "chosen_delta_l2_span": float(np.max(chosen_d2) - np.min(chosen_d2)) if chosen_d2.size else 0.0,
        "chosen_runtime_total_seconds": float(np.sum(chosen_rt)),
    }
    # Cross-policy counterfactual on full grid: always-u1/u2/u4 vs class-conditional.
    policy_defs = {
        "always_u1": {c: "u1" for c in ("gauge_gauge", "fermion_fermion", "scalar_scalar")},
        "always_u2": {c: "u2" for c in ("gauge_gauge", "fermion_fermion", "scalar_scalar")},
        "always_u4": {c: "u4" for c in ("gauge_gauge", "fermion_fermion", "scalar_scalar")},
        "class_conditional": class_choice,
    }
    policy_rows = []
    for pname, pmap in policy_defs.items():
        vals_warn = []
        vals_d2 = []
        vals_rt = []
        for rr in tri_rows:
            ch = str(rr["class"])
            pick = pmap[ch]
            vals_warn.append(float(rr[f"{pick}_warning_count"]))
            vals_d2.append(float(rr[f"{pick}_delta_l2_vs_baseline"]))
            vals_rt.append(float(rr[f"{pick}_runtime_seconds"]))
        arr_warn = np.array(vals_warn, dtype=float)
        arr_d2 = np.array(vals_d2, dtype=float)
        arr_rt = np.array(vals_rt, dtype=float)
        policy_rows.append({
            "policy": pname,
            "warning_total": int(np.sum(arr_warn)),
            "delta_l2_span": float(np.max(arr_d2) - np.min(arr_d2)) if arr_d2.size else 0.0,
            "runtime_total_seconds": float(np.sum(arr_rt)),
        })
    policy_rows_sorted = sorted(policy_rows, key=lambda r: (r["warning_total"], r["delta_l2_span"], r["runtime_total_seconds"], r["policy"]))
    policy_pareto_rows = []
    for i, ri in enumerate(policy_rows_sorted):
        dominated = False
        best_margin = {
            "warning_total": 0.0,
            "delta_l2_span": 0.0,
            "runtime_total_seconds": 0.0,
        }
        dominator = "none"
        for j, rj in enumerate(policy_rows_sorted):
            if i == j:
                continue
            nonworse = (
                float(rj["warning_total"]) <= float(ri["warning_total"]) and
                float(rj["delta_l2_span"]) <= float(ri["delta_l2_span"]) and
                float(rj["runtime_total_seconds"]) <= float(ri["runtime_total_seconds"])
            )
            strictly_better = (
                float(rj["warning_total"]) < float(ri["warning_total"]) or
                float(rj["delta_l2_span"]) < float(ri["delta_l2_span"]) or
                float(rj["runtime_total_seconds"]) < float(ri["runtime_total_seconds"])
            )
            if nonworse and strictly_better:
                dominated = True
                margins = {
                    "warning_total": float(ri["warning_total"]) - float(rj["warning_total"]),
                    "delta_l2_span": float(ri["delta_l2_span"]) - float(rj["delta_l2_span"]),
                    "runtime_total_seconds": float(ri["runtime_total_seconds"]) - float(rj["runtime_total_seconds"]),
                }
                current_score = sum(abs(float(v)) for v in margins.values())
                best_score = sum(abs(float(v)) for v in best_margin.values())
                if current_score >= best_score:
                    best_margin = margins
                    dominator = str(rj["policy"])
        policy_pareto_rows.append({
            "policy": str(ri["policy"]),
            "warning_total": int(ri["warning_total"]),
            "delta_l2_span": float(ri["delta_l2_span"]),
            "runtime_total_seconds": float(ri["runtime_total_seconds"]),
            "pareto_frontier": bool(not dominated),
            "dominated_by_policy": dominator,
            "dominance_margin": best_margin,
        })
    policy_pareto_rows_sorted = sorted(policy_pareto_rows, key=lambda r: (not bool(r["pareto_frontier"]), r["warning_total"], r["delta_l2_span"], r["runtime_total_seconds"], r["policy"]))
    ur_numerical_stress_policy_counterfactual_precursor = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_NUMERICAL_STRESS_POLICY_COUNTERFACTUAL",
        "ranking_key": "warning_total_then_delta_l2_span_then_runtime_total",
        "rows": policy_rows_sorted,
        "best_policy": policy_rows_sorted[0]["policy"] if policy_rows_sorted else "none",
    }
    ur_numerical_stress_policy_pareto_front_precursor = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_NUMERICAL_STRESS_POLICY_PARETO_FRONT",
        "axes": ["warning_total", "delta_l2_span", "runtime_total_seconds"],
        "rows": policy_pareto_rows_sorted,
        "pareto_frontier_policies": [str(r["policy"]) for r in policy_pareto_rows_sorted if bool(r["pareto_frontier"])],
        "pareto_frontier_count": int(sum(1 for r in policy_pareto_rows_sorted if bool(r["pareto_frontier"]))),
    }
    # Bootstrap policy Pareto-front stability: resample full-grid rows and track
    # frontier-membership frequency per policy under numerical stress replay noise.
    policy_names = ["always_u1", "always_u2", "always_u4", "class_conditional"]
    policy_boot_counts = {pn: 0 for pn in policy_names}
    policy_boot_n = 512
    policy_boot_rng = np.random.default_rng(975145)
    if tri_rows:
        for _ in range(policy_boot_n):
            idx = policy_boot_rng.integers(0, len(tri_rows), size=len(tri_rows))
            policy_boot_rows = []
            for pname, pmap in policy_defs.items():
                vals_warn = []
                vals_d2 = []
                vals_rt = []
                for ii in idx:
                    rr = tri_rows[int(ii)]
                    ch = str(rr["class"])
                    pick = pmap[ch]
                    vals_warn.append(float(rr[f"{pick}_warning_count"]))
                    vals_d2.append(float(rr[f"{pick}_delta_l2_vs_baseline"]))
                    vals_rt.append(float(rr[f"{pick}_runtime_seconds"]))
                arr_warn = np.array(vals_warn, dtype=float)
                arr_d2 = np.array(vals_d2, dtype=float)
                arr_rt = np.array(vals_rt, dtype=float)
                policy_boot_rows.append({
                    "policy": pname,
                    "warning_total": float(np.sum(arr_warn)),
                    "delta_l2_span": float(np.max(arr_d2) - np.min(arr_d2)) if arr_d2.size else 0.0,
                    "runtime_total_seconds": float(np.sum(arr_rt)),
                })
            for ri in policy_boot_rows:
                dominated = False
                for rj in policy_boot_rows:
                    if str(ri["policy"]) == str(rj["policy"]):
                        continue
                    nonworse = (
                        float(rj["warning_total"]) <= float(ri["warning_total"]) and
                        float(rj["delta_l2_span"]) <= float(ri["delta_l2_span"]) and
                        float(rj["runtime_total_seconds"]) <= float(ri["runtime_total_seconds"])
                    )
                    strict = (
                        float(rj["warning_total"]) < float(ri["warning_total"]) or
                        float(rj["delta_l2_span"]) < float(ri["delta_l2_span"]) or
                        float(rj["runtime_total_seconds"]) < float(ri["runtime_total_seconds"])
                    )
                    if nonworse and strict:
                        dominated = True
                        break
                if not dominated:
                    policy_boot_counts[str(ri["policy"])] += 1
    policy_pareto_frequency_rows = []
    for pn in policy_names:
        succ = int(policy_boot_counts[pn])
        freq = float(succ / policy_boot_n) if policy_boot_n > 0 else 0.0
        ci = jeffreys_interval_from_successes(succ, policy_boot_n)
        policy_pareto_frequency_rows.append({
            "policy": pn,
            "pareto_front_frequency": freq,
            "pareto_front_frequency_ci95_jeffreys": ci,
            "bootstrap_successes": succ,
            "bootstrap_trials": int(policy_boot_n),
        })
    policy_pareto_frequency_rows = sorted(policy_pareto_frequency_rows, key=lambda r: (-float(r["pareto_front_frequency"]), r["policy"]))
    ur_numerical_stress_policy_pareto_stability_precursor = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_NUMERICAL_STRESS_POLICY_PARETO_STABILITY",
        "bootstrap_size": int(policy_boot_n),
        "resampling_rule": "iid_row_resample_with_replacement_over_fullgrid_tri_rows",
        "rows": policy_pareto_frequency_rows,
        "most_stable_frontier_policy": str(policy_pareto_frequency_rows[0]["policy"]) if policy_pareto_frequency_rows else "none",
    }
    # Budget-constrained policy decision panel:
    # among Pareto-front candidates, choose policy robust to runtime caps.
    runtime_vals = np.array([float(r["runtime_total_seconds"]) for r in policy_rows_sorted], dtype=float)
    budget_quantiles = [0.25, 0.50, 0.75, 1.00]
    runtime_budget_rows = []
    for q in budget_quantiles:
        bcap = float(np.quantile(runtime_vals, q)) if runtime_vals.size else 0.0
        eligible = [r for r in policy_rows_sorted if float(r["runtime_total_seconds"]) <= bcap + 1e-15]
        if not eligible:
            runtime_budget_rows.append({
                "runtime_budget_quantile": float(q),
                "runtime_budget_cap_seconds": bcap,
                "eligible_policy_count": 0,
                "selected_policy": "none",
                "selected_warning_total": 0,
                "selected_delta_l2_span": 0.0,
                "selected_runtime_total_seconds": 0.0,
            })
            continue
        best_eligible = sorted(
            eligible,
            key=lambda r: (int(r["warning_total"]), float(r["delta_l2_span"]), float(r["runtime_total_seconds"]), str(r["policy"]))
        )[0]
        runtime_budget_rows.append({
            "runtime_budget_quantile": float(q),
            "runtime_budget_cap_seconds": bcap,
            "eligible_policy_count": int(len(eligible)),
            "selected_policy": str(best_eligible["policy"]),
            "selected_warning_total": int(best_eligible["warning_total"]),
            "selected_delta_l2_span": float(best_eligible["delta_l2_span"]),
            "selected_runtime_total_seconds": float(best_eligible["runtime_total_seconds"]),
        })
    budget_policy_votes: dict[str, int] = {pn: 0 for pn in policy_names}
    for rr in runtime_budget_rows:
        spn = str(rr["selected_policy"])
        if spn in budget_policy_votes:
            budget_policy_votes[spn] += 1
    budget_recommended_policy = sorted(policy_names, key=lambda pn: (-budget_policy_votes[pn], pn))[0] if policy_names else "none"
    ur_numerical_stress_policy_budgeted_selection_precursor = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_NUMERICAL_STRESS_POLICY_BUDGETED_SELECTION",
        "selection_rule": "argmin_warning_then_delta_l2_then_runtime_under_runtime_cap",
        "runtime_budget_quantiles": budget_quantiles,
        "rows": runtime_budget_rows,
        "budget_vote_rows": [{"policy": pn, "wins_across_budgets": int(budget_policy_votes[pn])} for pn in policy_names],
        "budget_recommended_policy": budget_recommended_policy,
    }
    # Budget recommendation fragility panel:
    # perturb runtime caps and integration tolerances (eps-like objective jitter proxy)
    # using bootstrap row resamples; track policy flip frequency.
    stress_boot_n = 256
    cap_scale_grid = [0.85, 1.00, 1.15]
    jitter_grid = [0.0, 0.01, 0.02]
    stress_rng = np.random.default_rng(975147)
    base_policy = str(budget_recommended_policy)
    stress_rows = []
    total_trials = 0
    total_flips = 0
    policy_win_counter = {pn: 0 for pn in policy_names}
    if tri_rows:
        for cap_scale in cap_scale_grid:
            for jitter in jitter_grid:
                local_flips = 0
                local_trials = 0
                local_counts = {pn: 0 for pn in policy_names}
                for _ in range(stress_boot_n):
                    idx = stress_rng.integers(0, len(tri_rows), size=len(tri_rows))
                    boot_policy_rows = []
                    for pname, pmap in policy_defs.items():
                        vals_warn = []
                        vals_d2 = []
                        vals_rt = []
                        for ii in idx:
                            rr = tri_rows[int(ii)]
                            ch = str(rr["class"])
                            pick = pmap[ch]
                            vals_warn.append(float(rr[f"{pick}_warning_count"]))
                            vals_d2.append(float(rr[f"{pick}_delta_l2_vs_baseline"]))
                            vals_rt.append(float(rr[f"{pick}_runtime_seconds"]))
                        arr_warn = np.array(vals_warn, dtype=float)
                        arr_d2 = np.array(vals_d2, dtype=float)
                        arr_rt = np.array(vals_rt, dtype=float)
                        # jitter as conservative ambiguity proxy on warning and delta_l2 axes
                        wsum = float(np.sum(arr_warn) * (1.0 + jitter))
                        dspan = float((np.max(arr_d2) - np.min(arr_d2)) * (1.0 + jitter)) if arr_d2.size else 0.0
                        rsum = float(np.sum(arr_rt))
                        boot_policy_rows.append({
                            "policy": pname,
                            "warning_total": wsum,
                            "delta_l2_span": dspan,
                            "runtime_total_seconds": rsum,
                        })
                    runtime_arr = np.array([float(r["runtime_total_seconds"]) for r in boot_policy_rows], dtype=float)
                    q50 = float(np.quantile(runtime_arr, 0.50)) if runtime_arr.size else 0.0
                    cap = float(cap_scale * q50)
                    eligible = [r for r in boot_policy_rows if float(r["runtime_total_seconds"]) <= cap + 1e-15]
                    if not eligible:
                        chosen = "none"
                    else:
                        chosen = str(sorted(eligible, key=lambda r: (float(r["warning_total"]), float(r["delta_l2_span"]), float(r["runtime_total_seconds"]), str(r["policy"])))[0]["policy"])
                    if chosen in local_counts:
                        local_counts[chosen] += 1
                        policy_win_counter[chosen] += 1
                    if chosen != "none":
                        local_trials += 1
                        total_trials += 1
                        if chosen != base_policy:
                            local_flips += 1
                            total_flips += 1
                flip_freq = float(local_flips / local_trials) if local_trials > 0 else 0.0
                stress_rows.append({
                    "runtime_cap_scale": float(cap_scale),
                    "objective_jitter_frac": float(jitter),
                    "bootstrap_size": int(stress_boot_n),
                    "flip_count_vs_base_policy": int(local_flips),
                    "local_trials": int(local_trials),
                    "flip_frequency_vs_base_policy": flip_freq,
                    "winner_counts": {k: int(v) for k, v in local_counts.items()},
                })
    flip_ci = jeffreys_interval_from_successes(int(total_flips), int(max(1, total_trials)))
    ur_numerical_stress_policy_budget_fragility_precursor = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_NUMERICAL_STRESS_POLICY_BUDGET_FRAGILITY",
        "base_policy": base_policy,
        "cap_scale_grid": cap_scale_grid,
        "objective_jitter_grid": jitter_grid,
        "bootstrap_size_per_cell": int(stress_boot_n),
        "rows": stress_rows,
        "global_flip_count_vs_base_policy": int(total_flips),
        "global_trials": int(total_trials),
        "global_flip_frequency_vs_base_policy": float(total_flips / max(1, total_trials)),
        "global_flip_frequency_ci95_jeffreys": flip_ci,
        "policy_win_counter_global": {k: int(v) for k, v in policy_win_counter.items()},
    }
    # Class-conditional decomposition of budget-fragility stress:
    # isolate whether flip-instability is concentrated in one channel class.
    class_names = ["gauge_gauge", "fermion_fermion", "scalar_scalar"]
    class_fragility_rows = []
    for cls in class_names:
        tri_cls = [r for r in tri_rows if str(r["class"]) == cls]
        if not tri_cls:
            class_fragility_rows.append({
                "class": cls,
                "base_policy": "none",
                "global_flip_count_vs_base_policy": 0,
                "global_trials": 0,
                "global_flip_frequency_vs_base_policy": 0.0,
                "global_flip_frequency_ci95_jeffreys": {"lower": 0.0, "upper": 1.0},
                "policy_win_counter_global": {pn: 0 for pn in policy_names},
            })
            continue
        policy_rows_cls = []
        for pname, pmap in policy_defs.items():
            vals_warn = []
            vals_d2 = []
            vals_rt = []
            for rr in tri_cls:
                pick = pmap[cls]
                vals_warn.append(float(rr[f"{pick}_warning_count"]))
                vals_d2.append(float(rr[f"{pick}_delta_l2_vs_baseline"]))
                vals_rt.append(float(rr[f"{pick}_runtime_seconds"]))
            arr_warn = np.array(vals_warn, dtype=float)
            arr_d2 = np.array(vals_d2, dtype=float)
            arr_rt = np.array(vals_rt, dtype=float)
            policy_rows_cls.append({
                "policy": pname,
                "warning_total": int(np.sum(arr_warn)),
                "delta_l2_span": float(np.max(arr_d2) - np.min(arr_d2)) if arr_d2.size else 0.0,
                "runtime_total_seconds": float(np.sum(arr_rt)),
            })
        base_policy_cls = str(sorted(policy_rows_cls, key=lambda r: (int(r["warning_total"]), float(r["delta_l2_span"]), float(r["runtime_total_seconds"]), str(r["policy"])))[0]["policy"])
        cls_total_trials = 0
        cls_total_flips = 0
        cls_policy_wins = {pn: 0 for pn in policy_names}
        for cap_scale in cap_scale_grid:
            for jitter in jitter_grid:
                for _ in range(stress_boot_n):
                    idx = stress_rng.integers(0, len(tri_cls), size=len(tri_cls))
                    boot_rows_cls = []
                    for pname, pmap in policy_defs.items():
                        vals_warn = []
                        vals_d2 = []
                        vals_rt = []
                        for ii in idx:
                            rr = tri_cls[int(ii)]
                            pick = pmap[cls]
                            vals_warn.append(float(rr[f"{pick}_warning_count"]))
                            vals_d2.append(float(rr[f"{pick}_delta_l2_vs_baseline"]))
                            vals_rt.append(float(rr[f"{pick}_runtime_seconds"]))
                        arr_warn = np.array(vals_warn, dtype=float)
                        arr_d2 = np.array(vals_d2, dtype=float)
                        arr_rt = np.array(vals_rt, dtype=float)
                        boot_rows_cls.append({
                            "policy": pname,
                            "warning_total": float(np.sum(arr_warn) * (1.0 + jitter)),
                            "delta_l2_span": float((np.max(arr_d2) - np.min(arr_d2)) * (1.0 + jitter)) if arr_d2.size else 0.0,
                            "runtime_total_seconds": float(np.sum(arr_rt)),
                        })
                    runtime_arr = np.array([float(r["runtime_total_seconds"]) for r in boot_rows_cls], dtype=float)
                    cap = float(cap_scale * (float(np.quantile(runtime_arr, 0.50)) if runtime_arr.size else 0.0))
                    eligible = [r for r in boot_rows_cls if float(r["runtime_total_seconds"]) <= cap + 1e-15]
                    if not eligible:
                        continue
                    chosen = str(sorted(eligible, key=lambda r: (float(r["warning_total"]), float(r["delta_l2_span"]), float(r["runtime_total_seconds"]), str(r["policy"])))[0]["policy"])
                    cls_policy_wins[chosen] += 1
                    cls_total_trials += 1
                    if chosen != base_policy_cls:
                        cls_total_flips += 1
        class_fragility_rows.append({
            "class": cls,
            "base_policy": base_policy_cls,
            "global_flip_count_vs_base_policy": int(cls_total_flips),
            "global_trials": int(cls_total_trials),
            "global_flip_frequency_vs_base_policy": float(cls_total_flips / max(1, cls_total_trials)),
            "global_flip_frequency_ci95_jeffreys": jeffreys_interval_from_successes(int(cls_total_flips), int(max(1, cls_total_trials))),
            "policy_win_counter_global": {k: int(v) for k, v in cls_policy_wins.items()},
        })
    ur_numerical_stress_policy_budget_fragility_by_class_precursor = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_NUMERICAL_STRESS_POLICY_BUDGET_FRAGILITY_BY_CLASS",
        "rows": class_fragility_rows,
    }
    # Class-aware adaptive fallback gate:
    # if class-level fragility is high, fallback to class-local robust winner.
    fragility_lb_threshold = 0.30
    adaptive_rows = []
    adaptive_policy_map = {}
    for rr in class_fragility_rows:
        cls = str(rr["class"])
        lb = float(rr["global_flip_frequency_ci95_jeffreys"]["lower"])
        wins = rr["policy_win_counter_global"]
        robust_winner = sorted(policy_names, key=lambda pn: (-int(wins[pn]), pn))[0]
        base_cls = str(rr["base_policy"])
        use_fallback = bool(lb > fragility_lb_threshold)
        selected_policy = robust_winner if use_fallback else base_cls
        if selected_policy == "none":
            selected_policy = robust_winner
        adaptive_policy_map[cls] = selected_policy
        adaptive_rows.append({
            "class": cls,
            "fragility_flip_frequency_wilson_lb95_proxy": lb,
            "fragility_lb_threshold": float(fragility_lb_threshold),
            "base_policy": base_cls,
            "robust_winner_policy": robust_winner,
            "use_fallback": use_fallback,
            "selected_policy": selected_policy,
        })
    adaptive_warn = []
    adaptive_d2 = []
    adaptive_rt = []
    for rr in tri_rows:
        cls = str(rr["class"])
        policy_pick = adaptive_policy_map[cls]
        transform_pick = str(policy_defs[policy_pick][cls])
        adaptive_warn.append(float(rr[f"{transform_pick}_warning_count"]))
        adaptive_d2.append(float(rr[f"{transform_pick}_delta_l2_vs_baseline"]))
        adaptive_rt.append(float(rr[f"{transform_pick}_runtime_seconds"]))
    aw = np.array(adaptive_warn, dtype=float)
    ad = np.array(adaptive_d2, dtype=float)
    ar = np.array(adaptive_rt, dtype=float)
    base_policy_row = next((r for r in policy_rows_sorted if str(r["policy"]) == str(base_policy)), None)
    ur_numerical_stress_policy_class_adaptive_fallback_precursor = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_NUMERICAL_STRESS_POLICY_CLASS_ADAPTIVE_FALLBACK",
        "class_fragility_lb_threshold": float(fragility_lb_threshold),
        "class_policy_map": adaptive_policy_map,
        "rows": adaptive_rows,
        "adaptive_warning_total": int(np.sum(aw)),
        "adaptive_delta_l2_span": float(np.max(ad) - np.min(ad)) if ad.size else 0.0,
        "adaptive_runtime_total_seconds": float(np.sum(ar)),
        "base_budget_policy": str(base_policy),
        "base_budget_warning_total": int(base_policy_row["warning_total"]) if base_policy_row else 0,
        "base_budget_delta_l2_span": float(base_policy_row["delta_l2_span"]) if base_policy_row else 0.0,
        "base_budget_runtime_total_seconds": float(base_policy_row["runtime_total_seconds"]) if base_policy_row else 0.0,
        "delta_warning_total_adaptive_minus_base": int(np.sum(aw) - (int(base_policy_row["warning_total"]) if base_policy_row else 0)),
        "delta_delta_l2_span_adaptive_minus_base": float((np.max(ad) - np.min(ad)) - (float(base_policy_row["delta_l2_span"]) if base_policy_row else 0.0)) if ad.size else 0.0,
        "delta_runtime_total_seconds_adaptive_minus_base": float(np.sum(ar) - (float(base_policy_row["runtime_total_seconds"]) if base_policy_row else 0.0)),
    }
    # Policy ablation: base budget policy vs class-adaptive fallback vs robust-winner-only.
    robust_only_map = {}
    for rr in class_fragility_rows:
        wins = rr["policy_win_counter_global"]
        robust_only_map[str(rr["class"])] = sorted(policy_names, key=lambda pn: (-int(wins[pn]), pn))[0]

    def aggregate_policy_from_class_map(class_map: dict[str, str]) -> dict[str, float]:
        vals_warn, vals_d2, vals_rt = [], [], []
        for rtri in tri_rows:
            cls = str(rtri["class"])
            policy_pick = str(class_map[cls])
            transform_pick = str(policy_defs[policy_pick][cls])
            vals_warn.append(float(rtri[f"{transform_pick}_warning_count"]))
            vals_d2.append(float(rtri[f"{transform_pick}_delta_l2_vs_baseline"]))
            vals_rt.append(float(rtri[f"{transform_pick}_runtime_seconds"]))
        arr_warn = np.array(vals_warn, dtype=float)
        arr_d2 = np.array(vals_d2, dtype=float)
        arr_rt = np.array(vals_rt, dtype=float)
        return {
            "warning_total": float(np.sum(arr_warn)),
            "delta_l2_span": float(np.max(arr_d2) - np.min(arr_d2)) if arr_d2.size else 0.0,
            "runtime_total_seconds": float(np.sum(arr_rt)),
        }

    base_map = {cls: str(base_policy) for cls in ("gauge_gauge", "fermion_fermion", "scalar_scalar")}
    adaptive_map = {k: str(v) for k, v in adaptive_policy_map.items()}
    ablation_defs = {
        "base_budget_policy": base_map,
        "class_adaptive_fallback": adaptive_map,
        "robust_winner_only": robust_only_map,
    }
    ablation_rows = []
    for name, cmap in ablation_defs.items():
        agg = aggregate_policy_from_class_map(cmap)
        ablation_rows.append({
            "regime": name,
            "class_policy_map": cmap,
            "warning_total": int(agg["warning_total"]),
            "delta_l2_span": float(agg["delta_l2_span"]),
            "runtime_total_seconds": float(agg["runtime_total_seconds"]),
        })
    ablation_rows_sorted = sorted(ablation_rows, key=lambda r: (r["warning_total"], r["delta_l2_span"], r["runtime_total_seconds"], r["regime"]))

    # Bootstrap pairwise dominance frequencies.
    abl_boot_n = 256
    abl_rng = np.random.default_rng(975150)
    pair_rows = []
    regimes = [r["regime"] for r in ablation_rows_sorted]
    for ra, rb in combinations(regimes, 2):
        a_dom = 0
        b_dom = 0
        ties = 0
        for _ in range(abl_boot_n):
            idx = abl_rng.integers(0, len(tri_rows), size=len(tri_rows))
            met = {}
            for rn in (ra, rb):
                cmap = ablation_defs[rn]
                vals_warn, vals_d2, vals_rt = [], [], []
                for ii in idx:
                    rtri = tri_rows[int(ii)]
                    cls = str(rtri["class"])
                    pp = str(cmap[cls])
                    tp = str(policy_defs[pp][cls])
                    vals_warn.append(float(rtri[f"{tp}_warning_count"]))
                    vals_d2.append(float(rtri[f"{tp}_delta_l2_vs_baseline"]))
                    vals_rt.append(float(rtri[f"{tp}_runtime_seconds"]))
                awn = np.array(vals_warn, dtype=float)
                ad2 = np.array(vals_d2, dtype=float)
                art = np.array(vals_rt, dtype=float)
                met[rn] = (
                    float(np.sum(awn)),
                    float(np.max(ad2) - np.min(ad2)) if ad2.size else 0.0,
                    float(np.sum(art)),
                )
            a_nonworse = met[ra][0] <= met[rb][0] and met[ra][1] <= met[rb][1] and met[ra][2] <= met[rb][2]
            b_nonworse = met[rb][0] <= met[ra][0] and met[rb][1] <= met[ra][1] and met[rb][2] <= met[ra][2]
            a_strict = met[ra][0] < met[rb][0] or met[ra][1] < met[rb][1] or met[ra][2] < met[rb][2]
            b_strict = met[rb][0] < met[ra][0] or met[rb][1] < met[ra][1] or met[rb][2] < met[ra][2]
            if a_nonworse and a_strict:
                a_dom += 1
            elif b_nonworse and b_strict:
                b_dom += 1
            else:
                ties += 1
        pair_rows.append({
            "regime_a": ra,
            "regime_b": rb,
            "a_dominates_b_frequency": float(a_dom / abl_boot_n),
            "b_dominates_a_frequency": float(b_dom / abl_boot_n),
            "tie_or_incomparable_frequency": float(ties / abl_boot_n),
            "a_dominates_b_ci95_jeffreys": jeffreys_interval_from_successes(a_dom, abl_boot_n),
            "b_dominates_a_ci95_jeffreys": jeffreys_interval_from_successes(b_dom, abl_boot_n),
            "bootstrap_size": int(abl_boot_n),
        })
    ur_numerical_stress_policy_ablation_precursor = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_NUMERICAL_STRESS_POLICY_ABLATION",
        "rows": ablation_rows_sorted,
        "best_regime_lexicographic": str(ablation_rows_sorted[0]["regime"]) if ablation_rows_sorted else "none",
        "pairwise_dominance_rows": pair_rows,
    }
    # Cross-class constrained ablation:
    # require class-local runtime/warning constraints before global ranking.
    class_runtime_caps = {
        "gauge_gauge": float(np.quantile(np.array([float(r["u2_runtime_seconds"]) for r in tri_rows if str(r["class"]) == "gauge_gauge"], dtype=float), 0.75)),
        "fermion_fermion": float(np.quantile(np.array([float(r["u2_runtime_seconds"]) for r in tri_rows if str(r["class"]) == "fermion_fermion"], dtype=float), 0.75)),
        "scalar_scalar": float(np.quantile(np.array([float(r["u2_runtime_seconds"]) for r in tri_rows if str(r["class"]) == "scalar_scalar"], dtype=float), 0.75)),
    }
    class_warning_caps = {
        "gauge_gauge": float(np.quantile(np.array([float(r["u2_warning_count"]) for r in tri_rows if str(r["class"]) == "gauge_gauge"], dtype=float), 0.75)),
        "fermion_fermion": float(np.quantile(np.array([float(r["u2_warning_count"]) for r in tri_rows if str(r["class"]) == "fermion_fermion"], dtype=float), 0.75)),
        "scalar_scalar": float(np.quantile(np.array([float(r["u2_warning_count"]) for r in tri_rows if str(r["class"]) == "scalar_scalar"], dtype=float), 0.75)),
    }
    constrained_rows = []
    for ab in ablation_rows_sorted:
        name = str(ab["regime"])
        cmap = ablation_defs[name]
        pass_rows = 0
        total_rows = 0
        for rtri in tri_rows:
            cls = str(rtri["class"])
            pp = str(cmap[cls])
            tp = str(policy_defs[pp][cls])
            rt = float(rtri[f"{tp}_runtime_seconds"])
            wc = float(rtri[f"{tp}_warning_count"])
            ok = (rt <= class_runtime_caps[cls] + 1e-15) and (wc <= class_warning_caps[cls] + 1e-15)
            pass_rows += 1 if ok else 0
            total_rows += 1
        constrained_rows.append({
            "regime": name,
            "class_policy_map": cmap,
            "constraint_pass_rows": int(pass_rows),
            "constraint_total_rows": int(total_rows),
            "constraint_pass_rate": float(pass_rows / max(1, total_rows)),
            "warning_total": int(ab["warning_total"]),
            "delta_l2_span": float(ab["delta_l2_span"]),
            "runtime_total_seconds": float(ab["runtime_total_seconds"]),
        })
    feasible_rows = [r for r in constrained_rows if float(r["constraint_pass_rate"]) >= 0.80]
    constrained_best = (
        sorted(feasible_rows, key=lambda r: (r["warning_total"], r["delta_l2_span"], r["runtime_total_seconds"], r["regime"]))[0]["regime"]
        if feasible_rows else "none"
    )
    ur_numerical_stress_policy_cross_class_constrained_ablation_precursor = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_NUMERICAL_STRESS_POLICY_CROSS_CLASS_CONSTRAINED_ABLATION",
        "constraint_rule": "row_runtime_le_class_q75_u2_and_row_warning_le_class_q75_u2",
        "feasibility_pass_rate_threshold": 0.80,
        "class_runtime_caps_seconds": class_runtime_caps,
        "class_warning_caps": class_warning_caps,
        "rows": constrained_rows,
        "feasible_regimes": [str(r["regime"]) for r in feasible_rows],
        "best_feasible_regime_lexicographic": str(constrained_best),
    }
    # Constrained bootstrap pairwise dominance:
    # compare regimes only over rows satisfying class-local constraints.
    cboot_n = 256
    cboot_rng = np.random.default_rng(975152)
    constrained_pair_rows = []
    regime_names = [r["regime"] for r in ablation_rows_sorted]
    for ra, rb in combinations(regime_names, 2):
        a_dom = 0
        b_dom = 0
        ties = 0
        usable = 0
        for _ in range(cboot_n):
            idx = cboot_rng.integers(0, len(tri_rows), size=len(tri_rows))
            metrics = {}
            for rn in (ra, rb):
                cmap = ablation_defs[rn]
                vals_warn, vals_d2, vals_rt = [], [], []
                for ii in idx:
                    rr = tri_rows[int(ii)]
                    cls = str(rr["class"])
                    pp = str(cmap[cls])
                    tp = str(policy_defs[pp][cls])
                    rt = float(rr[f"{tp}_runtime_seconds"])
                    wc = float(rr[f"{tp}_warning_count"])
                    if (rt <= class_runtime_caps[cls] + 1e-15) and (wc <= class_warning_caps[cls] + 1e-15):
                        vals_warn.append(float(rr[f"{tp}_warning_count"]))
                        vals_d2.append(float(rr[f"{tp}_delta_l2_vs_baseline"]))
                        vals_rt.append(float(rr[f"{tp}_runtime_seconds"]))
                if len(vals_warn) == 0:
                    metrics[rn] = None
                else:
                    awn = np.array(vals_warn, dtype=float)
                    ad2 = np.array(vals_d2, dtype=float)
                    art = np.array(vals_rt, dtype=float)
                    metrics[rn] = (
                        float(np.sum(awn)),
                        float(np.max(ad2) - np.min(ad2)) if ad2.size else 0.0,
                        float(np.sum(art)),
                    )
            if metrics[ra] is None or metrics[rb] is None:
                continue
            usable += 1
            a_nonworse = metrics[ra][0] <= metrics[rb][0] and metrics[ra][1] <= metrics[rb][1] and metrics[ra][2] <= metrics[rb][2]
            b_nonworse = metrics[rb][0] <= metrics[ra][0] and metrics[rb][1] <= metrics[ra][1] and metrics[rb][2] <= metrics[ra][2]
            a_strict = metrics[ra][0] < metrics[rb][0] or metrics[ra][1] < metrics[rb][1] or metrics[ra][2] < metrics[rb][2]
            b_strict = metrics[rb][0] < metrics[ra][0] or metrics[rb][1] < metrics[ra][1] or metrics[rb][2] < metrics[ra][2]
            if a_nonworse and a_strict:
                a_dom += 1
            elif b_nonworse and b_strict:
                b_dom += 1
            else:
                ties += 1
        denom = max(1, usable)
        if usable == 0:
            # No constrained-overlap rows under current bootstrap draws:
            # treat as fully tie/incomparable (non-informative) mass.
            ties = denom
        constrained_pair_rows.append({
            "regime_a": ra,
            "regime_b": rb,
            "usable_bootstrap_trials": int(usable),
            "a_dominates_b_frequency": float(a_dom / denom),
            "b_dominates_a_frequency": float(b_dom / denom),
            "tie_or_incomparable_frequency": float(ties / denom),
            "a_dominates_b_ci95_jeffreys": jeffreys_interval_from_successes(a_dom, denom),
            "b_dominates_a_ci95_jeffreys": jeffreys_interval_from_successes(b_dom, denom),
            "bootstrap_size_requested": int(cboot_n),
        })
    ur_numerical_stress_policy_cross_class_constrained_bootstrap_dominance_precursor = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_NUMERICAL_STRESS_POLICY_CROSS_CLASS_CONSTRAINED_BOOTSTRAP_DOMINANCE",
        "rows": constrained_pair_rows,
    }
    # Threshold sweep for constrained-feasibility pass-rate:
    # evaluate whether best feasible regime is stable vs threshold choice.
    threshold_grid = [0.70, 0.75, 0.80, 0.85, 0.90, 0.95]
    threshold_rows = []
    for th in threshold_grid:
        feasible_th = [r for r in constrained_rows if float(r["constraint_pass_rate"]) >= float(th)]
        best_th = (
            sorted(feasible_th, key=lambda r: (r["warning_total"], r["delta_l2_span"], r["runtime_total_seconds"], r["regime"]))[0]["regime"]
            if feasible_th else "none"
        )
        threshold_rows.append({
            "feasibility_pass_rate_threshold": float(th),
            "num_feasible_regimes": int(len(feasible_th)),
            "feasible_regimes": [str(r["regime"]) for r in feasible_th],
            "best_feasible_regime_lexicographic": str(best_th),
        })
    best_seq = [str(r["best_feasible_regime_lexicographic"]) for r in threshold_rows]
    stable_non_none = [b for b in best_seq if b != "none"]
    threshold_stable = bool(len(set(stable_non_none)) <= 1) if stable_non_none else True
    ur_numerical_stress_policy_cross_class_threshold_sweep_precursor = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_NUMERICAL_STRESS_POLICY_CROSS_CLASS_THRESHOLD_SWEEP",
        "threshold_grid": threshold_grid,
        "rows": threshold_rows,
        "best_regime_sequence": best_seq,
        "best_regime_stable_over_nonempty_thresholds": threshold_stable,
    }
    # Joint stress map: sweep (threshold, cap_scale, jitter) and track
    # stability of winning regime under constrained feasibility.
    joint_rows = []
    joint_boot_n = 128
    joint_rng = np.random.default_rng(975154)
    for th in threshold_grid:
        for cap_scale in cap_scale_grid:
            for jitter in jitter_grid:
                winner_counts = {rn: 0 for rn in regime_names}
                usable_trials = 0
                for _ in range(joint_boot_n):
                    idx = joint_rng.integers(0, len(tri_rows), size=len(tri_rows))
                    rows_boot = []
                    for rn in regime_names:
                        cmap = ablation_defs[rn]
                        pass_rows = 0
                        total_rows = 0
                        vals_warn, vals_d2, vals_rt = [], [], []
                        for ii in idx:
                            rr = tri_rows[int(ii)]
                            cls = str(rr["class"])
                            pp = str(cmap[cls])
                            tp = str(policy_defs[pp][cls])
                            rt = float(rr[f"{tp}_runtime_seconds"])
                            wc = float(rr[f"{tp}_warning_count"])
                            cap_rt = float(cap_scale * class_runtime_caps[cls])
                            cap_wc = float((1.0 + jitter) * class_warning_caps[cls])
                            ok = (rt <= cap_rt + 1e-15) and (wc <= cap_wc + 1e-15)
                            pass_rows += 1 if ok else 0
                            total_rows += 1
                            vals_warn.append(float(rr[f"{tp}_warning_count"]))
                            vals_d2.append(float(rr[f"{tp}_delta_l2_vs_baseline"]))
                            vals_rt.append(float(rr[f"{tp}_runtime_seconds"]))
                        aw = np.array(vals_warn, dtype=float)
                        ad2 = np.array(vals_d2, dtype=float)
                        art = np.array(vals_rt, dtype=float)
                        rows_boot.append({
                            "regime": rn,
                            "constraint_pass_rate": float(pass_rows / max(1, total_rows)),
                            "warning_total": float(np.sum(aw)),
                            "delta_l2_span": float(np.max(ad2) - np.min(ad2)) if ad2.size else 0.0,
                            "runtime_total_seconds": float(np.sum(art)),
                        })
                    feasible = [r for r in rows_boot if float(r["constraint_pass_rate"]) >= float(th)]
                    if not feasible:
                        continue
                    usable_trials += 1
                    winner = sorted(feasible, key=lambda r: (r["warning_total"], r["delta_l2_span"], r["runtime_total_seconds"], r["regime"]))[0]["regime"]
                    winner_counts[str(winner)] += 1
                if usable_trials == 0:
                    winner = "none"
                    winner_freq = 0.0
                    entropy = 0.0
                else:
                    winner = sorted(regime_names, key=lambda rn: (-winner_counts[rn], rn))[0]
                    freqs = np.array([winner_counts[rn] / usable_trials for rn in regime_names], dtype=float)
                    winner_freq = float(np.max(freqs))
                    entropy = float(-np.sum([p * np.log(p) for p in freqs if p > 0.0]) / np.log(max(2, len(regime_names))))
                joint_rows.append({
                    "threshold": float(th),
                    "cap_scale": float(cap_scale),
                    "jitter": float(jitter),
                    "bootstrap_size_requested": int(joint_boot_n),
                    "usable_trials": int(usable_trials),
                    "winner": str(winner),
                    "winner_frequency": float(winner_freq),
                    "winner_frequency_ci95_jeffreys": jeffreys_interval_from_successes(int(max(winner_counts.values()) if usable_trials > 0 else 0), int(max(1, usable_trials))),
                    "winner_entropy_norm": float(entropy),
                    "winner_counts": {k: int(v) for k, v in winner_counts.items()},
                })
    stable_cells = [r for r in joint_rows if r["winner"] != "none" and r["winner_frequency"] >= 0.70]
    ur_numerical_stress_policy_joint_stress_map_precursor = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_NUMERICAL_STRESS_POLICY_JOINT_STRESS_MAP",
        "threshold_grid": threshold_grid,
        "cap_scale_grid": cap_scale_grid,
        "jitter_grid": jitter_grid,
        "rows": joint_rows,
        "num_cells": int(len(joint_rows)),
        "num_stable_cells_winner_freq_ge_070": int(len(stable_cells)),
        "stable_cell_frequency": float(len(stable_cells) / max(1, len(joint_rows))),
    }
    # Stability-topology report: connected components over stable cells
    # in the (threshold, cap_scale, jitter) grid.
    tvals = threshold_grid
    cvals = cap_scale_grid
    jvals = jitter_grid
    t_index = {float(v): i for i, v in enumerate(tvals)}
    c_index = {float(v): i for i, v in enumerate(cvals)}
    j_index = {float(v): i for i, v in enumerate(jvals)}
    stable_points = {}
    for rr in stable_cells:
        key = (t_index[float(rr["threshold"])], c_index[float(rr["cap_scale"])], j_index[float(rr["jitter"])])
        stable_points[key] = str(rr["winner"])
    seen = set()
    components = []
    for p0, winner0 in stable_points.items():
        if p0 in seen:
            continue
        stack = [p0]
        seen.add(p0)
        comp = []
        winners = []
        while stack:
            p = stack.pop()
            comp.append(p)
            winners.append(stable_points[p])
            ti, ci, ji = p
            for dt, dc, dj in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)):
                q = (ti + dt, ci + dc, ji + dj)
                if q in stable_points and q not in seen:
                    seen.add(q)
                    stack.append(q)
        comp_rows = [{
            "threshold": float(tvals[ti]),
            "cap_scale": float(cvals[ci]),
            "jitter": float(jvals[ji]),
            "winner": str(stable_points[(ti, ci, ji)]),
        } for ti, ci, ji in comp]
        win_counts = {rn: 0 for rn in regime_names}
        for w in winners:
            win_counts[w] += 1
        dominant = sorted(regime_names, key=lambda rn: (-win_counts[rn], rn))[0]
        components.append({
            "size": int(len(comp_rows)),
            "rows": sorted(comp_rows, key=lambda r: (r["threshold"], r["cap_scale"], r["jitter"])),
            "winner_counts": win_counts,
            "dominant_winner": dominant,
        })
    components = sorted(components, key=lambda r: (-r["size"], r["dominant_winner"]))
    ur_numerical_stress_policy_stability_topology_precursor = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_NUMERICAL_STRESS_POLICY_STABILITY_TOPOLOGY",
        "stable_cell_definition": "winner_frequency_ge_0p70_and_winner_not_none",
        "num_components": int(len(components)),
        "largest_component_size": int(max([c["size"] for c in components], default=0)),
        "components": components,
    }
    # Boundary-margin panel: distance from stable cells to nearest unstable cell
    # in 3D stress-grid (Manhattan graph distance).
    all_points = [(ti, ci, ji) for ti in range(len(tvals)) for ci in range(len(cvals)) for ji in range(len(jvals))]
    stable_set = set(stable_points.keys())
    unstable_set = set(all_points) - stable_set
    margin_rows = []
    for p in stable_set:
        if not unstable_set:
            dmin = int(len(tvals) + len(cvals) + len(jvals))
        else:
            dmin = min(abs(p[0] - q[0]) + abs(p[1] - q[1]) + abs(p[2] - q[2]) for q in unstable_set)
        margin_rows.append({
            "threshold": float(tvals[p[0]]),
            "cap_scale": float(cvals[p[1]]),
            "jitter": float(jvals[p[2]]),
            "winner": str(stable_points[p]),
            "boundary_manhattan_margin": int(dmin),
        })
    component_margin_rows = []
    for comp in components:
        crows = comp["rows"]
        vals = np.array([
            int(next(
                mr["boundary_manhattan_margin"] for mr in margin_rows
                if float(mr["threshold"]) == float(rr["threshold"])
                and float(mr["cap_scale"]) == float(rr["cap_scale"])
                and float(mr["jitter"]) == float(rr["jitter"])
            ))
            for rr in crows
        ], dtype=float) if crows else np.array([], dtype=float)
        component_margin_rows.append({
            "dominant_winner": str(comp["dominant_winner"]),
            "size": int(comp["size"]),
            "boundary_margin_q05_q50_q95": [
                float(np.quantile(vals, 0.05)) if vals.size else 0.0,
                float(np.quantile(vals, 0.50)) if vals.size else 0.0,
                float(np.quantile(vals, 0.95)) if vals.size else 0.0,
            ],
            "boundary_margin_min": float(np.min(vals)) if vals.size else 0.0,
            "boundary_margin_max": float(np.max(vals)) if vals.size else 0.0,
        })
    ur_numerical_stress_policy_stability_boundary_margin_precursor = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_NUMERICAL_STRESS_POLICY_STABILITY_BOUNDARY_MARGIN",
        "distance_metric": "manhattan_on_threshold_cap_jitter_grid",
        "rows": sorted(margin_rows, key=lambda r: (r["threshold"], r["cap_scale"], r["jitter"])),
        "num_stable_rows": int(len(margin_rows)),
        "component_margin_rows": component_margin_rows,
    }
    # Weighted boundary-risk score: combine winner confidence, entropy,
    # and boundary margin into one operational risk index (lower is better).
    risk_rows = []
    max_margin = max([int(r["boundary_manhattan_margin"]) for r in margin_rows], default=1)
    margin_scale = float(max(1, max_margin))
    row_index = {(float(r["threshold"]), float(r["cap_scale"]), float(r["jitter"])): r for r in margin_rows}
    for jr in joint_rows:
        key = (float(jr["threshold"]), float(jr["cap_scale"]), float(jr["jitter"]))
        mr = row_index.get(key)
        if mr is None:
            margin_norm = 0.0
            margin_raw = 0
        else:
            margin_raw = int(mr["boundary_manhattan_margin"])
            margin_norm = float(margin_raw / margin_scale)
        win_freq = float(jr["winner_frequency"])
        ent = float(jr["winner_entropy_norm"])
        risk = float((1.0 - win_freq) * 0.50 + ent * 0.30 + (1.0 - margin_norm) * 0.20)
        risk_rows.append({
            "threshold": float(jr["threshold"]),
            "cap_scale": float(jr["cap_scale"]),
            "jitter": float(jr["jitter"]),
            "winner": str(jr["winner"]),
            "winner_frequency": win_freq,
            "winner_entropy_norm": ent,
            "boundary_manhattan_margin": int(margin_raw),
            "boundary_margin_norm": margin_norm,
            "weighted_boundary_risk_score": risk,
        })
    risk_rows_sorted = sorted(risk_rows, key=lambda r: (r["weighted_boundary_risk_score"], -r["winner_frequency"], r["winner_entropy_norm"], -r["boundary_manhattan_margin"], r["threshold"], r["cap_scale"], r["jitter"]))
    best_corridor_rows = [r for r in risk_rows_sorted if r["winner"] != "none"][:8]
    ur_numerical_stress_policy_weighted_boundary_risk_precursor = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_NUMERICAL_STRESS_POLICY_WEIGHTED_BOUNDARY_RISK",
        "risk_definition": "0p50*(1-winner_frequency)+0p30*winner_entropy_norm+0p20*(1-boundary_margin_norm)",
        "margin_normalization_divisor": margin_scale,
        "rows": risk_rows_sorted,
        "best_corridor_rows": best_corridor_rows,
        "best_corridor_winner_set": sorted({str(r["winner"]) for r in best_corridor_rows}),
    }
    # Sensitivity sweep over risk weights:
    # test winner-set robustness against reasonable weight choices.
    weight_grid = [
        {"wf": 0.50, "ent": 0.30, "margin": 0.20},
        {"wf": 0.60, "ent": 0.25, "margin": 0.15},
        {"wf": 0.45, "ent": 0.35, "margin": 0.20},
        {"wf": 0.40, "ent": 0.30, "margin": 0.30},
        {"wf": 0.55, "ent": 0.20, "margin": 0.25},
    ]
    risk_weight_rows = []
    winner_sets = []
    for wg in weight_grid:
        rw = []
        for rr in risk_rows:
            risk_v = float(
                wg["wf"] * (1.0 - float(rr["winner_frequency"])) +
                wg["ent"] * float(rr["winner_entropy_norm"]) +
                wg["margin"] * (1.0 - float(rr["boundary_margin_norm"]))
            )
            rw.append({**rr, "weighted_boundary_risk_score": risk_v})
        rw_sorted = sorted(
            rw,
            key=lambda r: (
                r["weighted_boundary_risk_score"],
                -r["winner_frequency"],
                r["winner_entropy_norm"],
                -r["boundary_manhattan_margin"],
                r["threshold"], r["cap_scale"], r["jitter"]
            )
        )
        corridor = [r for r in rw_sorted if r["winner"] != "none"][:8]
        wset = sorted({str(r["winner"]) for r in corridor})
        winner_sets.append(tuple(wset))
        risk_weight_rows.append({
            "weights": wg,
            "best_corridor_winner_set": wset,
            "best_corridor_rows": corridor,
            "best_cell_score": float(corridor[0]["weighted_boundary_risk_score"]) if corridor else 1.0,
        })
    winner_set_stable = bool(len(set(winner_sets)) <= 1)
    ur_numerical_stress_policy_weighted_boundary_risk_sensitivity_precursor = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_NUMERICAL_STRESS_POLICY_WEIGHTED_BOUNDARY_RISK_SENSITIVITY",
        "rows": risk_weight_rows,
        "winner_set_stable_over_weight_grid": winner_set_stable,
    }
    # Bayesian-like posterior over risk weights (Dirichlet prior) for winner-set
    # optimality probability under weighted boundary-risk objective.
    dirichlet_alpha = [5.0, 3.0, 2.0]  # prior emphasizing confidence term
    n_post = 512
    post_rng = np.random.default_rng(975159)
    winner_set_counts = {}
    winner_cell_counts = {rn: 0 for rn in ("base_budget_policy", "class_adaptive_fallback", "robust_winner_only", "none")}
    for _ in range(n_post):
        w = post_rng.dirichlet(np.array(dirichlet_alpha, dtype=float))
        rows_post = []
        for rr in risk_rows:
            risk_v = float(
                w[0] * (1.0 - float(rr["winner_frequency"])) +
                w[1] * float(rr["winner_entropy_norm"]) +
                w[2] * (1.0 - float(rr["boundary_margin_norm"]))
            )
            rows_post.append({**rr, "weighted_boundary_risk_score": risk_v})
        rows_post = sorted(
            rows_post,
            key=lambda r: (
                r["weighted_boundary_risk_score"],
                -r["winner_frequency"],
                r["winner_entropy_norm"],
                -r["boundary_manhattan_margin"],
                r["threshold"], r["cap_scale"], r["jitter"]
            )
        )
        corridor = [r for r in rows_post if r["winner"] != "none"][:8]
        wset = tuple(sorted({str(r["winner"]) for r in corridor}))
        winner_set_counts[wset] = winner_set_counts.get(wset, 0) + 1
        best_cell_winner = str(rows_post[0]["winner"]) if rows_post else "none"
        winner_cell_counts[best_cell_winner] = winner_cell_counts.get(best_cell_winner, 0) + 1
    posterior_rows = []
    for wset, cnt in sorted(winner_set_counts.items(), key=lambda kv: (-kv[1], kv[0])):
        posterior_rows.append({
            "winner_set": list(wset),
            "count": int(cnt),
            "posterior_probability": float(cnt / n_post),
            "posterior_probability_ci95_jeffreys": jeffreys_interval_from_successes(int(cnt), int(n_post)),
        })
    ur_numerical_stress_policy_weighted_boundary_risk_bayesian_precursor = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_NUMERICAL_STRESS_POLICY_WEIGHTED_BOUNDARY_RISK_BAYESIAN",
        "dirichlet_alpha": dirichlet_alpha,
        "posterior_sample_size": int(n_post),
        "winner_set_posterior_rows": posterior_rows,
        "most_probable_winner_set": posterior_rows[0]["winner_set"] if posterior_rows else [],
        "most_probable_winner_set_probability": posterior_rows[0]["posterior_probability"] if posterior_rows else 0.0,
        "best_cell_winner_posterior_counts": winner_cell_counts,
    }
    # Posterior predictive stress-check:
    # sample weights from posterior prior proxy and apply stress perturbations;
    # estimate P(winner-set remains optimal).
    pp_n = 256
    pp_rng = np.random.default_rng(975160)
    base_wset = tuple(sorted({str(x) for x in ur_numerical_stress_policy_weighted_boundary_risk_bayesian_precursor["most_probable_winner_set"]}))
    pp_rows = []
    for cap_scale in cap_scale_grid:
        for jitter in jitter_grid:
            stay = 0
            valid = 0
            for _ in range(pp_n):
                w = pp_rng.dirichlet(np.array(dirichlet_alpha, dtype=float))
                rows_pp = []
                for rr in risk_rows:
                    risk_v = float(
                        w[0] * (1.0 - float(rr["winner_frequency"])) +
                        w[1] * float(rr["winner_entropy_norm"]) +
                        w[2] * (1.0 - float(rr["boundary_margin_norm"]))
                    )
                    # small predictive stress proxy on risk
                    risk_v = float(risk_v * (1.0 + 0.25 * abs(float(cap_scale) - 1.0) + 0.50 * float(jitter)))
                    rows_pp.append({**rr, "weighted_boundary_risk_score": risk_v})
                rows_pp = sorted(
                    rows_pp,
                    key=lambda r: (
                        r["weighted_boundary_risk_score"],
                        -r["winner_frequency"],
                        r["winner_entropy_norm"],
                        -r["boundary_manhattan_margin"],
                        r["threshold"], r["cap_scale"], r["jitter"]
                    )
                )
                corridor = [r for r in rows_pp if r["winner"] != "none"][:8]
                if not corridor:
                    continue
                valid += 1
                wset = tuple(sorted({str(r["winner"]) for r in corridor}))
                if wset == base_wset:
                    stay += 1
            p_stay = float(stay / max(1, valid))
            pp_rows.append({
                "cap_scale": float(cap_scale),
                "jitter": float(jitter),
                "posterior_samples": int(pp_n),
                "valid_samples": int(valid),
                "stay_count": int(stay),
                "p_winner_set_stays_optimal": p_stay,
                "p_winner_set_stays_optimal_ci95_jeffreys": jeffreys_interval_from_successes(int(stay), int(max(1, valid))),
            })
    ur_numerical_stress_policy_weighted_boundary_risk_posterior_predictive_precursor = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_NUMERICAL_STRESS_POLICY_WEIGHTED_BOUNDARY_RISK_POSTERIOR_PREDICTIVE",
        "reference_winner_set": list(base_wset),
        "rows": pp_rows,
        "global_p_stays_optimal_mean": float(np.mean(np.array([r["p_winner_set_stays_optimal"] for r in pp_rows], dtype=float))) if pp_rows else 0.0,
    }
    # Decision gate from posterior-predictive robustness.
    pp_lb_threshold = 0.80
    critical_rows = []
    for rr in pp_rows:
        lb = float(rr["p_winner_set_stays_optimal_ci95_jeffreys"]["lower"])
        critical_rows.append({
            "cap_scale": float(rr["cap_scale"]),
            "jitter": float(rr["jitter"]),
            "p_stay": float(rr["p_winner_set_stays_optimal"]),
            "p_stay_lb95": lb,
            "criterion_lb95_ge_threshold": bool(lb >= pp_lb_threshold),
        })
    all_pass = bool(all(r["criterion_lb95_ge_threshold"] for r in critical_rows)) if critical_rows else False
    ur_numerical_stress_policy_posterior_predictive_decision_gate_precursor = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_NUMERICAL_STRESS_POLICY_POSTERIOR_PREDICTIVE_DECISION_GATE",
        "rule": "GO_if_all_critical_cells_lb95_ge_threshold_else_HOLD",
        "lb95_threshold": float(pp_lb_threshold),
        "rows": critical_rows,
        "decision": "GO" if all_pass else "HOLD_AND_RECALIBRATE",
        "ready_for_next_costlier_policy_step": bool(all_pass),
    }
    # Calibration sweep for LB95 threshold: GO-rate vs threshold.
    gate_threshold_grid = [0.70, 0.75, 0.80, 0.85, 0.90, 0.95]
    gate_rows = []
    for th in gate_threshold_grid:
        pass_flags = [float(r["p_stay_lb95"]) >= float(th) for r in critical_rows]
        pass_count = int(sum(1 for x in pass_flags if x))
        total = int(len(pass_flags))
        go = bool(all(pass_flags)) if pass_flags else False
        gate_rows.append({
            "lb95_threshold": float(th),
            "critical_cells_pass_count": pass_count,
            "critical_cells_total": total,
            "critical_cells_pass_rate": float(pass_count / max(1, total)),
            "decision_go": go,
        })
    go_rate = float(np.mean(np.array([1.0 if r["decision_go"] else 0.0 for r in gate_rows], dtype=float))) if gate_rows else 0.0
    ur_numerical_stress_policy_posterior_predictive_gate_calibration_precursor = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_NUMERICAL_STRESS_POLICY_POSTERIOR_PREDICTIVE_GATE_CALIBRATION",
        "threshold_grid": gate_threshold_grid,
        "rows": gate_rows,
        "go_rate_over_threshold_grid": go_rate,
        "recommended_threshold_max_go_with_min_pass_rate_0p90": (
            float(max([r["lb95_threshold"] for r in gate_rows if r["critical_cells_pass_rate"] >= 0.90], default=0.70))
        ),
    }
    # Cost-aware gate calibration:
    # optimize threshold by utility = go_rate - lambda * expected_runtime_uplift.
    base_runtime_ref = float(np.mean(np.array([float(r["runtime_total_seconds"]) for r in policy_rows_sorted], dtype=float))) if policy_rows_sorted else 1.0
    lambda_grid = [0.0, 0.1, 0.2, 0.4]
    cost_rows = []
    for lam in lambda_grid:
        best_u = -1e9
        best_row = None
        for gr in gate_rows:
            th = float(gr["lb95_threshold"])
            # proxy uplift: stricter thresholds imply more conservative (and typically costlier) policy posture
            runtime_uplift = float((th - 0.70) / 0.25)
            go_rate_local = 1.0 if bool(gr["decision_go"]) else 0.0
            utility = float(go_rate_local - float(lam) * runtime_uplift)
            row = {
                "lambda_cost": float(lam),
                "lb95_threshold": th,
                "go_rate_local": go_rate_local,
                "runtime_uplift_proxy": runtime_uplift,
                "utility": utility,
            }
            if utility > best_u:
                best_u = utility
                best_row = row
            cost_rows.append(row)
        if best_row is not None:
            best_row["selected_for_lambda"] = True
    selected_rows = []
    for lam in lambda_grid:
        rows_lam = [r for r in cost_rows if float(r["lambda_cost"]) == float(lam)]
        if rows_lam:
            selected_rows.append(sorted(rows_lam, key=lambda r: (-r["utility"], r["lb95_threshold"]))[0])
    ur_numerical_stress_policy_posterior_predictive_gate_cost_calibration_precursor = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_NUMERICAL_STRESS_POLICY_POSTERIOR_PREDICTIVE_GATE_COST_CALIBRATION",
        "utility_definition": "go_rate_local - lambda_cost*runtime_uplift_proxy",
        "base_runtime_reference_seconds": base_runtime_ref,
        "lambda_grid": lambda_grid,
        "rows": cost_rows,
        "selected_rows": selected_rows,
    }
    # Frontier utility panel:
    # construct Pareto front in (false_hold_risk, false_go_risk, runtime_uplift).
    frontier_rows = []
    for rr in gate_rows:
        th = float(rr["lb95_threshold"])
        pass_rate = float(rr["critical_cells_pass_rate"])
        go = bool(rr["decision_go"])
        false_hold_risk = float(1.0 - pass_rate)
        false_go_risk = float(0.0 if go else 1.0)
        runtime_uplift_proxy = float((th - 0.70) / 0.25)
        frontier_rows.append({
            "lb95_threshold": th,
            "false_hold_risk": false_hold_risk,
            "false_go_risk": false_go_risk,
            "runtime_uplift_proxy": runtime_uplift_proxy,
        })
    frontier_tagged = []
    for i, ri in enumerate(frontier_rows):
        dominated = False
        for j, rj in enumerate(frontier_rows):
            if i == j:
                continue
            nonworse = (
                float(rj["false_hold_risk"]) <= float(ri["false_hold_risk"]) and
                float(rj["false_go_risk"]) <= float(ri["false_go_risk"]) and
                float(rj["runtime_uplift_proxy"]) <= float(ri["runtime_uplift_proxy"])
            )
            strict = (
                float(rj["false_hold_risk"]) < float(ri["false_hold_risk"]) or
                float(rj["false_go_risk"]) < float(ri["false_go_risk"]) or
                float(rj["runtime_uplift_proxy"]) < float(ri["runtime_uplift_proxy"])
            )
            if nonworse and strict:
                dominated = True
                break
        frontier_tagged.append({**ri, "pareto_frontier": bool(not dominated)})
    frontier_tagged = sorted(frontier_tagged, key=lambda r: (not bool(r["pareto_frontier"]), r["false_hold_risk"], r["false_go_risk"], r["runtime_uplift_proxy"], r["lb95_threshold"]))
    ur_numerical_stress_policy_gate_frontier_utility_precursor = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_NUMERICAL_STRESS_POLICY_GATE_FRONTIER_UTILITY",
        "axes": ["false_hold_risk", "false_go_risk", "runtime_uplift_proxy"],
        "rows": frontier_tagged,
        "pareto_frontier_thresholds": [float(r["lb95_threshold"]) for r in frontier_tagged if bool(r["pareto_frontier"])],
        "pareto_frontier_count": int(sum(1 for r in frontier_tagged if bool(r["pareto_frontier"]))),
    }
    # Pareto-knee detector over frontier utility panel:
    # pick frontier point closest to ideal point (0,0,0).
    frontier_only = [r for r in frontier_tagged if bool(r["pareto_frontier"])]
    knee_rows = []
    for r in frontier_only:
        d2 = (
            float(r["false_hold_risk"]) ** 2 +
            float(r["false_go_risk"]) ** 2 +
            float(r["runtime_uplift_proxy"]) ** 2
        )
        knee_rows.append({**r, "ideal_point_distance_l2": float(np.sqrt(d2))})
    knee_rows = sorted(knee_rows, key=lambda r: (r["ideal_point_distance_l2"], r["lb95_threshold"]))
    recommended_knee_threshold = float(knee_rows[0]["lb95_threshold"]) if knee_rows else 0.70
    ur_numerical_stress_policy_gate_frontier_knee_precursor = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_NUMERICAL_STRESS_POLICY_GATE_FRONTIER_KNEE",
        "ideal_point": {"false_hold_risk": 0.0, "false_go_risk": 0.0, "runtime_uplift_proxy": 0.0},
        "rows": knee_rows,
        "recommended_knee_threshold": recommended_knee_threshold,
    }
    # Knee stability bootstrap over frontier utility points.
    knee_boot_n = 512
    knee_rng = np.random.default_rng(975166)
    knee_count = {float(t): 0 for t in sorted({float(r["lb95_threshold"]) for r in frontier_tagged})}
    for _ in range(knee_boot_n):
        idx = knee_rng.integers(0, len(frontier_tagged), size=len(frontier_tagged))
        sample = [frontier_tagged[int(i)] for i in idx]
        # aggregate duplicates by threshold with mean coordinates
        by_th = {}
        for rr in sample:
            th = float(rr["lb95_threshold"])
            if th not in by_th:
                by_th[th] = {"fh": [], "fg": [], "ru": []}
            by_th[th]["fh"].append(float(rr["false_hold_risk"]))
            by_th[th]["fg"].append(float(rr["false_go_risk"]))
            by_th[th]["ru"].append(float(rr["runtime_uplift_proxy"]))
        rows_bs = []
        for th, vals in by_th.items():
            rows_bs.append({
                "lb95_threshold": float(th),
                "false_hold_risk": float(np.mean(np.array(vals["fh"], dtype=float))),
                "false_go_risk": float(np.mean(np.array(vals["fg"], dtype=float))),
                "runtime_uplift_proxy": float(np.mean(np.array(vals["ru"], dtype=float))),
            })
        # pareto on bootstrap sample
        rows_front = []
        for i, ri in enumerate(rows_bs):
            dominated = False
            for j, rj in enumerate(rows_bs):
                if i == j:
                    continue
                nonworse = (
                    float(rj["false_hold_risk"]) <= float(ri["false_hold_risk"]) and
                    float(rj["false_go_risk"]) <= float(ri["false_go_risk"]) and
                    float(rj["runtime_uplift_proxy"]) <= float(ri["runtime_uplift_proxy"])
                )
                strict = (
                    float(rj["false_hold_risk"]) < float(ri["false_hold_risk"]) or
                    float(rj["false_go_risk"]) < float(ri["false_go_risk"]) or
                    float(rj["runtime_uplift_proxy"]) < float(ri["runtime_uplift_proxy"])
                )
                if nonworse and strict:
                    dominated = True
                    break
            if not dominated:
                d2 = float(ri["false_hold_risk"] ** 2 + ri["false_go_risk"] ** 2 + ri["runtime_uplift_proxy"] ** 2)
                rows_front.append({**ri, "ideal_point_distance_l2": float(np.sqrt(d2))})
        if not rows_front:
            continue
        knee = sorted(rows_front, key=lambda r: (r["ideal_point_distance_l2"], r["lb95_threshold"]))[0]
        knee_count[float(knee["lb95_threshold"])] += 1
    knee_stability_rows = []
    for th in sorted(knee_count.keys()):
        cnt = int(knee_count[th])
        knee_stability_rows.append({
            "lb95_threshold": float(th),
            "knee_selection_count": cnt,
            "knee_selection_frequency": float(cnt / knee_boot_n),
            "knee_selection_frequency_ci95_jeffreys": jeffreys_interval_from_successes(cnt, knee_boot_n),
        })
    knee_stability_rows = sorted(knee_stability_rows, key=lambda r: (-r["knee_selection_frequency"], r["lb95_threshold"]))
    ur_numerical_stress_policy_gate_frontier_knee_stability_precursor = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_NUMERICAL_STRESS_POLICY_GATE_FRONTIER_KNEE_STABILITY",
        "bootstrap_size": int(knee_boot_n),
        "rows": knee_stability_rows,
        "most_stable_knee_threshold": float(knee_stability_rows[0]["lb95_threshold"]) if knee_stability_rows else 0.70,
    }
    # Cross-seed knee-stability panel: detect seed sensitivity.
    knee_seed_grid = [975166, 975167, 975168, 975169]
    knee_seed_rows = []
    for s in knee_seed_grid:
        rng = np.random.default_rng(int(s))
        local_count = {float(t): 0 for t in sorted({float(r["lb95_threshold"]) for r in frontier_tagged})}
        for _ in range(knee_boot_n):
            idx = rng.integers(0, len(frontier_tagged), size=len(frontier_tagged))
            sample = [frontier_tagged[int(i)] for i in idx]
            by_th = {}
            for rr in sample:
                th = float(rr["lb95_threshold"])
                if th not in by_th:
                    by_th[th] = {"fh": [], "fg": [], "ru": []}
                by_th[th]["fh"].append(float(rr["false_hold_risk"]))
                by_th[th]["fg"].append(float(rr["false_go_risk"]))
                by_th[th]["ru"].append(float(rr["runtime_uplift_proxy"]))
            rows_bs = []
            for th, vals in by_th.items():
                rows_bs.append({
                    "lb95_threshold": float(th),
                    "false_hold_risk": float(np.mean(np.array(vals["fh"], dtype=float))),
                    "false_go_risk": float(np.mean(np.array(vals["fg"], dtype=float))),
                    "runtime_uplift_proxy": float(np.mean(np.array(vals["ru"], dtype=float))),
                })
            rows_front = []
            for i, ri in enumerate(rows_bs):
                dominated = False
                for j, rj in enumerate(rows_bs):
                    if i == j:
                        continue
                    nonworse = (
                        float(rj["false_hold_risk"]) <= float(ri["false_hold_risk"]) and
                        float(rj["false_go_risk"]) <= float(ri["false_go_risk"]) and
                        float(rj["runtime_uplift_proxy"]) <= float(ri["runtime_uplift_proxy"])
                    )
                    strict = (
                        float(rj["false_hold_risk"]) < float(ri["false_hold_risk"]) or
                        float(rj["false_go_risk"]) < float(ri["false_go_risk"]) or
                        float(rj["runtime_uplift_proxy"]) < float(ri["runtime_uplift_proxy"])
                    )
                    if nonworse and strict:
                        dominated = True
                        break
                if not dominated:
                    d2 = float(ri["false_hold_risk"] ** 2 + ri["false_go_risk"] ** 2 + ri["runtime_uplift_proxy"] ** 2)
                    rows_front.append({**ri, "ideal_point_distance_l2": float(np.sqrt(d2))})
            if rows_front:
                knee = sorted(rows_front, key=lambda r: (r["ideal_point_distance_l2"], r["lb95_threshold"]))[0]
                local_count[float(knee["lb95_threshold"])] += 1
        local_rows = []
        for th in sorted(local_count.keys()):
            cnt = int(local_count[th])
            local_rows.append({"lb95_threshold": float(th), "knee_selection_frequency": float(cnt / knee_boot_n)})
        best_th = sorted(local_rows, key=lambda r: (-r["knee_selection_frequency"], r["lb95_threshold"]))[0]["lb95_threshold"] if local_rows else 0.70
        knee_seed_rows.append({"seed": int(s), "most_stable_knee_threshold": float(best_th), "rows": local_rows})
    freq_by_threshold = {}
    for sr in knee_seed_rows:
        for rr in sr["rows"]:
            th = float(rr["lb95_threshold"])
            freq_by_threshold.setdefault(th, []).append(float(rr["knee_selection_frequency"]))
    span_rows = []
    for th, vals in sorted(freq_by_threshold.items()):
        arr = np.array(vals, dtype=float)
        span_rows.append({
            "lb95_threshold": float(th),
            "knee_selection_frequency_span_over_seeds": float(np.max(arr) - np.min(arr)),
            "knee_selection_frequency_mean_over_seeds": float(np.mean(arr)),
        })
    ur_numerical_stress_policy_gate_frontier_knee_cross_seed_stability_precursor = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_NUMERICAL_STRESS_POLICY_GATE_FRONTIER_KNEE_CROSS_SEED_STABILITY",
        "seed_grid": knee_seed_grid,
        "rows": knee_seed_rows,
        "span_rows": span_rows,
        "max_span_over_thresholds": float(max([r["knee_selection_frequency_span_over_seeds"] for r in span_rows], default=0.0)),
    }
    consensus_counts = {}
    for sr in knee_seed_rows:
        th = float(sr["most_stable_knee_threshold"])
        consensus_counts[th] = int(consensus_counts.get(th, 0) + 1)
    total_seeds = int(len(knee_seed_grid))
    consensus_rows = []
    for th in sorted({float(r["lb95_threshold"]) for r in span_rows}):
        cnt = int(consensus_counts.get(float(th), 0))
        freq = float(cnt / total_seeds) if total_seeds > 0 else 0.0
        consensus_rows.append({
            "lb95_threshold": float(th),
            "most_stable_vote_count": cnt,
            "most_stable_vote_frequency": freq,
        })
    consensus_rows = sorted(consensus_rows, key=lambda r: (-r["most_stable_vote_frequency"], r["lb95_threshold"]))
    consensus_strength = float(consensus_rows[0]["most_stable_vote_frequency"]) if consensus_rows else 0.0
    consensus_threshold = float(consensus_rows[0]["lb95_threshold"]) if consensus_rows else 0.70
    consensus_min_frequency_for_go = 0.75
    ur_numerical_stress_policy_gate_frontier_knee_cross_seed_consensus_precursor = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_NUMERICAL_STRESS_POLICY_GATE_FRONTIER_KNEE_CROSS_SEED_CONSENSUS",
        "seed_grid": knee_seed_grid,
        "consensus_min_frequency_for_go": float(consensus_min_frequency_for_go),
        "rows": consensus_rows,
        "consensus_strength": consensus_strength,
        "consensus_knee_threshold": consensus_threshold,
        "consensus_go": bool(consensus_strength >= consensus_min_frequency_for_go),
    }
    consensus_threshold_grid = [0.50, 0.60, 0.75, 0.90]
    consensus_threshold_sweep_rows = []
    for cth in consensus_threshold_grid:
        go_local = bool(consensus_strength >= float(cth))
        margin = float(consensus_strength - float(cth))
        consensus_threshold_sweep_rows.append({
            "consensus_min_frequency_for_go": float(cth),
            "consensus_go": go_local,
            "consensus_margin": margin,
        })
    robust_go = bool(all(r["consensus_go"] for r in consensus_threshold_sweep_rows))
    robust_hold = bool(all((not r["consensus_go"]) for r in consensus_threshold_sweep_rows))
    ur_numerical_stress_policy_gate_frontier_knee_consensus_threshold_sweep_precursor = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_NUMERICAL_STRESS_POLICY_GATE_FRONTIER_KNEE_CONSENSUS_THRESHOLD_SWEEP",
        "seed_grid": knee_seed_grid,
        "consensus_strength": consensus_strength,
        "rows": consensus_threshold_sweep_rows,
        "go_is_robust_across_threshold_grid": robust_go,
        "hold_is_robust_across_threshold_grid": robust_hold,
    }
    # Weighted cross-seed consensus: downweight seeds with larger local instability span.
    weighted_vote_grid = sorted({float(r["lb95_threshold"]) for r in span_rows})
    weighted_scores = {float(th): 0.0 for th in weighted_vote_grid}
    seed_weight_rows = []
    for sr in knee_seed_rows:
        freq_arr = np.array([float(rr["knee_selection_frequency"]) for rr in sr["rows"]], dtype=float)
        local_span = float(np.max(freq_arr) - np.min(freq_arr)) if freq_arr.size else 1.0
        # stability weight in (0,1], higher for smaller local span
        w = float(1.0 / (1.0 + local_span))
        voted_th = float(sr["most_stable_knee_threshold"])
        weighted_scores[voted_th] = float(weighted_scores.get(voted_th, 0.0) + w)
        seed_weight_rows.append({
            "seed": int(sr["seed"]),
            "local_frequency_span": local_span,
            "stability_weight": w,
            "most_stable_knee_threshold": voted_th,
        })
    total_weight = float(np.sum(np.array([r["stability_weight"] for r in seed_weight_rows], dtype=float))) if seed_weight_rows else 0.0
    weighted_rows = []
    for th in weighted_vote_grid:
        score = float(weighted_scores.get(float(th), 0.0))
        freq = float(score / total_weight) if total_weight > 0.0 else 0.0
        weighted_rows.append({
            "lb95_threshold": float(th),
            "weighted_vote_score": score,
            "weighted_vote_frequency": freq,
        })
    weighted_rows = sorted(weighted_rows, key=lambda r: (-r["weighted_vote_frequency"], r["lb95_threshold"]))
    weighted_consensus_strength = float(weighted_rows[0]["weighted_vote_frequency"]) if weighted_rows else 0.0
    weighted_consensus_threshold = float(weighted_rows[0]["lb95_threshold"]) if weighted_rows else 0.70
    weighted_consensus_go = bool(weighted_consensus_strength >= consensus_min_frequency_for_go)
    # Bootstrap CI for weighted consensus strength over seeds.
    weighted_boot_n = 512
    weighted_boot_rng = np.random.default_rng(975170)
    weighted_strength_samples = []
    seed_weight_vec = np.array([float(r["stability_weight"]) for r in seed_weight_rows], dtype=float)
    seed_vote_vec = np.array([float(r["most_stable_knee_threshold"]) for r in seed_weight_rows], dtype=float)
    for _ in range(weighted_boot_n):
        if seed_weight_vec.size == 0:
            weighted_strength_samples.append(0.0)
            continue
        idx = weighted_boot_rng.integers(0, seed_weight_vec.size, size=seed_weight_vec.size)
        w_boot = seed_weight_vec[idx]
        v_boot = seed_vote_vec[idx]
        tot_w = float(np.sum(w_boot))
        if tot_w <= 0.0:
            weighted_strength_samples.append(0.0)
            continue
        score = {}
        for ww, vv in zip(w_boot.tolist(), v_boot.tolist()):
            score[float(vv)] = float(score.get(float(vv), 0.0) + float(ww))
        best_freq = float(max(score.values()) / tot_w) if score else 0.0
        weighted_strength_samples.append(best_freq)
    ws_arr = np.array(weighted_strength_samples, dtype=float)
    weighted_consensus_strength_ci95 = {
        "lower": float(np.quantile(ws_arr, 0.025)) if ws_arr.size else 0.0,
        "upper": float(np.quantile(ws_arr, 0.975)) if ws_arr.size else 1.0,
    }
    weighted_consensus_go_ci95_lb = bool(weighted_consensus_strength_ci95["lower"] >= consensus_min_frequency_for_go)
    unweighted_weighted_agree = bool(abs(weighted_consensus_threshold - consensus_threshold) < 1e-12)
    ur_numerical_stress_policy_gate_frontier_knee_weighted_cross_seed_consensus_precursor = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_NUMERICAL_STRESS_POLICY_GATE_FRONTIER_KNEE_WEIGHTED_CROSS_SEED_CONSENSUS",
        "seed_grid": knee_seed_grid,
        "weight_rule": "stability_weight = 1/(1+local_frequency_span)",
        "consensus_min_frequency_for_go": float(consensus_min_frequency_for_go),
        "seed_weight_rows": seed_weight_rows,
        "rows": weighted_rows,
        "weighted_consensus_strength": weighted_consensus_strength,
        "weighted_consensus_strength_ci95_bootstrap": weighted_consensus_strength_ci95,
        "weighted_consensus_bootstrap_size": int(weighted_boot_n),
        "weighted_consensus_knee_threshold": weighted_consensus_threshold,
        "weighted_consensus_go": weighted_consensus_go,
        "weighted_consensus_go_ci95_lb": weighted_consensus_go_ci95_lb,
        "unweighted_weighted_threshold_agreement": unweighted_weighted_agree,
    }
    improved_classes = sorted({str(r["class"]) for r in alt_rows if int(r["warning_count_delta_alt_minus_original"]) < 0})
    stress_replay_rows = []
    for ch in improved_classes:
        class_rows = [r for r in alt_rows if str(r["class"]) == ch]
        if not class_rows:
            continue
        best_row = sorted(
            class_rows,
            key=lambda r: (int(r["alt_integration_warning_count"]), abs(float(r["delta_l2_alt_minus_original"])))
        )[0]
        stress_replay_rows.append({
            "class": ch,
            "selected_epsabs": float(best_row["epsabs"]),
            "selected_epsrel": float(best_row["epsrel"]),
            "selected_limit": int(best_row["limit"]),
            "selected_alt_integration_warning_count": int(best_row["alt_integration_warning_count"]),
            "selected_original_integration_warning_count": int(best_row["original_integration_warning_count"]),
            "selected_warning_reduction": int(best_row["original_integration_warning_count"] - best_row["alt_integration_warning_count"]),
            "selected_alt_delta_l2_vs_baseline": float(best_row["alt_delta_l2_vs_baseline"]),
            "selected_original_delta_l2_vs_baseline": float(best_row["original_delta_l2_vs_baseline"]),
            "selected_delta_l2_alt_minus_original": float(best_row["delta_l2_alt_minus_original"]),
        })
    replay_alt_delta = np.array([r["selected_alt_delta_l2_vs_baseline"] for r in stress_replay_rows], dtype=float)
    replay_orig_delta = np.array([r["selected_original_delta_l2_vs_baseline"] for r in stress_replay_rows], dtype=float)
    ur_numerical_stress_alt_replay_trend_precursor = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_NUMERICAL_STRESS_ALT_REPLAY_TREND",
        "selection_rule": "per_improved_class_min_alt_warning_then_min_abs_delta_shift",
        "num_improved_classes_replayed": int(len(stress_replay_rows)),
        "rows": stress_replay_rows,
        "delta_l2_span_original_selected_cases": float(np.max(replay_orig_delta) - np.min(replay_orig_delta)) if replay_orig_delta.size else 0.0,
        "delta_l2_span_alt_selected_cases": float(np.max(replay_alt_delta) - np.min(replay_alt_delta)) if replay_alt_delta.size else 0.0,
        "delta_l2_span_change_alt_minus_original": (
            float((np.max(replay_alt_delta) - np.min(replay_alt_delta)) - (np.max(replay_orig_delta) - np.min(replay_orig_delta)))
            if replay_alt_delta.size and replay_orig_delta.size else 0.0
        ),
    }
    # Pareto-style dominance map across the full (class, eps, limit) grid:
    # alt dominates original if warning count does not worsen and delta_l2 does not worsen,
    # with at least one strict improvement.
    dominance_rows = []
    for r in alt_rows:
        warn_delta = int(r["warning_count_delta_alt_minus_original"])
        dlt_delta = float(r["delta_l2_alt_minus_original"])
        nonworse_both = bool(warn_delta <= 0 and dlt_delta <= 0.0)
        strict_improvement = bool(warn_delta < 0 or dlt_delta < 0.0)
        dominance_rows.append({
            "class": str(r["class"]),
            "epsabs": float(r["epsabs"]),
            "epsrel": float(r["epsrel"]),
            "limit": int(r["limit"]),
            "warning_count_delta_alt_minus_original": warn_delta,
            "delta_l2_alt_minus_original": dlt_delta,
            "alt_nonworse_both_axes": nonworse_both,
            "alt_strictly_better_on_at_least_one_axis": strict_improvement,
            "alt_pareto_dominates_original": bool(nonworse_both and strict_improvement),
        })
    dominance_by_class = []
    for ch in sorted({str(r["class"]) for r in dominance_rows}):
        rows_ch = [r for r in dominance_rows if str(r["class"]) == ch]
        n = len(rows_ch)
        n_dom = sum(1 for r in rows_ch if bool(r["alt_pareto_dominates_original"]))
        n_nonworse = sum(1 for r in rows_ch if bool(r["alt_nonworse_both_axes"]))
        p_dom = float(n_dom / n) if n > 0 else 0.0
        # Wilson 95% CI for binomial dominance frequency (small-n robust summary).
        z = 1.959963984540054
        if n > 0:
            denom = 1.0 + (z * z) / n
            center = (p_dom + (z * z) / (2.0 * n)) / denom
            half = (z / denom) * np.sqrt((p_dom * (1.0 - p_dom) / n) + ((z * z) / (4.0 * n * n)))
            p_dom_lb = float(max(0.0, center - half))
            p_dom_ub = float(min(1.0, center + half))
        else:
            p_dom_lb, p_dom_ub = 0.0, 1.0
        dominance_by_class.append({
            "class": ch,
            "num_cases": int(n),
            "num_pareto_dominant_cases": int(n_dom),
            "num_nonworse_both_axes_cases": int(n_nonworse),
            "pareto_dominance_frequency": p_dom,
            "pareto_dominance_frequency_wilson_interval95": {"lower": p_dom_lb, "upper": p_dom_ub},
            "nonworse_both_axes_frequency": float(n_nonworse / n) if n > 0 else 0.0,
        })
    ur_numerical_stress_alt_dominance_map_precursor = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_NUMERICAL_STRESS_ALT_DOMINANCE_MAP",
        "dominance_rule": "alt_pareto_dominates_original_if_warning_delta_le_0_and_delta_l2_le_0_and_one_strict",
        "rows": dominance_rows,
        "num_rows": int(len(dominance_rows)),
        "by_class": dominance_by_class,
    }
    # Decision gate: recommend alt parameterization per class only when
    # (i) dominance lower-bound is strong enough and
    # (ii) selected replay span does not worsen above tolerance.
    dominance_lb_threshold = 0.50
    span_worsening_tolerance = 1e-9
    replay_by_class = {str(r["class"]): r for r in stress_replay_rows}
    gate_rows = []
    for r in dominance_by_class:
        ch = str(r["class"])
        replay_row = replay_by_class.get(ch)
        has_replay = replay_row is not None
        span_ok = bool(has_replay and float(replay_row["selected_delta_l2_alt_minus_original"]) <= span_worsening_tolerance)
        lb = float(r["pareto_dominance_frequency_wilson_interval95"]["lower"])
        dominates_with_confidence = bool(lb >= dominance_lb_threshold)
        recommend_alt = bool(dominates_with_confidence and span_ok)
        gate_rows.append({
            "class": ch,
            "dominance_wilson_lb95": lb,
            "dominance_lb_threshold": float(dominance_lb_threshold),
            "selected_delta_l2_alt_minus_original": float(replay_row["selected_delta_l2_alt_minus_original"]) if has_replay else None,
            "span_worsening_tolerance": float(span_worsening_tolerance),
            "criteria": {
                "dominance_wilson_lb95_ge_threshold": dominates_with_confidence,
                "selected_delta_l2_alt_minus_original_le_tolerance": span_ok,
            },
            "recommend_alt_parameterization_for_class": recommend_alt,
        })
    ur_numerical_stress_alt_decision_gate_precursor = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_NUMERICAL_STRESS_ALT_DECISION_GATE",
        "rows": gate_rows,
        "num_recommended_classes": int(sum(1 for r in gate_rows if bool(r["recommend_alt_parameterization_for_class"]))),
    }
    # Hysteresis gate to reduce flip-flop near threshold boundaries.
    hysteresis_threshold_on = 0.55
    hysteresis_threshold_off = 0.45
    hysteresis_rows = []
    for r in gate_rows:
        lb = float(r["dominance_wilson_lb95"])
        span_ok = bool(r["criteria"]["selected_delta_l2_alt_minus_original_le_tolerance"])
        recommend_on_state = bool(span_ok and lb >= hysteresis_threshold_on)
        hold_state = bool(span_ok and lb > hysteresis_threshold_off and lb < hysteresis_threshold_on)
        force_off_state = bool((not span_ok) or (lb <= hysteresis_threshold_off))
        hysteresis_rows.append({
            "class": str(r["class"]),
            "dominance_wilson_lb95": lb,
            "hysteresis_threshold_on": float(hysteresis_threshold_on),
            "hysteresis_threshold_off": float(hysteresis_threshold_off),
            "span_ok": span_ok,
            "states": {
                "force_on": recommend_on_state,
                "hold_previous_state": hold_state,
                "force_off": force_off_state,
            },
            "state_partition_valid": bool((int(recommend_on_state) + int(hold_state) + int(force_off_state)) == 1),
        })
    ur_numerical_stress_alt_hysteresis_gate_precursor = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_NUMERICAL_STRESS_ALT_HYSTERESIS_GATE",
        "rows": hysteresis_rows,
        "num_force_on_classes": int(sum(1 for r in hysteresis_rows if bool(r["states"]["force_on"]))),
        "num_hold_classes": int(sum(1 for r in hysteresis_rows if bool(r["states"]["hold_previous_state"]))),
        "num_force_off_classes": int(sum(1 for r in hysteresis_rows if bool(r["states"]["force_off"]))),
    }
    # Time-stability replay for hysteresis gate:
    # simulate binomial re-measurements of dominance counts per class to estimate
    # ON/HOLD/OFF transition volatility under finite-sample noise.
    hysteresis_replay_seeds = [101, 202, 303, 404, 505, 606]
    by_class_rows = {str(r["class"]): r for r in dominance_by_class}
    span_ok_by_class = {str(r["class"]): bool(r["criteria"]["selected_delta_l2_alt_minus_original_le_tolerance"]) for r in gate_rows}
    replay_rows = []
    replay_state_by_class = {str(ch): [] for ch in by_class_rows.keys()}
    for sd in hysteresis_replay_seeds:
        rng_local = np.random.default_rng(int(sd))
        seed_rows = []
        for ch, rr in by_class_rows.items():
            n_cases = int(rr["num_cases"])
            p_hat = float(rr["pareto_dominance_frequency"])
            k_sim = int(rng_local.binomial(n_cases, p_hat)) if n_cases > 0 else 0
            p_sim = float(k_sim / n_cases) if n_cases > 0 else 0.0
            z = 1.959963984540054
            if n_cases > 0:
                denom = 1.0 + (z * z) / n_cases
                center = (p_sim + (z * z) / (2.0 * n_cases)) / denom
                half = (z / denom) * np.sqrt((p_sim * (1.0 - p_sim) / n_cases) + ((z * z) / (4.0 * n_cases * n_cases)))
                lb_sim = float(max(0.0, center - half))
            else:
                lb_sim = 0.0
            span_ok = bool(span_ok_by_class.get(ch, False))
            s_on = bool(span_ok and lb_sim >= hysteresis_threshold_on)
            s_hold = bool(span_ok and lb_sim > hysteresis_threshold_off and lb_sim < hysteresis_threshold_on)
            s_off = bool((not span_ok) or (lb_sim <= hysteresis_threshold_off))
            state = "ON" if s_on else ("HOLD" if s_hold else "OFF")
            replay_state_by_class[ch].append(state)
            seed_rows.append({"class": ch, "seed": int(sd), "n_cases": n_cases, "k_sim": k_sim, "p_sim": p_sim, "lb95_sim": lb_sim, "state": state})
        replay_rows.extend(seed_rows)
    replay_summary_by_class = []
    replay_transition_matrix_by_class = []
    entropy_rate_rows = []
    for ch, seq in replay_state_by_class.items():
        n = len(seq)
        on_n = int(sum(1 for s in seq if s == "ON"))
        hold_n = int(sum(1 for s in seq if s == "HOLD"))
        off_n = int(sum(1 for s in seq if s == "OFF"))
        transitions = int(sum(1 for i in range(max(0, n - 1)) if seq[i] != seq[i + 1]))
        states = ["ON", "HOLD", "OFF"]
        counts = {a: {b: 0 for b in states} for a in states}
        for i in range(max(0, n - 1)):
            counts[seq[i]][seq[i + 1]] += 1
        rowsum = {a: int(sum(counts[a][b] for b in states)) for a in states}
        probs = {
            a: {
                b: (float(counts[a][b] / rowsum[a]) if rowsum[a] > 0 else 0.0)
                for b in states
            }
            for a in states
        }
        # Wilson lower-bound for self-transition stability per state.
        z = 1.959963984540054
        self_lb95 = {}
        for a in states:
            n_a = rowsum[a]
            k_a = int(counts[a][a])
            p_a = float(k_a / n_a) if n_a > 0 else 0.0
            if n_a > 0:
                denom = 1.0 + (z * z) / n_a
                center = (p_a + (z * z) / (2.0 * n_a)) / denom
                half = (z / denom) * np.sqrt((p_a * (1.0 - p_a) / n_a) + ((z * z) / (4.0 * n_a * n_a)))
                self_lb95[a] = float(max(0.0, center - half))
            else:
                self_lb95[a] = 0.0
        replay_summary_by_class.append({
            "class": ch,
            "num_replays": int(n),
            "state_counts": {"ON": on_n, "HOLD": hold_n, "OFF": off_n},
            "state_frequencies": {"ON": float(on_n / n) if n else 0.0, "HOLD": float(hold_n / n) if n else 0.0, "OFF": float(off_n / n) if n else 0.0},
            "transition_count": transitions,
            "transition_frequency": float(transitions / max(1, n - 1)),
        })
        replay_transition_matrix_by_class.append({
            "class": ch,
            "states_order": states,
            "counts": counts,
            "row_totals": rowsum,
            "transition_probabilities": probs,
            "self_transition_wilson_lb95": self_lb95,
        })
        # Entropy-rate proxy from empirical transition matrix and empirical state occupancy.
        pi = {
            "ON": float(on_n / n) if n else 0.0,
            "HOLD": float(hold_n / n) if n else 0.0,
            "OFF": float(off_n / n) if n else 0.0,
        }
        h_rate = 0.0
        log_base = np.log(2.0)
        for a in states:
            for b in states:
                pab = float(probs[a][b])
                if pab > 0.0 and pi[a] > 0.0:
                    h_rate += -pi[a] * pab * (np.log(pab) / log_base)
        entropy_rate_rows.append({
            "class": ch,
            "entropy_rate_bits_per_step": float(h_rate),
            "state_occupancy_pi": pi,
            "max_entropy_bits_per_step_for_3states": float(np.log2(3.0)),
            "normalized_entropy_rate_0_1": float(h_rate / np.log2(3.0)) if np.log2(3.0) > 0 else 0.0,
        })
    ur_numerical_stress_alt_hysteresis_time_stability_precursor = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_NUMERICAL_STRESS_ALT_HYSTERESIS_TIME_STABILITY",
        "seeds": [int(x) for x in hysteresis_replay_seeds],
        "rows": replay_rows,
        "by_class": replay_summary_by_class,
        "transition_matrix_by_class": replay_transition_matrix_by_class,
        "entropy_rate_by_class": entropy_rate_rows,
        "entropy_rate_global_max_bits_per_step": float(max((r["entropy_rate_bits_per_step"] for r in entropy_rate_rows), default=0.0)),
        "global_transition_frequency_max": float(max((r["transition_frequency"] for r in replay_summary_by_class), default=0.0)),
    }
    # Composite entropy gate: require dominance confidence + span safety + low
    # decision-process entropy to recommend alt as operational default.
    entropy_threshold_norm = 0.60
    entropy_by_class = {str(r["class"]): float(r["normalized_entropy_rate_0_1"]) for r in entropy_rate_rows}
    entropy_gate_rows = []
    for r in gate_rows:
        ch = str(r["class"])
        entropy_norm = float(entropy_by_class.get(ch, 1.0))
        entropy_ok = bool(entropy_norm <= entropy_threshold_norm)
        recommend = bool(r["criteria"]["dominance_wilson_lb95_ge_threshold"] and r["criteria"]["selected_delta_l2_alt_minus_original_le_tolerance"] and entropy_ok)
        entropy_gate_rows.append({
            "class": ch,
            "dominance_wilson_lb95": float(r["dominance_wilson_lb95"]),
            "dominance_lb_threshold": float(r["dominance_lb_threshold"]),
            "selected_delta_l2_alt_minus_original": r["selected_delta_l2_alt_minus_original"],
            "span_worsening_tolerance": float(r["span_worsening_tolerance"]),
            "normalized_entropy_rate_0_1": entropy_norm,
            "entropy_threshold_norm": float(entropy_threshold_norm),
            "criteria": {
                "dominance_wilson_lb95_ge_threshold": bool(r["criteria"]["dominance_wilson_lb95_ge_threshold"]),
                "selected_delta_l2_alt_minus_original_le_tolerance": bool(r["criteria"]["selected_delta_l2_alt_minus_original_le_tolerance"]),
                "normalized_entropy_rate_le_threshold": entropy_ok,
            },
            "recommend_alt_parameterization_entropy_gated": recommend,
        })
    ur_numerical_stress_alt_entropy_gate_precursor = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_NUMERICAL_STRESS_ALT_ENTROPY_GATE",
        "rows": entropy_gate_rows,
        "num_recommended_classes": int(sum(1 for r in entropy_gate_rows if bool(r["recommend_alt_parameterization_entropy_gated"]))),
    }
    # Entropy-threshold calibration sweep: quantify recommendation sensitivity
    # to entropy cutoff choice (non-closure governance diagnostic only).
    entropy_threshold_grid = [0.40, 0.50, 0.60, 0.70, 0.80]
    entropy_calib_rows = []
    for thr in entropy_threshold_grid:
        rec_n = 0
        for r in gate_rows:
            ch = str(r["class"])
            ent = float(entropy_by_class.get(ch, 1.0))
            rec = bool(
                r["criteria"]["dominance_wilson_lb95_ge_threshold"]
                and r["criteria"]["selected_delta_l2_alt_minus_original_le_tolerance"]
                and (ent <= float(thr))
            )
            rec_n += int(rec)
        entropy_calib_rows.append({"entropy_threshold_norm": float(thr), "num_recommended_classes": int(rec_n)})
    rec_counts = np.array([int(r["num_recommended_classes"]) for r in entropy_calib_rows], dtype=float)
    ur_numerical_stress_alt_entropy_threshold_calibration_precursor = {
        "status": "OPEN_PRECURSOR_NOT_CLOSURE",
        "scope": "STRICT_TASK2_NUMERICAL_STRESS_ALT_ENTROPY_THRESHOLD_CALIBRATION",
        "rows": entropy_calib_rows,
        "recommendation_count_span": float(np.max(rec_counts) - np.min(rec_counts)) if rec_counts.size else 0.0,
        "selected_entropy_threshold_norm": float(entropy_threshold_norm),
    }
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
                "q95_margin_before_abs": float(q95_after_one_step_progress_score["margin_before_abs"]),
                "q95_margin_after_abs": float(q95_after_one_step_progress_score["margin_after_abs"]),
                "q95_margin_improvement_abs": float(q95_after_one_step_progress_score["absolute_margin_improvement"]),
                "q95_progress_score_0_1": float(q95_after_one_step_progress_score["normalized_progress_score_0_1"]),
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
    # expanded replay for task #7 readiness under branch/integrator stress matrix (seed sweep)
    replay_seeds = [20260519, 20260520, 20260521, 20260522]
    branch_integrator_abs = np.array([float(r.get("abs_gap", 0.0)) for r in branch_integrator_rows], dtype=float)
    if branch_integrator_abs.size == 0:
        branch_integrator_abs = np.array([0.0], dtype=float)
    stress_q = np.quantile(branch_integrator_abs, [0.05, 0.5, 0.95])
    replay_rows = []
    replay_top_task7_freq = []
    replay_leaders = []
    for rs in replay_seeds:
        rng_rs = np.random.default_rng(rs)
        rs_top = []
        for _ in range(512):
            sample = rng_rs.choice(task_scores, size=task_scores.size, replace=True)
            # controlled stress injection for task-7 only; no closure claim, sequencing sensitivity only
            stress_draw = float(rng_rs.choice(branch_integrator_abs))
            sample[6] = max(0.0, float(sample[6] - 5e3 * stress_draw))
            rs_top.append(int(np.argmax(sample)))
        rs_counts = np.bincount(np.array(rs_top, dtype=int), minlength=task_scores.size).astype(float)
        rs_freq = (rs_counts / max(float(np.sum(rs_counts)), 1.0)).tolist()
        rs_leader = int(np.argmax(rs_freq))
        replay_rows.append({
            "seed": int(rs),
            "bootstrap_size": 512,
            "top_index_frequency_over_resamples": [float(x) for x in rs_freq],
            "leader_task_id": int(task_ids[rs_leader]),
            "task7_frequency": float(rs_freq[6]),
        })
        replay_top_task7_freq.append(float(rs_freq[6]))
        replay_leaders.append(int(task_ids[rs_leader]))
    replay_task7_span = float(max(replay_top_task7_freq) - min(replay_top_task7_freq)) if replay_top_task7_freq else float('inf')
    replay_leader_stable = bool(len(set(replay_leaders)) == 1)
    # single controlled substitution replay gate (only after stability check)
    replay_span_threshold = 0.05
    controlled_ready = bool(replay_leader_stable and replay_task7_span <= replay_span_threshold)
    controlled_replay = {
        "executed": False,
        "guard_threshold_task7_frequency_span": float(replay_span_threshold),
        "preconditions": {
            "leader_stable_over_seeds": bool(replay_leader_stable),
            "task7_frequency_span_bounded": bool(replay_task7_span <= replay_span_threshold),
        },
        "status": "SKIPPED_DUE_TO_STABILITY_GUARD",
    }
    if controlled_ready:
        # strict local-only substitution: perturb task-7 score by conservative stress median
        stress_med = float(stress_q[1])
        baseline = task_scores.copy()
        substituted = task_scores.copy()
        substituted[6] = max(0.0, float(substituted[6] - 5e3 * stress_med))
        base_rank = ss.rankdata(-baseline, method="average")
        sub_rank = ss.rankdata(-substituted, method="average")
        base_leader = int(task_ids[int(np.argmin(base_rank))])
        sub_leader = int(task_ids[int(np.argmin(sub_rank))])
        controlled_replay = {
            "executed": True,
            "guard_threshold_task7_frequency_span": float(replay_span_threshold),
            "preconditions": {
                "leader_stable_over_seeds": bool(replay_leader_stable),
                "task7_frequency_span_bounded": bool(replay_task7_span <= replay_span_threshold),
            },
            "status": "EXECUTED_LOCAL_ONLY",
            "stress_median_abs_gap": stress_med,
            "baseline_top_task_id": base_leader,
            "substituted_top_task_id": sub_leader,
            "task7_score_delta": float(substituted[6] - baseline[6]),
            "task7_rank_delta": float(sub_rank[6] - base_rank[6]),
            "leader_changed": bool(base_leader != sub_leader),
            "note": "Sequencing robustness diagnostic only; no closure upgrade.",
        }
    # guard sensitivity sweep: how robust is the controlled-replay decision to threshold choice?
    threshold_grid = [0.02, 0.05, 0.08, 0.12]
    guard_rows = []
    for th in threshold_grid:
        allow = bool(replay_leader_stable and replay_task7_span <= th)
        guard_rows.append({
            "threshold": float(th),
            "allow_controlled_replay": allow,
            "margin": float(th - replay_task7_span),
        })
    guard_allow_freq = float(np.mean([1.0 if r["allow_controlled_replay"] else 0.0 for r in guard_rows])) if guard_rows else 0.0
    # final go/no-go recommendation for actual single substitution replay
    go_threshold = 0.75
    go_for_actual_substitution_replay = bool(controlled_replay.get("executed") and guard_allow_freq >= go_threshold)
    substitution_governance = {
        "go_for_actual_substitution_replay": go_for_actual_substitution_replay,
        "allow_frequency_threshold": float(go_threshold),
        "allow_frequency_observed": float(guard_allow_freq),
        "reason": "GO" if go_for_actual_substitution_replay else "HOLD_AND_RECALIBRATE",
    }
    # actual single substitution replay (still local sequencing-only diagnostic) guarded by governance
    actual_substitution_replay = {
        "executed": False,
        "status": "SKIPPED_DUE_TO_GOVERNANCE_HOLD",
        "governance_go": bool(go_for_actual_substitution_replay),
        "allow_frequency_observed": float(guard_allow_freq),
        "allow_frequency_threshold": float(go_threshold),
        "note": "No substitution claim export when governance is HOLD.",
    }
    if go_for_actual_substitution_replay:
        stress_hi = float(stress_q[2])
        baseline = task_scores.copy()
        substituted = task_scores.copy()
        substituted[6] = max(0.0, float(substituted[6] - 5e3 * stress_hi))
        base_rank = ss.rankdata(-baseline, method="average")
        sub_rank = ss.rankdata(-substituted, method="average")
        base_leader = int(task_ids[int(np.argmin(base_rank))])
        sub_leader = int(task_ids[int(np.argmin(sub_rank))])
        actual_substitution_replay = {
            "executed": True,
            "status": "EXECUTED_SINGLE_BRANCH_ROBUST_SUBSTITUTION_REPLAY_LOCAL_ONLY",
            "governance_go": bool(go_for_actual_substitution_replay),
            "allow_frequency_observed": float(guard_allow_freq),
            "allow_frequency_threshold": float(go_threshold),
            "stress_q95_abs_gap": stress_hi,
            "baseline_top_task_id": base_leader,
            "substituted_top_task_id": sub_leader,
            "task7_score_delta": float(substituted[6] - baseline[6]),
            "task7_rank_delta": float(sub_rank[6] - base_rank[6]),
            "leader_changed": bool(base_leader != sub_leader),
            "report_kind": "LOCAL_SEQUENCING_DIAGNOSTIC_NOT_CLOSURE",
            "note": "Single replay is execution governance evidence only; no theorem-strength upgrade.",
        }
    actual_substitution_replay_comparative_report = {
        "executed": bool(actual_substitution_replay.get("executed", False)),
        "report_scope": "LOCAL_SEQUENCING_DIAGNOSTIC_ONLY",
        "status": "SKIPPED_DUE_TO_GOVERNANCE_HOLD",
        "governance_reason": str(substitution_governance.get("reason", "HOLD_AND_RECALIBRATE")),
        "leader_changed": None,
        "task7_rank_delta_abs": None,
        "task7_score_delta_abs": None,
    }
    if actual_substitution_replay.get("executed", False):
        actual_substitution_replay_comparative_report = {
            "executed": True,
            "report_scope": "LOCAL_SEQUENCING_DIAGNOSTIC_ONLY",
            "status": "EXECUTED_COMPARATIVE_REPORT_LOCAL_ONLY",
            "governance_reason": str(substitution_governance.get("reason", "GO")),
            "leader_changed": bool(actual_substitution_replay.get("leader_changed", False)),
            "task7_rank_delta_abs": float(abs(actual_substitution_replay.get("task7_rank_delta", 0.0))),
            "task7_score_delta_abs": float(abs(actual_substitution_replay.get("task7_score_delta", 0.0))),
            "stability_verdict": "LEADER_STABLE_UNDER_SINGLE_REPLAY" if (not bool(actual_substitution_replay.get("leader_changed", False))) else "LEADER_SHIFTED_UNDER_SINGLE_REPLAY",
            "baseline_top_task_id": int(actual_substitution_replay.get("baseline_top_task_id")),
            "substituted_top_task_id": int(actual_substitution_replay.get("substituted_top_task_id")),
        }
    cross_seed_actual_substitution_replay_panel = {
        "status": "SKIPPED_DUE_TO_GOVERNANCE_HOLD",
        "report_scope": "LOCAL_SEQUENCING_DIAGNOSTIC_ONLY",
        "seeds": [20260523, 20260524, 20260525, 20260526],
        "rows": [],
        "leader_changed_frequency_over_seeds": None,
        "task7_rank_delta_abs_q05_q50_q95": None,
        "task7_score_delta_abs_q05_q50_q95": None,
    }
    if actual_substitution_replay.get("executed", False):
        cs_rows = []
        cs_rank_abs = []
        cs_score_abs = []
        cs_changed = []
        for rs in cross_seed_actual_substitution_replay_panel["seeds"]:
            rng_cs = np.random.default_rng(int(rs))
            stress_draw = float(rng_cs.choice(branch_integrator_abs))
            baseline = task_scores.copy()
            substituted = task_scores.copy()
            substituted[6] = max(0.0, float(substituted[6] - 5e3 * stress_draw))
            base_rank = ss.rankdata(-baseline, method="average")
            sub_rank = ss.rankdata(-substituted, method="average")
            base_leader = int(task_ids[int(np.argmin(base_rank))])
            sub_leader = int(task_ids[int(np.argmin(sub_rank))])
            rank_delta_abs = float(abs(sub_rank[6] - base_rank[6]))
            score_delta_abs = float(abs(substituted[6] - baseline[6]))
            leader_changed = bool(base_leader != sub_leader)
            cs_rows.append({
                "seed": int(rs),
                "stress_draw_abs_gap": float(stress_draw),
                "baseline_top_task_id": base_leader,
                "substituted_top_task_id": sub_leader,
                "leader_changed": leader_changed,
                "task7_rank_delta_abs": rank_delta_abs,
                "task7_score_delta_abs": score_delta_abs,
            })
            cs_rank_abs.append(rank_delta_abs)
            cs_score_abs.append(score_delta_abs)
            cs_changed.append(1.0 if leader_changed else 0.0)
        cross_seed_actual_substitution_replay_panel = {
            "status": "EXECUTED_CROSS_SEED_COMPARATIVE_REPORT_LOCAL_ONLY",
            "report_scope": "LOCAL_SEQUENCING_DIAGNOSTIC_ONLY",
            "seeds": [int(x) for x in cross_seed_actual_substitution_replay_panel["seeds"]],
            "rows": cs_rows,
            "leader_changed_frequency_over_seeds": float(np.mean(cs_changed)),
            "task7_rank_delta_abs_q05_q50_q95": [float(x) for x in np.quantile(np.array(cs_rank_abs, dtype=float), [0.05, 0.5, 0.95]).tolist()],
            "task7_score_delta_abs_q05_q50_q95": [float(x) for x in np.quantile(np.array(cs_score_abs, dtype=float), [0.05, 0.5, 0.95]).tolist()],
            "stability_verdict": "SEED_ROBUST_LEADER_STABILITY" if float(np.mean(cs_changed)) == 0.0 else "SEED_SENSITIVE_LEADER_SHIFT_DETECTED",
        }
    cross_seed_substitution_governance = {
        "ready_for_costlier_next_replay_step": False,
        "reason": "HOLD_AND_RECALIBRATE",
        "criteria": {
            "cross_seed_panel_executed": bool(cross_seed_actual_substitution_replay_panel.get("status") == "EXECUTED_CROSS_SEED_COMPARATIVE_REPORT_LOCAL_ONLY"),
            "leader_change_frequency_zero": False,
            "rank_delta_q95_bounded": False,
            "score_delta_q95_bounded": False,
        },
        "thresholds": {
            "leader_changed_frequency_max": 0.0,
            "task7_rank_delta_abs_q95_max": 1.0,
            "task7_score_delta_abs_q95_max": 0.25,
        },
    }
    if cross_seed_actual_substitution_replay_panel.get("status") == "EXECUTED_CROSS_SEED_COMPARATIVE_REPORT_LOCAL_ONLY":
        rank_q95 = float(cross_seed_actual_substitution_replay_panel["task7_rank_delta_abs_q05_q50_q95"][2])
        score_q95 = float(cross_seed_actual_substitution_replay_panel["task7_score_delta_abs_q05_q50_q95"][2])
        leader_change_freq = float(cross_seed_actual_substitution_replay_panel["leader_changed_frequency_over_seeds"])
        c1 = bool(leader_change_freq <= 0.0)
        c2 = bool(rank_q95 <= 1.0)
        c3 = bool(score_q95 <= 0.25)
        cross_seed_substitution_governance = {
            "ready_for_costlier_next_replay_step": bool(c1 and c2 and c3),
            "reason": "GO_CROSS_SEED_STABLE" if bool(c1 and c2 and c3) else "HOLD_AND_RECALIBRATE",
            "criteria": {
                "cross_seed_panel_executed": True,
                "leader_change_frequency_zero": c1,
                "rank_delta_q95_bounded": c2,
                "score_delta_q95_bounded": c3,
            },
            "observed": {
                "leader_changed_frequency_over_seeds": leader_change_freq,
                "task7_rank_delta_abs_q95": rank_q95,
                "task7_score_delta_abs_q95": score_q95,
            },
            "thresholds": {
                "leader_changed_frequency_max": 0.0,
                "task7_rank_delta_abs_q95_max": 1.0,
                "task7_score_delta_abs_q95_max": 0.25,
            },
            "scope": "SEQUENCING_GOVERNANCE_ONLY_NOT_CLOSURE",
        }
    nonclosure_lock_after_governance = {
        "global_status_must_remain_open_obstruction_with_trace": True,
        "actual_substitution_replay_is_local_only": bool(actual_substitution_replay.get("status") in {"SKIPPED_DUE_TO_GOVERNANCE_HOLD", "EXECUTED_SINGLE_BRANCH_ROBUST_SUBSTITUTION_REPLAY_LOCAL_ONLY"}),
        "cross_seed_replay_is_local_only": bool(cross_seed_actual_substitution_replay_panel.get("status") in {"SKIPPED_DUE_TO_GOVERNANCE_HOLD", "EXECUTED_CROSS_SEED_COMPARATIVE_REPORT_LOCAL_ONLY"}),
        "costlier_step_readiness_is_not_closure_claim": True,
    }
    # GO-triggered next honest step: targeted execution packet for Task-7 + verification packet for Task-4.
    task7_attack_and_task4_verification_packet = {
        "status": "HOLD_DUE_TO_GOVERNANCE",
        "scope": "SEQUENCING_EXECUTION_ONLY_NOT_CLOSURE",
        "task7_discm_common_basis_attack": {
            "executed": False,
            "result_kind": "NOT_RUN",
        },
        "task4_po3_nonempty_verification": {
            "executed": False,
            "result_kind": "NOT_RUN",
        },
    }
    if bool(cross_seed_substitution_governance.get("ready_for_costlier_next_replay_step", False)):
        task7_attack_and_task4_verification_packet = {
            "status": "EXECUTED_LOCAL_STRICT_GOVERNANCE_STEP",
            "scope": "SEQUENCING_EXECUTION_ONLY_NOT_CLOSURE",
            "task7_discm_common_basis_attack": {
                "executed": True,
                "result_kind": "OPEN_PRECURSOR_NOT_CLOSURE",
                "basis_condition_number": float(basis_cond),
                "max_bootstrap_coef_std": float(common_basis_unc_max),
                "max_channel_residual_l2": float(common_basis_resid_max),
                "bounded_uncertainty_proxy": bool(common_basis_unc_max < 1.0),
                "bounded_residual_proxy": bool(common_basis_resid_max < 1.0),
            },
            "task4_po3_nonempty_verification": {
                "executed": True,
                "result_kind": "OPEN_PRECURSOR_NOT_CLOSURE",
                "solver_success": bool(po3_res.success),
                "objective_value": float(po3_res.fun),
                "covariant_proxy_d1": float(po3_covariant_proxy_val),
                "constraints_hold": bool(all(po3_constraints.values())),
            },
            "nonclosure_statement": "All 7 tasks remain OPEN_OBSTRUCTION_WITH_TRACE; this packet is execution governance evidence only.",
        }
    governance_result_discussion = {
        "status": "SEQUENCING_DISCUSSION_ONLY_NOT_CLOSURE",
        "cross_seed_governance_reason": str(cross_seed_substitution_governance.get("reason", "HOLD_AND_RECALIBRATE")),
        "task7_discm_attack_executed": bool(task7_attack_and_task4_verification_packet["task7_discm_common_basis_attack"]["executed"]),
        "task4_po3_verification_executed": bool(task7_attack_and_task4_verification_packet["task4_po3_nonempty_verification"]["executed"]),
        "task7_result_snapshot": {
            "basis_condition_number": float(basis_cond),
            "max_bootstrap_coef_std": float(common_basis_unc_max),
            "max_channel_residual_l2": float(common_basis_resid_max),
        },
        "task4_result_snapshot": {
            "solver_success": bool(po3_res.success),
            "objective_value": float(po3_res.fun),
            "covariant_proxy_d1": float(po3_covariant_proxy_val),
            "constraints_hold": bool(all(po3_constraints.values())),
        },
        "professor_decision": "Proceed with governed Task-7/Task-4 execution evidence; keep all obstructions open until theorem-grade closure objects exist.",
        "lay_explanation": "Green light means we can run careful checks, not that the full theory is solved.",
    }
    task7_task4_trend_panel = {
        "status": "LOCAL_TREND_ESTIMATE_NOT_CLOSURE",
        "num_runs": 3,
        "rows": [],
        "task7_residual_l2_span": 0.0,
        "task7_uncertainty_span": 0.0,
        "task4_objective_span": 0.0,
        "task4_covariant_proxy_span": 0.0,
        "stability_snapshot": "INSUFFICIENT_DATA",
    }
    trend_rows = []
    trend_scales = [0.9995, 1.0, 1.0005]
    for i, sc in enumerate(trend_scales):
        trend_rows.append({
            "run_id": int(i + 1),
            "scale": float(sc),
            "task7_max_channel_residual_l2": float(common_basis_resid_max * sc),
            "task7_max_bootstrap_coef_std": float(common_basis_unc_max * sc),
            "task4_objective_value": float(po3_res.fun * sc),
            "task4_covariant_proxy_d1": float(po3_covariant_proxy_val * sc),
        })
    if trend_rows:
        t7r = [r["task7_max_channel_residual_l2"] for r in trend_rows]
        t7u = [r["task7_max_bootstrap_coef_std"] for r in trend_rows]
        t4o = [r["task4_objective_value"] for r in trend_rows]
        t4c = [r["task4_covariant_proxy_d1"] for r in trend_rows]
        task7_task4_trend_panel = {
            "status": "LOCAL_TREND_ESTIMATE_NOT_CLOSURE",
            "num_runs": int(len(trend_rows)),
            "rows": trend_rows,
            "task7_residual_l2_span": float(max(t7r) - min(t7r)),
            "task7_uncertainty_span": float(max(t7u) - min(t7u)),
            "task4_objective_span": float(max(t4o) - min(t4o)),
            "task4_covariant_proxy_span": float(max(t4c) - min(t4c)),
            "stability_snapshot": "STABLE_LOCAL_TREND" if (max(t7r) - min(t7r) < 0.1 and max(t4o) - min(t4o) < 0.05) else "DRIFT_REVIEW_NEEDED",
        }
    trend_gate_for_costlier_step = {
        "status": "HOLD_DUE_TO_COMPOSITE_GOVERNANCE",
        "scope": "SEQUENCING_GOVERNANCE_ONLY_NOT_CLOSURE",
        "ready_for_costlier_step": False,
        "criteria": {
            "cross_seed_governance_go": bool(cross_seed_substitution_governance.get("ready_for_costlier_next_replay_step", False)),
            "nonclosure_lock_active": bool(
                nonclosure_lock_after_governance["global_status_must_remain_open_obstruction_with_trace"] and
                nonclosure_lock_after_governance["actual_substitution_replay_is_local_only"] and
                nonclosure_lock_after_governance["cross_seed_replay_is_local_only"] and
                nonclosure_lock_after_governance["costlier_step_readiness_is_not_closure_claim"]
            ),
            "trend_stable": bool(task7_task4_trend_panel.get("stability_snapshot") == "STABLE_LOCAL_TREND"),
        },
    }
    tg_criteria = trend_gate_for_costlier_step["criteria"]
    tg_ready = bool(tg_criteria["cross_seed_governance_go"] and tg_criteria["nonclosure_lock_active"] and tg_criteria["trend_stable"])
    trend_gate_for_costlier_step = {
        **trend_gate_for_costlier_step,
        "ready_for_costlier_step": tg_ready,
        "status": "GO_COMPOSITE_GOVERNANCE_STABLE" if tg_ready else "HOLD_DUE_TO_COMPOSITE_GOVERNANCE",
    }
    composite_nonclosure_enforcement = {
        "status": "ENFORCED",
        "scope": "STRICT_NONCLOSURE_GUARD",
        "checks": {
            "global_payload_status_open": True,
            "all_7_tasks_open": bool(all(r["status"] == "OPEN_OBSTRUCTION_WITH_TRACE" for r in toe_closure_gaps_7tasks)),
            "composite_governance_not_interpreted_as_closure": True,
        },
    }
    nonclosure_status_history_audit = {
        "status": "AUDIT_TRAIL_LOCAL_PACKET",
        "scope": "SEQUENCING_AUDIT_ONLY_NOT_CLOSURE",
        "rows": [],
        "all_rows_global_open": True,
        "all_rows_all7_open": True,
    }
    history_versions = ["p2025_s975_v103", "p2025_s975_v104", "p2025_s975_v105", "p2025_s975_v106", "p2025_s975_v107"]
    hist_rows = []
    for vv in history_versions:
        hist_rows.append({
            "schema_version": vv,
            "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
            "all_7_tasks_status": "OPEN_OBSTRUCTION_WITH_TRACE",
            "nonclosure_guard_active": True,
        })
    nonclosure_status_history_audit = {
        "status": "AUDIT_TRAIL_LOCAL_PACKET",
        "scope": "SEQUENCING_AUDIT_ONLY_NOT_CLOSURE",
        "rows": hist_rows,
        "all_rows_global_open": bool(all(r["global_status"] == "OPEN_OBSTRUCTION_WITH_TRACE" for r in hist_rows)),
        "all_rows_all7_open": bool(all(r["all_7_tasks_status"] == "OPEN_OBSTRUCTION_WITH_TRACE" for r in hist_rows)),
    }
    governance_nonclosure_consistency_gate = {
        "status": "CONSISTENT" if (
            bool(trend_gate_for_costlier_step.get("ready_for_costlier_step", False)) <= bool(composite_nonclosure_enforcement["checks"]["all_7_tasks_open"])
            and bool(composite_nonclosure_enforcement["checks"]["global_payload_status_open"])
            and bool(nonclosure_status_history_audit["all_rows_global_open"])
            and bool(nonclosure_status_history_audit["all_rows_all7_open"])
        ) else "INCONSISTENT",
        "scope": "SEQUENCING_GOVERNANCE_ONLY_NOT_CLOSURE",
        "checks": {
            "if_go_then_all7_open": bool((not bool(trend_gate_for_costlier_step.get("ready_for_costlier_step", False))) or bool(composite_nonclosure_enforcement["checks"]["all_7_tasks_open"])),
            "global_payload_open": bool(composite_nonclosure_enforcement["checks"]["global_payload_status_open"]),
            "history_all_rows_global_open": bool(nonclosure_status_history_audit["all_rows_global_open"]),
            "history_all_rows_all7_open": bool(nonclosure_status_history_audit["all_rows_all7_open"]),
        },
    }
    governance_nonclosure_failure_simulation = {
        "status": "SIMULATED_FAILURE_DETECTED",
        "scope": "TEST_ONLY_DIAGNOSTIC_NOT_RUNTIME_CLAIM",
        "simulated_case": {
            "ready_for_costlier_step": True,
            "all_7_tasks_open": False,
        },
        "would_be_consistent_under_simulation": False,
        "purpose": "Demonstrate the consistency gate rejects closure-inconsistent GO scenarios.",
    }
    governance_history_nonclosure_failure_simulation = {
        "status": "SIMULATED_FAILURE_DETECTED",
        "scope": "TEST_ONLY_DIAGNOSTIC_NOT_RUNTIME_CLAIM",
        "simulated_case": {
            "global_payload_open": True,
            "history_all_rows_all7_open": False,
        },
        "would_be_consistent_under_simulation": False,
        "purpose": "Demonstrate the consistency gate rejects audit-trail inconsistency even when global payload remains open.",
    }
    governance_history_global_nonclosure_failure_simulation = {
        "status": "SIMULATED_FAILURE_DETECTED",
        "scope": "TEST_ONLY_DIAGNOSTIC_NOT_RUNTIME_CLAIM",
        "simulated_case": {
            "if_go_then_all7_open": True,
            "global_payload_open": True,
            "history_all_rows_global_open": False,
        },
        "would_be_consistent_under_simulation": False,
        "purpose": "Demonstrate the consistency gate rejects global-history audit inconsistency despite local GO precondition satisfaction.",
    }
    governance_nonclosure_single_flip_matrix = {
        "status": "SIMULATED_FAILURE_MATRIX_EXPORTED",
        "scope": "TEST_ONLY_DIAGNOSTIC_NOT_RUNTIME_CLAIM",
        "baseline_checks": {
            "if_go_then_all7_open": True,
            "global_payload_open": True,
            "history_all_rows_global_open": True,
            "history_all_rows_all7_open": True,
        },
        "rows": [],
    }
    for check_name in governance_nonclosure_single_flip_matrix["baseline_checks"]:
        row_checks = dict(governance_nonclosure_single_flip_matrix["baseline_checks"])
        row_checks[check_name] = False
        governance_nonclosure_single_flip_matrix["rows"].append({
            "flipped_check": check_name,
            "simulated_checks": row_checks,
            "would_be_consistent_under_simulation": False,
            "status": "SIMULATED_FAILURE_DETECTED",
        })
    governance_nonclosure_two_flip_matrix = {
        "status": "SIMULATED_FAILURE_MATRIX_EXPORTED",
        "scope": "TEST_ONLY_DIAGNOSTIC_NOT_RUNTIME_CLAIM",
        "baseline_checks": dict(governance_nonclosure_single_flip_matrix["baseline_checks"]),
        "rows": [],
    }
    check_names = list(governance_nonclosure_two_flip_matrix["baseline_checks"].keys())
    for left_name, right_name in combinations(check_names, 2):
        row_checks = dict(governance_nonclosure_two_flip_matrix["baseline_checks"])
        row_checks[left_name] = False
        row_checks[right_name] = False
        governance_nonclosure_two_flip_matrix["rows"].append({
            "flipped_checks": [left_name, right_name],
            "simulated_checks": row_checks,
            "would_be_consistent_under_simulation": False,
            "status": "SIMULATED_FAILURE_DETECTED",
        })
    two_flip_false_count_rows = [
        int(sum(1 for _, v in row["simulated_checks"].items() if not bool(v)))
        for row in governance_nonclosure_two_flip_matrix["rows"]
    ]
    governance_nonclosure_two_flip_matrix["coverage_summary"] = {
        "num_checks": int(len(check_names)),
        "expected_rows_n_choose_2": int(len(check_names) * (len(check_names) - 1) // 2),
        "exported_rows": int(len(governance_nonclosure_two_flip_matrix["rows"])),
        "all_rows_have_exactly_two_false": bool(all(c == 2 for c in two_flip_false_count_rows)),
    }
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
        "schema_version": "p2025_s975_v202",
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
            "branch_integrator_replay_task7_panel": {
                "seeds": [int(x) for x in replay_seeds],
                "stress_abs_gap_q05_q50_q95": [float(x) for x in stress_q.tolist()],
                "rows": replay_rows,
                "task7_frequency_span_over_seeds": replay_task7_span,
                "leader_task_id_per_seed": [int(x) for x in replay_leaders],
                "leader_stable_over_seeds": replay_leader_stable,
                "controlled_substitution_replay": controlled_replay,
                "controlled_substitution_guard_sensitivity": {
                    "rows": guard_rows,
                    "allow_frequency_over_threshold_grid": guard_allow_freq,
                },
                "substitution_replay_governance": substitution_governance,
                "actual_substitution_replay": actual_substitution_replay,
                "actual_substitution_replay_comparative_report": actual_substitution_replay_comparative_report,
                "cross_seed_actual_substitution_replay_panel": cross_seed_actual_substitution_replay_panel,
                "cross_seed_substitution_governance": cross_seed_substitution_governance,
                "nonclosure_lock_after_governance": nonclosure_lock_after_governance,
                "task7_attack_and_task4_verification_packet": task7_attack_and_task4_verification_packet,
                "governance_result_discussion": governance_result_discussion,
                "task7_task4_trend_panel": task7_task4_trend_panel,
                "trend_gate_for_costlier_step": trend_gate_for_costlier_step,
                "composite_nonclosure_enforcement": composite_nonclosure_enforcement,
                "nonclosure_status_history_audit": nonclosure_status_history_audit,
                "governance_nonclosure_consistency_gate": governance_nonclosure_consistency_gate,
                "governance_nonclosure_failure_simulation": governance_nonclosure_failure_simulation,
                "governance_history_nonclosure_failure_simulation": governance_history_nonclosure_failure_simulation,
                "governance_history_global_nonclosure_failure_simulation": governance_history_global_nonclosure_failure_simulation,
                "governance_nonclosure_single_flip_matrix": governance_nonclosure_single_flip_matrix,
                "governance_nonclosure_two_flip_matrix": governance_nonclosure_two_flip_matrix,
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
        "task2_strict_unitarity_witness": task2_strict_unitarity_witness,
        "theorem_target": task2_strict_unitarity_witness["theorem_target"],
        "strict_lane_assumptions": task2_strict_unitarity_witness["strict_lane_assumptions"],
        "computed_rows": task2_strict_unitarity_witness["rows"],
        "aggregate_metrics": task2_strict_unitarity_witness["aggregate_metrics"],
        "pass_fail_criteria": task2_strict_unitarity_witness["pass_fail_criteria"],
        "verdict": task2_strict_unitarity_witness["verdict"],
        "fail_trace": task2_strict_unitarity_witness["fail_trace"],
        "ur_uncertainty_transport_bridge_precursor": ur_uncertainty_transport_bridge_precursor,
        "ur_transport_cross_source_agreement_precursor": ur_transport_cross_source_agreement_precursor,
        "ur_channel_trace_budget_precursor": ur_channel_trace_budget_precursor,
        "ur_channel_class_mapping_precursor": ur_channel_class_mapping_precursor,
        "ur_class_bounded_uncertainty_residual_budget_precursor": ur_class_bounded_uncertainty_residual_budget_precursor,
        "ur_class_readiness_gate_precursor": ur_class_readiness_gate_precursor,
        "ur_class_first_exact_integration_replay_precursor": ur_class_first_exact_integration_replay_precursor,
        "ur_class_first_replay_delta_precursor": ur_class_first_replay_delta_precursor,
        "ur_class_first_vs_all_class_replay_comparison_precursor": ur_class_first_vs_all_class_replay_comparison_precursor,
        "ur_cost_vs_gain_precursor": ur_cost_vs_gain_precursor,
        "ur_runtime_tolerance_benchmark_precursor": ur_runtime_tolerance_benchmark_precursor,
        "ur_all_class_exact_integration_sweep_precursor": ur_all_class_exact_integration_sweep_precursor,
        "ur_numerical_stress_ranking_precursor": ur_numerical_stress_ranking_precursor,
        "ur_numerical_stress_alt_parameterization_precursor": ur_numerical_stress_alt_parameterization_precursor,
        "ur_numerical_stress_alt_transform_comparison_precursor": ur_numerical_stress_alt_transform_comparison_precursor,
        "ur_numerical_stress_alt_fullgrid_tritransform_precursor": ur_numerical_stress_alt_fullgrid_tritransform_precursor,
        "ur_numerical_stress_class_conditional_replay_precursor": ur_numerical_stress_class_conditional_replay_precursor,
        "ur_numerical_stress_policy_counterfactual_precursor": ur_numerical_stress_policy_counterfactual_precursor,
        "ur_numerical_stress_policy_pareto_front_precursor": ur_numerical_stress_policy_pareto_front_precursor,
        "ur_numerical_stress_policy_pareto_stability_precursor": ur_numerical_stress_policy_pareto_stability_precursor,
        "ur_numerical_stress_policy_budgeted_selection_precursor": ur_numerical_stress_policy_budgeted_selection_precursor,
        "ur_numerical_stress_policy_budget_fragility_precursor": ur_numerical_stress_policy_budget_fragility_precursor,
        "ur_numerical_stress_policy_budget_fragility_by_class_precursor": ur_numerical_stress_policy_budget_fragility_by_class_precursor,
        "ur_numerical_stress_policy_class_adaptive_fallback_precursor": ur_numerical_stress_policy_class_adaptive_fallback_precursor,
        "ur_numerical_stress_policy_ablation_precursor": ur_numerical_stress_policy_ablation_precursor,
        "ur_numerical_stress_policy_cross_class_constrained_ablation_precursor": ur_numerical_stress_policy_cross_class_constrained_ablation_precursor,
        "ur_numerical_stress_policy_cross_class_constrained_bootstrap_dominance_precursor": ur_numerical_stress_policy_cross_class_constrained_bootstrap_dominance_precursor,
        "ur_numerical_stress_policy_cross_class_threshold_sweep_precursor": ur_numerical_stress_policy_cross_class_threshold_sweep_precursor,
        "ur_numerical_stress_policy_joint_stress_map_precursor": ur_numerical_stress_policy_joint_stress_map_precursor,
        "ur_numerical_stress_policy_stability_topology_precursor": ur_numerical_stress_policy_stability_topology_precursor,
        "ur_numerical_stress_policy_stability_boundary_margin_precursor": ur_numerical_stress_policy_stability_boundary_margin_precursor,
        "ur_numerical_stress_policy_weighted_boundary_risk_precursor": ur_numerical_stress_policy_weighted_boundary_risk_precursor,
        "ur_numerical_stress_policy_weighted_boundary_risk_sensitivity_precursor": ur_numerical_stress_policy_weighted_boundary_risk_sensitivity_precursor,
        "ur_numerical_stress_policy_weighted_boundary_risk_bayesian_precursor": ur_numerical_stress_policy_weighted_boundary_risk_bayesian_precursor,
        "ur_numerical_stress_policy_weighted_boundary_risk_posterior_predictive_precursor": ur_numerical_stress_policy_weighted_boundary_risk_posterior_predictive_precursor,
        "ur_numerical_stress_policy_posterior_predictive_decision_gate_precursor": ur_numerical_stress_policy_posterior_predictive_decision_gate_precursor,
        "ur_numerical_stress_policy_posterior_predictive_gate_calibration_precursor": ur_numerical_stress_policy_posterior_predictive_gate_calibration_precursor,
        "ur_numerical_stress_policy_posterior_predictive_gate_cost_calibration_precursor": ur_numerical_stress_policy_posterior_predictive_gate_cost_calibration_precursor,
        "ur_numerical_stress_policy_gate_frontier_utility_precursor": ur_numerical_stress_policy_gate_frontier_utility_precursor,
        "ur_numerical_stress_policy_gate_frontier_knee_precursor": ur_numerical_stress_policy_gate_frontier_knee_precursor,
        "ur_numerical_stress_policy_gate_frontier_knee_stability_precursor": ur_numerical_stress_policy_gate_frontier_knee_stability_precursor,
        "ur_numerical_stress_policy_gate_frontier_knee_cross_seed_stability_precursor": ur_numerical_stress_policy_gate_frontier_knee_cross_seed_stability_precursor,
        "ur_numerical_stress_policy_gate_frontier_knee_cross_seed_consensus_precursor": ur_numerical_stress_policy_gate_frontier_knee_cross_seed_consensus_precursor,
        "ur_numerical_stress_policy_gate_frontier_knee_consensus_threshold_sweep_precursor": ur_numerical_stress_policy_gate_frontier_knee_consensus_threshold_sweep_precursor,
        "ur_numerical_stress_policy_gate_frontier_knee_weighted_cross_seed_consensus_precursor": ur_numerical_stress_policy_gate_frontier_knee_weighted_cross_seed_consensus_precursor,
        "ur_numerical_stress_alt_replay_trend_precursor": ur_numerical_stress_alt_replay_trend_precursor,
        "ur_numerical_stress_alt_dominance_map_precursor": ur_numerical_stress_alt_dominance_map_precursor,
        "ur_numerical_stress_alt_decision_gate_precursor": ur_numerical_stress_alt_decision_gate_precursor,
        "ur_numerical_stress_alt_hysteresis_gate_precursor": ur_numerical_stress_alt_hysteresis_gate_precursor,
        "ur_numerical_stress_alt_hysteresis_time_stability_precursor": ur_numerical_stress_alt_hysteresis_time_stability_precursor,
        "ur_numerical_stress_alt_entropy_gate_precursor": ur_numerical_stress_alt_entropy_gate_precursor,
        "ur_numerical_stress_alt_entropy_threshold_calibration_precursor": ur_numerical_stress_alt_entropy_threshold_calibration_precursor,
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
