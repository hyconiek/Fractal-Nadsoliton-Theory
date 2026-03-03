#!/usr/bin/env python3
"""
QW-2044: Pointwise micro-derivation of Z_beta(d) and delta_eta(d).

Objective:
- derive local damping parameters directly from microstate profile distributions,
- without sector-specific retuning and without using external mass/flavor/GW targets
  in the fit stage.

Model used on local envelope windows:
    m(d) ~ A / (1 + beta * d**eta)

Then:
    Z_beta(d) = beta(d) / 0.01
    delta_eta(d) = eta(d) - 1
"""

from __future__ import annotations

import importlib.util
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2044_pointwise_micro_zbeta_deltaeta_derivation.json"
OUT_MD = ROOT / "RAPORT_QW2044_POINTWISE_MICRO_ZBETA_DELTAETA_DERIVATION.md"


def load_qw1739_module():
    path = ROOT / "QW_1739_SIGNED_DYNAMIC_MICROMODEL_DERIVATION.py"
    spec = importlib.util.spec_from_file_location("qw1739_mod_2044", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["qw1739_mod_2044"] = mod
    spec.loader.exec_module(mod)
    return mod


def smooth_abs(y: np.ndarray) -> np.ndarray:
    a = np.abs(y)
    if len(a) < 3:
        return a
    out = a.copy()
    out[1:-1] = 0.25 * a[:-2] + 0.5 * a[1:-1] + 0.25 * a[2:]
    return out


def extract_peak_indices(y: np.ndarray) -> np.ndarray:
    a = smooth_abs(y)
    n = len(a)
    idx = []
    for i in range(n):
        l = a[i - 1] if i - 1 >= 0 else -np.inf
        r = a[i + 1] if i + 1 < n else -np.inf
        if a[i] >= l and a[i] >= r:
            idx.append(i)

    idx = np.array(sorted(set(idx)), dtype=int)

    # remove adjacent duplicates by amplitude ordering
    if len(idx) > 1:
        keep = []
        last = -10
        for i in idx:
            if i - last <= 1:
                if a[i] > a[last]:
                    if keep:
                        keep.pop()
                    keep.append(i)
                    last = i
            else:
                keep.append(i)
                last = i
        idx = np.array(keep, dtype=int)

    # if too few peaks, fallback to top amplitudes spread over d
    if len(idx) < 4:
        top = np.argsort(-a)[: min(8, len(a))]
        top = np.sort(top)
        idx = top.astype(int)

    return idx


def fit_local_envelope(d: np.ndarray, m: np.ndarray) -> Dict[str, float]:
    # Linearization for each eta:
    # 1/m = u + v d^eta, u=1/A, v=beta/A => A=1/u, beta=v/u
    eta_grid = np.linspace(1.0, 3.0, 161)
    best = None

    inv_m = 1.0 / np.clip(m, 1e-9, None)

    for eta in eta_grid:
        x = d**eta
        mat = np.column_stack([np.ones(len(d), dtype=float), x])
        coef, *_ = np.linalg.lstsq(mat, inv_m, rcond=None)
        u, v = float(coef[0]), float(coef[1])
        if u <= 0.0 or v <= 0.0:
            continue

        a_hat = 1.0 / u
        beta_hat = v / u
        pred = a_hat / (1.0 + beta_hat * (d**eta))
        rmse = float(np.sqrt(np.mean((m - pred) ** 2)))

        if best is None or rmse < best["rmse"]:
            best = {
                "A": float(a_hat),
                "beta": float(beta_hat),
                "eta": float(eta),
                "rmse": rmse,
            }

    if best is None:
        return {"A": float("nan"), "beta": float("nan"), "eta": float("nan"), "rmse": float("inf")}
    return best


def derive_pointwise_windows(y: np.ndarray) -> List[Dict[str, float]]:
    idx = extract_peak_indices(y)
    if len(idx) < 4:
        return []

    a = smooth_abs(y)
    rows: List[Dict[str, float]] = []

    # Sliding windows over peaks to get local beta(d), eta(d)
    for j in range(0, len(idx) - 3):
        ids = idx[j : j + 4]
        d = (ids + 1).astype(float)
        m = a[ids].astype(float)

        # basic quality filters
        if bool(np.any(m <= 1e-6)):
            continue
        if float(np.std(d)) <= 0.25:
            continue

        fit = fit_local_envelope(d, m)
        if not np.isfinite(float(fit["beta"])) or not np.isfinite(float(fit["eta"])):
            continue

        rows.append(
            {
                "d_center": float(np.mean(d)),
                "d_min": float(np.min(d)),
                "d_max": float(np.max(d)),
                "beta": float(fit["beta"]),
                "eta": float(fit["eta"]),
                "A": float(fit["A"]),
                "rmse": float(fit["rmse"]),
            }
        )

    return rows


def build_profiles(n_grid: List[int], seeds_per_n: int, dmax: int) -> List[Dict[str, object]]:
    mod = load_qw1739_module()
    out: List[Dict[str, object]] = []

    for n in n_grid:
        for sidx in range(seeds_per_n):
            seed = 204400 + 1000 * n + sidx
            cfg = mod.MicroCfg(n_nodes=n, seed=seed)
            rng = np.random.default_rng(seed)

            theta, q = mod.build_micro_state(cfg, rng)
            w = mod.build_signed_dynamic_matrix(cfg, theta, q, rng)
            g = mod.effective_path_sum(w, lam=cfg.dyn_lambda, kmax=cfg.dyn_kmax)
            y = mod.profile_from_matrix(g, dmax=min(int(dmax), n // 2))

            rows = derive_pointwise_windows(np.array(y, dtype=float))
            out.append(
                {
                    "n_nodes": int(n),
                    "seed": int(seed),
                    "n_windows": int(len(rows)),
                    "rows": rows,
                }
            )
    return out


def summarize_bins(rows: List[Dict[str, float]], beta_target: float, eta_target: float) -> Dict[str, object]:
    bins: Dict[int, List[Dict[str, float]]] = {}
    for r in rows:
        b = int(np.clip(round(float(r["d_center"])), 1, 24))
        bins.setdefault(b, []).append(r)

    out_bins = []
    for b in sorted(bins):
        arr = bins[b]
        if len(arr) < 6:
            continue

        beta_vals = np.array([float(x["beta"]) for x in arr], dtype=float)
        eta_vals = np.array([float(x["eta"]) for x in arr], dtype=float)
        rmse_vals = np.array([float(x["rmse"]) for x in arr], dtype=float)

        beta_ci = [float(np.quantile(beta_vals, 0.025)), float(np.quantile(beta_vals, 0.975))]
        eta_ci = [float(np.quantile(eta_vals, 0.025)), float(np.quantile(eta_vals, 0.975))]

        z_vals = beta_vals / 0.01
        de_vals = eta_vals - 1.0

        out_bins.append(
            {
                "d_bin": int(b),
                "n": int(len(arr)),
                "beta_median": float(np.median(beta_vals)),
                "beta_ci95": beta_ci,
                "eta_median": float(np.median(eta_vals)),
                "eta_ci95": eta_ci,
                "z_beta_median": float(np.median(z_vals)),
                "delta_eta_median": float(np.median(de_vals)),
                "rmse_median": float(np.median(rmse_vals)),
                "target_beta_in_ci95": bool(beta_ci[0] <= beta_target <= beta_ci[1]),
                "target_eta_in_ci95": bool(eta_ci[0] <= eta_target <= eta_ci[1]),
            }
        )

    return {"bins": out_bins, "n_bins": int(len(out_bins))}


def main() -> None:
    d2039 = json.loads((ROOT / "report_qw2039_derivation_compatible_refrozen_kernel_gate.json").read_text(encoding="utf-8"))
    beta_target = float(d2039["selected_kernel"]["beta"])
    eta_target = float(d2039["selected_kernel"]["eta"])

    profiles = build_profiles(n_grid=[64, 96, 128], seeds_per_n=18, dmax=24)

    all_rows: List[Dict[str, float]] = []
    for p in profiles:
        all_rows.extend(p["rows"])

    if len(all_rows) < 60:
        raise RuntimeError("Too few pointwise windows for stable inference.")

    beta_all = np.array([float(r["beta"]) for r in all_rows], dtype=float)
    eta_all = np.array([float(r["eta"]) for r in all_rows], dtype=float)
    rmse_all = np.array([float(r["rmse"]) for r in all_rows], dtype=float)

    global_est = {
        "n_profiles": int(len(profiles)),
        "n_windows_total": int(len(all_rows)),
        "beta_median": float(np.median(beta_all)),
        "beta_ci95": [float(np.quantile(beta_all, 0.025)), float(np.quantile(beta_all, 0.975))],
        "eta_median": float(np.median(eta_all)),
        "eta_ci95": [float(np.quantile(eta_all, 0.025)), float(np.quantile(eta_all, 0.975))],
        "z_beta_median": float(np.median(beta_all / 0.01)),
        "delta_eta_median": float(np.median(eta_all - 1.0)),
        "rmse_median": float(np.median(rmse_all)),
    }

    binned = summarize_bins(all_rows, beta_target=beta_target, eta_target=eta_target)
    bins = binned["bins"]

    frac_beta = float(np.mean([1.0 if b["target_beta_in_ci95"] else 0.0 for b in bins])) if bins else 0.0
    frac_eta = float(np.mean([1.0 if b["target_eta_in_ci95"] else 0.0 for b in bins])) if bins else 0.0
    frac_joint = (
        float(np.mean([1.0 if (b["target_beta_in_ci95"] and b["target_eta_in_ci95"]) else 0.0 for b in bins]))
        if bins
        else 0.0
    )

    flags = {
        "enough_pointwise_bins_ge_6": bool(int(binned["n_bins"]) >= 6),
        "global_beta_target_inside_ci95": bool(global_est["beta_ci95"][0] <= beta_target <= global_est["beta_ci95"][1]),
        "global_eta_target_inside_ci95": bool(global_est["eta_ci95"][0] <= eta_target <= global_est["eta_ci95"][1]),
        "bin_beta_target_coverage_ge_0p50": bool(frac_beta >= 0.50),
        "bin_eta_target_coverage_ge_0p50": bool(frac_eta >= 0.50),
        "bin_joint_target_coverage_ge_0p35": bool(frac_joint >= 0.35),
        "median_window_rmse_le_0p12": bool(global_est["rmse_median"] <= 0.12),
    }

    pass_count = int(sum(1 for v in flags.values() if v))
    total_flags = int(len(flags))

    if pass_count == total_flags:
        verdict = "POINTWISE_MICRO_DERIVATION_PASS"
        readiness = "POINTWISE_MICRO_ZBETA_DELTAETA_SUPPORTED"
    elif pass_count >= 5:
        verdict = "POINTWISE_MICRO_DERIVATION_PARTIAL"
        readiness = "PARTIAL_POINTWISE_SUPPORT"
    else:
        verdict = "POINTWISE_MICRO_DERIVATION_FAIL"
        readiness = "POINTWISE_SUPPORT_NOT_ESTABLISHED"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "target_refrozen_kernel": {
            "beta": beta_target,
            "eta": eta_target,
            "z_beta_target": float(beta_target / 0.01),
            "delta_eta_target": float(eta_target - 1.0),
        },
        "global_estimates": global_est,
        "pointwise_bins": binned,
        "coverage": {
            "beta_target_in_ci95_fraction": frac_beta,
            "eta_target_in_ci95_fraction": frac_eta,
            "joint_target_in_ci95_fraction": frac_joint,
        },
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "readiness": readiness,
        "required_next_step": (
            "RUN_MICRO_COMPATIBILITY_INTERSECTION_GATE_WITH_STAGEC_PASS_SET"
            if pass_count >= 5
            else "IMPROVE_POINTWISE_IDENTIFIABILITY_WITH_PHASE_CONDITIONED_WINDOWS"
        ),
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2044: POINTWISE MICRO Z_BETA(d) / DELTA_ETA(d) DERIVATION",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- Readiness: **{readiness}**",
        f"- pass_count: {pass_count}/{total_flags}",
        "",
        "## Target (Refrozen QW-2039)",
        f"- beta_target: {beta_target:.6f}",
        f"- eta_target: {eta_target:.6f}",
        f"- Z_beta_target: {beta_target/0.01:.6f}",
        f"- delta_eta_target: {eta_target-1.0:.6f}",
        "",
        "## Global Micro Estimates",
        f"- n_profiles: {global_est['n_profiles']}",
        f"- n_windows_total: {global_est['n_windows_total']}",
        f"- beta median/CI95: {global_est['beta_median']:.6f} / [{global_est['beta_ci95'][0]:.6f}, {global_est['beta_ci95'][1]:.6f}]",
        f"- eta median/CI95: {global_est['eta_median']:.6f} / [{global_est['eta_ci95'][0]:.6f}, {global_est['eta_ci95'][1]:.6f}]",
        f"- Z_beta median: {global_est['z_beta_median']:.6f}",
        f"- delta_eta median: {global_est['delta_eta_median']:.6f}",
        f"- median window rmse: {global_est['rmse_median']:.6f}",
        "",
        "## Pointwise Coverage",
        f"- n_bins: {binned['n_bins']}",
        f"- beta target in CI95 fraction: {frac_beta:.4f}",
        f"- eta target in CI95 fraction: {frac_eta:.4f}",
        f"- joint target in CI95 fraction: {frac_joint:.4f}",
        "",
        "## Flags",
    ]
    for k, v in flags.items():
        lines.append(f"- {k}: {v}")

    top_bins = sorted(
        bins,
        key=lambda b: (
            -(1 if b["target_beta_in_ci95"] and b["target_eta_in_ci95"] else 0),
            b["rmse_median"],
        ),
    )[:8]
    lines.extend(["", "## Representative Bins (top-8)"])
    for b in top_bins:
        lines.append(
            "- d={d} n={n} | beta={bm:.4f} [{b0:.4f},{b1:.4f}] | "
            "eta={em:.4f} [{e0:.4f},{e1:.4f}] | rmse={r:.4f} | joint={j}".format(
                d=b["d_bin"],
                n=b["n"],
                bm=b["beta_median"],
                b0=b["beta_ci95"][0],
                b1=b["beta_ci95"][1],
                em=b["eta_median"],
                e0=b["eta_ci95"][0],
                e1=b["eta_ci95"][1],
                r=b["rmse_median"],
                j=(b["target_beta_in_ci95"] and b["target_eta_in_ci95"]),
            )
        )

    lines.extend(
        [
            "",
            "## Required Next Step",
            f"- {out['required_next_step']}",
            "",
            "## Artifacts",
            f"- JSON: `{OUT_JSON.name}`",
        ]
    )

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2044] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2044] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2044] verdict={verdict} pass_count={pass_count}/{total_flags}")


if __name__ == "__main__":
    main()
