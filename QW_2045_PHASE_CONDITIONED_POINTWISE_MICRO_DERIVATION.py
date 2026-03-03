#!/usr/bin/env python3
"""
QW-2045: Phase-conditioned pointwise micro-derivation of Z_beta(d), delta_eta(d).

Improvement over QW-2044:
- explicitly conditions local windows on high phase observability
  (|cos(omega d + phi)| near maxima),
- deconvolves phase magnitude from |y(d)| before damping fit.
"""

from __future__ import annotations

import importlib.util
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2045_phase_conditioned_pointwise_micro_derivation.json"
OUT_MD = ROOT / "RAPORT_QW2045_PHASE_CONDITIONED_POINTWISE_MICRO_DERIVATION.md"


def load_qw1739_module():
    path = ROOT / "QW_1739_SIGNED_DYNAMIC_MICROMODEL_DERIVATION.py"
    spec = importlib.util.spec_from_file_location("qw1739_mod_2045", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["qw1739_mod_2045"] = mod
    spec.loader.exec_module(mod)
    return mod


def fit_local_envelope(d: np.ndarray, m: np.ndarray) -> Dict[str, float]:
    eta_grid = np.linspace(1.0, 3.0, 161)
    inv_m = 1.0 / np.clip(m, 1e-9, None)
    best = None

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


def derive_phase_conditioned_rows(y: np.ndarray, mod) -> List[Dict[str, float]]:
    d = np.arange(1, len(y) + 1, dtype=float)

    # phase proxy from QW-1739 (no external targets used)
    est = mod.derive_omega_phi_beta(np.array(y, dtype=float))
    om = float(est["omega_hat"])
    ph = float(est["phi_hat"])

    phase_mag = np.abs(np.cos(om * d + ph))
    idx = np.array([], dtype=int)
    # Adaptive phase gate: keep strongest conditioning first, relax only if needed.
    for thr in [0.80, 0.75, 0.70, 0.65, 0.60]:
        sel = np.where(phase_mag >= thr)[0]
        if len(sel) >= 5:
            idx = sel
            break
    if len(idx) < 5:
        thr_q = float(np.quantile(phase_mag, 0.55))
        idx = np.where(phase_mag >= max(thr_q, 0.50))[0]
    if len(idx) < 5:
        return []

    amp = np.abs(y) / np.clip(phase_mag, 0.20, None)

    rows: List[Dict[str, float]] = []
    win = 4
    for j in range(0, len(idx) - win + 1):
        ids = idx[j : j + win]
        d_w = d[ids]
        m_w = amp[ids]

        # robustness filters
        if bool(np.any(m_w <= 1e-6)):
            continue
        if float(np.std(d_w)) <= 0.4:
            continue

        fit = fit_local_envelope(d_w, m_w)
        if not np.isfinite(float(fit["beta"])) or not np.isfinite(float(fit["eta"])):
            continue

        rows.append(
            {
                "d_center": float(np.mean(d_w)),
                "d_min": float(np.min(d_w)),
                "d_max": float(np.max(d_w)),
                "beta": float(fit["beta"]),
                "eta": float(fit["eta"]),
                "A": float(fit["A"]),
                "rmse": float(fit["rmse"]),
                "omega_hat": om,
                "phi_hat": ph,
                "phase_min": float(np.min(phase_mag[ids])),
            }
        )

    return rows


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
        phm_vals = np.array([float(x["phase_min"]) for x in arr], dtype=float)

        beta_ci = [float(np.quantile(beta_vals, 0.025)), float(np.quantile(beta_vals, 0.975))]
        eta_ci = [float(np.quantile(eta_vals, 0.025)), float(np.quantile(eta_vals, 0.975))]

        out_bins.append(
            {
                "d_bin": int(b),
                "n": int(len(arr)),
                "beta_median": float(np.median(beta_vals)),
                "beta_ci95": beta_ci,
                "eta_median": float(np.median(eta_vals)),
                "eta_ci95": eta_ci,
                "z_beta_median": float(np.median(beta_vals / 0.01)),
                "delta_eta_median": float(np.median(eta_vals - 1.0)),
                "rmse_median": float(np.median(rmse_vals)),
                "phase_min_median": float(np.median(phm_vals)),
                "target_beta_in_ci95": bool(beta_ci[0] <= beta_target <= beta_ci[1]),
                "target_eta_in_ci95": bool(eta_ci[0] <= eta_target <= eta_ci[1]),
            }
        )

    return {"bins": out_bins, "n_bins": int(len(out_bins))}


def main() -> None:
    mod = load_qw1739_module()
    d2039 = json.loads((ROOT / "report_qw2039_derivation_compatible_refrozen_kernel_gate.json").read_text(encoding="utf-8"))
    beta_target = float(d2039["selected_kernel"]["beta"])
    eta_target = float(d2039["selected_kernel"]["eta"])

    profiles = []
    all_rows: List[Dict[str, float]] = []

    for n in [64, 96, 128]:
        for sidx in range(22):
            seed = 204500 + 1000 * n + sidx
            cfg = mod.MicroCfg(n_nodes=n, seed=seed)
            rng = np.random.default_rng(seed)
            theta, q = mod.build_micro_state(cfg, rng)
            w = mod.build_signed_dynamic_matrix(cfg, theta, q, rng)
            g = mod.effective_path_sum(w, lam=cfg.dyn_lambda, kmax=cfg.dyn_kmax)
            y = mod.profile_from_matrix(g, dmax=min(24, n // 2))

            rows = derive_phase_conditioned_rows(np.array(y, dtype=float), mod)
            profiles.append({"n_nodes": int(n), "seed": int(seed), "n_rows": int(len(rows))})
            all_rows.extend(rows)

    if len(all_rows) < 60:
        raise RuntimeError("Too few phase-conditioned rows for robust pointwise inference.")

    beta_all = np.array([float(r["beta"]) for r in all_rows], dtype=float)
    eta_all = np.array([float(r["eta"]) for r in all_rows], dtype=float)
    rmse_all = np.array([float(r["rmse"]) for r in all_rows], dtype=float)
    phase_all = np.array([float(r["phase_min"]) for r in all_rows], dtype=float)

    global_est = {
        "n_profiles": int(len(profiles)),
        "n_rows_total": int(len(all_rows)),
        "beta_median": float(np.median(beta_all)),
        "beta_ci95": [float(np.quantile(beta_all, 0.025)), float(np.quantile(beta_all, 0.975))],
        "eta_median": float(np.median(eta_all)),
        "eta_ci95": [float(np.quantile(eta_all, 0.025)), float(np.quantile(eta_all, 0.975))],
        "z_beta_median": float(np.median(beta_all / 0.01)),
        "delta_eta_median": float(np.median(eta_all - 1.0)),
        "rmse_median": float(np.median(rmse_all)),
        "phase_min_median": float(np.median(phase_all)),
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
        "median_window_rmse_le_0p10": bool(global_est["rmse_median"] <= 0.10),
        "phase_condition_strength_ge_0p75": bool(global_est["phase_min_median"] >= 0.75),
    }

    pass_count = int(sum(1 for v in flags.values() if v))
    total_flags = int(len(flags))

    if pass_count == total_flags:
        verdict = "PHASE_CONDITIONED_POINTWISE_DERIVATION_PASS"
        readiness = "POINTWISE_IDENTIFIABILITY_REPAIRED"
    elif pass_count >= 6:
        verdict = "PHASE_CONDITIONED_POINTWISE_DERIVATION_PARTIAL"
        readiness = "PARTIAL_IDENTIFIABILITY_REPAIR"
    else:
        verdict = "PHASE_CONDITIONED_POINTWISE_DERIVATION_FAIL"
        readiness = "POINTWISE_IDENTIFIABILITY_STILL_INSUFFICIENT"

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
            "RUN_MICRO_STAGEC_INTERSECTION_GATE"
            if pass_count >= 6
            else "ADD_SIGNED_PHASE_TORSION_OBSERVABLE_TO_BREAK_REMAINING_DEGENERACY"
        ),
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2045: PHASE-CONDITIONED POINTWISE MICRO DERIVATION",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- Readiness: **{readiness}**",
        f"- pass_count: {pass_count}/{total_flags}",
        "",
        "## Target (Refrozen QW-2039)",
        f"- beta_target: {beta_target:.6f}",
        f"- eta_target: {eta_target:.6f}",
        "",
        "## Global Estimates",
        f"- n_profiles: {global_est['n_profiles']}",
        f"- n_rows_total: {global_est['n_rows_total']}",
        f"- beta median/CI95: {global_est['beta_median']:.6f} / [{global_est['beta_ci95'][0]:.6f}, {global_est['beta_ci95'][1]:.6f}]",
        f"- eta median/CI95: {global_est['eta_median']:.6f} / [{global_est['eta_ci95'][0]:.6f}, {global_est['eta_ci95'][1]:.6f}]",
        f"- rmse median: {global_est['rmse_median']:.6f}",
        f"- phase_min median: {global_est['phase_min_median']:.6f}",
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
            "eta={em:.4f} [{e0:.4f},{e1:.4f}] | rmse={r:.4f} | phase={ph:.3f} | joint={j}".format(
                d=b["d_bin"],
                n=b["n"],
                bm=b["beta_median"],
                b0=b["beta_ci95"][0],
                b1=b["beta_ci95"][1],
                em=b["eta_median"],
                e0=b["eta_ci95"][0],
                e1=b["eta_ci95"][1],
                r=b["rmse_median"],
                ph=b["phase_min_median"],
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

    print(f"[QW-2045] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2045] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2045] verdict={verdict} pass_count={pass_count}/{total_flags}")


if __name__ == "__main__":
    main()
