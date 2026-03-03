#!/usr/bin/env python3
"""
QW-2047: Signed phase-torsion observable scan for pointwise identifiability repair.

Adds an internal observable built from:
- phase-aligned micro profile,
- antisymmetry/torsion index of the signed dynamic micro matrix.

Scans lambda_tau to test whether this observable improves pointwise
Z_beta(d), delta_eta(d) identifiability without sector retune.
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
OUT_JSON = ROOT / "report_qw2047_signed_phase_torsion_observable_scan.json"
OUT_MD = ROOT / "RAPORT_QW2047_SIGNED_PHASE_TORSION_OBSERVABLE_SCAN.md"


def load_qw1739_module():
    path = ROOT / "QW_1739_SIGNED_DYNAMIC_MICROMODEL_DERIVATION.py"
    spec = importlib.util.spec_from_file_location("qw1739_mod_2047", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["qw1739_mod_2047"] = mod
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
            best = {"A": float(a_hat), "beta": float(beta_hat), "eta": float(eta), "rmse": rmse}

    if best is None:
        return {"A": float("nan"), "beta": float("nan"), "eta": float("nan"), "rmse": float("inf")}
    return best


def torsion_index_from_matrix(w: np.ndarray) -> float:
    ws = 0.5 * (w + w.T)
    wa = 0.5 * (w - w.T)
    nrm_s = float(np.linalg.norm(ws))
    nrm_a = float(np.linalg.norm(wa))
    return float(nrm_a / max(nrm_s + nrm_a, 1e-12))


def derive_rows_for_lambda(y: np.ndarray, omega_hat: float, phi_hat: float, tau_idx: float, lam_tau: float) -> List[Dict[str, float]]:
    d = np.arange(1, len(y) + 1, dtype=float)
    phase = np.cos(omega_hat * d + phi_hat)
    phase_mag = np.abs(phase)

    # signed phase alignment
    aligned = y * np.sign(phase + 1e-12)
    # torsion-dependent shaping (internal, profile-level)
    torsion_shape = np.exp(float(lam_tau) * tau_idx * (d / max(float(np.max(d)), 1e-12) - 0.5))
    # Use signed alignment to reduce phase flips, then positive envelope for local damping fit.
    obs = np.abs(aligned) * torsion_shape / np.clip(phase_mag, 0.25, None)

    idx = np.array([], dtype=int)
    for thr in [0.70, 0.60, 0.50, 0.40]:
        sel = np.where(phase_mag >= thr)[0]
        if len(sel) >= 5:
            idx = sel
            break
    if len(idx) < 5:
        thr_q = float(np.quantile(phase_mag, 0.55))
        idx = np.where(phase_mag >= max(thr_q, 0.35))[0]
    if len(idx) < 5:
        return []

    win = 4
    rows: List[Dict[str, float]] = []
    for j in range(0, len(idx) - win + 1):
        ids = idx[j : j + win]
        d_w = d[ids]
        m_w = obs[ids]

        if bool(np.any(m_w <= 1e-8)):
            continue
        if float(np.std(d_w)) <= 0.35:
            continue

        fit = fit_local_envelope(d_w, m_w)
        if not np.isfinite(float(fit["beta"])) or not np.isfinite(float(fit["eta"])):
            continue

        rows.append(
            {
                "d_center": float(np.mean(d_w)),
                "beta": float(fit["beta"]),
                "eta": float(fit["eta"]),
                "rmse": float(fit["rmse"]),
                "phase_min": float(np.min(phase_mag[ids])),
            }
        )

    return rows


def summarize(rows: List[Dict[str, float]], beta_target: float, eta_target: float) -> Dict[str, object]:
    bins = {}
    for r in rows:
        b = int(np.clip(round(float(r["d_center"])), 1, 24))
        bins.setdefault(b, []).append(r)

    out_bins = []
    for b in sorted(bins):
        arr = bins[b]
        if len(arr) < 6:
            continue
        beta_vals = np.array([x["beta"] for x in arr], dtype=float)
        eta_vals = np.array([x["eta"] for x in arr], dtype=float)
        rmse_vals = np.array([x["rmse"] for x in arr], dtype=float)
        phase_vals = np.array([x["phase_min"] for x in arr], dtype=float)

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
                "rmse_median": float(np.median(rmse_vals)),
                "phase_min_median": float(np.median(phase_vals)),
                "target_beta_in_ci95": bool(beta_ci[0] <= beta_target <= beta_ci[1]),
                "target_eta_in_ci95": bool(eta_ci[0] <= eta_target <= eta_ci[1]),
            }
        )

    n_bins = int(len(out_bins))
    frac_beta = float(np.mean([1.0 if b["target_beta_in_ci95"] else 0.0 for b in out_bins])) if out_bins else 0.0
    frac_eta = float(np.mean([1.0 if b["target_eta_in_ci95"] else 0.0 for b in out_bins])) if out_bins else 0.0
    frac_joint = (
        float(np.mean([1.0 if (b["target_beta_in_ci95"] and b["target_eta_in_ci95"]) else 0.0 for b in out_bins]))
        if out_bins
        else 0.0
    )

    beta_all = np.array([r["beta"] for r in rows], dtype=float)
    eta_all = np.array([r["eta"] for r in rows], dtype=float)
    rmse_all = np.array([r["rmse"] for r in rows], dtype=float)
    phase_all = np.array([r["phase_min"] for r in rows], dtype=float)

    global_est = {
        "n_rows_total": int(len(rows)),
        "beta_median": float(np.median(beta_all)) if len(beta_all) else float("nan"),
        "beta_ci95": [float(np.quantile(beta_all, 0.025)), float(np.quantile(beta_all, 0.975))] if len(beta_all) else [float("nan"), float("nan")],
        "eta_median": float(np.median(eta_all)) if len(eta_all) else float("nan"),
        "eta_ci95": [float(np.quantile(eta_all, 0.025)), float(np.quantile(eta_all, 0.975))] if len(eta_all) else [float("nan"), float("nan")],
        "rmse_median": float(np.median(rmse_all)) if len(rmse_all) else float("inf"),
        "phase_min_median": float(np.median(phase_all)) if len(phase_all) else 0.0,
    }

    flags = {
        "enough_pointwise_bins_ge_6": bool(n_bins >= 6),
        "global_beta_target_inside_ci95": bool(global_est["beta_ci95"][0] <= beta_target <= global_est["beta_ci95"][1]) if len(beta_all) else False,
        "global_eta_target_inside_ci95": bool(global_est["eta_ci95"][0] <= eta_target <= global_est["eta_ci95"][1]) if len(eta_all) else False,
        "bin_beta_target_coverage_ge_0p50": bool(frac_beta >= 0.50),
        "bin_eta_target_coverage_ge_0p50": bool(frac_eta >= 0.50),
        "bin_joint_target_coverage_ge_0p35": bool(frac_joint >= 0.35),
        "median_window_rmse_le_0p10": bool(global_est["rmse_median"] <= 0.10),
        "phase_condition_strength_ge_0p75": bool(global_est["phase_min_median"] >= 0.75),
    }

    pass_count = int(sum(1 for v in flags.values() if v))

    return {
        "n_bins": n_bins,
        "bins": out_bins,
        "coverage": {
            "beta": frac_beta,
            "eta": frac_eta,
            "joint": frac_joint,
        },
        "global": global_est,
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": int(len(flags)),
    }


def main() -> None:
    mod = load_qw1739_module()
    d2039 = json.loads((ROOT / "report_qw2039_derivation_compatible_refrozen_kernel_gate.json").read_text(encoding="utf-8"))
    beta_target = float(d2039["selected_kernel"]["beta"])
    eta_target = float(d2039["selected_kernel"]["eta"])

    lam_grid = [0.0, 0.4, 0.8, 1.2, 1.6]

    micro_profiles = []
    for n in [64, 96, 128]:
        for sidx in range(22):
            seed = 204700 + 1000 * n + sidx
            cfg = mod.MicroCfg(n_nodes=n, seed=seed)
            rng = np.random.default_rng(seed)
            theta, q = mod.build_micro_state(cfg, rng)
            w = mod.build_signed_dynamic_matrix(cfg, theta, q, rng)
            g = mod.effective_path_sum(w, lam=cfg.dyn_lambda, kmax=cfg.dyn_kmax)
            y = mod.profile_from_matrix(g, dmax=min(24, n // 2))
            est = mod.derive_omega_phi_beta(np.array(y, dtype=float))
            tau_idx = torsion_index_from_matrix(np.array(w, dtype=float))
            micro_profiles.append(
                {
                    "n": int(n),
                    "seed": int(seed),
                    "y": np.array(y, dtype=float),
                    "omega_hat": float(est["omega_hat"]),
                    "phi_hat": float(est["phi_hat"]),
                    "tau_idx": float(tau_idx),
                }
            )

    scan_rows = []
    for lam in lam_grid:
        all_rows = []
        for p in micro_profiles:
            rows = derive_rows_for_lambda(
                y=p["y"],
                omega_hat=float(p["omega_hat"]),
                phi_hat=float(p["phi_hat"]),
                tau_idx=float(p["tau_idx"]),
                lam_tau=float(lam),
            )
            all_rows.extend(rows)

        if len(all_rows) < 40:
            summary = {
                "lambda_tau": float(lam),
                "n_rows_total": int(len(all_rows)),
                "n_bins": 0,
                "pass_count": 0,
                "total_flags": 8,
                "flags": {},
                "coverage": {"beta": 0.0, "eta": 0.0, "joint": 0.0},
                "global": {
                    "n_rows_total": int(len(all_rows)),
                    "beta_median": float("nan"),
                    "beta_ci95": [float("nan"), float("nan")],
                    "eta_median": float("nan"),
                    "eta_ci95": [float("nan"), float("nan")],
                    "rmse_median": float("inf"),
                    "phase_min_median": 0.0,
                },
                "bins": [],
            }
        else:
            s = summarize(all_rows, beta_target=beta_target, eta_target=eta_target)
            summary = {"lambda_tau": float(lam), **s}

        scan_rows.append(summary)

    # Rank: pass_count -> n_bins -> phase strength -> rmse
    best = sorted(
        scan_rows,
        key=lambda r: (
            -int(r["pass_count"]),
            -int(r["n_bins"]),
            -float(r["global"]["phase_min_median"]),
            float(r["global"]["rmse_median"]),
        ),
    )[0]
    baseline = next(r for r in scan_rows if abs(float(r["lambda_tau"]) - 0.0) < 1e-12)

    improved = bool(int(best["pass_count"]) > int(baseline["pass_count"]))

    if int(best["pass_count"]) == int(best["total_flags"]):
        verdict = "SIGNED_PHASE_TORSION_OBSERVABLE_PASS"
        readiness = "POINTWISE_IDENTIFIABILITY_REPAIRED"
    elif int(best["pass_count"]) >= 6:
        verdict = "SIGNED_PHASE_TORSION_OBSERVABLE_PARTIAL"
        readiness = "PARTIAL_REPAIR_RETAINED"
    else:
        verdict = "SIGNED_PHASE_TORSION_OBSERVABLE_FAIL"
        readiness = "NO_SUFFICIENT_IDENTIFIABILITY_REPAIR"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "target_refrozen": {
            "beta": beta_target,
            "eta": eta_target,
        },
        "scan": scan_rows,
        "baseline_lambda0": baseline,
        "selected": best,
        "improved_vs_lambda0": improved,
        "verdict": verdict,
        "readiness": readiness,
        "required_next_step": (
            "PROMOTE_LAMBDA_TAU_SELECTED_INTO_MICRO_STAGEC_INTERSECTION_GATE"
            if improved
            else "REBUILD_MICRO_OBSERVABLE_BASIS_WITH_ADDITIONAL_TORSION_PHASE_CHANNEL"
        ),
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2047: SIGNED PHASE-TORSION OBSERVABLE SCAN",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- Readiness: **{readiness}**",
        f"- improved_vs_lambda0: {improved}",
        "",
        "## Selected lambda_tau",
        f"- lambda_tau: {best['lambda_tau']:.3f}",
        f"- pass_count: {best['pass_count']}/{best['total_flags']}",
        f"- n_bins: {best['n_bins']}",
        f"- coverage beta/eta/joint: {best['coverage']['beta']:.4f}/{best['coverage']['eta']:.4f}/{best['coverage']['joint']:.4f}",
        f"- phase_min_median: {best['global']['phase_min_median']:.6f}",
        f"- rmse_median: {best['global']['rmse_median']:.6f}",
        "",
        "## Baseline lambda_tau=0",
        f"- pass_count: {baseline['pass_count']}/{baseline['total_flags']}",
        f"- n_bins: {baseline['n_bins']}",
        f"- phase_min_median: {baseline['global']['phase_min_median']:.6f}",
        f"- rmse_median: {baseline['global']['rmse_median']:.6f}",
        "",
        "## Lambda Scan Summary",
    ]

    for r in scan_rows:
        lines.append(
            "- lambda={l:.2f} | pass={p}/{t} | bins={b} | cov(beta/eta/joint)={cb:.3f}/{ce:.3f}/{cj:.3f} | phase={ph:.3f} | rmse={rm:.4f}".format(
                l=float(r["lambda_tau"]),
                p=int(r["pass_count"]),
                t=int(r["total_flags"]),
                b=int(r["n_bins"]),
                cb=float(r["coverage"]["beta"]),
                ce=float(r["coverage"]["eta"]),
                cj=float(r["coverage"]["joint"]),
                ph=float(r["global"]["phase_min_median"]),
                rm=float(r["global"]["rmse_median"]),
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

    print(f"[QW-2047] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2047] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2047] verdict={verdict} selected_lambda={best['lambda_tau']} pass={best['pass_count']}/{best['total_flags']}")


if __name__ == "__main__":
    main()
