#!/usr/bin/env python3
"""
QW-1813: Sequence-window dynamic model (phase-2 entry).

Uses windowed cross-DFA trajectories per pulsar pair, then builds
sequence-level features:
  f_trend   - slope over window index
  f_var     - std over windows
  f_switch  - regime-switch proxy from sign changes

Model:
  M2S(theta, f) = A * HD(theta)^q + b1*f_trend + b2*f_var + b3*f_switch + C
"""

from __future__ import annotations

import importlib.util
import json
from itertools import combinations
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np
from scipy.special import logsumexp


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1813_sequence_window_model.json"
OUT_MD = ROOT / "RAPORT_QW1813_SEQUENCE_WINDOW_MODEL.md"


def load_helper():
    path = ROOT / "QW_1787_COVERAGE_FRACTION_SWEEP.py"
    spec = importlib.util.spec_from_file_location("qw1787_helper", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def center(v: np.ndarray) -> np.ndarray:
    return v - float(np.mean(v))


def window_trajectory(helper, x: np.ndarray, y: np.ndarray) -> np.ndarray:
    n = min(len(x), len(y))
    if n < 140:
        return np.array([], dtype=float)
    w = max(70, int(round(0.6 * n)))
    step = max(12, int(round(0.1 * n)))
    vals = []
    for start in range(0, n - w + 1, step):
        xx = x[start : start + w]
        yy = y[start : start + w]
        h = helper.cross_dfa(xx, yy, min_scale=10)
        if np.isfinite(h):
            vals.append(float(h))
    return np.array(vals, dtype=float)


def sequence_features(traj: np.ndarray) -> Dict[str, float] | None:
    if len(traj) < 4:
        return None
    t = np.arange(len(traj), dtype=float)
    slope = float(np.polyfit(t, traj, 1)[0])
    var = float(np.std(traj))
    signs = np.sign(traj)
    sign_changes = np.sum(signs[1:] != signs[:-1])
    switch = float(sign_changes / max(len(traj) - 1, 1))
    return {"f_trend": slope, "f_var": var, "f_switch": switch, "n_windows": int(len(traj))}


def loglike(H: np.ndarray, sigma: float, model: np.ndarray) -> float:
    z = (H - model) / sigma
    return float(-0.5 * np.sum(z * z))


def evidence_flat(H: np.ndarray, n_mc: int, seed: int) -> float:
    rng = np.random.default_rng(seed)
    C = rng.uniform(-1.0, 2.0, n_mc)
    sigma = max(float(np.std(H)), 1e-6)
    ll = np.array([loglike(H, sigma, np.full_like(H, c)) for c in C], dtype=float)
    return float(logsumexp(ll) - np.log(len(ll)))


def evidence_m2(helper, theta: np.ndarray, H: np.ndarray, q_center: float, q_width: float, n_mc: int, seed: int) -> float:
    rng = np.random.default_rng(seed)
    A = rng.uniform(-2.0, 2.0, n_mc)
    C = rng.uniform(-1.0, 2.0, n_mc)
    q = np.clip(rng.normal(loc=q_center, scale=q_width, size=n_mc), 0.8, 2.4)
    sigma = max(float(np.std(H)), 1e-6)
    hd0 = np.clip(helper.hellings_downs(theta), 1e-9, None)
    ll = np.array([loglike(H, sigma, a * (hd0 ** qq) + c) for a, c, qq in zip(A, C, q)], dtype=float)
    return float(logsumexp(ll) - np.log(len(ll)))


def evidence_m2s(helper, theta: np.ndarray, H: np.ndarray, F1: np.ndarray, F2: np.ndarray, F3: np.ndarray, q_center: float, q_width: float, n_mc: int, seed: int) -> float:
    rng = np.random.default_rng(seed)
    A = rng.uniform(-2.0, 2.0, n_mc)
    B1 = rng.uniform(-1.2, 1.2, n_mc)
    B2 = rng.uniform(-1.2, 1.2, n_mc)
    B3 = rng.uniform(-1.2, 1.2, n_mc)
    C = rng.uniform(-1.0, 2.0, n_mc)
    q = np.clip(rng.normal(loc=q_center, scale=q_width, size=n_mc), 0.8, 2.4)
    sigma = max(float(np.std(H)), 1e-6)
    hd0 = np.clip(helper.hellings_downs(theta), 1e-9, None)
    ll = np.array(
        [loglike(H, sigma, a * (hd0 ** qq) + b1 * F1 + b2 * F2 + b3 * F3 + c) for a, b1, b2, b3, c, qq in zip(A, B1, B2, B3, C, q)],
        dtype=float,
    )
    return float(logsumexp(ll) - np.log(len(ll)))


def main() -> None:
    helper = load_helper()
    d1773 = json.loads((ROOT / "report_qw1773_omega_suppressed_legacy_projection.json").read_text(encoding="utf-8"))
    d1793 = json.loads((ROOT / "report_qw1793_model_lockin_gate.json").read_text(encoding="utf-8"))

    q_center = float(d1773["projection"]["p"])
    q_width = float(d1793["operational_protocol"]["q_width"])
    frac = float(d1793["operational_protocol"]["fraction"])
    cohort = d1793["operational_protocol"]["cohort"]
    n_match_min = int(cohort["n_match_min"])
    stability_max = float(cohort["stability_max"])

    residuals = helper.load_residuals(ROOT / "nano15/residuals/NANOGrav15yr_PulsarTiming_v2.1.0/residuals", max_psr=34)
    positions = helper.load_positions(ROOT / "nano15/parfiles")

    rows: List[Dict[str, float]] = []
    psr_list = list(residuals.keys())
    for p1, p2 in combinations(psr_list, 2):
        sep = helper.angular_sep(p1, p2, positions)
        if sep is None:
            continue
        x, y = helper.match_epochs(residuals[p1], residuals[p2], tol_days=30.0)
        if x is None:
            continue
        hxy = helper.cross_dfa(x, y, min_scale=15)
        if not np.isfinite(hxy):
            continue
        stab = helper.split_half_stability(x, y)
        if len(x) < n_match_min or stab > stability_max:
            continue
        traj = window_trajectory(helper, x, y)
        feats = sequence_features(traj)
        if feats is None:
            continue
        rows.append(
            {
                "theta_deg": float(sep),
                "hxy": float(hxy),
                "f_trend": feats["f_trend"],
                "f_var": feats["f_var"],
                "f_switch": feats["f_switch"],
                "n_windows": feats["n_windows"],
            }
        )

    if len(rows) < 80:
        raise RuntimeError(f"Sequence cohort too small: {len(rows)}")

    theta_all = np.array([r["theta_deg"] for r in rows], dtype=float)
    H_all = np.array([r["hxy"] for r in rows], dtype=float)
    F1_all = center(np.array([r["f_trend"] for r in rows], dtype=float))
    F2_all = center(np.array([r["f_var"] for r in rows], dtype=float))
    F3_all = center(np.array([r["f_switch"] for r in rows], dtype=float))

    z0 = evidence_flat(H_all, n_mc=8500, seed=19201)
    z2 = evidence_m2(helper, theta_all, H_all, q_center=q_center, q_width=q_width, n_mc=13000, seed=19202)
    zs = evidence_m2s(helper, theta_all, H_all, F1_all, F2_all, F3_all, q_center=q_center, q_width=q_width, n_mc=19000, seed=19203)
    full_m2 = float(z2 - z0)
    full_m2s = float(zs - z0)
    full_delta = float(full_m2s - full_m2)

    rng = np.random.default_rng(19204)
    rep_rows = []
    for i in range(12):
        idx = helper.stratified_subsample_indices(theta_all, frac=frac, rng=rng, n_bins=8)
        if len(idx) < 75:
            continue
        th = theta_all[idx]
        hh = H_all[idx]
        f1 = F1_all[idx]
        f2 = F2_all[idx]
        f3 = F3_all[idx]
        b0 = evidence_flat(hh, n_mc=3900, seed=19210 + 20 * i + 1)
        b2 = evidence_m2(helper, th, hh, q_center=q_center, q_width=q_width, n_mc=6200, seed=19210 + 20 * i + 2)
        bs = evidence_m2s(helper, th, hh, f1, f2, f3, q_center=q_center, q_width=q_width, n_mc=8600, seed=19210 + 20 * i + 3)
        l2 = float(b2 - b0)
        ls = float(bs - b0)
        rep_rows.append({"rep": i, "logB_m2_vs_flat": l2, "logB_m2s_vs_flat": ls, "delta_m2s_vs_m2": float(ls - l2)})

    arr_d = np.array([r["delta_m2s_vs_m2"] for r in rep_rows], dtype=float)
    arr_s = np.array([r["logB_m2s_vs_flat"] for r in rep_rows], dtype=float)

    summary = {
        "n_pairs_sequence_cohort": len(rows),
        "mean_windows_per_pair": float(np.mean([r["n_windows"] for r in rows])),
        "fraction": frac,
        "q_width": q_width,
        "full_logB_m2_vs_flat": full_m2,
        "full_logB_m2s_vs_flat": full_m2s,
        "full_delta_m2s_vs_m2": full_delta,
        "replications": len(rep_rows),
        "prob_m2s_gt_m2": float(np.mean(arr_d > 0.0)),
        "prob_m2s_gt_flat": float(np.mean(arr_s > 0.0)),
        "median_delta_m2s_vs_m2": float(np.median(arr_d)),
        "std_delta_m2s_vs_m2": float(np.std(arr_d)),
    }

    pass_full = summary["full_delta_m2s_vs_m2"] > 0.0
    pass_rep = summary["prob_m2s_gt_m2"] >= 0.80 and summary["prob_m2s_gt_flat"] >= 0.95
    pass_disp = summary["std_delta_m2s_vs_m2"] <= 0.30

    if pass_full and pass_rep and pass_disp:
        verdict = "SEQUENCE_WINDOW_MODEL_SUPPORTED"
    elif pass_full and pass_rep:
        verdict = "SEQUENCE_WINDOW_MODEL_PARTIAL"
    else:
        verdict = "SEQUENCE_WINDOW_MODEL_WEAK"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "summary": summary,
        "pass_flags": {
            "full_gain": bool(pass_full),
            "replication_gain": bool(pass_rep),
            "dispersion_control": bool(pass_disp),
        },
        "verdict": verdict,
        "replications": rep_rows,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1813: SEQUENCE WINDOW MODEL",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Sequence cohort size: {len(rows)}",
        f"- Mean windows per pair: {summary['mean_windows_per_pair']:.2f}",
        f"- Full logB M2/M2S: {full_m2:.4f} / {full_m2s:.4f}",
        f"- Full delta M2S-M2: {full_delta:.4f}",
        f"- P(M2S>M2): {summary['prob_m2s_gt_m2']:.3f}",
        f"- P(M2S>flat): {summary['prob_m2s_gt_flat']:.3f}",
        f"- Std delta M2S-M2: {summary['std_delta_m2s_vs_m2']:.3f}",
        f"- Verdict: **{verdict}**",
        "",
        "## Pass Flags",
        f"- full_gain: {pass_full}",
        f"- replication_gain: {pass_rep}",
        f"- dispersion_control: {pass_disp}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1813] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1813] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
