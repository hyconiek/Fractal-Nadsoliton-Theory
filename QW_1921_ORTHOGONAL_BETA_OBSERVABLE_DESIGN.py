#!/usr/bin/env python3
"""
QW-1921: Orthogonal beta observable design on blind external data.

Purpose:
- execute next step from QW-1920:
  DESIGN_ORTHOGONAL_BETA_OBSERVABLE_ON_BLIND_EXTERNAL_WITH_INTERVENTION_PROTOCOL

Protocol:
1) deterministic discovery/holdout split by pair_id hash,
2) define omega-proxy from oscillatory residual structure,
3) evaluate beta-candidate observables for:
   - low coupling to omega-proxy,
   - high response to beta-like damping intervention,
   - low leakage under omega-like phase intervention,
4) select best candidate on discovery, validate on holdout.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1921_orthogonal_beta_observable_design.json"
OUT_MD = ROOT / "RAPORT_QW1921_ORTHOGONAL_BETA_OBSERVABLE_DESIGN.md"


def split_index(pair_id: str, k: int = 2) -> int:
    h = hashlib.sha256(pair_id.encode("utf-8")).hexdigest()
    return int(h[-8:], 16) % k


def moving_average(x: np.ndarray, win: int) -> np.ndarray:
    w = max(3, int(win) | 1)
    pad = w // 2
    xp = np.pad(x, (pad, pad), mode="edge")
    ker = np.ones(w, dtype=float) / float(w)
    y = np.convolve(xp, ker, mode="same")
    return y[pad : pad + len(x)]


def robust_scale(x: np.ndarray) -> float:
    med = float(np.median(x))
    mad = float(np.median(np.abs(x - med)))
    return max(1e-9, 1.4826 * mad)


def omega_proxy(theta: np.ndarray, hxy: np.ndarray) -> np.ndarray:
    order = np.argsort(theta)
    y = hxy[order]
    trend = moving_average(y, win=41)
    resid = y - trend

    d1 = np.diff(resid, n=1, prepend=resid[0])
    d2 = np.diff(resid, n=2, prepend=[resid[0], resid[1] if len(resid) > 1 else resid[0]])
    sign_flip = (np.sign(d1[1:]) * np.sign(d1[:-1]) < 0.0).astype(float)
    sign_flip = np.concatenate([[0.0], sign_flip])

    proxy = 0.65 * np.abs(d2) + 0.35 * sign_flip

    out = np.zeros_like(proxy)
    out[order] = proxy
    return out


def recompute_with_hxy(df: pd.DataFrame, hxy_new: np.ndarray) -> pd.DataFrame:
    out = df.copy()
    out["hxy"] = hxy_new

    theta = out["theta_deg"].to_numpy(dtype=float)
    order = np.argsort(theta)
    y = hxy_new[order]

    trend = moving_average(y, win=41)
    resid = y - trend
    d1 = np.diff(resid, n=1, prepend=resid[0])
    d2 = np.diff(resid, n=2, prepend=[resid[0], resid[1] if len(resid) > 1 else resid[0]])

    abs_resid = np.abs(resid)
    abs_d1 = np.abs(d1)
    abs_d2 = np.abs(d2)

    inv = np.zeros_like(y)
    inv[order] = np.arange(len(y))
    inv = inv.astype(int)

    out["c_abs_resid"] = abs_resid[inv]
    out["c_abs_d1"] = abs_d1[inv]
    out["c_abs_d2"] = abs_d2[inv]
    out["c_resid_std_local"] = moving_average(abs_resid, win=31)[inv]

    return out


def candidate_map(df: pd.DataFrame) -> Dict[str, np.ndarray]:
    theta = df["theta_deg"].to_numpy(dtype=float)
    t = (theta - float(np.min(theta))) / max(float(np.max(theta) - np.min(theta)), 1e-12)

    f_std = df["f_std"].to_numpy(dtype=float)
    f_autoc1 = np.clip(df["f_autoc1"].to_numpy(dtype=float), -1.0, 1.0)
    f_switch = df["f_switch"].to_numpy(dtype=float)
    f_slope = np.abs(df["f_slope"].to_numpy(dtype=float))

    c_abs_resid = df["c_abs_resid"].to_numpy(dtype=float)
    c_abs_d1 = df["c_abs_d1"].to_numpy(dtype=float)
    c_abs_d2 = df["c_abs_d2"].to_numpy(dtype=float)
    c_resid_std_local = df["c_resid_std_local"].to_numpy(dtype=float)

    cands = {
        "B1_tail_weighted_resid": c_abs_resid * (0.2 + t),
        "B2_tail_weighted_d2": c_abs_d2 * (0.2 + t),
        "B3_std_feature": f_std,
        "B4_one_minus_autoc1": 1.0 - f_autoc1,
        "B5_switch_rate": f_switch,
        "B6_abs_slope": f_slope,
        "B7_local_resid_std": c_resid_std_local,
        "B8_hybrid": 0.5 * c_abs_d1 + 0.35 * c_abs_d2 + 0.15 * (1.0 - f_autoc1),
    }
    return cands


def median_shift_effect(x0: np.ndarray, x1: np.ndarray) -> float:
    s = robust_scale(x0)
    return float(abs(np.median(x1) - np.median(x0)) / s)


def evaluate_candidates(df: pd.DataFrame) -> Tuple[Dict[str, Dict[str, float]], str]:
    theta = df["theta_deg"].to_numpy(dtype=float)
    hxy = df["hxy"].to_numpy(dtype=float)

    op = omega_proxy(theta, hxy)

    # Interventions (designed, deterministic, no threshold fitting):
    # beta-like damping: stronger attenuation toward large theta.
    t = (theta - float(np.min(theta))) / max(float(np.max(theta) - np.min(theta)), 1e-12)
    hxy_beta = np.clip(hxy / (1.0 + 0.80 * t), 0.0, 1.0)

    # omega-like phase intervention: cyclic shift in theta-ordered residual component.
    order = np.argsort(theta)
    y = hxy[order]
    trend = moving_average(y, win=41)
    resid = y - trend
    shift = max(5, int(0.09 * len(resid)))
    resid_shift = np.roll(resid, shift)
    hxy_omega_ordered = np.clip(trend + resid_shift, 0.0, 1.0)
    hxy_omega = np.zeros_like(hxy_omega_ordered)
    hxy_omega[order] = hxy_omega_ordered

    base = recompute_with_hxy(df, hxy)
    beta_i = recompute_with_hxy(df, hxy_beta)
    omega_i = recompute_with_hxy(df, hxy_omega)

    c0 = candidate_map(base)
    cb = candidate_map(beta_i)
    co = candidate_map(omega_i)

    out: Dict[str, Dict[str, float]] = {}
    for k in sorted(c0.keys()):
        x0 = c0[k]
        xb = cb[k]
        xo = co[k]

        if np.std(x0) <= 1e-12 or np.std(op) <= 1e-12:
            corr_abs = 0.0
        else:
            corr_abs = float(abs(np.corrcoef(x0, op)[0, 1]))

        beta_sens = median_shift_effect(x0, xb)
        omega_leak = median_shift_effect(x0, xo)

        score = float(beta_sens - omega_leak - corr_abs)

        out[k] = {
            "corr_abs_with_omega_proxy": corr_abs,
            "beta_sensitivity": beta_sens,
            "omega_leakage": omega_leak,
            "orthogonal_score": score,
        }

    best = max(out.items(), key=lambda kv: kv[1]["orthogonal_score"])[0]
    return out, best


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--pta-path",
        type=str,
        default="external_confirmatory_v2/confirmatory_dataset_external_source_rebuild_v2_1831cfg/pta_v2_pairs.csv",
    )
    args = ap.parse_args()

    pta_path = (ROOT / args.pta_path).resolve()
    if not pta_path.exists():
        raise RuntimeError(f"PTA file not found: {pta_path}")

    df = pd.read_csv(pta_path)
    req = ["pair_id", "theta_deg", "hxy", "f_std", "f_autoc1", "f_switch", "f_slope"]
    miss = [c for c in req if c not in df.columns]
    if miss:
        raise RuntimeError(f"Missing required columns: {miss}")

    split = np.array([split_index(str(x), k=2) for x in df["pair_id"].astype(str)], dtype=int)
    disc = df.loc[split == 0].reset_index(drop=True)
    hold = df.loc[split == 1].reset_index(drop=True)

    if len(disc) < 300 or len(hold) < 300:
        raise RuntimeError("Discovery or holdout split too small.")

    disc_eval, selected = evaluate_candidates(disc)
    hold_eval, _ = evaluate_candidates(hold)

    dsel = disc_eval[selected]
    hsel = hold_eval[selected]

    flags = {
        "holdout_corr_abs_le_0p20": bool(hsel["corr_abs_with_omega_proxy"] <= 0.20),
        "holdout_beta_sensitivity_ge_0p35": bool(hsel["beta_sensitivity"] >= 0.35),
        "holdout_omega_leakage_le_0p25": bool(hsel["omega_leakage"] <= 0.25),
        "holdout_orthogonal_score_positive": bool(hsel["orthogonal_score"] > 0.0),
    }

    if all(flags.values()):
        verdict = "ORTHOGONAL_BETA_OBSERVABLE_DESIGN_PASS"
    elif sum(1 for v in flags.values() if v) >= 3:
        verdict = "ORTHOGONAL_BETA_OBSERVABLE_DESIGN_PARTIAL"
    else:
        verdict = "ORTHOGONAL_BETA_OBSERVABLE_DESIGN_WEAK"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "pta_path": str(pta_path),
        "split": {
            "method": "sha256(pair_id) parity",
            "n_discovery": int(len(disc)),
            "n_holdout": int(len(hold)),
        },
        "discovery": {
            "candidate_metrics": disc_eval,
            "selected_candidate": selected,
            "selected_metrics": dsel,
        },
        "holdout": {
            "candidate_metrics": hold_eval,
            "selected_candidate_metrics": hsel,
        },
        "pass_flags": flags,
        "verdict": verdict,
        "required_next_step": (
            "IMPLEMENT_SELECTED_BETA_OBSERVABLE_IN_BLIND_EXTERNAL_INTERVENTION_RUN"
            if verdict != "ORTHOGONAL_BETA_OBSERVABLE_DESIGN_WEAK"
            else "REVISE_CANDIDATE_LIBRARY_AND_REDESIGN_INTERVENTION_PROTOCOL"
        ),
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1921: ORTHOGONAL BETA OBSERVABLE DESIGN",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- Selected candidate (discovery): `{selected}`",
        "",
        "## Split",
        f"- n_discovery: {len(disc)}",
        f"- n_holdout: {len(hold)}",
        "",
        "## Selected Candidate Metrics",
        f"- discovery corr_abs_with_omega_proxy: {dsel['corr_abs_with_omega_proxy']:.4f}",
        f"- discovery beta_sensitivity: {dsel['beta_sensitivity']:.4f}",
        f"- discovery omega_leakage: {dsel['omega_leakage']:.4f}",
        f"- discovery orthogonal_score: {dsel['orthogonal_score']:.4f}",
        f"- holdout corr_abs_with_omega_proxy: {hsel['corr_abs_with_omega_proxy']:.4f}",
        f"- holdout beta_sensitivity: {hsel['beta_sensitivity']:.4f}",
        f"- holdout omega_leakage: {hsel['omega_leakage']:.4f}",
        f"- holdout orthogonal_score: {hsel['orthogonal_score']:.4f}",
        "",
        "## Pass Flags",
        f"- holdout_corr_abs_le_0p20: {flags['holdout_corr_abs_le_0p20']}",
        f"- holdout_beta_sensitivity_ge_0p35: {flags['holdout_beta_sensitivity_ge_0p35']}",
        f"- holdout_omega_leakage_le_0p25: {flags['holdout_omega_leakage_le_0p25']}",
        f"- holdout_orthogonal_score_positive: {flags['holdout_orthogonal_score_positive']}",
        "",
        "## Required Next Step",
        f"- {out['required_next_step']}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1921] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1921] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1921] verdict={verdict} selected={selected}")


if __name__ == "__main__":
    main()
