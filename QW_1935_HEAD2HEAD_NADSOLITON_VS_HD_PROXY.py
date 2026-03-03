#!/usr/bin/env python3
"""
QW-1935: Post-closure head-to-head predictive audit vs SM+GR proxy.

Models:
- Legacy nadsoliton kernel (eta=1 baseline from QW-1932).
- Reparameterized nadsoliton kernel (selected eta branch from QW-1932).
- SM+GR proxy in PTA space: Hellings-Downs angular template (HD) + affine nuisance.

Note:
- This is an internal predictive comparison on currently available true-external package.
- It is not a final independent confirmatory publication-grade protocol.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1935_head2head_nadsoliton_vs_hd_proxy.json"
OUT_MD = ROOT / "RAPORT_QW1935_HEAD2HEAD_NADSOLITON_VS_HD_PROXY.md"


def split_index(pair_id: str, salt: int, k: int = 2) -> int:
    h = hashlib.sha256(f"{salt}|{pair_id}".encode("utf-8")).hexdigest()
    return int(h[-8:], 16) % k


def hd_curve(theta_deg: np.ndarray) -> np.ndarray:
    th = np.deg2rad(theta_deg)
    x = (1.0 - np.cos(th)) / 2.0
    x = np.clip(x, 1e-12, 1.0)
    return 0.5 + 1.5 * x * np.log(x) - 0.25 * x


def nad_kernel(theta_deg: np.ndarray, omega: float, phi: float, beta: float, eta: float) -> np.ndarray:
    tmin, tmax = float(np.min(theta_deg)), float(np.max(theta_deg))
    d_eff = 1.0 + 11.0 * (theta_deg - tmin) / max(tmax - tmin, 1e-12)
    return np.cos(omega * d_eff + phi) / (1.0 + beta * (d_eff ** eta))


def fit_affine_on_discovery(feature: np.ndarray, y: np.ndarray, disc: np.ndarray) -> np.ndarray:
    X = np.column_stack([feature[disc], np.ones(int(np.sum(disc)), dtype=float)])
    coef, *_ = np.linalg.lstsq(X, y[disc], rcond=None)
    a, b = float(coef[0]), float(coef[1])
    return a * feature + b


def metrics_holdout(y_true: np.ndarray, y_pred: np.ndarray, y_disc: np.ndarray) -> Dict[str, float]:
    corr = float(np.corrcoef(y_true, y_pred)[0, 1])
    rmse = float(np.sqrt(np.mean((y_true - y_pred) ** 2)))
    base = float(np.mean(y_disc))
    rmse_base = float(np.sqrt(np.mean((y_true - base) ** 2)))
    gain = float((rmse_base - rmse) / max(rmse_base, 1e-12))
    return {
        "pearson_corr": corr,
        "rmse": rmse,
        "rmse_base": rmse_base,
        "rmse_gain_ratio": gain,
    }


def sign_test_pvalue(win_count: int, n: int) -> float:
    # One-sided exact binomial tail under p=0.5.
    tail = 0.0
    for k in range(win_count, n + 1):
        tail += math.comb(n, k)
    return float(tail / (2**n))


def summarize_model(rows: List[Dict[str, object]], key: str) -> Dict[str, float]:
    corr = np.array([float(r[key]["pearson_corr"]) for r in rows], dtype=float)
    gain = np.array([float(r[key]["rmse_gain_ratio"]) for r in rows], dtype=float)
    rmse = np.array([float(r[key]["rmse"]) for r in rows], dtype=float)
    return {
        "corr_median": float(np.median(corr)),
        "gain_median": float(np.median(gain)),
        "rmse_median": float(np.median(rmse)),
    }


def summarize_pairwise(rows: List[Dict[str, object]], a: str, b: str) -> Dict[str, float]:
    # Positive delta_rmse means model a has lower RMSE than model b.
    d_rmse = np.array([float(r[b]["rmse"] - r[a]["rmse"]) for r in rows], dtype=float)
    d_corr = np.array([float(r[a]["pearson_corr"] - r[b]["pearson_corr"]) for r in rows], dtype=float)
    wins = int(np.sum(d_rmse > 0.0))
    n = int(len(d_rmse))
    return {
        "delta_rmse_median": float(np.median(d_rmse)),
        "delta_rmse_q10": float(np.quantile(d_rmse, 0.10)),
        "delta_rmse_q90": float(np.quantile(d_rmse, 0.90)),
        "delta_corr_median": float(np.median(d_corr)),
        "win_rate": float(wins / max(n, 1)),
        "win_count": wins,
        "n_salts": n,
        "sign_test_p_one_sided": sign_test_pvalue(wins, n) if n > 0 else 1.0,
    }


def eval_dataset(df: pd.DataFrame, triads: Dict[str, Dict[str, float]], salts: List[int]) -> Dict[str, object]:
    req = {"pair_id", "theta_deg", "hxy"}
    miss = sorted(req - set(df.columns))
    if miss:
        raise RuntimeError(f"Missing columns: {miss}")

    pid = df["pair_id"].astype(str).to_numpy()
    theta = df["theta_deg"].to_numpy(dtype=float)
    y = df["hxy"].to_numpy(dtype=float)

    feat_legacy = nad_kernel(
        theta,
        omega=float(triads["legacy"]["omega"]),
        phi=float(triads["legacy"]["phi"]),
        beta=float(triads["legacy"]["beta"]),
        eta=1.0,
    )
    feat_reparam = nad_kernel(
        theta,
        omega=float(triads["reparam"]["omega"]),
        phi=float(triads["reparam"]["phi"]),
        beta=float(triads["reparam"]["beta"]),
        eta=float(triads["reparam"]["eta"]),
    )
    feat_hd = hd_curve(theta)

    rows = []
    for salt in salts:
        fold = np.array([split_index(x, salt=salt, k=2) for x in pid], dtype=int)
        disc = fold == 0
        hold = fold == 1

        if int(np.sum(disc)) < 300 or int(np.sum(hold)) < 300:
            raise RuntimeError(f"Salt {salt}: too small split.")

        yhat_legacy = fit_affine_on_discovery(feat_legacy, y, disc)
        yhat_reparam = fit_affine_on_discovery(feat_reparam, y, disc)
        yhat_hd = fit_affine_on_discovery(feat_hd, y, disc)

        m_legacy = metrics_holdout(y[hold], yhat_legacy[hold], y[disc])
        m_reparam = metrics_holdout(y[hold], yhat_reparam[hold], y[disc])
        m_hd = metrics_holdout(y[hold], yhat_hd[hold], y[disc])

        rows.append(
            {
                "salt": int(salt),
                "n_discovery": int(np.sum(disc)),
                "n_holdout": int(np.sum(hold)),
                "legacy": m_legacy,
                "reparam": m_reparam,
                "hd_proxy": m_hd,
            }
        )

    summary = {
        "legacy": summarize_model(rows, key="legacy"),
        "reparam": summarize_model(rows, key="reparam"),
        "hd_proxy": summarize_model(rows, key="hd_proxy"),
        "reparam_vs_hd_proxy": summarize_pairwise(rows, a="reparam", b="hd_proxy"),
        "legacy_vs_hd_proxy": summarize_pairwise(rows, a="legacy", b="hd_proxy"),
        "reparam_vs_legacy": summarize_pairwise(rows, a="reparam", b="legacy"),
    }
    return {
        "rows": rows,
        "summary": summary,
    }


def load_qw1932_triads() -> Dict[str, Dict[str, float]]:
    d = json.loads((ROOT / "report_qw1932_physical_reparameterization_eta_scan.json").read_text(encoding="utf-8"))
    base = d.get("eta1_baseline", {}).get("fit", {})
    sel = d.get("selected", {})
    fit = sel.get("fit", {})
    return {
        "legacy": {
            "omega": float(base["omega"]),
            "phi": float(base["phi"]),
            "beta": float(base["beta"]),
        },
        "reparam": {
            "omega": float(fit["omega"]),
            "phi": float(fit["phi"]),
            "beta": float(fit["beta"]),
            "eta": float(sel["eta"]),
        },
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--primary-pta",
        type=str,
        default="external_confirmatory_v2/beta_channel_true_external/beta_channel_pairs.csv",
    )
    ap.add_argument(
        "--stress-pta",
        type=str,
        default="external_confirmatory_v2/confirmatory_dataset_external_source_alpha6_1831cfg/pta_v2_pairs.csv",
    )
    ap.add_argument("--n-salts", type=int, default=24)
    args = ap.parse_args()

    triads = load_qw1932_triads()
    p_primary = (ROOT / args.primary_pta).resolve()
    p_stress = (ROOT / args.stress_pta).resolve()
    if not p_primary.exists():
        raise RuntimeError(f"Primary PTA not found: {p_primary}")
    if not p_stress.exists():
        raise RuntimeError(f"Stress PTA not found: {p_stress}")

    salts = list(range(int(args.n_salts)))
    d_primary = pd.read_csv(p_primary)
    d_stress = pd.read_csv(p_stress)

    primary = eval_dataset(d_primary, triads=triads, salts=salts)
    stress = eval_dataset(d_stress, triads=triads, salts=salts)

    p_cmp = primary["summary"]["reparam_vs_hd_proxy"]
    s_cmp = stress["summary"]["reparam_vs_hd_proxy"]

    flags = {
        "primary_reparam_beats_hd_win_rate_ge_0p80": bool(p_cmp["win_rate"] >= 0.80),
        "primary_reparam_beats_hd_sign_p_le_0p05": bool(p_cmp["sign_test_p_one_sided"] <= 0.05),
        "stress_reparam_beats_hd_win_rate_ge_0p80": bool(s_cmp["win_rate"] >= 0.80),
        "stress_reparam_beats_hd_sign_p_le_0p05": bool(s_cmp["sign_test_p_one_sided"] <= 0.05),
    }

    if all(flags.values()):
        verdict = "HEAD2HEAD_REPARAM_OUTPERFORMS_HD_PROXY_BOTH_DATASETS"
        required_next = "FREEZE_PUBLIC_CONFIRMATORY_HEAD2HEAD_PROTOCOL"
    elif flags["primary_reparam_beats_hd_win_rate_ge_0p80"] and flags["primary_reparam_beats_hd_sign_p_le_0p05"]:
        verdict = "HEAD2HEAD_MIXED_PRIMARY_ONLY_REPARAM_BETTER"
        required_next = "CHARACTERIZE_REGIME_SPLIT_AND_DEFINE_DOMAIN_OF_VALIDITY_BEFORE_STRONG_TOE_CLAIM"
    else:
        verdict = "HEAD2HEAD_REPARAM_NOT_BETTER_THAN_HD_PROXY"
        required_next = "REWORK_PREDICTIVE_CHANNEL_BEFORE_SM_GR_SUPERIORITY_CLAIM"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "triads": triads,
        "config": {
            "n_salts": int(args.n_salts),
            "datasets": {
                "primary": str(p_primary),
                "stress": str(p_stress),
            },
        },
        "results": {
            "primary": primary,
            "stress": stress,
        },
        "flags": flags,
        "verdict": verdict,
        "required_next_step": required_next,
        "caveat": (
            "HD here is a pragmatic SM+GR PTA proxy with affine nuisance, "
            "not a full detector-noise Bayesian pipeline."
        ),
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    p_sm = primary["summary"]
    s_sm = stress["summary"]
    lines = [
        "# RAPORT QW-1935: HEAD-TO-HEAD NADSOLITON VS HD PROXY",
        "",
        f"- Data UTC: {out['generated_utc']}",
        (
            "- Reparam triad: "
            f"omega={triads['reparam']['omega']:.6f}, phi={triads['reparam']['phi']:.6f}, "
            f"beta={triads['reparam']['beta']:.6f}, eta={triads['reparam']['eta']:.2f}"
        ),
        f"- Verdict: **{verdict}**",
        "",
        "## Primary Summary",
        f"- reparam corr/gain/rmse medians: {p_sm['reparam']['corr_median']:.4f}/{p_sm['reparam']['gain_median']:.4f}/{p_sm['reparam']['rmse_median']:.4f}",
        f"- hd corr/gain/rmse medians: {p_sm['hd_proxy']['corr_median']:.4f}/{p_sm['hd_proxy']['gain_median']:.4f}/{p_sm['hd_proxy']['rmse_median']:.4f}",
        f"- reparam_vs_hd delta_rmse median: {p_sm['reparam_vs_hd_proxy']['delta_rmse_median']:.6f}",
        f"- reparam_vs_hd win_rate/sign_p: {p_sm['reparam_vs_hd_proxy']['win_rate']:.4f}/{p_sm['reparam_vs_hd_proxy']['sign_test_p_one_sided']:.4g}",
        "",
        "## Stress Summary",
        f"- reparam corr/gain/rmse medians: {s_sm['reparam']['corr_median']:.4f}/{s_sm['reparam']['gain_median']:.4f}/{s_sm['reparam']['rmse_median']:.4f}",
        f"- hd corr/gain/rmse medians: {s_sm['hd_proxy']['corr_median']:.4f}/{s_sm['hd_proxy']['gain_median']:.4f}/{s_sm['hd_proxy']['rmse_median']:.4f}",
        f"- reparam_vs_hd delta_rmse median: {s_sm['reparam_vs_hd_proxy']['delta_rmse_median']:.6f}",
        f"- reparam_vs_hd win_rate/sign_p: {s_sm['reparam_vs_hd_proxy']['win_rate']:.4f}/{s_sm['reparam_vs_hd_proxy']['sign_test_p_one_sided']:.4g}",
        "",
        "## Flags",
    ]
    for k, v in flags.items():
        lines.append(f"- {k}: {v}")
    lines.extend(
        [
            "",
            "## Caveat",
            f"- {out['caveat']}",
            "",
            "## Required Next Step",
            f"- {required_next}",
            "",
            "## Artifacts",
            f"- JSON: `{OUT_JSON.name}`",
        ]
    )
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1935] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1935] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1935] verdict={verdict}")


if __name__ == "__main__":
    main()

