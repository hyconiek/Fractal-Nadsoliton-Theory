#!/usr/bin/env python3
"""
QW-2035: Multifield stability audit of external eta-kernel signal.

Purpose:
- test stability of blind external signal across multiple deterministic folds,
- quantify p-value robustness and practical effect consistency.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Tuple

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2035_eta_primary_signal_multifold_stability.json"
OUT_MD = ROOT / "RAPORT_QW2035_ETA_PRIMARY_SIGNAL_MULTIFOLD_STABILITY.md"


def assign_fold(pair_id: str, salt: str, k_folds: int) -> int:
    h = hashlib.sha256((salt + "|" + pair_id).encode("utf-8")).hexdigest()
    return int(h[-8:], 16) % k_folds


def safe_corr(a: np.ndarray, b: np.ndarray) -> float:
    sa = float(np.std(a))
    sb = float(np.std(b))
    if sa <= 1e-12 or sb <= 1e-12:
        return 0.0
    return float(np.corrcoef(a, b)[0, 1])


def kernel_eta(d_eff: np.ndarray, omega: float, phi: float, beta: float, eta: float) -> np.ndarray:
    return np.cos(omega * d_eff + phi) / (1.0 + beta * (d_eff**eta))


def eval_fold(
    pair_id: np.ndarray,
    theta: np.ndarray,
    y: np.ndarray,
    fold_id: int,
    k_folds: int,
    salt: str,
    kernel: Dict[str, float],
    n_perm: int,
    rng_seed: int,
) -> Dict[str, float]:
    fold = np.array([assign_fold(str(x), salt=salt, k_folds=k_folds) for x in pair_id], dtype=int)
    hold = fold == fold_id
    disc = ~hold
    if int(np.sum(hold)) < 120 or int(np.sum(disc)) < 300:
        raise RuntimeError("Fold too small for stable audit.")

    tmin = float(np.min(theta))
    tmax = float(np.max(theta))
    d_eff = 1.0 + 11.0 * (theta - tmin) / max(tmax - tmin, 1e-12)
    k = kernel_eta(d_eff, kernel["omega"], kernel["phi"], kernel["beta"], kernel["eta"])

    x = np.column_stack([k[disc], np.ones(int(np.sum(disc)), dtype=float)])
    coef, *_ = np.linalg.lstsq(x, y[disc], rcond=None)
    a_hat, b_hat = float(coef[0]), float(coef[1])
    yhat = a_hat * k + b_hat

    y_hold = y[hold]
    yh_hold = yhat[hold]
    corr = safe_corr(y_hold, yh_hold)
    rmse = float(np.sqrt(np.mean((y_hold - yh_hold) ** 2)))
    base = float(np.mean(y[disc]))
    rmse_base = float(np.sqrt(np.mean((y_hold - base) ** 2)))
    gain = float((rmse_base - rmse) / max(rmse_base, 1e-12))

    rng = np.random.default_rng(rng_seed)
    corr_null = np.empty(n_perm, dtype=float)
    gain_null = np.empty(n_perm, dtype=float)
    for i in range(n_perm):
        yp = rng.permutation(y_hold)
        corr_null[i] = safe_corr(yp, yh_hold)
        rmse_p = float(np.sqrt(np.mean((yp - yh_hold) ** 2)))
        rmse_b = float(np.sqrt(np.mean((yp - base) ** 2)))
        gain_null[i] = float((rmse_b - rmse_p) / max(rmse_b, 1e-12))

    q95_corr = float(np.quantile(corr_null, 0.95))
    q95_gain = float(np.quantile(gain_null, 0.95))
    p_corr = float((1 + np.sum(corr_null >= corr)) / (n_perm + 1))
    p_gain = float((1 + np.sum(gain_null >= gain)) / (n_perm + 1))

    return {
        "n_holdout": int(np.sum(hold)),
        "n_discovery": int(np.sum(disc)),
        "corr": corr,
        "gain": gain,
        "q95_corr_null": q95_corr,
        "q95_gain_null": q95_gain,
        "p_corr": p_corr,
        "p_gain": p_gain,
    }


def eval_dataset(path: Path, kernel: Dict[str, float], salt: str, k_folds: int, n_perm: int, seed0: int) -> Dict:
    df = pd.read_csv(path)
    req = {"pair_id", "theta_deg", "hxy"}
    if not req.issubset(set(df.columns)):
        missing = sorted(req - set(df.columns))
        raise RuntimeError(f"Missing required columns in {path}: {missing}")

    pair_id = df["pair_id"].astype(str).to_numpy()
    theta = df["theta_deg"].to_numpy(dtype=float)
    y = df["hxy"].to_numpy(dtype=float)

    folds = []
    for f in range(k_folds):
        folds.append(
            eval_fold(
                pair_id=pair_id,
                theta=theta,
                y=y,
                fold_id=f,
                k_folds=k_folds,
                salt=salt,
                kernel=kernel,
                n_perm=n_perm,
                rng_seed=seed0 + 17 * f,
            )
        )

    corr = np.array([x["corr"] for x in folds], dtype=float)
    gain = np.array([x["gain"] for x in folds], dtype=float)
    p_corr = np.array([x["p_corr"] for x in folds], dtype=float)
    p_gain = np.array([x["p_gain"] for x in folds], dtype=float)
    q95_corr = np.array([x["q95_corr_null"] for x in folds], dtype=float)

    summary = {
        "n_rows": int(len(df)),
        "k_folds": int(k_folds),
        "median_corr": float(np.median(corr)),
        "median_gain": float(np.median(gain)),
        "min_corr": float(np.min(corr)),
        "min_gain": float(np.min(gain)),
        "median_p_corr": float(np.median(p_corr)),
        "median_p_gain": float(np.median(p_gain)),
        "frac_p_corr_le_0p05": float(np.mean(p_corr <= 0.05)),
        "frac_p_corr_le_0p01": float(np.mean(p_corr <= 0.01)),
        "frac_gain_positive": float(np.mean(gain > 0.0)),
        "median_q95_corr_null": float(np.median(q95_corr)),
    }

    flags = {
        "median_p_corr_le_0p01": bool(summary["median_p_corr"] <= 0.01),
        "frac_p_corr_le_0p05_ge_0p80": bool(summary["frac_p_corr_le_0p05"] >= 0.80),
        "frac_gain_positive_ge_0p80": bool(summary["frac_gain_positive"] >= 0.80),
        "median_corr_gt_median_q95_null": bool(summary["median_corr"] > summary["median_q95_corr_null"]),
    }

    return {
        "path": str(path),
        "folds": folds,
        "summary": summary,
        "flags": flags,
        "pass_count": int(sum(1 for v in flags.values() if v)),
        "total_flags": int(len(flags)),
    }


def load_kernel() -> Dict[str, float]:
    d = json.loads((ROOT / "report_qw2030_final_stage_c_gate_combined_branch.json").read_text(encoding="utf-8"))
    k = d["kernel"]
    return {"omega": float(k["omega"]), "phi": float(k["phi"]), "beta": float(k["beta"]), "eta": float(k["eta"])}


def main() -> None:
    kernel = load_kernel()
    k_folds = 7
    n_perm = 2000

    primary = eval_dataset(
        path=ROOT / "external_confirmatory_v2/beta_channel_true_external_v2/beta_channel_pairs.csv",
        kernel=kernel,
        salt="QW2035_PRIMARY",
        k_folds=k_folds,
        n_perm=n_perm,
        seed0=203510,
    )
    stress = eval_dataset(
        path=ROOT / "external_confirmatory_v2/confirmatory_dataset_external_source_alpha6_1831cfg/pta_v2_pairs.csv",
        kernel=kernel,
        salt="QW2035_STRESS",
        k_folds=k_folds,
        n_perm=n_perm,
        seed0=203570,
    )

    global_flags = {
        "primary_pass_ge_3of4": bool(primary["pass_count"] >= 3),
        "stress_pass_ge_3of4": bool(stress["pass_count"] >= 3),
    }
    pass_count = int(sum(1 for v in global_flags.values() if v))
    total_flags = int(len(global_flags))

    if pass_count == total_flags:
        verdict = "ETA_SIGNAL_MULTIFOLD_STABILITY_PASS"
        readiness = "SIGNAL_STABILITY_ACCEPTABLE"
    elif pass_count == 1:
        verdict = "ETA_SIGNAL_MULTIFOLD_STABILITY_PARTIAL"
        readiness = "SIGNAL_STABILITY_PARTIAL"
    else:
        verdict = "ETA_SIGNAL_MULTIFOLD_STABILITY_FAIL"
        readiness = "SIGNAL_STABILITY_NOT_ACCEPTABLE"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "kernel_source": "report_qw2030_final_stage_c_gate_combined_branch.json:kernel",
        "kernel": kernel,
        "config": {"k_folds": k_folds, "n_perm_per_fold": n_perm},
        "datasets": {"primary": primary, "stress": stress},
        "global_flags": global_flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "readiness": readiness,
        "verdict": verdict,
        "required_next_step": (
            "INCLUDE_MULTIFOLD_STABILITY_IN_RELEASE_CANDIDATE"
            if verdict == "ETA_SIGNAL_MULTIFOLD_STABILITY_PASS"
            else "REASSESS_EXTERNAL_SIGNAL_GENERALIZATION"
        ),
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    p = primary["summary"]
    s = stress["summary"]
    lines = [
        "# RAPORT QW-2035: ETA PRIMARY SIGNAL MULTIFOLD STABILITY",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Readiness: **{readiness}**",
        f"- Verdict: **{verdict}**",
        f"- global pass_count: {pass_count}/{total_flags}",
        "",
        "## Primary Summary",
        (
            f"- median corr/gain: {p['median_corr']:.6f}/{p['median_gain']:.6f}, "
            f"median p_corr/p_gain: {p['median_p_corr']:.6f}/{p['median_p_gain']:.6f}"
        ),
        (
            f"- frac p_corr<=0.05 / <=0.01: {p['frac_p_corr_le_0p05']:.3f} / "
            f"{p['frac_p_corr_le_0p01']:.3f}"
        ),
        f"- frac gain>0: {p['frac_gain_positive']:.3f}",
        "",
        "## Stress Summary",
        (
            f"- median corr/gain: {s['median_corr']:.6f}/{s['median_gain']:.6f}, "
            f"median p_corr/p_gain: {s['median_p_corr']:.6f}/{s['median_p_gain']:.6f}"
        ),
        (
            f"- frac p_corr<=0.05 / <=0.01: {s['frac_p_corr_le_0p05']:.3f} / "
            f"{s['frac_p_corr_le_0p01']:.3f}"
        ),
        f"- frac gain>0: {s['frac_gain_positive']:.3f}",
        "",
        "## Global Flags",
    ]
    for k, v in global_flags.items():
        lines.append(f"- {k}: {v}")
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

    print(f"[QW-2035] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2035] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2035] readiness={readiness} verdict={verdict} pass={pass_count}/{total_flags}")


if __name__ == "__main__":
    main()
