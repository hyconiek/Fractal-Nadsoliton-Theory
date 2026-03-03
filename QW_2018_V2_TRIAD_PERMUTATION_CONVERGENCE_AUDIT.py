#!/usr/bin/env python3
"""
QW-2018: Permutation convergence audit for v2 triad blind validation.

Goal:
- quantify stability of holdout p-values around strict alpha=0.01,
- avoid one-run threshold artifacts.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2018_v2_triad_permutation_convergence_audit.json"
OUT_MD = ROOT / "RAPORT_QW2018_V2_TRIAD_PERMUTATION_CONVERGENCE_AUDIT.md"


def split_index(pair_id: str, k: int = 2) -> int:
    h = hashlib.sha256(pair_id.encode("utf-8")).hexdigest()
    return int(h[-8:], 16) % k


def rankdata(a: np.ndarray) -> np.ndarray:
    order = np.argsort(a)
    ranks = np.empty_like(order, dtype=float)
    n = len(a)
    i = 0
    while i < n:
        j = i
        while j + 1 < n and a[order[j + 1]] == a[order[i]]:
            j += 1
        r = 0.5 * (i + j) + 1.0
        ranks[order[i : j + 1]] = r
        i = j + 1
    return ranks


def spearmanr(a: np.ndarray, b: np.ndarray) -> float:
    ra = rankdata(a)
    rb = rankdata(b)
    if np.std(ra) <= 1e-12 or np.std(rb) <= 1e-12:
        return 0.0
    return float(np.corrcoef(ra, rb)[0, 1])


def build_holdout_vectors(df: pd.DataFrame, omega: float, phi: float, beta: float) -> Tuple[np.ndarray, np.ndarray, Dict[str, float]]:
    pair_id = df["pair_id"].astype(str).to_numpy()
    theta = df["theta_deg"].to_numpy(dtype=float)
    y = df["hxy"].to_numpy(dtype=float)

    tmin, tmax = float(np.min(theta)), float(np.max(theta))
    d_eff = 1.0 + 11.0 * (theta - tmin) / max(tmax - tmin, 1e-12)

    k = np.cos(omega * d_eff + phi) / (1.0 + beta * d_eff)

    fold = np.array([split_index(x, k=2) for x in pair_id], dtype=int)
    disc = fold == 0
    hold = fold == 1

    X = np.column_stack([k[disc], np.ones(int(np.sum(disc)), dtype=float)])
    coef, *_ = np.linalg.lstsq(X, y[disc], rcond=None)
    a_hat, b_hat = float(coef[0]), float(coef[1])

    yhat = a_hat * k + b_hat
    y_hold = y[hold]
    yhat_hold = yhat[hold]

    corr_hold = float(np.corrcoef(y_hold, yhat_hold)[0, 1])
    spearman_hold = float(spearmanr(y_hold, yhat_hold))
    rmse_hold = float(np.sqrt(np.mean((y_hold - yhat_hold) ** 2)))
    base = float(np.mean(y[disc]))
    rmse_base = float(np.sqrt(np.mean((y_hold - base) ** 2)))
    rmse_gain = float((rmse_base - rmse_hold) / max(rmse_base, 1e-12))

    metrics = {
        "pearson_corr": corr_hold,
        "spearman_corr": spearman_hold,
        "rmse": rmse_hold,
        "rmse_base": rmse_base,
        "rmse_gain_ratio": rmse_gain,
    }
    return y_hold, yhat_hold, metrics


def permutation_pvals(y_hold: np.ndarray, yhat_hold: np.ndarray, rmse_gain_ref: float, n_perm: int, seed: int) -> Dict[str, float]:
    rng = np.random.default_rng(seed)

    corr_hold = float(np.corrcoef(y_hold, yhat_hold)[0, 1])
    corr_null = np.empty(n_perm, dtype=float)
    gain_null = np.empty(n_perm, dtype=float)

    base = float(np.mean(y_hold))
    for i in range(n_perm):
        yp = rng.permutation(y_hold)
        corr_null[i] = float(np.corrcoef(yp, yhat_hold)[0, 1])
        rmse_p = float(np.sqrt(np.mean((yp - yhat_hold) ** 2)))
        rmse_b = float(np.sqrt(np.mean((yp - base) ** 2)))
        gain_null[i] = float((rmse_b - rmse_p) / max(rmse_b, 1e-12))

    p_corr = float((1 + np.sum(corr_null >= corr_hold)) / (n_perm + 1))
    p_gain = float((1 + np.sum(gain_null >= rmse_gain_ref)) / (n_perm + 1))
    q95_corr = float(np.quantile(corr_null, 0.95))
    q95_gain = float(np.quantile(gain_null, 0.95))

    return {
        "p_corr": p_corr,
        "p_rmse_gain": p_gain,
        "q95_corr": q95_corr,
        "q95_rmse_gain": q95_gain,
    }


def main() -> None:
    d1917 = json.loads((ROOT / "report_qw1917_triad_derivation_no_ansatz_profile.json").read_text(encoding="utf-8"))
    triad = d1917.get("optimum", {})
    omega = float(triad.get("omega"))
    phi = float(triad.get("phi"))
    beta = float(triad.get("beta"))

    p_primary = ROOT / "external_confirmatory_v2" / "beta_channel_true_external_v2" / "beta_channel_pairs.csv"
    if not p_primary.exists():
        raise RuntimeError(f"Dataset not found: {p_primary}")

    df = pd.read_csv(p_primary)
    y_hold, yhat_hold, base_metrics = build_holdout_vectors(df, omega=omega, phi=phi, beta=beta)

    n_perm_grid = [1000, 3000, 5000, 10000, 20000]
    n_seeds = 12

    rows: List[Dict[str, object]] = []
    for n_perm in n_perm_grid:
        p_corr_vals = []
        p_gain_vals = []
        for s in range(n_seeds):
            res = permutation_pvals(
                y_hold=y_hold,
                yhat_hold=yhat_hold,
                rmse_gain_ref=float(base_metrics["rmse_gain_ratio"]),
                n_perm=int(n_perm),
                seed=201800 + 1000 * n_perm + s,
            )
            p_corr_vals.append(float(res["p_corr"]))
            p_gain_vals.append(float(res["p_rmse_gain"]))

        p_corr_arr = np.array(p_corr_vals, dtype=float)
        p_gain_arr = np.array(p_gain_vals, dtype=float)

        row = {
            "n_perm": int(n_perm),
            "n_seeds": int(n_seeds),
            "p_corr_median": float(np.median(p_corr_arr)),
            "p_corr_q10": float(np.quantile(p_corr_arr, 0.10)),
            "p_corr_q90": float(np.quantile(p_corr_arr, 0.90)),
            "p_gain_median": float(np.median(p_gain_arr)),
            "p_gain_q10": float(np.quantile(p_gain_arr, 0.10)),
            "p_gain_q90": float(np.quantile(p_gain_arr, 0.90)),
            "frac_corr_le_0p01": float(np.mean(p_corr_arr <= 0.01)),
            "frac_gain_le_0p01": float(np.mean(p_gain_arr <= 0.01)),
        }
        rows.append(row)

    last = rows[-1]
    robust_corr = bool(last["frac_corr_le_0p01"] >= 0.90)
    robust_gain = bool(last["frac_gain_le_0p01"] >= 0.90)

    if robust_corr and robust_gain:
        verdict = "V2_TRIAD_PERMUTATION_CONVERGENCE_ROBUST_PASS"
        required_next = "PROCEED_TO_TRIAD_STAGE_WITH_V2_PACKAGE"
    elif (last["frac_corr_le_0p01"] >= 0.50) and (last["frac_gain_le_0p01"] >= 0.50):
        verdict = "V2_TRIAD_PERMUTATION_CONVERGENCE_BORDERLINE"
        required_next = "INCREASE_HOLDOUT_ROWS_AND_RECHECK"
    else:
        verdict = "V2_TRIAD_PERMUTATION_CONVERGENCE_FAIL"
        required_next = "REWORK_TRIAD_MAP_OR_COLLECT_MORE_EXTERNAL_ROWS"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "triad": {"omega": omega, "phi": phi, "beta": beta},
        "dataset_path": str(p_primary),
        "holdout_base_metrics": base_metrics,
        "grid_results": rows,
        "robust_flags": {
            "robust_corr_p_le_0p01": robust_corr,
            "robust_gain_p_le_0p01": robust_gain,
        },
        "verdict": verdict,
        "required_next_step": required_next,
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2018: V2 TRIAD PERMUTATION CONVERGENCE AUDIT",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- Required next step: `{required_next}`",
        "",
        "## Holdout Base Metrics",
        f"- pearson_corr: {base_metrics['pearson_corr']:.6f}",
        f"- spearman_corr: {base_metrics['spearman_corr']:.6f}",
        f"- rmse_gain_ratio: {base_metrics['rmse_gain_ratio']:.6f}",
        "",
        "## Permutation Grid",
    ]

    for r in rows:
        lines.append(
            f"- n_perm={r['n_perm']}: "
            f"p_corr_med={r['p_corr_median']:.5f} [q10={r['p_corr_q10']:.5f}, q90={r['p_corr_q90']:.5f}], "
            f"p_gain_med={r['p_gain_median']:.5f} [q10={r['p_gain_q10']:.5f}, q90={r['p_gain_q90']:.5f}], "
            f"frac_corr<=0.01={r['frac_corr_le_0p01']:.3f}, "
            f"frac_gain<=0.01={r['frac_gain_le_0p01']:.3f}"
        )

    lines.extend(
        [
            "",
            "## Robust Flags",
            f"- robust_corr_p_le_0p01: {robust_corr}",
            f"- robust_gain_p_le_0p01: {robust_gain}",
            "",
            "## Artifacts",
            f"- JSON: `{OUT_JSON.name}`",
        ]
    )

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2018] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2018] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2018] verdict={verdict}")


if __name__ == "__main__":
    main()
