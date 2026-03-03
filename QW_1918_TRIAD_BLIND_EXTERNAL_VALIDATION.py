#!/usr/bin/env python3
"""
QW-1918: Blind external validation of triad (omega, phi, beta) from QW-1917.

Protocol:
- no tuning of (omega, phi, beta) on external data,
- deterministic discovery/holdout split by pair_id hash,
- only nuisance affine calibration (a, b) fitted on discovery,
- holdout tested against permutation null.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1918_triad_blind_external_validation.json"
OUT_MD = ROOT / "RAPORT_QW1918_TRIAD_BLIND_EXTERNAL_VALIDATION.md"


def load_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


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


def eval_one_dataset(
    df: pd.DataFrame,
    omega: float,
    phi: float,
    beta: float,
    n_perm: int,
    rng_seed: int,
) -> Dict:
    req = {"pair_id", "theta_deg", "hxy"}
    miss = sorted(req - set(df.columns))
    if miss:
        raise RuntimeError(f"Missing columns: {miss}")

    pair_id = df["pair_id"].astype(str).to_numpy()
    theta = df["theta_deg"].to_numpy(dtype=float)
    y = df["hxy"].to_numpy(dtype=float)

    # Map theta to an effective lattice distance in [1,12].
    tmin, tmax = float(np.min(theta)), float(np.max(theta))
    d_eff = 1.0 + 11.0 * (theta - tmin) / max(tmax - tmin, 1e-12)

    k = np.cos(omega * d_eff + phi) / (1.0 + beta * d_eff)

    fold = np.array([split_index(x, k=2) for x in pair_id], dtype=int)
    disc = fold == 0
    hold = fold == 1

    if int(np.sum(disc)) < 300 or int(np.sum(hold)) < 300:
        raise RuntimeError("Discovery/holdout split too small.")

    # Fit nuisance affine map on discovery only.
    X = np.column_stack([k[disc], np.ones(int(np.sum(disc)), dtype=float)])
    coef, *_ = np.linalg.lstsq(X, y[disc], rcond=None)
    a_hat, b_hat = float(coef[0]), float(coef[1])

    yhat = a_hat * k + b_hat
    y_hold = y[hold]
    yhat_hold = yhat[hold]

    corr_hold = float(np.corrcoef(y_hold, yhat_hold)[0, 1])
    spear_hold = float(spearmanr(y_hold, yhat_hold))

    rmse_hold = float(np.sqrt(np.mean((y_hold - yhat_hold) ** 2)))
    base = float(np.mean(y[disc]))
    rmse_base = float(np.sqrt(np.mean((y_hold - base) ** 2)))
    rmse_gain = float((rmse_base - rmse_hold) / max(rmse_base, 1e-12))

    # Permutation null on holdout.
    rng = np.random.default_rng(rng_seed)
    corr_null = np.empty(n_perm, dtype=float)
    gain_null = np.empty(n_perm, dtype=float)

    for i in range(n_perm):
        yp = rng.permutation(y_hold)
        corr_null[i] = float(np.corrcoef(yp, yhat_hold)[0, 1])
        rmse_p = float(np.sqrt(np.mean((yp - yhat_hold) ** 2)))
        rmse_b = float(np.sqrt(np.mean((yp - base) ** 2)))
        gain_null[i] = float((rmse_b - rmse_p) / max(rmse_b, 1e-12))

    p_corr = float((1 + np.sum(corr_null >= corr_hold)) / (n_perm + 1))
    p_gain = float((1 + np.sum(gain_null >= rmse_gain)) / (n_perm + 1))

    q95_corr = float(np.quantile(corr_null, 0.95))
    q95_gain = float(np.quantile(gain_null, 0.95))

    flags = {
        "corr_gt_perm_q95": bool(corr_hold > q95_corr),
        "rmse_gain_gt_perm_q95": bool(rmse_gain > q95_gain),
        "p_corr_le_0p01": bool(p_corr <= 0.01),
        "p_gain_le_0p01": bool(p_gain <= 0.01),
        "rmse_gain_positive": bool(rmse_gain > 0.0),
    }

    return {
        "n_rows": int(len(df)),
        "n_discovery": int(np.sum(disc)),
        "n_holdout": int(np.sum(hold)),
        "nuisance_affine": {"a": a_hat, "b": b_hat},
        "holdout_metrics": {
            "pearson_corr": corr_hold,
            "spearman_corr": spear_hold,
            "rmse": rmse_hold,
            "rmse_base": rmse_base,
            "rmse_gain_ratio": rmse_gain,
        },
        "permutation": {
            "n_perm": int(n_perm),
            "q95_corr": q95_corr,
            "q95_rmse_gain": q95_gain,
            "p_corr": p_corr,
            "p_rmse_gain": p_gain,
        },
        "pass_flags": flags,
        "all_pass": bool(all(flags.values())),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--primary-pta",
        type=str,
        default="external_confirmatory_v2/confirmatory_dataset_external_source_rebuild_v2_1831cfg/pta_v2_pairs.csv",
    )
    ap.add_argument(
        "--stress-pta",
        type=str,
        default="external_confirmatory_v2/confirmatory_dataset_external_source_alpha6_1831cfg/pta_v2_pairs.csv",
    )
    ap.add_argument("--n-perm", type=int, default=5000)
    args = ap.parse_args()

    d1917 = load_json("report_qw1917_triad_derivation_no_ansatz_profile.json")
    triad = d1917.get("optimum", {})
    omega = float(triad.get("omega"))
    phi = float(triad.get("phi"))
    beta = float(triad.get("beta"))

    p_primary = (ROOT / args.primary_pta).resolve()
    p_stress = (ROOT / args.stress_pta).resolve()
    if not p_primary.exists():
        raise RuntimeError(f"Primary PTA not found: {p_primary}")
    if not p_stress.exists():
        raise RuntimeError(f"Stress PTA not found: {p_stress}")

    df_primary = pd.read_csv(p_primary)
    df_stress = pd.read_csv(p_stress)

    primary = eval_one_dataset(df_primary, omega=omega, phi=phi, beta=beta, n_perm=args.n_perm, rng_seed=191801)
    stress = eval_one_dataset(df_stress, omega=omega, phi=phi, beta=beta, n_perm=args.n_perm, rng_seed=191802)

    primary_pass = bool(primary["all_pass"])
    stress_soft_pass = bool(
        stress["holdout_metrics"]["pearson_corr"] > stress["permutation"]["q95_corr"]
        and stress["holdout_metrics"]["rmse_gain_ratio"] > 0.0
        and stress["permutation"]["p_corr"] <= 0.05
    )

    if primary_pass and stress_soft_pass:
        verdict = "TRIAD_BLIND_EXTERNAL_VALIDATION_PASS_STRONG"
    elif primary_pass:
        verdict = "TRIAD_BLIND_EXTERNAL_VALIDATION_PASS_PRIMARY_ONLY"
    else:
        verdict = "TRIAD_BLIND_EXTERNAL_VALIDATION_FAIL"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "triad_source": "report_qw1917_triad_derivation_no_ansatz_profile.json",
        "triad": {"omega": omega, "phi": phi, "beta": beta},
        "datasets": {
            "primary": {
                "path": str(p_primary),
                **primary,
            },
            "stress": {
                "path": str(p_stress),
                **stress,
            },
        },
        "pass_flags": {
            "primary_all_pass": primary_pass,
            "stress_soft_pass": stress_soft_pass,
        },
        "verdict": verdict,
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    pm = primary["holdout_metrics"]
    pp = primary["permutation"]
    sm = stress["holdout_metrics"]
    sp = stress["permutation"]

    lines = [
        "# RAPORT QW-1918: TRIAD BLIND EXTERNAL VALIDATION",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Triad: omega={omega:.6f}, phi={phi:.6f}, beta={beta:.6f}",
        f"- Verdict: **{verdict}**",
        "",
        "## Primary Dataset (blind holdout)",
        f"- path: `{p_primary}`",
        f"- pearson_corr: {pm['pearson_corr']:.4f} (q95 null: {pp['q95_corr']:.4f}, p={pp['p_corr']:.4g})",
        f"- spearman_corr: {pm['spearman_corr']:.4f}",
        f"- rmse_gain_ratio: {pm['rmse_gain_ratio']:.4f} (q95 null: {pp['q95_rmse_gain']:.4f}, p={pp['p_rmse_gain']:.4g})",
        f"- all_pass: {primary['all_pass']}",
        "",
        "## Stress Dataset (blind holdout)",
        f"- path: `{p_stress}`",
        f"- pearson_corr: {sm['pearson_corr']:.4f} (q95 null: {sp['q95_corr']:.4f}, p={sp['p_corr']:.4g})",
        f"- spearman_corr: {sm['spearman_corr']:.4f}",
        f"- rmse_gain_ratio: {sm['rmse_gain_ratio']:.4f} (q95 null: {sp['q95_rmse_gain']:.4f}, p={sp['p_rmse_gain']:.4g})",
        f"- all_pass: {stress['all_pass']}",
        "",
        "## Gate Flags",
        f"- primary_all_pass: {primary_pass}",
        f"- stress_soft_pass: {stress_soft_pass}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1918] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1918] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1918] verdict={verdict}")


if __name__ == "__main__":
    main()
