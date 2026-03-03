#!/usr/bin/env python3
"""
QW-1933: Hash-split robustness audit for reparameterized triad (eta branch).

Goal:
- verify that QW-1932 selected candidate is not a single-split artifact,
- run blind discovery/holdout evaluation across multiple deterministic hash salts,
- quantify pass-rate under permutation-null criteria.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1933_reparam_hash_split_robustness.json"
OUT_MD = ROOT / "RAPORT_QW1933_REPARAM_HASH_SPLIT_ROBUSTNESS.md"


def load_qw1932_selected() -> Dict[str, float]:
    p = ROOT / "report_qw1932_physical_reparameterization_eta_scan.json"
    d = json.loads(p.read_text(encoding="utf-8"))
    sel = d.get("selected", {})
    fit = sel.get("fit", {})
    return {
        "omega": float(fit["omega"]),
        "phi": float(fit["phi"]),
        "beta": float(fit["beta"]),
        "eta": float(sel["eta"]),
    }


def split_index(pair_id: str, salt: int, k: int = 2) -> int:
    h = hashlib.sha256(f"{salt}|{pair_id}".encode("utf-8")).hexdigest()
    return int(h[-8:], 16) % k


def eval_one_salt(
    df: pd.DataFrame,
    triad: Dict[str, float],
    salt: int,
    n_perm: int,
    seed_base: int,
) -> Dict[str, object]:
    req = {"pair_id", "theta_deg", "hxy"}
    miss = sorted(req - set(df.columns))
    if miss:
        raise RuntimeError(f"Missing columns: {miss}")

    pair_id = df["pair_id"].astype(str).to_numpy()
    theta = df["theta_deg"].to_numpy(dtype=float)
    y = df["hxy"].to_numpy(dtype=float)

    tmin, tmax = float(np.min(theta)), float(np.max(theta))
    d_eff = 1.0 + 11.0 * (theta - tmin) / max(tmax - tmin, 1e-12)

    omega = float(triad["omega"])
    phi = float(triad["phi"])
    beta = float(triad["beta"])
    eta = float(triad["eta"])
    k = np.cos(omega * d_eff + phi) / (1.0 + beta * (d_eff ** eta))

    fold = np.array([split_index(x, salt=salt, k=2) for x in pair_id], dtype=int)
    disc = fold == 0
    hold = fold == 1

    n_disc = int(np.sum(disc))
    n_hold = int(np.sum(hold))
    if n_disc < 300 or n_hold < 300:
        raise RuntimeError(f"Salt {salt}: too small split n_disc={n_disc} n_hold={n_hold}")

    X = np.column_stack([k[disc], np.ones(n_disc, dtype=float)])
    coef, *_ = np.linalg.lstsq(X, y[disc], rcond=None)
    a_hat, b_hat = float(coef[0]), float(coef[1])

    yhat = a_hat * k + b_hat
    y_hold = y[hold]
    yhat_hold = yhat[hold]

    corr_hold = float(np.corrcoef(y_hold, yhat_hold)[0, 1])
    rmse_hold = float(np.sqrt(np.mean((y_hold - yhat_hold) ** 2)))
    base = float(np.mean(y[disc]))
    rmse_base = float(np.sqrt(np.mean((y_hold - base) ** 2)))
    rmse_gain = float((rmse_base - rmse_hold) / max(rmse_base, 1e-12))

    rng = np.random.default_rng(seed_base + 97 * salt)
    corr_null = np.empty(n_perm, dtype=float)
    gain_null = np.empty(n_perm, dtype=float)

    for i in range(n_perm):
        yp = rng.permutation(y_hold)
        corr_null[i] = float(np.corrcoef(yp, yhat_hold)[0, 1])
        rmse_p = float(np.sqrt(np.mean((yp - yhat_hold) ** 2)))
        rmse_b = float(np.sqrt(np.mean((yp - base) ** 2)))
        gain_null[i] = float((rmse_b - rmse_p) / max(rmse_b, 1e-12))

    q95_corr = float(np.quantile(corr_null, 0.95))
    q95_gain = float(np.quantile(gain_null, 0.95))
    p_corr = float((1 + np.sum(corr_null >= corr_hold)) / (n_perm + 1))
    p_gain = float((1 + np.sum(gain_null >= rmse_gain)) / (n_perm + 1))

    pass_flags = {
        "corr_gt_perm_q95": bool(corr_hold > q95_corr),
        "rmse_gain_gt_perm_q95": bool(rmse_gain > q95_gain),
        "p_corr_le_0p01": bool(p_corr <= 0.01),
        "p_gain_le_0p01": bool(p_gain <= 0.01),
        "rmse_gain_positive": bool(rmse_gain > 0.0),
    }

    return {
        "salt": int(salt),
        "n_discovery": n_disc,
        "n_holdout": n_hold,
        "nuisance_affine": {"a": a_hat, "b": b_hat},
        "holdout": {
            "pearson_corr": corr_hold,
            "rmse_gain_ratio": rmse_gain,
            "rmse": rmse_hold,
            "rmse_base": rmse_base,
        },
        "permutation": {
            "n_perm": int(n_perm),
            "q95_corr": q95_corr,
            "q95_rmse_gain": q95_gain,
            "p_corr": p_corr,
            "p_rmse_gain": p_gain,
        },
        "pass_flags": pass_flags,
        "all_pass": bool(all(pass_flags.values())),
    }


def summarize_salts(rows: List[Dict[str, object]]) -> Dict[str, object]:
    all_pass = np.array([bool(r["all_pass"]) for r in rows], dtype=bool)
    corr = np.array([float(r["holdout"]["pearson_corr"]) for r in rows], dtype=float)
    gain = np.array([float(r["holdout"]["rmse_gain_ratio"]) for r in rows], dtype=float)

    out = {
        "n_salts": int(len(rows)),
        "pass_rate": float(np.mean(all_pass)),
        "corr_median": float(np.median(corr)),
        "corr_q10": float(np.quantile(corr, 0.10)),
        "corr_q90": float(np.quantile(corr, 0.90)),
        "gain_median": float(np.median(gain)),
        "gain_q10": float(np.quantile(gain, 0.10)),
        "gain_q90": float(np.quantile(gain, 0.90)),
        "all_pass_count": int(np.sum(all_pass)),
    }
    return out


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
    ap.add_argument("--n-perm", type=int, default=800)
    args = ap.parse_args()

    triad = load_qw1932_selected()
    p_primary = (ROOT / args.primary_pta).resolve()
    p_stress = (ROOT / args.stress_pta).resolve()
    if not p_primary.exists():
        raise RuntimeError(f"Primary PTA not found: {p_primary}")
    if not p_stress.exists():
        raise RuntimeError(f"Stress PTA not found: {p_stress}")

    d_primary = pd.read_csv(p_primary)
    d_stress = pd.read_csv(p_stress)

    salts = list(range(int(args.n_salts)))
    rows_primary = [
        eval_one_salt(d_primary, triad=triad, salt=s, n_perm=int(args.n_perm), seed_base=193301) for s in salts
    ]
    rows_stress = [
        eval_one_salt(d_stress, triad=triad, salt=s, n_perm=int(args.n_perm), seed_base=193302) for s in salts
    ]

    sum_primary = summarize_salts(rows_primary)
    sum_stress = summarize_salts(rows_stress)

    # Strict robustness targets:
    # - primary: at least 70% full-pass across deterministic splits (weak channel),
    # - stress:  at least 90% full-pass across deterministic splits.
    flags = {
        "primary_pass_rate_ge_0p70": bool(sum_primary["pass_rate"] >= 0.70),
        "stress_pass_rate_ge_0p90": bool(sum_stress["pass_rate"] >= 0.90),
        "primary_corr_median_positive": bool(sum_primary["corr_median"] > 0.0),
        "primary_gain_median_positive": bool(sum_primary["gain_median"] > 0.0),
        "stress_corr_median_positive": bool(sum_stress["corr_median"] > 0.0),
        "stress_gain_median_positive": bool(sum_stress["gain_median"] > 0.0),
    }

    if all(flags.values()):
        verdict = "REPARAM_HASH_SPLIT_ROBUSTNESS_PASS"
        required_next = "INTEGRATE_QW1932_QW1933_IN_STRICT_CLOSURE_GATE"
    else:
        verdict = "REPARAM_HASH_SPLIT_ROBUSTNESS_PARTIAL_OR_FAIL"
        required_next = "INCREASE_SIGNAL_CHANNEL_OR_REDEFINE_OBSERVABLE_BEFORE_CLOSURE_CLAIM"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "triad_source": "report_qw1932_physical_reparameterization_eta_scan.json:selected",
        "triad": triad,
        "config": {
            "n_salts": int(args.n_salts),
            "n_perm_per_salt": int(args.n_perm),
        },
        "datasets": {
            "primary": {
                "path": str(p_primary),
                "summary": sum_primary,
                "rows": rows_primary,
            },
            "stress": {
                "path": str(p_stress),
                "summary": sum_stress,
                "rows": rows_stress,
            },
        },
        "flags": flags,
        "verdict": verdict,
        "required_next_step": required_next,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1933: REPARAM HASH-SPLIT ROBUSTNESS",
        "",
        f"- Data UTC: {out['generated_utc']}",
        (
            "- Triad (eta branch): "
            f"omega={triad['omega']:.6f}, phi={triad['phi']:.6f}, beta={triad['beta']:.6f}, eta={triad['eta']:.2f}"
        ),
        f"- Verdict: **{verdict}**",
        "",
        "## Primary Summary",
        f"- pass_rate: {sum_primary['pass_rate']:.4f} ({sum_primary['all_pass_count']}/{sum_primary['n_salts']})",
        f"- corr median [q10, q90]: {sum_primary['corr_median']:.4f} [{sum_primary['corr_q10']:.4f}, {sum_primary['corr_q90']:.4f}]",
        f"- gain median [q10, q90]: {sum_primary['gain_median']:.4f} [{sum_primary['gain_q10']:.4f}, {sum_primary['gain_q90']:.4f}]",
        "",
        "## Stress Summary",
        f"- pass_rate: {sum_stress['pass_rate']:.4f} ({sum_stress['all_pass_count']}/{sum_stress['n_salts']})",
        f"- corr median [q10, q90]: {sum_stress['corr_median']:.4f} [{sum_stress['corr_q10']:.4f}, {sum_stress['corr_q90']:.4f}]",
        f"- gain median [q10, q90]: {sum_stress['gain_median']:.4f} [{sum_stress['gain_q10']:.4f}, {sum_stress['gain_q90']:.4f}]",
        "",
        "## Strict Flags",
    ]
    for k, v in flags.items():
        lines.append(f"- {k}: {v}")
    lines.extend(
        [
            "",
            "## Required Next Step",
            f"- {required_next}",
            "",
            "## Artifacts",
            f"- JSON: `{OUT_JSON.name}`",
        ]
    )
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1933] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1933] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1933] verdict={verdict}")


if __name__ == "__main__":
    main()

