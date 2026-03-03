#!/usr/bin/env python3
"""
QW-2036: Primary multifold stability repair via stratified fold protocol.

Goal:
- reduce fold-instability on primary external dataset without retuning kernel,
- enforce deterministic stratified folds by intervention/regime/theta-bin.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2036_primary_stratified_multifold_repair.json"
OUT_MD = ROOT / "RAPORT_QW2036_PRIMARY_STRATIFIED_MULTIFOLD_REPAIR.md"


def hash_u64(s: str) -> int:
    d = hashlib.sha256(s.encode("utf-8")).digest()
    return int.from_bytes(d[:8], "big", signed=False)


def safe_corr(a: np.ndarray, b: np.ndarray) -> float:
    sa = float(np.std(a))
    sb = float(np.std(b))
    if sa <= 1e-12 or sb <= 1e-12:
        return 0.0
    return float(np.corrcoef(a, b)[0, 1])


def kernel_eta(d_eff: np.ndarray, omega: float, phi: float, beta: float, eta: float) -> np.ndarray:
    return np.cos(omega * d_eff + phi) / (1.0 + beta * (d_eff**eta))


def assign_fold_hash(pair_id: Iterable[str], k: int, salt: str) -> np.ndarray:
    return np.array([hash_u64(f"{salt}|{x}") % k for x in pair_id], dtype=int)


def assign_fold_stratified(
    df: pd.DataFrame,
    k: int,
    salt: str,
    strata_cols: list[str],
) -> np.ndarray:
    fold = np.full(len(df), -1, dtype=int)
    work = df.copy()
    work["_idx"] = np.arange(len(work), dtype=int)
    work["_pair_hash"] = work["pair_id"].astype(str).map(lambda x: hash_u64(f"{salt}|{x}"))

    grouped = work.groupby(strata_cols, sort=True, dropna=False)
    for _, g in grouped:
        gg = g.sort_values("_pair_hash")
        idx = gg["_idx"].to_numpy(dtype=int)
        for i, ii in enumerate(idx):
            fold[ii] = int(i % k)
    if np.any(fold < 0):
        raise RuntimeError("Stratified fold assignment failed.")
    return fold


def eval_folds(
    df: pd.DataFrame,
    fold: np.ndarray,
    k_folds: int,
    kernel: Dict[str, float],
    n_perm: int,
    seed0: int,
) -> Dict:
    pair_id = df["pair_id"].astype(str).to_numpy()
    theta = df["theta_deg"].to_numpy(dtype=float)
    y = df["hxy"].to_numpy(dtype=float)

    tmin = float(np.min(theta))
    tmax = float(np.max(theta))
    d_eff = 1.0 + 11.0 * (theta - tmin) / max(tmax - tmin, 1e-12)
    k = kernel_eta(d_eff, kernel["omega"], kernel["phi"], kernel["beta"], kernel["eta"])

    rows = []
    for f in range(k_folds):
        hold = fold == f
        disc = ~hold
        if int(np.sum(hold)) < 120 or int(np.sum(disc)) < 300:
            raise RuntimeError(f"Fold {f} too small.")

        x = np.column_stack([k[disc], np.ones(int(np.sum(disc)), dtype=float)])
        coef, *_ = np.linalg.lstsq(x, y[disc], rcond=None)
        a_hat, b_hat = float(coef[0]), float(coef[1])
        yhat = a_hat * k + b_hat

        y_hold = y[hold]
        yp_hold = yhat[hold]
        corr = safe_corr(y_hold, yp_hold)

        rmse = float(np.sqrt(np.mean((y_hold - yp_hold) ** 2)))
        base = float(np.mean(y[disc]))
        rmse0 = float(np.sqrt(np.mean((y_hold - base) ** 2)))
        gain = float((rmse0 - rmse) / max(rmse0, 1e-12))

        rng = np.random.default_rng(seed0 + 13 * f)
        corr_null = np.empty(n_perm, dtype=float)
        gain_null = np.empty(n_perm, dtype=float)
        for i in range(n_perm):
            yp = rng.permutation(y_hold)
            corr_null[i] = safe_corr(yp, yp_hold)
            rmse_p = float(np.sqrt(np.mean((yp - yp_hold) ** 2)))
            rmse_b = float(np.sqrt(np.mean((yp - base) ** 2)))
            gain_null[i] = float((rmse_b - rmse_p) / max(rmse_b, 1e-12))

        q95_corr = float(np.quantile(corr_null, 0.95))
        q95_gain = float(np.quantile(gain_null, 0.95))
        p_corr = float((1 + np.sum(corr_null >= corr)) / (n_perm + 1))
        p_gain = float((1 + np.sum(gain_null >= gain)) / (n_perm + 1))

        rows.append(
            {
                "fold": int(f),
                "n_holdout": int(np.sum(hold)),
                "n_discovery": int(np.sum(disc)),
                "corr": corr,
                "gain": gain,
                "q95_corr_null": q95_corr,
                "q95_gain_null": q95_gain,
                "p_corr": p_corr,
                "p_gain": p_gain,
            }
        )

    arr_corr = np.array([r["corr"] for r in rows], dtype=float)
    arr_gain = np.array([r["gain"] for r in rows], dtype=float)
    arr_pc = np.array([r["p_corr"] for r in rows], dtype=float)
    arr_pg = np.array([r["p_gain"] for r in rows], dtype=float)
    arr_q95 = np.array([r["q95_corr_null"] for r in rows], dtype=float)

    summary = {
        "n_rows": int(len(df)),
        "k_folds": int(k_folds),
        "median_corr": float(np.median(arr_corr)),
        "median_gain": float(np.median(arr_gain)),
        "min_corr": float(np.min(arr_corr)),
        "min_gain": float(np.min(arr_gain)),
        "median_p_corr": float(np.median(arr_pc)),
        "median_p_gain": float(np.median(arr_pg)),
        "frac_p_corr_le_0p05": float(np.mean(arr_pc <= 0.05)),
        "frac_p_corr_le_0p01": float(np.mean(arr_pc <= 0.01)),
        "frac_gain_positive": float(np.mean(arr_gain > 0.0)),
        "median_q95_corr_null": float(np.median(arr_q95)),
    }
    flags = {
        "median_p_corr_le_0p05": bool(summary["median_p_corr"] <= 0.05),
        "frac_p_corr_le_0p05_ge_0p80": bool(summary["frac_p_corr_le_0p05"] >= 0.80),
        "frac_gain_positive_ge_0p80": bool(summary["frac_gain_positive"] >= 0.80),
        "median_corr_gt_median_q95_null": bool(summary["median_corr"] > summary["median_q95_corr_null"]),
    }
    return {
        "fold_rows": rows,
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
    d2035 = json.loads((ROOT / "report_qw2035_eta_primary_signal_multifold_stability.json").read_text(encoding="utf-8"))
    baseline_primary = d2035["datasets"]["primary"]["summary"]

    p_primary = ROOT / "external_confirmatory_v2/beta_channel_true_external_v2/beta_channel_pairs.csv"
    p_stress = ROOT / "external_confirmatory_v2/confirmatory_dataset_external_source_alpha6_1831cfg/pta_v2_pairs.csv"
    dfp = pd.read_csv(p_primary)
    dfs = pd.read_csv(p_stress)

    # Deterministic theta stratification bins for better distributional balance.
    qbins = pd.qcut(dfp["theta_deg"], q=5, labels=False, duplicates="drop")
    dfp = dfp.copy()
    dfp["theta_bin"] = qbins.astype(int)

    k_folds = 5
    n_perm = 2000

    fold_primary = assign_fold_stratified(
        df=dfp,
        k=k_folds,
        salt="QW2036_PRIMARY_STRAT",
        strata_cols=["intervention_id", "regime", "theta_bin"],
    )
    fold_stress = assign_fold_hash(dfs["pair_id"].astype(str), k=k_folds, salt="QW2036_STRESS_HASH")

    primary = eval_folds(dfp, fold_primary, k_folds=k_folds, kernel=kernel, n_perm=n_perm, seed0=203610)
    stress = eval_folds(dfs, fold_stress, k_folds=k_folds, kernel=kernel, n_perm=n_perm, seed0=203650)

    improvement = {
        "primary_median_p_corr_drop": float(baseline_primary["median_p_corr"] - primary["summary"]["median_p_corr"]),
        "primary_frac_p_corr_le_0p05_gain": float(
            primary["summary"]["frac_p_corr_le_0p05"] - baseline_primary["frac_p_corr_le_0p05"]
        ),
        "primary_median_corr_gain": float(primary["summary"]["median_corr"] - baseline_primary["median_corr"]),
    }
    improvement_flags = {
        "median_p_corr_improved": bool(improvement["primary_median_p_corr_drop"] > 0.0),
        "frac_sig_improved": bool(improvement["primary_frac_p_corr_le_0p05_gain"] > 0.0),
    }

    global_flags = {
        "primary_pass_ge_3of4": bool(primary["pass_count"] >= 3),
        "stress_pass_ge_3of4": bool(stress["pass_count"] >= 3),
        "improvement_flags_all": bool(all(improvement_flags.values())),
    }
    pass_count = int(sum(1 for v in global_flags.values() if v))
    total_flags = int(len(global_flags))

    if pass_count == total_flags:
        verdict = "PRIMARY_STRATIFIED_MULTIFOLD_REPAIR_PASS"
        readiness = "SIGNAL_STABILITY_REPAIRED"
    elif pass_count >= 2:
        verdict = "PRIMARY_STRATIFIED_MULTIFOLD_REPAIR_PARTIAL"
        readiness = "SIGNAL_STABILITY_PARTIAL"
    else:
        verdict = "PRIMARY_STRATIFIED_MULTIFOLD_REPAIR_FAIL"
        readiness = "SIGNAL_STABILITY_NOT_REPAIRED"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "kernel_source": "report_qw2030_final_stage_c_gate_combined_branch.json:kernel",
        "baseline_source": "report_qw2035_eta_primary_signal_multifold_stability.json:datasets.primary.summary",
        "config": {"k_folds": k_folds, "n_perm_per_fold": n_perm, "primary_strata": ["intervention_id", "regime", "theta_bin"]},
        "primary": primary,
        "stress": stress,
        "baseline_primary_summary_qw2035": baseline_primary,
        "improvement": improvement,
        "improvement_flags": improvement_flags,
        "global_flags": global_flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "readiness": readiness,
        "verdict": verdict,
        "required_next_step": (
            "USE_STRATIFIED_PROTOCOL_AS_PRIMARY_CONFIRMATORY_TEMPLATE"
            if verdict == "PRIMARY_STRATIFIED_MULTIFOLD_REPAIR_PASS"
            else "CONTINUE_PRIMARY_SIGNAL_REPAIR"
        ),
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    p = primary["summary"]
    s = stress["summary"]
    lines = [
        "# RAPORT QW-2036: PRIMARY STRATIFIED MULTIFOLD REPAIR",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Readiness: **{readiness}**",
        f"- Verdict: **{verdict}**",
        f"- global pass_count: {pass_count}/{total_flags}",
        "",
        "## Primary (Stratified)",
        (
            f"- median corr/gain: {p['median_corr']:.6f}/{p['median_gain']:.6f}, "
            f"median p_corr: {p['median_p_corr']:.6f}, frac p<=0.05: {p['frac_p_corr_le_0p05']:.3f}"
        ),
        "",
        "## Stress",
        (
            f"- median corr/gain: {s['median_corr']:.6f}/{s['median_gain']:.6f}, "
            f"median p_corr: {s['median_p_corr']:.6f}, frac p<=0.05: {s['frac_p_corr_le_0p05']:.3f}"
        ),
        "",
        "## Improvement vs QW-2035 Primary",
        f"- median_p_corr_drop: {improvement['primary_median_p_corr_drop']:.6f}",
        f"- frac_p_corr_le_0p05_gain: {improvement['primary_frac_p_corr_le_0p05_gain']:.6f}",
        f"- median_corr_gain: {improvement['primary_median_corr_gain']:.6f}",
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

    print(f"[QW-2036] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2036] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2036] readiness={readiness} verdict={verdict} pass={pass_count}/{total_flags}")


if __name__ == "__main__":
    main()
