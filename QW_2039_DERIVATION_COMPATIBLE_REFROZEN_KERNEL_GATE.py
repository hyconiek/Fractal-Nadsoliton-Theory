#!/usr/bin/env python3
"""
QW-2039: Gate for derivation-compatible refrozen kernel candidate.

Steps:
1) select a QW-2038 pass candidate with beta inside QW-2034 beta CI95,
2) verify Stage-C hard flags for this candidate,
3) run blind external transfer check (QW-2031 style) on primary + stress.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2039_derivation_compatible_refrozen_kernel_gate.json"
OUT_MD = ROOT / "RAPORT_QW2039_DERIVATION_COMPATIBLE_REFROZEN_KERNEL_GATE.md"


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


def safe_corr(a: np.ndarray, b: np.ndarray) -> float:
    sa = float(np.std(a))
    sb = float(np.std(b))
    if sa <= 1e-12 or sb <= 1e-12:
        return 0.0
    return float(np.corrcoef(a, b)[0, 1])


def spearmanr(a: np.ndarray, b: np.ndarray) -> float:
    ra = rankdata(a)
    rb = rankdata(b)
    return safe_corr(ra, rb)


def kernel_eta(d_eff: np.ndarray, omega: float, phi: float, beta: float, eta: float) -> np.ndarray:
    return np.cos(omega * d_eff + phi) / (1.0 + beta * (d_eff**eta))


def eval_one_dataset(
    df: pd.DataFrame,
    kernel: Dict[str, float],
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

    tmin = float(np.min(theta))
    tmax = float(np.max(theta))
    d_eff = 1.0 + 11.0 * (theta - tmin) / max(tmax - tmin, 1e-12)
    k = kernel_eta(d_eff, kernel["omega"], kernel["phi"], kernel["beta"], kernel["eta"])

    fold = np.array([split_index(x, k=2) for x in pair_id], dtype=int)
    disc = fold == 0
    hold = fold == 1
    if int(np.sum(disc)) < 300 or int(np.sum(hold)) < 300:
        raise RuntimeError("Discovery/holdout split too small.")

    x = np.column_stack([k[disc], np.ones(int(np.sum(disc)), dtype=float)])
    coef, *_ = np.linalg.lstsq(x, y[disc], rcond=None)
    a_hat, b_hat = float(coef[0]), float(coef[1])
    yhat = a_hat * k + b_hat

    y_hold = y[hold]
    yhat_hold = yhat[hold]
    corr_hold = safe_corr(y_hold, yhat_hold)
    spear_hold = spearmanr(y_hold, yhat_hold)

    rmse_hold = float(np.sqrt(np.mean((y_hold - yhat_hold) ** 2)))
    base = float(np.mean(y[disc]))
    rmse_base = float(np.sqrt(np.mean((y_hold - base) ** 2)))
    rmse_gain = float((rmse_base - rmse_hold) / max(rmse_base, 1e-12))

    rng = np.random.default_rng(rng_seed)
    corr_null = np.empty(n_perm, dtype=float)
    gain_null = np.empty(n_perm, dtype=float)
    for i in range(n_perm):
        yp = rng.permutation(y_hold)
        corr_null[i] = safe_corr(yp, yhat_hold)
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
    d2034 = json.loads((ROOT / "report_qw2034_eta_kernel_derivational_stability_audit.json").read_text(encoding="utf-8"))
    d2038 = json.loads((ROOT / "report_qw2038_derivation_compatible_kernel_refreeze_scan.json").read_text(encoding="utf-8"))

    beta_upper = float(d2034["bootstrap_stats"]["beta"]["q97p5"])
    candidates = [r for r in d2038.get("pass_examples", []) if float(r["kernel"]["beta"]) <= beta_upper]
    if not candidates:
        raise RuntimeError("No QW-2038 pass candidate with beta inside QW-2034 CI95 upper bound.")

    selected = min(candidates, key=lambda r: float(r["score"]))
    kernel = {
        "omega": float(selected["kernel"]["omega"]),
        "phi": float(selected["kernel"]["phi"]),
        "beta": float(selected["kernel"]["beta"]),
        "eta": float(selected["kernel"]["eta"]),
    }

    # Stage-C hard flags reused from candidate.
    stagec_flags = selected["flags"]
    stagec_all_pass = bool(all(stagec_flags.values()))

    p_primary = ROOT / "external_confirmatory_v2/beta_channel_true_external_v2/beta_channel_pairs.csv"
    p_stress = ROOT / "external_confirmatory_v2/confirmatory_dataset_external_source_alpha6_1831cfg/pta_v2_pairs.csv"
    primary = eval_one_dataset(pd.read_csv(p_primary), kernel=kernel, n_perm=5000, rng_seed=203901)
    stress = eval_one_dataset(pd.read_csv(p_stress), kernel=kernel, n_perm=5000, rng_seed=203902)

    primary_pass = bool(primary["all_pass"])
    stress_soft_pass = bool(
        stress["holdout_metrics"]["pearson_corr"] > stress["permutation"]["q95_corr"]
        and stress["holdout_metrics"]["rmse_gain_ratio"] > 0.0
        and stress["permutation"]["p_corr"] <= 0.05
    )

    global_flags = {
        "beta_inside_derivational_ci95": bool(kernel["beta"] <= beta_upper),
        "stagec_candidate_all_pass": stagec_all_pass,
        "external_primary_all_pass": primary_pass,
        "external_stress_soft_pass": stress_soft_pass,
    }
    pass_count = int(sum(1 for v in global_flags.values() if v))
    total_flags = int(len(global_flags))

    if pass_count == total_flags:
        verdict = "DERIVATION_COMPATIBLE_REFROZEN_KERNEL_GATE_PASS"
        readiness = "TOE_INTERNAL_GAPS_MINIMIZED_PENDING_EXTERNAL_MULTITEAM_AUDIT"
    elif pass_count >= 3:
        verdict = "DERIVATION_COMPATIBLE_REFROZEN_KERNEL_GATE_PARTIAL"
        readiness = "PARTIAL"
    else:
        verdict = "DERIVATION_COMPATIBLE_REFROZEN_KERNEL_GATE_FAIL"
        readiness = "NOT_READY"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "candidate_pool": "report_qw2038_derivation_compatible_kernel_refreeze_scan.json:pass_examples",
            "beta_ci": "report_qw2034_eta_kernel_derivational_stability_audit.json:bootstrap_stats.beta",
        },
        "selected_candidate": selected,
        "selected_kernel": kernel,
        "beta_ci95_upper_from_qw2034": beta_upper,
        "stagec_flags": stagec_flags,
        "external": {
            "primary": primary,
            "stress": stress,
            "pass_flags": {
                "primary_all_pass": primary_pass,
                "stress_soft_pass": stress_soft_pass,
            },
        },
        "global_flags": global_flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "readiness": readiness,
        "verdict": verdict,
        "required_next_step": (
            "PROMOTE_SELECTED_KERNEL_TO_NEW_FROZEN_BASELINE_AND_PREPARE_EXTERNAL_MULTITEAM_AUDIT"
            if verdict == "DERIVATION_COMPATIBLE_REFROZEN_KERNEL_GATE_PASS"
            else "CONTINUE_BETA_IDENTIFIABILITY_REPAIR"
        ),
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    pm = primary["holdout_metrics"]
    pp = primary["permutation"]
    sm = stress["holdout_metrics"]
    sp = stress["permutation"]
    lines = [
        "# RAPORT QW-2039: DERIVATION-COMPATIBLE REFROZEN KERNEL GATE",
        "",
        f"- Data UTC: {out['generated_utc']}",
        (
            f"- Selected kernel omega/phi/beta/eta: {kernel['omega']:.6f} / "
            f"{kernel['phi']:.6f} / {kernel['beta']:.6f} / {kernel['eta']:.6f}"
        ),
        f"- QW-2034 beta CI95 upper: {beta_upper:.6f}",
        f"- Readiness: **{readiness}**",
        f"- Verdict: **{verdict}**",
        f"- pass_count: {pass_count}/{total_flags}",
        "",
        "## External Primary",
        f"- corr: {pm['pearson_corr']:.6f} (q95: {pp['q95_corr']:.6f}, p={pp['p_corr']:.6f})",
        f"- gain: {pm['rmse_gain_ratio']:.6f} (q95: {pp['q95_rmse_gain']:.6f}, p={pp['p_rmse_gain']:.6f})",
        "",
        "## External Stress",
        f"- corr: {sm['pearson_corr']:.6f} (q95: {sp['q95_corr']:.6f}, p={sp['p_corr']:.6f})",
        f"- gain: {sm['rmse_gain_ratio']:.6f} (q95: {sp['q95_rmse_gain']:.6f}, p={sp['p_rmse_gain']:.6f})",
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

    print(f"[QW-2039] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2039] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2039] readiness={readiness} verdict={verdict} pass={pass_count}/{total_flags}")


if __name__ == "__main__":
    main()
