#!/usr/bin/env python3
"""
QW-2046: Micro-derived -> Stage-C intersection gate.

Combines:
- micro-derived support from QW-2045,
- Stage-C pass-candidate pool from QW-2038,
- external blind validation protocol from QW-2039.

No sector retune in this gate.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2046_micro_stagec_intersection_gate.json"
OUT_MD = ROOT / "RAPORT_QW2046_MICRO_STAGEC_INTERSECTION_GATE.md"


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
    return safe_corr(rankdata(a), rankdata(b))


def kernel_eta(d_eff: np.ndarray, omega: float, phi: float, beta: float, eta: float) -> np.ndarray:
    return np.cos(omega * d_eff + phi) / (1.0 + beta * (d_eff**eta))


def eval_one_dataset(df: pd.DataFrame, kernel: Dict[str, float], n_perm: int, rng_seed: int) -> Dict:
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
    d2045 = json.loads((ROOT / "report_qw2045_phase_conditioned_pointwise_micro_derivation.json").read_text(encoding="utf-8"))
    d2038 = json.loads((ROOT / "report_qw2038_derivation_compatible_kernel_refreeze_scan.json").read_text(encoding="utf-8"))

    beta_ci = [float(x) for x in d2045["global_estimates"]["beta_ci95"]]
    eta_ci = [float(x) for x in d2045["global_estimates"]["eta_ci95"]]

    pass_examples: List[Dict[str, object]] = d2038.get("pass_examples", [])

    inter = []
    for r in pass_examples:
        k = r["kernel"]
        b = float(k["beta"])
        e = float(k["eta"])
        if beta_ci[0] <= b <= beta_ci[1] and eta_ci[0] <= e <= eta_ci[1]:
            inter.append(r)

    if not inter:
        # fallback nearest by normalized distance
        b_scale = max(beta_ci[1] - beta_ci[0], 1e-9)
        e_scale = max(eta_ci[1] - eta_ci[0], 1e-9)

        def dist(r):
            k = r["kernel"]
            db = (float(k["beta"]) - 0.5 * (beta_ci[0] + beta_ci[1])) / b_scale
            de = (float(k["eta"]) - 0.5 * (eta_ci[0] + eta_ci[1])) / e_scale
            return float(np.sqrt(db * db + de * de))

        selected = min(pass_examples, key=dist)
        from_intersection = False
    else:
        selected = min(inter, key=lambda r: float(r["score"]))
        from_intersection = True

    kernel = {
        "omega": float(selected["kernel"]["omega"]),
        "phi": float(selected["kernel"]["phi"]),
        "beta": float(selected["kernel"]["beta"]),
        "eta": float(selected["kernel"]["eta"]),
    }

    p_primary = ROOT / "external_confirmatory_v2/beta_channel_true_external_v2/beta_channel_pairs.csv"
    p_stress = ROOT / "external_confirmatory_v2/confirmatory_dataset_external_source_alpha6_1831cfg/pta_v2_pairs.csv"
    primary = eval_one_dataset(pd.read_csv(p_primary), kernel=kernel, n_perm=5000, rng_seed=204601)
    stress = eval_one_dataset(pd.read_csv(p_stress), kernel=kernel, n_perm=5000, rng_seed=204602)

    primary_pass = bool(primary["all_pass"])
    stress_soft_pass = bool(
        stress["holdout_metrics"]["pearson_corr"] > stress["permutation"]["q95_corr"]
        and stress["holdout_metrics"]["rmse_gain_ratio"] > 0.0
        and stress["permutation"]["p_corr"] <= 0.05
    )

    beta_overlap = bool(beta_ci[0] <= kernel["beta"] <= beta_ci[1])
    eta_overlap = bool(eta_ci[0] <= kernel["eta"] <= eta_ci[1])

    flags = {
        "stagec_candidate_all_pass": bool(all(selected["flags"].values())),
        "micro_beta_overlap": beta_overlap,
        "micro_eta_overlap": eta_overlap,
        "external_primary_all_pass": primary_pass,
        "external_stress_soft_pass": stress_soft_pass,
        "micro_pointwise_bins_ge_6": bool(d2045["flags"]["enough_pointwise_bins_ge_6"]),
        "micro_phase_condition_strength_ge_0p75": bool(d2045["flags"]["phase_condition_strength_ge_0p75"]),
    }

    pass_count = int(sum(1 for v in flags.values() if v))
    total_flags = int(len(flags))

    if pass_count == total_flags:
        verdict = "MICRO_STAGEC_INTERSECTION_GATE_PASS"
        readiness = "MICRO_TO_STAGEC_BRIDGE_STRONG"
    elif pass_count >= 5:
        verdict = "MICRO_STAGEC_INTERSECTION_GATE_PARTIAL"
        readiness = "MICRO_TO_STAGEC_BRIDGE_PARTIAL"
    else:
        verdict = "MICRO_STAGEC_INTERSECTION_GATE_FAIL"
        readiness = "MICRO_TO_STAGEC_BRIDGE_NOT_ESTABLISHED"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "micro_support": {
            "beta_ci95": beta_ci,
            "eta_ci95": eta_ci,
            "source": "report_qw2045_phase_conditioned_pointwise_micro_derivation.json",
        },
        "stagec_pool": {
            "pass_examples_count": int(len(pass_examples)),
            "intersection_count": int(len(inter)),
            "selected_from_intersection": bool(from_intersection),
            "selected_candidate": selected,
            "selected_kernel": kernel,
        },
        "external": {
            "primary": primary,
            "stress": stress,
            "pass_flags": {
                "primary_all_pass": primary_pass,
                "stress_soft_pass": stress_soft_pass,
            },
        },
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "readiness": readiness,
        "required_next_step": (
            "ADD_SIGNED_PHASE_TORSION_OBSERVABLE_AND_RE-RUN_POINTWISE_DERIVATION"
            if verdict != "MICRO_STAGEC_INTERSECTION_GATE_PASS"
            else "PROMOTE_MICRO_BRIDGE_TO_STRICT_INTERNAL_BASELINE"
        ),
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    pm = primary["holdout_metrics"]
    pp = primary["permutation"]
    sm = stress["holdout_metrics"]
    sp = stress["permutation"]

    lines = [
        "# RAPORT QW-2046: MICRO STAGE-C INTERSECTION GATE",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- Readiness: **{readiness}**",
        f"- pass_count: {pass_count}/{total_flags}",
        "",
        "## Micro Support (from QW-2045)",
        f"- beta CI95: [{beta_ci[0]:.6f}, {beta_ci[1]:.6f}]",
        f"- eta CI95: [{eta_ci[0]:.6f}, {eta_ci[1]:.6f}]",
        "",
        "## Stage-C Intersection",
        f"- pass_examples_count: {len(pass_examples)}",
        f"- intersection_count: {len(inter)}",
        f"- selected_from_intersection: {from_intersection}",
        f"- selected kernel omega/phi/beta/eta: {kernel['omega']:.6f} / {kernel['phi']:.6f} / {kernel['beta']:.6f} / {kernel['eta']:.6f}",
        "",
        "## External Primary",
        f"- corr: {pm['pearson_corr']:.6f} (q95={pp['q95_corr']:.6f}, p={pp['p_corr']:.6f})",
        f"- gain: {pm['rmse_gain_ratio']:.6f} (q95={pp['q95_rmse_gain']:.6f}, p={pp['p_rmse_gain']:.6f})",
        "",
        "## External Stress",
        f"- corr: {sm['pearson_corr']:.6f} (q95={sp['q95_corr']:.6f}, p={sp['p_corr']:.6f})",
        f"- gain: {sm['rmse_gain_ratio']:.6f} (q95={sp['q95_rmse_gain']:.6f}, p={sp['p_rmse_gain']:.6f})",
        "",
        "## Flags",
    ]
    for k, v in flags.items():
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

    print(f"[QW-2046] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2046] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2046] verdict={verdict} pass_count={pass_count}/{total_flags}")


if __name__ == "__main__":
    main()
