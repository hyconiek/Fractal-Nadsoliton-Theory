#!/usr/bin/env python3
"""
QW-1853: Joint external confirmatory execution (PTA_V2 + GW).

Execution is blocked unless QW-1852 reports a validated external dataset.
"""

from __future__ import annotations

import importlib.util
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd
from scipy.stats import beta


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1853_joint_external_confirmatory_v2.json"
OUT_MD = ROOT / "RAPORT_QW1853_JOINT_EXTERNAL_CONFIRMATORY_V2.md"


def load_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def pinball_single(y: float, qhat: float, tau: float) -> float:
    e = y - qhat
    return float(tau * e if e >= 0.0 else (tau - 1.0) * e)


def quantile_score_single(y: float, mu: float, resid_train: np.ndarray, taus: np.ndarray) -> float:
    vals = []
    for tau in taus:
        q_res = float(np.quantile(resid_train, tau))
        vals.append(pinball_single(y, mu + q_res, float(tau)))
    return float(np.mean(vals))


def auc_pos_gt_neg(pos: np.ndarray, neg: np.ndarray) -> float:
    y = np.concatenate([np.ones(len(pos), dtype=int), np.zeros(len(neg), dtype=int)])
    s = np.concatenate([pos, neg])
    n1 = len(pos)
    n0 = len(neg)
    order = np.argsort(s)
    ranks = np.empty_like(order, dtype=float)
    ranks[order] = np.arange(1, len(s) + 1, dtype=float)
    rs = float(np.sum(ranks[y == 1]))
    return float((rs - n1 * (n1 + 1) / 2.0) / (n1 * n0))


def score_raw(df: pd.DataFrame) -> np.ndarray:
    return (
        0.55 * df["max_abs_corr"].to_numpy(dtype=float)
        + 0.25 * df["mean_abs_corr"].to_numpy(dtype=float)
        + 0.10 * df["corr_at_0ms"].to_numpy(dtype=float)
        + 0.10 * df["corr_at_10ms"].to_numpy(dtype=float)
    )


def metrics_from_scores(pair: np.ndarray, s: np.ndarray) -> Dict[str, float]:
    sh = s[pair == "H1-L1"]
    hv = s[pair == "H1-V1"]
    lv = s[pair == "L1-V1"]
    ctrl = np.concatenate([hv, lv])

    auc = auc_pos_gt_neg(sh, ctrl)
    q90 = float(np.quantile(ctrl, 0.90))
    p_shared = float(np.mean(sh > q90))
    p_ctrl = float(np.mean(ctrl > q90))
    adv = float(p_shared - p_ctrl)

    return {
        "auc_h1l1_vs_ctrl": float(auc),
        "adv_shared_minus_ctrl_q90": float(adv),
        "control_gap": float(abs(np.median(hv) - np.median(lv))),
    }


def eval_pta_v2(df_pta: pd.DataFrame, thresholds: Dict[str, float]) -> Dict:
    m1817 = load_module(ROOT / "QW_1817_SEQUENCE_OOS_VALIDATION.py", "qw1817_oos_1853")
    d1773 = load_json("report_qw1773_omega_suppressed_legacy_projection.json")

    q_center = float(d1773["projection"]["p"])

    req = [
        "theta_deg",
        "hxy",
        "f_mean",
        "f_std",
        "f_slope",
        "f_quad",
        "f_spread",
        "f_autoc1",
        "f_switch",
    ]
    miss = [c for c in req if c not in df_pta.columns]
    if miss:
        raise RuntimeError(f"PTA missing columns: {miss}")

    theta = df_pta["theta_deg"].to_numpy(dtype=float)
    y = df_pta["hxy"].to_numpy(dtype=float)

    helper = load_module(ROOT / "QW_1787_COVERAGE_FRACTION_SWEEP.py", "qw1787_helper_1853")
    hd = np.clip(helper.hellings_downs(theta), 1e-9, None)
    E_raw = np.column_stack([df_pta[c].to_numpy(dtype=float) for c in req[2:]])

    q_grid = np.linspace(max(0.9, q_center - 0.25), min(2.3, q_center + 0.25), 17)
    alpha_grid = np.array([0.08, 0.15, 0.30, 0.60, 1.20, 2.40], dtype=float)
    taus = np.array([0.10, 0.50, 0.90], dtype=float)

    n_rep = 60
    rng = np.random.default_rng(21931)

    pair_gains: List[List[float]] = [[] for _ in range(len(df_pta))]
    rep_rows = []

    for rep in range(n_rep):
        tr_idx, te_idx = m1817.stratified_split_indices(theta, train_frac=0.72, rng=rng, n_bins=8)
        if len(tr_idx) < 50 or len(te_idx) < 20:
            continue

        hd_tr, hd_te = hd[tr_idx], hd[te_idx]
        y_tr, y_te = y[tr_idx], y[te_idx]
        E_tr_raw, E_te_raw = E_raw[tr_idx], E_raw[te_idx]
        E_tr, E_te = m1817.standardize_train_test(E_tr_raw, E_te_raw)

        q2, beta2, _ = m1817.best_q_m2(hd_tr, y_tr, q_grid=q_grid)
        X2_tr = np.column_stack([hd_tr ** q2, np.ones_like(hd_tr)])
        X2_te = np.column_stack([hd_te ** q2, np.ones_like(hd_te)])
        mu2_tr = X2_tr @ beta2
        mu2_te = X2_te @ beta2
        r2_tr = y_tr - mu2_tr

        qE, alphaE, betaE, _ = m1817.best_q_alpha_m2e(
            hd_tr,
            y_tr,
            E_tr,
            q_grid=q_grid,
            alpha_grid=alpha_grid,
            rng=np.random.default_rng(61000 + rep),
        )
        XE_tr = np.column_stack([hd_tr ** qE, E_tr, np.ones_like(hd_tr)])
        XE_te = np.column_stack([hd_te ** qE, E_te, np.ones_like(hd_te)])
        muE_tr = XE_tr @ betaE
        muE_te = XE_te @ betaE
        rE_tr = y_tr - muE_tr

        gains = []
        for j, idx in enumerate(te_idx):
            qs2 = quantile_score_single(float(y_te[j]), float(mu2_te[j]), r2_tr, taus=taus)
            qsE = quantile_score_single(float(y_te[j]), float(muE_te[j]), rE_tr, taus=taus)
            g = float(qs2 - qsE)
            gains.append(g)
            pair_gains[int(idx)].append(g)

        rep_rows.append(
            {
                "rep": int(rep),
                "n_test": int(len(te_idx)),
                "q_m2": float(q2),
                "q_m2e": float(qE),
                "alpha_m2e": float(alphaE),
                "mean_pair_gain": float(np.mean(gains)),
                "prob_pair_gain_positive": float(np.mean(np.array(gains) > 0.0)),
            }
        )

    if len(rep_rows) < 50:
        raise RuntimeError(f"PTA too few valid reps: {len(rep_rows)}")

    pair_mean_gain = np.array([float(np.mean(v)) for v in pair_gains if len(v) > 0], dtype=float)
    n_pair = int(len(pair_mean_gain))
    k_pos = int(np.sum(pair_mean_gain > 0.0))

    mean_pair = float(np.mean(pair_mean_gain))
    prob_pair = float(k_pos / n_pair)
    lower95_prob = float(beta.ppf(0.05, k_pos, n_pair - k_pos + 1)) if k_pos > 0 else 0.0

    rng_bs = np.random.default_rng(21932)
    idx = rng_bs.integers(0, n_pair, size=(30000, n_pair))
    boot_means = np.mean(pair_mean_gain[idx], axis=1)
    lower95_mean = float(np.quantile(boot_means, 0.05))

    pass_flags = {
        "mean_pair_mean_gain": mean_pair >= float(thresholds["mean_pair_mean_gain_min"]),
        "bootstrap_lower95_mean_pair_mean_gain": lower95_mean >= float(thresholds["bootstrap_lower95_mean_pair_mean_gain_min"]),
        "prob_pair_mean_gain_positive": prob_pair >= float(thresholds["prob_pair_mean_gain_positive_min"]),
        "one_sided_lower95_prob_pair_mean_gain_positive": lower95_prob >= float(thresholds["one_sided_lower95_prob_pair_mean_gain_positive_min"]),
    }

    return {
        "summary": {
            "n_pairs": n_pair,
            "n_replications": int(len(rep_rows)),
            "mean_pair_mean_gain": mean_pair,
            "bootstrap_lower95_mean_pair_mean_gain": lower95_mean,
            "prob_pair_mean_gain_positive": prob_pair,
            "one_sided_lower95_prob_pair_mean_gain_positive": lower95_prob,
        },
        "pass_flags": pass_flags,
        "replications": rep_rows,
    }


def eval_gw(df_gw: pd.DataFrame, thresholds: Dict[str, float]) -> Dict:
    req = ["pair", "window_idx", "max_abs_corr", "mean_abs_corr", "corr_at_0ms", "corr_at_10ms"]
    miss = [c for c in req if c not in df_gw.columns]
    if miss:
        raise RuntimeError(f"GW missing columns: {miss}")

    d = df_gw.dropna(subset=req).copy()
    d["score_raw"] = score_raw(d)

    fold_rows = []
    n_folds = 5
    for fold in range(n_folds):
        tr = d[(d["window_idx"].astype(int) % n_folds) != fold].copy()
        te = d[(d["window_idx"].astype(int) % n_folds) == fold].copy()

        med_hv = float(np.median(tr.loc[tr["pair"] == "H1-V1", "score_raw"]))
        med_lv = float(np.median(tr.loc[tr["pair"] == "L1-V1", "score_raw"]))
        med_ctrl = float(np.median(tr.loc[tr["pair"] != "H1-L1", "score_raw"]))

        off_hv = med_hv - med_ctrl
        off_lv = med_lv - med_ctrl

        te["score_cal"] = te["score_raw"].to_numpy(dtype=float)
        te.loc[te["pair"] == "H1-V1", "score_cal"] = te.loc[te["pair"] == "H1-V1", "score_cal"] - off_hv
        te.loc[te["pair"] == "L1-V1", "score_cal"] = te.loc[te["pair"] == "L1-V1", "score_cal"] - off_lv

        m = metrics_from_scores(te["pair"].to_numpy(), te["score_cal"].to_numpy(dtype=float))
        fold_rows.append(
            {
                "fold": int(fold),
                "n_test": int(len(te)),
                "cal_auc": float(m["auc_h1l1_vs_ctrl"]),
                "cal_adv": float(m["adv_shared_minus_ctrl_q90"]),
                "cal_control_gap": float(m["control_gap"]),
            }
        )

    auc = np.array([r["cal_auc"] for r in fold_rows], dtype=float)
    adv = np.array([r["cal_adv"] for r in fold_rows], dtype=float)
    gap = np.array([r["cal_control_gap"] for r in fold_rows], dtype=float)

    summary = {
        "n_rows": int(len(d)),
        "n_folds": n_folds,
        "calibrated_mean_auc": float(np.mean(auc)),
        "calibrated_mean_adv": float(np.mean(adv)),
        "calibrated_mean_control_gap": float(np.mean(gap)),
        "calibrated_prob_adv_positive": float(np.mean(adv > 0.0)),
    }

    pass_flags = {
        "calibrated_mean_auc": summary["calibrated_mean_auc"] >= float(thresholds["calibrated_mean_auc_min"]),
        "calibrated_mean_adv": summary["calibrated_mean_adv"] >= float(thresholds["calibrated_mean_adv_min"]),
        "calibrated_mean_control_gap": summary["calibrated_mean_control_gap"] <= float(thresholds["calibrated_mean_control_gap_max"]),
        "calibrated_prob_adv_positive": summary["calibrated_prob_adv_positive"] >= float(thresholds["calibrated_prob_adv_positive_min"]),
    }

    return {
        "summary": summary,
        "pass_flags": pass_flags,
        "folds": fold_rows,
    }


def main() -> None:
    d1839 = load_json("report_qw1839_joint_confirmatory_prereg_protocol.json")
    d1850 = load_json("report_qw1850_pta_v2_prereg_protocol.json")
    d1852 = load_json("report_qw1852_external_confirmatory_data_precheck.json")

    precheck_ready = d1852.get("readiness") == "EXTERNAL_DATASET_READY_FOR_QW1853"
    candidate_dir = Path(d1852.get("validation", {}).get("candidate_dir", ""))

    if not precheck_ready:
        out = {
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "precheck_readiness": d1852.get("readiness"),
            "hard_gate": "PARTIAL",
            "readiness": "BLOCKED_WAITING_VALID_EXTERNAL_DATASET",
            "verdict": "JOINT_EXTERNAL_CONFIRMATORY_V2_NOT_EXECUTED",
        }
        OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

        lines = [
            "# RAPORT QW-1853: JOINT EXTERNAL CONFIRMATORY V2",
            "",
            f"- Data UTC: {out['generated_utc']}",
            f"- Hard gate: **{out['hard_gate']}**",
            f"- Readiness: **{out['readiness']}**",
            f"- Precheck readiness: `{out['precheck_readiness']}`",
            "",
            "## Status",
            "- Execution blocked until QW-1852 validates external candidate dataset.",
            "",
            "## Artifacts",
            f"- JSON: `{OUT_JSON.name}`",
        ]
        OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

        print(f"[QW-1853] Saved JSON: {OUT_JSON.name}")
        print(f"[QW-1853] Saved MD:   {OUT_MD.name}")
        return

    manifest = json.loads((candidate_dir / "manifest.json").read_text(encoding="utf-8"))
    role_map = {x.get("role"): x for x in manifest.get("files", []) if isinstance(x, dict)}

    pta_path = candidate_dir / str(role_map["pta_pairs"]["path"])
    gw_path = candidate_dir / str(role_map["gw_windows"]["path"])

    df_pta = pd.read_csv(pta_path)
    df_gw = pd.read_csv(gw_path)

    pta_thresholds = d1850["protocol"]["pta_v2_protocol"]["thresholds"]
    gw_thresholds = d1839["protocol"]["gw_protocol"]["thresholds"]

    pta_eval = eval_pta_v2(df_pta, pta_thresholds)
    gw_eval = eval_gw(df_gw, gw_thresholds)

    pass_pta = all(bool(v) for v in pta_eval["pass_flags"].values())
    pass_gw = all(bool(v) for v in gw_eval["pass_flags"].values())

    if pass_pta and pass_gw:
        hard_gate = "PASS"
        readiness = "JOINT_EXTERNAL_CONFIRMATORY_V2_PASS"
    elif pass_pta or pass_gw:
        hard_gate = "PARTIAL"
        readiness = "JOINT_EXTERNAL_CONFIRMATORY_V2_PARTIAL"
    else:
        hard_gate = "FAIL"
        readiness = "JOINT_EXTERNAL_CONFIRMATORY_V2_FAIL"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "dataset": {
            "candidate_dir": str(candidate_dir),
            "dataset_id": manifest.get("dataset", {}).get("dataset_id"),
        },
        "pta_v2": pta_eval,
        "gw": gw_eval,
        "joint_pass_flags": {
            "pta_v2_all": bool(pass_pta),
            "gw_all": bool(pass_gw),
        },
        "hard_gate": hard_gate,
        "readiness": readiness,
        "verdict": "JOINT_EXTERNAL_CONFIRMATORY_V2_EXECUTED",
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    ps = pta_eval["summary"]
    gs = gw_eval["summary"]
    lines = [
        "# RAPORT QW-1853: JOINT EXTERNAL CONFIRMATORY V2",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Dataset ID: {out['dataset']['dataset_id']}",
        f"- Hard gate: **{hard_gate}**",
        f"- Readiness: **{readiness}**",
        "",
        "## PTA V2",
        f"- mean_pair_mean_gain: {ps['mean_pair_mean_gain']:.6f}",
        f"- bootstrap_lower95_mean_pair_mean_gain: {ps['bootstrap_lower95_mean_pair_mean_gain']:.6f}",
        f"- prob_pair_mean_gain_positive: {ps['prob_pair_mean_gain_positive']:.3f}",
        f"- one_sided_lower95_prob_pair_mean_gain_positive: {ps['one_sided_lower95_prob_pair_mean_gain_positive']:.3f}",
        f"- all thresholds pass: {pass_pta}",
        "",
        "## GW",
        f"- calibrated_mean_auc: {gs['calibrated_mean_auc']:.6f}",
        f"- calibrated_mean_adv: {gs['calibrated_mean_adv']:.6f}",
        f"- calibrated_mean_control_gap: {gs['calibrated_mean_control_gap']:.6f}",
        f"- calibrated_prob_adv_positive: {gs['calibrated_prob_adv_positive']:.3f}",
        f"- all thresholds pass: {pass_gw}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1853] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1853] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
