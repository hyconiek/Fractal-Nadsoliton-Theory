#!/usr/bin/env python3
"""
QW-1989: Constrained-adversarial hard audit.
Adversarial null is restricted by sign-complexity (flip budget) and exact class balance.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

from QW_1983_FOLD_ROBUST_OPERATOR_SEARCH import ROOT, bootstrap_pass_rate
from QW_1986_TRI_BASIS_STRICT_5OF5_ATTEMPT import build_fold_channels_tri, fold_null_stats_tri
from QW_1988_TRI_BASIS_HARD_LOCK_STRESS import bootstrap_pass_rate_block


IN_QW1970 = ROOT / "report_qw1970_structural_gw_control_term_gate.json"
IN_QW1969 = ROOT / "report_qw1969_bootstrap_robust_recenter_search.json"
IN_QW1987 = ROOT / "report_qw1987_tri_basis_fold2_targeted_strict_push.json"
OUT_JSON = ROOT / "report_qw1989_constrained_adversarial_audit.json"
OUT_MD = ROOT / "RAPORT_QW1989_CONSTRAINED_ADVERSARIAL_AUDIT.md"

REAL_IID_BOOT = 3000
REAL_BLOCK_BOOT = 1800
REAL_BLOCK_LEN = 10

NULL_RANDOM_TRIALS = 32
NULL_RANDOM_BOOT = 120

ADV_FLIP_BUDGETS = [4, 6, 8]
ADV_BOOT = 500


def optimal_sign_sequence_with_flip_budget(
    t: np.ndarray,
    n_plus: int,
    max_flips: int,
) -> np.ndarray:
    """
    Maximize sum_i sign_i * t_i
    with sign_i in {-1,+1}, exact n_plus positives, and <= max_flips sign changes.
    """
    n = len(t)
    neg_inf = -1e100

    # dp[p, f, s] where s=0 => -1, s=1 => +1
    dp_prev = np.full((n_plus + 1, max_flips + 1, 2), neg_inf, dtype=float)
    back_p = np.full((n, n_plus + 1, max_flips + 1, 2), -1, dtype=np.int16)
    back_f = np.full((n, n_plus + 1, max_flips + 1, 2), -1, dtype=np.int16)
    back_s = np.full((n, n_plus + 1, max_flips + 1, 2), -1, dtype=np.int8)

    # i=0 initialization
    # sign -1
    dp_prev[0, 0, 0] = -float(t[0])
    # sign +1
    if n_plus >= 1:
        dp_prev[1, 0, 1] = float(t[0])

    for i in range(1, n):
        dp_curr = np.full((n_plus + 1, max_flips + 1, 2), neg_inf, dtype=float)
        p_max_prev = min(i, n_plus)
        for p_prev in range(p_max_prev + 1):
            for f_prev in range(max_flips + 1):
                for s_prev in (0, 1):
                    v = dp_prev[p_prev, f_prev, s_prev]
                    if v <= neg_inf / 2:
                        continue
                    for s_new in (0, 1):
                        p_new = p_prev + (1 if s_new == 1 else 0)
                        if p_new > n_plus:
                            continue
                        f_new = f_prev + (1 if s_new != s_prev else 0)
                        if f_new > max_flips:
                            continue
                        contrib = float(t[i]) if s_new == 1 else -float(t[i])
                        v_new = v + contrib
                        if v_new > dp_curr[p_new, f_new, s_new]:
                            dp_curr[p_new, f_new, s_new] = v_new
                            back_p[i, p_new, f_new, s_new] = p_prev
                            back_f[i, p_new, f_new, s_new] = f_prev
                            back_s[i, p_new, f_new, s_new] = s_prev
        dp_prev = dp_curr

    # choose best terminal state with exact plus count
    best_v = neg_inf
    best_state = (-1, -1)
    for f in range(max_flips + 1):
        for s in (0, 1):
            v = dp_prev[n_plus, f, s]
            if v > best_v:
                best_v = v
                best_state = (f, s)
    f_cur, s_cur = best_state
    if f_cur < 0:
        raise RuntimeError("No feasible constrained-adversarial sign sequence found")

    signs = np.empty(n, dtype=float)
    p_cur = n_plus
    for i in range(n - 1, -1, -1):
        signs[i] = 1.0 if s_cur == 1 else -1.0
        if i == 0:
            break
        p_prev = int(back_p[i, p_cur, f_cur, s_cur])
        f_prev = int(back_f[i, p_cur, f_cur, s_cur])
        s_prev = int(back_s[i, p_cur, f_cur, s_cur])
        p_cur, f_cur, s_cur = p_prev, f_prev, s_prev
    return signs


def constrained_adversarial_rate(
    s_base: np.ndarray,
    pairs: np.ndarray,
    c1: np.ndarray,
    c3: np.ndarray,
    c4: np.ndarray,
    xi1: float,
    xi3: float,
    xi4: float,
    thr: Dict[str, float],
    control_order: np.ndarray,
    max_flips: int,
    seed: int,
) -> float:
    ctrl_idx = np.where(pairs != 0)[0]
    n_plus = int(np.sum(pairs == 1))

    t = xi1 * c1 + xi3 * c3 + xi4 * c4
    ordered_ctrl = ctrl_idx[control_order]
    t_ord = t[ordered_ctrl]
    signs_ord = optimal_sign_sequence_with_flip_budget(t_ord, n_plus=n_plus, max_flips=max_flips)

    rand_sign = np.zeros_like(pairs, dtype=float)
    rand_sign[ordered_ctrl] = signs_ord

    c1n_raw = rand_sign * c1
    c3n_raw = rand_sign * c3
    c4n_raw = rand_sign * c4
    c1n = c1n_raw / max(float(np.std(c1n_raw)), 1e-12)
    c3n = c3n_raw / max(float(np.std(c3n_raw)), 1e-12)
    c4n = c4n_raw / max(float(np.std(c4n_raw)), 1e-12)
    s = s_base + xi1 * c1n + xi3 * c3n + xi4 * c4n
    s_hl, s_hv, s_lv = s[pairs == 0], s[pairs == 1], s[pairs == 2]
    return float(bootstrap_pass_rate(s_hl, s_hv, s_lv, thr, n_boot=ADV_BOOT, seed=seed))


def main() -> None:
    r1970 = json.loads(IN_QW1970.read_text(encoding="utf-8"))
    r1969 = json.loads(IN_QW1969.read_text(encoding="utf-8"))
    r1987 = json.loads(IN_QW1987.read_text(encoding="utf-8"))
    kernel = r1970["fixed_components"]["kernel"]
    params = r1970["fixed_components"]["params"]
    thr = r1969["thresholds"]
    xi1 = float(r1987["best"]["xi1"])
    xi3 = float(r1987["best"]["xi3"])
    xi4 = float(r1987["best"]["xi4"])

    df = pd.read_csv(ROOT / "gw1831_window_features.csv")
    df = df.copy()
    df["fold"] = df["window_idx"].astype(int) % 5
    fold_dfs = [df[df["fold"] == k].reset_index(drop=True) for k in range(5)]

    fold_rows = []
    for f, dff in enumerate(fold_dfs):
        s_hl, s_hv, s_lv, pairs, c1, c3, c4 = build_fold_channels_tri(dff, kernel, params, xi1, xi3, xi4)

        real_iid = bootstrap_pass_rate(s_hl, s_hv, s_lv, thr, n_boot=REAL_IID_BOOT, seed=191000 + f)
        real_block = bootstrap_pass_rate_block(
            s_hl=s_hl,
            s_hv=s_hv,
            s_lv=s_lv,
            thr=thr,
            n_boot=REAL_BLOCK_BOOT,
            seed=192000 + f,
            block_len=REAL_BLOCK_LEN,
        )

        pair_map = {"H1-L1": 0, "H1-V1": 1, "L1-V1": 2}
        pairs_vec = dff["pair"].map(pair_map).to_numpy(dtype=int)
        s_full = np.zeros(len(pairs_vec), dtype=float)
        s_full[pairs_vec == 0] = s_hl
        s_full[pairs_vec == 1] = s_hv
        s_full[pairs_vec == 2] = s_lv
        s_base = s_full - xi1 * c1 - xi3 * c3 - xi4 * c4

        null_mean, null_p90 = fold_null_stats_tri(
            s_base=s_base,
            pairs=pairs_vec,
            c1=c1,
            c3=c3,
            c4=c4,
            xi1=xi1,
            xi3=xi3,
            xi4=xi4,
            thr=thr,
            seed=193000 + f,
            n_trials=NULL_RANDOM_TRIALS,
            n_boot=NULL_RANDOM_BOOT,
        )

        ctrl_idx = np.where(pairs_vec != 0)[0]
        # order controls by window index for temporal complexity constraint
        control_order = np.argsort(dff.iloc[ctrl_idx]["window_idx"].to_numpy(dtype=int))

        adv_by_budget = {}
        for k in ADV_FLIP_BUDGETS:
            adv_rate = constrained_adversarial_rate(
                s_base=s_base,
                pairs=pairs_vec,
                c1=c1,
                c3=c3,
                c4=c4,
                xi1=xi1,
                xi3=xi3,
                xi4=xi4,
                thr=thr,
                control_order=control_order,
                max_flips=k,
                seed=194000 + f * 50 + k,
            )
            adv_by_budget[str(k)] = adv_rate

        adv_constrained_worst = float(max(adv_by_budget.values()))
        fold_rows.append(
            {
                "fold": f,
                "real_iid": real_iid,
                "real_block": real_block,
                "null_random_mean": null_mean,
                "null_random_p90": null_p90,
                "adv_constrained_by_flip_budget": adv_by_budget,
                "adv_constrained_worst": adv_constrained_worst,
            }
        )
        print(f"[QW-1989] fold {f} done", flush=True)

    agg = {
        "real_iid_min": float(min(r["real_iid"] for r in fold_rows)),
        "real_block_min": float(min(r["real_block"] for r in fold_rows)),
        "null_random_p90_max": float(max(r["null_random_p90"] for r in fold_rows)),
        "adv_constrained_worst_max": float(max(r["adv_constrained_worst"] for r in fold_rows)),
    }

    strict_gate = bool(
        agg["real_iid_min"] >= 0.95
        and agg["real_block_min"] >= 0.90
        and agg["null_random_p90_max"] <= 0.40
        and agg["adv_constrained_worst_max"] <= 0.45
    )
    pragmatic_gate = bool(
        agg["real_iid_min"] >= 0.94
        and agg["real_block_min"] >= 0.90
        and agg["null_random_p90_max"] <= 0.41
        and agg["adv_constrained_worst_max"] <= 0.60
    )
    verdict = (
        "CONSTRAINED_ADV_AUDIT_STRICT_PASS"
        if strict_gate
        else "CONSTRAINED_ADV_AUDIT_PRAGMATIC_PASS"
        if pragmatic_gate
        else "CONSTRAINED_ADV_AUDIT_FAIL"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source_reports": [IN_QW1970.name, IN_QW1969.name, IN_QW1987.name],
        "candidate": {"xi1": xi1, "xi3": xi3, "xi4": xi4},
        "config": {
            "real_iid_boot": REAL_IID_BOOT,
            "real_block_boot": REAL_BLOCK_BOOT,
            "real_block_len": REAL_BLOCK_LEN,
            "null_random_trials": NULL_RANDOM_TRIALS,
            "null_random_boot": NULL_RANDOM_BOOT,
            "adv_flip_budgets": ADV_FLIP_BUDGETS,
            "adv_boot": ADV_BOOT,
        },
        "fold_results": fold_rows,
        "aggregate": agg,
        "strict_gate": strict_gate,
        "pragmatic_gate": pragmatic_gate,
        "verdict": verdict,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1989: CONSTRAINED ADVERSARIAL AUDIT",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        "",
        "## Aggregate",
        f"- real_iid_min: {100.0 * agg['real_iid_min']:.2f}%",
        f"- real_block_min: {100.0 * agg['real_block_min']:.2f}%",
        f"- null_random_p90_max: {100.0 * agg['null_random_p90_max']:.2f}%",
        f"- adv_constrained_worst_max: {100.0 * agg['adv_constrained_worst_max']:.2f}%",
        "",
        "## Gates",
        f"- strict_gate: {strict_gate}",
        f"- pragmatic_gate: {pragmatic_gate}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1989] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1989] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1989] verdict={verdict}")


if __name__ == "__main__":
    main()

