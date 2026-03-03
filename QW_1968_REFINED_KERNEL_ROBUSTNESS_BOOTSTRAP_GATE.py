#!/usr/bin/env python3
"""
QW-1968: Robustness gate for the single-pass point found in QW-1967.

Purpose:
- test whether the QW-1967 triad pass has finite local stability volume,
- quantify GW-channel statistical robustness via bootstrap,
- keep one frozen kernel, fixed mass branch, and fixed q_nu ordering.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
IN_QW1967 = ROOT / "report_qw1967_isospin_split_local_refinement_gate.json"
OUT_JSON = ROOT / "report_qw1968_refined_kernel_robustness_bootstrap_gate.json"
OUT_MD = ROOT / "RAPORT_QW1968_REFINED_KERNEL_ROBUSTNESS_BOOTSTRAP_GATE.md"


CKM_REF = np.array(
    [
        [0.97401, 0.22650, 0.00361],
        [0.22636, 0.97320, 0.04053],
        [0.00854, 0.03978, 0.99917],
    ],
    dtype=float,
)
PMNS_REF = np.array(
    [
        [0.821, 0.550, 0.150],
        [0.432, 0.582, 0.693],
        [0.378, 0.598, 0.707],
    ],
    dtype=float,
)


def kernel_fn(d: np.ndarray, omega: float, phi: float, beta: float, eta: float) -> np.ndarray:
    return np.cos(omega * d + phi) / (1.0 + beta * (d**eta))


def cyclic_distance_matrix(q_left: np.ndarray, q_right: np.ndarray, modulus: int = 24) -> np.ndarray:
    dq = np.abs(q_left[:, None] - q_right[None, :]) % float(modulus)
    return np.minimum(dq, float(modulus) - dq)


def rank_auc_pos_gt_neg(pos: np.ndarray, neg: np.ndarray) -> float:
    y = np.concatenate([np.ones(len(pos), dtype=int), np.zeros(len(neg), dtype=int)])
    s = np.concatenate([pos, neg])
    n1 = len(pos)
    n0 = len(neg)
    order = np.argsort(s)
    ranks = np.empty_like(order, dtype=float)
    ranks[order] = np.arange(1, len(s) + 1, dtype=float)
    rs = float(np.sum(ranks[y == 1]))
    return float((rs - n1 * (n1 + 1) / 2.0) / (n1 * n0))


def flavor_hamiltonian(
    q: np.ndarray,
    iso_tag: float,
    sector_tag: float,
    params: Dict[str, float],
    kernel: Dict[str, float],
) -> np.ndarray:
    n = len(q)
    d = 1.0 + cyclic_distance_matrix(q, q, modulus=24)
    k = kernel_fn(d, kernel["omega"], kernel["phi"], kernel["beta"], kernel["eta"])
    base = np.sign(k) * (np.abs(k) ** params["p_amp"]) * (d**params["r_dist"])

    idx = np.arange(n, dtype=float)
    i_minus_j = idx[:, None] - idx[None, :]
    gap = np.abs(i_minus_j)
    q_diff = q[:, None] - q[None, :]
    near = np.exp(-params["rho_gap"] * gap)

    phase = (
        params["phase_q"] * q_diff
        + iso_tag * params["theta_iso"] * i_minus_j
        + sector_tag * params["theta_sector"] * i_minus_j
    )
    strength = params["lambda_mix"] * base * near

    re_raw = strength * np.cos(phase)
    im_raw = strength * np.sin(phase)

    re = 0.5 * (re_raw + re_raw.T)
    im_asym = 0.5 * (im_raw - im_raw.T)
    np.fill_diagonal(re, 0.0)
    np.fill_diagonal(im_asym, 0.0)

    idx_centered = idx - np.mean(idx)
    q_centered = q - np.mean(q)
    diag = (
        params["diag_q_coeff"] * q_centered
        + iso_tag * params["diag_iso"] * idx_centered
        + sector_tag * params["diag_sector"] * idx_centered
    )

    h = re + 1j * params["chi_im"] * im_asym + np.diag(diag)
    return 0.5 * (h + h.conj().T)


def flavor_metrics(
    q_up: np.ndarray,
    q_down: np.ndarray,
    q_nu: np.ndarray,
    q_lep: np.ndarray,
    params: Dict[str, float],
    kernel: Dict[str, float],
) -> Dict[str, float]:
    hu = flavor_hamiltonian(q_up, iso_tag=+1.0, sector_tag=+1.0, params=params, kernel=kernel)
    hd = flavor_hamiltonian(q_down, iso_tag=-1.0, sector_tag=+1.0, params=params, kernel=kernel)
    hn = flavor_hamiltonian(q_nu, iso_tag=+1.0, sector_tag=-1.0, params=params, kernel=kernel)
    hl = flavor_hamiltonian(q_lep, iso_tag=-1.0, sector_tag=-1.0, params=params, kernel=kernel)

    _, uu = np.linalg.eigh(hu)
    _, ud = np.linalg.eigh(hd)
    _, un = np.linalg.eigh(hn)
    _, ul = np.linalg.eigh(hl)

    ckm = np.abs(uu.conj().T @ ud)
    pmns = np.abs(un.conj().T @ ul)
    ckm_rel = np.abs(ckm - CKM_REF) / np.clip(CKM_REF, 1e-12, None)
    pmns_rel = np.abs(pmns - PMNS_REF) / np.clip(PMNS_REF, 1e-12, None)
    return {
        "ckm_mean_rel_pct": float(100.0 * np.mean(ckm_rel)),
        "pmns_mean_rel_pct": float(100.0 * np.mean(pmns_rel)),
    }


def gw_metrics(p_amp: float, r_dist: float, kernel: Dict[str, float], df_gw: pd.DataFrame) -> Dict[str, float]:
    d = np.array([1.0, 2.0, 3.0, 4.0], dtype=float)
    raw_w = (np.abs(kernel_fn(d, kernel["omega"], kernel["phi"], kernel["beta"], kernel["eta"])) ** p_amp) * (d**r_dist)
    w = raw_w / np.sum(raw_w)

    score = (
        w[0] * df_gw["max_abs_corr"].to_numpy(dtype=float)
        + w[1] * df_gw["mean_abs_corr"].to_numpy(dtype=float)
        + w[2] * df_gw["corr_at_0ms"].to_numpy(dtype=float)
        + w[3] * df_gw["corr_at_10ms"].to_numpy(dtype=float)
    )
    pair = df_gw["pair"].astype(str).to_numpy()
    s_hl = score[pair == "H1-L1"]
    s_hv = score[pair == "H1-V1"]
    s_lv = score[pair == "L1-V1"]
    s_ctrl = np.concatenate([s_hv, s_lv])

    q90 = float(np.quantile(s_ctrl, 0.90))
    auc = rank_auc_pos_gt_neg(s_hl, s_ctrl)
    adv = float(np.mean(s_hl > q90) - np.mean(s_ctrl > q90))
    sep = float(np.median(s_hl) - np.median(s_ctrl))
    ctrl_gap = float(abs(np.median(s_hv) - np.median(s_lv)))
    return {
        "auc_h1l1_vs_ctrl": auc,
        "adv_shared_minus_ctrl_q90": adv,
        "sep_median_h1l1_minus_ctrl": sep,
        "control_median_gap": ctrl_gap,
    }


def eval_flags(
    mass_flags: Dict[str, bool],
    flavor: Dict[str, float],
    gw: Dict[str, float],
    thresholds: Dict[str, float],
) -> Dict[str, bool]:
    return {
        **mass_flags,
        "ckm_mean_rel_pct_le_max": bool(flavor["ckm_mean_rel_pct"] <= thresholds["ckm_mean_rel_pct_max"]),
        "pmns_mean_rel_pct_le_max": bool(flavor["pmns_mean_rel_pct"] <= thresholds["pmns_mean_rel_pct_max"]),
        "gw_sep_ge_min": bool(gw["sep_median_h1l1_minus_ctrl"] >= thresholds["gw_sep_min"]),
        "gw_adv_ge_min": bool(gw["adv_shared_minus_ctrl_q90"] >= thresholds["gw_adv_min"]),
        "gw_auc_ge_min": bool(gw["auc_h1l1_vs_ctrl"] >= thresholds["gw_auc_min"]),
        "gw_control_gap_le_max": bool(gw["control_median_gap"] <= thresholds["gw_control_gap_max"]),
    }


def violation_score(flavor: Dict[str, float], gw: Dict[str, float], thresholds: Dict[str, float]) -> float:
    # 0.0 means full pass; positive values measure normalized violation.
    return float(
        max(0.0, flavor["ckm_mean_rel_pct"] / thresholds["ckm_mean_rel_pct_max"] - 1.0)
        + max(0.0, flavor["pmns_mean_rel_pct"] / thresholds["pmns_mean_rel_pct_max"] - 1.0)
        + 40.0 * max(0.0, thresholds["gw_auc_min"] - gw["auc_h1l1_vs_ctrl"])
        + 45.0 * max(0.0, thresholds["gw_adv_min"] - gw["adv_shared_minus_ctrl_q90"])
        + 250.0 * max(0.0, thresholds["gw_sep_min"] - gw["sep_median_h1l1_minus_ctrl"])
        + 30.0 * max(0.0, gw["control_median_gap"] - thresholds["gw_control_gap_max"])
    )


def clip_to_range(v: float, lo: float, hi: float) -> float:
    return float(min(hi, max(lo, v)))


def local_jitter_samples(
    best_params: Dict[str, float],
    ranges: Dict[str, List[float]],
    radius: float,
    n_samples: int,
    rng: np.random.Generator,
) -> List[Dict[str, float]]:
    keys = sorted(best_params.keys())
    out: List[Dict[str, float]] = []
    for _ in range(n_samples):
        p: Dict[str, float] = {}
        for k in keys:
            lo, hi = float(ranges[k][0]), float(ranges[k][1])
            span = hi - lo
            delta = rng.uniform(-radius, radius) * span
            p[k] = clip_to_range(float(best_params[k]) + float(delta), lo, hi)
        out.append(p)
    return out


def local_robustness_scan(
    best_params: Dict[str, float],
    ranges: Dict[str, List[float]],
    radii: List[float],
    n_samples_per_radius: int,
    kernel: Dict[str, float],
    q_up: np.ndarray,
    q_down: np.ndarray,
    q_nu: np.ndarray,
    q_lep: np.ndarray,
    df_gw: pd.DataFrame,
    mass_flags: Dict[str, bool],
    thresholds: Dict[str, float],
    seed: int,
) -> Tuple[List[Dict[str, object]], Dict[str, object]]:
    rng = np.random.default_rng(seed)
    rows: List[Dict[str, object]] = []
    global_best = None

    for radius in radii:
        samples = local_jitter_samples(best_params, ranges, radius, n_samples_per_radius, rng)
        pass_count = 0
        fail_breakdown = {
            "flavor_only_fail": 0,
            "gw_only_fail": 0,
            "both_flavor_gw_fail": 0,
            "mass_fail": 0,
        }
        best_row = None

        for params in samples:
            f = flavor_metrics(q_up, q_down, q_nu, q_lep, params, kernel)
            g = gw_metrics(params["p_amp"], params["r_dist"], kernel, df_gw)
            flags = eval_flags(mass_flags, f, g, thresholds)
            all_pass = bool(all(flags.values()))
            if all_pass:
                pass_count += 1
            else:
                if not (
                    flags["mass_mean_rel_pct_le_max"]
                    and flags["mass_max_rel_pct_le_max"]
                    and flags["mass_tau_charm_ratio_err_le_max"]
                ):
                    fail_breakdown["mass_fail"] += 1
                flavor_ok = flags["ckm_mean_rel_pct_le_max"] and flags["pmns_mean_rel_pct_le_max"]
                gw_ok = (
                    flags["gw_sep_ge_min"]
                    and flags["gw_adv_ge_min"]
                    and flags["gw_auc_ge_min"]
                    and flags["gw_control_gap_le_max"]
                )
                if (not flavor_ok) and gw_ok:
                    fail_breakdown["flavor_only_fail"] += 1
                elif flavor_ok and (not gw_ok):
                    fail_breakdown["gw_only_fail"] += 1
                elif (not flavor_ok) and (not gw_ok):
                    fail_breakdown["both_flavor_gw_fail"] += 1

            v = violation_score(f, g, thresholds)
            row = {
                "radius_fraction_of_range": float(radius),
                "params": params,
                "flavor": f,
                "gw": g,
                "flags": flags,
                "all_pass": all_pass,
                "violation_score": float(v),
            }
            if best_row is None or row["violation_score"] < best_row["violation_score"]:
                best_row = row
            if global_best is None or row["violation_score"] < global_best["violation_score"]:
                global_best = row

        rows.append(
            {
                "radius_fraction_of_range": float(radius),
                "n_samples": int(n_samples_per_radius),
                "pass_count": int(pass_count),
                "pass_rate": float(pass_count / n_samples_per_radius),
                "best": best_row,
                "fail_breakdown": fail_breakdown,
            }
        )
    return rows, global_best


def bootstrap_gw(
    best_params: Dict[str, float],
    kernel: Dict[str, float],
    df_gw: pd.DataFrame,
    thresholds: Dict[str, float],
    mass_flags: Dict[str, bool],
    flavor_at_best: Dict[str, float],
    n_boot: int,
    seed: int,
) -> Dict[str, object]:
    rng = np.random.default_rng(seed)
    by_pair = {
        "H1-L1": df_gw[df_gw["pair"] == "H1-L1"].reset_index(drop=True),
        "H1-V1": df_gw[df_gw["pair"] == "H1-V1"].reset_index(drop=True),
        "L1-V1": df_gw[df_gw["pair"] == "L1-V1"].reset_index(drop=True),
    }
    n_per_pair = {k: len(v) for k, v in by_pair.items()}

    gw_rows: List[Dict[str, float]] = []
    gw_pass_count = 0
    triad_pass_count = 0

    flavor_flags = {
        "ckm_mean_rel_pct_le_max": bool(flavor_at_best["ckm_mean_rel_pct"] <= thresholds["ckm_mean_rel_pct_max"]),
        "pmns_mean_rel_pct_le_max": bool(flavor_at_best["pmns_mean_rel_pct"] <= thresholds["pmns_mean_rel_pct_max"]),
    }

    for _ in range(n_boot):
        parts = []
        for pair_name, df_pair in by_pair.items():
            idx = rng.integers(0, len(df_pair), size=n_per_pair[pair_name], endpoint=False)
            parts.append(df_pair.iloc[idx].assign(pair=pair_name))
        boot_df = pd.concat(parts, axis=0, ignore_index=True)

        g = gw_metrics(best_params["p_amp"], best_params["r_dist"], kernel, boot_df)
        gw_flags = {
            "gw_sep_ge_min": bool(g["sep_median_h1l1_minus_ctrl"] >= thresholds["gw_sep_min"]),
            "gw_adv_ge_min": bool(g["adv_shared_minus_ctrl_q90"] >= thresholds["gw_adv_min"]),
            "gw_auc_ge_min": bool(g["auc_h1l1_vs_ctrl"] >= thresholds["gw_auc_min"]),
            "gw_control_gap_le_max": bool(g["control_median_gap"] <= thresholds["gw_control_gap_max"]),
        }
        gw_pass = bool(all(gw_flags.values()))
        if gw_pass:
            gw_pass_count += 1

        triad_flags = {**mass_flags, **flavor_flags, **gw_flags}
        triad_pass = bool(all(triad_flags.values()))
        if triad_pass:
            triad_pass_count += 1

        gw_rows.append(g)

    df = pd.DataFrame(gw_rows)
    ci = {}
    for col in df.columns:
        q025, q50, q975 = np.quantile(df[col].to_numpy(dtype=float), [0.025, 0.5, 0.975])
        ci[col] = {"q025": float(q025), "q50": float(q50), "q975": float(q975)}

    return {
        "seed": int(seed),
        "n_boot": int(n_boot),
        "n_per_pair": n_per_pair,
        "gw_pass_rate": float(gw_pass_count / n_boot),
        "triad_pass_rate": float(triad_pass_count / n_boot),
        "gw_metric_ci_95": ci,
    }


def one_at_a_time_sensitivity(
    best_params: Dict[str, float],
    ranges: Dict[str, List[float]],
    kernel: Dict[str, float],
    q_up: np.ndarray,
    q_down: np.ndarray,
    q_nu: np.ndarray,
    q_lep: np.ndarray,
    df_gw: pd.DataFrame,
    thresholds: Dict[str, float],
) -> List[Dict[str, object]]:
    base_f = flavor_metrics(q_up, q_down, q_nu, q_lep, best_params, kernel)
    base_g = gw_metrics(best_params["p_amp"], best_params["r_dist"], kernel, df_gw)
    base_v = violation_score(base_f, base_g, thresholds)

    rows = []
    for k in sorted(best_params.keys()):
        lo, hi = float(ranges[k][0]), float(ranges[k][1])
        step = 0.005 * (hi - lo)
        p_plus = dict(best_params)
        p_minus = dict(best_params)
        p_plus[k] = clip_to_range(best_params[k] + step, lo, hi)
        p_minus[k] = clip_to_range(best_params[k] - step, lo, hi)

        f_plus = flavor_metrics(q_up, q_down, q_nu, q_lep, p_plus, kernel)
        g_plus = gw_metrics(p_plus["p_amp"], p_plus["r_dist"], kernel, df_gw)
        v_plus = violation_score(f_plus, g_plus, thresholds)

        f_minus = flavor_metrics(q_up, q_down, q_nu, q_lep, p_minus, kernel)
        g_minus = gw_metrics(p_minus["p_amp"], p_minus["r_dist"], kernel, df_gw)
        v_minus = violation_score(f_minus, g_minus, thresholds)

        rows.append(
            {
                "param": k,
                "step_abs": float(step),
                "base_violation_score": float(base_v),
                "plus_violation_score": float(v_plus),
                "minus_violation_score": float(v_minus),
                "sensitivity_abs": float(abs(v_plus - v_minus) / (2.0 * max(step, 1e-12))),
            }
        )

    rows.sort(key=lambda x: x["sensitivity_abs"], reverse=True)
    return rows


def verdict_from_robustness(local_rows: List[Dict[str, object]], triad_boot_pass_rate: float) -> Tuple[str, str]:
    by_radius = {float(r["radius_fraction_of_range"]): float(r["pass_rate"]) for r in local_rows}
    r005 = by_radius.get(0.005, 0.0)
    r010 = by_radius.get(0.01, 0.0)
    r020 = by_radius.get(0.02, 0.0)

    if r005 >= 0.01 and r010 >= 0.005 and r020 >= 0.001 and triad_boot_pass_rate >= 0.95:
        return (
            "ROBUST_LOCKABLE_TRIAD_PASS",
            "LOCK_SHARED_OPERATOR_AND_RUN_FULLY_EXTERNAL_CONFIRMATORY",
        )
    if (max(by_radius.values()) > 0.0) and triad_boot_pass_rate >= 0.80:
        return (
            "NARROW_BUT_STATISTICALLY_SUPPORTED_PASS",
            "REGULARIZE_OPERATOR_AND_RETEST_ROBUSTNESS_BEFORE_EXTERNAL",
        )
    return (
        "FRAGILE_PASS_NOT_YET_LOCKABLE",
        "INCREASE_PASS_VOLUME_WITHOUT_KERNEL_RETUNE_THEN_REPEAT_QW1968",
    )


def main() -> None:
    r1967 = json.loads(IN_QW1967.read_text(encoding="utf-8"))
    kernel = r1967["kernel"]
    thresholds = r1967["thresholds"]
    mass_flags = r1967["source_mass_branch"]["mass_flags_fixed"]
    best = r1967["summary"]["best"]
    best_params = best["params"]
    best_q_nu = np.array(best["q_nu"], dtype=float)

    ranges = r1967["search"]["ranges"]

    r1962 = json.loads((ROOT / "report_qw1962_noncircular_branch_unified_triad_gate.json").read_text(encoding="utf-8"))
    q_assign = r1962["q_assignment_used"]
    q_up = np.array([q_assign["Top"], q_assign["Charm"], q_assign["Muon"]], dtype=float)
    q_down = np.array([q_assign["Bottom"], q_assign["Tau"], q_assign["Muon"]], dtype=float)
    q_lep = np.array([q_assign["Electron"], q_assign["Muon"], q_assign["Tau"]], dtype=float)

    df_gw = pd.read_csv(ROOT / "gw1831_window_features.csv")

    flavor_best = flavor_metrics(q_up, q_down, best_q_nu, q_lep, best_params, kernel)
    gw_best = gw_metrics(best_params["p_amp"], best_params["r_dist"], kernel, df_gw)
    best_flags = eval_flags(mass_flags, flavor_best, gw_best, thresholds)

    local_radii = [0.0025, 0.005, 0.01, 0.02, 0.05]
    n_samples_per_radius = 12000
    local_rows, local_global_best = local_robustness_scan(
        best_params=best_params,
        ranges=ranges,
        radii=local_radii,
        n_samples_per_radius=n_samples_per_radius,
        kernel=kernel,
        q_up=q_up,
        q_down=q_down,
        q_nu=best_q_nu,
        q_lep=q_lep,
        df_gw=df_gw,
        mass_flags=mass_flags,
        thresholds=thresholds,
        seed=1968,
    )

    boot = bootstrap_gw(
        best_params=best_params,
        kernel=kernel,
        df_gw=df_gw,
        thresholds=thresholds,
        mass_flags=mass_flags,
        flavor_at_best=flavor_best,
        n_boot=5000,
        seed=2968,
    )

    sensitivity = one_at_a_time_sensitivity(
        best_params=best_params,
        ranges=ranges,
        kernel=kernel,
        q_up=q_up,
        q_down=q_down,
        q_nu=best_q_nu,
        q_lep=q_lep,
        df_gw=df_gw,
        thresholds=thresholds,
    )

    verdict, required_next = verdict_from_robustness(local_rows, boot["triad_pass_rate"])

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source_report": IN_QW1967.name,
        "fixed_branch": {
            "kernel": kernel,
            "best_q_nu": [int(x) for x in best["q_nu"]],
            "best_params": best_params,
            "best_metrics": {"flavor": flavor_best, "gw": gw_best, "flags": best_flags},
            "mass_flags_fixed": mass_flags,
        },
        "thresholds": thresholds,
        "local_robustness_scan": {
            "seed": 1968,
            "radii_fraction_of_range": local_radii,
            "n_samples_per_radius": int(n_samples_per_radius),
            "results_by_radius": local_rows,
            "global_best_in_scan": local_global_best,
        },
        "bootstrap_gw_and_triad": boot,
        "one_at_a_time_sensitivity_top10": sensitivity[:10],
        "verdict": verdict,
        "required_next_step": required_next,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    radius_lines = []
    for row in local_rows:
        radius_lines.append(
            (
                f"- radius={row['radius_fraction_of_range']:.4f}: "
                f"pass_rate={100.0 * row['pass_rate']:.3f}% "
                f"({row['pass_count']}/{row['n_samples']})"
            )
        )

    sens_lines = []
    for item in sensitivity[:5]:
        sens_lines.append(
            f"- {item['param']}: sensitivity_abs={item['sensitivity_abs']:.3f}, step={item['step_abs']:.6f}"
        )

    lines = [
        "# RAPORT QW-1968: REFINED KERNEL ROBUSTNESS BOOTSTRAP GATE",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        "",
        "## Baseline (from QW-1967 best point)",
        (
            f"- flavor CKM/PMNS mean rel%: "
            f"{flavor_best['ckm_mean_rel_pct']:.3f}/{flavor_best['pmns_mean_rel_pct']:.3f}"
        ),
        (
            f"- GW auc/adv/sep/gap: "
            f"{gw_best['auc_h1l1_vs_ctrl']:.4f}/"
            f"{gw_best['adv_shared_minus_ctrl_q90']:.4f}/"
            f"{gw_best['sep_median_h1l1_minus_ctrl']:.6f}/"
            f"{gw_best['control_median_gap']:.6f}"
        ),
        f"- all_pass: {all(best_flags.values())}",
        "",
        "## Local Robustness (fixed q_nu, no sector retune)",
        *radius_lines,
        "",
        "## GW Bootstrap (n=5000)",
        f"- GW pass rate: {100.0 * boot['gw_pass_rate']:.2f}%",
        f"- Triad pass rate (mass+flavor fixed + GW bootstrap): {100.0 * boot['triad_pass_rate']:.2f}%",
        "",
        "## Sensitivity Top 5 (one-at-a-time)",
        *sens_lines,
        "",
        "## Required Next Step",
        f"- {required_next}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1968] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1968] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1968] verdict={verdict}")


if __name__ == "__main__":
    main()
