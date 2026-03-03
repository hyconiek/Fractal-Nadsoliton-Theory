#!/usr/bin/env python3
"""
QW-1969: Bootstrap-robust recentering search around QW-1967/1968 pass point.

Goal:
- keep frozen kernel and fixed mass branch,
- keep fixed q_nu (no sector retune),
- find nearby shared-operator parameters with improved GW bootstrap triad pass-rate.
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
IN_QW1968 = ROOT / "report_qw1968_refined_kernel_robustness_bootstrap_gate.json"
OUT_JSON = ROOT / "report_qw1969_bootstrap_robust_recenter_search.json"
OUT_MD = ROOT / "RAPORT_QW1969_BOOTSTRAP_ROBUST_RECENTER_SEARCH.md"


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


def gw_weights(params: Dict[str, float], kernel: Dict[str, float]) -> np.ndarray:
    d = np.array([1.0, 2.0, 3.0, 4.0], dtype=float)
    raw = (np.abs(kernel_fn(d, kernel["omega"], kernel["phi"], kernel["beta"], kernel["eta"])) ** params["p_amp"]) * (
        d**params["r_dist"]
    )
    return raw / np.sum(raw)


def gw_metrics_from_score_arrays(s_hl: np.ndarray, s_hv: np.ndarray, s_lv: np.ndarray) -> Dict[str, float]:
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


def gw_metrics(params: Dict[str, float], kernel: Dict[str, float], gw_features: np.ndarray, pairs: np.ndarray) -> Dict[str, float]:
    w = gw_weights(params, kernel)
    score = gw_features @ w
    s_hl = score[pairs == 0]
    s_hv = score[pairs == 1]
    s_lv = score[pairs == 2]
    return gw_metrics_from_score_arrays(s_hl, s_hv, s_lv)


def triad_flags(
    mass_flags: Dict[str, bool],
    flavor: Dict[str, float],
    gw: Dict[str, float],
    thr: Dict[str, float],
) -> Dict[str, bool]:
    return {
        **mass_flags,
        "ckm_mean_rel_pct_le_max": bool(flavor["ckm_mean_rel_pct"] <= thr["ckm_mean_rel_pct_max"]),
        "pmns_mean_rel_pct_le_max": bool(flavor["pmns_mean_rel_pct"] <= thr["pmns_mean_rel_pct_max"]),
        "gw_sep_ge_min": bool(gw["sep_median_h1l1_minus_ctrl"] >= thr["gw_sep_min"]),
        "gw_adv_ge_min": bool(gw["adv_shared_minus_ctrl_q90"] >= thr["gw_adv_min"]),
        "gw_auc_ge_min": bool(gw["auc_h1l1_vs_ctrl"] >= thr["gw_auc_min"]),
        "gw_control_gap_le_max": bool(gw["control_median_gap"] <= thr["gw_control_gap_max"]),
    }


def normalized_margins(flavor: Dict[str, float], gw: Dict[str, float], thr: Dict[str, float]) -> Dict[str, float]:
    return {
        "ckm_margin": (thr["ckm_mean_rel_pct_max"] - flavor["ckm_mean_rel_pct"]) / thr["ckm_mean_rel_pct_max"],
        "pmns_margin": (thr["pmns_mean_rel_pct_max"] - flavor["pmns_mean_rel_pct"]) / thr["pmns_mean_rel_pct_max"],
        "auc_margin": (gw["auc_h1l1_vs_ctrl"] - thr["gw_auc_min"]) / thr["gw_auc_min"],
        "adv_margin": (gw["adv_shared_minus_ctrl_q90"] - thr["gw_adv_min"]) / max(thr["gw_adv_min"], 1e-12),
        "sep_margin": (gw["sep_median_h1l1_minus_ctrl"] - thr["gw_sep_min"]) / max(thr["gw_sep_min"], 1e-12),
        "gap_margin": (thr["gw_control_gap_max"] - gw["control_median_gap"]) / max(thr["gw_control_gap_max"], 1e-12),
    }


def robust_margin_score(margins: Dict[str, float]) -> float:
    vals = np.array(list(margins.values()), dtype=float)
    return float(np.min(vals) + 0.25 * np.mean(vals))


def clip_to_range(v: float, lo: float, hi: float) -> float:
    return float(min(hi, max(lo, v)))


def sample_local_params(
    best_params: Dict[str, float],
    ranges: Dict[str, List[float]],
    radius: float,
    rng: np.random.Generator,
) -> Dict[str, float]:
    p = {}
    for k in sorted(best_params.keys()):
        lo, hi = float(ranges[k][0]), float(ranges[k][1])
        span = hi - lo
        delta = rng.uniform(-radius, radius) * span
        p[k] = clip_to_range(float(best_params[k]) + float(delta), lo, hi)
    return p


def bootstrap_triad_pass_rate(
    params: Dict[str, float],
    kernel: Dict[str, float],
    thr: Dict[str, float],
    mass_flags: Dict[str, bool],
    flavor: Dict[str, float],
    s_hl: np.ndarray,
    s_hv: np.ndarray,
    s_lv: np.ndarray,
    n_boot: int,
    rng: np.random.Generator,
) -> float:
    flavor_flags = {
        "ckm_mean_rel_pct_le_max": bool(flavor["ckm_mean_rel_pct"] <= thr["ckm_mean_rel_pct_max"]),
        "pmns_mean_rel_pct_le_max": bool(flavor["pmns_mean_rel_pct"] <= thr["pmns_mean_rel_pct_max"]),
    }

    n_hl, n_hv, n_lv = len(s_hl), len(s_hv), len(s_lv)
    pass_count = 0
    for _ in range(n_boot):
        b_hl = s_hl[rng.integers(0, n_hl, size=n_hl, endpoint=False)]
        b_hv = s_hv[rng.integers(0, n_hv, size=n_hv, endpoint=False)]
        b_lv = s_lv[rng.integers(0, n_lv, size=n_lv, endpoint=False)]
        g = gw_metrics_from_score_arrays(b_hl, b_hv, b_lv)
        flags = {
            **mass_flags,
            **flavor_flags,
            "gw_sep_ge_min": bool(g["sep_median_h1l1_minus_ctrl"] >= thr["gw_sep_min"]),
            "gw_adv_ge_min": bool(g["adv_shared_minus_ctrl_q90"] >= thr["gw_adv_min"]),
            "gw_auc_ge_min": bool(g["auc_h1l1_vs_ctrl"] >= thr["gw_auc_min"]),
            "gw_control_gap_le_max": bool(g["control_median_gap"] <= thr["gw_control_gap_max"]),
        }
        if all(flags.values()):
            pass_count += 1
    return float(pass_count / n_boot)


def deterministic_local_pass_rate(
    center_params: Dict[str, float],
    ranges: Dict[str, List[float]],
    radii: List[float],
    n_per_radius: int,
    kernel: Dict[str, float],
    q_up: np.ndarray,
    q_down: np.ndarray,
    q_nu: np.ndarray,
    q_lep: np.ndarray,
    gw_features: np.ndarray,
    pairs: np.ndarray,
    thr: Dict[str, float],
    mass_flags: Dict[str, bool],
    seed: int,
) -> List[Dict[str, object]]:
    rng = np.random.default_rng(seed)
    rows = []
    for radius in radii:
        pass_count = 0
        for _ in range(n_per_radius):
            p = sample_local_params(center_params, ranges, radius, rng)
            f = flavor_metrics(q_up, q_down, q_nu, q_lep, p, kernel)
            g = gw_metrics(p, kernel, gw_features, pairs)
            flags = triad_flags(mass_flags, f, g, thr)
            if all(flags.values()):
                pass_count += 1
        rows.append(
            {
                "radius_fraction_of_range": float(radius),
                "n_samples": int(n_per_radius),
                "pass_count": int(pass_count),
                "pass_rate": float(pass_count / n_per_radius),
            }
        )
    return rows


def verdict(best_boot_rate: float, local_rows: List[Dict[str, object]]) -> Tuple[str, str]:
    by_r = {float(r["radius_fraction_of_range"]): float(r["pass_rate"]) for r in local_rows}
    r01 = by_r.get(0.01, 0.0)
    r02 = by_r.get(0.02, 0.0)
    r05 = by_r.get(0.05, 0.0)

    if best_boot_rate >= 0.95 and r01 >= 0.95 and r02 >= 0.90 and r05 >= 0.70:
        return (
            "LOCKABLE_BOOTSTRAP_STABLE_TRIAD",
            "FREEZE_RECENTERED_OPERATOR_AND_RUN_TRUE_EXTERNAL_CONFIRMATORY",
        )
    if best_boot_rate >= 0.80 and r01 >= 0.90 and r02 >= 0.80:
        return (
            "PARTIAL_BOOTSTRAP_ROBUSTNESS_IMPROVEMENT",
            "RUN_TARGETED_GW_CONTROL_GAP_REGULARIZATION_THEN_REPEAT_QW1969",
        )
    return (
        "INSUFFICIENT_BOOTSTRAP_ROBUSTNESS",
        "ADD_STRUCTURAL_GW_CONTROL_TERM_WITH_SHARED_ORIGIN_AND_REPEAT",
    )


def main() -> None:
    r1967 = json.loads(IN_QW1967.read_text(encoding="utf-8"))
    r1968 = json.loads(IN_QW1968.read_text(encoding="utf-8"))
    kernel = r1967["kernel"]
    thr = r1967["thresholds"]
    mass_flags = r1967["source_mass_branch"]["mass_flags_fixed"]
    best_params = r1967["summary"]["best"]["params"]
    ranges = r1967["search"]["ranges"]
    q_nu = np.array(r1967["summary"]["best"]["q_nu"], dtype=float)

    r1962 = json.loads((ROOT / "report_qw1962_noncircular_branch_unified_triad_gate.json").read_text(encoding="utf-8"))
    q_assign = r1962["q_assignment_used"]
    q_up = np.array([q_assign["Top"], q_assign["Charm"], q_assign["Muon"]], dtype=float)
    q_down = np.array([q_assign["Bottom"], q_assign["Tau"], q_assign["Muon"]], dtype=float)
    q_lep = np.array([q_assign["Electron"], q_assign["Muon"], q_assign["Tau"]], dtype=float)

    df_gw = pd.read_csv(ROOT / "gw1831_window_features.csv")
    gw_features = df_gw[["max_abs_corr", "mean_abs_corr", "corr_at_0ms", "corr_at_10ms"]].to_numpy(dtype=float)
    pair_map = {"H1-L1": 0, "H1-V1": 1, "L1-V1": 2}
    pairs = df_gw["pair"].map(pair_map).to_numpy(dtype=int)

    # Baseline metrics from QW-1967 point.
    base_f = flavor_metrics(q_up, q_down, q_nu, q_lep, best_params, kernel)
    base_g = gw_metrics(best_params, kernel, gw_features, pairs)
    base_margins = normalized_margins(base_f, base_g, thr)
    base_score = robust_margin_score(base_margins)
    base_flags = triad_flags(mass_flags, base_f, base_g, thr)

    w_base = gw_weights(best_params, kernel)
    s_base = gw_features @ w_base
    s_base_hl = s_base[pairs == 0]
    s_base_hv = s_base[pairs == 1]
    s_base_lv = s_base[pairs == 2]

    baseline_boot = {
        "triad_pass_rate_screen_1200": bootstrap_triad_pass_rate(
            params=best_params,
            kernel=kernel,
            thr=thr,
            mass_flags=mass_flags,
            flavor=base_f,
            s_hl=s_base_hl,
            s_hv=s_base_hv,
            s_lv=s_base_lv,
            n_boot=1200,
            rng=np.random.default_rng(3969),
        ),
        "triad_pass_rate_from_qw1968_5000": float(r1968["bootstrap_gw_and_triad"]["triad_pass_rate"]),
    }

    # Stage A: local deterministic pass search.
    rng = np.random.default_rng(1969)
    n_search = 60000
    search_radii = [0.01, 0.02, 0.05]
    search_probs = [0.45, 0.35, 0.20]
    pass_candidates: List[Dict[str, object]] = []

    for _ in range(n_search):
        radius = float(rng.choice(search_radii, p=search_probs))
        p = sample_local_params(best_params, ranges, radius, rng)
        f = flavor_metrics(q_up, q_down, q_nu, q_lep, p, kernel)
        g = gw_metrics(p, kernel, gw_features, pairs)
        flags = triad_flags(mass_flags, f, g, thr)
        if all(flags.values()):
            margins = normalized_margins(f, g, thr)
            pass_candidates.append(
                {
                    "radius_drawn": radius,
                    "params": p,
                    "flavor": f,
                    "gw": g,
                    "margins": margins,
                    "robust_margin_score": robust_margin_score(margins),
                }
            )

    pass_candidates.append(
        {
            "radius_drawn": 0.0,
            "params": best_params,
            "flavor": base_f,
            "gw": base_g,
            "margins": base_margins,
            "robust_margin_score": base_score,
        }
    )
    pass_candidates.sort(key=lambda x: x["robust_margin_score"], reverse=True)

    top_k = min(100, len(pass_candidates))
    top_candidates = pass_candidates[:top_k]

    # Stage B: bootstrap screening for the best deterministic candidates.
    screened = []
    for i, c in enumerate(top_candidates):
        p = c["params"]
        w = gw_weights(p, kernel)
        s = gw_features @ w
        s_hl = s[pairs == 0]
        s_hv = s[pairs == 1]
        s_lv = s[pairs == 2]
        boot_rate = bootstrap_triad_pass_rate(
            params=p,
            kernel=kernel,
            thr=thr,
            mass_flags=mass_flags,
            flavor=c["flavor"],
            s_hl=s_hl,
            s_hv=s_hv,
            s_lv=s_lv,
            n_boot=1200,
            rng=np.random.default_rng(5000 + i),
        )
        screened.append({**c, "triad_boot_pass_rate_screen_1200": float(boot_rate)})

    screened.sort(
        key=lambda x: (
            x["triad_boot_pass_rate_screen_1200"],
            x["robust_margin_score"],
        ),
        reverse=True,
    )

    finalists = screened[:5]
    finals = []
    for i, c in enumerate(finalists):
        p = c["params"]
        w = gw_weights(p, kernel)
        s = gw_features @ w
        s_hl = s[pairs == 0]
        s_hv = s[pairs == 1]
        s_lv = s[pairs == 2]
        boot_rate_5000 = bootstrap_triad_pass_rate(
            params=p,
            kernel=kernel,
            thr=thr,
            mass_flags=mass_flags,
            flavor=c["flavor"],
            s_hl=s_hl,
            s_hv=s_hv,
            s_lv=s_lv,
            n_boot=5000,
            rng=np.random.default_rng(7000 + i),
        )
        finals.append({**c, "triad_boot_pass_rate_final_5000": float(boot_rate_5000)})

    finals.sort(
        key=lambda x: (
            x["triad_boot_pass_rate_final_5000"],
            x["robust_margin_score"],
        ),
        reverse=True,
    )
    best_final = finals[0]

    local_rows = deterministic_local_pass_rate(
        center_params=best_final["params"],
        ranges=ranges,
        radii=[0.01, 0.02, 0.05],
        n_per_radius=8000,
        kernel=kernel,
        q_up=q_up,
        q_down=q_down,
        q_nu=q_nu,
        q_lep=q_lep,
        gw_features=gw_features,
        pairs=pairs,
        thr=thr,
        mass_flags=mass_flags,
        seed=8969,
    )

    v, nxt = verdict(best_final["triad_boot_pass_rate_final_5000"], local_rows)

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "qw1967": IN_QW1967.name,
            "qw1968": IN_QW1968.name,
        },
        "frozen_kernel": kernel,
        "fixed_q_nu": [int(x) for x in q_nu.tolist()],
        "thresholds": thr,
        "baseline": {
            "params": best_params,
            "flavor": base_f,
            "gw": base_g,
            "flags": base_flags,
            "margins": base_margins,
            "robust_margin_score": base_score,
            "bootstrap": baseline_boot,
        },
        "search": {
            "seed": 1969,
            "n_search": int(n_search),
            "search_radii": search_radii,
            "search_probs": search_probs,
            "deterministic_pass_candidates_count": int(len(pass_candidates)),
            "top_k_for_bootstrap_screen": int(top_k),
        },
        "screened_top10": screened[:10],
        "finalists": finals,
        "best_recentered_candidate": best_final,
        "best_recentered_local_pass_rates": local_rows,
        "verdict": v,
        "required_next_step": nxt,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1969: BOOTSTRAP ROBUST RECENTER SEARCH",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{v}**",
        "",
        "## Baseline vs Best Recentered",
        (
            f"- baseline triad bootstrap pass (5000 from QW-1968): "
            f"{100.0 * baseline_boot['triad_pass_rate_from_qw1968_5000']:.2f}%"
        ),
        (
            f"- best recentered triad bootstrap pass (5000): "
            f"{100.0 * best_final['triad_boot_pass_rate_final_5000']:.2f}%"
        ),
        (
            f"- deterministic candidate pool (triad pass): "
            f"{len(pass_candidates)} / {n_search + 1}"
        ),
        "",
        "## Best Recentered Deterministic Metrics",
        (
            f"- CKM/PMNS mean rel%: "
            f"{best_final['flavor']['ckm_mean_rel_pct']:.3f}/"
            f"{best_final['flavor']['pmns_mean_rel_pct']:.3f}"
        ),
        (
            f"- GW auc/adv/sep/gap: "
            f"{best_final['gw']['auc_h1l1_vs_ctrl']:.4f}/"
            f"{best_final['gw']['adv_shared_minus_ctrl_q90']:.4f}/"
            f"{best_final['gw']['sep_median_h1l1_minus_ctrl']:.6f}/"
            f"{best_final['gw']['control_median_gap']:.6f}"
        ),
        "",
        "## Local Deterministic Pass Around Best Recentered",
        *[
            (
                f"- radius={r['radius_fraction_of_range']:.3f}: "
                f"{100.0 * r['pass_rate']:.2f}% ({r['pass_count']}/{r['n_samples']})"
            )
            for r in local_rows
        ],
        "",
        "## Required Next Step",
        f"- {nxt}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1969] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1969] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1969] verdict={v}")


if __name__ == "__main__":
    main()
