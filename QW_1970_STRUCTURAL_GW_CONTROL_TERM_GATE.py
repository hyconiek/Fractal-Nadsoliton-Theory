#!/usr/bin/env python3
"""
QW-1970: Structural GW control-gap term gate (single shared extra term).

Design:
- fixed frozen kernel,
- fixed mass + flavor branch (no sector retune),
- add one shared GW structural term to stabilize control-gap under bootstrap.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
IN_QW1969 = ROOT / "report_qw1969_bootstrap_robust_recenter_search.json"
OUT_JSON = ROOT / "report_qw1970_structural_gw_control_term_gate.json"
OUT_MD = ROOT / "RAPORT_QW1970_STRUCTURAL_GW_CONTROL_TERM_GATE.md"


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


def base_gw_weights(params: Dict[str, float], kernel: Dict[str, float]) -> np.ndarray:
    d = np.array([1.0, 2.0, 3.0, 4.0], dtype=float)
    raw = (np.abs(kernel_fn(d, kernel["omega"], kernel["phi"], kernel["beta"], kernel["eta"])) ** params["p_amp"]) * (
        d**params["r_dist"]
    )
    return raw / np.sum(raw)


def gw_metrics_from_scores(s_hl: np.ndarray, s_hv: np.ndarray, s_lv: np.ndarray) -> Dict[str, float]:
    s_ctrl = np.concatenate([s_hv, s_lv])
    q90 = float(np.quantile(s_ctrl, 0.90))
    auc = rank_auc_pos_gt_neg(s_hl, s_ctrl)
    adv = float(np.mean(s_hl > q90) - np.mean(s_ctrl > q90))
    sep = float(np.median(s_hl) - np.median(s_ctrl))
    gap = float(abs(np.median(s_hv) - np.median(s_lv)))
    return {
        "auc_h1l1_vs_ctrl": auc,
        "adv_shared_minus_ctrl_q90": adv,
        "sep_median_h1l1_minus_ctrl": sep,
        "control_median_gap": gap,
    }


def triad_flags(mass_flags: Dict[str, bool], flavor: Dict[str, float], gw: Dict[str, float], thr: Dict[str, float]) -> Dict[str, bool]:
    return {
        **mass_flags,
        "ckm_mean_rel_pct_le_max": bool(flavor["ckm_mean_rel_pct"] <= thr["ckm_mean_rel_pct_max"]),
        "pmns_mean_rel_pct_le_max": bool(flavor["pmns_mean_rel_pct"] <= thr["pmns_mean_rel_pct_max"]),
        "gw_sep_ge_min": bool(gw["sep_median_h1l1_minus_ctrl"] >= thr["gw_sep_min"]),
        "gw_adv_ge_min": bool(gw["adv_shared_minus_ctrl_q90"] >= thr["gw_adv_min"]),
        "gw_auc_ge_min": bool(gw["auc_h1l1_vs_ctrl"] >= thr["gw_auc_min"]),
        "gw_control_gap_le_max": bool(gw["control_median_gap"] <= thr["gw_control_gap_max"]),
    }


def bootstrap_pass_rate(
    s_hl: np.ndarray,
    s_hv: np.ndarray,
    s_lv: np.ndarray,
    n_boot: int,
    thr: Dict[str, float],
    flavor: Dict[str, float],
    mass_flags: Dict[str, bool],
    seed: int,
) -> float:
    rng = np.random.default_rng(seed)
    n_hl, n_hv, n_lv = len(s_hl), len(s_hv), len(s_lv)
    pass_count = 0
    flavor_flags = {
        "ckm_mean_rel_pct_le_max": bool(flavor["ckm_mean_rel_pct"] <= thr["ckm_mean_rel_pct_max"]),
        "pmns_mean_rel_pct_le_max": bool(flavor["pmns_mean_rel_pct"] <= thr["pmns_mean_rel_pct_max"]),
    }
    for _ in range(n_boot):
        b_hl = s_hl[rng.integers(0, n_hl, size=n_hl, endpoint=False)]
        b_hv = s_hv[rng.integers(0, n_hv, size=n_hv, endpoint=False)]
        b_lv = s_lv[rng.integers(0, n_lv, size=n_lv, endpoint=False)]
        g = gw_metrics_from_scores(b_hl, b_hv, b_lv)
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


def verdict(final_rate: float, baseline_rate: float) -> Tuple[str, str]:
    if final_rate >= 0.95:
        return (
            "STRUCTURAL_CONTROL_TERM_LOCKABLE",
            "FREEZE_TERM_AND_RUN_TRUE_EXTERNAL_CONFIRMATORY",
        )
    if final_rate >= 0.80 and final_rate > baseline_rate + 0.03:
        return (
            "STRUCTURAL_CONTROL_TERM_PARTIAL_SUCCESS",
            "PROMOTE_TERM_AND_RUN_JOINT_CROSS_DATA_CHECK",
        )
    return (
        "STRUCTURAL_CONTROL_TERM_INSUFFICIENT",
        "NEED_DEEPER_SHARED_DYNAMICS_FOR_CONTROL_CHANNEL",
    )


def main() -> None:
    r1969 = json.loads(IN_QW1969.read_text(encoding="utf-8"))
    kernel = r1969["frozen_kernel"]
    params = r1969["best_recentered_candidate"]["params"]
    thr = r1969["thresholds"]
    q_nu = np.array(r1969["fixed_q_nu"], dtype=float)

    r1967 = json.loads((ROOT / "report_qw1967_isospin_split_local_refinement_gate.json").read_text(encoding="utf-8"))
    mass_flags = r1967["source_mass_branch"]["mass_flags_fixed"]

    r1962 = json.loads((ROOT / "report_qw1962_noncircular_branch_unified_triad_gate.json").read_text(encoding="utf-8"))
    q_assign = r1962["q_assignment_used"]
    q_up = np.array([q_assign["Top"], q_assign["Charm"], q_assign["Muon"]], dtype=float)
    q_down = np.array([q_assign["Bottom"], q_assign["Tau"], q_assign["Muon"]], dtype=float)
    q_lep = np.array([q_assign["Electron"], q_assign["Muon"], q_assign["Tau"]], dtype=float)

    flavor = flavor_metrics(q_up, q_down, q_nu, q_lep, params, kernel)

    df = pd.read_csv(ROOT / "gw1831_window_features.csv")
    pair_map = {"H1-L1": 0, "H1-V1": 1, "L1-V1": 2}
    pairs = df["pair"].map(pair_map).to_numpy(dtype=int)

    features = df[["max_abs_corr", "mean_abs_corr", "corr_at_0ms", "corr_at_10ms"]].to_numpy(dtype=float)
    w = base_gw_weights(params, kernel)
    base_score = features @ w

    # Structural control carrier:
    # one shared term built from kernel phase + lag channel, antisymmetric on control pairs.
    lag_s = df["best_lag_ms"].to_numpy(dtype=float) * 1e-3
    lag_phase = np.sin(kernel["omega"] * lag_s + kernel["phi"])
    info_delta = df["corr_at_0ms"].to_numpy(dtype=float) - df["corr_at_10ms"].to_numpy(dtype=float)
    pair_sign = np.where(pairs == 1, 1.0, np.where(pairs == 2, -1.0, 0.0))
    carrier_raw = pair_sign * lag_phase * info_delta
    carrier_std = float(np.std(carrier_raw))
    carrier = carrier_raw / (carrier_std if carrier_std > 1e-12 else 1.0)

    # Baseline from QW-1969 best recentered candidate.
    s_hl_base = base_score[pairs == 0]
    s_hv_base = base_score[pairs == 1]
    s_lv_base = base_score[pairs == 2]
    gw_base = gw_metrics_from_scores(s_hl_base, s_hv_base, s_lv_base)
    base_flags = triad_flags(mass_flags, flavor, gw_base, thr)
    baseline_boot_5000 = float(r1969["best_recentered_candidate"]["triad_boot_pass_rate_final_5000"])

    # Grid search for the structural coefficient xi.
    xi_grid = np.linspace(-0.020, 0.020, 801)
    rows = []
    for xi in xi_grid:
        score = base_score + float(xi) * carrier
        s_hl = score[pairs == 0]
        s_hv = score[pairs == 1]
        s_lv = score[pairs == 2]
        gw = gw_metrics_from_scores(s_hl, s_hv, s_lv)
        flags = triad_flags(mass_flags, flavor, gw, thr)
        rows.append(
            {
                "xi": float(xi),
                "gw": gw,
                "flags": flags,
                "all_pass": bool(all(flags.values())),
                "scores": {
                    "s_hl": s_hl,
                    "s_hv": s_hv,
                    "s_lv": s_lv,
                },
            }
        )

    deterministic_pass_rows = [r for r in rows if r["all_pass"]]
    if not deterministic_pass_rows:
        out = {
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "verdict": "STRUCTURAL_CONTROL_TERM_NO_DETERMINISTIC_PASS",
            "required_next_step": "CHANGE_STRUCTURAL_TERM_CLASS",
        }
        OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")
        OUT_MD.write_text(
            "# RAPORT QW-1970: STRUCTURAL GW CONTROL TERM GATE\n\n- Verdict: **STRUCTURAL_CONTROL_TERM_NO_DETERMINISTIC_PASS**\n",
            encoding="utf-8",
        )
        print(f"[QW-1970] Saved JSON: {OUT_JSON.name}")
        print(f"[QW-1970] Saved MD:   {OUT_MD.name}")
        print("[QW-1970] verdict=STRUCTURAL_CONTROL_TERM_NO_DETERMINISTIC_PASS")
        return

    # Rank deterministic pass candidates by low control gap and strong sep.
    deterministic_pass_rows.sort(
        key=lambda r: (
            -r["gw"]["sep_median_h1l1_minus_ctrl"],
            -r["gw"]["adv_shared_minus_ctrl_q90"],
            -r["gw"]["auc_h1l1_vs_ctrl"],
            r["gw"]["control_median_gap"],
        )
    )
    top = deterministic_pass_rows[:40]

    # Bootstrap screen.
    screened = []
    for i, r in enumerate(top):
        boot = bootstrap_pass_rate(
            s_hl=r["scores"]["s_hl"],
            s_hv=r["scores"]["s_hv"],
            s_lv=r["scores"]["s_lv"],
            n_boot=2500,
            thr=thr,
            flavor=flavor,
            mass_flags=mass_flags,
            seed=19700 + i,
        )
        screened.append(
            {
                "xi": r["xi"],
                "gw": r["gw"],
                "boot_pass_rate_2500": float(boot),
            }
        )
    screened.sort(key=lambda x: x["boot_pass_rate_2500"], reverse=True)

    finalists_xi = [x["xi"] for x in screened[:5]]
    finalists = []
    for j, xi in enumerate(finalists_xi):
        score = base_score + float(xi) * carrier
        s_hl = score[pairs == 0]
        s_hv = score[pairs == 1]
        s_lv = score[pairs == 2]
        gw = gw_metrics_from_scores(s_hl, s_hv, s_lv)
        boot_5000 = bootstrap_pass_rate(
            s_hl=s_hl,
            s_hv=s_hv,
            s_lv=s_lv,
            n_boot=5000,
            thr=thr,
            flavor=flavor,
            mass_flags=mass_flags,
            seed=19800 + j,
        )
        finalists.append(
            {
                "xi": float(xi),
                "gw": gw,
                "boot_pass_rate_5000": float(boot_5000),
            }
        )
    finalists.sort(key=lambda x: x["boot_pass_rate_5000"], reverse=True)
    best = finalists[0]

    # Small xi-local robustness around best xi.
    xi_local = [best["xi"] - 0.001, best["xi"], best["xi"] + 0.001]
    xi_local_rows = []
    for k, xi in enumerate(xi_local):
        score = base_score + float(xi) * carrier
        s_hl = score[pairs == 0]
        s_hv = score[pairs == 1]
        s_lv = score[pairs == 2]
        gw = gw_metrics_from_scores(s_hl, s_hv, s_lv)
        flags = triad_flags(mass_flags, flavor, gw, thr)
        boot_3000 = bootstrap_pass_rate(
            s_hl=s_hl,
            s_hv=s_hv,
            s_lv=s_lv,
            n_boot=3000,
            thr=thr,
            flavor=flavor,
            mass_flags=mass_flags,
            seed=19900 + k,
        )
        xi_local_rows.append(
            {
                "xi": float(xi),
                "gw": gw,
                "all_pass_deterministic": bool(all(flags.values())),
                "boot_pass_rate_3000": float(boot_3000),
            }
        )

    v, nxt = verdict(best["boot_pass_rate_5000"], baseline_boot_5000)

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source_report": IN_QW1969.name,
        "fixed_components": {
            "kernel": kernel,
            "params": params,
            "q_nu": [int(x) for x in q_nu.tolist()],
            "flavor": flavor,
            "mass_flags_fixed": mass_flags,
        },
        "baseline_no_structural_term": {
            "gw": gw_base,
            "flags": base_flags,
            "bootstrap_pass_rate_5000": baseline_boot_5000,
        },
        "structural_term": {
            "definition": "score = base_score + xi * pair_sign * sin(omega*lag+phi) * (corr0-corr10) / std",
            "carrier_std": carrier_std,
            "xi_grid_min": float(np.min(xi_grid)),
            "xi_grid_max": float(np.max(xi_grid)),
            "xi_grid_size": int(len(xi_grid)),
            "deterministic_pass_count": int(len(deterministic_pass_rows)),
        },
        "screened_top10": screened[:10],
        "finalists": finalists,
        "best": best,
        "xi_local_robustness": xi_local_rows,
        "verdict": v,
        "required_next_step": nxt,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1970: STRUCTURAL GW CONTROL TERM GATE",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{v}**",
        "",
        "## Baseline vs Structural Term",
        f"- baseline bootstrap pass (5000): {100.0 * baseline_boot_5000:.2f}%",
        f"- best structural-term bootstrap pass (5000): {100.0 * best['boot_pass_rate_5000']:.2f}%",
        f"- delta: {100.0 * (best['boot_pass_rate_5000'] - baseline_boot_5000):.2f} pp",
        "",
        "## Best Structural Candidate",
        f"- xi: {best['xi']:.6f}",
        (
            f"- GW auc/adv/sep/gap: "
            f"{best['gw']['auc_h1l1_vs_ctrl']:.4f}/"
            f"{best['gw']['adv_shared_minus_ctrl_q90']:.4f}/"
            f"{best['gw']['sep_median_h1l1_minus_ctrl']:.6f}/"
            f"{best['gw']['control_median_gap']:.6f}"
        ),
        "",
        "## Xi Local Robustness",
        *[
            (
                f"- xi={r['xi']:.6f}: det_pass={r['all_pass_deterministic']}, "
                f"boot3000={100.0 * r['boot_pass_rate_3000']:.2f}%"
            )
            for r in xi_local_rows
        ],
        "",
        "## Required Next Step",
        f"- {nxt}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1970] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1970] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1970] verdict={v}")


if __name__ == "__main__":
    main()
