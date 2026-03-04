#!/usr/bin/env python3
"""
QW-2059: Grepped dedup + historical flavor transfer audit.

Goal:
- verify by grep which flavor/generator paths were already executed,
- avoid duplicate reruns of same methodological family,
- test transfer of strongest historical operators onto current frozen kernel (QW-2049).
"""

from __future__ import annotations

import json
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2059_grepped_dedup_historical_flavor_transfer_audit.json"
OUT_MD = ROOT / "RAPORT_QW2059_GREPPED_DEDUP_HISTORICAL_FLAVOR_TRANSFER_AUDIT.md"


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


def run_rg_files() -> List[str]:
    proc = subprocess.run(["rg", "--files"], cwd=ROOT, capture_output=True, text=True, check=False)
    if proc.returncode not in (0, 1):
        return []
    return [x.strip() for x in proc.stdout.splitlines() if x.strip()]


def grep_dedup_index() -> Dict[str, object]:
    files = run_rg_files()
    pattern = re.compile(r"(QW_|RAPORT_QW|report_qw).*(flavor|ckm|pmns|su3|nonabel|generator)", re.IGNORECASE)
    indexed = sorted([f for f in files if pattern.search(f)])

    key_methods = {
        "QW_1966_isospin_split_scan": "QW_1966_ISOSPIN_SPLIT_SHARED_FLAVOR_DYNAMICS_SCAN.py",
        "QW_2029_shared_flavor_basis_scan": "QW_2029_CKM_BLOCKER_SHARED_FLAVOR_BASIS_SCAN.py",
        "QW_2012_no_ansatz_no_fit": "QW_2012_SINGLE_KERNEL_FLAVOR_NO_ANSATZ_STRICT.py",
        "QW_2056_operator_family_frontier": "QW_2056_FIRST_PRINCIPLES_FLAVOR_OPERATOR_FAMILY_FRONTIER.py",
        "QW_2057_su3_rotation_frontier": "QW_2057_SU3_ROTATION_FLAVOR_FRONTIER_NO_FIT.py",
        "QW_2058_nonabelian_no_fit": "QW_2058_NONABELIAN_FLAVOR_GENERATOR_NO_FIT.py",
    }
    exists = {k: bool((ROOT / v).exists()) for k, v in key_methods.items()}

    return {
        "n_indexed_files": int(len(indexed)),
        "indexed_files_top50": indexed[:50],
        "key_methods": key_methods,
        "key_methods_exists": exists,
    }


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


def gw_metrics(kernel: Dict[str, float], p_amp: float, r_dist: float) -> Dict[str, float]:
    d = np.array([1.0, 2.0, 3.0, 4.0], dtype=float)
    raw = (np.abs(kernel_fn(d, kernel["omega"], kernel["phi"], kernel["beta"], kernel["eta"])) ** p_amp) * (d**r_dist)
    w = raw / np.sum(raw)

    df = pd.read_csv(ROOT / "gw1831_window_features.csv")
    score = (
        w[0] * df["max_abs_corr"].to_numpy(dtype=float)
        + w[1] * df["mean_abs_corr"].to_numpy(dtype=float)
        + w[2] * df["corr_at_0ms"].to_numpy(dtype=float)
        + w[3] * df["corr_at_10ms"].to_numpy(dtype=float)
    )
    pair = df["pair"].astype(str).to_numpy()
    s_hl = score[pair == "H1-L1"]
    s_hv = score[pair == "H1-V1"]
    s_lv = score[pair == "L1-V1"]
    s_ctrl = np.concatenate([s_hv, s_lv])
    q90 = float(np.quantile(s_ctrl, 0.90))
    return {
        "auc_h1l1_vs_ctrl": float(rank_auc_pos_gt_neg(s_hl, s_ctrl)),
        "adv_shared_minus_ctrl_q90": float(np.mean(s_hl > q90) - np.mean(s_ctrl > q90)),
        "sep_median_h1l1_minus_ctrl": float(np.median(s_hl) - np.median(s_ctrl)),
        "control_median_gap": float(abs(np.median(s_hv) - np.median(s_lv))),
    }


def flavor_ham_1966(
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


def flavor_ham_2029(
    q: np.ndarray,
    iso_tag: float,
    sector_tag: float,
    p_amp: float,
    r_dist: float,
    params: Dict[str, float],
    kernel: Dict[str, float],
) -> np.ndarray:
    n = len(q)
    d = 1.0 + cyclic_distance_matrix(q, q, modulus=24)
    k = kernel_fn(d, kernel["omega"], kernel["phi"], kernel["beta"], kernel["eta"])

    q_diff = q[:, None] - q[None, :]
    q_abs = np.abs(q_diff) / 24.0

    base_amp = np.sign(k) * (np.abs(k) ** p_amp) * (d**r_dist)
    amp_mod = 1.0 + params["amp_qbias"] * (q_abs - float(np.mean(q_abs)))
    base = base_amp * amp_mod

    idx = np.arange(n, dtype=float)
    i_minus_j = idx[:, None] - idx[None, :]
    gap = np.abs(i_minus_j)
    near = np.exp(-params["rho_gap"] * gap)

    phase = (
        params["phase_q"] * q_diff
        + params["phase_q3"] * ((q_diff / 24.0) ** 3)
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
        + iso_tag * params.get("diag_iso", 0.0) * idx_centered
        + sector_tag * params.get("diag_sector", 0.0) * idx_centered
    )

    h = re + 1j * params["chi_im"] * im_asym + np.diag(diag)
    return 0.5 * (h + h.conj().T)


def flavor_err_from_evec(hu: np.ndarray, hd: np.ndarray, hn: np.ndarray, hl: np.ndarray) -> Dict[str, float]:
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


def main() -> None:
    grep_index = grep_dedup_index()

    d2049 = json.loads((ROOT / "report_qw2049_spectral_micro_stagec_intersection_gate.json").read_text(encoding="utf-8"))
    kernel = {k: float(v) for k, v in d2049["stagec_pool"]["selected_kernel"].items()}

    r1961 = json.loads((ROOT / "report_qw1961_noncircular_gamma_q_derivation_matrix.json").read_text(encoding="utf-8"))
    q_name = str(r1961["summary"]["best_noncircular"]["q_assignment"])
    q_map = r1961["inputs"]["q_assignments"][q_name]
    q_up = np.array([q_map["Top"], q_map["Charm"], q_map["Muon"]], dtype=float)
    q_down = np.array([q_map["Bottom"], q_map["Tau"], q_map["Muon"]], dtype=float)
    q_lep = np.array([q_map["Electron"], q_map["Muon"], q_map["Tau"]], dtype=float)

    # Case A: QW-1966 best row transferred to current kernel.
    r1966 = json.loads((ROOT / "report_qw1966_isospin_split_shared_flavor_dynamics_scan.json").read_text(encoding="utf-8"))
    c1966 = r1966["summary"]["best"]
    p1966 = c1966["params"]
    q_nu_1966 = np.array(c1966["q_nu"], dtype=float)
    hu = flavor_ham_1966(q_up, +1.0, +1.0, p1966, kernel)
    hd = flavor_ham_1966(q_down, -1.0, +1.0, p1966, kernel)
    hn = flavor_ham_1966(q_nu_1966, +1.0, -1.0, p1966, kernel)
    hl = flavor_ham_1966(q_lep, -1.0, -1.0, p1966, kernel)
    f1966 = flavor_err_from_evec(hu, hd, hn, hl)
    g1966 = gw_metrics(kernel, p_amp=float(p1966["p_amp"]), r_dist=float(p1966["r_dist"]))

    # Case B: QW-2029 best row transferred to current kernel.
    r2029 = json.loads((ROOT / "report_qw2029_ckm_blocker_shared_flavor_basis_scan.json").read_text(encoding="utf-8"))
    c2029 = r2029["best_row"]
    p2029 = c2029["params"]
    q_nu_2029 = np.array(c2029["q_nu"], dtype=float)
    p_amp_2029 = float(r2029["fixed_context"]["p_amp"])
    r_dist_2029 = float(r2029["fixed_context"]["r_dist"])
    hu2 = flavor_ham_2029(q_up, +1.0, +1.0, p_amp_2029, r_dist_2029, p2029, kernel)
    hd2 = flavor_ham_2029(q_down, -1.0, +1.0, p_amp_2029, r_dist_2029, p2029, kernel)
    hn2 = flavor_ham_2029(q_nu_2029, +1.0, -1.0, p_amp_2029, r_dist_2029, p2029, kernel)
    hl2 = flavor_ham_2029(q_lep, -1.0, -1.0, p_amp_2029, r_dist_2029, p2029, kernel)
    f2029 = flavor_err_from_evec(hu2, hd2, hn2, hl2)
    g2029 = gw_metrics(kernel, p_amp=p_amp_2029, r_dist=r_dist_2029)

    # Reference current strict nonabelian no-fit result.
    r2058 = json.loads((ROOT / "report_qw2058_nonabelian_flavor_generator_no_fit.json").read_text(encoding="utf-8"))
    f2058 = {
        "ckm_mean_rel_pct": float(r2058["metrics"]["flavor_nonabelian"]["ckm_mean_rel_pct"]),
        "pmns_mean_rel_pct": float(r2058["metrics"]["flavor_nonabelian"]["pmns_mean_rel_pct"]),
    }
    g2058 = r2058["metrics"]["gw"]

    thresholds = {
        "ckm_mean_rel_pct_max": 15.0,
        "pmns_mean_rel_pct_max": 15.0,
        "gw_sep_min": 0.0020,
        "gw_adv_min": 0.30,
        "gw_auc_min": 0.75,
        "gw_control_gap_max": 0.0025,
    }

    def flags(f: Dict[str, float], g: Dict[str, float], first_principles: bool) -> Dict[str, bool]:
        return {
            "ckm_mean_rel_pct_le_max": bool(f["ckm_mean_rel_pct"] <= thresholds["ckm_mean_rel_pct_max"]),
            "pmns_mean_rel_pct_le_max": bool(f["pmns_mean_rel_pct"] <= thresholds["pmns_mean_rel_pct_max"]),
            "gw_sep_ge_min": bool(g["sep_median_h1l1_minus_ctrl"] >= thresholds["gw_sep_min"]),
            "gw_adv_ge_min": bool(g["adv_shared_minus_ctrl_q90"] >= thresholds["gw_adv_min"]),
            "gw_auc_ge_min": bool(g["auc_h1l1_vs_ctrl"] >= thresholds["gw_auc_min"]),
            "gw_control_gap_le_max": bool(g["control_median_gap"] <= thresholds["gw_control_gap_max"]),
            "first_principles_no_fit_protocol": bool(first_principles),
        }

    eval_rows = {
        "transfer_qw1966_best_to_qw2049": {
            "origin": "QW-1966 summary.best",
            "method_family": "isospin_split_shared_flavor_dynamics_scan",
            "flavor": f1966,
            "gw": g1966,
            "flags": flags(f1966, g1966, first_principles=False),
            "notes": "historical scan-derived parameters, not no-fit first-principles",
        },
        "transfer_qw2029_best_to_qw2049": {
            "origin": "QW-2029 best_row + fixed_context",
            "method_family": "shared_flavor_basis_scan",
            "flavor": f2029,
            "gw": g2029,
            "flags": flags(f2029, g2029, first_principles=False),
            "notes": "historical scan-derived parameters, not no-fit first-principles",
        },
        "reference_qw2058_nonabelian_no_fit": {
            "origin": "QW-2058 report",
            "method_family": "nonabelian_no_fit",
            "flavor": f2058,
            "gw": g2058,
            "flags": flags(f2058, g2058, first_principles=True),
            "notes": "current strict first-principles baseline",
        },
    }

    for _, row in eval_rows.items():
        row["pass_count"] = int(sum(1 for v in row["flags"].values() if v))
        row["total_flags"] = int(len(row["flags"]))
        row["all_pass"] = bool(all(row["flags"].values()))

    any_pass = any(bool(row["all_pass"]) for row in eval_rows.values())
    verdict = (
        "DEDUP_AUDIT_IDENTIFIES_EXISTING_METHODS_AND_NO_STRICT_PASS_UNDER_CURRENT_KERNEL"
        if not any_pass
        else "DEDUP_AUDIT_FOUND_STRICT_PASS"
    )
    required_next = (
        "BUILD_QW2060_FIRST_PRINCIPLES_RECONSTRUCTION_OF_QW2029_FLAVOR_BASIS_FROM_KERNEL_INVARIANTS_ONLY"
        if not any_pass
        else "FREEZE_AND_EXTERNALIZE_STRICT_PASS_BRANCH"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "grep_dedup_index": grep_index,
        "kernel_source": "report_qw2049_spectral_micro_stagec_intersection_gate.json:stagec_pool.selected_kernel",
        "kernel": kernel,
        "q_assignment_source": {
            "file": "report_qw1961_noncircular_gamma_q_derivation_matrix.json",
            "q_assignment_name": q_name,
        },
        "thresholds": thresholds,
        "evaluations": eval_rows,
        "verdict": verdict,
        "required_next_step": required_next,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2059: GREPPED DEDUP + HISTORICAL FLAVOR TRANSFER AUDIT",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        "",
        "## Grep Dedup Index",
        f"- indexed files count: {grep_index['n_indexed_files']}",
    ]
    for k, v in grep_index["key_methods_exists"].items():
        lines.append(f"- {k}: {v}")

    lines.append("")
    lines.append("## Evaluations")
    for name, row in eval_rows.items():
        f = row["flavor"]
        g = row["gw"]
        lines.extend(
            [
                f"### {name}",
                f"- method_family: {row['method_family']}",
                f"- flavor CKM/PMNS mean rel%: {f['ckm_mean_rel_pct']:.3f}/{f['pmns_mean_rel_pct']:.3f}",
                (
                    f"- GW auc/adv/sep/gap: "
                    f"{g['auc_h1l1_vs_ctrl']:.4f}/{g['adv_shared_minus_ctrl_q90']:.4f}/"
                    f"{g['sep_median_h1l1_minus_ctrl']:.6f}/{g['control_median_gap']:.6f}"
                ),
                f"- pass_count: {row['pass_count']}/{row['total_flags']}",
                f"- all_pass: {row['all_pass']}",
                f"- notes: {row['notes']}",
                "",
            ]
        )

    lines.extend(
        [
            "## Required Next Step",
            f"- {required_next}",
            "",
            "## Artifacts",
            f"- JSON: `{OUT_JSON.name}`",
        ]
    )
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2059] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2059] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2059] verdict={verdict}")


if __name__ == "__main__":
    main()
