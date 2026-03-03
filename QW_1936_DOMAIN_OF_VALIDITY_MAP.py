#!/usr/bin/env python3
"""
QW-1936: Domain-of-validity map after mixed head-to-head result.

Builds angular regime diagnostics for reparameterized nadsoliton vs HD proxy:
- per-theta-bin MSE delta (HD - reparam),
- win fractions across deterministic hash splits,
- dataset-level domain statement.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1936_domain_of_validity_map.json"
OUT_MD = ROOT / "RAPORT_QW1936_DOMAIN_OF_VALIDITY_MAP.md"


def split_index(pair_id: str, salt: int, k: int = 2) -> int:
    h = hashlib.sha256(f"{salt}|{pair_id}".encode("utf-8")).hexdigest()
    return int(h[-8:], 16) % k


def hd_curve(theta_deg: np.ndarray) -> np.ndarray:
    th = np.deg2rad(theta_deg)
    x = (1.0 - np.cos(th)) / 2.0
    x = np.clip(x, 1e-12, 1.0)
    return 0.5 + 1.5 * x * np.log(x) - 0.25 * x


def nad_kernel(theta_deg: np.ndarray, omega: float, phi: float, beta: float, eta: float) -> np.ndarray:
    tmin, tmax = float(np.min(theta_deg)), float(np.max(theta_deg))
    d_eff = 1.0 + 11.0 * (theta_deg - tmin) / max(tmax - tmin, 1e-12)
    return np.cos(omega * d_eff + phi) / (1.0 + beta * (d_eff ** eta))


def fit_affine(x: np.ndarray, y: np.ndarray, disc: np.ndarray) -> np.ndarray:
    X = np.column_stack([x[disc], np.ones(int(np.sum(disc)), dtype=float)])
    coef, *_ = np.linalg.lstsq(X, y[disc], rcond=None)
    a, b = float(coef[0]), float(coef[1])
    return a * x + b


def load_reparam_triad() -> Dict[str, float]:
    d = json.loads((ROOT / "report_qw1932_physical_reparameterization_eta_scan.json").read_text(encoding="utf-8"))
    sel = d.get("selected", {})
    fit = sel.get("fit", {})
    return {
        "omega": float(fit["omega"]),
        "phi": float(fit["phi"]),
        "beta": float(fit["beta"]),
        "eta": float(sel["eta"]),
    }


def theta_bins() -> List[Tuple[float, float]]:
    return [(0.0, 30.0), (30.0, 60.0), (60.0, 90.0), (90.0, 120.0), (120.0, 150.0), (150.0, 180.0)]


def evaluate_dataset(df: pd.DataFrame, triad: Dict[str, float], salts: List[int]) -> Dict[str, object]:
    req = {"pair_id", "theta_deg", "hxy"}
    miss = sorted(req - set(df.columns))
    if miss:
        raise RuntimeError(f"Missing columns: {miss}")

    pair_id = df["pair_id"].astype(str).to_numpy()
    theta = df["theta_deg"].to_numpy(dtype=float)
    y = df["hxy"].to_numpy(dtype=float)

    feat_reparam = nad_kernel(
        theta,
        omega=float(triad["omega"]),
        phi=float(triad["phi"]),
        beta=float(triad["beta"]),
        eta=float(triad["eta"]),
    )
    feat_hd = hd_curve(theta)

    bins = theta_bins()
    per_bin = {f"{lo:.0f}-{hi:.0f}": [] for lo, hi in bins}
    global_delta = []

    for salt in salts:
        fold = np.array([split_index(x, salt=salt, k=2) for x in pair_id], dtype=int)
        disc = fold == 0
        hold = fold == 1
        if int(np.sum(disc)) < 300 or int(np.sum(hold)) < 300:
            raise RuntimeError(f"Salt {salt}: too small split.")

        yhat_r = fit_affine(feat_reparam, y, disc)
        yhat_h = fit_affine(feat_hd, y, disc)

        th_h = theta[hold]
        y_h = y[hold]
        e2_r = (y_h - yhat_r[hold]) ** 2
        e2_h = (y_h - yhat_h[hold]) ** 2

        d_global = float(np.mean(e2_h) - np.mean(e2_r))
        global_delta.append(d_global)

        for lo, hi in bins:
            m = (th_h >= lo) & (th_h < hi)
            key = f"{lo:.0f}-{hi:.0f}"
            if int(np.sum(m)) < 20:
                per_bin[key].append(None)
            else:
                per_bin[key].append(float(np.mean(e2_h[m]) - np.mean(e2_r[m])))

    bin_summary = {}
    for key, vals in per_bin.items():
        arr = np.array([v for v in vals if v is not None], dtype=float)
        if arr.size == 0:
            bin_summary[key] = {
                "n_valid_salts": 0,
                "delta_mse_median": None,
                "delta_mse_q10": None,
                "delta_mse_q90": None,
                "reparam_win_rate": None,
            }
            continue
        bin_summary[key] = {
            "n_valid_salts": int(arr.size),
            "delta_mse_median": float(np.median(arr)),
            "delta_mse_q10": float(np.quantile(arr, 0.10)),
            "delta_mse_q90": float(np.quantile(arr, 0.90)),
            "reparam_win_rate": float(np.mean(arr > 0.0)),
        }

    g = np.array(global_delta, dtype=float)
    dataset_summary = {
        "delta_mse_global_median": float(np.median(g)),
        "delta_mse_global_q10": float(np.quantile(g, 0.10)),
        "delta_mse_global_q90": float(np.quantile(g, 0.90)),
        "reparam_global_win_rate": float(np.mean(g > 0.0)),
    }
    return {
        "summary": dataset_summary,
        "bin_summary": bin_summary,
    }


def make_domain_statement(s: Dict[str, object], min_win: float = 0.70) -> Dict[str, object]:
    bins_reparam = []
    bins_hd = []
    uncertain = []
    for key, row in s["bin_summary"].items():
        wr = row["reparam_win_rate"]
        if wr is None:
            uncertain.append(key)
        elif wr >= min_win:
            bins_reparam.append(key)
        elif wr <= (1.0 - min_win):
            bins_hd.append(key)
        else:
            uncertain.append(key)
    return {
        "reparam_dominant_bins": bins_reparam,
        "hd_dominant_bins": bins_hd,
        "uncertain_bins": uncertain,
        "global_reparam_better": bool(s["summary"]["delta_mse_global_median"] > 0.0),
    }


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
    args = ap.parse_args()

    triad = load_reparam_triad()
    p_primary = (ROOT / args.primary_pta).resolve()
    p_stress = (ROOT / args.stress_pta).resolve()
    if not p_primary.exists():
        raise RuntimeError(f"Primary PTA not found: {p_primary}")
    if not p_stress.exists():
        raise RuntimeError(f"Stress PTA not found: {p_stress}")

    salts = list(range(int(args.n_salts)))
    d_primary = pd.read_csv(p_primary)
    d_stress = pd.read_csv(p_stress)

    primary = evaluate_dataset(d_primary, triad=triad, salts=salts)
    stress = evaluate_dataset(d_stress, triad=triad, salts=salts)

    domain_primary = make_domain_statement(primary, min_win=0.70)
    domain_stress = make_domain_statement(stress, min_win=0.70)

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "triad_source": "report_qw1932_physical_reparameterization_eta_scan.json:selected",
        "triad": triad,
        "config": {
            "n_salts": int(args.n_salts),
            "dominance_win_rate_threshold": 0.70,
        },
        "datasets": {
            "primary": {
                "path": str(p_primary),
                **primary,
                "domain_statement": domain_primary,
            },
            "stress": {
                "path": str(p_stress),
                **stress,
                "domain_statement": domain_stress,
            },
        },
        "verdict": "DOMAIN_OF_VALIDITY_MAP_COMPLETE",
        "required_next_step": "FREEZE_DOMAIN_CONDITIONAL_PREDICTIONS_AND_PREPARE_EXTERNAL_CONFIRMATORY_PROTOCOL",
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    ps = primary["summary"]
    ss = stress["summary"]
    lines = [
        "# RAPORT QW-1936: DOMAIN OF VALIDITY MAP",
        "",
        f"- Data UTC: {out['generated_utc']}",
        (
            "- Reparam triad: "
            f"omega={triad['omega']:.6f}, phi={triad['phi']:.6f}, beta={triad['beta']:.6f}, eta={triad['eta']:.2f}"
        ),
        f"- Verdict: **{out['verdict']}**",
        "",
        "## Global Delta MSE (HD - Reparam)",
        f"- Primary median [q10, q90]: {ps['delta_mse_global_median']:.6f} [{ps['delta_mse_global_q10']:.6f}, {ps['delta_mse_global_q90']:.6f}]",
        f"- Primary global win_rate(reparam): {ps['reparam_global_win_rate']:.4f}",
        f"- Stress median [q10, q90]: {ss['delta_mse_global_median']:.6f} [{ss['delta_mse_global_q10']:.6f}, {ss['delta_mse_global_q90']:.6f}]",
        f"- Stress global win_rate(reparam): {ss['reparam_global_win_rate']:.4f}",
        "",
        "## Domain Statement: Primary",
        f"- reparam_dominant_bins: {domain_primary['reparam_dominant_bins']}",
        f"- hd_dominant_bins: {domain_primary['hd_dominant_bins']}",
        f"- uncertain_bins: {domain_primary['uncertain_bins']}",
        "",
        "## Domain Statement: Stress",
        f"- reparam_dominant_bins: {domain_stress['reparam_dominant_bins']}",
        f"- hd_dominant_bins: {domain_stress['hd_dominant_bins']}",
        f"- uncertain_bins: {domain_stress['uncertain_bins']}",
        "",
        "## Required Next Step",
        f"- {out['required_next_step']}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1936] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1936] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1936] verdict={out['verdict']}")


if __name__ == "__main__":
    main()

