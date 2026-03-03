#!/usr/bin/env python3
"""
QW-1922: Selected beta observable blind external intervention run.

Implements next step from QW-1921:
IMPLEMENT_SELECTED_BETA_OBSERVABLE_IN_BLIND_EXTERNAL_INTERVENTION_RUN

Selected observable:
- B7_local_resid_std (from QW-1921 discovery/holdout pass)

Datasets:
- primary external rebuild v2
- stress external alpha6

Evaluation:
- deterministic discovery/holdout split by pair_id hash,
- beta-like intervention vs omega-like intervention,
- contrast = beta_effect - omega_leakage,
- bootstrap lower bound on holdout contrast.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1922_beta_observable_blind_external_intervention.json"
OUT_MD = ROOT / "RAPORT_QW1922_BETA_OBSERVABLE_BLIND_EXTERNAL_INTERVENTION.md"


def split_index(pair_id: str, k: int = 2) -> int:
    h = hashlib.sha256(pair_id.encode("utf-8")).hexdigest()
    return int(h[-8:], 16) % k


def moving_average(x: np.ndarray, win: int) -> np.ndarray:
    w = max(3, int(win) | 1)
    pad = w // 2
    xp = np.pad(x, (pad, pad), mode="edge")
    ker = np.ones(w, dtype=float) / float(w)
    y = np.convolve(xp, ker, mode="same")
    return y[pad : pad + len(x)]


def robust_scale(x: np.ndarray) -> float:
    med = float(np.median(x))
    mad = float(np.median(np.abs(x - med)))
    return max(1e-9, 1.4826 * mad)


def recompute_components(theta: np.ndarray, hxy: np.ndarray) -> Dict[str, np.ndarray]:
    order = np.argsort(theta)
    y = hxy[order]
    trend = moving_average(y, win=41)
    resid = y - trend

    abs_resid = np.abs(resid)
    c_resid_std_local = moving_average(abs_resid, win=31)

    inv = np.zeros_like(y)
    inv[order] = np.arange(len(y))
    inv = inv.astype(int)

    return {
        "B7_local_resid_std": c_resid_std_local[inv],
        "trend": trend,
        "resid": resid,
        "order": order,
    }


def interventions(theta: np.ndarray, hxy: np.ndarray) -> Dict[str, np.ndarray]:
    t = (theta - float(np.min(theta))) / max(float(np.max(theta) - np.min(theta)), 1e-12)

    # beta-like attenuation toward large theta.
    hxy_beta = np.clip(hxy / (1.0 + 0.80 * t), 0.0, 1.0)

    # omega-like phase scramble preserving smooth trend.
    c = recompute_components(theta, hxy)
    order = c["order"]
    trend = c["trend"]
    resid = c["resid"]
    shift = max(5, int(0.09 * len(resid)))
    resid_shift = np.roll(resid, shift)
    hxy_omega_ord = np.clip(trend + resid_shift, 0.0, 1.0)
    hxy_omega = np.zeros_like(hxy_omega_ord)
    hxy_omega[order] = hxy_omega_ord

    return {
        "beta": hxy_beta,
        "omega": hxy_omega,
    }


def effect(base: np.ndarray, alt: np.ndarray) -> float:
    s = robust_scale(base)
    return float(abs(np.median(alt) - np.median(base)) / s)


def bootstrap_contrast(base: np.ndarray, beta_alt: np.ndarray, omega_alt: np.ndarray, n_boot: int, seed: int) -> Dict[str, float]:
    rng = np.random.default_rng(seed)
    n = len(base)
    vals = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        b0 = base[idx]
        bb = beta_alt[idx]
        bo = omega_alt[idx]
        e_beta = effect(b0, bb)
        e_omega = effect(b0, bo)
        vals.append(e_beta - e_omega)

    arr = np.array(vals, dtype=float)
    return {
        "median": float(np.median(arr)),
        "q05": float(np.quantile(arr, 0.05)),
        "q95": float(np.quantile(arr, 0.95)),
    }


def eval_dataset(df: pd.DataFrame, n_boot: int, seed: int) -> Dict[str, object]:
    req = ["pair_id", "theta_deg", "hxy"]
    miss = [c for c in req if c not in df.columns]
    if miss:
        raise RuntimeError(f"Missing columns: {miss}")

    split = np.array([split_index(str(x), k=2) for x in df["pair_id"].astype(str)], dtype=int)
    disc = df.loc[split == 0].reset_index(drop=True)
    hold = df.loc[split == 1].reset_index(drop=True)

    def one(part: pd.DataFrame, offset: int) -> Dict[str, object]:
        theta = part["theta_deg"].to_numpy(dtype=float)
        hxy = part["hxy"].to_numpy(dtype=float)

        iv = interventions(theta, hxy)

        b0 = recompute_components(theta, hxy)["B7_local_resid_std"]
        bb = recompute_components(theta, iv["beta"])["B7_local_resid_std"]
        bo = recompute_components(theta, iv["omega"])["B7_local_resid_std"]

        e_beta = effect(b0, bb)
        e_omega = effect(b0, bo)
        contrast = float(e_beta - e_omega)

        ci = bootstrap_contrast(b0, bb, bo, n_boot=n_boot, seed=seed + offset)

        return {
            "n": int(len(part)),
            "effect_beta": e_beta,
            "effect_omega": e_omega,
            "contrast": contrast,
            "contrast_bootstrap": ci,
        }

    d_res = one(disc, offset=1)
    h_res = one(hold, offset=2)

    flags = {
        "holdout_effect_beta_ge_0p60": bool(h_res["effect_beta"] >= 0.60),
        "holdout_effect_omega_le_0p25": bool(h_res["effect_omega"] <= 0.25),
        "holdout_contrast_ge_0p35": bool(h_res["contrast"] >= 0.35),
        "holdout_contrast_boot_q05_ge_0p20": bool(h_res["contrast_bootstrap"]["q05"] >= 0.20),
    }

    return {
        "split": {
            "n_discovery": int(len(disc)),
            "n_holdout": int(len(hold)),
        },
        "discovery": d_res,
        "holdout": h_res,
        "pass_flags": flags,
        "all_pass": bool(all(flags.values())),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--primary-pta",
        type=str,
        default="external_confirmatory_v2/confirmatory_dataset_external_source_rebuild_v2_1831cfg/pta_v2_pairs.csv",
    )
    ap.add_argument(
        "--stress-pta",
        type=str,
        default="external_confirmatory_v2/confirmatory_dataset_external_source_alpha6_1831cfg/pta_v2_pairs.csv",
    )
    ap.add_argument("--n-boot", type=int, default=1200)
    args = ap.parse_args()

    p_primary = (ROOT / args.primary_pta).resolve()
    p_stress = (ROOT / args.stress_pta).resolve()

    if not p_primary.exists():
        raise RuntimeError(f"Primary PTA not found: {p_primary}")
    if not p_stress.exists():
        raise RuntimeError(f"Stress PTA not found: {p_stress}")

    d_primary = pd.read_csv(p_primary)
    d_stress = pd.read_csv(p_stress)

    primary = eval_dataset(d_primary, n_boot=int(args.n_boot), seed=192201)
    stress = eval_dataset(d_stress, n_boot=int(args.n_boot), seed=192202)

    primary_pass = bool(primary["all_pass"])
    stress_pass = bool(stress["all_pass"])

    if primary_pass and stress_pass:
        verdict = "BETA_OBSERVABLE_BLIND_EXTERNAL_INTERVENTION_PASS"
    elif primary_pass:
        verdict = "BETA_OBSERVABLE_BLIND_EXTERNAL_INTERVENTION_PARTIAL"
    else:
        verdict = "BETA_OBSERVABLE_BLIND_EXTERNAL_INTERVENTION_FAIL"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "selected_observable": "B7_local_resid_std",
        "datasets": {
            "primary": {
                "path": str(p_primary),
                **primary,
            },
            "stress": {
                "path": str(p_stress),
                **stress,
            },
        },
        "pass_flags": {
            "primary_all_pass": primary_pass,
            "stress_all_pass": stress_pass,
        },
        "verdict": verdict,
        "required_next_step": (
            "INTEGRATE_BETA_OBSERVABLE_AS_EXPLICIT_CONSTRAINT_IN_TRIAD_DERIVATION"
            if verdict != "BETA_OBSERVABLE_BLIND_EXTERNAL_INTERVENTION_FAIL"
            else "REWORK_INTERVENTION_DEFINITION_AND_REPEAT"
        ),
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    p_h = primary["holdout"]
    s_h = stress["holdout"]

    lines = [
        "# RAPORT QW-1922: BETA OBSERVABLE BLIND EXTERNAL INTERVENTION",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Selected observable: `{out['selected_observable']}`",
        f"- Verdict: **{verdict}**",
        "",
        "## Primary (holdout)",
        f"- effect_beta: {p_h['effect_beta']:.4f}",
        f"- effect_omega: {p_h['effect_omega']:.4f}",
        f"- contrast: {p_h['contrast']:.4f}",
        f"- contrast_boot_q05: {p_h['contrast_bootstrap']['q05']:.4f}",
        f"- all_pass: {primary['all_pass']}",
        "",
        "## Stress (holdout)",
        f"- effect_beta: {s_h['effect_beta']:.4f}",
        f"- effect_omega: {s_h['effect_omega']:.4f}",
        f"- contrast: {s_h['contrast']:.4f}",
        f"- contrast_boot_q05: {s_h['contrast_bootstrap']['q05']:.4f}",
        f"- all_pass: {stress['all_pass']}",
        "",
        "## Global Flags",
        f"- primary_all_pass: {primary_pass}",
        f"- stress_all_pass: {stress_pass}",
        "",
        "## Required Next Step",
        f"- {out['required_next_step']}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1922] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1922] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1922] verdict={verdict}")


if __name__ == "__main__":
    main()
