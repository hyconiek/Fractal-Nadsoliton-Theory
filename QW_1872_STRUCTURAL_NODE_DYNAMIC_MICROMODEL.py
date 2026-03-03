#!/usr/bin/env python3
"""
QW-1872: Structural node-dynamic micromodel audit.

Compares baseline micromodel fit versus node-augmented structural model on
profiles from QW-1739 and checks compatibility with canonical tuple.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1872_structural_node_dynamic_micromodel.json"
OUT_MD = ROOT / "RAPORT_QW1872_STRUCTURAL_NODE_DYNAMIC_MICROMODEL.md"


def read_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def circular_diff(a: float, b: float) -> float:
    d = (a - b + math.pi) % (2.0 * math.pi) - math.pi
    return abs(d)


def estimate_beta_from_envelope(y: np.ndarray) -> float:
    d = np.arange(1, len(y) + 1, dtype=float)
    a = np.abs(y)
    q = float(np.quantile(a, 0.35))
    mask = a > max(q, 1e-4)
    x = d[mask]
    z = 1.0 / np.clip(a[mask], 1e-6, None)

    beta = np.nan
    if len(x) >= 3:
        m, c = np.polyfit(x, z, deg=1)
        if c > 0 and m > 0:
            beta = float(m / c)

    if not np.isfinite(beta) or beta <= 0:
        beta_grid = np.linspace(0.001, 0.30, 1200)
        best_rmse = float("inf")
        best_beta = 0.05
        for b in beta_grid:
            env = 1.0 / (1.0 + b * d)
            scale = float(np.dot(a, env) / max(np.dot(env, env), 1e-12))
            pred = scale * env
            rmse = float(np.sqrt(np.mean((a - pred) ** 2)))
            if rmse < best_rmse:
                best_rmse = rmse
                best_beta = float(b)
        beta = best_beta

    return float(np.clip(beta, 0.001, 0.30))


def fit_phase_with_fixed_omega(y_flat: np.ndarray, d: np.ndarray, omega: float) -> Tuple[float, float, float]:
    c = np.cos(omega * d)
    s = np.sin(omega * d)
    b = np.column_stack([c, s])
    coeff, *_ = np.linalg.lstsq(b, y_flat, rcond=None)
    c0, s0 = float(coeff[0]), float(coeff[1])
    amp = float(np.hypot(c0, s0))
    phi = float(np.arctan2(-s0, c0))
    pred = amp * np.cos(omega * d + phi)
    sse = float(np.sum((y_flat - pred) ** 2))
    return amp, phi, sse


def derive_omega_phi_beta(y: np.ndarray) -> Dict[str, float]:
    d = np.arange(1, len(y) + 1, dtype=float)
    beta = estimate_beta_from_envelope(y)
    y_flat = y * (1.0 + beta * d)
    y_ctr = y_flat - np.mean(y_flat)

    n_fft = 256
    spec = np.fft.rfft(y_ctr, n=n_fft)
    freqs = np.fft.rfftfreq(n_fft, d=1.0)
    power = np.abs(spec) ** 2
    if len(power) > 1:
        idx = int(np.argmax(power[1:]) + 1)
        omega0 = float(2.0 * math.pi * freqs[idx])
    else:
        omega0 = 0.5
    omega0 = float(np.clip(omega0, 0.05, 1.8))

    grid = np.linspace(max(0.05, omega0 - 0.28), min(1.8, omega0 + 0.28), 621)
    best = {"omega": omega0, "amp": 0.0, "phi": 0.0, "sse": float("inf")}
    for om in grid:
        amp, phi, sse = fit_phase_with_fixed_omega(y_flat, d, float(om))
        if sse < best["sse"]:
            best = {"omega": float(om), "amp": amp, "phi": float(phi), "sse": sse}

    omega = best["omega"]
    amp = best["amp"]
    phi = best["phi"]

    y_hat = amp * np.cos(omega * d + phi) / (1.0 + beta * d)
    ss_res = float(np.sum((y - y_hat) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 if ss_tot <= 1e-15 else float(1.0 - ss_res / ss_tot)
    rmse = float(np.sqrt(np.mean((y - y_hat) ** 2)))

    return {
        "omega": omega,
        "phi": phi,
        "beta": beta,
        "amp": amp,
        "rmse": rmse,
        "r2": r2,
        "y_hat": y_hat,
    }


def node_field(d: np.ndarray, nodes: List[int], sigma: float) -> np.ndarray:
    if not nodes:
        return np.zeros_like(d, dtype=float)
    terms = []
    for n in nodes:
        terms.append(np.exp(-0.5 * ((d - float(n)) / sigma) ** 2))
    m = np.max(np.vstack(terms), axis=0)
    return np.clip(m, 0.0, 1.0)


def canonical_score(omega: float, phi: float, beta: float, t: Dict[str, float]) -> float:
    z_o = abs(omega - t["omega"]) / 0.20
    z_p = circular_diff(phi, t["phi"]) / 0.30
    z_b = abs(beta - t["beta"]) / 0.030
    z = min(20.0, math.sqrt((z_o * z_o + z_p * z_p + z_b * z_b) / 3.0))
    return float(math.exp(-0.5 * z * z))


def main() -> None:
    d1739 = read_json("report_qw1739_signed_dynamic_micromodel_derivation.json")
    d1862 = read_json("report_qw1862_micromodel_canonical_compatibility_audit.json")
    d1871 = read_json("report_qw1871_primary_node_evidence_corpus.json")

    rows = d1739.get("rows", [])
    if not rows:
        raise RuntimeError("No rows in report_qw1739_signed_dynamic_micromodel_derivation.json")

    target = d1862.get("target_tuple", {})
    t = {
        "omega": float(target.get("omega", math.pi / 4.0)),
        "phi": float(target.get("phi", math.pi / 6.0)),
        "beta": float(target.get("beta", 0.01)),
    }

    post = d1871.get("model_posterior", {})
    pA = float(post.get("M_A_2_5_8_11_or_2plus3n", 1.0 / 3.0))
    pB = float(post.get("M_B_2_8_14", 1.0 / 3.0))
    pM = float(post.get("M_MIX_A_and_B", 1.0 / 3.0))

    model_defs = {
        "M_A": {"nodes": [2, 5, 8, 11], "prior": pA},
        "M_B": {"nodes": [2, 8, 14], "prior": pB},
        "M_MIX": {"nodes": [2, 5, 8, 11, 14], "prior": pM},
    }

    kappa_grid = [0.10, 0.20, 0.35, 0.50, 0.70, 0.90]
    sigma_grid = [0.35, 0.55, 0.80, 1.10, 1.50]

    out_rows = []

    for r in rows:
        y = np.array([float(r["profile"][str(i)]) for i in range(1, 13)], dtype=float)
        d = np.arange(1, len(y) + 1, dtype=float)

        base = derive_omega_phi_beta(y)
        base_score = canonical_score(base["omega"], base["phi"], base["beta"], t)

        best = {
            "model": "BASELINE",
            "kappa": 0.0,
            "sigma": None,
            "omega": base["omega"],
            "phi": base["phi"],
            "beta": base["beta"],
            "amp": base["amp"],
            "rmse": base["rmse"],
            "r2": base["r2"],
            "canon_score": base_score,
            "objective": base["rmse"] + 0.18 * (1.0 - base_score),
        }

        for mid, md in model_defs.items():
            nodes = [n for n in md["nodes"] if n <= len(y) + 2]
            prior = max(1e-6, md["prior"])

            for kappa in kappa_grid:
                for sigma in sigma_grid:
                    nf = node_field(d, nodes, sigma)
                    mod = np.clip(1.0 - kappa * nf, 0.12, 1.0)

                    y_demod = y / mod
                    est = derive_omega_phi_beta(y_demod)

                    y_hat = est["amp"] * np.cos(est["omega"] * d + est["phi"]) / (1.0 + est["beta"] * d)
                    y_hat *= mod

                    rmse = float(np.sqrt(np.mean((y - y_hat) ** 2)))
                    ss_res = float(np.sum((y - y_hat) ** 2))
                    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
                    r2 = 1.0 if ss_tot <= 1e-15 else float(1.0 - ss_res / ss_tot)

                    cscore = canonical_score(est["omega"], est["phi"], est["beta"], t)
                    obj = rmse + 0.18 * (1.0 - cscore) + 0.03 * (-math.log(prior)) + 0.015 * kappa

                    if obj < best["objective"]:
                        best = {
                            "model": mid,
                            "kappa": kappa,
                            "sigma": sigma,
                            "omega": est["omega"],
                            "phi": est["phi"],
                            "beta": est["beta"],
                            "amp": est["amp"],
                            "rmse": rmse,
                            "r2": r2,
                            "canon_score": cscore,
                            "objective": obj,
                        }

        out_rows.append(
            {
                "seed": int(r["seed"]),
                "n_nodes": int(r["n_nodes"]),
                "baseline": {
                    "omega": base["omega"],
                    "phi": base["phi"],
                    "beta": base["beta"],
                    "rmse": base["rmse"],
                    "r2": base["r2"],
                    "canon_score": base_score,
                },
                "best_node_model": best,
                "gains": {
                    "rmse_gain": base["rmse"] - best["rmse"],
                    "canon_gain": best["canon_score"] - base_score,
                },
            }
        )

    rmse_base = np.array([x["baseline"]["rmse"] for x in out_rows], dtype=float)
    rmse_best = np.array([x["best_node_model"]["rmse"] for x in out_rows], dtype=float)
    c_base = np.array([x["baseline"]["canon_score"] for x in out_rows], dtype=float)
    c_best = np.array([x["best_node_model"]["canon_score"] for x in out_rows], dtype=float)

    model_counts: Dict[str, int] = {}
    for x in out_rows:
        m = x["best_node_model"]["model"]
        model_counts[m] = model_counts.get(m, 0) + 1

    summary = {
        "n_profiles": len(out_rows),
        "baseline_rmse_median": float(np.median(rmse_base)),
        "node_rmse_median": float(np.median(rmse_best)),
        "baseline_canon_median": float(np.median(c_base)),
        "node_canon_median": float(np.median(c_best)),
        "rmse_improved_fraction": float(np.mean(rmse_best < rmse_base)),
        "canon_improved_fraction": float(np.mean(c_best > c_base)),
        "median_canon_gain": float(np.median(c_best - c_base)),
        "best_model_counts": model_counts,
    }

    if (
        summary["node_canon_median"] >= 0.20
        and summary["median_canon_gain"] >= 0.05
        and summary["rmse_improved_fraction"] >= 0.60
    ):
        verdict = "STRUCTURAL_NODE_MODEL_PARTIAL_BRIDGE"
    elif summary["median_canon_gain"] > 0.0:
        verdict = "STRUCTURAL_NODE_MODEL_WEAK_BRIDGE"
    else:
        verdict = "STRUCTURAL_NODE_MODEL_NOT_BRIDGING"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "target_tuple": t,
        "node_model_prior": {
            "A": pA,
            "B": pB,
            "MIX": pM,
        },
        "summary": summary,
        "rows": out_rows,
        "verdict": verdict,
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1872: STRUCTURAL NODE-DYNAMIC MICROMODEL",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        "",
        "## Summary",
        f"- n_profiles: {summary['n_profiles']}",
        f"- baseline rmse median: {summary['baseline_rmse_median']:.4f}",
        f"- node rmse median: {summary['node_rmse_median']:.4f}",
        f"- baseline canon median: {summary['baseline_canon_median']:.4e}",
        f"- node canon median: {summary['node_canon_median']:.4e}",
        f"- rmse improved fraction: {summary['rmse_improved_fraction']:.3f}",
        f"- canon improved fraction: {summary['canon_improved_fraction']:.3f}",
        f"- median canon gain: {summary['median_canon_gain']:.4e}",
        "",
        "## Best Model Counts",
    ]

    for k, v in sorted(model_counts.items(), key=lambda x: x[1], reverse=True):
        lines.append(f"- {k}: {v}")

    lines += ["", "## Artifacts", f"- JSON: `{OUT_JSON.name}`"]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1872] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1872] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
