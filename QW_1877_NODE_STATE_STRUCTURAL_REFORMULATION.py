#!/usr/bin/env python3
"""
QW-1877: Structural reformulation with explicit node-state dynamics.

Introduces a separate node-state recurrence instead of static node amplitude
modulation and tests fit/compatibility on QW-1739 profiles.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1877_node_state_structural_reformulation.json"
OUT_MD = ROOT / "RAPORT_QW1877_NODE_STATE_STRUCTURAL_REFORMULATION.md"


def read_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def circular_diff(a: float, b: float) -> float:
    d = (a - b + math.pi) % (2.0 * math.pi) - math.pi
    return abs(d)


def canon_score(omega: float, phi: float, beta: float, t: Dict[str, float]) -> float:
    z_o = abs(omega - t["omega"]) / 0.20
    z_p = circular_diff(phi, t["phi"]) / 0.30
    z_b = abs(beta - t["beta"]) / 0.030
    z = min(20.0, math.sqrt((z_o * z_o + z_p * z_p + z_b * z_b) / 3.0))
    return float(math.exp(-0.5 * z * z))


def node_indicator(d: np.ndarray, nodes: List[int]) -> np.ndarray:
    s = np.zeros_like(d, dtype=float)
    node_set = set(int(x) for x in nodes)
    for i, x in enumerate(d):
        s[i] = 1.0 if int(x) in node_set else 0.0
    return s


def simulate_node_state(d: np.ndarray, omega: float, phi: float, rho: float, xi: float, zeta: float, nodes: List[int]) -> np.ndarray:
    ind = node_indicator(d, nodes)
    st = np.zeros_like(d, dtype=float)
    prev = 0.0
    for i, di in enumerate(d):
        drive = xi * math.sin(omega * float(di) + phi)
        node_term = -zeta * ind[i]
        cur = rho * prev + drive + node_term
        st[i] = cur
        prev = cur
    return st


def model_signal(d: np.ndarray, omega: float, phi: float, beta: float, amp: float, gamma: float, rho: float, xi: float, zeta: float, nodes: List[int]) -> np.ndarray:
    env = 1.0 / (1.0 + beta * d)
    base = amp * np.cos(omega * d + phi) * env
    st = simulate_node_state(d, omega, phi, rho, xi, zeta, nodes)
    return base + gamma * st


def estimate_amp_gamma(y: np.ndarray, d: np.ndarray, omega: float, phi: float, beta: float, rho: float, xi: float, zeta: float, nodes: List[int]) -> Dict[str, float]:
    env_basis = np.cos(omega * d + phi) / (1.0 + beta * d)
    st_basis = simulate_node_state(d, omega, phi, rho, xi, zeta, nodes)

    X = np.column_stack([env_basis, st_basis])
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    amp = float(coef[0])
    gamma = float(coef[1])
    yhat = X @ coef

    rmse = float(np.sqrt(np.mean((y - yhat) ** 2)))
    ss_res = float(np.sum((y - yhat) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 if ss_tot <= 1e-15 else float(1.0 - ss_res / ss_tot)

    return {
        "amp": amp,
        "gamma": gamma,
        "rmse": rmse,
        "r2": r2,
        "yhat": yhat,
    }


def main() -> None:
    d1739 = read_json("report_qw1739_signed_dynamic_micromodel_derivation.json")
    d1862 = read_json("report_qw1862_micromodel_canonical_compatibility_audit.json")
    d1871 = read_json("report_qw1871_primary_node_evidence_corpus.json")
    d1874 = read_json("report_qw1874_beta_omega_orthogonal_forcing.json")

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

    models = {
        "M_A": {"nodes": [2, 5, 8, 11], "prior": pA},
        "M_B": {"nodes": [2, 8, 14], "prior": pB},
        "M_MIX": {"nodes": [2, 5, 8, 11, 14], "prior": pM},
    }

    base_rows = {int(x["seed"]): x for x in d1874.get("rows", [])}

    rng = np.random.default_rng(187700)

    rows_out = []
    for r in d1739.get("rows", []):
        seed = int(r["seed"])
        y = np.array([float(r["profile"][str(i)]) for i in range(1, 13)], dtype=float)
        d = np.arange(1, len(y) + 1, dtype=float)

        base = base_rows.get(seed, {}).get("best", None)
        if base is None:
            base = {
                "omega": float(r["omega_hat"]),
                "phi": float(r["phi_hat"]),
                "beta": float(r["beta_hat"]),
                "rmse": float(r["rmse"]),
                "canon_score": canon_score(float(r["omega_hat"]), float(r["phi_hat"]), float(r["beta_hat"]), t),
            }

        best = {
            "model": "BASELINE",
            "omega": float(base["omega"]),
            "phi": float(base["phi"]),
            "beta": float(base["beta"]),
            "rho": 0.0,
            "xi": 0.0,
            "zeta": 0.0,
            "amp": 0.0,
            "gamma": 0.0,
            "rmse": float(base["rmse"]),
            "r2": 0.0,
            "canon_score": float(base["canon_score"]),
            "objective": float(base["rmse"] + 0.15 * (1.0 - float(base["canon_score"]))),
        }

        for mid, md in models.items():
            nodes = [n for n in md["nodes"] if n <= len(y) + 2]
            prior = max(1e-6, md["prior"])

            # Canonical anchors + shifted anchors
            anchors = [
                (0.78, t["phi"], 0.02, 0.60, 0.25, 0.30),
                (0.20, 0.48, 0.25, 0.70, 0.12, 0.15),
                (float(base["omega"]), float(base["phi"]), float(base["beta"]), 0.55, 0.20, 0.20),
            ]
            for om, ph, be, rho, xi, zeta in anchors:
                est = estimate_amp_gamma(y, d, om, ph, be, rho, xi, zeta, nodes)
                cscore = canon_score(om, ph, be, t)
                obj = est["rmse"] + 0.15 * (1.0 - cscore) + 0.03 * (-math.log(prior))
                if obj < best["objective"]:
                    best = {
                        "model": mid,
                        "omega": om,
                        "phi": ph,
                        "beta": be,
                        "rho": rho,
                        "xi": xi,
                        "zeta": zeta,
                        "amp": est["amp"],
                        "gamma": est["gamma"],
                        "rmse": est["rmse"],
                        "r2": est["r2"],
                        "canon_score": cscore,
                        "objective": obj,
                    }

            for _ in range(420):
                coin = rng.random()
                if coin < 0.45:
                    om = float(rng.uniform(0.60, 0.95))
                    be = float(rng.uniform(0.005, 0.06))
                elif coin < 0.80:
                    om = float(rng.uniform(0.08, 0.30))
                    be = float(rng.uniform(0.10, 0.30))
                else:
                    om = float(rng.uniform(0.08, 0.95))
                    be = float(rng.uniform(0.005, 0.30))

                ph = float(rng.uniform(-1.1, 1.1))
                rho = float(rng.uniform(0.30, 0.95))
                xi = float(rng.uniform(0.02, 0.55))
                zeta = float(rng.uniform(0.02, 0.60))

                est = estimate_amp_gamma(y, d, om, ph, be, rho, xi, zeta, nodes)
                cscore = canon_score(om, ph, be, t)
                obj = est["rmse"] + 0.15 * (1.0 - cscore) + 0.03 * (-math.log(prior))

                if obj < best["objective"]:
                    best = {
                        "model": mid,
                        "omega": om,
                        "phi": ph,
                        "beta": be,
                        "rho": rho,
                        "xi": xi,
                        "zeta": zeta,
                        "amp": est["amp"],
                        "gamma": est["gamma"],
                        "rmse": est["rmse"],
                        "r2": est["r2"],
                        "canon_score": cscore,
                        "objective": obj,
                    }

        rows_out.append(
            {
                "seed": seed,
                "baseline": {
                    "omega": float(base["omega"]),
                    "phi": float(base["phi"]),
                    "beta": float(base["beta"]),
                    "rmse": float(base["rmse"]),
                    "canon_score": float(base["canon_score"]),
                },
                "reform_best": best,
                "gains": {
                    "rmse_gain": float(base["rmse"]) - best["rmse"],
                    "canon_gain": best["canon_score"] - float(base["canon_score"]),
                },
            }
        )

    rmse_b = np.array([x["baseline"]["rmse"] for x in rows_out], dtype=float)
    rmse_r = np.array([x["reform_best"]["rmse"] for x in rows_out], dtype=float)
    cs_b = np.array([x["baseline"]["canon_score"] for x in rows_out], dtype=float)
    cs_r = np.array([x["reform_best"]["canon_score"] for x in rows_out], dtype=float)

    cnt = {}
    for x in rows_out:
        m = x["reform_best"]["model"]
        cnt[m] = cnt.get(m, 0) + 1

    summary = {
        "n_profiles": len(rows_out),
        "baseline_rmse_median": float(np.median(rmse_b)),
        "reform_rmse_median": float(np.median(rmse_r)),
        "baseline_canon_median": float(np.median(cs_b)),
        "reform_canon_median": float(np.median(cs_r)),
        "rmse_improved_fraction": float(np.mean(rmse_r < rmse_b)),
        "canon_improved_fraction": float(np.mean(cs_r > cs_b)),
        "canon_gain_median": float(np.median(cs_r - cs_b)),
        "best_model_counts": cnt,
    }

    if (
        summary["rmse_improved_fraction"] >= 0.60
        and summary["canon_improved_fraction"] >= 0.60
        and summary["canon_gain_median"] >= 1e-4
    ):
        verdict = "NODE_STATE_REFORMULATION_PARTIAL_PROGRESS"
    elif summary["canon_gain_median"] > 0:
        verdict = "NODE_STATE_REFORMULATION_WEAK_PROGRESS"
    else:
        verdict = "NODE_STATE_REFORMULATION_NO_PROGRESS"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "target_tuple": t,
        "node_model_prior": {"A": pA, "B": pB, "MIX": pM},
        "summary": summary,
        "rows": rows_out,
        "verdict": verdict,
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1877: NODE-STATE STRUCTURAL REFORMULATION",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        "",
        "## Summary",
        f"- n_profiles: {summary['n_profiles']}",
        f"- rmse median: {summary['baseline_rmse_median']:.4f} -> {summary['reform_rmse_median']:.4f}",
        f"- canon median: {summary['baseline_canon_median']:.4e} -> {summary['reform_canon_median']:.4e}",
        f"- rmse improved fraction: {summary['rmse_improved_fraction']:.3f}",
        f"- canon improved fraction: {summary['canon_improved_fraction']:.3f}",
        f"- canon gain median: {summary['canon_gain_median']:.4e}",
        "",
        "## Best Model Counts",
    ]

    for k, v in sorted(cnt.items(), key=lambda x: x[1], reverse=True):
        lines.append(f"- {k}: {v}")

    lines += ["", "## Artifacts", f"- JSON: `{OUT_JSON.name}`"]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1877] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1877] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
