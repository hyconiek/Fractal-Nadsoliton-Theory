#!/usr/bin/env python3
"""
QW-1720: Extended flavor operator closure test.

Goal:
1) Keep unitarity exactly.
2) Use richer topological invariants than QW-1719.
3) Test single shared parameter family for CKM + PMNS.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.linalg import expm


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1720_flavor_extended_operator.json"
OUT_MD = ROOT / "RAPORT_QW1720_FLAVOR_EXTENDED_OPERATOR.md"

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


def matrix_error(pred_abs: np.ndarray, ref: np.ndarray) -> dict:
    abs_err = np.abs(pred_abs - ref)
    rel = abs_err / np.clip(ref, 1e-12, None)
    return {
        "mae": float(np.mean(abs_err)),
        "max_abs": float(np.max(abs_err)),
        "mean_rel_pct": float(np.mean(rel) * 100.0),
        "max_rel_pct": float(np.max(rel) * 100.0),
    }


def antisym_part(m: np.ndarray) -> np.ndarray:
    return 0.5 * (m - m.T)


def build_generators(q_left: np.ndarray, q_right: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Return three antisymmetric generator templates:
    G1: inverse distance
    G2: inverse squared distance
    G3: generation-locality mask weighted by distance asymmetry
    """
    n = len(q_left)
    g1 = np.zeros((n, n), dtype=float)
    g2 = np.zeros((n, n), dtype=float)
    g3 = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            dq = abs(q_left[i] - q_right[j])
            gap = abs(i - j)
            sign = 1.0 if i < j else -1.0
            g1[i, j] = sign / max(dq, 1.0)
            g2[i, j] = sign / max(dq * dq, 1.0)
            locality = 1.0 / (1.0 + gap)
            asym = (q_left[i] - q_right[j]) / (1.0 + dq)
            g3[i, j] = sign * locality * asym
    return antisym_part(g1), antisym_part(g2), antisym_part(g3)


def build_unitary(
    g1: np.ndarray,
    g2: np.ndarray,
    g3: np.ndarray,
    params: dict,
) -> np.ndarray:
    a1 = params["a1"]
    a2 = params["a2"]
    a3 = params["a3"]
    a4 = params["a4"]
    delta = params["delta"]

    # Hermitian generator from antisymmetric real templates
    # H = i * (linear combination of antisymmetric matrices) + CP-phase rank-1 structure
    base = a1 * g1 + a2 * g2 + a3 * g3 + a4 * (g1 @ g2 - g2 @ g1)
    h = 1j * base

    # Add minimal CP-odd Hermitian term that preserves unitarity on exponentiation
    cp = np.zeros_like(h, dtype=complex)
    cp[0, 2] = np.exp(-1j * delta)
    cp[2, 0] = np.conj(cp[0, 2])
    h = h + 0.05 * cp

    return expm(1j * h)


def clamp(value: float, lo: float, hi: float) -> float:
    return float(min(max(value, lo), hi))


def evaluate_shared(
    params: dict,
    ckm_gens: tuple[np.ndarray, np.ndarray, np.ndarray],
    pmns_gens: tuple[np.ndarray, np.ndarray, np.ndarray],
) -> dict:
    u_ckm = np.abs(build_unitary(*ckm_gens, params))
    u_pmns = np.abs(build_unitary(*pmns_gens, params))
    e_ckm = matrix_error(u_ckm, CKM_REF)
    e_pmns = matrix_error(u_pmns, PMNS_REF)
    reg = 0.03 * (
        abs(params["a1"])
        + abs(params["a2"])
        + abs(params["a3"])
        + abs(params["a4"])
    )
    score = (
        e_ckm["mean_rel_pct"]
        + e_pmns["mean_rel_pct"]
        + 0.10 * (e_ckm["max_rel_pct"] + e_pmns["max_rel_pct"])
        + reg
    )
    return {
        "score": float(score),
        "params": dict(params),
        "ckm_error": e_ckm,
        "pmns_error": e_pmns,
        "ckm_pred_abs": u_ckm,
        "pmns_pred_abs": u_pmns,
    }


def search_shared() -> dict:
    q_ckm_l = np.array([0.0, 9.0, 14.0], dtype=float)
    q_ckm_r = np.array([7.0, 9.0, 14.0], dtype=float)
    q_pmns_l = np.array([0.0, 1.0, 2.0], dtype=float)
    q_pmns_r = np.array([24.0, 14.0, 9.0], dtype=float)

    ckm_gens = build_generators(q_ckm_l, q_ckm_r)
    pmns_gens = build_generators(q_pmns_l, q_pmns_r)

    bounds = {
        "a1": (-2.0, 2.0),
        "a2": (-1.2, 1.2),
        "a3": (-1.2, 1.2),
        "a4": (-0.5, 0.5),
        "delta": (0.0, math.pi),
    }
    coarse_grids = {
        "a1": np.linspace(-2.0, 2.0, 9),
        "a2": np.linspace(-1.2, 1.2, 9),
        "a3": np.linspace(-1.2, 1.2, 9),
        "a4": np.linspace(-0.5, 0.5, 7),
        "delta": np.linspace(0.0, math.pi, 9),
    }

    best = None
    eval_count = 0

    # Coarse global scan
    for a1 in coarse_grids["a1"]:
        for a2 in coarse_grids["a2"]:
            for a3 in coarse_grids["a3"]:
                for a4 in coarse_grids["a4"]:
                    for delta in coarse_grids["delta"]:
                        params = {
                            "a1": float(a1),
                            "a2": float(a2),
                            "a3": float(a3),
                            "a4": float(a4),
                            "delta": float(delta),
                        }
                        cand = evaluate_shared(params, ckm_gens, pmns_gens)
                        eval_count += 1
                        if best is None or cand["score"] < best["score"]:
                            best = cand

    # Deterministic coordinate refinement around best coarse point
    steps = {
        "a1": 0.35,
        "a2": 0.22,
        "a3": 0.22,
        "a4": 0.09,
        "delta": 0.28,
    }
    keys = ["a1", "a2", "a3", "a4", "delta"]
    for shrink in [1.0, 0.65, 0.42, 0.28, 0.18]:
        improved = True
        while improved:
            improved = False
            center = best["params"]
            for k in keys:
                for sign in (-1.0, 1.0):
                    trial = dict(center)
                    lo, hi = bounds[k]
                    trial[k] = clamp(center[k] + sign * steps[k] * shrink, lo, hi)
                    cand = evaluate_shared(trial, ckm_gens, pmns_gens)
                    eval_count += 1
                    if cand["score"] + 1e-12 < best["score"]:
                        best = cand
                        improved = True

    best["eval_count"] = int(eval_count)
    return best


def main() -> None:
    best = search_shared()
    threshold = 15.0
    passed = (
        best["ckm_error"]["mean_rel_pct"] < threshold
        and best["pmns_error"]["mean_rel_pct"] < threshold
    )
    verdict = "FLAVOR_EXTENDED_OPERATOR_CLOSED" if passed else "FLAVOR_EXTENDED_OPERATOR_NOT_CLOSED"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "threshold_mean_rel_pct": threshold,
        "search_protocol": {
            "coarse_grid_sizes": {"a1": 9, "a2": 9, "a3": 9, "a4": 7, "delta": 9},
            "refinement": "deterministic_coordinate_descent",
        },
        "best": {
            "score": best["score"],
            "eval_count": best["eval_count"],
            "params": best["params"],
            "ckm_error": best["ckm_error"],
            "pmns_error": best["pmns_error"],
            "ckm_pred_abs": best["ckm_pred_abs"].tolist(),
            "pmns_pred_abs": best["pmns_pred_abs"].tolist(),
        },
        "verdict": verdict,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    p = best["params"]
    lines = [
        "# RAPORT QW-1720: FLAVOR EXTENDED OPERATOR",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Werdykt: **{verdict}**",
        "",
        "## Best shared params",
        f"- a1={p['a1']:.4f}, a2={p['a2']:.4f}, a3={p['a3']:.4f}, a4={p['a4']:.4f}, delta={p['delta']:.4f}",
        f"- CKM mean rel. error: {best['ckm_error']['mean_rel_pct']:.2f}%",
        f"- PMNS mean rel. error: {best['pmns_error']['mean_rel_pct']:.2f}%",
        f"- Liczba ewaluacji: {best['eval_count']}",
        "",
        "## Artefakty",
        f"- JSON szczegółowy: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"[QW-1720] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1720] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
