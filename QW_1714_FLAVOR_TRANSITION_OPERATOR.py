#!/usr/bin/env python3
"""
QW-1714: Flavor closure attempt with a structured transition operator.

Idea:
Build a generator from topological distances and map to a unitary matrix:
U = exp(i * G), with G = a1*A + a2*A^2 + a3*A^3
Then compare |U| to CKM/PMNS benchmarks.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.linalg import expm


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1714_flavor_transition_operator.json"
OUT_MD = ROOT / "RAPORT_QW1714_FLAVOR_TRANSITION_OPERATOR.md"

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
    rel_err = abs_err / np.clip(ref, 1e-12, None)
    return {
        "mae": float(np.mean(abs_err)),
        "max_abs": float(np.max(abs_err)),
        "mean_rel_pct": float(np.mean(rel_err) * 100.0),
        "max_rel_pct": float(np.max(rel_err) * 100.0),
    }


def build_generator(q_left: np.ndarray, q_right: np.ndarray) -> np.ndarray:
    n = len(q_left)
    a = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            dq = abs(q_left[i] - q_right[j])
            # signed antisymmetric pattern by generation ordering
            sign = 1.0 if i < j else -1.0
            a[i, j] = sign / max(dq, 1.0)
    # Make exactly antisymmetric real matrix
    a = 0.5 * (a - a.T)
    return a


def transition_unitary(a_mat: np.ndarray, a1: float, a2: float, a3: float) -> np.ndarray:
    g = a1 * a_mat + a2 * (a_mat @ a_mat) + a3 * (a_mat @ a_mat @ a_mat)
    return expm(1j * g)


def search_best(q_left: np.ndarray, q_right: np.ndarray, ref: np.ndarray) -> dict:
    a_mat = build_generator(q_left, q_right)
    best = None
    grid = np.linspace(-1.5, 1.5, 31)
    for a1 in grid:
        for a2 in np.linspace(-0.8, 0.8, 21):
            for a3 in np.linspace(-0.4, 0.4, 11):
                u = transition_unitary(a_mat, float(a1), float(a2), float(a3))
                pred = np.abs(u)
                err = matrix_error(pred, ref)
                reg = 0.10 * (abs(a1) + abs(a2) + abs(a3))
                score = err["mean_rel_pct"] + 0.25 * err["max_rel_pct"] + reg
                if best is None or score < best["score"]:
                    best = {
                        "score": float(score),
                        "a1": float(a1),
                        "a2": float(a2),
                        "a3": float(a3),
                        "pred_abs": pred,
                        "error": err,
                    }
    return best


def main() -> None:
    q_up = np.array([0.0, 9.0, 14.0], dtype=float)
    q_down = np.array([7.0, 9.0, 14.0], dtype=float)
    q_nu = np.array([0.0, 1.0, 2.0], dtype=float)
    q_lep = np.array([24.0, 14.0, 9.0], dtype=float)

    ckm_best = search_best(q_up, q_down, CKM_REF)
    pmns_best = search_best(q_nu, q_lep, PMNS_REF)

    # single shared coefficients test
    a_mat_ckm = build_generator(q_up, q_down)
    a_mat_pmns = build_generator(q_nu, q_lep)
    shared_best = None
    for a1 in np.linspace(-1.5, 1.5, 31):
        for a2 in np.linspace(-0.8, 0.8, 21):
            for a3 in np.linspace(-0.4, 0.4, 11):
                u1 = transition_unitary(a_mat_ckm, float(a1), float(a2), float(a3))
                u2 = transition_unitary(a_mat_pmns, float(a1), float(a2), float(a3))
                e1 = matrix_error(np.abs(u1), CKM_REF)
                e2 = matrix_error(np.abs(u2), PMNS_REF)
                reg = 0.10 * (abs(a1) + abs(a2) + abs(a3))
                score = e1["mean_rel_pct"] + e2["mean_rel_pct"] + 0.20 * (e1["max_rel_pct"] + e2["max_rel_pct"]) + reg
                if shared_best is None or score < shared_best["score"]:
                    shared_best = {
                        "score": float(score),
                        "a1": float(a1),
                        "a2": float(a2),
                        "a3": float(a3),
                        "ckm_error": e1,
                        "pmns_error": e2,
                        "ckm_pred_abs": np.abs(u1),
                        "pmns_pred_abs": np.abs(u2),
                    }

    threshold = 12.0
    closure_shared = (
        shared_best["ckm_error"]["mean_rel_pct"] < threshold
        and shared_best["pmns_error"]["mean_rel_pct"] < threshold
    )
    verdict = "FLAVOR_CLOSED_WITH_STRUCTURED_OPERATOR" if closure_shared else "FLAVOR_STILL_OPEN_AFTER_STRUCTURED_OPERATOR"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "independent_best": {
            "ckm": {
                "a1": ckm_best["a1"],
                "a2": ckm_best["a2"],
                "a3": ckm_best["a3"],
                "error": ckm_best["error"],
            },
            "pmns": {
                "a1": pmns_best["a1"],
                "a2": pmns_best["a2"],
                "a3": pmns_best["a3"],
                "error": pmns_best["error"],
            },
        },
        "shared_best": {
            "a1": shared_best["a1"],
            "a2": shared_best["a2"],
            "a3": shared_best["a3"],
            "ckm_error": shared_best["ckm_error"],
            "pmns_error": shared_best["pmns_error"],
            "ckm_pred_abs": shared_best["ckm_pred_abs"].tolist(),
            "pmns_pred_abs": shared_best["pmns_pred_abs"].tolist(),
        },
        "threshold_mean_rel_pct": threshold,
        "verdict": verdict,
    }

    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1714: FLAVOR TRANSITION OPERATOR",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Werdykt: **{verdict}**",
        "",
        "## 1) Shared operator test",
        f"- a1={shared_best['a1']:.4f}, a2={shared_best['a2']:.4f}, a3={shared_best['a3']:.4f}",
        f"- CKM mean rel. error: {shared_best['ckm_error']['mean_rel_pct']:.2f}%",
        f"- PMNS mean rel. error: {shared_best['pmns_error']['mean_rel_pct']:.2f}%",
        "",
        "## 2) Interpretation",
        "- Jeśli shared operator nadal nie przechodzi progu, flavor wymaga nowej struktury niż sam wielomian przejść topologicznych 3. rzędu.",
        "",
        "## Artefakty",
        f"- JSON szczegółowy: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"[QW-1714] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1714] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
