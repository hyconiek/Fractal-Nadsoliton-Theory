#!/usr/bin/env python3
"""
QW-1717: Parameter identifiability audit for FIN effective closure models.

Tests:
1) Mass correction model identifiability (l1,l2,l3) from QW-1711.
2) Flavor transition model identifiability (a1,a2,a3) from QW-1714.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.linalg import expm


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1717_identifiability_test.json"
OUT_MD = ROOT / "RAPORT_QW1717_IDENTIFIABILITY_TEST.md"


CKM_REF = np.array(
    [
        [0.97401, 0.22650, 0.00361],
        [0.22636, 0.97320, 0.04053],
        [0.00854, 0.03978, 0.99917],
    ],
    dtype=float,
)


def matrix_error(pred_abs: np.ndarray, ref: np.ndarray) -> float:
    rel = np.abs(pred_abs - ref) / np.clip(ref, 1e-12, None)
    return float(np.mean(rel) * 100.0)


def build_generator(q_left: np.ndarray, q_right: np.ndarray) -> np.ndarray:
    n = len(q_left)
    a = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            dq = abs(q_left[i] - q_right[j])
            sign = 1.0 if i < j else -1.0
            a[i, j] = sign / max(dq, 1.0)
    return 0.5 * (a - a.T)


def main() -> None:
    # --- Mass identifiability ---
    gamma = 1.52
    m_top = 173_000.0
    particles = {
        "Top": {"Q": 0, "mass_mev": 173_000.0, "sector": 0.0, "gen": 3.0},
        "Bottom": {"Q": 7, "mass_mev": 4_180.0, "sector": 1.0, "gen": 3.0},
        "Tau": {"Q": 9, "mass_mev": 1_776.9, "sector": -1.0, "gen": 3.0},
        "Charm": {"Q": 9, "mass_mev": 1_270.0, "sector": 1.0, "gen": 2.0},
        "Muon": {"Q": 14, "mass_mev": 105.7, "sector": -1.0, "gen": 2.0},
        "Electron": {"Q": 24, "mass_mev": 0.511, "sector": -1.0, "gen": 1.0},
    }
    names = list(particles.keys())
    q = np.array([particles[n]["Q"] / 24.0 for n in names], dtype=float)
    s = np.array([particles[n]["sector"] for n in names], dtype=float)
    g = np.array([particles[n]["gen"] - 2.0 for n in names], dtype=float)
    y = []
    for n in names:
        p = particles[n]
        base = m_top * (4.0 ** (-(gamma * p["Q"] / 4.0)))
        y.append(math.log(p["mass_mev"] / base))
    y = np.array(y, dtype=float)
    x = np.column_stack([q, s, g])

    # OLS solution + condition number
    beta_hat, *_ = np.linalg.lstsq(x, y, rcond=None)
    xtx = x.T @ x
    cond_mass = float(np.linalg.cond(xtx))

    # Near-optimal parameter volume (score <= min + tol)
    def mass_score(b):
        y_pred = x @ b
        return float(np.mean((y_pred - y) ** 2))

    min_score = mass_score(beta_hat)
    tol = min_score + 1e-4
    candidates_mass = []
    for l1 in np.linspace(beta_hat[0] - 0.25, beta_hat[0] + 0.25, 41):
        for l2 in np.linspace(beta_hat[1] - 0.25, beta_hat[1] + 0.25, 41):
            for l3 in np.linspace(beta_hat[2] - 0.25, beta_hat[2] + 0.25, 41):
                b = np.array([l1, l2, l3], dtype=float)
                sc = mass_score(b)
                if sc <= tol:
                    candidates_mass.append((l1, l2, l3, sc))
    mass_non_identifiable = len(candidates_mass) > 500

    # --- Flavor identifiability (CKM side only, same structure as QW-1714) ---
    q_up = np.array([0.0, 9.0, 14.0], dtype=float)
    q_down = np.array([7.0, 9.0, 14.0], dtype=float)
    a_mat = build_generator(q_up, q_down)

    def ckm_score(a1, a2, a3):
        gmat = a1 * a_mat + a2 * (a_mat @ a_mat) + a3 * (a_mat @ a_mat @ a_mat)
        u = expm(1j * gmat)
        return matrix_error(np.abs(u), CKM_REF)

    # coarse best
    best = None
    for a1 in np.linspace(-1.5, 1.5, 31):
        for a2 in np.linspace(-0.8, 0.8, 21):
            for a3 in np.linspace(-0.4, 0.4, 11):
                sc = ckm_score(a1, a2, a3)
                if best is None or sc < best[0]:
                    best = (sc, a1, a2, a3)
    best_sc, b1, b2, b3 = best

    # near-optimal volume
    tol_ckm = best_sc + 1.0  # 1 percentage point envelope
    candidates_ckm = []
    for a1 in np.linspace(b1 - 0.4, b1 + 0.4, 31):
        for a2 in np.linspace(b2 - 0.3, b2 + 0.3, 21):
            for a3 in np.linspace(b3 - 0.2, b3 + 0.2, 11):
                sc = ckm_score(a1, a2, a3)
                if sc <= tol_ckm:
                    candidates_ckm.append((a1, a2, a3, sc))
    flavor_non_identifiable = len(candidates_ckm) > 800

    # overall verdict
    if mass_non_identifiable and flavor_non_identifiable:
        verdict = "NON_IDENTIFIABLE_BOTH_SECTORS"
    elif mass_non_identifiable or flavor_non_identifiable:
        verdict = "PARTIALLY_NON_IDENTIFIABLE"
    else:
        verdict = "IDENTIFIABLE_WITHIN_TEST_GRID"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "mass_identifiability": {
            "beta_hat": {"l1": float(beta_hat[0]), "l2": float(beta_hat[1]), "l3": float(beta_hat[2])},
            "min_score_mse": float(min_score),
            "condition_number_xtx": cond_mass,
            "n_near_optimal_candidates": len(candidates_mass),
            "non_identifiable_flag": mass_non_identifiable,
        },
        "flavor_identifiability": {
            "best_ckm_mean_rel_pct": float(best_sc),
            "best_params": {"a1": float(b1), "a2": float(b2), "a3": float(b3)},
            "n_near_optimal_candidates": len(candidates_ckm),
            "non_identifiable_flag": flavor_non_identifiable,
        },
        "verdict": verdict,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1717: IDENTIFIABILITY TEST",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Werdykt: **{verdict}**",
        "",
        "## 1) Mass model",
        f"- cond(X^T X) = {cond_mass:.3e}",
        f"- near-optimal candidates = {len(candidates_mass)}",
        f"- non-identifiable flag = {mass_non_identifiable}",
        "",
        "## 2) Flavor model (CKM)",
        f"- best mean rel error = {best_sc:.2f}%",
        f"- near-optimal candidates = {len(candidates_ckm)}",
        f"- non-identifiable flag = {flavor_non_identifiable}",
        "",
        "## Interpretacja",
        "- Duża liczba prawie-równoważnych parametrów oznacza słabą identyfikowalność i ryzyko nadinterpretacji dopasowań.",
        "",
        "## Artefakty",
        f"- JSON szczegółowy: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1717] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1717] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
