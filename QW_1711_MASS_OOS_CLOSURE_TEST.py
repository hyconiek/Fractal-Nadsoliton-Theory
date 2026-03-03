#!/usr/bin/env python3
"""
QW-1711: Out-of-sample closure test for FIN mass sector.

Objective:
1) Use frozen core formula as baseline.
2) Fit minimal effective correction on train set only.
3) Evaluate generalization on holdout particles.
4) Check if sector can be considered closed without leakage.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1711_mass_oos_closure_test.json"
OUT_MD = ROOT / "RAPORT_QW1711_MASS_OOS_CLOSURE_TEST.md"


def rel_err_pct(pred: float, target: float) -> float:
    return abs(pred - target) / abs(target) * 100.0


def main() -> None:
    # Frozen FIN core
    gamma = 1.52
    m_top = 173_000.0  # MeV

    # Particle table: Q assignments from current FIN narrative
    particles = {
        "Top": {"Q": 0, "mass_mev": 173_000.0, "sector": 0.0, "gen": 3.0},
        "Bottom": {"Q": 7, "mass_mev": 4_180.0, "sector": 1.0, "gen": 3.0},
        "Tau": {"Q": 9, "mass_mev": 1_776.9, "sector": -1.0, "gen": 3.0},
        "Charm": {"Q": 9, "mass_mev": 1_270.0, "sector": 1.0, "gen": 2.0},
        "Muon": {"Q": 14, "mass_mev": 105.7, "sector": -1.0, "gen": 2.0},
        "Electron": {"Q": 24, "mass_mev": 0.511, "sector": -1.0, "gen": 1.0},
    }

    # Baseline prediction
    for p in particles.values():
        q = p["Q"]
        p["base_pred"] = m_top * (4.0 ** (-(gamma * q / 4.0)))
        p["base_log_delta"] = math.log(p["mass_mev"] / p["base_pred"])

    # Feature model for correction term:
    # Delta_a = l1*(Q/24) + l2*sector + l3*(gen-2)
    names = list(particles.keys())
    q_feat = np.array([particles[n]["Q"] / 24.0 for n in names], dtype=float)
    sec_feat = np.array([particles[n]["sector"] for n in names], dtype=float)
    gen_feat = np.array([particles[n]["gen"] - 2.0 for n in names], dtype=float)
    y = np.array([particles[n]["base_log_delta"] for n in names], dtype=float)

    # Fixed split (transparent, deterministic)
    train_names = ["Bottom", "Muon", "Electron"]
    test_names = ["Tau", "Charm", "Top"]
    idx_train = np.array([names.index(n) for n in train_names], dtype=int)
    idx_test = np.array([names.index(n) for n in test_names], dtype=int)

    x_train = np.column_stack([q_feat[idx_train], sec_feat[idx_train], gen_feat[idx_train]])
    y_train = y[idx_train]

    # OLS for lambda parameters (minimal correction model)
    lam, *_ = np.linalg.lstsq(x_train, y_train, rcond=None)

    def corrected_pred(i: int) -> float:
        delta = lam[0] * q_feat[i] + lam[1] * sec_feat[i] + lam[2] * gen_feat[i]
        return particles[names[i]]["base_pred"] * math.exp(delta)

    train_rows = []
    test_rows = []
    for i, n in enumerate(names):
        pred = corrected_pred(i)
        row = {
            "particle": n,
            "exp_mev": particles[n]["mass_mev"],
            "base_pred_mev": particles[n]["base_pred"],
            "corr_pred_mev": pred,
            "base_err_pct": rel_err_pct(particles[n]["base_pred"], particles[n]["mass_mev"]),
            "corr_err_pct": rel_err_pct(pred, particles[n]["mass_mev"]),
        }
        if i in idx_train:
            train_rows.append(row)
        else:
            test_rows.append(row)

    train_mae = float(np.mean([r["corr_err_pct"] for r in train_rows]))
    test_mae = float(np.mean([r["corr_err_pct"] for r in test_rows]))
    generalization_gap = test_mae - train_mae

    # Closure criterion: holdout error and gap must be both small
    closure_test_threshold = 10.0
    gap_threshold = 7.0
    closed_oos = (test_mae < closure_test_threshold) and (generalization_gap < gap_threshold)

    verdict = "MASS_SECTOR_OOS_CLOSED" if closed_oos else "MASS_SECTOR_NOT_OOS_CLOSED"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "frozen_core": {
            "gamma": gamma,
            "m_top_mev": m_top,
        },
        "model": {
            "delta_formula": "Delta = l1*(Q/24) + l2*sector + l3*(gen-2)",
            "lambda": {
                "l1": float(lam[0]),
                "l2": float(lam[1]),
                "l3": float(lam[2]),
            },
        },
        "split": {
            "train": train_names,
            "test": test_names,
        },
        "train_rows": train_rows,
        "test_rows": test_rows,
        "metrics": {
            "train_mean_err_pct": train_mae,
            "test_mean_err_pct": test_mae,
            "generalization_gap_pct": generalization_gap,
            "closure_test_threshold_pct": closure_test_threshold,
            "gap_threshold_pct": gap_threshold,
        },
        "verdict": verdict,
    }

    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1711: MASS OOS CLOSURE TEST",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Werdykt: **{verdict}**",
        "",
        "## 1) Model korekcyjny",
        "- Delta = l1*(Q/24) + l2*sector + l3*(gen-2)",
        f"- l1={lam[0]:+.6f}, l2={lam[1]:+.6f}, l3={lam[2]:+.6f}",
        "",
        "## 2) Wynik out-of-sample",
        f"- Train mean error: {train_mae:.2f}%",
        f"- Test mean error: {test_mae:.2f}%",
        f"- Generalization gap: {generalization_gap:.2f} p.p.",
        "",
        "## 3) Test holdout",
    ]
    for r in test_rows:
        lines.append(
            f"- {r['particle']}: base={r['base_err_pct']:.2f}%, corrected={r['corr_err_pct']:.2f}%"
        )
    lines.extend(
        [
            "",
            "## Interpretacja",
            "- Jeśli test nie przechodzi, obecny sektor mas nie jest domknięty predykcyjnie i wymaga bogatszych niezmienników topologicznych.",
            "",
            "## Artefakty",
            f"- JSON szczegółowy: `{OUT_JSON.name}`",
        ]
    )
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1711] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1711] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
