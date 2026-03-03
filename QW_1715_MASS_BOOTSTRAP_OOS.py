#!/usr/bin/env python3
"""
QW-1715: Bootstrap out-of-sample stress test for FIN mass closure.

Runs many train/test splits for the correction model from QW-1711 and reports
distribution of holdout errors and generalization gaps.
"""

from __future__ import annotations

import itertools
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1715_mass_bootstrap_oos.json"
OUT_MD = ROOT / "RAPORT_QW1715_MASS_BOOTSTRAP_OOS.md"


def rel_err_pct(pred: float, target: float) -> float:
    return abs(pred - target) / abs(target) * 100.0


def main() -> None:
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

    for p in particles.values():
        p["base_pred"] = m_top * (4.0 ** (-(gamma * p["Q"] / 4.0)))
        p["base_log_delta"] = math.log(p["mass_mev"] / p["base_pred"])

    q_feat = np.array([particles[n]["Q"] / 24.0 for n in names], dtype=float)
    sec_feat = np.array([particles[n]["sector"] for n in names], dtype=float)
    gen_feat = np.array([particles[n]["gen"] - 2.0 for n in names], dtype=float)
    y = np.array([particles[n]["base_log_delta"] for n in names], dtype=float)

    idx_all = np.arange(len(names))
    split_rows = []

    # All train subsets of size 3 (combinatorial bootstrap-style stress test)
    for train_idx in itertools.combinations(idx_all, 3):
        train_idx = np.array(train_idx, dtype=int)
        test_idx = np.array([i for i in idx_all if i not in set(train_idx)], dtype=int)

        x_train = np.column_stack([q_feat[train_idx], sec_feat[train_idx], gen_feat[train_idx]])
        y_train = y[train_idx]
        lam, *_ = np.linalg.lstsq(x_train, y_train, rcond=None)

        def pred(i: int) -> float:
            delta = lam[0] * q_feat[i] + lam[1] * sec_feat[i] + lam[2] * gen_feat[i]
            return particles[names[i]]["base_pred"] * math.exp(delta)

        train_errs = [rel_err_pct(pred(i), particles[names[i]]["mass_mev"]) for i in train_idx]
        test_errs = [rel_err_pct(pred(i), particles[names[i]]["mass_mev"]) for i in test_idx]
        split_rows.append(
            {
                "train": [names[i] for i in train_idx],
                "test": [names[i] for i in test_idx],
                "lambda": {"l1": float(lam[0]), "l2": float(lam[1]), "l3": float(lam[2])},
                "train_mean_err_pct": float(np.mean(train_errs)),
                "test_mean_err_pct": float(np.mean(test_errs)),
                "gap_pct": float(np.mean(test_errs) - np.mean(train_errs)),
            }
        )

    train_means = np.array([r["train_mean_err_pct"] for r in split_rows], dtype=float)
    test_means = np.array([r["test_mean_err_pct"] for r in split_rows], dtype=float)
    gaps = np.array([r["gap_pct"] for r in split_rows], dtype=float)

    metrics = {
        "n_splits": int(len(split_rows)),
        "train_mean_pct": float(np.mean(train_means)),
        "train_median_pct": float(np.median(train_means)),
        "test_mean_pct": float(np.mean(test_means)),
        "test_median_pct": float(np.median(test_means)),
        "test_p90_pct": float(np.quantile(test_means, 0.90)),
        "gap_mean_pct": float(np.mean(gaps)),
        "gap_median_pct": float(np.median(gaps)),
    }

    threshold_test_median = 10.0
    threshold_gap_mean = 7.0
    robust_oos = (metrics["test_median_pct"] < threshold_test_median) and (metrics["gap_mean_pct"] < threshold_gap_mean)
    verdict = "MASS_OOS_ROBUST" if robust_oos else "MASS_OOS_NOT_ROBUST"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "frozen_core": {"gamma": gamma, "m_top_mev": m_top},
        "metrics": metrics,
        "thresholds": {
            "test_median_pct": threshold_test_median,
            "gap_mean_pct": threshold_gap_mean,
        },
        "verdict": verdict,
        "splits_top_40": split_rows[:40],
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1715: MASS BOOTSTRAP OOS",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Werdykt: **{verdict}**",
        "",
        "## Agregat splitów",
        f"- n_splits: {metrics['n_splits']}",
        f"- train mean/median: {metrics['train_mean_pct']:.2f}% / {metrics['train_median_pct']:.2f}%",
        f"- test mean/median/p90: {metrics['test_mean_pct']:.2f}% / {metrics['test_median_pct']:.2f}% / {metrics['test_p90_pct']:.2f}%",
        f"- gap mean/median: {metrics['gap_mean_pct']:.2f} / {metrics['gap_median_pct']:.2f} p.p.",
        "",
        "## Interpretacja",
        "- Jeśli bootstrap OOS nie jest robust, sektor mas wymaga nowego modelu cech topologicznych i/lub innej struktury korekt.",
        "",
        "## Artefakty",
        f"- JSON szczegółowy: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"[QW-1715] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1715] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
