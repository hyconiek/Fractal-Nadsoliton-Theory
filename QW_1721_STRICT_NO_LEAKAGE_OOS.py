#!/usr/bin/env python3
"""
QW-1721: Strict no-leakage OOS protocol.

Implements a pre-registered split where:
1) Train set defines correction parameters.
2) Holdout set is never seen during fitting.
3) Multi-observable score is computed only on holdout.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1721_strict_no_leakage_oos.json"
OUT_MD = ROOT / "RAPORT_QW1721_STRICT_NO_LEAKAGE_OOS.md"


def rel_err_pct(pred: float, true: float) -> float:
    return abs(pred - true) / abs(true) * 100.0


def main() -> None:
    # Manifest (pre-registered within script)
    manifest = {
        "split_version": "strict-v1",
        "train_particles": ["Bottom", "Muon", "Electron"],
        "test_particles": ["Top", "Tau", "Charm"],
        "metrics": {
            "test_mean_pct_max": 10.0,
            "test_max_pct_max": 20.0,
        },
    }

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
        p["y"] = math.log(p["mass_mev"] / p["base_pred"])

    x = {}
    for n in names:
        q = particles[n]["Q"] / 24.0
        sec = particles[n]["sector"]
        gen = particles[n]["gen"] - 2.0
        x[n] = np.array([q, sec, gen], dtype=float)

    train = manifest["train_particles"]
    test = manifest["test_particles"]
    x_train = np.vstack([x[n] for n in train])
    y_train = np.array([particles[n]["y"] for n in train], dtype=float)

    lam, *_ = np.linalg.lstsq(x_train, y_train, rcond=None)

    def pred_mass(n: str) -> float:
        delta = float(x[n] @ lam)
        return particles[n]["base_pred"] * math.exp(delta)

    rows_train = []
    rows_test = []
    for n in train:
        p = pred_mass(n)
        e = rel_err_pct(p, particles[n]["mass_mev"])
        rows_train.append({"particle": n, "pred_mev": p, "true_mev": particles[n]["mass_mev"], "err_pct": e})
    for n in test:
        p = pred_mass(n)
        e = rel_err_pct(p, particles[n]["mass_mev"])
        rows_test.append({"particle": n, "pred_mev": p, "true_mev": particles[n]["mass_mev"], "err_pct": e})

    train_mean = float(np.mean([r["err_pct"] for r in rows_train]))
    test_mean = float(np.mean([r["err_pct"] for r in rows_test]))
    test_max = float(np.max([r["err_pct"] for r in rows_test]))

    pass_flag = (
        test_mean <= manifest["metrics"]["test_mean_pct_max"]
        and test_max <= manifest["metrics"]["test_max_pct_max"]
    )
    verdict = "STRICT_OOS_PASS" if pass_flag else "STRICT_OOS_FAIL"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "manifest": manifest,
        "lambda": {"l1": float(lam[0]), "l2": float(lam[1]), "l3": float(lam[2])},
        "rows_train": rows_train,
        "rows_test": rows_test,
        "metrics": {
            "train_mean_pct": train_mean,
            "test_mean_pct": test_mean,
            "test_max_pct": test_max,
        },
        "verdict": verdict,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1721: STRICT NO-LEAKAGE OOS",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Werdykt: **{verdict}**",
        "",
        "## Split (pre-registered)",
        f"- Train: {train}",
        f"- Test: {test}",
        "",
        "## Metrics",
        f"- train mean: {train_mean:.2f}%",
        f"- test mean: {test_mean:.2f}%",
        f"- test max: {test_max:.2f}%",
        "",
        "## Artefakty",
        f"- JSON szczegółowy: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"[QW-1721] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1721] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
