#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

FAR = Path(__file__).resolve().parent
GEN = FAR / "generated"


def k_strict(d: float, omega: float, phi: float, beta: float, eta: float) -> float:
    return math.cos(omega * d + phi) / (1.0 + beta * (d ** eta))


def moment(order: int, n: int, omega: float, phi: float, beta: float, eta: float) -> float:
    return sum(((i / n) ** order) * k_strict(i / n, omega, phi, beta, eta) for i in range(1, n + 1)) / n


def map_values(m0: float, m1: float, m2: float, m3: float, gain: float, suppress: float) -> dict[str, float]:
    return {
        "g1": abs(m0 * gain),
        "g2": abs(m1 * gain) * 3.0,
        "g3": abs(m2 * suppress) * 12.0,
        "gravity_effective_observable_1": abs(m3 * suppress) * 40.0,
    }


def load_ref(path: Path) -> dict[str, tuple[float, float]]:
    out = {}
    with path.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f)
        for row in r:
            out[row["observable"]] = (float(row["observed"]), float(row["sigma"]))
    return out


def score(pred: dict[str, float], ref: dict[str, tuple[float, float]], gain: float, suppress: float) -> dict:
    zvals = [abs(pred[k] - ref[k][0]) / ref[k][1] for k in pred]
    max_abs_z = max(zvals)
    mean_abs_z = sum(zvals) / len(zvals)
    reg = abs(gain - 1.0) + abs(suppress - 1.0)
    objective = mean_abs_z + 0.1 * reg
    return {"max_abs_z": max_abs_z, "mean_abs_z": mean_abs_z, "regularization": reg, "objective": objective}


def main() -> None:
    ap = argparse.ArgumentParser(description="P1368 sweep gain/suppress and report holdout-readiness")
    ap.add_argument("--ref", type=Path, default=GEN / "p1358_reference_observed_values.csv")
    ap.add_argument("--n-grid", type=int, default=4096)
    ap.add_argument("--omega", type=float, default=0.18575)
    ap.add_argument("--phi", type=float, default=0.1625)
    ap.add_argument("--beta", type=float, default=1.0)
    ap.add_argument("--eta", type=float, default=1.8)
    ap.add_argument("--out", type=Path, default=GEN / "p1368_gain_suppress_sweep_summary.json")
    a = ap.parse_args()

    ref = load_ref(a.ref)
    m0, m1, m2, m3 = [moment(i, a.n_grid, a.omega, a.phi, a.beta, a.eta) for i in range(4)]

    gains = [0.95, 1.00, 1.05, 1.10, 1.15]
    suppresses = [0.75, 0.80, 0.85, 0.90, 0.95]

    trials = []
    for g in gains:
        for s in suppresses:
            pred = map_values(m0, m1, m2, m3, g, s)
            sc = score(pred, ref, g, s)
            trials.append({"gain": g, "suppress": s, **sc})

    best = min(trials, key=lambda x: x["objective"])
    holdout_ready = best["max_abs_z"] < 5.0

    out = {
        "packet": "P1368",
        "as_of": "2026-05-12",
        "trial_count": len(trials),
        "best": best,
        "holdout_ready": holdout_ready,
        "next_priority": "P1369_HOLDOUT_PROTOCOL" if holdout_ready else "P1369_MODEL_CLASS_REFINEMENT_BEFORE_HOLDOUT",
        "trials": trials,
    }
    a.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1368] wrote {a.out}; holdout_ready={holdout_ready}; best_max_abs_z={best['max_abs_z']:.3f}")


if __name__ == "__main__":
    main()
