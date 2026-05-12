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


def _moment(order: int, n: int, omega: float, phi: float, beta: float, eta: float) -> float:
    vals = []
    for i in range(1, n + 1):
        d = i / n
        vals.append((d ** order) * k_strict(d, omega, phi, beta, eta))
    return sum(vals) / n


def main() -> None:
    ap = argparse.ArgumentParser(description="P1358: provenance-locked candidate value generator from strict kernel only")
    ap.add_argument("--omega", type=float, default=0.18575)
    ap.add_argument("--phi", type=float, default=0.16250)
    ap.add_argument("--beta", type=float, default=1.0)
    ap.add_argument("--eta", type=float, default=1.8)
    ap.add_argument("--n-grid", type=int, default=4096)
    ap.add_argument("--out-csv", type=Path, default=GEN / "p1358_kernel_predicted_values.csv")
    ap.add_argument("--out-json", type=Path, default=GEN / "p1358_kernel_value_generator_summary.json")
    a = ap.parse_args()

    if a.n_grid <= 16:
        raise SystemExit("n-grid too small for stable candidate moments")

    m0 = _moment(0, a.n_grid, a.omega, a.phi, a.beta, a.eta)
    m1 = _moment(1, a.n_grid, a.omega, a.phi, a.beta, a.eta)
    m2 = _moment(2, a.n_grid, a.omega, a.phi, a.beta, a.eta)
    m3 = _moment(3, a.n_grid, a.omega, a.phi, a.beta, a.eta)

    # Candidate, dimensionless, strict-kernel-only map (no observational injection).
    pred = {
        "g1": abs(m0),
        "g2": abs(m1) * 3.0,
        "g3": abs(m2) * 12.0,
        "gravity_effective_observable_1": abs(m3) * 40.0,
    }

    with a.out_csv.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["observable", "predicted"])
        for k, v in pred.items():
            w.writerow([k, f"{v:.12f}"])

    out = {
        "packet": "P1358",
        "as_of": "2026-05-12",
        "strict_kernel": {
            "omega": a.omega,
            "phi": a.phi,
            "beta": a.beta,
            "eta": a.eta,
            "n_grid": a.n_grid,
        },
        "moments": {"m0": m0, "m1": m1, "m2": m2, "m3": m3},
        "predicted_values_dimensionless": pred,
        "observational_injection": False,
        "status": "CANDIDATE_NONCALIBRATED",
        "next_priority": "P1358B_JOIN_WITH_REFERENCE_DATA_AND_RESIDUAL_EVALUATION",
    }
    a.out_json.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1358] wrote {a.out_csv} and {a.out_json}")


if __name__ == "__main__":
    main()
