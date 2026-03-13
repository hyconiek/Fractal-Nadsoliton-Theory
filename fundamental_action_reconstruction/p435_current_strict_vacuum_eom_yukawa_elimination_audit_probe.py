#!/usr/bin/env python3
from __future__ import annotations

import json
import random
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"

OUT_JSON = GENERATED / "p435_current_strict_vacuum_eom_yukawa_elimination_audit_probe.json"
OUT_SUMMARY = GENERATED / "p435_current_strict_vacuum_eom_yukawa_elimination_audit_probe_summary.json"


def max_abs(values: list[float]) -> float:
    return float(max(abs(float(v)) for v in values)) if values else 0.0


def gen_nonzero_uniform(rng: random.Random, *, low: float, high: float) -> float:
    x = float(rng.uniform(low, high))
    if rng.random() < 0.5:
        x = -x
    if x == 0.0:
        return float(low)
    return x


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    seed = 435_20260313
    rng = random.Random(seed)

    n = 12
    tol = 1e-12

    vpsi = [gen_nonzero_uniform(rng, low=0.1, high=1.0) for _ in range(n)]
    vphi = gen_nonzero_uniform(rng, low=0.1, high=1.0)

    g4 = [float(rng.uniform(-1.0, 1.0)) for _ in range(n)]
    g6 = [float(rng.uniform(-0.5, 0.5)) for _ in range(n)]
    gY = [float(rng.uniform(-1.0, 1.0)) for _ in range(n)]

    # Ksym(i,j) represents the symmetrized mixing coefficient (K_{i,j}+K_{j,i})/2.
    Ksym: list[list[float]] = [[0.0 for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            val = float(rng.uniform(-0.5, 0.5))
            Ksym[i][j] = val
            Ksym[j][i] = val

    def mixing_row_dot(i: int) -> float:
        return float(sum(Ksym[i][j] * vpsi[j] for j in range(n) if j != i))

    # Construct m2_i so that the constant-vacuum EoM residual is exactly 0 (up to float).
    m2: list[float] = []
    eom_residual: list[float] = []
    for i in range(n):
        mix = mixing_row_dot(i)
        m2_i = float(
            -(mix / vpsi[i])
            - g4[i] * (vpsi[i] ** 2)
            - g6[i] * (vpsi[i] ** 4)
            - 2.0 * gY[i] * (vphi**2)
        )
        m2.append(m2_i)
        res = float(mix + g4[i] * (vpsi[i] ** 3) + g6[i] * (vpsi[i] ** 5) + 2.0 * gY[i] * (vphi**2) * vpsi[i] + m2_i * vpsi[i])
        eom_residual.append(res)

    # Compare diagonal Hessian entry (original) vs Yukawa-free rewrite (N474).
    h_diag: list[float] = []
    h_diag_yukawa_free: list[float] = []
    elim_error: list[float] = []
    for i in range(n):
        h_i = float(3.0 * g4[i] * (vpsi[i] ** 2) + 5.0 * g6[i] * (vpsi[i] ** 4) + 2.0 * gY[i] * (vphi**2) + m2[i])
        mix = mixing_row_dot(i)
        h_i_free = float(-(mix / vpsi[i]) + 2.0 * g4[i] * (vpsi[i] ** 2) + 4.0 * g6[i] * (vpsi[i] ** 4))
        h_diag.append(h_i)
        h_diag_yukawa_free.append(h_i_free)
        elim_error.append(float(h_i - h_i_free))

    artifact: dict[str, Any] = {
        "probe_id": "P435_CURRENT_STRICT_VACUUM_EOM_YUKAWA_ELIMINATION_AUDIT_PROBE",
        "as_of": "2026-03-13",
        "no_false_pass": True,
        "seed": seed,
        "n_sites": n,
        "tolerance": tol,
        "toy_instantiation": {
            "vpsi": vpsi,
            "vphi": vphi,
            "g4": g4,
            "g6": g6,
            "gY": gY,
            "Ksym": Ksym,
            "m2_solved_from_constant_vacuum_eom": m2,
        },
        "checks": {
            "max_abs_constant_vacuum_eom_residual": max_abs(eom_residual),
            "max_abs_yukawa_elimination_error": max_abs(elim_error),
            "passes_at_tolerance": bool(max_abs(eom_residual) <= tol and max_abs(elim_error) <= tol),
        },
        "theorem_pointer": "N474",
        "notes": [
            "Toy audit only: this does not claim any physical or strict-derived instantiation.",
            "m2_i is chosen to satisfy the constant-vacuum EoM residual = 0 for each i.",
            "The check verifies that the diagonal Hessian entry matches the Yukawa-free rewrite from N474.",
        ],
    }

    summary = {
        "probe_id": artifact["probe_id"],
        "as_of": artifact["as_of"],
        "status": "EXECUTED_TOY_AUDIT_NO_FALSE_PASS",
        "n_sites": n,
        "tolerance": tol,
        "max_abs_constant_vacuum_eom_residual": artifact["checks"]["max_abs_constant_vacuum_eom_residual"],
        "max_abs_yukawa_elimination_error": artifact["checks"]["max_abs_yukawa_elimination_error"],
        "passes_at_tolerance": artifact["checks"]["passes_at_tolerance"],
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()

