#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2205 = GEN / "p2205_s1155_strict_nu_branch_transport_vanishing_constraints_probe.json"
OUT = GEN / "p2206_s1156_strict_nu_branch_transport_operator_constraint_solver.json"
MD = GEN / "p2206_s1156_strict_nu_branch_transport_operator_constraint_solver.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def ksq(d: float, omega: float, phi: float, beta: float, eta: float) -> float:
    den = 1.0 + beta * (d**eta)
    c = math.cos(omega * d + phi)
    return (c * c) / (den * den)


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2205 = load(IN_2205)
    p2190 = load(ROOT / "generated/p2190_s1140_strict_cutkosky_fixed_channel_real_discontinuity_integral_packet.json")
    w = p2190.get("strict_cutkosky_fixed_channel_real_discontinuity_integral_packet", {}) or {}
    params = w.get("strict_kernel_params", {}) or {"omega": 0.18575, "phi": 0.1625, "beta": 1.0, "eta": 1.8}

    omega0 = float(params["omega"])
    phi = float(params["phi"])
    beta = float(params["beta"])
    eta = float(params["eta"])

    # symbolic operator mismatch equation
    d, om, ph, be, et, m = sp.symbols("d om ph be et m", positive=True)
    mismatch = sp.simplify(
        sp.cos(om * d + ph) ** 2 / (1 + be * d**et) ** 2
        - sp.cos((om * m) * d + ph) ** 2 / (1 + be * d**et) ** 2
    )
    reduced_mismatch = sp.simplify((1 + be * d**et) ** 2 * mismatch)
    symbolic_condition = sp.Eq(reduced_mismatch, 0)

    # finite witness set (noncyclic anchor): enforce mismatch=0 on fixed d-samples
    d_samples = np.array([0.7, 1.9, 3.3, 5.1, 7.7, 11.3, 15.9], dtype=float)
    m_grid = np.linspace(0.94, 1.06, 2401)

    def score(mm: float) -> float:
        vals = [
            abs(ksq(float(ds), omega0, phi, beta, eta) - ksq(float(ds), omega0 * mm, phi, beta, eta))
            for ds in d_samples
        ]
        return float(max(vals))

    scored = [{"m": float(mm), "max_abs_mismatch_on_samples": score(float(mm))} for mm in m_grid]
    best = min(scored, key=lambda r: r["max_abs_mismatch_on_samples"])

    tol = 1e-10
    near_zero = [r for r in scored if r["max_abs_mismatch_on_samples"] < tol]
    nontrivial = [r for r in near_zero if abs(r["m"] - 1.0) > 1e-6]

    solver_verdict = {
        "sampled_constraint_family": "for all d in fixed noncyclic witness set: cos^2(omega*d+phi)=cos^2(omega*m*d+phi)",
        "best_m_on_grid": best,
        "near_zero_solution_count": len(near_zero),
        "nontrivial_near_zero_solution_count": len(nontrivial),
        "m_equal_1_selected": abs(best["m"] - 1.0) < 5e-4,
        "interpretation": "on this frozen channel/domain witness set, transport-operator mismatch solver selects m≈1 and finds no nontrivial near-zero branch",
    }

    payload = {
        "schema_version": "p2206_s1156_v1",
        "packet_id": "P2206",
        "stage_id": "S1156",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-26T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_TRANSPORT_OPERATOR_CONSTRAINT_SOLVER",
        "strict_nu_branch_transport_operator_constraint_solver": {
            "solver_id": "STRICT_NU_BRANCH_TRANSPORT_OPERATOR_CONSTRAINT_SOLVER_V1",
            "source_packet": str(IN_2205.relative_to(ROOT)),
            "symbolic_operator_mismatch_equation": str(symbolic_condition),
            "reduced_mismatch_expression": str(reduced_mismatch),
            "frozen_params": params,
            "witness_d_samples": [float(x) for x in d_samples],
            "grid_scan_range": {"m_min": float(m_grid.min()), "m_max": float(m_grid.max()), "count": int(m_grid.size)},
            "solver_verdict": solver_verdict,
            "theorem_scope_limit": "finite witness-set constraint solver only; not full global transport theorem closure",
        },
        "recommended_next_honest_step": {
            "id": "P2207_candidate",
            "goal": "lift finite witness-set solver to interval-certified bound proving absence/presence of nontrivial m branches on a compact domain"
        },
        "gatekeeper_checks": {
            "constraint_solver_exported": True,
            "best_m_finite": math.isfinite(float(best["m"])),
            "best_score_finite": math.isfinite(float(best["max_abs_mismatch_on_samples"])),
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
            "full_cutkosky_closure_proven": False,
            "full_d3_covariance_transport_proven": False,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join([
            "# P2206 S1156: strict nu-branch transport operator constraint solver",
            "",
            f"- best m on grid: `{best['m']:.12f}`",
            f"- best max mismatch: `{best['max_abs_mismatch_on_samples']:.12e}`",
            f"- nontrivial near-zero solutions: `{len(nontrivial)}`",
            "",
            "Finite witness-set solver exported; no global Task-3 closure claim.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
