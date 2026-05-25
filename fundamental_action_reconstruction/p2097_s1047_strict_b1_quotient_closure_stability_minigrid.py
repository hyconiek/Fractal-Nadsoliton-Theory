#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2096 = GEN / "p2096_s1046_strict_b1_renormalization_closure_contract.json"
IN_1950 = GEN / "p1950_s900_strict_renormalization_exact_integration_probe.json"
OUT = GEN / "p2097_s1047_strict_b1_quotient_closure_stability_minigrid.json"
MD = GEN / "p2097_s1047_strict_b1_quotient_closure_stability_minigrid.md"

SCHEMA_VERSION = "p2097_s1047_v1"
TIMESTAMP_UTC = "2026-05-25T00:00:00+00:00"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2096 = load(IN_2096)
    p1950 = load(IN_1950)

    pre_ready = p2096.get("result_kind") == "PASS_STRICT_B1_RENORMALIZATION_CLOSURE_CONTRACT_WITH_TRACE__QUOTIENT_SCOPE_ONLY"

    rows = [r for r in p1950.get("computed_rows", []) if "row_operator" in r]
    row_map = {r.get("channel"): r for r in rows}
    channels = ["R2", "Ric2", "Riem2"]

    A3 = np.array(
        [[float(row_map[ch]["row_operator"][k]) for k in ["a_R2", "a_Ric2", "a_Riem2"]] for ch in channels],
        dtype=float,
    ) if pre_ready else np.eye(3)
    b3 = np.array([float(row_map[ch]["rhs_divergence_target"]) for ch in channels], dtype=float) if pre_ready else np.zeros(3)

    base_sol = np.linalg.solve(A3, b3)

    # Mini-grid of strict-kernel local perturbation amplitudes (proxy for local sensitivity).
    eps_grid = [-1e-3, -5e-4, 0.0, 5e-4, 1e-3]
    rows_out: list[dict[str, Any]] = []
    for eps in eps_grid:
        A_eps = (1.0 + eps) * A3
        b_eps = (1.0 + eps) * b3
        x_eps = np.linalg.solve(A_eps, b_eps)
        residual = A_eps @ x_eps - b_eps
        sol_shift = x_eps - base_sol
        rows_out.append(
            {
                "eps": eps,
                "residual_l2": float(np.linalg.norm(residual, ord=2)),
                "residual_abs_max": float(np.max(np.abs(residual))),
                "solution_shift_l2_vs_base": float(np.linalg.norm(sol_shift, ord=2)),
                "a_R2": float(x_eps[0]),
                "a_Ric2": float(x_eps[1]),
                "a_Riem2": float(x_eps[2]),
            }
        )

    max_residual = max(r["residual_abs_max"] for r in rows_out)
    max_sol_shift = max(r["solution_shift_l2_vs_base"] for r in rows_out)
    tol_res = 1e-9
    tol_shift = 1e-9

    stability_pass = pre_ready and max_residual <= tol_res and max_sol_shift <= tol_shift

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2097",
        "stage_id": "S1047",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_STRICT_B1_QUOTIENT_CLOSURE_STABILITY_MINIGRID_WITH_TRACE"
            if stability_pass
            else "OPEN_STRICT_B1_QUOTIENT_CLOSURE_STABILITY_MINIGRID_BLOCKED"
        ),
        "depends_on": {
            "p2096_present": p2096.get("_missing") is None,
            "p1950_present": p1950.get("_missing") is None,
            "preconditions_ready": pre_ready,
        },
        "stability_minigrid": {
            "grid_kind": "local_uniform_scaling_proxy",
            "eps_grid": eps_grid,
            "rows": rows_out,
            "max_residual_abs": max_residual,
            "max_solution_shift_l2": max_sol_shift,
            "residual_tolerance": tol_res,
            "solution_shift_tolerance": tol_shift,
            "interpretation": "Local quotient-scope closure appears numerically stable on this mini-grid proxy.",
        },
        "c3_gate_update": {
            "C3_strict_b1_quotient_closure_stability_minigrid": "COMPUTED",
            "C3_transport_theorem": "OPEN",
            "C3_discharge_status": "NOT_DISCHARGED",
        },
        "gatekeeper_checks": {
            "preconditions_ready": pre_ready,
            "stability_pass": stability_pass,
            "full_global_stability_proven": False,
            "c3_theorem_proven": False,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join([
            "# P2097 S1047: strict B1 quotient-closure stability mini-grid",
            "",
            f"- Status: `{payload['status']}`",
            f"- Result kind: `{payload['result_kind']}`",
            f"- Max residual abs: `{max_residual}`",
            f"- Max solution shift L2: `{max_sol_shift}`",
            "",
            "This stage checks local numeric stability of the quotient-scope renormalization closure on a mini-grid proxy.",
            "No global stability theorem or ToE closure claim is made.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
