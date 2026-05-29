#!/usr/bin/env python3
"""Scratch probe: global admissibility stress test for the local Puiseux inverse.

The formal inverse packet showed that local coefficient matching is easy because
normalized legacy has a nonzero linear term.  This packet asks the next harder
question: can the finite formal inverse be promoted even heuristically to a
real, monotone, nonnegative Z12 distance map?

It cannot.  The cutoff-5 Puiseux branch is excellent near d=0 but turns negative
before d=1 and then blows up on the Z12 range.  This is an honest obstruction to
promoting the local bridge-prep trace into a global bridge without a new
strict-side transport law and admissibility mechanism.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np

HERE = Path(__file__).resolve().parent
FORMAL_INVERSE = HERE / "bridge_phase_normalized_formal_inverse_report.json"
OUT_JSON = HERE / "bridge_phase_normalized_global_admissibility_report.json"
OUT_MD = HERE / "bridge_phase_normalized_global_admissibility_report.md"

LEGACY = {"omega": math.pi / 4.0, "phi": math.pi / 6.0, "beta_tors": 0.01}
STRICT = {"omega": 0.18575, "phi": 0.16250, "beta": 1.0, "eta": 9.0 / 5.0}
LOCAL_D = np.array([1e-4, 1e-3, 1e-2, 5e-2, 1e-1], dtype=float)
Z12_D = np.arange(0, 12, dtype=float)
GRID = np.linspace(0.0, 11.0, 2201)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def parse_exp(label: str) -> float:
    if "/" in label:
        num, den = [int(part) for part in label.split("/")]
        return num / den
    return float(label)


def legacy_normalized(x: np.ndarray) -> np.ndarray:
    return (
        np.cos(LEGACY["omega"] * x + LEGACY["phi"])
        / math.cos(LEGACY["phi"])
        / (1.0 + LEGACY["beta_tors"] * x)
    )


def strict_normalized(d: np.ndarray) -> np.ndarray:
    return (
        np.cos(STRICT["omega"] * d + STRICT["phi"])
        / math.cos(STRICT["phi"])
        / (1.0 + STRICT["beta"] * d ** STRICT["eta"])
    )


def truncated_x(d: np.ndarray, coeffs: dict[str, float]) -> np.ndarray:
    out = np.zeros_like(d, dtype=float)
    for label, coeff in coeffs.items():
        out += float(coeff) * d ** parse_exp(label)
    return out


def residual_rows(d_values: np.ndarray, coeffs: dict[str, float]) -> list[dict[str, float]]:
    x = truncated_x(d_values, coeffs)
    residual = legacy_normalized(x) - strict_normalized(d_values)
    rows = []
    for idx, d in enumerate(d_values):
        rows.append(
            {
                "d": float(d),
                "x_truncated": float(x[idx]),
                "legacy_norm_at_x": float(legacy_normalized(np.array([x[idx]]))[0]),
                "strict_norm": float(strict_normalized(np.array([d]))[0]),
                "residual": float(residual[idx]),
                "abs_residual": float(abs(residual[idx])),
                "nonnegative_x": bool(x[idx] >= 0.0),
            }
        )
    return rows


def admissibility_metrics(coeffs: dict[str, float]) -> dict[str, Any]:
    x_grid = truncated_x(GRID, coeffs)
    dx = np.diff(x_grid)
    first_negative = None
    negative_indices = np.where(x_grid < 0.0)[0]
    if len(negative_indices):
        first_negative = float(GRID[int(negative_indices[0])])
    first_decrease = None
    decrease_indices = np.where(dx < -1e-8)[0]
    if len(decrease_indices):
        first_decrease = float(GRID[int(decrease_indices[0])])
    z_rows = residual_rows(Z12_D, coeffs)
    local_rows = residual_rows(LOCAL_D, coeffs)
    return {
        "grid_size": int(len(GRID)),
        "x_min_on_0_11": float(np.min(x_grid)),
        "x_max_on_0_11": float(np.max(x_grid)),
        "first_grid_d_with_negative_x": first_negative,
        "first_grid_d_with_decreasing_x": first_decrease,
        "is_nonnegative_on_grid": bool(np.all(x_grid >= 0.0)),
        "is_monotone_non_decreasing_on_grid": bool(np.all(dx >= -1e-8)),
        "z12_rows": z_rows,
        "local_rows": local_rows,
        "max_abs_z12_residual": float(max(row["abs_residual"] for row in z_rows)),
        "max_abs_local_residual": float(max(row["abs_residual"] for row in local_rows)),
        "z12_negative_x_count": int(sum(not row["nonnegative_x"] for row in z_rows)),
    }


def main() -> None:
    formal = load_json(FORMAL_INVERSE)
    coeffs = formal["formal_inverse_recurrence"]["x_series_coefficients"]
    metrics = admissibility_metrics(coeffs)
    global_admissible = (
        metrics["is_nonnegative_on_grid"]
        and metrics["is_monotone_non_decreasing_on_grid"]
        and metrics["max_abs_z12_residual"] < 1e-6
    )

    report = {
        "status": "OPEN_GLOBAL_ADMISSIBILITY_STRESS_OBSTRUCTION_NO_BRIDGE_THEOREM",
        "result_kind": "SCRATCH_PHASE_NORMALIZED_GLOBAL_ADMISSIBILITY_PROBE__NOT_A_THEOREM",
        "source_reports": {"formal_inverse": str(FORMAL_INVERSE.relative_to(HERE.parents[1]))},
        "truncated_branch": {
            "max_exp": formal["constants"]["max_exp"],
            "coefficient_count": len(coeffs),
            "coefficients": coeffs,
        },
        "admissibility_metrics": metrics,
        "global_admissible_candidate": global_admissible,
        "restricted_obstruction": {
            "claim": "The cutoff-5 local formal inverse is not a nonnegative monotone global Z12 distance map and does not reproduce the strict normalized kernel on d=0..11.",
            "is_obstructed": not global_admissible,
            "first_negative_d": metrics["first_grid_d_with_negative_x"],
            "first_decreasing_d": metrics["first_grid_d_with_decreasing_x"],
            "max_abs_z12_residual": metrics["max_abs_z12_residual"],
        },
        "honest_interpretation": [
            "The formal inverse remains a valid local asymptotic object near d=0.",
            "The same finite branch fails simple global admissibility checks on the Z12 range.",
            "A future bridge would need a derived global continuation/admissibility rule, not just more local coefficients.",
        ],
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is exported.",
            "No global Z12 distance map is derived.",
            "No strict-side continuation law for the Puiseux branch is exported.",
            "No legacy physical-role transfer is licensed.",
            "No QW-2191 selector discharge and no ToE closure are claimed.",
        ],
        "next_honest_step": "Replace unconstrained local inversion by a strict-side global admissibility problem: prove or falsify a monotone nonnegative continuation law for x(d) on the finite Z12 distances.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch phase-normalized global admissibility probe\n\n"
        "Status: cutoff formal inverse fails global Z12 admissibility; no bridge theorem.\n\n"
        f"- Local rows remain small with max abs residual `{metrics['max_abs_local_residual']:.3e}`.\n"
        f"- On grid [0,11], first negative `x(d)` occurs at `d={metrics['first_grid_d_with_negative_x']}` and first decrease at `d={metrics['first_grid_d_with_decreasing_x']}`.\n"
        f"- Z12 stress: negative x count `{metrics['z12_negative_x_count']}` and max abs residual `{metrics['max_abs_z12_residual']:.3e}`.\n"
        "- Honest read: local formal inversion is not enough; a global strict-side continuation/admissibility law is the real proof object.\n"
        "- No false pass: no kernel identity, no global Z12 map, no strict continuation law, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
