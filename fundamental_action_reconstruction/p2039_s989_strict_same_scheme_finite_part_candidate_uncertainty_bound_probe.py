#!/usr/bin/env python3
"""P2039 S989: uncertainty-bound probe for P2038 finite-part candidate.

Next honest step after P2038: compute a first uncertainty bound for the
nonzero finite-part correction candidate matrix on a small strict-kernel grid,
for one controlled background pair.

This packet remains probe-level: no theorem-grade finite-part transport claim,
no C3 discharge, no background globalization.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2039_s989_strict_same_scheme_finite_part_candidate_uncertainty_bound_probe.json"
MD = GEN / "p2039_s989_strict_same_scheme_finite_part_candidate_uncertainty_bound_probe.md"

SCHEMA_VERSION = "p2039_s989_v1"
TIMESTAMP_UTC = "2026-05-25T00:00:00+00:00"
BASIS = ["R2_bar", "Ric2_bar", "Riem2_bar"]


def load(name: str) -> dict[str, Any]:
    path = GEN / name
    if not path.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def file_sha256(path: Path) -> str:
    if not path.exists():
        return "MISSING"
    return hashlib.sha256(path.read_bytes()).hexdigest()


def as_bool(x: Any) -> bool:
    return bool(x is True)


def max_abs_matrix(m: list[list[float]]) -> float:
    return max(abs(v) for row in m for v in row)


def mat_scale(m: list[list[float]], s: float) -> list[list[float]]:
    return [[float(v * s) for v in row] for row in m]


def mat_sub(a: list[list[float]], b: list[list[float]]) -> list[list[float]]:
    return [[a[i][j] - b[i][j] for j in range(3)] for i in range(3)]


def linf(m: list[list[float]]) -> float:
    return max(abs(v) for row in m for v in row)


def strict_kernel(d: np.ndarray, omega: float, phi: float, beta: float, eta: float) -> np.ndarray:
    return np.cos(omega * d + phi) / (1.0 + beta * np.power(d, eta))


def strict_kernel_scale(omega: float, phi: float, beta: float, eta: float) -> float:
    # Data-derived same-scheme sensitivity scale from strict kernel profile change
    # on a bounded d-grid, normalized to the baseline profile.
    d = np.linspace(0.05, 6.0, 240)
    k0 = strict_kernel(d, 0.18575, 0.16250, 1.0, 1.8)
    k1 = strict_kernel(d, omega, phi, beta, eta)
    n0 = float(np.trapezoid(np.abs(k0), d))
    diff = float(np.trapezoid(np.abs(k1 - k0), d))
    rel = diff / n0 if n0 > 0.0 else 0.0
    # Keep scaling close to 1 on the local grid; this is still a probe-level map.
    return 1.0 + rel


def main() -> None:
    GEN.mkdir(exist_ok=True)

    p2038 = load("p2038_s988_strict_same_scheme_finite_part_map_data_import_audit.json")
    checks_2038 = p2038.get("gatekeeper_checks") or {}
    import_data = p2038.get("candidate_data_import") or {}

    preconditions_ready = (
        p2038.get("result_kind")
        == "PASS_FIRST_NONZERO_SAME_SCHEME_FINITE_PART_CANDIDATE_IMPORTED_WITH_TRACE__C3_STILL_OPEN"
        and as_bool(checks_2038.get("nonzero_candidate_present"))
    )

    delta_base = import_data.get("delta_finite_part_candidate_matrix") or [[0.0, 0.0, 0.0]] * 3
    delta_base = [[float(x) for x in row] for row in delta_base]

    # Small strict-kernel grid around current working tuple.
    omega_grid = [0.18565, 0.18575, 0.18585]
    phi_grid = [0.16240, 0.16250, 0.16260]
    beta_grid = [0.99, 1.00, 1.01]
    eta_grid = [1.79, 1.80, 1.81]

    points: list[dict[str, Any]] = []
    deviations: list[float] = []

    for omega in omega_grid:
        for phi in phi_grid:
            for beta in beta_grid:
                for eta in eta_grid:
                    s = strict_kernel_scale(omega, phi, beta, eta)
                    delta_here = mat_scale(delta_base, s)
                    d = linf(mat_sub(delta_here, delta_base))
                    deviations.append(d)
                    points.append(
                        {
                            "omega": omega,
                            "phi": phi,
                            "beta": beta,
                            "eta": eta,
                            "scale": s,
                            "linf_deviation_from_base": d,
                        }
                    )

    max_dev = max(deviations) if deviations else math.inf
    min_dev = min(deviations) if deviations else math.inf
    mean_dev = sum(deviations) / len(deviations) if deviations else math.inf
    base_linf = max_abs_matrix(delta_base)
    rel_bound = (max_dev / base_linf) if base_linf > 0 else math.inf

    uncertainty_bound = {
        "bound_type": "small_grid_linf_uncertainty_bound",
        "absolute_linf_bound": max_dev,
        "relative_linf_bound_vs_base": rel_bound,
        "base_candidate_linf_norm": base_linf,
        "min_grid_deviation_linf": min_dev,
        "mean_grid_deviation_linf": mean_dev,
        "grid_size": len(points),
        "stability_pass_under_threshold": bool(max_dev <= 2.0e-5),
        "threshold_linf": 2.0e-5,
    }

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2039",
        "stage_id": "S989",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_FIRST_COMPUTED_UNCERTAINTY_BOUND_FOR_NONZERO_FINITE_PART_CANDIDATE_WITH_TRACE__C3_STILL_OPEN"
            if preconditions_ready and math.isfinite(max_dev)
            else "OPEN_UNCERTAINTY_BOUND_PROBE_BLOCKED"
        ),
        "route": "strict_only_small_grid_stability_probe",
        "depends_on": {
            "p2038_present": p2038.get("_missing") is None,
            "p2038_preconditions_ready": preconditions_ready,
        },
        "input_hashes": {
            "p2038_json_sha256": file_sha256(GEN / "p2038_s988_strict_same_scheme_finite_part_map_data_import_audit.json"),
        },
        "probe_scope": {
            "controlled_background_pair": import_data.get("controlled_background_pair", "UNKNOWN"),
            "basis": BASIS,
            "strict_kernel_grid": {
                "omega": omega_grid,
                "phi": phi_grid,
                "beta": beta_grid,
                "eta": eta_grid,
            },
            "modeling_note": "Grid-scanned strict-kernel profile sensitivity scaling on bounded d-grid for first uncertainty bound.",
        },
        "uncertainty_bound": uncertainty_bound,
        "grid_samples": points,
        "c3_gate_update": {
            "C3_uncertainty_bound_probe": "COMPUTED_FOR_ONE_CONTROLLED_PAIR",
            "C3_transport_theorem": "OPEN",
            "C3_discharge_status": "NOT_DISCHARGED",
            "missing_for_c3_theorem": [
                "proof of same-scheme subtraction compatibility",
                "operator-level transport identity theorem",
                "cross-background finite-part lock theorem",
            ],
        },
        "gatekeeper_checks": {
            "preconditions_ready": preconditions_ready,
            "grid_size_81": len(points) == 81,
            "base_candidate_nonzero": base_linf > 0.0,
            "uncertainty_bound_finite": math.isfinite(max_dev),
            "stability_threshold_checked": True,
            "c3_theorem_proven": False,
            "no_background_globalization_claimed": True,
            "no_tensor_component_claimed": True,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    md = [
        "# P2039 S989: uncertainty bound probe for nonzero finite-part candidate",
        "",
        f"- Status: `{payload['status']}`",
        f"- Result kind: `{payload['result_kind']}`",
        "",
        "## Export",
        "",
        "Computed first small-grid uncertainty bound (linf) for the nonzero finite-part",
        "candidate matrix imported in P2038, on basis `R2_bar/Ric2_bar/Riem2_bar`.",
        "",
        "## Gate update",
        "",
        "- `C3`: uncertainty bound probe computed for one controlled pair.",
        "- `C3`: theorem remains open (not discharged).",
        "",
        "## Discipline",
        "",
        "No theorem-grade transport claim, no background globalization claim, no ToE closure claim.",
    ]
    MD.write_text("\n".join(md) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
