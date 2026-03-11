#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_SIGMA_INT = GENERATED / "sigma_int_strict_derived_v1.json"
IN_EPS = GENERATED / "eps_sigma_int_E_pair_amplitude_strict_provenance_v1.json"
IN_SIGN_MASK = GENERATED / "b_sigma_int_E_pair_sign_mask_strict_provenance_v1.json"
IN_QW2190 = REPO / "report_qw2190_kernel_mode_representation_emergence_gate.json"

OUT_JSON = (
    GENERATED / "p403_current_strict_sigma_int_theta_pair_delta_d_sensitivity_audit_probe.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p403_current_strict_sigma_int_theta_pair_delta_d_sensitivity_audit_probe_summary.json"
)


def load_json_path(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


@dataclass(frozen=True)
class KernelTuple:
    omega: float
    phi: float
    beta: float
    eta: float


def corridor_local_window(k: KernelTuple) -> dict[str, float]:
    delta_barrier = (math.pi / 2.0) - abs(k.phi)
    eps_local = 0.5 * delta_barrier
    d_local = eps_local / k.omega
    delta_max = d_local / 11.0
    return {
        "delta_barrier": float(delta_barrier),
        "eps_local": float(eps_local),
        "d_local": float(d_local),
        "delta_max": float(delta_max),
    }


def phasor_pair_theta(
    *,
    k: KernelTuple,
    sigma_int: float,
    eps: float,
    b_pair: list[int],
    delta_d: float,
) -> dict[str, float]:
    x = 0.0
    y = 0.0
    for idx_k in range(12):
        d = float(idx_k) * float(delta_d)
        w = (1.0 + float(sigma_int) * float(eps) * float(b_pair[idx_k])) / 12.0
        theta_d = k.omega * d + k.phi
        denom = 1.0 + k.beta * (d**k.eta)
        x += w * (math.cos(theta_d) / denom)
        y += w * (math.sin(theta_d) / denom)
    theta = math.atan2(y, x)
    r = math.sqrt(x * x + y * y)
    return {"X": float(x), "Y": float(y), "R": float(r), "theta": float(theta)}


def summarize_theta(rows: list[dict[str, Any]]) -> dict[str, Any]:
    t1 = np.array([float(r["pair1"]["theta"]) for r in rows], dtype=float)
    t2 = np.array([float(r["pair2"]["theta"]) for r in rows], dtype=float)
    t1u = np.unwrap(t1)
    t2u = np.unwrap(t2)
    return {
        "pair1": {
            "theta_min": float(np.min(t1)),
            "theta_max": float(np.max(t1)),
            "theta_range": float(np.max(t1) - np.min(t1)),
            "theta_std": float(np.std(t1)),
            "theta_unwrapped_min": float(np.min(t1u)),
            "theta_unwrapped_max": float(np.max(t1u)),
            "theta_unwrapped_range": float(np.max(t1u) - np.min(t1u)),
            "theta_unwrapped_std": float(np.std(t1u)),
        },
        "pair2": {
            "theta_min": float(np.min(t2)),
            "theta_max": float(np.max(t2)),
            "theta_range": float(np.max(t2) - np.min(t2)),
            "theta_std": float(np.std(t2)),
            "theta_unwrapped_min": float(np.min(t2u)),
            "theta_unwrapped_max": float(np.max(t2u)),
            "theta_unwrapped_range": float(np.max(t2u) - np.min(t2u)),
            "theta_unwrapped_std": float(np.std(t2u)),
        },
    }


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    sigma_int_obj = load_json_path(IN_SIGMA_INT)
    eps_obj = load_json_path(IN_EPS)
    sign_mask_obj = load_json_path(IN_SIGN_MASK)
    q2190 = load_json_path(IN_QW2190)

    sigma_int = float(sigma_int_obj["value"])
    eps = float(eps_obj["value"])
    b1 = [int(x) for x in sign_mask_obj["value"]["pair1"]]
    b2 = [int(x) for x in sign_mask_obj["value"]["pair2"]]

    kpars = q2190["kernel"]
    k = KernelTuple(
        omega=float(kpars["omega"]),
        phi=float(kpars["phi"]),
        beta=float(kpars["beta"]),
        eta=float(kpars["eta"]),
    )

    corridor = corridor_local_window(k)
    delta_max = float(corridor["delta_max"])

    # Grid: strictly inside the corridor, plus the historically used 0.25 instance.
    grid = list(np.linspace(delta_max / 12.0, delta_max, 12, dtype=float))
    grid.append(0.25)
    delta_grid = sorted({float(x) for x in grid if 0.0 < float(x) <= delta_max})

    rows: list[dict[str, Any]] = []
    for delta_d in delta_grid:
        p1 = phasor_pair_theta(
            k=k,
            sigma_int=sigma_int,
            eps=eps,
            b_pair=b1,
            delta_d=delta_d,
        )
        p2 = phasor_pair_theta(
            k=k,
            sigma_int=sigma_int,
            eps=eps,
            b_pair=b2,
            delta_d=delta_d,
        )
        rows.append(
            {
                "delta_d": float(delta_d),
                "pair1": p1,
                "pair2": p2,
                "theta_diff_pair1_minus_pair2": float(p1["theta"] - p2["theta"]),
            }
        )

    stats = summarize_theta(rows)
    tol = 1e-12
    depends = bool(
        (float(stats["pair1"]["theta_unwrapped_range"]) > tol)
        or (float(stats["pair2"]["theta_unwrapped_range"]) > tol)
    )

    report = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "stage": "P403",
        "goal": "audit_delta_d_sensitivity_of_sigma_int_strict_theta_pair_inside_positive_window_corridor",
        "inputs": {
            "sigma_int_strict_derived_v1": str(IN_SIGMA_INT.relative_to(REPO)),
            "eps_sigma_int_E_pair_amplitude_strict_provenance_v1": str(IN_EPS.relative_to(REPO)),
            "b_sigma_int_E_pair_sign_mask_strict_provenance_v1": str(IN_SIGN_MASK.relative_to(REPO)),
            "qw2190_kernel_mode_report": str(IN_QW2190.relative_to(REPO)),
        },
        "kernel": {
            "omega": k.omega,
            "phi": k.phi,
            "beta": k.beta,
            "eta": k.eta,
        },
        "corridor": corridor,
        "delta_d_grid": delta_grid,
        "rows": rows,
        "summary_stats": stats,
        "verdict": {
            "tolerance_for_dependence_test": tol,
            "theta_pair_depends_on_delta_d_choice": depends,
            "statement": (
                (
                    "Across admissible delta_d choices in the positive-window corridor, the computed theta-pair varies "
                    f"(unwrapped theta range exceeds tolerance {tol:g}). Therefore, without an explicit delta_d selector "
                    "premise, the sigma_int→theta construction remains candidate-only and cannot be promoted to strict-core "
                    "theta export or to actual object-support discharge."
                )
                if depends
                else (
                    "On the sampled admissible delta_d grid, the theta-pair variation is within numeric tolerance "
                    f"{tol:g} after unwrapping. This probe therefore does not establish delta_d dependence (nor "
                    "independence) as a strict statement; it only reports the sampled-grid outcome without promotion."
                )
            ),
        },
        "hard_limits": [
            "Does not promote any theta value to strict-core theta export.",
            "Does not discharge N302/N395 (object-support above the export-map object).",
            "Does not claim admissible S_sel_int nor strict-core selector closure.",
            "Does not claim QW-2191 discharge.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": report["stage"],
        "theta_pair_depends_on_delta_d_choice": bool(depends),
        "tolerance_for_dependence_test": float(tol),
        "corridor_delta_max": float(delta_max),
        "pair1_theta_unwrapped_range": float(stats["pair1"]["theta_unwrapped_range"]),
        "pair2_theta_unwrapped_range": float(stats["pair2"]["theta_unwrapped_range"]),
        "strict_core_promotion": False,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(json.dumps(summary, ensure_ascii=True))


if __name__ == "__main__":
    main()
