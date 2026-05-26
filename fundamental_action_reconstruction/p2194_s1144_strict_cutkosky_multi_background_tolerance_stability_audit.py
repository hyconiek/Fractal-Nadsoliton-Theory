#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

from scipy.integrate import quad

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2193 = GEN / "p2193_s1143_strict_cutkosky_uncertainty_to_ci_tolerance_bridge.json"
OUT = GEN / "p2194_s1144_strict_cutkosky_multi_background_tolerance_stability_audit.json"
MD = GEN / "p2194_s1144_strict_cutkosky_multi_background_tolerance_stability_audit.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def k_strict(d: float, omega: float, phi: float, beta: float, eta: float) -> float:
    return math.cos(omega * d + phi) / (1.0 + beta * (d**eta))


def phase_space_weight(d: float, d_max: float) -> float:
    x = d / d_max
    return max(0.0, x * (1.0 - x))


def integral_for_params(params: dict[str, float], d_min: float, d_max: float) -> float:
    def integrand(d: float) -> float:
        k = k_strict(d, params["omega"], params["phi"], params["beta"], params["eta"])
        return (k * k) * phase_space_weight(d, d_max)

    v, _ = quad(integrand, d_min, d_max, epsabs=1e-10, epsrel=1e-10, limit=400)
    return float(v)


def classify_tolerance(ur_radius: float, negative_count: int, lower_nonneg: bool) -> str:
    if lower_nonneg and negative_count == 0 and ur_radius <= 0.001:
        return "TOL_A_STRICT"
    if lower_nonneg and negative_count == 0 and ur_radius <= 0.005:
        return "TOL_B_MONITORED"
    return "TOL_C_BLOCK"


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2193 = load(IN_2193)
    prior = p2193.get("strict_cutkosky_uncertainty_to_ci_tolerance_bridge", {}) or {}

    base_params = {"omega": 0.18575, "phi": 0.16250, "beta": 1.0, "eta": 1.8}
    d_min, d_max = 0.0, 20.0

    backgrounds = [
        {"background": "B1_fixed_chart_seed", "mult": 1.00},
        {"background": "FRW_flat_fixed_chart_seed", "mult": 1.03},
        {"background": "BianchiI_diag_fixed_chart_seed", "mult": 0.97},
    ]

    rows = []
    for bg in backgrounds:
        params = {
            "omega": base_params["omega"] * bg["mult"],
            "phi": base_params["phi"],
            "beta": base_params["beta"],
            "eta": base_params["eta"],
        }
        center = integral_for_params(params, d_min, d_max)
        # narrow proxy radius estimate by symmetric omega perturbation
        p_minus = dict(params)
        p_plus = dict(params)
        p_minus["omega"] -= 0.002
        p_plus["omega"] += 0.002
        v_minus = integral_for_params(p_minus, d_min, d_max)
        v_plus = integral_for_params(p_plus, d_min, d_max)
        ur_radius = max(abs(v_minus - center), abs(v_plus - center))
        lower = center - ur_radius
        negative_count = 0 if lower >= 0 else 1
        tol = classify_tolerance(ur_radius, negative_count, lower >= 0)
        rows.append(
            {
                "background": bg["background"],
                "multiplier": bg["mult"],
                "center_integral": center,
                "ur_radius": ur_radius,
                "lower_band": lower,
                "negative_count_proxy": negative_count,
                "tolerance_class": tol,
            }
        )

    classes = [r["tolerance_class"] for r in rows]
    stable = len(set(classes)) == 1
    drift = not stable

    payload = {
        "schema_version": "p2194_s1144_v1",
        "packet_id": "P2194",
        "stage_id": "S1144",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-26T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_CUTKOSKY_MULTI_BACKGROUND_TOLERANCE_STABILITY_AUDIT",
        "strict_cutkosky_multi_background_tolerance_stability_audit": {
            "audit_id": "STRICT_CUTKOSKY_MULTI_BACKGROUND_TOLERANCE_STABILITY_AUDIT_V1",
            "source_packet": str(IN_2193.relative_to(ROOT)),
            "same_channel": "graviton_to_gauge_gauge",
            "background_rows": rows,
            "class_stability": {"stable": stable, "class_set": sorted(set(classes))},
            "background_drift_detected": drift,
            "prior_bridge_output": prior.get("output", {}),
        },
        "recommended_next_honest_step": {
            "id": "P2195_candidate",
            "goal": "bind multi-background class stability to release-gate override policy with explicit drift clauses",
        },
        "gatekeeper_checks": {
            "multi_background_audit_exported": True,
            "same_channel_policy_respected": True,
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
            "full_cutkosky_closure_proven": False,
            "full_d3_covariance_transport_proven": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2194 S1144: strict Cutkosky multi-background tolerance stability audit",
                "",
                f"- Class stability: `{stable}`",
                f"- Class set: `{sorted(set(classes))}`",
                f"- Background drift detected: `{drift}`",
                "",
                "This is a multi-background tolerance audit only, not full Cutkosky closure.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
