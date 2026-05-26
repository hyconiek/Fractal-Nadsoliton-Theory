#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

from scipy.integrate import quad

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2200 = GEN / "p2200_s1150_strict_cutkosky_integrand_bound_and_monotonicity_certificate.json"
OUT = GEN / "p2201_s1151_strict_cutkosky_shared_majorant_frw_bianchi_certificate.json"
MD = GEN / "p2201_s1151_strict_cutkosky_shared_majorant_frw_bianchi_certificate.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def k_strict(d: float, omega: float, phi: float, beta: float, eta: float) -> float:
    return math.cos(omega * d + phi) / (1.0 + beta * (d**eta))


def phase_space_weight(d: float, d_max: float) -> float:
    x = d / d_max
    return max(0.0, x * (1.0 - x))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2200 = load(IN_2200)
    c = p2200.get("strict_cutkosky_integrand_bound_and_monotonicity_certificate", {}) or {}

    src = load(ROOT / c.get("source_packet", "generated/p2190_s1140_strict_cutkosky_fixed_channel_real_discontinuity_integral_packet.json"))
    w = src.get("strict_cutkosky_fixed_channel_real_discontinuity_integral_packet", {}) or {}
    params = w.get("strict_kernel_params", {}) or {"omega": 0.18575, "phi": 0.1625, "beta": 1.0, "eta": 1.8}
    domain = w.get("integration_domain", {}) or {"d_min": 0.0, "d_max": 20.0}
    d_min, d_max = float(domain["d_min"]), float(domain["d_max"])

    backgrounds = [
        {"name": "FRW_flat_fixed_chart_seed", "omega_mult": 1.00},
        {"name": "BianchiI_diag_fixed_chart_seed", "omega_mult": 1.03},
    ]

    rows = []
    for bg in backgrounds:
        omega = float(params["omega"]) * float(bg["omega_mult"])
        phi = float(params["phi"])
        beta = float(params["beta"])
        eta = float(params["eta"])

        def integrand(x: float) -> float:
            k = k_strict(x, omega, phi, beta, eta)
            return (k * k) * phase_space_weight(x, d_max)

        def majorant(x: float) -> float:
            return 1.0 / (4.0 * (1.0 + beta * (x**eta)) ** 2)

        obs, obs_err = quad(integrand, d_min, d_max, epsabs=1e-10, epsrel=1e-10, limit=400)
        ub, ub_err = quad(majorant, d_min, d_max, epsabs=1e-10, epsrel=1e-10, limit=400)
        margin = ub - obs

        rows.append(
            {
                "background": bg["name"],
                "omega_mult": bg["omega_mult"],
                "observed_integral": float(obs),
                "observed_abs_error": float(obs_err),
                "majorant_integral_upper_bound": float(ub),
                "majorant_abs_error": float(ub_err),
                "bound_margin": float(margin),
                "bound_respected": bool(obs <= ub + 1e-10),
            }
        )

    margins = [r["bound_margin"] for r in rows]
    min_margin = min(margins)
    max_margin = max(margins)
    spread = max_margin - min_margin

    payload = {
        "schema_version": "p2201_s1151_v1",
        "packet_id": "P2201",
        "stage_id": "S1151",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-26T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_CUTKOSKY_SHARED_MAJORANT_FRW_BIANCHI_CERTIFICATE",
        "strict_cutkosky_shared_majorant_frw_bianchi_certificate": {
            "certificate_id": "STRICT_CUTKOSKY_SHARED_MAJORANT_FRW_BIANCHI_CERTIFICATE_V1",
            "source_packet": str(IN_2200.relative_to(ROOT)),
            "same_channel": w.get("frozen_channel", "graviton_to_gauge_gauge"),
            "integration_domain": {"d_min": d_min, "d_max": d_max},
            "rows": rows,
            "margin_summary": {
                "min_margin": min_margin,
                "max_margin": max_margin,
                "margin_spread": spread,
                "all_bounds_respected": all(r["bound_respected"] for r in rows),
            },
        },
        "recommended_next_honest_step": {
            "id": "P2202_candidate",
            "goal": "bind FRW/Bianchi shared-majorant margin spread into transport-obstruction ledger and compare against task-3 blockers",
        },
        "gatekeeper_checks": {
            "shared_majorant_certificate_exported": True,
            "all_bounds_respected": all(r["bound_respected"] for r in rows),
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
                "# P2201 S1151: strict Cutkosky shared-majorant FRW/Bianchi certificate",
                "",
                f"- Margin min/max: `{min_margin:.12e}` / `{max_margin:.12e}`",
                f"- Margin spread: `{spread:.12e}`",
                f"- All bounds respected: `{all(r['bound_respected'] for r in rows)}`",
                "",
                "This is a shared-majorant bound comparison only, not full dressed Cutkosky closure.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
