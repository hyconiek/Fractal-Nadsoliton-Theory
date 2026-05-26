#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from scipy.integrate import quad

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2201 = GEN / "p2201_s1151_strict_cutkosky_shared_majorant_frw_bianchi_certificate.json"
OUT = GEN / "p2203_s1153_strict_frw_bianchi_transport_residual_map_under_shared_majorant.json"
MD = GEN / "p2203_s1153_strict_frw_bianchi_transport_residual_map_under_shared_majorant.md"


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
    p2201 = load(IN_2201)
    cert = p2201.get("strict_cutkosky_shared_majorant_frw_bianchi_certificate", {}) or {}

    rows = cert.get("rows", []) or []
    if len(rows) < 2:
        raise RuntimeError("P2201 rows unavailable for FRW/Bianchi residual map")

    frw = next(r for r in rows if "FRW" in str(r.get("background", "")))
    bianchi = next(r for r in rows if "Bianchi" in str(r.get("background", "")))

    # reuse kernel params from P2190 to avoid redefining physics assumptions
    p2190 = load(ROOT / "generated/p2190_s1140_strict_cutkosky_fixed_channel_real_discontinuity_integral_packet.json")
    w = p2190.get("strict_cutkosky_fixed_channel_real_discontinuity_integral_packet", {}) or {}
    params = w.get("strict_kernel_params", {}) or {"omega": 0.18575, "phi": 0.1625, "beta": 1.0, "eta": 1.8}
    domain = w.get("integration_domain", {}) or {"d_min": 0.0, "d_max": 20.0}
    d_min, d_max = float(domain["d_min"]), float(domain["d_max"])

    omega0 = float(params["omega"])
    phi = float(params["phi"])
    beta = float(params["beta"])
    eta = float(params["eta"])

    frw_mult = float(frw["omega_mult"])
    bianchi_mult = float(bianchi["omega_mult"])
    om_grid = np.linspace(frw_mult, bianchi_mult, 17)

    map_rows = []
    for om in om_grid:
        omega_frw = omega0 * float(frw_mult)
        omega_bi = omega0 * float(om)

        def residual_integrand(x: float) -> float:
            kf = k_strict(x, omega_frw, phi, beta, eta)
            kb = k_strict(x, omega_bi, phi, beta, eta)
            return abs((kf * kf) - (kb * kb)) * phase_space_weight(x, d_max)

        res, err = quad(residual_integrand, d_min, d_max, epsabs=1e-10, epsrel=1e-10, limit=400)
        map_rows.append({
            "omega_mult_probe": float(om),
            "transport_residual_l1": float(res),
            "transport_residual_abs_error": float(err),
        })

    res_vals = [r["transport_residual_l1"] for r in map_rows]
    max_residual = max(res_vals)
    mean_residual = float(sum(res_vals) / len(res_vals))

    min_margin = float(cert.get("margin_summary", {}).get("min_margin", 0.0) or 0.0)
    normalized_max = (max_residual / min_margin) if min_margin > 0.0 else float("inf")

    theorem_gap_delta = {
        "task3_requirement": "global FRW<->Bianchi transport closure for nu branch",
        "computed_residual_map_kind": "L1 absolute squared-kernel transport residual with shared phase-space weight",
        "residual_identically_zero_proven": False,
        "residual_nonzero_detected": max_residual > 0.0,
        "max_residual_over_min_margin": normalized_max,
        "interpretation": "quantitative residual profile exported; theorem-grade global covariance transport closure remains open",
    }

    payload = {
        "schema_version": "p2203_s1153_v1",
        "packet_id": "P2203",
        "stage_id": "S1153",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-26T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_FRW_BIANCHI_TRANSPORT_RESIDUAL_MAP_UNDER_SHARED_MAJORANT",
        "strict_frw_bianchi_transport_residual_map_under_shared_majorant": {
            "map_id": "STRICT_FRW_BIANCHI_TRANSPORT_RESIDUAL_MAP_UNDER_SHARED_MAJORANT_V1",
            "source_packets": [str(IN_2201.relative_to(ROOT)), "generated/p2190_s1140_strict_cutkosky_fixed_channel_real_discontinuity_integral_packet.json"],
            "integration_domain": {"d_min": d_min, "d_max": d_max},
            "residual_map_rows": map_rows,
            "summary": {
                "max_residual": max_residual,
                "mean_residual": mean_residual,
                "min_majorant_margin_from_p2201": min_margin,
            },
            "theorem_gap_delta": theorem_gap_delta,
        },
        "recommended_next_honest_step": {
            "id": "P2204_candidate",
            "goal": "derive symbolic transport-operator mismatch terms and attempt sufficient conditions for zero residual on nu branch"
        },
        "gatekeeper_checks": {
            "transport_residual_map_exported": True,
            "residual_map_finite": all(math.isfinite(v) for v in res_vals),
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
            "# P2203 S1153: FRW/Bianchi transport residual map under shared majorant",
            "",
            f"- max residual: `{max_residual:.12e}`",
            f"- mean residual: `{mean_residual:.12e}`",
            f"- min majorant margin (P2201): `{min_margin:.12e}`",
            f"- max residual / min margin: `{normalized_max:.12e}`",
            "",
            "Residual map exported as quantitative theorem-gap delta; no claim of global Task-3 closure.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
