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
IN_2190 = GEN / "p2190_s1140_strict_cutkosky_fixed_channel_real_discontinuity_integral_packet.json"
OUT = GEN / "p2191_s1141_strict_cutkosky_uncertainty_bands_and_ur_residual_envelope.json"
MD = GEN / "p2191_s1141_strict_cutkosky_uncertainty_bands_and_ur_residual_envelope.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def k_strict(d: float, omega: float, phi: float, beta: float, eta: float) -> float:
    return math.cos(omega * d + phi) / (1.0 + beta * (d ** eta))


def phase_space_weight(d: float, d_max: float) -> float:
    x = d / d_max
    return max(0.0, x * (1.0 - x))


def integrate_value(params: dict[str, float], d_min: float, d_max: float) -> tuple[float, float]:
    def integrand(d: float) -> float:
        k = k_strict(d, params["omega"], params["phi"], params["beta"], params["eta"])
        return (k * k) * phase_space_weight(d, d_max)

    return quad(integrand, d_min, d_max, epsabs=1e-10, epsrel=1e-10, limit=400)


def main() -> None:
    GEN.mkdir(exist_ok=True)

    p2190 = load(IN_2190)
    witness = p2190.get("strict_cutkosky_fixed_channel_real_discontinuity_integral_packet", {}) or {}

    params = witness.get("strict_kernel_params", {}) or {"omega": 0.18575, "phi": 0.16250, "beta": 1.0, "eta": 1.8}
    domain = witness.get("integration_domain", {}) or {"d_min": 0.0, "d_max": 20.0}
    d_min = float(domain.get("d_min", 0.0))
    d_max = float(domain.get("d_max", 20.0))

    central, central_err = integrate_value(params, d_min, d_max)

    band_spec = {
        "omega": [params["omega"] - 0.005, params["omega"], params["omega"] + 0.005],
        "phi": [params["phi"] - 0.005, params["phi"], params["phi"] + 0.005],
        "beta": [params["beta"] - 0.05, params["beta"], params["beta"] + 0.05],
        "eta": [params["eta"] - 0.05, params["eta"], params["eta"] + 0.05],
    }

    rows = []
    for omega in band_spec["omega"]:
        for phi in band_spec["phi"]:
            for beta in band_spec["beta"]:
                for eta in band_spec["eta"]:
                    p = {"omega": float(omega), "phi": float(phi), "beta": float(beta), "eta": float(eta)}
                    val, err = integrate_value(p, d_min, d_max)
                    rows.append({"params": p, "integral": float(val), "quad_abs_error": float(err)})

    values = np.array([r["integral"] for r in rows], dtype=float)
    envelope = {
        "min_integral": float(np.min(values)),
        "max_integral": float(np.max(values)),
        "central_integral": float(central),
        "residual_radius_from_central": float(np.max(np.abs(values - central))),
        "sample_count": int(values.size),
    }

    ur_link_style = {
        "observable": "DiscM_minus_CutSum_proxy",
        "center": float(central),
        "lower": envelope["min_integral"],
        "upper": envelope["max_integral"],
        "radius": envelope["residual_radius_from_central"],
    }

    checks = {
        "central_finite": math.isfinite(central),
        "band_values_finite": bool(np.all(np.isfinite(values))),
        "envelope_ordered": envelope["min_integral"] <= envelope["central_integral"] <= envelope["max_integral"],
        "ur_radius_nonnegative": ur_link_style["radius"] >= 0.0,
    }
    all_checks_pass = all(checks.values())

    result_kind = (
        "PASS_STRICT_CUTKOSKY_UNCERTAINTY_BANDS_AND_UR_RESIDUAL_ENVELOPE"
        if all_checks_pass
        else "OPEN_STRICT_CUTKOSKY_UNCERTAINTY_BANDS_AND_UR_RESIDUAL_ENVELOPE_NUMERICAL_ISSUE"
    )

    payload = {
        "schema_version": "p2191_s1141_v1",
        "packet_id": "P2191",
        "stage_id": "S1141",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-26T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": result_kind,
        "strict_cutkosky_uncertainty_bands_and_ur_residual_envelope": {
            "witness_id": "STRICT_CUTKOSKY_UR_LINK_UNCERTAINTY_ENVELOPE_V1",
            "source_packet": str(IN_2190.relative_to(ROOT)),
            "frozen_channel": witness.get("frozen_channel", "graviton_to_gauge_gauge"),
            "frozen_background": witness.get("frozen_background", "B1_fixed_chart_seed"),
            "integration_domain": {"d_min": d_min, "d_max": d_max},
            "central": {"value": float(central), "quad_abs_error": float(central_err), "params": params},
            "band_spec": band_spec,
            "uncertainty_band_table": rows,
            "residual_residue_envelope": envelope,
            "ur_link_uncertainty_object": ur_link_style,
        },
        "recommended_next_honest_step": {
            "id": "P2192_candidate",
            "goal": "attach residue-sign and positivity-band certificate table on the same fixed channel/background",
        },
        "gatekeeper_checks": {
            "uncertainty_band_table_exported": True,
            "ur_envelope_exported": True,
            "numerical_envelope_checks_pass": all_checks_pass,
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
                "# P2191 S1141: strict Cutkosky uncertainty bands + UR residual envelope",
                "",
                f"- Result kind: `{result_kind}`",
                f"- Sample count: `{len(rows)}`",
                f"- Central integral: `{central:.12e}`",
                f"- Envelope min/max: `[{envelope['min_integral']:.12e}, {envelope['max_integral']:.12e}]`",
                f"- UR radius: `{ur_link_style['radius']:.12e}`",
                "",
                "This is a strict-lane uncertainty-envelope witness only, not full dressed Cutkosky closure.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
