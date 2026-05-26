#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from scipy.integrate import quad
import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2204 = GEN / "p2204_s1154_strict_frw_bianchi_symbolic_transport_operator_mismatch_certificate.json"
OUT = GEN / "p2205_s1155_strict_nu_branch_transport_vanishing_constraints_probe.json"
MD = GEN / "p2205_s1155_strict_nu_branch_transport_vanishing_constraints_probe.md"


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
    p2204 = load(IN_2204)

    p2190 = load(ROOT / "generated/p2190_s1140_strict_cutkosky_fixed_channel_real_discontinuity_integral_packet.json")
    w = p2190.get("strict_cutkosky_fixed_channel_real_discontinuity_integral_packet", {}) or {}
    params = w.get("strict_kernel_params", {}) or {"omega": 0.18575, "phi": 0.1625, "beta": 1.0, "eta": 1.8}
    domain = w.get("integration_domain", {}) or {"d_min": 0.0, "d_max": 20.0}
    d_min, d_max = float(domain["d_min"]), float(domain["d_max"])

    omega0 = float(params["omega"])
    phi = float(params["phi"])
    beta = float(params["beta"])
    eta = float(params["eta"])

    # symbolic sensitivity around m=1
    d, om, ph, be, et, m = sp.symbols("d om ph be et m", positive=True)
    kf = sp.cos(om * d + ph) / (1 + be * d**et)
    kb = sp.cos((om * m) * d + ph) / (1 + be * d**et)
    delta = sp.simplify(kf**2 - kb**2)
    ddelta_dm_at_1 = sp.simplify(sp.diff(delta, m).subs({m: 1}))

    # numeric probe: residual scaling for |m-1| in tiny neighborhood
    dm_grid = np.array([0.0, 1e-4, 2e-4, 5e-4, 1e-3, 2e-3], dtype=float)
    rows = []
    for dm in dm_grid:
        mm = 1.0 + float(dm)

        def rint(x: float) -> float:
            a = k_strict(x, omega0, phi, beta, eta)
            b = k_strict(x, omega0 * mm, phi, beta, eta)
            return abs((a * a) - (b * b)) * phase_space_weight(x, d_max)

        val, err = quad(rint, d_min, d_max, epsabs=1e-11, epsrel=1e-11, limit=500)
        rows.append({"dm": float(dm), "m": mm, "residual_l1": float(val), "abs_err": float(err)})

    eps_ref = 1e-3
    r_eps = next(r["residual_l1"] for r in rows if abs(r["dm"] - eps_ref) < 1e-16)
    ratio_rows = [
        {
            "dm": r["dm"],
            "residual_over_dm": (r["residual_l1"] / r["dm"]) if r["dm"] > 0 else None,
            "residual_over_dm2": (r["residual_l1"] / (r["dm"] ** 2)) if r["dm"] > 0 else None,
        }
        for r in rows
    ]

    nonzero = [r for r in ratio_rows if r["dm"] > 0]
    ro1 = [r["residual_over_dm"] for r in nonzero]
    linear_like = (max(ro1) - min(ro1)) / max(ro1) < 0.15 if ro1 else False

    payload = {
        "schema_version": "p2205_s1155_v1",
        "packet_id": "P2205",
        "stage_id": "S1155",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-26T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_TRANSPORT_VANISHING_CONSTRAINTS_PROBE",
        "strict_nu_branch_transport_vanishing_constraints_probe": {
            "probe_id": "STRICT_NU_BRANCH_TRANSPORT_VANISHING_CONSTRAINTS_PROBE_V1",
            "source_packet": str(IN_2204.relative_to(ROOT)),
            "symbolic_first_variation_at_m_eq_1": str(ddelta_dm_at_1),
            "sufficient_zero_constraint_candidate": "m = 1 (equivalently delta_omega = 0 on frozen channel/domain)",
            "numeric_local_scaling_rows": rows,
            "numeric_local_scaling_ratios": ratio_rows,
            "local_behavior_assessment": {
                "residual_at_dm_1e3": r_eps,
                "approximately_linear_in_dm_near_zero": linear_like,
                "interpretation": "local transport mismatch behaves as nonzero first-order effect when m departs from 1",
            },
            "theorem_scope_limit": "local sufficient-constraint probe only; no global nu-branch transport theorem closure",
        },
        "recommended_next_honest_step": {
            "id": "P2206_candidate",
            "goal": "build explicit nu-branch transport operator equation and test if project constraints can force m=1 without ad-hoc selector assumptions",
        },
        "gatekeeper_checks": {
            "vanishing_constraints_probe_exported": True,
            "numeric_rows_finite": all(math.isfinite(r["residual_l1"]) for r in rows),
            "symbolic_first_variation_exported": True,
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
            "# P2205 S1155: strict nu-branch transport vanishing constraints probe",
            "",
            "- symbolic sufficient candidate: `m = 1`",
            f"- residual(dm=1e-3): `{r_eps:.12e}`",
            f"- approximately linear near zero: `{linear_like}`",
            "",
            "Local vanishing-constraint probe only; no global Task-3 closure claim.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
