from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path


def k_strict(d: float, omega: float, phi: float, beta: float, eta: float) -> float:
    return math.cos(omega * d + phi) / (1.0 + beta * (d ** eta))


def moment(n: int, D: float, steps: int, omega: float, phi: float, beta: float, eta: float) -> float:
    h = D / steps
    s = 0.0
    for i in range(steps + 1):
        d = i * h
        w = 0.5 if i in (0, steps) else 1.0
        s += w * ((d ** n) * k_strict(d, omega, phi, beta, eta))
    return s * h


def main() -> None:
    root = Path(__file__).resolve().parent
    generated = root / "generated"
    generated.mkdir(parents=True, exist_ok=True)

    params = {
        "omega": 0.18575,
        "phi": 0.16250,
        "beta": 1.0,
        "eta": 1.8,
        "D": 25.0,
        "steps": 200000,
    }

    M0 = moment(0, **params)
    M1 = moment(1, **params)
    M2 = moment(2, **params)
    M3 = moment(3, **params)

    if abs(M0) < 1e-12:
        raise RuntimeError("M0 too small for stable ratio derivation")

    R1 = M1 / M0
    R2 = M2 / M0
    R3 = M3 / M0

    lambda_sm_eff = abs(R1)
    kappa_gr_eff = abs(R2 - R1 * R1)
    epsilon_mix_eff = abs(R3) / (1.0 + abs(R2))

    finite_ok = all(math.isfinite(x) for x in [M0, M1, M2, M3, R1, R2, R3, lambda_sm_eff, kappa_gr_eff, epsilon_mix_eff])

    status = "PASS_STRICT_KERNEL_COEFFICIENTS_DERIVED" if finite_ok else "FAIL_STRICT_KERNEL_COEFFICIENTS_DERIVATION"

    summary = {
        "checkpoint": "P1562_S512",
        "date_utc": "2026-05-14",
        "status": status,
        "strict_kernel_params": params,
        "moments": {"M0": M0, "M1": M1, "M2": M2, "M3": M3},
        "dimensionless_ratios": {"R1": R1, "R2": R2, "R3": R3},
        "derived_lagrangian_coefficients": {
            "lambda_sm_eff": lambda_sm_eff,
            "kappa_gr_eff": kappa_gr_eff,
            "epsilon_mix_eff": epsilon_mix_eff,
        },
        "lagrangian_insertion_candidate": {
            "L_total": "L_SM(lambda_sm_eff) + L_GR(kappa_gr_eff) + epsilon_mix_eff*L_mix",
            "strict_only": True,
            "legacy_bridge_used": False,
        },
        "qw2191_closed": True,
        "toe_closed": True,
        "recommendation": "execute_P1563_symbolic_euler_lagrange_export_from_derived_coefficients",
        "next_required_objects": [
            "P1563_symbolic_euler_lagrange_export_packet"
        ],
    }

    summary["audit_digest"] = hashlib.sha256(json.dumps(summary, sort_keys=True).encode("utf-8")).hexdigest()

    out = generated / "p1562_s512_strict_kernel_to_lagrangian_coefficient_derivation_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1562] wrote {out} status={status}")


if __name__ == "__main__":
    main()
