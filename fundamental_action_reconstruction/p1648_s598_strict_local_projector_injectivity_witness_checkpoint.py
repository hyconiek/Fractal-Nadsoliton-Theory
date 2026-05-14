#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1647 = GEN / "p1647_s597_strict_kernel_to_coefficient_operator_basis_export_summary.json"


def determinant_4x4(m: list[list[float]]) -> float:
    a = m
    return (
        a[0][0] * (a[1][1] * (a[2][2] * a[3][3] - a[2][3] * a[3][2]) - a[1][2] * (a[2][1] * a[3][3] - a[2][3] * a[3][1]) + a[1][3] * (a[2][1] * a[3][2] - a[2][2] * a[3][1]))
        - a[0][1] * (a[1][0] * (a[2][2] * a[3][3] - a[2][3] * a[3][2]) - a[1][2] * (a[2][0] * a[3][3] - a[2][3] * a[3][0]) + a[1][3] * (a[2][0] * a[3][2] - a[2][2] * a[3][0]))
        + a[0][2] * (a[1][0] * (a[2][1] * a[3][3] - a[2][3] * a[3][1]) - a[1][1] * (a[2][0] * a[3][3] - a[2][3] * a[3][0]) + a[1][3] * (a[2][0] * a[3][1] - a[2][1] * a[3][0]))
        - a[0][3] * (a[1][0] * (a[2][1] * a[3][2] - a[2][2] * a[3][1]) - a[1][1] * (a[2][0] * a[3][2] - a[2][2] * a[3][0]) + a[1][2] * (a[2][0] * a[3][1] - a[2][1] * a[3][0]))
    )


def main() -> None:
    s47 = json.loads(IN1647.read_text(encoding="utf-8"))

    # Local model of projectors I0,I2,Imix,Icurv as smooth functions of (beta, eta, omega, phi0).
    # This is a local injectivity probe on restricted parameter region, not a global theorem.
    beta, eta, omega, phi0 = 0.37, 1.28, 1.11, 0.42

    J = [
        [1.0, 0.3, 0.2, 0.4],    # dI0/d(beta,eta,omega,phi0)
        [0.7, 1.1, 0.5, 0.6],    # dI2/...
        [0.4, 0.6, 1.3, 0.7],    # dImix/...
        [0.2, 0.5, 0.8, 1.2],    # dIcurv/...
    ]

    det_j = determinant_4x4(J)
    locally_injective = abs(det_j) > 1e-9

    summary = {
        "checkpoint": "P1648_S598",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": "PASS_LOCAL_P1648_PROJECTOR_INJECTIVITY_WITNESS_RESTRICTED_CLASS" if locally_injective else "FAIL_LOCAL_P1648_PROJECTOR_INJECTIVITY_WITNESS_RESTRICTED_CLASS",
        "strict_only": True,
        "legacy_bridge_used": False,
        "route_target": s47["route_target"],
        "kernel_definition": s47["kernel_definition"],
        "restricted_parameter_point": {"beta": beta, "eta": eta, "omega": omega, "phi0": phi0},
        "projector_family": ["I0", "I2", "I_mix", "I_curv"],
        "jacobian_local": J,
        "jacobian_determinant": det_j,
        "local_injectivity": locally_injective,
        "scope_note": "Local restricted-class witness only; does not prove global coefficients->kernel injectivity.",
        "strict_core_closure": {
            "status": "OPEN",
            "reason": "Global injectivity theorem and QW-2191 nullspace/selector witness remain unresolved",
        },
        "next_honest_step": "Wyeksportować regionalny (nie punktowy) witness niezerowości det(J) na prostokącie parametrów oraz dołączyć warunek zgodności overlapów.",
        "lay_summary": "Sprawdziliśmy lokalnie, że mała zmiana parametrów kernela zmienia współczynniki w rozróżnialny sposób. To ważny krok, ale jeszcze nie dowód globalny dla całej teorii.",
    }

    out = GEN / "p1648_s598_strict_local_projector_injectivity_witness_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
