#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1646 = GEN / "p1646_s596_strict_gauge_metric_cross_variation_scaffold_summary.json"


def main() -> None:
    s46 = json.loads(IN1646.read_text(encoding="utf-8"))

    kernel_def = "K_strict(d) = cos(omega*d + phi0)/(1 + beta*d^eta)"

    coeff_map = {
        "c_kin_phi": "I0[K_strict]",
        "lambda_phi": "I2[K_strict]",
        "lambda_phiH": "I_mix[K_strict, H]",
        "xi_phiR": "I_curv[K_strict]",
        "Mpl2_half": "I_grav0[K_strict]",
        "c1_R2": "I_grav2a[K_strict]",
        "c2_Ricci2": "I_grav2b[K_strict]",
        "c3_Riem2": "I_grav2c[K_strict]",
    }

    integral_semantics = {
        "I0": "normalized even-moment projection onto scalar kinetic operator",
        "I2": "quartic scalar self-interaction moment",
        "I_mix": "mixed projector onto phi^2 H^†H sector",
        "I_curv": "scalar-curvature coupling projector",
        "I_grav0": "Einstein-Hilbert-scale projector",
        "I_grav2a/b/c": "independent quadratic-curvature projectors",
    }

    reverse_obligation = {
        "target": "coefficients -> K_strict",
        "status": "OPEN",
        "missing": [
            "injective theorem for projector family {I*}",
            "nullspace elimination witness compatible with QW-2191",
            "global overlap consistency of recovered kernel parameters",
        ],
    }

    summary = {
        "checkpoint": "P1647_S597",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": "KEEP_OPEN_P1647_KERNEL_TO_COEFFICIENT_OPERATOR_BASIS_EXPORTED",
        "strict_only": True,
        "legacy_bridge_used": False,
        "route_target": s46["route_target"],
        "kernel_definition": kernel_def,
        "coefficient_operator_basis_map": coeff_map,
        "projector_semantics": integral_semantics,
        "full_lagrangian_reference": s46["full_lagrangian_density"],
        "reverse_obligation_coeff_to_kernel": reverse_obligation,
        "strict_core_closure": {
            "status": "OPEN",
            "reason": "Coefficient->kernel injectivity and selector uniqueness witnesses still missing",
        },
        "next_honest_step": "Wyeksportować lokalny witness iniektywności dla podrodziny {I0, I2, I_mix, I_curv} na ograniczonej klasie parametrów (beta, eta, omega, phi0).",
        "lay_summary": "Dodaliśmy mapę: jak z kernela strict wyprowadzać współczynniki do pełnego lagranżianu. Nadal brakuje formalnego dowodu, że tę mapę da się odwrócić jednoznacznie w całej teorii.",
    }

    out = GEN / "p1647_s597_strict_kernel_to_coefficient_operator_basis_export_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
