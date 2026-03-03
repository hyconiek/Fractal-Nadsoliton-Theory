#!/usr/bin/env python3
"""
QW-1712: IR Lorentz recovery test for FIN effective dispersion.

Defines a minimal effective dispersion model and quantifies:
Delta_Lorentz(mu) = |c_phase(mu)/c - 1|
across IR/UV scales, including sensitivity to beta_tors.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1712_ir_lorentz_recovery_test.json"
OUT_MD = ROOT / "RAPORT_QW1712_IR_LORENTZ_RECOVERY_TEST.md"


def c_phase_ratio(mu: np.ndarray, beta_tors: float, alpha_geo: float) -> np.ndarray:
    """
    Effective phase velocity ratio c_phase/c from a FIN-style EFT expansion:
    w^2 = c^2 k^2 [1 + xi2*mu^2 + xi4*mu^4],  mu = k/Lambda
    """
    xi2 = beta_tors * alpha_geo / 12.0
    xi4 = (beta_tors ** 2) * alpha_geo / 4.0
    return np.sqrt(np.maximum(1.0 + xi2 * mu**2 + xi4 * mu**4, 0.0))


def main() -> None:
    alpha_geo = 4.0 * math.log(2.0)
    beta0 = 0.01

    mu = np.logspace(-6, 0, 2000)  # IR -> UV proxy

    # Baseline
    c_ratio = c_phase_ratio(mu, beta0, alpha_geo)
    delta = np.abs(c_ratio - 1.0)

    # IR criterion window
    mu_ir = 1e-3
    ir_mask = mu <= mu_ir
    delta_ir_max = float(np.max(delta[ir_mask]))
    delta_ir_mean = float(np.mean(delta[ir_mask]))

    # Sensitivity in beta_tors
    beta_scan = [0.005, 0.01, 0.02, 0.05]
    sensitivity = []
    for b in beta_scan:
        d = np.abs(c_phase_ratio(mu, b, alpha_geo) - 1.0)
        sensitivity.append(
            {
                "beta_tors": b,
                "delta_ir_max": float(np.max(d[ir_mask])),
                "delta_ir_mean": float(np.mean(d[ir_mask])),
                "delta_uv_at_mu1": float(d[-1]),
            }
        )

    # Closure criterion for IR compatibility
    threshold_ir = 1e-3
    ir_ok = delta_ir_max < threshold_ir
    verdict = "IR_LORENTZ_RECOVERY_OK" if ir_ok else "IR_LORENTZ_RECOVERY_NOT_OK"

    # Additional diagnostics
    sample_points = [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0]
    point_rows = []
    for m in sample_points:
        d = abs(c_phase_ratio(np.array([m]), beta0, alpha_geo)[0] - 1.0)
        point_rows.append({"mu": m, "delta_lorentz": float(d)})

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "params": {
            "alpha_geo": alpha_geo,
            "beta_tors_baseline": beta0,
        },
        "definition": {
            "delta_lorentz": "|c_phase/c - 1|",
            "mu": "k/Lambda",
            "mu_ir_threshold": mu_ir,
        },
        "baseline": {
            "delta_ir_max": delta_ir_max,
            "delta_ir_mean": delta_ir_mean,
            "delta_uv_at_mu1": float(delta[-1]),
        },
        "sensitivity_beta_scan": sensitivity,
        "sample_points": point_rows,
        "thresholds": {
            "ir_max_threshold": threshold_ir,
        },
        "verdict": verdict,
    }

    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1712: IR LORENTZ RECOVERY TEST",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Werdykt: **{verdict}**",
        "",
        "## 1) Baseline (beta_tors = 0.01)",
        f"- max Delta_Lorentz (mu <= 1e-3): {delta_ir_max:.6e}",
        f"- mean Delta_Lorentz (mu <= 1e-3): {delta_ir_mean:.6e}",
        f"- Delta_Lorentz (mu=1): {delta[-1]:.6e}",
        "",
        "## 2) Sensitivity beta_tors",
    ]
    for row in sensitivity:
        lines.append(
            f"- beta={row['beta_tors']:.3f}: delta_ir_max={row['delta_ir_max']:.6e}, delta_ir_mean={row['delta_ir_mean']:.6e}, delta_uv(mu=1)={row['delta_uv_at_mu1']:.6e}"
        )
    lines.extend(
        [
            "",
            "## 3) Interpretacja",
            "- Test mówi, czy model EFT zachowuje zgodność Lorentzowską w limicie IR.",
            "- Nie jest to pełny dowód relatywistycznej kompletności; to test konieczny, nie wystarczający.",
            "",
            "## Artefakty",
            f"- JSON szczegółowy: `{OUT_JSON.name}`",
        ]
    )

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"[QW-1712] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1712] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
