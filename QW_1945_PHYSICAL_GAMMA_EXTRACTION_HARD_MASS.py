#!/usr/bin/env python3
"""
QW-1945: Physical gamma extraction from kernel for hard mass formula.

Mass law tested:
    m(Q) = m_top * 4^(-(gamma * Q / 4))

No correction terms, no fitting to masses.
Gamma is extracted from kernel decay with fixed physical recipes.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1945_physical_gamma_extraction_hard_mass.json"
OUT_MD = ROOT / "RAPORT_QW1945_PHYSICAL_GAMMA_EXTRACTION_HARD_MASS.md"


PARTICLES = [
    ("Top", 0.0, 173_000.0),
    ("Bottom", 7.0, 4_180.0),
    ("Tau", 9.0, 1_776.9),
    ("Charm", 9.0, 1_270.0),
    ("Muon", 14.0, 105.7),
    ("Electron", 24.0, 0.511),
]


def kernel_fn(d: np.ndarray, omega: float, phi: float, beta: float, eta: float) -> np.ndarray:
    return np.cos(omega * d + phi) / (1.0 + beta * (d**eta))


def gamma_from_two_points(k_a: float, k_b: float, d_a: float, d_b: float) -> float:
    slope = math.log(max(k_b, 1e-15) / max(k_a, 1e-15)) / max(d_b - d_a, 1e-15)
    return float(-4.0 * slope / math.log(4.0))


def gamma_wls_short_range(d: np.ndarray, k_abs: np.ndarray, d_max: float = 6.0) -> Dict[str, float]:
    m = d <= d_max
    x = d[m]
    y = np.log(np.clip(k_abs[m], 1e-15, None))
    w = np.clip(k_abs[m], 1e-12, None)

    x_bar = float(np.sum(w * x) / np.sum(w))
    y_bar = float(np.sum(w * y) / np.sum(w))
    num = float(np.sum(w * (x - x_bar) * (y - y_bar)))
    den = float(np.sum(w * (x - x_bar) ** 2))
    slope = num / max(den, 1e-15)
    intercept = y_bar - slope * x_bar

    y_hat = intercept + slope * x
    ss_res = float(np.sum(w * (y - y_hat) ** 2))
    ss_tot = float(np.sum(w * (y - y_bar) ** 2))
    r2 = float(1.0 - ss_res / max(ss_tot, 1e-15))
    gamma = float(-4.0 * slope / math.log(4.0))
    return {
        "gamma": gamma,
        "slope": slope,
        "intercept": intercept,
        "r2_weighted": r2,
    }


def mass_errors(gamma: float) -> Dict[str, object]:
    rows = []
    errs = []
    for name, q, m_exp in PARTICLES:
        pred = 173_000.0 * (4.0 ** (-(gamma * q / 4.0)))
        err = abs(pred - m_exp) / max(m_exp, 1e-15) * 100.0
        errs.append(err)
        rows.append(
            {
                "particle": name,
                "q": q,
                "exp_mev": m_exp,
                "pred_mev": float(pred),
                "rel_err_pct": float(err),
            }
        )
    return {
        "rows": rows,
        "mean_rel_err_pct": float(np.mean(errs)),
        "max_rel_err_pct": float(np.max(errs)),
    }


def main() -> None:
    d1932 = json.loads((ROOT / "report_qw1932_physical_reparameterization_eta_scan.json").read_text(encoding="utf-8"))
    sel = d1932["selected"]
    omega = float(sel["fit"]["omega"])
    phi = float(sel["fit"]["phi"])
    beta = float(sel["fit"]["beta"])
    eta = float(sel["eta"])

    d = np.arange(1.0, 13.0, dtype=float)
    k_abs = np.abs(kernel_fn(d, omega, phi, beta, eta))

    g12 = gamma_from_two_points(float(k_abs[0]), float(k_abs[1]), 1.0, 2.0)
    g14 = gamma_from_two_points(float(k_abs[0]), float(k_abs[3]), 1.0, 4.0)
    gwls = gamma_wls_short_range(d, k_abs, d_max=6.0)

    variants = [
        {"label": "local_1_to_2", "gamma": g12, "fit_quality": None},
        {"label": "local_1_to_4", "gamma": g14, "fit_quality": None},
        {"label": "wls_short_range_1_to_6", "gamma": gwls["gamma"], "fit_quality": gwls},
    ]

    for v in variants:
        m = mass_errors(v["gamma"])
        v["mass"] = m

    primary = next(v for v in variants if v["label"] == "wls_short_range_1_to_6")
    thresholds = {
        "mean_rel_err_pct_max": 15.0,
        "max_rel_err_pct_max": 75.0,
    }
    primary_pass = bool(
        primary["mass"]["mean_rel_err_pct"] <= thresholds["mean_rel_err_pct_max"]
        and primary["mass"]["max_rel_err_pct"] <= thresholds["max_rel_err_pct_max"]
    )

    verdict = "PHYSICAL_GAMMA_EXTRACTION_HARD_MASS_PASS" if primary_pass else "PHYSICAL_GAMMA_EXTRACTION_HARD_MASS_FAIL"
    required_next = (
        "USE_PRIMARY_GAMMA_IN_FINAL_SINGLE_KERNEL_GATE"
        if primary_pass
        else "MASS_STILL_OPEN_UNDER_PHYSICAL_GAMMA_EXTRACTION"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "kernel_source": "report_qw1932_physical_reparameterization_eta_scan.json:selected",
        "kernel": {"omega": omega, "phi": phi, "beta": beta, "eta": eta},
        "k_abs_d1_to_d12": [float(x) for x in k_abs],
        "variants": variants,
        "primary_variant": primary["label"],
        "thresholds": thresholds,
        "primary_pass": primary_pass,
        "verdict": verdict,
        "required_next_step": required_next,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1945: PHYSICAL GAMMA EXTRACTION HARD MASS",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        "",
        "## Primary Rule",
        "- weighted log-decay fit on d=1..6",
        f"- gamma: {primary['gamma']:.6f}",
    ]
    if primary.get("fit_quality") is not None:
        lines.append(f"- weighted R2: {primary['fit_quality']['r2_weighted']:.6f}")
    lines.extend(
        [
            (
                "- primary mass mean/max rel%: "
                f"{primary['mass']['mean_rel_err_pct']:.3f}/{primary['mass']['max_rel_err_pct']:.3f}"
            ),
            "",
            "## All Variants",
        ]
    )
    for v in variants:
        lines.append(
            f"- {v['label']}: gamma={v['gamma']:.6f}, "
            f"mean/max rel%={v['mass']['mean_rel_err_pct']:.3f}/{v['mass']['max_rel_err_pct']:.3f}"
        )
    lines.extend(
        [
            "",
            "## Required Next Step",
            f"- {required_next}",
            "",
            "## Artifacts",
            f"- JSON: `{OUT_JSON.name}`",
        ]
    )
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1945] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1945] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1945] verdict={verdict}")


if __name__ == "__main__":
    main()

