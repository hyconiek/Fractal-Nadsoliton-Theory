#!/usr/bin/env python3
"""
QW-1939: Hard mass-formula baseline (no correction terms, no fitting).

Exact tested law:
    m(Q) = m_top * 4^(-(gamma * Q / 4))

Protocol:
- no Delta correction,
- no train/test fitting,
- gamma either fixed canonical (1.52) or derived directly from frozen kernel.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1939_hard_mass_formula_baseline.json"
OUT_MD = ROOT / "RAPORT_QW1939_HARD_MASS_FORMULA_BASELINE.md"


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


def rel_err_pct(pred: float, exp: float) -> float:
    return float(abs(pred - exp) / max(abs(exp), 1e-15) * 100.0)


def run_variant(gamma: float, label: str) -> Dict[str, object]:
    rows: List[Dict[str, float]] = []
    errs = []
    for name, q, m_exp in PARTICLES:
        pred = 173_000.0 * (4.0 ** (-(gamma * q / 4.0)))
        err = rel_err_pct(pred, m_exp)
        errs.append(err)
        rows.append(
            {
                "particle": name,
                "q": float(q),
                "exp_mev": float(m_exp),
                "pred_mev": float(pred),
                "rel_err_pct": float(err),
            }
        )
    return {
        "label": label,
        "gamma": float(gamma),
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

    k1 = float(abs(kernel_fn(np.array([1.0]), omega, phi, beta, eta)[0]))
    k2 = float(abs(kernel_fn(np.array([2.0]), omega, phi, beta, eta)[0]))
    k4 = float(abs(kernel_fn(np.array([4.0]), omega, phi, beta, eta)[0]))

    # Direct kernel-derived exponents (no fit).
    gamma_12 = float(-4.0 * math.log(max(k2 / max(k1, 1e-15), 1e-15), 4.0))
    gamma_14 = float(-4.0 * math.log(max(k4 / max(k1, 1e-15), 1e-15), 4.0) / 3.0)

    variants = [
        run_variant(1.52, "canonical_gamma_1p52"),
        run_variant(gamma_12, "kernel_derived_gamma_1to2"),
        run_variant(gamma_14, "kernel_derived_gamma_1to4"),
    ]

    # Primary hard-baseline criterion: kernel-derived gamma from d=1..4.
    primary = next(v for v in variants if v["label"] == "kernel_derived_gamma_1to4")
    thresholds = {
        "mean_rel_err_pct_max": 15.0,
        "max_rel_err_pct_max": 75.0,
    }
    hard_pass = bool(
        primary["mean_rel_err_pct"] <= thresholds["mean_rel_err_pct_max"]
        and primary["max_rel_err_pct"] <= thresholds["max_rel_err_pct_max"]
    )

    verdict = "HARD_MASS_FORMULA_BASELINE_PASS" if hard_pass else "HARD_MASS_FORMULA_BASELINE_FAIL"
    required_next = (
        "USE_HARD_MASS_FORMULA_AS_LOCKED_SECTOR_INPUT"
        if hard_pass
        else "MASS_SECTOR_REMAINS_OPEN_UNDER_EXACT_FORMULA"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "kernel_source": "report_qw1932_physical_reparameterization_eta_scan.json:selected",
        "kernel": {"omega": omega, "phi": phi, "beta": beta, "eta": eta},
        "kernel_values": {"abs_k1": k1, "abs_k2": k2, "abs_k4": k4},
        "variants": variants,
        "primary_variant": primary["label"],
        "thresholds": thresholds,
        "hard_pass": hard_pass,
        "verdict": verdict,
        "required_next_step": required_next,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1939: HARD MASS FORMULA BASELINE",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Kernel: omega={omega:.6f}, phi={phi:.6f}, beta={beta:.6f}, eta={eta:.2f}",
        f"- Verdict: **{verdict}**",
        "",
        "## Exact Formula",
        "- m(Q) = m_top * 4^(-(gamma * Q / 4))",
        "- No Delta correction, no fitting.",
        "",
        "## Variants",
    ]
    for v in variants:
        lines.append(
            f"- {v['label']}: gamma={v['gamma']:.6f}, "
            f"mean/max rel%={v['mean_rel_err_pct']:.3f}/{v['max_rel_err_pct']:.3f}"
        )
    lines.extend(
        [
            "",
            "## Primary Hard Baseline",
            f"- variant: {primary['label']}",
            f"- hard_pass: {hard_pass}",
            "",
            "## Required Next Step",
            f"- {required_next}",
            "",
            "## Artifacts",
            f"- JSON: `{OUT_JSON.name}`",
        ]
    )
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1939] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1939] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1939] verdict={verdict}")


if __name__ == "__main__":
    main()
