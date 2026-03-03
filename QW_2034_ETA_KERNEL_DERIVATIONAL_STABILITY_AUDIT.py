#!/usr/bin/env python3
"""
QW-2034: Eta-kernel derivational stability audit.

Purpose:
- quantify how stable (omega, phi, beta, eta) are when re-derived from microprofiles,
- test whether the QW-2030 frozen kernel lies inside bootstrap derivational intervals.
"""

from __future__ import annotations

import importlib.util
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2034_eta_kernel_derivational_stability_audit.json"
OUT_MD = ROOT / "RAPORT_QW2034_ETA_KERNEL_DERIVATIONAL_STABILITY_AUDIT.md"


def load_qw2021_module():
    path = ROOT / "QW_2021_V2_ETA_OPERATOR_BETA_CONSTRAINT_SCAN.py"
    spec = importlib.util.spec_from_file_location("qw2021_mod_2034", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["qw2021_mod_2034"] = mod
    spec.loader.exec_module(mod)
    return mod


def load_target_kernel() -> Dict[str, float]:
    p = ROOT / "report_qw2030_final_stage_c_gate_combined_branch.json"
    d = json.loads(p.read_text(encoding="utf-8"))
    return {k: float(v) for k, v in d["kernel"].items() if k in {"omega", "phi", "beta", "eta"}}


def ci95(x: np.ndarray) -> Dict[str, float]:
    return {
        "q02p5": float(np.quantile(x, 0.025)),
        "q50": float(np.quantile(x, 0.50)),
        "q97p5": float(np.quantile(x, 0.975)),
        "mean": float(np.mean(x)),
        "std": float(np.std(x)),
    }


def main() -> None:
    mod = load_qw2021_module()
    target = load_target_kernel()

    d, y, w = mod.build_profiles(n_grid=[48, 64, 80, 96], seeds_per_n=6, dmax=24)
    n_profiles = int(y.shape[0])

    rng = np.random.default_rng(203401)
    n_boot = 32
    fits: List[Dict[str, float]] = []

    for _ in range(n_boot):
        idx = rng.integers(0, n_profiles, size=n_profiles)
        yb = y[idx, :]
        wb = w[idx]
        fit = mod.fit_global(d, yb, wb, beta_target=None, beta_scale=None, lambda_beta=0.0)["optimum"]
        fits.append(
            {
                "omega": float(fit["omega"]),
                "phi": float(fit["phi"]),
                "beta": float(fit["beta"]),
                "eta": float(fit["eta"]),
                "objective": float(fit["objective"]),
            }
        )

    arr = {k: np.array([f[k] for f in fits], dtype=float) for k in ["omega", "phi", "beta", "eta", "objective"]}
    stats = {k: ci95(v) for k, v in arr.items()}

    flags = {
        "target_omega_in_bootstrap_ci95": bool(stats["omega"]["q02p5"] <= target["omega"] <= stats["omega"]["q97p5"]),
        "target_phi_in_bootstrap_ci95": bool(stats["phi"]["q02p5"] <= target["phi"] <= stats["phi"]["q97p5"]),
        "target_beta_in_bootstrap_ci95": bool(stats["beta"]["q02p5"] <= target["beta"] <= stats["beta"]["q97p5"]),
        "target_eta_in_bootstrap_ci95": bool(stats["eta"]["q02p5"] <= target["eta"] <= stats["eta"]["q97p5"]),
        "beta_std_le_0p60": bool(stats["beta"]["std"] <= 0.60),
        "eta_std_le_0p50": bool(stats["eta"]["std"] <= 0.50),
    }
    pass_count = int(sum(1 for v in flags.values() if v))
    total_flags = int(len(flags))

    if pass_count >= 5:
        verdict = "ETA_KERNEL_DERIVATIONAL_STABILITY_PASS"
        readiness = "DERIVATIONAL_STABILITY_ACCEPTABLE"
    elif pass_count >= 4:
        verdict = "ETA_KERNEL_DERIVATIONAL_STABILITY_PARTIAL"
        readiness = "DERIVATIONAL_STABILITY_PARTIAL"
    else:
        verdict = "ETA_KERNEL_DERIVATIONAL_STABILITY_FAIL"
        readiness = "DERIVATIONAL_STABILITY_NOT_ACCEPTABLE"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source": {
            "microprofiles": "QW_2021_V2_ETA_OPERATOR_BETA_CONSTRAINT_SCAN.py:build_profiles",
            "fit_method": "QW_2021_V2_ETA_OPERATOR_BETA_CONSTRAINT_SCAN.py:fit_global",
            "target_kernel": "report_qw2030_final_stage_c_gate_combined_branch.json:kernel",
        },
        "config": {
            "n_grid": [48, 64, 80, 96],
            "seeds_per_n": 6,
            "dmax": 24,
            "n_bootstrap": n_boot,
            "lambda_beta": 0.0,
        },
        "target_kernel": target,
        "bootstrap_stats": stats,
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "readiness": readiness,
        "verdict": verdict,
        "required_next_step": (
            "LOCK_DERIVATIONAL_APPENDIX_FOR_FROZEN_KERNEL"
            if verdict == "ETA_KERNEL_DERIVATIONAL_STABILITY_PASS"
            else "TIGHTEN_MICRO_DERIVATION_AND_REPEAT"
        ),
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2034: ETA KERNEL DERIVATIONAL STABILITY AUDIT",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Readiness: **{readiness}**",
        f"- Verdict: **{verdict}**",
        f"- pass_count: {pass_count}/{total_flags}",
        "",
        "## Target Kernel (QW-2030)",
        f"- omega={target['omega']:.6f}, phi={target['phi']:.6f}, beta={target['beta']:.6f}, eta={target['eta']:.6f}",
        "",
        "## Bootstrap CI95",
        (
            f"- omega q02.5/q50/q97.5: {stats['omega']['q02p5']:.6f} / "
            f"{stats['omega']['q50']:.6f} / {stats['omega']['q97p5']:.6f}"
        ),
        (
            f"- phi q02.5/q50/q97.5: {stats['phi']['q02p5']:.6f} / "
            f"{stats['phi']['q50']:.6f} / {stats['phi']['q97p5']:.6f}"
        ),
        (
            f"- beta q02.5/q50/q97.5 (std): {stats['beta']['q02p5']:.6f} / "
            f"{stats['beta']['q50']:.6f} / {stats['beta']['q97p5']:.6f} ({stats['beta']['std']:.6f})"
        ),
        (
            f"- eta q02.5/q50/q97.5 (std): {stats['eta']['q02p5']:.6f} / "
            f"{stats['eta']['q50']:.6f} / {stats['eta']['q97p5']:.6f} ({stats['eta']['std']:.6f})"
        ),
        "",
        "## Flags",
    ]
    for k, v in flags.items():
        lines.append(f"- {k}: {v}")
    lines.extend(
        [
            "",
            "## Required Next Step",
            f"- {out['required_next_step']}",
            "",
            "## Artifacts",
            f"- JSON: `{OUT_JSON.name}`",
        ]
    )
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2034] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2034] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2034] readiness={readiness} verdict={verdict} pass={pass_count}/{total_flags}")


if __name__ == "__main__":
    main()
