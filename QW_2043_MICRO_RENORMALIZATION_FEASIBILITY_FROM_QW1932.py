#!/usr/bin/env python3
"""
QW-2043: Feasibility audit of strong renormalization from existing micro branch.

Uses QW-1932 strict-pass envelope to test whether required matching constants
(from QW-2042) are already supported by pre-existing micro-derivation branch,
without introducing a new sector retune.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2043_micro_renormalization_feasibility_from_qw1932.json"
OUT_MD = ROOT / "RAPORT_QW2043_MICRO_RENORMALIZATION_FEASIBILITY_FROM_QW1932.md"


def rel_dist(a: float, b: float, ref: float) -> float:
    return float(abs(a - b) / max(abs(ref), 1e-15))


def main() -> None:
    d1932 = json.loads((ROOT / "report_qw1932_physical_reparameterization_eta_scan.json").read_text(encoding="utf-8"))
    d2039 = json.loads((ROOT / "report_qw2039_derivation_compatible_refrozen_kernel_gate.json").read_text(encoding="utf-8"))
    d2042 = json.loads((ROOT / "report_qw2042_eft_matching_naturalness_audit.json").read_text(encoding="utf-8"))

    rows: List[Dict[str, object]] = d1932.get("results", [])
    strict_rows = [r for r in rows if bool(r.get("comparison", {}).get("strict_all_pass", False))]
    if not strict_rows:
        raise RuntimeError("No strict_all_pass rows in QW-1932 report.")

    beta0 = 0.01
    eta0 = 1.0

    target_beta = float(d2039["selected_kernel"]["beta"])
    target_eta = float(d2039["selected_kernel"]["eta"])
    target_z_beta = float(target_beta / beta0)
    target_delta_eta = float(target_eta - eta0)

    vals = []
    for r in strict_rows:
        fit = r["fit"]
        beta = float(fit["beta"])
        eta = float(r["eta"])
        z_beta = float(beta / beta0)
        d_eta = float(eta - eta0)
        vals.append(
            {
                "eta": eta,
                "beta": beta,
                "z_beta": z_beta,
                "delta_eta": d_eta,
                "objective": float(fit["objective"]),
                "primary_corr_ratio": float(r["comparison"]["primary_corr_ratio"]),
                "primary_gain_ratio": float(r["comparison"]["primary_gain_ratio"]),
                "stress_corr_ratio": float(r["comparison"]["stress_corr_ratio"]),
                "stress_gain_ratio": float(r["comparison"]["stress_gain_ratio"]),
            }
        )

    z_arr = np.array([v["z_beta"] for v in vals], dtype=float)
    de_arr = np.array([v["delta_eta"] for v in vals], dtype=float)

    z_ci95 = [float(np.quantile(z_arr, 0.025)), float(np.quantile(z_arr, 0.975))]
    de_ci95 = [float(np.quantile(de_arr, 0.025)), float(np.quantile(de_arr, 0.975))]

    inside_z = bool(z_ci95[0] <= target_z_beta <= z_ci95[1])
    inside_de = bool(de_ci95[0] <= target_delta_eta <= de_ci95[1])

    # Joint nearest strict point in normalized coordinates.
    z_scale = max(float(np.std(z_arr)), 1e-12)
    de_scale = max(float(np.std(de_arr)), 1e-12)
    for v in vals:
        dz = (float(v["z_beta"]) - target_z_beta) / z_scale
        dd = (float(v["delta_eta"]) - target_delta_eta) / de_scale
        v["joint_norm_distance"] = float(math.sqrt(dz * dz + dd * dd))

    nearest = min(vals, key=lambda v: float(v["joint_norm_distance"]))

    flags = {
        "target_Zbeta_inside_strict_ci95": inside_z,
        "target_delta_eta_inside_strict_ci95": inside_de,
        "nearest_joint_distance_le_0p75": bool(float(nearest["joint_norm_distance"]) <= 0.75),
        "strict_pass_count_ge_4": bool(len(strict_rows) >= 4),
    }
    pass_count = int(sum(1 for v in flags.values() if v))
    total_flags = int(len(flags))

    if pass_count == total_flags:
        verdict = "MICRO_RENORMALIZATION_FEASIBILITY_PASS"
        readiness = "STRONG_RENORMALIZATION_MICRO_FEASIBLE_WITHIN_EXISTING_BRANCH"
    elif pass_count >= 3:
        verdict = "MICRO_RENORMALIZATION_FEASIBILITY_PARTIAL"
        readiness = "PARTIAL_FEASIBILITY"
    else:
        verdict = "MICRO_RENORMALIZATION_FEASIBILITY_FAIL"
        readiness = "MICRO_BRANCH_DOES_NOT_SUPPORT_REQUIRED_RENORMALIZATION"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source_reports": {
            "qw1932": "report_qw1932_physical_reparameterization_eta_scan.json",
            "qw2039": "report_qw2039_derivation_compatible_refrozen_kernel_gate.json",
            "qw2042": "report_qw2042_eft_matching_naturalness_audit.json",
        },
        "target_from_refrozen": {
            "beta": target_beta,
            "eta": target_eta,
            "Z_beta_target": target_z_beta,
            "delta_eta_target": target_delta_eta,
        },
        "strict_pass_envelope_qw1932": {
            "n_rows": int(len(strict_rows)),
            "Z_beta_ci95": z_ci95,
            "delta_eta_ci95": de_ci95,
            "Z_beta_range": [float(np.min(z_arr)), float(np.max(z_arr))],
            "delta_eta_range": [float(np.min(de_arr)), float(np.max(de_arr))],
        },
        "nearest_strict_point": nearest,
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "readiness": readiness,
        "required_next_step": (
            "DERIVE_THE_SAME_ZBETA_AND_DELTA_ETA_POINTWISE_FROM_MICROSTATE_DISTRIBUTIONS"
            if verdict != "MICRO_RENORMALIZATION_FEASIBILITY_PASS"
            else "LOCK_MICRODERIVATION_PROTOCOL_AND_RUN_UNIFIED_NO_RETUNE_GATE"
        ),
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2043: MICRO RENORMALIZATION FEASIBILITY FROM QW-1932",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- Readiness: **{readiness}**",
        f"- pass_count: {pass_count}/{total_flags}",
        "",
        "## Target (from QW-2039 refrozen)",
        f"- Z_beta_target: {target_z_beta:.6f}",
        f"- delta_eta_target: {target_delta_eta:.6f}",
        "",
        "## Strict Envelope (QW-1932 strict_all_pass)",
        f"- n_rows: {len(strict_rows)}",
        f"- Z_beta CI95: [{z_ci95[0]:.6f}, {z_ci95[1]:.6f}]",
        f"- delta_eta CI95: [{de_ci95[0]:.6f}, {de_ci95[1]:.6f}]",
        "",
        "## Nearest Strict Point",
        f"- eta: {nearest['eta']:.6f}",
        f"- beta: {nearest['beta']:.6f}",
        f"- Z_beta: {nearest['z_beta']:.6f}",
        f"- delta_eta: {nearest['delta_eta']:.6f}",
        f"- joint_norm_distance: {nearest['joint_norm_distance']:.6f}",
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

    print(f"[QW-2043] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2043] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2043] verdict={verdict} pass_count={pass_count}/{total_flags}")


if __name__ == "__main__":
    main()
