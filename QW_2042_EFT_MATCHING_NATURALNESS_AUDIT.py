#!/usr/bin/env python3
"""
QW-2042: EFT-style matching and naturalness audit for canonical -> refrozen kernel.

Purpose:
Estimate what renormalization/matching constants are required to map
pre-1700 canonical kernel semantics into QW-2039 refrozen kernel semantics,
and classify whether the required bridge looks perturbative or strongly
non-perturbative.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2042_eft_matching_naturalness_audit.json"
OUT_MD = ROOT / "RAPORT_QW2042_EFT_MATCHING_NATURALNESS_AUDIT.md"


def min_orders_for_factor(target_factor: float, per_order_factor: float = 2.0) -> int:
    if target_factor <= 1.0:
        return 0
    return int(math.ceil(math.log(target_factor) / math.log(per_order_factor)))


def main() -> None:
    d2039 = json.loads((ROOT / "report_qw2039_derivation_compatible_refrozen_kernel_gate.json").read_text(encoding="utf-8"))
    d2041 = json.loads((ROOT / "report_qw2041_canonical_refrozen_reparameterization_audit.json").read_text(encoding="utf-8"))

    canonical = {
        "omega": float(math.pi / 4.0),
        "phi": float(math.pi / 6.0),
        "beta": 0.01,
        "eta": 1.0,
    }
    refrozen = {
        "omega": float(d2039["selected_kernel"]["omega"]),
        "phi": float(d2039["selected_kernel"]["phi"]),
        "beta": float(d2039["selected_kernel"]["beta"]),
        "eta": float(d2039["selected_kernel"]["eta"]),
    }

    z_omega = float(refrozen["omega"] / max(canonical["omega"], 1e-15))
    delta_phi = float(refrozen["phi"] - canonical["phi"])
    z_beta = float(refrozen["beta"] / max(canonical["beta"], 1e-15))
    delta_eta = float(refrozen["eta"] - canonical["eta"])

    ln_z_omega = float(math.log(max(z_omega, 1e-15)))
    ln_z_beta = float(math.log(max(z_beta, 1e-15)))

    naturalness_flags = {
        "abs_lnZomega_le_1": bool(abs(ln_z_omega) <= 1.0),
        "abs_lnZbeta_le_1": bool(abs(ln_z_beta) <= 1.0),
        "abs_delta_eta_le_0p30": bool(abs(delta_eta) <= 0.30),
        "abs_delta_phi_le_pi_over_6": bool(abs(delta_phi) <= (math.pi / 6.0)),
    }

    n_orders_omega = min_orders_for_factor(max(z_omega, 1.0 / max(z_omega, 1e-15)), per_order_factor=2.0)
    n_orders_beta = min_orders_for_factor(max(z_beta, 1.0 / max(z_beta, 1e-15)), per_order_factor=2.0)

    bridge_class = {
        "beta_channel": "NONPERTURBATIVE_REQUIRED" if n_orders_beta >= 5 else "PERTURBATIVE_PLAUSIBLE",
        "omega_channel": "NONPERTURBATIVE_REQUIRED" if n_orders_omega >= 5 else "PERTURBATIVE_PLAUSIBLE",
        "eta_channel": "ANOMALOUS_DIMENSION_ACTIVE" if abs(delta_eta) >= 0.5 else "WEAK_ANOMALOUS_DIMENSION",
    }

    naturalness_pass_count = int(sum(1 for v in naturalness_flags.values() if v))
    naturalness_total = int(len(naturalness_flags))

    reparam_fail = bool(d2041.get("verdict") == "CANONICAL_REFROZEN_REPARAMETERIZATION_FAIL")
    if naturalness_pass_count == naturalness_total and not reparam_fail:
        verdict = "EFT_MATCHING_NATURAL"
        readiness = "CANONICAL_TO_REFROZEN_BRIDGE_SIMPLE"
    elif naturalness_pass_count >= 2:
        verdict = "EFT_MATCHING_TENSIONED"
        readiness = "BRIDGE_REQUIRES_EXTENDED_OPERATOR_BASIS"
    else:
        verdict = "EFT_MATCHING_STRONGLY_NONNATURAL"
        readiness = "MICRODERIVATION_OF_STRONG_RENORMALIZATION_REQUIRED"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "canonical_kernel": canonical,
        "refrozen_kernel": refrozen,
        "matching_constants": {
            "Z_omega": z_omega,
            "delta_phi": delta_phi,
            "Z_beta": z_beta,
            "delta_eta": delta_eta,
            "ln_Z_omega": ln_z_omega,
            "ln_Z_beta": ln_z_beta,
        },
        "naturalness_flags": naturalness_flags,
        "naturalness_pass_count": naturalness_pass_count,
        "naturalness_total": naturalness_total,
        "minimal_orders_if_factor_per_order_le_2": {
            "omega_channel": n_orders_omega,
            "beta_channel": n_orders_beta,
        },
        "bridge_classification": bridge_class,
        "dependency": {
            "qw2041_verdict": d2041.get("verdict"),
            "qw2041_readiness": d2041.get("readiness"),
        },
        "verdict": verdict,
        "readiness": readiness,
        "required_next_step": (
            "DERIVE_Z_BETA_AND_DELTA_ETA_FROM_NADSOLITON_MICRODYNAMICS_WITHOUT_SECTOR_RETUNE"
            if verdict != "EFT_MATCHING_NATURAL"
            else "FORMALIZE_MATCHING_LAYER_IN_MAIN_THEORY"
        ),
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2042: EFT MATCHING NATURALNESS AUDIT",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- Readiness: **{readiness}**",
        "",
        "## Matching Constants (canonical -> refrozen)",
        f"- Z_omega = {z_omega:.6f}",
        f"- delta_phi = {delta_phi:.6f}",
        f"- Z_beta = {z_beta:.6f}",
        f"- delta_eta = {delta_eta:.6f}",
        f"- ln(Z_omega) = {ln_z_omega:.6f}",
        f"- ln(Z_beta) = {ln_z_beta:.6f}",
        "",
        "## Naturalness Flags",
    ]
    for k, v in naturalness_flags.items():
        lines.append(f"- {k}: {v}")

    lines.extend(
        [
            "",
            "## Minimal Orders (if |factor_per_order| <= 2)",
            f"- omega channel: {n_orders_omega}",
            f"- beta channel: {n_orders_beta}",
            "",
            "## Bridge Classification",
            f"- beta channel: {bridge_class['beta_channel']}",
            f"- omega channel: {bridge_class['omega_channel']}",
            f"- eta channel: {bridge_class['eta_channel']}",
            "",
            "## Dependency",
            f"- QW-2041 verdict: {d2041.get('verdict')}",
            "",
            "## Required Next Step",
            f"- {out['required_next_step']}",
            "",
            "## Artifacts",
            f"- JSON: `{OUT_JSON.name}`",
        ]
    )

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2042] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2042] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2042] verdict={verdict} naturalness={naturalness_pass_count}/{naturalness_total}")


if __name__ == "__main__":
    main()
