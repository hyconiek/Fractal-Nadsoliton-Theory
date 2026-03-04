#!/usr/bin/env python3
"""
QW-2070: Full radiative program baseline for SM+GR closure package.

Purpose:
- define an explicit radiative-workflow map for SM+GR closure,
- execute implemented one-loop sectors (QED and QCD),
- expose unresolved radiative sectors as formal blockers.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2070_full_radiative_program_baseline.json"
OUT_MD = ROOT / "RAPORT_QW2070_FULL_RADIATIVE_PROGRAM_BASELINE.md"


def load_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def load_optional_json(name: str) -> Dict | None:
    path = ROOT / name
    if not path.exists():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def is_closure_resolved_entry(e: Dict) -> bool:
    lvl = str(e.get("strict_level", ""))
    st = str(e.get("status", ""))
    if lvl == "si_definition" and st == "definition_constant":
        return True
    if lvl == "strict_internal_gate" and st.startswith("derived") and st != "derived_strict_target_miss":
        within = e.get("within_tolerance")
        if isinstance(within, bool):
            return within
        return True
    return False


def active_qed_sum_q2(mu_gev: float) -> float:
    # Charged fermions; Nc * Q^2 contributes to the one-loop QED beta coefficient.
    fermions = [
        ("e", 0.000511, -1.0, 1.0),
        ("mu", 0.10566, -1.0, 1.0),
        ("tau", 1.77686, -1.0, 1.0),
        ("u", 0.00216, +2.0 / 3.0, 3.0),
        ("d", 0.00467, -1.0 / 3.0, 3.0),
        ("s", 0.093, -1.0 / 3.0, 3.0),
        ("c", 1.27, +2.0 / 3.0, 3.0),
        ("b", 4.18, -1.0 / 3.0, 3.0),
        ("t", 173.0, +2.0 / 3.0, 3.0),
    ]
    return float(sum(nc * (q * q) for _, m, q, nc in fermions if mu_gev >= m))


def run_alpha_qed_one_loop(mu0: float, alpha0: float, mu_target: float) -> float:
    if mu0 <= 0.0 or mu_target <= 0.0 or alpha0 <= 0.0:
        raise ValueError("Scales and alpha0 must be > 0.")
    if abs(mu_target - mu0) < 1e-15:
        return float(alpha0)

    fermion_thresholds = sorted([0.000511, 0.00216, 0.00467, 0.093, 0.10566, 1.27, 1.77686, 4.18, 173.0])
    forward = mu_target > mu0
    cuts = [x for x in fermion_thresholds if mu0 < x < mu_target] if forward else [x for x in fermion_thresholds if mu_target < x < mu0]
    cuts = sorted(cuts)
    if not forward:
        cuts = list(reversed(cuts))

    boundaries = [mu0] + cuts + [mu_target]
    inv_alpha = 1.0 / alpha0

    for i in range(len(boundaries) - 1):
        a = boundaries[i]
        b = boundaries[i + 1]
        mu_mid = math.sqrt(a * b)
        sum_q2 = active_qed_sum_q2(mu_mid)
        # 1/alpha(mu2) = 1/alpha(mu1) - (2/3pi) * sum(Q_f^2) * ln(mu2/mu1)
        inv_alpha = inv_alpha - (2.0 / (3.0 * math.pi)) * sum_q2 * math.log(b / a)

    if inv_alpha <= 0.0:
        raise ValueError("Unphysical QED running: inverse alpha <= 0.")
    return float(1.0 / inv_alpha)


def active_nf_qcd(mu_gev: float) -> int:
    quark_masses = [0.00216, 0.00467, 0.093, 1.27, 4.18, 173.0]
    nf = sum(1 for m in quark_masses if mu_gev >= m)
    return int(max(3, min(6, nf)))


def run_alpha_s_one_loop(mu0: float, alpha0: float, mu_target: float) -> float:
    if mu0 <= 0.0 or mu_target <= 0.0 or alpha0 <= 0.0:
        raise ValueError("Scales and alpha0 must be > 0.")
    if abs(mu_target - mu0) < 1e-15:
        return float(alpha0)

    thresholds = sorted([1.27, 4.18, 173.0])
    forward = mu_target > mu0
    cuts = [x for x in thresholds if mu0 < x < mu_target] if forward else [x for x in thresholds if mu_target < x < mu0]
    cuts = sorted(cuts)
    if not forward:
        cuts = list(reversed(cuts))

    boundaries = [mu0] + cuts + [mu_target]
    inv_alpha = 1.0 / alpha0

    for i in range(len(boundaries) - 1):
        a = boundaries[i]
        b = boundaries[i + 1]
        mu_mid = math.sqrt(a * b)
        nf = active_nf_qcd(mu_mid)
        beta0 = 11.0 - (2.0 / 3.0) * nf
        # 1/alpha_s(mu2) = 1/alpha_s(mu1) + beta0/(2pi) * ln(mu2/mu1)
        inv_alpha = inv_alpha + (beta0 / (2.0 * math.pi)) * math.log(b / a)

    if inv_alpha <= 0.0:
        raise ValueError("Unphysical QCD running: inverse alpha_s <= 0.")
    return float(1.0 / inv_alpha)


def main() -> None:
    r2068 = load_json("report_qw2068_sm_gr_parameter_registry.json")
    r2069 = load_json("report_qw2069_full_sm_gr_derivation_package.json")
    r2072 = load_optional_json("report_qw2072_ew_yukawa_flavor_radiative_baselines.json")
    r2073 = load_optional_json("report_qw2073_radiative_channels_closure_upgrade.json")

    registry = r2068["groups"]
    alpha_inv_ref = float(next(x for x in registry["gauge_and_electroweak"] if x["id"] == "alpha_em_inv_mz")["value"])
    alpha_s_ref = float(next(x for x in registry["gauge_and_electroweak"] if x["id"] == "alpha_s_mz")["value"])
    mz_ref = float(next(x for x in registry["gauge_and_electroweak"] if x["id"] == "m_z")["value"])

    alpha0 = 1.0 / alpha_inv_ref
    alpha_qed_runs = []
    for mu in [0.000511, 0.10566, 1.77686, mz_ref]:
        alpha_mu = run_alpha_qed_one_loop(mu0=0.000511, alpha0=alpha0, mu_target=mu)
        alpha_qed_runs.append(
            {
                "mu_gev": float(mu),
                "alpha_em": alpha_mu,
                "alpha_em_inv": float(1.0 / alpha_mu),
            }
        )

    alpha_s_runs = []
    for mu in [2.0, 10.0, mz_ref, 173.0, 1000.0]:
        alpha_mu = run_alpha_s_one_loop(mu0=mz_ref, alpha0=alpha_s_ref, mu_target=mu)
        alpha_s_runs.append(
            {
                "mu_gev": float(mu),
                "alpha_s": alpha_mu,
                "n_f_active": active_nf_qcd(mu),
            }
        )

    update_map = {}
    if r2072 is not None:
        for row in r2072.get("channel_updates", []):
            if isinstance(row, dict) and "id" in row:
                update_map[row["id"]] = row
    if r2073 is not None:
        for row in r2073.get("channel_updates", []):
            if isinstance(row, dict) and "id" in row:
                update_map[row["id"]] = row

    def channel_from_update_or_missing(channel_id: str, missing_note: str) -> Dict:
        upd = update_map.get(channel_id)
        if upd is None:
            return {
                "id": channel_id,
                "status": "missing",
                "strict_level": "not_derived",
                "closure_ready": False,
                "notes": missing_note,
            }
        return {
            "id": channel_id,
            "status": upd.get("status", "implemented_baseline_nonclosing"),
            "strict_level": upd.get("strict_level", "model_radiative_baseline"),
            "closure_ready": bool(upd.get("closure_ready", False)),
            "notes": upd.get("notes", "Updated from QW-2072 radiative baseline report."),
        }

    channels: List[Dict] = [
        {
            "id": "qed_running_one_loop",
            "status": "implemented_one_loop_baseline",
            "strict_level": "model_radiative_baseline",
            "closure_ready": True,
            "notes": "Piecewise one-loop QED running across charged-fermion thresholds.",
        },
        {
            "id": "qcd_running_one_loop",
            "status": "implemented_one_loop_baseline",
            "strict_level": "model_radiative_baseline",
            "closure_ready": True,
            "notes": "Piecewise one-loop QCD running with nf threshold switching.",
        },
        channel_from_update_or_missing(
            "electroweak_oblique_and_delta_r",
            "Full electroweak precision loop closure is not yet implemented in this chain.",
        ),
        channel_from_update_or_missing(
            "yukawa_running_and_threshold_matching",
            "Full Yukawa matrix RG transport and threshold matching absent.",
        ),
        channel_from_update_or_missing(
            "ckm_pmns_rge_transport",
            "No complete CKM/PMNS radiative transport chain in strict package.",
        ),
        channel_from_update_or_missing(
            "gr_eft_running",
            "No full EFT curvature/torsion running closure for GR constants in package.",
        ),
        channel_from_update_or_missing(
            "cosmological_backreaction_radiative_chain",
            "No full radiative+backreaction closure for Lambda/H0 in strict package.",
        ),
    ]

    implemented_count = sum(1 for c in channels if c["status"].startswith("implemented"))
    closure_ready_count = sum(1 for c in channels if bool(c.get("closure_ready", False)))
    missing_count = sum(1 for c in channels if c["status"] == "missing")

    derivation_missing_ids = [e["id"] for e in r2069["entries"] if e.get("status") == "missing"]
    derivation_unresolved_ids = [
        e["id"] for e in r2069["entries"] if not is_closure_resolved_entry(e)
    ]
    radiative_sensitive_ids = {
        "alpha_s_mz",
        "g_f",
        "m_w",
        "m_z",
        "m_h",
        "v_higgs",
        "m_up",
        "m_down",
        "m_strange",
        "m_nu1",
        "m_nu2",
        "m_nu3",
        "delta_cp_ckm",
        "delta_cp_pmns",
        "G_newton",
        "lambda_cosmological",
        "h0",
    }
    radiative_sensitive_missing_direct = sorted([x for x in derivation_missing_ids if x in radiative_sensitive_ids])
    radiative_sensitive_unresolved = sorted([x for x in derivation_unresolved_ids if x in radiative_sensitive_ids])

    if (
        missing_count == 0
        and closure_ready_count == len(channels)
        and len(radiative_sensitive_unresolved) == 0
    ):
        verdict = "FULL_RADIATIVE_PROGRAM_PASS"
    elif implemented_count >= 2:
        verdict = "FULL_RADIATIVE_PROGRAM_PARTIAL_BASELINE"
    else:
        verdict = "FULL_RADIATIVE_PROGRAM_FAIL"

    if missing_count > 0:
        required_next_step = "IMPLEMENT_MISSING_RADIATIVE_CHANNELS"
    elif closure_ready_count < len(channels):
        required_next_step = "UPGRADE_NONCLOSING_RADIATIVE_CHANNELS_TO_CLOSURE_READY"
    elif len(radiative_sensitive_unresolved) > 0:
        required_next_step = "DERIVE_REMAINING_RADIATIVE_SENSITIVE_SM_GR_PARAMETERS"
    else:
        required_next_step = "MERGE_WITH_FULL_PRECISION_CLOSURE_GATE"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "registry": "report_qw2068_sm_gr_parameter_registry.json",
            "full_derivation_package": "report_qw2069_full_sm_gr_derivation_package.json",
            "ew_yukawa_flavor_baselines": (
                "report_qw2072_ew_yukawa_flavor_radiative_baselines.json"
                if r2072 is not None
                else None
            ),
            "radiative_channels_closure_upgrade": (
                "report_qw2073_radiative_channels_closure_upgrade.json"
                if r2073 is not None
                else None
            ),
        },
        "channels": channels,
        "qed_one_loop_baseline": {
            "mu0_gev": 0.000511,
            "alpha0": alpha0,
            "samples": alpha_qed_runs,
        },
        "qcd_one_loop_baseline": {
            "mu0_gev": mz_ref,
            "alpha_s_mu0": alpha_s_ref,
            "samples": alpha_s_runs,
        },
        "coverage": {
            "n_channels_total": len(channels),
            "n_channels_implemented": implemented_count,
            "n_channels_closure_ready": closure_ready_count,
            "n_channels_missing": missing_count,
            "n_radiative_sensitive_missing_direct_parameters": len(radiative_sensitive_missing_direct),
            "radiative_sensitive_missing_direct_parameters": radiative_sensitive_missing_direct,
            "n_radiative_sensitive_unresolved_parameters": len(radiative_sensitive_unresolved),
            "radiative_sensitive_unresolved_parameters": radiative_sensitive_unresolved,
        },
        "verdict": verdict,
        "required_next_step": required_next_step,
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2070: FULL RADIATIVE PROGRAM BASELINE",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        "",
        "## Channel Coverage",
        f"- total channels: {len(channels)}",
        f"- implemented: {implemented_count}",
        f"- closure-ready: {closure_ready_count}",
        f"- missing: {missing_count}",
        f"- radiative-sensitive missing-direct parameters from QW-2069: {len(radiative_sensitive_missing_direct)}",
        f"- radiative-sensitive strict-unresolved parameters from QW-2069: {len(radiative_sensitive_unresolved)}",
        "",
        "## Implemented Baselines",
        "- QED one-loop running: implemented",
        "- QCD one-loop running: implemented",
        "- EW/Yukawa/CKM-PMNS baseline blocks: imported from QW-2072 when available",
        "",
        "## QED Sample",
    ]
    for row in alpha_qed_runs:
        lines.append(
            f"- mu={row['mu_gev']:.6g} GeV -> alpha^-1={row['alpha_em_inv']:.6f}"
        )

    lines.extend(["", "## QCD Sample"])
    for row in alpha_s_runs:
        lines.append(
            f"- mu={row['mu_gev']:.6g} GeV -> alpha_s={row['alpha_s']:.6f} (n_f={row['n_f_active']})"
        )

    lines.extend(["", "## Artifact", f"- JSON: `{OUT_JSON.name}`"])
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2070] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2070] Saved MD:   {OUT_MD.name}")
    print(
        f"[QW-2070] verdict={verdict} implemented={implemented_count}/{len(channels)} "
        f"closure_ready={closure_ready_count}/{len(channels)} "
        f"radiative_unresolved={len(radiative_sensitive_unresolved)}"
    )


if __name__ == "__main__":
    main()
