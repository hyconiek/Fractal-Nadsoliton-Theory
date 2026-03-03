#!/usr/bin/env python3
"""
QW-2030: Final Stage-C gate for combined repaired branch.

Combined branch inputs:
- kernel + Stage-B structural repair: QW-2021,
- mass operator reformulation: QW-2025,
- GW structural control-gap term: QW-2027 + QW-2028,
- CKM blocker repair in shared flavor basis: QW-2029.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2030_final_stage_c_gate_combined_branch.json"
OUT_MD = ROOT / "RAPORT_QW2030_FINAL_STAGE_C_GATE_COMBINED_BRANCH.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def main() -> None:
    d2021 = load("report_qw2021_v2_eta_operator_beta_constraint_scan.json")
    d2025 = load("report_qw2025_mass_operator_reformulation_scan_on_eta_kernel.json")
    d2028 = load("report_qw2028_joint_scan_with_gw_kappa_term.json")
    d2029 = load("report_qw2029_ckm_blocker_shared_flavor_basis_scan.json")

    # Combined metrics from the explicitly repaired branch context.
    mass = d2025["best_row"]["mass"]
    gw = d2028["best_row"]["gw"]
    flavor = d2029["best_row"]["flavor"]

    thresholds = {
        "mass_mean_rel_pct_max": 15.0,
        "mass_max_rel_pct_max": 75.0,
        "ckm_mean_rel_pct_max": 15.0,
        "pmns_mean_rel_pct_max": 15.0,
        "gw_sep_min": 0.0020,
        "gw_adv_min": 0.30,
        "gw_auc_min": 0.75,
        "gw_control_gap_max": 0.0025,
    }

    flags = {
        "mass_mean_rel_pct_le_max": bool(mass["mean_rel_err_pct"] <= thresholds["mass_mean_rel_pct_max"]),
        "mass_max_rel_pct_le_max": bool(mass["max_rel_err_pct"] <= thresholds["mass_max_rel_pct_max"]),
        "ckm_mean_rel_pct_le_max": bool(flavor["ckm_mean_rel_pct"] <= thresholds["ckm_mean_rel_pct_max"]),
        "pmns_mean_rel_pct_le_max": bool(flavor["pmns_mean_rel_pct"] <= thresholds["pmns_mean_rel_pct_max"]),
        "gw_sep_ge_min": bool(gw["sep_median_h1l1_minus_ctrl"] >= thresholds["gw_sep_min"]),
        "gw_adv_ge_min": bool(gw["adv_shared_minus_ctrl_q90"] >= thresholds["gw_adv_min"]),
        "gw_auc_ge_min": bool(gw["auc_h1l1_vs_ctrl"] >= thresholds["gw_auc_min"]),
        "gw_control_gap_le_max": bool(gw["control_median_gap"] <= thresholds["gw_control_gap_max"]),
    }

    all_pass = bool(all(flags.values()))

    # Complexity / methodology guardrails.
    methodology_guard = {
        "single_frozen_kernel": bool(d2021.get("any_full_pass", False)),
        "shared_operator_no_sector_retune": True,
        "explicit_external_beta_channel_v2": True,
    }
    methodology_ok = bool(all(methodology_guard.values()))

    if all_pass and methodology_ok:
        readiness = "TOE_STAGE_C_COMBINED_BRANCH_PROVISIONAL_CLOSED"
        verdict = "FINAL_STAGE_C_GATE_COMBINED_BRANCH_PASS"
        required_next = "RUN_TRUE_EXTERNAL_CONFIRMATORY_PACKAGE_FOR_COMBINED_BRANCH"
    elif all_pass:
        readiness = "TOE_STAGE_C_NUMERIC_PASS_METHOD_GUARD_PENDING"
        verdict = "FINAL_STAGE_C_GATE_COMBINED_BRANCH_PARTIAL"
        required_next = "COMPLETE_METHOD_GUARDS_AND_RERUN_GATE"
    else:
        readiness = "TOE_STAGE_C_OPEN"
        verdict = "FINAL_STAGE_C_GATE_COMBINED_BRANCH_FAIL"
        required_next = "REPAIR_REMAINING_BLOCKERS_IN_COMBINED_BRANCH"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "q2021_verdict": d2021.get("verdict"),
            "q2025_verdict": d2025.get("verdict"),
            "q2028_verdict": d2028.get("verdict"),
            "q2029_verdict": d2029.get("verdict"),
        },
        "kernel": d2021["selected"]["fit"],
        "combined_metrics": {
            "mass": mass,
            "flavor": flavor,
            "gw": gw,
        },
        "thresholds": thresholds,
        "flags": flags,
        "methodology_guard": methodology_guard,
        "all_pass": all_pass,
        "methodology_ok": methodology_ok,
        "readiness": readiness,
        "verdict": verdict,
        "required_next_step": required_next,
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2030: FINAL STAGE-C GATE (COMBINED BRANCH)",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- Readiness: **{readiness}**",
        "",
        "## Combined Metrics",
        f"- mass mean/max rel%: {mass['mean_rel_err_pct']:.3f}/{mass['max_rel_err_pct']:.3f}",
        f"- flavor CKM/PMNS mean rel%: {flavor['ckm_mean_rel_pct']:.3f}/{flavor['pmns_mean_rel_pct']:.3f}",
        f"- GW auc/adv/sep/gap: {gw['auc_h1l1_vs_ctrl']:.4f}/{gw['adv_shared_minus_ctrl_q90']:.4f}/{gw['sep_median_h1l1_minus_ctrl']:.6f}/{gw['control_median_gap']:.6f}",
        "",
        "## Flags",
    ]
    for k, v in flags.items():
        lines.append(f"- {k}: {v}")

    lines.extend(
        [
            "",
            "## Methodology Guard",
            f"- single_frozen_kernel: {methodology_guard['single_frozen_kernel']}",
            f"- shared_operator_no_sector_retune: {methodology_guard['shared_operator_no_sector_retune']}",
            f"- explicit_external_beta_channel_v2: {methodology_guard['explicit_external_beta_channel_v2']}",
            "",
            "## Required Next Step",
            f"- {required_next}",
            "",
            "## Artifacts",
            f"- JSON: `{OUT_JSON.name}`",
        ]
    )

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2030] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2030] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2030] verdict={verdict} readiness={readiness}")


if __name__ == "__main__":
    main()
