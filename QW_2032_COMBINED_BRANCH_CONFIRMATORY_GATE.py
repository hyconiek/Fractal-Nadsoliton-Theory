#!/usr/bin/env python3
"""
QW-2032: Combined-branch confirmatory gate.

Aggregates strict prerequisites into one status gate:
- QW-2030 combined Stage-C pass,
- QW-2015 true external package readiness,
- QW-2017 blind intervention pass,
- QW-2031 blind eta-kernel transfer pass.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2032_combined_branch_confirmatory_gate.json"
OUT_MD = ROOT / "RAPORT_QW2032_COMBINED_BRANCH_CONFIRMATORY_GATE.md"


def load_json(path: Path) -> Dict:
    if not path.exists():
        raise RuntimeError(f"Missing required artifact: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    d2030 = load_json(ROOT / "report_qw2030_final_stage_c_gate_combined_branch.json")
    d2015 = load_json(ROOT / "report_qw2015_true_external_beta_channel_v2_readiness_gate.json")
    d2017 = load_json(ROOT / "report_qw2017_v2_beta_observable_blind_external_intervention.json")
    d2031 = load_json(ROOT / "report_qw2031_v2_eta_triad_blind_external_validation.json")

    flags = {
        "stage_c_combined_pass": bool(d2030.get("verdict") == "FINAL_STAGE_C_GATE_COMBINED_BRANCH_PASS"),
        "methodology_guard_no_sector_retune": bool(
            all(bool(x) for x in d2030.get("methodology_guard", {}).values())
        ),
        "external_package_ready": bool(d2015.get("readiness") == "READY"),
        "beta_intervention_pass": bool(d2017.get("verdict") == "BETA_OBSERVABLE_BLIND_EXTERNAL_INTERVENTION_PASS"),
        "eta_external_transfer_pass": bool(
            d2031.get("verdict") in {
                "ETA_TRIAD_BLIND_EXTERNAL_VALIDATION_PASS_STRONG",
                "ETA_TRIAD_BLIND_EXTERNAL_VALIDATION_PASS_PRIMARY_ONLY",
            }
        ),
    }

    pass_count = int(sum(1 for x in flags.values() if x))
    total_flags = int(len(flags))

    if pass_count == total_flags and d2031.get("verdict") == "ETA_TRIAD_BLIND_EXTERNAL_VALIDATION_PASS_STRONG":
        verdict = "COMBINED_BRANCH_CONFIRMATORY_GATE_PASS_STRONG"
        readiness = "STAGE_C_PLUS_EXTERNAL_PRECONFIRMATORY_CLOSED"
        next_step = "RUN_TRULY_INDEPENDENT_MULTITEAM_CONFIRMATORY_PACKAGE"
    elif pass_count == total_flags:
        verdict = "COMBINED_BRANCH_CONFIRMATORY_GATE_PASS_PARTIAL"
        readiness = "STAGE_C_PLUS_EXTERNAL_PRECONFIRMATORY_PARTIAL"
        next_step = "UPGRADE_ETA_TRANSFER_TO_STRONG_AND_RUN_INDEPENDENT_PACKAGE"
    else:
        verdict = "COMBINED_BRANCH_CONFIRMATORY_GATE_FAIL"
        readiness = "NOT_READY"
        next_step = "REPAIR_FAILED_PREREQUISITES"

    summary = {
        "mass_mean_rel_pct": float(d2030["combined_metrics"]["mass"]["mean_rel_err_pct"]),
        "mass_max_rel_pct": float(d2030["combined_metrics"]["mass"]["max_rel_err_pct"]),
        "ckm_mean_rel_pct": float(d2030["combined_metrics"]["flavor"]["ckm_mean_rel_pct"]),
        "pmns_mean_rel_pct": float(d2030["combined_metrics"]["flavor"]["pmns_mean_rel_pct"]),
        "gw_auc": float(d2030["combined_metrics"]["gw"]["auc_h1l1_vs_ctrl"]),
        "gw_adv": float(d2030["combined_metrics"]["gw"]["adv_shared_minus_ctrl_q90"]),
        "gw_sep": float(d2030["combined_metrics"]["gw"]["sep_median_h1l1_minus_ctrl"]),
        "gw_gap": float(d2030["combined_metrics"]["gw"]["control_median_gap"]),
        "ext_pairs_n": int(d2015["counts"]["n_pairs_total"]),
        "ext_events_n": int(d2015["counts"]["n_events_total"]),
        "eta_primary_corr": float(d2031["datasets"]["primary"]["holdout_metrics"]["pearson_corr"]),
        "eta_primary_gain": float(d2031["datasets"]["primary"]["holdout_metrics"]["rmse_gain_ratio"]),
        "eta_stress_corr": float(d2031["datasets"]["stress"]["holdout_metrics"]["pearson_corr"]),
        "eta_stress_gain": float(d2031["datasets"]["stress"]["holdout_metrics"]["rmse_gain_ratio"]),
    }

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "qw2030": "report_qw2030_final_stage_c_gate_combined_branch.json",
            "qw2015": "report_qw2015_true_external_beta_channel_v2_readiness_gate.json",
            "qw2017": "report_qw2017_v2_beta_observable_blind_external_intervention.json",
            "qw2031": "report_qw2031_v2_eta_triad_blind_external_validation.json",
        },
        "summary": summary,
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "readiness": readiness,
        "verdict": verdict,
        "required_next_step": next_step,
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2032: COMBINED BRANCH CONFIRMATORY GATE",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Readiness: **{readiness}**",
        f"- Verdict: **{verdict}**",
        f"- pass_count: {pass_count}/{total_flags}",
        "",
        "## Stage-C Combined Metrics (QW-2030)",
        f"- mass mean/max rel%: {summary['mass_mean_rel_pct']:.3f}/{summary['mass_max_rel_pct']:.3f}",
        f"- flavor CKM/PMNS mean rel%: {summary['ckm_mean_rel_pct']:.3f}/{summary['pmns_mean_rel_pct']:.3f}",
        f"- GW auc/adv/sep/gap: {summary['gw_auc']:.4f}/{summary['gw_adv']:.4f}/{summary['gw_sep']:.6f}/{summary['gw_gap']:.6f}",
        "",
        "## External Package (QW-2015)",
        f"- n_pairs/n_events: {summary['ext_pairs_n']}/{summary['ext_events_n']}",
        "",
        "## Eta Transfer (QW-2031)",
        f"- primary corr/gain: {summary['eta_primary_corr']:.4f}/{summary['eta_primary_gain']:.6f}",
        f"- stress corr/gain: {summary['eta_stress_corr']:.4f}/{summary['eta_stress_gain']:.6f}",
        "",
        "## Flags",
    ]
    for k, v in flags.items():
        lines.append(f"- {k}: {v}")

    lines.extend(
        [
            "",
            "## Required Next Step",
            f"- {next_step}",
            "",
            "## Artifacts",
            f"- JSON: `{OUT_JSON.name}`",
        ]
    )

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2032] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2032] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2032] readiness={readiness} verdict={verdict} pass={pass_count}/{total_flags}")


if __name__ == "__main__":
    main()
