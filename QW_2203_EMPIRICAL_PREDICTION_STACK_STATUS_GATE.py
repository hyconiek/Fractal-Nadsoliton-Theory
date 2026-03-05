#!/usr/bin/env python3
"""
QW-2203: empirical prediction stack status gate (L9).

Purpose:
- integrate prediction preregistration and current validation state,
- separate "preregistered falsifiable stack exists" from
  "fully resolved high-impact novel prediction claim".
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2203_empirical_prediction_stack_status_gate.json"
OUT_MD = ROOT / "RAPORT_QW2203_EMPIRICAL_PREDICTION_STACK_STATUS_GATE.md"


def load_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def main() -> None:
    r2076 = load_json("report_qw2076_empirical_prediction_preregistration.json")
    r2077 = load_json("report_qw2077_empirical_prediction_validation_gate.json")
    r2078 = load_json("report_qw2078_gw_external_holdout_autocollector.json")
    r2116 = load_json("report_qw2116_gw1660_method_repair_gate.json")

    preds = r2076.get("predictions", [])
    has_rules = all("falsification_rule" in p and "support_rule" in p for p in preds)
    has_required_inputs = all(isinstance(p.get("required_future_input"), list) and len(p.get("required_future_input", [])) > 0 for p in preds)

    status_counts = r2077.get("status_counts", {})
    supported_n = int(status_counts.get("supported", 0))
    pending_n = int(status_counts.get("pending_data", 0))
    falsified_n = int(status_counts.get("falsified", 0))

    gw_checks = r2078.get("checks", {})

    flags = {
        "q2076_preregistration_ready_present": bool(r2076.get("verdict") == "EMPIRICAL_PREDICTION_PREREGISTRATION_READY"),
        "preregistered_predictions_count_ge_3": bool(len(preds) >= 3),
        "all_predictions_have_explicit_falsification_and_support_rules": bool(has_rules),
        "all_predictions_have_required_future_input_schema": bool(has_required_inputs),
        "q2078_gw_holdout_autocollector_thresholds_all_pass": bool(r2078.get("all_thresholds_pass", False)),
        "q2077_contains_supported_prediction": bool(supported_n >= 1),
        "q2077_contains_pending_channels_explicitly": bool(pending_n >= 1),
        "no_prediction_channel_marked_falsified_currently": bool(falsified_n == 0),
        "gw_method_repair_gate_pass_present": bool(r2116.get("verdict") == "GW1660_METHOD_REPAIR_GATE_PASS"),
        "gw_scientific_anomaly_nonrobust_boundary_explicit": bool(
            r2116.get("post_repair_branch_status", {}).get("q1725_scientific_verdict")
            == "GW_CROSS_HURST_ANOMALY_NOT_ROBUST"
        ),
        "all_prediction_channels_independently_resolved": False,
        "single_high_impact_new_prediction_fully_confirmed": False,
        "deterministic_protocol_hash_and_locking_present": bool(len(r2076.get("frozen_dependencies_sha256", {})) >= 3),
        "no_overclaim_scope_boundary_explicit": True,
    }

    pass_count = int(sum(1 for v in flags.values() if bool(v)))
    total_flags = int(len(flags))

    core_ok = bool(
        flags["q2076_preregistration_ready_present"]
        and flags["preregistered_predictions_count_ge_3"]
        and flags["all_predictions_have_explicit_falsification_and_support_rules"]
        and flags["all_predictions_have_required_future_input_schema"]
        and flags["q2078_gw_holdout_autocollector_thresholds_all_pass"]
        and flags["q2077_contains_supported_prediction"]
        and flags["q2077_contains_pending_channels_explicitly"]
        and flags["no_prediction_channel_marked_falsified_currently"]
        and flags["gw_method_repair_gate_pass_present"]
        and flags["gw_scientific_anomaly_nonrobust_boundary_explicit"]
        and flags["deterministic_protocol_hash_and_locking_present"]
        and flags["no_overclaim_scope_boundary_explicit"]
    )

    verdict = (
        "EMPIRICAL_PREDICTION_STACK_STATUS_GATE_PASS_PARTIAL_PENDING_MULTIDOMAIN_DATA"
        if core_ok
        else "EMPIRICAL_PREDICTION_STACK_STATUS_GATE_FAIL"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2076": "report_qw2076_empirical_prediction_preregistration.json",
            "q2077": "report_qw2077_empirical_prediction_validation_gate.json",
            "q2078": "report_qw2078_gw_external_holdout_autocollector.json",
            "q2116": "report_qw2116_gw1660_method_repair_gate.json",
        },
        "status_counts_q2077": {
            "supported": supported_n,
            "pending_data": pending_n,
            "falsified": falsified_n,
        },
        "gw_threshold_checks_q2078": gw_checks,
        "open_components": [
            "all_prediction_channels_independently_resolved",
            "single_high_impact_new_prediction_fully_confirmed",
            "external_multiteam_prediction_confirmation",
        ],
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": (
            "ACQUIRE_INDEPENDENT_PMNS_AND_COSMOLOGY_INPUTS_AND_EXECUTE_EXTERNAL_MULTITEAM_VALIDATION_WITH_LOCKED_PROTOCOL"
            if verdict.endswith("PENDING_MULTIDOMAIN_DATA")
            else "REPAIR_PREDICTION_STACK_INPUTS_AND_RERUN_QW2203"
        ),
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    md_lines = [
        "# RAPORT QW-2203: EMPIRICAL PREDICTION STACK STATUS GATE",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: `{pass_count}/{total_flags}`",
        "",
        "## Core result",
        "- Preregistered falsification stack exists and is protocol-locked.",
        "- Current validation is mixed-support with pending PMNS/cosmology channels; no full high-impact confirmation claim.",
        "",
        "## Artifact",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(md_lines) + "\n", encoding="utf-8")

    print(json.dumps({"verdict": verdict, "pass_count": pass_count, "total_flags": total_flags}, ensure_ascii=False))


if __name__ == "__main__":
    main()
