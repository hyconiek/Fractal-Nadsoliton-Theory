#!/usr/bin/env python3
from __future__ import annotations
import json
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT_POLICY = ROOT / "generated" / "p1797_s747_strict_sv1_to_sv5_run_manifest_binding_and_gate_decision_checkpoint.json"
OUT_TEMPLATE = ROOT / "generated" / "p1797_s747_sv1_sv5_run_manifest_template.json"

policy = {
    "packet_id": "P1797_S747",
    "as_of": str(date(2026, 5, 15)),
    "bw_pass_requirements": [
        "all_sv_bindings_present",
        "common_freeze_id",
        "p1793_intake_pass",
        "bw_residual_explicit_zero_with_digest",
    ],
    "bw_fail_verdict": "OPEN_OBSTRUCTION_WITH_TRACE",
    "downstream_locks_if_bw_not_pass": {"G_BRST": "LOCKED", "G_CUT": "LOCKED"},
}

template = {
    "freeze_id_common": "FREEZE_ID_REQUIRED",
    "sv_bindings": {
        "SV1": "EA_covariant_nonproxy_artifact_id",
        "SV2": "EH_covariant_nonproxy_artifact_id",
        "SV3": "ELg_nonproxy_artifact_id",
        "SV4": "boundary_control_witness_id",
        "SV5": "H1_4D_weak_form_witness_id",
    },
    "validation": {
        "intake_validation_result_id": "p1793_result_id",
        "state_update_protocol_id": "p1794_result_id",
    },
    "gate_decision_inputs": {
        "bw_residual_report_id": "REPORT_ID",
        "bw_check_digest": "HASH_REQUIRED"
    }
}

for path, data in ((OUT_POLICY, policy), (OUT_TEMPLATE, template)):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
        f.write("\n")
    print(path)
