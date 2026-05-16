#!/usr/bin/env python3
"""P1821 S771 strict priority lane sequencer checkpoint."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"


def load(name: str) -> dict:
    return json.loads((GEN / name).read_text(encoding="utf-8"))


def main() -> None:
    p1820 = load("p1820_s770_strict_current_priority_bottleneck_to_execution_contract_checkpoint.json")
    p1819 = load("p1819_s769_strict_pre_tg_verdict_statevector_lock_binding_checkpoint.json")

    snap = p1820.get("current_priority_snapshot", {})

    lane = [
        {
            "step": "S1_EXPORT_EA_EH_ELG_BOUNDARY_H1_BW_INPUTS",
            "targets": [
                "explicit_covariant_nonproxy_E_A_mu",
                "explicit_covariant_nonproxy_E_H",
                "metric_EL_g_export",
                "boundary_term_control",
                "H1_4D_weak_form_readiness",
                "Bianchi_Ward_consistency",
            ],
            "status_vector": {k: snap.get(k, "OPEN_UNKNOWN") for k in [
                "explicit_covariant_nonproxy_E_A_mu",
                "explicit_covariant_nonproxy_E_H",
                "metric_EL_g_export",
                "boundary_term_control",
                "H1_4D_weak_form_readiness",
                "Bianchi_Ward_consistency",
            ]},
            "completion_rule": "all targets must have residual/witness traces attached",
        },
        {
            "step": "S2_EXECUTE_TG1_UNIFIED_NONPROXY_RESIDUAL_RUN",
            "precondition": "P1819 status == LOCK_CONSISTENT",
            "precondition_observed": p1819.get("status", "OPEN_UNKNOWN"),
            "required_output": ["binary TG1 verdict", "full residual trace", "no-autopromotion state sync"],
            "completion_rule": "TG1 verdict exported with witness trace",
        },
        {
            "step": "S3_GATE_FORWARD_TO_TG2_TG3_ONLY_IF_TG1_PASS_ZERO",
            "targets": ["TG2_BRST", "TG3_CUT"],
            "status_vector": {
                "TG1_BW": snap.get("TG1_BW", "OPEN_UNKNOWN"),
                "TG2_BRST": snap.get("TG2_BRST", "OPEN_UNKNOWN"),
                "TG3_CUT": snap.get("TG3_CUT", "OPEN_UNKNOWN"),
            },
            "completion_rule": "maintain OPEN/LOCKED unless strict witness closes predecessor gate",
        },
    ]

    out = {
        "packet_id": "P1821",
        "stage_id": "S771",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "execution_lane": lane,
        "technical_progress": "Priority obligations converted into an ordered execution lane with strict completion rules.",
        "proven": "Current pipeline ordering is explicit: exports -> TG1 witness -> conditional TG2/TG3 progression.",
        "open": "TG1 residual witness and theorem-level gate closures remain open.",
        "false_pass_risk": "Skipping lane order could reintroduce passive gate promotion without witness trace.",
        "next_honest_step": "Run S1 evidence packaging, then execute S2 TG1 runpack with full trace.",
        "lay_explanation": "Mamy teraz kolejkę kroków: najpierw komplet dowodów, potem test TG1, a dopiero na końcu kolejne bramki.",
    }

    path = GEN / "p1821_s771_strict_priority_lane_sequencer_checkpoint.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
