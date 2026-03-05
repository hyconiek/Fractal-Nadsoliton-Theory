#!/usr/bin/env python3
"""
QW-2204: external multiteam execution status gate (L10).

Purpose:
- integrate external replication readiness chain (bundle/rehearsal/governance/lock),
- keep explicit boundary between "ready packet" and "truly independent executed".
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2204_external_multiteam_execution_status_gate.json"
OUT_MD = ROOT / "RAPORT_QW2204_EXTERNAL_MULTITEAM_EXECUTION_STATUS_GATE.md"


def load_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def main() -> None:
    r2032 = load_json("report_qw2032_combined_branch_confirmatory_gate.json")
    r2033 = load_json("report_qw2033_independent_confirmatory_freeze_bundle.json")
    r2050 = load_json("report_qw2050_spectral_micro_bridge_freeze_bundle.json")
    r2051 = load_json("report_qw2051_independent_rehearsal_gate.json")
    r2052 = load_json("report_qw2052_external_source_only_governance_gate.json")
    r2053 = load_json("report_qw2053_independent_multiteam_protocol_lock.json")
    r2016 = load_json("report_qw2016_v2_triad_blind_external_validation.json")
    r2017 = load_json("report_qw2017_v2_beta_observable_blind_external_intervention.json")

    flags = {
        "q2032_combined_confirmatory_gate_pass_strong": bool(
            r2032.get("verdict") == "COMBINED_BRANCH_CONFIRMATORY_GATE_PASS_STRONG"
        ),
        "q2016_external_blind_validation_pass_strong": bool(
            r2016.get("verdict") == "TRIAD_BLIND_EXTERNAL_VALIDATION_PASS_STRONG"
        ),
        "q2017_external_intervention_pass": bool(
            r2017.get("verdict") == "BETA_OBSERVABLE_BLIND_EXTERNAL_INTERVENTION_PASS"
        ),
        "q2033_independent_freeze_bundle_ready": bool(
            r2033.get("verdict") == "INDEPENDENT_CONFIRMATORY_FREEZE_BUNDLE_READY"
        ),
        "q2050_spectral_bundle_ready": bool(
            r2050.get("verdict") == "SPECTRAL_MICRO_BRIDGE_FREEZE_BUNDLE_READY"
        ),
        "q2051_independent_rehearsal_pass": bool(
            r2051.get("verdict") == "INDEPENDENT_REHEARSAL_GATE_PASS"
        ),
        "q2052_external_source_only_governance_pass": bool(
            r2052.get("verdict") == "EXTERNAL_SOURCE_ONLY_GOVERNANCE_PASS"
        ),
        "q2053_independent_multiteam_protocol_lock_ready": bool(
            r2053.get("verdict") == "INDEPENDENT_MULTITEAM_PROTOCOL_LOCK_READY"
        ),
        "locked_package_handoff_instructions_present": bool(
            r2053["flags"].get("runbook_written", False) and bool(r2053.get("lock_file"))
        ),
        "hash_locked_protocol_manifest_present": bool(bool(r2053.get("lock_sha256"))),
        "truly_independent_multiteam_execution_completed": False,
        "at_least_two_external_teams_completed_and_reported": False,
        "independent_team_reports_public_and_signed": False,
        "no_overclaim_scope_boundary_explicit": True,
    }

    pass_count = int(sum(1 for v in flags.values() if bool(v)))
    total_flags = int(len(flags))

    core_ok = bool(
        flags["q2032_combined_confirmatory_gate_pass_strong"]
        and flags["q2016_external_blind_validation_pass_strong"]
        and flags["q2017_external_intervention_pass"]
        and flags["q2033_independent_freeze_bundle_ready"]
        and flags["q2050_spectral_bundle_ready"]
        and flags["q2051_independent_rehearsal_pass"]
        and flags["q2052_external_source_only_governance_pass"]
        and flags["q2053_independent_multiteam_protocol_lock_ready"]
        and flags["locked_package_handoff_instructions_present"]
        and flags["hash_locked_protocol_manifest_present"]
        and flags["no_overclaim_scope_boundary_explicit"]
    )

    verdict = (
        "EXTERNAL_MULTITEAM_EXECUTION_STATUS_GATE_PASS_PARTIAL_PACKET_READY_EXECUTION_PENDING"
        if core_ok
        else "EXTERNAL_MULTITEAM_EXECUTION_STATUS_GATE_FAIL"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2032": "report_qw2032_combined_branch_confirmatory_gate.json",
            "q2033": "report_qw2033_independent_confirmatory_freeze_bundle.json",
            "q2050": "report_qw2050_spectral_micro_bridge_freeze_bundle.json",
            "q2051": "report_qw2051_independent_rehearsal_gate.json",
            "q2052": "report_qw2052_external_source_only_governance_gate.json",
            "q2053": "report_qw2053_independent_multiteam_protocol_lock.json",
            "q2016": "report_qw2016_v2_triad_blind_external_validation.json",
            "q2017": "report_qw2017_v2_beta_observable_blind_external_intervention.json",
        },
        "open_components": [
            "truly_independent_multiteam_execution_completed",
            "at_least_two_external_teams_completed_and_reported",
            "independent_team_reports_public_and_signed",
        ],
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": (
            "EXECUTE_LOCKED_PROTOCOL_BY_TRULY_INDEPENDENT_MULTITEAM_AND_PUBLISH_SIGNED_REPORTS"
            if verdict.endswith("EXECUTION_PENDING")
            else "REPAIR_EXTERNAL_PACKET_CHAIN_AND_RERUN_QW2204"
        ),
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    md_lines = [
        "# RAPORT QW-2204: EXTERNAL MULTITEAM EXECUTION STATUS GATE",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: `{pass_count}/{total_flags}`",
        "",
        "## Core result",
        "- External replication packet chain is locked and ready for handoff.",
        "- Truly independent multiteam execution/public signed reports remain pending.",
        "",
        "## Artifact",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(md_lines) + "\n", encoding="utf-8")

    print(json.dumps({"verdict": verdict, "pass_count": pass_count, "total_flags": total_flags}, ensure_ascii=False))


if __name__ == "__main__":
    main()
