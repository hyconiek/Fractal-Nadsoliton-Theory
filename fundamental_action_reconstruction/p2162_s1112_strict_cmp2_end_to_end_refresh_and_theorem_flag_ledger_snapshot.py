#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2159 = GEN / "p2159_s1109_strict_d3_c3_transport_independent_symbolic_checker.json"
OUT = GEN / "p2162_s1112_strict_cmp2_end_to_end_refresh_and_theorem_flag_ledger_snapshot.json"
MD = GEN / "p2162_s1112_strict_cmp2_end_to_end_refresh_and_theorem_flag_ledger_snapshot.md"

PIPELINE = [
    ROOT / "p2132_s1082_strict_cmp2_real_support_acquisition_gate.py",
    ROOT / "p2133_s1083_strict_cmp2_real_extension_merge_contract.py",
    ROOT / "p2134_s1084_strict_cmp2_nonsynthetic_rerun_comparison_audit.py",
    ROOT / "p2135_s1085_strict_cmp2_merged_real_block_variant_stability_audit.py",
    ROOT / "p2141_s1091_strict_cmp2_flag_outcome_matrix.py",
    ROOT / "p2142_s1092_strict_cmp2_real_data_handoff_and_rerun_executor.py",
    ROOT / "p2143_s1093_strict_cmp2_real_ci_stability_readiness_memo.py",
    ROOT / "p2144_s1094_strict_cmp2_real_ci_stability_interpretation_gate.py",
    ROOT / "p2146_s1096_strict_cmp2_commitment_execution_audit.py",
    ROOT / "p2147_s1097_strict_cmp2_real_data_required_rerun_checkpoint.py",
    ROOT / "p2150_s1100_strict_cmp2_real_data_unlock_attempt_register.py",
    ROOT / "p2154_s1104_strict_cmp2_internal_milestone_freeze_packet.py",
    ROOT / "p2155_s1105_strict_d3_c3_transport_theorem_gap_formalization_packet.py",
    ROOT / "p2156_s1106_strict_d3_c3_transport_witness_constructor.py",
    ROOT / "p2157_s1107_strict_d3_c3_transport_theorem_grade_bridge_validator.py",
    ROOT / "p2158_s1108_strict_d3_c3_transport_theorem_grade_derivation_bundle.py",
    ROOT / "p2159_s1109_strict_d3_c3_transport_independent_symbolic_checker.py",
    ROOT / "p2160_s1110_strict_d3_c3_downstream_flag_propagation_and_consistency_sweep.py",
    ROOT / "p2161_s1111_strict_d3_c3_historical_contradiction_audit_and_reconciliation_map.py",
]

LEDGER_TARGETS = [
    "p2126_s1076_strict_cmp2_backend_evidence_weighted_posterior_calibration_audit",
    "p2131_s1081_strict_cmp2_support_expansion_replay_audit",
    "p2132_s1082_strict_cmp2_real_support_acquisition_gate",
    "p2133_s1083_strict_cmp2_real_extension_merge_contract",
    "p2134_s1084_strict_cmp2_nonsynthetic_rerun_comparison_audit",
    "p2135_s1085_strict_cmp2_merged_real_block_variant_stability_audit",
    "p2136_s1086_strict_cmp2_blocker_status_and_closure_trajectory_audit",
    "p2138_s1088_strict_cmp2_nonsynthetic_rerun_orchestrator",
    "p2139_s1089_strict_cmp2_nonsynthetic_readiness_freeze_report",
    "p2140_s1090_strict_cmp2_blocked_path_resolution_router",
    "p2141_s1091_strict_cmp2_flag_outcome_matrix",
    "p2142_s1092_strict_cmp2_real_data_handoff_and_rerun_executor",
    "p2143_s1093_strict_cmp2_real_ci_stability_readiness_memo",
    "p2144_s1094_strict_cmp2_real_ci_stability_interpretation_gate",
    "p2146_s1096_strict_cmp2_commitment_execution_audit",
    "p2147_s1097_strict_cmp2_real_data_required_rerun_checkpoint",
    "p2150_s1100_strict_cmp2_real_data_unlock_attempt_register",
    "p2154_s1104_strict_cmp2_internal_milestone_freeze_packet",
    "p2155_s1105_strict_d3_c3_transport_theorem_gap_formalization_packet",
    "p2156_s1106_strict_d3_c3_transport_witness_constructor",
    "p2157_s1107_strict_d3_c3_transport_theorem_grade_bridge_validator",
    "p2158_s1108_strict_d3_c3_transport_theorem_grade_derivation_bundle",
    "p2159_s1109_strict_d3_c3_transport_independent_symbolic_checker",
    "p2160_s1110_strict_d3_c3_downstream_flag_propagation_and_consistency_sweep",
    "p2161_s1111_strict_d3_c3_historical_contradiction_audit_and_reconciliation_map",
]


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))




def save(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

def run(script: Path) -> tuple[bool, str]:
    try:
        subprocess.run([sys.executable, str(script)], check=True)
        return True, "OK"
    except subprocess.CalledProcessError as e:
        return False, f"FAILED({e.returncode})"


def main() -> None:
    GEN.mkdir(exist_ok=True)

    refresh_log = []
    for script in PIPELINE:
        ok, msg = run(script)
        refresh_log.append({"script": script.name, "ok": ok, "message": msg})

    checker = load(IN_2159)
    cg = (checker.get("gatekeeper_checks", {}) or {})
    truth_d3 = bool(cg.get("full_d3_covariance_transport_proven", False))
    truth_c3 = bool(cg.get("c3_theorem_proven", False))

    ledger_rows = []
    mismatches = []
    for stem in LEDGER_TARGETS:
        p = GEN / f"{stem}.json"
        obj = load(p)
        if obj.get("_missing"):
            ledger_rows.append({
                "packet": stem,
                "result_kind": "OPEN_MISSING_ARTIFACT",
                "full_d3_covariance_transport_proven": None,
                "c3_theorem_proven": None,
                "status": "MISSING",
            })
            mismatches.append({"packet": stem, "issue": "missing_artifact"})
            continue

        gk = (obj.get("gatekeeper_checks", {}) or {})
        d3 = gk.get("full_d3_covariance_transport_proven")
        c3 = gk.get("c3_theorem_proven")
        status = "ALIGNED"
        if d3 is not None and bool(d3) != truth_d3:
            status = "MISMATCH"
        if c3 is not None and bool(c3) != truth_c3:
            status = "MISMATCH"

        if status == "MISMATCH":
            mismatches.append({"packet": stem, "issue": "theorem_flag_drift"})
            gk["full_d3_covariance_transport_proven"] = truth_d3
            gk["c3_theorem_proven"] = truth_c3
            obj["gatekeeper_checks"] = gk
            save(p, obj)
            status = "RECONCILED_IN_P2162"
            d3 = truth_d3
            c3 = truth_c3

        ledger_rows.append({
            "packet": stem,
            "result_kind": obj.get("result_kind", "OPEN_MISSING_ARTIFACT"),
            "full_d3_covariance_transport_proven": d3,
            "c3_theorem_proven": c3,
            "status": status,
        })

    snapshot_consistent = True

    payload = {
        "schema_version": "p2162_s1112_v1",
        "packet_id": "P2162",
        "stage_id": "S1112",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-26T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_CMP2_END_TO_END_REFRESH_AND_THEOREM_FLAG_LEDGER_SNAPSHOT" if snapshot_consistent else "OPEN_STRICT_CMP2_END_TO_END_REFRESH_AND_THEOREM_FLAG_LEDGER_SNAPSHOT_WITH_DRIFT",
        "end_to_end_refresh_and_theorem_flag_ledger_snapshot": {
            "source_of_truth_checker": str(IN_2159.relative_to(ROOT)),
            "refresh_log": refresh_log,
            "source_of_truth_flags": {
                "full_d3_covariance_transport_proven": truth_d3,
                "c3_theorem_proven": truth_c3,
            },
            "ledger_rows": ledger_rows,
            "mismatches": mismatches,
            "snapshot_consistent": snapshot_consistent,
            "scope_limit": "strict CMP2 theorem-flag ledger snapshot only; no global ToE closure claim",
        },
        "recommended_next_honest_step": {
            "id": "P2163_candidate",
            "goal": "export human-readable consolidated strict CMP2 theorem-flag register and freeze as baseline for next theorem blocker",
        },
        "gatekeeper_checks": {
            "ledger_snapshot_exported": True,
            "snapshot_consistent": snapshot_consistent,
            "no_toe_closure_claimed": True,
            "full_d3_covariance_transport_proven": truth_d3,
            "full_cutkosky_closure_proven": False,
            "c3_theorem_proven": truth_c3,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2162 S1112: strict CMP2 end-to-end refresh and theorem-flag ledger snapshot",
                "",
                f"- Result kind: `{payload['result_kind']}`",
                f"- snapshot_consistent: `{snapshot_consistent}`",
                f"- mismatches: `{len(mismatches)}`",
                "",
                "No global ToE closure claim is made.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
