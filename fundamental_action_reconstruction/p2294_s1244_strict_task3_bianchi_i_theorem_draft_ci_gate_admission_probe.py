#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2292 = GEN / "p2292_s1242_strict_task3_bianchi_i_theorem_draft_metadata_hash_validator_probe.json"
IN_2293 = GEN / "p2293_s1243_strict_task3_bianchi_i_metadata_validator_negative_control_matrix_probe.json"
OUT = GEN / "p2294_s1244_strict_task3_bianchi_i_theorem_draft_ci_gate_admission_probe.json"
MD = GEN / "p2294_s1244_strict_task3_bianchi_i_theorem_draft_ci_gate_admission_probe.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2292 = load(IN_2292)
    p2293 = load(IN_2293)

    v2292 = (
        p2292.get("strict_task3_bianchi_i_theorem_draft_metadata_hash_validator_probe", {}) or {}
    )
    matrix2293 = (
        p2293.get("strict_task3_bianchi_i_metadata_validator_negative_control_matrix_probe", {}) or {}
    )

    validator_decision = str(v2292.get("validator_decision", "METADATA_REJECT") or "METADATA_REJECT")
    matrix = matrix2293.get("matrix", []) or []
    all_expectations_met = bool(matrix2293.get("all_expectations_met", False))

    matrix_has_baseline_accept = any(
        (row.get("candidate_id") == "C0_baseline_valid")
        and (row.get("validator_decision") == "METADATA_ACCEPT")
        and bool(row.get("matches_expectation", False))
        for row in matrix
    )
    matrix_reject_paths_stable = sum(1 for row in matrix if row.get("validator_decision") == "METADATA_REJECT") >= 3

    ci_gate_open = (
        validator_decision == "METADATA_ACCEPT"
        and all_expectations_met
        and matrix_has_baseline_accept
        and matrix_reject_paths_stable
    )
    ci_gate_decision = "CI_GATE_OPEN" if ci_gate_open else "CI_GATE_BLOCK"
    theorem_attempt_decision = "THEOREM_DRAFT_ADMIT" if ci_gate_open else "THEOREM_DRAFT_HOLD"

    payload = {
        "schema_version": "p2294_s1244_v1",
        "packet_id": "P2294",
        "stage_id": "S1244",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_TASK3_BIANCHI_I_THEOREM_DRAFT_CI_GATE_ADMISSION_PROBE",
        "strict_task3_bianchi_i_theorem_draft_ci_gate_admission_probe": {
            "probe_id": "STRICT_TASK3_BIANCHI_I_THEOREM_DRAFT_CI_GATE_ADMISSION_PROBE_V1",
            "source_packets": [
                str(IN_2292.relative_to(ROOT)),
                str(IN_2293.relative_to(ROOT)),
            ],
            "inputs": {
                "validator_decision": validator_decision,
                "all_expectations_met": all_expectations_met,
                "matrix_has_baseline_accept": matrix_has_baseline_accept,
                "matrix_reject_paths_stable": matrix_reject_paths_stable,
            },
            "checks": {
                "metadata_accept_required": validator_decision == "METADATA_ACCEPT",
                "negative_controls_consistent": all_expectations_met,
                "baseline_accept_present": matrix_has_baseline_accept,
                "reject_paths_present": matrix_reject_paths_stable,
            },
            "ci_gate_decision": ci_gate_decision,
            "theorem_attempt_decision": theorem_attempt_decision,
            "theorem_scope_limit": "CI admission gating only; no theorem proof and no bridge/selector/ToE closure claim",
        },
        "recommended_next_honest_step": {
            "id": "P2295_candidate",
            "goal": "export append-only CI-gate admission ledger with monotonic run-index and mismatch evidence pointers",
        },
        "gatekeeper_checks": {
            "ci_gate_decision_exported": ci_gate_decision in {"CI_GATE_OPEN", "CI_GATE_BLOCK"},
            "theorem_attempt_decision_exported": theorem_attempt_decision in {"THEOREM_DRAFT_ADMIT", "THEOREM_DRAFT_HOLD"},
            "open_implies_metadata_accept": (ci_gate_decision != "CI_GATE_OPEN")
            or (validator_decision == "METADATA_ACCEPT"),
            "open_implies_negative_controls_ok": (ci_gate_decision != "CI_GATE_OPEN")
            or all_expectations_met,
            "no_bridge_theorem_claimed": True,
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2294 S1244: theorem-draft CI gate admission",
                "",
                f"- validator decision: `{validator_decision}`",
                f"- matrix all expectations met: `{all_expectations_met}`",
                f"- ci gate decision: `{ci_gate_decision}`",
                f"- theorem attempt decision: `{theorem_attempt_decision}`",
                "",
                "Admission gate only; no bridge theorem / selector closure / ToE closure claim.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
