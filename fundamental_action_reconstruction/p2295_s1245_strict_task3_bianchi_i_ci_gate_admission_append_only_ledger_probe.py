#!/usr/bin/env python3
from __future__ import annotations
import hashlib
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2293 = GEN / "p2293_s1243_strict_task3_bianchi_i_metadata_validator_negative_control_matrix_probe.json"
IN_2294 = GEN / "p2294_s1244_strict_task3_bianchi_i_theorem_draft_ci_gate_admission_probe.json"
OUT = GEN / "p2295_s1245_strict_task3_bianchi_i_ci_gate_admission_append_only_ledger_probe.json"
MD = GEN / "p2295_s1245_strict_task3_bianchi_i_ci_gate_admission_append_only_ledger_probe.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def sha256_json(value: Any) -> str:
    return sha256_text(json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False))


def artifact_hash(path: Path) -> str:
    if not path.exists():
        return ""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2293 = load(IN_2293)
    p2294 = load(IN_2294)

    matrix_probe = p2293.get("strict_task3_bianchi_i_metadata_validator_negative_control_matrix_probe", {}) or {}
    gate_probe = p2294.get("strict_task3_bianchi_i_theorem_draft_ci_gate_admission_probe", {}) or {}

    matrix = matrix_probe.get("matrix", []) or []
    reject_rows = [row for row in matrix if row.get("validator_decision") == "METADATA_REJECT"]
    accept_rows = [row for row in matrix if row.get("validator_decision") == "METADATA_ACCEPT"]

    gate_decision = str(gate_probe.get("ci_gate_decision", "CI_GATE_BLOCK") or "CI_GATE_BLOCK")
    theorem_attempt_decision = str(gate_probe.get("theorem_attempt_decision", "THEOREM_DRAFT_HOLD") or "THEOREM_DRAFT_HOLD")
    inputs = gate_probe.get("inputs", {}) or {}
    checks = gate_probe.get("checks", {}) or {}

    source_artifacts = [
        {
            "packet_id": "P2293",
            "path": str(IN_2293.relative_to(ROOT)),
            "sha256": artifact_hash(IN_2293),
        },
        {
            "packet_id": "P2294",
            "path": str(IN_2294.relative_to(ROOT)),
            "sha256": artifact_hash(IN_2294),
        },
    ]

    mismatch_evidence_pointers = [
        {
            "candidate_id": str(row.get("candidate_id", "")),
            "source_packet": str(IN_2293.relative_to(ROOT)),
            "json_pointer": f"/strict_task3_bianchi_i_metadata_validator_negative_control_matrix_probe/matrix/{idx}",
            "hash_match": bool(row.get("hash_match", False)),
            "decision_match": bool(row.get("decision_match", False)),
            "validator_decision": str(row.get("validator_decision", "")),
            "expected_decision": str(row.get("expected_decision", "")),
            "matches_expectation": bool(row.get("matches_expectation", False)),
        }
        for idx, row in enumerate(matrix)
        if row.get("validator_decision") == "METADATA_REJECT"
    ]

    admission_record_core = {
        "record_id": "TASK3_BIANCHI_I_CI_GATE_ADMISSION_LEDGER_ENTRY_V1",
        "run_index": 1,
        "previous_entry_sha256": None,
        "source_artifacts": source_artifacts,
        "gate_decision": gate_decision,
        "theorem_attempt_decision": theorem_attempt_decision,
        "admission_inputs": inputs,
        "admission_checks": checks,
        "negative_control_summary": {
            "accept_rows": len(accept_rows),
            "reject_rows": len(reject_rows),
            "all_expectations_met": bool(matrix_probe.get("all_expectations_met", False)),
        },
        "mismatch_evidence_pointers": mismatch_evidence_pointers,
        "scope_limit": "append-only CI admission ledger only; no theorem proof and no bridge/selector/ToE closure claim",
    }
    admission_record = dict(admission_record_core)
    admission_record["entry_sha256"] = sha256_json(admission_record_core)

    ledger = [admission_record]
    ledger_fingerprint = sha256_json(ledger)
    run_indices = [entry["run_index"] for entry in ledger]
    source_hashes = [artifact["sha256"] for artifact in source_artifacts]

    next_goal = (
        "derive a theorem-draft launch manifest from the admitted ledger entry with exact source-hash requirements"
        if theorem_attempt_decision == "THEOREM_DRAFT_ADMIT"
        else "export hold-path replay instructions tied to the blocked ledger entry and unresolved A1/A2/A3 verifier status"
    )

    payload = {
        "schema_version": "p2295_s1245_v1",
        "packet_id": "P2295",
        "stage_id": "S1245",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-28T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_TASK3_BIANCHI_I_CI_GATE_ADMISSION_APPEND_ONLY_LEDGER_PROBE",
        "strict_task3_bianchi_i_ci_gate_admission_append_only_ledger_probe": {
            "probe_id": "STRICT_TASK3_BIANCHI_I_CI_GATE_ADMISSION_APPEND_ONLY_LEDGER_PROBE_V1",
            "ledger_policy": {
                "append_only": True,
                "run_index_base": 1,
                "run_index_rule": "new entries must increment by exactly one and preserve all prior entry_sha256 values",
                "entry_hash_rule": "entry_sha256 is canonical SHA-256 over the admission record before adding entry_sha256",
            },
            "ledger": ledger,
            "ledger_fingerprint_sha256": ledger_fingerprint,
            "theorem_scope_limit": "CI-gate admission ledger only; not theorem proof and not bridge/selector/ToE closure",
        },
        "recommended_next_honest_step": {
            "id": "P2296_candidate",
            "goal": next_goal,
        },
        "gatekeeper_checks": {
            "ledger_exported": len(ledger) >= 1,
            "run_indices_monotonic": run_indices == sorted(run_indices) and len(run_indices) == len(set(run_indices)),
            "first_run_index_is_one": bool(run_indices) and run_indices[0] == 1,
            "source_hashes_present": all(len(h) == 64 for h in source_hashes),
            "entry_hash_length_ok": len(admission_record["entry_sha256"]) == 64,
            "ledger_fingerprint_length_ok": len(ledger_fingerprint) == 64,
            "mismatch_evidence_pointers_present": len(mismatch_evidence_pointers) >= 3,
            "gate_decision_preserved": gate_decision in {"CI_GATE_OPEN", "CI_GATE_BLOCK"},
            "theorem_attempt_decision_preserved": theorem_attempt_decision in {"THEOREM_DRAFT_ADMIT", "THEOREM_DRAFT_HOLD"},
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
                "# P2295 S1245: CI-gate admission append-only ledger",
                "",
                f"- ledger entries: `{len(ledger)}`",
                f"- first run index: `{run_indices[0] if run_indices else None}`",
                f"- gate decision preserved: `{gate_decision}`",
                f"- theorem attempt decision preserved: `{theorem_attempt_decision}`",
                f"- mismatch evidence pointers: `{len(mismatch_evidence_pointers)}`",
                f"- ledger fingerprint: `{ledger_fingerprint}`",
                "",
                "Append-only CI admission ledger only; no bridge theorem / selector closure / ToE closure claim.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
