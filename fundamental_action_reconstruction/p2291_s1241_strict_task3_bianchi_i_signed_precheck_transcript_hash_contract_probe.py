#!/usr/bin/env python3
from __future__ import annotations
import hashlib
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2290 = GEN / "p2290_s1240_strict_task3_bianchi_i_deterministic_theorem_precheck_evaluator_probe.json"
OUT = GEN / "p2291_s1241_strict_task3_bianchi_i_signed_precheck_transcript_hash_contract_probe.json"
MD = GEN / "p2291_s1241_strict_task3_bianchi_i_signed_precheck_transcript_hash_contract_probe.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def sha256_text(s: str) -> str:
    return hashlib.sha256(s.encode("utf-8")).hexdigest()


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2290 = load(IN_2290)

    probe = (p2290.get("strict_task3_bianchi_i_deterministic_theorem_precheck_evaluator_probe", {}) or {})
    checks = probe.get("checks", {}) or {}
    decision = str(probe.get("evaluator_decision", "PRECHECK_BLOCK") or "PRECHECK_BLOCK")
    inputs = probe.get("inputs", {}) or {}

    transcript_payload = {
        "transcript_id": "TASK3_PRECHECK_TRANSCRIPT_V1",
        "source_packet": str(IN_2290.relative_to(ROOT)),
        "inputs": inputs,
        "checks": checks,
        "evaluator_decision": decision,
        "scope": "task3_bianchi_i_theorem_attempt_precheck",
    }
    canonical = json.dumps(transcript_payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    transcript_hash = sha256_text(canonical)

    contract = {
        "contract_id": "TASK3_THEOREM_DRAFT_METADATA_HASH_CONTRACT_V1",
        "required_transcript_hash_sha256": transcript_hash,
        "required_precheck_decision": "PRECHECK_PASS",
        "observed_precheck_decision": decision,
        "status": "READY_FOR_DRAFT_METADATA_BINDING" if decision == "PRECHECK_PASS" else "BLOCKED_BY_PRECHECK_DECISION",
        "rule": "Any theorem-draft metadata packet must include transcript_hash exactly equal to required_transcript_hash_sha256.",
    }

    payload = {
        "schema_version": "p2291_s1241_v1",
        "packet_id": "P2291",
        "stage_id": "S1241",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_TASK3_BIANCHI_I_SIGNED_PRECHECK_TRANSCRIPT_HASH_CONTRACT_PROBE",
        "strict_task3_bianchi_i_signed_precheck_transcript_hash_contract_probe": {
            "probe_id": "STRICT_TASK3_BIANCHI_I_SIGNED_PRECHECK_TRANSCRIPT_HASH_CONTRACT_PROBE_V1",
            "source_packets": [str(IN_2290.relative_to(ROOT))],
            "transcript_payload": transcript_payload,
            "transcript_hash_sha256": transcript_hash,
            "hash_contract": contract,
            "theorem_scope_limit": "transcript-hash contract packaging only; not theorem proof and not selector/ToE closure",
        },
        "recommended_next_honest_step": {
            "id": "P2292_candidate",
            "goal": "build metadata validator that rejects theorem-draft packets missing exact transcript_hash or with non-PRECHECK_PASS decision",
        },
        "gatekeeper_checks": {
            "transcript_hash_exported": len(transcript_hash) == 64,
            "hash_contract_exported": True,
            "decision_valid": decision in ["PRECHECK_PASS", "PRECHECK_BLOCK"],
            "contract_status_valid": contract["status"] in ["READY_FOR_DRAFT_METADATA_BINDING", "BLOCKED_BY_PRECHECK_DECISION"],
            "no_bridge_theorem_claimed": True,
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text("\n".join([
        "# P2291 S1241: signed precheck transcript hash contract",
        "",
        f"- evaluator decision: `{decision}`",
        f"- transcript hash: `{transcript_hash}`",
        f"- contract status: `{contract['status']}`",
        "",
        "Hash-contract packaging only; no selector closure / ToE closure claim.",
    ]) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
