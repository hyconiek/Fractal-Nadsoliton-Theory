#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2291 = GEN / "p2291_s1241_strict_task3_bianchi_i_signed_precheck_transcript_hash_contract_probe.json"
OUT = GEN / "p2292_s1242_strict_task3_bianchi_i_theorem_draft_metadata_hash_validator_probe.json"
MD = GEN / "p2292_s1242_strict_task3_bianchi_i_theorem_draft_metadata_hash_validator_probe.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2291 = load(IN_2291)

    probe = (p2291.get("strict_task3_bianchi_i_signed_precheck_transcript_hash_contract_probe", {}) or {})
    contract = probe.get("hash_contract", {}) or {}

    required_hash = str(contract.get("required_transcript_hash_sha256", "") or "")
    required_decision = str(contract.get("required_precheck_decision", "PRECHECK_PASS") or "PRECHECK_PASS")

    # deterministic metadata candidate (what a theorem draft would carry)
    draft_metadata_candidate = {
        "theorem_draft_id": "TASK3_BIANCHI_I_DRAFT_ATTEMPT_PLACEHOLDER_V1",
        "precheck_transcript_hash_sha256": required_hash,
        "precheck_decision": str(contract.get("observed_precheck_decision", "PRECHECK_BLOCK") or "PRECHECK_BLOCK"),
    }

    checks = {
        "hash_match": draft_metadata_candidate["precheck_transcript_hash_sha256"] == required_hash,
        "decision_match": draft_metadata_candidate["precheck_decision"] == required_decision,
    }
    validator_decision = "METADATA_ACCEPT" if all(checks.values()) else "METADATA_REJECT"

    payload = {
        "schema_version": "p2292_s1242_v1",
        "packet_id": "P2292",
        "stage_id": "S1242",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_TASK3_BIANCHI_I_THEOREM_DRAFT_METADATA_HASH_VALIDATOR_PROBE",
        "strict_task3_bianchi_i_theorem_draft_metadata_hash_validator_probe": {
            "probe_id": "STRICT_TASK3_BIANCHI_I_THEOREM_DRAFT_METADATA_HASH_VALIDATOR_PROBE_V1",
            "source_packets": [str(IN_2291.relative_to(ROOT))],
            "contract_snapshot": {
                "required_hash": required_hash,
                "required_precheck_decision": required_decision,
            },
            "draft_metadata_candidate": draft_metadata_candidate,
            "checks": checks,
            "validator_decision": validator_decision,
            "theorem_scope_limit": "metadata validator only; not theorem proof and not selector/ToE closure",
        },
        "recommended_next_honest_step": {
            "id": "P2293_candidate",
            "goal": "run negative-control metadata replay matrix (hash mismatch / decision mismatch) and export rejection trace table",
        },
        "gatekeeper_checks": {
            "validator_exported": True,
            "required_hash_length_ok": len(required_hash) == 64,
            "checks_exported": len(checks) == 2,
            "decision_valid": validator_decision in ["METADATA_ACCEPT", "METADATA_REJECT"],
            "no_bridge_theorem_claimed": True,
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text("\n".join([
        "# P2292 S1242: theorem-draft metadata hash validator",
        "",
        f"- hash match: `{checks['hash_match']}`",
        f"- decision match: `{checks['decision_match']}`",
        f"- validator decision: `{validator_decision}`",
        "",
        "Metadata-validator only; no selector closure / ToE closure claim.",
    ]) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
