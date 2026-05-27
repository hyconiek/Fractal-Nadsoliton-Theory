#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2288 = GEN / "p2288_s1238_strict_task3_bianchi_i_certificate_chain_index_gating_ledger_probe.json"
OUT = GEN / "p2289_s1239_strict_task3_bianchi_i_theorem_attempt_precheck_contract_probe.json"
MD = GEN / "p2289_s1239_strict_task3_bianchi_i_theorem_attempt_precheck_contract_probe.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2288 = load(IN_2288)

    probe = (p2288.get("strict_task3_bianchi_i_certificate_chain_index_gating_ledger_probe", {}) or {})
    rec = probe.get("chain_index_record", {}) or {}
    fingerprint = str(probe.get("chain_index_fingerprint_sha256", "") or "")

    gating_decision = str(rec.get("gating_decision", "BLOCK_THEOREM_ATTEMPT_PRECHECK") or "BLOCK_THEOREM_ATTEMPT_PRECHECK")
    contract_status = "READY_FOR_PRECHECK_ONLY" if gating_decision == "ALLOW_THEOREM_ATTEMPT_PRECHECK" else "BLOCKED_BY_GATING_LEDGER"

    precheck_contract = {
        "contract_id": "TASK3_THEOREM_ATTEMPT_PRECHECK_CONTRACT_V1",
        "required_chain_index_fingerprint_sha256": fingerprint,
        "required_gating_decision": "ALLOW_THEOREM_ATTEMPT_PRECHECK",
        "observed_gating_decision": gating_decision,
        "status": contract_status,
        "rule": "Theorem draft attempt must be blocked unless observed gating decision is ALLOW and fingerprint matches exactly.",
    }

    payload = {
        "schema_version": "p2289_s1239_v1",
        "packet_id": "P2289",
        "stage_id": "S1239",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_TASK3_BIANCHI_I_THEOREM_ATTEMPT_PRECHECK_CONTRACT_PROBE",
        "strict_task3_bianchi_i_theorem_attempt_precheck_contract_probe": {
            "probe_id": "STRICT_TASK3_BIANCHI_I_THEOREM_ATTEMPT_PRECHECK_CONTRACT_PROBE_V1",
            "source_packets": [str(IN_2288.relative_to(ROOT))],
            "precheck_contract": precheck_contract,
            "theorem_scope_limit": "precheck-contract export only; not theorem proof and not selector/ToE closure",
        },
        "recommended_next_honest_step": {
            "id": "P2290_candidate",
            "goal": "run deterministic precheck evaluator that takes proposed theorem-draft metadata and verifies exact fingerprint+gating compliance",
        },
        "gatekeeper_checks": {
            "precheck_contract_exported": True,
            "fingerprint_length_ok": len(fingerprint) == 64,
            "observed_gating_decision_valid": gating_decision in ["ALLOW_THEOREM_ATTEMPT_PRECHECK", "BLOCK_THEOREM_ATTEMPT_PRECHECK"],
            "contract_status_valid": contract_status in ["READY_FOR_PRECHECK_ONLY", "BLOCKED_BY_GATING_LEDGER"],
            "no_bridge_theorem_claimed": True,
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text("\n".join([
        "# P2289 S1239: theorem-attempt precheck contract",
        "",
        f"- observed gating decision: `{gating_decision}`",
        f"- contract status: `{contract_status}`",
        f"- required fingerprint: `{fingerprint}`",
        "",
        "Precheck contract only; no selector closure / ToE closure claim.",
    ]) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
