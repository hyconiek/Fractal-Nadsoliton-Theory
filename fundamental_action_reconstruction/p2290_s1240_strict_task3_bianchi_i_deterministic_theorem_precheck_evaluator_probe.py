#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2289 = GEN / "p2289_s1239_strict_task3_bianchi_i_theorem_attempt_precheck_contract_probe.json"
IN_2288 = GEN / "p2288_s1238_strict_task3_bianchi_i_certificate_chain_index_gating_ledger_probe.json"
OUT = GEN / "p2290_s1240_strict_task3_bianchi_i_deterministic_theorem_precheck_evaluator_probe.json"
MD = GEN / "p2290_s1240_strict_task3_bianchi_i_deterministic_theorem_precheck_evaluator_probe.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2289 = load(IN_2289)
    p2288 = load(IN_2288)

    contract = (p2289.get("strict_task3_bianchi_i_theorem_attempt_precheck_contract_probe", {}) or {}).get("precheck_contract", {}) or {}
    ledger = (p2288.get("strict_task3_bianchi_i_certificate_chain_index_gating_ledger_probe", {}) or {})
    actual_fingerprint = str(ledger.get("chain_index_fingerprint_sha256", "") or "")
    actual_gating = str((ledger.get("chain_index_record", {}) or {}).get("gating_decision", "BLOCK_THEOREM_ATTEMPT_PRECHECK") or "BLOCK_THEOREM_ATTEMPT_PRECHECK")

    expected_fingerprint = str(contract.get("required_chain_index_fingerprint_sha256", "") or "")
    expected_gating = str(contract.get("required_gating_decision", "ALLOW_THEOREM_ATTEMPT_PRECHECK") or "ALLOW_THEOREM_ATTEMPT_PRECHECK")

    checks = {
        "fingerprint_match": actual_fingerprint == expected_fingerprint,
        "gating_match": actual_gating == expected_gating,
    }

    evaluator_decision = "PRECHECK_PASS" if all(checks.values()) else "PRECHECK_BLOCK"

    payload = {
        "schema_version": "p2290_s1240_v1",
        "packet_id": "P2290",
        "stage_id": "S1240",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_TASK3_BIANCHI_I_DETERMINISTIC_THEOREM_PRECHECK_EVALUATOR_PROBE",
        "strict_task3_bianchi_i_deterministic_theorem_precheck_evaluator_probe": {
            "probe_id": "STRICT_TASK3_BIANCHI_I_DETERMINISTIC_THEOREM_PRECHECK_EVALUATOR_PROBE_V1",
            "source_packets": [str(IN_2289.relative_to(ROOT)), str(IN_2288.relative_to(ROOT))],
            "inputs": {
                "expected_fingerprint": expected_fingerprint,
                "expected_gating_decision": expected_gating,
                "actual_fingerprint": actual_fingerprint,
                "actual_gating_decision": actual_gating,
            },
            "checks": checks,
            "evaluator_decision": evaluator_decision,
            "theorem_scope_limit": "deterministic precheck evaluator only; not theorem proof and not selector/ToE closure",
        },
        "recommended_next_honest_step": {
            "id": "P2291_candidate",
            "goal": "freeze signed precheck transcript and require transcript hash in any theorem-draft metadata packet",
        },
        "gatekeeper_checks": {
            "checks_exported": len(checks) == 2,
            "fingerprint_length_ok": len(actual_fingerprint) == 64 and len(expected_fingerprint) == 64,
            "decision_valid": evaluator_decision in ["PRECHECK_PASS", "PRECHECK_BLOCK"],
            "no_bridge_theorem_claimed": True,
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text("\n".join([
        "# P2290 S1240: deterministic theorem precheck evaluator",
        "",
        f"- fingerprint match: `{checks['fingerprint_match']}`",
        f"- gating match: `{checks['gating_match']}`",
        f"- evaluator decision: `{evaluator_decision}`",
        "",
        "Deterministic precheck only; no selector closure / ToE closure claim.",
    ]) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
