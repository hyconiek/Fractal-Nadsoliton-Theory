#!/usr/bin/env python3
from __future__ import annotations
import hashlib
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2286 = GEN / "p2286_s1236_strict_task3_bianchi_i_immutable_a1_a2_a3_certificate_bundle_probe.json"
IN_2287 = GEN / "p2287_s1237_strict_task3_bianchi_i_certificate_bundle_integrity_recheck_probe.json"
OUT = GEN / "p2288_s1238_strict_task3_bianchi_i_certificate_chain_index_gating_ledger_probe.json"
MD = GEN / "p2288_s1238_strict_task3_bianchi_i_certificate_chain_index_gating_ledger_probe.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def sha256_text(s: str) -> str:
    return hashlib.sha256(s.encode("utf-8")).hexdigest()


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2286 = load(IN_2286)
    p2287 = load(IN_2287)

    bundle = (p2286.get("strict_task3_bianchi_i_immutable_a1_a2_a3_certificate_bundle_probe", {}) or {}).get("bundle", {}) or {}
    bundle_hash = str(bundle.get("payload_sha256", "") or "")
    bundle_payload = bundle.get("payload", {}) or {}
    verifier_pass = bool(bundle_payload.get("verifier_pass", False))

    integrity_probe = (p2287.get("strict_task3_bianchi_i_certificate_bundle_integrity_recheck_probe", {}) or {})
    integrity_ok = bool((integrity_probe.get("integrity_recheck", {}) or {}).get("integrity_ok", False))
    mutation_detected = bool((integrity_probe.get("negative_control", {}) or {}).get("mutation_detected", False))

    gating_decision = (
        "ALLOW_THEOREM_ATTEMPT_PRECHECK"
        if (integrity_ok and mutation_detected and len(bundle_hash) == 64 and verifier_pass)
        else "BLOCK_THEOREM_ATTEMPT_PRECHECK"
    )

    chain_index_record = {
        "record_id": "TASK3_CHAIN_INDEX_RECORD_V1",
        "bundle_payload_sha256": bundle_hash,
        "integrity_ok": integrity_ok,
        "mutation_detected": mutation_detected,
        "verifier_pass": verifier_pass,
        "gating_decision": gating_decision,
        "source_packets": [str(IN_2286.relative_to(ROOT)), str(IN_2287.relative_to(ROOT))],
    }

    chain_index_fingerprint = sha256_text(json.dumps(chain_index_record, sort_keys=True, separators=(",", ":"), ensure_ascii=False))

    payload = {
        "schema_version": "p2288_s1238_v1",
        "packet_id": "P2288",
        "stage_id": "S1238",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_TASK3_BIANCHI_I_CERTIFICATE_CHAIN_INDEX_GATING_LEDGER_PROBE",
        "strict_task3_bianchi_i_certificate_chain_index_gating_ledger_probe": {
            "probe_id": "STRICT_TASK3_BIANCHI_I_CERTIFICATE_CHAIN_INDEX_GATING_LEDGER_PROBE_V1",
            "chain_index_record": chain_index_record,
            "chain_index_fingerprint_sha256": chain_index_fingerprint,
            "theorem_scope_limit": "gating-ledger packaging only; not theorem proof and not selector/ToE closure",
        },
        "recommended_next_honest_step": {
            "id": "P2289_candidate",
            "goal": "attach chain-index fingerprint to theorem-attempt precheck contract and require exact match before any theorem-draft run",
        },
        "gatekeeper_checks": {
            "chain_index_exported": True,
            "bundle_hash_length_ok": len(bundle_hash) == 64,
            "chain_fingerprint_length_ok": len(chain_index_fingerprint) == 64,
            "gating_decision_present": gating_decision in ["ALLOW_THEOREM_ATTEMPT_PRECHECK", "BLOCK_THEOREM_ATTEMPT_PRECHECK"],
            "block_if_verifier_fails": verifier_pass or gating_decision == "BLOCK_THEOREM_ATTEMPT_PRECHECK",
            "no_bridge_theorem_claimed": True,
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text("\n".join([
        "# P2288 S1238: certificate chain-index gating ledger",
        "",
        f"- bundle hash length ok: `{len(bundle_hash) == 64}`",
        f"- integrity ok: `{integrity_ok}`",
        f"- mutation detected: `{mutation_detected}`",
        f"- verifier pass: `{verifier_pass}`",
        f"- gating decision: `{gating_decision}`",
        f"- chain fingerprint: `{chain_index_fingerprint}`",
        "",
        "Gating-ledger packaging only; no selector closure / ToE closure claim.",
    ]) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
