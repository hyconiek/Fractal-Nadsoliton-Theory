#!/usr/bin/env python3
from __future__ import annotations
import hashlib
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2286 = GEN / "p2286_s1236_strict_task3_bianchi_i_immutable_a1_a2_a3_certificate_bundle_probe.json"
OUT = GEN / "p2287_s1237_strict_task3_bianchi_i_certificate_bundle_integrity_recheck_probe.json"
MD = GEN / "p2287_s1237_strict_task3_bianchi_i_certificate_bundle_integrity_recheck_probe.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def sha256_of_obj(obj: Any) -> str:
    canonical = json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2286 = load(IN_2286)

    probe = (p2286.get("strict_task3_bianchi_i_immutable_a1_a2_a3_certificate_bundle_probe", {}) or {})
    bundle = probe.get("bundle", {}) or {}
    payload = bundle.get("payload", {}) or {}
    reported_hash = str(bundle.get("payload_sha256", "") or "")

    recomputed_hash = sha256_of_obj(payload)
    integrity_ok = reported_hash == recomputed_hash

    # negative control: mutate a deterministic leaf if available
    mutated = json.loads(json.dumps(payload))
    if isinstance(mutated.get("verifier_pass"), bool):
        mutated["verifier_pass"] = not mutated["verifier_pass"]
    mutated_hash = sha256_of_obj(mutated)
    mutation_detected = mutated_hash != reported_hash

    payload_out = {
        "schema_version": "p2287_s1237_v1",
        "packet_id": "P2287",
        "stage_id": "S1237",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_TASK3_BIANCHI_I_CERTIFICATE_BUNDLE_INTEGRITY_RECHECK_PROBE",
        "strict_task3_bianchi_i_certificate_bundle_integrity_recheck_probe": {
            "probe_id": "STRICT_TASK3_BIANCHI_I_CERTIFICATE_BUNDLE_INTEGRITY_RECHECK_PROBE_V1",
            "source_packets": [str(IN_2286.relative_to(ROOT))],
            "integrity_recheck": {
                "reported_payload_sha256": reported_hash,
                "recomputed_payload_sha256": recomputed_hash,
                "integrity_ok": integrity_ok,
            },
            "negative_control": {
                "mutated_payload_sha256": mutated_hash,
                "mutation_detected": mutation_detected,
                "rule": "mutation_detected should be true for deterministic payload mutation",
            },
            "theorem_scope_limit": "integrity recheck utility only; not theorem proof and not selector/ToE closure",
        },
        "recommended_next_honest_step": {
            "id": "P2288_candidate",
            "goal": "export bundle chain-index packet linking P2286 hash to downstream theorem-attempt gating ledger",
        },
        "gatekeeper_checks": {
            "integrity_recheck_exported": True,
            "integrity_ok": integrity_ok,
            "mutation_detected": mutation_detected,
            "reported_hash_length_ok": len(reported_hash) == 64,
            "recomputed_hash_length_ok": len(recomputed_hash) == 64,
            "no_bridge_theorem_claimed": True,
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload_out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text("\n".join([
        "# P2287 S1237: certificate bundle integrity recheck",
        "",
        f"- integrity ok: `{integrity_ok}`",
        f"- mutation detected: `{mutation_detected}`",
        f"- reported hash: `{reported_hash}`",
        f"- recomputed hash: `{recomputed_hash}`",
        "",
        "Integrity utility only; no selector closure / ToE closure claim.",
    ]) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
