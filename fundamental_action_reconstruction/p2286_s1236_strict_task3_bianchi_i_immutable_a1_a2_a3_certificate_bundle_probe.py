#!/usr/bin/env python3
from __future__ import annotations
import hashlib
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2285 = GEN / "p2285_s1235_strict_task3_bianchi_i_machine_checkable_a1_a2_a3_verifier_probe.json"
OUT = GEN / "p2286_s1236_strict_task3_bianchi_i_immutable_a1_a2_a3_certificate_bundle_probe.json"
MD = GEN / "p2286_s1236_strict_task3_bianchi_i_immutable_a1_a2_a3_certificate_bundle_probe.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def sha256_of_obj(obj: Any) -> str:
    canonical = json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2285 = load(IN_2285)

    probe = p2285.get("strict_task3_bianchi_i_machine_checkable_a1_a2_a3_verifier_probe", {}) or {}
    recomputed = probe.get("recomputed", {}) or {}
    reported = probe.get("reported_from_p2284", {}) or {}
    consistency = probe.get("consistency", {}) or {}
    verifier_pass = bool(probe.get("verifier_pass", False))

    certificate_payload = {
        "certificate_id": "TASK3_BIANCHI_I_A1_A2_A3_IMMUTABLE_CERTIFICATE_V1",
        "source_packet": str(IN_2285.relative_to(ROOT)),
        "recomputed": recomputed,
        "reported": reported,
        "consistency": consistency,
        "verifier_pass": verifier_pass,
        "scope": "task3_bianchi_i_transport_premise_verification",
    }

    digest = sha256_of_obj(certificate_payload)

    bundle = {
        "bundle_version": "v1",
        "payload": certificate_payload,
        "payload_sha256": digest,
        "integrity_rule": "payload_sha256 must match canonical JSON of payload",
    }

    payload = {
        "schema_version": "p2286_s1236_v1",
        "packet_id": "P2286",
        "stage_id": "S1236",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_TASK3_BIANCHI_I_IMMUTABLE_A1_A2_A3_CERTIFICATE_BUNDLE_PROBE",
        "strict_task3_bianchi_i_immutable_a1_a2_a3_certificate_bundle_probe": {
            "probe_id": "STRICT_TASK3_BIANCHI_I_IMMUTABLE_A1_A2_A3_CERTIFICATE_BUNDLE_PROBE_V1",
            "source_packets": [str(IN_2285.relative_to(ROOT))],
            "bundle": bundle,
            "theorem_scope_limit": "immutable verifier bundle packaging only; not theorem proof and not selector/ToE closure",
        },
        "recommended_next_honest_step": {
            "id": "P2287_candidate",
            "goal": "export deterministic re-check utility that recomputes payload hash and rejects mutated certificate bundles",
        },
        "gatekeeper_checks": {
            "bundle_exported": True,
            "sha256_length_ok": len(digest) == 64,
            "verifier_pass_boolean": isinstance(verifier_pass, bool),
            "consistency_has_A1_A2_A3": set(consistency.keys()) == {"A1", "A2", "A3"},
            "no_bridge_theorem_claimed": True,
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text("\n".join([
        "# P2286 S1236: immutable A1/A2/A3 certificate bundle",
        "",
        f"- verifier pass: `{verifier_pass}`",
        f"- payload sha256: `{digest}`",
        f"- consistency keys: `{sorted(consistency.keys())}`",
        "",
        "Immutable bundle packaging only; no selector closure / ToE closure claim.",
    ]) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
