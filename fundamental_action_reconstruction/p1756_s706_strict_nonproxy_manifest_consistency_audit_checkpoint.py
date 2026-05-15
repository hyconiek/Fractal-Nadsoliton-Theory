#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1752 = GEN / "p1752_s702_strict_nonproxy_h1_4d_execution_trigger_checkpoint.json"
IN1754 = GEN / "p1754_s704_strict_nonproxy_minimal_delivery_manifest_checkpoint.json"
OUT = GEN / "p1756_s706_strict_nonproxy_manifest_consistency_audit_checkpoint.json"


def main() -> None:
    p1752 = json.loads(IN1752.read_text(encoding="utf-8"))
    p1754 = json.loads(IN1754.read_text(encoding="utf-8"))

    missing_trigger = set(p1752.get("nonproxy_h1_4d_trigger_gate", {}).get("missing", []))
    manifest = p1754.get("minimal_nonproxy_delivery_manifest", {})

    key_map = {
        "M1_explicit_covariant_E_A_mu_expression_nonproxy": "explicit_covariant_E_A_mu_expression_nonproxy",
        "M2_explicit_covariant_E_H_expression_nonproxy": "explicit_covariant_E_H_expression_nonproxy",
        "M3_shared_background_family_contract": "shared_background_family_contract",
        "M4_boundary_term_control_clause": "boundary_term_control_clause",
        "M5_boundary_control_contract": "boundary_control_contract",
    }

    mismatches = []
    for m_key, t_key in key_map.items():
        m_status = manifest.get(m_key, {}).get("status")
        expected = "MISSING" if t_key in missing_trigger else "EXPORTED"
        if m_status != expected:
            mismatches.append({"manifest_key": m_key, "trigger_key": t_key, "manifest_status": m_status, "expected": expected})

    payload = {
        "checkpoint": "P1756_S706",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "strict_only": True,
        "legacy_bridge_used": False,
        "chain": "K_strict -> coefficients -> full non-skeleton L_total -> EOM -> trigger/manifest consistency audit",
        "input_anchors": ["p1752_s702", "p1754_s704"],
        "consistency_result": {
            "manifest_keys_checked": list(key_map.keys()),
            "trigger_missing_keys": sorted(missing_trigger),
            "mismatches": mismatches,
            "status": "PASS_CONSISTENT" if not mismatches else "OBSTRUCTION_INCONSISTENT_MANIFEST",
        },
        "status": "KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED",
        "no_false_pass_claim": True,
        "next_honest_step": "Jeśli status PASS_CONSISTENT, dostarczyć M1..M5; jeśli OBSTRUCTION, najpierw zsynchronizować manifest z trigger gate i dopiero kontynuować 4D H1.",
        "lay_summary": "Sprawdziliśmy, czy lista braków i bramka uruchomienia mówią to samo. Dzięki temu nie powstaje chaos proceduralny przed krytycznym testem 4D.",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
