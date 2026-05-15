#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1752 = GEN / "p1752_s702_strict_nonproxy_h1_4d_execution_trigger_checkpoint.json"
IN1753 = GEN / "p1753_s703_strict_full_chain_forward_reverse_state_vector_checkpoint.json"
OUT = GEN / "p1754_s704_strict_nonproxy_minimal_delivery_manifest_checkpoint.json"


def main() -> None:
    p1752 = json.loads(IN1752.read_text(encoding="utf-8"))
    p1753 = json.loads(IN1753.read_text(encoding="utf-8"))

    missing = p1752.get("nonproxy_h1_4d_trigger_gate", {}).get("missing", [])

    manifest = {
        "M1_explicit_covariant_E_A_mu_expression_nonproxy": {
            "required_for": ["R3_nonproxy_h1_4d_trigger", "H1_strict_local"],
            "status": "MISSING" if "explicit_covariant_E_A_mu_expression_nonproxy" in missing else "EXPORTED",
        },
        "M2_explicit_covariant_E_H_expression_nonproxy": {
            "required_for": ["R3_nonproxy_h1_4d_trigger", "H1_strict_local"],
            "status": "MISSING" if "explicit_covariant_E_H_expression_nonproxy" in missing else "EXPORTED",
        },
        "M3_shared_background_family_contract": {
            "required_for": ["H1_strict_local", "H1_weak_form", "metric_residual"],
            "status": "MISSING" if "shared_background_family_contract" in missing else "EXPORTED",
        },
        "M4_boundary_term_control_clause": {
            "required_for": ["H1_weak_form", "promotion_safety"],
            "status": "MISSING" if "boundary_term_control_clause" in missing else "EXPORTED",
        },
        "M5_boundary_control_contract": {
            "required_for": ["H1_4D_trigger_gate"],
            "status": "MISSING" if "boundary_control_contract" in missing else "EXPORTED",
        },
    }

    ready = all(v["status"] == "EXPORTED" for v in manifest.values())

    payload = {
        "checkpoint": "P1754_S704",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "strict_only": True,
        "legacy_bridge_used": False,
        "chain": "K_strict -> coefficients -> full non-skeleton L_total -> EOM -> reverse minimal nonproxy delivery manifest",
        "input_anchors": ["p1752_s702", "p1753_s703"],
        "forward_reverse_state_vector": p1753.get("forward_reverse_state_vector", {}),
        "minimal_nonproxy_delivery_manifest": manifest,
        "manifest_ready_for_h1_4d": ready,
        "status": "READY_FOR_NONPROXY_H1_4D" if ready else "KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED",
        "no_false_pass_claim": True,
        "next_honest_step": "Uzupełnić wszystkie pozycje M1..M5 i dopiero uruchomić nonproxy H1 4D; następnie przejść do metrycznego EL_g-E_munu w tej samej rodzinie teł.",
        "lay_summary": "To jest lista kontrolna minimalnych braków. Dopóki choć jedna pozycja jest pusta, teoria nie może uczciwie przejść do kluczowego testu 4D.",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
