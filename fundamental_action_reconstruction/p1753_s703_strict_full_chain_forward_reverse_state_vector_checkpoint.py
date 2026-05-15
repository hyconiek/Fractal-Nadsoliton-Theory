#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1728 = GEN / "p1728_s678_strict_full_lagrangian_non_skeleton_and_bidirectional_closure_gap_checkpoint.json"
IN1752 = GEN / "p1752_s702_strict_nonproxy_h1_4d_execution_trigger_checkpoint.json"
OUT = GEN / "p1753_s703_strict_full_chain_forward_reverse_state_vector_checkpoint.json"


def main() -> None:
    p1728 = json.loads(IN1728.read_text(encoding="utf-8"))
    p1752 = json.loads(IN1752.read_text(encoding="utf-8"))

    forward = p1728.get("bidirectional_witness_map", {})
    h1_gate = p1752.get("nonproxy_h1_4d_trigger_gate", {})

    state_vector = {
        "S1_kernel_to_coefficients": forward.get("forward_kernel_to_coefficients", "UNKNOWN"),
        "S2_coefficients_to_full_lagrangian": forward.get("forward_coefficients_to_full_lagrangian", "UNKNOWN"),
        "S3_lagrangian_to_eom": forward.get("forward_lagrangian_to_eom", "UNKNOWN"),
        "R1_eom_to_variational_origin": forward.get("reverse_eom_to_variational_origin", "UNKNOWN"),
        "R2_variational_origin_to_kernel_identifiability": forward.get(
            "reverse_variational_origin_to_kernel_identifiability", "UNKNOWN"
        ),
        "R3_nonproxy_h1_4d_trigger": "READY" if h1_gate.get("trigger_ready") else "BLOCKED",
    }

    blockers = h1_gate.get("missing", [])
    closure_ready = state_vector["R3_nonproxy_h1_4d_trigger"] == "READY"

    payload = {
        "checkpoint": "P1753_S703",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "strict_only": True,
        "legacy_bridge_used": False,
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "chain": "K_strict -> coefficients -> full non-skeleton L_total -> EOM -> reverse-state-vector",
        "full_lagrangian_density_nonskeleton_instantiated": p1728.get(
            "full_lagrangian_density_nonskeleton_instantiated", {}
        ),
        "forward_reverse_state_vector": state_vector,
        "closure_readiness": {
            "strict_core_ready_for_nonproxy_h1_4d": closure_ready,
            "blocking_items": blockers,
            "status": "KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED" if not closure_ready else "READY_FOR_H1_4D_RUN_ONLY",
        },
        "no_false_pass_claim": True,
        "next_honest_step": "Zrealizować brakujące eksporty blokujące R3 i uruchomić pierwszy nonproxy H1 4D; dopiero potem przejść do sprzężenia z metrycznym EL_g-E_munu theorem gate.",
        "lay_summary": "Łańcuch do równań ruchu jest mocny, ale część odwrotna nadal czeka na kilka brakujących elementów. Ten raport zbiera to w jedną tablicę stanu, żeby nie tracić czasu na powtarzanie tych samych kroków.",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
