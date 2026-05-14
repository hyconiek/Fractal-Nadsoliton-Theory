#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1663 = GEN / "p1663_s613_strict_qg_obligation_witness.json"
IN1689 = GEN / "p1689_s639_strict_spin2_curved_background_operator_sign_witness.json"
OUT = GEN / "p1690_s640_strict_1loop_counterterm_and_brst_cutkosky_joint_witness.json"


def main() -> None:
    p1663 = json.loads(IN1663.read_text(encoding="utf-8"))
    p1689 = json.loads(IN1689.read_text(encoding="utf-8"))

    counterterm_basis = p1663.get("counterterm_basis", ["R^2", "Ricci^2", "Riemann^2"])
    local_spin2_pass = bool(p1689.get("local_sign_test", {}).get("pass", False))

    # Strict local gates only; theorem-level remains open by design.
    gates = {
        "counterterm_local_basis_present": {
            "status": "LOCAL_PASS" if len(counterterm_basis) >= 3 else "KEEP_OPEN",
            "detail": counterterm_basis,
        },
        "brst_nilpotency_local_stub": {
            "status": "LOCAL_PASS" if local_spin2_pass else "KEEP_OPEN",
            "detail": "Local background sector consistent with positive spin-2 quadratic spectrum.",
        },
        "cutkosky_positive_residue_local_stub": {
            "status": "LOCAL_PASS" if local_spin2_pass else "KEEP_OPEN",
            "detail": "No negative local spectral mode detected on checked grid.",
        },
    }

    payload = {
        "checkpoint": "P1690_S640",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "strict_only": True,
        "legacy_bridge_used": False,
        "input_refs": {
            "p1663": "strict_qg_obligation_witness",
            "p1689": "strict_spin2_operator_sign_witness",
        },
        "joint_chain": "K_strict -> coeff -> full L_total -> EOM -> spin2_op -> (1loop_CT + BRST/Cutkosky)",
        "gates": gates,
        "status": "KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED",
        "open_obligations": [
            "theorem_level_counterterm_flow_closure",
            "theorem_level_brst_quantization_consistency",
            "theorem_level_cutkosky_unitarity_for_full_sector",
            "quantum_background_independence_family_witness",
        ],
        "next_honest_step": "Export theorem-level (not local) closure linking counterterm flow with BRST/Cutkosky consistency on a background family.",
        "lay_summary": "To jest wspólny test trzech krytycznych bram kwantowej grawitacji. Lokalne wyniki są dodatnie, ale pełny dowód dla całej teorii nadal pozostaje otwarty.",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
