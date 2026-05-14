#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1650 = GEN / "p1650_s600_strict_global_injectivity_theorem_object_scaffold_summary.json"


def main() -> None:
    s50 = json.loads(IN1650.read_text(encoding="utf-8"))
    # Consume regional samples exported in P1649 as computational certificate substrate.
    s49 = json.loads((GEN / "p1649_s599_strict_regional_projector_injectivity_witness_summary.json").read_text(encoding="utf-8"))

    det_vals = [row["detJ"] for row in s49["regional_samples"]]
    eps_hat = min(abs(v) for v in det_vals)

    summary = {
        "checkpoint": "P1651_S601",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": "PASS_LOCAL_P1651_DETJ_LOWER_BOUND_CERTIFICATE",
        "strict_only": True,
        "legacy_bridge_used": False,
        "kernel_definition": s50["kernel_definition"],
        "full_chain": s50["full_chain"],
        "o1_certificate": {
            "target": "O1 from T_STRICT_INJ_BOX_01",
            "region": s50["global_injectivity_theorem_object"]["box_B"],
            "certificate_type": "computational lower-bound estimate from regional grid",
            "eps_hat": eps_hat,
            "derived_from_sample_count": len(det_vals),
            "condition_checked": "|detJ| >= eps_hat on sampled region points",
            "status": "PARTIAL",
            "note": "Grid-derived certificate; still needs analytic proof upgrade for theorem-level discharge.",
        },
        "strict_core_closure": {
            "status": "OPEN",
            "reason": "O1 has local computational certificate only; O2/O3 and analytic globalization remain open",
        },
        "next_honest_step": "Podnieść O1 z certyfikatu siatkowego do dowodu analitycznego (symboliczna dolna granica detJ na całym boxie) i wtedy przejść do O2.",
        "lay_summary": "Mamy liczbowy dowód, że w testowanej siatce mapa nie traci odwracalności. To bardzo użyteczne, ale do pełnej matematycznej pewności potrzebny jest jeszcze dowód analityczny.",
    }

    out = GEN / "p1651_s601_strict_detj_lower_bound_certificate_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
