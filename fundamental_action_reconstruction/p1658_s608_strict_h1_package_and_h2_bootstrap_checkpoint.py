#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1655 = GEN / "p1655_s605_strict_helmholtz_h1_local_witness_summary.json"
IN1656 = GEN / "p1656_s606_strict_helmholtz_h1_gauge_scalar_covariant_summary.json"
IN1657 = GEN / "p1657_s607_strict_helmholtz_h1_gauge_metric_covariant_summary.json"


def main() -> None:
    s55 = json.loads(IN1655.read_text(encoding="utf-8"))
    s56 = json.loads(IN1656.read_text(encoding="utf-8"))
    s57 = json.loads(IN1657.read_text(encoding="utf-8"))

    summary = {
        "checkpoint": "P1658_S608",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": "KEEP_OPEN_P1658_H1_PACKAGE_INTEGRATED_H2_BOOTSTRAP_EXPORTED",
        "strict_only": True,
        "legacy_bridge_used": False,
        "kernel_definition": s57["kernel_definition"],
        "forward_chain_context": "K_strict -> coefficients -> full L_total -> EOM",
        "reverse_gate_target": "EOM -> L_total Helmholtz H1..H4 discharge",
        "h1_package": {
            "phi_h_local": s55["helmholtz_h1_local"]["status"],
            "gauge_scalar_covariant": s56["h1_local_covariant_check"]["status"],
            "gauge_metric_covariant": s57["h1_local_gauge_metric_check"]["status"],
            "integration_status": "PARTIAL",
        },
        "h2_bootstrap": {
            "condition": "self-adjointness of linearized EOM operator on coupled field multiplet",
            "exported_schema": [
                "define linearized operator L_E around background-agnostic chart",
                "check <u, L_E v> - <L_E u, v> boundary/covariant terms",
                "separate gauge-fixing dependent and invariant pieces",
            ],
            "status": "PARTIAL",
        },
        "qg_gates": {
            "renormalization": "OPEN",
            "unitarity": "OPEN",
            "background_independence": "OPEN",
        },
        "final_strict_core_closure": {
            "status": "OPEN",
            "reason": "H2/H3/H4 and QG theorem-level exports missing",
        },
        "next_honest_step": "Wyeksportować pierwszy jawny rachunek dla H2: bilans <u,L_E v>-<L_E u,v> w bloku gauge-scalar z rozdzieleniem członów brzegowych.",
        "lay_summary": "Połączyliśmy trzy lokalne testy zgodności H1 w jedną całość i przygotowaliśmy plan następnego testu H2. To ważny postęp, ale pełne domknięcie teorii nadal wymaga dalszych formalnych dowodów.",
    }

    out = GEN / "p1658_s608_strict_h1_package_and_h2_bootstrap_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
