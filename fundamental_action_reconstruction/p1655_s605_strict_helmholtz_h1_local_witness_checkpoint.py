#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1654 = GEN / "p1654_s604_strict_bidirectional_theorem_requirement_matrix_with_qg_summary.json"


def main() -> None:
    s54 = json.loads(IN1654.read_text(encoding="utf-8"))

    # V(phi,h)=mphi2/2*phi^2 + lam4/24*phi^4 + muH2*h^2 + lamH*h^4 + lamPH*phi^2*h^2
    # dV/dphi = mphi2*phi + lam4/6*phi^3 + 2*lamPH*phi*h^2
    # dV/dh   = 2*muH2*h + 4*lamH*h^3 + 2*lamPH*phi^2*h
    # d/dh(dV/dphi) = 4*lamPH*phi*h
    # d/dphi(dV/dh) = 4*lamPH*phi*h
    cross_1 = "4*lamPH*phi*h"
    cross_2 = "4*lamPH*phi*h"
    residual = "0"

    summary = {
        "checkpoint": "P1655_S605",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": "KEEP_OPEN_P1655_LOCAL_H1_WITNESS_EXPORTED",
        "strict_only": True,
        "legacy_bridge_used": False,
        "inherits_from": s54["checkpoint"],
        "kernel_definition": s54["kernel_definition"],
        "reverse_gate_target": "EOM -> L_total Helmholtz H1..H4 discharge",
        "local_sector": "scalar phi + Higgs-radial h with lambda_phiH mixing",
        "helmholtz_h1_local": {
            "cross_derivative_dphi_dh": cross_1,
            "cross_derivative_dh_dphi": cross_2,
            "residual": residual,
            "status": "PARTIAL",
        },
        "global_reverse_chain_status": "OPEN",
        "open_items": [
            "Extend H1 witness to full gauge/fermion/metric coupled system",
            "Export H2/H3/H4 theorem-level witnesses",
            "Integrate Helmholtz package into EOM->L_total reconstruction theorem",
        ],
        "next_honest_step": "Rozszerzyć H1 z sektora phi-H do sektora gauge+scalar i dodać jawne operatory dywergencji kowariantnej.",
        "lay_summary": "Pokazaliśmy pierwszy matematyczny test cofania teorii: lokalnie równania zachowują symetrię pochodnych wymaganą do odzyskiwania lagranżianu. To dopiero początek — pełny układ nadal jest otwarty.",
    }

    out = GEN / "p1655_s605_strict_helmholtz_h1_local_witness_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
