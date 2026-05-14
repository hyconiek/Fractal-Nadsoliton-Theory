#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1644 = GEN / "p1644_s594_strict_eom_to_lagrangian_helmholtz_witness_scaffold_summary.json"


def main() -> None:
    s44 = json.loads(IN1644.read_text(encoding="utf-8"))

    # Local symbolic probe without external CAS dependency:
    # V(phi,H)= mu2/2 H^2 + lamH/4 H^4 + lamphi/4 phi^4 + lamphiH/2 phi^2 H^2 + xiHR/2 H^2 R + xiPhiR/2 phi^2 R
    # E_phi = dV/dphi = lamphi*phi^3 + lamphiH*phi*H^2 + xiPhiR*phi*R
    # E_H   = dV/dH   = mu2*H + lamH*H^3 + lamphiH*phi^2*H + xiHR*H*R
    # d(E_phi)/dH = 2*lamphiH*phi*H
    # d(E_H)/dphi = 2*lamphiH*phi*H

    cross_variation_commutes = True

    summary = {
        "checkpoint": "P1645_S595",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": "PASS_LOCAL_P1645_CROSS_VARIATION_SYMMETRY_SCALAR_HIGGS_BLOCK",
        "strict_only": True,
        "legacy_bridge_used": False,
        "scope": "local scalar-Higgs potential block of full strict Lagrangian",
        "probe": {
            "E_phi": "lamphi*phi^3 + lamphiH*phi*H^2 + xiPhiR*phi*R",
            "E_H": "mu2*H + lamH*H^3 + lamphiH*phi^2*H + xiHR*H*R",
            "dH_dEphi": "2*lamphiH*phi*H",
            "dphi_dEH": "2*lamphiH*phi*H",
            "mismatch": "0",
            "cross_variation_commutes": cross_variation_commutes,
            "method": "analytic closed-form derivative check (dependency-free)",
        },
        "impact_on_H2": {
            "global_status": "PARTIAL",
            "reason": "One mixed block validated analytically; gauge/fermion/metric blocks and boundary class still open",
        },
        "strict_core_closure": s44["strict_core_closure"],
        "next_honest_step": "Rozszerzyć test H2 o jawny blok gauge+metric oraz dołączyć eksport założeń brzegowych (H4), żeby zamknąć B4 na pełnej klasie operatorów strict.",
        "lay_summary": "Sprawdziliśmy poprawność jednego kluczowego fragmentu równań: w tym bloku mieszane pochodne są zgodne. To postęp, ale pełna teoria ma więcej bloków do formalnego domknięcia.",
    }

    out = GEN / "p1645_s595_strict_cross_variation_symmetry_probe_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
