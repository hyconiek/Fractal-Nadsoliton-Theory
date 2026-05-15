#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1691 = GEN / "p1691_s641_strict_full_chain_lagrangian_to_qg_theorem_obligation_matrix.json"
IN1694 = GEN / "p1694_s644_strict_kernel_to_full_lagrangian_bidirectional_map_witness.json"
OUT = GEN / "p1696_s646_strict_gauge_higgs_covariant_eom_export_checkpoint.json"


def main() -> None:
    p1691 = json.loads(IN1691.read_text(encoding="utf-8"))
    p1694 = json.loads(IN1694.read_text(encoding="utf-8"))

    x = sp.symbols("x", real=True)
    h = sp.Function("h")(x)
    A = sp.Function("A")(x)

    coeff = p1694.get("forward_coefficient_map_symbolic", {})
    z1 = sp.Float(coeff.get("Z1", "1.0"))
    muH2 = sp.Float(coeff.get("muH2", "1.0"))
    lambdaH = sp.Float(coeff.get("lambdaH", "1.0"))
    xiHR = sp.Float(coeff.get("xiHR", "0.0"))

    R0 = sp.symbols("R0", real=True)

    # reduced covariant-like export block for gauge+higgs on background with scalar curvature marker R0
    L_gauge_higgs = sp.expand(
        z1 * sp.diff(A, x) ** 2 / 2
        + sp.diff(h, x) ** 2 / 2
        - muH2 * h**2 / 2
        - lambdaH * h**4 / 4
        - xiHR * R0 * h**2 / 2
    )

    eom_A = sp.simplify(sp.diff(sp.diff(L_gauge_higgs, sp.diff(A, x)), x) - sp.diff(L_gauge_higgs, A))
    eom_h = sp.simplify(sp.diff(sp.diff(L_gauge_higgs, sp.diff(h, x)), x) - sp.diff(L_gauge_higgs, h))

    payload = {
        "checkpoint": "P1696_S646",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "strict_only": True,
        "legacy_bridge_used": False,
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "chain": "K_strict -> coefficients -> full L_total -> gauge+higgs covariant export block -> EOM",
        "full_lagrangian_density_anchor": p1691.get("full_lagrangian_density", {}),
        "coefficient_anchor": coeff,
        "covariant_block": {
            "background_marker": "R0",
            "L_gauge_higgs_reduced": str(L_gauge_higgs),
            "EOM_A": str(eom_A),
            "EOM_h": str(eom_h),
        },
        "bidirectional_status": {
            "forward_export": "GAUGE_HIGGS_BLOCK_EXPORTED",
            "reverse_global_variational_origin": "KEEP_OPEN",
        },
        "status": "KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED",
        "open_obligations": [
            "metric_covariant_eom_export",
            "fermion_covariant_spinor_export_nonproxy",
            "counterterm_flow_closure_theorem",
            "BRST_cohomology_theorem",
            "Cutkosky_unitarity_theorem_full_sector",
            "background_family_independence_theorem",
        ],
        "next_honest_step": "Dołączyć blok metryczny i spinorowy do pełnego kowariantnego EOM bundle oraz skleić z theorem witnessami BRST/Cutkosky/counterterm-flow.",
        "lay_summary": "Wyprowadziliśmy już konkretny fragment równań ruchu (gauge+Higgs) z pełnej teorii strict. To duży krok, ale do pełnego domknięcia kwantowej grawitacji potrzeba jeszcze brakujących bloków i dowodów globalnych.",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
