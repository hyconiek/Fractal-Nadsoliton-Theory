#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1696 = GEN / "p1696_s646_strict_gauge_higgs_covariant_eom_export_checkpoint.json"
OUT = GEN / "p1697_s647_strict_metric_spinor_covariant_eom_bundle_checkpoint.json"


def main() -> None:
    p1696 = json.loads(IN1696.read_text(encoding="utf-8"))
    coeff = p1696.get("coefficient_anchor", {})

    x = sp.symbols("x", real=True)
    g = sp.Function("g")(x)
    psi = sp.Function("psi")(x)

    cR = sp.Float(coeff.get("cR", "1.0"))
    xiGR = sp.Float(coeff.get("xiGR", "1.0"))
    ypsi = sp.Float(coeff.get("ypsi", "1.0"))
    mpsi = sp.Float(coeff.get("mpsi", "1.0"))

    R1 = sp.symbols("R1", real=True)

    # Reduced covariant-like metric + spinor export block (strict operational lane).
    L_metric_spinor = sp.expand(
        cR * sp.diff(g, x) ** 2 / 2
        - xiGR * R1 * g**2 / 2
        + sp.I * psi * sp.diff(psi, x)
        - mpsi * psi**2
        - ypsi * g * psi**2
    )

    eom_g = sp.simplify(sp.diff(sp.diff(L_metric_spinor, sp.diff(g, x)), x) - sp.diff(L_metric_spinor, g))
    eom_psi = sp.simplify(sp.diff(sp.diff(L_metric_spinor, sp.diff(psi, x)), x) - sp.diff(L_metric_spinor, psi))

    payload = {
        "checkpoint": "P1697_S647",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "strict_only": True,
        "legacy_bridge_used": False,
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "chain": "K_strict -> coefficients -> full L_total -> gauge+higgs block + metric+spinor block -> EOM bundle",
        "full_lagrangian_density_anchor": p1696.get("full_lagrangian_density_anchor", {}),
        "coefficient_anchor": coeff,
        "covariant_bundle": {
            "gauge_higgs_from_p1696": p1696.get("covariant_block", {}),
            "metric_spinor_block": {
                "background_marker": "R1",
                "L_metric_spinor_reduced": str(L_metric_spinor),
                "EOM_g": str(eom_g),
                "EOM_psi": str(eom_psi),
            },
        },
        "bidirectional_status": {
            "forward_export": "GAUGE_HIGGS_METRIC_SPINOR_BUNDLE_EXPORTED",
            "reverse_global_variational_origin": "KEEP_OPEN",
        },
        "status": "KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED",
        "open_obligations": [
            "full_metric_tensor_covariant_export_nonproxy",
            "full_spinor_bundle_gamma_connection_export_nonproxy",
            "counterterm_flow_closure_theorem",
            "BRST_cohomology_theorem",
            "Cutkosky_unitarity_theorem_full_sector",
            "background_family_independence_theorem",
            "strict_selector_source_or_symmetry_breaking_bridge_for_QW2191",
        ],
        "next_honest_step": "Podnieść metrykę i spinory z modelu reduced do pełnych obiektów tensorowo-spinorowych (bez proxy), a następnie skleić theorem-level pakiet QG: renormalizacja+unitarność+background-independence.",
        "lay_summary": "Dołożyliśmy kolejne fragmenty równań ruchu (metryka i spinory) do wcześniejszego bloku gauge+Higgs. To zwiększa kompletność łańcucha strict, ale wciąż brakuje pełnych dowodów kwantowej grawitacji i pełnej wersji kowariantnej bez uproszczeń.",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
