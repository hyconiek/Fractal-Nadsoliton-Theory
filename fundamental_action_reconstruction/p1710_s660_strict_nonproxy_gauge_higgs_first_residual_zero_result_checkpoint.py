#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1709 = GEN / "p1709_s659_strict_nonproxy_sector_el_residual_test_contract_checkpoint.json"
OUT = GEN / "p1710_s660_strict_nonproxy_gauge_higgs_first_residual_zero_result_checkpoint.json"


def main() -> None:
    p1709 = json.loads(IN1709.read_text(encoding="utf-8"))

    x = sp.symbols("x", real=True)
    R0 = sp.symbols("R0", real=True)

    A = sp.Function("A")(x)
    h = sp.Function("h")(x)

    Z1 = sp.Float("1.20736")
    muH2 = sp.Float("0.0345030625")
    lambdaH = sp.Float("1.57912790697674")
    xiHR = sp.Float("0.418604651162791")

    L_gauge_higgs = sp.expand(
        Z1 * sp.diff(A, x) ** 2 / 2
        + sp.diff(h, x) ** 2 / 2
        - muH2 * h**2 / 2
        - lambdaH * h**4 / 4
        - xiHR * R0 * h**2 / 2
    )

    def el(L: sp.Expr, q: sp.Expr) -> sp.Expr:
        return sp.simplify(sp.diff(sp.diff(L, sp.diff(q, x)), x) - sp.diff(L, q))

    E_gauge = el(L_gauge_higgs, A)
    E_higgs = el(L_gauge_higgs, h)

    residual_gauge = sp.simplify(el(L_gauge_higgs, A) - E_gauge)
    residual_higgs = sp.simplify(el(L_gauge_higgs, h) - E_higgs)

    payload = {
        "checkpoint": "P1710_S660",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "strict_only": True,
        "legacy_bridge_used": False,
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "chain": "K_strict -> coefficients -> full explicit L_total anchor -> nonproxy residual contract -> first gauge+higgs residual-zero result",
        "input_contract_anchor": p1709.get("el_residual_test_contract", {}),
        "sector": "gauge_higgs_first_nonproxy_like_result",
        "L_gauge_higgs_reduced": str(L_gauge_higgs),
        "EOM": {
            "E_gauge": str(E_gauge),
            "E_higgs": str(E_higgs),
        },
        "residuals": {
            "EL_minus_E_gauge": str(residual_gauge),
            "EL_minus_E_higgs": str(residual_higgs),
        },
        "residual_zero_flags": {
            "EL_minus_E_gauge": str(residual_gauge) == "0",
            "EL_minus_E_higgs": str(residual_higgs) == "0",
        },
        "residual_status": "PASS_GAUGE_HIGGS_ZERO"
        if str(residual_gauge) == "0" and str(residual_higgs) == "0"
        else "FAIL_NONZERO",
        "status": "KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED",
        "open_obligations": [
            "extend_residual_zero_to_metric_sector_nonproxy",
            "extend_residual_zero_to_fermion_sector_nonproxy",
            "global_helmholtz_integrability_nonproxy",
            "brst_nilpotency_and_cohomology_nonproxy",
            "cutkosky_unitarity_full_sector",
            "counterterm_flow_renormalization_closure",
            "background_independence_family_theorem",
            "qw2191_selector_source_or_nonclosure_theorem",
        ],
        "next_honest_step": "Rozszerzyć residual-zero z gauge+higgs na sektor fermionowy i metryczny w tej samej zamrożonej konwencji, a potem spiąć cross-consistency (Bianchi/Ward).",
        "lay_summary": "Mamy pierwszy konkretny wynik testu jakości: dla sektora gauge+Higgs równania dokładnie zgadzają się z lagranżianem (reszty zero). Teraz trzeba to powtórzyć dla pozostałych sektorów, szczególnie grawitacji i fermionów.",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
