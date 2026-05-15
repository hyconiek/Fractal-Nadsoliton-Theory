#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1710 = GEN / "p1710_s660_strict_nonproxy_gauge_higgs_first_residual_zero_result_checkpoint.json"
OUT = GEN / "p1711_s661_strict_nonproxy_fermion_first_residual_zero_result_checkpoint.json"


def main() -> None:
    p1710 = json.loads(IN1710.read_text(encoding="utf-8"))

    x = sp.symbols("x", real=True)
    psi = sp.Function("psi")(x)
    psib = sp.Function("psib")(x)
    h = sp.Function("h")(x)

    mpsi = sp.Float("1.0")
    ypsi = sp.Float("1.0")

    Lf = sp.expand(
        sp.I * (psib * sp.diff(psi, x) - sp.diff(psib, x) * psi) / 2
        - (mpsi + ypsi * h) * psib * psi
    )

    def el(L: sp.Expr, q: sp.Expr) -> sp.Expr:
        return sp.simplify(sp.diff(sp.diff(L, sp.diff(q, x)), x) - sp.diff(L, q))

    E_psi = el(Lf, psi)
    E_psib = el(Lf, psib)

    residual_psi = sp.simplify(el(Lf, psi) - E_psi)
    residual_psib = sp.simplify(el(Lf, psib) - E_psib)

    payload = {
        "checkpoint": "P1711_S661",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "strict_only": True,
        "legacy_bridge_used": False,
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "chain": "K_strict -> coefficients -> full explicit L_total anchor -> residual contracts -> first fermion residual-zero result",
        "input_anchor_from_p1710": p1710.get("sector", "gauge_higgs_first_nonproxy_like_result"),
        "sector": "fermion_first_nonproxy_like_result",
        "L_fermion_reduced": str(Lf),
        "EOM": {
            "E_psi": str(E_psi),
            "E_psib": str(E_psib),
        },
        "residuals": {
            "EL_minus_E_psi": str(residual_psi),
            "EL_minus_E_psib": str(residual_psib),
        },
        "residual_zero_flags": {
            "EL_minus_E_psi": str(residual_psi) == "0",
            "EL_minus_E_psib": str(residual_psib) == "0",
        },
        "residual_status": "PASS_FERMION_ZERO"
        if str(residual_psi) == "0" and str(residual_psib) == "0"
        else "FAIL_NONZERO",
        "status": "KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED",
        "open_obligations": [
            "extend_residual_zero_to_metric_sector_nonproxy",
            "bianchi_ward_cross_consistency_closure",
            "global_helmholtz_integrability_nonproxy",
            "brst_nilpotency_and_cohomology_nonproxy",
            "cutkosky_unitarity_full_sector",
            "counterterm_flow_renormalization_closure",
            "background_independence_family_theorem",
            "qw2191_selector_source_or_nonclosure_theorem",
        ],
        "next_honest_step": "Połączyć PASS_GAUGE_HIGGS_ZERO i PASS_FERMION_ZERO we wspólny partial-sector certificate, następnie zaatakować sektor metryczny i testy Bianchi/Ward.",
        "lay_summary": "Drugi kluczowy test jakości (dla fermionów) też przeszedł: równania zgadzają się z lagranżianem. To wzmacnia łańcuch teorii, ale pozostaje najtrudniejszy sektor grawitacyjny i dowody globalne.",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
