#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1704 = GEN / "p1704_s654_strict_nonproxy_metric_spinor_eom_bundle_export_scaffold_checkpoint.json"
OUT = GEN / "p1705_s655_strict_nonproxy_metric_spinor_variational_objects_export_checkpoint.json"


def main() -> None:
    p1704 = json.loads(IN1704.read_text(encoding="utf-8"))

    # Symbolic proxy objects for nonproxy-target variational export contract
    d = sp.symbols("d", positive=True, real=True)
    omega, phi0, beta, eta = sp.symbols("omega phi0 beta eta", real=True)
    K_strict = sp.cos(omega * d + phi0) / (1 + beta * d**eta)

    g = sp.Symbol("g")
    R = sp.Symbol("R")
    Ric2 = sp.Symbol("Ric2")
    Riem2 = sp.Symbol("Riem2")
    H2 = sp.Symbol("H2")
    phi2 = sp.Symbol("phi2")

    Mpl2, c1, c2, c3, xiHR, xiphiR = sp.symbols("Mpl2 c1 c2 c3 xiHR xiphiR", real=True)

    # Computable variation templates (contract-level exports, still nonproxy-open)
    dL_dR = sp.simplify(Mpl2 / 2 + 2 * c1 * R - xiHR * H2 - xiphiR * phi2)
    dL_dRic2 = sp.simplify(c2)
    dL_dRiem2 = sp.simplify(c3)

    psi = sp.Function("psi")
    psib = sp.Function("psib")
    x = sp.symbols("x", real=True)
    h = sp.Function("h")
    mpsi, ypsi = sp.symbols("mpsi ypsi", real=True)
    Lf = sp.I * (psib(x) * sp.diff(psi(x), x) - sp.diff(psib(x), x) * psi(x)) / 2 - (mpsi + ypsi * h(x)) * psib(x) * psi(x)
    eom_psi = sp.simplify(sp.diff(sp.diff(Lf, sp.diff(psi(x), x)), x) - sp.diff(Lf, psi(x)))
    eom_psib = sp.simplify(sp.diff(sp.diff(Lf, sp.diff(psib(x), x)), x) - sp.diff(Lf, psib(x)))

    payload = {
        "checkpoint": "P1705_S655",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "strict_only": True,
        "legacy_bridge_used": False,
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "chain": "K_strict -> coefficients -> full explicit L_total -> variational objects export (metric+spinor) -> EOM readiness",
        "strict_kernel_symbolic": str(K_strict),
        "source_scaffold": p1704.get("nonproxy_export_scaffold", {}),
        "metric_variational_objects_export": {
            "dL_dR_contract": str(dL_dR),
            "dL_dRicci2_contract": str(dL_dRic2),
            "dL_dRiemann2_contract": str(dL_dRiem2),
            "notes": "Contract-level symbolic objects for nonproxy metric variation assembly.",
        },
        "spinor_variational_objects_export": {
            "L_fermion_reduced_contract": str(sp.expand(Lf)),
            "EOM_psi_contract": str(eom_psi),
            "EOM_psib_contract": str(eom_psib),
            "notes": "Computable spinor-side variation templates; full curved gamma/tetrad export remains OPEN.",
        },
        "status": "KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED",
        "open_obligations": {
            "full_tensor_metric_variation_nonproxy": "OPEN",
            "full_curved_gamma_tetrad_spin_connection_export": "OPEN",
            "global_helmholtz_integrability_nonproxy": "OPEN",
            "brst_nilpotency_nonproxy_proof": "OPEN",
            "cutkosky_unitarity_full_sector": "OPEN",
            "counterterm_flow_renormalization_closure": "OPEN",
            "background_independence_family_theorem": "OPEN",
        },
        "next_honest_step": "Połączyć te obiekty kontraktowe w jeden jawny nonproxy tensor+spinor bundle i wyeksportować pełny covariant EOM set jako wejście do Helmholtz/BRST theorem proofs.",
        "lay_summary": "Przeszliśmy z samego planu do policzalnych klocków wariacyjnych dla grawitacji i fermionów. To jeszcze nie pełna finalna matematyka, ale ważny krok techniczny potrzebny do dalszych dowodów.",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
