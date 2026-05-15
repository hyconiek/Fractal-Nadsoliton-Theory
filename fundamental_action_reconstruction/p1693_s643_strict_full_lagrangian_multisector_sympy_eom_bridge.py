#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1691 = GEN / "p1691_s641_strict_full_chain_lagrangian_to_qg_theorem_obligation_matrix.json"
IN1664 = GEN / "p1664_s614_strict_full_lagrangian_manifest_and_inversion.json"
OUT = GEN / "p1693_s643_strict_full_lagrangian_multisector_sympy_eom_bridge.json"


def main() -> None:
    p1691 = json.loads(IN1691.read_text(encoding="utf-8"))
    p1664 = json.loads(IN1664.read_text(encoding="utf-8"))
    c = p1664["coefficient_map"]

    x = sp.symbols("x", real=True)
    phi = sp.Function("phi")(x)
    h = sp.Function("h")(x)
    A = sp.Function("A")(x)
    psi = sp.Function("psi")(x)

    # strict-derived coefficient anchors (no legacy bridge)
    mphi2 = sp.Float(c["muH2"]) * (1 + sp.Float(p1664["kernel_input"]["beta"]))
    lam3 = sp.Float(c["cR2"])
    lam4 = sp.Float(c["cRic2"])
    lamh = sp.Float(c["lambdaH"])
    lamph = sp.Float(c["xiHR"]) * lamh
    mH2 = sp.Float(c["muH2"])
    ZA = sp.Float(c["Z1"])
    y = sp.Float(c["ye"])

    # reduced multisector strict density approximating L_SM + L_GR carrier pieces
    L_phi = sp.Rational(1, 2) * sp.diff(phi, x) ** 2 - mphi2 * phi**2 / 2 - lam3 * phi**3 / 6 - lam4 * phi**4 / 24
    L_h = sp.Rational(1, 2) * sp.diff(h, x) ** 2 - mH2 * h**2 / 2 - lamh * h**4 / 4
    L_mix = -lamph * phi**2 * h**2 / 2
    L_gauge = ZA * sp.diff(A, x) ** 2 / 2
    L_yuk = -y * phi * psi**2
    L_total_reduced = sp.expand(L_phi + L_h + L_mix + L_gauge + L_yuk)

    eom_phi = sp.simplify(sp.diff(sp.diff(L_total_reduced, sp.diff(phi, x)), x) - sp.diff(L_total_reduced, phi))
    eom_h = sp.simplify(sp.diff(sp.diff(L_total_reduced, sp.diff(h, x)), x) - sp.diff(L_total_reduced, h))
    eom_A = sp.simplify(sp.diff(sp.diff(L_total_reduced, sp.diff(A, x)), x) - sp.diff(L_total_reduced, A))

    # fermion proxy in reduced scalarized form
    eom_psi_proxy = sp.simplify(sp.diff(L_total_reduced, psi))

    payload = {
        "checkpoint": "P1693_S643",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "strict_only": True,
        "legacy_bridge_used": False,
        "forward_chain": "K_strict -> coefficients -> full nonskeleton L_total anchor -> multisector symbolic replay -> EOM witnesses",
        "full_lagrangian_density_anchor": p1691.get("full_lagrangian_density", {}),
        "reduced_multisector_density": str(L_total_reduced),
        "symbolic_eom_multisector": {
            "phi": str(eom_phi),
            "h": str(eom_h),
            "A": str(eom_A),
            "psi_proxy": str(eom_psi_proxy),
        },
        "theory_chain_bidirectional_note": "Forward replay exported; reverse global variational+QG theorem closure remains open.",
        "status": "KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED",
        "open_obligations": [
            "covariant_spinor_gauge_EL_export_nonproxy",
            "full_counterterm_1loop_flow_closure_theorem",
            "BRST_nilpotency_cohomology_theorem_full_sector",
            "Cutkosky_unitarity_theorem_spin2_plus_SM_mix",
            "background_family_independence_quantum_theorem",
        ],
        "next_honest_step": "Zastąpić proxy fermionowy pełnym kowariantnym operatorem spinorowym i dołączyć BRST/Cutkosky theorem witnesses na rodzinie teł.",
        "lay_summary": "Mamy teraz automatycznie policzony szkic równań ruchu dla wielu sektorów naraz, ale pełny dowód bezpieczeństwa kwantowej grawitacji nadal jest otwarty.",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
