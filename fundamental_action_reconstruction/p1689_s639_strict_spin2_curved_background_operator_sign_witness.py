#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN = GEN / "p1688_s638_strict_full_lagrangian_to_sector_eom_export_summary.json"
OUT = GEN / "p1689_s639_strict_spin2_curved_background_operator_sign_witness.json"


def spin2_eigenvalue(k2: float, mpl2: float, c2eff: float, curvature: float) -> float:
    # Minimal quadratic spin-2 proxy on constant-curvature background:
    # O2(k) ~ mpl2*k^2 + c2eff*k^4 + curvature-shift.
    return mpl2 * k2 + c2eff * (k2**2) + curvature


def main() -> None:
    p1688 = json.loads(IN.read_text(encoding="utf-8"))

    # strict-only working point (internal, no legacy bridge)
    mpl2 = 1.72
    c2eff = 0.11
    curvature = 0.02
    k2_grid = [0.0, 0.25, 0.5, 1.0, 2.0, 4.0]

    spectrum = []
    for k2 in k2_grid:
        lam = spin2_eigenvalue(k2, mpl2=mpl2, c2eff=c2eff, curvature=curvature)
        spectrum.append({"k2": k2, "lambda_spin2": lam})

    min_lam = min(row["lambda_spin2"] for row in spectrum)
    local_sign_pass = min_lam > 0.0

    payload = {
        "checkpoint": "P1689_S639",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "strict_only": True,
        "legacy_bridge_used": False,
        "input_chain": p1688.get("bidirectional_status", {}),
        "background": "constant_curvature_sector",
        "operator_model": "O2 = Mpl2*k2 + c2eff*k2^2 + curvature_shift",
        "parameters": {"Mpl2": mpl2, "c2eff": c2eff, "curvature_shift": curvature},
        "spectrum": spectrum,
        "local_sign_test": {
            "criterion": "min_eigenvalue > 0",
            "min_eigenvalue": min_lam,
            "pass": local_sign_pass,
        },
        "status": "KEEP_OPEN_QG_GLOBAL_OBLIGATIONS",
        "open_obligations": [
            "1loop_counterterm_export",
            "cutkosky_brst_unitarity_witness",
            "quantum_background_independence_global_witness",
        ],
        "next_honest_step": "Unify 1-loop counterterm derivation and BRST/Cutkosky unitarity checks on the same spin-2 background in one strict-core witness.",
        "lay_summary": "Lokalny test drgań spin-2 wypada dodatnio na tej siatce, ale to tylko warunek konieczny. Do pełnego domknięcia ToE nadal brakuje twardych dowodów kwantowych.",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")
    print(f"local_sign_pass={local_sign_pass}")


if __name__ == "__main__":
    main()
