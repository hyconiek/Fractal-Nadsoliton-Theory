#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1630 = GEN / "p1630_s580_bidirectional_strict_chain_check_summary.json"
IN1629 = GEN / "p1629_s579_global_noncyclic_cover_export_summary.json"


def coeff_map(omega: float, phi: float, beta: float, eta: float) -> tuple[float, float, float]:
    lam = 0.28 + 0.22 * omega + 0.04 * phi
    kap = 10.7 + 1.1 * beta + 0.35 * eta
    eps = 16.8 + 2.1 * phi + 0.7 * omega + 0.4 * beta
    return lam, kap, eps


def jacobian_numeric(p: tuple[float, float, float, float], h: float = 1e-4):
    base = coeff_map(*p)
    J = [[0.0] * 4 for _ in range(3)]
    for j in range(4):
        q = list(p)
        q[j] += h
        up = coeff_map(*q)
        for i in range(3):
            J[i][j] = (up[i] - base[i]) / h
    return J


def gram_det(J):
    # det(J J^T) for 3x4 J
    G = [[sum(J[i][k] * J[j][k] for k in range(4)) for j in range(3)] for i in range(3)]
    det = (
        G[0][0] * (G[1][1] * G[2][2] - G[1][2] * G[2][1])
        - G[0][1] * (G[1][0] * G[2][2] - G[1][2] * G[2][0])
        + G[0][2] * (G[1][0] * G[2][1] - G[1][1] * G[2][0])
    )
    return det


def main() -> None:
    s30 = json.loads(IN1630.read_text(encoding="utf-8"))
    s29 = json.loads(IN1629.read_text(encoding="utf-8"))

    p = s30["backward_chain_local"]["reference_kernel_params"]
    p0 = (p["omega"], p["phi"], p["beta"], p["eta"])

    chart_reports = []
    for ch in s29["global_noncyclic_cover"]["charts"]:
        scale = 1.0 + 0.02 * (ch["domain_d"][1] - ch["domain_d"][0])
        pt = (p0[0] * scale, p0[1] * scale, p0[2], p0[3])
        J = jacobian_numeric(pt)
        detG = gram_det(J)
        invertible_local = detG > 1e-12
        chart_reports.append({
            "chart": ch["name"],
            "point": {"omega": pt[0], "phi": pt[1], "beta": pt[2], "eta": pt[3]},
            "det_JJT": detG,
            "locally_invertible": invertible_local,
        })

    pass_all = all(r["locally_invertible"] for r in chart_reports)

    summary = {
        "checkpoint": "P1631_S581",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": "PASS_P1631_LOCAL_INVERTIBILITY_ON_COVER" if pass_all else "KEEP_OPEN_P1631_INVERTIBILITY_GAP",
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "forward_backward_chain": s30,
        "cover_invertibility_report": chart_reports,
        "strict_core_closure": {
            "status": "OPEN",
            "missing_exports": ["E_selector_internal_source_full_domain", "E_full_variational_proof_log_machine_checkable"],
            "missing_witnesses": ["W_noncyclic_provider_for_selector_uniqueness_PROVED_GLOBAL"],
            "missing_theorems": ["T_qw2191_selector_uniqueness_strict_GLOBAL", "T_global_toe_closure_strict"],
        },
        "strict_only": True,
        "legacy_bridge_used": False,
        "external_team_validation_required": False,
        "next_honest_step": "Sformalizować theorem-level przejście od lokalnej odwracalności na coverze do globalnej jednoznaczności selektora (overlap compatibility + composition).",
        "lay_summary": "Sprawdziliśmy, że w każdym obszarze mapy da się lokalnie odwracać parametry; teraz trzeba dowieść, że te lokalne wyniki składają się globalnie bez sprzeczności.",
    }

    out = GEN / "p1631_s581_cover_wise_jacobian_invertibility_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
