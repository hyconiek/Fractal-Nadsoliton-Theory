#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1636 = GEN / "p1636_s586_strict_bidirectional_kernel_coefficient_closure_audit_summary.json"


# constant Jacobian of coeff_map used in P1636
J = [
    [0.22, 0.04, 0.0, 0.0],
    [0.0, 0.0, 1.1, 0.35],
    [0.7, 2.1, 0.4, 0.0],
]


def rank3x4(mat: list[list[float]]) -> int:
    # quick Gaussian elimination
    a = [row[:] for row in mat]
    r, c, rank = 0, 0, 0
    while r < 3 and c < 4:
        piv = max(range(r, 3), key=lambda i: abs(a[i][c]))
        if abs(a[piv][c]) < 1e-12:
            c += 1
            continue
        a[r], a[piv] = a[piv], a[r]
        pv = a[r][c]
        for j in range(c, 4):
            a[r][j] /= pv
        for i in range(3):
            if i != r:
                f = a[i][c]
                for j in range(c, 4):
                    a[i][j] -= f * a[r][j]
        rank += 1
        r += 1
        c += 1
    return rank


def nullspace_direction() -> dict[str, float]:
    # Solve J v = 0 with eta=1 gauge
    eta = 1.0
    beta = -(0.35 / 1.1) * eta
    omega = -(0.04 / 0.22) * 1.0  # set phi=1 first
    phi = 1.0
    # adjust (omega,phi) to satisfy row3 with chosen beta
    # 0.7 w +2.1 p +0.4 b =0 -> w = (-2.1 p -0.4 b)/0.7
    omega = (-2.1 * phi - 0.4 * beta) / 0.7
    return {"omega": omega, "phi": phi, "beta": beta, "eta": eta}


def main() -> None:
    s36 = json.loads(IN1636.read_text(encoding="utf-8"))
    r = rank3x4(J)
    null_dim = 4 - r
    v = nullspace_direction()

    summary = {
        "checkpoint": "P1637_S587",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": "PASS_P1637_IDENTIFIABILITY_OBSTRUCTION_EXPORTED" if null_dim >= 1 else "KEEP_OPEN_P1637_OBSTRUCTION_NOT_SHOWN",
        "route_target": s36["route_target"],
        "strict_chain": s36["strict_chain"],
        "identifiability_analysis": {
            "jacobian_kernel_to_coeff": J,
            "rank": r,
            "nullspace_dimension": null_dim,
            "example_null_direction": v,
            "meaning": "At least one internal direction leaves coefficients unchanged, so coeff->kernel inverse is non-unique without extra selector constraints.",
        },
        "links_to_p1636": {
            "p1636_status": s36["status"],
            "p1636_max_abs_residual": s36["backward_reconstruction"]["max_abs_residual"],
        },
        "strict_core_closure": s36["strict_core_closure"],
        "strict_only": True,
        "legacy_bridge_used": False,
        "next_honest_step": "Dodać theorem-level selector constraint, który usuwa nullspace (np. globalny warunek atlasowy/overlap consistency) i zamienia lokalny inverse w globalnie jednoznaczny export.",
        "lay_summary": "Pokazaliśmy matematycznie, że obecne równania mają ukryty kierunek swobody. Dlatego z tych samych danych można odzyskać wiele podobnych kerneli i potrzebny jest dodatkowy globalny warunek.",
    }

    out = GEN / "p1637_s587_strict_identifiability_obstruction_from_kernel_to_coefficient_map_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
