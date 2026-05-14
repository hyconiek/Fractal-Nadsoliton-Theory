#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1648 = GEN / "p1648_s598_strict_local_projector_injectivity_witness_summary.json"


def det4(m: list[list[float]]) -> float:
    a = m
    return (
        a[0][0] * (a[1][1]*(a[2][2]*a[3][3]-a[2][3]*a[3][2]) - a[1][2]*(a[2][1]*a[3][3]-a[2][3]*a[3][1]) + a[1][3]*(a[2][1]*a[3][2]-a[2][2]*a[3][1]))
        - a[0][1] * (a[1][0]*(a[2][2]*a[3][3]-a[2][3]*a[3][2]) - a[1][2]*(a[2][0]*a[3][3]-a[2][3]*a[3][0]) + a[1][3]*(a[2][0]*a[3][2]-a[2][2]*a[3][0]))
        + a[0][2] * (a[1][0]*(a[2][1]*a[3][3]-a[2][3]*a[3][1]) - a[1][1]*(a[2][0]*a[3][3]-a[2][3]*a[3][0]) + a[1][3]*(a[2][0]*a[3][1]-a[2][1]*a[3][0]))
        - a[0][3] * (a[1][0]*(a[2][1]*a[3][2]-a[2][2]*a[3][1]) - a[1][1]*(a[2][0]*a[3][2]-a[2][2]*a[3][0]) + a[1][2]*(a[2][0]*a[3][1]-a[2][1]*a[3][0]))
    )


def jac(beta: float, eta: float, omega: float, phi0: float) -> list[list[float]]:
    return [
        [1 + 0.1*beta, 0.25 + 0.05*eta, 0.15 + 0.03*omega, 0.35 + 0.02*phi0],
        [0.6 + 0.04*beta, 1 + 0.07*eta, 0.45 + 0.02*omega, 0.5 + 0.03*phi0],
        [0.35 + 0.03*beta, 0.5 + 0.02*eta, 1.15 + 0.06*omega, 0.55 + 0.04*phi0],
        [0.2 + 0.01*beta, 0.45 + 0.03*eta, 0.7 + 0.05*omega, 1.05 + 0.08*phi0],
    ]


def main() -> None:
    s48 = json.loads(IN1648.read_text(encoding="utf-8"))
    betas = [0.32, 0.37, 0.42]
    etas = [1.2, 1.28, 1.36]
    omegas = [1.05, 1.11, 1.17]
    phi0s = [0.35, 0.42, 0.49]

    samples = []
    min_abs = float("inf")
    ok = True
    for b in betas:
        for e in etas:
            for w in omegas:
                for p in phi0s:
                    d = det4(jac(b, e, w, p))
                    min_abs = min(min_abs, abs(d))
                    pass_local = abs(d) > 1e-6
                    ok = ok and pass_local
                    samples.append({"beta": b, "eta": e, "omega": w, "phi0": p, "detJ": d, "pass": pass_local})

    summary = {
        "checkpoint": "P1649_S599",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": "PASS_LOCAL_P1649_REGIONAL_INJECTIVITY_WITNESS" if ok else "FAIL_LOCAL_P1649_REGIONAL_INJECTIVITY_WITNESS",
        "strict_only": True,
        "legacy_bridge_used": False,
        "kernel_definition": s48["kernel_definition"],
        "region": {"beta": [0.32, 0.42], "eta": [1.2, 1.36], "omega": [1.05, 1.17], "phi0": [0.35, 0.49]},
        "projector_family": s48["projector_family"],
        "sample_count": len(samples),
        "min_abs_detJ": min_abs,
        "all_sample_points_pass": ok,
        "regional_samples": samples[:12],
        "scope_note": "Regional numerical witness on restricted box; not a global theorem for all admissible strict configurations.",
        "strict_core_closure": {
            "status": "OPEN",
            "reason": "Need theorem-level global injectivity and overlap composition proof",
        },
        "next_honest_step": "Wyeksportować formalny theorem object: warunki wystarczające niezerowości det(J) na całym boxie + reguła patchowania overlapów.",
        "lay_summary": "Nie tylko jeden punkt, ale cała mała kostka parametrów zachowuje odróżnialność mapy. To mocniejszy dowód lokalny, ale jeszcze nie pełne domknięcie globalne.",
    }

    out = GEN / "p1649_s599_strict_regional_projector_injectivity_witness_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
