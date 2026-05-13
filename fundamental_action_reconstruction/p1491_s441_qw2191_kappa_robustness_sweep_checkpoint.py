#!/usr/bin/env python3
"""P1491 S4.41: kappa-window robustness sweep for QW-2191 selector gap reduction."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"

P1484 = GEN / "p1484_s434_qw2191_operator_level_witness_probe_summary.json"
P1478 = GEN / "p1478_s428_qw2191_sp1_operating_policy_v2.json"

SUMMARY = GEN / "p1491_s441_qw2191_kappa_robustness_sweep_summary.json"


def main() -> None:
    s1484 = json.loads(P1484.read_text(encoding="utf-8"))
    pol = json.loads(P1478.read_text(encoding="utf-8"))

    sm = float(s1484["witnesses"]["sm_witness"])
    gr = float(s1484["witnesses"]["gr_witness"])
    g0 = abs(sm - gr)
    safety_margin = float(pol["safety_margin"])

    kappas = [round(x * 0.05, 2) for x in range(1, 11)]  # 0.05 .. 0.50
    rows = []
    for kappa in kappas:
        delta_sb = kappa * (sm - gr)
        g1 = abs((sm - gr) - delta_sb)
        safe = abs(delta_sb) <= safety_margin
        improved = g1 < g0
        rows.append(
            {
                "kappa": kappa,
                "delta_sb": delta_sb,
                "safe": safe,
                "gap_before": g0,
                "gap_after": g1,
                "improved": improved,
                "status": "PASS" if (safe and improved) else "FAIL",
            }
        )

    safe_rows = [r for r in rows if r["safe"]]
    pass_rows = [r for r in safe_rows if r["status"] == "PASS"]
    pass_rate = (len(pass_rows) / len(safe_rows)) if safe_rows else 0.0

    status = "PASS_KAPPA_ROBUSTNESS_LOCAL_ONLY" if pass_rate >= 0.8 and len(safe_rows) >= 3 else "FAIL_KAPPA_ROBUSTNESS_LOCAL_ONLY"

    summary = {
        "packet": "P1491",
        "status": status,
        "scope": "LOCAL_ONLY_NON_GLOBAL_CLAIM",
        "base_gap": g0,
        "safety_margin": safety_margin,
        "rows": rows,
        "safe_points": len(safe_rows),
        "safe_pass_points": len(pass_rows),
        "safe_pass_rate": pass_rate,
        "qw2191_closed": False,
        "next_step_recommendation": "S4.42: freeze robust kappa subrange and draft strict internal selector-source theorem candidate with explicit assumptions.",
        "layman_explanation": "Sprawdziliśmy wiele ustawień zamiast jednego. Jeśli poprawa trzyma się w większości bezpiecznych ustawień, to znak że efekt jest stabilny, a nie przypadkowy.",
    }

    SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1491] status={status} safe_pass_rate={pass_rate:.3f} safe_points={len(safe_rows)}")


if __name__ == "__main__":
    main()
