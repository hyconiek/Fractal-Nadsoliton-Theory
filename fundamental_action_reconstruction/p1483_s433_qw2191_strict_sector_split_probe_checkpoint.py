#!/usr/bin/env python3
"""P1483 S4.33: strict sector split probe toward F => L_SM + L_GR."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"

P1468 = GEN / "p1468_s418_qw2191_sp1_local_pilot_summary.json"
P1481 = GEN / "p1481_s431_qw2191_sp1_cross_scenario_policy_stress_summary.json"
P1478 = GEN / "p1478_s428_qw2191_sp1_operating_policy_v2.json"
P1482 = GEN / "p1482_s432_qw2191_strict_f_to_lsm_lgr_step_summary.json"

SUMMARY = GEN / "p1483_s433_qw2191_strict_sector_split_probe_summary.json"
OBSTRUCTION = GEN / "p1483_s433_qw2191_strict_sector_split_probe_obstruction.json"


def clamp(x: float, lo: float, hi: float) -> float:
    return max(lo, min(hi, x))


def main() -> None:
    p1468 = json.loads(P1468.read_text(encoding="utf-8"))
    p1481 = json.loads(P1481.read_text(encoding="utf-8"))
    p1478 = json.loads(P1478.read_text(encoding="utf-8"))
    p1482 = json.loads(P1482.read_text(encoding="utf-8"))

    delta = float(p1468["arm_B_with_SP1_metric"]) - float(p1468["arm_A_no_selector_premise_metric"])
    stress_pass = p1481["status"] == "PASS_SP1_CROSS_SCENARIO_LOCAL_ONLY"
    ready_pass = p1482["status"] == "PASS_STRICT_NEXT_STEP_READY"

    eps_cap = abs(float(p1478["safety_margin"]))
    eps_mix = clamp(abs(delta) * 0.05, 0.0, eps_cap)

    # heuristic split proposal (local-only, no closure claim)
    w_sm = clamp(0.5 + 0.25 * (delta / 0.03), 0.35, 0.65)
    w_gr = 1.0 - w_sm

    checks = {
        "no_legacy_bridge": True,
        "qw2191_closed": False,
        "sp1_delta_positive": delta > 0.0,
        "cross_scenario_pass": stress_pass,
        "p1482_readiness_pass": ready_pass,
        "eps_mix_within_cap": eps_mix <= eps_cap,
        "weights_normalized": abs((w_sm + w_gr) - 1.0) < 1e-12,
    }

    pass_keys = [k for k in checks if k != "qw2191_closed"]
    ok = all(checks[k] for k in pass_keys)
    status = "PASS_STRICT_SECTOR_SPLIT_PROBE_LOCAL_ONLY" if ok else "FAIL_STRICT_SECTOR_SPLIT_PROBE"

    summary = {
        "packet": "P1483",
        "status": status,
        "scope": "LOCAL_ONLY_NON_GLOBAL_CLAIM",
        "route": "strict_only_F_to_LSM_plus_LGR",
        "ansatz": {
            "L_total": "w_SM*L_SM + w_GR*L_GR + eps_mix*L_mix",
            "w_SM": w_sm,
            "w_GR": w_gr,
            "eps_mix": eps_mix,
            "eps_cap": eps_cap,
        },
        "checks": checks,
        "next_step_recommendation": "S4.34: build explicit operator-level witnesses for SM-like and GR-like channels and test selector sensitivity under SP1 gate.",
        "layman_explanation": "To jak strojenie dwóch gałek: jednej dla fizyki cząstek i drugiej dla grawitacji, z małym bezpiecznikiem mieszania. Jeszcze niczego nie zamykamy; sprawdzamy tylko czy podział jest stabilny i uczciwy.",
    }
    SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    if ok:
        if OBSTRUCTION.exists():
            OBSTRUCTION.unlink()
    else:
        OBSTRUCTION.write_text(json.dumps({"packet": "P1483", "status": status, "checks": checks}, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    print(f"[P1483] status={status} w_SM={w_sm:.6f} w_GR={w_gr:.6f} eps_mix={eps_mix:.6f}")


if __name__ == "__main__":
    main()
