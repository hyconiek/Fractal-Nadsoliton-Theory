from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

TOP_POINTS = [
    (0.36, 0.08, 1.15, 1.55),
    (0.36, 0.08, 1.05, 1.55),
    (0.36, 0.08, 0.95, 1.55),
    (0.36, 0.08, 0.85, 1.55),
    (0.36, 0.08, 1.15, 1.70),
]
NOISE_LEVELS = [1e-4, 3e-4, 1e-3]
TOL_SM = 0.10
TOL_GR = 0.08


def coeffs(o, ph, b, e):
    lam = abs(math.cos(ph)) / (1.0 + 0.3 * b)
    kap = (o * o + 0.5 * e) / (0.2 + abs(math.cos(ph)))
    eps = (b * e + 0.1 * o) / (1.0 + abs(math.sin(ph)))
    return lam, kap, eps


def main() -> None:
    root = Path(__file__).resolve().parent
    out_dir = root / "generated"
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    worst = {"R_bundle": -1.0}
    for p in TOP_POINTS:
        lam0, kap0, eps0 = coeffs(*p)
        for nl in NOISE_LEVELS:
            lam1 = lam0 * (1.0 + 0.35 * nl / 1e-3)
            kap1 = kap0 * (1.0 + 0.95 * nl / 1e-3)
            drift_sm = abs(lam1 - lam0) / abs(lam0)
            drift_gr = abs(kap1 - kap0) / abs(kap0)
            r_bundle = max(drift_sm / TOL_SM, drift_gr / TOL_GR)
            row = {
                "point": p,
                "noise": nl,
                "drift_sm": drift_sm,
                "drift_gr": drift_gr,
                "R_bundle": r_bundle,
            }
            rows.append(row)
            if r_bundle > worst["R_bundle"]:
                worst = row

    status = "PASS_W1576A_CANDIDATE" if worst["R_bundle"] <= 1.0 else "FAIL_W1576A_CANDIDATE"

    summary = {
        "checkpoint": "P1577_S527",
        "date_utc": "2026-05-14",
        "status": status,
        "strict_only": True,
        "legacy_bridge_used": False,
        "qw2191_closed": False,
        "toe_closed": False,
        "num_tested_points": len(TOP_POINTS),
        "noise_levels": NOISE_LEVELS,
        "tolerances": {"sm": TOL_SM, "gr": TOL_GR},
        "worst_case": worst,
        "missing_exports_witnesses_theorems": [
            "T1577A_gr_channel_drift_reduction_theorem",
            "T1577B_joint_sm_gr_stabilization_theorem",
            "T1577C_final_strict_core_bundle_closure_theorem",
        ],
        "next_honest_step": "P1578_apply_local_gr_channel_renormalization_and_rerun_witness",
        "lay_summary": "Kanał GR dominuje wrażliwość na szum i blokuje pełny witness odporności bundle'u.",
    }
    summary["audit_digest"] = hashlib.sha256(json.dumps(summary, sort_keys=True).encode("utf-8")).hexdigest()

    out = out_dir / "p1577_s527_strict_sm_gr_bundle_robustness_witness_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1577] wrote {out} status={status}")


if __name__ == "__main__":
    main()
