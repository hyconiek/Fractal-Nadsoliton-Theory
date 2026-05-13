#!/usr/bin/env python3
"""P1460 S4.10: export h2 safety-band and verify h1/h3 non-degradation in local replay."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
P1452 = GEN / "p1452_s43_holdout_input.json"
P1458 = GEN / "p1458_s48_h2_robustness_sweep_summary.json"
P1459 = GEN / "p1459_s49_h2_stress_edge_summary.json"
SUMMARY = GEN / "p1460_s410_h2_delta_safety_band_summary.json"


def readj(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def gain(case: dict) -> float:
    return float(case["margin_after"]) - float(case["margin_before"])


def main() -> None:
    h1452 = readj(P1452)
    h1458 = readj(P1458)
    h1459 = readj(P1459)

    lower = float(h1459["last_pass"]["delta"])
    upper = max(float(v["delta_margin_boost_h2"]) for v in h1458["variants"] if bool(v["all_pass"]))

    min_gain = float(h1452["min_gain"])
    replay_tol = float(h1452["replay_tol"])

    h1 = next(c for c in h1452["cases"] if c["id"] == "h1")
    h3 = next(c for c in h1452["cases"] if c["id"] == "h3")

    controls = []
    for case in (h1, h3):
        g = gain(case)
        controls.append(
            {
                "id": case["id"],
                "gain": g,
                "gain_ok": g >= min_gain,
                "replay_gap": float(case["replay_gap"]),
                "replay_ok": float(case["replay_gap"]) <= replay_tol,
                "status": "PASS" if (g >= min_gain and float(case["replay_gap"]) <= replay_tol) else "FAIL",
            }
        )

    control_fail = next((c for c in controls if c["status"] == "FAIL"), None)
    status = "PASS_SAFETY_BAND_LOCAL_ONLY" if control_fail is None and lower <= upper else "FAIL_SAFETY_BAND_LOCAL"

    summary = {
        "packet": "P1460",
        "status": status,
        "scope": "LOCAL_ONLY_NON_GLOBAL_CLAIM",
        "strict_core_qw2191_closed": False,
        "legacy_bridge_used": False,
        "h2_safety_band": {"delta_min": lower, "delta_max": upper},
        "source": {
            "lower_from_p1459_last_pass": lower,
            "upper_from_p1458_max_pass_variant": upper,
        },
        "control_replay_h1_h3": controls,
    }

    SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1460] status={status} band=[{lower:.4f},{upper:.4f}] controls={len(controls)}")


if __name__ == "__main__":
    main()
