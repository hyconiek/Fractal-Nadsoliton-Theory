#!/usr/bin/env python3
"""P1461 S4.11: verify safety-band edges and near-outside points for h2."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
HOLDOUT = GEN / "p1452_s43_holdout_input.json"
BAND = GEN / "p1460_s410_h2_delta_safety_band_summary.json"
SUMMARY = GEN / "p1461_s411_band_edge_verification_summary.json"
OBSTRUCTION = GEN / "p1461_s411_band_edge_verification_obstruction.json"


def readj(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def eval_h2(h2: dict, delta: float, min_gain: float, replay_tol: float) -> dict:
    before = float(h2["margin_before"])
    after = float(h2["margin_after"]) + delta
    g = after - before
    gain_ok = g >= min_gain
    replay_gap = float(h2["replay_gap"])
    replay_ok = replay_gap <= replay_tol
    return {
        "delta": delta,
        "gain": g,
        "gain_ok": gain_ok,
        "replay_gap": replay_gap,
        "replay_ok": replay_ok,
        "status": "PASS" if gain_ok and replay_ok else "FAIL",
    }


def main() -> None:
    holdout = readj(HOLDOUT)
    band = readj(BAND)["h2_safety_band"]
    min_gain = float(holdout["min_gain"])
    replay_tol = float(holdout["replay_tol"])
    h2 = next(c for c in holdout["cases"] if c["id"] == "h2")

    dmin = float(band["delta_min"])
    dmax = float(band["delta_max"])
    points = [
        {"label": "edge_min", "delta": dmin},
        {"label": "edge_max", "delta": dmax},
        {"label": "outside_below", "delta": round(dmin - 0.0001, 6)},
        {"label": "outside_above", "delta": round(dmax + 0.0001, 6)},
    ]

    rows = []
    first_fail = None
    for point in points:
        res = eval_h2(h2, point["delta"], min_gain, replay_tol)
        row = {"label": point["label"], **res}
        rows.append(row)
        if row["status"] == "FAIL" and first_fail is None:
            first_fail = row

    expected = {
        "edge_min": "PASS",
        "edge_max": "PASS",
        "outside_below": "FAIL",
        "outside_above": "PASS",
    }
    expectation_mismatch = [r for r in rows if expected[r["label"]] != r["status"]]
    status = "PASS_EDGE_VERIFICATION_LOCAL_ONLY" if not expectation_mismatch else "FAIL_EDGE_VERIFICATION"

    summary = {
        "packet": "P1461",
        "status": status,
        "scope": "LOCAL_ONLY_NON_GLOBAL_CLAIM",
        "strict_core_qw2191_closed": False,
        "legacy_bridge_used": False,
        "safety_band": {"delta_min": dmin, "delta_max": dmax},
        "expected_pattern": expected,
        "results": rows,
        "expectation_mismatch_count": len(expectation_mismatch),
    }
    SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    if expectation_mismatch:
        obstruction = {
            "packet": "P1461",
            "status": "FAIL_EDGE_PATTERN_MISMATCH",
            "first_fail": first_fail,
            "mismatch": expectation_mismatch,
        }
        OBSTRUCTION.write_text(json.dumps(obstruction, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        print(f"[P1461] status={status} mismatch={len(expectation_mismatch)}")
    else:
        if OBSTRUCTION.exists():
            OBSTRUCTION.unlink()
        print(f"[P1461] status={status} results={len(rows)}")


if __name__ == "__main__":
    main()
