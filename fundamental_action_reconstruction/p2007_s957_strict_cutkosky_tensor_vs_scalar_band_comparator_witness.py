#!/usr/bin/env python3
"""P2007 S957 strict Cutkosky tensor-vs-scalar band comparator witness.

Next honest step after P2006: quantitatively compare P2005 (scalar backend band)
vs P2006 (tensor-surrogate covariance band) to measure whether tensorization
stabilizes classifier behavior.
"""
from __future__ import annotations
import json, platform
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2007_s957_strict_cutkosky_tensor_vs_scalar_band_comparator_witness.json"
TS = "2026-05-18T00:00:00+00:00"


def load(name: str) -> dict[str, Any]:
    p = GEN / name
    if not p.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(p.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2005 = load("p2005_s955_strict_cutkosky_gx_backend_amplitude_bound_and_classifier_witness.json")
    p2006 = load("p2006_s956_strict_cutkosky_gx_backend_tensor_surrogate_and_covariance_classifier_witness.json")

    s_scan = p2005.get("classifier_scan", {})
    t_scan = p2006.get("classifier_scan", {})

    s_dom = s_scan.get("dominant_classifier", "MIXED_OR_INCONCLUSIVE")
    t_dom = t_scan.get("dominant_classifier", "MIXED_OR_INCONCLUSIVE")
    s_share = float(s_scan.get("dominant_share", 0.0))
    t_share = float(t_scan.get("dominant_share", 0.0))

    s_count = int(s_scan.get("count", 0))
    t_count = int(t_scan.get("count", 0))

    s_counts = s_scan.get("counts", {})
    t_counts = t_scan.get("counts", {})
    s_missing = int(s_counts.get("MISSING_CHANNEL_PRESSURE_SUPPORTED", 0))
    t_missing = int(t_counts.get("MISSING_CHANNEL_PRESSURE_SUPPORTED", 0))

    share_shift = t_share - s_share
    missing_frac_s = s_missing / s_count if s_count > 0 else 0.0
    missing_frac_t = t_missing / t_count if t_count > 0 else 0.0
    missing_frac_shift = missing_frac_t - missing_frac_s

    if share_shift > 0.05 and missing_frac_shift > 0.05:
        verdict = "TENSORIZATION_IMPROVES_MISSING_CHANNEL_SIGNAL"
    elif share_shift < -0.05 and missing_frac_shift < -0.05:
        verdict = "TENSORIZATION_WEAKENS_MISSING_CHANNEL_SIGNAL"
    else:
        verdict = "TENSORIZATION_NEUTRAL_OR_MIXED"

    gate = {
        "p2005_present": p2005.get("result_kind") == "PASS_GX_BACKEND_BOUND_CLASSIFIER_WITNESS",
        "p2006_present": p2006.get("result_kind") == "PASS_GX_TENSOR_COVARIANCE_CLASSIFIER_WITNESS",
        "scalar_scan_nonempty": s_count > 0,
        "tensor_scan_nonempty": t_count > 0,
        "shares_bounded": 0.0 <= s_share <= 1.0 and 0.0 <= t_share <= 1.0,
    }

    out = {
        "ledger_id": "P2007_S957_STRICT_CUTKOSKY_TENSOR_VS_SCALAR_BAND_COMPARATOR_WITNESS",
        "packet_id": "P2007",
        "stage_id": "S957",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TS,
        "route": "strict_only",
        "depends_on": {"p2005": p2005.get("ledger_id"), "p2006": p2006.get("ledger_id")},
        "scalar_band_summary": {
            "dominant_classifier": s_dom,
            "dominant_share": s_share,
            "scan_count": s_count,
            "missing_fraction": missing_frac_s,
        },
        "tensor_band_summary": {
            "dominant_classifier": t_dom,
            "dominant_share": t_share,
            "scan_count": t_count,
            "missing_fraction": missing_frac_t,
        },
        "comparison": {
            "dominant_share_shift_tensor_minus_scalar": share_shift,
            "missing_fraction_shift_tensor_minus_scalar": missing_frac_shift,
            "verdict": verdict,
        },
        "gatekeeper_checks": gate,
        "result_kind": "PASS_TENSOR_VS_SCALAR_COMPARATOR_WITNESS" if all(gate.values()) else "OPEN_OBSTRUCTION_WITH_TRACE",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "false_pass_guard": "Comparator summarizes diagnostics only; no theorem-grade closure claim.",
        "next_honest_step": "Export explicit backend loop-amplitude tensor object and rerun comparator with direct tensor source (without surrogate).",
        "lay_explanation": "Porównaliśmy dwa podejścia: prosty zakres skali i zakres z macierzy backendu. Sprawdzamy, czy wersja tensorowa daje stabilniejszy sygnał klasyfikacji.",
        "environment": {"python": platform.python_version()},
    }

    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P2007] wrote witness: {OUT}")


if __name__ == "__main__":
    main()
