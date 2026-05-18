#!/usr/bin/env python3
"""P2003 S953 strict Cutkosky classifier robustness-band witness.

Next honest step after P2002: robustness check for missing-channel-vs-structural
classification by scanning uncertainty bands over gx weight and kernel amplitude.
"""
from __future__ import annotations
import json, platform
from pathlib import Path
from typing import Any
import numpy as np
import scipy.linalg as la

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2003_s953_strict_cutkosky_classifier_robustness_band_witness.json"
TS = "2026-05-18T00:00:00+00:00"


def load(name: str) -> dict[str, Any]:
    p = GEN / name
    if not p.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(p.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2001 = load("p2001_s951_strict_cutkosky_full_three_loop_kernel_plus_extra_channel_witness.json")
    p2002 = load("p2002_s952_strict_cutkosky_missing_channel_vs_structural_classifier_witness.json")

    rows = p2001.get("grid_table", [])
    deltas = np.array([float(r.get("Delta_opt", 0.0)) for r in rows], dtype=float)
    gx_vals = np.array([float(r.get("Cut_channels", {}).get("gx", 0.0)) for r in rows], dtype=float)
    base_l2 = float(la.norm(deltas, 2)) if deltas.size else float("nan")

    # Uncertainty scan: gx weight ±40%, gx kernel amplitude ±30%.
    weight_scales = [0.6, 0.8, 1.0, 1.2, 1.4]
    amp_scales = [0.7, 0.85, 1.0, 1.15, 1.3]

    scenarios = []
    l2_values = []
    base_classifier = p2002.get("classifier", "MIXED_OR_INCONCLUSIVE")

    # Reconstruct proxy effect on delta when gx term is scaled by factor f
    # delta_new = delta_old - (f-1)*gx_term
    for ws in weight_scales:
        for a in amp_scales:
            f = ws * a
            d_new = deltas - (f - 1.0) * gx_vals
            l2_new = float(la.norm(d_new, 2))
            l2_values.append(l2_new)
            ratio = l2_new / base_l2 if base_l2 > 0 else float("inf")
            if ratio < 0.95:
                cls = "MISSING_CHANNEL_PRESSURE_SUPPORTED"
            elif ratio > 1.05:
                cls = "STRUCTURAL_OBSTRUCTION_PRESSURE_SUPPORTED"
            else:
                cls = "MIXED_OR_INCONCLUSIVE"
            scenarios.append({
                "weight_scale": ws,
                "amplitude_scale": a,
                "combined_scale": f,
                "l2_ratio_vs_p2001": ratio,
                "classifier": cls,
            })

    l2_min = float(np.min(l2_values)) if l2_values else float("nan")
    l2_max = float(np.max(l2_values)) if l2_values else float("nan")

    counts = {
        "MISSING_CHANNEL_PRESSURE_SUPPORTED": sum(1 for s in scenarios if s["classifier"] == "MISSING_CHANNEL_PRESSURE_SUPPORTED"),
        "STRUCTURAL_OBSTRUCTION_PRESSURE_SUPPORTED": sum(1 for s in scenarios if s["classifier"] == "STRUCTURAL_OBSTRUCTION_PRESSURE_SUPPORTED"),
        "MIXED_OR_INCONCLUSIVE": sum(1 for s in scenarios if s["classifier"] == "MIXED_OR_INCONCLUSIVE"),
    }
    total = max(1, len(scenarios))
    robustness_score = counts[base_classifier] / total if base_classifier in counts else 0.0

    robust = robustness_score >= 0.7

    gate = {
        "p2001_present": p2001.get("result_kind") == "PASS_FULL_THREE_LOOP_PLUS_EXTRA_CHANNEL_WITNESS",
        "p2002_present": p2002.get("result_kind") == "PASS_DELTA_NORM_CLASSIFIER_WITNESS",
        "nonempty_scan": len(scenarios) > 0,
        "finite_band": bool(np.isfinite(l2_min) and np.isfinite(l2_max)),
        "robustness_score_bounded": 0.0 <= robustness_score <= 1.0,
    }

    out = {
        "ledger_id": "P2003_S953_STRICT_CUTKOSKY_CLASSIFIER_ROBUSTNESS_BAND_WITNESS",
        "packet_id": "P2003",
        "stage_id": "S953",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TS,
        "route": "strict_only",
        "depends_on": {
            "p2001": p2001.get("ledger_id"),
            "p2002": p2002.get("ledger_id"),
        },
        "base_classifier_from_p2002": base_classifier,
        "uncertainty_band": {
            "gx_weight_scale": [min(weight_scales), max(weight_scales)],
            "gx_amplitude_scale": [min(amp_scales), max(amp_scales)],
            "l2_band": [l2_min, l2_max],
        },
        "classifier_counts": counts,
        "robustness_score": robustness_score,
        "robust_under_band": robust,
        "scan_table": scenarios,
        "gatekeeper_checks": gate,
        "result_kind": "PASS_CLASSIFIER_ROBUSTNESS_BAND_WITNESS" if all(gate.values()) else "OPEN_OBSTRUCTION_WITH_TRACE",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "false_pass_guard": "Robustness scan is sensitivity diagnostics only; not theorem-grade unitarity closure.",
        "next_honest_step": "Replace gx proxy with explicit loop-derived amplitude and rerun robustness band with channelwise uncertainty propagated from backend loop data.",
        "lay_explanation": "Sprawdziliśmy, czy wniosek z P2002 utrzyma się przy rozsądnych zmianach niepewności nowego kanału. To test stabilności wniosku, nie końcowy dowód.",
        "environment": {"python": platform.python_version(), "numpy": np.__version__, "scipy": __import__('scipy').__version__},
    }

    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P2003] wrote witness: {OUT}")


if __name__ == "__main__":
    main()
