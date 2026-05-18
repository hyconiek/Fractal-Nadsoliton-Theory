#!/usr/bin/env python3
"""P2002 S952 strict Cutkosky missing-channel vs structural classifier witness.

Reads P2000 and P2001 and classifies whether adding the extra channel class gx
contracts or expands optical-defect norms; exports a strict-lane diagnostic only.
"""
from __future__ import annotations
import json, platform
from pathlib import Path
from typing import Any
import numpy as np
import scipy.linalg as la

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2002_s952_strict_cutkosky_missing_channel_vs_structural_classifier_witness.json"
TS = "2026-05-18T00:00:00+00:00"


def load(name: str) -> dict[str, Any]:
    p = GEN / name
    if not p.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(p.read_text(encoding="utf-8"))


def delta_vector(packet: dict[str, Any]) -> np.ndarray:
    rows = packet.get("grid_table", [])
    return np.array([float(r.get("Delta_opt", 0.0)) for r in rows], dtype=float)


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2000 = load("p2000_s950_strict_cutkosky_loop_kernel_two_channel_witness.json")
    p2001 = load("p2001_s951_strict_cutkosky_full_three_loop_kernel_plus_extra_channel_witness.json")

    d0 = delta_vector(p2000)
    d1 = delta_vector(p2001)

    l2_2000 = float(la.norm(d0, 2)) if d0.size else float("nan")
    l2_2001 = float(la.norm(d1, 2)) if d1.size else float("nan")
    max_2000 = float(np.max(np.abs(d0))) if d0.size else float("nan")
    max_2001 = float(np.max(np.abs(d1))) if d1.size else float("nan")

    l2_ratio = l2_2001 / l2_2000 if l2_2000 > 0 else float("inf")
    max_ratio = max_2001 / max_2000 if max_2000 > 0 else float("inf")

    # Conservative classifier:
    # ratio < 0.95 => contraction (missing-channel pressure supported)
    # ratio > 1.05 => expansion (structural-obstruction pressure supported)
    # else mixed/inconclusive.
    if l2_ratio < 0.95 and max_ratio < 0.95:
        cls = "MISSING_CHANNEL_PRESSURE_SUPPORTED"
    elif l2_ratio > 1.05 and max_ratio > 1.05:
        cls = "STRUCTURAL_OBSTRUCTION_PRESSURE_SUPPORTED"
    else:
        cls = "MIXED_OR_INCONCLUSIVE"

    gate = {
        "p2000_present": p2000.get("result_kind") == "PASS_LOOP_KERNEL_TWO_CHANNEL_WITNESS",
        "p2001_present": p2001.get("result_kind") == "PASS_FULL_THREE_LOOP_PLUS_EXTRA_CHANNEL_WITNESS",
        "shared_grid_length": d0.size == d1.size and d0.size > 0,
        "l2_ratio_finite": bool(np.isfinite(l2_ratio)),
        "max_ratio_finite": bool(np.isfinite(max_ratio)),
    }

    out = {
        "ledger_id": "P2002_S952_STRICT_CUTKOSKY_MISSING_CHANNEL_VS_STRUCTURAL_CLASSIFIER_WITNESS",
        "packet_id": "P2002",
        "stage_id": "S952",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TS,
        "route": "strict_only",
        "inputs": {
            "p2000": p2000.get("ledger_id"),
            "p2001": p2001.get("ledger_id"),
        },
        "delta_norm_comparison": {
            "l2_p2000": l2_2000,
            "l2_p2001": l2_2001,
            "l2_ratio_p2001_over_p2000": l2_ratio,
            "max_abs_p2000": max_2000,
            "max_abs_p2001": max_2001,
            "max_ratio_p2001_over_p2000": max_ratio,
        },
        "classifier": cls,
        "gatekeeper_checks": gate,
        "result_kind": "PASS_DELTA_NORM_CLASSIFIER_WITNESS" if all(gate.values()) else "OPEN_OBSTRUCTION_WITH_TRACE",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "false_pass_guard": "Classifier is comparative diagnostic only; it is not a proof of full unitarity closure.",
        "next_honest_step": "Replace gx by explicit loop-derived amplitude and rerun classifier with uncertainty bands to test robustness of obstruction class.",
        "lay_explanation": "Porównaliśmy dwa kolejne modele: jeśli po dodaniu kanału błąd maleje, brakowało kanału; jeśli rośnie, bardziej prawdopodobna jest obstrukcja strukturalna.",
        "environment": {"python": platform.python_version(), "numpy": np.__version__, "scipy": __import__('scipy').__version__},
    }

    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P2002] wrote witness: {OUT}")


if __name__ == "__main__":
    main()
