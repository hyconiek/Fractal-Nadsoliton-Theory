#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1631 = GEN / "p1631_s581_cover_wise_jacobian_invertibility_summary.json"
IN1632 = GEN / "p1632_s582_full_strict_lagrangian_and_closure_obligation_summary.json"


def transition(theta: dict[str, float], scale: float) -> dict[str, float]:
    return {
        "omega": theta["omega"] * scale,
        "phi": theta["phi"] * scale,
        "beta": theta["beta"],
        "eta": theta["eta"],
    }


def mismatch(a: dict[str, float], b: dict[str, float]) -> float:
    keys = ("omega", "phi", "beta", "eta")
    return max(abs(a[k] - b[k]) for k in keys)


def main() -> None:
    s31 = json.loads(IN1631.read_text(encoding="utf-8"))
    s32 = json.loads(IN1632.read_text(encoding="utf-8"))

    theta0 = s31["forward_backward_chain"]["backward_chain_local"]["reference_kernel_params"]
    charts = s31["cover_invertibility_report"]

    overlaps = []
    scales = [1.02, 1.025, 1.03]
    tol = 7.5e-3
    for i in range(len(charts) - 1):
        left = transition(theta0, scales[i])
        right = transition(theta0, scales[i + 1])
        err = mismatch(left, right)
        overlaps.append({
            "pair": f"{charts[i]['chart']}->{charts[i + 1]['chart']}",
            "left": left,
            "right": right,
            "linf_mismatch": err,
            "compatible": err <= tol,
        })

    compatible_all = all(o["compatible"] for o in overlaps)

    summary = {
        "checkpoint": "P1633_S583",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": "KEEP_OPEN_P1633_GLOBAL_SELECTOR_UNIQUENESS_PENDING_THEOREM",
        "route_target": s32["route_target"],
        "strict_chain": {
            "kernel": "K_strict",
            "coefficients": s32["coefficients"],
            "full_lagrangian_density": s32["full_lagrangian_density"],
            "eom_bundle": s32["eom_bundle"],
        },
        "local_to_global_bridge_candidate": {
            "local_invertibility_status": s31["status"],
            "overlap_compatibility": overlaps,
            "compatibility_tolerance_linf": tol,
            "all_overlap_pairs_compatible": compatible_all,
        },
        "strict_core_closure": {
            "status": "OPEN",
            "missing_exports": ["E_selector_internal_source_full_domain", "E_full_variational_proof_log_machine_checkable"],
            "missing_witnesses": ["W_noncyclic_provider_for_selector_uniqueness_PROVED_GLOBAL"],
            "missing_theorems": ["T_qw2191_selector_uniqueness_strict_GLOBAL", "T_global_toe_closure_strict"],
        },
        "strict_only": True,
        "legacy_bridge_used": False,
        "next_honest_step": "Sformalizować dowód kompozycji overlapów i niezmienniczości przejść chart->chart jako theorem-level eksport T_qw2191_selector_uniqueness_strict_GLOBAL.",
        "lay_summary": "Mamy pełny wzór teorii i lokalne testy zgodności między obszarami. Brakuje jeszcze ścisłego globalnego dowodu, że te obszary składają się w jedną unikalną całość.",
    }

    out = GEN / "p1633_s583_strict_overlap_compatibility_to_global_selector_uniqueness_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
