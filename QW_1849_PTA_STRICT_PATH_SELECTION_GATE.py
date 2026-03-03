#!/usr/bin/env python3
"""
QW-1849: Strict path selection gate for PTA criterion closure.

Integrates:
- QW-1842 execution readiness,
- QW-1844 strict gate,
- QW-1847 replication targets (V1 criterion),
- QW-1848 unit-of-analysis audit.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1849_pta_strict_path_selection_gate.json"
OUT_MD = ROOT / "RAPORT_QW1849_PTA_STRICT_PATH_SELECTION_GATE.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def clamp01(x: float) -> float:
    return max(0.0, min(1.0, x))


def main() -> None:
    d1842 = load("report_qw1842_joint_confirmatory_execution_gate.json")
    d1844 = load("report_qw1844_strict_joint_rigor_gate.json")
    d1847 = load("report_qw1847_pta_replication_target_proposal.json")
    d1848 = load("report_qw1848_pta_unit_of_analysis_audit.json")

    ops_ready = d1842.get("hard_gate") == "PASS"
    strict_ready_v1 = d1844.get("hard_gate") == "PASS"

    s1848 = d1848.get("summary", {})
    pair_prob = float(s1848.get("pair_level", {}).get("prob_pair_mean_gain_positive", 0.0))
    split_prob = float(s1848.get("split_level", {}).get("prob_rep_mean_gain_positive", 0.0))
    compression_gap = float(s1848.get("compression_gap", 0.0))
    pval_prob_0p9 = float(s1848.get("inference", {}).get("pvalue_h0_prob_le_0p90", 1.0))

    unit_mismatch = d1848.get("verdict") == "PTA_UNIT_MISMATCH_REQUIRES_CRITERION_REDESIGN"
    severe_gap = compression_gap >= 0.10
    threshold_0p9_not_supported = (pair_prob < 0.90) or (pval_prob_0p9 > 0.05)

    n_target = d1847.get("targets", {}).get("recommended_target_n")
    n_add = d1847.get("targets", {}).get("additional_needed")

    if unit_mismatch and threshold_0p9_not_supported:
        best_path = "PATH_B_VERSIONED_CRITERION_REPARAM_WITH_EXTERNAL_CONFIRMATORY"
        rationale = (
            "Current PTA V1 probability criterion (0.9) is not supported under pair-level analysis; "
            "split-level positivity appears inflated by averaging compression."
        )
    elif (not strict_ready_v1) and (n_add is not None) and (int(n_add) > 0):
        best_path = "PATH_A_KEEP_V1_AND_EXPAND_REPLICATIONS"
        rationale = "Strict gate misses only due to underpower; increase independent replications under frozen V1 criterion."
    else:
        best_path = "PATH_UNCHANGED"
        rationale = "No path switch required from current strict evaluation."

    # weighted strictness score (high means methodologically consistent)
    score = clamp01(
        0.25 * float(ops_ready)
        + 0.20 * float(strict_ready_v1)
        + 0.25 * float(not unit_mismatch)
        + 0.15 * float(not severe_gap)
        + 0.15 * float(not threshold_0p9_not_supported)
    )

    if best_path.startswith("PATH_B"):
        hard_gate = "PARTIAL"
        readiness = "REPARAM_PROTOCOL_REQUIRED_BEFORE_STRICT_CONFIRMATORY"
    elif best_path.startswith("PATH_A"):
        hard_gate = "PARTIAL"
        readiness = "DATA_EXPANSION_REQUIRED_FOR_V1"
    else:
        hard_gate = "PASS" if strict_ready_v1 else "PARTIAL"
        readiness = "STRICT_READY" if strict_ready_v1 else "STRICT_PARTIAL"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "qw1842_hard_gate": d1842.get("hard_gate"),
            "qw1844_hard_gate": d1844.get("hard_gate"),
            "qw1847_recommended_target_n": n_target,
            "qw1847_additional_needed": n_add,
            "qw1848_verdict": d1848.get("verdict"),
        },
        "diagnostics": {
            "split_level_prob_positive": split_prob,
            "pair_level_prob_positive": pair_prob,
            "compression_gap": compression_gap,
            "pvalue_h0_prob_le_0p90": pval_prob_0p9,
        },
        "flags": {
            "ops_ready": bool(ops_ready),
            "strict_ready_v1": bool(strict_ready_v1),
            "unit_mismatch": bool(unit_mismatch),
            "severe_gap": bool(severe_gap),
            "threshold_0p9_not_supported": bool(threshold_0p9_not_supported),
        },
        "decision": {
            "best_path": best_path,
            "rationale": rationale,
            "alternative_path_a_keep_v1": {
                "target_n": n_target,
                "additional_needed": n_add,
            },
            "alternative_path_b_reparam": {
                "requires_new_prereg_version": True,
                "requires_external_confirmatory_dataset": True,
            },
        },
        "strict_score": score,
        "hard_gate": hard_gate,
        "readiness": readiness,
        "verdict": "PTA_STRICT_PATH_SELECTION_COMPLETE",
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1849: PTA STRICT PATH SELECTION GATE",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Strict score: {score:.3f}",
        f"- Hard gate: **{hard_gate}**",
        f"- Readiness: **{readiness}**",
        f"- Best path: **{best_path}**",
        "",
        "## Diagnostics",
        f"- split-level prob positive: {split_prob:.3f}",
        f"- pair-level prob positive: {pair_prob:.3f}",
        f"- compression gap: {compression_gap:.3f}",
        f"- p-value H0(prob<=0.90): {pval_prob_0p9:.6f}",
        "",
        "## Decision Rationale",
        f"- {rationale}",
        "",
        "## Alternative A (keep V1)",
        f"- target_n: {n_target}",
        f"- additional_needed: {n_add}",
        "",
        "## Alternative B (reparam)",
        "- requires new prereg version: True",
        "- requires external confirmatory dataset: True",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1849] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1849] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
