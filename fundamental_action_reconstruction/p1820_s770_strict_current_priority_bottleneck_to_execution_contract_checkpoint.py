#!/usr/bin/env python3
"""P1820 S770 strict current-priority bottleneck to execution contract checkpoint."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"


def load(name: str) -> dict:
    return json.loads((GEN / name).read_text(encoding="utf-8"))


def rank(status: str) -> int:
    if status.startswith("PASS_ZERO"):
        return 0
    if status.startswith("LOCK_CONSISTENT"):
        return 1
    if status.startswith("OPEN"):
        return 3
    if status.startswith("KEEP_OPEN"):
        return 3
    if status.startswith("LOCKED"):
        return 2
    return 5




def is_active_blocker(status: str) -> bool:
    return status.startswith("OPEN") or status.startswith("KEEP_OPEN")

def main() -> None:
    p1764 = load("p1764_s714_strict_covariant_nonproxy_ea_eh_explicit_export_checkpoint.json")
    p1765 = load("p1765_s715_strict_nonproxy_metric_elg_explicit_export_checkpoint.json")
    p1762 = load("p1762_s712_strict_boundary_control_contract_finalization_checkpoint.json")
    p1763 = load("p1763_s713_strict_nonproxy_h1_4d_first_execution_attempt_checkpoint.json")
    p1766 = load("p1766_s716_strict_forward_reverse_state_vector_update_with_bianchi_ward_gate_checkpoint.json")
    p1819 = load("p1819_s769_strict_pre_tg_verdict_statevector_lock_binding_checkpoint.json")
    p1812 = load("p1812_s762_strict_canonical_gate_status_source_checkpoint.json")

    items = {
        "explicit_covariant_nonproxy_E_A_mu": p1764.get("status", "OPEN_UNKNOWN"),
        "explicit_covariant_nonproxy_E_H": p1764.get("status", "OPEN_UNKNOWN"),
        "metric_EL_g_export": p1765.get("status", "OPEN_UNKNOWN"),
        "boundary_term_control": p1762.get("status", "OPEN_UNKNOWN"),
        "H1_4D_weak_form_readiness": p1763.get("status", "OPEN_UNKNOWN"),
        "Bianchi_Ward_consistency": p1766.get("status", "OPEN_UNKNOWN"),
        "pre_TG_lock": p1819.get("status", "OPEN_UNKNOWN"),
        "TG1_BW": p1812.get("canonical_gate_status", {}).get("TG1_BW", "OPEN_UNKNOWN"),
        "TG2_BRST": p1812.get("canonical_gate_status", {}).get("TG2_BRST", "OPEN_UNKNOWN"),
        "TG3_CUT": p1812.get("canonical_gate_status", {}).get("TG3_CUT", "OPEN_UNKNOWN"),
    }

    ranked = sorted(
        [{"item": k, "status": v, "severity": rank(v)} for k, v in items.items()],
        key=lambda x: x["severity"],
        reverse=True,
    )

    active_ranked = [x for x in ranked if is_active_blocker(x["status"])]
    ranked_top = active_ranked[:3]

    out = {
        "packet_id": "P1820",
        "stage_id": "S770",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "current_priority_snapshot": items,
        "bottleneck_ranking": ranked,
        "top_bottlenecks": ranked_top,
        "top_bottlenecks_policy": "ACTIVE_OPEN_ONLY",
        "execution_contract": {
            "hard_precondition": "P1819 must be LOCK_CONSISTENT",
            "active_bottleneck": "TG1_BW unified nonproxy residual witness",
            "required_inputs": [
                "explicit covariant nonproxy E_A_mu export",
                "explicit covariant nonproxy E_H export",
                "metric EL_g export",
                "boundary-term control clauses",
                "H1 4D weak-form residual pack",
                "Bianchi/Ward divergence trace"
            ],
            "required_output": [
                "binary TG1 verdict",
                "residual witness trace",
                "state-vector sync without PASS promotion"
            ]
        },
        "technical_progress": "Priority bottlenecks normalized with dependency-aware severity and active-open-only top list for TG1 execution focus.",
        "proven": "Pre-TG consistency lock is satisfied at lock level, but theorem gates remain open.",
        "open": "No global nonproxy residual witness exported yet for TG1/TG2/TG3 closure.",
        "false_pass_risk": "Confusing lock consistency with theorem closure would overstate progress.",
        "next_honest_step": "Execute unified TG1 nonproxy residual run with full witness trace and replay P1810/P1811/P1812.",
        "lay_explanation": "Mamy uporządkowaną listę braków: dane są spójne, ale nadal brakuje kluczowego testu fizycznego TG1.",
    }

    path = GEN / "p1820_s770_strict_current_priority_bottleneck_to_execution_contract_checkpoint.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
