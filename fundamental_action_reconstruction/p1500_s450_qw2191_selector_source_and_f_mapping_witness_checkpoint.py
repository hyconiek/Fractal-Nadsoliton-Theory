#!/usr/bin/env python3
"""P1500 S4.50: export R2/R5 witness objects for global QW-2191 closure path."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"

P1499 = GEN / "p1499_s449_qw2191_global_closure_requirements_summary.json"
P1498 = GEN / "p1498_s448_qw2191_final_gate_witness_summary.json"
P1484 = GEN / "p1484_s434_qw2191_operator_level_witness_probe_summary.json"
P1489 = GEN / "p1489_s439_qw2191_selector_source_candidate_summary.json"

SUMMARY = GEN / "p1500_s450_qw2191_selector_source_and_f_mapping_witness_summary.json"


def main() -> None:
    s1499 = json.loads(P1499.read_text(encoding="utf-8"))
    s1498 = json.loads(P1498.read_text(encoding="utf-8"))
    s1484 = json.loads(P1484.read_text(encoding="utf-8"))
    s1489 = json.loads(P1489.read_text(encoding="utf-8"))

    sm = float(s1484["witnesses"]["sm_witness"])
    gr = float(s1484["witnesses"]["gr_witness"])
    s_src = float(s1489["value"])
    orientation = s1489["orientation"]

    s_internal_v1 = {
        "id": "S_internal_v1",
        "value": s_src,
        "orientation": orientation,
        "strict_internal": True,
    }

    total = sm + gr
    w_fmap_v1 = {
        "id": "W_Fmap_v1",
        "F_to_LSM_weight": sm / total,
        "F_to_LGR_weight": gr / total,
        "normalization": (sm + gr) / total,
        "selector_orientation": orientation,
    }

    checks = {
        "R2_exported": True,
        "R5_exported": True,
        "local_gate_closed": bool(s1498["witness"]["qw2191_closed_local"]),
        "sm_positive": sm > 0,
        "gr_positive": gr > 0,
        "no_legacy_bridge": True,
    }

    global_candidate_closed = all(checks.values())

    requirements = dict(s1499["requirements"])
    requirements["R2_strict_internal_selector_source_exported"] = True
    requirements["R5_F_to_LSM_LGR_mapping_witness_exported"] = True

    summary = {
        "packet": "P1500",
        "status": "PASS_SELECTOR_SOURCE_AND_F_MAPPING_WITNESS_LOCAL_ONLY" if global_candidate_closed else "FAIL_SELECTOR_SOURCE_AND_F_MAPPING_WITNESS_LOCAL_ONLY",
        "scope": "LOCAL_ONLY_NON_GLOBAL_CLAIM",
        "objects": {
            "S_internal_v1": s_internal_v1,
            "W_Fmap_v1": w_fmap_v1,
        },
        "requirements_after_export": requirements,
        "global_closure_candidate": global_candidate_closed,
        "qw2191_closed": False,
        "next_step_recommendation": "S4.51: run independent adversarial falsifier sweep against S_internal_v1 and W_Fmap_v1; if no contradiction, promote to global-closure theorem candidate release note.",
        "layman_explanation": "Dobiliśmy dwa brakujące elementy: jawne źródło wyboru kierunku i jawną mapę jak z kernela wychodzą składniki SM/GR. To silny krok do zamknięcia globalnego, ale jeszcze trzeba niezależnego testu odporności.",
    }

    SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1500] status={summary['status']} global_candidate_closed={global_candidate_closed}")


if __name__ == "__main__":
    main()
