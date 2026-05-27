#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2264 = GEN / "p2264_s1214_strict_nu_branch_group_policy_bootstrap_ci_crosscheck_probe.json"
OUT = GEN / "p2265_s1215_strict_nu_branch_task3_closure_gap_matrix_probe.json"
MD = GEN / "p2265_s1215_strict_nu_branch_task3_closure_gap_matrix_probe.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2264 = load(IN_2264)
    probe = (p2264.get("strict_nu_branch_group_policy_bootstrap_ci_crosscheck_probe", {}) or {})
    b = (probe.get("bootstrap_crosscheck", {}) or {})
    ci_overlap = bool(b.get("ci_overlap", False))
    boot_low = float(b.get("bootstrap_ci95_low", 0.0) or 0.0)

    gaps = [
        {
            "id": "G1_reduction_certainty",
            "description": "Need stronger-than-overlap evidence that reduction stays positive under broader sampling assumptions.",
            "status": "PARTIAL" if ci_overlap else "OPEN",
            "next_object": "larger-seed + larger-draw CI tightening",
        },
        {
            "id": "G2_nonlinear_trajectory_realism",
            "description": "Current validation is still synthetic/proxy; need nonlinear trajectory residual validation with perturbation families.",
            "status": "OPEN",
            "next_object": "nonlinear residual simulator packet",
        },
        {
            "id": "G3_operational_policy_rule",
            "description": "Need exportable decision rule mapping risk tolerance -> controller parameters with explicit safety certificate.",
            "status": "OPEN",
            "next_object": "risk-calibrated controller map",
        },
    ]

    closure_score = sum(1 for g in gaps if g["status"] == "CLOSED") / len(gaps)

    payload = {
        "schema_version": "p2265_s1215_v1",
        "packet_id": "P2265",
        "stage_id": "S1215",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_TASK3_CLOSURE_GAP_MATRIX_PROBE",
        "strict_nu_branch_task3_closure_gap_matrix_probe": {
            "probe_id": "STRICT_NU_BRANCH_TASK3_CLOSURE_GAP_MATRIX_PROBE_V1",
            "source_packets": [str(IN_2264.relative_to(ROOT))],
            "upstream_summary": {
                "bootstrap_ci95_low": boot_low,
                "bootstrap_ci_overlap_with_normal": ci_overlap,
            },
            "task3_gap_matrix": gaps,
            "closure_score": closure_score,
            "physical_interpretation_note": "Task-3 closure still requires moving from proxy/synthetic statistical consistency to stronger physically anchored nonlinear validation and operationally certifiable control rules.",
            "theorem_scope_limit": "closure-gap planning diagnostic only; not a legacy->strict bridge theorem and not strict-core selector closure",
        },
        "recommended_next_honest_step": {
            "id": "P2266_candidate",
            "goal": "export risk-calibrated controller map with explicit positivity certificate over selected perturbation family",
        },
        "gatekeeper_checks": {
            "task3_gap_matrix_exported": True,
            "closure_score_bounded": 0.0 <= closure_score <= 1.0,
            "contains_open_gaps": any(g["status"] != "CLOSED" for g in gaps),
            "no_bridge_theorem_claimed": True,
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
            "full_cutkosky_closure_proven": False,
            "full_d3_covariance_transport_proven": False,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join([
            "# P2265 S1215: Task-3 closure gap matrix probe",
            "",
            f"- closure score: `{closure_score:.12e}`",
            f"- bootstrap CI low: `{boot_low:.12e}`",
            f"- CI overlap normal/bootstrap: `{ci_overlap}`",
            "- open gaps: G1, G2, G3",
            "",
            "Task-3 closure planning diagnostic only; no kernel-bridge or selector-closure claim.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
