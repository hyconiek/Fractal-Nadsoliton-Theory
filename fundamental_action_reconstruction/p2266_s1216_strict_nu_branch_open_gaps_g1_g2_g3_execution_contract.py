#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2265 = GEN / "p2265_s1215_strict_nu_branch_task3_closure_gap_matrix_probe.json"
OUT = GEN / "p2266_s1216_strict_nu_branch_open_gaps_g1_g2_g3_execution_contract.json"
MD = GEN / "p2266_s1216_strict_nu_branch_open_gaps_g1_g2_g3_execution_contract.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2265 = load(IN_2265)
    probe = p2265.get("strict_nu_branch_task3_closure_gap_matrix_probe", {}) or {}
    gaps = probe.get("task3_gap_matrix", []) or []

    execution_items = {
        "G1_reduction_certainty": {
            "objective": "Tighten reduction certainty beyond overlap-only evidence.",
            "deliverable": "seed_draw_sensitivity_table + one-sided positive-lift certificate",
            "acceptance_rule": "bootstrap_ci95_low > 0 for every declared seed/draw slice",
            "status_after_contract": "PARTIAL_READY_FOR_NUMERIC_REPLAY",
        },
        "G2_nonlinear_trajectory_realism": {
            "objective": "Move from proxy consistency to nonlinear trajectory residual validation.",
            "deliverable": "nonlinear trajectory family residual ledger (bounded perturbation classes)",
            "acceptance_rule": "max_abs_residual <= declared epsilon on each perturbation family",
            "status_after_contract": "PARTIAL_READY_FOR_RESIDUAL_REPLAY",
        },
        "G3_operational_policy_rule": {
            "objective": "Export risk->controller policy rule with explicit no-false-pass safety bounds.",
            "deliverable": "risk-calibrated parameter map + safety certificate",
            "acceptance_rule": "every exported policy row includes bounded-risk domain + positivity_guard=true",
            "status_after_contract": "PARTIAL_READY_FOR_POLICY_REPLAY",
        },
    }

    contracted = []
    for gap in gaps:
        gid = gap.get("id", "")
        contract = execution_items.get(gid, {})
        contracted.append({
            "id": gid,
            "prior_status": gap.get("status", "OPEN"),
            "description": gap.get("description", ""),
            "execution_contract": contract,
            "post_contract_status": contract.get("status_after_contract", "OPEN"),
        })

    closure_score = 0.0

    payload = {
        "schema_version": "p2266_s1216_v1",
        "packet_id": "P2266",
        "stage_id": "S1216",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_OPEN_GAPS_G1_G2_G3_EXECUTION_CONTRACT",
        "strict_nu_branch_open_gaps_execution_contract": {
            "contract_id": "STRICT_NU_BRANCH_OPEN_GAPS_CONTRACT_V1",
            "source_packets": [str(IN_2265.relative_to(ROOT))],
            "gap_contract_rows": contracted,
            "closure_score": closure_score,
            "closure_score_interpretation": "Execution contracts exported for G1/G2/G3, but no theorem-grade gap discharge claimed yet.",
            "theorem_scope_limit": "execution-contract freeze only; not strict-core selector closure and not ToE closure",
        },
        "recommended_next_honest_step": {
            "id": "P2267_candidate",
            "goal": "run G1 contract replay and export one-sided positivity lift certificate table",
        },
        "gatekeeper_checks": {
            "all_three_gaps_contracted": len(contracted) == 3,
            "g1_contracted": any(r.get("id") == "G1_reduction_certainty" for r in contracted),
            "g2_contracted": any(r.get("id") == "G2_nonlinear_trajectory_realism" for r in contracted),
            "g3_contracted": any(r.get("id") == "G3_operational_policy_rule" for r in contracted),
            "closure_score_bounded": 0.0 <= closure_score <= 1.0,
            "no_bridge_theorem_claimed": True,
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join([
            "# P2266 S1216: Open gaps G1/G2/G3 execution contract",
            "",
            "- Action taken: exported explicit execution contracts for G1, G2, G3.",
            f"- closure score: `{closure_score:.12e}`",
            "- status: OPEN (contracts ready, no theorem-grade discharge claimed).",
            "",
            "No selector closure / ToE closure claim.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
