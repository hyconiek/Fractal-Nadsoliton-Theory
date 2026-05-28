#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2281 = GEN / "p2281_s1231_strict_nu_branch_group_policy_minimal_config_fresh_replay_validation_probe.json"
IN_2282 = GEN / "p2282_s1232_strict_task3_global_bianchi_i_g1_g2_g3_closure_matrix_probe.json"
IN_2283 = GEN / "p2283_s1233_strict_task3_bianchi_i_sufficient_condition_implication_stub_probe.json"
OUT = GEN / "p2284_s1234_strict_task3_bianchi_i_quantified_premise_table_and_theorem_bridge_draft_probe.json"
MD = GEN / "p2284_s1234_strict_task3_bianchi_i_quantified_premise_table_and_theorem_bridge_draft_probe.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2281 = load(IN_2281)
    p2282 = load(IN_2282)
    p2283 = load(IN_2283)

    rows2281 = (p2281.get("strict_nu_branch_group_policy_minimal_config_fresh_replay_validation_probe", {}) or {}).get("rows", []) or []
    lock = (p2281.get("strict_nu_branch_group_policy_minimal_config_fresh_replay_validation_probe", {}) or {}).get("locked_minimal_config", {}) or {}
    closure_rows = (p2282.get("strict_task3_global_bianchi_i_g1_g2_g3_closure_matrix_probe", {}) or {}).get("gap_rows", []) or []
    assumptions = (p2283.get("strict_task3_bianchi_i_sufficient_condition_implication_stub_probe", {}) or {}).get("assumptions", {}) or {}

    g2_row = next((r for r in closure_rows if r.get("id") == "G2_nonlinear_trajectory_realism"), {})
    g3_row = next((r for r in closure_rows if r.get("id") == "G3_operational_policy_rule"), {})
    g2_metric = float(g2_row.get("metric", 1.0) or 1.0)

    a1_min_margin = min((float(r.get("margin_to_target", -1.0)) for r in rows2281), default=-1.0)
    a2_max_residual = g2_metric
    a3_lock_cost = float(g3_row.get("metric", -1.0) or -1.0)
    a3_policy_lock_closed = g3_row.get("status") == "CLOSED"

    premise_table = [
        {
            "premise_id": "A1",
            "symbol": "min_margin_to_target",
            "value": a1_min_margin,
            "required_relation": ">= 0",
            "satisfied": a1_min_margin >= -1e-12,
        },
        {
            "premise_id": "A2",
            "symbol": "max_transport_residual_l1",
            "value": a2_max_residual,
            "required_relation": "<= 5e-5",
            "satisfied": a2_max_residual <= 5e-5,
        },
        {
            "premise_id": "A3",
            "symbol": "policy_lock_cost_proxy",
            "value": a3_lock_cost,
            "required_relation": "G3 status == CLOSED with real P2280 feasible policy lock",
            "satisfied": a3_policy_lock_closed,
        },
    ]

    all_quantified_satisfied = all(p["satisfied"] for p in premise_table)
    stub_ready = bool(assumptions.get("A1_minimal_replay_margin_nonnegative", False)) and bool(assumptions.get("A2_bianchi_transport_residual_threshold", False)) and bool(assumptions.get("A3_operational_policy_lock_exists", False))

    theorem_bridge_draft = {
        "draft_id": "TASK3_BIANCHI_I_THEOREM_BRIDGE_DRAFT_V1",
        "premise_layer": "quantified machine-checkable A1/A2/A3 table",
        "logical_form": "(A1_quant & A2_quant & A3_quant) => TASK3_TRANSPORT_SUFFICIENT_CONDITION_CANDIDATE",
        "compatibility_with_stub": stub_ready,
        "quantified_ready": all_quantified_satisfied,
        "note": "draft object only; no theorem proof export",
    }

    payload = {
        "schema_version": "p2284_s1234_v1",
        "packet_id": "P2284",
        "stage_id": "S1234",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_TASK3_BIANCHI_I_QUANTIFIED_PREMISE_TABLE_AND_THEOREM_BRIDGE_DRAFT_PROBE",
        "strict_task3_bianchi_i_quantified_premise_table_and_theorem_bridge_draft_probe": {
            "probe_id": "STRICT_TASK3_BIANCHI_I_QUANTIFIED_PREMISE_TABLE_AND_THEOREM_BRIDGE_DRAFT_PROBE_V1",
            "source_packets": [str(IN_2281.relative_to(ROOT)), str(IN_2282.relative_to(ROOT)), str(IN_2283.relative_to(ROOT))],
            "premise_table": premise_table,
            "theorem_bridge_draft": theorem_bridge_draft,
            "theorem_scope_limit": "quantified draft packaging only; not theorem proof and not selector/ToE closure",
        },
        "recommended_next_honest_step": {
            "id": "P2285_candidate",
            "goal": "export machine-checkable verifier that re-evaluates A1/A2/A3 from source packets and produces deterministic PASS/OPEN certificate",
        },
        "gatekeeper_checks": {
            "premise_table_exported": len(premise_table) == 3,
            "all_premise_ids_present": {"A1", "A2", "A3"} == {p["premise_id"] for p in premise_table},
            "draft_exported": True,
            "no_bridge_theorem_claimed": True,
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text("\n".join([
        "# P2284 S1234: quantified premise table + theorem bridge draft",
        "",
        f"- A1 min margin: `{a1_min_margin:.12e}`",
        f"- A2 max residual: `{a2_max_residual:.12e}`",
        f"- A3 locked trials budget: `{a3_lock_cost:.12e}`",
        f"- quantified premises satisfied: `{all_quantified_satisfied}`",
        f"- stub compatibility: `{stub_ready}`",
        "",
        "Draft packaging only; no selector closure / ToE closure claim.",
    ]) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
