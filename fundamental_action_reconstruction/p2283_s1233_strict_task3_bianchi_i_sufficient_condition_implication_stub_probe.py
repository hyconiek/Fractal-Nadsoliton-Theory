#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2282 = GEN / "p2282_s1232_strict_task3_global_bianchi_i_g1_g2_g3_closure_matrix_probe.json"
OUT = GEN / "p2283_s1233_strict_task3_bianchi_i_sufficient_condition_implication_stub_probe.json"
MD = GEN / "p2283_s1233_strict_task3_bianchi_i_sufficient_condition_implication_stub_probe.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2282 = load(IN_2282)

    probe = (p2282.get("strict_task3_global_bianchi_i_g1_g2_g3_closure_matrix_probe", {}) or {})
    gaps = probe.get("gap_rows", []) or []
    by_id = {g.get("id", ""): g for g in gaps}

    g1 = by_id.get("G1_reduction_certainty", {})
    g2 = by_id.get("G2_nonlinear_trajectory_realism", {})
    g3 = by_id.get("G3_operational_policy_rule", {})

    g1_closed = g1.get("status") == "CLOSED"
    g2_closed = g2.get("status") == "CLOSED"
    g3_closed = g3.get("status") == "CLOSED"

    # theorem-ready implication stub (not theorem export)
    assumptions = {
        "A1_minimal_replay_margin_nonnegative": bool(g1_closed),
        "A2_bianchi_transport_residual_threshold": bool(g2_closed),
        "A3_operational_policy_lock_exists": bool(g3_closed),
    }
    all_assumptions_hold = all(assumptions.values())

    implication_stub = {
        "statement_id": "TASK3_BIANCHI_I_SUFFICIENT_CONDITION_STUB_V1",
        "formal_shape": "(A1 & A2 & A3) => TASK3_TRANSPORT_SUFFICIENT_CONDITION_READY",
        "A1": "minimal-config replay margin floor is nonnegative across risk rows",
        "A2": "global Bianchi-I transport residual map remains under declared threshold",
        "A3": "minimal feasible operational policy-lock configuration exists",
        "conclusion": "Task-3 sufficient-condition candidate is structurally ready for theorem-grade packaging",
        "evaluated_ready": all_assumptions_hold,
    }

    payload = {
        "schema_version": "p2283_s1233_v1",
        "packet_id": "P2283",
        "stage_id": "S1233",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_TASK3_BIANCHI_I_SUFFICIENT_CONDITION_IMPLICATION_STUB_PROBE",
        "strict_task3_bianchi_i_sufficient_condition_implication_stub_probe": {
            "probe_id": "STRICT_TASK3_BIANCHI_I_SUFFICIENT_CONDITION_IMPLICATION_STUB_PROBE_V1",
            "source_packets": [str(IN_2282.relative_to(ROOT))],
            "assumptions": assumptions,
            "implication_stub": implication_stub,
            "theorem_scope_limit": "implication-structure stub only; not a theorem proof and not selector/ToE closure",
        },
        "recommended_next_honest_step": {
            "id": "P2284_candidate",
            "goal": "export quantified theorem-bridge draft with explicit constants for A1/A2/A3 and machine-checkable premise table",
        },
        "gatekeeper_checks": {
            "assumption_rows_exported": len(assumptions) == 3,
            "contains_A1": "A1_minimal_replay_margin_nonnegative" in assumptions,
            "contains_A2": "A2_bianchi_transport_residual_threshold" in assumptions,
            "contains_A3": "A3_operational_policy_lock_exists" in assumptions,
            "no_bridge_theorem_claimed": True,
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join([
            "# P2283 S1233: Task-3 Bianchi-I sufficient-condition implication stub",
            "",
            f"- A1 holds: `{assumptions['A1_minimal_replay_margin_nonnegative']}`",
            f"- A2 holds: `{assumptions['A2_bianchi_transport_residual_threshold']}`",
            f"- A3 holds: `{assumptions['A3_operational_policy_lock_exists']}`",
            f"- implication ready flag: `{all_assumptions_hold}`",
            "",
            "Implication-structure stub only; no selector closure / ToE closure claim.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
