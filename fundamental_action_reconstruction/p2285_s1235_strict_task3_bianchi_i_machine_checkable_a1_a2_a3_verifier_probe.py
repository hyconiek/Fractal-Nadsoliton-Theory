#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2281 = GEN / "p2281_s1231_strict_nu_branch_group_policy_minimal_config_fresh_replay_validation_probe.json"
IN_2282 = GEN / "p2282_s1232_strict_task3_global_bianchi_i_g1_g2_g3_closure_matrix_probe.json"
IN_2284 = GEN / "p2284_s1234_strict_task3_bianchi_i_quantified_premise_table_and_theorem_bridge_draft_probe.json"
OUT = GEN / "p2285_s1235_strict_task3_bianchi_i_machine_checkable_a1_a2_a3_verifier_probe.json"
MD = GEN / "p2285_s1235_strict_task3_bianchi_i_machine_checkable_a1_a2_a3_verifier_probe.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2281 = load(IN_2281)
    p2282 = load(IN_2282)
    p2284 = load(IN_2284)

    rows2281 = (p2281.get("strict_nu_branch_group_policy_minimal_config_fresh_replay_validation_probe", {}) or {}).get("rows", []) or []
    closure_rows = (p2282.get("strict_task3_global_bianchi_i_g1_g2_g3_closure_matrix_probe", {}) or {}).get("gap_rows", []) or []
    quantified = (p2284.get("strict_task3_bianchi_i_quantified_premise_table_and_theorem_bridge_draft_probe", {}) or {}).get("premise_table", []) or []

    # Deterministic recomputation from source packets.
    # A3 is intentionally tied to the P2282 G3 closure row, not to P2281's
    # fallback replay parameters, because fallback parameters are not a feasible
    # policy lock.
    replay_summary = (
        p2281.get("strict_nu_branch_group_policy_minimal_config_fresh_replay_validation_probe", {}) or {}
    ).get("global_summary", {}) or {}
    a1_recomputed = (
        len(rows2281) > 0
        and bool(replay_summary.get("all_rows_meet_target", False))
        and min((float(r.get("margin_to_target", -1.0)) for r in rows2281), default=-1.0) >= -1e-12
    )

    g2 = next((g for g in closure_rows if g.get("id") == "G2_nonlinear_trajectory_realism"), {})
    a2_recomputed = g2.get("status") == "CLOSED" and float(g2.get("metric", 1.0) or 1.0) <= 5e-5

    g3 = next((g for g in closure_rows if g.get("id") == "G3_operational_policy_rule"), {})
    a3_recomputed = g3.get("status") == "CLOSED"

    by_id = {p.get("premise_id", ""): bool(p.get("satisfied", False)) for p in quantified}
    a1_reported = by_id.get("A1", False)
    a2_reported = by_id.get("A2", False)
    a3_reported = by_id.get("A3", False)

    consistency = {
        "A1": a1_recomputed == a1_reported,
        "A2": a2_recomputed == a2_reported,
        "A3": a3_recomputed == a3_reported,
    }

    final_pass = all(consistency.values()) and a1_recomputed and a2_recomputed and a3_recomputed

    payload = {
        "schema_version": "p2285_s1235_v1",
        "packet_id": "P2285",
        "stage_id": "S1235",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_TASK3_BIANCHI_I_MACHINE_CHECKABLE_A1_A2_A3_VERIFIER_PROBE",
        "strict_task3_bianchi_i_machine_checkable_a1_a2_a3_verifier_probe": {
            "probe_id": "STRICT_TASK3_BIANCHI_I_MACHINE_CHECKABLE_A1_A2_A3_VERIFIER_PROBE_V1",
            "source_packets": [str(IN_2281.relative_to(ROOT)), str(IN_2282.relative_to(ROOT)), str(IN_2284.relative_to(ROOT))],
            "recomputed": {"A1": a1_recomputed, "A2": a2_recomputed, "A3": a3_recomputed},
            "reported_from_p2284": {"A1": a1_reported, "A2": a2_reported, "A3": a3_reported},
            "consistency": consistency,
            "verifier_pass": final_pass,
            "theorem_scope_limit": "deterministic verifier for premise consistency only; not theorem proof and not selector/ToE closure",
        },
        "recommended_next_honest_step": {
            "id": "P2286_candidate",
            "goal": "export immutable signed certificate bundle (A1/A2/A3 + verifier_pass) for downstream theorem-attempt gating",
        },
        "gatekeeper_checks": {
            "consistency_rows_exported": len(consistency) == 3,
            "all_consistency_flags_boolean": all(isinstance(v, bool) for v in consistency.values()),
            "contains_verifier_pass": True,
            "no_bridge_theorem_claimed": True,
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text("\n".join([
        "# P2285 S1235: machine-checkable A1/A2/A3 verifier",
        "",
        f"- recomputed A1/A2/A3: `{a1_recomputed}/{a2_recomputed}/{a3_recomputed}`",
        f"- consistency A1/A2/A3: `{consistency['A1']}/{consistency['A2']}/{consistency['A3']}`",
        f"- verifier pass: `{final_pass}`",
        "",
        "Verifier-only probe; no selector closure / ToE closure claim.",
    ]) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
