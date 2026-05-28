#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2203 = GEN / "p2203_s1153_strict_frw_bianchi_transport_residual_map_under_shared_majorant.json"
IN_2280 = GEN / "p2280_s1230_strict_nu_branch_group_policy_minimal_lock_criterion_grid_probe.json"
IN_2281 = GEN / "p2281_s1231_strict_nu_branch_group_policy_minimal_config_fresh_replay_validation_probe.json"
OUT = GEN / "p2282_s1232_strict_task3_global_bianchi_i_g1_g2_g3_closure_matrix_probe.json"
MD = GEN / "p2282_s1232_strict_task3_global_bianchi_i_g1_g2_g3_closure_matrix_probe.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2203 = load(IN_2203)
    p2280 = load(IN_2280)
    p2281 = load(IN_2281)

    # G1: reduction certainty from fresh minimal-config replay margins.
    # Require both row-level nonnegative margins and the upstream all-rows summary,
    # so an empty/default replay cannot silently close G1.
    replay_probe = p2281.get("strict_nu_branch_group_policy_minimal_config_fresh_replay_validation_probe", {}) or {}
    rows2281 = replay_probe.get("rows", []) or []
    replay_summary = replay_probe.get("global_summary", {}) or {}
    g1_metric = min((float(r.get("margin_to_target", -1.0)) for r in rows2281), default=-1.0)
    g1_pass = (
        len(rows2281) > 0
        and bool(replay_summary.get("all_rows_meet_target", False))
        and all(float(r.get("margin_to_target", -1.0)) >= -1e-12 for r in rows2281)
    )

    # G2: nonlinear trajectory realism anchored to Bianchi-I transport residual map amplitude control.
    # P2203 exports residual_map_rows; keep residual_map as backward-compatible fallback only.
    transport_probe = p2203.get("strict_frw_bianchi_transport_residual_map_under_shared_majorant", {}) or {}
    map2203 = transport_probe.get("residual_map_rows", []) or transport_probe.get("residual_map", []) or []
    max_transport_residual = max((float(r.get("transport_residual_l1", 1.0)) for r in map2203), default=1.0)
    g2_threshold = 5e-5
    g2_pass = max_transport_residual <= g2_threshold

    # G3: operational policy rule requires an actual feasible minimal lock config from P2280.
    # Do not count P2281's deterministic fallback replay parameters as a real policy lock.
    lock_probe = p2280.get("strict_nu_branch_group_policy_minimal_lock_criterion_grid_probe", {}) or {}
    minimal = lock_probe.get("minimal_feasible_config", {}) or {}
    feasible_count = int(lock_probe.get("feasible_count", 0) or 0)
    g3_pass = feasible_count > 0 and bool(minimal)
    g3_metric = float(minimal.get("cost_proxy", -1.0) or -1.0)

    gaps = [
        {
            "id": "G1_reduction_certainty",
            "status": "CLOSED" if g1_pass else "OPEN",
            "metric": g1_metric,
            "criterion": "min margin_to_target from P2281 >= 0",
        },
        {
            "id": "G2_nonlinear_trajectory_realism",
            "status": "CLOSED" if g2_pass else "OPEN",
            "metric": max_transport_residual,
            "criterion": f"max Bianchi-I transport_residual_l1 <= {g2_threshold}",
        },
        {
            "id": "G3_operational_policy_rule",
            "status": "CLOSED" if g3_pass else "OPEN",
            "metric": g3_metric,
            "criterion": "P2280 feasible_count > 0 and minimal feasible policy-lock config exists",
        },
    ]

    closure_score = sum(1 for g in gaps if g["status"] == "CLOSED") / len(gaps)

    payload = {
        "schema_version": "p2282_s1232_v1",
        "packet_id": "P2282",
        "stage_id": "S1232",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE" if closure_score < 1.0 else "PASS_ZERO",
        "result_kind": "PASS_STRICT_TASK3_GLOBAL_BIANCHI_I_G1_G2_G3_CLOSURE_MATRIX_PROBE",
        "strict_task3_global_bianchi_i_g1_g2_g3_closure_matrix_probe": {
            "probe_id": "STRICT_TASK3_GLOBAL_BIANCHI_I_G1_G2_G3_CLOSURE_MATRIX_PROBE_V1",
            "source_packets": [
                str(IN_2203.relative_to(ROOT)),
                str(IN_2280.relative_to(ROOT)),
                str(IN_2281.relative_to(ROOT)),
            ],
            "task3_scope": "global Bianchi-I transport calibration lane with G1/G2/G3 closure diagnostics",
            "gap_rows": gaps,
            "closure_score": closure_score,
            "theorem_scope_limit": "closure-matrix diagnostic only; not strict-core selector closure and not ToE closure",
        },
        "recommended_next_honest_step": {
            "id": "P2283_candidate",
            "goal": "export theorem-ready implication stub from Bianchi-I residual threshold + minimal-config replay margin to Task-3 transport closure sufficient condition",
        },
        "gatekeeper_checks": {
            "gap_rows_exported": len(gaps) == 3,
            "contains_G1": any(g["id"] == "G1_reduction_certainty" for g in gaps),
            "contains_G2": any(g["id"] == "G2_nonlinear_trajectory_realism" for g in gaps),
            "contains_G3": any(g["id"] == "G3_operational_policy_rule" for g in gaps),
            "closure_score_bounded": 0.0 <= closure_score <= 1.0,
            "no_bridge_theorem_claimed": True,
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE" if closure_score < 1.0 else "PASS_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join([
            "# P2282 S1232: Task-3 global Bianchi-I G1/G2/G3 closure matrix",
            "",
            f"- closure score: `{closure_score:.12e}`",
            f"- G1 closed: `{g1_pass}` (min margin: `{g1_metric:.12e}`)",
            f"- G2 closed: `{g2_pass}` (max residual L1: `{max_transport_residual:.12e}`)",
            f"- G3 closed: `{g3_pass}` (cost proxy: `{g3_metric:.12e}`)",
            "",
            "Closure-matrix diagnostic only; no selector closure / ToE closure claim.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
